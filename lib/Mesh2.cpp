
#include "Mesh.h"
#include "Plot.h"
#include "TimeFunction.h"
#include "PropertyTypeEnum.h"
#include "ComputerTime.h"


extern Plot               plot;
extern List<TimeFunction> timeFunction;
extern ComputerTime       computerTime;



static const unsigned int pow2[8] = {1,2,4,8,16,32,64,128};


using namespace std;







void Mesh::plotBoun_hlp(double *x, double &dpt, int &symb)
{
  switch (symb)
  {
    case XDASH   : plot.dash(x,dpt,0); break;

    case YDASH   : plot.dash(x,dpt,1); break;

    case ZDASH   : plot.dash(x,dpt,2); break;

    case CIRCLE  : plot.circle(x,dpt); break;
		   
    case TRIANGLE: plot.triangle(x,dpt); break;
		   
    case SQUARE  : plot.square(x,dpt); break;
		   
    default      : prgWarning(1,"Mesh::plotBoun_hlp","ignoring j > 6 !");
  }
  return;
}










bool Mesh::isBoundaryNode2D(int nd, int &nxtNd)
{
  int ne = nodeElem[nd-1].n;

  if (ne == 0) return false;

  int e, j = -1, m = 0, m0, np, *ndEl = &(nodeElem[nd-1][0]);

  Vector<int>  plotSeq;
  Vector<bool> virgin;
  
  for (e=0; e<ne; e++) 
  { 
    if (((ElementGroup*)elem[ndEl[e]]->elemGrp)->closed) 
      { virgin.append(true);  if (j < 0) j = e; }
    else 
      { virgin.append(false); m++; }
  }
  if (m == ne) return false;
 
  e = j;
  elem[ndEl[e]]->givePlotSequence2D(plotSeq);
  np = plotSeq.n;
  plotSeq.append(plotSeq[0]);
  j = 0; 
  while (plotSeq[j] != nd && j < np) j++;
  if (j == np) return false; 
  m = plotSeq[++j];
  plotSeq.free();
  m0 = m;
  
  while (1)
  {
    //cout << nd << "," << m << "\n\n";	  
    e = 0;
    while (e < ne)
    {
      if (virgin[e])
      {
        elem[ndEl[e]]->givePlotSequence2D(plotSeq);
        np = plotSeq.n;
        plotSeq.append(plotSeq[0]);
        j = 0;
        while (plotSeq[j] != m && j < np) j++;
        if (j < np && plotSeq[j+1] == nd) 
	{
	  j += 2; if (j > np) j = 1;
          m = plotSeq[j];
	  if (m == m0) return false;
          virgin[e] = false;
          plotSeq.free();
          break;
	}
      }
      e++;
    }
    if (e == ne) { nxtNd = m; return true; }
  }
}








    
void Mesh::initialiseUDep(List<DependentDoF> *UD)
{
  int c, i, ii, j, m, nn = ndf;
  
  Vector<int>        *UBCInt = &ubcIntTmp; 
  Vector<double>     *UBCDbl = &ubcDblTmp;
  MatrixFullArray<int> *ID = &idu;
 
  if (UD == &xDep)
  {
    UBCInt = &xbcIntTmp;
    UBCDbl = &xbcDblTmp;
    ID     = &idx;
    nn     = ndm;
  }
  
  List<DependentDoF>     &uD = *UD;
  Vector<int>        &ubcInt = *UBCInt;
  Vector<double>     &ubcDbl = *UBCDbl;
  MatrixFullArray<int> &id = *ID;  
	
  c = uD.n;
 
  m = intDiv(ubcInt.n,1+nn);
  
  for (i=0; i<m; i++)
  {
    ii = ubcInt[i*(1+nn)];

    if (ii<1 || ii>numnp) prgError(1,"Mesh::initialiseUDep","invalid node number!");
    
    for (j=1; j<nn+1; j++)
    {
      if (id(ii,j) != 0 && abs(ubcDbl[i*nn+j-1]) > 1.e-10) 
      {
        uD.add(new DependentDoF);
	      
	uD[c].nd      = ii;
	uD[c].dof     = j;
	uD[c].ucBase  = ubcDbl[i*nn+j-1];
	uD[c++].tmFct = ubcInt[i*(1+nn)+j];
      }
    }
  }
  ubcInt.free();
  ubcDbl.free();

  return;
}





  

void Mesh::prepareAndCheckUDep(List<DependentDoF> *UD)
{
  char fct[] = "Mesh::prepareAndCheckUDep";
  
  int i, j, ii, jj, m, nn = ndf;
  MatrixFullArray<int> *ID = &idu;
  if (UD == &xDep) { nn = ndm; ID = &idx; }
  List<DependentDoF>     &uD = *UD;
  MatrixFullArray<int> &id = *ID;  
 
  // set master index
  
  for (i=0; i<uD.n; i++)
  {
    uD[i].master.setDim(uD[i].masterNd.n);

    for (j=0; j<uD[i].masterNd.n; j++)
    {
      ii = uD[i].masterNd[j];
      jj = uD[i].masterDoF[j];
      if (ii < 1 || ii > numnp) prgError(1,fct,"dependent displacements: invalid node number!");
      if (jj < 1 || jj > nn)   
	prgError(1,fct,"dependent displacements: invalid degree of freedom!");
      m = id(ii,jj);
      if (m == 0) 
        prgError(1,fct,"dependent / master displacement inconsistency!");
      if (m < 0) if (uD[-1-m].master.n > 0) 
        prgError(2,fct,"dependent / master displacement inconsistency!");
      uD[i].master[j] = m;
    }
  }

  // check distinctness of slaves

  return;

  for (i=0; i<uD.n-1; i++)
    for (j=i+1; j<uD.n; j++)
      if (uD[i].nd == uD[j].nd && uD[i].dof == uD[j].dof)
      {
	//cout << "   " << uD[i].nd << "," << uD[i].dof << "    ";
	prgError(2,fct,"double dependent displacement!");
      }

  return;
}










void Mesh::updateUDepInc_hlp(List<DependentDoF> &uD, MatrixFullArray<int> &id)
{
  int i, j, ii;

  for (i=0; i<uD.n; i++) uD[i].timeUpdate();

  for (i=0; i<uD.n; i++) 
  {
    for (j=0; j<uD[i].master.n; j++)
    {
      ii = id(uD[i].masterNd[j],uD[i].masterDoF[j]);

      if (ii < 0)
      {
        ii = - 1 - ii;

	uD[i].duc += uD[ii].duc * uD[i].alpha[j];
      }
    }
  }
  return;
}












void Mesh::replaceUDepTmFct(List<DependentDoF> *UD)
{
  int i, j;
  
  List<DependentDoF> &uD = *UD;
  
  for (i=0; i<uD.n; i++)
  {
    if (uD[i].tmFct != 0)	  
    {
      j = 0; while (j < timeFunction.n && uD[i].tmFct != timeFunction[j].id) j++;

      if (j == timeFunction.n) prgError(1,"Mesh::replaceUDepTmFct",
		               "invalid time function id in uDep or xDep data!");
      uD[i].tmFct = j;
    }
    else uD[i].tmFct = -1;

    uD[i].init();
  }
  return;
}
 













void Mesh::replaceFrcTmFct(void)
{
  int i, j;
  
  for (i=0; i<frcTmFct.n; i++)
  {
    if (frcTmFct[i] != 0)	  
    {
      j = 0; while (j<timeFunction.n && frcTmFct[i]!=timeFunction[j].id) j++;

      if (j == timeFunction.n) prgError(1,"Mesh::replaceFrcTmFct",
		               "invalid time function id in 'forces'!");
      frcTmFct[i] = j;
    }
    else frcTmFct[i] = -1;
  }

  return;
}











bool Mesh::sameFace(Vector<int> &face1, Vector<int> &face2)
{
  if (face1.n != face2.n) return false;
  
  int i, j;

  for (i=0; i<face1.n; i++)
  {
    j = 0;
    while (j<face2.n && face1[i] != face2[j]) j++;

    if (j == face2.n) return false;
  }

  return true;
}












void Mesh::separateBndNodes(Vector<int> &queue, Vector<int> &tmpBndNd, VectorArray<int> &nd2bnd,
                                                VectorArray<unsigned int> &faceBits)
{
  int e, i, j, k, q = queue[0];

  bool addFlag;

  queue.del(0);

  Vector<int> face;

  nodeFlag[q] = false;

  nd2bnd[q] = tmpBndNd.n;

  tmpBndNd.append(q+1);

  for (i=0; i<nodeElem[q].n; i++)  // loop over all elements attached to node q
  {
    e = nodeElem[q][i];

    for (j=0; j<elem[e]->nFaces(); j++) // loop over all faces of the current element
    {
      if (faceBits[e] & pow2[j])  // if the face is on the surface and active
      {
        elem[e]->giveFace3D(j+1,face);

        k=0; while (k < face.n) if (face[k]!=q+1) k++; else break;

        if (k<face.n) // if node q is one of the nodes defining the face then ...
        {
          addFlag = true;

          for (k=0; k<face.n; k++) // loop over the nodes of the face
          {
            if (nodeFlag[face[k]-1]) // if the node has not been visited before ...
            {
              addFlag = false;
              if (!queue.contains(face[k]-1)) queue.append(face[k]-1); // ... add it to the queue
            }
          }
          if (addFlag)
          {
            surf3D->addElementFace(e,j+1,face,nd2bnd); // ... add it to the object surface

            faceBits[e] = faceBits[e] - pow2[j];       // and make it inactive
          }
        }
      }
    } 
  }
  return;
}












void Mesh::generateNodeNodeConnectivity(void)
{
  int e, i, ii, j, jj, nen1;

  Vector<int> *iTmpPtr;

  if (nodeNode != NULL) delete [] nodeNode;

  nodeNode = new VectorArray<int> [numnp];
  iTmpPtr  = new Vector<int>      [numnp];
  
  for (e=0; e<numel; e++)
  {
    nen1 = elem[e]->nen();

    for (i=0; i<nen1; i++)
    {
      ii = elem[e]->ix[i];
      
      for (j=0; j<nen1; j++)
        if (j!=i)
	{
          jj = elem[e]->ix[j];

          if (!iTmpPtr[ii-1].contains(jj)) iTmpPtr[ii-1].append(jj);
	}
    }
  }
  for (i=0; i<numnp; i++) nodeNode[i] = iTmpPtr[i];

  delete [] iTmpPtr;

  return;
}











void Mesh::generateNodeElemConnectivity(void)
{
  int e, i, ii, j, nen1;

  Vector<int> *iTmpPtr;

  if (nodeElem != NULL) delete [] nodeElem;

  nodeElem = new VectorArray<int> [numnp];
  iTmpPtr  = new Vector<int>        [numnp];
  
  for (e=0; e<numel; e++)
  {
    nen1 = elem[e]->nen();

    for (i=0; i<nen1; i++)
    {
      ii = elem[e]->ix[i] - 1;
      iTmpPtr[ii].append(e);
    }
  } 
  for (i=0; i<numnp; i++) 
  { 
    nodeElem[i].setDim(iTmpPtr[i].n); 
    for (j=0; j<iTmpPtr[i].n; j++) nodeElem[i][j] = iTmpPtr[i][j];
  }
  delete [] iTmpPtr;

  return;
}












void Mesh::generateBoundaryData2D(Vector<int> &tmpBndNd, Vector<int> &tmpBnd)
{
  VectorArray<bool> isBnd;

  double *X = x.x;

  int i, ii, j, iindm;

  // find one node on outer boundary (smallest first coordinate)

  j = 0; while (j < numnp) { if (isBoundaryNode2D(j+1,ii)) break; j++; }

  if (j == numnp) return;  // -> all elements might be 'open' (e.g. trusses or beams)

  isBnd.setDim(numnp);

  for (i=0; i<numnp; i++)  isBnd[i] = false;

  tmpBndNd.free();

  tmpBndNd.append(j); j *= ndm;
  
  iindm = j;

  for (i=tmpBndNd[0]+1; i<numnp; i++) 
  {  
    iindm += ndm;
    if (X[iindm] < X[j]) if (isBoundaryNode2D(i+1,ii)) { tmpBndNd[0] = i; j = iindm; }
  }

  tmpBndNd[0]++;
  
  isBnd[tmpBndNd[0]-1] = isBoundaryNode2D(tmpBndNd[0],ii);

  // get ordered list of boundary nodes

  j = 0;

  tmpBnd.free();

  tmpBnd.append(0);

  while (j < numnp)
  {
    tmpBndNd.append(ii);
    while (1)
    {
      i = ii;
      isBnd[i-1] = isBoundaryNode2D(ii,ii);
      if (ii == tmpBndNd[tmpBnd.lastCoeff()]) break;
      else tmpBndNd.append(ii);
    }
    ii = 0;
    while (ii == 0 && j < numnp)
    {
      while (j < numnp) if (isBnd[j]) j++; else break;
      if (j < numnp)
      { 
        if (isBoundaryNode2D(j+1,ii)) 
        { 
          tmpBnd.append(tmpBndNd.n);
          tmpBndNd.append(j+1);
          isBnd[j] = true;
        }
        else j++;
      }
    }
  }
  return;
}













void Mesh::generateBoundaryData3D(Vector<int> &tmpBndNd, Vector<int> &tmpBnd)
{
  char fct[] = "Mesh::generateBoundaryData3D";

  VectorArray<unsigned int> faceBits;

  Vector<int> face, face1, face2;

  double *X = x.x;

  int e, e1, e2, i, i1, i2, ii, j, jj, k, m, n, el1, el2;

  // find all boundary faces

  faceBits.setDim(numel);
  for (e=0; e<numel; e++)
  {
    switch (elem[e]->nFaces())
    {
      case  4: faceBits[e] = 15; break;
      case  6: faceBits[e] = 63; break;
      default: prgError(1,fct,"invalid number of element faces!");
    }
  }

  for (i=0; i<numnp; i++) nodeFlag[i] = false;

  for (i=0; i<numnp; i++)
  {
    for (el1=0; el1<nodeElem[i].n-1; el1++)
    {
      e1 = nodeElem[i][el1];

      for (i1=0; i1<elem[e1]->nFaces(); i1++)
      {
        elem[e1]->giveFace3D(i1+1,face1);

        el2 = el1 + 1;
        while (el2<nodeElem[i].n && faceBits[e1] & pow2[i1])
        {
          e2 = nodeElem[i][el2];

          i2 = 0;
          while (i2<elem[e2]->nFaces())
          {
            if (faceBits[e2] & pow2[i2])
            {
              elem[e2]->giveFace3D(i2+1,face2);
              if (sameFace(face1,face2))
              {
                faceBits[e1] = faceBits[e1] ^ pow2[i1];
                faceBits[e2] = faceBits[e2] ^ pow2[i2];

                break;
              }
            }
            i2++;
          }
          el2++;
        }
      }
    }
  }

  // mark all boundary nodes

  for (e1=0; e1<numel; e1++)
  {
    for (i1=0; i1<elem[e1]->nFaces(); i1++)
    {
      if (faceBits[e1] & pow2[i1])
      {
        elem[e1]->giveFace3D(i1+1,face1);

        for (j=0; j<face1.n; j++) nodeFlag[face1[j]-1] = true;
      }
    }
  }

  // find boundary node on external boundary

  i = 0; while (!nodeFlag[i]) i++;

  j = i; while (j < numnp) { if (nodeFlag[j]) if (x.x[j*3] < x.x[i*3]) i = j; j++; }

  // separate nodes on different boundary surfaces

  surf3D = new ObjectSurface(this);

  VectorArray<int> nd2bnd; 

  nd2bnd.setDim(numnp);

  Vector<int> queue;

  queue.append(i);

  while (1)
  {
    tmpBnd.append(tmpBndNd.n);

    while (queue.n > 0) separateBndNodes(queue,tmpBndNd,nd2bnd,faceBits);

    i = 0; while (i<numnp && !nodeFlag[i]) i++;

    if (i == numnp) break; else queue.append(i);
  }

  surf3D->finalise(tmpBndNd);

  return;
}









