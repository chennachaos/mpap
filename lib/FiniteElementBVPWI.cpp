
#include <iostream>

#include "FiniteElementBVPWI.h"
#include "DomainTree.h"
#include "DataBlockTemplate.h"
#include "Plot.h"
#include "ComputerTime.h"
#include "SolverMA41.h"
#include "MathGeom.h"
#include "InterfaceMatch.h"



extern ComputerTime computerTime;
extern DomainTree   domain;
extern Plot         plot;



using namespace std;



FiniteElementBVPWI::FiniteElementBVPWI(void)                       
{ 
  show = false;

  membraneFlag = false;

  isConnected = false;

  dxdu = .5; // to be overwritten by interface routines

  ctimDeriv1 = 0.;  doneDeriv1 = false;
  ctimDeriv2 = 0.;  doneDeriv2 = false;
  ctimDeriv3 = 0.;  doneDeriv3 = false;
  ctimDeriv4 = 0.;  doneDeriv4 = false;
  ctimDeriv5 = 0.;  doneDeriv5 = false;
  ctimDeriv6 = 0.;  doneDeriv6 = false;
  ctimDeriv7 = 0.;  doneDeriv7 = false;
  ctimDeriv8 = 0.;  doneDeriv8 = false;

  ctimBack3  = 0.;
  ctimBack5  = 0.;
  ctimBack8  = 0.;

  ctimElim   = 0.;  doneElim   = false;
  ctimPrepExt= 0.;  donePrepExt= false;

  // add new type
  
  DomainType *finiteElementBVPWI = domain.newType(FINITEELEMENTBVPWI,FINITEELEMENTBVP);

  if (finiteElementBVPWI == NULL) return;  // domain type exists already

  finiteElementBVPWI->key.addNew("interface mesh nodes only",
                                 "interface nodes",
                                 "interface interpolations",
                                 "rigid body interface",
                                 "show",
                                 "setMembraneFlag");
  return;
}







	                                          
FiniteElementBVPWI::~FiniteElementBVPWI(void)                     
{
  for (int i=0; i<n2bl.n; i++) if (n2bl[i] != NULL) delete n2bl[i];

  return;
}










void FiniteElementBVPWI::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString tmpl, *word;
 
  char tmp[30], fct[] = "FiniteElementBVPWI::readInputData";

  int fD, mD, nw, e, i, j, k, m, n, ku, kx, ndfRBI = 3;

  if (ndm == 3) ndfRBI = 6; 

  double xTmp[10];

  Vector<double> dTmp;
  Vector<int>    iTmp, lTmp;
  Vector<void*>  vTmp, v2Tmp;
  MyStringList   sTmp;

  List< Vector<int> > lviTmp;

  DataBlockTemplate t1, t2;

  switch (domain[FINITEELEMENTBVPWI].key.whichBegins(line))
  {
    case  0: cout << "     FINITEELEMENTBVPWI: reading interface mesh nodes only ...\n\n";

             if (numnp < 1) prgError(3,fct,"'coor' must precede 'interface mesh nodes only'!");

             if (!prgReadColonSepListVectorInt(Ifile,line,lviTmp))

               prgError(1,fct,"invalid node number or invalid keyword!");

             if (lviTmp.n != 3) prgError(6,fct,
                "specify three parts of 'interface mesh nodes only' separated by ':'!");

             if (lviTmp[0].n != ndf) prgError(6,fct,"specify interface degree of freedoms!");

             if (lviTmp[1].n != 1 && lviTmp[1].n != 2) 

               prgError(6,fct,"specify interface impact depth!");
             
             if (lviTmp[2].n < 1)

               if (lviTmp.n == 3) prgError(6,fct,"specify interface free node numbers!");

             j = 0; for (i=0; i<ndf; i++) if (lviTmp[0][i] != 0) j++;
             n = dofU.n;
             dofU.add(new VectorArray<int>); 
             dofU[n].setDim(j);
             j = 0; for (i=0; i<ndf; i++) if (lviTmp[0][i] != 0) dofU[n][j++] = i;
    
             if (dofX.n != n) prgError(10,fct,"dofr.n != dofU.n");
             dofX.add(new VectorArray<int>); 
             dofX[n].setDim(ndm);
             for (i=0; i<ndm; i++) dofX[n][i] = i;

             if (lviTmp[1].n == 1) lviTmp[1].append(-1);
             fD = lviTmp[1][0];
             mD = lviTmp[1][1];
             if (fD < 1) prgError(1,fct,"minimum freeDepth is 1!");
             if (mD < -1) prgError(1,fct,"minimum meshDepth is -1!");

             v2Tmp.free();

             for (i=0; i<lviTmp[2].n; i++)
             {
               j = lviTmp[2][i];
               if (j < 1 || j > numnp) prgError(1,fct,"invalid mesh node number!");
               freeNodeTmp.add(new FreeNode(INTERPOLATION_FN,x.x+j*ndm-ndm,fD,mD,
                                            &(dofU[n]),&(dofX[n])));
               freeNodeIsConnected.append(true);
               vTmp.free(); vTmp.append((void*)(&freeNodeTmp[freeNodeTmp.n-1]));
               dTmp.free(); dTmp.append(1.);
               if (isALE()) bndNodeTmp.add(new BndNode(j,x.x+(j-1)*ndm,vTmp,vTmp,dTmp));
               else         bndNodeTmp.add(new BndNode(j,x.x+(j-1)*ndm,vTmp,v2Tmp,dTmp));
             }

	     break;

    case  1: cout << "     FINITEELEMENTBVPWI: reading interface nodes ...\n\n";

             if (freeNodeTmp.n > 0) prgError(1,fct,"interface nodes have already been generated!");

	     sprintf(tmp,"123 %di %df",ndm+ndf+2,ndm+3*ndf);  
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	     
             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
		     prgError(2,fct,"data error in 'interface nodes'!");
	     
             n = dTmp.n / (ndm+3*ndf);

             for (i=0; i<n; i++)
             {
               for (j=0; j<ndm; j++) xTmp[j] = dTmp[i*(ndm+3*ndf)+j];

               k = 0; for (j=0; j<ndm; j++) if (iTmp[i*(ndm+ndf+2)+j] != 0) k++;
               m = dofX.n;
               dofX.add(new VectorArray<int>); 
               dofX[m].setDim(k);
               k = 0; for (j=0; j<ndm; j++) if (iTmp[i*(ndm+ndf+2)+j] != 0) dofX[m][k++] = j;
               kx = 0; while (kx < m) if (dofX[kx] == dofX[m]) break; else kx++;
               if (kx != m) dofX.del(m);

               k = 0; for (j=0; j<ndf; j++) if (iTmp[i*(ndm+ndf+2)+ndm+j] != 0) k++;
               m = dofU.n;
               dofU.add(new VectorArray<int>); 
               dofU[m].setDim(k);
               k = 0; for (j=0; j<ndf; j++) if (iTmp[i*(ndm+ndf+2)+ndm+j] != 0) dofU[m][k++] = j;
               ku = 0; while (ku < m) if (dofU[ku] == dofU[m]) break; else ku++;
               if (ku != m) dofU.del(m);

               freeNodeTmp.add(new FreeNode(INTERPOLATION_FN,
                                            xTmp,
                                            iTmp[i*(ndm+ndf+2)+ndm+ndf],
                                            iTmp[i*(ndm+ndf+2)+ndm+ndf+1],
                                            &(dofU[ku]),
                                            &(dofX[kx])));

               freeNodeIsConnected.append(false);
             }	     

	     break;

    case  2: cout << "     FINITEELEMENTBVPWI: reading interface interpolations ...\n\n";
	   
             if (intpTmp.n > 0) prgError(2,fct,"interpolation already specified!");

             m = 9;

	     sprintf(tmp,"%di",m);
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	     
             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
		     prgError(1,fct,"data error in 'interface interpolations'!");
	
             n = intDiv(iTmp.n,m);
	     
             for (i=0; i<n; i++)
             {
               j = m - 1; while (j > 0 && iTmp[i*m+j] == 0) j--;
 
               if (j < ndm-1) prgError(2,fct,"data error in 'interface interpolations'!");

               intpTmp.add(new Vector<int>);

               for (e=0; e<=j; e++) intpTmp[i].append(iTmp[i*m+e]);
             }

	     break;

    case  3: cout << "     FINITEELEMENTBVPWI: reading rigid body interface ...\n\n";
	   
             if (numnp < 1) prgError(3,fct,"'coor' must precede 'rigid body interface'!");

             line.getNextLine(Ifile);

             nw = line.split(&word);
            
             if (nw < ndm) prgError(1,fct,"specify rigid body coordinates!");
		     
             for (i=0; i<ndm; i++) if (!word[i].toDbl(&xTmp[i])) 
               prgError(2,fct,"input error in coordinates of 'rigid body interface'!");
		     
             for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     
             if (!prgReadColonSepListVectorInt(Ifile,line,lviTmp))

               prgError(100,fct,"invalid node number or invalid keyword!");

             if (lviTmp.n != 3 && lviTmp.n != 4)
               prgError(60,fct,"data error in 'rigid body interface'!");

             if (lviTmp[0].n != ndfRBI) prgError(7,fct,"specify rigid body degrees of freedom!");

             if (lviTmp[1].n != 2) prgError(7,fct,"specify free node and mesh impact depths!");

             if (lviTmp[2].n < 1) prgError(6,fct,"specify boundary mesh nodes!");

             j = 0; for (i=0; i<ndfRBI; i++) if (lviTmp[0][i] != 0) j++;
             m = dofU.n;
             dofU.add(new VectorArray<int>); 
             dofU[m].setDim(j);
             j = 0; for (i=0; i<ndfRBI; i++) if (lviTmp[0][i] != 0) dofU[m][j++] = i;
             ku = 0; while (ku < m) if (dofU[ku] == dofU[m]) break; else ku++;
             if (ku != m) dofU.del(m);
 
             if (dofX.n != m) prgError(10,fct,"dofr.n != dofU.n");
             kx = 0; while (kx < m) if (dofX[kx] == dofU[ku]) break; else kx++;
             if (kx == m)
             {
               dofX.add(new VectorArray<int>);
               dofX[m].setDim(dofU[m].n);
               for (i=0; i<dofU[m].n; i++) dofX[kx][i] = dofU[ku][i];
             }

             if (lviTmp[1].n == 1) lviTmp[1].append(-1);
             fD = lviTmp[1][0];
             mD = lviTmp[1][1];
             if (fD < 1) prgError(1,fct,"minimum freeDepth is 1!");
             if (mD < -1) prgError(1,fct,"minimum meshDepth is -1!");

             freeNodeTmp.add(new FreeNode(RIGIDBODY_FN,xTmp,fD,mD,&dofU[ku],&dofX[kx]));
             freeNodeIsConnected.append(true);

             vTmp.free(); vTmp.append((void*)(&freeNodeTmp[freeNodeTmp.n-1]));
             dTmp.free();

             for (i=0; i<lviTmp[2].n; i++)
             {
               j = lviTmp[2][i];
               if (j < 1 || j > numnp) prgError(10,fct,"invalid mesh node number!");
               bndNodeTmp.add(new BndNode(j,x.x+(j-1)*ndm,vTmp,vTmp,dTmp));
             }

             if (lviTmp.n == 3) break;

             if (lviTmp[3].n == 0) break;

             for (i=0; i<numnp; i++) nodeFlag[i] = false;

             for (i=0; i<bndNodeTmp.n; i++) nodeFlag[bndNodeTmp[i].nd-1] = true;

             v2Tmp.free();

             for (i=0; i<lviTmp[3].n; i++)
             {
               j = lviTmp[3][i];
               if (j < 1 || j > numnp) prgError(11,fct,"invalid mesh node number!");
               if (!nodeFlag[j-1])
                 bndNodeTmp.add(new BndNode(j,x.x+(j-1)*ndm,v2Tmp,vTmp,dTmp));
             }
             
	     break;

    case  4: show = true;

             line.getNextLine(Ifile);

             break;

    case  5: membraneFlag = true;

             line.getNextLine(Ifile);

             break;

    case -1: // go and inherit from FiniteElementBVP
	     
	     this->FiniteElementBVP::readInputData(Ifile, line); 
	     
	     break;
  }

  return;
}














void FiniteElementBVPWI::prepareInputData(void)
{
  char fct[] = "FiniteElementBVPWI::prepareInputData"; 

  // call ancestor function

  FiniteElementBVP::prepareInputData();
 
  cout << "     FINITEELEMENTBVPWI: prepare input data ...\n\n";

  if (freeNodeTmp.n < 1) return;

  int e, i, j, m, n, nenEl, *IX, *DOF;

  // copy freeNodeTmp to ListArray,
  // update bndNodeTmp.freeU and bndNodeTmp.freeX accordingly,
  // set freeNode identification numbers nd

  freeNode.setDim(freeNodeTmp.n);

  for (i=0; i<freeNodeTmp.n; i++) { freeNodeTmp[i].nd = i+1; freeNode[i] = freeNodeTmp[i]; }  

  for (i=0; i<bndNodeTmp.n; i++)
  {
    for (j=0; j<bndNodeTmp[i].freeU.n; j++)
      bndNodeTmp[i].freeU[j] = freeNode.x + bndNodeTmp[i].freeU[j]->nd - 1;

    for (j=0; j<bndNodeTmp[i].freeX.n; j++)
      bndNodeTmp[i].freeX[j] = freeNode.x + bndNodeTmp[i].freeX[j]->nd - 1;
  }
  freeNodeTmp.free();

  // connect interpolation free nodes and check for unconnected free nodes

  if (freeNode.n != freeNodeIsConnected.n) prgError(1,fct,"freeNode.n != freeNodeIsConnected.n");

  for (i=0; i<intpTmp.n; i++)
    for (j=0; j<intpTmp[i].n; j++)
      freeNodeIsConnected[intpTmp[i][j]-1] = true;

  for (i=0; i<freeNode.n; i++) if (!freeNodeIsConnected[i])
    prgError(2,fct,"error in input data: unconnected interface free nodes detected!");

  // generate interpolation data if required

  init1InterpolationData();

  // bndNodeTmp is now complete. Check whether it's ok and copy to ListArray.

  for (i=0; i<bndNodeTmp.n; i++)
    if (!bndNodeTmp[i].checkOK()) prgError(20,fct,"problem in bndNode!");

  bndNode.setDim(bndNodeTmp.n);
  for (i=0; i<bndNodeTmp.n; i++) bndNode[i] = bndNodeTmp[i];

  if (show)
  {
    cout << "bndNode:\n\n  bndNode.n = " << bndNode.n << "\n\n";
    for (i=0; i<bndNode.n; i++) bndNodeTmp[i].print(ndm); cout << "\n";
  }

  bndNodeTmp.free();

  // recalculate idu (generate and take into account uDep data for bnd nodes) and
  // recalculate idx (generate and take into account xDep data for bnd nodes)

  for (i=0; i<numnp*ndf; i++) if (idu.x[i] < 1) idu.x[i] = 1; else idu.x[i] = 0;

  if (isALE()) for (i=0; i<numnp*ndm; i++) if (idx.x[i] < 1) idx.x[i] = 1; else idx.x[i] = 0;

  m = uDep.n;
  n = xDep.n;

  for (i=0; i<m; i++) idu(uDep[i].nd,uDep[i].dof) = - 1 - i;
  for (i=0; i<n; i++) idx(xDep[i].nd,xDep[i].dof) = - 1 - i;

  for (i=0; i<bndNode.n; i++)
  {
    bndNode[i].generateDeps(ndm,isALE());

    for (j=0; j<bndNode[i].uDep.n; j++)
    {
      uDep.add(bndNode[i].uDep[j]);
      idu(uDep[m].nd,uDep[m].dof) = - 1 - m++;
    }
    for (j=0; j<bndNode[i].xDep.n; j++)
    {
      xDep.add(bndNode[i].xDep[j]);
      idx(xDep[n].nd,xDep[n].dof) = - 1 - n++;
    }
  }

  nequ = 0;
  for (i=1; i<numnp+1; i++)
    for (j=1; j<ndf+1; j++)
      if (idu(i,j) == 0) idu(i,j) = ++nequ; else if (idu(i,j) > 0) idu(i,j) = 0;

  prepareAndCheckUDep(&uDep);

  if (isALE())
  {
    neqx = 0;
    for (i=1; i<numnp+1; i++)
      for (j=1; j<ndm+1; j++)
        if (idx(i,j) == 0) idx(i,j) = ++neqx; else if (idx(i,j) > 0) idx(i,j) = 0;

    prepareAndCheckUDep(&xDep);
  }

  // generate nodeHasU and nodeHasX

  nodeHasU.setDim(numnp);
  nodeHasX.setDim(numnp);

  for (i=0; i<numnp; i++) { nodeHasU[i] = nodeHasDoFU(i+1); nodeHasX[i] = nodeHasDoFX(i+1); }

  // check for elements with no degrees of freedom

  fixedU = true;

  if (nequ > 0)
  {
    e = 0; while (e < numel)
    {
      nenEl = elem[e]->nen();
      IX    = elem[e]->ix;
      j = 0; while (j < nenEl) if (nodeHasU[IX[j]-1]) break; else j++;
      if (j < nenEl) break; else e++;
    }
    if (e == numel) fixedU = false;
  }

  fixedX = true;

  if (neqx > 0)
  {
    e = 0; while (e < numel)
    {
      nenEl = elem[e]->nen();
      IX    = elem[e]->ix;
      j = 0; while (j < nenEl) if (nodeHasX[IX[j]-1]) break; else j++;
      if (j < nenEl) break; else e++;
    }
    if (e == numel) fixedX = false;
  }

  // generate extended node to node connectivity,
  // extended nodal connectivity is required for 'init2LayerData'

  COUT << "generateExtNodeNode ...\n\n";

  ListArray< Vector<int> > extNodeNode;

  generateExtNodeNode(extNodeNode);

  if (show)
  {
    for (i=0; i<numnp; i++) cout << i+1 << " -> " << extNodeNode[i] << "\n";
    for (i=0; i<freeNode.n; i++)
      cout << i+numnp+1 << "(" << i+1 << ") -> " << extNodeNode[i+numnp] << "\n"; cout << "\n";
  }

  // generate base tree for impact connectivities and delete temporary node node connectivities

  ListArray< List< Vector<int> > > base;

  COUT << "generateBaseTree ...\n\n";

  generateBaseTree(base,extNodeNode);

  extNodeNode.free();

  if (show)
  {
    for (i=0; i<freeNode.n; i++)
    {
      cout << freeNode[i].nd+numnp << "(" << i << ") -> " 
           << freeNode[i].freeDepth << ", " << freeNode[i].meshDepth << "\n\n  "; 
      for (j=0; j<base[i].n; j++) cout << base[i][j] << "\n  "; cout << "\n";
    }
  }

  // calculate impact connectivities, assembly data, etc

  COUT << "init2LayerData ...\n\n";

  init2LayerData(base);

  COUT << "init3AssemblyData ...\n\n";

  init3AssemblyData();

  COUT << "init4Finalise ...\n\n";

  init4finalise();

  return;
}








void FiniteElementBVPWI::prepareInteractions(void)
{
  //cout << "     FINITEELEMENTBVPWI: preparing interactions ...\n\n"; 

  // go and inherit from ancestors

  FiniteElementBVP::prepareInteractions();

  return;
}







void FiniteElementBVPWI::generateExtNodeNode(ListArray< Vector<int> > &extNodeNode)
{
  char fct[] = "FiniteElementBVPWI::generateNodeNodeTmp";

  int b, c, i, j, k, l, m, n;

  ListArray< Vector<int> > bndToFree;

  bndToFree.setDim(bndNode.n);

  VectorArray<int> nodeToBnd;

  // copy nodeNode to extNodeNode

  extNodeNode.setDim(numnp + freeNode.n);

  for (i=0; i<numnp; i++)

    for (j=0; j<nodeNode[i].n; j++) extNodeNode[i].append(nodeNode[i][j]);

  // generate nodeToBnd

  nodeToBnd.setDim(numnp+freeNode.n);

  for (i=0; i<nodeToBnd.n; i++) nodeToBnd[i] = -1;

  for (i=0; i<bndNode.n; i++) nodeToBnd[bndNode[i].nd-1] = i;

  // generate bndToFree based on freeU //and freeX

  for (i=0; i<bndNode.n; i++)
  {
    for (j=0; j<bndNode[i].freeU.n; j++) bndToFree[i].append(bndNode[i].freeU[j]->nd + numnp);

    /*if (bndNode[i].xFlag == INDX)
    {
      for (j=0; j<bndNode[i].freeX.n; j++) 
      {
        k = bndNode[i].freeX[j]->nd + numnp;

        if (!bndToFree[i].contains(k)) bndToFree[i].append(k);
      }
    }*/
  }

  // extend extNodeNode

  for (i=0; i<bndNode.n; i++)
  {
    b = bndNode[i].nd - 1;

    // link boundary node to free nodes and vice versa

    for (k=0; k<bndToFree[i].n; k++) 
    {
      m = bndToFree[i][k] - 1;

      if (!extNodeNode[b].contains(m+1)) extNodeNode[b].append(m+1);
      if (!extNodeNode[m].contains(b+1)) extNodeNode[m].append(b+1);
    }

    // make all nodes that are linked with the boundary node point to its master free nodes
    // and vice versa

    for (j=0; j<nodeNode[b].n; j++)
    {
      c = nodeNode[b][j] - 1;

      if (!extNodeNode[c].contains(b+1)) prgError(1,fct,"fatal error!");

      for (k=0; k<bndToFree[i].n; k++)
      {
        m = bndToFree[i][k] - 1;

        if (!extNodeNode[m].contains(c+1)) extNodeNode[m].append(c+1);
        if (!extNodeNode[c].contains(m+1)) extNodeNode[c].append(m+1);
      }

      // if node c is a boundary node, interconnect the master free nodes

      c = nodeToBnd[c];

      if (c > -1)
      {
        for (k=0; k<bndToFree[i].n; k++)
        {
          m = bndToFree[i][k] - 1;

          for (l=0; l<bndToFree[c].n; l++)
          {
            n = bndToFree[c][l] - 1;

            if (!extNodeNode[m].contains(n+1)) extNodeNode[m].append(n+1);
          }
        }
      }
    }
  }
  return;
}










void FiniteElementBVPWI::generateBaseTree(ListArray< List< Vector<int> > > &base,
                                          ListArray< Vector<int> > &extNodeNode)
{
  int i, j, d;

  base.setDim(freeNode.n);

  for (i=0; i<freeNode.n; i++)
  {
    base[i].add(new Vector<int>);
    base[i][0].append(freeNode[i].nd + numnp);
  }

  nodeFlag.setDim(numnp + freeNode.n);

  for (i=0; i<nodeFlag.n; i++) nodeFlag[i] = false;

  cout << "                        ";

  for (i=0; i<freeNode.n; i++) 
  {
    printf("\b\b\b\b\b\b\b\b\b\b\b\b\b%5d / %5d",i+1,freeNode.n);  

    nodeFlag[freeNode[i].nd + numnp - 1] = true;

    d = max(freeNode[i].freeDepth,freeNode[i].meshDepth);

    generateBaseTreeHelp(base[i],extNodeNode,d);
   
    nodeFlag[freeNode[i].nd + numnp - 1] = false;
 
    while (base[i][base[i].n-1].n == 0) base[i].del(base[i].n-1);

    base[i].del(0);

    base[i][0].append(freeNode[i].nd + numnp);

    freeNode[i].freeDepth = min(freeNode[i].freeDepth,base[i].n);
    freeNode[i].meshDepth = min(freeNode[i].meshDepth,base[i].n);
  }
  cout << "\n\n";

  nodeFlag.setDim(numnp);

  return;
}










void FiniteElementBVPWI::generateBaseTreeHelp(List< Vector<int> > &tmp, 
                                              ListArray< Vector<int> > &extNodeNode,
                                              int depth)
{
  int i, j, ii, jj, l = 0;

  Vector<int> *tmp1 = &(tmp[0]), *tmp2, flagChanged;

  while (l < depth+1 && tmp1->n > 0)
  {
    tmp2 = new Vector<int>;

    tmp.add(tmp2);
  
    for (i=0; i<tmp1->n; i++) 
    { 
      ii = (*tmp1)[i];

      for (j=0; j<extNodeNode[ii-1].n; j++)
      {
	jj = extNodeNode[ii-1][j];
	
        if (!nodeFlag[jj-1])
        {
          tmp2->append(jj);
          nodeFlag[jj-1] = true;
          flagChanged.append(jj-1);
        }
      }
    }
    l++;
    tmp1 = tmp2;
  }
  for (i=0; i<flagChanged.n; i++) nodeFlag[flagChanged[i]] = false;

  return;
}








void FiniteElementBVPWI::init1InterpolationData(void)
{
  if (intpTmp.n == 0) return;

  char fct[] = "FiniteElementBVPWI::init1InterpolationData";

  if (ndm != 2 && ndm != 3) prgError(1,fct,"only for 2D or 3D!");

  int e, i;

  VectorArray<bool> nodeIsBndNode;

  // check order of interpolation

  for (e=0; e<intpTmp.n; e++)
    if (intpTmp[e].n != ndm) prgError(3,fct,"so far, linear interpolations only (n = ndm)!");

  // check consistency of dofU and dofX for interface interpolation elements

  for (e=0; e<intpTmp.n; e++)
    for (i=1; i<intpTmp[e].n; i++)
    {
      if (freeNode[intpTmp[e][0]-1].dofU != freeNode[intpTmp[e][i]-1].dofU)
        prgError(4,fct,"dofU inconsistency in interface interpolations!");
      if (freeNode[intpTmp[e][0]-1].dofX != freeNode[intpTmp[e][i]-1].dofX)
        prgError(4,fct,"dofX inconsistency in interface interpolations!");
    }

  // prepare nodeFlag (true for all boundary nodes) and initialise nodeIsBndNode,

  nodeIsBndNode.setDim(numnp);

  for (i=0; i<numnp; i++) { nodeFlag[i] = false; nodeIsBndNode[i] = false; }

  for (i=0; i<bndNodeTmp.n; i++) nodeIsBndNode[bndNodeTmp[i].nd-1] = true;

  if (nen == 2) for (i=0; i<numnp; i++) nodeFlag[i] = true;  // for trusses & beams

  for (i=0; i<nBndNd; i++) nodeFlag[bndNd[i]-1] = true;

  // generate interpolation bnd nodes

  findBndNodesForInterpolations(nodeIsBndNode,true);

  if (membraneFlag)
  {
    for (i=0; i<numnp; i++) if (nodeIsBndNode[i]) nodeFlag[i] = false;

    findBndNodesForInterpolations(nodeIsBndNode,false);
  }

  return;
}









void FiniteElementBVPWI::init2LayerData(ListArray< List< Vector<int> > > &base)
{
  char fct[] = "FiniteElementBVPWI::init2LayerData"; 

  int d, e, e2, i, j, k, l, m, n, q,
      *IX, nenEl;

// generate nodeToDepth

  if (show) cout << " ---- generate nodeToDepth ----\n\n";

  VectorArray<int> nodeToDepth;

  nodeToDepth.setDim(numnp+freeNode.n);

  for (i=0; i<numnp+freeNode.n; i++) nodeToDepth[i] = -1;

  m = 1;

  for (i=0; i<freeNode.n; i++) m = max(m,max(freeNode[i].meshDepth,freeNode[i].freeDepth));

  for (d=0; d<m; d++)
  {
    for (i=0; i<freeNode.n; i++)
    {
      if (base[i].n > d)
      {
        for (j=0; j<base[i][d].n; j++)
        {
          k = base[i][d][j] - 1;

          if (nodeToDepth[k] < 0) nodeToDepth[k] = d + 1;
        }
      }
    }
  }

  if (show) cout << nodeToDepth << "\n\n";

  // generate layElem

  if (show) cout << " ---- generate layElem ----\n\n";

  ListArray< Vector<int> > layElemTmp;

  layElemTmp.setDim(m);

  for (e=0; e<numel; e++)
  {
    IX    = elem[e]->ix;
    nenEl = elem[e]->nen();

    d = - 1;
    j = 0;
    while (j < nenEl)
    {
      if (nodeToDepth[IX[j]-1] == -1) { d = -1; break; }
      else if (d < nodeToDepth[IX[j]-1]) d = nodeToDepth[IX[j]-1];
      j++;
    }

    if (d > 0 && d < m+1) layElemTmp[d-1].append(e+1);
  }

  while (m > 0 && layElemTmp[m-1].n == 0) m--;
  if (m == 0) prgError(1,fct,"layElem is empty!");

  // sort layElem[0]; elements without free mesh motion come first

  nLayElem1 = 0;
  nLayElem2 = layElemTmp[0].n;
  nLayElem3 = nLayElem2;

  for (e2=0; e2<layElemTmp[0].n; e2++)
  {
    e = layElemTmp[0][e2] - 1;

    IX    = elem[e]->ix;
    nenEl = elem[e]->nen();

    j = 0; while (j < nenEl) if (nodeHasX[IX[j]-1]) break; else j++;

    if (j == nenEl) { nLayElem1++; layElemTmp[0].move(e2,0); }
  }

  // extend layElem[0] by all elements containing bndNode[IX].freeX.n > 0

  for (e=0; e<numel; e++) elem[e]->flag = false;

  for (e=0; e<layElemTmp[0].n; e++) elem[layElemTmp[0][e]-1]->flag = true;

  for (i=0; i<bndNode.n; i++)
  {
    if (bndNode[i].freeX.n > 0)
    {
      m = bndNode[i].nd - 1;
    
      if (nodeToDepth[m] > -1)
      {
        for (k=0; k<nodeElem[m].n; k++)
        {
          e = nodeElem[m][k];
    
          if (!elem[e]->flag)
          {
            IX    = elem[e]->ix;
            nenEl = elem[e]->nen();

            j = 0; while (j < nenEl) if (nodeHasX[IX[j]-1]) break; else j++;

            if (j == nenEl) layElemTmp[0].append(e+1);

            else layElemTmp[0].insert(e+1,nLayElem3++);

            elem[e]->flag = true; 
          }
        }
      }
    }
  }
  layElem = layElemTmp; layElemTmp.free();

  if (show) { for (d=0; d<layElem.n; d++) cout << layElem[d] << "\n"; cout << "\n"; }
  if (show) { cout << nLayElem1 << ", " << nLayElem2 << ", " << nLayElem3 << "\n\n"; }

  // generate nodeToMesh

  if (show) cout << "---- nodeToMesh ----\n\n";

  Vector<int> nodeToMeshTmp;

  nodeToMesh.setDim(numnp);

  for (i=0; i<numnp; i++)
  {
    if (nodeToDepth[i] != -1 && nodeHasX[i])
    {
      for (j=0; j<nodeNode[i].n; j++)
      { 
        m = nodeNode[i][j] - 1;

        if (nodeToDepth[m] != -1 && nodeHasU[m]) nodeToMeshTmp.append(m);
      }
      if (nodeHasU[i]) nodeToMeshTmp.append(i);

      nodeToMesh[i].setDim(nodeToMeshTmp.n);

      for (j=0; j<nodeToMeshTmp.n; j++) nodeToMesh[i][j].nd = nodeToMeshTmp[j];

      nodeToMeshTmp.free();
    }
  }

  if (show)
  {
    for (i=0; i<numnp; i++)
    {
      if (nodeToMesh[i].n > 0) cout << i+1 << ": " << nodeToMesh[i][0].nd+1;
      for (j=1; j<nodeToMesh[i].n; j++) cout << ", " << nodeToMesh[i][j].nd+1; cout << "\n";
    }
  }

  // generate ixbx, ixbu, ixlx, ixlu

  if (show) cout << " ---- ixbx, ixbu, ixlx, ixlu ----\n\n";

  ixbu.setDim(nLayElem2);
  ixlu.setDim(layElem[0].n);
  ixbx.setDim(layElem[0].n);
  ixlx.setDim(nLayElem3);

  VectorArray<int> nodeToBndNodeX, nodeToLayNodeX,
                   nodeToBndNodeU, nodeToLayNodeU;

  nodeToBndNodeX.setDim(numnp);
  nodeToLayNodeX.setDim(numnp);
  nodeToBndNodeU.setDim(numnp);
  nodeToLayNodeU.setDim(numnp);

  for (i=0; i<numnp; i++) { nodeToBndNodeU[i] = -1; nodeToBndNodeX[i] = -1; }

  for (i=0; i<bndNode.n; i++)
  {
    if (bndNode[i].freeU.n > 0) nodeToBndNodeU[bndNode[i].nd-1] = i;
    if (bndNode[i].freeX.n > 0) nodeToBndNodeX[bndNode[i].nd-1] = i;
  }

  for (i=0; i<numnp; i++)
  {
    if (nodeHasU[i]) nodeToLayNodeU[i] = i; else nodeToLayNodeU[i] = -1;
    if (nodeHasX[i]) nodeToLayNodeX[i] = i; else nodeToLayNodeX[i] = -1;
  }

  for (e2=0; e2<layElem[0].n; e2++)
  {
    e = layElem[0][e2] - 1;

    nenEl = elem[e]->nen();
    IX    = elem[e]->ix;

    if (show) { cout << e + 1 << ": "; for (j=0; j<nenEl; j++) cout << IX[j] << " "; }

    if (e2 < nLayElem2)
    {
      ixbu[e2].setDim(nenEl);

      for (j=0; j<nenEl; j++) { k = IX[j] - 1; ixbu[e2][j] = nodeToBndNodeU[k]; }

      if (show) { cout << " ixbu -> "; for (j=0; j<nenEl; j++) cout << ixbu[e2][j] << " "; }
    }
    if (e2 < nLayElem3)
    {
      ixlx[e2].setDim(nenEl);

      for (j=0; j<nenEl; j++) { k = IX[j] - 1; ixlx[e2][j] = nodeToLayNodeX[k]; }

      if (show) { cout << " ixlx -> "; for (j=0; j<nenEl; j++) cout << ixlx[e2][j] << " "; }
    }

    ixlu[e2].setDim(nenEl);
    ixbx[e2].setDim(nenEl);

    for (j=0; j<nenEl; j++)
    {
      k = IX[j] - 1;
      ixbx[e2][j] = nodeToBndNodeX[k];
      ixlu[e2][j] = nodeToLayNodeU[k];
    }

    if (show)
    {
      cout << " ixbx -> "; for (j=0; j<nenEl; j++) cout << ixbx[e2][j] << " ";
      cout << " ixlu -> "; for (j=0; j<nenEl; j++) cout << ixlu[e2][j] << " ";
      cout << "\n";
    }
  }
  if (show) cout << "\n";

  nodeToBndNodeX.free();
  nodeToBndNodeU.free();
  nodeToLayNodeU.free();

  nodeToDepth.free();

  // generate freeNode.toFree

  if (show) cout << " ---- freeNode.toFree ----\n\n";

  List< Vector<int> > deriv;

  //if (nequ > 0 || neqx > 0)
  {
    for (i=0; i<freeNode.n; i++)
    {
      for (j=0; j<freeNode[i].freeDepth; j++)
      {
        deriv.add(new Vector<int>);
        for (k=0; k<base[i][j].n; k++)
        {
          m = base[i][j][k] - 1 - numnp;
          if (m > -1) deriv[j].append(m);
        }
      }
      while (freeNode[i].freeDepth > 0)
      {
        if (deriv[deriv.n-1].n == 0) { deriv.del(deriv.n-1); freeNode[i].freeDepth--; } else break;
      }
      generateDependencies(freeNode[i].toFree,deriv); deriv.free();
    }
  }
  /*else
  {
    for (i=0; i<freeNode.n; i++)
    {
      deriv.add(new Vector<int>);
      for (k=0; k<base[i][0].n; k++)
      {
        m = base[i][0][k] - 1 - numnp;
        if (m > -1) deriv[0].append(m);
      }
      freeNode[i].freeDepth = 1;

      generateDependencies(freeNode[i].toFree,deriv); deriv.free();
    }
  }*/

  maxFreeDepth = freeNode[0].freeDepth;

  for (i=1; i<freeNode.n; i++) maxFreeDepth = max(maxFreeDepth,freeNode[i].freeDepth);

  if (show)
  {
    for (i=0;i<freeNode.n;i++) 
      { cout << i << "\n"; printDependencies(freeNode[i].toFree); cout << "\n"; }
    cout <<"\n";
  }

  // generate freeNode.toMesh

  if (show) cout << " ---- generate freeNode.toMesh ----\n\n";

  for (i=0; i<numnp; i++) nodeFlag[i] = nodeHasX[i];

  for (i=0; i<freeNode.n; i++)
  {
    if (freeNode[i].meshDepth > 0)
    {
      for (j=0; j<freeNode[i].meshDepth; j++)
      {
        deriv.add(new Vector<int>);
        for (k=0; k<base[i][j].n; k++)
        {
          m = base[i][j][k] - 1;
          if (m < numnp) if (nodeFlag[m]) deriv[j].append(m);
        }
      }
      while (deriv.n > 0) if (deriv[deriv.n-1].n == 0) deriv.del(deriv.n-1); else break;

      k = 0; while (k < deriv.n) if (deriv[k].n > 0) k++; else deriv.trunc(k);

      freeNode[i].meshDepth = deriv.n;

      generateDependencies(freeNode[i].toMesh,deriv); deriv.free();
    }
  }

  maxMeshDepth = freeNode[0].meshDepth;

  for (i=1; i<freeNode.n; i++) maxMeshDepth = max(maxMeshDepth,freeNode[i].meshDepth);

  if (show)
  {
    for (i=0;i<freeNode.n;i++)
      { cout << i << "\n"; printDependencies(freeNode[i].toMesh); cout << "\n";}
    cout <<"\n";
  }

  // generate freeNode.toBLay

  if (show) cout << " ---- freeNode.toBLay ----\n\n";

  for (i=0; i<numnp; i++) if (nodeHasU[i]) nodeFlag[i] = true;
  
  for (i=0; i<freeNode.n; i++)
  {
    deriv.add(new Vector<int>);
    for (k=0; k<base[i][0].n; k++)
    {
      m = base[i][0][k];
      if (m < numnp + 1) if (nodeFlag[m-1]) deriv[0].append(m);
      //if (m < numnp + 1) deriv[0].append(m);
    }
    freeNode[i].toBLay.setDim(deriv[0].n);
    for (k=0; k<deriv[0].n; k++) freeNode[i].toBLay[k].setData(deriv[0][k]);
    deriv.free();
  }

  base.free();

  if (show)
  {
    for (i=0;i<freeNode.n;i++) 
    {
      if (freeNode[i].toBLay.n < 1) cout << i << ": { ";
      else 
      {
        cout << i << ": {" << freeNode[i].toBLay[0].nd;
        for (j=1; j<freeNode[i].toBLay.n; j++) cout << "," << freeNode[i].toBLay[j].nd;
      }
      cout <<"}\n";
    }
    cout <<"\n";
  }

  // generate freeNode.toBLayMesh

  if (show) cout << " ---- freeNode.toBLayMesh ----\n\n";

  for (i=0; i<numnp; i++) nodeFlag[i] = nodeHasX[i];
  
  deriv.free(); for (i=0; i<freeNode.n; i++) deriv.add(new Vector<int>);

  for (i=0; i<bndNode.n; i++)
  {
    if (bndNode[i].freeX.n > 0)
    {
      m = bndNode[i].nd - 1;

      for (j=0; j<nodeNode[m].n; j++)
      {
        n = nodeNode[m][j];

        if (nodeFlag[n-1])
        {
          for (k=0; k<bndNode[i].freeX.n; k++)
          {
            l = bndNode[i].freeX[k]->nd - 1;

            if (!deriv[l].contains(n)) deriv[l].append(n);
          }
        }
      }
    }
  }

  for (i=0; i<freeNode.n; i++)
  {
    freeNode[i].toBLayMesh.setDim(deriv[i].n);
    for (k=0; k<deriv[i].n; k++) freeNode[i].toBLayMesh[k].setData(deriv[i][k]);
  }
  deriv.free();

  if (show)
  {
    for (i=0;i<freeNode.n;i++) 
    {
      if (freeNode[i].toBLayMesh.n < 1) cout << i << ": { ";
      else 
      {
        cout << i << ": {" << freeNode[i].toBLayMesh[0].nd;
        for (j=1; j<freeNode[i].toBLayMesh.n; j++) cout << "," << freeNode[i].toBLayMesh[j].nd;
      }
      cout <<"}\n";
    }
    cout <<"\n";
  }

  // set freeNode.toFree.invPtr

  if (show) cout << " ---- set freeNode.toFree.invPtr ----\n\n";

  for (d=0; d<maxFreeDepth; d++)
  {
    for (i=0; i<freeNode.n; i++)
    {
      if (d<freeNode[i].freeDepth)
      {
        for (j=0; j<freeNode[i].toFree[d].n; j++)
        {
          if (freeNode[i].toFree[d][j].invPtr == NULL)
          {
            k = freeNode[i].toFree[d][j].nd;
          
            q = 0;
            while (q<freeNode[k].toFree[d].n)
              if (freeNode[k].toFree[d][q].nd == i) break; else q++;

            if (q == freeNode[k].toFree[d].n)
            {
              freeNode[i].toFree[d][j].invPtr = &(freeNode[k].toFree[d][q]);
              freeNode[k].toFree[d][q].invPtr = NULL;
              prgWarning(1,fct,"no symmetry in freeNode.toFree!");
            }
            else
            {
              freeNode[i].toFree[d][j].invPtr = &(freeNode[k].toFree[d][q]);
              freeNode[k].toFree[d][q].invPtr = &(freeNode[i].toFree[d][j]);
            }
          }
        }
      }
    }
  }

  // some stuff related to calcDerivatives6

  VectorArray<int> countThis;

  if (isALE())
  {
    // sort freeNode[].toMesh (in each layer, BLay nodes with free X come first)  
    
    nodeToLayNodeX.zero();
    
    for (i=0; i<freeNode.n; i++)
      for (j=0; j<freeNode[i].toBLay.n; j++)
      {
        k = freeNode[i].toBLay[j].nd - 1;
        if (nodeHasX[k]) nodeToLayNodeX[k] = 1;
      }
    
    for (i=0; i<freeNode.n; i++)
    {
      for (d=0; d<freeNode[i].toMesh.n; d++)
      {
        j = 0;
        while (j < freeNode[i].toMesh[d].n) 
          if (nodeToLayNodeX[freeNode[i].toMesh[d][j].nd] != 0) j++; else break;
    
        k = j++;
    
        while (j < freeNode[i].toMesh[d].n)
        {
          if (nodeToLayNodeX[freeNode[i].toMesh[d][j].nd] != 0) 
    
            freeNode[i].toMesh[d][j].swap(freeNode[i].toMesh[d][k++]);
    
          j++;
        }
      }
    }
    if (show)
    {
      for (i=0;i<freeNode.n;i++)
        { cout << i << "\n"; printDependencies(freeNode[i].toMesh); cout << "\n";}
      cout <<"\n";
    }

    // generate f2f

    f2f.setDim(freeNode.n);

    for (i=0; i<f2f.n; i++) f2f[i] = NULL;

    // allocate memory for n2bl
    
    countThis.setDim(numnp);

    countThis.zero();

    for (i=0; i<freeNode.n; i++)

      for (j=0; j<freeNode[i].toBLay.n; j++)

        countThis[freeNode[i].toBLay[j].nd-1]++;

    n2bl.setDim(numnp);

    for (i=0; i<numnp; i++)

      if (countThis[i] == 0) n2bl[i] = NULL;

      else n2bl[i] = new SomeData(countThis[i]);

    // set pointers in n2bl
 
    countThis.zero();

    for (i=0; i<freeNode.n; i++)
    {
      for (j=0; j<freeNode[i].toBLay.n; j++)
      {
        m = freeNode[i].toBLay[j].nd - 1;

        if (show) cout << i+1 << " : " << m+1 << "\n";

        n2bl[m]->bl[countThis[m]] = &(freeNode[i].toBLay[j]);

        n2bl[m]->f[countThis[m]++] = i;
      }
    }
  }

  nodeToLayNodeX.free();

  return;
}







void FiniteElementBVPWI::init3AssemblyData(void)
{
  char fct[] = "FiniteElementBVPWI::init3AssemblyData";

  int c, d, e, e2, e3, i, ir, ic, nr, nc, jr, jc, fr, fc, nenEl, *IX;

  Vector<DataToFree*>     pffuTmp, pffxTmp;
  Vector<DataToBLay*>     pfluTmp, pflxTmp, plfxTmp;
  Vector<DataToBLayMesh*> plfmTmp;

  Vector<DataToMesh*>     pmmTmp;

  // generate pffu, pflu, pffx, pflx, plfx

  pffu.setDim(nLayElem2);
  pflu.setDim(nLayElem2);
  pffx.setDim(nLayElem2);
  pflx.setDim(nLayElem2);

  plfx.setDim(layElem[0].n);

  for (e2=0; e2<layElem[0].n; e2++)
  {
    if (e2 < nLayElem2)
    {
      nenEl = ixbu[e2].n;
    
      for (ir=0; ir<nenEl; ir++)
      {
        nr = ixbu[e2][ir];
    
        if (nr > -1) // if row node is a boundary node (u)
        {
          for (jr=0; jr<bndNode[nr].freeU.n; jr++)
          {
            fr = bndNode[nr].freeU[jr]->nd-1;
 
            for (ic=0; ic<nenEl; ic++)
            {
              nc = ixbu[e2][ic];
    
              if (nc > -1) // if column node is a boundary node (u)
              {
                for (jc=0; jc<bndNode[nc].freeU.n; jc++)
                {
                  fc = bndNode[nc].freeU[jc]->nd-1;
   
                  c = 0; 
                  while (c < freeNode[fr].toFree[0].n)
    
                    if (freeNode[fr].toFree[0][c].nd != fc) c++; else break;
     
                  if (c == freeNode[fr].toFree[0].n)
    
                    prgError(1,fct,"no pair (fr,fc) in first layer of freeToFree!");
    
                  pffuTmp.append(&(freeNode[fr].toFree[0][c]));
                }
              }
    
              nc = ixbx[e2][ic];
    
              if (nc > -1) // if column node is a boundary node (x)
              {
                for (jc=0; jc<bndNode[nc].freeX.n; jc++)
                {
                  fc = bndNode[nc].freeX[jc]->nd-1;
    
                  c = 0; 
                  while (c < freeNode[fr].toFree[0].n)
    
                    if (freeNode[fr].toFree[0][c].nd != fc) c++; else break;
     
                  if (c == freeNode[fr].toFree[0].n)
    
                    prgError(2,fct,"no pair (fr,fc) in first layer of freeToFree!");
    
                  pffxTmp.append(&(freeNode[fr].toFree[0][c]));
                }
              }
    
              nc = ixlu[e2][ic];
    
              if (nc > -1) // if column node is a layer node (u)
              {
                c = 0;
                while (c < freeNode[fr].toBLay.n)
    
                  if (freeNode[fr].toBLay[c].nd != nc+1) c++; else break;
    
                if (c == freeNode[fr].toBLay.n)
    
                  prgError(1,fct,"no pair (fr,nc) in freeToBLay!");
    
                pfluTmp.append(&(freeNode[fr].toBLay[c]));
              }

              nc = ixlx[e2][ic];
    
              if (nc > -1) // if column node is a layer node (x)
              {
                c = 0;
                while (c < freeNode[fr].toBLay.n)
    
                  if (freeNode[fr].toBLay[c].nd != nc+1) c++; else break;
    
                if (c == freeNode[fr].toBLay.n)
    
                  prgError(2,fct,"no pair (fr,nc) in freeToBLay!");
    
                pflxTmp.append(&(freeNode[fr].toBLay[c]));
              }
            }
          }
        }
      }
      pffu[e2] = pffuTmp; pffuTmp.free();
      pflu[e2] = pfluTmp; pfluTmp.free();
      pffx[e2] = pffxTmp; pffxTmp.free();
      pflx[e2] = pflxTmp; pflxTmp.free();
    }

    nenEl = ixlu[e2].n;
    
    for (ir=0; ir<nenEl; ir++)
    {
      nr = ixlu[e2][ir];
    
      if (nr > -1) // if row node is a layer node (u)
      {
        for (ic=0; ic<nenEl; ic++)
        {
          nc = ixbx[e2][ic];
    
          if (nc > -1) // if column node is a boundary node (x)
          {
            for (jc=0; jc<bndNode[nc].freeX.n; jc++)
            {
              fc = bndNode[nc].freeX[jc]->nd-1;
 
              c = 0; 
              while (c < freeNode[fc].toBLay.n)
    
                if (freeNode[fc].toBLay[c].nd != nr+1) c++; else break;
    
              if (c == freeNode[fc].toBLay.n)
    
                prgError(3,fct,"no pair (fc,nr) in freeNode.toBLay!");
    
              plfxTmp.append(&(freeNode[fc].toBLay[c]));
            }
          }
        }
      }
    }
    plfx[e2] = plfxTmp; plfxTmp.free();
  }

  // generate plfm

  plfm.setDim(nLayElem3-nLayElem1);

  for (e3=0; e3<nLayElem3-nLayElem1; e3++)
  {
    e2 = e3 + nLayElem1;

    nenEl = ixbx[e2].n;
    
    for (ic=0; ic<nenEl; ic++)
    {
      nc = ixbx[e2][ic];

      if (nc > -1) // if column node is a boundary node (x)
      {
        for (jc=0; jc<bndNode[nc].freeX.n; jc++)
        {
          fc = bndNode[nc].freeX[jc]->nd-1;

          for (ir=0; ir<nenEl; ir++)
          {
            nr = ixlx[e2][ir];
       
            if (nr > -1) // if row node is a layer node (m)
            {
              c = 0; 
              while (c < freeNode[fc].toBLayMesh.n)
    
                if (freeNode[fc].toBLayMesh[c].nd != nr+1) c++; else break;
    
              if (c == freeNode[fc].toBLayMesh.n)
    
                prgError(1,fct,"no pair (fc,nr) in first layer of freeToBLayMesh!");
    
              plfmTmp.append(&(freeNode[fc].toBLayMesh[c]));
            }
          }
        }
      }
    }
    plfm[e3] = plfmTmp; plfmTmp.free();
  }

  // generate pmm

  for (d=0; d<layElem.n; d++)
  {
    for (e2=0; e2<layElem[d].n; e2++)
    {
      e = layElem[d][e2] - 1;

      nenEl = elem[e]->nen();

      IX    = elem[e]->ix;

      for (ir=0; ir<nenEl; ir++)
      {
        nr = IX[ir] - 1;

        if (nodeHasU[nr])
        {
          for (ic=0; ic<nenEl; ic++)
          {
            nc = IX[ic] - 1;

            if (nodeHasX[nc])
            {
              c = 0; 
              while (c < nodeToMesh[nc].n)

                { if (nodeToMesh[nc][c].nd == nr) break; else c++; }

              if (c == nodeToMesh[nc].n) prgError(1,fct,"fatal error!");

              pmmTmp.append(&(nodeToMesh[nc][c]));
            }
          }
        }
      }
    }
  }
  pmm = pmmTmp;

  return;
}







void FiniteElementBVPWI::init4finalise(void)
{
  int d, i, j, k, ndmi, ndmj, ndfi, ndfj;

  // allocate memory and calculate nAllDoFU and nAllDoFX

  nAllDoFU = 0;
  nAllDoFX = 0;

  for (i=0; i<freeNode.n; i++)
  {
    ndfi = freeNode[i].dofU->n;
    ndmi = freeNode[i].dofX->n;

    nAllDoFU += ndfi;
    nAllDoFX += ndmi;

    for (d=0; d<freeNode[i].toFree.n; d++)
    {
      for (j=0; j<freeNode[i].toFree[d].n; j++)
      {
        k = freeNode[i].toFree[d][j].nd;

        ndfj = freeNode[k].dofU->n;
        ndmj = freeNode[k].dofX->n;

        freeNode[i].toFree[d][j].allocMem(ndfi,ndmi,ndfj,ndmj);
      }
    }

    for (d=0; d<freeNode[i].toMesh.n; d++)
    {
      for (j=0; j<freeNode[i].toMesh[d].n; j++)
      {
        freeNode[i].toMesh[d][j].allocMem(ndmi,ndm);
      }
    }
 
    for (j=0; j<freeNode[i].toBLay.n; j++)
    {
      freeNode[i].toBLay[j].allocMem(ndf,ndfi,ndmi,ndm);
    }

    for (j=0; j<freeNode[i].toBLayMesh.n; j++)
    {
      freeNode[i].toBLayMesh[j].allocMem(ndm);
    }
  }

  for (i=0; i<nodeToMesh.n; i++)
  {
    for (j=0; j<nodeToMesh[i].n; j++)

      nodeToMesh[i][j].allocMem(ndf,ndm);
  }

  col.setDim(max(nequ,neqx));

  // generate profile vectors for u and x degrees of freedom

  profU.setDim(freeNode.n);
  profX.setDim(freeNode.n);

  profU[0] = 0;
  profX[0] = 0;

  for (i=1; i<freeNode.n; i++)
  {
    profU[i] = profU[i-1] + freeNode[i-1].dofU->n;
    profX[i] = profX[i-1] + freeNode[i-1].dofX->n;
  }

  return;
}











void FiniteElementBVPWI::setLagrangianFreeNodes(VectorArray<int> &lagrFlag)
{
  int dp, i, j;

  if (!isALE(CHECKNEQX)) return;
 
  for (i=0; i<freeNode.n; i++)
  {
    for (dp=0; dp<freeNode[i].toFree.n; dp++)

      for (j=0; j<freeNode[i].toFree[dp].n; j++)
 
        if (lagrFlag[freeNode[i].toFree[dp][j].nd] == 1)
        {
          freeNode[i].toFree[dp][j].makeLagr();
        }

    if (lagrFlag[i] == 1)
    {
      for (j=0; j<freeNode[i].toBLay.n; j++)

        freeNode[i].toBLay[j].makeLagr();

      freeNode[i].isLagr = true;
    }
  }
  return;
}











void FiniteElementBVPWI::plotInterfaceNodes(int num, int flag, int fnd, bool defFlg)
{
  char fct[] = "FiniteElementBVPWI::plotInterfaceNodes";

  int i, j, k, dd;

  double *X = x0.x, 
         d = (plot.dAct[0] + plot.dAct[1]) * .0055;

  if (defFlg) X = x.x;

  if (flag == 1)
  {
    COUT << "number of free nodes: " << freeNode.n << "\n\n";

    if (ndm < 3)
    {
      for (i=0; i<freeNode.n; i++) 
      {
        if (num) plot.point(freeNode[i].x.x,d,i+1); 
        else     plot.point(freeNode[i].x.x,d);
      }
    }
    else  prgWarning(1,fct,"no free node plotting in 3D!");
  }
  else if (flag == 2)
  {
    COUT << "number of boundary nodes: " << bndNode.n << "\n\n";

    if (ndm < 3)
    {
      for (i=0; i<bndNode.n; i++) 
      {
        j = bndNode[i].nd - 1;
        if (num) plot.point(X+j*ndm,d,j+1); 
        else     plot.point(X+j*ndm,d);
      }
    }
    else
    {
      for (i=0; i<nBndNd; i++) nodeFlag[bndNd[i]-1] = false;
      for (i=0; i<bndNode.n; i++) nodeFlag[bndNode[i].nd-1] = true;
      for (i=0; i<nBndNd; i++)
      {
        if (nodeFlag[bndNd[i]-1])
        {
          if (num) surf3D->plotNodePoint(i,d,bndNd[i]);
          else     surf3D->plotNodePoint(i,d);
        }
      }
    }
  }
  /*else if (flag == 3)
  {
    i = 0; while (i<freeNode.n) if (freeNode[i] == fnd) break; else i++;
    if (i==freeNode.n) { COUT << "The specified node is not a free node!\n\n"; return; }
    for (dd=0; dd<freeToBLay[i].n; dd++)  
      for (j=0; j<freeToBLay[i].m(dd+1); j++) 
      {
        k = freeToBLay[i][dd][j].nd - 1;
        if (num) plot.point(&(X[k*ndm]),d,k+1); 
        else     plot.point(&(X[k*ndm]),d);
      }
  }*/
  else prgError(1,fct,"invalid option flag!");

  return;
}









void FiniteElementBVPWI::eliminate(bool deriv, bool printRes)
{
  char fct[] = "FiniteElementBVPWI::eliminate";

  int i, maxIter = 20, pR;

  if (printRes) pR = 3; else pR = 1;

  //cout << fct << "\n";

  computerTime.go(fct);

  if (deriv) for (i=0; i<freeNode.n; i++) freeNode[i].zeroDerivativeMtx();

  // update uDep to provide correct boundary conditions

  for (i=0; i<bndNode.n; i++) bndNode[i].getUX(ndf,ndm,u.x,x.x,xn.x,x0.x);

  // update mesh and calculate mesh derivative

  updateMesh(1, printRes);

  if (deriv) { calculateDerivatives4(); calculateDerivatives5(); }

  // Newton-Raphson loop to solve for system response

  if (nequ > 0) if (!solverOK) prgError(2,fct,"some error in macro sequence!");

  firstIter = true;

  localStiffnessError = 0;

  i = 0;

  while (true)
  {
    if (calcStiffnessAndResidual(1+2*(pR==3 && !(nequ==0 && i==1))) != 0) 
      // argument evaluates to 1 or 3 and prevents double output of "nothing to solve for!"

      prgError(3,fct,"problem in calcStiffnessAndResidual!");

    if (factoriseSolveAndUpdate() != 0) prgError(4,fct,"problem in factoriseSolveAndUpdate!");
 
    updateIterStep();

    if (++i > maxIter) break;

    if (fixedU) { if (converged() && i > 1) break; } else if (converged()) break;
    //
    // why this fixedU stuff?? apparently required for some updates
    //
  }
  if (printRes) cout << "\n";

  if (!converged()) prgWarning(1,fct,"NO CONVERGENCE!");

  // calculate reaction forces in free nodes

  for (i=0; i<freeNode.n; i++) freeNode[i].reac.zero();

  for (i=0; i<bndNode.n; i++) bndNode[i].giveReactions(ndf,ndm,reac.x,u.x,x.x,x0.x);

  // generate stiffness matrices  d gf / d uf  and  d gf / d xf

  if (!deriv) { ctimElim += computerTime.stop(fct); doneElim = true; return; }

  calculateDerivatives1();

  calculateDerivatives2();

  calculateDerivatives6();

  calculateDerivatives7();

  calculateDerivatives3();

  calculateDerivatives8();

  ctimElim += computerTime.stop(fct);

  doneElim = true;

  return;
}






void FiniteElementBVPWI::eliminateDiffTest(double maxU, double maxX,
                                           double ddd, int dig, int dig2, bool gfrmt)
{
  char fct[] = "FiniteElementBVPWI::eliminateDiffTest";

  //cout << fct << "\n\n"; 

  if (solver == NULL) { prgWarning(1,fct,"solver == NULL!"); return; }

  int i, ii, j, ic, jc, ir, jr, cu, cx, k, m = max(nAllDoFU,nAllDoFX);

  double dd[6]   = {-3.*ddd, -2.*ddd, -ddd, +ddd, +2.*ddd, +3.*ddd };

  VectorArray<double> R, Rx, kdiff, kanly, kxdiff, kxanly, random;

  random.setDim(m);

  generateRandomDbl(random.x,nAllDoFU,-1.,+1.);

  cu = 0;
  for (i=0; i<freeNode.n; i++)
    for (j=0; j<freeNode[i].dofU->n; j++)
      freeNode[i].u[j] += random[cu++] * maxU;

  generateRandomDbl(random.x,nAllDoFX,-1.,+1.);

  cx = 0;
  for (i=0; i<freeNode.n; i++)
    for (j=0; j<freeNode[i].dofX->n; j++)
      freeNode[i].x[j] += random[cx++] * maxX;

  random.free();

// diff test for stiffness with respect to intfU

      R.setDim(6*m);
  kdiff.setDim(m*m);
  kanly.setDim(m*m);

  j = 0;

  for (ic=0; ic<freeNode.n; ic++)
  {
    for (jc=0; jc<freeNode[ic].dofU->n; jc++) // loop over columns of stiffness
    {
      for (k=0; k<6; k++) // loop over perturbations
      {
        // apply pertubation

        freeNode[ic].u[jc] += dd[k];

        if (freeNode[ic].isLagr) freeNode[ic].x[jc] += dd[k] * dxdu;

        // calculate residual

        eliminate(false,true);

        // remove pertubation
  	
        freeNode[ic].u[jc] -= dd[k];

        if (freeNode[ic].isLagr) freeNode[ic].x[jc] -= dd[k] * dxdu;

        // store residual

        i = 0;

        for (ir=0; ir<freeNode.n; ir++)

          for (jr=0; jr<freeNode[ir].dofU->n; jr++)

            R[k+6*i++] = freeNode[ir].reac[jr];
      }

      for (i=0; i<nAllDoFU; i++)  // loop over rows of stiffness
  		
        kdiff[j*nAllDoFU+i] = (+       R[i*6+0]
                               -  9. * R[i*6+1]
                               + 45. * R[i*6+2]
                               - 45. * R[i*6+3]
                               +  9. * R[i*6+4]
                               -       R[i*6+5] ) / (60. * ddd);
      j++;
    }
  }

  // calculate stiffness

  eliminate(true,true);

  generateSimpleStiffnessMatrixU(&(kanly.x));

  prgCompareTwoSimpleMatrices(kdiff.x,                       // matrix 1
		              kanly.x,                       // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      nAllDoFU,nAllDoFU,             // matrix dimension
			      dig,dig2,gfrmt,                // format
			      0,                             // indentation
			      true,                          // interactive
			      true);                         // row/column numbers

  if (!isALE()) return;

// diff test for stiffness with respect to intfX

  kdiff.zero();

      Rx.setDim(6*neqx);
  kxdiff.setDim(m*neqx);
  kxanly.setDim(m*neqx);

  j = 0;

  for (ic=0; ic<freeNode.n; ic++)
  {
    for (jc=0; jc<freeNode[ic].dofX->n; jc++) // loop over columns of stiffness
    {
      for (k=0; k<6; k++) // loop over perturbations
      {
        // apply pertubation
  
        freeNode[ic].x[jc] += dd[k];
  
        // calculate residual
  
        eliminate(false,true);
  
        // remove pertubation
  	
        freeNode[ic].x[jc] -= dd[k];
  
        // store residual
  
        i = 0;
  
        for (ir=0; ir<freeNode.n; ir++)
  
          for (jr=0; jr<freeNode[ir].dofX->n; jr++)
  
            R[k+6*i++] = freeNode[ir].reac[jr];
  
        for (i=0; i<numnp*ndm; i++) { ii = idx.x[i]-1; if (ii > -1) Rx[ii*6+k] = x.x[i]; }
      }
  
      if (!freeNode[ic].isLagr)
      {
        for (i=0; i<nAllDoFU; i++)  // loop over rows of stiffness
  		
          kdiff[j*nAllDoFU+i] = (+       R[i*6+0]
                                 -  9. * R[i*6+1]
                                 + 45. * R[i*6+2]
                                 - 45. * R[i*6+3]
                                 +  9. * R[i*6+4]
                                 -       R[i*6+5] ) / (60. * ddd);
      }
  
      for (i=0; i<neqx; i++)  // loop over rows of stiffness
  
        kxdiff[j*neqx+i] = ( +       Rx[i*6+0]
                             -  9. * Rx[i*6+1]
                             + 45. * Rx[i*6+2]
                             - 45. * Rx[i*6+3]
                             +  9. * Rx[i*6+4]
                             -       Rx[i*6+5] ) / (60. * ddd);
      j++;
    }
  }

  // calculate stiffness

  eliminate(true,true);

  generateSimpleStiffnessMatrixXX(&(kxanly.x));

  prgCompareTwoSimpleMatrices(kxdiff.x,                      // matrix 1
		              kxanly.x,                      // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      neqx,nAllDoFX,                 // matrix dimension
			      dig,dig2,gfrmt,                // format
			      0,                             // indentation
			      true,                          // interactive
			      true);                         // row/column numbers

  generateSimpleStiffnessMatrixX(&(kanly.x));

  prgCompareTwoSimpleMatrices(kdiff.x,                       // matrix 1
		              kanly.x,                       // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      nAllDoFU,nAllDoFX,             // matrix dimension
			      dig,dig2,gfrmt,                // format
			      0,                             // indentation
			      true,                          // interactive
			      true);                         // row/column numbers
  return;
}






void FiniteElementBVPWI::generateSimpleStiffnessMatrixU(double **k)
{
  int d, i, j, l, m, o, nl, nm;

  for (i=0; i<nAllDoFU*nAllDoFU; i++) (*k)[i] = 0.;

  for (i=0; i<freeNode.n; i++)
  {
    nl = freeNode[i].dofU->n;
  
    for (d=0; d<freeNode[i].toFree.n; d++)
    {
      for (j=0; j<freeNode[i].toFree[d].n; j++)
      {
        o  = freeNode[i].toFree[d][j].nd;

        nm = freeNode[o].dofU->n;

        for (l=0; l<nl; l++)
        {
          for (m=0; m<nm; m++)
          {
            (*k)[(profU[o]+m)*nAllDoFU+profU[i]+l] = freeNode[i].toFree[d][j].dgfduf[m*nl+l];
          }
        }
      }
    }
  }

  return;
}





void FiniteElementBVPWI::generateSimpleStiffnessMatrixX(double **k)
{
  int d, i, j, l, m, o, nl, nm;

  for (i=0; i<nAllDoFU*nAllDoFX; i++) (*k)[i] = 0.;

  for (i=0; i<freeNode.n; i++)
  {
    nl = freeNode[i].dofU->n;
  
    for (d=0; d<freeNode[i].toFree.n; d++)
    {
      for (j=0; j<freeNode[i].toFree[d].n; j++)
      {
        o  = freeNode[i].toFree[d][j].nd;

        if (!freeNode[o].isLagr)
        {
          nm = freeNode[o].dofX->n;

          for (l=0; l<nl; l++)
          {
            for (m=0; m<nm; m++)
            {
              (*k)[(profX[o]+m)*nAllDoFU+profU[i]+l] = freeNode[i].toFree[d][j].dgfdxf[m*nl+l];
            }
          }
        }
      }
    }
  }

  return;
}





void FiniteElementBVPWI::generateSimpleStiffnessMatrixXX(double **k)
{
  int d, i, j, l, ll, m, o, nm;

  double fact;

  for (i=0; i<neqx*nAllDoFX; i++) (*k)[i] = 0.;

  for (i=0; i<freeNode.n; i++)
  {
    nm = freeNode[i].dofX->n;

    if (freeNode[i].isLagr) fact = 1./dxdu; else fact = 1.; 
 
    for (d=0; d<freeNode[i].toMesh.n; d++)
    {
      for (j=0; j<freeNode[i].toMesh[d].n; j++)
      {
        o  = freeNode[i].toMesh[d][j].nd;

        for (l=0; l<ndm; l++)
        {
          ll = idx.x[o*ndm+l] - 1;

          if (ll > -1)
          {
            for (m=0; m<nm; m++)
            {
              (*k)[(profX[i]+m)*neqx+ll] = freeNode[i].toMesh[d][j].xx[m*ndm+l] * fact;
            }
          }
        }
      }
    }
  }

  return;
}








void FiniteElementBVPWI::prepareForExternalSolver(void *interfaceSolver,
                                                  double *rMesh,
                                                  bool printMeshRes)
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::
// This member function is tailored to be evoked from a 
// child of InterfaceMatch, specifically InterfaceN.
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::
{
  char fct[] = "FiniteElementBVPWI::prepareForExternalSolver";

  //cout << fct << "\n";

  computerTime.go(fct);

  int e, i, j, k, n;

  //cout << firstIter << "\n";

  bool meshFlag = (isALE() && neqx > 0), 
       elimMeshFlag = (rMesh == NULL);

  double fact;

  for (i=0; i<freeNode.n; i++) freeNode[i].zeroDerivativeMtx();

  // update uDep to provide correct boundary conditions

  for (i=0; i<bndNode.n; i++) 
  {
    bndNode[i].getUX(ndf,ndm,u.x,x.x,xn.x,x0.x);

    for (j=0; j<bndNode[i].uDep.n; j++)
    {
      k = bndNode[i].uDep[j]->dof- 1;
      n = bndNode[i].uDep[j]->nd - 1;

      u.x[n*ndf+k] = bndNode[i].uDep[j]->uc;
    }

    if (meshFlag && !elimMeshFlag)

      for (j=0; j<bndNode[i].xDep.n; j++)
      {
        k = bndNode[i].xDep[j]->dof- 1;
        n = bndNode[i].xDep[j]->nd - 1;

        d.x[n*ndm+k] = bndNode[i].xDep[j]->uc;

        x.x[n*ndm+k] = xn.x[n*ndm+k] + d.x[n*ndm+k];
        v.x[n*ndm+k] = td[22] * vn.x[n*ndm+k] + td[20] * d.x[n*ndm+k];
      }
  }

  updateIterStep();

  // update mesh and calculate mesh derivative

  if (meshFlag)
  {
    if (elimMeshFlag)
    {
      updateMesh(1, printMeshRes); if (printMeshRes) cout << "\n";

      calculateDerivatives4();
      calculateDerivatives5();
    }
    else
    {
      solverMesh = (Solver*)interfaceSolver;

      solverMesh->whatToSolveFor = MSH;

      solverMesh->dom = (Domain*)this;

      r.zero();

      for (e=0; e<numel; e++)
      {
        if (elem[e]->calcStiffnessAndResidualMesh() != 0)

          prgWarning(1,fct,"Error in calcStiffnessAndResidualMesh");

        solverMesh->assembleElemMtx(elem[e]);

        solverMesh->assembleElemVec(elem[e],firstIter);
      }
      solverMesh = NULL;

      for (i=0; i<neqx; i++) rMesh[i] = r[i];

      calculateDerivatives4();
    }
  }

  localStiffnessError = 0;

  if (solver != NULL) prgError(1,fct,"solver != NULL");

  solver = (Solver*)interfaceSolver;

  solver->whatToSolveFor = DOF;

  solver->dom = (Domain*)this;

  if (calcStiffnessAndResidual(false,false) != 0)

    prgError(3,fct,"problem in calcStiffnessAndResidual!");

  solver = NULL;

  // calculate reaction forces in free nodes

  for (i=0; i<freeNode.n; i++) freeNode[i].reac.zero();

  for (i=0; i<bndNode.n; i++) bndNode[i].giveReactions(ndf,ndm,reac.x,u.x,x.x,x0.x);

  // calculate derivatives in layer adjacent to the interface

  calculateDerivatives1();

  calculateDerivatives2();

  if (meshFlag && elimMeshFlag) calculateDerivatives6();

  calculateDerivatives7();

  ctimPrepExt += computerTime.stop(fct);

  donePrepExt = true;

  return;
}










void FiniteElementBVPWI::printComputerTime(bool reset, int detailFlg)
{
  FiniteElementBVP::printComputerTime(reset,detailFlg);

  if (doneElim || donePrepExt || doneDeriv1 || doneDeriv2 || doneDeriv3 || doneDeriv4
                              || doneDeriv5 || doneDeriv6 || doneDeriv7 || doneDeriv8)

    COUT << "----------------------------------------------------\n";

  if (doneElim) { COUT;
        printf("FiniteElementBVPWI::eliminate: %7.3f sec ->%5.1f %\n",
               ctimElim, ctimElim/ctimSinceLastCall*100.); }

  if (donePrepExt) { COUT;
        printf("FiniteElementBVPWI::prepExtSlv:%7.3f sec ->%5.1f %\n",
               ctimPrepExt, ctimPrepExt/ctimSinceLastCall*100.); }

  if (doneDeriv1) { COUT;
        printf("FiniteElementBVPWI::calcDeriv1:%7.3f sec ->%5.1f %\n",
               ctimDeriv1, ctimDeriv1/ctimSinceLastCall*100.); }

  if (doneDeriv2) { COUT;
        printf("FiniteElementBVPWI::calcDeriv2:%7.3f sec ->%5.1f %\n",
               ctimDeriv2, ctimDeriv2/ctimSinceLastCall*100.); }

  if (doneDeriv3) { COUT;
        printf("FiniteElementBVPWI::calcDeriv3:%7.3f sec ->%5.1f %\n",
               ctimDeriv3, ctimDeriv3/ctimSinceLastCall*100.); }

  if (doneDeriv4) { COUT;
        printf("FiniteElementBVPWI::calcDeriv4:%7.3f sec ->%5.1f %\n",
               ctimDeriv4, ctimDeriv4/ctimSinceLastCall*100.); }

  if (doneDeriv5) { COUT;
        printf("FiniteElementBVPWI::calcDeriv5:%7.3f sec ->%5.1f %\n",
               ctimDeriv5, ctimDeriv5/ctimSinceLastCall*100.); }

  if (doneDeriv6) { COUT;
        printf("FiniteElementBVPWI::calcDeriv6:%7.3f sec ->%5.1f %\n",
               ctimDeriv6, ctimDeriv6/ctimSinceLastCall*100.); }

  if (doneDeriv7) { COUT;
        printf("FiniteElementBVPWI::calcDeriv7:%7.3f sec ->%5.1f %\n",
               ctimDeriv7, ctimDeriv7/ctimSinceLastCall*100.); }

  if (doneDeriv8) { COUT;
        printf("FiniteElementBVPWI::calcDeriv8:%7.3f sec ->%5.1f %\n",
               ctimDeriv8, ctimDeriv8/ctimSinceLastCall*100.); }

  if (doneDeriv3) { COUT;
        printf("FiniteElementBVPWI/backSubst3: %7.3f sec ->%5.1f %\n",
               ctimBack3, ctimBack3/ctimSinceLastCall*100.); }

  if (doneDeriv5) { COUT;
        printf("FiniteElementBVPWI/backSubst5: %7.3f sec ->%5.1f %\n",
               ctimBack5, ctimBack5/ctimSinceLastCall*100.); }

  if (doneDeriv8) { COUT;
        printf("FiniteElementBVPWI/backSubst8: %7.3f sec ->%5.1f %\n",
               ctimBack8, ctimBack8/ctimSinceLastCall*100.); }

  if (reset)
  {
    ctimDeriv1 = 0.;  doneDeriv1 = false;
    ctimDeriv2 = 0.;  doneDeriv2 = false;
    ctimDeriv3 = 0.;  doneDeriv3 = false;
    ctimDeriv4 = 0.;  doneDeriv4 = false;
    ctimDeriv5 = 0.;  doneDeriv5 = false;
    ctimDeriv6 = 0.;  doneDeriv6 = false;
    ctimDeriv7 = 0.;  doneDeriv7 = false;
    ctimDeriv8 = 0.;  doneDeriv8 = false;

    ctimBack3  = 0.;
    ctimBack5  = 0.;
    ctimBack8  = 0.;

    ctimElim   = 0.;  doneElim   = false;
    ctimPrepExt= 0.;  donePrepExt= false;
  }

  return;
}

