
#include <iostream>

#include "SolverSparse.h"
#include "Mesh.h"
#include "Debug.h"
#include "MathMatrix.h"
#include "FunctionsProgram.h"
#include "ComputerTime.h"
#include "DependentDoF.h"


extern ComputerTime computerTime;


using namespace std;



SolverSparse::SolverSparse(void)
{
  if (debug) cout << " SolverSparse constructor\n\n";

  comprMtxFlg = false;

  return;
}




SolverSparse::~SolverSparse()
{
  free();
	
  if (debug)  cout << " SolverSparse destructor\n\n";

  return;
}





void SolverSparse::prepareMatrixPattern(Domain *domain)
{
  char fct[] = "SolverSparse::prepareMatrixPattern";
  
  if (currentStatus != EMPTY)
    { prgWarning(1,fct,"solver is not EMPTY! but I go on ..."); }
    
  computerTime.go(fct);
 
  int ii1, ii2, i1, i2, j1, j2, k1, k2, e, ee, e1, flg,
      r, c, nen, nen1, ndf, ndf1, ndfMin, nst, nst1, nel,
      *sparse, *sparse1, *id1, *id2, *nodeElem, *ix1, 
      numel, numnp, neq;

  Element *elem, *elem1;

  MatrixSparse<double> mtxs;

  // get stuff from the mesh .....................................................
  
  dom = domain;

  if (dom == NULL) prgError(1,fct,"invalid domain pointer!");

  Mesh *mesh = (Mesh*) dom; 
 
  flg = DOF; if (this == mesh->solverMesh) flg = MSH;

  whatToSolveFor = flg;
  
  Element **elemList = mesh->elem;

  VectorArray<int> *nodeElement = mesh->nodeElem;
  
  MatrixFullArray<int> *ID = &(mesh->idu); if (flg == MSH) ID = &(mesh->idx);
  MatrixFullArray<int> &id = *ID;

  List<DependentDoF> *UD = &(mesh->uDep); if (flg == MSH) UD = &(mesh->xDep);
  List<DependentDoF> &uDep = *UD;
  
  numel = mesh->numel;
  numnp = mesh->numnp;
  
  neq   = mesh->nequ; if (flg == MSH) neq = mesh->neqx;

  // calculate standard connectivities ...........................................
  
  for (e=0; e<numel; e++)
  {
    elem = elemList[e];
	  
    ndf = elem->ndfm(flg);
    
    nst  = elem->nen() * ndf;

    nst1 = nst * nst;
    
    if (elem->sparse[flg] != NULL) delete [] elem->sparse[flg];
	  
    elem->sparse[flg] = new int [nst1];

    for (r=0; r<nst1; r++) elem->sparse[flg][r] = -1;
  }
  
//  cout << " numnp " << numnp << endl;
  
  for (ii1=0; ii1<numnp; ii1++)
  {
    //cout << ii1 << "\n";

    nel      = nodeElement[ii1].n;

    if (nel == 0) goto jump1;
    
    nodeElem = &(nodeElement[ii1][0]);
    
//    cout << " nel " << nel << " ......... " ;
//    for(e=0;e<nel;e++)
//      cout << '\t' << nodeElem[e];
//    cout << endl;
    
    id1 = &(id(ii1+1,1));
    
    for (e=0; e<nel; e++)
    {
      ee     = nodeElem[e];

      elem   = elemList[ee];
      
      nen    = elem->nen();
      
      ndf    = elem->ndfm(flg);

      nst    = nen * ndf;
      
      sparse = elem->sparse[flg];
    
      i1 = 0; while (elem->ix[i1] != ii1+1) i1++;
      
//      cout << '\t' << " i1 ... : " << i1 << endl;
      
      for (i2=0; i2<nen; i2++)
      {
        ii2 = elem->ix[i2];

        id2 = &(id(ii2,1));
	
	k2 = -1;
	
        e1 = 0;

//        cout << '\t' << " elem1 ..... " ;

	while (e1 < e)
	{
          elem1 = elemList[nodeElem[e1]];
          
//          cout << '\t' << nodeElem[e1];
		
	  nen1  = elem1->nen();
		  
          k2++;

          ix1 = elem1->ix;
	  
	  while (k2 < nen1 && ix1[k2] != ii2) k2++;
	  
	  if (k2 == nen1) { e1++; k2 = -1; } else break;
	}

//        cout << endl;	
          
        if (k2 > -1)
	{

          //cout << " exists!  " << ee+1 << "," << e1+1 << ": " << ii1+1 << "," << ii2 << "\n";

	  ndf1    = elem1->ndfm(flg);
	  
          nst1    = ndf1 * nen1;
	  
          sparse1 = elem1->sparse[flg];
	  
          k1 = 0; while (ix1[k1] != ii1+1) k1++;
	  
          if (ndf1 > ndf) ndfMin = ndf; else ndfMin = ndf1;
	  
          for (j1=0; j1<ndfMin; j1++)
	  		  
	    for (j2=0; j2<ndfMin; j2++)

              sparse[nst*(i2*ndf+j2)+i1*ndf+j1] = sparse1[nst1*(k2*ndf1+j2)+k1*ndf1+j1];
	    
          // if  ndf > ndf1  generate remaining matrix coefficients
	  
	  for (j1=ndfMin; j1<ndf; j1++)
	  {
            r = id1[j1];
	    
            if (r > 0)
	    {
	      for (j2=0; j2<ndf; j2++)
	      {
		c = id2[j2];

		if (c > 0)  sparse[nst*(i2*ndf+j2)+i1*ndf+j1] = mtxs.append(r,c); 
	      }
	    }

	    c = id2[j1];

	    if (c > 0)
	    {
	      for (j2=0; j2<ndfMin; j2++)
	      {
	        r = id1[j2];

		if (r > 0)  sparse[nst*(i2*ndf+j1)+i1*ndf+j2] = mtxs.append(r,c); 
	      }
	    }	        
	  }
	  // .....................................................
	}
	else
	{
          //cout << " new!     " << ee+1 << "  : " << ii1+1 << "," << ii2 << "\n";

          for (j1=0; j1<ndf; j1++)
	  {		  
            r = id1[j1];

	    if (r > 0)
	    {
	      for (j2=0; j2<ndf; j2++)
	      {
                c = id2[j2];

		if (c > 0)  sparse[nst*(i2*ndf+j2)+i1*ndf+j1] = mtxs.append(r,c); 
	      }
	    }
	  }
	}
      }
    }
    jump1: continue;
  }

  // account for dependent degrees of freedom .........................................

  int m, i, j, k, hr, hc, ne = mtxs.x.n, pos, mstr, mstr2;

  List< Vector<int> > *sMasterRowAll = new List< Vector<int> > [numel], 
                      *sMasterColAll = new List< Vector<int> > [numel],
                      *sMasterRowElm, *sMasterColElm;

  ListArray< List< Vector<int> > >    elemSMaster;
  ListArray< List< Vector<double> > > elemSFact;
  
  elemSMaster.setDim(numel);
  elemSFact.setDim(numel);

  Vector<int>    *sm, *sr, *sc;
  Vector<double> *sf;
  
  DependentDoF *uDep1, *uDep2;

  for (e1=0; e1<numel; e1++)
  {
    elem1  = elemList[e1];
      
    ix1    = elem1->ix;

    nen    = elem1->nen();
    
    ndf    = elem1->ndfm(flg);
      
    nst    = nen * ndf;
    
    sparse = elem1->sparse[flg];
    
    m = 0;

    sMasterRowElm = &(sMasterRowAll[e1]);
    sMasterColElm = &(sMasterColAll[e1]);
    
    for (i1=0; i1<nen; i1++)
    {
      ii1 = ix1[i1];

      id1 = &(id(ii1,1));
      
      for (j1=0; j1<ndf; j1++)
      {
        r = id1[j1];

	if (r < 0)
	{
          uDep1 = &(uDep[-1-r]);

	  if (uDep1->master.n > 0)
	  {
	    j = 0; while (uDep1->master[j] < 1 && j < uDep1->master.n) j++;

            if (j < uDep1->master.n)
	    {	
              // degree of freedom j1 of node i1/ii1 of element e1
	      // is a dependent with 'real' masters
		    
              hr = i1 * ndf + j1;
	  
              hc = (i1*ndf+j1)*nst;
	  
              for (i2=0; i2<nen; i2++)
	      {
	        ii2 = ix1[i2];
	      
	        id2 = &(id(ii2,1));
            
	        for (j2=0; j2<ndf; j2++)
	        {
	          c = id2[j2];	  
	  
	          if (c != 0)
	          { 
                    if (c < 0)
		    {			    
                      uDep2 = &(uDep[-1-c]);

		      if (uDep2->master.n == 0) goto jump2;

                      j = 0; while (uDep2->master[j] < 1 && j < uDep2->master.n) j++;

                      if (j == uDep2->master.n) goto jump2;
		    }
		
                    // degree of freedom j2 of node i2/ii2 of element e1 
		    // is free or is a dependent of 'real' masters
		    
                    j = i2 * ndf + j2;
	        
	            if (sparse[j*nst+hr] == -1)
	            {
                      elemSMaster[e1].add(new Vector<int>);
		      sMasterRowElm->add (new Vector<int>);
		      sMasterColElm->add (new Vector<int>);
                      elemSFact[e1].add  (new Vector<double>);
                      sparse[j*nst+hr] = - 2 - m++;
	            }
	            i  = -2-sparse[j*nst+hr];
                    sm = &(elemSMaster[e1][i]);
		    sr = &((*sMasterRowElm)[i]);
		    sc = &((*sMasterColElm)[i]);
                    sf = &(elemSFact[e1][i]);
	          
                    j = sm->n; 
	          
                    if (c > 0)
	            { 
	              for (i=0; i<uDep1->master.n; i++)
	              { 
			mstr = uDep1->master[i];
			if (mstr > 0)
                        {
                          if (uDepPosExists(elemList, nodeElement, elemSMaster,
					    sMasterRowAll,sMasterColAll,
				            uDep1,i,mstr,c,true,&pos))
			    //(*sm)[j] = pos;
			    sm->append(pos);
			  else
                          {
	                    //(*sm)[j] = mtxs.append(mstr,c);
	                    sm->append(mtxs.append(mstr,c));
                          }
			  //(*sr)[j] = mstr;
			  //(*sc)[j] = c;
			  
	                  //(*sf)[j++] = uDep1->beta[i];

                          sr->append(mstr);
                          sc->append(c);
                          sf->append(uDep1->beta[i]);
			}
	              }
	            }
		    else
	            {
	              for (i=0; i<uDep1->master.n; i++)
	              {
			mstr = uDep1->master[i];
                        if (mstr > 0)
			{
	                  for (k=0; k<uDep2->master.n; k++)
  	                  {
		            mstr2 = uDep2->master[k];
			    if (mstr2 > 0)
			    {
			      if (uDepPosExists(elemList, nodeElement, elemSMaster,
					        sMasterRowAll,sMasterColAll,
						uDep1,i,mstr,mstr2,true,&pos))
				//(*sm)[j] = pos;
				sm->append(pos);
			      else
			      {
			        //(*sm)[j] = mtxs.append(mstr,mstr2);
                                sm->append(mtxs.append(mstr,mstr2));
			      }
                              //(*sr)[j] = mstr;
			      //(*sc)[j] = mstr2;

	                      //(*sf)[j++] = uDep1->beta[i] * uDep2->alpha[k];
	                      sr->append(mstr);
                              sc->append(mstr2);
                              sf->append(uDep1->beta[i] * uDep2->alpha[k]);
			    }
	                  }
			}
	              }
	            }
		    
                    j = i2 * ndf + j2;
  	      
  	            if (sparse[j+hc] == -1)
  	            {
                      elemSMaster[e1].add(new Vector<int>);
		      sMasterRowElm->add (new Vector<int>);
		      sMasterColElm->add (new Vector<int>);
                      elemSFact[e1].add  (new Vector<double>);
                      sparse[j+hc] = - 2 - m++;
  	            }
                    if (c > 0) // within this if-clause, 'c' stands for 'row'
  	            {
		      i  = -2-sparse[j+hc];
                      sm = &(elemSMaster[e1][i]);
                      sr = &((*sMasterRowElm)[i]);
		      sc = &((*sMasterColElm)[i]);
                      sf = &(elemSFact[e1][i]);
               
	              j = sm->n; 
	          
	              for (i=0; i<uDep1->master.n; i++)
	              {		
			mstr = uDep1->master[i];
                        if (mstr > 0) 
			{
		          if (uDepPosExists(elemList, nodeElement, elemSMaster,
					    sMasterRowAll,sMasterColAll,
				            uDep1,i,c,mstr,false,&pos))
			    //(*sm)[j] = pos;
			    sm->append(pos);
			  else
			  {
                            //(*sm)[j] = mtxs.append(c,mstr);
                            sm->append(mtxs.append(c,mstr));
			  }
                          //(*sr)[j] = c;
			  //(*sc)[j] = mstr;
			  
	                  //(*sf)[j++] = uDep1->alpha[i];
	                  sr->append(c);
                          sc->append(mstr);
                          sf->append(uDep1->alpha[i]);
			}
		      }
		    }
		  }
                  jump2: continue;
	        }
	      }
	    }
	  }
	}
      }
    }
  }

  delete [] sMasterRowAll;
  delete [] sMasterColAll;

  for (e=0; e<numel; e++)
  {
    elem = elemList[e];
    elem->sMaster[flg] = elemSMaster[e];
    elem->sFact[flg] = elemSFact[e];
  }

  VectorArray<int> perm, tmp;

  VectorArray<int> *sMaster;

  if (comprMtxFlg)
  {
    // sort matrix entries

    COUT << "sorting matrix entries for compressed row format ...\n\n";

    computerTime.go("sorting matrix entries");

    mtxs.quickSort(true,true,perm);

    computerTime.stopAndPrint("sorting matrix entries");

    // check for double entries

    for (i=1; i<mtxs.row.n; i++)
    {
      if (mtxs.row[i] == mtxs.row[i-1])
        if (mtxs.col[i] == mtxs.col[i-1])
          prgError(100,fct,"double entries! Could the solver cope?!");
    }

    // adjust elem->sparse and elem->sMaster

    tmp.setDim(perm.n);

    for (i=0; i<perm.n; i++) tmp[perm[i]] = i;

    for (e=0; e<numel; e++)
    {
      elem = elemList[e];
    
      nst  = elem->nen() * elem->ndfm(flg);

      nst1 = nst * nst;

      for (i=0; i<nst1; i++)
        if (elem->sparse[flg][i] > -1)
          elem->sparse[flg][i] = tmp.x[elem->sparse[flg][i]];

      sMaster = elem->sMaster[flg].x;
      k       = elem->sMaster[flg].n;
      for (i=0; i<k; i++)
        for (j=0; j<sMaster[i].n; j++)
          sMaster[i][j] = tmp.x[sMaster[i][j]];
    }

    perm.free();
    tmp.free();

    // generate compr

    compr.setDim(mtxs.nRow+1);
    j = 0;
    compr[j++] = 1;
    for (i=1; i<mtxs.row.n; i++)
      if (mtxs.row[i] != mtxs.row[i-1]) compr[j++] = i+1;

    compr[j] = mtxs.col.n+1;

    if (j != mtxs.nRow) prgError(200,fct,"fatal error!");

    //cout << mtxs.col << "\n\n" << mtxs.row << "\n\n" << compr << "\n\n";

    COUT << "matrix compression done.\n\n";
  }

  mtx = mtxs;

  mtxs.free();

/*  for (e=0; e<dom->numel; e++)
  {  
    sparse = dom->elem[e]->sparse;

    nst = dom->elem[e]->ndf() * dom->elem[e]->nen();
	  
    for (i=0; i<nst; i++) 
      { for (j=0; j<nst; j++) printf("%5d",sparse[j*nst+i]); cout << "\n"; }

    cout << "\n";

    for (i=0; i<dom->elem[e]->sMaster[flg].n; i++)

	    cout << -2-i << ": " << dom->elem[e]->sMaster[flg][i] 
		              << dom->elem[e]->sFact[flg][i] << "\n";

    cout << "\n\n";
    
  }
  cout << " 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 \n";
  cout << mtx.msl.row << "\n" << mtx.msl.col << "\n"; */

  if (neq != mtx.nRow) 
  {
    cout << neq << " != " << mtx.nRow << "\n";
    prgError(2,fct,"neq inconsistency! undefined dof ?");
  }

  if (ne != mtx.x.n) 
    COUT << "number of coefficients due to dependencies: " << mtx.x.n - ne << "\n\n";

  currentStatus = PATTERN_OK;
  
/*  
    cout << endl;
    cout << endl;    
  for (e=0; e<numel; e++)
  {
    elem = elemList[e];
	  
    ndf = elem->ndfm(flg);
    nst  = elem->nen() * ndf;

    i1=0;
    for(r=0; r<nst; r++)
    {
      for(c=0; c<nst; c++)
      {
         cout << '\t' << elem->sparse[flg][i1] ;
         i1++;
      }
      cout << endl;
    }
    cout << endl;
    cout << endl;    
  }
*/    

  computerTime.stopAndPrint(fct);

  return;
}










bool SolverSparse::uDepPosExists(Element **elemList,
		                 VectorArray<int> *nodeElement,
                                 ListArray< List< Vector<int> > > &elemSMaster,
		                 List< Vector<int> > *sMRA, 
		                 List< Vector<int> > *sMCA,
				 DependentDoF* uDep, int mi, int r, int c, 
		                 bool rowFlg, int* pos)
{
  int e, ee, ii1, i1, j1, l1, k, l, sp, ndf, nst, *sparse, *nodeElem, nel, flg = whatToSolveFor;

  Element *elem;

  List< Vector<int> > *sMRE, *sMCE;
  
  Vector<int> *sm, *sr, *sc;
      
  ii1 = uDep->nd - 1;
  j1  = uDep->dof - 1;
	
  nel = nodeElement[ii1].n;

  if (nel == 0) return false;
    
  nodeElem = &(nodeElement[ii1][0]);
  
  for (e=0; e<nel; e++)
  {
    ee     = nodeElem[e];

    elem   = elemList[ee];
     
    ndf    = elem->ndfm(flg);
    
    nst    = elem->nen() * ndf;
      
    sparse = elem->sparse[flg];
    
    i1 = 0; while (elem->ix[i1] != ii1+1) i1++;

    l1 = i1*ndf + j1;
    
    sMRE = &(sMRA[ee]);
    sMCE = &(sMCA[ee]);
    
    for (k=0; k<nst; k++)
    {
      if (rowFlg) sp = sparse[k*nst+l1];
      else        sp = sparse[l1*nst+k];

      if (sp < -1)
      {
        sm = &(elemSMaster[ee][-2-sp]);
	sr = &((*sMRE)      [-2-sp]);
	sc = &((*sMCE)      [-2-sp]);

        for (l=0; l<sm->n; l++)
	{
          if ((*sr)[l] == r && (*sc)[l] == c) { *pos = (*sm)[l]; return true; }
	}
      }
    }
  }
  
  return false;	
}









void SolverSparse::zeroMtx(void)

{
  mtx.zero();

  return;
}






	
int SolverSparse::assembleElemMtx(Element* elem)
{
  int    flg = whatToSolveFor, nst2 = elem->ndfm(flg) * elem->nen(),
         *sp  = elem->sparse[flg], i, d, m;

  double *s = dom->s;

  //prgPrintSimpleMatrix(s,6,6,12,5,true,5,false);
 
  nst2 *= nst2;
  
  for (i=0; i<nst2; i++)    
  {
    //if (sp[i] > -1)   { cout << sp[i] << "<" << mtx.x.n << "\n";     mtx[sp[i]] += s[i]; }

    if (sp[i] > -1)    mtx[sp[i]] += s[i];
  
    else if (sp[i] < -1)   
    {
       d = - 2 - sp[i];
       
       for (m=0; m<elem->sMaster[flg][d].n; m++)
       
         mtx[elem->sMaster[flg][d][m]] += elem->sFact[flg][d][m] * s[i];
    }
  }

  return 0;
}
    




int SolverSparse::assembleElemVec(Element* elem, bool firstIter)
{
  int    flg     = whatToSolveFor, 
         nen     = elem->nen(),
         ndf     = elem->ndfm(flg),
	 nst     = ndf * nen,
	 *ix     = elem->ix,
	 indfDom = 0, 
	 indf    = 0,
	 i, ii, j, jj, m, n, k, *l = new int[nst];
 
  // get stuff from the mesh
  
  Mesh *mesh = (Mesh*) dom; 
 
  List<DependentDoF> *UD = &(mesh->uDep); 
  int   *ID = mesh->idu.x,
         ndfDom  = mesh->ndf;
  double *s = mesh->s,
	 *p = mesh->p,
      *reac = mesh->reac.x,
         *r = mesh->r.x;
  
  if (flg == MSH) 
  {
    UD = &(mesh->xDep);
    ID = mesh->idx.x;  
    reac = mesh->reacMesh.x;
    ndfDom = dom->ndm;
  }
 
  List<DependentDoF> &uDep = *UD;
 
  DependentDoF *uDep1, *uDep2;
 
  // get local idu and add up reactions
  
  for (i=0; i<nen; i++)  
  { 
    indfDom = (ix[i]-1) * ndfDom;
    for (j=0; j<ndf; j++)  
    {
      l[indf+j]       = ID[indfDom+j] - 1;
      reac[indfDom+j] += p[indf+j];
    }
    indf += ndf;
  } 

  // add up residuals

  //cout << '\t' << "  nst  : " << nst << endl;
 
  for (i=0; i<nst; i++) 
  {	  
    ii = l[i];
	  
    if (ii > -1)  r[ii] += p[i];    // standard residuals

    else if (ii < -1) 
    {
      uDep1 = &(uDep[-2-ii]);
      
      if (firstIter)  // eliminate prescribed displacements
      { 
        for (j=0; j<nst; j++)
	{
          jj = l[j];
          
          cout << uDep1->duc << endl;
	 
	  if (jj > -1)  r[jj] -= s[i*nst+j] * uDep1->duc; 

	  else if (jj < -1)
	  {
            uDep2 = &(uDep[-2-jj]);

	    for (k=0; k<uDep2->master.n; k++)
	    {
              m = uDep2->master[k] - 1;
	      
	      if (m > -1)  r[m] -= s[i*nst+j] * uDep1->duc * uDep2->beta[k];
	    }
	  }
	}
      }

      for (k=0; k<uDep1->master.n; k++)  // add up equations according dependent dof
      {
        m = uDep1->master[k] - 1;

        if (m > -1)  r[m] += uDep1->beta[k] * p[i];
      }      
    }
  }

  delete [] l;
  
  return 0;
}






void SolverSparse::free(void)
{
  mtx.free();
	
  return;
}
    






void SolverSparse::printInfo(void)
{
  COUT << "sparse solver:  neq = " << mtx.nRow << "\n";
  COUT << "                ne  = " << mtx.x.n << "\n\n"; 

  return;
}







void SolverSparse::printMatrix(int dig, int dig2, bool gfrmt, int indent, bool interactive)
{
  int i, j, dim = mtx.nRow;

  char frmt[30];
  
  if (gfrmt) sprintf(frmt," %%%d.%dg",dig,dig2);
  else       sprintf(frmt," %%%d.%df",dig,dig2);
 
  for (i=0; i<dim; i++)
  { 
    for (j=0; j<indent; j++) printf(" ");
    for (j=0; j<dim; j++)  printf(frmt,mtx(i+1,j+1,true)); printf("\n");
  }
  printf("\n");
	
  return;
}
    





double SolverSparse::giveMatrixCoefficient(int row, int col)
{ 
  return mtx(row,col,true);
}




void SolverSparse::copyToSimpleMatrix(double* &sMtx)
{
  mtx.copyToSimpleMatrix(sMtx);

  return;
}

