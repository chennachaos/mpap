
#include <iostream>

#include "FunctionsProgram.h"
#include "InterfaceN.h"
#include "DomainTree.h"
#include "SolverMA41.h"
//#include "SolverPARDISO.h"
#include "Solid.h"
#include "Fluid.h"
#include "FreeSurface.h"
#include "ComputerTime.h"
#include "MpapTime.h"


extern DomainTree   domain;
extern ComputerTime computerTime;
extern MpapTime     mpapTime;




InterfaceN::InterfaceN(void)                       
{                                                  
  ctimCalcSub1     = 0.;  doneCalcSub1     = false;
  ctimCalcSub2     = 0.;  doneCalcSub2     = false;
  ctimCalcSub3     = 0.;  doneCalcSub3     = false;
  ctimCalcFreeSurf = 0.;  doneCalcFreeSurf = false;

  ctimFactSolvUpdt = 0.;
  ctimCalcStiffRes = 0.;

  symFlag = false;

  // add new type
  
  DomainType *interfaceN = domain.newType(INTERFACEN,INTERFACEMATCH);
  
  if (interfaceN == NULL) return;  // domain type already exists

  //interfaceN->key.addNew("");

  return;
}




	                                          
InterfaceN::~InterfaceN(void)                     
{                                                 
  return;
}







void InterfaceN::readInputData(std::ifstream &Ifile, MyString &line)
{
  char fct[] = "InterfaceN::readInputData";

  switch (domain[INTERFACEN].key.whichBegins(line))
  {
    case -1: // go and inherit from INTERFACEMATCH
	     
	     InterfaceMatch::readInputData(Ifile,line); 
	     
	     break;
  }
 
  return;
}








void InterfaceN::prepareInputData(void)
{
  char fct[] = "InterfaceN::prepareInputData";

  // go and inherit from ancestors

  InterfaceMatch::prepareInputData();
 
  
  std::cout << "     INTERFACEN: prepare input data ...\n\n";

  // check that all domMode values admissible

  int admissibleDomMode[] = { ELIMINATE, ELIMINATE_MESH, RESOLVE, DUMMY };

  for (int dm=0; dm<domMode.n; dm++)
  {
    if (!isElementOf(domMode[dm], admissibleDomMode, 4))

      prgError(1,fct,"invalid domain mode!");
  }

  return;
}









void InterfaceN::prepareInteractions(void)
{
  // go and inherit from ancestors

  InterfaceMatch::prepareInteractions();

  cout << "     INTERFACEN: preparing interactions ...\n\n"; 

  char fct[] = "InterfaceN::prepareInteractions";
 
  return;
}









void InterfaceN::printInfo(void)
{ 
  COUT << "ndm ....... = " << ndm   << "\n";
  COUT << "ndf ....... = " << ndf   << "\n";
  COUT << "numnp ..... = " << numnp << "\n";
  COUT << "nequ        = " << nequ  << "\n";
  COUT << "neqx        = " << neqx  << "\n";
  COUT << "neqTot      = " << neqTot << "\n";
  COUT << "ne (global) = " << ((SolverSparse*)solver)->mtx.x.n << "\n";

  if (uType == VELOCITY) COUT << "uType       = VELOCITY\n";
  else                   COUT << "uType       = DISPLACEMENT\n";

  if (initGuess == KEEP_VELOCITY) COUT << "initGuess   = KEEP_VELOCITY\n";
  else                            COUT << "initGuess   = KEEP_DISPLACEMENT\n";

  cout << "\n";

  return;
}








void InterfaceN::setSolver(int slv, int *parm, bool cIO)
{
  char fct[] = "InterfaceN::setSolver";

  int numProc;

  symFlag = false;

  switch (slv) 
  {
    case  1: // MA41 ..........................

             solver = (Solver*) new SolverMA41;	     

             prepareMatrixPattern();

             solver->printInfo();

             if (solver->initialise() != 0) return; 

             break;
/*	     
    case  2: // PARDISO(sym)

             symFlag = true;

    case  3: // PARDISO(unsym)

#ifndef WINDOWS

             if (parm != NULL) numProc = parm[0]; else numProc = 1;

             numProc = min(MAX_PROCESSORS,numProc);

             solver = (Solver*) new SolverPARDISO;

             prepareMatrixPattern();

             solver->printInfo();

             if (slv == 2) if (solver->initialise(numProc,PARDISO_STRUCT_SYM) != 0) return;
             if (slv == 3) if (solver->initialise(numProc,PARDISO_UNSYM) != 0) return;

             break;

#else
             prgError(1,fct,"sparse solver PARDISO not yet available under windows!");
#endif
*/
    default: // invalid slv ...................
	     
	     cout << " this solver has not been implemented yet!\n\n"; 
	     
	     break;
  }

  solverOK = true;
 
  if (solver != NULL) solver->checkIO = cIO;

  return;
}












void InterfaceN::prepareMatrixPattern(void)
{
  char fct[] = "InterfaceN::prepareMatrixPattern";

  int dm, i, j, k, ku, kx, n;

  MatrixSparse<double> mtx;

  VectorArray<int> nidu, pidu, nidx, pidx;

  computerTime.go(fct);

  // calculate offsets for degree of freedom numbering of individual domains

  nOff.setDim(domType.n);

  nOff.zero();

  nOff[0] = nequ + neqx + 1;

  for (dm=1; dm<domType.n; dm++)
  {
    nOff[dm] = nOff[dm-1];

    if (domMode[dm-1] != ELIMINATE)
    {
      nOff[dm] += domPtr[dm-1]->nequ;

      if (domMode[dm-1] == RESOLVE && domPtr[dm-1]->isALE(CHECKNEQX))

        nOff[dm] += domPtr[dm-1]->neqx;
    }
  }

  // calculate neqTot and dimension residual vector

  dm = domType.n - 1;

  neqTot = nOff[dm] - 1;

  if (domMode[dm] != ELIMINATE)
  {
    neqTot += domPtr[dm]->nequ;

    if (domMode[dm]==RESOLVE && domPtr[dm]->isALE(CHECKNEQX)) neqTot += domPtr[dm]->neqx;
  }

  r.setDim(neqTot);

  // generate nidu, pidu and nidx, pidx

  nidu.setDim(numnp);
  pidu.setDim(numnp);
  nidx.setDim(numnp);
  pidx.setDim(numnp);

  pidu[0] = 1;
  pidx[0] = 1+nequ;
  ku      = 0;
  kx      = 0;

  for (i=0; i<numnp; i++)
  {
    nidu[i] = 0; for (j=0; j<ndf; j++) if (idu.x[ku++] > 0) nidu[i]++;
    nidx[i] = 0; for (j=0; j<ndm; j++) if (idx.x[kx++] > 0) nidx[i]++;

    if (i > 0)
    {
      pidu[i] = pidu[i-1] + nidu[i-1];
      pidx[i] = pidx[i-1] + nidx[i-1];
    }
  }

  // assemble pattern of sub-matrices

  ListArray< Vector<int> > i2I, i2IPos;

  i2I.setDim(numnp);

  i2IPos.setDim(numnp);

  posMtxSub2.setDim(domType.n);
  posMtxSub3.setDim(domType.n);
  posMtxSub4.setDim(domType.n);

  colDblPos.setDim(domType.n);

  for (dm=0; dm<domType.n; dm++)
  {
    F2I = freeNodeToIntfNode[dm].x;
    
    FN  = freeNode(dm).x;
    nFN = freeNode(dm).n;

    if (!isFreeSurface(*domPtr[dm]))
    {
      if (isFluid(*domPtr[dm]) && domMode[dm] != ELIMINATE && uType == DISPLACEMENT)
    
        prgError(1,fct,"use uType = VELOCITY if fluid domain is not eliminated!");
    
      prepareMatrixPatternSub1(dm, mtx, nidu.x, pidu.x, nidx.x, pidx.x, i2I, i2IPos);

      posMtxSub2[dm] = mtx.x.n;
                   
      prepareMatrixPatternSub2(dm, mtx, nidu.x, pidu.x, nidx.x, pidx.x);
      
      posMtxSub3[dm] = mtx.x.n;
        
      prepareMatrixPatternSub3(dm, mtx, nidu.x, pidu.x, nidx.x, pidx.x);
      
      posMtxSub4[dm] = mtx.x.n;
    }
    else
    {
      prepareMatrixPatternSub1(dm, mtx, nidu.x, pidu.x, nidx.x, pidx.x, i2I, i2IPos);
      
      posMtxSub2[dm] = mtx.x.n;

      prepareMatrixPatternFreeSurface(dm, mtx, nidu.x, pidu.x, nidx.x, pidx.x);
    }
  }

  // if required rearrange matrix in compressed row format
  // if required make matrix symmetric

  VectorArray<int> perm, tmp, &compr = ((SolverSparse*)solver)->compr;

  if (((SolverSparse*)solver)->comprMtxFlg)
  {
    // sort matrix entries

    COUT << "sorting matrix entries for compressed row format ...\n\n";

    computerTime.go("sorting matrix entries");

    //mtx.quickSort(true,true,perm);
    //if (symFlag) mtx.makeSymmetric(true,perm);

    mtx.copySort(true,true,symFlag,perm);

    computerTime.stopAndPrint("sorting matrix entries");

    // check for double entries

    for (i=1; i<mtx.row.n; i++)
      if (mtx.row[i] == mtx.row[i-1])
        if (mtx.col[i] == mtx.col[i-1])
          prgError(100,fct,"double entries! Can the solver cope?!");

    // generate compr

    compr.setDim(mtx.nRow+1);
    j = 0;
    compr[j++] = 1;
    for (i=1; i<mtx.row.n; i++)
      if (mtx.row[i] != mtx.row[i-1]) compr[j++] = i+1;

    compr[j] = mtx.col.n+1;
    
    if (j != mtx.nRow) prgError(200,fct,"fatal error!");

    COUT << "matrix compression done.\n\n";

    // prepare vector tmp for calculation of posSwap

    tmp.setDim(perm.n);

    for (i=0; i<perm.n; i++) tmp.x[perm.x[i]] = i;

    perm.free();
  }
  else 
  {
    tmp.setDim(mtx.x.n);

    for (i=0; i<tmp.n; i++) tmp[i] = i;
  }

  // copy to MatrixSparseArray, print pattern if matrix is small

  MatrixSparseArray<double> &mtx2 = ((SolverSparse*)solver)->mtx;

  mtx2 = mtx;

  if (neqTot < 31) mtx2.printPattern();

  if (((SolverSparse*)solver)->comprMtxFlg)
    prgCheckCompressedMatrix(mtx2.nRow, mtx2.nCol, mtx2.x, mtx2.row, mtx2.col, compr);

  // if required, prepare posSwap

  List< Vector<int> > posSwapTmp;

  VectorArray<bool> doneFlag;

  doneFlag.setDim(mtx.x.n);

  for (i=0; i<tmp.n; i++) doneFlag[i] = false;

  n = -1;
  j = 0;

  while (1)
  {
    while (j < tmp.n) if (tmp[j] == j || doneFlag[j]) j++; else break;

    if (j == tmp.n) break;

    posSwapTmp.add(new Vector<int>); n++;

    posSwapTmp[n].append(j);

    i = j;

    k = 0;

    do
    {
      posSwapTmp[n].append(tmp[i] - i);

      doneFlag[i] = true;

      i = tmp[i];

    } while (i != j);  //} while (!doneFlag[i]); is equivalent
  }
  mtx.free();

  tmp.free();

  posSwap = posSwapTmp;

  //for (i=0; i<posSwap.n; i++) cout << posSwap[i] << "\n";

  // finalise

  ((SolverSparse*)solver)->currentStatus = PATTERN_OK;

  for (dm=0; dm<domType.n; dm++)

    if (domMode[dm] != ELIMINATE) { delete domPtr[dm]->solver; domPtr[dm]->solver = NULL; }

  computerTime.stopAndPrint(fct);

  return;
}












int InterfaceN::calcStiffnessAndResidual(int printRes, bool zeroMtx, bool zeroRes)
{
  char fct[] = "InterfaceN::calcStiffnessAndResidual";

  computerTime.go(fct);

  int dm, i, j, k, n, numnpDom, ndfDom, ndmDom;

  double fact[7], dmy0, dmy1,
         *X, *Xn, *D, *V, *Vn, *U, *Un, *dU, *dUn, *ddU, *ddUn;

  if (zeroMtx) ((SolverSparse*)solver)->mtx.zero();

  if (zeroRes) r.zero();

  s = ((SolverSparse*)solver)->mtx.x.x;

  for (dm=0; dm<domType.n; dm++)
  {
    // provide correct initial guess for resolved domain interiors
    // (Should this be turned into member functions for Fluid, Solid, etc.?)

    if (firstIter && domMode[dm] != ELIMINATE)
    {
      numnpDom = domPtr[dm]->numnp;
      ndfDom   = domPtr[dm]->ndf;
      ndmDom   = domPtr[dm]->ndm;

      X  = domPtr[dm]->x.x;
      Xn = domPtr[dm]->xn.x;
      D  = domPtr[dm]->d.x;
      V  = domPtr[dm]->v.x;
      Vn = domPtr[dm]->vn.x;

      U    = domPtr[dm]->u.x;
      Un   = domPtr[dm]->un.x;
      dU   = domPtr[dm]->u3.x;
      dUn  = domPtr[dm]->u4.x;
      ddU  = domPtr[dm]->u5.x;
      ddUn = domPtr[dm]->u6.x;

      if (isFluid(*domPtr[dm]) && domMode[dm] == RESOLVE)
      {
        // sort out initial guees for mesh motion

        n = numnpDom * ndmDom;

        if (initGuess == KEEP_DISPLACEMENT)
        {
          fact[0] = domPtr[dm]->td[22];

          for (i=0; i<n; i++) V[i] = fact[0] * Vn[i];
        }
        else if (initGuess == KEEP_VELOCITY)
        {
          fact[0] = (1. - domPtr[dm]->td[22]) / domPtr[dm]->td[20];
      
          for (i=0; i<n; i++) { D[i] = fact[0] * Vn[i]; X[i] = Xn[i] + D[i]; }
        }
      }
      else if (isSolid(*domPtr[dm]) && initGuess == KEEP_VELOCITY)
      {
        // sort out initial guess for solid if KEEP_VELOCITY

        fact[0] = - domPtr[dm]->td[11] / domPtr[dm]->td[10];
        fact[1] = (1. - domPtr[dm]->td[12]) / domPtr[dm]->td[10];
        fact[2] = - domPtr[dm]->td[13] / domPtr[dm]->td[10];
        fact[3] = domPtr[dm]->td[14];
        fact[4] = domPtr[dm]->td[15];
        fact[5] = domPtr[dm]->td[16];
        fact[6] = domPtr[dm]->td[17];

        n = 0;

        for (i=0; i<numnpDom; i++) for (j=0; j<ndmDom; j++)
        {
          k = i * ndfDom + j;

          dU[k] = dUn[k];
           U[k] = fact[0] * Un[k] + fact[1] * dUn[k] + fact[2] * ddUn[k];
           X[n] = Xn[n++] + U[k];
         ddU[i] = fact[3] * U[k] + fact[4] * Un[k] + fact[5] * dUn[k] + fact[6] * ddUn[k];
        }
      }
    }

    F2I = freeNodeToIntfNode[dm].x;

    FN  = freeNode(dm).x;
    nFN = freeNode(dm).n;

    // prepare derivative factors and set freeNode degrees of freedom

    fact[0] = 1. / domReacTime[dm];
    fact[1] = - (1. - domReacTime[dm]) * fact[0];
    fact[2] = fact[0];
    fact[3] = 1.;
    fact[4] = 1. / domPtr[dm]->td[20];
    fact[5] = fact[0] * fact[4];

    if (isSolid(*domPtr[dm]))
    {
      if (uType == VELOCITY) { fact[2] *= td[9]; fact[3] = td[9]; }
    
      for (i=0; i<nFN; i++)
      {
        for (j=0; j<FN[i].dofU->n; j++) FN[i].u[j] = u.x[F2I[i]*ndf+j];
      }
    }
    else if (isFluid(*domPtr[dm]) || isFreeSurface(*domPtr[dm]))
    { 
      if (uType == DISPLACEMENT) fact[2] *= td[9];
    
      for (i=0; i<nFN; i++)
      {
        for (j=0; j<FN[i].dofU->n; j++) FN[i].u[j] = du.x[F2I[i]*ndf+j];
        for (j=0; j<FN[i].dofX->n; j++) FN[i].x[j] =  x.x[F2I[i]*ndm+j];
      }
    }
    else prgError(10,fct,"not yet implemented!");
   
    // calculate domain contributions to residual and stiffness
 
    if (!isFreeSurface(*domPtr[dm]))
    {
      if (domMode[dm] == ELIMINATE) domPtr[dm]->eliminate(true,(printRes==3));
      
      else if (domMode[dm] == RESOLVE)
      
        domPtr[dm]->prepareForExternalSolver((void*)solver,r.x+nOff[dm]-1+domPtr[dm]->nequ,false);
      
      else  domPtr[dm]->prepareForExternalSolver((void*)solver,NULL,(printRes==3));
      
      calcStiffnessAndResidualSub1(dm,fact);
      
      calcStiffnessAndResidualSub2(dm,fact);
      
      calcStiffnessAndResidualSub3(dm,fact);
    }
    else
    {
      domPtr[dm]->prepareForExternalSolver(NULL,NULL,false);
      
      calcStiffnessAndResidualSub1(dm,fact);
      
      calcStiffnessAndResidualFreeSurface(dm,fact);
    }
  }

  // permute matrix coefficients if required

  for (k=0; k<posSwap.n; k++)
  {
    i = posSwap[k][0];

    dmy0 = s[i];

    for (j=1; j<posSwap[k].n; j++)
    {
      i += posSwap[k][j];

      dmy1 = s[i];

      s[i] = dmy0;

      dmy0 = dmy1;
    }
  }

  //((SolverSparse*)solver)->mtx.print(6,3);

  firstIter = false;

  // print residual

  rNormPrev = rNorm; 
  rNorm     = r.norm();

  if (printRes > 1) 
  {
    COUT << domain.name(this);
    printf("  %11.4e\n",rNorm);
    if (printRes == 3)
    {
      i=0; while (i<domMode.n)
        if (domMode[i] != RESOLVE && domMode[i] != DUMMY)
          { cout << "\n"; break; } else i++;
    }
  }

  // finalise

  solver->currentStatus = ASSEMBLY_OK;

  s = NULL;

  ctimCalcStiffRes += computerTime.stop(fct);

  return 0;

  MatrixSparseArray<double> &mtx2 = ((SolverSparse*)solver)->mtx;

  if (((SolverSparse*)solver)->comprMtxFlg)
    prgCheckCompressedMatrix(mtx2.nRow, mtx2.nCol, mtx2.x, mtx2.row, mtx2.col, 
                             ((SolverSparse*)solver)->compr);

  return 0;
}













int InterfaceN::factoriseSolveAndUpdate(void)
{
  char fct[] = "InterfaceN::factoriseSolveAndUpdate";

  computerTime.go(fct);

  int i, ii, j, k, l, m, dm, indf, indm, n, ndfDom, ndmDom, *ID, *IDU, *IDX;

  double *DD, *UPDT, *D, *X, *V, *Xn, *Vn, fact1, fact2;

  BndNode *BNDNODE;

  // factorise and solve

  DD = solver->factoriseAndSolve(r.x);

  // update freeNode degrees of freedom
 
  if (uType == DISPLACEMENT) UPDT = u.x;

    else if (uType == VELOCITY) UPDT = du.x;

      else prgError(1,fct,"what on earth!!??");

  ID = idu.x;

  n = numnp * ndf;

  indf = -1;

  for (i=0; i<n; i++)
  {
    ii = ID[++indf] - 1;

    if (ii > -1) UPDT[indf] += DD[ii];
  }

  // update non-Lagrangian freeNode coordinate displacements (if required)

  if (neqx > 0)
  {
    UPDT = u.x;

    ID = idx.x;

    n = numnp * ndm;

    indm = -1;

    for (i=0; i<n; i++)
    {
      ii = ID[++indm] - 1;

      if (ii > -1) UPDT[indm] += DD[nequ+ii] * mpapTime.dt;
    }
  }

  // update domain data

  for (dm=0; dm<domType.n; dm++)
  {
    F2I = freeNodeToIntfNode[dm].x;

    if (domMode[dm] != ELIMINATE)
    {
      // update degrees of freedom

      ndfDom = domPtr[dm]->ndf;
      IDU    = domPtr[dm]->idu.x;
      UPDT   = domPtr[dm]->u.x;

      if (domPtr[dm]->nequ > 0)
      {
        n      = domPtr[dm]->numnp * ndfDom;
        m      = nOff[dm] - 2;
        indf   = -1;

        for (i=0; i<n; i++) { ii = IDU[++indf]; if (ii > 0) UPDT[indf] += DD[m+ii]; } 
      }

      for (i=0; i<domPtr[dm]->uDep.n; i++) domPtr[dm]->uDep[i].update(UPDT,ndfDom);

      // update mesh

      if (domMode[dm] == RESOLVE && domPtr[dm]->isALE(CHECKNEQX))
      {
        ndmDom = domPtr[dm]->ndm;
        IDX    = domPtr[dm]->idx.x;
        UPDT   = domPtr[dm]->d.x;
        X      = domPtr[dm]->x.x;
        Xn     = domPtr[dm]->xn.x;
        D      = domPtr[dm]->d.x;
        V      = domPtr[dm]->v.x;
        Vn     = domPtr[dm]->vn.x;
        fact1  = 1. / domPtr[dm]->td[20];
        fact2  = - domPtr[dm]->td[22] * fact1;
        n      = domPtr[dm]->numnp * ndmDom;
        m      = nOff[dm] - 2 + domPtr[dm]->nequ;
        indm   = -1;

        for (i=0; i<n; i++)
        {
          ii = IDX[++indm];
          if (ii > 0)
          {
            V[indm] += DD[m+ii];

            D[indm] = fact1 * V[indm] + fact2 * Vn[indm];

            X[indm] = Xn[indm] + D[indm]; 
          }
        }

        for (i=0; i<domPtr[dm]->xDep.n; i++)  domPtr[dm]->xDep[i].updateWithIncrement(UPDT,ndmDom);
      }

      domPtr[dm]->updateIterStep();
    }
  }
  ctimFactSolvUpdt += computerTime.stop(fct);

  return 0;
}












void InterfaceN::globalDiffStiffTest(double ddd, int dig, int dig2, bool gfrmt)
{
  char fct[] = "InterfaceN::globalDiffStiffTest";

  int    dm, d, e, jc, i, ii, j, k, *IDU = idu.x, *IDX = idx.x, o, q, n;

  double dd[6] = {-3.*ddd, -2.*ddd, -ddd, +ddd, +2.*ddd, +3.*ddd },
         fact1, fact2, *U, *V, *Vn, *D, *X, *Xn,
	 *R     = new double[6*neqTot],
	 *kdiff = new double[neqTot*neqTot],
	 *kanly;

  COUT << "\n";
  COUT << "diff-test for " << domain.name(this) << ", global stiffness matrix\n";
  COUT << "-----------------------------------------------------------\n\n";

  // loop over columns of stiffness

  // freeNode degrees of freedom

  n = numnp * ndf;

  if (uType == VELOCITY) U = du.x; else U = u.x;

  for (jc=0; jc<n; jc++) // nodes
  {
    j = IDU[jc] - 1;

    if (j > -1)
    {
      // loop over perturbations
          
      for (k=0; k<6; k++)
      {
        // apply pertubation
            
        U[jc] += dd[k];

        updateIterStep();
  
        // calculate residual

        calcStiffnessAndResidual(0);

        // remove pertubation
      
        U[jc] -= dd[k];
      
        updateIterStep();
      
        // loop over rows of s, store residual

        for (i=0; i<neqTot; i++) R[i*6+k] = r[i];
      }

      // loop over rows of s

      for (i=0; i<neqTot; i++)
      	
        kdiff[j*neqTot+i] = ( +       R[i*6+0]
                              -  9. * R[i*6+1]
                              + 45. * R[i*6+2]
                              - 45. * R[i*6+3]
                              +  9. * R[i*6+4]
                              -       R[i*6+5] ) / (60. * ddd);
    }
  }

  n = numnp * ndm;

  U = x.x;

  q = nequ;

  for (jc=0; jc<n; jc++) // nodes
  {
    j = IDX[jc] - 1;

    if (j > -1)
    {
      // loop over perturbations
          
      for (k=0; k<6; k++)
      {
        // apply pertubation
            
        U[jc] += dd[k] * mpapTime.dt;

        //updateIterStep();
  
        // calculate residual

        calcStiffnessAndResidual(0);

        // remove pertubation
      
        U[jc] -= dd[k] * mpapTime.dt;
      
        //updateIterStep();
      
        // loop over rows of s, store residual

        for (i=0; i<neqTot; i++) R[i*6+k] = r[i];
      }

      // loop over rows of s

      for (i=0; i<neqTot; i++)

        kdiff[(j+q)*neqTot+i] = ( +       R[i*6+0]
                                  -  9. * R[i*6+1]
                                  + 45. * R[i*6+2]
                                  - 45. * R[i*6+3]
                                  +  9. * R[i*6+4]
                                  -       R[i*6+5] ) / (60. * ddd);
    }
  }

  for (dm=0; dm<domType.n; dm++)
  {
    if (domMode[dm] != ELIMINATE)
    {
      U   = domPtr[dm]->u.x;
      n   = domPtr[dm]->numnp * domPtr[dm]->ndf;
      IDU = domPtr[dm]->idu.x;

      o = nOff[dm] - 1;

      for (jc=0; jc<n; jc++)
      {
        j = IDU[jc] - 1;

        if (j > -1)
        {
          // loop over perturbations
      	  
          for (k=0; k<6; k++)
          {
            // apply pertubation
      	    
            U[jc] += dd[k];

            updateIterStep();
      
            // calculate residual

            calcStiffnessAndResidual(0);

            // remove pertubation
      	
            U[jc] -= dd[k];

            updateIterStep();
      	
            // loop over rows of s, store residual

            for (i=0; i<neqTot; i++) R[i*6+k] = r[i];
          }

          // loop over rows of s

          for (i=0; i<neqTot; i++)

            kdiff[(j+o)*neqTot+i] = ( +       R[i*6+0]
                                      -  9. * R[i*6+1]
                                      + 45. * R[i*6+2]
                                      - 45. * R[i*6+3]
                                      +  9. * R[i*6+4]
                                      -       R[i*6+5] ) / (60. * ddd);
        }
      }

      if (domMode[dm] == RESOLVE && domPtr[dm]->isALE(CHECKNEQX))
      {
        n   = domPtr[dm]->ndm * domPtr[dm]->numnp;
        IDX = domPtr[dm]->idx.x;

        q = o + domPtr[dm]->nequ;
 
        fact1 = 1. / domPtr[dm]->td[20];
        fact2 = domPtr[dm]->td[22];
 
        V  = domPtr[dm]->v.x;
        Vn = domPtr[dm]->vn.x;
        D  = domPtr[dm]->d.x;
        X  = domPtr[dm]->x.x;
        Xn = domPtr[dm]->xn.x;

        for (jc=0; jc<n; jc++) // nodes
        {
          j = IDX[jc] - 1;

          if (j > -1)
          {
            // loop over perturbations
                  
            for (k=0; k<6; k++)
            {
              // apply pertubation
                    
              V[jc] += dd[k];
              D[jc] = (V[jc] - fact2 * Vn[jc]) * fact1;
              X[jc] = Xn[jc] + D[jc];

              // calculate residual
      
              calcStiffnessAndResidual(0);
      
              // remove pertubation
              
              V[jc] -= dd[k];
              D[jc] = (V[jc] - fact2 * Vn[jc]) * fact1;
              X[jc] = Xn[jc] + D[jc];
      
              // loop over rows of s, store residual
      
              for (i=0; i<neqTot; i++) R[i*6+k] = r[i];
            }
      
            // loop over rows of s
      
            for (i=0; i<neqTot; i++)
      
              kdiff[(j+q)*neqTot+i] = ( +       R[i*6+0]
                                        -  9. * R[i*6+1]
                                        + 45. * R[i*6+2]
                                        - 45. * R[i*6+3]
                                        +  9. * R[i*6+4]
                                        -       R[i*6+5] ) / (60. * ddd);
          }
        }
      }
    }
  }
 
  // calculate stiffness

  calcStiffnessAndResidual(0);

  // compare the matrices

  solver->copyToSimpleMatrix(kanly);
  
  prgCompareTwoSimpleMatrices(kdiff,                         // matrix 1
		              kanly,                         // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      neqTot,neqTot,                 // matrix dimension
			      dig,dig2,gfrmt,                // format
			      0,                             // indentation
			      true,                          // interactive
			      true);                         // row/column numbers

  delete [] kdiff;
  delete [] kanly;
  delete [] R;
 
  return;
}








void InterfaceN::plotSolverMatrixPattern(char *fileName)
{
  if (solver != NULL) ((SolverSparse*)solver)->mtx.plotPattern(fileName);
 
  return;
}










void InterfaceN::printComputerTime(bool reset, int detailFlg)
{
  InterfaceMatch::printComputerTime(reset,detailFlg);

  COUT << "----------------------------------------------------\n";

  if (doneCalcSub1) { COUT;
        printf("InterfaceN::calcStiffResSub1:  %7.3f sec ->%5.1f %\n",
               ctimCalcSub1, ctimCalcSub1/ctimSinceLastCall*100.); }

  if (doneCalcSub2) { COUT;
        printf("InterfaceN::calcStiffResSub2:  %7.3f sec ->%5.1f %\n",
               ctimCalcSub2, ctimCalcSub2/ctimSinceLastCall*100.); }

  if (doneCalcSub3) { COUT;
        printf("InterfaceN::calcStiffResSub3:  %7.3f sec ->%5.1f %\n",
               ctimCalcSub3, ctimCalcSub3/ctimSinceLastCall*100.); }

  if (doneCalcFreeSurf) { COUT;
        printf("InterfaceN::calcStiffResFree:  %7.3f sec ->%5.1f %\n",
               ctimCalcFreeSurf, ctimCalcFreeSurf/ctimSinceLastCall*100.); }

  COUT; printf("InterfaceN::calcStiffnessRes:  %7.3f sec ->%5.1f %\n",
               ctimCalcStiffRes, ctimCalcStiffRes/ctimSinceLastCall*100.);

  COUT; printf("InterfaceN::factoriseSolvUpdt: %7.3f sec ->%5.1f %\n",
               ctimFactSolvUpdt, ctimFactSolvUpdt/ctimSinceLastCall*100.);

  if (reset)
  {
    ctimCalcSub1     = 0.;  doneCalcSub1     = false;
    ctimCalcSub2     = 0.;  doneCalcSub2     = false;
    ctimCalcSub3     = 0.;  doneCalcSub3     = false;
    ctimCalcFreeSurf = 0.;  doneCalcFreeSurf = false;

    ctimFactSolvUpdt = 0.;
    ctimCalcStiffRes = 0.;
  }

  return;
}

