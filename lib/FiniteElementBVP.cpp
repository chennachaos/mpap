
#include <iostream>

#include "FunctionsProgram.h"
#include "FiniteElementBVP.h"
#include "DomainTree.h"
#include "DomainType.h"
#include "Debug.h"
#include "List.h"
#include "TimeFunction.h"
#include "ComputerTime.h"
#include "SolverMA41.h"

//#ifndef WINDOWS
//  #include "SolverPARDISO.h"
//#endif

#include "PropertyTypeEnum.h"
#include "MathGeom.h"

extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern ComputerTime       computerTime;


using namespace std;



FiniteElementBVP::FiniteElementBVP(void)                       
{                                                  
  solver = NULL;

  localStiffnessError = 0;

  ctimFactSolvUpdt = 0.;
  ctimCalcStiffRes = 0.;

  tol = -2.;

  freeCirculationFlag = false;

  checkElementOutput = false;
  checkNodalData     = false;

  // add new type
  
  DomainType *finiteElementBVP = domain.newType(FINITEELEMENTBVP,MESH);
  
  if (finiteElementBVP == NULL) return;  // domain type already exists

  finiteElementBVP->key.addNew("control",
                               "free circulation",
                               "check element output",
                               "check nodal data");
  return;
}




	                                          
FiniteElementBVP::~FiniteElementBVP(void)                     
{               
  if (solver != NULL) delete solver;
	
  if (debug) std::cout << " FiniteElementBVP destructor\n\n";

  return;
}






void FiniteElementBVP::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString *word;
 
  char fct[] = "FiniteElementBVP::readInputData";

  int nw, i;

  switch (domain[FINITEELEMENTBVP].key.whichBegins(line))
  {
    case  0: cout << "     FINITEELEMENTBVP: reading control ...\n\n";

             if (tol > -1)         prgError(1,fct,"'control' has already been read!");
	     
	     line.getNextLine(Ifile);
	     
	     nw = line.split(&word);
            
             if (!word[0].toDbl(&tol))      prgError(2,fct,"input error in 'control'!");

             for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     
       	     line.getNextLine(Ifile);

	     break;

    case  1: freeCirculationFlag = true;

             line.getNextLine(Ifile);

             break;

    case  2: checkElementOutput = true;

             line.getNextLine(Ifile);

             break;

    case  3: checkNodalData = true;

             line.getNextLine(Ifile);

             break;

    case -1: // go and inherit from MESH
	     
	     this->Mesh::readInputData(Ifile, line); 
	     
	     break;
  }
 
  return;
}





void FiniteElementBVP::prepareInputData(void)
{
  char fct[] = "FiniteElementBVP::prepareInputData";

  // call ancestor function

  Mesh::prepareInputData();

  
  std::cout << "     FINITEELEMENTBVP: prepare input data ...\n\n";
 
  if (tol < -1.) prgError(1,fct,"control data is missing!");


  return;
}





void FiniteElementBVP::printInfo(void)
{
  this->Mesh::printInfo();
	
  return;
}







void FiniteElementBVP::timeUpdate(void)
{
  int e, i, j, k, indf = 0, indm = 0, *idu1 = &(idu(1,1)), numnpXndf = numnp * ndf, nivEl;

  double *U = u.x, *Un = un.x, *dU = u3.x, 
	 *dUn = u4.x, *ddU = u5.x, *ddUn = u6.x;

  // basic time update only.
  
  // un <- u
  // dun <- du
  // ddun <- ddu
  
  for (i=0; i<numnpXndf; i++)
  {
      Un[i] =   U[i];
     dUn[i] =  dU[i];
    ddUn[i] = ddU[i];
  }
	     
  // prescribed displacements (set increments only)

  updateUDepIncrements();

  // set iteration flag
  
  firstIter = true;
  
  localStiffnessError = 0;

  // update dU, ddU if required
  
  updateIterStep();
  
  return;   
}






void FiniteElementBVP::updateIterStep(void)
{

  // We don't do anything else here. Things have to be taken care of by the elements.
      
  return;
}







int FiniteElementBVP::calcStiffnessAndResidual(int printRes, bool zeroMtx, bool zeroRes)
{
  char fct[] = "FiniteElementBVP::calcStiffnessAndResidual";

  computerTime.go(fct);

  int e, nst;

  bool fs, fp;

  reac.zero();

  //cout << u << "\n";

  if (nequ == 0) goto noDoF;

  if (solver->currentStatus < INIT_OK) 
  {
    prgWarning(1,"FiniteElementBVP::calcStiffnessAndResidual","initialise solver first!"); 
    return 2;
  }

  if (solver == NULL) 
  {
    COUT << "You need to select a solver first!\n\n";    return 1;
  }

  if (zeroMtx) solver->zeroMtx();

  if (firstIter) rNorm = -1.;
  
  if (zeroRes) r.zero();

  if (checkNodalData)
  {
    if (prgNAN(x.x,numnp*ndm)) prgError(2,fct,"nan entry in x!");
    if (prgNAN(u.x,numnp*ndf)) prgError(2,fct,"nan entry in u!");
  }
 
  for (e=0; e<numel; e++)
  {
    localStiffnessError = elem[e]->calcStiffnessAndResidual();

    if (localStiffnessError != 0) 
    {
      cout; COUT << "local element failure!\n\n";
      return localStiffnessError;
    }

    if (checkElementOutput)
    {
      nst = elem[e]->nen() * elem[e]->ndf();

      fs = prgNAN(s,nst*nst);
      fp = prgNAN(p,nst);
      
      if (fs || fp)
      {
        cout << "\n";
        cout << "  " << domain.name(this) << ":  failure in element " << e + 1 << "\n\n";
        if (!fp) prgError(1,fct,"nan entry in s!");
        if (!fs) prgError(1,fct,"nan entry in p!");
        prgError(1,fct,"nan entry in p and s!");
      }
    }

    solver->assembleElemMtx(elem[e]);

    solver->assembleElemVec(elem[e],firstIter);
  }


 // cout << x << endl;   cout << endl;   cout << endl;
//
  cout << "      r Vector ...:  " << endl;
  cout << endl;
  for(int ii=0;ii<r.n;ii++)
  {
    cout << fixed << '\t' << ii << '\t' << r[ii] << endl;
  }
  cout << endl;
  cout << endl;
//

  addForces();

//
  cout << "      r Vector ...:  " << endl;
  cout << endl;
  for(int ii=0;ii<r.n;ii++)
  {
    cout << fixed << '\t' << ii << '\t' << r[ii] << endl;
  }
  cout << endl;
  cout << endl;
//

  if (freeCirculationFlag) applyCirculationCutFor2DPotentialFlow();

  firstIter = false;
  rNormPrev = rNorm; 
  rNorm     = r.norm();
  
  if (printRes > 1) { COUT << domain.name(this); printf("  %11.4e\n",rNorm); }
  
  solver->currentStatus = ASSEMBLY_OK;

  ctimCalcStiffRes += computerTime.stop(fct);

  //cout << r << "\n";

  //((SolverSparse*)solver)->mtx.printPattern();
  //((SolverSparse*)solver)->mtx.print(8,3);

  //exit(0);

  return 0;

  noDoF:  // this might sometimes be the case for truss or membrane structures
          // in the context of domain interaction

  for (e=0; e<numel; e++)
  {
    localStiffnessError = elem[e]->calcStiffnessAndResidual();

    if (localStiffnessError != 0) 
    {
      cout; COUT << "local element failure!\n\n";
      return localStiffnessError;
    }
    elem[e]->assembleReactionForces();
  }
  addForces();

  rNorm = 0.;

  if (printRes > 1) { COUT << domain.name(this); printf("  %11.4e (nothing to solve for!)\n",0.); }

  ctimCalcStiffRes += computerTime.stop(fct);

  return 0;
}









void FiniteElementBVP::setSolver(int slv, int *parm, bool cIO)
{
  char fct[] = "FiniteElementBVP::setSolver";

  if (nequ == 0 && neqx == 0) 
  { 
    prgWarning(1,fct,"no degrees of freedom! -> nothing to be done!");
    solverOK = true;
    return; 
  }

  int numProc;

  switch (slv) 
  {
    case  1: // MA41 ..........................

             solver = (Solver*) new SolverMA41;

             solver->prepareMatrixPattern(this);

             solver->printInfo();

	     if (solver->initialise() != 0) return; 
	    
             if (!isALE() || neqx == 0) break;

             solverMesh = (Solver*) new SolverMA41;	     

             solverMesh->prepareMatrixPattern(this);

             COUT << "mesh solver:\n\n";

             solverMesh->printInfo();

	     if (solverMesh->initialise() != 0) return; 
	    
	     break;

/*	     
    case  2: // PARDISO(sym)
    case  3: // PARDISO(unsym)

#ifndef WINDOWS

             if (parm != NULL) numProc = parm[0]; else numProc = 1;

             numProc = min(MAX_PROCESSORS,numProc);

             solver = (Solver*) new SolverPARDISO;

             solver->prepareMatrixPattern(this);

             solver->printInfo();

             if (slv == 2) if (solver->initialise(numProc,PARDISO_STRUCT_SYM) != 0) return;
             if (slv == 3) if (solver->initialise(numProc,PARDISO_UNSYM) != 0) return;

             if (!isALE() || neqx == 0) break;

             solverMesh = (Solver*) new SolverPARDISO;

             solverMesh->prepareMatrixPattern(this);

             COUT << "mesh solver:\n\n";

             solverMesh->printInfo();

             if (slv == 2) if (solverMesh->initialise(numProc,PARDISO_STRUCT_SYM) != 0) return;
             if (slv == 3) if (solverMesh->initialise(numProc,PARDISO_UNSYM) != 0) return;

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
  
  if (solver != NULL)      solver    ->checkIO = cIO;
  if (solverMesh != NULL)  solverMesh->checkIO = cIO;

/*
  //  for(int ii=0;ii<uDep.n;ii++)
    cout << '\t' << " uDep.n " << '\t' << uDep.n << endl;
    cout << endl;
    cout << '\t' << uDep[0].uc << '\t' << uDep[0].duc << '\t' << uDep[0].ucBase << endl;
*/
  return;
}





int FiniteElementBVP::factoriseSolveAndUpdate(void)
{
  char fct[] = "FiniteElementBVP::factoriseSolveAndUpdate";

  computerTime.go(fct);

  int i, j, ii, *IDU = idu.x, indf = 0, numnpXndf = numnp * ndf;

  double *U  = u.x, *du;


  // factorise and solve

  if (nequ > 0) 
  {
    du = solver->factoriseAndSolve(r.x);
  
    //cout << r << "\n";
 
    if (du == NULL) 
      { cout << "\n"; COUT << "Error! solver (substitution) has failed!\n\n"; return 2; }
   
    // update solution vector

    for (i=0; i<numnp; i++)  
    {	  	  
      for (j=0; j<ndf; j++)  {  ii = IDU[indf+j];  if (ii > 0)  U[indf+j] += du[ii-1]; }

      indf += ndf; 
    }
  }

  // update dependent displacements

  for (i=0; i<uDep.n; i++)  uDep[i].update(U,ndf);


  // cout << u << endl;   cout << endl;   cout << endl;

  //cout << u << " B\n";

  ctimFactSolvUpdt += computerTime.stop(fct);

  return 0;
}









bool FiniteElementBVP::converged(void)
{
  if (rNorm < tol && localStiffnessError == 0) return true;  

  return false;
}






bool FiniteElementBVP::diverging(double factor)
{
  if (rNormPrev > -0.1 && rNorm / rNormPrev > factor) return true;  

  if (localStiffnessError != 0) return true;

  if (prgNAN(rNorm)) return true;
  
  return false;
}






void FiniteElementBVP::addForces(void)
{
  int i, j, k, *idu1;

  double f, *reac1;
 
  for (k=0; k<frcNd.n; k++)
  {
    i = frcNd[k];
	  
    idu1  = &(idu(i,1));

    reac1 = &(reac(i,1));
    
    for (j=0; j<ndf; j++)
    {
      f = frc[k*ndf+j];

      if (frcTmFct[k*ndf+j] > -1) f *= timeFunction[frcTmFct[k*ndf+j]].prop;
	    
      if (idu1[j] > 0) { r[idu1[j]-1] += f; } //std::cout << f << "\n\n";  }

      reac1[j] += f;
    }
  }

  return;
}







void FiniteElementBVP::elementDiffStiffTest(double ddd, int el, int dig, int dig2, bool gfrmt)
{
  int e = el - 1;

  if (e < 0)      e = 0; 
  if (e >= numel) e = numel - 1;

  COUT << "\n";
  COUT << "diff-test for " << domain.name(this) << ", element " << e+1 << "\n";
  COUT << "---------------------------------------------------------\n";
  
  elem[e]->diffStiffTest(ddd,dig,dig2,gfrmt);
  
  return;
}







void FiniteElementBVP::elementDiffMeshDerivTest(double ddd, int el, int dig, int dig2, bool gfrmt)
{
  if (!isALE()) { COUT << "no ALE domain !!\n\n"; return; }

  int e = el - 1;

  if (e < 0)      e = 0; 
  if (e >= numel) e = numel - 1;

  COUT << "\n";
  COUT << "mesh derivative diff-test for " << domain.name(this) << ", element " << e+1 << "\n";
  COUT << "---------------------------------------------------------\n";
  
  elem[e]->diffMeshDerivTest(ddd,dig,dig2,gfrmt);
  
  return;
}








void FiniteElementBVP::globalDiffStiffTest(double ddd, int dig, int dig2, bool gfrmt)
{
  int    d, e, jc = 0, i, ii, jnd, jj, j, k, *IDU = idu.x, nivEl, row0 = 1, col0 = 1;

  bool   keep = firstIter; firstIter = false;
  
  double dd[6] = {-3.*ddd, -2.*ddd, -ddd, +ddd, +2.*ddd, +3.*ddd },
	 *R     = new double[6*nequ],
	 *kdiff = new double[nequ*nequ],
	 *kanly,// = new double[nequ*nequ],
	 *intVar1, *intVar2, *U = u.x;

  COUT << "\n";
  COUT << "diff-test for " << domain.name(this) << ", global stiffness matrix\n";
  COUT << "-----------------------------------------------------------\n\n";

  // loop over columns of stiffness

  for (jnd=0; jnd<numnp; jnd++) // nodes
  {
    for (jj=0; jj<ndf; jj++)  // dof
    {
      j = IDU[jc++];

      if (j > 0)
      {
        j--;
	      
        // loop over perturbations
  	    
        for (k=0; k<6; k++)
        {
          // apply pertubation
  	      
          u(jnd+1,jj+1) += dd[k];

          for (d=0; d<uDep.n; d++) uDep[d].update(U,ndf);
	  
          updateIterStep();
  
          // calculate residual

	  calcStiffnessAndResidual(0);

          // remove pertubation
  	
          u(jnd+1,jj+1) -= dd[k];
  	
          for (d=0; d<uDep.n; d++) uDep[d].update(U,ndf);
	  
          updateIterStep();
  	
          // loop over rows of s, store residual

          for (i=0; i<nequ; i++) R[i*6+k] = r[i];
        }

        // loop over rows of s

        for (i=0; i<nequ; i++)
  		
          kdiff[j*nequ+i] = ( +       R[i*6+0]
                              -  9. * R[i*6+1]
                              + 45. * R[i*6+2]
                              - 45. * R[i*6+3]
                              +  9. * R[i*6+4]
                              -       R[i*6+5] ) / (60. * ddd);
      }
    }
  }
 
  // calculate stiffness

  calcStiffnessAndResidual(0);

  for (e=0; e<numel; e++)
  {
    intVar1 = elem[e]->intVar1;
    intVar2 = elem[e]->intVar2;
	  
    nivEl = elem[e]->nivGP() * elem[e]->nGaussPoints();

    for (i=0; i<nivEl; i++) intVar2[i] = intVar1[i];
  }	

  firstIter = keep;
  
  // compare the matrices

  solver->copyToSimpleMatrix(kanly);
  
  prgCompareTwoSimpleMatrices(kdiff,                         // matrix 1
		              kanly,                         // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      nequ,nequ,                     // matrix dimension
			      dig,dig2,gfrmt,                // format
			      0,                             // indentation
			      true,                          // interactive
			      true);                         // row/column numbers

  delete [] kdiff;
  delete [] kanly;
  delete [] R;
  
  return;
}






void FiniteElementBVP::reset(void)
{
  int i, numnpXndf = numnp * ndf;

  double *U = u.x, *Un = un.x, *dU = u3.x, 
	 *dUn = u4.x, *ddU = u5.x, *ddUn = u6.x;

  for (i=0; i<numnpXndf; i++)
  {
      U[i] =   Un[i];
     dU[i] =  dUn[i];
    ddU[i] = ddUn[i];
  }
  
  for (i=0; i<uDep.n; i++) uDep[i].timeUpdate();
 
  return;
}









void FiniteElementBVP::printComputerTime(bool reset, int detailFlg)
{
  Mesh::printComputerTime(reset,detailFlg);

  COUT << "----------------------------------------------------\n";

  COUT; printf("FiniteElementBVP::calcStiffRes:%7.3f sec ->%5.1f %\n",
               ctimCalcStiffRes, ctimCalcStiffRes/ctimSinceLastCall*100.);

  COUT; printf("FiniteElementBVP::factSolvUpdt:%7.3f sec ->%5.1f %\n",
               ctimFactSolvUpdt, ctimFactSolvUpdt/ctimSinceLastCall*100.);

  if (reset)
  {
    ctimFactSolvUpdt = 0.;
    ctimCalcStiffRes = 0.;
  }

  return;
}







void FiniteElementBVP::applyCirculationCutFor2DPotentialFlow(void)

{
  char fct[] = "FiniteElementBVP::applyCirculationCutFor2DPotentialFlow";

  if (nBnd > 1) { prgWarning(1,fct,"Invalid mesh topology! No cut? -> abort"); return; }

  if (elemGrp.n != 1) { prgWarning(2,fct,"multiple element groups? -> abort"); return; }

  if (ndm != 2) { prgWarning(2,fct,"use 2D with Element2D3nodedLinearPoisson"); return; }

  if (elemGrp[0].elemProp[ELEMENTTYPE]->name != "2D3nodedLinearPoisson")

    { prgWarning(2,fct,"wrong element type! -> abort"); return; }

  if (!solverOK) { prgWarning(2,fct,"initialise solver first"); return; }

  int row, col, e, i, j, k, l;

  Vector<int> nodesOnCut, elementsAttached;

  List< Vector<int> > nodesInElement;

  // generate nodesOnCut

  if (nBnd > 1) prgError(1,fct,"multiply connected domain!");

  for (i=0; i<numnp; i++)  nodeFlag[i]          = false;
  for (i=0; i<nBndNd; i++) nodeFlag[bndNd[i]-1] = true;

  for (i=0; i<uDep.n; i++)
  {
    if (uDep[i].master.n > 0) 
    {
      if (!nodeFlag[uDep[i].nd-1]) prgError(1,fct,"uDep nodes need to be boundary nodes!");

      else nodeFlag[uDep[i].nd-1] = false;

      //if (uDep[i].master.n != 2) prgError(1,fct,"uDep nodes need to have 2 masters exactly!");

      k = uDep[i].masterNd[0]; if (k == numnp) k = uDep[i].masterNd[1];

      l = uDep[i].nd;

      k = k + k - 2;
      l = l + l - 2;

      if (abs(x.x[k] - x.x[l]) > 1.e-8 || abs(x.x[k+1] - x.x[l+1]) > 1.e-8)
      {
        //cout << k/2+1 << " = " << l/2+1 << "\n";

        x.x[k]   = x.x[l];
        x.x[k+1] = x.x[l+1];

        prgError(1,fct,"nodes do not match!");
      }
    }
  }

  return; 
  i = 0; while (i < nBndNd)
  {
    if (!nodeFlag[bndNd[i]-1]) 
    {
      nodesOnCut.append(bndNd[i]);

      if (nodesOnCut.n > 1)

        if (nodeFlag[bndNd[i-1]-1]) prgError(2,fct,"error in uDep data; cut not continuous?!");
    }
    i++;
  }

  // generate elementsAttached and nodesInElement

  for (i=0; i<nodesOnCut.n; i++)
  {
    for (j=0; j<nodeElem[nodesOnCut[i]-1].n; j++) 
    {
      e = nodeElem[nodesOnCut[i]-1][j];

      if (!elementsAttached.contains(e,&l))
      {
        l = elementsAttached.n;
        elementsAttached.append(e);
        nodesInElement.add(new Vector<int>);
      }

      k = 0; while (k < 3) if (elem[e]->ix[k] == nodesOnCut[i]) break; else k++;

      if (k == 3) prgError(1,fct,"fatal error!");

      nodesInElement[l].append(k);
    }
  }

  // loop over elementsAttached

  //cout << nodesOnCut << "\n";
  //cout << elementsAttached << "\n";
  //for (i=0; i<nodesInElement.n; i++)
  //  cout << elementsAttached[i]+1 << ":" << nodesInElement[i] << "\n";

  MatrixSparseArray<double> &mtx = ((SolverSparse*)solver)->mtx;

  j = idu.x[numnp-1];

  i = 0; while (i < mtx.x.n) if (mtx.row[i] == j && mtx.col[i] == j) break; else i++;

  double *MTX = mtx.x.x + i,
         *R   = r.x + j - 1,
         Circ = u.x[numnp-1];

  for (i=0; i<elementsAttached.n; i++)
  {
    elem[elementsAttached[i]]->calcStiffnessAndResidual();
 
    for (j=0; j<nodesInElement[i].n; j++)
    {
      row = nodesInElement[i][j];

      for (k=0; k<nodesInElement[i].n; k++)
      {
        col = nodesInElement[i][k];
        
        *MTX -= 2. * s[col*3+row];

        *R   += Circ * 2. * s[col*3+row];
      }
    }
  }

  return;
}










void FiniteElementBVP::adjustElementAssembly(int p0, int fl)
{
  char fct[] = "FiniteElementBVP::adjustElementAssembly";

  int e;

  for (e=0; e<numel; e++) elem[e]->adjustAssembly(p0,fl);

  return;
}








void FiniteElementBVP::plotSolverMatrixPattern(char *fileName)
{
  if (solver != NULL) ((SolverSparse*)solver)->mtx.plotPattern(fileName);
 
  return;
}






