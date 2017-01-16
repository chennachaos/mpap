
#include <iostream>

#include "FunctionsProgram.h"
#include "FreeSurface.h"
#include "DomainTree.h"
#include "DomainType.h"
#include "Fluid.h"
#include "DataBlockTemplate.h"
#include "PropertyTypeEnum.h"


 


extern DomainTree domain;




FreeSurface::FreeSurface(void)                       
{                                                  
  // add new type
  
  DomainType *freeSurface = domain.newType(FREESURFACE,FLUID);
  
  if (freeSurface == NULL) return;  // domain type already exists

  freeSurface->key.addNew("dimensions");

  return;
}




	                                          
FreeSurface::~FreeSurface(void)                     
{                                                 
  return;
}







void FreeSurface::readInputData(std::ifstream &Ifile, MyString &line)
{
  int i, nw;
	
  MyString *word;
 
  char fct[] = "FreeSurface::readInputData";
 
  switch (domain[FREESURFACE].key.whichBegins(line))
  {
    case 0: cout << "     FREESURFACE: reading dimensions ...\n\n";
	   
	    line.getNextLine(Ifile);
	    
	    nw = line.split(&word);
            
	    if (nw < 1)                prgError(1,fct,"input error in 'dimensions'!");
	            
            if (!word[0].toInt(&ndf))  prgError(2,fct,"input error in 'dimensions'!");
	            
            if (nw>1) { if (!word[1].toInt(&ndm)) prgError(3,fct,"input error in 'dimensions'!"); }

            else ndm = ndf;

            for (i=0; i<nw; i++) word[i].free(); delete [] word;

            if (ndf != ndm || (ndf != 2 && ndf != 3))

              prgError(4,fct,"invalid values in 'dimensions'!");
	    
       	    line.getNextLine(Ifile);

	    break;

    case -1: // go and inherit from MESH
	     
	     this->Fluid::readInputData(Ifile, line); 
	     
	     break;
  }
  return;
}








void FreeSurface::prepareInputData(void)
{
  char fct[] = "FreeSurface::prepareInputData";

  int e, i, j, n;

  // issue warning if freeNodes have been set in input file

  if (freeNodeTmp.n > 0)

    prgWarning(1,fct,"all nodes will be used as interface nodes for this domain type!");

  // define all nodes as freeNodes
 
  Vector<void*> vTmp, v2Tmp;

  Vector<double> dTmp;

  dTmp.append(1.);

  freeNodeTmp.free();
  bndNodeTmp.free();
  freeNodeIsConnected.free();

  n = dofU.n;
  dofU.add(new VectorArray<int>); 
  dofU[n].setDim(ndf);
  for (i=0; i<ndf; i++) dofU[n][i] = i;
  
  if (dofX.n != n) prgError(10,fct,"dofX.n != dofU.n");
  dofX.add(new VectorArray<int>); 
  dofX[n].setDim(ndm);
  for (i=0; i<ndm; i++) dofX[n][i] = i;

  for (i=0; i<numnp; i++)
  {
    freeNodeTmp.add(new FreeNode(INTERPOLATION_FN,x.x+i*ndm,1,1,&(dofU[n]),&(dofX[n])));

    vTmp.free();

    vTmp.append((void*)&(freeNodeTmp.lastItem()));

    v2Tmp = vTmp;

    bndNodeTmp.add(new BndNode(i+1,x.x+i*ndm,vTmp,v2Tmp,dTmp));

    freeNodeIsConnected.append(true);
  }

  tol = 1.;

  // go and inherit from ancestors

  Fluid::prepareInputData();
  
  std::cout << "     FREESURFACE: prepare input data ...\n\n";

  // for two dimensional case, set boundary nodes and
  // check orientation of elements

  VectorArray<int> tmp;

  if (nen == 2 && ndm == 2)
  {
    delete [] bnd;   bnd   = NULL;
    delete [] bndNd; bndNd = NULL;

    nBnd   = 1;
    nBndNd = 2;

    bnd   = new int* [nBnd+1];
    bndNd = new int  [nBndNd+1];

    bnd[0] = bndNd;
    bnd[1] = bndNd + 2;

    tmp.setDim(numnp);

    j = 0;

    tmp.zero(); for (e=0; e<numel; e++) tmp[elem[e]->ix[0]-1]++;

    for (i=0; i<numnp; i++)
    {
      if (tmp[i] == 0)
      {
        if (j != 0) prgError(1,fct,"invalid mesh!"); else bndNd[j++] = i + 1;
      }
      else if (tmp[i] != 1) prgError(2,fct,"invalid mesh!");
    }

    tmp.zero(); for (e=0; e<numel; e++) tmp[elem[e]->ix[1]-1]++;

    for (i=0; i<numnp; i++)
    {
      if (tmp[i] == 0)
      {
        if (j != 1) prgError(3,fct,"invalid mesh!"); else bndNd[j++] = i + 1;
      }
      else if (tmp[i] != 1) prgError(4,fct,"invalid mesh!");
    }

    if (j != 2) prgError(5,fct,"invalid mesh!");
  }

  // resize s

  int nst = nen * ndm;

  delete [] s; s = new double [nst * nst * 2];


  

  return;
}









void FreeSurface::prepareInteractions(void)
{
  // go and inherit from ancestors

  Fluid::prepareInteractions();

  cout << "     FREESURFACE: preparing interactions ...\n\n"; 

  char fct[] = "FreeSurface::prepareInteractions";
 
  return;
}









void FreeSurface::printInfo(void)
{ 
  COUT << "What do you want to know?\n\n";

  return;
}






void FreeSurface::setSolver(int, int *, bool)
{
  COUT << "FreeSurface::setSolver: nothing to be done here!\n\n";

  return;
}






void FreeSurface::prepareForExternalSolver(void *interfaceSolver,
                                                  double *rMesh,
                                                  bool printMeshRes)
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::
// This member function is tailored to be evoked from a 
// child of InterfaceMatch, specifically InterfaceN.
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::
{
  char fct[] = "FreeSurface::prepareForExternalSolver";

  //cout << fct << "\n";

  int e, i, j, k, n;

  //cout << firstIter << "\n";

  bool meshFlag = (isALE() && neqx > 0), 
       elimMeshFlag = (rMesh == NULL);

  //computerTime.go(fct);

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

    for (j=0; j<bndNode[i].xDep.n; j++)
    {
      k = bndNode[i].xDep[j]->dof- 1;
      n = bndNode[i].xDep[j]->nd - 1;

      x.x[n*ndm+k] = xn.x[n*ndm+k] + bndNode[i].xDep[j]->uc;
      v.x[n*ndm+k] = td[22] * vn.x[n*ndm+k] + td[20] * bndNode[i].xDep[j]->uc;
    }
  }

  //updateIterStep();

  return;

/*
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

  //ctimElim += computerTime.stop(fct);
*/
  return;
}









void FreeSurface::setLagrangianFreeNodes(VectorArray<int> &)
{
  // dummy function to hide inherited member function

  // nothing to be done here!

  return;
}


