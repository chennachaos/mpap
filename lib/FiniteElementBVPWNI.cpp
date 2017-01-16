
#include "FiniteElementBVPWNI.h"
#include "InterfaceGS.h"
#include "Fluid.h"
#include "Solid.h"
#include "ComputerTime.h"


extern ComputerTime computerTime;



FiniteElementBVPWNI::FiniteElementBVPWNI(void)
{ 
  isNeumannProblem = false;

  ctimSolve = 0.;  doneSolve = false;

  // add new type
  
  DomainType *finiteElementBVPWNI = domain.newType(FINITEELEMENTBVPWNI,FINITEELEMENTBVPWI);

  if (finiteElementBVPWNI == NULL) return;  // domain type exists already

  //finiteElementBVPWNI->key.addNew();

  return;
}







	                                          
FiniteElementBVPWNI::~FiniteElementBVPWNI(void)                     
{

  return;
}









void FiniteElementBVPWNI::readInputData(std::ifstream &Ifile, MyString &line)
{
  char fct[] = "FiniteElementBVPWNI::readInputData";

  switch (domain[FINITEELEMENTBVPWNI].key.whichBegins(line))
  {
    case -1: // go and inherit from FiniteElementBVPWI
	     
	     this->FiniteElementBVPWI::readInputData(Ifile, line); 
	     
	     break;
  }

  return;
}











void FiniteElementBVPWNI::prepareInputData(void)
{
  char fct[] = "FiniteElementBVPWNI::prepareInputData"; 

  // call ancestor function

  FiniteElementBVPWI::prepareInputData();
 
  //cout << "     FINITEELEMENTBVPWNI: prepare input data ...\n\n";

  return;
}










void FiniteElementBVPWNI::prepareInteractions(void)
{
  // go and inherit from ancestors

  FiniteElementBVPWI::prepareInteractions();

  cout << "     FINITEELEMENTBVPWNI: preparing interactions ...\n\n"; 

  // check whether there is an InterfaceGS domain which 
  // requires this domain to be a Neumann problem

  int i = 0, j = -1;

  if (freeNode.n > 0)
  {
    while (i < domain.nDomainOfType(INTERFACEGS))
    {
      j = domain(INTERFACEGS,i++).doIfThisIsNeumannProb(this);

      if (j == 1) isNeumannProblem = true;

      if (j != -1) break;
    }

    // transform into Neumann problem if required

    if (isNeumannProblem) transformIntoNeumannProblem();
  }

  return;
}










void FiniteElementBVPWNI::printInfo(void)
{
  this->FiniteElementBVPWI::printInfo();
	
  return;
}










void FiniteElementBVPWNI::solve(int DNFlag, bool printRes)
{
  char fct[] = "FiniteElementBVPWNI::solve";

  computerTime.go(fct);

  if (DNFlag == RETURN_FORCES)
  {
    if (isNeumannProblem) prgError(1,fct,"this domain is set-up as a Neumann problem");

    eliminate(false,printRes);

    ctimSolve += computerTime.stop(fct);

    doneSolve = true;

    return;
  }

  if (DNFlag != RETURN_DISPLACEMENTS) prgError(2,fct,"invalid DNFlag");
  
  if (!isNeumannProblem) prgError(2,fct,"this domain is not set-up as a Neumann problem");
    
  // solve Neumann problem

  if (!isSolid(*this)) prgError(2,fct,"so far this has been tested only for Neumann Solids");

  int i, j, k, maxIter = 30, *IDU = idu.x, *DOFU, pR;

  if (printRes) pR = 2; else pR = 0;

  k = frcNd.n - freeNode.n;

  for (i=0; i<freeNode.n; i++)
  {
    DOFU = freeNode[i].dofU->x;

    for (j=0; j<freeNode[i].dofU->n; j++)
    {
      frc[k*ndf+DOFU[j]] = freeNode[i].reac[j];
    }
    k++;
  }

  i = 0;

  while (1)
  {
    calcStiffnessAndResidual(pR);

    if (converged()) break;

    factoriseSolveAndUpdate();

    updateIterStep();

    if (++i > maxIter) break;
  }

  if (!converged()) prgWarning(1,fct,"NO CONVERGENCE!");

  if (pR) cout << "\n";

  // retrieve displacements for freeNodes

  k = numnp - freeNode.n;
 
  for (i=0; i<freeNode.n; i++)
  {
    DOFU = freeNode[i].dofU->x;

    for (j=0; j<freeNode[i].dofU->n; j++)
    {
      freeNode[i].u[j] = u.x[k*ndf+DOFU[j]];
    }
    k++;
  }

  ctimSolve += computerTime.stop(fct);

  doneSolve = true;

  return;
}











void FiniteElementBVPWNI::transformIntoNeumannProblem(void)
{
  char fct[] = "FiniteElementBVPWNI::transformIntoNeumannProblem";  

  if (isFluid(*this) && neqx > 0)

    prgError(1,fct,"not possible for domain with independent mesh motion!");

  //cout << fct << "\n";

  int i, j, k, m, n, m1, n1, *IDU, *DOFU;

  //cout << "\n" << idu << "\n";

  //for (i=0; i<uDep.n; i++) uDep[i].printInfo(); cout << "\n";

  MatrixFullArray<double> dtmp;
  MatrixFullArray<int>    itmp;

  // delete Dirichlet solver

  if (solver != NULL) { solver->free(); solver = NULL; }

  // extend node data by freeNodes

  n  = numnp * ndm;
  m  = numnp * ndf;
  n1 = (numnp + freeNode.n) * ndm;
  m1 = (numnp + freeNode.n) * ndf;

  dtmp = x;
  x.setDim(numnp + freeNode.n, ndm, true);
  for (i=0; i<n; i++) x.x[i] = dtmp.x[i];
  for (i=n; i<n1; i++) x.x[i] = 0.;

  dtmp = x0;
  x0.setDim(numnp + freeNode.n, ndm, true);
  for (i=0; i<n; i++) x0.x[i] = dtmp.x[i];
  for (i=n; i<n1; i++) x0.x[i] = 0.;

  dtmp = xn;
  xn.setDim(numnp + freeNode.n, ndm, true);
  for (i=0; i<n; i++) xn.x[i] = dtmp.x[i];
  for (i=n; i<n1; i++) xn.x[i] = 0.;

   d.setDim(numnp + freeNode.n, ndm, true);
  d0.setDim(numnp + freeNode.n, ndm, true);
   v.setDim(numnp + freeNode.n, ndm, true);
  vn.setDim(numnp + freeNode.n, ndm, true);

  reacMesh.setDim(numnp + freeNode.n, ndm, true);

  dtmp = u;
  u.setDim(numnp + freeNode.n, ndf, true);
  for (i=0; i<m; i++) u.x[i] = dtmp.x[i];
  for (i=m; i<m1; i++) u.x[i] = 0.;

  dtmp = un;
  un.setDim(numnp + freeNode.n, ndf, true);
  for (i=0; i<m; i++) un.x[i] = dtmp.x[i];
  for (i=m; i<m1; i++) un.x[i] = 0.;

  dtmp = u3;
  u3.setDim(numnp + freeNode.n, ndf, true);
  for (i=0; i<m; i++) u3.x[i] = dtmp.x[i];
  for (i=m; i<m1; i++) u3.x[i] = 0.;

  dtmp = u4;
  u4.setDim(numnp + freeNode.n, ndf, true);
  for (i=0; i<m; i++) u4.x[i] = dtmp.x[i];
  for (i=m; i<m1; i++) u4.x[i] = 0.;

  dtmp = u5;
  u5.setDim(numnp + freeNode.n, ndf, true);
  for (i=0; i<m; i++) u5.x[i] = dtmp.x[i];
  for (i=m; i<m1; i++) u5.x[i] = 0.;

  dtmp = u6;
  u6.setDim(numnp + freeNode.n, ndf, true);
  for (i=0; i<m; i++) u6.x[i] = dtmp.x[i];
  for (i=m; i<m1; i++) u6.x[i] = 0.;

  dtmp = reac;
  reac.setDim(numnp + freeNode.n, ndf, true);
  for (i=0; i<m; i++) reac.x[i] = dtmp.x[i];
  for (i=m; i<m1; i++) reac.x[i] = 0.;

  nodeFlag.setDim(numnp + freeNode.n);

  r.setDim(max(n,m));

  outp.setDim(max(max(n,m),max(numel*ndf,numel*ndm))); 

  // extend point load data

  VectorArray<int>    itmpv;
  VectorArray<double> dtmpv;

  itmpv = frcNd;
  frcNd.setDim(itmpv.n + freeNode.n);
  for (i=0; i<itmpv.n; i++) frcNd[i] = itmpv[i];
  for (i=0; i<freeNode.n; i++) frcNd[itmpv.n+i] = numnp + i + 1;

  itmpv = frcTmFct;
  frcTmFct.setDim(frcNd.n * ndf);
  for (i=0; i<itmpv.n; i++) frcTmFct[i] = itmpv[i];
  for (i=itmpv.n; i<frcTmFct.n; i++) frcTmFct[i] = -1;

  dtmpv = frc;
  frc.setDim(frcNd.n * ndf);
  for (i=0; i<dtmpv.n; i++) frc[i] = dtmpv[i];
  for (i=dtmpv.n; i<frc.n; i++) frc[i] = 0.;
 
  // regenerate idu and idx

  itmp = idu;

  idu.setDim(numnp + freeNode.n, ndf, true);

  IDU = idu.x;
  for (i=0; i<m; i++) IDU[i] = itmp.x[i];

  for (i=0; i<freeNode.n; i++)
  {
    if (freeNode[i].idu == NULL) prgError(1,fct,"freeNode[i].idu == NULL");

    for (j=0; j<ndf; j++) IDU[m+j] = 0;

    DOFU = freeNode[i].dofU->x;
    n    = freeNode[i].dofU->n;

    for (j=0; j<n; j++)
    {
      if (freeNode[i].idu[j] != 0)
      {
        IDU[m+DOFU[j]] = 0;     // modify this to account for i
                                // non-zero prescribed displacements
        freeNode[i].idu[j] = 0; //
      }
      else
      {
        IDU[m+DOFU[j]] = ++nequ;

        freeNode[i].idu[j] = nequ;
      }
    }
    m += ndf;
  }

  idx.setDim(numnp + freeNode.n, ndm, true);
  idx.zero();
  
  // adjust uDep  

  DependentDoF *uDepj;

  for (i=0; i<bndNode.n; i++)
  {
    for (j=0; j<bndNode[i].uDep.n; j++)
    {
      uDepj = bndNode[i].uDep[j];

      uDepj->ucBase = 0.;

      n = 0;
      for (k=0; k<bndNode[i].freeU.n; k++)
      {
        if (IDU[(bndNode[i].freeU[k]->nd + numnp - 1)*ndf+uDepj->dof-1] > 0) n++;
      }

      uDepj->masterNd.setDim(n);
      uDepj->masterDoF.setDim(n);
      uDepj->alpha.setDim(n);
      uDepj->beta.setDim(n);

      n = 0;
      for (k=0; k<bndNode[i].freeU.n; k++)
      { 
        if (IDU[(bndNode[i].freeU[k]->nd + numnp - 1)*ndf+uDepj->dof-1] > 0)
        {
          uDepj->masterNd[n]  = bndNode[i].freeU[k]->nd + numnp;
          uDepj->masterDoF[n] = uDepj->dof;
          uDepj->alpha[n]     = bndNode[i].dat[k];
          uDepj->beta[n++]    = bndNode[i].dat[k];
        }
      }
      //cout << uDepj->alpha << ", " << uDepj->beta << "\n";
    }
  }

  // finalise

  VectorArray<int> *tmp;

  tmp = new VectorArray<int> [numnp + freeNode.n];

  for (i=0; i<numnp; i++) tmp[i] = nodeElem[i];

  delete [] nodeElem;

  nodeElem = tmp;

  tmp = new VectorArray<int> [numnp + freeNode.n];

  for (i=0; i<numnp; i++) tmp[i] = nodeNode[i];

  delete [] nodeNode;

  nodeNode = tmp;

  numnp += freeNode.n;

  prepareAndCheckUDep(&uDep);

  //cout << "\n" << idu << "\n";

  //for (i=0; i<uDep.n; i++) uDep[i].printInfo(); cout << "\n";

  return;
}








void FiniteElementBVPWNI::printComputerTime(bool reset, int detailFlg)
{
  FiniteElementBVPWI::printComputerTime(reset,detailFlg);

  if (doneSolve)
  {
    COUT << "----------------------------------------------------\n";

    COUT; printf("FiniteElementBVPWNI::solve:    %7.3f sec ->%5.1f %\n",
                 ctimSolve, ctimSolve/ctimSinceLastCall*100.);
  }

  if (reset)
  {
    ctimSolve = 0.;  doneSolve = false;
  }

  return;
}

