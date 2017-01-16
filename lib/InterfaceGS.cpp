
#include <iostream>

#include "FunctionsProgram.h"
#include "InterfaceGS.h"
#include "DomainTree.h"
#include "SolverMA41.h"
//#include "SolverPARDISO.h"
#include "Solid.h"
#include "Fluid.h"
#include "FreeSurface.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "FunctionsSupport.h"


extern DomainTree   domain;
extern ComputerTime computerTime;
extern MpapTime     mpapTime;




InterfaceGS::InterfaceGS(void)                       
{                                                  
  solverOK = true;

  relaxMode = 0;

  relaxParam[0] = 1.;

  // add new type
  
  DomainType *interfaceGS = domain.newType(INTERFACEGS,INTERFACEMATCH);
  
  if (interfaceGS == NULL) return;  // domain type already exists

  interfaceGS->key.addNew("relaxation");

  return;
}




	                                          
InterfaceGS::~InterfaceGS(void)                     
{                                                 
  return;
}







void InterfaceGS::readInputData(std::ifstream &Ifile, MyString &line)
{
  char fct[] = "InterfaceGS::readInputData";

  char *relaxModeName[] = { "NONE", "FIXED", "AITKEN", "MYAITKEN", NULL };

  int i, nw;

  MyString *word;

  switch (domain[INTERFACEGS].key.whichBegins(line))
  {
    case  0: cout << "     INTERFACEGS: reading relaxation ...\n\n";

             if (relaxMode > 0) prgError(1,fct,"'relaxation' has already been read!");
	     
             line.getNextLine(Ifile);
	     
             nw = line.split(&word);
            
             if (nw < 1) prgError(1,fct,"input error in 'relaxation'!");
		     
             relaxMode = word[0].which(relaxModeName);

             if (relaxMode < 0) prgError(3,fct,"invalid relaxMode name in 'relaxation'!");

             if (relaxMode != NONE)
             {
               if (nw < 2) prgError(1,fct,"provide parameter for relaxation method!");

               if (!word[1].toDbl(relaxParam)) prgError(5,fct,"input error in 'relaxation'!");
             }
             relaxParam[1] = relaxParam[0];
 
             for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     
       	     line.getNextLine(Ifile);

             break;

    case -1: // go and inherit from INTERFACEMATCH
	     
	     InterfaceMatch::readInputData(Ifile,line); 
	     
	     break;
  }
 
  return;
}








void InterfaceGS::prepareInputData(void)
{
  char fct[] = "InterfaceGS::prepareInputData";

  // go and inherit from ancestors

  InterfaceMatch::prepareInputData();
 
  
  std::cout << "     INTERFACEGS: prepare input data ...\n\n";

  
// intialise whatever


  return;
}









void InterfaceGS::prepareInteractions(void)
{
  // go and inherit from ancestors

  InterfaceMatch::prepareInteractions();

  cout << "     INTERFACEGS: preparing interactions ...\n\n"; 

  char fct[] = "InterfaceGS::prepareInteractions";
 
  int admissibleDomMode[] = { RETURN_DISPLACEMENTS, RETURN_FORCES },
      aitkenRelax[] = { AITKEN, MYAITKEN };

  neumannDom = -1;

  for (int dm=0; dm<domMode.n; dm++)
  {
    if (!isElementOf(domMode[dm], admissibleDomMode, 2))

      prgError(1,fct,"invalid domain mode!");

    if (domMode[dm] == RETURN_DISPLACEMENTS)
    {
      if (neumannDom != -1)

        prgError(1,fct,"set one Neumann (RETURN_DISPLACEMENTS) subproblem only!");

      neumannDom = dm;
    }
  }

  if (uType != DISPLACEMENT) prgError(2,fct,"set uType = DISPLACEMENT!");

  // allocate some memory

  r.setDim(numnp * ndf);

  trac.setDim(numnp * ndf);

  tracPrev.setDim(numnp * ndf);

  uInc.setDim(numnp * ndf);

  if (isElementOf(relaxMode, aitkenRelax, 2))
  {
    DDui.setDim(numnp * ndf);
    DDuj.setDim(numnp * ndf);
  }

  return;
}









void InterfaceGS::printInfo(void)
{ 
  COUT << "ndm ....... = " << ndm   << "\n";
  COUT << "ndf ....... = " << ndf   << "\n";
  COUT << "numnp ..... = " << numnp << "\n";
  COUT << "nequ        = " << nequ  << "\n";
  COUT << "neqx        = " << neqx  << "\n";

  if (uType == VELOCITY) COUT << "uType       = VELOCITY\n";
  else                   COUT << "uType       = DISPLACEMENT\n";

  if (initGuess == KEEP_VELOCITY) COUT << "initGuess   = KEEP_VELOCITY\n";
  else                            COUT << "initGuess   = KEEP_DISPLACEMENT\n";

  cout << "\n";

  return;
}








int InterfaceGS::doIfThisIsNeumannProb(Domain *dom)
{
  int i, j, dm = 0, *IDU = idu.x;

  while (dm < domType.n)
  {
    if (domPtr[dm] == dom) break; else dm++;
  }

  if (dm == domType.n) return -1; // dom is not connected to this interface

  if (domMode[dm] != RETURN_DISPLACEMENTS) return 0; // dom is connected as Dirichlet problem

  F2I = freeNodeToIntfNode[dm].x;

  for (i=0; i<freeNode(dm).n; i++)
  {
    freeNode(dm)[i].idu = new int [ndf];

    for (j=0; j<ndf; j++)

      if (IDU[F2I[i]*ndf+j] > 0) freeNode(dm)[i].idu[j] = 0; else freeNode(dm)[i].idu[j] = 1;
  }

  return 1; // dom is connected as Neumann problem and dom->freeNode[].idu have been set
}










void InterfaceGS::setSolver(int, int*, bool)
{
  char fct[] = "InterfaceGS::setSolver";

  COUT << fct << ": nothing to be done here!\n\n";

  return;
}









int InterfaceGS::calcStiffnessAndResidual(int printRes, bool, bool)
{
  char fct[] = "InterfaceGS::calcStiffnessAndResidual";

  int dm, i, j, k, l, n, *IDU;

  double fact, fact1, fact2;

  // solve Dirichlet problems

  tracPrev = trac;

  trac.zero();

  for (dm=0; dm<domType.n; dm++)
  {
    if (domMode[dm] == RETURN_FORCES)
    {
      F2I = freeNodeToIntfNode[dm].x;
      
      FN  = freeNode(dm).x;
      nFN = freeNode(dm).n;

      // set freeNode degrees of freedom
      
      if (isSolid(*domPtr[dm]))
      {
        for (i=0; i<nFN; i++)
        {
          for (j=0; j<FN[i].dofU->n; j++) FN[i].u[j] = u.x[F2I[i]*ndf+j];
        }
      }
      else if (isFluid(*domPtr[dm]))
      { 
        for (i=0; i<nFN; i++)
        {
          for (j=0; j<FN[i].dofU->n; j++) FN[i].u[j] = du.x[F2I[i]*ndf+j];
          for (j=0; j<FN[i].dofX->n; j++) FN[i].x[j] =  x.x[F2I[i]*ndm+j];
        }
      }
      else prgError(10,fct,"not yet implemented!");

      // solve
 
      domPtr[dm]->solve(RETURN_FORCES,(printRes==3));

      // assemble reaction forces

      fact1 = 1. / domReacTime[dm];
      fact2 = - (1. - domReacTime[dm]) * fact1;

      for (i=0; i<nFN; i++)

        for (j=0; j<FN[i].dofU->n; j++)
        {
          trac.x[F2I[i]*ndf+j] += fact1 * FN[i].reac[j] + fact2 * FN[i].reacn[j];
        }
    }
  }

  // check for convergence based on traction force residual

  if (!firstIter) for (i=0; i<numnp*ndf; i++) r[i] = trac.x[i] - tracPrev.x[i];

  else { fact1 = tol + tol; for (i=0; i<numnp*ndf; i++) r[i] = fact1; }

  if (printRes > 1)
  {
    COUT << domain.name(this);
    if (!firstIter) 
    {
      if (relaxMode == MYAITKEN)    printf("  %11.4e%3d\n",r.norm(),nMyAitken-1);
      else if (relaxMode == AITKEN) printf("  %11.4e%7.3g\n",r.norm(),relaxParam[1]);
      else                          printf("  %11.4e\n",r.norm());
    }
    else printf("   **********\n");
    if (printRes == 3) cout << "\n";
  }

  if (r.norm() < tol) return 0;

  // solve Neumann problems

  dm = neumannDom;

  // set freeNode loads      

  fact1 = domReacTime[dm];
  fact2 = 1. - domReacTime[dm];

  F2I = freeNodeToIntfNode[dm].x;
      
  FN  = freeNode(dm).x;
  nFN = freeNode(dm).n;

  for (i=0; i<nFN; i++)
  {
    for (j=0; j<FN[i].dofU->n; j++)
    {
      FN[i].reac[j] = fact1 * trac.x[F2I[i]*ndf+j] + fact2 * FN[i].reacn[j];
    }
  }

  // solve

  domPtr[dm]->solve(RETURN_DISPLACEMENTS,(printRes==3));

  // calculate displacement increments

  for (i=0; i<nFN; i++)
  {
    for (j=0; j<ndf; j++)

      if (idu.x[F2I[i]*ndf+j] > 0) uInc.x[F2I[i]*ndf+j] = FN[i].u[j] - u.x[F2I[i]*ndf+j];

      else uInc.x[F2I[i]*ndf+j] = 0.;
  }

  // update interface displacements

  k = ndf * numnp;

  if (relaxMode == NONE || relaxMode == FIXED) // fixed relaxation
  {
    fact = relaxParam[0];

    for (i=0; i<k; i++) u.x[i] += fact * uInc.x[i];
  }
  else if (relaxMode == AITKEN) // standard Aitken
  {
    if (firstIter) fact = min(relaxParam[0],relaxParam[1]);
    else
    {
      for (l=0; l<k; l++) DDui.x[l] = DDuj.x[l] - uInc.x[l];

      fact = relaxParam[1] * dot_(DDuj.x,DDui.x,&k) / dot_(DDui.x,DDui.x,&k);
    }

    for (l=0; l<k; l++) 
    {
      u.x[l]   += fact * uInc.x[l];
      DDuj.x[l] = uInc.x[l];
    }
    relaxParam[1] = fact;
  }
  else if (relaxMode == MYAITKEN) // my Aitken
  {
    if (nMyAitken > roundToInt(relaxParam[0]))
    {
      ua.add(ua.takeOut(0));
      Du.add(Du.takeOut(0));
      nMyAitken--;
    }

    n = nMyAitken;

    if (ua.n < n+1)
    {
      if (ua.n != n) prgError(1,fct,"impossible error!");
      ua.add(new VectorArray<double>);
      Du.add(new VectorArray<double>);
      ua[n].setDim(k);
      Du[n].setDim(k);
    }
    fact = 10000.;
    for (l=0; l<k; l++)
    {
      ua[n].x[l] = .5 * (u.x[l] + u.x[l] + uInc.x[l]);
      uInc.x [l]*= fact;
      Du[n].x[l] =                         uInc.x[l];
    }
    n = getAlphaForMyAitken();

    //for (i=0; i<n; i++) cout << alp[i] << ", "; cout << "\n";

    fact = 1.;

    for (i=0; i<n; i++) fact -= alp[i];

    for (l=0; l<k; l++) u.x[l] = fact * ua[n].x[l];

    for (i=0; i<n; i++)

      for (l=0; l<k; l++) u.x[l] += alp[i] * ua[i].x[l];

    nMyAitken = n + 1;
  }
  else prgError(2,fct,"invalid relaxation mode!");

  // finalise

  firstIter = false;

  return 0;
}













int InterfaceGS::getAlphaForMyAitken(void)
{
  char fct[] = "InterfaceGS::getAlphaForMyAitken";

  if (nMyAitken == 0) return 0;

  int i, j, k, l, isw = -1, n0 = nMyAitken, n = n0, d = 0;

  k = numnp * ndf;

  if (rhs.n < n) 
  {
    rhs.setDim(n);
    alp.setDim(n);
    pos.setDim(n);
    mtx.setDim(n*n);
    dmx.setDim(n*n);
  }

  for (i=0; i<n0; i++)
  {
    for (l=0; l<k; l++) DDui.x[l] = Du[n0].x[l] - Du[i].x[l];
 
    for (j=0; j<n0; j++)
    {
      for (l=0; l<k; l++) DDuj.x[l] = Du[n0].x[l] - Du[j].x[l];

      mtx[j*n0+i] = dot_(DDuj.x,DDui.x,&k);
    }

    rhs[i] = dot_(Du[n0].x,DDui.x,&k);
  }

  while (1)
  {
    for (i=0; i<n; i++)
      for (j=0; j<n; j++)
        dmx.x[j*n+i] = mtx.x[(d+j)*n0+d+i];

    //prgPrintSimpleMatrix(dmx.x,n,n,8,4,true,3,false);

    if (decomplr_matrix_(dmx.x,pos.x,&n,&isw) == 0) break; else { d++; n--; }

    ua.del(0);
    Du.del(0);

    if (n == 0) prgError(1,fct,"impossible error!");
  }

  //cout << "\n" << rhs << "\n";

  for (i=0; i<n; i++) rhs.x[i] = rhs.x[i+d];

  //cout << rhs << "\n\n";

  solve_matrix_(dmx.x,pos.x,rhs.x,alp.x,&n);

  return n;
}









int InterfaceGS::factoriseSolveAndUpdate(void)
{
  char fct[] = "InterfaceGS::factoriseAndUpdate";

  //COUT << fct << ": nothing to be done here!\n\n";

  return 0;
}






bool InterfaceGS::converged(void)
{
  if (r.norm() < tol) return true; else return false;
}





void InterfaceGS::timeUpdate(void)
{
  nMyAitken = 0;

  InterfaceMatch::timeUpdate();

  return;
}


