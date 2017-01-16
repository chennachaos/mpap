
#include <iostream>

#include "FunctionsProgram.h"
#include "FunctionsMaterial.h"
#include "PropertyTypeEnum.h"
#include "DomainTypeEnum.h"
#include "Solid.h"
#include "DomainTree.h"
#include "DomainType.h"
#include "Debug.h"
#include "Definitions.h"
#include "List.h"
#include "TimeFunction.h"
#include "MpapTime.h"
#include "WorkSpace.h"


extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern MpapTime           mpapTime;
extern WorkSpace          workSpace;


using namespace std;



Solid::Solid(void)                       
{                                                  
  if (debug) std::cout << " Solid constructor\n\n";

  addElemGrpProp(MATERIAL);

  // add new type
  
  DomainType *solid = domain.newType(SOLID,FINITEELEMENTBVPWNI);
  
  if (solid == NULL) return;  // domain type already exists

  solid->key.addNew("material", "control");

  return;
}




	                                          
Solid::~Solid(void)                     
{                                                 
  if (debug) std::cout << " Solid destructor\n\n";

  return;
}









void Solid::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString *word;
 
  char fct[] = "Solid::readInputData";

  int nw, i, nn;
 
  switch (domain[SOLID].key.whichBegins(line))
  {
    case  0: cout << "     SOLID: reading material ...\n\n";

             elemProp.add(new PropertyItem(MATERIAL));

	     elemProp[elemProp.n-1].readInputData(Ifile,line,"input error in 'material'!");

	     break;
	     
    case  1: cout << "     SOLID: reading control ...\n\n";

             if (tis > -1)         prgError(1,fct,"'control' has already been read!");

             line.getNextLine(Ifile);

             nw = line.split(&word);
            
             if (nw < 3)                    prgError(1,fct,"input error in 'control'!");
		     
             if (!word[0].toDbl(&tol))      prgError(2,fct,"input error in 'control'!");
		     
             if (!word[1].toInt(&lumpType)) prgError(3,fct,"input error in 'control'!");
	     if (lumpType < 0) lumpType = 0;

	     if (!word[2].toInt(&tis))      prgError(4,fct,"input error in 'control'!");
	     if (tis < 0) tis = 0;

             for (i=0; i<10; i++) td[i] = 0.;
	     
             nn = nw; if (nn > 13) nn = 13;
	     
             for (i=3; i<nn; i++) 
	       if (!word[i].toDbl(&(td[i-3]))) prgError(5,fct,"input error in 'control'!");

             for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     
       	     line.getNextLine(Ifile);

	     break;
	     
    case -1: // go and inherit from FINITEELEMENTBVPWNI
	     
	     this->FiniteElementBVPWNI::readInputData(Ifile, line); 
	     
	     break;
  }
 
  return;
}





void Solid::prepareInputData(void)
{
  char fct[] = "Solid::prepareInputData";

  int intf, e;

  // go and inherit from ancestors

  FiniteElementBVPWNI::prepareInputData();
 
  
  std::cout << "     SOLID: prepare input data ...\n\n";

  
  if (lumpType == -1) prgError(1,"Solid::prepareInputData","'control' data undefined!");

  if (ndf < ndm)      prgError(2,"Solid::prepareInputData","ndf < ndm not admissible!");
  
// intialise materials

  char *matTypeNames[] = MATERIAL_TYPE_NAMES;
  
  prepareElemProp(MATERIAL,matTypeNames);

  
// initialise internal variables

  for (e=0; e<numel; e++)  elem[e]->initialiseIntVar(); 


  return;
}





void Solid::prepareInteractions(void)
{
  // go and inherit from ancestors

  FiniteElementBVPWNI::prepareInteractions();


  cout << "     SOLID: preparing interactions ...\n\n"; 
  
  int i, j, e;
 

  // set Gauss point data id for multiscaleMaterial2D/3D
  
  char *matName[] = MATERIAL_TYPE_NAMES;
  
  for (int i=0; i<elemGrp.n; i++)
  {
    j = elemGrp[i].elemProp[MATERIAL]->id;
    
    if (strcmp(matName[j],"multiscaleMaterial2D") == 0
     || strcmp(matName[j],"multiscaleMaterial3D") == 0)
    {
      for (e=0; e<elemGrp[i].elem.n; e++) elemGrp[i].elem[e].setGaussPointDataId();
    }
  }
  
  return;
}





void Solid::printInfo(void)
{
  this->FiniteElementBVPWNI::printInfo();
	
  return;
}





void Solid::setTimeParam(void)
{
  double dt = mpapTime.dt, alpf, alpm, beta, gamm, rho;
	
  for (int i=10; i<TD_DIM; i++) td[i] = 0.;

  td[5-1] = dt;

  switch (tis)
  {
    case  0: // quasi static

             td[6-1] = 1.0;
	    
	     break;

    case  1: // generalised midpoint rule

             gamm = td[0];

	     td[6-1] = gamm;
	     td[7-1] = gamm;
	     td[8-1] = gamm;
	     
             td[10] = - (1. - gamm) / gamm;         
             td[11] = 1. / (gamm * dt);

             td[9-1] = td[11] * td[11];  // coefficient of u_n+1 in ddu_n+1 
	     
	     break;

    case  2: // generalised alpha-method

             rho  = td[0];
	     
             alpf = 1. / (1. + rho);
	     alpm = (2. - rho) / (1. + rho);
	     beta = .25 * (1. + alpm - alpf)*(1. + alpm - alpf);
             gamm = .5 + alpm - alpf;
	     
	     td[6-1] = alpf;
	     td[7-1] = alpf;
	     td[8-1] = alpm;
	     
             td[10] = gamm / (beta * dt);             // U_n+1   in  dU_n+1
             td[11] = - td[10];                       // U_n
             td[12] = 1. - gamm / beta;               // dU_n    
             td[13] = dt * (1. - gamm / (2. * beta)); // ddU_n
	     
             td[14] = 1. / (beta * dt * dt);          // U_n+1   in  ddU_n+1
             td[15] = - td[14];                       // U_n
             td[16] = - 1. / (beta * dt);             // dU_n
             td[17] = 1. - 1. / (2. * beta);          // ddU_n

             td[9-1] = td[14];  // coefficient of U_n+1 in ddU_n+1 
	     
	     break;
	     
    default: prgError(1,"Solid::setTimeParam","invalid value of tis!");
  }
  return;
}






void Solid::timeUpdate(void)
{
  int e, i, j, k, indf = 0, indm = 0, numnpXndf = numnp * ndf, 
                                      numnpXndm = numnp * ndm, nivEl;

  double *X = x.x, *Xn = xn.x, *X0 = x0.x,
	 *U = u.x, *Un = un.x, *dU = u3.x, 
	 *dUn = u4.x, *ddU = u5.x, *ddUn = u6.x,
	 *intVar1, *intVar2;

  // xn <- x
 
  for (i=0; i<numnpXndm; i++) Xn[i] =  X[i];

  // un <- u
  
  for (i=0; i<numnpXndf; i++)   Un[i] = U[i];
  
  // velocities, accelerations
  
  switch (tis)
  {
    case  0: // quasi static

	     // nothing to be done here

	     break;
	     
    case  1: // generalised midpoint rule

    case  2: // generalised alpha-method

	     for (i=0; i<numnpXndf; i++)
	     {
	        dUn[i] =  dU[i];
	       ddUn[i] = ddU[i];
	     }
	     
	     break;

    default: prgError(1,"Solid::timeUpdate","invalid value of tis!");
  }
  
  // update internal variables

  for (e=0; e<numel; e++)
  {
    intVar1 = elem[e]->intVar1;
    intVar2 = elem[e]->intVar2;

    nivEl = elem[e]->nivGP() * elem[e]->nGaussPoints();

    for (i=0; i<nivEl; i++) intVar1[i] = intVar2[i];
  }

  // prescribed displacements  (set increments only)

  updateUDepIncrements();
  
  // set iteration flag
  
  firstIter = true;
 
  localStiffnessError = 0;

  // update X, dU, ddU    !important!
  
  updateIterStep();

 // cout << u << endl;   cout << endl;   cout << endl;

  // cout << x << endl;   cout << endl;   cout << endl;	 
  return;   
}






void Solid::updateIterStep(void)
{
  int i, j, indf = 0, indm = 0, numnpXndf = numnp * ndf;

  double *X = x.x, *X0 = x0.x, *U = u.x, 
	 *Un = un.x, *dU = u3.x, *dUn = u4.x, *ddU = u5.x, *ddUn = u6.x;

  // update coordinates
 
  for (i=0; i<numnp; i++)
  {
    for (j=0; j<ndm; j++) {  X[indm+j] = X0[indm+j] + U[indf+j]; }
    
    indf += ndf;
    indm += ndm;
  }

  switch (tis)  // note that changing this stuff requires
	         // adjusting of ElementSolid::diffStiffTest,
		  // otherwise, 'diff' command for elements is in trouble
  {
    case  0: // quasi static

	     // nothing to be done here

	     break;

    case  1: // generalised midpoint rule

	     for (i=0; i<numnpXndf; i++)
	     {
                dU[i] = td[10]* dUn[i] + td[11]*( U[i]- Un[i]);
               ddU[i] = td[10]*ddUn[i] + td[11]*(dU[i]-dUn[i]);
	     }
	     
	     break;

    case  2: // generalised alpha-method

	     for (i=0; i<numnpXndf; i++)
	     {
                dU[i] = td[10]*(U[i]-Un[i]) + td[12]*dUn[i] + td[13]*ddUn[i];
               ddU[i] = td[14]*(U[i]-Un[i]) + td[16]*dUn[i] + td[17]*ddUn[i];
	     }
	     
	     break;
	     
    default: prgError(1,"Solid::updateIterStep","invalid value of tis!");
  }

  if (ndm == 3)
  {
    surf3D->updtFlagX = true;
    surf3D->updtFlagU = true;
  }

  return;
}





int Solid::calcStiffnessAndResidual(int printRes, bool zeroMtx, bool zeroRes)
{
  // copy internal variables intVar1 to intVar2
  // (basically forget values at t_n+1 from non-converged Newton step)

  int e, i, nivEl;
	
  double *intVar1, *intVar2;
  
  for (e=0; e<numel; e++)
  {
    intVar1 = elem[e]->intVar1;
    intVar2 = elem[e]->intVar2;
	  
    nivEl = elem[e]->nivGP() * elem[e]->nGaussPoints();

    for (i=0; i<nivEl; i++) intVar2[i] = intVar1[i]; 
  }	

  // then call the usual stuff
  
  return FiniteElementBVPWNI::calcStiffnessAndResidual(printRes, zeroMtx, zeroRes);
}








void Solid::reset(void)
{
  int e, i, numnpXndf = numnp * ndf, numnpXndm = numnp * ndm, nivEl;

  double *X = x.x, *Xn = xn.x, *U = u.x, 
	 *Un = un.x, *dU = u3.x, *dUn = u4.x, *ddU = u5.x, *ddUn = u6.x,
	 *intVar1, *intVar2;

  for (i=0; i<numnpXndm; i++) X[i] = Xn[i];
  
  for (i=0; i<numnpXndf; i++) {  U[i] =   Un[i];
                                dU[i] =  dUn[i];
                               ddU[i] = ddUn[i]; }

  for (e=0; e<numel; e++)
  {
    intVar1 = elem[e]->intVar1;
    intVar2 = elem[e]->intVar2;
	  
    nivEl = elem[e]->nivGP() * elem[e]->nGaussPoints();

    for (i=0; i<nivEl; i++) intVar2[i] = intVar1[i];
  }	

  for (i=0; i<uDep.n; i++) uDep[i].timeUpdate();
  
  return;
}












void Solid::domainSpecificNodalDataTransfer(int numnpNewAll)
{
  int i;

  for (i=0; i<numnpNewAll*ndm; i++) u.x[i] = x.x[i] - x0.x[i];

  return;
}





