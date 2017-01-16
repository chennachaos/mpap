
#include <iostream>

#include "FunctionsProgram.h"
#include "DomainTypeEnum.h"
#include "DomainTree.h"
#include "DomainType.h"
#include "Debug.h"
#include "Fluid.h"
#include "MpapTime.h"
#include "PropertyTypeEnum.h"


extern DomainTree domain;
extern MpapTime   mpapTime;


Fluid::Fluid(void)                       
{                                                  
  if (debug) std::cout << " Fluid constructor\n\n";

  //std::cout << ++nObj << "\n\n";

  addElemGrpProp(ALETYPE);
  
  
  // add new type

  DomainType *fluid = domain.newType(FLUID,FINITEELEMENTBVPWNI);
 
  if (fluid == NULL) return;  // domain type already exists
  
  fluid->key.addNew("control","time derivative dof");

  return;
}




	                                          
Fluid::~Fluid(void)                     
{                                                 
  if (debug) std::cout << " Fluid destructor\n\n";
  
  //std::cout << --nObj << "\n\n";

  return;
}





void Fluid::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString *word;
 
  char fct[] = "Fluid::readInputData";

  int nw, i, nn;

  VectorArray<int> tmp;

  switch (domain[FLUID].key.whichBegins(line))
  {
    case  0: cout << "     FLUID: reading control ...\n\n";

             if (tis > -1)         prgError(1,fct,"'control' has already been read!");
	     
	     line.getNextLine(Ifile);
	     
	     nw = line.split(&word);
            
	     if (nw < 3)                    prgError(1,fct,"input error in 'control'!");
		     
             if (!word[0].toDbl(&tol))      prgError(2,fct,"input error in 'control'!");
		     
             if (!word[1].toDbl(&tolMesh))  prgError(3,fct,"input error in 'control'!");
		     
	     if (!word[2].toInt(&tis))      prgError(4,fct,"input error in 'control'!");
	     if (tis < 0) tis = 0;

             for (i=0; i<10; i++) td[i] = 0.;
	    
             nn = nw; if (nn > 13) nn = 13;
	     
             for (i=3; i<nn; i++) 
	       if (!word[i].toDbl(&(td[i-3]))) prgError(5,fct,"input error in 'control'!");

             for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     
       	     line.getNextLine(Ifile);

	     break;

    case  1: // cout << "     FLUID: reading time derivative dof ...\n\n";

             tmp.setDim(ndf);

	     line.getNextLine(Ifile);
	     
	     nw = line.split(&word);
            
	     if (nw != ndf) prgError(1,fct,"input error in 'time derivative dof'!");

             for (i=0; i<ndf; i++)

               if (!word[i].toInt(tmp.x+i))

                 prgError(2,fct,"input error in 'time derivative dof'!");
		     
             for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     
             nn = 0; for (i=0; i<ndf; i++) if (tmp[i] == 1) nn++; timeIntSwap.setDim(nn);

             nn = 0; for (i=0; i<ndf; i++) if (tmp[i] == 1) timeIntSwap[nn++] = i;

       	     line.getNextLine(Ifile);

             break;

    case -1: // go and inherit from FINITEELEMENTBVPWNI
	     
	     this->FiniteElementBVPWNI::readInputData(Ifile,line); 
	     
	     break;
  }
 
  return;
}





void Fluid::prepareInputData(void)
{
  char fct[] = "Fluid::prepareInputData";

  // go and inherit from ancestors

  FiniteElementBVPWNI::prepareInputData();
 
  
  std::cout << "     FLUID: prepare input data ...\n\n";

  
// intialise whatever


  return;
}







void Fluid::prepareInteractions(void)
{
  // go and inherit from ancestors

  FiniteElementBVPWNI::prepareInteractions();


  cout << "     FLUID: preparing interactions ...\n\n"; 
  
 
  return;
}








void Fluid::setTimeParam(void)
{
  double dt = mpapTime.dt, alpf, alpm, gamm, rho;
	
  for (int i=10; i<TD_DIM; i++) td[i] = 0.;

  switch (tis)
  {
    case  1: // generalised midpoint rule

             gamm = td[0];

	     td[6-1] = gamm;
	     td[7-1] = gamm;
	     td[8-1] = gamm;
	     
             td[10] = 1. / (gamm * dt);      // U_n+1   in  dU_n+1
             td[11] = - td[10];              // U_n
             td[12] = - (1. - gamm) / gamm;  // dU_n    
	     
             td[9-1] = td[10];    // coefficient of U_n+1 in dU_n+1 

             td[20]  = td[10];   // X_n+1  in  V_n+1
             td[21]  = td[11];   // X_n
             td[22]  = td[12];   // V_n
	     
             td[30]  = 1. / td[10];        // dU_n+1   in  U_n+1
             td[31]  = - td[30] * td[11];  // U_n
             td[32]  = - td[30] * td[12];  // dU_n

             break;

    case  2: // generalised alpha-method

             rho  = td[0];
	     
             alpf = 1. / (1. + rho);
	     alpm = .5 * (3. - rho) / (1. + rho);
             gamm = .5 + alpm - alpf;
	     
	     td[6-1] = alpf;
	     td[7-1] = alpf;
	     td[8-1] = alpm;
	     
             td[10] = 1. / (gamm * dt);      // U_n+1   in  dU_n+1
             td[11] = - td[10];              // U_n
             td[12] = - (1. - gamm) / gamm;  // dU_n    
	     
             td[9-1] = td[10];    // coefficient of U_n+1 in dU_n+1 

             td[20]  = td[10];   // X_n+1  in  V_n+1
             td[21]  = td[11];   // X_n
             td[22]  = td[12];   // V_n
	     
             td[30]  = 1. / td[10];        // dU_n+1   in  U_n+1
             td[31]  = - td[30] * td[11];  // U_n
             td[32]  = - td[30] * td[12];  // dU_n

	     break;
	     
    default: prgError(1,"Fluid::setTimeParam","invalid value of tis!");
  }
  return;
}






void Fluid::timeUpdate(void)
{
  int e, i, j, k, indf = 0, indm = 0, *idu1 = &(idu(1,1)), 
      numnpXndf = numnp * ndf, numnpXndm = numnp * ndm, nivEl;

  double *X = x.x, *Xn = xn.x, *X0 = x0.x,
	 *U = u.x, *Un = un.x, *dU = u3.x, 
	 *V = v.x, *Vn = vn.x, 
	 *dUn = u4.x, *ddU = u5.x, *ddUn = u6.x,
	 *intVar1, *intVar2;

  // un <- u
  
  for (i=0; i<numnpXndf; i++)   Un[i] = U[i];
  
  // ALE stuff
  
  for (i=0; i<numnpXndm; i++)  
  { 
    Vn[i] = V[i]; 
    Xn[i] = X[i];
  }
  d.zero();
  d0.zero();

  // velocities, accelerations
  
  switch (tis)
  {
    case  0: // quasi static

	     // nothing to be done here

	     break;
	     
    case  1: // generalised midpoint rule

    case  2: // generalised alpha-method

	     for (i=0; i<numnpXndf; i++)  dUn[i] =  dU[i];
	     
	     break;

    default: prgError(1,"Fluid::timeUpdate","invalid value of tis!");
  }
  
  // prescribed displacements (set increments only)

  updateUDepIncrements();

  // set iteration flag

  firstIter = true;
  
  localStiffnessError = 0;

  // update dU  !important!
  
  updateIterStep();
  
  return;   
}






void Fluid::updateIterStep(void)
{
  int i, j, k, *idu1 = &(idu(1,1)), indf = 0, indm = 0, 
      numnpXndf = numnp * ndf;

  double *U = u.x, *Un = un.x, *dU = u3.x, *dUn = u4.x, *ddU = u5.x, *ddUn = u6.x;

  switch (tis)
  {
    case  0: // quasi static

	     // nothing to be done here

	     break;

    case  1: // generalised midpoint rule

    case  2: // generalised alpha-method

	     for (i=0; i<numnpXndf; i++)
		     
                dU[i] = td[10]*(U[i]-Un[i]) + td[12]*dUn[i];
	     
             for (k=0; k<timeIntSwap.n; k++)
             {
               j = timeIntSwap[k];

               for (i=0; i<numnp; i++)

                 dU[i*ndf+j] = td[30]*U[i*ndf+j] + td[31]*dUn[i*ndf+j] + td[32]*Un[i*ndf+j]; 
             }
	     break;
	     
    default: prgError(1,"Fluid::updateIterStep","invalid value of tis!");
  }

  if (ndm == 3) surf3D->updtFlagU = true;

  return;
}






void Fluid::printInfo(void)
{
  this->FiniteElementBVPWNI::printInfo();
	
  return;
}







void Fluid::reset(void)
{
  int i, j, numnpXndf = numnp * ndf;

  double *U = u.x, *Un = un.x, *dU = u3.x, *dUn = u4.x, *ddU = u5.x, *ddUn = u6.x;

  for (i=0; i<numnpXndf; i++) {  U[i] =  Un[i];
	                        dU[i] = dUn[i]; }	     

  for (i=0; i<uDep.n; i++) uDep[i].timeUpdate();

  return;
}


