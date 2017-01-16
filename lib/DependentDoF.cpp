
#include <iostream>

#include "DependentDoF.h"
#include "TimeFunction.h"
#include "List.h"
#include "Debug.h"



extern List<TimeFunction> timeFunction;



DependentDoF::DependentDoF(void)
{
  if (debug) std::cout << " DependentDoF constructor\n\n";

  uc = 0.;
  duc = 0.;
  ucBase = 0.;

  return;
}




DependentDoF::~DependentDoF()
{
  if (debug) std::cout << " DependentDoF destructor\n\n";
 
  return;	
}






void DependentDoF::init(void)
{
  duc = 0.;
	
  uc = 0.;
      
  return;  
}






void DependentDoF::timeUpdate(void)
{
  duc = - uc;
	
  uc = ucBase;

  if (tmFct > -1) uc *= timeFunction[tmFct].prop;

  duc += uc;
  
  return;  
}






void DependentDoF::update(double *u, int &ndf)
{
  double *U = &(u[(nd-1)*ndf+dof-1]);

  //std::cout << *U << " <- " << uc << "\n";

  *U = uc;

  for (int i=0; i<masterNd.n; i++)  *U += alpha[i] * u[(masterNd[i]-1)*ndf+masterDoF[i]-1];

  return;
}






void DependentDoF::updateWithIncrement(double *u, int &ndf)
{
  double *U = &(u[(nd-1)*ndf+dof-1]);

  *U += duc;

  return;
}





void DependentDoF::printInfo(void)
{
  cout << " nd = " << nd << ", dof = " << dof << ", masterNd = " 
       << masterNd << ", masterDoF = " << masterDoF << "\n";

  return;
}

