
#include "ElementSolid.h"
#include "Debug.h"
#include "ElementGroup.h"
#include "FunctionsMaterial.h"
#include "FunctionsProgram.h"
#include "PropertyTypeEnum.h"


using namespace std;



ElementSolid::ElementSolid(void)
{
  if (debug) cout << " constructor ElementSolid\n\n";

  return;
}




ElementSolid::~ElementSolid()
{
  if (debug) cout << " destructor ElementSolid\n\n";

  return;
}





int ElementSolid::finiteStrain(void)
{
  return 0;
}





int ElementSolid::stressStrainState(void)
{
  return 0;
}





int ElementSolid::nGaussPoints(void)
{
  return 0;
}





int ElementSolid::nivGP(void)
{
  // get number of internal variable per Gauss point

  double dmy[10], *matDat = &(((ElementGroup*) elemGrp)->elemProp[MATERIAL]->data[0]);

  int    isw = 1, n, dmyI[3],
         matId  = ((ElementGroup*) elemGrp)->elemProp[MATERIAL]->id + 1,
  	 mDim   = matdim_(&matId),
         finite = finiteStrain(),
         sss    = stressStrainState();

  //cout << matId << " " << finite << " " << sss << " " << ndm() << "\n\n";
  
  if (mDim == 1) matlib1d_(matDat,dmy,dmy,dmy,dmy,dmy,dmy,dmy,
                           &matId,&n,&finite,&sss,&isw,dmyI);
  
    else if (mDim == 2) matlib2d_(matDat,dmy,dmy,dmy,dmy,dmy,dmy,dmy,
                                  &matId,&n,&finite,&sss,&isw,dmyI);
    
      else if (mDim == 3) matlib3d_(matDat,dmy,dmy,dmy,dmy,dmy,dmy,
                                    &matId,&n,&finite,&isw,dmyI);
      
        else prgError(1,"ElementSolid::nivGP","invalid value of ndm!");

  return n;
}






void ElementSolid::initialiseIntVar(void)
{
  // check or set nivGP of element group
	
  int n = nivGP(), *eGnivGP = &(((ElementGroup*) elemGrp)->nivGP);

  if (*eGnivGP == 0) *eGnivGP = n;
  
    else if (*eGnivGP != n) 
      prgError(1,"ElementSolid::initialiseIntVar","nivGP inconsistency!");
	
  // check or set nGaussPoints of element group
	
  n = nGaussPoints();
   
  int *eGnGP = &(((ElementGroup*) elemGrp)->nGP);

  if (*eGnGP == 0) *eGnGP = n;
  
    else if (*eGnGP != n) 
      prgError(2,"ElementSolid::initialiseIntVar","nGP inconsistency!");
	
  // allocate memory

  n = (*eGnivGP) * (*eGnGP);
    
  intVar1 = new double [n];
  intVar2 = new double [n];	
 
  // initialise values

  double dmy[10], *matDat = &(((ElementGroup*) elemGrp)->elemProp[MATERIAL]->data[0]);

  int    l, ll = 0, isw = 2, dmyI[3],
         matId  = ((ElementGroup*) elemGrp)->elemProp[MATERIAL]->id + 1,
         mDim   = matdim_(&matId),
         finite = finiteStrain(),
         sss    = stressStrainState();

  //cout << '\t' << "  mDim  : " << mDim << endl;

  for (l=0; l<(*eGnGP); l++)
  {
    if (mDim == 1) matlib1d_(matDat,dmy,dmy,dmy,dmy,&(intVar1[ll]),&(intVar2[ll]),dmy,
	                     &matId,eGnivGP,&finite,&sss,&isw,dmyI);
  
      else if (mDim == 2) matlib2d_(matDat,dmy,dmy,dmy,dmy,&(intVar1[ll]),&(intVar2[ll]),dmy,
  		                    &matId,eGnivGP,&finite,&sss,&isw,dmyI);
    
        else if (mDim == 3) matlib3d_(matDat,dmy,dmy,dmy,&(intVar1[ll]),&(intVar2[ll]),dmy,
			              &matId,eGnivGP,&finite,&isw,dmyI);
      
          else prgError(1,"ElementSolid::initialiseIntVar","invalid value of ndm!");
	  
    ll += *eGnivGP;
  }

  //for (int i=0; i<n; i++)  cout << intVar1[i] << ","<< intVar2[i] << "+\n"; 

  return;
}








void ElementSolid::setGaussPointDataId(void)
{
  int ngp = nGaussPoints(),
      counter = GpDatId[0]-1, i;
  
  for (i=0; i<ngp; i++) GpDatId[i] = counter*ngp+i;

  return;
}

