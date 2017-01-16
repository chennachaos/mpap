
#include <math.h>
#include "Debug.h"
#include "MpapTime.h"
#include "Plot.h"
#include "NurbsElemKirchhoffPlate.h"
#include "NurbsShapeFunctions.h"
#include <assert.h>
#include "ComputerTime.h"

using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern Plot plot;


NurbsElemKirchhoffPlate::NurbsElemKirchhoffPlate()
{
  if (debug) cout << " constructor NurbsElemKirchhoffPlate\n\n";
}



NurbsElemKirchhoffPlate::~NurbsElemKirchhoffPlate()
{
  if (debug) cout << " destructor NurbsElemKirchhoffPlate\n\n";
}



int NurbsElemKirchhoffPlate::calcStiffnessAndResidual()
{

//  cout << '\t' << " Total Volume for element # " << elenum << " is = " << totvol << endl; cout << endl;

//   computerTime.stopAndPrint(fct);

  return 0;
}






void NurbsElemKirchhoffPlate::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{

  return;
}



void NurbsElemKirchhoffPlate::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{


  return;
}




void NurbsElemKirchhoffPlate::projectStress(int varindex, double* outval)
{
 
  return;
}




void NurbsElemKirchhoffPlate::projectStrain(int vartype, int varindex, double* outval)
{

  return;
}



void NurbsElemKirchhoffPlate::projectIntVar(int index, double* outval)
{
   int ind1, ii, jj;

   ind1 = 0;
   for(jj=0;jj<nGP2;jj++)
   {
       for(ii=0;ii<nGP1;ii++)
       {
           outval[ind1] = intVar2[ind1*nivGP+index];
           ind1++;
       }
   }


   return;
}



int NurbsElemKirchhoffPlate::calcStiffnessMatrix(double dt)
{

  return 0;
}





int NurbsElemKirchhoffPlate::calcInternalForces()
{
 

  return 0;
}





int NurbsElemKirchhoffPlate::calcMassMatrix(int lumpInd, double dt)
{
 
  return 0;
}





int NurbsElemKirchhoffPlate::calcOutput(double u1, double v1)
{

  return 0;
}

