
#include "Debug.h"
#include "FunctionsProgram.h"
#include "PropertyTypeEnum.h"
#include "NurbsElement.h"
#include <iostream>
#include <iomanip>
#include <stdio.h>
////#include "Plot.h"
#include "MpapTime.h"
#include "NurbsShapeFunctions.h"
#include "TimeFunction.h"

extern MpapTime           mpapTime;
extern List<TimeFunction> timeFunction;

using namespace std;



NurbsElement::NurbsElement(void)
{
  if (debug) cout << " constructor NurbsElement\n\n";

  // cout << "     NurbsElement: constructor ...\n\n";

  nlbf = ndof = nsize = nivGP = nGP = elenum = patchnum = counter = 0;
  
  tracflag = false;

  curve0   = NULL;
  curve1   = NULL;

  surf0    = NULL;
  surf1    = NULL;
  surf2    = NULL;

  solid0   = NULL;
  solid1   = NULL;
  solid2   = NULL;

  intVar1  = NULL;
  intVar2  = NULL;
  elmDat   = NULL;
  matDat   = NULL;

}





NurbsElement::~NurbsElement()
{
  if (debug)   cout << " destructor NurbsElement\n\n";

  if (intVar1!=NULL) delete [] intVar1;   intVar1 = NULL;
  if (intVar2!=NULL) delete [] intVar2;   intVar2 = NULL;

  if (elmDat!=NULL) elmDat = NULL;
  if (matDat!=NULL) matDat = NULL;

  if(curve0 != NULL)  // dont delete this pointer as it is just a reference
     curve0 = NULL;

  if(curve1 != NULL)  // dont delete this pointer as it is just a reference
     curve1 = NULL;

  if(surf0 != NULL)  // dont delete this pointer as it is just a reference
     surf0 = NULL;

  if(surf1 != NULL)  // dont delete this pointer as it is just a reference
     surf1 = NULL;

  if(surf2 != NULL)  // dont delete this pointer as it is just a reference
     surf2 = NULL;

  if(solid0 != NULL)
     solid0 = NULL;

  if(solid1 != NULL)
    solid1 = NULL;

  if(solid2 != NULL)
    solid2 = NULL;

//  cout << "     NurbsElement: destructor ...\n\n";
}



void NurbsElement::prepareElemData()
{
    // no. of local basis functions
    if(curve0 != NULL)
       nlbf = curve0->p + 1;

    if(surf0 != NULL)
       nlbf = surf0->nlbf;

    if(solid0 != NULL)
       nlbf = solid0->nlbf;


    nsize = nlbf*ndof;

    int ii, jj;

    forassy.setDim(nsize);

    forassembly.setDim(nsize);
    for(ii=0;ii<nsize;ii++)
    {
      forassembly[ii].setDim(nsize);
      for(jj=0;jj<nsize;jj++)
        forassembly[ii][jj] = -1;
    }

    primvar.setDim(nsize);

    stiffness_local.setDim(nsize);
    for(ii=0;ii<nsize;ii++)
    {
      stiffness_local[ii].setDim(nsize);
      stiffness_local[ii].zero();
    }

    //resi.setDim(nsize);
    
    Klocal.resize(nsize, nsize);
    Flocal.resize(nsize);

  return;
}



void NurbsElement::AssembleElementMatrix(int index, MatrixSparseArray<double>& mtx)
{
   int nn=0, aa, bb;
     
   for(aa=0;aa<nsize;aa++)
   {
     for(bb=0;bb<nsize;bb++)
     {
       nn = forassembly[aa][bb];
       if(nn != -1)
       {
         mtx.x[nn-1] += Klocal(aa,bb);
       }
     }
   }

   return;
}




void NurbsElement::AssembleElementMatrix(int index, SparseMatrixXd& mtx)
{
   cout << " Need to be implemented " << endl;
   int nn=0, aa, bb;
     
   for(aa=0;aa<nsize;aa++)
   {
     for(bb=0;bb<nsize;bb++)
     {
       nn = forassembly[aa][bb];
       if(nn != -1)
       {
         //mtx.coeffRef(aa,bb) += Klocal(aa,bb);
       }
     }
   }

   return;
}




void NurbsElement::AssembleElementMatrix3(int index, double fact, SparseMatrixXd& mtx)
{
  cout << " NurbsElement::AssembleElementMatrix3 ... need to correct this subroutine " << endl;

   int nn=0, aa, bb;
     
   for(aa=0;aa<nsize;aa++)
   {
     for(bb=0;bb<nsize;bb++)
     {
       nn = forassembly[aa][bb];
       if(nn != -1)
         mtx.coeffRef(aa,bb) += (fact * stiffness_local[aa][bb]);
     }
   }

  return;
}



void NurbsElement::AssembleElementVector(bool firstIter, bool flag, double* rhs, double* reac, int start1, int start2)
{
   // flag == true  ---> just external force vector
   // flag == false ---> internal load vector + contributions from nodes with specified displacement BCs

  //cout << " primvar " << endl;

   int *tt;

   if(curve0 != NULL)
     tt = &(curve0->LM[elenum][0]);

   if(surf0 != NULL)
     tt = &(surf0->LM[elenum][0]);

   if(solid0 != NULL)
     tt = &(solid0->LM[elenum][0]);
   
   //cout << surf0->LM[elenum] << endl;

   if(flag)
   {
      for(int aa=0;aa<nsize;aa++)
      {
         if(tt[aa] != -1)
            rhs[tt[aa]] += Flocal[aa];
      }
   }
   else
   {
      int aa, bb;
      for(aa=0;aa<nsize;aa++)
      {
         if(tt[aa] != -1)
           rhs[tt[aa]] += Flocal[aa];

         // add up reaction forces
         reac[forassy[aa]] += Flocal[aa];
      }
      if(firstIter)
      {
         double fact1;

         //printPrimVariable();
         //printf("\n\n");

         //fact1 = timeFunction[0].prop;
         fact1 = mpapTime.dt;
         //fact1 = 1.0;

         double fact;
         for(aa=0;aa<nsize;aa++)
         {
            if(tt[aa] == -1)
            {
               fact = fact1 * primvar[aa];
               //cout << " fact " << fact << endl;
               for(bb=0;bb<nsize;bb++)
               {
                  if(tt[bb] != -1)
                  {
                    //cout << " fact " << Klocal(bb, aa) * fact << endl;
                    rhs[tt[bb]] -= Klocal(bb, aa) * fact;
                  }
               }
            }
         }
      }
   }

//  cout << " Flocal " << Flocal << endl;

  tt= NULL;

  return;
}



void NurbsElement::printStiffnessMatrix()
{
   int ii, jj;
   printf(" Stiffness Matrix for element : %5d\n\n",elenum);
   for(ii=0;ii<nsize;ii++)
   {
      for(jj=0;jj<nsize;jj++)
      {
         printf("%12.6f\t", Klocal(ii,jj)) ;
      }
      printf("\n");
   }
   printf("\n");
   printf("\n");

   return;
}




void NurbsElement::printForceVector()
{
   printf(" Residue Vector for element : %5d\n\n",elenum);
   for(int ii=0;ii<nsize;ii++)
   {
     printf("\t%5d\t%14.12f\n",ii, Flocal[ii]);
   }
   printf("\n");
   printf("\n");

  return;
}






void NurbsElement::printPrimVariable()
{
  int index, ii, dof;
  for(ii=0;ii<nlbf;ii++)
  {
    index = ndof*ii;
    printf("\t %3d \t", ii);
    for(dof=0;dof<ndof;dof++)
      printf(" %12.8f \t", primvar[index+dof] );
    printf("\n");
  }

  return;
}




