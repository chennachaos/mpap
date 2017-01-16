
#include "LagrangeElement.h"

#include "Debug.h"
#include "PropertyTypeEnum.h"
#include "MpapTime.h"
#include "TimeFunction.h"

#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "Functions.h"
#include "QuadratureUtil.h"
#include "FunctionsMaterial.h"

#include "KimMoinFlow.h"

extern MpapTime           mpapTime;
extern List<TimeFunction> timeFunction;

using namespace std;



LagrangeElement::LagrangeElement()
{
  if (debug) cout << " constructor LagrangeElement\n\n";

  // cout << "     LagrangeElement: constructor ...\n\n";

  //nlbf = ndof = nsize = nivGP = nGP = elenum = patchnum = counter = 0;
  subdomId = 0;

  tracflag = false;

  intVar1  = NULL;
  intVar2  = NULL;
  elmDat   = NULL;
  matDat   = NULL;

}





LagrangeElement::~LagrangeElement()
{
  if (debug)   cout << " destructor LagrangeElement\n\n";

  if (intVar1!=NULL) delete [] intVar1;   intVar1 = NULL;
  if (intVar2!=NULL) delete [] intVar2;   intVar2 = NULL;

  if (elmDat!=NULL) elmDat = NULL;
  if (matDat!=NULL) matDat = NULL;

//  cout << "     LagrangeElement: destructor ...\n\n";
}




void LagrangeElement::prepareElemData()
{
    int ii, jj, kk, ind;

    nlbf = nodeNums.size();
    nsize = nlbf*ndof;

    //printVector(nodeNums);
    globalDOFnums.resize(nsize);

    kk=0;
    for(ii=0;ii<nodeNums.size();ii++)
    {
      ind = nodeNums[ii]*ndof;
      for(jj=0;jj<ndof;jj++)
        globalDOFnums[kk++] = ind+jj;
    }

    //primvar.setDim(nsize);
    //resi.setDim(nsize);
    //printVector(globalDOFnums);

  return;
}


double  LagrangeElement::computeGeomOrig(int dir, VectorXd& NN)
{
  double val=0.0;

  for(int ii=0;ii<nlbf;ii++)
    val += GeomData->NodePosOrig[nodeNums[ii]][dir] * NN[ii];

  return  val;
}


double  LagrangeElement::computeGeomNew(int dir, VectorXd& NN)
{
  double val=0.0;

  for(int ii=0;ii<nlbf;ii++)
    val += GeomData->NodePosNew[nodeNums[ii]][dir] * NN[ii];

  return  val;
}



double  LagrangeElement::computeGeomCur(int dir, VectorXd& NN)
{
  double val=0.0;

  for(int ii=0;ii<nlbf;ii++)
    val += GeomData->NodePosCur[nodeNums[ii]][dir] * NN[ii];

  return  val;
}



double  LagrangeElement::computeValue(int dir, VectorXd& NN)
{
   double val=0.0;

   for(int ii=0;ii<nlbf;ii++)
     val += SolnData->var1[nodeNums[ii]*ndof+dir] * NN[ii];
   
   return  val;
}



double  LagrangeElement::computeValueExtrap(int dir, VectorXd& NN)
{
  double val=0.0;

  for(int ii=0;ii<nlbf;ii++)
    val += SolnData->var1Extrap[nodeNums[ii]*ndof+dir] * NN[ii];

  return  val;
}




double  LagrangeElement::computeValuePrev(int dir, VectorXd& NN)
{
   double val=0.0;

   for(int ii=0;ii<nlbf;ii++)
     val += SolnData->var1Prev[nodeNums[ii]*ndof+dir] * NN[ii];
   
   return  val;
}




double  LagrangeElement::computeValuePrev2(int dir, VectorXd& NN)
{
   double val=0.0;

   for(int ii=0;ii<nlbf;ii++)
     val += SolnData->var1Prev2[nodeNums[ii]*ndof+dir] * NN[ii];
   
   return  val;
}






double  LagrangeElement::computeValuePrev3(int dir, VectorXd& NN)
{
  double val=0.0;

  for(int ii=0;ii<nlbf;ii++)
    val += SolnData->var1Prev3[nodeNums[ii]*ndof+dir] * NN[ii];

  return  val;
}






double  LagrangeElement::computeValuePrev4(int dir, VectorXd& NN)
{
  double val=0.0;

  for(int ii=0;ii<nlbf;ii++)
    val += SolnData->var1Prev4[nodeNums[ii]*ndof+dir] * NN[ii];

  return  val;
}






double  LagrangeElement::computeValueCur(int dir, VectorXd& NN)
{
   double val=0.0;

   for(int ii=0;ii<nlbf;ii++)
     val += SolnData->var1Cur[nodeNums[ii]*ndof+dir] * NN[ii];
   
   return  val;
}




double  LagrangeElement::computeValueDot(int dir, VectorXd& NN)
{
   double val=0.0;

   for(int ii=0;ii<nlbf;ii++)
     val += SolnData->var1Dot[nodeNums[ii]*ndof+dir] * NN[ii];
   
   return  val;
}



double  LagrangeElement::computeValueDotDot(int dir, VectorXd& NN)
{
   double val=0.0;

   for(int ii=0;ii<nlbf;ii++)
     val += SolnData->var1DotDot[nodeNums[ii]*ndof+dir] * NN[ii];
   
   return  val;
}





double  LagrangeElement::computeValueDotCur(int dir, VectorXd& NN)
{
   double val=0.0;

   for(int ii=0;ii<nlbf;ii++)
     val += SolnData->var1DotCur[nodeNums[ii]*ndof+dir] * NN[ii];
   
   return  val;
}


double  LagrangeElement::computeValueDotDotCur(int dir, VectorXd& NN)
{
   double val=0.0;

   for(int ii=0;ii<nlbf;ii++)
     val += SolnData->var1DotDotCur[nodeNums[ii]*ndof+dir] * NN[ii];
   
   return  val;
}





double  LagrangeElement::computeValue2(int dir, VectorXd& NN)
{
   double val=0.0;

   for(int ii=0;ii<nlbf;ii++)
     val += SolnData->var2[nodeNums[ii]*ndof+dir] * NN[ii];
   
   return  val;
}



double  LagrangeElement::computeValue2Prev(int dir, VectorXd& NN)
{
   double val=0.0;

   for(int ii=0;ii<nlbf;ii++)
     val += SolnData->var2Prev[nodeNums[ii]*ndof+dir] * NN[ii];
   
   return  val;
}



double  LagrangeElement::computeValue2Cur(int dir, VectorXd& NN)
{
   double val=0.0;

   for(int ii=0;ii<nlbf;ii++)
     val += SolnData->var2Cur[nodeNums[ii]*ndof+dir] * NN[ii];
   
   return  val;
}


double  LagrangeElement::computeForce(int dir, VectorXd& NN)
{
  double val=0.0;

  return  val;
}



double  LagrangeElement::computeForcePrev(int dir, VectorXd& NN)
{
  double val=0.0;

  return  val;
}



double  LagrangeElement::computeForceCur(int dir, VectorXd& NN)
{
  double val=0.0;

  return  val;
}



/*
void  LagrangeElement::AssembleElementVector(int ind, bool flag, double* rhs, double* reac, int start1, int start2)
{
  for(int ii=0;ii<nsize;ii++)
  {
    if(forAssyVec[ii] != -1)
      rhs[start1+forAssyVec[ii]] += Flocal(ii);
  }

  return;
}
*/


void  LagrangeElement::AssembleElementMatrix(int start, SparseMatrixXd& mtx)
{
/*
  int aa, bb, ii, jj, r;

  for(ii=0;ii<nsize;ii++)
  {
    aa = forAssyVec[ii];
    if( aa != -1 )
    {
      r = start + aa;
      for(jj=0;jj<nsize;jj++)
      {
        bb = forAssyVec[jj];
        if( bb != -1 )
          mtx.coeffRef(r, start+bb) += Klocal(ii,jj);
      }
    }
  }
*/
  return;
}


void  LagrangeElement::AssembleElementMatrix(int start, Mat mtx)
{
/*
  PetscErrorCode ierr;
  int aa, bb, ii, jj, r;

  for(ii=0;ii<nsize;ii++)
  {
    aa = forAssyVec[ii];
    if( aa != -1 )
    {
      r = start + aa;
      for(jj=0;jj<nsize;jj++)
      {
        bb = forAssyVec[jj];
        if( bb != -1 )
          ierr = MatSetValue(mtx, aa, bb, Klocal(ii, jj), ADD_VALUES);
      }
    }
  }
*/
  return;
}



void LagrangeElement::AssembleElementMatrixAndVector(int start, SparseMatrixXd& mtx, double* rhs)
{
/*
  int aa, bb, ii, jj, r;

  for(ii=0;ii<nsize;ii++)
  {
    aa = forAssyVec[ii];
    if( aa != -1 )
    {
      r = start + aa;
      rhs[r] += Flocal(ii);
      for(jj=0;jj<nsize;jj++)
      {
        bb = forAssyVec[jj];
        if( bb != -1 )
          mtx.coeffRef(r, start+bb) += Klocal(ii,jj);
      }
    }
  }
*/
  return;
}





void LagrangeElement::AssembleElementVector(bool firstIter, bool flag, double* rhs, double* reac, int start1, int start2)
{
  // flag == true  ---> just external force vector
  // flag == false ---> internal load vector + contributions from nodes with specified displacement BCs
/*
  int aa, bb, ii, jj, r, c;

  if(flag)
  {
      for(ii=0;ii<nsize;ii++)
      {
        aa = forAssyVec[ii];
        r = start1 + aa;
        if( aa != -1)
          rhs[r] += Flocal[ii];
      }
  }
  else
  {
      //printVector(forAssyVec);
      for(ii=0;ii<nsize;ii++)
      {
        aa = forAssyVec[ii];
        r = start1 + aa;

        if(aa != -1)
          rhs[r] += Flocal[ii];

        //cout << " vvvvvvvvvvvvvvvv " << endl;
        // add up reaction forces
        reac[globalDOFnums[ii]] += Flocal[ii];
      }
      if(firstIter)
      {
         //cout << " vvvvvvvvvvvvvvvv " << endl;
         double fact1, fact;
         
         fact1 = timeFunction[0].prop;
         //fact1 = mpapTime.dt;
         fact1 = 1.0;

         //if(mpapTime.cur <= 5.0e-3)
         //{
           //fact1 = 0.5*( 1.0-cos(628.3185*mpapTime.cur));
           //fact1 = 0.5*( 1.0-cos(20.0*mpapTime.cur));
           //fact1 = 0.5*( 1.0-cos(20.0*mpapTime.cur)) - 0.5*( 1.0-cos(20.0*(mpapTime.cur-mpapTime.dt)));
           //fact = sin(628.3185*mpapTime.cur);
           //fact1 = 1.0;
           //fact = sin(20.0*mpapTime.cur);
         //}
         //else
           //fact1 = 0.0;

         for(ii=0;ii<nsize;ii++)
         {
            aa = forAssyVec[ii];
            if(aa == -1)
            {
               //cout << " aa " << aa << endl;
               fact = fact1 * SolnData->var1applied[globalDOFnums[ii]];
               //cout << " fact " << fact << endl;

               for(jj=0;jj<nsize;jj++)
               {
                  bb = forAssyVec[jj];
                  r = start1 + bb;
                  if( bb != -1 )
                    rhs[r] -= Klocal(jj, ii) * fact;
               }
            }
         }
      }
  }
*/
//  cout << " resi " << resi << endl;

  return;
}



void LagrangeElement::AssembleElementVector(bool firstIter, bool flag, Vec rhs, Vec reac, int start1, int start2)
{
/*
  // flag == true  ---> just external force vector
  // flag == false ---> internal load vector + contributions from nodes with specified displacement BCs

  int aa, bb, ii, jj, r, c;

  if(flag)
  {
      for(ii=0;ii<nsize;ii++)
      {
        aa = forAssyVec[ii];
        r = start1 + aa;
        if( aa != -1)
          VecSetValue(rhs, r, Flocal(ii), ADD_VALUES);
      }
  }
  else
  {
      //printVector(forAssyVec);
      for(ii=0;ii<nsize;ii++)
      {
        aa = forAssyVec[ii];
        r = start1 + aa;

        if(aa != -1)
          VecSetValue(rhs, r, Flocal(ii), ADD_VALUES);

        // add up reaction forces
        //reac[globalDOFnums[ii]] += Flocal[ii];
        VecSetValue(reac, globalDOFnums[ii], Flocal(ii), ADD_VALUES);
      }
      if(firstIter)
      {
         //cout << " aaaaaaaaaaaaaa " << endl;
         double fact1=0.0, val=0.0;

         fact1 = timeFunction[0].prop;
         //fact1 = mpapTime.dt;

         double fact;
         for(ii=0;ii<nsize;ii++)
         {
            aa = forAssyVec[ii];
            if(aa == -1)
            {
               //cout << " aa " << aa << endl;
               fact = fact1 * SolnData->var1applied[globalDOFnums[ii]];
               //cout << " fact " << fact << endl;

               for(jj=0;jj<nsize;jj++)
               {
                  bb = forAssyVec[jj];
                  r = start1 + bb;
                  if( bb != -1 )
		  {
                    val = -Klocal(jj, ii) * fact;
                    VecSetValue(rhs, r, val, ADD_VALUES);
		  }
               }
            }
         }
      }
  }
*/
//  cout << " resi " << resi << endl;

  return;
}



void LagrangeElement::printStiffnessMatrix()
{
/*
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
*/

  return;
}




void LagrangeElement::printForceVector()
{
   printf(" Residue Vector for element : %5d\n\n",elenum);
   for(int ii=0;ii<nsize;ii++)
   {
     printf("\t%5d\t%12.8f\n",ii,resi[ii]);
   }
   printf("\n");
   printf("\n");

  return;
}



void LagrangeElement::printPrimVariable()
{
     int index, ii, dof;
     for(ii=0;ii<nlbf;ii++)
     {
        index = ndof*ii;
        cout << '\t' << ii << '\t' ;
        for(dof=0;dof<ndof;dof++)
           cout << primvar[index+dof] << '\t' ;
        cout << endl;
     }

  return;
}



int  LagrangeElement::calcError(int index)
{
    //cout << " Error ... " << endl;

    // computing error for Navier-Stokes
    ///////////////////////////////////////////////////////////
    
    //Stokes2DEx1  analy;
    //Stokes2DEx2  analy;

    //Stokes2DEx3  analy;
    //Kovasznay  analy;
    //analy.SetPressure(0.0);

    //PearsonVortex analy;

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    //matDat = &(SolnData->MatlProp[matType].data[0]);

    double rho = elmDat[4];
    double mu  = elmDat[5];

    KimMoinFlow  analy(rho, mu, 1.0);
    //Kovasznay  analy;

    int      ii, jj, gp1, gp2, TI;
    double   Jac, dvol, diff, fact, val;
    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);
    double  param[2], geom[2];

    elemError = 0.0;
    


 if(index < 3) // L2 norm in x-velocity (index=0), y-velocity (index=1) and pressure (index=2)
 {
    for(gp2=0;gp2<nGP;gp2++)
    {
        param[1] = GeomData->gausspoints2[gp2];

    for(gp1=0;gp1<nGP;gp1++)
    {
          param[0] = GeomData->gausspoints1[gp1];

          GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = GeomData->gaussweights2[gp2] * GeomData->gaussweights1[gp1] * Jac;

          geom[0] = computeGeomOrig(0, N);
          geom[1] = computeGeomOrig(1, N);

          val = analy.computeValue(index, geom[0], geom[1], mpapTime.cur);

          val -= computeValue(index, N);

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f  \t %12.8f   \t %12.8f \n", vx, vx2, vy, vy2, diff);
          
          elemError += ( (val*val) * dvol );
    }//gp1
    }//gp2
 }
 else // H1 semi-norm in velocity
 {
    MatrixXd  F(2,2);
    double  v[2];

    for(gp2=0;gp2<nGP;gp2++)
    {
        param[1] = GeomData->gausspoints2[gp2];

    for(gp1=0;gp1<nGP;gp1++)
    {
          param[0] = GeomData->gausspoints1[gp1];

          GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = GeomData->gaussweights2[gp2] * GeomData->gaussweights1[gp1] * Jac;


          geom[0] = computeGeomOrig(0, N);
          geom[1] = computeGeomOrig(1, N);

          v[0] = analy.computeValue(0, geom[0], geom[1], mpapTime.cur);
          v[1] = analy.computeValue(1, geom[0], geom[1], mpapTime.cur);

          v[0] -= computeValue(0, N);
          v[1] -= computeValue(1, N);

          F.setZero();
          analy.computeDerivatives(geom[0], geom[1], mpapTime.cur, &(F(0,0)));

          F(0,0) -= computeValue(0, dN_dx);
          F(0,1) -= computeValue(0, dN_dy);
          F(1,0) -= computeValue(1, dN_dx);
          F(1,1) -= computeValue(1, dN_dy);

          val = v[0]*v[0]+v[1]*v[1] + F(0,0)*F(0,0)+F(0,1)*F(0,1)+F(1,0)*F(1,0)+F(1,1)*F(1,1);

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", dx, dy, diff);
        
          elemError += ( val * dvol );
    }//gp1
    }//gp2
 }

  return 1;
}



void  LagrangeElement::computeMomentum(int index1, int index2, VectorXd& momentum)
{
  // compute total energy for structural dynamic problems
  ///////////////////////////////////////////////////////////

  elmDat = &(SolnData->ElemProp[elmType].data[0]);
  //matDat = &(SolnData->MatlProp[matType].data[0]);

  double  F[4], detF=0.0, F33, fact, dvol, dvol0, Jac, dt, JacMult;
  double  posi[2], disp[2], velo[2], acce[2];
  double  stre[4], cc[4][4], bc[2][4], param[2], bforce[2];

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);

  int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp1, gp2;
  int   ind1, ind2, kk, DIM=2;
  
  double  rho  = elmDat[5] ;

  count = 1;   ll = 0;   err = 0;   isw = 3;
  dt = mpapTime.dt;

  //cout << " finite = " << finite << endl;

  int nGP = (int) elmDat[0] ;

  vector<double>  gausspoints1, gaussweights1;

  getGaussPoints1D(nGP, gausspoints1, gaussweights1);

  momentum.setZero();

  for(gp2=0;gp2<nGP;gp2++)
  {
    JacMult = gaussweights1[gp2] * thick;

    param[1] = gausspoints1[gp2];

  for(gp1=0;gp1<nGP;gp1++)
  {
        param[0] = gausspoints1[gp1];

        GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        //for(ii=0; ii<nlbf; ii++)
          //cout << N[ii] << '\t' << dN_dx[ii] << '\t' << dN_dy[ii] << endl;

        fact = gaussweights1[gp1] * JacMult;

        dvol0 = Jac * fact;
        dvol = dvol0;

        //if(axsy)
          //dvol *= 2.0*PI*yy;

        //if(finite)
        //{
          //GeomData->computeBasisFunctions2D(1, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          //dvol = Jac * fact;
        //}

        disp[0] = computeValue(0, N);
        disp[1] = computeValue(1, N);

        posi[0] = computeGeomOrig(0, N);
        posi[1] = computeGeomOrig(1, N);

        posi[0] += disp[0];
        posi[1] += disp[1];

        velo[0] = computeValueDot(0, N);
        velo[1] = computeValueDot(1, N);

        acce[0] = computeValueDotDot(0, N);
        acce[1] = computeValueDotDot(1, N);

        //printf(" %14.12f \t %14.12f \n ", posi[0], posi[1]);
        //printf(" %14.12f \t %14.12f \n ", velo[0], velo[1]);

        fact = rho*dvol0;
        momentum[0] += fact*velo[0];
        momentum[1] += fact*velo[1];

        fact = rho*dvol0;
        momentum[5] += fact*(posi[0]*velo[1]-posi[1]*velo[0]);
  }//gp1
  }//gp2

  //printf(" \n \n ");

  return;
}

  
  
  
  

  
  
  
  
  
  
  
  