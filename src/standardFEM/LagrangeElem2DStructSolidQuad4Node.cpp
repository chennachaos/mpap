
#include "LagrangeElem2DStructSolidQuad4Node.h"
#include "Debug.h"
#include "MpapTime.h"
#include "ComputerTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "FunctionsMaterial.h"
#include "QuadratureUtil.h"
#include "TimeFunction.h"

using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;



LagrangeElem2DStructSolidQuad4Node::LagrangeElem2DStructSolidQuad4Node()
{
  if (debug) cout << " constructor LagrangeElem2DStructSolidQuad4Node\n\n";

  ndim   = 2;
  degree = 1;
  npElem = 4;
  nlbf   = 4;
  ndof   = 2;
  nsize  = npElem*ndof;
}


LagrangeElem2DStructSolidQuad4Node::~LagrangeElem2DStructSolidQuad4Node()
{
  if (debug) cout << " destructor LagrangeElem2DStructSolidQuad4Node\n\n";
}


void LagrangeElem2DStructSolidQuad4Node::prepareElemData()
{
  LagrangeElement::prepareElemData();

  return;
}



void LagrangeElem2DStructSolidQuad4Node::prepareElemData2()
{
  cout << " LagrangeElem2DStructSolidQuad4Node::prepareElemData2() " << endl;
  return;
}


int LagrangeElem2DStructSolidQuad4Node::calcLoadVector()
{
  return 0;
}



int LagrangeElem2DStructSolidQuad4Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
    double  F[4], detF, F33, fact, dvol, dvol0;
    double  Jac, bb1, bb2, bb3, cc1, cc2, cc3;

    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);
    VectorXd  dispC(nsize), velC(nsize), accC(nsize);
    MatrixXd  Mlocal(nsize, nsize);

    double  stre[4], cc[4][4], bc[2][4], param[2], bforce[2];

    int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp, TI, TIp1, TJ, TJp1;
    int   ind1, ind2, kk;

    double rho  = elmDat[5] ;
    double rho0 = rho ;
    //bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    //bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    bforce[0]   = rho0*elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = rho0*elmDat[7]*timeFunction[0].prop ;
    //bforce[0]   = rho0*elmDat[6] ;
    //bforce[1]   = rho0*elmDat[7] ;
    double af = SolnData->td(2);
    double d1 = SolnData->td(5);
    double acceFact = SolnData->td(10);
    double dt = mpapTime.dt;

    //cout << "finite = " << finite << '\t' << af << '\t' << d1 << '\t' << acceFact << endl;
    //cout << "timeFunction[0].prop=" << timeFunction[0].prop << endl;
    //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", BULK, mu, bforce[0], bforce[1]);
    //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", rho, af, d1, aa);
    //printf(" %5d \t %5d \t %5d \t %5d \n ", nodeNums[0], nodeNums[1], nodeNums[2], nodeNums[3]);

    double xNode[4], yNode[4], xx, yy;

    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    for(ii=0;ii<npElem;ii++)
    {
      ind1 = ndof*ii;
      ind2 = nodeNums[ii]*ndof;

      for(kk=0;kk<ndof;kk++)
      {
        dispC(ind1+kk)  =  SolnData->var1Cur[ind2+kk];
        velC(ind1+kk)   =  SolnData->var1DotCur[ind2+kk];
        accC(ind1+kk)   =  SolnData->var1DotDotCur[ind2+kk];
      }
    }

    //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", xNode[0], xNode[1], xNode[2], xNode[3]);
    //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n\n\n ", yNode[0], yNode[1], yNode[2], yNode[3]);

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();
    Mlocal.setZero();

    vector<double>  gausspoints1, gausspoints2, gaussweights;

    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(thick*Jac);
        dvol  = dvol0;


//         xx = yy= 0.0;
//         for(ii=0;ii<nlbf;ii++)
//         {
//           xx += N[ii]*xNode[ii];
//           yy += N[ii]*yNode[ii];
//         }


        //for(ii=0; ii<nlbf; ii++)
          //cout << N[ii] << '\t' << dN_dx[ii] << '\t' << dN_dy[ii] << endl;

        if(axsy)
          dvol *= 2.0*PI*yy;

        //GeomData->computeDeformationGradient(1, nodeNums, &dN_dx(0), &dN_dy(0), F, detF);

        F[0] = computeValueCur(0, dN_dx) + 1.0;
        F[2] = computeValueCur(0, dN_dy);
        F[1] = computeValueCur(1, dN_dx);
        F[3] = computeValueCur(1, dN_dy) + 1.0;

        detF =  F[0]*F[3] - F[1]*F[2];

        if(detF < 0.0)   return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(1, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          dvol = gaussweights[gp]*(thick*Jac);
        }

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
          F33 = 1.0;

        //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", F[0], F[1], F[2], F[3]);

        matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        //eps[0] = F[0] - 1.0;   eps[1] = F[3] - 1.0;  eps[2] = F33  - 1.0;   eps[3] = 0.5 * (F[1] + F[2]);
        //fact   = lamb * (eps[0] + eps[1] + eps[2]);        fact2  = 2.0 * mu;
        //stre[0] = fact + fact2 * eps[0];   stre[1] = fact + fact2 * eps[1];
        //stre[2] = fact + fact2 * eps[2];   stre[3] =        fact2 * eps[3];

        //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", eps[0], eps[1], eps[2], eps[3]);
        //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", stre[0], stre[1], stre[2], stre[3]);
        //printf("\n %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", cc[0][0], cc[0][1], cc[0][2], cc[0][3]);
        //printf("   %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", cc[1][0], cc[1][1], cc[1][2], cc[1][3]);
        //printf("   %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", cc[2][0], cc[2][1], cc[2][2], cc[2][3]);
        //printf("   %14.12f \t %14.12f \t %14.12f \t %14.12f \n\n ", cc[3][0], cc[3][1], cc[3][2], cc[3][3]);
        //printf("  dvol =  %14.12f \n\n ", dvol);
        //printf(" %14.12f \t %14.12f \n ", bforce[0], bforce[1]);

        if(err !=0)    return 1;

        if(finite)  dvol *= F33;

//         for(ii=0;ii<4;ii++)
//         {
//           stre[ii] *= dvol;
//           for(jj=0;jj<4;jj++)
//             cc[ii][jj] *= dvol;
//         }

        //==============================================
        // CALCULATE TANGENT STIFFNESS
        //==============================================

        //   part 1. -- material part (not necessarily symmetric!!)


        for(ii=0;ii<nlbf;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol0*(N[ii]*rho);

          bc[0][0] = bb1 * cc[0][0] + bb2 * cc[3][0];
          bc[0][1] = bb1 * cc[0][1] + bb2 * cc[3][1];
          bc[0][2] = bb1 * cc[0][3] + bb2 * cc[3][3];

          bc[1][0] = bb2 * cc[1][0] + bb1 * cc[3][0];
          bc[1][1] = bb2 * cc[1][1] + bb1 * cc[3][1];
          bc[1][2] = bb2 * cc[1][3] + bb1 * cc[3][3];

          TI   = 2*ii;
          TIp1 = TI+1;

          Flocal(TI)   += (bb3*bforce[0] - bb1*stre[0] - bb2*stre[3]) ;
          Flocal(TIp1) += (bb3*bforce[1] - bb1*stre[3] - bb2*stre[1]) ;

          for(jj=0; jj<nlbf; jj++)
          {
              cc1 = af*dN_dx[jj];
              cc2 = af*dN_dy[jj];

              TJ   = 2*jj;
              TJp1 = TJ+1;

              Klocal(TI,   TJ)    +=  (bc[0][0] * cc1 + bc[0][2] * cc2) ;
              Klocal(TI,   TJp1)  +=  (bc[0][1] * cc2 + bc[0][2] * cc1) ;
              Klocal(TIp1, TJ)    +=  (bc[1][0] * cc1 + bc[1][2] * cc2) ;
              Klocal(TIp1, TJp1)  +=  (bc[1][1] * cc2 + bc[1][2] * cc1) ;

              Mlocal(TI,   TJ)    +=  bb3*N[jj];
              Mlocal(TIp1, TJp1)  +=  bb3*N[jj];
          }
        }

        //   part 2. -- geometrical matrix  (if geometry nonlinear)

        if(finite)
        {
          for(ii=0; ii<nlbf; ii++)
          {
              bb1 = dvol*(dN_dx[ii]*stre[0] + dN_dy[ii]*stre[3]);
              bb2 = dvol*(dN_dx[ii]*stre[3] + dN_dy[ii]*stre[1]);

              TI   = 2*ii;
              TIp1 = TI+1;

              for(jj=0; jj<nlbf; jj++)
              {
                TJ   = 2*jj;
                TJp1 = TJ+1;

                fact = af*(bb1*dN_dx[jj] + bb2*dN_dy[jj] );

                Klocal(TI,TJ)     += fact ;
                Klocal(TIp1,TJp1) += fact ;
              }
          }
        }
    } //gp

    /*
    int vec[]={0,1,2,3,6,7,4,5};

    MatrixXd  K2(nsize, nsize);
    for(ii=0; ii<nsize; ii++)
    {
      for(jj=0; jj<nsize; jj++)
        K2(ii,jj) = Klocal(vec[ii], vec[jj]);
    }

    printMatrix(Klocal);
    printf("\n\n\n");
    printMatrix(K2);
    */

    Klocal +=  d1*Mlocal;
    Flocal -=  Mlocal*accC;

    //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);
    //printf("\n\n\n");  printVector(accC);

    return 0;
}


int LagrangeElem2DStructSolidQuad4Node::calcInternalForces()
{
  return 0;
}



void LagrangeElem2DStructSolidQuad4Node::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem2DStructSolidQuad4Node::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem2DStructSolidQuad4Node::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem2DStructSolidQuad4Node::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem2DStructSolidQuad4Node::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem2DStructSolidQuad4Node::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem2DStructSolidQuad4Node::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}




void  LagrangeElem2DStructSolidQuad4Node::computeEnergy(int index1, int index2, VectorXd& energy)
{
    // compute total energy for structural dynamic problems
    ///////////////////////////////////////////////////////////

    double  F[4], detF=0.0, F33, fact, dvol, dvol0, Jac;
    double  posi[2], disp[2], velo[2];
    double  stress[4], strain[4], cc[4][4], param[2], bforce[2];

    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);

    int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp, TI, TIp1, TJ, TJp1;
    int   ind1, ind2, kk;
  
    double  rho  = elmDat[5] ;

    count = 1;   ll = 0;   err = 0;   isw = 3;
    double dt = mpapTime.dt;


    vector<double>  gausspoints1, gausspoints2, gaussweights;

    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    energy.setZero();

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        fact  = gaussweights[gp] * thick;
        dvol0 = Jac * fact;
        dvol  = dvol0;

        //if(axsy)
          //dvol *= 2.0*PI*yy;

        //GeomData->computeDeformationGradient(1, nodeNums, &dN_dx(0), &dN_dy(0), F, detF);

        F[0] = computeValue(0, dN_dx) + 1.0;
        F[2] = computeValue(0, dN_dy);
        F[1] = computeValue(1, dN_dx);
        F[3] = computeValue(1, dN_dy) + 1.0;

        detF =  F[0]*F[3] - F[1]*F[2];

        //if(detF < 0.0)   return 1;

        //if(finite)
        //{
          //GeomData->computeBasisFunctions2D(1, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          //dvol = Jac * fact;
        //}

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
          F33 = 1.0;
        
        //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", F[0], F[1], F[2], F[3]);

        matlib2d_(matDat, F, &F33, stress, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        strain[0] = F[0] - 1.0;
        strain[1] = F[3] - 1.0;
        strain[2] = F33  - 1.0;
        strain[3] = 0.5 * (F[1] + F[2]);

        //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", strain[0], strain[1], strain[2], strain[3]);
        //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", stress[0], stress[1], stress[2], stress[3]);
        //printf("  dvol =  %14.12f \n\n ", dvol);

        velo[0] = computeValueDot(0, N);
        velo[1] = computeValueDot(1, N);

        // kinetic energy
        energy[0] += (0.5*dvol0*rho*(velo[0]*velo[0]+velo[1]*velo[1]));
        energy[1] += (0.5*dvol*(strain[0]*stress[0]+strain[1]*stress[1]+strain[2]*stress[2]+strain[3]*stress[3]));
    }//gp

    //printf("\n\n");
    return;
}




