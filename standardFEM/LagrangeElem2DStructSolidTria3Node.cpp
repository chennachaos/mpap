
#include "LagrangeElem2DStructSolidTria3Node.h"

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



LagrangeElem2DStructSolidTria3Node::LagrangeElem2DStructSolidTria3Node()
{
  if (debug) cout << " constructor LagrangeElem2DStructSolidTria3Node\n\n";

  degree = 1;
  npElem = 3;
  nlbf   = 3;
  ndof   = 2;
  nsize  = nlbf*ndof;

}

LagrangeElem2DStructSolidTria3Node::~LagrangeElem2DStructSolidTria3Node()
{
  if (debug) cout << " destructor LagrangeElem2DStructSolidTria3Node\n\n";
}


void LagrangeElem2DStructSolidTria3Node::prepareElemData()
{
  LagrangeElement::prepareElemData();

  int ii, jj;

  degree = 1;
  npElem = 3;
  nlbf   = 3;
  ndof   = 2;
  nsize  = nlbf*ndof;


   // set the element property variables

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    matDat = &(SolnData->MatlProp[matType].data[0]);

   int nGP1  = (int) elmDat[0] ;
   int nGP2  = nGP1;
   nGP   = nGP1 * nGP2;
   finiteInt  = (int) elmDat[2] ;
   sss        = (int) elmDat[3] ;
   thick      = elmDat[4] ;
   //rho        = elmDat[5] ;
   //bforce[0]  = elmDat[6] ;
   //bforce[1]  = elmDat[7] ;
   matId      = SolnData->MatlProp[matType].id + 1;
   finite     = (finiteInt >= 1) ;
   axsy       = (sss == 3);
   //followerLoadFlag = (elmDat[8] == 1);
   followerLoadFlag = false;

  if(sss != 1) thick = 1.0; // for plane strain and axisymmetric problems

  return;
}



void LagrangeElem2DStructSolidTria3Node::prepareElemData2()
{
  cout << " LagrangeElem2DStructSolidTria3Node::prepareElemData2() " << endl;
  return;
}



int LagrangeElem2DStructSolidTria3Node::calcLoadVector()
{
  return 0;
}


int LagrangeElem2DStructSolidTria3Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
//   char fct[] = "LagrangeElem2DStructSolidTria3Node::calcStiffnessAndResidual";
//   computerTime.go(fct);

  double  F[4], detF=0.0, F33, fact, fact1, fact2, dvol, dvol0;
  double  Jac, dt, totvol=0.0, bb1, bb2, bb3, JacMult;
  double  af, am, d1, aa;

  VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);
  VectorXd  dispC(nsize), velC(nsize), accC(nsize);
  MatrixXd  Mlocal(nsize, nsize);

  double  stre[4], cc[4][4], bc[2][4], param[2], bforce[2];

  int   err,  isw,  count,  count1, index, ll = 0, ii, jj, gp, TI, TIp1, TJ, TJp1;
  int   ind1, ind2, kk;
  
  double  BULK = matDat[0] ;
  double  mu   = matDat[1] ;
  double  lamb = BULK - 2.0*mu/3.0;

  double rho = elmDat[5] ;
  bforce[0]  = 0.0;
  bforce[1]  = 0.0;

  //bforce[0]  = elmDat[6]*timeFunction[0].prop ;
  //bforce[1]  = elmDat[7]*timeFunction[0].prop ;

  af = SolnData->td(2);
  d1 = SolnData->td(5);
  aa = SolnData->td(10);
  
  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", BULK, mu, bforce[0], bforce[1]);
  //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", rho, af, d1, aa);
  //printf(" %5d \t %5d \t %5d \t %5d \n ", nodeNums[0], nodeNums[1], nodeNums[2], nodeNums[3]);

  double xNode[3], yNode[3], xx, yy;

  for(ii=0;ii<npElem;ii++)
  {
    //xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
    //yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];

    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
  }

    cout << xNode[0] << endl;
    cout << yNode[0] << endl;
    cout << xNode[1] << endl;
    cout << yNode[1] << endl;
    cout << xNode[2] << endl;
    cout << yNode[2] << endl;
    

  double totVolume = 0.5*(xNode[0]*(yNode[1]-yNode[2]) + xNode[1]*(yNode[2]-yNode[0]) + xNode[2]*(yNode[0]-yNode[1]));

  cout << " volume = " << totVolume << endl;

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

  count = 1;   ll = 0;   err = 0;   isw = 3;
  dt = mpapTime.dt;

  //cout << " finite = " << finite << endl;

  nGP = (int) elmDat[0] ;

  vector<double>  gausspoints1, gausspoints2, gaussweights;

  getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);


    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];

          GeomData->computeBasisFunctions2D(0, 1, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          fact = gaussweights[gp];

          dvol0 = Jac * fact;
          dvol = dvol0;

          xx = yy= 0.0;
          for(ii=0;ii<nlbf;ii++)
          {
            xx += N[ii]*xNode[ii];
            yy += N[ii]*yNode[ii];
          }

        //for(ii=0; ii<nlbf; ii++)
          //cout << N[ii] << '\t' << dN_dx[ii] << '\t' << dN_dy[ii] << endl;


        //if(axsy)
          //dvol *= 2.0*PI*yy;

        //GeomData->computeDeformationGradient(1, nodeNums, &dN_dx(0), &dN_dy(0), F, detF);

        F[0] = computeValueCur(0, dN_dx) + 1.0;
        F[2] = computeValueCur(0, dN_dy);
        F[1] = computeValueCur(1, dN_dx);
        F[3] = computeValueCur(1, dN_dy) + 1.0;

        detF =  F[0]*F[3] - F[1]*F[2];

        //if(detF < 0.0)   return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(1, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          dvol = Jac * fact;
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

        if(err !=0)    return 1;

        if(finite)
          dvol *= F33;

        for(ii=0;ii<4;ii++)
        {
          stre[ii] *= dvol;
          for(jj=0;jj<4;jj++)
            cc[ii][jj] *= dvol;
        }

        //==============================================
        // CALCULATE TANGENT STIFFNESS
        //==============================================

        //   part 1. -- material part (not necessarily symmetric!!)

        for(ii=0;ii<nlbf;ii++)
        {
          bb1 = dN_dx[ii];
          bb2 = dN_dy[ii];

          bb3 = N[ii]*dvol0*rho;

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
              bb1 = af*dN_dx[jj];
              bb2 = af*dN_dy[jj];

              TJ   = 2*jj;
              TJp1 = TJ+1;

              Klocal(TI,TJ)     +=  (bc[0][0] * bb1 + bc[0][2] * bb2) ;
              Klocal(TI,TJp1)   +=  (bc[0][1] * bb2 + bc[0][2] * bb1) ;
              Klocal(TIp1,TJ)   +=  (bc[1][0] * bb1 + bc[1][2] * bb2) ;
              Klocal(TIp1,TJp1) +=  (bc[1][1] * bb2 + bc[1][2] * bb1) ;
              
              Mlocal(TI,  TJ)   +=  bb3*N[jj];
              Mlocal(TIp1,TJp1) +=  bb3*N[jj];
          }
        }

        //   part 2. -- geometrical matrix  (if geometry nonlinear)

        if(finite)
        {
          for(ii=0;ii<nlbf;ii++)
          {
              fact1 = dN_dx[ii] * stre[0] + dN_dy[ii] * stre[3] ;
              fact2 = dN_dx[ii] * stre[3] + dN_dy[ii] * stre[1] ;

              TI   = 2*ii;
              TIp1 = TI+1;

              for(jj=0; jj<nlbf; jj++)
              {
                TJ   = 2*jj;
                TJp1 = TJ+1;

                fact = af*(fact1 * dN_dx[jj] + fact2 * dN_dy[jj] );

                Klocal(TI,TJ)     += fact ;
                Klocal(TIp1,TJp1) += fact ;
              }
          }
        }
  }//gp

  //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);

  Klocal +=  d1*Mlocal;
  Flocal -=  Mlocal*accC;

  //Klocal /= aa;
  
  //printVector(Flocal);

  return 0;
}






int LagrangeElem2DStructSolidTria3Node::calcInternalForces()
{
  return 0;
}



void LagrangeElem2DStructSolidTria3Node::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem2DStructSolidTria3Node::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem2DStructSolidTria3Node::projectStress(int varindex, double* outval)
{
  return;
}



void LagrangeElem2DStructSolidTria3Node::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem2DStructSolidTria3Node::projectIntVar(int index, double* outval)
{
  return;
}


int LagrangeElem2DStructSolidTria3Node::calcOutput(double u1, double v1)
{
  return 0;
}



void LagrangeElem2DStructSolidTria3Node::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}


