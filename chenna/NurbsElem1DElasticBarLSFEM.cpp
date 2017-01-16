
#include "Debug.h"
//#include "Plot.h"
//#include "FunctionsElement.h"
#include "NurbsElem1DElasticBarLSFEM.h"
#include "NurbsShapeFunctions.h"
#include "ComputerTime.h"
#include "TimeFunction.h"


using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
//extern Plot plot;
extern List<TimeFunction> timeFunction;






NurbsElem1DElasticBarLSFEM::NurbsElem1DElasticBarLSFEM(void)
{
  if (debug) cout << " constructor NurbsElem1DElasticBarLSFEM\n\n";

  return;
}



NurbsElem1DElasticBarLSFEM::~NurbsElem1DElasticBarLSFEM()
{
  if (debug) cout << " destructor NurbsElem1DElasticBarLSFEM\n\n";

  return;
}



void NurbsElem1DElasticBarLSFEM::prepareElemData()
{
   // set the data in the base class "NurbsElement"

   startindex.setDim(1);

   NurbsElement::prepareElemData();

   elmDat = &(curve0->ElemProp.data[0]);
   matDat = &(curve0->MatlProp.data[0]);

   nGP = (int) elmDat[0];

   int ff = (int) elmDat[1];

   finite = (ff == 1);

   //cout << " finite " << finite << endl;

   a = elmDat[2];
   c = elmDat[3];
   rho0 = elmDat[4];
   bforce = elmDat[5];

   return;
}




void NurbsElem1DElasticBarLSFEM::initialiseDOFvalues()
{
    int ind1, ind2, ind4, ii, jj;

    int *tt;
    tt = &(curve0->IEN[elenum][0]);

    ind1 = 0;
    for(ii=0;ii<nlbf;ii++)
    {
       ind2 = ndof * tt[ii];
       for(jj=0;jj<ndof;jj++)
       {
          ind4 = ind2+jj;

          forassy[ind1] = ind4;
          primvar[ind1] = curve0->Uinit[ind4];

          ind1++;
       }
    }

  return;
}




void NurbsElem1DElasticBarLSFEM::initialiseKnotsAtGPs()
{
   startindex[0] = curve0->INC[0][curve0->IEN[elenum][0]];

   knotsAtGPs.setDim(nGP);

   int ind1 = startindex[0]+curve0->p, ind2 = ind1+1;

   double *gausspoints = &(curve0->gausspoints[0]);

   double  val1 = 0.5*(curve0->U[ind2] - curve0->U[ind1]);
   double  val2 = 0.5*(curve0->U[ind2] + curve0->U[ind1]);

   // loop over Gauss points
   for(int gp=0;gp<nGP;gp++)
       knotsAtGPs[gp] = val1 * gausspoints[gp] + val2;

//   cout << '\t' << knotsAtGPs << endl;

    return;
}


void NurbsElem1DElasticBarLSFEM::createTractionDataVariable()
{
    assert(tracflag == true);
    
    if( !(tracdata.n > 0) )
    {
       tracdata.setDim(2);
       for(int ii=0;ii<2;ii++)
       {
         tracdata[ii].setDim(2);
         tracdata[ii][0] = tracdata[ii][1] = 7777.0;
       }
    }

    cout << " NurbsElem1DElasticBarLSFEM::createTractionDataVariable() .... DONE " << endl;

   return;
}

int NurbsElem1DElasticBarLSFEM::calcStiffnessAndResidual()
{
  if(finite)
    calcStiffnessAndResidual2();
  else
    calcStiffnessAndResidual1();

   return 0;
}


int NurbsElem1DElasticBarLSFEM::calcLoadVector()
{
  if(finite)
    calcLoadVector2();
  else
    calcLoadVector1();

   return 0;
}



int NurbsElem1DElasticBarLSFEM::calcStiffnessAndResidual1()
{
    int ii, jj, gp;

    double  fact, dvol0, Jac, dt, ff1, ff2, lambda, BULK, mu, b1, b2, b3, Jmod, D2u, f;
    vector<double>  N(nlbf), dN_dx(nlbf), d2N_dx2(nlbf);

    VectorXd  D(nsize);
   
    double *gaussweights = &(curve0->gaussweights[0]);
    double  *values1 = &(curve1->Values[0][0]);
   
    int *tt = &(curve0->IEN[elenum][0]);

    for(ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();
   
    resi.zero();

    //a = 240.565;

    for(gp=0;gp<nGP;gp++)   // loop over Gauss points
    {
        curve1->ShapeFunsAndDerivatives2(startindex[0], knotsAtGPs[gp], &N[0], &dN_dx[0], &d2N_dx2[0], Jac);

        //for(ii=0;ii<nlbf;ii++)
          //printf(" shape functions ... %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", N[ii], dN_dx[ii], d2N_dx2[ii], Jac);
        //printf("\n\n");

        Jmod = Jac * gaussweights[gp];

        D2u = 0.0;
        for(ii=0;ii<nlbf;ii++)
        {
           D2u += values1[tt[ii]]*d2N_dx2[ii];

           D(ii)  = a*d2N_dx2[ii];
        }

        //printf(" %14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f \t%14.8f\n\n", ff, Du(0), Du(1), dp(0), dp(1), f(0), f(1), pres);

        f = -a*D2u - bforce * rho0 * timeFunction[0].prop;

        //printf(" f ... %12.6f  \t %12.6f  \t %12.6f \n", D2u, f, Jmod);

        for(ii=0;ii<nsize;ii++)
        {
           resi[ii] += Jmod*(D(ii)*f);

           for(jj=0;jj<nsize;jj++)
             stiffness_local[ii][jj]  +=  Jmod*(D(ii)*D(jj));
        }
    }

//printStiffnessMatrix();
//printf("\n\n");
//printForceVector();

   return 0;
}





int NurbsElem1DElasticBarLSFEM::calcLoadVector1()
{
   if(tracflag)
   {
      for(int ii=0;ii<nsize;ii++)
        stiffness_local[ii].zero();

      resi.zero();

      double  *values1 = &(curve1->Values[0][0]);
   
      int *tt = &(curve0->IEN[elenum][0]);

      int p = curve0->p, ii, jj, gp, index;
      vector<double>  N(nlbf), dN_dx(nlbf);
      double  Jac, fact, ALPHA, ALPHA1, res, J, Jmod, b1, fact1;

      ALPHA1 = ALPHA = 1.0;

      // side #1
      if(!CompareDoubles(tracdata[0][0], 7777))
      {
          curve1->ShapeFunsAndDerivatives(startindex[0], 0.0, &N[0], &dN_dx[0], J);

          Jmod = J;

          res = tracdata[0][0] ;
          for(ii=0;ii<nlbf;ii++)
          {
              b1 = values1[tt[ii]];

              res -= b1*N[ii];
          }

          Jmod *= ALPHA;
          for(ii=0;ii<nlbf;ii++)
          {
              fact1 = N[ii]*Jmod;

              resi[ii] += fact1*res;

              for(jj=0;jj<nlbf;jj++)
                stiffness_local[ii][jj]  +=  fact1 * N[jj];
          }           
          //cout << " side1 done " << endl;
      }
      // side #2
      if(!CompareDoubles(tracdata[1][0],7777))
      {
          curve1->ShapeFunsAndDerivatives(startindex[0], 1.0, &N[0], &dN_dx[0], J);

          Jmod = J;

          res = tracdata[1][0] ;
          for(ii=0;ii<nlbf;ii++)
          {
              b1 = values1[tt[ii]];

              res -= b1*N[ii];
          }

          Jmod *= ALPHA;
          for(ii=0;ii<nlbf;ii++)
          {
              fact1 = N[ii]*Jmod;

              resi[ii] += fact1*res;

              for(jj=0;jj<nlbf;jj++)
                stiffness_local[ii][jj]  +=  fact1 * N[jj];
          }            
          //cout << " side2 done " << endl;
        }
    }
//printStiffnessMatrix();
//printf("\n\n");
//printForceVector();

  return 0;
}



/*
int NurbsElem1DElasticBarLSFEM::calcStiffnessAndResidual2()
{
    // in the current configuration

    int ii, jj, gp;

    double  fact, dvol, dvol0, Jac, dt, lambda, BULK, mu, J, fact1, fact2, b1, rho;
    double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], d2N_dx2[nlbf], d2N_dy2[nlbf], d2N_dxy[nlbf];
    double  D2u, DF, Db, f, DJ, F, G, b;

    VectorXd  D(nsize), ff(3);

    BULK = matDat[0];
    mu = matDat[1];
    lambda = BULK - 2.0*mu/3.0;

    double *gaussweights = &(curve0->gaussweights[0]);
    double  *values1 = &(curve1->Values[0][0]);
   
    int *tt = &(curve0->IEN[elenum][0]);

    for(ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();
   
    resi.zero();

    for(gp=0;gp<nGP;gp++)   // loop over Gauss points
    {
        curve0->ShapeFunsAndDerivatives2(startindex[0], knotsAtGPs[gp], N, dN_dx, d2N_dx2, Jac);

        //for(ii=0;ii<nlbf;ii++)
          //printf(" shape functions ... %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", N[ii], dN_dx[ii], d2N_dx2[ii], Jac);
        //printf("\n\n");

        D2u = 0.0;
        F = 1.0;
        for(ii=0;ii<nlbf;ii++)
        {
           b1 = values1[tt[ii]];
           
           F += b1*dN_dx[ii];

           D2u += b1*d2N_dx2[ii];
        }

        b = F*F;
        G = 1.0/F;

        DF = D2u * G;
        Db = 2.0 * F * DF;
        J  = F;
        DJ = DF;

        b1 = lambda*(1.0-log(J));

        ff(0) = mu*Db - (lambda/J)*DJ;
        ff(1) = mu*(1.0+b)+ b1;
        ff(2) = (mu*(1.0+b)+ b1)*DJ;

        fact1 = 1.0/J;
        fact2 = -fact1/J;

        ff(0) *= fact1;
        ff(1) *= fact1;
        ff(2) *= fact2;

        curve1->ShapeFunsAndDerivatives2(startindex[0], knotsAtGPs[gp], N, dN_dx, d2N_dx2, Jac);

        dvol = Jac * gaussweights[gp];

        rho = rho0*timeFunction[0].prop/J;

        //fact1 = bforce*rho;
        fact1 = 0.0;

        for(ii=0;ii<nlbf;ii++)
           D(ii) = (ff(0)+ ff(2)-fact1)*dN_dx[ii]+ff(1)*d2N_dx2[ii];

        f = -bforce * rho;

        fact1 = mu/J;
        fact2 = 1.0/J/J;

        f -= (fact1*Db + fact2*(lambda*(1.0-log(J))-mu*(b-1.0))*DJ);

        //printf("\t f vector  = %14.8f \t%14.8f \n", dp(0), dp(1));
        //printf("\t f vector  = %14.8f \t%14.8f \t%14.8f \n", f(0), f(1), f(2));

        //dvol *= J;

        for(ii=0;ii<nsize;ii++)
        {
           resi[ii] += dvol*(D(ii)*f);

           for(jj=0;jj<nsize;jj++)
             stiffness_local[ii][jj]  +=  dvol*(D(ii)*D(jj)) ;
        }
    }//gp

//printStiffnessMatrix();
//printf("\n\n");
//printForceVector();

   return 0;
}
*/


//
int NurbsElem1DElasticBarLSFEM::calcStiffnessAndResidual2()
{
    // in reference configuration

    int ii, jj, gp;

    double  fact, dvol0, Jac, dt, lambda, BULK, mu, J, b1, rho;
    vector<double>  N(nlbf), dN_dx(nlbf), d2N_dx2(nlbf);
    double  D2u, DF, f, DJ, F;

    VectorXd  D(nsize), ff(2);

    BULK = matDat[0];
    mu = matDat[1];
    lambda = BULK - 2.0*mu/3.0;

    double *gaussweights = &(curve0->gaussweights[0]);
    double  *values1 = &(curve1->Values[0][0]);
   
    int *tt = &(curve0->IEN[elenum][0]);

    for(ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();
   
    resi.zero();

    for(gp=0;gp<nGP;gp++)   // loop over Gauss points
    {
        curve0->ShapeFunsAndDerivatives2(startindex[0], knotsAtGPs[gp], &N[0], &dN_dx[0], &d2N_dx2[0], Jac);

        dvol0 = Jac * gaussweights[gp];

        //for(ii=0;ii<nlbf;ii++)
          //printf(" shape functions ... %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n", N[ii], dN_dx[ii], d2N_dx2[ii], Jac);
        //printf("\n\n");

        D2u = 0.0;
        F = 1.0;
        for(ii=0;ii<nlbf;ii++)
        {
           b1 = values1[tt[ii]];
           
           F += b1*dN_dx[ii];

           D2u += b1*d2N_dx2[ii];
        }

        DF = D2u;
        J  = F;
        DJ = DF;

        ff(0) = (-3.0*lambda+2.0*lambda*log(J)-2.0*mu)*DJ/J/J/J;
        ff(1) = mu + (lambda*(1.0-log(J))+mu)/J/J;

        //curve1->ShapeFunsAndDerivatives2(startindex[0], knotsAtGPs[gp], N, dN_dx, d2N_dx2, Jac);

        for(ii=0;ii<nlbf;ii++)
           D(ii) = ff(0)*dN_dx[ii]+ff(1)*d2N_dx2[ii];

        rho = rho0*timeFunction[0].prop;

        f = -bforce * rho;

        f -= (mu*D2u + (lambda*(1.0-log(J))+mu)*DJ/J/J);

        //printf("\t f vector  = %14.8f \t%14.8f \t%14.8f \n", f(0), f(1), f(2));

        for(ii=0;ii<nsize;ii++)
        {
           resi[ii] += dvol0*(D(ii)*f);

           for(jj=0;jj<nsize;jj++)
             stiffness_local[ii][jj]  +=  dvol0*(D(ii)*D(jj)) ;
        }
    }//gp

   return 0;
}
//



//
int NurbsElem1DElasticBarLSFEM::calcLoadVector2()
{
   if(tracflag)
   {
      for(int ii=0;ii<nsize;ii++)
        stiffness_local[ii].zero();

      resi.zero();

      double  *values1 = &(curve1->Values[0][0]);
   
      int *tt = &(curve0->IEN[elenum][0]);

      int p = curve0->p, ii, jj, gp, index;
      vector<double>  N(nlbf), dN_dx(nlbf);
      double Jac, fact, ALPHA, ALPHA1, res, J, Jmod, b1, fact1;

      ALPHA1 = ALPHA = 1.0;

      // side #1
      if(!CompareDoubles(tracdata[0][0], 7777))
      {
          curve1->ShapeFunsAndDerivatives(startindex[0], 0.0, &N[0], &dN_dx[0], J);

          Jmod = J;

          res = tracdata[0][0] ;
          for(ii=0;ii<nlbf;ii++)
          {
              b1 = values1[tt[ii]];

              res -= b1*N[ii];
          }

          Jmod *= ALPHA;
          for(ii=0;ii<nlbf;ii++)
          {
              fact1 = N[ii]*Jmod;

              resi[ii] += fact1*res;

              for(jj=0;jj<nlbf;jj++)
                stiffness_local[ii][jj]  +=  fact1 * N[jj];
          }           
          //cout << " side1 done " << endl;
      }
      // side #2
      if(!CompareDoubles(tracdata[1][0],7777))
      {
          curve1->ShapeFunsAndDerivatives(startindex[0], 1.0, &N[0], &dN_dx[0], J);

          Jmod = J;

          res = tracdata[1][0] ;
          for(ii=0;ii<nlbf;ii++)
          {
              b1 = values1[tt[ii]];

              res -= b1*N[ii];
          }

          Jmod *= ALPHA;
          for(ii=0;ii<nlbf;ii++)
          {
              fact1 = N[ii]*Jmod;

              resi[ii] += fact1*res;

              for(jj=0;jj<nlbf;jj++)
                stiffness_local[ii][jj]  +=  fact1 * N[jj];
          }            
          //cout << " side2 done " << endl;
        }
    }
//printStiffnessMatrix();
//printf("\n\n");
//printForceVector();
  return 0;
}





int NurbsElem1DElasticBarLSFEM::calcMassMatrix(int lumpInd, double dt)
{
   /*
   *  lumpInd = 1 --> consistent Mass matrix

   *          = 2 --> Row-sum Mass lumping
   *          = 3 --> proportional Mass lumping
   */

/*
    int ii, jj, gp;

    double  Jmod, fact, Jac, N[nlbf], dN_dx[nlbf];

    mass_local.setDim(nsize);
    for(ii=0;ii<nsize;ii++)
    {
       mass_local[ii].setDim(nsize);
       mass_local[ii].zero();
    }

    double *gaussweights = &(curve0->gaussweights[0]);

    for(gp=0;gp<nGP;gp++)   // loop over Gauss points
    {
        curve1->ShapeFunsAndDerivatives(startindex[0], knotsAtGPs[gp], N, dN_dx, Jac);

        Jmod = Jac * gaussweights[gp];

        for(ii=0;ii<nlbf;ii++)
        {
            fact = N[ii] * Jmod;

            for(jj=0;jj<nlbf;jj++)
               mass_local[ii][jj]  +=  fact * N[jj];
        }
    }
*/
/*
    for(ii=0;ii<nsize;ii++)
    {
        for(jj=0;jj<nsize;jj++)
           cout << '\t' << mass_local[ii][jj];
        cout << endl;
    }
*/
  return 0;
}




void NurbsElem1DElasticBarLSFEM::AssembleElementMatrix2(int index, MatrixXd& stiff, MatrixXd& MM)
{
    int  row, aa, bb;

    for(aa=0;aa<nsize;aa++)
    {
        row = forassy[aa];
        for(bb=0;bb<nsize;bb++)
        {
            stiff(row, forassy[bb]) += stiffness_local[aa][bb];
            MM(row, forassy[bb])  += mass_local[aa][bb];
        }
    }

  return;
}





int NurbsElem1DElasticBarLSFEM::calcOutput(double u1, double v1)
{
  return 0;
}


void NurbsElem1DElasticBarLSFEM::discreteContourplot(int, int, int, int, double, double)
{
  return;
}

void NurbsElem1DElasticBarLSFEM::projectToKnots(bool, int, int, int)
{
  return;
}

void NurbsElem1DElasticBarLSFEM::projectStrain(int, int, double*)
{
  return;
}

void NurbsElem1DElasticBarLSFEM::projectStress(int, double*)
{
  return;
}

void NurbsElem1DElasticBarLSFEM::projectIntVar(int, double*)
{
  return;
}

void NurbsElem1DElasticBarLSFEM::toPostprocess(int, int, int,  SparseMatrixXd&, VectorXd&)
{
  return;
}

int NurbsElem1DElasticBarLSFEM::calcStiffnessMatrix(double tt)
{
  return 0;
}

int NurbsElem1DElasticBarLSFEM::calcInternalForces()
{
  return 0;
}

int NurbsElem1DElasticBarLSFEM::applyDirichletBCs()
{
  return 0;
}
