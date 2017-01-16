
#include "Debug.h"
//#include "Plot.h"
//#include "FunctionsElement.h"
#include "NurbsElem1DElasticBar.h"
#include "NurbsShapeFunctions.h"
#include "ComputerTime.h"
#include "TimeFunction.h"

using namespace std;

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
//extern Plot plot;
extern List<TimeFunction> timeFunction;


NurbsElem1DElasticBar::NurbsElem1DElasticBar()
{
  if (debug) cout << " constructor NurbsElem1DElasticBar\n\n";

  return;
}



NurbsElem1DElasticBar::~NurbsElem1DElasticBar()
{
  if (debug) cout << " destructor NurbsElem1DElasticBar\n\n";

  return;
}



void NurbsElem1DElasticBar::prepareElemData()
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




void NurbsElem1DElasticBar::initialiseDOFvalues()
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




void NurbsElem1DElasticBar::initialiseKnotsAtGPs()
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


/*
int NurbsElem1DElasticBar::calcStiffnessAndResidual()
{
    int ii, jj, gp;

    //set the dimentions of stiffness matrix and force vector
    for(ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();

    resi.zero();

    double *gaussweights = &(curve0->gaussweights[0]);

    double  Jmod, fact, Jac, N[nlbf], dN_dx[nlbf];

    for(gp=0;gp<nGP;gp++)   // loop over Gauss points
    {
       curve0->ShapeFunsAndDerivatives(startindex[0], knotsAtGPs[gp], N, dN_dx, Jac);

       Jmod = Jac * gaussweights[gp];

       for(ii=0;ii<nsize;ii++)
       {
          for(jj=0;jj<nsize;jj++)
             stiffness_local[ii][jj] += (a*dN_dx[ii]*dN_dx[jj] + c*N[ii]*N[jj] )*Jmod;
       }
    }

//    printStiffnessMatrix();

   return 0;
}
*/




int NurbsElem1DElasticBar::calcStiffnessAndResidual()
{
    int ii, jj, gp;

    //set the dimentions of stiffness matrix and force vector
    for(ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();

    resi.zero();

    double *gaussweights = &(curve0->gaussweights[0]);

    double  dvol0, dvol, fact, Jac, F, J, lambda, BULK, mu, b, rho, D, sig;
    vector<double>  N(nlbf), dN_dx(nlbf);

    BULK = matDat[0];
    mu = matDat[1];
    lambda = BULK - 2.0*mu/3.0;

    for(gp=0;gp<nGP;gp++)   // loop over Gauss points
    {
       curve0->ShapeFunsAndDerivatives(startindex[0], knotsAtGPs[gp], &N[0], &dN_dx[0], Jac);

       dvol0 = Jac * gaussweights[gp];
       dvol  = dvol0;

       curve1->deformationGradient(startindex[0], 1, &dN_dx[0], F);

       //for(ii=0;ii<nlbf;ii++)
          //printf("\t F ...  = %14.8f \t%14.8f \t%14.8f \n", F, N[ii], dN_dx[ii]);

       J = F;
       b = F*F;

       D = lambda/J - 2.0*(lambda*log(J)-mu)/J;
       sig = ((lambda*log(J)-mu)+mu*b)/J;

       if(finite)
       {
          curve1->ShapeFunsAndDerivatives(startindex[0], knotsAtGPs[gp], &N[0], &dN_dx[0], Jac);
          dvol = Jac * gaussweights[gp];
       }

       //printf("\t F ...  = %14.8f \t%14.8f \t%14.8f \n", F, Jac, dvol);

       D *= dvol;
       sig *= dvol;

       rho = bforce * rho0 * timeFunction[0].prop*dvol0;

       for(ii=0;ii<nlbf;ii++)
       {
          resi[ii] += (-dN_dx[ii]*sig + rho*N[ii]);

          fact = D * dN_dx[ii];
          for(jj=0;jj<nlbf;jj++)
             stiffness_local[ii][jj] += (fact*dN_dx[jj]);
       }
       if(finite)
       {
          for(ii=0;ii<nlbf;ii++)
          {
             fact = dN_dx[ii]*sig;
             for(jj=0;jj<nlbf;jj++)
               stiffness_local[ii][jj] += (fact*dN_dx[jj]);
          }
       }
    }

//    printStiffnessMatrix();

   return 0;
}





int NurbsElem1DElasticBar::calcMassMatrix(int lumpInd, double dt)
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





int NurbsElem1DElasticBar::calcLoadVector()
{
    // initialize resi vector
    resi.zero();


  return 0;
}




void NurbsElem1DElasticBar::AssembleElementMatrix2(int index, MatrixXd& stiff, MatrixXd& MM)
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



