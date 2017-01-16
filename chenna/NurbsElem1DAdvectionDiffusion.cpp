
#include "Debug.h"
#include "Plot.h"
#include "FunctionsElement.h"
#include "NurbsElem1DAdvectionDiffusion.h"
#include "NurbsShapeFunctions.h"
using namespace std;

extern Plot     plot;
//extern MpapTime mpapTime;






NurbsElem1DAdvectionDiffusion::NurbsElem1DAdvectionDiffusion(void)
{
  if (debug) cout << " constructor NurbsElem1DAdvectionDiffusion\n\n";

  startindex.setDim(1);

  return;
}



NurbsElem1DAdvectionDiffusion::~NurbsElem1DAdvectionDiffusion()
{
  if (debug) cout << " destructor NurbsElem1DAdvectionDiffusion\n\n";

  return;
}


void NurbsElem1DAdvectionDiffusion::prepareElemData()
{
   // set the data in the base class "NurbsElement"
   NurbsElement::prepareElemData();

   elmDat = &(curve0->ElemProp.data[0]);

   nGP  =  (int) elmDat[0];
     a  = elmDat[1];
     mu = elmDat[2];
     s  = elmDat[3];

    tau = 1.0/(2*curve0->nelem*a);
    
    cout << uvalues[1] << '\t' << uvalues[0] << endl;
    
    //tau = (uvalues[1] - uvalues[0])/(2.0*a);
    
    if( mu < 0 )
      tau = - tau;
   return;
}




void NurbsElem1DAdvectionDiffusion::initialiseDOFvalues()
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




void NurbsElem1DAdvectionDiffusion::initialiseKnotsAtGPs()
{
   startindex[0] = curve0->INC[0][curve0->IEN[elenum][0]];

   knotsAtGPs.setDim(nGP);

   int ind1 = startindex[0]+curve0->p, ind2 = ind1+1;

   double *gausspoints = &(curve0->gausspoints[0]);

   double  val1 = 0.5*(curve0->U[ind2] - curve0->U[ind1]);
   double  val2 = 0.5*(curve0->U[ind2] + curve0->U[ind1]);

   // loop over Gauss points
   for(int gp=0;gp<nGP;gp++)
   {
       knotsAtGPs[gp] = val1 * gausspoints[gp] + val2;
   }

//   cout << '\t' << knotsAtGPs << endl;

  return;
}


/*
int NurbsElem1DAdvectionDiffusion::calcStiffnessAndResidual()
{
    int ii, jj, gp;

    //zero stiffness matrix and residue vector
    for(ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();

    resi.zero();

    double *gausspoints  = &(curve0->gausspoints[0]);
    double *gaussweights = &(curve0->gaussweights[0]);

    double  Jac, Jmod, N[nlbf], dN_dx[nlbf], d2N_dx2[nlbf];

    double  temp[nsize]; // temporary vector to store terms due to LS stabilization

    for(gp=0;gp<nGP;gp++)   // loop over Gauss points
    {
       // Generate Shape Functions at the Gauss Point

       curve1->ShapeFunsAndDerivatives2(startindex[0], knotsAtGPs[gp], N, dN_dx, d2N_dx2, Jac);

       Jmod = Jac * gaussweights[gp];

       // Least Square Stabilization term
       for(ii=0;ii<nsize;ii++)
         temp[ii] = a * dN_dx[ii] - mu * d2N_dx2[ii];

       //stiffness matrix computations
       for(ii=0;ii<nsize;ii++)
       {
          for(jj=0;jj<nsize;jj++)
          {
              stiffness_local[ii][jj] += (mu * dN_dx[ii] * dN_dx[jj] + a * N[ii] * dN_dx[jj] + tau*temp[ii]*temp[jj])*Jmod;
          }
       }

       //force vector computations
       if(!CompareDoubles(s, 0.0))
       {
          for(ii=0;ii<nsize;ii++)
          {
             resi[ii] += s*N[ii]*Jmod;
          }
       }
    }

  return 0;
}
*/



int NurbsElem1DAdvectionDiffusion::calcStiffnessAndResidual()
{
    int ii, jj, gp;

    //zero stiffness matrix and residue vector
    for(ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();

    resi.zero();

    double *gausspoints  = &(curve0->gausspoints[0]);
    double *gaussweights = &(curve0->gaussweights[0]);

    double  Jac, Jmod, x0, xcoord, L, h, Pe, xi;

    vector<double>  N(nlbf), dN_dx(nlbf), d2N_dx2(nlbf), temp(nsize), temp2(nsize); // temporary vector to store terms due to LS stabilization
    
    x0 = -1.0;
    
    //cout << tau << endl;

    //tau *= 100.0 ;
    
    L = 2.0;
    h = L/curve0->nelem;

    for(gp=0;gp<nGP;gp++)   // loop over Gauss points
    {
       // Generate Shape Functions at the Gauss Point

       curve1->ShapeFunsAndDerivatives2(startindex[0], knotsAtGPs[gp], &N[0], &dN_dx[0], &d2N_dx2[0], Jac);

       Jmod = Jac * gaussweights[gp];

       xcoord = x0 + 2.0 * knotsAtGPs[gp];
       a = xcoord;
       
       // Least Square Stabilization term
       for(ii=0;ii<nsize;ii++)
       {
         temp[ii] = a * dN_dx[ii];
         temp2[ii] = a * dN_dx[ii] - mu * d2N_dx2[ii];
       }

       //tau = 2.0/(curve0->nelem*abs(a));
       //tau = -tau;

       Pe = a*h*L/mu/2.0;
       xi = 1.0/tanh(Pe) - 1.0/Pe;
       tau = h * xi*L/a/2.0;

       //stiffness matrix computations
       for(ii=0;ii<nsize;ii++)
       {
          for(jj=0;jj<nsize;jj++)
          {
              //stiffness_local[ii][jj] += (mu * dN_dx[ii] * dN_dx[jj] + a * N[ii] * dN_dx[jj] + tau*temp[ii]*temp[jj])*Jmod;

              stiffness_local[ii][jj] += (mu * dN_dx[ii] * dN_dx[jj] + a * N[ii] * dN_dx[jj] + tau*temp2[ii]*temp2[jj])*Jmod;
              
              //stiffness_local[ii][jj] += (mu * dN_dx[ii] * dN_dx[jj] + a * N[jj] * dN_dx[ii] + tau*temp2[ii]*temp[jj])*Jmod;
          }
       }

       s = (mu*PI*PI*cos(PI*xcoord) - PI*xcoord*sin(PI*xcoord) ) * Jmod;

       for(ii=0;ii<nsize;ii++)
         //resi[ii] += (s*N[ii] );
         resi[ii] += (s*(N[ii] + tau * temp2[ii]));
    }

  return 0;
}





int NurbsElem1DAdvectionDiffusion::calcLoadVector()
{
    // initialize resi vector
    resi.zero();


  return 0;
}
