
#include "Debug.h"
#include "Plot.h"
#include "FunctionsElement.h"
#include "NurbsElem1DEulerBeam.h"
#include "NurbsShapeFunctions.h"
using namespace std;

extern Plot     plot;
//extern MpapTime mpapTime;






NurbsElem1DEulerBeam::NurbsElem1DEulerBeam(void)
{
  if (debug) cout << " constructor NurbsElem1DEulerBeam\n\n";

  return;
}



NurbsElem1DEulerBeam::~NurbsElem1DEulerBeam()
{
  if (debug) cout << " destructor NurbsElem1DEulerBeam\n\n";

  return;
}



void   NurbsElem1DEulerBeam::prepareElemData()
{
   // set the data in the base class "NurbsElement"

   startindex.setDim(1);

   NurbsElement::prepareElemData();

   elmDat = &(curve0->ElemProp.data[0]);

   nGP = (int) elmDat[0];
   thick = elmDat[1];
     EI  = elmDat[2] * elmDat[3];
      cf = elmDat[4];
    pres = elmDat[5];

   return;
}




void   NurbsElem1DEulerBeam::initialiseDOFvalues()
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




void   NurbsElem1DEulerBeam::initialiseKnotsAtGPs()
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



int NurbsElem1DEulerBeam::calcStiffnessAndResidual()
{
    int   ii, jj,  gp;

    double   Jmod, Jac;
    vector<double>  dN_dx(nlbf), d2N_dx2(nlbf), N(nlbf);

    for(ii=0;ii<nsize;ii++)
      stiffness_local[ii].zero();

    resi.zero();

    for(gp=0;gp<nGP;gp++)
    {
        curve1->ShapeFunsAndDerivatives2(startindex[0], knotsAtGPs[gp], &N[0], &dN_dx[0], &d2N_dx2[0], Jac);
        Jmod = Jac*gaussweights[gp];

       for(ii=0;ii<nlbf;ii++)
       {
          for(jj=0;jj<nlbf;jj++)
          {
              stiffness_local[ii][jj] += (EI * d2N_dx2[ii] * d2N_dx2[jj] )*Jmod;
          }
       }
    }

  return 0;
}




int NurbsElem1DEulerBeam::calcMassMatrix(int lumpInd, double dt)
{
  return 0;
}





int NurbsElem1DEulerBeam::calcLoadVector()
{
    // initialize resi vector
    resi.zero();


  return 0;
}



void   NurbsElem1DEulerBeam::AssembleElementMatrix2(int index, MatrixXd& stiff, MatrixXd& MM)
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

