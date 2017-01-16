
#include "Debug.h"
#include "FunctionsProgram.h"
#include "PropertyTypeEnum.h"
#include "Plot.h"
#include "NurbsShapeFunctions.h"
#include "NurbsElementSolid.h"
#include <assert.h>

#include "ExactSolutionsElasticity.h"

using namespace std;

extern Plot     plot;
extern MpapTime mpapTime;




NurbsElementSolid::NurbsElementSolid()
{
  if (debug) cout << " constructor NurbsElementSolid\n\n";

  // cout << "     NurbsElementSolid: constructor ...\n\n";
}





NurbsElementSolid::~NurbsElementSolid()
{
  if (debug)   cout << " destructor NurbsElementSolid\n\n";

//NurbsElement::~NurbsElement();

//    cout << "     NurbsElementSolid: destructor ...\n\n";
}






void NurbsElementSolid::initialiseDOFvalues()
{
    int ind1, ind2, ind4, ii, jj;

    int *tt;

    tt = &(surf0->IEN[elenum][0]);
    ind1 = 0;
    for(ii=0;ii<nlbf;ii++)
    {
       ind2 = ndof * tt[ii];
       for(jj=0;jj<ndof;jj++)
       {
          ind4 = ind2+jj;

          forassy[ind1] = ind4;
          primvar[ind1] = surf0->Uinit[ind4];

          ind1++;
       }
    }

  return;
}



void NurbsElementSolid::createTractionDataVariable()
{
    assert(tracflag == true);
    
    if( !(tracdata.n > 0) )
    {
       tracdata.setDim(4);
       for(int ii=0;ii<4;ii++)
       {
         tracdata[ii].setDim(2);
         tracdata[ii].zero();
         //tracdata[ii][0] = tracdata[ii][1] = 7777.0;
       }
    }

   return;
}



void NurbsElementSolid::prepareElemData()
{
   // set the data in the base class "NurbsElement"
   NurbsElement::prepareElemData();

   uvalues[2] = uvalues[1] - uvalues[0];
   vvalues[2] = vvalues[1] - vvalues[0];
   
   JacMultFact = 0.25 * uvalues[2] * vvalues[2];

   vals2project.setDim(4);

   startindex.setDim(2);

   // set the element property variables

   elmDat = &(surf0->ElemProp.data[0]);
   matDat = &(surf0->MatlProp.data[0]);

         nGP1    		= (int) elmDat[0] ;
         nGP2    		= (int) elmDat[1] ;
         nGP    		= nGP1 * nGP2;
         finiteInt		= (int) elmDat[2] ;
         sss		 	= (int) elmDat[3] ;
         thick  		= elmDat[4] ;
         bforce[0]		= elmDat[5] ;
         bforce[1]		= elmDat[6] ;
         rho0   		= elmDat[7] ;
         followerLoadFlag	= (elmDat[8] == 1);
         matId  		= surf0->MatlProp.id + 1;
         finite		= (finiteInt >= 1) ;
         axsy           = (sss == 3);

         if(sss != 1) thick = 1.0; // for plane strain and axisymmetric problems

    if(surf0->ElemProp.id == 7)
    {
       int ii, jj, sizep;

       sizep = surf2->nsize;

       forassyKup.setDim(nsize);
       for(ii=0;ii<nsize;ii++)
          forassyKup[ii].setDim(sizep);

       forassyKpu.setDim(sizep);

       for(ii=0;ii<sizep;ii++)
       {
          forassyKpu[ii].setDim(nsize);
       }

       for(ii=0;ii<nsize;ii++)
       {
          for(jj=0;jj<sizep;jj++)
          {
             forassyKup[ii][jj]  = -1;
             forassyKpu[jj][ii]  = -1;
          }
       }

       forassyKtt.setDim(sizep);
       for(ii=0;ii<sizep;ii++)
       {
          forassyKtt[ii].setDim(sizep);
          for(jj=0;jj<sizep;jj++)
            forassyKtt[ii][jj]  = -1;
       }

       resi2.setDim(sizep);
    }

    if(surf0->ElemProp.id == 8)
    {
       int ii, jj, sizep;

       sizep = surf2->nsize;

       resi2.setDim(2*sizep);

       forassyKut.setDim(nsize);
       forassyKup.setDim(nsize);

       for(ii=0;ii<nsize;ii++)
       {
          forassyKut[ii].setDim(sizep);
          forassyKup[ii].setDim(sizep);
       }

       forassyKtu.setDim(sizep);
       forassyKpu.setDim(sizep);

       for(ii=0;ii<sizep;ii++)
       {
          forassyKtu[ii].setDim(nsize);
          forassyKpu[ii].setDim(nsize);
       }

       for(ii=0;ii<nsize;ii++)
       {
          for(jj=0;jj<sizep;jj++)
          {
             forassyKut[ii][jj]  = -1;
             forassyKup[ii][jj]  = -1;
             forassyKtu[jj][ii]  = -1;
             forassyKpu[jj][ii]  = -1;
          }
       }

       forassyKtt.setDim(sizep);
       forassyKtp.setDim(sizep);
       forassyKpt.setDim(sizep);

       for(ii=0;ii<sizep;ii++)
       {
          forassyKtt[ii].setDim(sizep);
          forassyKtp[ii].setDim(sizep);
          forassyKpt[ii].setDim(sizep);

          for(jj=0;jj<sizep;jj++)
          {
             forassyKtt[ii][jj]  = -1;
             forassyKtp[ii][jj]  = -1;
             forassyKpt[ii][jj]  = -1;
          }
       }
    }

    //cout << " yyyyyyyyyyy " << endl;
    // Initialize element internal variables, if any are there
    NurbsElementSolid::initialiseIntVar();
    //cout << " yyyyyyyyyyy " << endl;

  return;
}





void NurbsElementSolid::setnivGP()
{
  // get number of internal variable per Gauss point

   double dmy[10];

  int  dmyI[3],
       isw   = 1,
  	mDim  = matdim_(&matId);

  if (mDim == 1) matlib1d_(matDat,dmy,dmy,dmy,dmy,dmy,dmy,dmy,
                           &matId,&nivGP,&finiteInt,&sss,&isw,dmyI);

  else if (mDim == 2) matlib2d_(matDat,dmy,dmy,dmy,dmy,dmy,dmy,dmy,
                                  &matId,&nivGP,&finiteInt,&sss,&isw,dmyI);

  else if (mDim == 3) matlib3d_(matDat,dmy,dmy,dmy,dmy,dmy,dmy,
                                    &matId,&nivGP,&finiteInt,&isw,dmyI);

  else prgError(1,"NurbsElementSolid::nivGP","invalid value of ndm!");

  return;
}




void NurbsElementSolid::initialiseIntVar()
{
  // set nivGP value
  setnivGP();

  if(nivGP > 0)
  {
      // allocate memory
      int n = nivGP * nGP;

      intVar1 = new double [n];
      intVar2 = new double [n];

      // initialise values

      double dmy[10];

      int   l, dmyI[3],
            ll   = 0,
            isw  = 2,
            mDim = matdim_(&matId);

     for(int gp2=0;gp2<nGP2;gp2++)
     {
        for(int gp1=0;gp1<nGP1;gp1++)
        {
           if (mDim == 1) matlib1d_(matDat,dmy,dmy,dmy,dmy,&(intVar1[ll]),&(intVar2[ll]),dmy,
   	                         &matId,&nivGP,&finiteInt,&sss,&isw,dmyI);

           else if (mDim == 2) matlib2d_(matDat,dmy,dmy,dmy,dmy,&(intVar1[ll]),&(intVar2[ll]),dmy,
                                &matId,&nivGP,&finiteInt,&sss,&isw,dmyI);

           else if (mDim == 3) matlib3d_(matDat,dmy,dmy,dmy,&(intVar1[ll]),&(intVar2[ll]),dmy,
	                         &matId,&nivGP,&finiteInt,&isw,dmyI);

           else prgError(1,"NurbsElementSolid::initialiseIntVar","invalid value of ndm!");

           ll += nivGP;
        }
     }

   //   for (int i=0; i<n; i++)  cout << '\t' << i << '\t' << intVar1[i] << '\t' << intVar2[i] << endl;
  }
  return;
}




void NurbsElementSolid::initialiseKnotsAtGPs()
{
   int  ind1, ind2, ind3, ind4, count=0, index, gp1, gp2;

   knotsAtGPs.setDim(nGP*2);

   double val1, val2, val3, val4, val5;

   ind1 = surf0->IEN[elenum][0];

   startindex[0] = surf0->INC[0][ind1];
   startindex[1] = surf0->INC[1][ind1];


   ind1 = startindex[0]+surf0->p;
   ind2 = ind1+1;
   ind3 = startindex[1]+surf0->q;
   ind4 = ind3+1;


   double *gausspoints1 = &(surf0->gausspoints1[0]);
   double *gausspoints2 = &(surf0->gausspoints2[0]);

   val1 = 0.5*uvalues[2];
   val2 = 0.5*(uvalues[1]+uvalues[0]);

   val3 = 0.5*vvalues[2];
   val4 = 0.5*(vvalues[1]+vvalues[0]);

   // loop over Gauss points
   for(gp2=0;gp2<nGP2;gp2++)
   {
      val5 = val3 * gausspoints2[gp2] + val4;
      for(gp1=0;gp1<nGP1;gp1++)
      {
         index = count*2;
         knotsAtGPs[index]   = val1 * gausspoints1[gp1] + val2;
         knotsAtGPs[index+1] = val5;

         count++;
      }
   }

  return;
}






void NurbsElementSolid::plotGaussPoints(int num, bool defFlg)
{
   int ind1 = startindex[0]+surf0->p, ind2 = ind1+1;
   int ind3 = startindex[1]+surf0->q, ind4 = ind3+1;

   double d = (plot.dAct[0] + plot.dAct[1]) * .0025, X[2] = {0.0, 0.0};
   int index;
   EPOINT EP;

   // loop over Gauss points
   for(int gp=0;gp<nGP;gp++)
   {
         index = gp*2;
         if(defFlg)
            EP = surf1->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();
         else
            EP = surf0->SurfacePoint(knotsAtGPs[index], knotsAtGPs[index+1]).CalcEuclid();

         X[0] = EP.x;
         X[1] = EP.y;

         plot.point(X,d);
   }

return;
}




//
int NurbsElementSolid::calcLoadVector()
{

//cout << " AAAAAAAAAA " << endl;

    // initialize Flocal vector
    Flocal.setZero();

   double *gausspoints1 = &(surf0->gausspoints1[0]);
   double *gausspoints2 = &(surf0->gausspoints2[0]);
   double *gaussweights1 = &(surf0->gaussweights1[0]);
   double *gaussweights2 = &(surf0->gaussweights2[0]);

     //   FORCES DUE TO SURFACE TRACTION
     //============================================
     // normal traction into the surface is taken as positive

     if(tracflag)
     {
        //cout << "       elem... : " << elenum << endl;
        //cout << tracdata[0][0] << '\t' << tracdata[0][1] << endl;
        //cout << tracdata[1][0] << '\t' << tracdata[1][1] << endl;
        //cout << tracdata[2][0] << '\t' << tracdata[2][1] << endl;
        //cout << tracdata[3][0] << '\t' << tracdata[3][1] << endl;

        double dvol0=0.0, fact1, J, Jmod, rad=0.0, dircos[2], tracX, tracY, param, theta;
        int index=0, ngbf1, ngbf2;
        double stress[2][2];

        int p = surf0->p, q = surf0->q, ind1, ind2, ii, gp;

        ListArray<CPOINT>  Pw1;
        ListArray<EPOINT>  SKL;

        EPOINT  EP1, EP2, Normal, EP;

        ngbf1 = surf0->ngbf1;
        ngbf2 = surf0->ngbf2;
        
        PlateWithHole  pwhole(sss, matDat[1], matDat[2]);
        
        double  len=0.0;

        // side #1
        if(!CompareDoubles(tracdata[0][0],0.0) || !CompareDoubles(tracdata[0][1],0.0))
        {
            double NN[p+1];

            Pw1.setDim(ngbf1);

            for(ii=0;ii<ngbf1;ii++)
               Pw1[ii] = surf0->Pw[ii][0];

           NurbsCURVE   curve_temp(Pw1, surf0->U, p);

           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              if( finite && followerLoadFlag )
                 NurbsShapeFunctions2DAlg2(surf1, startindex[0], startindex[1], gausspoints1[gp], -1.0, NN, J, dircos);
              else
                 NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], -1.0, NN, J, dircos);

              Jmod = J * gaussweights1[gp] * thick;
              
              len += Jmod;

              //   for axisymmetric problems compute radius
              if(axsy)
              {
                 if(finite && followerLoadFlag )
                    rad  = surf1->SurfacePoint(knotsAtGPs[2*gp], vvalues[0]).CalcEuclid().x;
                 else
                    rad  = surf0->SurfacePoint(knotsAtGPs[2*gp], vvalues[0]).CalcEuclid().x;

                 Jmod *= (twoPI * rad);
              }

              tracX = tracdata[0][0] * (-dircos[0]) + tracdata[0][1] * (-dircos[1]);
              tracY = tracdata[0][0] * (-dircos[1]) + tracdata[0][1] * (dircos[0]);

              ////////////////////////
              if(2<1)
	      {
              ind1 = startindex[0]+p; ind2 = ind1+1;

              param = ((surf1->U[ind2]-surf1->U[ind1])*gausspoints1[gp] + (surf1->U[ind2]+surf1->U[ind1])) / 2;

              curve_temp.CurveDerPointRat(param, 1, SKL);

              Normal = SKL[1]/SKL[1].Norm();

              //Normal.print2screen();

              Jmod =  SKL[1].Norm() * gaussweights1[gp] * thick * 0.5 * uvalues[2];

              //tracX = tracdata[0][0] * -Normal.y;
              //tracY = tracdata[0][0] * Normal.x;
              }

              ////////////////////

              //cout << tracX << '\t' << tracY <<  '\t' << J << '\t' << Jmod << endl;
              //cout << dircos[0] << '\t' << dircos[1] << endl;
              //cout << endl;

              for(ii=0;ii<=p;ii++)
              {
                 index = ndof*ii;
                 fact1 = NN[ii] * Jmod;
                 Flocal[index]   += tracX * fact1 ;
                 Flocal[index+1] += tracY * fact1 ;
              }
           }
        }
        //cout << " length = " << len << endl;

        // side #2
        if(!CompareDoubles(tracdata[1][0],0.0) || !CompareDoubles(tracdata[1][1],0.0))
        {
           double NN[q+1];

           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              if( finite && followerLoadFlag )
                 NurbsShapeFunctions2DAlg2(surf1, startindex[0], startindex[1], 1.0, gausspoints2[gp], NN, J, dircos);
              else
                 NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], 1.0, gausspoints2[gp], NN, J, dircos);

              Jmod = J * gaussweights2[gp] * thick;

              //printf("\t J, Jmod  = %12.6f \t %12.6f \n",J, Jmod);

              //   for axisymmetric problems compute radius
              if(axsy)
              {
                 if(finite && followerLoadFlag )
                    rad  = surf1->SurfacePoint(uvalues[1], knotsAtGPs[2*(nGP1*(gp+1)-1)+1]).CalcEuclid().x;
                 else
                    rad  = surf0->SurfacePoint(uvalues[1], knotsAtGPs[2*(nGP1*(gp+1)-1)+1]).CalcEuclid().x;

                 Jmod *= (twoPI * rad);
              }

              tracX = tracdata[1][0] * (-dircos[0]) + tracdata[1][1] * (-dircos[1]);
              tracY = tracdata[1][0] * (-dircos[1]) + tracdata[1][1] * (dircos[0]);

              for(ii=0;ii<=q;ii++)
              {
                 index = ndof * ((p + 1) * (ii+1) - 1);
                 fact1 = NN[ii] * Jmod;
                 Flocal[index]   += tracX * fact1 ;
                 Flocal[index+1] += tracY * fact1 ;
              }
           }
        }

        // side #3
        if(!CompareDoubles(tracdata[2][0],0.0) || !CompareDoubles(tracdata[2][1],0.0))
        {
           double NN[p+1];

           for(gp=0;gp<nGP1;gp++)   // loop over Gauss points
           {
              //if( finite && followerLoadFlag )
                // NurbsShapeFunctions2DAlg2(surf1, startindex[0], startindex[1], gausspoints1[gp], 1.0, NN, J, dircos);
              //else
                 NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], gausspoints1[gp], 1.0, NN, J, dircos);

              Jmod = J * gaussweights1[gp] * thick;

              //   for axisymmetric problems compute radius
              if(axsy)
              {
                 if(finite && followerLoadFlag )
                    rad  = surf1->SurfacePoint(knotsAtGPs[2*((nGP2-1)*nGP1+gp)], vvalues[1]).CalcEuclid().x;
                 else
                    rad  = surf0->SurfacePoint(knotsAtGPs[2*((nGP2-1)*nGP1+gp)], vvalues[1]).CalcEuclid().x;

                 Jmod *= (twoPI * rad);
              }

              EP  = surf0->SurfacePoint(knotsAtGPs[2*((nGP2-1)*nGP1+gp)], vvalues[1]).CalcEuclid() ;

              rad = sqrt(EP.x*EP.x + EP.y*EP.y) ;

              theta = atan2(EP.y, EP.x);

              //printf(" dircos ... %14.10f \t %14.10f \n", dircos[0], dircos[1]);

              //tracX = tracdata[2][0] * (-dircos[0]) + tracdata[2][1] * (-dircos[1]);
              //tracY = tracdata[2][0] * (-dircos[1]) + tracdata[2][1] * (dircos[0]);

              tracX = tracdata[2][0] ;
              tracY = tracdata[2][1] ;

              //stress[0][0] = pwhole.stressXX(rad, theta);
              //stress[0][1] = pwhole.stressXY(rad, theta);
              //stress[1][0] = stress[0][1];
              //stress[1][1] = pwhole.stressYY(rad, theta);

              //tracX = -( stress[0][0] * dircos[0] + stress[0][1] * dircos[1] );
              //tracY = -( stress[1][0] * dircos[0] + stress[1][1] * dircos[1] );

              //printf(" tracX ... %12.6f \t %12.6f \n", tracX, tracY);

              for(ii=0;ii<=p;ii++)
              {
                index = ndof * ((p+1) * q + ii);
                fact1 = NN[ii] * Jmod;
                Flocal[index]   += tracX * fact1 ;
                Flocal[index+1] += tracY * fact1 ;
              }
           }
        }

        // side #4
        if(!CompareDoubles(tracdata[3][0],0.0) || !CompareDoubles(tracdata[3][1],0.0))
        {
           double NN[q+1];

           for(gp=0;gp<nGP2;gp++)   // loop over Gauss points
           {
              if( finite && followerLoadFlag )
                 NurbsShapeFunctions2DAlg2(surf1, startindex[0], startindex[1], -1.0, gausspoints2[gp],  NN, J, dircos);
              else
                 NurbsShapeFunctions2DAlg2(surf0, startindex[0], startindex[1], -1.0, gausspoints2[gp],  NN, J, dircos);

              Jmod = J * gaussweights2[gp] * thick;

              //   for axisymmetric problems compute radius
              if(axsy)
              {
                 if(finite && followerLoadFlag )
                    rad  = surf1->SurfacePoint(uvalues[0], knotsAtGPs[gp*nGP1*2+1]).CalcEuclid().x;
                 else
                    rad  = surf0->SurfacePoint(uvalues[0], knotsAtGPs[gp*nGP1*2+1]).CalcEuclid().x;
                 
                 Jmod *= (twoPI * rad);
              }

              tracX = -tracdata[3][0] * (-dircos[0]) + -tracdata[3][1] * (-dircos[1]);
              tracY = -tracdata[3][0] * (-dircos[1]) + -tracdata[3][1] * (dircos[0]);
              
              //cout << tracX << '\t' << tracY << '\t' << rad << '\t' << J << '\t' << Jmod << endl;

              for(ii=0;ii<=q;ii++)
              {
                 index = ndof * ((p+1) * ii);
                 fact1 = NN[ii] * Jmod;
                 Flocal[index]   += tracX * fact1;
                 Flocal[index+1] += tracY * fact1;
              }

           }
        }
    }
//cout << " AAAAAAAAAA " << endl;
//printForceVector();
//
  return 0;
}
//





void NurbsElementSolid::diffStiffTest(double ddd, int dig, int dig2, bool gfrmt)
{
  int    i, jnd, jj, j, k, ii;

  double dd[6]  = {-3.*ddd, -2.*ddd, -ddd, +ddd, +2.*ddd, +3.*ddd },
         *r     = new double[6*nsize];

  ListArray<VectorArray<double> > sdiff, sdiff2;

  sdiff.setDim(nsize);
  sdiff2.setDim(nsize);
    for(int ii=0;ii<nsize;ii++)
    {
      stiffness_local[ii].zero();
      sdiff[ii].setDim(nsize);
      sdiff2[ii].setDim(nsize);
      sdiff[ii].zero();
      sdiff2[ii].zero();
    }

  int loc_num=0;
  vector<int>  uindex(nlbf);

  for(jj=0; jj<=surf0->q; jj++)
  {
    for(ii=0; ii<=surf0->p; ii++)
    {
      uindex[loc_num] = (startindex[0]+ii) + (surf0->ngbf1 * (startindex[1]+jj));
      loc_num++;
    }
  }
    for(ii=0; ii<nlbf; ii++)
       cout << uindex[ii] << '\t';
  cout << endl;


  // loop over columns of s

  for (jnd=0; jnd<nlbf; jnd++) // nodes
  {
    for (jj=0; jj<ndof; jj++)  // dof
    {
      j = jnd*ndof + jj;

      // loop over perturbations

      for (k=0; k<6; k++)
      {
        // apply pertubation
        surf1->updateCoordsSingleCP(uindex[jnd], jj, dd[k]);

        surf1->computeNET();

        // calculate residual
        calcStiffnessAndResidual();

        // remove pertubation

        surf1->updateCoordsSingleCP(uindex[jnd], jj, -dd[k]);

        // loop over rows of s, store residual

        for (i=0; i<nsize; i++) r[i*6+k] = Flocal[i];
      }

      // loop over rows of s
      for (i=0; i<nsize; i++)
        sdiff[j][i] = ( +       r[i*6+0]
                           -  9. * r[i*6+1]
                           + 45. * r[i*6+2]
                           - 45. * r[i*6+3]
                           +  9. * r[i*6+4]
                           -       r[i*6+5] ) / (60. * ddd);

    }
  }


  // calculate stiffness
  calcStiffnessAndResidual();


  cout.setf(ios::fixed);
  cout.setf(ios::showpoint);
  cout.precision(5);

  cout << "    Analytical Stiffness Matrix   " << endl;
  cout << endl;
  cout << endl;
  for(ii=0;ii<nsize;ii++)
  {
    cout << '\t' ;
    for(jj=0;jj<nsize;jj++)
    {
      cout  <<  fixed << stiffness_local[ii][jj] << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;


  cout << "    diff Stiffness Matrix   " << endl;
  cout << endl;
  cout << endl;
  for(ii=0;ii<nsize;ii++)
  {
    cout << '\t' ;
    for(jj=0;jj<nsize;jj++)
    {
      sdiff2[ii][jj] = stiffness_local[ii][jj] - sdiff[ii][jj];
      cout  <<  fixed <<  sdiff[ii][jj] << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;

  cout << "    difference of Stiffness Matrices   " << endl;
  cout << endl;
  cout << endl;
   for(ii=0;ii<nsize;ii++)
  {
    cout << '\t' ;
    for(jj=0;jj<nsize;jj++)
    {
      cout  <<  fixed << sdiff2[ii][jj] << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;


/*
  prgCompareTwoSimpleMatrices(sdiff,                         // matrix 1
		              s,                             // matrix 2
		              "numerical differentiation",   // title matrix 1
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      nsize, nsize,                      // matrix dimension
			      dig,dig2,gfrmt,                // format
			      0,                             // indentation
			      false,                         // interactive
			      false);                        // row/column numbers
*/
  delete [] r;
  //delete [] sdiff;

  return;
}





double NurbsElementSolid::volume(bool flag)
{
  double  Jac, dvol0, volume=0.0;
//
   int  count1, ii, jj, index, gp1, gp2;

   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf];

   double *gaussweights = &(surf0->gaussweights[0]);

   count1 = 0;
   volume = 0.0;
   for(gp2=0;gp2<nGP2;gp2++)
   {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*2;

          surf0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, Jac);

          dvol0 = Jac * gaussweights[count1] * JacMultFact;

          volume += dvol0;

          count1++;
     }
   }
//
  return volume;
}

