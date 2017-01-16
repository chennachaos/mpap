
#include "Debug.h"
#include "FunctionsProgram.h"
#include "PropertyTypeEnum.h"
//#include "Plot.h"
#include "NurbsShapeFunctions.h"
#include "NurbsElementSolid3D.h"
#include <assert.h>

using namespace std;

//extern Plot     plot;
extern MpapTime mpapTime;




NurbsElementSolid3D::NurbsElementSolid3D()
{
  if (debug) cout << " constructor NurbsElementSolid3D\n\n";

  // cout << "     NurbsElementSolid3D: constructor ...\n\n";
}





NurbsElementSolid3D::~NurbsElementSolid3D()
{
  if (debug)   cout << " destructor NurbsElementSolid3D\n\n";

//NurbsElement::~NurbsElement();

//  cout << "     NurbsElementSolid3D: destructor ...\n\n";
}






void NurbsElementSolid3D::initialiseDOFvalues()
{
    int ind1, ind2, ind4, ii, jj;

    int *tt;

    tt = &(solid0->IEN[elenum][0]);
    ind1 = 0;
    for(ii=0;ii<nlbf;ii++)
    {
       ind2 = ndof * tt[ii];
       for(jj=0;jj<ndof;jj++)
       {
          ind4 = ind2+jj;

          forassy[ind1] = ind4;
          primvar[ind1] = solid1->Uinit[ind4];

          ind1++;
       }
    }

  return;
}


void NurbsElementSolid3D::createTractionDataVariable()
{
    assert(tracflag == true);
    
    if(!(tracdata.n > 0))
    {
       tracdata.setDim(18);
       for(int ii=0;ii<18;ii++)
       {
          tracdata[ii].setDim(3);
          tracdata[ii].zero();
       }
    }

   return;
}

void NurbsElementSolid3D::prepareElemData()
{
    // set the data in the base class "NurbsElement"
    NurbsElement::prepareElemData();

    uvalues[2] = uvalues[1] - uvalues[0];
    vvalues[2] = vvalues[1] - vvalues[0];
    wvalues[2] = wvalues[1] - wvalues[0];

    JacMultFact = 0.125 * uvalues[2] * vvalues[2] * wvalues[2];

    vals2project.setDim(8);
    startindex.setDim(3);

    elmDat = &(solid0->ElemProp.data[0]);
    matDat = &(solid0->MatlProp.data[0]);

    nGP1               = (int) elmDat[0] ;
    nGP2               = (int) elmDat[1] ;
    nGP3               = (int) elmDat[2] ;
    nGP                = nGP1 * nGP2 * nGP3;
    finiteInt          = (int) elmDat[3] ;

    bforce[0]          = elmDat[4] ;
    bforce[1]          = elmDat[5] ;
    bforce[2]          = elmDat[6] ;
    rho0               = elmDat[7] ;
    followerLoadFlag   = (elmDat[8] == 1);
    matId              = solid0->MatlProp.id + 1;
    finite             = (finiteInt >= 1) ;

    if(solid0->ElemProp.id == 12)
    {
       int ii, jj, sizep;

       sizep = solid2->nsize;

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

    // Initialize element internal variables, if any are there
    NurbsElementSolid3D::initialiseIntVar();

  return;
}





void NurbsElementSolid3D::setnivGP()
{
  // get number of internal variable per Gauss point

   double dmy[10];

  int  dmyI[3],
       isw   = 1,
  	mDim  = matdim_(&matId);

  if(mDim != 3) prgError(1,"NurbsElementSolid3D::nivGP","invalid value of ndm!");

  matlib3d_(matDat,dmy,dmy,dmy,dmy,dmy,dmy,
                               &matId,&nivGP,&finiteInt,&isw,dmyI);

  return;
}




void NurbsElementSolid3D::initialiseIntVar()
{
  // set nivGP value
  setnivGP();

  if(nivGP > 0)
  {
    // allocate memory
    int n = nivGP * nGP;

    intVar1 = new double [n];
    intVar2 = new double [n];

    double dmy[10];

    int   l, dmyI[3], ll   = 0, isw  = 2, mDim = matdim_(&matId), gp1, gp2, gp3;

    for(gp3=0;gp3<nGP3;gp3++)
    {
      for(gp2=0;gp2<nGP2;gp2++)
      {
        for(gp1=0;gp1<nGP1;gp1++)
        {
          matlib3d_(matDat,dmy,dmy,dmy,&(intVar1[ll]),&(intVar2[ll]),dmy, &matId,&nivGP,&finiteInt,&isw,dmyI);

          ll += nivGP;
        }
      }
    }
  }

  return;
}




void NurbsElementSolid3D::initialiseKnotsAtGPs()
{
   int  ind1, ind2, ind3, ind4, ind5, ind6, count=0, index, gp1, gp2, gp3;

   knotsAtGPs.setDim(nGP*3);

   double  val1, val2, val3, val4, val5, val6, yc, zc;

   ind1 = solid0->IEN[elenum][0];

   startindex[0] = solid0->INC[0][ind1];
   startindex[1] = solid0->INC[1][ind1];
   startindex[2] = solid0->INC[2][ind1];

   ind1 = startindex[0] + solid0->p;
   ind2 = ind1+1;
   ind3 = startindex[1] + solid0->q;
   ind4 = ind3+1;
   ind5 = startindex[2] + solid0->r;
   ind6 = ind5+1;


   double *gausspoints1 = &(solid0->gausspoints1[0]);
   double *gausspoints2 = &(solid0->gausspoints2[0]);
   double *gausspoints3 = &(solid0->gausspoints3[0]);

   val1 = 0.5*uvalues[2];
   val2 = 0.5*(uvalues[1]+uvalues[0]);

   val3 = 0.5*vvalues[2];
   val4 = 0.5*(vvalues[1]+vvalues[0]);

   val5 = 0.5*wvalues[2];
   val6 = 0.5*(wvalues[1]+wvalues[0]);

   // loop over Gauss points
   for(gp3=0;gp3<nGP3;gp3++)
   {
      zc = val5 * gausspoints3[gp3] + val6;
      for(gp2=0;gp2<nGP2;gp2++)
      {
         yc = val3 * gausspoints2[gp2] + val4;
         for(gp1=0;gp1<nGP1;gp1++)
         {
            index = count*3;
            knotsAtGPs[index]   = val1 * gausspoints1[gp1] + val2;
            knotsAtGPs[index+1] = yc;
            knotsAtGPs[index+2] = zc;

            count++;
         }
      }
   }

//   cout << '\t' << " elenum # " << elenum << '\t' << val1 << '\t' << val2 << '\t' << val3 << '\t' << val4 << endl;
//   cout << '\t' << knotsAtGPs << endl;

  return;
}






void NurbsElementSolid3D::plotGaussPoints(int num, bool defFlg)
{
   cout << "  NurbsElementSolid3D::plotGaussPoints  ... Not implemented yet  " << endl;
return;
}





int NurbsElementSolid3D::calcLoadVector()
{
    Flocal.setZero();

    if(tracflag == true)
    {
      double *gausspoints1  = &(solid0->gausspoints1[0]);
      double *gausspoints2  = &(solid0->gausspoints2[0]);
      double *gausspoints3  = &(solid0->gausspoints3[0]);
      double *gaussweights1 = &(solid0->gaussweights1[0]);
      double *gaussweights2 = &(solid0->gaussweights2[0]);
      double *gaussweights3 = &(solid0->gaussweights3[0]);

      double  *gaussweights = &(solid0->gaussweights[0]);

      double   dvol0=0.0, fact, fact1, fact2, fact3, Jac, Jmod, dircos[3], tracX, tracY, tracZ;
      double   val1, val2, val3, multfact;

      int index=0, ind, gp1, gp2, gp3, ii, jj, kk, size_local, count;
      int  p = solid0->p, q = solid0->q, r = solid0->r, nlbf1, nlbf2, nlbf3;
      int  ngbf1, ngbf2, ngbf3, ngbf, ngbf1m2, nlbf1m2, ind1, ind2, ind3, local1, local2, nGP1m2;
      bool  type[18];
        
      for(ii=0;ii<18;ii++)
      {
        type[ii] = false;
        type[ii] = (!CompareDoubles(tracdata[ii][0],0.0) || !CompareDoubles(tracdata[ii][1],0.0) || !CompareDoubles(tracdata[ii][2],0.0));

        //cout << (ii+1) << '\t' << tracdata[ii][0] << '\t' << tracdata[ii][1] << '\t' << tracdata[ii][2] << '\t' << type[ii] << endl;
      }
        
      ngbf1 = solid0->ngbf1;
      ngbf2 = solid0->ngbf2;
      ngbf3 = solid0->ngbf3;

      nlbf1 = solid0->nlbf1;
      nlbf2 = solid0->nlbf2;
      nlbf3 = solid0->nlbf3;

      ngbf1m2 = solid0->ngbf1m2;
      nGP1m2  = solid0->nGP1m2;
      nlbf1m2 = solid0->nlbf1m2;
        
        
      CNET Pw1;
      ListArray<ListArray<EPOINT> >  SKL;
      EPOINT  EP, Normal, Su, Sv;

      // face #1
      if(type[0])
      {
        ind1 = 0;
        fact1 = tracdata[0][0];    fact2 = tracdata[0][1];     fact3 = tracdata[0][2];

        Pw1.setDim(ngbf1);
        for(ii=0;ii<ngbf1;ii++)
        {
          Pw1[ii].setDim(ngbf2);        
          for(jj=0;jj<ngbf2;jj++)
            Pw1[ii][jj] = solid0->Pw[ind1][ii][jj];
        }
        
        NurbsSURFACE  surf_temp(Pw1, solid0->U, solid0->V, p, q);
        surf_temp.ndof = 2;
        surf_temp.initializeBCdata();
        surf_temp.computeNET();
        
        size_local = nlbf1m2;
        vector<double>  NN(size_local);

        multfact = 0.25 * uvalues[2] * vvalues[2];

        count = 0;
        for(gp2=0;gp2<nGP2;gp2++)
        {
        for(gp1=0;gp1<nGP1;gp1++)
        {
            index = 3*count;
            count++;

            surf_temp.ShapeFunctions(knotsAtGPs[index], knotsAtGPs[index+1], &NN[0]);

            surf_temp.SurfDerPointRat(knotsAtGPs[index], knotsAtGPs[index+1], 1, SKL);

            CrossProduct(SKL[0][1], SKL[1][0], EP);

            Jmod = EP.Norm() * gaussweights1[gp1] * gaussweights2[gp2] * multfact;

            Normal = EP / EP.Norm();
            Sv = SKL[0][1] / SKL[0][1].Norm();
            Su = SKL[1][0] / SKL[1][0].Norm();

            //Normal.print2screen();

            tracX = fact1 * Normal.x + fact2 * Su.x + fact3 * Sv.x;
            tracY = fact1 * Normal.y + fact2 * Su.y + fact3 * Sv.y;
            tracZ = fact1 * Normal.z + fact2 * Su.z + fact3 * Sv.z;

            local2=0;
            for(jj=0;jj<=q;jj++)
            {
              for(ii=0;ii<=p;ii++)
              {
                ind2 = 3*local2;
                fact = NN[local2] * Jmod;
                Flocal[ind2]   += tracX * fact ;
                Flocal[ind2+1] += tracY * fact ;
                Flocal[ind2+2] += tracZ * fact ;
                local2++;
              }
            }
        }//gp1
        }//gp2
      }
      // face #2
      if(type[1])
      {
        ind1 = ngbf3-1;
        fact1 = tracdata[1][0];    fact2 = tracdata[1][1];     fact3 = tracdata[1][2];

        Pw1.setDim(ngbf1);
        for(ii=0;ii<ngbf1;ii++)
        {
          Pw1[ii].setDim(ngbf2);        
          for(jj=0;jj<ngbf2;jj++)
            Pw1[ii][jj] = solid0->Pw[ind1][ii][jj];
        }
        
        NurbsSURFACE  surf_temp(Pw1, solid0->U, solid0->V, p, q);
        surf_temp.ndof = 2;
        surf_temp.initializeBCdata();
        surf_temp.computeNET();
        
        size_local = nlbf1m2;
        vector<double>  NN(size_local);

        multfact = 0.25 * uvalues[2] * vvalues[2];

        count = nGP1m2 * (nGP3-1);
        for(gp2=0;gp2<nGP2;gp2++)
        {
          for(gp1=0;gp1<nGP1;gp1++)
          {
            index = 3*count;
            count++;

            surf_temp.ShapeFunctions(knotsAtGPs[index], knotsAtGPs[index+1], &NN[0]);

            surf_temp.SurfDerPointRat(knotsAtGPs[index], knotsAtGPs[index+1], 1, SKL);

            CrossProduct(SKL[0][1], SKL[1][0], EP);

            Jmod = EP.Norm() * gaussweights1[gp1] * gaussweights2[gp2] * multfact;

            Normal = EP / EP.Norm();
            Sv = SKL[0][1] / SKL[0][1].Norm();
            Su = SKL[1][0] / SKL[1][0].Norm();

            //Normal.print2screen();

            tracX = fact1 * Normal.x + fact2 * Su.x + fact3 * Sv.x;
            tracY = fact1 * Normal.y + fact2 * Su.y + fact3 * Sv.y;
            tracZ = fact1 * Normal.z + fact2 * Su.z + fact3 * Sv.z;

            //cout << tracX << '\t' << tracY << '\t' << tracZ << endl;

            local1 = nlbf1m2*r;
            local2=0;
            for(jj=0;jj<=q;jj++)
            {
              for(ii=0;ii<=p;ii++)
              {
                ind2 = 3*local1;
                fact = NN[local2++] * Jmod;
                Flocal[ind2]   += tracX * fact ;
                Flocal[ind2+1] += tracY * fact ;
                Flocal[ind2+2] += tracZ * fact ;
                local1++;
              }
            }
        }//gp1
        }//gp2
      }
      // face #3
      if( type[2] )
      {
        ind1 = 0;
        fact1 = tracdata[2][0];    fact2 = tracdata[2][1];     fact3 = tracdata[2][2];

        Pw1.setDim(ngbf1);
        for(ii=0;ii<ngbf1;ii++)
        {
          Pw1[ii].setDim(ngbf3);        
          for(jj=0;jj<ngbf3;jj++)
            Pw1[ii][jj] = solid0->Pw[jj][ii][ind1];
        }

        NurbsSURFACE  surf_temp(Pw1, solid0->U, solid0->W, p, r);
        surf_temp.ndof = 2;
        surf_temp.initializeBCdata();
        surf_temp.computeNET();

        size_local = nlbf1*nlbf3;
        vector<double>  NN(size_local);

        multfact = 0.25 * uvalues[2] * wvalues[2];

        count = 0;
        for(gp3=0;gp3<nGP3;gp3++)
        {
          count = nGP1m2*gp3;
          for(gp1=0;gp1<nGP1;gp1++)
          {
            index = 3*(count + gp1);

            surf_temp.ShapeFunctions(knotsAtGPs[index], knotsAtGPs[index+2], &NN[0]);

            surf_temp.SurfDerPointRat(knotsAtGPs[index], knotsAtGPs[index+2], 1, SKL);

            CrossProduct(SKL[0][1], SKL[1][0], EP);

            Jmod = EP.Norm() * gaussweights1[gp1] * gaussweights3[gp3] * multfact;

            //Normal = EP / sqrt(val1*val1*val2*val2 - val3*val3);
            Normal = EP / EP.Norm();
            Sv = SKL[0][1] / SKL[0][1].Norm();
            Su = SKL[1][0] / SKL[1][0].Norm();

            //Normal.print2screen();
            //Su.print2screen();
            //Sv.print2screen();

            tracX = fact1 * Normal.x + fact2 * Su.x + fact3 * Sv.x;
            tracY = fact1 * Normal.y + fact2 * Su.y + fact3 * Sv.y;
            tracZ = fact1 * Normal.z + fact2 * Su.z + fact3 * Sv.z;

            local2 = 0;
            for(jj=0;jj<=r;jj++)
            {
              ind1 = nlbf1m2*jj;
              for(ii=0;ii<=p;ii++)
              {
                ind2 = 3*(ind1 + ii);
                fact = NN[local2++] * Jmod;
                Flocal[ind2]   += tracX * fact ;
                Flocal[ind2+1] += tracY * fact ;
                Flocal[ind2+2] += tracZ * fact ;
              }
            }
        }//gp1
        }//gp2
      }

      // face  #4
      if(type[3])
      {
        ind1 = ngbf2-1;
        fact1 = tracdata[3][0];    fact2 = tracdata[3][1];     fact3 = tracdata[3][2];

        Pw1.setDim(ngbf1);
        for(ii=0;ii<ngbf1;ii++)
        {
          Pw1[ii].setDim(ngbf3);        
          for(jj=0;jj<ngbf3;jj++)
            Pw1[ii][jj] = solid0->Pw[jj][ii][ind1];
        }

        NurbsSURFACE  surf_temp(Pw1, solid0->U, solid0->W, solid0->p, solid0->r);
        surf_temp.ndof = 2;
        surf_temp.initializeBCdata();
        surf_temp.computeNET();
        
        size_local = nlbf1*nlbf3;
        vector<double>  NN(size_local);

        multfact = 0.25 * uvalues[2] * wvalues[2];

        count = 0;
        for(gp2=0;gp2<nGP3;gp2++)
        {
          count = (nGP1m2)*(gp2+1) - nGP1;
          for(gp1=0;gp1<nGP1;gp1++)
          {
            index = 3*(count + gp1);

            surf_temp.ShapeFunctions(knotsAtGPs[index], knotsAtGPs[index+2], &NN[0]);

            surf_temp.SurfDerPointRat(knotsAtGPs[index], knotsAtGPs[index+2], 1, SKL);

            CrossProduct(SKL[0][1], SKL[1][0], EP);

            Jmod = EP.Norm() * gaussweights1[gp1] * gaussweights3[gp2] * multfact;

            Normal = EP / EP.Norm();
            Sv = SKL[0][1] / SKL[0][1].Norm();
            Su = SKL[1][0] / SKL[1][0].Norm();

            tracX = fact1 * Normal.x + fact2 * Su.x + fact3 * Sv.x;
            tracY = fact1 * Normal.y + fact2 * Su.y + fact3 * Sv.y;
            tracZ = fact1 * Normal.z + fact2 * Su.z + fact3 * Sv.z;

            //cout << tracX << '\t' << tracY << '\t' << tracZ << '\t' << Jac << '\t' << Jmod << endl;  cout << endl;

            local2 = 0;
            for(jj=0;jj<=r;jj++)
            {
              ind1 = nlbf1m2*(jj+1) - p - 1;
              for(ii=0;ii<=p;ii++)
              {
                ind2 = 3*(ind1 + ii);
                fact = NN[local2++] * Jmod;
                Flocal[ind2]   += tracX * fact ;
                Flocal[ind2+1] += tracY * fact ;
                Flocal[ind2+2] += tracZ * fact ;
              }
            }
        }//gp1
        }//gp2
      }

      // face #5
      if( type[4] )
      {
        ind1 = 0;
        fact1 = tracdata[4][0];    fact2 = tracdata[4][1];     fact3 = tracdata[4][2];

        Pw1.setDim(ngbf3);
        for(ii=0;ii<ngbf3;ii++)
        {
          Pw1[ii].setDim(ngbf2);        
          for(jj=0;jj<ngbf2;jj++)
            Pw1[ii][jj] = solid0->Pw[ii][ind1][jj];
        }

        NurbsSURFACE  surf_temp(Pw1, solid0->W, solid0->V, r, q);
        surf_temp.ndof = 2;
        surf_temp.initializeBCdata();
        surf_temp.computeNET();
        
        size_local = (q+1)*(r+1);
        vector<double>  NN(size_local);

        multfact = 0.25 * vvalues[2] * wvalues[2];

        count = 0;
        for(gp2=0;gp2<nGP2;gp2++)
        {
          count = nGP1*gp2;
          for(gp3=0;gp3<nGP3;gp3++)
          {
            index = 3*(count + nGP1m2*gp3);

            surf_temp.ShapeFunctions(knotsAtGPs[index+2], knotsAtGPs[index+1], &NN[0]);

            surf_temp.SurfDerPointRat(knotsAtGPs[index+2], knotsAtGPs[index+1], 1, SKL);
                  
            CrossProduct(SKL[0][1], SKL[1][0], EP);
                    
            Jmod = EP.Norm() * gaussweights2[gp2] * gaussweights3[gp3] * multfact;
                  
            Normal = EP / EP.Norm();
            Sv = SKL[0][1] / SKL[0][1].Norm();
            Su = SKL[1][0] / SKL[1][0].Norm();
                  
            tracX = fact1 * Normal.x + fact2 * Su.x + fact3 * Sv.x;
            tracY = fact1 * Normal.y + fact2 * Su.y + fact3 * Sv.y;
            tracZ = fact1 * Normal.z + fact2 * Su.z + fact3 * Sv.z;
                  
            local2 = 0;
            for(jj=0;jj<=q;jj++)
            {
              ind1 = (p+1)*jj;
              for(ii=0;ii<=r;ii++)
              {
                ind2 = 3*(ind1 + nlbf1m2*ii);
                fact = NN[local2++] * Jmod;
                Flocal[ind2]   += tracX * fact ;
                Flocal[ind2+1] += tracY * fact ;
                Flocal[ind2+2] += tracZ * fact ;
              }
            }
        }//gp1
        }//gp2
      }
      // face #6
      if( type[5] )
      {
        ind1 = ngbf1-1;
        fact1 = tracdata[5][0];    fact2 = tracdata[5][1];     fact3 = tracdata[5][2];

        Pw1.setDim(ngbf3);
        for(ii=0;ii<ngbf3;ii++)
        {
          Pw1[ii].setDim(ngbf2);        
          for(jj=0;jj<ngbf2;jj++)
            Pw1[ii][jj] = solid0->Pw[ii][ind1][jj];
        }
        
        NurbsSURFACE  surf_temp(Pw1, solid0->W, solid0->V, r, q);
        surf_temp.ndof = 2;
        surf_temp.initializeBCdata();
        surf_temp.computeNET();
        
        size_local = (q+1)*(r+1);
        vector<double>  NN(size_local);

        multfact = 0.25 * vvalues[2] * wvalues[2];

        count = 0;
        for(gp2=0;gp2<nGP2;gp2++)
        {
          count = nGP1*(gp2+1) - 1;
          for(gp3=0;gp3<nGP3;gp3++)
          {
            index = 3*(count + nGP1m2*gp3);

            surf_temp.ShapeFunctions(knotsAtGPs[index+2], knotsAtGPs[index+1], &NN[0]);

            surf_temp.SurfDerPointRat(knotsAtGPs[index+2], knotsAtGPs[index+1], 1, SKL);

            CrossProduct(SKL[0][1], SKL[1][0], EP);

            Jmod = EP.Norm() * gaussweights2[gp2] * gaussweights3[gp3] * multfact;

            Normal = EP / EP.Norm();
            Sv = SKL[0][1] / SKL[0][1].Norm();
            Su = SKL[1][0] / SKL[1][0].Norm();

            //Normal.print2screen();

            tracX = fact1 * Normal.x + fact2 * Su.x + fact3 * Sv.x;
            tracY = fact1 * Normal.y + fact2 * Su.y + fact3 * Sv.y;
            tracZ = fact1 * Normal.z + fact2 * Su.z + fact3 * Sv.z;

            //cout << tracX << '\t' << tracY << '\t' << tracZ << endl;

            local2 = 0;
            for(jj=0;jj<=q;jj++)
            {
              ind1 = (p+1)*(jj+1) - 1;
              for(ii=0;ii<=r;ii++)
              {
                ind2 = 3*(ind1 + nlbf1m2*ii);
                fact = NN[local2++] * Jmod;
                Flocal[ind2]   += tracX * fact ;
                Flocal[ind2+1] += tracY * fact ;
                Flocal[ind2+2] += tracZ * fact ;
              }
            }
        }//gp1
        }//gp2
      }

      // edge #4
      if( type[9] )
      {
        ind1 = ngbf2-1;
        fact1 = tracdata[9][0];    fact2 = tracdata[9][1];     fact3 = tracdata[9][2];

        Pw1.setDim(ngbf1);
        for(ii=0;ii<ngbf1;ii++)
        {
          Pw1[ii].setDim(ngbf3);        
          for(jj=0;jj<ngbf3;jj++)
            Pw1[ii][jj] = solid0->Pw[jj][ii][ind1];
        }

        NurbsSURFACE  surf_temp(Pw1, solid0->U, solid0->W, solid0->p, solid0->r);
        surf_temp.ndof = 2;
        surf_temp.initializeBCdata();
        surf_temp.computeNET();

        size_local = nlbf1;
        vector<double>  NN(size_local);

        multfact = 0.5 * uvalues[2];

        count = nGP1m2 - nGP1;
        for(gp1=0;gp1<nGP1;gp1++)
        {
          index = 3*(count + gp1);

          BasisFuns(&(solid0->U[0]), solid0->U.n, p, knotsAtGPs[index], &NN[0]) ;

          surf_temp.SurfDerPointRat(knotsAtGPs[index], 0.0, 1, SKL);

          CrossProduct(SKL[0][1], SKL[1][0], EP);                  

          //EP.print2screen();                    SKL[1][0].print2screen();                    SKL[0][1].print2screen();

          Jmod = SKL[1][0].Norm() * gaussweights1[gp1] * multfact;

          Normal = EP / EP.Norm();

          tracX = fact1 * Normal.x ;
          tracY = fact1 * Normal.y ;
          tracZ = fact1 * Normal.z ;

          //Normal.print2screen();

          //cout << tracX << '\t' << tracY << '\t' << tracZ << '\t' << Jmod << endl;  cout << endl;

          ind1 = nlbf1m2 - nlbf1;
          for(ii=0;ii<=p;ii++)
          {
            ind2 = 3*(ind1 + ii);
            fact = NN[ii] * Jmod;
            Flocal[ind2]   += tracX * fact ;
            Flocal[ind2+1] += tracY * fact ;
            Flocal[ind2+2] += tracZ * fact ;
          }
        }//gp1
      }

      // edge  #11
      if(type[16])
      {
        ind1 = ngbf2-1;
        fact1 = tracdata[16][0];    fact2 = tracdata[16][1];     fact3 = tracdata[16][2];

        Pw1.setDim(ngbf1);
        for(ii=0;ii<ngbf1;ii++)
        {
          Pw1[ii].setDim(ngbf3);        
          for(jj=0;jj<ngbf3;jj++)
            Pw1[ii][jj] = solid0->Pw[jj][ii][ind1];
        }

        NurbsSURFACE  surf_temp(Pw1, solid0->U, solid0->W, solid0->p, solid0->r);
        surf_temp.ndof = 2;
        surf_temp.initializeBCdata();
        surf_temp.computeNET();
        
        size_local = nlbf1*nlbf3;
        vector<double>  NN(size_local);

        multfact = 0.5 * wvalues[2];

        count = 0;
        for(gp3=0;gp3<nGP3;gp3++)
        {
          index = 3*( (nGP1m2)*(gp3+1) - nGP1 );

          BasisFuns(&(solid0->W[0]), solid0->W.n, r, knotsAtGPs[index+2], &NN[0]) ;

          surf_temp.SurfDerPointRat(0.0, knotsAtGPs[index+2], 1, SKL);

          CrossProduct(SKL[0][1], SKL[1][0], EP);                  

          //EP.print2screen();
          //SKL[1][0].print2screen();
          //SKL[0][1].print2screen();

          Jmod = SKL[1][0].Norm() * gaussweights3[gp3] * multfact;

          Normal = EP / EP.Norm();
          Sv = SKL[0][1] / SKL[0][1].Norm();
          Su = SKL[1][0] / SKL[1][0].Norm();

          tracX = fact1 * Normal.x + fact2 * Su.x + fact3 * Sv.x;
          tracY = fact1 * Normal.y + fact2 * Su.y + fact3 * Sv.y;
          tracZ = fact1 * Normal.z + fact2 * Su.z + fact3 * Sv.z;

          Normal.print2screen();

          //cout << tracX << '\t' << tracY << '\t' << tracZ << '\t' << Jac << '\t' << Jmod << endl;  cout << endl;

          local2 = 0;
          for(jj=0;jj<=r;jj++)
          {
            ind2 = 3*( nlbf1m2*(jj+1) - nlbf1 );

            fact = NN[local2++] * Jmod;
            Flocal[ind2]   += tracX * fact ;
            Flocal[ind2+1] += tracY * fact ;
            Flocal[ind2+2] += tracZ * fact ;
          }
        }//gp2
      }
    }

  return 0;
}




void NurbsElementSolid3D::diffStiffTest(double ddd, int dig, int dig2, bool gfrmt)
{
  int    i, jnd, ii, jj, j, k, kk;

  double dd[6]  = {-3.*ddd, -2.*ddd, -ddd, +ddd, +2.*ddd, +3.*ddd },
         *r     = new double[6*nsize];

  ListArray<VectorArray<double> > sdiff, sdiff2;

  sdiff.setDim(nsize);
  sdiff2.setDim(nsize);
   for(ii=0;ii<nsize;ii++)
   {
      stiffness_local[ii].zero();
      sdiff[ii].setDim(nsize);
      sdiff2[ii].setDim(nsize);
      sdiff[ii].zero();
      sdiff2[ii].zero();
   }

  int  loc_num=0, nn, ind1, ind2, ind3;
  vector<int>   uindex(nlbf);

  nn = solid0->ngbf1 * solid0->ngbf2;
  
  for(kk=0;kk<=solid0->r;kk++)
  {
     ind3 = nn*(startindex[2]+kk);
     for(jj=0;jj<=solid0->q;jj++)
     {
        ind2 = ind3 + solid0->ngbf1 * (startindex[1]+jj) + startindex[0];
        for(ii=0;ii<=solid0->p;ii++)
        {
           uindex[loc_num] = ind2 + ii ;
           loc_num++;
        }
     }
  }
  cout << " Global basis function numbers for this element " << endl;
    for(ii=0; ii<nlbf; ii++)
       cout << uindex[ii] << '\t';
  cout << endl;

  // loop over columns of s

  for(jnd=0;jnd<nlbf;jnd++) // nodes
  {
    for(jj=0;jj<ndof;jj++)  // dof
    {
      j = jnd*ndof + jj;
      
      // loop over perturbations
	    
      for(k=0; k<6; k++)
      {
        // apply pertubation
        solid1->updateCoordsSingleCP(uindex[jnd], jj, dd[k]);

        solid1->computeNET();

        // calculate residual
        calcStiffnessAndResidual();

        // remove pertubation

        solid1->updateCoordsSingleCP(uindex[jnd], jj, -dd[k]);

        // loop over rows of s, store residual

        for (i=0; i<nsize; i++)
          r[i*6+k] = resi[i];
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




double NurbsElementSolid3D::volume(bool flag)
{
  double  Jac, dvol0, volume=0.0;
//
   int  count1, ii, jj, index, gp1, gp2, gp3;

   double  N[nlbf], dN_dx[nlbf], dN_dy[nlbf], dN_dz[nlbf];

   double *gaussweights = &(solid0->gaussweights[0]);

  count1 = 0;
  volume = 0.0;
  for(gp2=0;gp2<nGP2;gp2++)
  {
    for(gp2=0;gp2<nGP2;gp2++)
    {
      for(gp1=0;gp1<nGP1;gp1++)
      {
          index = count1*3;

          solid0->ShapeFunDerivatives(&(startindex[0]), &(knotsAtGPs[index]), N, dN_dx, dN_dy, dN_dz, Jac);

          dvol0 = Jac * gaussweights[count1] * JacMultFact;

          volume += dvol0;

          count1++;
      }
    }
  }
//
  return volume;
}


