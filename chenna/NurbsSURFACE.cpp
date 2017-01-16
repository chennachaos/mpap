/*=============================================================================
        File: NurbsSURFACE.cpp
  Created by: Chennakesava Kadapa          (08 Jan 2011)
 Purpose    : Implementation file for the definitions of NURBS SURFACE Class

 ============================================================================*/

#include "NurbsSURFACE.h"
#include "DataBlockTemplate.h"
#include "Plot.h"
#include <iomanip>
#include "QuadratureUtil.h"

using namespace std;

extern Plot plot;
extern MpapTime mpapTime;



NurbsSURFACE::NurbsSURFACE(CNET& Pw1, KNOTVECTOR& U1, KNOTVECTOR& V1, DEGREE p1, DEGREE q1)
{
  Pw = Pw1;
  U  = U1;
  V  = V1;
  p  = p1;
  q  = q1;
}




NurbsSURFACE::NurbsSURFACE(const NurbsSURFACE& surf1)
{
  Pw = surf1.Pw;
  U  = surf1.U;
  V  = surf1.V;
  p  = surf1.p;
  q  = surf1.q;
}





NurbsSURFACE::~NurbsSURFACE()
{
   //cout << "    NurbsSURFACE: destructor ...\n\n";
}



NurbsSURFACE& NurbsSURFACE::operator=(const NurbsSURFACE &rhs)
{
  this->Pw = rhs.Pw;
  this->U  = rhs.U;
  this->V  = rhs.V;
  this->p  = rhs.p;
  this->q  = rhs.q;

  return *this;
}


void NurbsSURFACE::initializeBCdata()
{
   int ii, jj;

   PLATE_BENDING = false;

   tracflag.setDim(4);
   for(ii=0;ii<4;ii++)
   tracflag[ii] = false;

   VectorArray<double> X1, X2;
   findunique(U, X1);
   findunique(V, X2);

   nelem1 = X1.n-1; //no. of elements in first direction
   nelem2 = X2.n-1; //no. of elements in second direction

   nelem = nelem1*nelem2;

   ngbf1 = U.n - p - 1;
   ngbf2 = V.n - q - 1;

   ngbf = ngbf1 * ngbf2;
   
   cout << nelem1 << '\t' << nelem2 << '\t' << ngbf1 << '\t' << ngbf2 << endl;

   nlbf1 = p + 1;
   nlbf2 = q + 1;
   nlbf  = nlbf1 * nlbf2;

   nsize = nlbf*ndof;

   dispBCs.setDim(ngbf);

   gbfnums.setDim(ngbf);

   for(ii=0;ii<ngbf;ii++)
   {
      gbfnums[ii] = -5555;

      dispBCs[ii].setDim(ndof);
      for(jj=0;jj<ndof;jj++)
         dispBCs[ii][jj] = -7777; // some fancy number
   }

   forceBCs = dispBCs;
   
   closed.setDim(2);
   closed[0] = closed[1] = false;

  // check if the surface is closed

/*
  closed1 = closed2 = false;
  for(ii=0;ii<ngbf2;ii++)
  {
    if( Pw[ngbf1-1][ii] == Pw[0][ii] )
      closed1 = true;
  }

  for(ii=0;ii<ngbf1;ii++)
  {
    if( Pw[ii][ngbf2-1] == Pw[ii][0] )
      closed2 = true;
  }

  // cout << '\t' << " closed1 && closed2 " << closed1 << '\t' << closed2 << endl;
*/

  Uinit.setDim(ndof*ngbf);
  
  Uinit.zero();

  Values.setDim(ndof);
  ValuesPrev.setDim(ndof);
  ValuesPrev2.setDim(ndof);
  ValuesPrev3.setDim(ndof);
  ValuesPrev4.setDim(ndof);

  ValuesDot.setDim(ndof);
  ValuesDotDot.setDim(ndof);
  ValuesDotPrev.setDim(ndof);
  ValuesDotDotPrev.setDim(ndof);

  ValuesCur.setDim(ndof);
  ValuesDotCur.setDim(ndof);
  ValuesDotDotCur.setDim(ndof);


  for(ii=0;ii<ndof;ii++)
  {
    Values[ii].setDim(ngbf);           Values[ii].zero();
    ValuesPrev[ii].setDim(ngbf);       ValuesPrev[ii].zero();
    ValuesPrev2[ii].setDim(ngbf);      ValuesPrev2[ii].zero();
    ValuesPrev3[ii].setDim(ngbf);      ValuesPrev3[ii].zero();
    ValuesPrev4[ii].setDim(ngbf);      ValuesPrev4[ii].zero();

    ValuesDot[ii].setDim(ngbf);        ValuesDot[ii].zero();
    ValuesDotPrev[ii].setDim(ngbf);    ValuesDotPrev[ii].zero();
    ValuesDotDot[ii].setDim(ngbf);     ValuesDotDot[ii].zero();
    ValuesDotDotPrev[ii].setDim(ngbf); ValuesDotDotPrev[ii].zero();


    ValuesCur[ii].setDim(ngbf);        ValuesCur[ii].zero();
    ValuesDotCur[ii].setDim(ngbf);     ValuesDotCur[ii].zero();
    ValuesDotDotCur[ii].setDim(ngbf);  ValuesDotDotCur[ii].zero();
  }


  //for(ii=0;ii<ngbf;ii++)
    //ValuesDotPrev[0][ii] = 10.0;


  toUfull.setDim(ndof*ngbf);
  toUfull.zero();

  dispBCdata.setDim(8);
  dispBCdata.zero();

  edgedata.setDim(4);
  for(ii=0;ii<4;ii++)
     edgedata[ii] = -1;
  intfdata.setDim(8);
  for(ii=0;ii<8;ii++)
      intfdata[ii] = -1;

   PP.setDim(ngbf1);
   for(ii=0;ii<ngbf1;ii++)
     PP[ii].setDim(ngbf2);


   dispdata.setDim(4);
   for(ii=0;ii<4;ii++)
   {
     dispdata[ii].setDim(2);

     dispdata[ii][0] = 7777;
     dispdata[ii][1] = 7777;
   }

   dispdata[0][0] = 0.0;
   dispdata[0][1] = 0.0;
   //dispdata[1][0] = -7777;
   //dispdata[1][1] = -7777;
   dispdata[2][0] =  0.0;
   dispdata[2][1] = -1.0;
   //dispdata[3][0] = -7777;
   //dispdata[4][1] = -7777;

   //cout << " AAAAAAAAAAAA " << endl;

   //    computeNET();

  return;
}


void NurbsSURFACE::updateValues(int ind, double* uu)
{
  ind -= 1;

  if(ind > ndof)
    cerr << " ERROR! in  NurbsSURFACE::updateValues .... " << endl;

  for(int ii=0;ii<ngbf;ii++)
    Values[ind][ii] += uu[toUfull[ii]];

   return;
}




void NurbsSURFACE::updateCoordinates(double* uu)
{
  int index, count, ii, jj;

  count=0;
  for(jj=0;jj<ngbf2;jj++)
  {
    for(ii=0;ii<ngbf1;ii++)
    {
      index = gbfnums[count]*ndof;

      //cout << jj << '\t' << ii << '\t' << count << '\t' << index << '\t' << Pw[ii][jj].x << '\t' << Pw[ii][jj].y << '\t' << Pw[ii][jj].w << endl;

      Pw[ii][jj].x += (uu[index]   * Pw[ii][jj].w);
      Pw[ii][jj].y += (uu[index+1] * Pw[ii][jj].w);

      //PP[ii][jj].x += (uu[index]   );
      //PP[ii][jj].y += (uu[index+1] );

      count++;
    }
  }

  if( closed[0] )
  {
    int  ind = ngbf1-1;
    for(jj=0;jj<ngbf2;jj++)
      Pw[ind][jj] = Pw[0][jj];
  }

  if( closed[1] )
  {
    int  ind = ngbf2-1;
    for(ii=0;ii<ngbf1;ii++)
      Pw[ii][ind] = Pw[ii][0];
  }

  computeNET();

  return;
}




void NurbsSURFACE::geomToVector(double* outp)
{
  int  ind1, ind3, ii, jj;
  EPOINT EP;

  ind1=0;
  for(jj=0;jj<ngbf2;jj++)
  {
    for(ii=0;ii<ngbf1;ii++)
    {
      EP = Pw[ii][jj].CalcEuclid();

      ind3 = ndof*gbfnums[ind1];
               
      outp[ind3]   = EP.x;
      outp[ind3+1] = EP.y;
      ind1++;
    }
  }
  return;
}




void NurbsSURFACE::resetGeometry(double* uu)
{
  int index, count, ii, jj;

  double  temp;

  count=0;
  for(jj=0;jj<ngbf2;jj++)
  {
    for(ii=0;ii<ngbf1;ii++)
    {
      index = gbfnums[count]*ndof;

      temp = Pw[ii][jj].w;

      Pw[ii][jj].x = (uu[index]   * temp);
      Pw[ii][jj].y = (uu[index+1] * temp);
      count++;
    }
  }

  computeNET();
  return;
}




void NurbsSURFACE::addInitDOFvalues()
{
  int  ii, jj, index, count=0;

  double temp;

  for(jj=0;jj<ngbf2;jj++)
  {
    for(ii=0;ii<ngbf1;ii++)
    {
      index = count*ndof;

      temp = mpapTime.dt * Pw[ii][jj].w;

      Pw[ii][jj].x += (temp * Uinit[index]);
      Pw[ii][jj].y += (temp * Uinit[index+1]);

      count++;
    }
  }

  return;
}



void NurbsSURFACE::computeNET()
{
  for(int jj=0;jj<ngbf2;jj++)
  {
    for(int ii=0;ii<ngbf1;ii++)
       PP[ii][jj] = Pw[ii][jj].CalcEuclid();
  }

  return;
}






int  NurbsSURFACE::GenerateConnectivityArrays1(int& ntoteqs)
{
    nGP1  = (int) ElemProp.data[0] ;
    nGP2  = (int) ElemProp.data[1] ;

    int e=0, a, ni, nj, ii, jj, kk, ll, ind, col;

    // get the Gauss Points and Weights

    getGaussPoints1D(nGP1, gausspoints1, gaussweights1);
    getGaussPoints1D(nGP2, gausspoints2, gaussweights2);

    gaussweights.resize(nGP1*nGP2);

    ind=0;
    for(jj=0;jj<nGP2;jj++)
    {
        for(ii=0;ii<nGP1;ii++)
           gaussweights[ind++] = gaussweights2[jj] * gaussweights1[ii];
    }

    // INC Array

    INC.setDim(2);
    for(ii=0;ii<2;ii++)
      INC[ii].setDim(ngbf);

    col=0;
    for(ii=0;ii<ngbf2;ii++)
    {
      for(jj=0;jj<ngbf1;jj++)
      {
         INC[0][col] = jj;
         INC[1][col] = ii;
         col++;
      }
    }

    // IEN Array

    IEN.setDim(nelem);
    for(ii=0;ii<nelem;ii++)
      IEN[ii].setDim(nlbf);

    for(ll=q;ll<ngbf2;ll++)
    {
       if( !CompareDoubles(V[ll], V[ll+1]) )
       {
          nj = ll - q;
          for(kk=p;kk<ngbf1;kk++)
          {
              if( !CompareDoubles(U[kk], U[kk+1]) )
              {
                  a=0;
                  ni = kk - p;
                  for(jj=0;jj<=q;jj++)
                  {
                      ind = ngbf1*(nj+jj) + ni;
                      for(ii=0;ii<=p;ii++)
                      {
                         IEN[e][a] = ind+ii;
                         a++;
                      }
                  }
                  e++;
              }
          }
       }
    }

    // ID Array

    ID.setDim(ndof);
    for(ii=0;ii<ndof;ii++)
      ID[ii].setDim(ngbf);
/*
    for(ii=0;ii<ngbf;ii++)
    {
        for(jj=0;jj<ndof;jj++)
        {
            if(dispBCs[ii][jj] == (int) -7777)
            {
               ID[jj][ii] = ntoteqs;
               ntoteqs++;
            }
            else
               ID[jj][ii] = -1;
        }
    }
*/
/*
    int  ind1, ind2;
    ind1 = ngbf1;
    ind2 = ngbf2;
    
    ind=0;
    for(jj=0;jj<ind2;jj++)
    {
        for(ii=0;ii<ind1;ii++)
        {
            for(kk=0;kk<ndof;kk++)
            {
                if(dispBCs[ind][kk] == (int) -7777)
                {
                   ID[kk][ind] = ntoteqs;
                   ntoteqs++;
                }
                else
                   ID[kk][ind] = -1;
            }
            ind++;
        }
    }
*/

     cout << closed[0] << '\t' << closed[1] << endl;

     int ii1, jj1;
     for(ii=0;ii<ngbf;ii++)
     {
         // check for the closedness in first parametric direction
         if( (ii+1)%ngbf1 == 0 )
         {
            ii1 = INC[0][ii];
            jj1 = INC[1][ii];

            if(closed[0])
            {
               for(kk=0;kk<ndof;kk++)
                  ID[kk][ii] = ID[kk][ii-ngbf1+1];
            }
            else  // check for the specified displacement BCs
            {
               for(jj=0;jj<ndof;jj++)
               {
                  if(dispBCs[ii][jj] == (int) -7777)
                  {
                     ID[jj][ii] = ntoteqs;
                     ntoteqs++;
                  }
                  else
                     ID[jj][ii] = -1;
               }
            }
         }
         else   // check for the specified displacement BCs
         {
            for(jj=0;jj<ndof;jj++)
            {
               if(dispBCs[ii][jj] == (int) -7777)
               {
                  ID[jj][ii] = ntoteqs;
                  ntoteqs++;
               }
               else
                  ID[jj][ii] = -1;
            }
         }
     }


    return  0;
}








void NurbsSURFACE::updateCoordsSingleCP(int num, int dir, double val)
{
   int ii, jj;

   ii = num%ngbf1;
   jj = num/ngbf1;

   if(dir == 0)
      Pw[ii][jj].x += (val * Pw[ii][jj].w) ;
   else if(dir == 1)
      Pw[ii][jj].y += (val * Pw[ii][jj].w) ;
   else
      Pw[ii][jj].z += (val * Pw[ii][jj].w) ;


return;
}





CPOINT NurbsSURFACE::SurfacePoint(double u, double v)
{
/*
  Evaluates the Surface point in 4D at parameters u and v
  return the 4D point at Sw(u,v)
  Algorithm A4.3 on pg#134 of the Nurbs Book
*/
      vector<double>  Nu(nlbf1), Nv(nlbf2);
      int row, col, ii, jj, ind;

      row = FindSpan(&(U[0]), U.n, p, u) - p;
      col = FindSpan(&(V[0]), V.n, q, v) - q;

      BasisFuns(&(U[0]), U.n, p, u, &Nu[0]) ;
      BasisFuns(&(V[0]), V.n, q, v, &Nv[0]) ;

      CPOINT  Sw, CPtemp; Sw = 0.0;
      for(jj=0;jj<=q;jj++)
      {
          CPtemp = 0.0;
          ind = col+jj;
          for(ii=0;ii<=p;ii++)
	     CPtemp = CPtemp + Nu[ii]*Pw[row+ii][ind];

          Sw = Sw + Nv[jj]*CPtemp;
      }

      return Sw;
}



EPOINT NurbsSURFACE::SurfacePoint2(double u, double v)
{
/*
  Evaluates the Surface point in 4D at parameters u and v
  return the 4D point at Sw(u,v)
  Algorithm A4.3 on pg#134 of the Nurbs Book
*/
      vector<double>  Nu(nlbf1), Nv(nlbf2);
      int row, col, ii, jj, ind;

      row = FindSpan(&(U[0]), U.n, p, u) - p;
      col = FindSpan(&(V[0]), V.n, q, v) - q;

      BasisFuns(&(U[0]), U.n, p, u, &Nu[0]) ;
      BasisFuns(&(V[0]), V.n, q, v, &Nv[0]) ;

      EPOINT  Sw, CPtemp; Sw = 0.0;
      for(jj=0;jj<=q;jj++)
      {
          CPtemp = 0.0;
          ind = col+jj;
          for(ii=0;ii<=p;ii++)
	     CPtemp = CPtemp + Nu[ii]*PP[row+ii][ind];

          Sw = Sw + Nv[jj]*CPtemp;
      }

      return Sw;
}



void NurbsSURFACE::CurveOnSurf(int dir, double value, NurbsCURVE* curve1)
{
//  Computes the curve on a Surface for a given u and direction

  if(dir == 1)
  {
    if(CompareDoubles(value,0.0))
    {
       curve1->Pw.setDim(ngbf2);
       curve1->U.setDim(V.n);

       for(int ii=0;ii<curve1->Pw.n;ii++)
         curve1->Pw[ii] = Pw[0][ii];

       curve1->U = V;
       curve1->p = q;
    }
    else if(CompareDoubles(value,1.0))
    {
       curve1->Pw.setDim(ngbf2);
       curve1->U.setDim(V.n);

       for(int ii=0;ii<curve1->Pw.n;ii++)
         curve1->Pw[ii] = Pw[ngbf1-1][ii];

       curve1->U = V;
       curve1->p = q;
    }
    else
    {
       double u = value;
       vector<double>  Nu(p+1);

       int uspan = FindSpan(&(U[0]), U.n, p, u) ;
       BasisFuns(&(U[0]), U.n, p, u, &Nu[0]) ;
       int row = uspan-p;

       ListArray<CPOINT> temp;
       temp.setDim(Pw[0].n);


       for(int ll=0;ll<Pw[0].n;ll++)
       {
          for(int kk=0;kk<=p;kk++)
	     temp[ll] = temp[ll] + Nu[kk]*Pw[row+kk][ll];
       }

       curve1->Pw.setDim(temp.n);
       curve1->U.setDim(V.n);

       for(int ii=0;ii<curve1->Pw.n;ii++)
         curve1->Pw[ii] = temp[ii];

       curve1->U = V;
       curve1->p = q;
    }

    return;
  }
  else
  {
    if(CompareDoubles(value,0.0))
    {
       curve1->Pw.setDim(ngbf1);
       curve1->U.setDim(U.n);

       for(int ii=0;ii<curve1->Pw.n;ii++)
         curve1->Pw[ii] = Pw[ii][0];

       curve1->U = U;
       curve1->p = p;
    }
    else if(CompareDoubles(value,1.0))
    {
       curve1->Pw.setDim(ngbf1);
       curve1->U.setDim(U.n);

       for(int ii=0;ii<curve1->Pw.n;ii++)
         curve1->Pw[ii] = Pw[ii][ngbf2-1];

       curve1->U = U;
       curve1->p = p;
    }
    else
    {
       double v = value;
       vector<double>  Nv(q+1) ;

       int vspan = FindSpan(&(V[0]), V.n, q, v) ;

       BasisFuns(&(V[0]), V.n, q, v, &Nv[0]) ;
       int coln = vspan-q;

       ListArray<CPOINT> temp;
       temp.setDim(p+1);

       for(int ll=0;ll<=Pw.n;ll++)
       {
         for(int kk=0;kk<=q;kk++)
           temp[ll] = temp[ll] + Pw[ll][coln+kk] * Nv[kk];
       }

       curve1->Pw.setDim(temp.n);
       curve1->U.setDim(U.n);

       for(int ii=0;ii<curve1->Pw.n;ii++)
         curve1->Pw[ii] = temp[ii];

       curve1->U = U;
       curve1->p = p;
    }
    return;
  }
}





void NurbsSURFACE::SurfDerPointHom(double u, double v, int d, ListArray<ListArray<CPOINT> >& SKL)
{

//  Evaluates the surface derivative in 4D at parameters u and v
//  return the 4D point at S(u,v)
//  Algorithm A3.6 on pg#111 of the Nurbs Book

  SKL.setDim(d+1);
  for(int ii=0;ii<=d;ii++)
    SKL[ii].setDim(d+1);


  int k, l, s, r, du, dv, dd, ii;
  du = min(d, p);
  for(k=p+1;k<=d;k++)
  {
    for(l=0;l<=d-k;l++)
      SKL[k][l] = 0.0;
  }

  dv = min(d, q);
  for(l=q+1;l<=d;l++)
  {
    for(k=0;k<=d-l;k++)
      SKL[k][l] = 0.0;
  }


  double** Nu = new double*[du+1];
  double** Nv = new double*[dv+1];

  for(ii=0;ii<=du;ii++)
     Nu[ii] = new double[nlbf1];

  for(ii=0;ii<=dv;ii++)
     Nv[ii] = new double[nlbf2];


  int uspan = FindSpan(&(U[0]), U.n, p, u) ;
  DersBasisFuns(&(U[0]), U.n, p, u, du, Nu);

  int vspan = FindSpan(&(V[0]), V.n, q, v) ;
  DersBasisFuns(&(V[0]), V.n, q, v, dv, Nv);


  ListArray<CPOINT> temp;
  temp.setDim(q+1);

  for(k=0;k<=du;k++)
  {
    for(s=0;s<=q;s++)
    {
      temp[s] = 0.0;
      for(r=0;r<=p;r++)
        temp[s] = temp[s] + Nu[k][r] * Pw[uspan-p+r][vspan-q+s];
    }
    dd = min(d-k, dv);
    for(l=0;l<=dd;l++)
    {
       SKL[k][l] = 0.0;
       for(s=0;s<=q;s++)
         SKL[k][l] = SKL[k][l] + Nv[l][s] * temp[s];
    }
  }


  for(ii=0;ii<=du;ii++)
     delete []  Nu[ii];

  for(ii=0;ii<=dv;ii++)
     delete []  Nv[ii];

  delete [] Nu;
  delete [] Nv;


  return;
}




void NurbsSURFACE::SurfDerPointRat(double u, double v, int d, ListArray<ListArray<EPOINT> >& SKL)
{

//  Evaluates the surface derivative in 3D at parameters u and v
//  return the 3D point at S(u,v)
//  Algorithm A4.4 on pg#137 of the Nurbs Book

  int k, l, i, j ;
  double temp=0.0;


  SKL.setDim(d+1);
  for(i=0;i<=d;i++)
    SKL[i].setDim(d+1);

  EPOINT pv, pv2 ;

  ListArray<ListArray<CPOINT> > ders;

  SurfDerPointHom(u,v,d,ders) ;

  for(k=0;k<=d;k++)
  {
     for(l=0;l<=d-k;l++)
     {
        pv.x = ders[k][l].x ;
        pv.y = ders[k][l].y ;
        pv.z = ders[k][l].z ;

        for(j=1;j<=l;j++)
	   pv = pv - Bin(l,j) * ders[0][j].w * SKL[k][l-j] ;

        for(i=1;i<=k;i++)
        {
           temp = Bin(k,i);

           pv = pv - temp * ders[i][0].w * SKL[k-i][l] ;

           pv2 = 0.0 ;

           for(j=1;j<=l;j++)
              pv2 = pv2 + Bin(l,j) * ders[i][j].w * SKL[k-i][l-j] ;

           pv = pv - temp*pv2 ;
        }

        SKL[k][l] = pv/ders[0][0].w ;
     }
  }


  return;
}



//
//
//
//
//
//

void NurbsSURFACE::SurfacePointInverse(CPOINT& CP, double u0, double v0, VectorArray<double>& value)
{

  EPOINT S, Su, Sv, Suv, Suu, Svv, res, CP_temp;
  ListArray<ListArray<EPOINT> > Sder;
  double num1, denom1, num2, denom2, val1, val2, u_new, v_new;
  double eps1 = 0.001, eps2 = 0.001;
  CPOINT dummy, temp;

  double temp1 = 0.0, temp2 = 0.0, temp3 = 0.0, detJ = 0.0, norm_res = 0.0;

//  u0 = 0.0; v0 = 0.05;

  double u_old, v_old; // initial guesses for Newton Iteration
  u_old = u0; v_old = v0;

  double J[2][2]= {{0.0, 0.0},{0.0, 0.0} }, invJ[2][2] = {{0.0, 0.0},{0.0, 0.0} };
  double f[2] = {0.0, 0.0};
  double delta[2] ={0.0, 0.0};

  int ii;

    cout << endl;
    cout << endl;

  for(ii=0;ii<100;ii++)
  {
    SurfDerPointRat(u_old, v_old, 2, Sder);

    Su = Sder[1][0];
    Sv = Sder[0][1];
    Suv = Sder[1][1];
    Suu = Sder[2][0];
    Svv = Sder[0][2];

    //   cout << "You are here CCC" << endl;

    //P_temp = Sder[0][0];

    S = SurfacePoint(u_old, v_old).CalcEuclid();

    CP_temp = CP.CalcEuclid();

    res = S - CP_temp;

  /*
    S.print2screen();
    Su.print2screen();
    Sv.print2screen();
    Suv.print2screen();
    Suu.print2screen();
    Svv.print2screen();
    res.print2screen();
  */

    temp1 = Su.Norm();
    temp2 = Sv.Norm();

    J[0][0] = temp1 * temp1 + DotProduct(res, Suu);
    J[0][1] = DotProduct(Su, Sv) + DotProduct(res, Suv);
    J[1][0] = J[0][1];
    J[1][1] = temp2 * temp2  + DotProduct(res, Svv);

    f[0] = -DotProduct(res, Su);
    f[1] = -DotProduct(res, Sv);

    detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];
    invJ[0][0] = J[1][1] / detJ;
    invJ[0][1] = - J[0][1] / detJ;
    invJ[1][0] = - J[1][0] / detJ;
    invJ[1][1] = J[0][0] / detJ;

    delta[0] = invJ[0][0] * f[0] + invJ[0][1] * f[1];
    delta[1] = invJ[1][0] * f[0] + invJ[1][1] * f[1];

    u_new = u_old + delta[0];
    v_new = v_old + delta[1];


    if(u_new < U[0])
      u_new = U[0];
    if(u_new > U[U.n-1])
      u_new = U[U.n-1];
    if(v_new < V[0])
      v_new = V[0];
    if(v_new > V[V.n-1])
      v_new = V[V.n-1];

    cout << "  u_new && v_new  " << u_new << "   " << v_new << endl;

    norm_res = res.Norm();
    num1 = abs(DotProduct(Su, res));
    denom1 = Su.Norm() * norm_res;
    num2 = abs(DotProduct(Sv, res));
    denom2 = Sv.Norm() * norm_res;

    val1 = num1/denom1;
    val2 = num2/denom2;

  /*
    cout << endl;
    cout << "    residue norm   " << norm_res << endl;
    cout << "    val1 & val2   " << val1 << '\t' << val2 << endl;
    cout << endl;
    cout << endl;
  */

    if( norm_res <= eps1 || val1 <= eps2 || val2 <= eps2)
      break;
    else
      u_old = u_new;  v_old = v_new;

  }//while(norm_res > eps1);


    if(u_new < U[0])
      u_new = U[0];
    if(u_new > U[U.n-1])
      u_new = U[U.n-1];
    if(v_new < V[0])
      v_new = V[0];
    if(v_new > V[V.n-1])
      v_new = V[V.n-1];

  value[0] = u_new;
  value[1] = v_new;

  return;
}





void NurbsSURFACE::PlotControlPoints(int col)
{
    EPOINT *PP1;

    double x1d[3];

    for(int ii=0;ii<ngbf1;ii++)
    {
        for(int jj=0;jj<ngbf2;jj++)
        {
           PP1 = &(PP[ii][jj]);

           x1d[0] = PP1->x; x1d[1] = PP1->y; x1d[2] = PP1->z;

           plot.point(x1d);
        }
    }

  return;
}





void NurbsSURFACE::writeToFile(MyString &fileName, int index)
{
    ofstream  fout;

    int ii, jj;
    char tmp[500];

    fout.open(fileName.asCharArray());

    if(!fout)
    {
       prgWarning(1,"IsogeometricFEM::::writeToFile","could not open file for writing!");
       return;
    }
    
    cout << index << endl;

    if(index == 1) // write mesh to file
    {
        VectorArray<double>  uu, vv;

        findunique(U, uu);
        findunique(V, vv);

        ListArray<ListArray<EPOINT> > S1;
        S1.setDim(uu.n);

        // Calculate the Points on the surface
        for(ii=0;ii<uu.n;ii++)
        {
           S1[ii].setDim(vv.n);
           for(jj=0;jj<vv.n;jj++)
              S1[ii][jj] = SurfacePoint(uu[ii], vv[jj]).CalcEuclid();
        }

        fout << '\t' << nelem2+1 << '\t' << nelem1+1 << '\t' << 0.0 << endl;

        // plot 'u' isolines

        for(ii=0;ii<uu.n;ii++)
        {
            if( IsKnot(&(U[0]), U.n ,uu[ii]) )
            {
               for(jj=0;jj<vv.n;jj++)
               {
                  sprintf(tmp," %14.8f", S1[ii][jj].x);
                  sprintf(&(tmp[strlen(tmp)])," %14.8f", S1[ii][jj].y);
                  sprintf(&(tmp[strlen(tmp)])," %14.8f", S1[ii][jj].z);
                  fout << tmp << "\n";

                  //fout << S1[ii][jj].x << setw(14) << S1[ii][jj].y << setw(14) << S1[ii][jj].z << endl ;
               }
	       fout << "\n\n";
            }
        }

        // plot 'v' isolines

        for(jj=0;jj<vv.n;jj++)
        {
           if( IsKnot(&(V[0]), V.n ,vv[jj]) )
           {
              for(ii=0;ii<uu.n;ii++)
              {
                 sprintf(tmp," %14.8f", S1[ii][jj].x);
                 sprintf(&(tmp[strlen(tmp)])," %14.8f", S1[ii][jj].y);
                 sprintf(&(tmp[strlen(tmp)])," %14.8f", S1[ii][jj].z);
                 fout << tmp << "\n";

                 //fout << S1[ii][jj].x << setw(14) << S1[ii][jj].y << setw(14) << S1[ii][jj].z << endl ;
              }
	       fout << "\n\n";
           }
        }
    }
    else if(index == 2) // write control points to file
    {
        for(jj=0;jj<ngbf2;jj++)
        {
            for(ii=0;ii<ngbf1;ii++)
            {
               sprintf(tmp," %12.8f", Pw[ii][jj].x);
               sprintf(&(tmp[strlen(tmp)])," %12.8f", Pw[ii][jj].y);
               sprintf(&(tmp[strlen(tmp)])," %12.8f", Pw[ii][jj].z);
               sprintf(&(tmp[strlen(tmp)])," %12.8f", Pw[ii][jj].w);
               fout << tmp << "\n";
            }
            fout << "\n";
            fout << "\n";
        }
  }
  else //write CPs, order, and KVs to a file
  {
     cout << " Not yet implemented ... " << endl;
  }

  fout.close();

  return;
}




void NurbsSURFACE::readSurfaceFromFile(MyString &fileName)
{
  char fct[] = "NurbsSURFACE::readSurfaceFromFile";

  ifstream Ffile;

  MyString line, tmpl;

  int nn, n1, n2, ind, ii, jj;

  Vector<double> dTmp;
  Vector<int>    iTmp, lTmp;
  MyStringList   sTmp;
  List<Vector<int> > lviTmp;
  List<Vector<double> > lviTmp2;

  DataBlockTemplate t1, t2;

  char tmp[30];

  MyStringList  names;
  names.addNew("ngbf1 and ngbf2", "Surface Coordinates" );

  Ffile.open(fileName.asCharArray());

  if (!Ffile) prgError(1,fct,"failed to open the input file.");

     line.read(Ffile).stripToMin();
  while( Ffile )
  {
     //line.read(Ffile).stripToMin();
     cout << '\t' << line << endl;

     switch(names.whichBegins(line))
     {
       case 0:  // read ngbf1 and ngbf2

                if (!prgReadLnBrkSepListVectorInt(Ffile,line,lviTmp))
                   prgError(2,fct," invalid number or invalid keyword! ");

                n1 = lviTmp[0][0];
                n2 = lviTmp[0][1];

                cout << '\t' << lviTmp.n << '\t' << lviTmp[0].n << endl;

                cout << '\t' << n1 << '\t' << n2 << endl;

                if (n1 != ngbf1 || n2 != ngbf2)
                   prgError(3,fct," speficied global basis functions do match ");

                break;

       case 1:  // read Surface Coordinates

                cout << '\t' << line << endl;

             if (!prgReadLnBrkSepListVectorDbl(Ffile,line,lviTmp2))
                prgError(4,fct," invalid number or invalid keyword! ");


                nn = lviTmp2.n;

                cout << '\t' << lviTmp2.n << '\t' << lviTmp2[0].n << endl;

                sprintf(tmp,"%df",4);

              //  if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);

                tmpl.append(tmp);

                t1.initialise(tmpl);
	         t2.initialise(tmp);
                t1.expandToMatch(t2);

                if (!t1.readBlock(Ffile,line,iTmp,dTmp,sTmp,lTmp))
	          prgError(4,fct,"data error in 'Surface coordinates'!");

                nn = dTmp.n / (4);

                cout << '\t' << nn << endl;

                if(nn != (ngbf1*ngbf2) )
                   prgError(5,fct," data mismatch in surface coordinates to be read ");

                for(int pp=0;pp<nn;pp++)
                {
                    ii = pp%ngbf1;
                    jj = pp/ngbf1;
                    Pw[ii][jj].x = dTmp[pp*4+0];
                    Pw[ii][jj].y = dTmp[pp*4+1];
                    Pw[ii][jj].z = dTmp[pp*4+2];
                    Pw[ii][jj].w = dTmp[pp*4+3];
                }

                break;

        default:

                   prgError(5,fct," invalid identifier in input file ");

	     break;
      }

   }

  Ffile.close();

  return;
}

