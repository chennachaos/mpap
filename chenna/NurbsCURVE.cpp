/*=============================================================================
        File: NurbsCURVE.cpp
  Created by: Chennakesava Kadapa          (08 Jan 2011)
 Purpose    : Implementation file for the definitions of NURBS CURVE Class

 ============================================================================*/

#include "NurbsCURVE.h"
#include "Plot.h"
#include <iomanip>
#include "QuadratureUtil.h"

using namespace std;

extern Plot plot;
extern MpapTime mpapTime;



NurbsCURVE::NurbsCURVE()
{
}



NurbsCURVE::~NurbsCURVE()
{
  // cout << "    NurbsCURVE: destructor ...\n\n";
}



NurbsCURVE::NurbsCURVE(CPOLYGON& Pw1, KNOTVECTOR& U1, DEGREE p1)
{
  Pw = Pw1;
  U = U1;
  p = p1;
}



NurbsCURVE::NurbsCURVE(const NurbsCURVE& curv1)
{
  Pw = curv1.Pw;
  U = curv1.U;
  p = curv1.p;
}



NurbsCURVE& NurbsCURVE::operator=(const NurbsCURVE &rhs)
{
  this->Pw = rhs.Pw;
  this->U = rhs.U;
  this->p = rhs.p;

  return *this;
}




void NurbsCURVE::geomToVector(double* outp)
{
    EPOINT EP;

    for(int ii=0;ii<ngbf;ii++)
    {
        EP = Pw[ii].CalcEuclid();
        outp[ii]  = EP.x;
    }
}




void NurbsCURVE::resetGeometry(double* uu)
{
    int index, count, ii, jj;
     
    double  temp;

    for(ii=0;ii<ngbf;ii++)
    {
        index = gbfnums[ii]*ndof;
           
        temp = Pw[ii].w;

        Pw[ii].x = (uu[index] * temp);
    }

    computeNET();   
}







//
//
//
//
//
void NurbsCURVE::initializeBCdata()
{
  ngbf = Pw.n;

  nlbf = p+1;

  dispBCs.setDim(ngbf);
  for(int ii=0;ii<ngbf;ii++)
    dispBCs[ii].setDim(ndof);

  for(int ii=0;ii<ngbf;ii++)
  {
     for(int jj=0;jj<ndof;jj++)
     {
        dispBCs[ii][jj] = -7777; // some fancy number
     }
  }

  forceBCs = dispBCs;

  VectorArray<double> X1;
  findunique(U, X1);

  nelem = X1.n-1;

  Uinit.setDim(ngbf*ndof);
  Uinit.zero();

  gbfnums.setDim(ngbf);

  toUfull.setDim(ngbf*ndof);

  Values.setDim(ndof); // data for constraint variable
  
  for(int ii=0;ii<ndof;ii++)
  {
     Values[ii].setDim(ngbf);
     Values[ii].zero();
  }

  PP.setDim(ngbf);

  computeNET();

  return;
}





void NurbsCURVE::addInitDOFvalues()
{
    int index, count=0;

    double temp;

    if(ndof == 1)
    {
       for(int ii=0;ii<ngbf;ii++)
       {
          Pw[ii].y += (mpapTime.dt * Pw[ii].w * Uinit[ii]);
       }
    }
    else
    {
       for(int ii=0;ii<ngbf;ii++)
       {
          index = count*ndof;

          temp = mpapTime.dt * Pw[ii].w;

          Pw[ii].x += (temp * Uinit[index]);
          Pw[ii].y += (temp * Uinit[index+1]);

          count++;
       }
    }

   return;
}



void NurbsCURVE::computeNET()
{
    for(int jj=0;jj<ngbf;jj++)
       PP[jj] = Pw[jj].CalcEuclid();

    return;
}





void NurbsCURVE::updateCoordinates(double* uu)
{
    // works only for ndof = 1

    for(int ii=0;ii<ngbf;ii++)
    {
      Pw[ii].x += (uu[gbfnums[ii]] * Pw[ii].w);
      //printf("%12.6f \n", Pw[ii].x);
    }

    computeNET();

    return;
}





double  NurbsCURVE::computeValue(int dof, double u)
{
    double  val;
    vector<double>  Nu(nlbf);

    int ni, nj, temp, ii;

    ni = FindSpan(&(U[0]), U.n, p, u) - p ;

    BasisFuns(&(U[0]), U.n, p, u, &Nu[0]) ;

    val=0.0;
    dof -= 1;
    for(ii=0; ii<=p; ii++)
      val += Nu[ii] * Values[dof][ni+ii];

    return val;
}


void NurbsCURVE::updateCoordsSingleCP(int num, int dir, double val)
{
      Pw[num].y += (val * Pw[num].w) ;

return;
}


// compute connectivity arrays INC, IEN and ID
int NurbsCURVE::GenerateConnectivityArrays1(int& ntoteqs)
{
    nGP  = (int) ElemProp.data[0] ;

    // get the Gauss Points and Weights in the first direction
    getGaussPoints1D(nGP, gausspoints, gaussweights);


    int ii, jj, kk, e, ni, val1;

    // INC Array
    //-------------------

    INC.setDim(1);
    INC[0].setDim(ngbf);

    for(ii=0;ii<ngbf;ii++)
      INC[0][ii] = ii;


    // IEN Array
    //-------------------

    IEN.setDim(nelem);
    for(ii=0;ii<nelem;ii++)
      IEN[ii].setDim(nlbf);


    val1 = U.n - p -1;

    e = 0;
    for(kk=p;kk<val1;kk++)
    {
      if( !CompareDoubles(U[kk], U[kk+1]) )
      {
         ni = kk - p;
         for(ii=0;ii<=p;ii++)
            IEN[e][ii] = ni+ii;
         e++;
      }
    }

    // ID Array
    //-------------------

    ID.setDim(ndof);
    for(ii=0;ii<ndof;ii++)
      ID[ii].setDim(ngbf);

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

  return 0;
}





void NurbsCURVE::ShapeFunctions(double u, double* N)
{
   BasisFuns(&(U[0]), U.n, p, u, N) ;

   return;
}





void NurbsCURVE::ShapeFunDerivatives(int ni, double u, double* dN_dx, double& Jac)
{
   return;
}




void NurbsCURVE::ShapeFunsAndDerivatives(int ni, double u, double* N, double* dN_dx, double& Jac)
{
   int ROWS = 2, ii;

   double** ders = new double*[ROWS];

   for(ii=0;ii<ROWS;ii++)
     ders[ii] = new double[p+1];

   DersBasisFuns(&(U[0]), U.n, p, u, 1, ders);

   // Gradient of mapping from parameter space to physical space
   double dx_du=0;
   for(ii=0; ii<=p; ii++)
   {
      N[ii] = ders[0][ii];
      dx_du += ( Pw[ni+ii].CalcEuclid().x * ders[1][ii] );
   }

  // Compute inverse of gradient
  double du_dx = 1.0/dx_du;

  // Compute derivatives of basis functions w.r.t physical coordinates

  for(ii=0; ii<=p; ii++)
    dN_dx[ii] = ders[1][ii] * du_dx;

  // Gradient of mapping from parent element to parameter space
  double du_dtildeu = 0.5*(U[ni+p+1]-U[ni+p]);

  //Calculate Jacobian
  Jac = dx_du * du_dtildeu;


  for(ii=0;ii<ROWS;ii++)
    delete [] ders[ii];

  delete [] ders;
   return;
}



void NurbsCURVE::deformationGradient(int ni, bool confflag, double* dN_dx, double& F)
{
      //   confflag = 1 --> coordiantes from current configuration
      //                    shape function derivatives w.r.t reference configuration

      //   confflag = 2 ==> coordiantes from reference configuration
      //                    shape function derivatives w.r.t current configuration

      int  ii;

      EPOINT *EP;

      F = 0.0;
      for(ii=0; ii<=p; ii++)
      {
          EP = &(PP[ni+ii]);
          F += EP->x * dN_dx[ii];
      }

      if(confflag == 1) // calculation of F
        F = F;
      else //confflag == 2.  calculation of F^-1
        F = 1.0/F;

  return;
}





void NurbsCURVE::ShapeFunsAndDerivatives2(int ni, double u, double* N, double* dN_dx, double* d2N_dx2, double& Jac)
{
    int ROWS = 3, ii;

    double** ders = new double*[ROWS];

    for(ii=0;ii<ROWS;ii++)
       ders[ii] = new double[p+1];

    DersBasisFuns(&(U[0]), U.n, p, u, 2, ders);

    double   dx_du=0.0, d2x_du2=0.0;

    EPOINT  EP;
    for(ii=0; ii<=p; ii++)
    {
        N[ii] = ders[0][ii];
        EP = Pw[ni+ii].CalcEuclid();
        dx_du   += ( EP.x * ders[1][ii] );
        d2x_du2 += ( EP.x * ders[2][ii] );
    }

    double   du_dx = 1.0/dx_du;
    double   d2u_dx2 = -d2x_du2/pow((dx_du),3);

    for(ii=0; ii<=p; ii++)
    {
        dN_dx[ii]   = ders[1][ii] * du_dx;
        d2N_dx2[ii] = ders[2][ii] * du_dx * du_dx + ders[1][ii] * d2u_dx2;
    }

    Jac = dx_du * 0.5*(U[ni+p+1] - U[ni+p]);

    for(ii=0;ii<ROWS;ii++)
      delete [] ders[ii];

    delete [] ders;

   return;
}



CPOINT NurbsCURVE::CurvePoint(double u)
{
// Evaluates the curve in 4D at parameter u
// Algorithm A4.1 on pg#124 of the Nurbs Book

   vector<double>  Nb(p+1);

   int span = FindSpan(&(U[0]), U.n, p, u) ;

   BasisFuns(&(U[0]), U.n, p, u, &Nb[0]) ;

   CPOINT Cw(0,0,0,0);
   for(int j=0;j<=p;j++)
     Cw = Cw + Nb[j] * Pw[span-p+j] ;

   return Cw;
}


//
//
//
//
//

void NurbsCURVE::CurveDerPointHom(double u, int d, ListArray<CPOINT>& Cder)
{
//  Evaluates the curve derivative in 4D(homogeneous coordinates) at parameter u and
//  return the 4D point at C(u)
//  Algorithm A3.2 on pg#93 of the Nurbs Book

   int k, j, du, span, ii;

   du = min(d, p);

   //resize Cder
   Cder.setDim(d+1);

   double** nders = new double*[du+1];

   for(ii=0;ii<=du;ii++)
      nders[ii] = new double[p+1];

   span = FindSpan(&(U[0]), U.n, p, u) ;

   DersBasisFuns(&(U[0]), U.n, p, u, du, nders);

   for(k=p+1;k<=d;k++)
     Cder[k] = 0.0;

   for(k=0;k<=du;k++)
   {
      Cder[k] = 0.0;
      for(j=0;j<=p;j++)
        Cder[k] = Cder[k] + nders[k][j] * Pw[span-p+j] ;
   }

  for(ii=0;ii<=du;ii++)
     delete [] nders[ii];

  delete [] nders;

 return;

}


//
//
//
//
void NurbsCURVE::CurveDerPointRat(double u, int d, ListArray<EPOINT>& Cder)
{
//  Evaluates the curve derivative in 3D(projected) at parameter u and
//  Algorithm A4.2 on pg#127 of the Nurbs Book

   //resize Cder
   Cder.setDim(d+1);

   EPOINT pv;

   ListArray<CPOINT> ders;

   CurveDerPointHom(u, d, ders) ;

/*
  cout << '\t' << "    Homogeneous Coords " << endl;
  for(int ii=0;ii<=d;ii++)
    ders[ii].print2screen();
  cout << endl;
*/
   int k, i ;

   for(k=0;k<=d;k++)
   {
      pv.x = ders[k].x ;
      pv.y = ders[k].y ;
      pv.z = ders[k].z ;

      for(i=1;i<=k;i++)
        pv = pv - Bin(k,i)*ders[i].w * Cder[k-i] ;

      Cder[k] = pv/ders[0].w ;
  }

 return;

}


//
//
//
//
double NurbsCURVE::CurvePointInverse(EPOINT EP, double u0)
{
  ListArray<EPOINT> Cder;
  double num, denom, val, val1, val2, u_new, u_old, resnorm;
  double eps1 = 0.0001;

  EPOINT C, Cu, Cuu, res;

  u_new = u_old = u0;  // initial guess for Newton Iteration

  if(u_old<U[0])
    u_old = U[0] ;
  if(u_old>U[U.n-1])
    u_old = U[U.n-1] ;

  //cout << "         u_old : " << u_old << endl;

  int t=0;

  while(1)
  {
    ++t ;
    if(t>100)
      return u_old;

    //Cu = CurvePoint(u_old).CalcEuclid();

    CurveDerPointRat(u_old, 2, Cder);

    C = Cder[0];
    Cu = Cder[1];
    Cuu = Cder[2];

    res = C - EP;
/*
    C.print2screen();
    Cu.print2screen();
    Cuu.print2screen();
*/
    resnorm = res.Norm();

    num = DotProduct(Cu, res);
    val1 = Cu.Norm();
    val2 = DotProduct(Cuu, res);

    val = abs(num / (val1 * resnorm));

    if( res.Norm() <= eps1 && val <= eps1)
    {
      if(u_old<U[0])
        u_old = U[0] ;
      if(u_old>U[U.n-1])
        u_old = U[U.n-1] ;

       return u_old;
    }
    else
    {
       denom = val1 * val1 + val2;

       u_new = u_old - (num/denom);

    //   cout << "         val1 : " << val1 << endl;
    //   cout << "         u_new : " << u_new << endl;

       if(u_new<U[0])
         u_new = U[0] ;
       if(u_new>U[U.n-1])
         u_new = U[U.n-1] ;

       if( abs((u_new - u_old) * val1 ) <= eps1 )
         return u_new;
    }

    u_old = u_new;
  }

  return u_new;
}





void NurbsCURVE::PlotElements(int col, bool PLOT_KNOT_LINES, int* resln)
{
      VectorArray<double> uu;

      findunique(U, uu);

//      create_vector2(U, 10, uu);

      ListArray<EPOINT> C1;
      C1.setDim(uu.n);

      for(int ii=0;ii<uu.n;ii++)
        C1[ii] = CurvePoint(uu[ii]).CalcEuclid();



      double x1d[3], x2d[3];
      EPOINT *PP1, *PP2;
      for(int ii=0;ii<uu.n-1;ii++)
      {
        PP1 = &(C1[ii]);
        PP2 = &(C1[ii+1]);

        x1d[0] = PP1->x; x1d[1] = PP1->y; x1d[2] = PP1->z;
        x2d[0] = PP2->x; x2d[1] = PP2->y; x2d[2] = PP2->z;

        plot.line(x1d,x2d);
        //plot.point(x1d);        plot.point(x2d);
      }

  return;
}




void NurbsCURVE::PlotCurve(int npoints)
{
      VectorArray<double> uu;

      create_vector2(U, npoints, uu);

      ListArray<EPOINT> C1;      C1.setDim(uu.n);

      for(int ii=0;ii<uu.n;ii++)
        C1[ii] = CurvePoint(uu[ii]).CalcEuclid();

      double x1d[3], x2d[3];
      EPOINT *PP1, *PP2;
      for(int ii=0;ii<uu.n-1;ii++)
      {
         PP1 = &(C1[ii]);
         PP2 = &(C1[ii+1]);

         x1d[0] = PP1->x; x1d[1] = PP1->y; x1d[2] = PP1->z;
         x2d[0] = PP2->x; x2d[1] = PP2->y; x2d[2] = PP2->z;

         plot.line(x1d,x2d);
      }

  return;
}





void NurbsCURVE::printGeomToFile()
{
    VectorArray<double> uu;

    if(U.n > 30)
       findunique(U, uu);
    else
       create_vector2(U, 10, uu);


    ListArray<EPOINT> C1;      C1.setDim(uu.n);

    for(int ii=0;ii<uu.n;ii++)
      C1[ii] = CurvePoint(uu[ii]).CalcEuclid();

    ofstream fout("NURBS_Curve_Geom.dat");

    if(fout.fail())
    {
      cout << " Could not open the Output file" << endl;
      exit(1);
    }

	fout.setf(ios::fixed);
	fout.setf(ios::showpoint);
	fout.precision(6);

	for(int ii=0;ii<uu.n;ii++)
	{
           fout << C1[ii].x << setw(14) << C1[ii].y << setw(14) << C1[ii].z << endl ;
	}

	fout.close();

  return;
}



void NurbsCURVE::PlotControlPoints(int col)
{
    EPOINT PP1;
    double x1d[3];
    for(int ii=0;ii<Pw.n;ii++)
    {
        PP1 = Pw[ii].CalcEuclid();
        x1d[0] = PP1.x; x1d[1] = PP1.y; x1d[2] = PP1.z ;
        cout << PP1.x << '\t' << PP1.y << endl;
        plot.point(x1d);
    }

  return;
}



void NurbsCURVE::createAndWriteBasisFunctions(int index)
{
    cout << " Knot Vector ... " << endl;
    cout << U << endl;
    cout << endl;
    cout << endl;

    VectorArray<double>  uu;

    int m = U.n-1;
    int n = m-p-1, ii, jj, kk, span;

    create_vector2(U, 20, uu);

    ListArray<VectorArray<double> > NN;
    NN.setDim(n+1);

    for(ii=0;ii<=n;ii++)
    {
        NN[ii].setDim(uu.n);
	NN[ii].zero();
    }

    cout << " AAAAAAAAAAAA " << endl;

    for(ii=0;ii<=n;ii++)
    {
        for(jj=0;jj<uu.n;jj++)
        {
            NN[ii][jj] = OneBasisFun_recurs(U, p, ii, uu[jj]);
        }
    }

        cout << " AAAAAAAAAAAA " << endl;

    ofstream fout("NURBS_Basis_Functions.dat");

    if(fout.fail())
    {
        cout << " Could not open the Output file" << endl;
        exit(1);
    }

    fout.setf(ios::fixed);
    fout.setf(ios::showpoint);
    fout.precision(6);

    for(jj=0;jj<uu.n;jj++)
    {
        fout << setw(10) << uu[jj];
        for(ii=0;ii<=n;ii++)
            fout << setw(14) << NN[ii][jj];
        fout << endl;
    }

    fout.close();

    return;
}

void NurbsCURVE::createAndWriteBasisFunctionsDerivatives(int, int)
{

  return;
}

