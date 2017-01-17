
#ifndef Functions_h
#define Functions_h


#include <iostream>
//#include <limits.h>
//#include <float.h>
//#include <math.h>
//#include <cmath>

#include "myConstants.h"

//using namespace Eigen;
//using namespace std;

using std::cout;
using std::endl;


//typedef  M_PI PI;


class  Function
{
  public:

    Function() {}

    virtual ~Function() {}

    virtual  double  computeValue(int dir, double xx = 0.0, double yy = 0.0, double zz = 0.0) { return -1111.0; };

    virtual  double  computeForce(int dir, double xx = 0.0, double yy = 0.0, double zz = 0.0) { return -1111.0; };

    virtual  void computeDerivatives(double xx, double yy, double* dv)
    { cout << " computeDerivatives() is not defined for this function " << endl; }

    virtual  void computeDerivatives2(double xx, double yy, double* dv)
    { cout << " computeDerivatives2() is not defined for this function " << endl; }

};


class ConstantFunction : public Function
{
  private: 
    double  val;
  
  public:

    ConstantFunction(double vv)
    { val = vv; }

    virtual ~ConstantFunction();
    
    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {return val; }
};



class ChannelFlowXaxis : public Function
{
  private: 
    double  ll, ul, fact;
  
  public:

    ChannelFlowXaxis(double aa, double bb, double cc)
    { ll = aa; ul = bb; fact = cc; }

    virtual ~ChannelFlowXaxis();
    
    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    { 
      if(dir == 0)
      {
        if(yy < ll || yy > ul)
          return 0.0;
        else
          return fact*(ul-yy)*(yy-ll);
      }
      else
        return 0.0;
    }
};



class ChannelFlowYaxis : public Function
{
  private: 
    double  ll, ul, fact;
  
  public:

    ChannelFlowYaxis(double aa, double bb, double cc)
    { ll = aa; ul = bb; fact = cc; }

    virtual ~ChannelFlowYaxis();
    
    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    { 
      if(dir == 1)
      {
        if(xx < ll || xx > ul)
          return 0.0;
        else
          return fact*(ul-xx)*(xx-ll);
      }
      else
        return 0.0;
    }
};




class ChannelFlowZaxis : public Function
{
  private: 
    double  ll, ul, fact;
  
  public:

    ChannelFlowZaxis(double aa, double bb, double cc)
    { ll = aa; ul = bb; fact = cc; }

    virtual ~ChannelFlowZaxis();
    
    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    { 
      if(dir == 2)
      {
        if(zz < ll || zz > ul)
          return 0.0;
        else
          return fact*(ul-zz)*(zz-ll);
      }
      else
        return 0.0;
    }
};

class PipeFlowXaxis : public Function
{
  private: 
    double  Rad, fact;
  
  public:

    PipeFlowXaxis(double aa, double bb)
    { Rad = aa; fact = bb;}

    virtual ~PipeFlowXaxis();

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    { 
      if(dir == 0)
      {
        double  ff = sqrt(yy*yy+zz*zz);
        if(ff >= Rad)
          return 0.0;
        else
          return fact*(1.0-ff*ff/Rad/Rad);
      }
      else
        return 0.0;
    }
};



class PipeFlowYaxis : public Function
{
  private: 
    double  Rad, fact;
  
  public:

    PipeFlowYaxis(double aa, double bb)
    { Rad = aa; fact = bb;}

    virtual ~PipeFlowYaxis();

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    { 
      if(dir == 1)
      {
        double  ff = sqrt(xx*xx+zz*zz);
        if(ff >= Rad)
          return 0.0;
        else
          return fact*(1.0-ff*ff/Rad/Rad);
      }
      else
        return 0.0;
    }
};


class PipeFlowZaxis : public Function
{
  private: 
    double  Rad, fact;
  
  public:

    PipeFlowZaxis(double aa, double bb)
    { Rad = aa; fact = bb;}

    virtual ~PipeFlowZaxis();

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    { 
      if(dir == 2)
      {
        double  ff = sqrt(xx*xx+yy*yy);
        if(ff >= Rad)
          return 0.0;
        else
          return fact*(1.0-ff*ff/Rad/Rad);
      }
      else
        return 0.0;
    }
};



class  AdvDiffExact1D : public Function
{
  public:
    double  L, a, mu, u0, u1;

    AdvDiffExact1D():L(1.0), a(1.0), mu(0.01), u0(0.0), u1(1.0) {}

    AdvDiffExact1D(double L1, double a1, double mu1, double v1, double v2)
    {
      L=L1; a=a1; mu=mu1; u0=v1; u1=v2;
    }

    virtual ~AdvDiffExact1D(){}

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      double  alpha = a/mu, fact;

      fact = (u1-u0)/(exp(alpha*L)-1.0);

      return  (u0+fact*(exp(alpha*xx)-1.0));
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      double  alpha = a/mu, fact;

      fact = (u1-u0)/(exp(alpha*L)-1.0);

      dv[0] = ((fact*a/mu)*exp(alpha*xx));

      return ;
    }
};




class  HemkerExact : public Function
{
  public:
    double  k, x0, L, fact;
   
    HemkerExact()
    {
      L = 2.0;
      k = 10.0e10;
      x0 = -1.0;
      fact = erf(1.0/sqrt(2.0*k));
    }

    virtual ~HemkerExact(){}

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  (cos(PI*xx) + erf(xx/sqrt(2.0*k))/fact);
    }

    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  (-k*PI*PI*cos(PI*xx) - PI*xx*sin(PI*xx) );
    }
};



class PoissonEx1 : public Function
{
  public:
    
    double  a, b;
    
    PoissonEx1() {a=2.0;b=2.0;}

    virtual ~PoissonEx1() {}

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  (sin(a*PI*xx) * sin(b*PI*yy));
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  ((a*a+b*b)*PI*PI*sin(a*PI*xx) * sin(b*PI*yy));
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = (a*PI*cos(a*PI*xx) * sin(b*PI*yy));
      dv[1] = (b*PI*sin(a*PI*xx) * cos(b*PI*yy));

      return ;
    }
};



class  PoissonEx2 : public Function
{
  public:
    PoissonEx2() {}

    virtual ~PoissonEx2() {}

    virtual double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  (xx*xx+yy*yy+zz*zz);
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      //return  (4.0*(xx*xx*xx*xx+yy*yy*yy*yy+zz*zz*zz*zz));
      return  6.0;
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = (2.0*xx);
      dv[1] = (2.0*yy);
      //dv[2] = (2.0*zz);

      return ;
    }

};




class  PoissonEx3 : public Function
{
  public:
    PoissonEx3() {}

    virtual ~PoissonEx3() {}

    virtual double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  (cosh(PI*yy) - sinh(PI*yy)/tanh(PI))*sin(PI*xx);
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  0.0;
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = (cosh(PI*yy) - sinh(PI*yy)/tanh(PI))*PI*cos(PI*xx);
      dv[1] = (PI*sinh(PI*yy) - PI*cosh(PI*yy)/tanh(PI))*sin(PI*xx);

      return ;
    }
};



class  PoissonInterfaceEx1 : public Function
{
  public:

    PoissonInterfaceEx1() {}

    virtual ~PoissonInterfaceEx1() {}

    virtual double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
        double  r = sqrt(xx*xx+yy*yy);
        if( r <= 0.5)
          return 1.0;
        else
          return  (1.0 + log(2.0*r));
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  0.0;
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      double  r = sqrt(xx*xx+yy*yy);

      if( r <= 0.5)
      {
        dv[0] = 0.0;
        dv[1] = 0.0;
      }
      else
      {
        dv[0] = xx/(r*r);
        dv[1] = yy/(r*r);
      }

      return ;
    }
};





class  PoissonInterfaceEx2 : public Function
{
  public:
    double  C, b;
    
    PoissonInterfaceEx2(double b1, double C1)
    {C=C1; b=b1;}

    virtual ~PoissonInterfaceEx2() {}

    virtual double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
        double  r = sqrt(xx*xx+yy*yy);
        if( r <= 0.5)
          return r*r;
        else
          return  0.25*(1.0-1.0/8.0/b-1.0/b) + (0.5*r*r*r*r+r*r)/b + C*log(2.0*r)/b;
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  8.0*(xx*xx+yy*yy)+4.0;
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      double  r = sqrt(xx*xx+yy*yy);

      if( r <= 0.5)
      {
        dv[0] = 2.0*xx;
        dv[1] = 2.0*yy;
      }
      else
      {
        dv[0] = 2.0*xx*(r*r+1.0)/b + C*xx/(b*r*r);
        dv[1] = 2.0*yy*(r*r+1.0)/b + C*yy/(b*r*r);
      }

      return ;
    }
};




class  PoissonInterfaceEx3 : public Function
{
  public:
  
    double  k1, k2, jumpi, val1, val2;

    PoissonInterfaceEx3(double a, double b, double c)
    {
      k1=a; k2=b; jumpi=c;

      val1 = 3.0*k1+k2;
      val2 = 4.0*k1*(k1+k2);
    }

    virtual ~PoissonInterfaceEx3() {}

    virtual double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if( xx <= 0.5)
        return (val1*xx/val2 - xx*xx/2.0/k1 );
      else
        return  ( (k2-k1+val1*xx)/val2 - xx*xx/2.0/k2 + jumpi );
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  1.0;
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      if( xx <= 0.5)
      {
        dv[0] = val1/val2 - xx/k1 ;
      }
      else
      {
        dv[0] = val1/val2 - xx/k2 ;
      }

      dv[1] = 0.0;

      return ;
    }
};



class  PoissonInterfaceEx4 : public Function
{
  public:
  
    double  R;

    PoissonInterfaceEx4(double a=0.5):R(a)
    {    }

    virtual ~PoissonInterfaceEx4() {}

    virtual double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      double rad = sqrt(xx*xx + yy*yy);
      if( rad >= R)
        return 0.0;
      else
        return  ( (rad*rad*rad-R*R*R)/9.0 );
    }

    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  (-sqrt(xx*xx + yy*yy) );
    }

    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      double rad = sqrt(xx*xx + yy*yy);
      if( rad >= R)
      {
        dv[0] = 0.0 ;
        dv[1] = 0.0 ;
      }
      else
      {
        dv[0] = xx*rad/3.0 ;
        dv[1] = yy*rad/3.0 ;
      }
      return ;
    }
};




class  AdvDiffExact2DEx1 : public Function
{
  public:

    double  a, theta, ax, ay;

    AdvDiffExact2DEx1():a(1.0), theta(35.0)
    {
      ax = 1.0;
      ay = tan(35.0*PI/180.0);
    }

    virtual ~AdvDiffExact2DEx1() {}

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if( yy > ay*xx )
        return 2.0;
      else
        return 1.0;
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  0.0;
    }
};




class  Stokes2DEx1 : public Function
{
  public:

    Stokes2DEx1() {}

    virtual ~Stokes2DEx1() {}

    double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if(dir == 0)
        return  ( (12.0*(1.0-2.0*yy)*(xx-2.0)*pow(xx,3.0)) + ((12.0-48.0*yy+72.0*yy*yy-48.0*yy*yy*yy)*pow(xx,2.0)) +  ((-2.0+24.0*yy-72*yy*yy+48*yy*yy*yy)*xx) + (1.0-4.0*yy+12.0*yy*yy-8.0*yy*yy*yy) );
      else //if(dir == 1)
        return  ( ((8.0-48.0*yy+48.0*yy*yy)*pow(xx,3.0)) + ((-12.0+72.0*yy-72.0*yy*yy)*pow(xx,2.0)) + ((4.0-24.0*yy+48*yy*yy-48*yy*yy*yy+24.0*yy*yy*yy*yy)*xx) + (-12.0*yy*yy+24.0*yy*yy*yy-12.0*yy*yy*yy*yy) );
    }

    double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if(dir == 0)
        return  ( pow((xx - xx*xx),2.0)*(2.0*yy - 6.0*yy*yy + 4.0*yy*yy*yy) );
      else if(dir == 1)
        return  ( -pow((yy - yy*yy),2.0)*(2.0*xx - 6.0*xx*xx + 4.0*xx*xx*xx) );
      else //if(dir == 2)
        return  (xx*(1.0-xx));
    }

    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = (2.0*xx-6.0*xx*xx+4.0*xx*xx*xx)*(2.0*yy-6.0*yy*yy+4.0*yy*yy*yy);
      dv[2] = (xx*xx-2.0*xx*xx*xx+xx*xx*xx*xx)*(2.0-12.0*yy+12.0*yy*yy);
      dv[1] = -(2.0-12.0*xx+12.0*xx*xx)*(yy*yy-2.0*yy*yy*yy+yy*yy*yy*yy);
      dv[3] = -(2.0*xx-6.0*xx*xx+4.0*xx*xx*xx)*(2.0*yy-6.0*yy*yy+4.0*yy*yy*yy);

      return ;
    }
};



class  Stokes2DEx2 : public Function
{
  public:

    double  nu, a, b;

    Stokes2DEx2(){nu=1.0;a=2.0*PI; b=3.0*PI;}

    virtual ~Stokes2DEx2() {}

    double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if(dir == 0)
        return  ( 2.0*nu*a*a*sin(a*xx)*cos(a*yy) + b*cos(b*xx)*cos(b*yy));
      else //if(dir == 1)
        return  (-2.0*nu*a*a*cos(a*xx)*sin(a*yy) - b*sin(b*xx)*sin(b*yy));
    }

    double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if(dir == 0)
        return  ( sin(a*xx)*cos(a*yy) );
      else if(dir == 1)
        return  ( -cos(a*xx)*sin(a*yy) );
      else //if(dir == 2)
        return  (sin(b*xx)*cos(b*yy));
    }

    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] =  a*cos(a*xx)*cos(a*yy); // du/dx
      dv[2] = -a*sin(a*xx)*sin(a*yy); // du/dy
      dv[1] =  a*sin(a*xx)*sin(a*yy); // dv/dx
      dv[3] = -a*cos(a*xx)*cos(a*yy); // dv/dy

      return ;
    }
};





class  Stokes2DEx3 : public Function
{
  public:

    Stokes2DEx3() {}

    virtual ~Stokes2DEx3() {}

    double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if(dir == 0)
        return  ( 3.0*(xx*xx)*(yy*yy) - yy - 1.0 );
      else //if(dir == 1)
        return  ( 2*(xx*xx*xx*yy) + 3.0*xx - 1.0 );
    }

    double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if(dir == 0)
        return  ( xx*(1.0 + xx - 2.0*yy + xx*xx - 3.0*yy*yy + xx*yy) );
      else if(dir == 1)
        return  ( yy*(-1.0 - 2.0*xx + yy - 3.0*xx*xx + yy*yy - xx*yy) );
      else //if(dir == 2)
        return  ( xx*yy + xx + yy + xx*xx*xx*yy*yy );
    }

    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      // du/dx
      dv[0] = 1.0 + 2.0*xx - 2.0*yy + 3.0*xx*xx - 3.0*yy*yy + 2.0*xx*yy;
      // du/dy
      dv[2] = -2.0*xx - 6.0*xx*yy + xx*xx ;
      // dv/dx
      dv[1] = -2.0*yy - 6.0*xx*yy - yy*yy ;
      // dv/dy
      dv[3] = -1.0 - 2.0*xx + 2.0*yy - 3.0*xx*xx + 3.0*yy*yy - 2.0*xx*yy;

      return ;
    }
};





class  PearsonVortex : public Function
{
   public:
   
    double  nu;
   
    PearsonVortex(double dd):nu(dd) {}

    PearsonVortex(){nu=0.0001;}

    virtual ~PearsonVortex(){}

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if(dir == 0)
        return  ( -cos(PI*xx)*sin(PI*yy) );
      else if(dir == 1)
        return  (  sin(PI*xx)*cos(PI*yy) );
      else //if(dir == 2)
        return  -0.25*(cos(2.0*PI*xx)+cos(2.0*PI*yy));
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  0.0;
    }
};



class  Kovasznay : public Function
{
  private:
    double  Re, lamda, p0;

  public:
    Kovasznay()
    {
      Re = 40.0;
      lamda = 0.5*Re - sqrt(0.25*Re*Re + 4.0*PI*PI);
      p0 = 0.5*exp(-lamda);
    }

    Kovasznay(double Re1)
    {
      Re = Re1;

      lamda = 0.5*Re - sqrt(0.25*Re*Re + 4.0*PI*PI);

      p0 = 0.5*exp(-lamda);
    }

    virtual ~Kovasznay(){}

    void SetPressure(double pp)
    {
      p0 = pp;
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  0.0;
    }     

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      if(dir == 0)
        return ( 1.0 - exp(lamda*xx) * cos(2.0*PI*yy) );
      else if(dir == 1)
        return ( (lamda/2.0/PI) * exp(lamda*xx) * sin(2.0*PI*yy) );
      else //if(dir == 2)
        return  (p0 - 0.5*exp(2.0*lamda*xx));
    }

    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      double  fact = exp(lamda*xx), fact2 = 2.0*PI*yy;;

      dv[0] = - lamda * fact * cos(fact2);
      dv[2] = 2.0*PI*fact* sin(fact2);
      dv[1] = (lamda*lamda/2.0/PI) * fact * sin(fact2);
      dv[3] = lamda*fact* cos(fact2);

      return;
    }
};




class  Biharmonic1DEx1 : public Function
{
  public:

    double  fact, a, aPI;

    Biharmonic1DEx1(){ a=1.0; aPI=a*PI; fact = aPI*aPI*aPI*aPI; }

    virtual ~Biharmonic1DEx1() {}

    double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  ( fact * sin(PI*xx) );
    }

    double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  ( sin(aPI*xx) );
    }

    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = aPI*cos(aPI*xx);

      return ;
    }

    virtual  void computeDerivatives2(double xx, double yy, double* dv)
    {
      dv[0] = -aPI*aPI*sin(aPI*xx); // d2u_dx2

      return ;
    }

};





class  Biharmonic2DEx1 : public Function
{
  public:

    double  fact, a, b;

    Biharmonic2DEx1(){fact = 4.0*(PI*PI*PI*PI); a=1.0;b=1.0;}

    virtual ~Biharmonic2DEx1() {}

    double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  ( fact * sin(PI*xx) * sin(PI*yy) );
    }

    double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  ( sin(PI*xx)*sin(PI*yy) );
    }

    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = PI*cos(PI*xx)*sin(PI*yy);
      dv[1] = PI*sin(PI*xx)*cos(PI*yy);

      return ;
    }

    virtual  void computeDerivatives2(double xx, double yy, double* dv)
    {
      dv[0] = -PI*PI*sin(PI*xx)*sin(PI*yy); // d2u_dx2
      dv[1] = -PI*PI*sin(PI*xx)*sin(PI*yy); // d2u_dy2
      dv[2] =  PI*PI*cos(PI*xx)*cos(PI*yy); // d2u_dxdy

      return ;
    }

};





#endif


