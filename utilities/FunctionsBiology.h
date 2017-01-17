
#ifndef Functions_Biology_h
#define Functions_Biology_h


#include <iostream>
//#include <limits.h>
//#include <float.h>
//#include <math.h>
#include <cmath>

#include "myConstants.h"
#include "Functions.h"

//using namespace Eigen;
//using namespace std;

using std::cout;
using std::endl;




//class  Poisson1DEx1 : public Function

class  Poisson1DEx1
{
  public:
    double  u0, u1, a, b, c;

    Poisson1DEx1()
    {
      u0=0.0; u1=0.0; a=1.0; b=-10.0; c=10.0;
    }

    virtual ~Poisson1DEx1() {}

    virtual double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  (u0 + (u1-u0-a/2.0-b/12.0-c/30.0)*xx + a*xx*xx/2.0 + b*xx*xx*xx*xx/12.0 + c*xx*xx*xx*xx*xx*xx/30.0);
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  (a + b*xx*xx + c*xx*xx*xx*xx);
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = (u1-u0-a/2.0-b/12.0-c/30.0) + a*xx + b*xx*xx*xx/3.0 + c*xx*xx*xx*xx*xx/5.0;

      return ;
    }
};



class  Fisher1DEx1
{
  public:
    double  mu, a, b, c, rho1, rho2;

    Fisher1DEx1()
    {
      mu=1.0; a=1.0; b=-10.0; c=10.0;
      rho1 = 1.0;
      rho2 = 1.0;
    }

    virtual ~Fisher1DEx1() {}

    virtual double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  ( a*xx*xx/2.0 + b*xx*xx*xx*xx/12.0 + c*xx*xx*xx*xx*xx*xx/30.0 );
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      double u = computeValue(0, xx);
      return  ( -mu*(a+b*xx*xx+c*xx*xx*xx*xx)-rho1*u+rho2*u*u);
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = a*xx + b*xx*xx*xx/3.0 + c*xx*xx*xx*xx*xx/5.0;

      return ;
    }
};


class  Fisher1DEx2
{
  public:
    double  mu, a, rho1, rho2;

    Fisher1DEx2()
    {
      mu=1.0; a=1.0;
      rho1 = 1.0;
      rho2 = 1.0;
    }

    virtual ~Fisher1DEx2() {}

    virtual double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  ( a*xx*sin(3.0*PI*xx) );
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      double u = computeValue(0, xx);
      
      return  ( -mu*(6.0*PI*a*cos(3.0*PI*xx)-9.0*PI*PI*xx*sin(3.0*PI*xx))-rho1*u+rho2*u*u);
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = a*sin(3.0*PI*xx) + 3.0*PI*a*xx*cos(3.0*PI*xx);

      return ;
    }
};




class  EFK1DEx1
{
  public:
    double  mu, a, rho1, rho2, gamm;

    EFK1DEx1()
    {
      mu=1.0; a=1.0;
      rho1 = 1.0;
      rho2 = 1.0;
      gamm = 0.01;
    }

    virtual ~EFK1DEx1() {}

    virtual double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  ( a*xx*sin(3.0*PI*xx) );
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      double u = computeValue(0, xx);
      
      return  ( -mu*(6.0*PI*a*cos(3.0*PI*xx)-9.0*PI*PI*xx*sin(3.0*PI*xx))-rho1*u+rho2*u*u);
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = a*sin(3.0*PI*xx) + 3.0*PI*a*xx*cos(3.0*PI*xx);

      return ;
    }
};





class  FK_unsteady_Ex1
{
  public:
    double  mu, rho;

    FK_unsteady_Ex1(double m1, double r1)
    {
      mu  = m1;
      rho = r1;
    }

    virtual ~FK_unsteady_Ex1() {}

    virtual double  computeValue(int dir, double tt=0.0, double xx=0.0)
    {
      double  val = 1.0+exp( (sqrt(rho/6.0))*xx - (5.0*rho/6.0)*tt);

      return  (1.0/val/val);
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  0.0;
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = 0.0;

      return ;
    }
};



class  FK_unsteady_Ex2
{
  public:
    double  mu, rho1, rho2;

    FK_unsteady_Ex2()
    {
      mu  = 1.0;
      rho1 = 0.5;
      rho2 = 1.0;
    }

    FK_unsteady_Ex2(double m1, double r1, double r2)
    {
      mu  = m1;
      rho1 = r1;
      rho2 = r2;
    }

    virtual ~FK_unsteady_Ex2() {}

    virtual double  computeValue(int dir, double tt=0.0, double xx=0.0)
    {
      double  val = -sqrt(rho1/24.0)*xx + 5.0*rho1*tt/12.0;

      return  (-rho1/4.0/rho2*(1.0/(cosh(val)*cosh(val)) - 2.0*tanh(val) - 2.0) );
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  0.0;
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = 0.0;

      return ;
    }
};




class FK2DsteadyEx1 : public Function
{
  public:
    
    double  mu, a, b, rho;
    
    FK2DsteadyEx1(double m1, double r1) {mu=m1; rho=r1; a=2.0; b=2.0;}

    virtual ~FK2DsteadyEx1() {}

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      return  (sin(a*PI*xx) * sin(b*PI*yy));
    }
    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double zz=0.0)
    {
      double val=sin(a*PI*xx) * sin(b*PI*yy);
      
      return  (val*(mu*(a*a+b*b)*PI*PI - rho) + rho*val*val*val);
    }
    virtual  void computeDerivatives(double xx, double yy, double* dv)
    {
      dv[0] = (a*PI*cos(a*PI*xx) * sin(b*PI*yy));
      dv[1] = (b*PI*sin(a*PI*xx) * cos(b*PI*yy));

      return ;
    }
};





#endif

