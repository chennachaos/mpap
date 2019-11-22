
#ifndef KimMoinFlow_h
#define KimMoinFlow_h


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



//class  KimMoinFlow : public Function

/*
class  KimMoinFlow
{
  private:
    double rho, mu, a, c;

  public:
    KimMoinFlow(double r, double m, double c1)
    {
      rho = r;
      mu  = m;
      c   = c1;
      a   = 2.0;
    }

    virtual ~KimMoinFlow(){}

    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      double cx = cos(a*PI*xx);
      double sx = sin(a*PI*xx);
      double cy = cos(a*PI*yy);
      double sy = sin(a*PI*yy);
      double gu = exp(-2.0*tt);
      double gp = exp(-4.0*tt);

      //if(dir == 0)
        //return  ( 2.0*cx*sy*gu*(rho-mu*a*a*PI*PI) + a*PI*sx*cx*gp*(1.0-rho) );
      //else //if(dir == 1)
        //return  (-2.0*sx*cy*gu*(rho-mu*a*a*PI*PI) + a*PI*sy*cy*gp*(1.0-rho) );

      if(dir == 0)
        return  ( 2.0*cx*sy*gu*(rho-mu*a*a*PI*PI) + a*PI*sx*cx*gp );
      else //if(dir == 1)
        return  (-2.0*sx*cy*gu*(rho-mu*a*a*PI*PI) + a*PI*sy*cy*gp );
    }

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      if(dir == 0)
        return ( -cos(a*PI*xx)*sin(a*PI*yy)*exp(-2.0*tt) );
      else if(dir == 1)
        return (  sin(a*PI*xx)*cos(a*PI*yy)*exp(-2.0*tt) );
      else //if(dir == 2)
        return  ( -0.25*(cos(2.0*a*PI*xx)+cos(2.0*a*PI*yy))*exp(-4.0*tt) );
    }

    virtual  void computeDerivatives(double xx, double yy, double tt, double* dv)
    {
      double c1 = exp(-2.0*tt);
      
      dv[0] =  a*PI*sin(a*PI*xx) * sin(a*PI*yy)*c1; // du/dx
      dv[2] = -a*PI*cos(a*PI*xx) * cos(a*PI*yy)*c1; // du/dy
      dv[1] =  a*PI*cos(a*PI*xx) * cos(a*PI*yy)*c1; // dv/dx
      dv[3] = -a*PI*sin(a*PI*xx) * sin(a*PI*yy)*c1; // dv/dy

      return;
    }
};
*/


/*
class  KimMoinFlow
{
  // Unsteady Stokes

  private:
    double rho, mu;

  public:
    KimMoinFlow(double r, double m, double c1)
    {
      rho = r;
      mu  = m;
    }

    virtual ~KimMoinFlow(){}

    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      double cx = cos(xx);
      double sx = sin(xx);
      double cy = cos(yy);
      double sy = sin(yy);

      if(dir == 0)
        return  (-2.0*cx*sy*(rho*cos(2.0*tt)+mu*sin(2.0*tt)) + sx*cx*sin(2.0*tt)*sin(2.0*tt) );
      else //if(dir == 1)
        return  ( 2.0*sx*cy*(rho*cos(2.0*tt)+mu*sin(2.0*tt)) + sy*cy*sin(2.0*tt)*sin(2.0*tt) );
    }

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      if(dir == 0)
        return ( -cos(xx)*sin(yy)*sin(2.0*tt) );
      else if(dir == 1)
        return (  sin(xx)*cos(yy)*sin(2.0*tt) );
      else //if(dir == 2)
        return  ( -0.25*(cos(2.0*xx)+cos(2.0*yy))*sin(2.0*tt)*sin(2.0*tt) );
    }

    virtual  void computeDerivatives(double xx, double yy, double tt, double* dv)
    {
      double c1 = sin(2.0*tt);
      
      dv[0] =  sin(xx)*sin(yy)*c1; // du/dx
      dv[2] = -cos(xx)*cos(yy)*c1; // du/dy
      dv[1] =  cos(xx)*cos(yy)*c1; // dv/dx
      dv[3] = -sin(xx)*sin(yy)*c1; // dv/dy

      return;
    }
};
*/


/*
class  KimMoinFlow
{
  // For Steady Stokes  
  
  private:
    double rho, mu, a, c;

  public:
    KimMoinFlow(double r, double m, double c1)
    {
      rho = r;
      mu  = m;
    }

    virtual ~KimMoinFlow(){}

    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      double cx = cos(xx);
      double sx = sin(xx);
      double cy = cos(yy);
      double sy = sin(yy);

      if(dir == 0)
        return  (-2.0*mu*cx*sy + sx*cx );
      else //if(dir == 1)
        return  ( 2.0*mu*sx*cy + sy*cy );
    }

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      if(dir == 0)
        return ( -cos(xx)*sin(yy) );
      else if(dir == 1)
        return (  sin(xx)*cos(yy) );
      else //if(dir == 2)
        return  ( -0.25*(cos(2.0*xx)+cos(2.0*yy)) );
    }

    virtual  void computeDerivatives(double xx, double yy, double tt, double* dv)
    {
      dv[0] =  sin(xx)*sin(yy); // du/dx
      dv[2] = -cos(xx)*cos(yy); // du/dy
      dv[1] =  cos(xx)*cos(yy); // dv/dx
      dv[3] = -sin(xx)*sin(yy); // dv/dy

      return;
    }
};
*/


/*
class  KimMoinFlow
{
  // For Steady Navier-Stokes  
  
  private:
    double rho, mu;

  public:
    KimMoinFlow(double r, double m, double c1)
    {
      rho = r;
      mu  = m;
    }

    virtual ~KimMoinFlow(){}

    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      double cx = cos(xx);
      double sx = sin(xx);
      double cy = cos(yy);
      double sy = sin(yy);

      if(dir == 0)
        return  (-rho*sx*cx - 2.0*mu*cx*sy + sx*cx );
      else //if(dir == 1)
        return  (-rho*sy*cy + 2.0*mu*sx*cy + sy*cy );
    }

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      if(dir == 0)
        return ( -cos(xx)*sin(yy) );
      else if(dir == 1)
        return (  sin(xx)*cos(yy) );
      else //if(dir == 2)
        return  ( -0.25*(cos(2.0*xx)+cos(2.0*yy)) );
    }

    virtual  void computeDerivatives(double xx, double yy, double tt, double* dv)
    {
      dv[0] =  sin(xx)*sin(yy); // du/dx
      dv[2] = -cos(xx)*cos(yy); // du/dy
      dv[1] =  cos(xx)*cos(yy); // dv/dx
      dv[3] = -sin(xx)*sin(yy); // dv/dy

      return;
    }
};
*/


class  KimMoinFlow
{
  // Unsteady Navier-Stokes

  private:
    double rho, mu;

  public:
    KimMoinFlow(double r, double m, double c1)
    {
      rho = r;
      mu  = m;
    }

    virtual ~KimMoinFlow(){}

    virtual  double  computeForce(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      double cx = cos(xx);
      double sx = sin(xx);
      double cy = cos(yy);
      double sy = sin(yy);

      if(dir == 0)
        return  (-2.0*cx*sy*(rho*cos(2.0*tt)+mu*sin(2.0*tt)) - rho*sx*cx*sin(2.0*tt)*sin(2.0*tt) + sx*cx*sin(2.0*tt)*sin(2.0*tt) );
      else //if(dir == 1)
        return  ( 2.0*sx*cy*(rho*cos(2.0*tt)+mu*sin(2.0*tt)) - rho*sy*cy*sin(2.0*tt)*sin(2.0*tt) + sy*cy*sin(2.0*tt)*sin(2.0*tt) );
    }

    virtual  double  computeValue(int dir, double xx=0.0, double yy=0.0, double tt=0.0)
    {
      if(dir == 0)
        return ( -cos(xx)*sin(yy)*sin(2.0*tt) );
      else if(dir == 1)
        return (  sin(xx)*cos(yy)*sin(2.0*tt) );
      else //if(dir == 2)
        return  ( -0.25*(cos(2.0*xx)+cos(2.0*yy))*sin(2.0*tt)*sin(2.0*tt) );
    }

    virtual  void computeDerivatives(double xx, double yy, double tt, double* dv)
    {
      double c1 = sin(2.0*tt);
      
      dv[0] =  sin(xx)*sin(yy)*c1; // du/dx
      dv[2] = -cos(xx)*cos(yy)*c1; // du/dy
      dv[1] =  cos(xx)*cos(yy)*c1; // dv/dx
      dv[3] = -sin(xx)*sin(yy)*c1; // dv/dy

      return;
    }
};



#endif


