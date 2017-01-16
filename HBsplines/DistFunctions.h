
#ifndef incl_DistanceFunction_h
#define incl_DistanceFunction_h


//#include "headersBasic.h"
//#include "util.h"
//#include "NurbsUtilitiesSOLID.h"

#include "myConstants.h"
#include <iostream>


using std::cout;
using std::endl;


//using namespace Eigen;


class DistanceFunction
{
  private:
    double dist;
  
  public:
    DistanceFunction() {}
    
    virtual ~DistanceFunction() {}

    virtual void computePoint(double param, double* coord)
    { cout << "   'computePoint(double param, double* coord)' is not defined for this distance function !\n\n"; return ; }

    virtual void computeNormal(double param, double* normal)
    { cout << "   'computeNormal(double param, double* normal)' is not defined for this distance function !\n\n"; return ; }

    virtual double ComputeDistance(double xx=0.0, double yy=0.0, double zz=0.0)
    { cout << "   'ComputeDistance(double xx=0.0, double yy=0.0, double zz=0.0)' is not defined for this distance function !\n\n"; return -9999.0; }

    virtual double ComputeDistance(myPoint& pt)
    { cout << "   'ComputeDistance(myPoint& pt)' is not defined for this distance function !\n\n"; return -9999.0; }
};


class  Circle: public DistanceFunction
{
   private:
     myPoint  center;
     double  radius;
   
   public:
     
     Circle()
     {
       center(0) = 0.0;
       center(1) = 0.0;
       radius    = 1.0;
     }

     Circle(double x1, double y1, double r)
     {
       center(0) = x1;
       center(1) = y1;
       radius    = r;
     }

     virtual ~Circle(){}

     virtual void computePoint(double param, double* coord)
     {
        // param is the sector angle starting from +ve x-axis
        coord[0] = center(0) + radius * cos(2.0*PI*param);
        coord[1] = center(1) + radius * sin(2.0*PI*param);
     }

     bool checkPointLocation(double xx, double yy)
     {
        double temp1 = xx - center(0);
        double temp2 = yy - center(1);
        
        return (sqrt(temp1 * temp1 + temp2 * temp2) < radius);
     }

     bool checkPointLocation(myPoint& pt)
     {
        myPoint  pp = pt - center;
        
        return (pp.norm() < radius);
     }

     virtual double ComputeDistance(myPoint& pt)
     {
        myPoint  pp = pt - center;
        
        return  (pp.norm() - radius);
     }

     virtual double ComputeDistance(double xx, double yy)
     {
        double temp1 = xx - center(0);
        double temp2 = yy - center(1);
        
        return (sqrt(temp1 * temp1 + temp2 * temp2) - radius);
     }

     void computeNormal(double param, double* normal)
     {
        // param is the sector angle starting from +ve x-axis
        normal[0] = cos(2.0*PI*param);
        normal[1] = sin(2.0*PI*param);
     }
};


class  Ellipse: public DistanceFunction
{
   private:     
     myPoint  center;
     double   radiusX, radiusY;
   
   public:

     Ellipse()
     {
       center(0) = 0.0;
       center(1) = 0.0;
       radiusX   = 1.0;
       radiusY   = 2.0;
     }

     Ellipse(double x1, double y1, double r1, double r2)
     {
       center(0) = x1;
       center(1) = y1;
       radiusX   = r1;
       radiusY   = r2;
     }

     virtual ~Ellipse(){}

     bool checkPointLocation(myPoint& pt)
     {
        myPoint  pp = pt - center;

        pp(0) /= radiusX;
        pp(1) /= radiusY;
        
        return  (pp.norm() < 1.0);
     }

     virtual double ComputeDistance(myPoint& pt)
     {
        myPoint  pp = pt - center;

        pp(0) /= radiusX;
        pp(1) /= radiusY;
        
        return  (pp.norm() - 1.0);
     }

     bool checkPointLocation(double xx, double yy)
     {
        double temp1 = (xx-center(0))/radiusX;
        double temp2 = (yy-center(1))/radiusY;
        
        return (sqrt(temp1*temp1 + temp2 * temp2) < 1.0);
     }
};




class  Rectangle: public DistanceFunction
{
   private:
     myPoint  pt1, pt2;
   
   public:

     Rectangle()
     {
       pt1(0) = 0.0;
       pt1(1) = 0.0;
       pt2(0) = 1.0;
       pt2(1) = 1.0;
     }

     Rectangle(double x1, double y1, double x2, double y2)
     {
       pt1(0) = x1;
       pt1(1) = y1;
       pt2(0) = x2;
       pt2(1) = y2;
     }

     Rectangle(myPoint& p1, myPoint& p2)
     {
       pt1 = p1;
       pt2 = p2;
     }
     
    virtual ~Rectangle(){}

     bool checkPointLocation(myPoint& pt)
     {
        return ( min(min(pt(0)-pt1(0), pt2(0)-pt(0)), min(pt(1)-pt1(1), pt2(1)-pt(1))) < 0.0);
     }

     virtual double ComputeDistance(myPoint& pt)
     {
       return min(min(pt(0)-pt1(0), pt2(0)-pt(0)), min(pt(1)-pt1(1), pt2(1)-pt(1)));
     }

     bool checkPointLocation(double xx, double yy)
     {
        return ( ComputeDistance(xx, yy) > 0.0 );
     }

     double ComputeDistance(double xx, double yy)
     {
       return min(min(xx-pt1(0), pt2(0)-xx), min(yy-pt1(1), pt2(1)-yy));
     }
};






#endif

