
#ifndef incl_DistanceFunction3D_h
#define incl_DistanceFunction3D_h


#include "headersBasic.h"

#include "DistFunctions.h"


using namespace Eigen;


class  dSphere: public DistanceFunction
{
   private:
     myPoint  center;
     double  radius;
   
   public:
     
     dSphere()
     {
       center(0) = 0.0;
       center(1) = 0.0;
       center(3) = 0.0;
       radius    = 1.0;
     }

     dSphere(double x1, double y1, double z1, double r)
     {
       center(0) = x1;
       center(1) = y1;
       center(2) = z1;
       radius    = r;
     }

     ~dSphere(){}
     
     void SetParameters(double x1, double y1, double z1, double r)
     {
       center(0) = x1;
       center(1) = y1;
       center(2) = z1;
       radius = r;
     }
     
     void SetCenter(double x1, double y1, double z1)
     {
       center(0) = x1;
       center(1) = y1;
       center(2) = z1;
     }

     void SetCenter(myPoint& pt)    {  center = pt;  }

     void SetRadius(double r)     {  radius = r;   }
     

     bool checkPointLocation(double xx, double yy, double zz)
     {
       double temp1 = xx - center(0);
       double temp2 = yy - center(1);
       double temp3 = zz - center(2);
        
        return (sqrt(temp1 * temp1 + temp2 * temp2 + temp3 * temp3) < radius);
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

};



class  dEllipsoid: public DistanceFunction
{
   private:     
     myPoint  center;
     double   radiusX, radiusY, radiusZ;
   
   public:

     dEllipsoid()
     {
       center(0) = 0.0;
       center(1) = 0.0;
       center(2) = 0.0;
       radiusX   = 1.0;
       radiusY   = 2.0;
       radiusZ   = 3.0;
     }

     dEllipsoid(double x1, double y1, double z1, double r1, double r2, double r3)
     {
       center(0) = x1;
       center(1) = y1;
       center(2) = z1;
       radiusX   = r1;
       radiusY   = r2;
       radiusZ   = r3;
     }

     ~dEllipsoid(){}

     void SetParameters(double x1, double y1, double z1, double r1, double r2, double r3)
     {
       center(0) = x1;
       center(1) = y1;
       center(2) = z1;
       radiusX   = r1;
       radiusY   = r2;
       radiusZ   = r3;
     }

     void SetCenter(double x1, double y1, double z1)
     {
       center(0) = x1;
       center(1) = y1;
       center(2) = z1;
     }

     void SetCenter(myPoint& pt)
     {
       center = pt;
     }

     void SetRadius(double r1, double r2, double r3)
     {
       radiusX = r1;
       radiusY = r2;
       radiusZ = r3;
     }

     bool checkPointLocation(myPoint& pt)
     {
        myPoint  pp = pt - center;

        pp(0) /= radiusX;
        pp(1) /= radiusY;
        pp(2) /= radiusZ;
        
        return  (pp.norm() < 1.0);
     }

     virtual double ComputeDistance(myPoint& pt)
     {
        myPoint  pp = pt - center;

        pp(0) /= radiusX;
        pp(1) /= radiusY;
        pp(2) /= radiusZ;
        
        return  (pp.norm() - 1.0);
     }

     bool checkPointLocation(double xx, double yy, double zz)
     {
        double temp1 = (xx-center(0))/radiusX;
        double temp2 = (yy-center(1))/radiusY;
        double temp3 = (zz-center(2))/radiusZ;
        
        return (sqrt(temp1*temp1 + temp2*temp2 + temp3*temp3) < 1.0);
     }
};





class  dCylinder: public DistanceFunction
{
   private:     
     myPoint  center;
     double   radiusX, radiusY, radiusZ;
   
   public:

     dCylinder()
     {
       center(0) = 0.0;
       center(1) = 0.0;
       center(2) = 0.0;
       radiusX   = 1.0;
       radiusY   = 2.0;
       radiusZ   = 3.0;
     }

     dCylinder(double x1, double y1, double z1, double r1, double r2, double r3)
     {
       center(0) = x1;
       center(1) = y1;
       center(2) = z1;
       radiusX   = r1;
       radiusY   = r2;
       radiusZ   = r3;
     }

     ~dCylinder(){}

     void SetParameters(double x1, double y1, double z1, double r1, double r2, double r3)
     {
       center(0) = x1;
       center(1) = y1;
       center(2) = z1;
       radiusX   = r1;
       radiusY   = r2;
       radiusZ   = r3;
     }

     void SetCenter(double x1, double y1, double z1)
     {
       center(0) = x1;
       center(1) = y1;
       center(2) = z1;
     }

     void SetCenter(myPoint& pt)
     {
       center = pt;
     }

     void SetRadius(double r1, double r2, double r3)
     {
       radiusX = r1;
       radiusY = r2;
       radiusZ = r3;
     }

     bool checkPointLocation(myPoint& pt)
     {
        myPoint  pp = pt - center;

        pp(0) /= radiusX;
        pp(1) /= radiusY;
        pp(2) /= radiusZ;
        
        return  (pp.norm() < 1.0);
     }

     virtual double ComputeDistance(myPoint& pt)
     {
        myPoint  pp = pt - center;

        pp(0) /= radiusX;
        pp(1) /= radiusY;
        pp(2) /= radiusZ;
        
        return  (pp.norm() - 1.0);
     }

     bool checkPointLocation(double xx, double yy, double zz)
     {
        double temp1 = (xx-center(0))/radiusX;
        double temp2 = (yy-center(1))/radiusY;
        double temp3 = (zz-center(2))/radiusZ;
        
        return (sqrt(temp1*temp1 + temp2*temp2 + temp3*temp3) < 1.0);
     }
};



class  dBox: public DistanceFunction
{
   private:
     myPoint  pt1, pt2;
   
   public:

     dBox()
     {
       pt1(0) = 0.0;
       pt1(1) = 0.0;
       pt1(2) = 0.0;
       pt2(0) = 1.0;
       pt2(1) = 2.0;
       pt2(2) = 3.0;
    }

     dBox(double x1, double y1, double z1, double x2, double y2, double z2)
     {
       pt1(0) = x1;
       pt1(1) = y1;
       pt1(2) = z1;
       pt2(0) = x2;
       pt2(1) = y2;
       pt2(2) = z2;
    }

     dBox(myPoint& p1, myPoint& p2)
     {
       pt1 = p1;
       pt2 = p2;
     }
     
     ~dBox(){}

     void SetParameters(double x1, double y1, double z1, double x2, double y2, double z2)
     {
       pt1(0) = x1;
       pt1(1) = y1;
       pt1(2) = z1;
       pt2(0) = x2;
       pt2(1) = y2;
       pt2(2) = z2;
     }

     bool checkPointLocation(myPoint& pt)
     {
        return (min(min(min(pt(0)-pt1(0), pt2(0)-pt(0)) , min(pt(1)-pt1(1), pt2(1)-pt(1))), min(pt(2)-pt1(2), pt2(2)-pt(2))) < 0.0);
     }

     virtual double ComputeDistance(myPoint& pt)
     {
       return (min(min(min(pt(0)-pt1(0), pt2(0)-pt(0)) , min(pt(1)-pt1(1), pt2(1)-pt(1))), min(pt(2)-pt1(2), pt2(2)-pt(2))));
     }

     virtual double ComputeDistance(double xx, double yy, double zz)
     {
       return (min(min(min(xx-pt1(0), pt2(0)-xx), min(yy-pt1(1), pt2(1)-yy)), min(zz-pt1(2), pt2(2)-zz)));
     }
};



#endif
