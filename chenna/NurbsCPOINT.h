/*=============================================================================
        File: NurbsCPOINT.h
  Created by: Chennakesava Kadapa          (21 Dec 2010)
 Purpose    : Header(interface) file for the definition of CPOINT class
 
============================================================================*/
#ifndef NurbsCPOINT_H
#define NurbsCPOINT_H


#include <iostream>
#include "NurbsEPOINT.h"


// Control Point
class CPOINT
{
  public:
    double   x, y, z, w;
 
    CPOINT() : x(0.0), y(0.0),  z(0.0), w(0.0) {}

    CPOINT(double xcoord, double ycoord, double zcoord, double wt) : x(xcoord), y(ycoord),  z(zcoord),  w(wt) { }

    CPOINT(const CPOINT& cp1) : x(cp1.x), y(cp1.y),  z(cp1.z),  w(cp1.w) { }

    ~CPOINT() { }


    // calculated the corresponding Euclidean Point
    EPOINT CalcEuclid()
    {   return EPOINT(x/w, y/w, z/w);  }

    // assignment operator
    CPOINT& operator = (const CPOINT &rhs)
    {
       this->x = rhs.x;   this->y = rhs.y;  this->z = rhs.z;  this->w = rhs.w;

       return *this;
    }

    CPOINT& operator = (const double value)
    {
       this->x = value;       this->y = value;       this->z = value;       this->w = value;

       return *this;
    }

    // operator +
    friend const CPOINT operator +(const CPOINT& cp1, const CPOINT& cp2)
    {
       return  CPOINT(cp1.x + cp2.x, cp1.y + cp2.y, cp1.z + cp2.z, cp1.w + cp2.w);
    }

    // operator -
    friend const CPOINT operator -(const CPOINT& cp1, const CPOINT& cp2)
    {
       return  CPOINT(cp1.x - cp2.x, cp1.y - cp2.y, cp1.z - cp2.z, cp1.w - cp2.w);
    }
    
    // post-multiplication by a constant
    friend const CPOINT operator *(const CPOINT& cp, double kk)
    {
       return CPOINT(cp.x * kk, cp.y * kk, cp.z * kk, cp.w * kk);
    }

    // pre-multiplication by a constant
    friend const CPOINT operator *(double kk, const CPOINT& cp)
    {
      return CPOINT(cp.x * kk, cp.y * kk, cp.z * kk, cp.w * kk);
    }

    // divide all entires by a constant
    friend const CPOINT operator /(const CPOINT& p, double kk)
    {
      return CPOINT(p.x/kk, p.y/kk, p.z/kk, p.w/kk);
    }


    // check for equality of two CPOINTs
    friend bool operator ==(const CPOINT& cp1, const CPOINT& cp2)
    {
      return (CompareDoubles(cp1.x,cp2.x) && CompareDoubles(cp1.y,cp2.y) && CompareDoubles(cp1.z,cp2.z) && CompareDoubles(cp1.w,cp2.w) );
    }
    
    // negative of a CPOINT
    friend const CPOINT operator -(const CPOINT& cp)
    {
      return CPOINT(-cp.x, -cp.y, -cp.z, -cp.w);
    }

    void print2screen()
    {
        printf("%12.8f\t%12.8f\t%12.8f\t%12.8f\n", x, y, z, w);
    }


};



#endif // NurbsCPOINT_H
