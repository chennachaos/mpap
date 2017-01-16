/*=============================================================================
        File: NurbsEPOINT.h
  Created by: Chennakesava Kadapa          (10 Dec 2010)
 Purpose    : Header file for the definition of EPOINT Class

 ============================================================================*/
#ifndef NurbsEPOINT_H
#define NurbsEPOINT_H

#include <iostream>
#include <math.h>
#include "MathMatrix.h"
#include "List.h"
#include "MathVector.h"
#include "Domain.h"
#include "myConstants.h"
#include "util.h"


// Euclidean Point 
class EPOINT
{
    public:
    double  x, y, z;

    EPOINT():x(0.0), y(0.0), z(0.0) {}

    ~EPOINT() { } ;

    EPOINT(double xcoord, double ycoord, double zcoord) : x(xcoord),   y(ycoord),   z(zcoord) { }

    EPOINT(const EPOINT& p1):  x(p1.x),  y(p1.y), z(p1.z) { }

    double Norm()
    {   return sqrt(x * x + y * y + z * z);  }
    

    EPOINT& operator = (const EPOINT &rhs)
    {
      this->x = rhs.x;      this->y = rhs.y;      this->z = rhs.z;

      return *this;
    }

    EPOINT& operator = (const double value)
    {
      this->x = value;      this->y = value;      this->z = value;

      return *this;
    }

    friend const EPOINT operator +(const EPOINT& p1, const EPOINT& p2)
    {
      return EPOINT(p1.x + p2.x, p1.y + p2.y, p1.z + p2.z);
    }

    friend const EPOINT operator -(const EPOINT& p1, const EPOINT& p2)
    {
      return EPOINT(p1.x - p2.x, p1.y - p2.y, p1.z - p2.z);
    }

    // post-multiplication by a constant
    friend const EPOINT operator *(const EPOINT& p, double kk)
    {
      return EPOINT(p.x * kk, p.y * kk, p.z * kk);
    }

    // pre-multiplication by a constant
    friend const EPOINT operator *(double kk, const EPOINT& p) 
    {
      return EPOINT(p.x * kk, p.y * kk, p.z * kk);
    }
    
    // divide all entires by a constant
    friend const EPOINT operator /(const EPOINT& p, double kk) 
    {
      return EPOINT(p.x/kk, p.y/kk, p.z/kk);
    }
    
    // check for equality of two CPOINTs
    friend bool operator ==(const EPOINT& p1, const EPOINT& p2)
    {
      return ( CompareDoubles(p1.x, p2.x) && CompareDoubles(p1.y, p2.y) && CompareDoubles(p1.z, p2.z) );
    }

    // negative of a CPOINT
    friend const EPOINT operator -(const EPOINT& p)
    {
      return EPOINT(-p.x, -p.y, -p.z);
    }

    void print2screen()
    {
        printf("%12.8f\t%12.8f\t%12.8f\n", x, y, z);
    }


};

#endif  //NurbsEPOINT_H












