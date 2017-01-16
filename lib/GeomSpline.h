
#ifndef incl_GeomSpline_h
#define incl_GeomSpline_h


#include "MathVector.h"
#include "List.h"
#include "GeomObject.h"
#include "GeomPoint.h"
#include "Domain.h"


class GeomSpline: public GeomObject
{ 
  public:

    GeomSpline(void);
    virtual ~GeomSpline();
    
    Vector<void*> suppPnt;

    virtual bool isType(GeomTypeEnum type) { return (type == SPLINE); }

    double xmin(void);

    double area(double*);

    void discretise(Domain*, int);

    void checkElemSizes(void);

    void draw(int);

  private:

    double tang1[2], tang2[2];

    void calcTangent(int, double*);

    void calcPoint(int, double, double*);

    void subdivide(int, ListInfinite<VectorFixed<double,5> >&,
                   double, double, double*, double*, double, double);
};

#endif




