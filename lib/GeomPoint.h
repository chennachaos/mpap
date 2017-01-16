
#ifndef incl_GeomPoint_h
#define incl_GeomPoint_h


#include "GeomObject.h"
#include "MathVector.h"



class GeomPoint: public GeomObject
{ 
  public:

    GeomPoint(void);

    GeomPoint(double*);

    virtual ~GeomPoint();
    
    double x[3];

    double h;

    Vector<double> transferFact;

    Vector<int>    transferNode;

    virtual bool isType(GeomTypeEnum type) { return (type == POINT); }

    void draw(int);

  private:

};

#endif




