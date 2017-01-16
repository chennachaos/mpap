
#ifndef incl_GeomSurface_h
#define incl_GeomSurface_h


#include "MathVector.h"
#include "MathMatrix.h"
#include "GeomObject.h"
#include "GeomSpline.h"

using namespace MatricesWulf;


enum OrientationEnum { FORWARD, BACKWARD };



class GeomSurface: public GeomObject
{ 
  public:

    GeomSurface(void);
    virtual ~GeomSurface();
    
    Vector<void*> suppSpln;

    Vector<int>  orientation, loop;

    bool hasBeenDiscretised;

    virtual bool isType(GeomTypeEnum type) { return (type == SURFACE); }

    bool initialise(void);

    void generateMesh(MatrixFullArray<double>&, MatrixFullArray<int>&, 
                      void*, int&, int&, int&,
                      Domain*, int, bool showFlg = false);

    void draw(int);

  private:

};

#endif




