
#ifndef incl_Perspective_h
#define incl_Perspective_h


#include "MathVector.h"
#include "MathMatrix.h"
#include "ObjectSurface.h"


enum RotModeEnum { ROTZ, ROTX, ROTY, ROTFREE, ROTOBS };

enum ZoomModeEnum { VIEW, OBSERVER };

using namespace MatricesWulf;



class Perspective
{ 
  public:

    Perspective(void);
    virtual ~Perspective();
    
    float observer[3], view[3], xCntr[3];

    RotModeEnum rotMode;

    ZoomModeEnum zoomMode;

    ObjectSurface *changePersObjSurfPtr;

    bool  setDesEqAct;

    float fit(float*, float*);

    void  reset(void);
    
    void  pictureCoor(double *x, float *picCoor)
    {
      float xf[3] = { x[0], x[1], x[2] };
      pictureCoor(xf,picCoor);
      return;
    }

    void  pictureCoor(float*, float*);

    float zDepth(float*);

    void  calcNew(float, float, float*, int);

    void  plotChangePersObjSurf(bool backFaceFlag, bool defmFlag);

    void  print(float*);

    void  prepare(Vector<double> &, int);

    void  searchCone(double*, double*, float*, float);

    void  generateBoundingBox(double*, double*);

    void  addObjSurf(ObjectSurface*);

    void  removeAllObjSurf(void);

    void  setNewPersFlag(void);

  private:

    ObjectSurface boundingBox, **objSurfPtr;

    int   nObjSurf;

    float uX[3], uY[3], viewDotView, viewDotObserver;

    MatrixFullArray<float> M1, M2, M; // keep here to avoid repeated memory allocation

    void drawCoorAxes(void);
};




#endif




