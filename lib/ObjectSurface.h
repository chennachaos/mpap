
#ifndef incl_ObjectSurface_h
#define incl_ObjectSurface_h


#include "MathMatrix.h"
#include "MathVector.h"
#include "BasicFace.h"
#include "List.h"
#include "Domain.h"





class SurfaceEdge
{
  public:

    SurfaceEdge(void) { return; }
    ~SurfaceEdge()    { return; }

    int j1, j2, f1, f2;

    void setVisibility(BasicFace**, bool);

    bool isKink(BasicFace**, float);
};





class ObjectSurface: public ListItem
{ 
  public:

    ObjectSurface(void);
    ObjectSurface(Domain*);
    virtual ~ObjectSurface();

    void free(void);

    bool updtFlagX, updtFlagU, newPersFlag, showBackFaces, showDeformed, carefulCheckDone;

    void addElementFace(int, int, Vector<int> &, VectorArray<int> &);

    void finalise(Vector<int> &);

    void generateBoundingBox(double*, double*);

    void draw(int, int, bool edgesOnly = false, bool backFaceSwap = false, bool defmFlag = false);

    bool getNeighboursOnSameGeometry(int, Vector<int>&, float);

    void plotNodePoint(int, float, int nb = -1);

    bool nodeIsReallyVisible(int);

    void updateU(int, int);

    void contourPlot(int);

    void writeInterfaceInterpolations(VectorArray<int> &, ofstream &);

  private:

    Domain *dom;

    float *x, *picx, *u, umn, umx, *view, *observer;

    int numfc, nPaint, numnp, numed, *scrx, *faceIsVisible, *perm;

    char *nodeIsVisible;   //  2 -> probably visible
                           //  1 -> visible
                           //  0 -> not visible

    Vector<int>   *bnd2face;

    VectorInfinite<int>   paintSequence;

    VectorInfinite<float> dblSortValue;

    List<BasicFace> faceList;

    BasicFace   **face;

    SurfaceEdge *edge;

    void updateVisible(bool backFaceSwap = false, bool defmFlag = false);

    void updateCoor(void);

    void carefulCheck(void);

    void setVisEdgeBits(int);
};


#endif




