
#ifndef incl_BasicFace_h
#define incl_BasicFace_h

#include "List.h"


class BasicFace: public ListItem
{ 
  public:

    BasicFace(void);
    virtual ~BasicFace();

    int ix[3];

    unsigned int edgeBits;

    float normal[3];

    void calcNormal(float *);

    bool facingOK(float*, float*);

    void simplePaint(int, int, int*, float*, bool edgesOnly = false);

    float calcZDepth(float*);

    bool frontContains(float*, float*);

    bool backContains(float*, float*);

    bool frontContains(int*, int*);

    bool backContains(int*, int*);

    void contourPlot(int, float*, float*, float, float, int*, float*, bool);

    //void backContourPlot(int, float*, float*, float, float, int*, float*);

    inline int n1(int n0)
      { int i=0; while (ix[i]!=n0) i++; if (i==2) return ix[0]; else return ix[i+1]; }

    inline int n2(int n0)
      { int i=0; while (ix[i]!=n0) i++; if (i==0) return ix[2]; else return ix[i-1]; }
};


#endif




