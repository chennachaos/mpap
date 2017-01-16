
#ifndef incl_FreeNodeInterpolation_h
#define incl_FreeNodeInterpolation_h


#include "FreeNode.h"
#include "ProfileListArray.h"
#include "MathVector.h"




class FreeNodeInterpolation: public FreeNode
{
  public:

    FreeNodeInterpolation(void);

    FreeNodeInterpolation(int, double*, int, int, VectorArray<int>*, VectorArray<int>*);

    virtual ~FreeNodeInterpolation();

    VectorArray<double> u, x, reac;

    int                 nd, typ, meshDepth, freeDepth;

    VectorArray<int>    *dofU, *dofX;

    ListArray< ListArray<DataDxdu> >   dxdu;
              
    ListArray< ListArray<DataToFree> > toFree;
              
    ListArray< ListArray<DataToMesh> > toMesh;
              
    ListArray<DataToBLay> toBLay;

    FreeNodeInterpolation &operator=(FreeNodeInterpolation &);

    void giveUX(double*, double*, double*, double*, double*, double*);
    void giveU (double*, double*, double*, double*, double*);
    void giveX (double*, double*, double*, double*, double*);

    void getReactions(double*, double*, double*, double*, double*);

    void print(void);      
};


#endif




