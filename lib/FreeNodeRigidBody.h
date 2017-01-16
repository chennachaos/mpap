
#ifndef incl_FreeNodeRigidBody_h
#define incl_FreeNodeRigidBody_h


#include "FreeNode.h"
#include "ProfileListArray.h"
#include "MathVector.h"




class FreeNodeRigidBody: public FreeNode
{
  public:

    FreeNodeRigidBody(void);

    FreeNodeRigidBody(int, double*, int, int, VectorArray<int>*, VectorArray<int>*);

    virtual ~FreeNodeRigidBody();

    VectorArray<double> u, x, reac;

    int                 nd, typ, meshDepth, freeDepth;

    VectorArray<int>    *dofU, *dofX;

    ListArray< ListArray<DataDxdu> >   dxdu;
              
    ListArray< ListArray<DataToFree> > toFree;
              
    ListArray< ListArray<DataToMesh> > toMesh;
              
    ListArray<DataToBLay> toBLay;

    FreeNodeRigidBody &operator=(FreeNodeRigidBody &);

    void giveUX(double*, double*, double*, double*, double*, double*);
    void giveU (double*, double*, double*, double*, double*);
    void giveX (double*, double*, double*, double*, double*);

    void getReactions(double*, double*, double*, double*, double*);

    void print(void);      
};


#endif




