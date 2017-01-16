
#ifndef incl_FreeNode_h
#define incl_FreeNode_h


#include "List.h"
#include "MathVector.h"



class DataToFree
{
  public:

    DataToFree(void) { invPtr = NULL; }

    virtual ~DataToFree() { }

    void setData(int i) { nd = i; } 

    void allocMem(int ndfi, int ndmi, int ndfj, int ndmj)

                        { dgfduf.setDim(ndfi*ndfj);
                          dgfdxf.setDim(ndmi*ndmj); } 

    void makeLagr(void) { if (dgfduf.n < dgfdxf.n) dgfduf.setDim(dgfdxf.n);
                          dgfdxf.free();
                          dgfdxf.x = dgfduf.x; }

    void zeroMtx(void)  { dgfduf.zero(); dgfdxf.zero(); }

    int nd; // smallest nd = 0

    DataToFree *invPtr;

    int pos;

    VectorArray<double> dgfduf, dgfdxf;
};






class DataToMesh
{
  public:

    DataToMesh(void) { }

    virtual ~DataToMesh() { }

    void setData(int i) { nd = i; }

    void allocMem(int ni, int nj) { xx.setDim(ni*nj); }

    void zeroMtx(void) { xx.zero(); }

    void swap(DataToMesh &b)

           { int m = nd; nd = b.nd; b.nd = m; VectorArray<double> a = xx; xx = b.xx; b.xx = a; }

    int nd; // smallest nd = 0
 
    VectorArray<double> xx;
};






class DataToBLay
{
  public:

    DataToBLay(void) { }

    virtual ~DataToBLay() { }

    void setData(int i) { nd = i; }

    void allocMem(int ndfEl, int ndfF, int ndmF, int ndm)
                        {
                          dgfdul.setDim(ndfF*ndfEl);
                          dgfdxl.setDim(ndfF*ndm);
                          dglduf.setDim(ndfEl*ndfF);
                          dgldxf.setDim(ndfEl*ndmF); }

    void makeLagr(void) {
                          if (dglduf.n < dgldxf.n) dglduf.setDim(dgldxf.n);
                          dgldxf.free();
                          dgldxf.x = dglduf.x; }

    void zeroMtx(void)  {
                          dgfdul.zero();
                          dgfdxl.zero();
                          dglduf.zero();
                          dgldxf.zero(); }

    int nd; // smallest nd = 1

    VectorArray<double> dgfdul, dgfdxl, dglduf, dgldxf;
};






class DataToBLayMesh
{
  public:

    DataToBLayMesh(void) { }

    virtual ~DataToBLayMesh() { }

    void setData(int i) { nd = i; }

    void allocMem(int ndm) { dmldxf.setDim(ndm*ndm); }

    void zeroMtx(void) { dmldxf.zero(); }

    int nd; // smallest nd = 1

    VectorArray<double> dmldxf;
};






class FreeNode: public ListItem
{
  public:

    FreeNode(void);

    FreeNode(int, double*, int, int, VectorArray<int>*, VectorArray<int>*);

    virtual ~FreeNode();

    VectorArray<double> u, x, reac, reacn;

    bool                isLagr;

    int                 typ, nd, meshDepth, freeDepth, *idu; // smallest nd = 1

    VectorArray<int>    *dofU, *dofX;

    ListArray< ListArray<DataToFree> > toFree;
              
    ListArray< ListArray<DataToMesh> > toMesh;
              
    ListArray<DataToBLay> toBLay;

    ListArray<DataToBLayMesh> toBLayMesh;

    FreeNode &operator=(FreeNode &);

    void makeLagrangian(void);

    void giveUX(double*, double*, double*, double*, double*, double*);
    void giveU (double*, double*, double*, double*, double*);
    void giveX (double*, double*, double*, double*, double*);

    void getReactions(double*, double*, double*, double*, double*);

    void print(void);

    void zeroDerivativeMtx(void);
};






template<class Type> void printDependencies(ListArray< ListArray<Type> > &d)
{
  int i, j;

  for (i=0; i<d.n; i++)
  {
    if (d[i].n < 1) cout << "{ ";
    else
    {
      cout << "{" << d[i][0].nd;
      for (j=1; j<d[i].n; j++)
        cout << "," << d[i][j].nd;
    }
    cout << "}\n";
  }
  return;
}






template<class Type> void generateDependencies(ListArray< ListArray<Type> > &d, 
                                                       List< Vector<int> > &tmp)
{
  int i, j;

  d.setDim(tmp.n);

  for (i=0; i<tmp.n; i++)
  {
    d[i].setDim(tmp[i].n);
    for (j=0; j<tmp[i].n; j++) d[i][j].setData(tmp[i][j]);
  }
  return;
}




enum { INTERPOLATION_FN, RIGIDBODY_FN, UNDEF_FN };


#endif




