
#ifndef incl_FiniteElementBVPWI_h
#define incl_FiniteElementBVPWI_h


#include "FiniteElementBVP.h"
#include "List.h"
#include "FreeNode.h"
#include "BndNode.h"



struct SomeData
{
  SomeData(int n) { bl.setDim(n); f.setDim(n); }

  VectorArray<DataToBLay*> bl; // pointers to freeNode[].toBLay

  VectorArray<int> f;
};



class FiniteElementBVPWI: public FiniteElementBVP
{
  public:

    FiniteElementBVPWI(void);

    virtual ~FiniteElementBVPWI();

    ListArray<FreeNode> freeNode;

    ListArray<BndNode> bndNode;

    List< VectorArray<int> > dofX, dofU;

    ListArray< ListArray<DataToMesh> > nodeToMesh;

    bool isConnected;

    double dxdu;

    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInteractions(void);

    virtual void setLagrangianFreeNodes(VectorArray<int> &);

    virtual void plotInterfaceNodes(int, int, int, bool);

    virtual void calculateDerivatives1(void);

    virtual void calculateDerivatives2(void);

    virtual void calculateDerivatives3(void);

    virtual void calculateDerivatives4(void);

    virtual void calculateDerivatives5(void);

    virtual void calculateDerivatives6(void);

    virtual void calculateDerivatives7(void);

    virtual void calculateDerivatives8(void);

    virtual void eliminate(bool, bool);

    virtual void eliminateDiffTest(double,double,double,int,int,bool);

    virtual void prepareForExternalSolver(void*,double*,bool);

    virtual void printComputerTime(bool reset = true, int detailFlg = 1);

  private:

    // permanent data

    int maxFreeDepth, maxMeshDepth, nAllDoFU, nAllDoFX,
        nLayElem1, nLayElem2, nLayElem3;

    bool fixedU, fixedX;

    ListArray< VectorArray<int> > layElem, ixbx, ixbu, ixlx, ixlu, ixlm;

    ListArray< VectorArray<DataToFree*> >     pffu, pffx;
    ListArray< VectorArray<DataToBLay*> >     pflu, pflx, plfx;
    ListArray< VectorArray<DataToBLayMesh*> > plfm;

    VectorArray<DataToMesh*> pmm;

    VectorArray<double> col;

    VectorArray<int> profU, profX;

    VectorArray<bool> nodeHasU, nodeHasX;

    VectorArray<DataToFree*>  f2f;

    VectorArray<SomeData*> n2bl;

    // data relevant only during initialisation phase, memory is then released

    bool show, membraneFlag;

    Vector<int> nodeFlagChanged;

    List< Vector<int> > intpTmp;

  protected:

    Vector<bool> freeNodeIsConnected;

    List<FreeNode> freeNodeTmp;

    List<BndNode>  bndNodeTmp;

  private:

    float ctimDeriv1, ctimDeriv2, ctimDeriv3, ctimDeriv4,
          ctimDeriv5, ctimDeriv6, ctimDeriv7, ctimDeriv8,
          ctimBack3,  ctimBack5,  ctimBack8,  ctimElim, ctimPrepExt;

    bool  doneDeriv1, doneDeriv2, doneDeriv3, doneDeriv4,
          doneDeriv5, doneDeriv6, doneDeriv7, doneDeriv8,
          doneElim,   donePrepExt;

    // private member functions

    void generateExtNodeNode(ListArray< Vector<int> > &);

    void generateBaseTree(ListArray< List< Vector<int> > > &, ListArray< Vector<int> > &);

    void generateBaseTreeHelp(List< Vector<int> > &, ListArray< Vector<int> > &, int);

    void init1InterpolationData(void);

    void init2LayerData(ListArray< List< Vector<int> > > &);

    void init3AssemblyData(void);

    void init4finalise(void);

    void findBndNodesForInterpolations(VectorArray<bool> &, bool);

    int  guessNearestMarkedNode(int, double*);

    int  findOneNodeInIntpElem(int, Vector<int> &, Vector<double> &f);

    void connectNodeToIntf(int, Vector<int> &, Vector<double> &);

    int  findAnotherNodeInIntpElem(int, Vector<int> &, Vector<int> &, Vector<double> &);

    bool nodeInIntpElem(int, Vector<int> &, Vector<double> &);

    void calcF2Fu(int, int, BndNode &, BndNode &, int, int, int, double *);

    void calcL2Fu(int, BndNode &, int, int, int, int, double *);

    void calcF2Lu(int, BndNode &, int, int, int, int, double *);

    void calcF2Fx(int, int, BndNode &, BndNode &, int, int, int, double *);

    void calcL2Fx(int, BndNode &, int, int, int, int, double *);

    void calcF2Lx(int, BndNode &, int, int, int, double *);

    void calcL2Fm(int, BndNode &, int, int, int, double *);

    void calcM2M(int, int, int, int, double *);

    void generateSimpleStiffnessMatrixU(double **);

    void generateSimpleStiffnessMatrixX(double **);

    void generateSimpleStiffnessMatrixXX(double **);
};




#include "DomainInlineFunctions.h"

define_reference_cast(finiteElementBVPWI,FiniteElementBVPWI)

define_isType(isFiniteElementBVPWI,FINITEELEMENTBVPWI)
	


#endif




