
#ifndef incl_InterfaceN_h
#define incl_InterfaceN_h


#include "InterfaceMatch.h"


class InterfaceN: public InterfaceMatch
{ 
  public:

    InterfaceN(void);

    virtual ~InterfaceN();

    Solver *solver;

    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInteractions(void);

    virtual void printInfo(void);

    virtual void setSolver(int, int *parm = NULL, bool cIO = false);

    virtual int  calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true);

    virtual int  factoriseSolveAndUpdate(void);

    virtual void globalDiffStiffTest(double,int,int,bool);

    virtual void plotSolverMatrixPattern(char*);

    virtual void printComputerTime(bool reset = true, int detailFlg = 1);

  private:

    int neqTot;

    double ctimCalcSub1, ctimCalcSub2, ctimCalcSub3, ctimCalcFreeSurf,
           ctimFactSolvUpdt, ctimCalcStiffRes;

    bool   doneCalcSub1, doneCalcSub2, doneCalcSub3, doneCalcFreeSurf,
           symFlag;

    ListArray< VectorArray<int> > posSwap;

    ListArray< ListArray< VectorArray<int> > > colDblPos;

    VectorArray<int> nOff, posMtxSub2, posMtxSub3, posMtxSub4;

    List< VectorArray<int> > posMtx;

    void prepareMatrixPattern(void);

    void prepareMatrixPatternSub1(int, MatrixSparse<double> &, int *, int *, int *, int *,
                                       ListArray< Vector<int> > &,
                                       ListArray< Vector<int> > &);

    void prepareMatrixPatternSub2(int, MatrixSparse<double> &, int *, int *, int *, int *);

    void prepareMatrixPatternSub3(int, MatrixSparse<double> &, int *, int *, int *, int *);

    void prepareMatrixPatternFreeSurface(int, MatrixSparse<double> &, int *, int *, int *, int *);

    void calcStiffnessAndResidualSub1(int, double *);

    void calcStiffnessAndResidualSub2(int, double *);

    void calcStiffnessAndResidualSub3(int, double *);

    void calcStiffnessAndResidualFreeSurface(int, double *);
};



#include "DomainInlineFunctions.h"

define_reference_cast(interfaceN,InterfaceN)

define_isType(isInterfaceN,INTERFACEN)

	
#endif


