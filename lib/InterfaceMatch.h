
#ifndef incl_InterfaceMatch_h
#define incl_InterfaceMatch_h

#include "MathVector.h"
#include "FiniteElementBVPWNI.h"
#include "SolverSparse.h"
#include "FreeNode.h"


using namespace MatricesWulf;


struct MoreData: public ListItem  // needed for temporary node matching data structure
{
  MoreData(int cc, FreeNode &f) { c = cc; x = f.x; fnd = false; chk = false; }
  Vector<MoreData*> depPtr;
  VectorArray<double> x;
  bool fnd, chk;
  int c;
};







class InterfaceMatch: public Domain
{ 
  public:

    InterfaceMatch(void);

    virtual ~InterfaceMatch();

    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInteractions(void);

    virtual void setTimeParam(void);

    virtual void printInfo(void);

    virtual bool converged(void);

    virtual bool diverging(double);

    virtual void timeUpdate(void);

    virtual void updateIterStep(void);

  protected:

    VectorArray<FiniteElementBVPWNI*> domPtr;

    VectorArray<int> domType, domId, domMode;

    VectorArray<double> domReacTime;

    ListArray< VectorArray<int> > intfNodeToFreeNode, freeNodeToIntfNode;

    int numnp, nequ, neqx, *F2I, uType, initGuess, nFN;

    FreeNode *FN;

    double rNorm, rNormPrev, dxdu;

    MatrixFullArray<int> idu, idx;

    VectorArray<bool> isLagr;

    VectorArray<double> r;

    MatrixFullArray<double> x, xn, x0, u, un, du, dun, ddu, ddun;

    ListArray<FreeNode> &freeNode(int dm) { return domPtr[dm]->freeNode; }

    ListArray<BndNode> &bndNode(int dm) { return domPtr[dm]->bndNode; }

  private:

    // data that becomes obsolete after initialisation stage

    VectorArray<int> bcTmp;

    List<MoreData> node;

    Vector<int> changed, searchFrom;

    int compCount;

    bool show;

    void getCoorAndMatchNodes(void);

    void generateIdux(void);

    void setLagrangianNodes(void);

    void initKinData(void);

    void findNodeMatches(int);

    void findNodeForMe(int, int);

    void addNonMatchingNodes(int);

    void extendConnectivity(int);
    
};



#include "DomainInlineFunctions.h"

define_reference_cast(interfaceMatch,InterfaceMatch)

define_isType(isInterfaceMatch,INTERFACEMATCH)



enum { RESOLVE, ELIMINATE, ELIMINATE_MESH, RETURN_FORCES, RETURN_DISPLACEMENTS, DUMMY };

enum { VELOCITY, DISPLACEMENT };

enum { KEEP_VELOCITY, KEEP_DISPLACEMENT };





	
#endif


