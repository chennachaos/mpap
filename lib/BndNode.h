
#ifndef incl_BndNode_h
#define incl_BndNode_h


#include "List.h"
#include "DependentDoF.h"
#include "FreeNode.h"
#include "MathVector.h"


class BndNode: public ListItem
{
  public:

    BndNode(void);

    BndNode(int, double*, Vector<void*>&, Vector<void*>&, Vector<double>&);

    virtual ~BndNode();

    double *x;

    int    nd, ndat, xFlag; // smallest nd = 1

    VectorArray<double> dat;

    ListArray<FreeNode*> freeU, freeX;

    ListArray<DependentDoF*> uDep, xDep;

    BndNode &operator=(BndNode &);

    void getUX(int, int, double*, double*, double*, double*);
    void getU (int, int, double*, double*, double*, double*);
    void getX (int, int, double*, double*, double*, double*);
 
    void giveReactions(int, int, double*, double*, double*, double*);

/*    void calcDerivatives();

    void getUXDeriv(int, int, double*, double*, double*, double*);
    void getUDeriv(int, int, double*, double*, double*, double*);
    void getXDeriv(int, int, double*, double*, double*, double*);*/

    void giveReactionsDeriv(int, int, double*, double*, double*, double*);

    bool checkOK(void);

    void print(int);

    void generateDeps(int, bool);

    int    giveDerivN(int);
    double *giveDeriv(int);

    int    getDerivN(int);
    double *getDeriv(int);
};


enum { NOX, XEQU, INDX };


#define MAX_DOF 10


#endif




