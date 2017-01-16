
#ifndef incl_BeamSection_h
#define incl_BeamSection_h


#include "List.h"
#include "Domain.h"
#include "MathVector.h"
#include "MathMatrix.h"
#include "Plot.h"


extern Plot plot;








class SecondMomentOfArea
{
  public:

    double yy, zz, yz;

    void zero(void) { yy = 0.; zz = 0.; yz = 0.; }
 
    double principalDirection(void);

    void rotate(double, SecondMomentOfArea &);

    void shift(double, SecondMomentOfArea &);
};







class Node: public ListItem
{
  public:

    Node(void) { }

    Node(double y, double z, double boom) { coor[0] = y; coor[1] = z; B = boom; return; } 

    virtual ~Node() { return; }

    double coor[2], B, q[4];

    void calcCoorG(double *, double &);

    void calcI(double *, SecondMomentOfArea &);

    void rotateCoor(double *);

    void calcShearFlow(double*, SecondMomentOfArea &);

    void draw(double, int, int);
};









class Flange: public ListItem
{
  public:

    Flange(void) { }

    Flange(int nd1, int nd2, double thick) 
    { 
      n1 = nd1-1; n2 = nd2-1; t = thick; 
      for (int i=0; i<10; i++) { q1[i] = 0.; qc[i] = 0.; q2[i] = 0.; Q[i] = 0.; } 
      return;  
    }

    virtual ~Flange() { return; }

    int n1, n2, l1, l2;

    double l, t, p, q1[10], qc[10], q2[10], Q[10];

    void calcGeom(List<Node> &);

    void calcCoorG(List<Node> &, double *, double &);

    void calcI(List<Node> &, double *, SecondMomentOfArea &);

    void calcShearFlow(List<Node> &, double *, SecondMomentOfArea &);

    void qForSandT(double *, double);

    void tauForSandT(void);

    void sigVonMises(void);

    void sigForMandN(List<Node> &, double*, double, SecondMomentOfArea &, double*, double);

    void draw(List<Node> &, double, int);

    void plotAndPrint(List<Node> &, int, char *, double, double, double, int);

    void plotLabels(List<Node> &, int, double, double);
};









class BSLoop: public ListItem
{
  public:

    BSLoop(void) { }

    BSLoop(Vector<int> &tmp) { flng = tmp; return; }

    virtual ~BSLoop() { return; }

    Vector<int> flng;

    double area;

    void sort(List<Node> &, List<Flange> &);
};










class BeamSection: public Domain
{ 
  public:

    BeamSection(void);

    virtual ~BeamSection();
      
    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInteractions(void);

    virtual void doForSection(bool, bool, bool, bool, bool, bool, bool, bool, int, double, int);
    
    virtual void findMinMaxX(double *xmn, double *xmx, bool defFlg = false);
 
  private:

    int    numnp, numfl;
 
    double mxThick, mxLength, mxBoom, area, J, alpha, diameter, 
           coorA[2], coorC[2], coorS[2], 
           S[2], M[2], T, N;

    SecondMomentOfArea I, pI;
 
    List<Node>   node;

    List<Flange> flange;

    List<BSLoop> loop;

    List< Vector<int> > nodeFlange, nodeNode;

    void analysis(void);
 
    void calcShearFlow(int);

    void calcBSLoops(void);

    bool isBSLoop(Vector<int> &);

    bool isNewBSLoop(Vector<int> &);

    void plotAndPrint(int, char*, double, int, bool);

    void checkAndCompleteBC(void);
};



#include "DomainInlineFunctions.h"

define_reference_cast(beamSection,BeamSection)

define_isType(isBeamSection,BEAMSECTION)
	


#endif




