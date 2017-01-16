
#ifndef incl_BeamBending_h
#define incl_BeamBending_h


#include "List.h"
#include "Domain.h"
#include "MathVector.h"
#include "MathMatrix.h"
#include "Plot.h"


extern Plot plot;







class BeamBending: public Domain
{ 
  public:

    BeamBending(void);
    virtual ~BeamBending();
      
    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInteractions(void);

    virtual void doForBending(bool, bool, bool, bool, bool, bool, bool, bool, int);
    
    virtual void findMinMaxXBB(bool, bool, bool, bool, bool, bool,  bool,
                               double *xmn, double *xmx);
  private:

    int    nSect, nIntC;
 
    double Iyy, Izz, Iyz,
           minl, refHdW, h0[10], dh[10], mn[8], mx[8];

    Vector<double> x, val, pz, py, length;
    Vector<int> typ;

    VectorArray<double> brhs, A;
    MatrixFullArray<double> bmtx;
 
    void checkAndCompleteBC(void);

    void analysis(void);

    void plotAndPrint(int, int, bool, bool, double, double);

    void plotSystem(int, int, double, double);

    void addW0(int, int, double, double);
                    
    void addW1(int, int, double, double);
                    
    void addV0(int, int, double, double);
                    
    void addV1(int, int, double, double);

    void addBz0(int, int, double, double);
                     
    void addBz1(int, int, double, double);
                     
    void addBy0(int, int, double, double);
                     
    void addBy1(int, int, double, double);
                     
    void addMy0(int, int, double, double);
                     
    void addMy1(int, int, double, double);
                     
    void addMz0(int, int, double, double);
                     
    void addMz1(int, int, double, double);
                     
    void addSz0(int, int, double, double);
                     
    void addSz1(int, int, double, double);
                     
    void addSy0(int, int, double, double);
                     
    void addSy1(int, int, double, double);

    double w(double, int);

    double v(double, int);

    double bz(double, int);

    double by(double, int);

    double My(double, int);

    double Mz(double, int);

    double Sz(double, int);

    double Sy(double, int);

    void  drawW(double, double *);

    void  drawV(double, double *);

    void drawBz(double, double *);

    void drawBy(double, double *);

    void drawMy(double, double *);

    void drawMz(double, double *);

    void drawSz(double, double *);

    void drawSy(double, double *);

    void drawCoorSys(double, double*, int);

    void drawClamp(double, double*);

    void drawSimpleSupp(double, double*);

    void drawMomHinge(double, double*);

    void drawShearHinge(double, double*);

    void drawPointLoad(double, double*, double);
};



#include "DomainInlineFunctions.h"

define_reference_cast(beamBending,BeamBending)

define_isType(isBeamBending,BEAMBENDING)
	











#endif




