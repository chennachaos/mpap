
#ifndef incl_Element_h
#define incl_Element_h


#include "List.h"
#include "Domain.h"
#include "MyString.h"
#include "DomainTypeEnum.h" // for derived elements
#include "MathVector.h"
#include "WhatToSolveForEnum.h"
#include "MathVector.h"


using namespace std;


class Element: public ListItem
{
  public:

    Element(void);

    virtual ~Element();
	  
    void   *elemGrp;
	 
    int    *ix, *sparse[2];

    double *intVar1, *intVar2;

    VectorArray<int> GpDatId;
    
    ListArray< VectorArray<int> > sMaster[2];
    
    ListArray< VectorArray<double> > sFact[2];
    
    bool flag;   // multipurpose flag, e.g. for searching
    int  iflg;   // multipurpose integer, e.g. for ordering and searching
    
    void   print(void);
    
    int   &idu(int,int);
    int   &idx(int,int);
    
    double &x   (int,int);
    double &xn  (int,int);
    double &x0  (int,int);
    double &u   (int,int);
    double &un  (int,int);
    double &u3  (int,int);
    double &u4  (int,int);
    double &u5  (int,int);
    double &u6  (int,int);
    double &d   (int,int);
    double &d0  (int,int);
    double &v   (int,int);
    double &vn  (int,int);
    double &outp(int,int j = 0);

    void   belongsTo(void*, void*);
    
    virtual int ndf(void) { return 0; }
    virtual int ndm(void) { return 0; }
    virtual int nen(void) { return 0; }
    
    virtual int ndfm(int) { return ndf(); }
    
    virtual int nivGP(void) { return 0; }
    virtual int nGaussPoints(void) { return 0; }
   
    virtual void diffStiffTest(double,int,int,bool);
    
    virtual bool forDomainType(int) { return false; }

    virtual double getElemSizeOpt(double*, double*);
    
    virtual void adjustAssembly(int, int);
 
    // declare all member functions of all derived elements
   
    virtual void assembleReactionForces(void);
 
    virtual int  calcStiffnessAndResidual(void)
      {	cout << "   'calcStiffnessAndResidual' is not defined for this element!\n\n";
	return 0; }
    
    virtual int  calcStiffnessAndResidualMesh(void)
      {	cout << "   'calcStiffnessAndResidualMesh' is not defined for this element!\n\n";
	return 0; }
    
    virtual int  calcMeshDerivatives(void)
      {	cout << "   'calcMeshDerivatives' is not defined for this element!\n\n";
	return 0; }
    
    virtual void plotOutline(bool defFlg = false)
      {	cout << "   'plotOutline' is not defined for this element!\n\n"; return; }
    
    virtual void paint(bool defFlg = false)
      {	cout << "   'paint' is not defined for this element!\n\n"; return; }
    
    virtual void initialiseIntVar(void)
      { cout << "   'initialiseIntVar' is not defined for this element!\n\n"; return; }

    virtual void putLabel(char*, bool defFlg = false)
      { cout << "   'putLabel' is not defined for this element!\n\n"; return; }

    virtual void contourPlot(int, int, int, double, double, bool defFlg = true)
      { cout << "   'contourPlot' is not defined for this element!\n\n"; return; }

    virtual void projectIntVar(int)
      { cout << "   'projectIntVar' is not defined for this element!\n\n"; return; }

    virtual void projectStress(int)
      { cout << "   'projectStress' is not defined for this element!\n\n"; return; }

    virtual void projectVorticity(int)
      { cout << "   'projectVorticity' is not defined for this element!\n\n"; return; }

    virtual void projectError(int)
      { cout << "   'projectError' is not defined for this element!\n\n"; return; }

    virtual double volume(bool init = false)
      { cout << "   'volume' is not defined for this element!\n\n"; return 0.; }

    virtual void givePlotSequence2D(Vector<int> &)
      { cout << "   'givePlotSequence2D' is not defined for this element!\n\n"; return; }

    virtual void giveFace3D(int, Vector<int> &)
      { cout << "   'giveFace3D' is not defined for this element!\n\n"; return; }

    virtual int nFaces(void)
      { cout << "   'nFaces' is not defined for this element!\n\n"; return 0; }

    virtual int nBasicFacesPerFace(void)
      { cout << "   'nBasicFacesPerFace' is not defined for this element!\n\n"; return 0; }

    virtual void defineBasicFace(int, int, int*, unsigned int*)
      { cout << "   'defineBasicFace' is not defined for this element!\n\n"; return; }

    virtual void plotGaussPoints(int,bool defFlg = false)
      { cout << "  'plotGaussPoints' is not available for this element!\n\n"; return; }
   
    virtual bool isClosed(void) { return true; }
    
    virtual void setGaussPointDataId(void)
      { cout << "  'setGaussPointDataId' is not available for this element!\n\n"; return; }

    virtual void diffAleTest(double,int,int,bool)
      { cout << "  'diffAleTest' is not available for this element!\n\n"; return; }

    virtual bool containsPoint(double*, double*)
      { cout << "  'containsPoint' is not available for this element!\n\n"; return false; }

    virtual double diameter(bool init = false)
      { cout << "  'diameter' is not available for this element!\n\n"; return -1.; }

    virtual void getDistLoadFact(Vector<double> &, Vector<int> &)
      { cout << "  'getDistLoadFact' is not available for this element!\n\n"; return; }

    virtual int finiteStrain(void)
      { cout << "  'finiteStrain' is not available for this element!\n\n"; return 0; }

    virtual void diffMeshDerivTest(double,int,int,bool)
      { cout << "  'diffMeshDerivTest' is not available for this element!\n\n"; return; }
    
    virtual void projectGradient(int, int)
      { cout << "  'projectGradient' is not available for this element!\n\n"; return; }
    
    virtual void projectNormOfGradient(int)
      { cout << "  'projectNormOfGradient' is not available for this element!\n\n"; return; }

    virtual void projectNormOfGradientSquared(int)
      { cout << "  'projectNormOfGradientSquared' is not available for this element!\n\n"; return; }

    virtual void getGradient(int,double*)
      { cout << "  'getGradient' is not available for this element!\n\n"; return; }
   
private:

};

#endif

