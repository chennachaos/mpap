
#ifndef incl_Mesh_h
#define incl_Mesh_h

#include "Geometry.h"
#include "MathMatrix.h"
#include "MathVector.h"
#include "Element.h"
#include "ElementGroup.h"
#include "MyStringList.h"
#include "DependentDoF.h"
#include "Solver.h"
#include "ObjectSurface.h"
#include "Definitions.h"
#include "ElementGeom2DEnum.h"
#include "NodeGeomLink.h"
#include "SelectNode.h"
#include "List.h"


using namespace MatricesWulf;


class PrescribedDisplacement: public ListItem
{
  public:
    int nd;
    double u[10];
};




class StructuredEdge2D: public ListItem
{
  public:
    int ix[3];
    double dat[5];
};


class StructuredQuad2D: public ListItem
{
  public:
    int ix[9];
    StructuredQuad2D* iq[4];
    StructuredEdge2D* ie[4];
    void addCentroidAndArea(double *, Vector<double> &);
    bool split(StructuredQuad2D *, int, int &, StructuredQuad2D *, StructuredEdge2D *,
               Vector<double> &, List<StructuredQuad2D> &,
               List<StructuredEdge2D> &, Vector<StructuredEdge2D*> &);
};








class Mesh: public Geometry
{
  public:

    Mesh(void);
    virtual ~Mesh();
      
    int        nen;
    int        numnp;
    int        numel;
    int        numElemGrp;
    int        nequ;
    int        neqx;

    int        lastSearchNode;
    
    double     tolMesh, ctimUpdateMesh;

    Vector<int> elemGrpProp, elemGrpToBeMeshed;
    
    List<PropertyItem> elemProp;

    MatrixFullArray<double> x0, x, xn, d, d0, v, vn, u, un, u3, u4, u5, u6, reac, reacMesh;

    VectorArray<double> r, frc, outp, wrndFact, elemSizeCurr, elemSizeOpt;

    VectorArray<int>    frcTmFct, frcNd, *nodeElem, *nodeNode, wrndType, wrndIndx;

    ListArray< VectorArray<int> >   wrndNode;

    List<DependentDoF> uDep, xDep;
    
    double rNormPrev, rNorm;
    
    MatrixFullArray<int>    idu, idx;

    List<ElementGroup>  elemGrp;

    Element **elem;

    int *bndNd, **bnd, nBndNd, nBnd;
	    
    Solver       *solverMesh;

    ObjectSurface *surf3D;

    SelectNode selectNode;

    VectorArray<bool> nodeFlag;

    MatrixFullArray<int> ixTmp;
    Vector<int>          ubcIntTmp, xbcIntTmp, elemGrpPropTmp;
    Vector<double>       ubcDblTmp, xbcDblTmp;

    void         addElemGrpProp(int);
    
    virtual bool isALE(bool flag = false);

    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInputData2(void) { return; }

    virtual void prepareInteractions(void);
    
    virtual void printInfo(void);
    
    virtual void findMinMaxX(double*, double*, bool defFlg = false);

    virtual void findMinMaxU(int, double&, double&);

    virtual void plotMesh(bool, bool);
    
    virtual void paintElemGrp(int, bool defFlg = false);

    virtual void plotNodes(int, int, bool defFlg = false);

    virtual void plotBoun(bool defFlg = false);
    
    virtual void plotFixed(bool defFlg = false);
    
    virtual void plotLoad(double, bool defFlg = false);
    
    virtual void plotReac(double, bool defFlg = false);
    
    virtual void plotElemNum(bool defFlg = false);

    virtual void prepareElemProp(int, char**);

    virtual void printNodalData(int);
   
    virtual void contourPlot(int, int, int, bool, double&, double&, bool);
   
    virtual void projectToNodes(MyString&, int, int);

    virtual void writeNodalData(void);
    
    virtual void plotGaussPoints(int, bool defFlg = false);

    virtual void plotU1D(int);

    virtual void updateUDepIncrements(void);
   
    virtual void reset(void);

    virtual int  updateMesh(int, bool printRes = true);
    
    virtual void elementDiffStiffTestMesh(double,int,int,int,bool);

    virtual void remeshElemGroups(bool showFlg = false);

    virtual void calcElemSizeOpt(double, double, double, bool);

    virtual void getInputElemSizeOpt(void);

    virtual double getElemSizeOptInternal(double*);

    virtual void calcElemSizeCurr(void);

    virtual void smoothElemSizeOptDist(int, int, double, double);

    virtual bool elemSizeRatioOK(double);

    virtual void transferNodalDataTest(double*);

    virtual void writeMeshToFile(MyString &, bool defm = false);

    virtual void interactiveNodeSelection(int, int, bool deselectFlag = false);

    virtual void simpleNodeSelection(void);

    inline  bool nodeHasDoFU (int i) { int j=0, n=i*ndf-ndf; 
                                       while (j<ndf) if (idu.x[n+j]>0) break; else j++; 
                                       if (j<ndf) return true; else return false; };

    inline  bool nodeHasDoFX (int i) { int j=0, n=i*ndm-ndm; 
                                       while (j<ndm) if (idx.x[n+j]>0) break; else j++; 
                                       if (j<ndm) return true; else return false; };

    virtual void printComputerTime(bool reset = true, int detailFlg = 1);

    bool    import(int, char*, int, int);

    bool    generateSimple2DMesh(int, int, int, int*, double*);

    bool    extrude(Domain*, double*, double*, int, int);
 
    void    setBoundaryConditions(bool*, bool remove = false);

    void    fixNodeMotion(bool*, bool remove = false);

    void    setFreeNodes(int);

    void    writeInputFile(char*, int, bool, bool, bool);

    void    setTmpPresDisp(int, double*, double*, double);

    int guessNearestNode(int, double *); // in Mesh3.cpp

  private:

    Vector<int> tmpFreeNode;

    List<PrescribedDisplacement> tmpPresDisp;

    ListInfinite<NodeGeomLink> ndGmLnk;

  // in Mesh2.cpp:

    void plotBoun_hlp(double*, double&, int&);
 
    bool isBoundaryNode2D(int, int&);
   
    void initialiseUDep(List<DependentDoF> *);
    
    void updateUDepInc_hlp(List<DependentDoF> &, MatrixFullArray<int> &);

    void replaceUDepTmFct(List<DependentDoF> *);

    void replaceFrcTmFct(void);

    bool sameFace(Vector<int> &, Vector<int> &);

    void separateBndNodes(Vector<int>&, Vector<int>&, VectorArray<int>&, 
                          VectorArray<unsigned int>&);
 
    void generateNodeNodeConnectivity(void);

    void generateNodeElemConnectivity(void);

    void generateBoundaryData2D(Vector<int>&, Vector<int>&);

    void generateBoundaryData3D(Vector<int>&, Vector<int>&);
 
  // in Mesh3.cpp:  stuff related to (re)meshing

    int  pointInElement(double*);

    int pointInElement_hlp(int, double*);

    void smoothElemSizeDist_hlp0(void);

    void smoothElemSizeDist_hlp1(int, double, double);

    void smoothElemSizeDist_hlp2(int, double, double);

    void smoothElemSizeDist_hlp3(double);

    void smoothElemSizeDist_hlp4(double);

    void reorganiseMesh(ListInfinite< MatrixFullArray<double> > &,
                        ListInfinite< MatrixFullArray<int> > &,
                        ListInfinite< ListInfinite<NodeGeomLink> > &,
                        Vector<int> &);

    void getTransferDataInternal(double *, Vector<int> &, Vector<double> &);

    void transfer_hlp(MatrixFullArray<double> &, List< Vector<int> > &, 
                      List< Vector<double> > &, int, int, int);

    void transferInternalVariables(Element**, int);

    void calcElemNodeDataAlongGeomObj(Vector<int> &, List< Vector<int> > &,
                                      List< Vector<int> > &, Vector<int> &, void *);

  // in Mesh4.cpp:  stuff related to preprocessing

    bool getNeighboursOnSameGeometry2D(int, Vector<int>&);

    void completeData(void);

  protected:  // in order to allow Deniz's microCells to generate elements

    Element *newElement(int);

    void prepareAndCheckUDep(List<DependentDoF> *); // in Mesh2.cpp, 
                                                    // also called from FiniteElementBVPWI.cpp
};



enum { XDASH, YDASH, ZDASH, SQUARE, TRIANGLE, CIRCLE };

#define CHECKNEQX true


#include "DomainInlineFunctions.h"

define_reference_cast(mesh,Mesh)

define_isType(isMesh,MESH)
	



#endif




