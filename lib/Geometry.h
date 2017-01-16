
#ifndef incl_Geometry_h
#define incl_Geometry_h


#include "Domain.h"
#include "GeomPoint.h"
#include "GeomSpline.h"
#include "GeomSurface.h"
//#include "GeomVolume.h"
#include "List.h"
#include "MathVector.h"




class HPointSource: public ListItem
{
  public:

    HPointSource(void) { }

    HPointSource(double *xx, double hh, double mxGrd, int ndm = 2)
    { 
      for (int i=0; i<ndm; i++) x[i] = xx[i];

      h = hh;
      maxGrad = mxGrd;
      return;
    }
    
    double x[3], h, maxGrad;

};









class GeomBndBase: public ListItem
{
  public:

    GeomBndBase(void) { geomObj = NULL; return; }
    
    void *geomObj;

    GeomSurface* whichSurface(void) { if (((GeomObject*)geomObj)->isType(SURFACE))
                                       return (GeomSurface*)geomObj; else return NULL; }

    GeomSpline* whichSpline(void) { if (((GeomObject*)geomObj)->isType(SPLINE))
                                     return (GeomSpline*)geomObj; else return NULL; }

    GeomPoint* whichPoint(void) { if (((GeomObject*)geomObj)->isType(POINT))
                                   return (GeomPoint*)geomObj; else return NULL; }
};










class GeomBndCond: public GeomBndBase
{
  public:

    Vector<int> idu, idx;
};



class GeomPresDisp: public GeomBndBase
{
  public:

    Vector<int>    tmFct;
    Vector<double> uBase;
};



class GeomDistLoad: public GeomBndBase
{
  public:

    Vector<int>    tmFct;
    Vector<double> fdl;
};












class GeomNodalDataOutput: public ListItem
{
  public:

    int    type, indx;
    double fact;

    Vector<void*> geomObj;
};











class Geometry: public Domain
{ 
  public:

    Geometry(void);
    virtual ~Geometry();

    List<GeomPoint> point;

    List<GeomSpline> spline;
    
    List<GeomSurface> surface;

    //List<GeomVolume> volume;

    List<HPointSource> hPntSrc;

    List<GeomBndCond> bndCond;

    List<GeomPresDisp> presDisp;

    List<GeomDistLoad> distLoad;

    List<GeomNodalDataOutput> wrndOutp;

    Vector<int> grpType;

    bool geometryDiscretisedAtLeastOnce;
 
    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void plotGeometry(unsigned int);

    virtual void findGeomMinMaxX(double*, double*);

    virtual void whichSurfaces(Vector<int>&, Vector<int>&);

    virtual void whichSplines(Vector<int>&, Vector<int>&);
    
    virtual void printComputerTime(bool reset = true, int detailFlg = 1);

  private:


};




#include "DomainInlineFunctions.h"

define_reference_cast(geometry,Geometry)

define_isType(isGeometry,GEOMETRY)
	



#endif




