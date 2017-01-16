/*=============================================================================
        File: NurbsSURFACE.h
  Created by: Chennakesava Kadapa          (09 Jan 2011)
 Purpose    : Header file for the definitions of NURBS SURFACE Class


 ============================================================================*/

#ifndef NurbsSURFACE_H
#define NurbsSURFACE_H

#include "NurbsUtilitiesCURVE.h"
#include "PropertyItem.h"
#include <Eigen/Dense>
//#include <Eigen/Sparse>

using namespace Eigen;

using namespace std;


// SURFACE Class
class NurbsSURFACE : public NurbsBASE
{
  public:
    NET PP;
    CNET Pw;
    KNOTVECTOR U, V;
    DEGREE  p, q;


    int   ngbf1, ngbf2, nelem1, nelem2, nGP1, nGP2, nlbf1, nlbf2;

    bool   PLATE_BENDING;

    vector<double>   gausspoints1, gausspoints2, gaussweights1, gaussweights2;

    ListArray<VectorArray<double> >  dispdata;


    NurbsSURFACE() {};

    NurbsSURFACE(CNET& Pw1, KNOTVECTOR& U1, KNOTVECTOR& V1, DEGREE p1, DEGREE q1);

    NurbsSURFACE(const NurbsSURFACE& surf1);

    virtual ~NurbsSURFACE();

    CPOINT SurfacePoint(double, double);
    
    EPOINT SurfacePoint2(double, double);

    void CurveOnSurf(int dir, double u, NurbsCURVE* curve1);

    NurbsSURFACE& operator = (const NurbsSURFACE &rhs);

    void SurfDerPointHom(double u, double v, int k, ListArray<ListArray<CPOINT> >&);

    void SurfDerPointRat(double u, double v, int k, ListArray<ListArray<EPOINT> >&);

    void SurfacePointInverse(CPOINT& CP, double u0, double v0, VectorArray<double>& );


    virtual void initializeBCdata();

    virtual int GenerateConnectivityArrays1(int& );

    virtual void addInitDOFvalues();

    virtual void updateCoordinates(double*);

    virtual void PlotElements(int, bool, int*);

    virtual void PlotControlPoints(int);

//    virtual void printGeomToFile();

    virtual void updateCoordsSingleCP(int num, int dir, double val);

    virtual void computeNET();

    virtual void writeToFile(MyString&, int);

    void readSurfaceFromFile(MyString&);

    void PlotElementsVTK(int, bool, int*);

    void PlotValues(int);

    virtual void updateValues(int, double*);

    double computeValue(int, double, double);

    void  ShapeFunctions(double, double, double*);

    void  ShapeFunDerivatives(int*, double*, double*, double*, double*, double&);

    void  ShapeFunDerivatives2(int*, double*, double*, double*, double*, double*, double*, double*, double*, double&);

    void  ShapeFunDerivatives5(int, double*, double*);

    void deformationGradient(int, int, bool, double*, double*, double*, double&);

    void deformationGradient2(int*, double*, double* , double* , MatrixXd&, MatrixXd&, double&);

    virtual void geomToVector(double*);

    virtual void resetGeometry(double*);

    double  computeValueAndShanpeFns(int, double, double, double*);
    
    double   computeValue2(int, int, double*);

};



#endif  //NurbsSURFACE_H
