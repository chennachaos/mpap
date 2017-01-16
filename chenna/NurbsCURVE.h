/*=============================================================================
        File: NurbsCURVE.h
  Created by: Chennakesava Kadapa          (08 Jan 2011)
 Purpose    : Header file for the definitions of NURBS CURVE Class

 ============================================================================*/
#ifndef NurbsCURVE_H
#define NurbsCURVE_H


#include <iostream>
#include <math.h>
#include "NURBS_1D.h"
#include "PropertyItem.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "NurbsBASE.h"


using namespace std;


// Curve Class
class NurbsCURVE: public NurbsBASE
{
  public:
    ListArray<EPOINT>  PP;
    CPOLYGON Pw;
    KNOTVECTOR U;
    DEGREE  p;



    NurbsCURVE() ;

    virtual ~NurbsCURVE();

    NurbsCURVE(CPOLYGON& Pw1, KNOTVECTOR& U1, DEGREE p1);

    NurbsCURVE(const NurbsCURVE& curv1);

    NurbsCURVE& operator = (const NurbsCURVE &rhs);

    CPOINT CurvePoint(double u);

    void CurveDerPointHom(double u, int d, ListArray<CPOINT>& Cder);

    void CurveDerPointRat(double u, int d, ListArray<EPOINT>& Cder);

    double CurvePointInverse(EPOINT EP, double u0);


    virtual void initializeBCdata();

    virtual int GenerateConnectivityArrays1(int& );

    virtual void addInitDOFvalues();

    virtual void updateCoordinates(double*);

    virtual void PlotElements(int, bool, int*);

    virtual void PlotControlPoints(int);

    virtual void printGeomToFile();

    virtual void updateCoordsSingleCP(int num, int dir, double val);

    virtual void computeNET();

    //virtual int GenerateConnectivityArrays2();

    void PlotCurve(int);

    void ShapeFunctions(double, double*);

    void ShapeFunDerivatives(int, double, double*, double&);

    void ShapeFunsAndDerivatives(int, double, double*, double*, double&);

    void ShapeFunsAndDerivatives2(int, double, double*, double*, double*, double&);

    void deformationGradient(int, bool, double*, double&);

    virtual void createAndWriteBasisFunctions(int);

    virtual void createAndWriteBasisFunctionsDerivatives(int, int);

    virtual void geomToVector(double*);

    virtual void resetGeometry(double*);

    double computeValue(int, double);

};



#endif  //NurbsCURVE_H
