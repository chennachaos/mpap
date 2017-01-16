/*=============================================================================
        File: NurbsSOLID.h
  Created by: Chennakesava Kadapa          (09 Jan 2011)
 Purpose    : Header file for the definitions of NURBS SOLID Class


 ============================================================================*/

#ifndef NurbsSOLID_H
#define NurbsSOLID_H

#include "NurbsUtilitiesSURFACE.h"

#include "PropertyItem.h"
using namespace std;


// SOLID Class
class NurbsSOLID : public NurbsBASE
{
  public:
    NET3D  PP;
    CNET3D  Pw;
    KNOTVECTOR  U, V, W;
    DEGREE  p, q, r;


    int   ngbf1, ngbf2, ngbf3, nelem1, nelem2, nelem3, nGP1, nGP2, nGP3, nlbf1, nlbf2, nlbf3, ngbf1m2, nlbf1m2, nelem1m2, nGP1m2;

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights1, gaussweights2, gaussweights3, face_rot_angle;
    
    ListArray<NurbsSURFACE>   boundary_faces;
    
    bool face_rotation_flag;


    NurbsSOLID() {};

    NurbsSOLID(CNET3D& Pw1, KNOTVECTOR& U1, KNOTVECTOR& V1, KNOTVECTOR& W1, DEGREE p1, DEGREE q1, DEGREE r1);

    NurbsSOLID(const NurbsSOLID&);

    virtual ~NurbsSOLID();

    CPOINT SolidPoint(double, double, double);

    void CurveOnSurf(int dir, double u, NurbsCURVE* curve1);

    NurbsSOLID& operator = (const NurbsSOLID &rhs);

    virtual void initializeBCdata();

    virtual int GenerateConnectivityArrays1(int& );

    virtual void addInitDOFvalues();

    virtual void updateCoordinates(double*);

    virtual void PlotElements(int, bool, int*);

    virtual void PlotControlPoints(int);

    virtual void geomToVector(double*);

    virtual void resetGeometry(double*);

//    virtual void printGeomToFile();

    virtual void updateCoordsSingleCP(int num, int dir, double val);

    virtual void computeNET();

    void PlotElementsVTK(int, bool, int*);

    void PlotValues(int);

    virtual void updateValues(int, double*);

    double computeValue(int, double, double, double);

    void ShapeFunctions(double, double, double, double*);
    
    double computeValueAndShanpeFns(int, double, double, double, double*);

    void ShapeFunDerivatives(int*, double*, double*, double*, double*, double*, double&);

    void ShapeFunDerivatives3(int, int*, double*, double*, double*, double&);

    void deformationGradient(int*, bool, double*, double*, double*, double*, double&);
    
    void print2screen();
    
    int  gbfNumFace1(int, int);
    int  gbfNumFace2(int, int);
    int  gbfNumFace3(int, int);
    int  gbfNumFace4(int, int);
    int  gbfNumFace5(int, int);
    int  gbfNumFace6(int, int);

    void  computeGBFnumbers(int, int&);
    
    void  computeToUfull();
        void New();

};



#endif  //NurbsSOLID_H
