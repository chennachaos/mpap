/*=============================================================================
        File: NurbsBASE.h
  Created by: Chennakesava Kadapa          (08 Jan 2011)
 Purpose    : Header file for the definitions of NURBS CURVE Class
 
 ============================================================================*/
#ifndef NurbsBASE_H
#define NurbsBASE_H


#include <iostream>
#include <math.h>
#include "NURBS_1D.h"
#include "PropertyItem.h"
#include "ComputerTime.h"


using namespace std;


class  NurbsBASE
{
  public:

    int  nelem, nlbf, ndof, ngbf, nGP, nsize, DIM, tis;

    double  rhoInfty;

    bool  PHYSICS_TYPE;

    PropertyItem  ElemProp, MatlProp;

    VectorArray<bool>  closed, tracflag;

    VectorArray<double> Uinit, disp, dispBCdata;
    
    vector<double>  gausspoints, gaussweights;

    VectorArray<int> edgedata, intfdata, gbfnums, toUfull;

    ListArray<VectorArray<double> > dispBCs, forceBCs;

    ListArray<VectorArray<double> >  Values, ValuesDot, ValuesDotDot;
    ListArray<VectorArray<double> >  ValuesDotPrev, ValuesDotDotPrev;
    ListArray<VectorArray<double> >  ValuesPrev, ValuesPrev2, ValuesPrev3, ValuesPrev4;
    ListArray<VectorArray<double> >  ValuesCur, ValuesDotCur, ValuesDotDotCur;

    ListArray<VectorArray<int> > INC, IEN, ID, LM;

    VectorXd  td;



    NurbsBASE() { };

    virtual ~NurbsBASE();

    virtual void initializeBCdata()
    { cout << "   'initializeBCdata' is not defined for this NURBS Object !\n\n"; }

    virtual int GenerateConnectivityArrays1(int& )
    { cout << "   'GenerateConnectivityArrays1' is not defined for this NURBS Object !\n\n"; return 0;}

    int GenerateConnectivityArrays2();

    void  printConnectivityArrays();

    void  setTimeIntegrationParameters(int ttt, double rho1);

    void  setTimeParam();

    void  updateIterStep();

    void  setSolidOrFluid(int ttt);

    virtual void addInitDOFvalues()
    { cout << "   'addInitDOFvalues' is not defined for this NURBS Object !\n\n"; }

    virtual void updateCoordinates(double*)
    { cout << "   'updateCoordinates' is not defined for this NURBS Object !\n\n"; }

    virtual void PlotElements(int, bool, int*)
    { cout << "   'PlotElements' is not defined for this NURBS Object !\n\n"; }

    virtual void PlotControlPoints(int)
    { cout << "   'PlotControlPoints' is not defined for this NURBS Object !\n\n"; }

    virtual void printGeomToFile()
    { cout << "   'printGeomToFile' is not defined for this NURBS Object !\n\n"; }

    virtual void updateCoordsSingleCP(int num, int dir, double val)
    { cout << "   'updateCoordsSingleCP' is not defined for this NURBS Object !\n\n"; }

    virtual void computeNET()
    { cout << "   'computeNET' is not defined for this NURBS Object !\n\n"; }

    virtual void updateValues(int, double*)
    { cout << "   'updateValues' is not defined for this NURBS Object !\n\n"; }

    virtual void geomToVector(double*)
    { cout << "   'geomToVector' is not defined for this NURBS Object !\n\n"; }

    virtual void resetGeometry(double*)
    { cout << "   'resetGeometry' is not defined for this NURBS Object !\n\n"; }

    virtual void writeToFile(MyString&, int)
    { cout << "   'writeToFile' is not defined for this NURBS Object !\n\n"; }

    virtual void createAndWriteBasisFunctions(int)
    { cout << "   'createAndWriteBasisFunctions' is not defined for this NURBS Object !\n\n"; }

    virtual void createAndWriteBasisFunctionsDerivatives(int, int)
    { cout << "   'createAndWriteBasisFunctionsDerivatives' is not defined for this NURBS Object !\n\n"; }


};



#endif  //NurbsBASE_H
