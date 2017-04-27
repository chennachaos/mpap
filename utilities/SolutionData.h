
#ifndef incl_SolutionData_h
#define incl_SolutionData_h

#include "headersBasic.h"
#include "headersEigen.h"
#include "PropertyItem.h"


class  SolutionData
{
  private:
    int  tis, size1, size2, size3, size4;
    double  rho, dt;

  public:

    bool firstIter, STAGGERED;

    PhysicsType  PHYSICS_TYPE;

    MatrixXd  solnCFM;

    VectorXd  var1, var1Prev, var1Cur, var1Dot, var1DotPrev, var1DotCur, var1DotDot, var1DotDotPrev, var1DotDotCur;
    VectorXd  var2, var2Prev, var2Cur, var2Dot, var2DotPrev, var2DotCur, var2DotDot, var2DotDotPrev, var2DotDotCur;
    VectorXd  var3, var3Prev, var3Cur, var4, var4Prev, var4Cur;
    VectorXd  dDot, dDotPrev, var1PrevIter;
    VectorXd  var1Init, var2Init, vorticity, var1applied, var1Prev2, var1Prev3, var1Prev4, var1Extrap;
    
    vector<double>  FluidProps, stagParams;

    VectorXd  td, rhsVec, forcePred, forceCur, reac;
    VectorXd  force, forcePrev, forcePrev2, forcePrev3, forcePrev4, forcePrev5, forcePrev6;
    VectorXd  soln, solnInit, forceTemp, forceTempPrev;
    VectorXd  bodyForce, bodyForceCur, bodyForcePrev, bodyForcePrev2, bodyForcePred, bodyForcePredPrev, bodyForcePredPrev2;

    //PropertyItem  ElemProp, MatlProp;
    List<PropertyItem>  ElemProp, MatlProp;

    vector<int>  node_map_new_to_old;
    vector<int>  node_map_old_to_new;


    SolutionData();

    ~SolutionData(){}

    void SetTimeIncrementType(int ttt)
    {  tis = ttt; }

    void SetRho(double ttt)
    {  rho = ttt; }

    void  SetStaggeredParams(vector<double>&  vecTemp);

    void  SetPhysicsTypetoSolid()
    {
      PHYSICS_TYPE = PHYSICS_TYPE_SOLID;
      return;
    }

    void  SetPhysicsTypetoFluid()
    { 
      PHYSICS_TYPE = PHYSICS_TYPE_FLUID;
      return;
    }

    void  printSelf();

    void  build();

    void  initialise(int size1=0, int size2=0, int size3=0, int size4=0);

    void  setTimeParam();

    void  timeUpdate();

    void  updateIterStep();

    void  reset();

    void  interpolateForce();

    void  perform_Aitken_accelerator_force();

    void  perform_Aitken_accelerator_displacement();

};





#endif


