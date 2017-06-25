
#include "SolutionData.h"
#include "ComputerTime.h"
#include "TimeFunction.h"
#include "MpapTime.h"
#include "Functions.h"
#include "PropertyItem.h"
#include "BasisFunctionsBSpline.h"
#include "util.h"


extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;



SolutionData::SolutionData()
{
  td.resize(100);
  
  STAGGERED = true;

  stagParams.resize(10);
  stagParams[0] = 1;
  stagParams[1] = 1;
  stagParams[2] = 1.0;

  size1 = size2 = size3 = size4 = 0;

  tis = 0;
  rho = 0.0;
  dt  = 1.0;

  return;
}




void SolutionData::SetStaggeredParams(vector<double>&  vecTemp)
{
  //cout << " SolutionData::SetStaggeredParams(vector<double>&  vecTemp) " << endl;
  //printVector(vecTemp);
  assert(vecTemp.size() >= 3);

  stagParams = vecTemp;
  
  STAGGERED = (stagParams[0] == 0);
  
  return;
}


void SolutionData::initialise(int s1, int s2, int s3, int s4)
{
  size1 = s1;
  size2 = s2;
  size3 = s3;
  size4 = s4;

    /////////////////////// 
    // for solid domain
    //
    // var1         --> disp
    // var1Dot      --> velocity
    // var1DotDot   --> acceleration
    // var2         --> pressure (for mixed formulations)
    //
    /////////////////////// 

    /////////////////////// 
    // for fluid domain
    //
    // var1     --> velocity
    // var1Dot  --> acceleration
    // var2     --> pressure
    // var3     --> Lagrange multipliers
    /////////////////////// 

    var1.resize(size1);
    var1.setZero();

    var1Prev    = var1;
    var1Prev2   = var1;
    var1Prev3   = var1;
    var1Prev4   = var1;
    var1Cur     = var1;
    
    var1Extrap     = var1;
    
    var1Dot     = var1;
    var1DotPrev = var1;
    var1DotCur  = var1;

    var1DotDot     = var1;
    var1DotDotPrev = var1;
    var1DotDotCur  = var1;

    var1Init    = var1;
    var1applied = var1;

    bodyForce         = var1;
    bodyForcePrev     = var1;
    bodyForcePred     = var1;
    bodyForcePredPrev = var1;
    bodyForceCur      = var1;

    force         = var1;
    forcePrev     = var1;
    forcePrev2    = var1;
    forcePrev3    = var1;
    forcePrev4    = var1;
    forcePrev5    = var1;
    
    forcePred     = var1;
    forceCur      = var1;
    forceTemp     = var1;
    forceTempPrev = var1;
    rhsVec        = var1;
    
    reac =  var1;
    
    dDot = var1;
    dDotPrev = var1;


    var2.resize(size2);
    var2.setZero();

    var2Prev = var2;
    var2Cur  = var2;

    var2Dot     = var2;
    var2DotPrev = var2;
    var2DotCur  = var2;

    var2DotDot     = var2;
    var2DotDotPrev = var2;
    var2DotDotCur  = var2;

    var3.resize(size3);
    var3.setZero();

    var3Prev = var3;
    var3Cur  = var3;

    var4.resize(size4);
    var4.setZero();

    var4Prev = var4;
    var4Cur  = var4;

  return;
}



void SolutionData::printSelf()
{
//   cout << " Degree and Jacobian " << degree[0] << '\t' << Jfull << endl;
   return;
}





void  SolutionData::setTimeParam()
{
  //cout << PHYSICS_TYPE << endl;
  //cout << tis << '\t' << rho << '\t' << mpapTime.dt << endl;

  if( PHYSICS_TYPE == PHYSICS_TYPE_SOLID )
    SetTimeParametersSolid(tis, rho, mpapTime.dt, td);
  else
    SetTimeParametersFluid(tis, rho, mpapTime.dt, td);
  
  //printVector(td);

   // solve for the initial profile

   return;
}



void  SolutionData::timeUpdate()
{
  // store the variables 

  var1Prev4 = var1Prev3;
  var1Prev3 = var1Prev2;
  var1Prev2 = var1Prev;

  var1Prev = var1;
  var2Prev = var2;
  var3Prev = var3;
  var4Prev = var4;
  
  var1Extrap = var1Prev;
  //var1Extrap = 2.0*var1Prev - var1Prev2;
  //var1Extrap = 3.0*var1Prev - 3.0*var1Prev2 + var1Prev3;

  var1DotPrev = var1Dot;
  var1DotDotPrev = var1DotDot;
  
  dDotPrev = dDot;

  var2DotPrev = var2Dot;
  var2DotDotPrev = var2DotDot;

  bodyForcePrev = bodyForce;

  // initialise the variables

  //if(!PHYSICS_TYPE)
  //{
    //var1.setZero();
    //var2.setZero();
    //var3.setZero();
    //var4.setZero();
  //}

  if(STAGGERED)
  {
    double  q1, q2, q3, q4, fact;

    int  predType = (int) stagParams[1];
    //cout << " predType = " << predType << endl;
    
    VectorXd  tempForce(force.rows());
    
    tempForce.setZero();

    //if( var1[0] < 1.0e-10 )
    //if( (var1[0] <= 1.0e-10) || (var1[0] >= 1.55) )
    //if( (var1[0] <= 1.0e-10) || (var1[0] > 0.5) )
    //if( (var1[1] > 1.0e-10) || (var1[2] > 1.0e-10) )
    //{
      //predType = 1;
      //var1DotPrev.setZero();
      //var1DotDotPrev.setZero();
    //}


    switch(predType)
    {
      case 1:
        q1 =  1.0;  q2 = 0.0;  q3 =  0.0;  q4 =  0.0;
        break;
        
      case 2:
        q1 =  2.0;  q2 = -1.0;  q3 =  0.0;  q4 =  0.0;
        //q1 =  4.0/3.0;  q2 = -1.0/3.0;  q3 =  0.0;  q4 =  0.0;
        break;
      
      case 3:
        q1 =  3.0;  q2 = -3.0;  q3 =  1.0;  q4 =  0.0;
        //q1 =  2.5;  q2 = -2.0;  q3 =  0.5;  q4 =  0.0;
        break;
      
      case 4:
        q1 =  4.0;  q2 = -6.0;  q3 =  4.0;  q4 = -1.0;
        //q1 =  104.0/35.0;  q2 = -114.0/35.0;  q3 =  56.0/35.0;  q4 = -11.0/35.0;
        break;
        
      default:
        cout << " unknown value for 'predType' " << endl;
        break;
    }

    forcePred = q1*force + q2*forcePrev + q3*forcePrev2 + q4*forcePrev3;

    forceCur  = td[2]*forcePred + (1.0-td[2])*force;

  }

  forcePrev5 = forcePrev4;
  forcePrev4 = forcePrev3;
  forcePrev3 = forcePrev2;
  forcePrev2 = forcePrev;
  forcePrev  = force;

  return;
}



void SolutionData::interpolateForce()
{
  //cout << " STAGGERED = " << STAGGERED << endl;
  if(STAGGERED)
  {
    //printVector(stagParams);
    double beta = stagParams[2];

    //if( var1[0] < 1.0e-10 )
    //if( (var1[0] <= 1.0e-10) || (var1[0] > 0.5) )
    //if( (var1[1] > 1.0e-10) || (var1[2] > 1.0e-10) )
      //beta = 1.0;

    //if( var1[1] > 1.0e-10 )
      //forceTemp[0] -= var1[1];
    //if( var1[2] > 1.0e-10 )
      //forceTemp[0] += var1[2];
    //else
    //{
    //}

    //printVector(forceTemp);

    //cout << " beta = " << beta << '\t' << var1[1] << endl;

    force     =  beta*forceTemp + (1.0-beta)*forcePred;

    forceCur  =  td[2]*force + (1.0-td[2])*forcePrev;
  }
  else
  {
    force     =  forceTemp;
    forceCur  =  td[2]*force + (1.0-td[2])*forcePrev;
  }
  //printf("\n\n");    printVector(forceCur);

  return;
}



void SolutionData::updateIterStep()
{
  char fct[] = "SolutionData::updateIterStep";
  
  if( PHYSICS_TYPE == PHYSICS_TYPE_SOLID )
  {
    //if(STAGGERED)
    //{
      // displacement as the primary variable
      //cout << "  SOLID SOLID SOLID " << endl;

      var1Dot     = td[10]*var1 + td[11]*var1Prev + td[12]*var1DotPrev + td[13]*var1DotDotPrev + td[14]*dDotPrev;
      var1DotDot  = td[15]*var1 + td[16]*var1Prev + td[17]*var1DotPrev + td[18]*var1DotDotPrev + td[19]*dDotPrev;
      dDot        = td[20]*var1 + td[21]*var1Prev + td[22]*var1DotPrev + td[23]*var1DotDotPrev + td[24]*dDotPrev;
    //}
    //else
    //{
      // velocity as the primary variable

      //var1       = td[40]*var1Dot + td[41]*var1Prev + td[42]*var1DotPrev + td[43]*var1DotDotPrev + td[44]*dDotPrev;
      //var1DotDot = td[45]*var1Dot + td[46]*var1Prev + td[47]*var1DotPrev + td[48]*var1DotDotPrev + td[49]*dDotPrev;
      //// ddot_{n+1} for modified state-space formulation
      //dDot       = td[50]*var1Dot + td[51]*var1Prev + td[52]*var1DotPrev + td[53]*var1DotDotPrev + td[54]*dDotPrev;
    //}

    // compute Current values

    var1Cur       = td[2]*var1       + (1.0-td[2])*var1Prev;
    var1DotCur    = td[2]*var1Dot    + (1.0-td[2])*var1DotPrev;
    var1DotDotCur = td[1]*var1DotDot + (1.0-td[1])*var1DotDotPrev;

    //printVector(var1Cur);
    //printVector(var1DotCur);
    //printVector(var1DotDotCur);

      //cout << " ggggggggggggg " << endl;
  }
  else
  {
    //cout << " uuuuuuuuuuuu  " << endl;
    //printVector(td);
    var1Dot    = td[9]*var1 + td[10]*var1Prev + td[11]*var1Prev2 + td[12]*var1Prev3 + td[13]*var1Prev4 + td[15]*var1DotPrev ;
    var2Dot    = td[9]*var2 + td[10]*var2Prev + td[15]*var2DotPrev ;

    var1Cur    = td[2]*var1    + (1.0-td[2])*var1Prev; // velocity
    var2Cur    = td[2]*var2    + (1.0-td[2])*var2Prev; // pressure
    var3Cur    = td[2]*var3    + (1.0-td[2])*var3Prev; // Lagrange parameters
    var4Cur    = td[2]*var4    + (1.0-td[2])*var4Prev; // Solid dof (for contact elements)

    var1DotCur = td[1]*var1Dot + (1.0-td[1])*var1DotPrev;
    var2DotCur = td[1]*var2Dot + (1.0-td[1])*var2DotPrev;
    //cout << " uuuuuuuuuuuu  " << endl;
  }

  return;
}


void  SolutionData::reset()
{
  var1 = var1Prev;
  var2 = var2Prev;
  var3 = var3Prev;
  var4 = var4Prev;

  var1Prev  = var1Prev2;
  var1Prev2 = var1Prev3;
  var1Prev3 = var1Prev4;

  var1Dot = var1DotPrev ;
  var1DotDot = var1DotDotPrev ;
  
  dDot = dDotPrev;

  var2Dot = var2DotPrev ;
  var2DotDot = var2DotDotPrev ;

  bodyForce = bodyForcePrev;
  bodyForcePrev = bodyForcePrev2;

  force = forcePrev;
  forcePrev = forcePrev2;
  forcePrev2 = forcePrev3;
  forcePrev3 = forcePrev4;
  forcePrev4 = forcePrev5;

  return;
}



void  SolutionData::perform_Aitken_accelerator_force()
{
  double  num, denom;

  for(int ii=0; ii<size1; ii++)
  {
    num = force[ii] - forcePrev[ii];
    num = num*num;

    denom = force[ii] - 2.0*forcePrev[ii] + forcePrev2[ii];

    force[ii] -= num/denom;
  }

  return;
}


void  SolutionData::perform_Aitken_accelerator_displacement()
{
  double  num, denom;

  for(int ii=0; ii<size1; ii++)
  {
    num = var1[ii] - var1Prev[ii];
    num = num*num;

    denom = var1[ii] - 2.0*var1Prev[ii] + var1Prev2[ii];

    force[ii] -= num/denom;
  }

  return;
}

    
    
    
    
    