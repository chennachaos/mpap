
#include "ImmersedRigidSolid.h"
#include "BasisFunctionsLagrange.h"
#include "myConstants.h"
#include "ImmersedIntegrationElement.h"
#include "FunctionsProgram.h"
#include "MpapTime.h"
#include "TimeFunction.h"


extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;


ImmersedRigidSolid::ImmersedRigidSolid()
{
}



ImmersedRigidSolid::ImmersedRigidSolid(int dd)
{
  DIM = dd;

  ndofRigidbody = 3*(DIM-1);

  matM.resize(ndofRigidbody, ndofRigidbody);
  matM.setZero();
  matC = matM;
  matK = matM;
  Klocal = matM;
  Flocal.resize(ndofRigidbody);

  dofData.resize(ndofRigidbody,-1);

  SolnData.initialise(3, 0, 0, 0);
  SolnData.setPhysicsTypetoSolid();

  preLoad.resize(ndofRigidbody, 0.0);
  initForcePred.resize(ndofRigidbody, 0.0);

  PRESC_MOTION = false;
  //PRESC_MOTION = true;

  PrescMotionTimeFuncs.resize(ndofRigidbody, -1);
}


ImmersedRigidSolid::~ImmersedRigidSolid()
{
}


void ImmersedRigidSolid::printSelf()
{
  printf("\n\n\n Rigid body %5d details \n\n", id);
  printf("--------------------------------------------------------------");
  printf("\n Mass matrix \n\n");
  printMatrix(matM);
  printf("\n Damping matrix \n\n");
  printMatrix(matC);
  printf("\n Stiffness matrix \n\n");
  printMatrix(matK);

  //if( std::any_of(PrescMotionTimeFuncs.begin(), PrescMotionTimeFuncs.end(), [](int i){return (i!=-1);}) )
  if(PRESC_MOTION)
  {
    printf(" Moves with a prescribed motion. Coupled FSI is ignored ... \n ");

    if(DIM==2)
    {
      if(PrescMotionTimeFuncs[0] == -1)
        printf("\t Is fixed in X-direction \n");
      else
        printf("\t Moves in X-direction as per time function %5d \n", PrescMotionTimeFuncs[0]);

      if(PrescMotionTimeFuncs[1] == -1)
        printf("\t Is fixed in Y-direction \n");
      else
        printf("\t Moves in Y-direction as per time function %5d \n", PrescMotionTimeFuncs[1]);

      if(PrescMotionTimeFuncs[2] == -1)
        printf("\t Is fixed in theta-direction \n");
      else
        printf("\t Moves in theta-direction as per time function %5d \n", PrescMotionTimeFuncs[2]);
    }
    if(DIM==3)
    {
      if(PrescMotionTimeFuncs[0] == -1)
        printf("\t Is fixed in X-direction \n");
      else
        printf("\t Moves in X-direction as per time function %5d \n", PrescMotionTimeFuncs[0]);

      if(PrescMotionTimeFuncs[1] == -1)
        printf("\t Is fixed in Y-direction \n");
      else
        printf("\t Moves in Y-direction as per time function %5d \n", PrescMotionTimeFuncs[1]);

      if(PrescMotionTimeFuncs[2] == -1)
        printf("\t Is fixed in Z-direction \n");
      else
        printf("\t Moves in Z-direction as per time function %5d \n", PrescMotionTimeFuncs[2]);

      if(PrescMotionTimeFuncs[3] == -1)
        printf("\t Is fixed in theta-direction \n");
      else
        printf("\t Moves in theta-direction as per time function %5d \n", PrescMotionTimeFuncs[3]);

      if(PrescMotionTimeFuncs[4] == -1)
        printf("\t Is fixed in phi-direction \n");
      else
        printf("\t Moves in phi-direction as per time function %5d \n", PrescMotionTimeFuncs[4]);

      if(PrescMotionTimeFuncs[5] == -1)
        printf("\t Is fixed in psi-direction \n");
      else
        printf("\t Moves in psi-direction as per time function %5d \n", PrescMotionTimeFuncs[5]);
    }
  }
  else
  {
    printf(" It has %5d free degree of freedom \n", totalDOF);
    printVector(dofData);

    if(DIM==2)
    {
      if(dofData[0] == -1)
        printf("\t Is fixed in X-direction \n");
      else
        printf("\t FSI in X-direction \n", dofData[0]);

      if(dofData[1] == -1)
        printf("\t Is fixed in Y-direction \n");
      else
        printf("\t FSI in Y-direction \n", dofData[1]);

      if(dofData[2] == -1)
        printf("\t Is fixed in theta-direction \n");
      else
        printf("\t FSI in theta-direction \n", dofData[2]);
    }
    if(DIM==3)
    {
      if(dofData[0] == -1)
        printf("\t Is fixed in X-direction \n");
      else
        printf("\t FSI in X-direction \n");

      if(dofData[1] == -1)
        printf("\t Is fixed in Y-direction \n");
      else
        printf("\t FSI in Y-direction \n");

      if(dofData[2] == -1)
        printf("\t Is fixed in Z-direction \n");
      else
        printf("\t FSI in Z-direction \n");

      if(dofData[3] == -1)
        printf("\t Is fixed in theta-direction \n");
      else
        printf("\t FSI in theta-direction \n");

      if(dofData[4] == -1)
        printf("\t Is fixed in phi-direction \n");
      else
        printf("\t FSI in phi-direction \n");

      if(dofData[5] == -1)
        printf("\t Is fixed in psi-direction \n");
      else
        printf("\t FSI in psi-direction \n");
    }
  }

  return;
}


void ImmersedRigidSolid::setMass(vector<double>& tempVec)
{
  assert( tempVec.size() >= ndofRigidbody*ndofRigidbody );
  
  int ii, jj, kk=0;

  for(ii=0;ii<ndofRigidbody;ii++)
  {
    for(jj=0;jj<ndofRigidbody;jj++)
      matM(ii, jj) = tempVec[kk++];
  }
  
  return;
}


void ImmersedRigidSolid::setDamping(vector<double>& tempVec)
{
  assert( tempVec.size() >= ndofRigidbody*ndofRigidbody );

  int ii, jj, kk=0;

  for(ii=0;ii<ndofRigidbody;ii++)
  {
    for(jj=0;jj<ndofRigidbody;jj++)
      matC(ii, jj) = tempVec[kk++];
  }
  
  return;
}


void ImmersedRigidSolid::setStiffness(vector<double>& tempVec)
{
  assert( tempVec.size() >= ndofRigidbody*ndofRigidbody );

  int ii, jj, kk=0;

  for(ii=0;ii<ndofRigidbody;ii++)
  {
    for(jj=0;jj<ndofRigidbody;jj++)
      matK(ii, jj) = tempVec[kk++];
  }
  
  return;
}


void ImmersedRigidSolid::setBoundaryConditions(vector<int>& vectemp)
{
  dofData = vectemp;

  //cout << " rigid body boundary conditions " << endl;  printVector(dofData);
  
  return;
}



void ImmersedRigidSolid::setPrescribedMotion(vector<int>& vectemp)
{
  //cout << "  vectemp.size() = " << vectemp.size() << endl;
  if(vectemp.size() > 0)
  {
    for(int ii=0;ii<vectemp.size()-1;ii++)
      PrescMotionTimeFuncs.push_back(vectemp[ii+1]);
  
    PRESC_MOTION = true;
  }
  else
    PRESC_MOTION = false;

  return;
}



void ImmersedRigidSolid::setPreload(vector<double>& tempVec)
{
  assert( tempVec.size() >= ndofRigidbody );

  preLoad = tempVec;

  return;
}

void  ImmersedRigidSolid::setInitialForcePredictor(vector<double>& tempVec)
{
  assert( tempVec.size() >= ndofRigidbody );

  initForcePred = tempVec;

  return;
}


void  ImmersedRigidSolid::setRigidBodyMotionLimits(vector<vector<double> >& tempVec)
{
  rigidBodyMotionLimits = tempVec;

  return;
}


void  ImmersedRigidSolid::setNodalPositions(vector<vector<double> >&  datatemp)
{
  nNode = datatemp.size();
  //cout << " nNode  = " << nNode << endl;

  GeomData.NodePosOrig.resize(nNode);
  GeomData.NodePosCur.resize(nNode);
  GeomData.NodePosNew.resize(nNode);  
  GeomData.NodePosPrev.resize(nNode);
  GeomData.NodePosPrevCur.resize(nNode);  
  GeomData.specValNew.resize(nNode);
  GeomData.specValCur.resize(nNode);

  GeomData.acceNew.resize(nNode);
  GeomData.acceCur.resize(nNode);

  int ii, jj;
  double val;

  for(ii=0;ii<nNode;ii++)
  {
    for(jj=0;jj<DIM;jj++)
    {
      val = datatemp[ii][jj];

      GeomData.NodePosOrig[ii][jj] = val;
      GeomData.NodePosCur[ii][jj]  = val;
      GeomData.NodePosNew[ii][jj]  = val;
      GeomData.NodePosPrev[ii][jj]  = val;
      GeomData.NodePosPrevCur[ii][jj]  = val;

      GeomData.specValNew[ii][jj]  = 0.0;
      GeomData.specValCur[ii][jj]  = 0.0;

      GeomData.acceNew[ii][jj]  = 0.0;
      GeomData.acceCur[ii][jj]  = 0.0;
    }
  }

  if(DIM == 2)
  {
    for(ii=0;ii<nNode;ii++)
    {
      GeomData.NodePosOrig[ii][2] = 0.0;
      GeomData.NodePosCur[ii][2]  = 0.0;
      GeomData.NodePosNew[ii][2]  = 0.0;
      GeomData.NodePosPrev[ii][2]  = 0.0;
      GeomData.NodePosPrevCur[ii][2]  = 0.0;

      GeomData.specValNew[ii][2]  = 0.0;
      GeomData.specValCur[ii][2]  = 0.0;

      GeomData.acceNew[ii][2]  = 0.0;
      GeomData.acceCur[ii][2]  = 0.0;
    }
  }

  return;
}






void ImmersedRigidSolid::initialise()
{
  SolnData.STAGGERED = STAGGERED;
  
  //printSelf();

  setSolver(1);
  setTimeParam();
  computeInitialAcceleration();

  SolnData.force[0] = initForcePred[1];

  SolnData.forcePrev = SolnData.force;

  return;
}



void ImmersedRigidSolid::setSolver(int slv, int *parm, bool cIO)
{
  prepareMatrixPattern();

  return;
}



void ImmersedRigidSolid::prepareMatrixPattern()
{
  //printf("\n     ImmersedRigidSolid::prepareMatrixPattern()  .... STARTED ...\n");

  //if( std::any_of(PrescMotionTimeFuncs.begin(), PrescMotionTimeFuncs.end(), [](int i){return (i!=-1);}) )
  if(PRESC_MOTION)
  {
    totalDOF = 0;
  }
  else
  {
    int  ii, jj, kk, ll, ind;
  
    //int size_temp = ndofRigidbody;
    
    // increase system size for contacts with Lagrange multipliers
    // contacts for limiting rigid-solid motion
    // two contacts for each DOF
    //   - one in negative direction
    //   - one is positive direction

    //size_temp += 2*ndofRigidbody;

    totalDOF = 0;
    for(ii=0;ii<dofData.size();ii++)
    {
      if( dofData[ii] == 1)
      {
        assy4r.push_back(ii);

        forAssyMat.push_back(vector<int>());
        kk=0;
        for(jj=0;jj<dofData.size();jj++)
        {
          if( dofData[jj] == 1)
            forAssyMat[totalDOF].push_back(kk++);
        }
        totalDOF++;
      }
    }

    //cout << " ndofRigidbody     " << ndofRigidbody << endl;
    //cout << " totalDOF " << totalDOF << endl;
    //cout << " assy4r " << endl;
    //printVector(assy4r);

    bool pp1=0;
    //pp1=1;
    if(pp1)
    {
      printf("   dof to dof connectivity ...:  \n\n");
      for(ii=0;ii<totalDOF;ii++)
      {
        cout << " dof # " << ii << " : ";
        for(jj=0;jj<totalDOF;jj++)
          cout << '\t' << forAssyMat[ii][jj];
        cout << endl;
        }
      printf("\n\n\n");
    }

    Kglobal.resize(totalDOF, totalDOF);
    rhsVec.resize(totalDOF);
    Kglobal.setZero();
    rhsVec.setZero();

    if(!STAGGERED)
    {
      forAssyCoupledHorz.resize(totalDOF);

      ll=0;
      for(ii=0;ii<dofData.size();ii++)
      {
        if( dofData[ii] == 1 )
        {
          for(jj=0; jj<nNode; jj++)
          {
            ind = jj*DIM;
            for(kk=0;kk<DIM;kk++)
              forAssyCoupledHorz[ll].push_back(ind+kk);
          }
          ll++;
        }
      }
      for(ii=0;ii<totalDOF;ii++)
        printVector(forAssyCoupledHorz[ii]);
    }
  }
  //printf("\n     ImmersedRigidSolid::prepareMatrixPattern()  .... FINISHED ...\n\n");

  return;
}


void  ImmersedRigidSolid::setBoundaryConditions(vector<vector<double> >& datatemp)
{
  return;
}


int  ImmersedRigidSolid::updatePointPositions()
{
  if(DIM == 1)
    updatePointPositions1D();
  else if(DIM == 2)
    updatePointPositions2D();
  else if(DIM == 3)
    updatePointPositions3D();

  if(PRESC_MOTION)
    return 1;
  else
    return (totalDOF>0);
}


void  ImmersedRigidSolid::updatePointPositions1D()
{
  return;
}



void  ImmersedRigidSolid::updatePointPositions2D()
{
    int  ee, elm, ii, ind, bb;
    double  thetaNew, thetaCur, omegaNew, omegaCur, fact;
    double  A, B, F, wn, ct, st, t0;

    myPoint  geom, param, rOrig, rNew, vNew, dCentNew, dCentCur, vCentNew, vCentCur;
    myPoint  centroidNew, centroidCur;
    MatrixXd  RotNew(3, 3), RotCur(3, 3);

    geom.setZero();
    param.setZero();
    rOrig.setZero();
    rNew.setZero();
    vNew.setZero();
    dCentNew.setZero();
    dCentCur.setZero();
    vCentNew.setZero();
    vCentCur.setZero();
    centroidNew.setZero();
    centroidCur.setZero();

    RotNew.setZero();
    RotCur.setZero();

    if(totalDOF > 0)
    {
      // values at t_{n+af}
      //dCentCur[0] = SolnData.var1Cur[0];
      //dCentCur[1] = SolnData.var1Cur[1];
      //thetaCur    = SolnData.var1Cur[2];

      //vCentCur[0] = SolnData.var1DotCur[0];
      //vCentCur[1] = SolnData.var1DotCur[1];
      //omegaCur    = SolnData.var1DotCur[2];

      dCentCur[0] = 0.0;
      dCentCur[1] = SolnData.var1Cur[0];
      dCentCur[2] = 0.0;
      thetaCur    = 0.0;

      vCentCur[0] = 0.0;
      vCentCur[1] = SolnData.var1DotCur[0];
      vCentCur[2] = 0.0;
      omegaCur    = 0.0;

      //dCentCur[0] = 0.0;
      //dCentCur[1] = 0.0;
      //thetaCur    = SolnData.var1Cur[0];

      //vCentCur[0] = 0.0;
      //vCentCur[1] = 0.0;
      //omegaCur    = SolnData.var1DotCur[0];

      ct = cos(thetaCur);
      st = sin(thetaCur);

      RotCur(0,0) =  ct; RotCur(0,1) = -st;
      RotCur(1,0) =  st; RotCur(1,1) = ct;

      // values at t_{n+1}
      //dCentNew[0] = SolnData.var1[0];
      //dCentNew[1] = SolnData.var1[1];
      //thetaNew    = SolnData.var1[2];

      //vCentNew[0] = SolnData.var1Dot[0];
      //vCentNew[1] = SolnData.var1Dot[1];
      //omegaNew    = SolnData.var1Dot[2];

      dCentNew[0] = 0.0;
      dCentNew[1] = SolnData.var1[0];
      dCentNew[2] = 0.0;
      thetaNew    = 0.0;

      vCentNew[0] = 0.0;
      vCentNew[1] = SolnData.var1Dot[0];
      vCentNew[2] = 0.0;
      omegaNew    = 0.0;

      //dCentNew[0] = 0.0;
      //dCentNew[1] = 0.0;
      //thetaNew    = SolnData.var1[0];

      //vCentNew[0] = 0.0;
      //vCentNew[1] = 0.0;
      //omegaNew    = SolnData.var1Dot[0];

      ct = cos(thetaNew);
      st = sin(thetaNew);

      RotNew(0,0) =  ct; RotNew(0,1) = -st;
      RotNew(1,0) =  st; RotNew(1,1) = ct;

      //cout << " disp    " << dCentCur[0] << '\t' << dCentCur[1] << '\t' << thetaCur << endl;
      //cout << " velo    " << vCentCur[0] << '\t' << vCentCur[1] << '\t' << omegaCur << endl;

      computeCentroid(0);
      //centroid[2] = 0.0;
      //centroidCur[2] = 0.0;
      //centroidNew[2] = 0.0;

      centroidCur[0] = centroid[0] + dCentCur[0];
      centroidCur[1] = centroid[1] + dCentCur[1];

      centroidNew[0] = centroid[0] + dCentNew[0];
      centroidNew[1] = centroid[1] + dCentNew[1];

      //cout << " centroid    " << centroid[0] << '\t' << centroid[1] << '\t' << centroid[2] << endl;
      //cout << " centroid    " << centroidNew[0] << '\t' << centroidNew[1] << '\t' << centroidNew[2] << endl;
      //cout << " centroid    " << centroidCur[0] << '\t' << centroidCur[1] << '\t' << centroidCur[2] << endl;
      //centroid[0] = 73.0;
      //centroid[1] = 60.0;

      //cout << " ssssssssssss " << DIM << endl;

      for(bb=0;bb<nNode;bb++)
      {
        for(ii=0;ii<DIM;ii++)
          rOrig(ii) = GeomData.NodePosOrig[bb][ii] ;

        // values at t_{n+1}
        rNew = centroid + dCentNew + RotNew*(rOrig - centroid);

        vNew[0] = vCentNew[0] - omegaNew*(rNew[1] - centroidNew[1]);
        vNew[1] = vCentNew[1] + omegaNew*(rNew[0] - centroidNew[0]);

        for(ii=0;ii<DIM;ii++)
        {
          GeomData.NodePosNew[bb][ii]  = rNew[ii];
          GeomData.specValNew[bb][ii]  = vNew[ii];
        }

        // values at t_{n+af}
        rNew = centroid + dCentCur + RotCur*(rOrig - centroid);

        vNew[0] = vCentCur[0] - omegaCur*(rNew[1] - centroidCur[1]);
        vNew[1] = vCentCur[1] + omegaCur*(rNew[0] - centroidCur[0]);

        for(ii=0;ii<DIM;ii++)
        {
          //cout << bb << '\t' << ii << '\t' << rNew[ii] << endl;
          GeomData.NodePosCur[bb][ii]  = rNew[ii];
          GeomData.specValCur[bb][ii]  = vNew[ii];
        }

      }
    }

    A = 80.0*PI/180.0;
    //A = 0.2;
    B = 5.0;
    F = 0.1666*1.1;
    //F = 1.1;

    t0 = 1.0/F;
    //t0 = 2.0*t0;
    t0 = 0.0;
    wn = 1.0;

    //cout << " PRESC_MOTION = " << PRESC_MOTION << endl;

    //if(PRESC_MOTION && (mpapTime.cur > t0)  && ( (id==0) || (id==1)) )
    if(PRESC_MOTION && (mpapTime.cur > t0)  && ( id==0) )
    {
      double  tNew, tCur, af, RPM, N, wn1;
      af = SolnData.td(2);

      tNew = mpapTime.cur;
      tCur = tNew - (1.0-af)*mpapTime.dt;

      /*
      N = 50.0;

      if(id == 0 )
        RPM = -N;
      else
        RPM = -N;

      double  RPS = RPM/60.0; // rotations per second

      double  RPMS = RPS;
      RPMS = RPS/1000.0; // Rotations per millisecond

      double  tOne = abs(1.0/RPMS); // time per one rotation in milliseconds

      //cout << " tOne = " << tOne << endl;

      wn  = -2.0*PI*RPMS; // radians per millisecond

      //if(tNew < tOne)
        //wn1 = wn*sin(0.5*PI*tNew/tOne);
      //else
        wn1 = wn;

      dCentNew[0] = 0.0;
      dCentNew[1] = 0.0;
      thetaNew    = wn1*tNew;
      //thetaNew    = -(0.25*PI)*(1.0+fact);

      vCentNew[0] = 0.0;
      vCentNew[1] = 0.0;
      omegaNew    = wn1;
      //omegaNew    = -(0.25*PI)*(1.0-fact*fact);

      //if(tCur < tOne)
        //wn1 = wn*sin(0.5*PI*tCur/tOne);
      //else
        wn1 = wn;

      dCentCur[0] = 0.0;
      dCentCur[1] = 0.0;
      thetaCur    = wn1*tCur;
      //thetaCur    = -(0.25*PI)*(1.0+fact);

      vCentCur[0] = 0.0;
      vCentCur[1] = 0.0;
      omegaCur    = wn1;
      //omegaCur    = -(0.25*PI)*(1.0-fact*fact);
      */

      /*
      wn = 2.0*PI*F;

      dCentNew[0] = 0.0;
      dCentNew[1] = A*sin(wn*(tNew-t0));
      thetaNew    = 0.0;

      vCentNew[0] = 0.0;
      vCentNew[1] = A*wn*cos(wn*(tNew-t0));
      omegaNew    = 0.0;

      dCentCur[0] = 0.0;
      dCentCur[1] = A*sin(wn*(tCur-t0));
      thetaCur    = 0.0;

      vCentCur[0] = 0.0;
      vCentCur[1] = A*wn*cos(wn*(tCur-t0));
      omegaCur    = 0.0;
      */

      wn = 2.0*PI*F;

      dCentNew[0] = 0.0;
      dCentNew[1] = 0.0;
      thetaNew    = A*sin(wn*(tNew-t0));

      vCentNew[0] = 0.0;
      vCentNew[1] = 0.0;
      omegaNew    = A*wn*cos(wn*(tNew-t0));

      dCentCur[0] = 0.0;
      dCentCur[1] = 0.0;
      thetaCur    = A*sin(wn*(tCur-t0));

      vCentCur[0] = 0.0;
      vCentCur[1] = 0.0;
      omegaCur    = A*wn*cos(wn*(tCur-t0));


      // from Fehmi's paper
      /*
      double  aa=0.2;
      double  d0=0.0;
      double  v0=0.0;
      double  ti=0.0;

      double  d1=0.5*aa*ti*ti + v0*ti + d0;

      if( tNew < ti)
      {
        dCentNew[0] = 0.5*aa*tNew*tNew + v0*tNew + d0;
        dCentNew[1] = 0.0;
        thetaNew    = 0.0;

        vCentNew[0] = aa*tNew + v0;
        vCentNew[1] = 0.0;
        omegaNew    = 0.0;

        dCentCur[0] = 0.5*aa*tCur*tCur + v0*tCur + d0;
        dCentCur[1] = 0.0;
        thetaCur    = 0.0;

        vCentCur[0] = aa*tCur + v0;
        vCentCur[1] = 0.0;
        omegaCur    = 0.0;
      }
      else
      {
        dCentNew[0] = 1.0*(tNew-ti) + d1;
        dCentNew[1] = 0.0;
        thetaNew    = 0.0;

        vCentNew[0] = 1.0;
        vCentNew[1] = 0.0;
        omegaNew    = 0.0;

        dCentCur[0] = 1.0*(tNew-ti) + d1;
        dCentCur[1] = 0.0;
        thetaCur    = 0.0;

        vCentCur[0] = 1.0;
        vCentCur[1] = 0.0;
        omegaCur    = 0.0;
      }
      */

      // hovering insect wing
      /* 
      double  A0=2.5, beta=PI/3.0, phi=0.0, alpha0=PI/4.0, T=0.025;
      wn = 2.0*PI/T;

      fact = 0.5*A0*(1.0+cos(wn*tNew));
      dCentNew[0] = fact*cos(beta);
      dCentNew[1] = fact*sin(beta);
      thetaNew    = alpha0*(1.0-sin(wn*tNew+phi));

      fact = 0.5*A0*(-sin(wn*tNew))*wn;
      vCentNew[0] = fact*cos(beta);
      vCentNew[1] = fact*sin(beta);
      omegaNew    = alpha0*(-cos(wn*tNew+phi)*wn);

      fact = 0.5*A0*(1.0+cos(wn*tCur));
      dCentCur[0] = fact*cos(beta);
      dCentCur[1] = fact*sin(beta);
      thetaCur    = alpha0*(1.0-sin(wn*tCur+phi));

      fact = 0.5*A0*(-sin(wn*tCur))*wn;
      vCentCur[0] = fact*cos(beta);
      vCentCur[1] = fact*sin(beta);
      omegaCur    = alpha0*(-cos(wn*tCur+phi)*wn);
      */

      //cout << " disp    " << dCentNew[0] << '\t' << dCentNew[1] << '\t' << thetaNew << endl;
      //cout << " velo    " << vCentNew[0] << '\t' << vCentNew[1] << '\t' << omegaNew << endl;

      //printf("\t %12.6f \t %12.6f \t %12.6f \n",tNew, tCur, wn);
      //printf("\t %12.6f \t %12.6f \t %12.6f \n",dCentCur[0], dCentCur[1], thetaCur);
      //printf("\t %12.6f \t %12.6f \t %12.6f \n",vCentCur[0], vCentCur[1], omegaCur);
      //cout << " disp    " << dCentCur[0] << '\t' << dCentCur[1] << '\t' << thetaCur << endl;
      //cout << " velo    " << vCentCur[0] << '\t' << vCentCur[1] << '\t' << omegaCur << endl;

      ct = cos(thetaNew);
      st = sin(thetaNew);

      RotNew(0,0) =  ct; RotNew(0,1) = -st;
      RotNew(1,0) =  st; RotNew(1,1) =  ct;

      ct = cos(thetaCur);
      st = sin(thetaCur);

      RotCur(0,0) =  ct; RotCur(0,1) = -st;
      RotCur(1,0) =  st; RotCur(1,1) =  ct;

      computeCentroid(0);

      centroidCur[0] = centroid[0] + dCentCur[0];
      centroidCur[1] = centroid[1] + dCentCur[1];

      centroidNew[0] = centroid[0] + dCentNew[0];
      centroidNew[1] = centroid[1] + dCentNew[1];

      //cout << " centroid    " << centroid[0] << '\t' << centroid[1] << '\t' << centroid[2] << endl;
      //cout << " centroid    " << centroidNew[0] << '\t' << centroidNew[1] << '\t' << centroidNew[2] << endl;
      //cout << " centroid    " << centroidCur[0] << '\t' << centroidCur[1] << '\t' << centroidCur[2] << endl;

      for(bb=0; bb<nNode; bb++)
      {
        for(ii=0;ii<DIM;ii++)
          rOrig(ii) = GeomData.NodePosOrig[bb][ii] ;

        // values at t_{n+1}
        rNew = centroid + dCentNew + RotNew*(rOrig - centroid);

        vNew[0] = vCentNew[0] - omegaNew*(rNew[1] - centroidNew[1]);
        vNew[1] = vCentNew[1] + omegaNew*(rNew[0] - centroidNew[0]);

        for(ii=0;ii<DIM;ii++)
        {
          GeomData.NodePosNew[bb][ii]  = rNew[ii];
          GeomData.specValNew[bb][ii]  = vNew[ii];
        }

        // values at t_{n+af}
        rNew = centroid + dCentCur + RotCur*(rOrig - centroid);

        vNew[0] = vCentCur[0] - omegaCur*(rNew[1] - centroidCur[1]);
        vNew[1] = vCentCur[1] + omegaCur*(rNew[0] - centroidCur[0]);

        for(ii=0;ii<DIM;ii++)
        {
          GeomData.NodePosCur[bb][ii]  = rNew[ii];
          GeomData.specValCur[bb][ii]  = vNew[ii];
        }
      }
    }

  return;
}




void  ImmersedRigidSolid::updatePointPositions3D()
{
    if(totalDOF == 0)
      return;

    int  ee, elm, ii, ind, bb, fact;
    double  thetaNew, thetaCur, omegaNew, omegaCur;
    double  A, B, F, wn, ct, st, t0;

    myPoint  geom, param, rOrig, rNew, vNew, dCentNew, dCentCur, vCentNew, vCentCur;
    MatrixXd  RotNew(3, 3), RotCur(3, 3);
    
    RotNew.setZero();
    RotCur.setZero();
    
    bool pres_motion=true;
    pres_motion = false;

    if(totalDOF > 0)
    {
      // values at t_{n+af}

      dCentCur[0] = 0.0;
      dCentCur[1] = SolnData.var1Cur[0];
      dCentCur[2] = 0.0;
      thetaCur    = 0.0;

      vCentCur[0] = 0.0;
      vCentCur[1] = SolnData.var1DotCur[0];
      vCentCur[2] = 0.0;
      omegaCur    = 0.0;

      ct = cos(thetaCur);
      st = sin(thetaCur);

      RotCur(0,0) =  ct; RotCur(0,1) = st;
      RotCur(1,0) = -st; RotCur(1,1) = ct;

      // values at t_{n+1}
      //dCentNew[0] = SolnData.var1[0];
      //dCentNew[1] = SolnData.var1[1];
      //thetaNew    = SolnData.var1[2];

      //vCentNew[0] = SolnData.var1Dot[0];
      //vCentNew[1] = SolnData.var1Dot[1];
      //omegaNew    = SolnData.var1Dot[2];

      dCentNew[0] = 0.0;
      dCentNew[1] = SolnData.var1[0];
      dCentNew[2] = 0.0;
      thetaNew    = 0.0;

      vCentNew[0] = 0.0;
      vCentNew[1] = SolnData.var1Dot[0];
      vCentNew[2] = 0.0;
      omegaNew    = 0.0;

      ct = cos(thetaNew);
      st = sin(thetaNew);

      RotNew(0,0) =  ct; RotNew(0,1) = st;
      RotNew(1,0) = -st; RotNew(1,1) = ct;

      //cout << " disp    " << dCentCur[0] << '\t' << dCentCur[1] << '\t' << thetaCur << endl;
      //cout << " velo    " << vCentCur[0] << '\t' << vCentCur[1] << '\t' << omegaCur << endl;

      computeCentroid(0);

      //cout << " ssssssssssss " << DIM << endl;

      for(bb=0;bb<nNode;bb++)
      {
        for(ii=0;ii<DIM;ii++)
          rOrig(ii) = GeomData.NodePosOrig[bb][ii] ;

        // values at t_{n+1}
        //rNew = centroid + dCentNew + RotNew*(rOrig - centroid);

        //vNew[0] = vCentNew[0] - omegaNew*(rNew[0] - centroid[0]);
        //vNew[1] = vCentNew[1] + omegaNew*(rNew[1] - centroid[1]);
        //
        rNew = dCentNew + rOrig;
        
        vNew = vCentNew;

        for(ii=0;ii<DIM;ii++)
        {
          GeomData.NodePosNew[bb][ii]  = rNew[ii];
          GeomData.specValNew[bb][ii]  = vNew[ii];
        }

        // values at t_{n+af}
        //rNew = centroid + dCentCur + RotCur*(rOrig - centroid);

        //vNew[0] = vCentCur[0] - omegaCur*(rNew[0] - centroid[0]);
        //vNew[1] = vCentCur[1] + omegaCur*(rNew[1] - centroid[1]);
        
        rNew = dCentCur + rOrig;

        vNew = vCentCur;

        for(ii=0;ii<DIM;ii++)
        {
          //cout << bb << '\t' << ii << endl;
          GeomData.NodePosCur[bb][ii]  = rNew[ii];
          GeomData.specValCur[bb][ii]  = vNew[ii];
        }

      }
    }

  return;
}





void ImmersedRigidSolid::computeInitialAcceleration()
{
  //for(int ii=0;ii<ndof;ii++)
    //SolidSolnData.acce[ii] = (SolidSolnData.force[ii] - C(ii,ii)*SolidSolnData.velo[ii] - K(ii,ii)*SolidSolnData.disp[ii])/M(ii,ii); // replace 0.0 with force if it is not zero
}


void ImmersedRigidSolid::updateForce()
{
  computeTotalForce();

  SolnData.forceTemp.setZero();

  //cout << totalForce[0] << '\t' << totalForce[1] << endl;

  //SolnData.forceTemp[0] = totalForce[0];  // Fx
  //SolnData.forceTemp[1] = totalForce[1];  // Fy
  //SolnData.forceTemp[2] = totalMoment[2]; // Mz

  SolnData.forceTemp[0] = totalForce[1];  // Fy
  SolnData.forceTemp[1] = 0.0;  // lamba - contact
  SolnData.forceTemp[2] = 0.0;

  SolnData.interpolateForce();

  //printf("\n\n");  printVector(SolidSolnData.forceTemp);  printf("\n\n");
  //printf("\n\n");  printVector(SolidSolnData.forceCur);  printf("\n\n");

  return;
}





void ImmersedRigidSolid::updateForce(double* data)
{
  //cout << totalForce[0] << '\t' << totalForce[1] << endl;

  if(DIM == 2)
  {
    totalForce[0]  = -data[0];
    totalForce[1]  = -data[1];
    totalForce[2]  = -data[5];
  }
  else if(DIM == 3)
  {
    totalForce[0]  = -data[0];
    totalForce[1]  = -data[1];
    totalForce[2]  = -data[2];
    totalForce[3]  = -data[3];
    totalForce[4]  = -data[4];
    totalForce[5]  = -data[5];
  }

  //PetscPrintf(MPI_COMM_WORLD, "   Forces for solid %5d =  %12.6f \t %12.6f \t %12.6f \n\n", id, totalForce[0], totalForce[1], totalForce[2]);

  //SolnData.forceTemp[0] = totalForce[0];  // Fx
  //SolnData.forceTemp[1] = totalForce[1];  // Fy
  //SolnData.forceTemp[2] = totalMoment[2]; // Mz

  SolnData.forceTemp[0] = totalForce[1];  // Fy
  //SolnData.forceTemp[0] = totalForce[2];  // Mz
  SolnData.forceTemp[1] = 0.0;  // lamba - contact
  SolnData.forceTemp[2] = 0.0;

  SolnData.interpolateForce();

  //printf("\n\n");  printVector(SolidSolnData.forceTemp);  printf("\n\n");
  //printf("\n\n");  printVector(SolidSolnData.forceCur);  printf("\n\n");

  return;
}


void ImmersedRigidSolid::updateDisplacement(double* data)
{
  // update solution vector     
  for(int ii=0;ii<totalDOF;ii++)
  {
    //SolnData.var1Dot[assy4r[ii]] = data[ii];
    //SolnData.var1Dot[ii] = data[ii];
    SolnData.var1[ii] += data[ii];
  }
  
  updateIterStep();

  return;
}


void ImmersedRigidSolid::solveTimeStep()
{
  if( totalDOF > 0 )
  {
    //printf("\n Solving Immersed Rigid Solid \n");
  
    tol = 1.0e-6;

    for(int iter=0;iter<10;iter++)
    {
      calcStiffnessAndResidual();

      cout << iter << '\t' << rNorm << endl;

      if(rNorm < tol)
        break;

      factoriseSolveAndUpdate();
    
      updateIterStep();

      //printVector(SolnData.var1Dot);
      //printVector(SolnData.var1);
    }
    //printf("\n Solving Immersed Rigid Solid ..... DONE  \n\n");
  }

  return;
}



int ImmersedRigidSolid::applyBoundaryConditions()
{   
  return 0;
}

int ImmersedRigidSolid::applyExternalForces()
{
  return 0;
}



void ImmersedRigidSolid::resetMatrixAndVector()
{
  return;
}

/*
  double F0=9.0;
  //double F0=0.0035;
  //F0 = 1.2*0.196349540849362*981.0;
  //F0 = 0.0072*1.7671e+03*9.81*0.001;
  //F0 = 0.001996*0.0314158*981.0;
  //F0 = 1.2*0.0314158*981.0;
*/

//
int ImmersedRigidSolid::calcStiffnessAndResidual(int solver_type, bool zeroMtx, bool zeroRes)
{
  if(firstIter)
    rNorm = -1.0;

  int ii, jj, ind;
  double  y1, y2, fact, tol=1.e-12, lamn, gn, af, cn, g0, disp, beta=1.0;
  
  Kglobal.setZero();
  rhsVec.setZero();
  
  //cout << " force = " << SolnData.forceCur(0) << '\t' << SolnData.forceCur(1) << endl;
  //cout << " values = " << SolnData.var1Cur(0) << '\t' << SolnData.forceCur(1) << endl;

  //cout << SolnData.td[5] << '\t' << SolnData.td[6] << '\t' << SolnData.td[7] << endl;
  
  //printMatrix(M);  printMatrix(C);  printMatrix(K);

  Kglobal(0,0) = SolnData.td[5]*matM(1,1) + SolnData.td[6]*matC(1,1) + SolnData.td[7]*matK(1,1);
  rhsVec(0)    = SolnData.forceCur(0) - matM(1,1)*SolnData.var1DotDotCur[0] - matC(1,1)*SolnData.var1DotCur[0] - matK(1,1)*SolnData.var1Cur[0];

  //Kglobal(0,0) = SolnData.td[5]*matM(2,2) + SolnData.td[6]*matC(2,2) + SolnData.td[7]*matK(2,2);
  //rhsVec(0)    = SolnData.forceCur(0) - matM(2,2)*SolnData.var1DotDotCur[0] - matC(2,2)*SolnData.var1DotCur[0] - matK(2,2)*SolnData.var1Cur[0];

  rhsVec(0)   += preLoad[1] ;

  //cout << " contact force = " << veloCur(1) << '\t' << velo(1) << '\t' << veloPrev(1) << endl;

  //Kglobal(0,0) /= SolnData.td[10];
  
  //cout << Kglobal(0,0) << '\t' << rhsVec(0) << endl;

  //////////////////////////////////////////
  // terms related to contact
  //////////////////////////////////////////

  disp = SolnData.var1Cur[0];

  y1 = 21.1; // top-most point on the lower block
  g0 = 0.05; // initial gap
  y2 = y1 + g0 + disp;
  cn = 1000.0;

  //gn = y2 - y1 - g0;
  gn = - disp; // gn is penetration variable

  //cout << " displacement = " << disp << '\t' << lamn << endl;

  //cout << " penetration = " << gn << '\t' << lamn << endl;

  af   = SolnData.td[2];
  fact = 1.0/SolnData.td[10];
  fact = 1.0;

if(totalDOF > 1)
{
  lamn = SolnData.var1Cur[1] ;
  //lamn = 0.0;
  //lamn = SolnData.var1DotCur[1];

  //cout << " penetration = " << gn << '\t' << lamn << endl;

  if( (lamn + cn*gn) > tol)
  {
    //cout << " Contact 1 is active " << endl;

    Kglobal(0,1) -= af*beta;
    Kglobal(1,0) -= af*beta*fact;

    //Kglobal(0,0) += 2.0*cn*af;

    rhsVec(0)   += lamn;
    //rhsVec(0)   += 2.0*cn*gn;
    rhsVec(1)   -= gn;  // residual: contact force
  }
  else
  {
    //cout << " Contact 1 is inactive " << endl;
    
    Kglobal(1,1) += af;
    rhsVec(1)    -= lamn;
  }

  disp = SolnData.var1Cur[0];

  //gn = disp-5.0;
  //gn = disp-1.55;
  gn = disp-20.4;
  //gn = disp-80.0*PI/180.0;

  double lamn2 = SolnData.var1Cur[2];
  //lamn2 = 0.0;

  //cout << " penetration = " << gn << '\t' << lamn2 << endl;

  if( (lamn2 + cn*gn) > tol)
  {
    //cout << " Contact 2 is active " << endl;

    Kglobal(0,2) += af*beta;
    Kglobal(2,0) += af*beta*fact;

    //Kglobal(0,0) += 2.0*cn*af;

    rhsVec(0)   -= lamn2 ;
    //rhsVec(0)   -= 2.0*cn*gn;
    rhsVec(2)   -= gn;  // residual: contact force
  }
  else
  {
    //cout << " Contact 2 is inactive " << endl;

    Kglobal(2,2) += af;
    rhsVec(2)    -= lamn2;
  }
}
  //printMatrix(Kglobal);
  //printf("\n\n");
  //printVector(rhsVec);

  //cout << " forceCur " << forceCur(0) << '\t' << Flocal(0) << '\t' << Klocal(0,0) << endl;

  firstIter = false;
  rNormPrev = rNorm;
  rNorm     = rhsVec.norm();
  //iterCount++;

  //printf(" ImmersedRigidSolid .... norm  %11.4e\n", rNorm);
  
  return 0;
}
//



/*
int ImmersedRigidSolid::calcStiffnessAndResidual(int solver_type, bool zeroMtx, bool zeroRes)
{
  if(firstIter)
    rNorm = -1.0;

  int ii, jj, ind;
  double  y1, y2, fact, tol=1.e-8, lamn, gn, af, cn, g0, disp;

  Kglobal.setZero();
  rhsVec.setZero();
  
  //cout << " force = " << SolnData.forceCur(0) << '\t' << SolnData.forceCur(1) << endl;
  //cout << " values = " << SolnData.var1Cur(0) << '\t' << SolnData.forceCur(1) << endl;

  //cout << SolnData.td[5] << '\t' << SolnData.td[6] << '\t' << SolnData.td[7] << endl;
  
  //printMatrix(M);  printMatrix(C);  printMatrix(K);

  //Kglobal(0,0) = SolnData.td[5]*matM(1,1) + SolnData.td[6]*matC(1,1) + SolnData.td[7]*matK(1,1);
  //rhsVec(0)    = SolnData.forceCur(0) - matM(1,1)*SolnData.var1DotDotCur[0] - matC(1,1)*SolnData.var1DotCur[0] - matK(1,1)*SolnData.var1Cur[0];

  Kglobal = SolnData.td[5]*matM + SolnData.td[6]*matC + SolnData.td[7]*matK;
  rhsVec  = SolnData.forceCur - matM*SolnData.var1DotDotCur - matC*SolnData.var1DotCur - matK*SolnData.var1Cur;

  double Rad=0.1, As=PI*Rad*Rad;

  double  F0 = 0.001996*As*981.0;
  
  rhsVec[1] -= F0;


if(totalDOF > 10)
{
  //////////////////////////////////////////
  // terms related to contact
  //////////////////////////////////////////

  disp = SolnData.var1Cur[0];

  cn = 1.0e6;

  gn = - disp; // gn is penetration variable

  af   = SolnData.td[2];

  lamn = 2.0*cn*gn ;

  cout << " penetration = " << gn << '\t' << lamn << endl;

  if( gn > tol)
  {
    cout << " Contact 1 is active " << endl;

    Kglobal(0,0) += 2.0*cn*af;

    rhsVec(0)   += lamn;
  }
  else
  {
    cout << " Contact 1 is inactive " << endl;
  }

  disp = SolnData.var1Cur[0];

  gn = disp-2.0;

  lamn = 2.0*cn*gn ;

  cout << " penetration = " << gn << '\t' << lamn << endl;

  if( gn > tol )
  {
    cout << " Contact 2 is active " << endl;

    Kglobal(0,0) += 2.0*cn*af;

    rhsVec(0)   -= lamn;
  }
  else
  {
    cout << " Contact 2 is inactive " << endl;
  }
}
  //printMatrix(Kglobal);
  //printf("\n\n");
  //printVector(rhsVec);

  //cout << " forceCur " << forceCur(0) << '\t' << Flocal(0) << '\t' << Klocal(0,0) << endl;

  firstIter = false;
  rNormPrev = rNorm;
  rNorm     = rhsVec.norm();
  //iterCount++;

  //printf(" ImmersedRigidSolid .... norm  %11.4e\n", rNorm);
  
  return 0;
}
*/




void ImmersedRigidSolid::calcCouplingMatrices()
{
  //////////////////////////////////////////
  // off-diagonal matrices
  //////////////////////////////////////////

    int  aa, nlb, ii, jj, kk, ind1, ind2;
    MatrixXd  Ktemp, Ktemp2;

    Khorz.resize(ndofRigidbody, nNode*DIM);
    Khorz.setZero();

    Kvert.resize(nNode*DIM, ndofRigidbody);
    Kvert.setZero();

    for(aa=0;aa<ImmIntgElems.size();aa++)
    {
      nlb = ImmIntgElems[aa]->pointNums.size();

      //cout << " aa = " << aa << endl;
      ImmIntgElems[aa]->computeKhorzKvertRigid(Ktemp, Ktemp2);
      //printMatrix(Ktemp);
      //printf("\n\n");
      //printMatrix(Ktemp2);
      //printf("\n\n");

      for(ii=0;ii<ndofRigidbody;ii++)
      {
        for(jj=0;jj<nlb;jj++)
        {
          ind1 = ImmIntgElems[aa]->pointNums[jj] * DIM;
          ind2 = jj*DIM;

          for(kk=0;kk<DIM;kk++)
          {
            Khorz(ii, ind1+kk) += Ktemp(ii, ind2+kk);
            Kvert(ind1+kk, ii) += Ktemp2(ind2+kk, ii);
          }
        }
      }
    }

    Khorz *= -SolnData.td[2];
    Kvert *= -SolnData.td[2];

  //printMatrix(Khorz);
}





int ImmersedRigidSolid::factoriseSolveAndUpdate()
{
  //cout << " totalDOF = " << totalDOF << endl;

  VectorXd  sln(totalDOF);
  
  sln = Kglobal.fullPivLu().solve(rhsVec);
  
  //printVector(sln);

  for(int ii=0;ii<sln.rows();ii++)
    //SolnData.var1Dot[assy4r[ii]] += sln[ii];
    //SolnData.var1Dot[ii] += sln[ii];
    SolnData.var1[ii] += sln[ii];

  return 0;
}



/*
int ImmersedRigidSolid::assembleGlobalMatrixAndVector(int start1, int start2, SparseMatrixXd& mtx, double* rhs)
{
  int ii, jj, r, c, kk, ll;

  if(STAGGERED)
  {
    for(ii=0;ii<totalDOF;ii++)
    {
      kk = assy4r[ii];

      r  = start2 + ii;

      rhs[r] += Flocal[kk];

      for(jj=0;jj<totalDOF;jj++)
        mtx.coeffRef(r, start2+jj) += Klocal(kk, assy4r[jj]);
    }
  }
  else
  {
    calcStiffnessAndResidual(1, 0, 0);
    calcCouplingMatrices();

    for(ii=0;ii<totalDOF;ii++)
    {
      r  = start2 + ii;
      kk = assy4r[ii];

      rhs[r] += Flocal[kk];

      for(jj=0;jj<totalDOF;jj++)
        mtx.coeffRef(r, start2+jj) += Klocal(kk, assy4r[jj]);

      for(jj=0;jj<forAssyCoupledHorz[ii].size();jj++)
      {
        c = start1 + forAssyCoupledHorz[ii][jj];

        mtx.coeffRef(r, c) += Khorz(kk, jj);
        mtx.coeffRef(c, r) += Kvert(jj, kk);
      }
    }
  }

  return 0;
}
*/



//
int ImmersedRigidSolid::assembleGlobalMatrixAndVector(int ind1, int ind2, SparseMatrixXd& mtx, double* rhs)
{
  if(totalDOF <= 0)
    return 1;


  calcStiffnessAndResidual(1, 0, 0);
  calcCouplingMatrices();

  int ii, jj, r, c, ll;
  
  for(ii=0;ii<totalDOF;ii++)
  {
    r  = ind2 + ii;

    rhs[r] += rhsVec[ii];

    for(jj=0;jj<totalDOF;jj++)
      mtx.coeffRef(r, ind2+jj) += Kglobal(ii, jj);
  }

    r  = ind2 + 0;

    for(jj=0;jj<forAssyCoupledHorz[0].size();jj++)
    {
      c = ind1 + forAssyCoupledHorz[0][jj];

      mtx.coeffRef(r, c) += Khorz(1, jj);
      mtx.coeffRef(c, r) += Kvert(jj, 1);
    }

  return 1;
}
//






int ImmersedRigidSolid::assembleGlobalMatrixAndVectorCutFEM(int start1, int start2, SolverPetsc* solverTemp)
{
  calcStiffnessAndResidual(1, 0, 0);
  //printMatrix(Klocal);
  //printVector(SolnData.forceCur);
  //printVector(Flocal);
  //cout << " totalDOF = " << totalDOF << endl;
  //printVector(assy4r);
  //cout << " disp = " << SolnData.var1Cur(1) << endl;

  cout << " ImmersedRigidSolid::AssembleGlobalMatrixAndVectorCutFEM ... needs to be corrected " << endl;

  int ii, jj, r, c, kk, ll;

  for(ii=0;ii<totalDOF;ii++)
  {
    r  = start1 + ii;
    kk = assy4r[ii];

    //rhs[r] += Flocal[kk];

    //for(jj=0;jj<totalDOF;jj++)
      //mtx.coeffRef(r, start2+jj) += Klocal(kk, assy4r[jj]);
  }

  return 0;
}



void ImmersedRigidSolid::assembleElementVector(int ind, bool flag, double* rhs)
{
  for(int ii=0;ii<ndof;ii++)
  {
    rhs[ii] += Flocal(ii);
  }

  return;
}


void ImmersedRigidSolid::writeOutput()
{
    if(BC_ENFORCE_TYPE == BC_ENFORCE_TYPE_PENALTY)
    {
      cout << " can't print the output data when PENALTY method is used for Immersed Bodies " << endl;
      return;
    }

    int ii, jj, bb, ee, id, type, nodenum, dof;

  for(bb=0; bb<OutputData.size(); bb++)
  {
    type    =  OutputData[bb][0];
    nodenum =  OutputData[bb][1] - 1;
    dof     =  OutputData[bb][2] - 1;

    char        tmp[100];
    MyString    tmpStr;

    //cout << id << '\t' << type << '\t' << nodenum << '\t' << dof << '\t' << fact << endl;
    //PetscPrintf(MPI_COMM_WORLD, "   Forces =  %12.6f \t %12.6f \t %12.6f \n\n", totalForce[0], totalForce[1], totalForce[2]);

    switch(type)
    {
      case  1 : // total force on the

            sprintf(tmp," \t %12.6E", totalForce[dof]);

      break;

      case  2 : // force on individual boundary point

             sprintf(tmp," \t %12.6E", SolnData.force[dof]);

      break;

      case  3 : // displacement

             sprintf(tmp," \t %12.6E", SolnData.var1[dof]);

      break;

      case  4 : // velocity

             sprintf(tmp," \t %12.6E", SolnData.var1Dot[dof]);

      break;

      case  5 : // acceleration

             sprintf(tmp," \t %12.6E", SolnData.var1DotDot[dof]);

      break;

      default : 

             prgError(1,"HBSplineFEM::writeImmersedSolidOutput()","invalid value of 'type'!");

      break;
    }

    tmpStr.append(tmp);

    prgWriteToTFile(tmpStr);
  }

  return;
}



void ImmersedRigidSolid::postProcess(int index)
{

  vtkSmartPointer<vtkXMLPolyDataWriter>  writerPolyData =   vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  char fname1[50];
  sprintf(fname1,"%s%d%s", "immersedpoly-", id,".vtp");

  writerPolyData->SetFileName(fname1);

#if VTK_MAJOR_VERSION == 5
    writerPolyData->SetInput(polyDataVTK);
#else
    writerPolyData->SetInputData(polyDataVTK);
#endif

  writerPolyData->Write();

  return;
}


