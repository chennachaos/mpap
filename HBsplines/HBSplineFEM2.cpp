
#include "HBSplineFEM.h"
#include "MpapTime.h"
#include "ImmersedIntegrationElement.h"
#include "myDataIntegrateCutFEM.h"
#include "QuadratureUtil.h"
#include "BasisFunctionsLagrange.h"

extern MpapTime mpapTime;


void  HBSplineFEM::computeTotalForce(int index)
{
  return;
}


void  HBSplineFEM::ImmersedBoundaryBodyForce1D()
{
    //
    //if(mpapTime.cur <= mpapTime.dt)
    //{
      cout << mpapTime.cur << '\t' << mpapTime.dt << endl;
      for(int ii=0;ii<fluidDOF;ii++)
        SolnData.bodyForce(ii) = 0.0;
      //return;
    //}
    // code for 1D Posson problem

    int  aa, bb, ii, jj, nlocal, ind, nodenum, size, ee, dir, numPoints;

    double   fact, val, xx, yy, incr1, incr2, *knots[2];
    double   t0, tf, tdiff, funcVal[2], temp, val1, val2, ds, h, C;
    myPoint  knotIncr, knotBegin;

    vector<int>  bfs;

    node*  ndtmp;
    numPoints = 1;
    geom[0] = 0.5;

    //C = LevelSetFunc[0][6];

    nlocal = degree[0] + 1;
    VectorXd  N(nlocal), NN, dN_dx(nlocal), d2N_dx2(nlocal);

        for(bb=0;bb<numPoints;bb++)
        {
            xx = geom[0];
            ee = findCellNumber(geom);

            //printf("xx = %12.6f,  cell = %5d, \n", geom[0], ee);

            ndtmp = elems[ee];

            knots[0] = ndtmp->GetKnots(0);

            knotBegin = ndtmp->GetKnotBegin();
            knotIncr  = ndtmp->GetKnotIncrement();


            bfs = ndtmp->GlobalBasisFuncs;

            geometryToParametric(geom, param);

            //printf("param ... = %12.6f \t %12.6f \t %12.6f, \n \n ", param, knots[0][0], knots[0][2]);

            GeomData.computeBasisFunctions1D(knotBegin, knotIncr, param, N, dN_dx, d2N_dx2);

            if(ndtmp->GetParent() == NULL )
              NN = N;
            else
            {
              NN = ndtmp->SubDivMat * N;
              //NN = NN/NN.sum();
            }

            //printf("BasisFuns \n");
            //for(ii=0;ii<bfs.size();ii++)
              //printf(" \t %5d \t %12.6f \t %12.6f \t %12.6f \n ", bfs[ii],  NN[ii], dN_dx[ii], d2N_dx2[ii]);

            //val2 = ndtmp->computeValueOld(0, N);
            //val1 = ndtmp->computeValueOld(0, d2N_dx2);
            //printf("\t val1  %12.8f \t   %12.8f \n", val1, ndtmp->computeValueOld(0, dN_dx));

            h = knots[0][2];
            //val1 = C/h;
            val1 = - GeomData.FluidProps[2]/h;
            val1 = 1.0/h;
            //val1 = (0.0-soln(totalDOF));
            //val1 = -val2/h;
            printf("\t val1  %12.8f \n", val1);

            for(ii=0;ii<bfs.size();ii++)
            {
               SolnData.bodyForce(bfs[ii]) += val1*NN(ii);
              //rhsVec(bfs[ii]) += val1*d2N_dx2(ii);
            }

            //rhsVec(totalDOF) -= ndtmp->computeValue(0, N);
            //rhsVec(totalDOF) -= FluidSolnData.soln(totalDOF);

            //printf("BasisFuns \n");
            //for(ii=0;ii<bfs.size();ii++)
               //FluidSolnData.bodyForce(bfs[ii]) -= val1*N(ii);
            //printf("\n\n");
        }

      //for(int ii=0;ii<totalDOF;ii++)
        //FluidSolnData.bodyForce(ii) -= 1.0;
    
    //printf("\n bodyForce \n");
    //printVector(FluidSolnData.bodyForce);

    return;
}
//


void  HBSplineFEM::computeTotalBodyForce(int index)
{
  double val, forceTotal[3];
  int dir, ee;
   
  for(dir=0;dir<DIM;dir++)
  {
    val = 0.0;
    for(ee=0;ee<activeElements.size();ee++)
      val += elems[activeElements[ee]]->computeTotalBodyForce(index, dir);

    forceTotal[dir] = val;
  }

  char        tmp[200];
  MyString    tmpStr;
    
  sprintf(tmp," \t %12.6E \t %12.6E ", forceTotal[0], forceTotal[1]);
  tmpStr.append(tmp);

  prgWriteToTFile(tmpStr);

  return;
}



void  HBSplineFEM::ImmersedBoundaryBodyForceLSFEM()
{
    int ee, elm, aa, bb;
    double  temp[2], fact=0.0, param[2], af, VelSolid[2], wn, t0, A, F;
    
    ImmersedIntegrationElement  *lme;
    
    t0 = 0.0;

    if(mpapTime.cur < 100000.0)
    {
      SolnData.bodyForce.setZero();

      VelSolid[0] = 0.0;
      VelSolid[1] = 0.0;

      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      {
        for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
        {
          // no need to update the positions and element numbers
       
          // compute the body force

          //cout << bb << '\t' << aa << endl;
          lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

          lme->computeBodyForce(1, VelSolid);
          //lme->mapDataToGlobalBodyForceVector(0, temp);
        }
      }

      af = SolnData.td(2); // alpha_f for the fluid

      SolnData.bodyForceCur = af*SolnData.bodyForce + (1.0-af)*SolnData.bodyForcePrev;

      computeTotalBodyForce(0);
      computeTotalBodyForce(2);
      //cout << " kkkkkkkkkkkk " << endl;
      //printVector(FluidSolnData.bodyForce);
    }
    /*
    else
    {
       bool flag1;
       flag1 = 0;
       
       if(flag1)
       {
         A = FluidSolnData.ElemProp.data[5];
         F = FluidSolnData.ElemProp.data[6];

         wn = 2.0*PI*0.16*F;
         fact = A*sin(wn*(mpapTime.cur-t0));
         
         VelSolid[1] = A*wn*cos(wn*(mpapTime.cur-t0));
       }
       else
       {       
         computeTotalBodyForce(0); // compute total body force
         
         //RigidBodyForce[0] *= -1.0;
         //RigidBodyForce[1] *= -1.0;
       
         ImmersedBodyObjects.setTimeParam();
         
         ImmersedBodyObjects.timeUpdate();

         ImmersedBodyObjects.SolveTimeStep(RigidBodyForce);

         fact = ImmersedBodyObjects.GetDisplacement(1);
       
         VelSolid[1] = ImmersedBodyObjects.GetVelocity(1);
       }

       FluidSolnData.bodyForce.setZero();

       for(ee=0;ee<IBpoints.size();ee++)
       {
          // update the positions and element numbers
       
          param[0] = IBpoints[ee]->GetPositionOrig(0);
          param[1] = IBpoints[ee]->GetPositionOrig(1);

          //param[0] = IBpoints[ee]->GetPositionOld(0);
          //param[1] = IBpoints[ee]->GetPositionOld(1);

          param[1] += fact;

          IBpoints[ee]->swapPositions();

          //IBpoints[ee]->SetPositionCur(0, param[0]);
          IBpoints[ee]->SetPositionCur(1, param[1]);

          elm = findCellNumber(param[0], param[1]);

          IBpoints[ee]->SetElementNum(elm);

          IBpoints[ee]->elem = elems[elm];

          //printf("xx = %12.6f, yy = %12.6f, cell = %5d, \n", param[0], param[1], ee);
          geometryToParametric(param[0], param[1], param);
          //printf("xx = %12.6f, yy = %12.6f, cell = %5d, \n", param[0], param[1], ee);

          //IBpoints[ee]->SetParam(0, param[0]);
          IBpoints[ee]->SetParam(1, param[1]);
       
          // compute the body force

          IBpoints[ee]->computeBodyForce(0, VelSolid);

          IBpoints[ee]->mapDataToGlobalBodyForceVector(0, temp);
       }
    }
    */
    //af = FluidSolnData.td(2); // alpha_f for the fluid

    //FluidSolnData.bodyForceCur = af*FluidSolnData.bodyForce + (1.0-af)*FluidSolnData.bodyForcePrev;
    //FluidSolnData.bodyForceCur = FluidSolnData.bodyForce ;

    //
    //printVector(FluidSolnData.bodyForce);
    
    //FluidSolnData.bodyForce = beta*FluidSolnData.bodyForce + (1.0-beta)*FluidSolnData.bodyForcePredPrev;

    //FluidSolnData.bodyForcePred = 2.0*FluidSolnData.bodyForce - FluidSolnData.bodyForcePrev;
    
    //FluidSolnData.bodyForceCur = af*FluidSolnData.bodyForcePred + (1.0-af)*FluidSolnData.bodyForcePredPrev;
    //    

    //computeTotalBodyForce(0); // compute total body force current

    //char        tmp[200];
    //MyString    tmpStr;
    
    //printf(" \t %12.6E \t %12.6E \t %12.6E \n", RigidBodyForce[0], RigidBodyForce[1], ImmersedBodyObjects.GetDisplacement(1));

    //sprintf(tmp," \t %12.6E \t %12.6E \t %12.6E \t %12.6E \t %12.6E ", RigidBodyForce[0], RigidBodyForce[1], fact, VelSolid[1], ImmersedBodyObjects.GetAcceleration(1));
    //sprintf(tmp," \t %12.6E \t %12.6E", RigidBodyForce(0,0), RigidBodyForce(0,1));
    //tmpStr.append(tmp);

    //prgWriteToTFile(tmpStr);

    return;
}





void  HBSplineFEM::ImmersedBoundaryBodyForce()
{
  return;
}

void  HBSplineFEM::ImmersedBoundaryBodyForce2D()
{
  return;
}


void  HBSplineFEM::solveSolidProblem()
{
  IB_MOVED = false;
  for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
  {
    //cout << " Lagrange multipliers " << endl;        printVector(&(FluidSolnData.var3(0)), IBDOF);
    //ImmersedBodyObjects[bb]->updateForce(&(totalForce(0)));

    //cout << " kkkkkkkkkkk " << endl;
    ImmersedBodyObjects[bb]->SolveTimeStep();

    IB_MOVED = (IB_MOVED || ImmersedBodyObjects[bb]->updatePointPositions() );
  }

  return;
}




void  HBSplineFEM::writeNodalData()
{
  writeFluidOutput();
  writeImmersedSolidOutput();

  return;
}


void  HBSplineFEM::writeFluidOutput()
{
  return;
}



void  HBSplineFEM::writeImmersedSolidOutput()
{
  for(int bb=0; bb<nImmSolids; bb++)
    ImmersedBodyObjects[bb]->writeOutput();

  return;
}



/*
void  HBSplineFEM::ImmersedBoundaryBodyForce2D()
{
    // code for 2D Poisson problem

    for(int ii=0;ii<fluidDOF;ii++)
      FluidSolnData.bodyForce(ii) = 0.0;

    int  aa, bb, ii, jj, nlocal, ind, nodenum, size, ee, dir, numPoints;

    double   fact, val, xx, yy, incr1, incr2, *knots[2], param[2], normal[2];
    double   t0, tf, tdiff, funcVal[2], temp, val1, val2, ds, h1, h2;

    vector<int>  bfs;

    node*  ndtmp;
    numPoints = 1;
    numPoints = IBpoints.size();

    nlocal = (degree[0] + 1)*(degree[1]+1);
    VectorXd  N(nlocal), NN;

        for(bb=0;bb<numPoints;bb++)
        {
            //param[0] = 0.5*gridLEN[0];
            //param[1] = 0.5*gridLEN[1];

            //ee = findCellNumber(param[0], param[1]);
            //
            ee = IBpoints[bb]->GetElementNum();


            //printf("xx = %12.6f, yy = %12.6f, cell = %5d, \n", param[0], param[1], ee);

            param[0] = IBpoints[bb]->GetParam(0);
            param[1] = IBpoints[bb]->GetParam(1);

            normal[0] = IBpoints[bb]->GetNormal(0);
            normal[1] = IBpoints[bb]->GetNormal(1);

            ndtmp = elems[ee];

            knots[0] = ndtmp->GetKnots(0);
            knots[1] = ndtmp->GetKnots(1);

            bfs = ndtmp->GlobalBasisFuncs;

            //geometryToParametric(param[0], param[1], param);

            //printf("param ... = %12.6f \t %12.6f, \n \n ", param[0], param[1]);
            //printf("param ... = %12.6f \t %12.6f, \n \n ", knots[0][0], knots[0][1]);
            //printf("param ... = %12.6f \t %12.6f, \n \n ", knots[1][0], knots[1][1]);

            FluidSolnData.computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], param[0], param[1], &N(0));

            if(ndtmp->GetParent() == NULL )
              NN = N;
            else
              NN = ndtmp->SubDivMat * N;

            h1 = gridLEN[0]*knots[0][2];
            h2 = gridLEN[1]*knots[1][2];
            fact = 1.0*IBpoints[bb]->GetArcLength()/h1/h2;

            //val1 = - FluidSolnData.ElemProp.data[4]*fact;
            val1 = fact;
            //val1 = (0.0-soln(totalDOF));
            //val1 = -val2/h;
            //printf("\t val1  %12.8f \n", val1);

            for(ii=0;ii<bfs.size();ii++)
            {
               //FluidSolnData.bodyForce(bfs[ii]) += val1*NN(ii);
               FluidSolnData.bodyForce(3*bfs[ii]) += (val1*NN(ii)*normal[0]);
               FluidSolnData.bodyForce(3*bfs[ii]+1) += (val1*NN(ii)*normal[1]);
            }

            //printf("BasisFuns \n");
            //for(ii=0;ii<bfs.size();ii++)
               //FluidSolnData.bodyForce(bfs[ii]) -= val1*N(ii);
            //printf("\n\n");
        }

    //printf("\n bodyForce \n");
    //printVector(FluidSolnData.bodyForce);

    return;
}
*/



void  HBSplineFEM::ImmersedBoundaryBodyForce3D()
{
    return;
}




void  HBSplineFEM::printResultAtPoint(int index, double u1, double v1, double w1)
{
  return;
}



void  HBSplineFEM::applyInterfaceTerms2D()
{
    ////////////////////////////////////////////////////////
    //
    // stiffness and residual for the immersed boundary points
    // Lagrange multipliers
    //  or
    // Penalty method
    //
    ////////////////////////////////////////////////////////

    int ii, jj, aa, bb, c, r, gp, nr;

    MatrixXd  Klocal, Klocal2;
    VectorXd  Flocal;
    node *nd;

    double PENALTY;
    ImmersedIntegrationElement  *lme;
    
    int nlf=(degree[0]+1)*(degree[0]+1);

    int nlb, ind1, ind2, nGauss;

      VectorXd  NN(nlf), dNN_dx(nlf), dNN_dy(nlf), dN_dx, dN_dy, Nf;
      VectorXd  Flocal2, vel(DIM), vel2(DIM), lagmults(DIM), Nb, dNb, xx, yy,  specValx, specValy, res(DIM);
      MatrixXd  Khorz;
      myPoint  knotIncr, knotBegin;

      double  detJ, *knots[2], af, dvol, fact, fact1, fact2;

      af = SolnData.td(2);
      
      bool axsy = (GeomData.FluidProps[2] == 1);
      double  rho = GeomData.FluidProps[3];
      double  mu = GeomData.FluidProps[4];


      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      {
        nlb = ImmersedBodyObjects[0]->ImmIntgElems[0]->pointNums.size();

        Nb.resize(nlb);
        dNb.resize(nlb);
        xx.resize(nlb);
        yy.resize(nlb);
        specValx.resize(nlb);
        specValy.resize(nlb);

        if(ImmersedBodyObjects[bb]->IsBoundaryConditionTypeLagrange())
        {
          for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
          {
            lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

            for(ii=0;ii<nlb;ii++)
            {
              xx[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][0];
              yy[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][1];
              
              specValx[ii] = lme->GeomDataLag->specValCur[lme->pointNums[ii]][0];
              specValy[ii] = lme->GeomDataLag->specValCur[lme->pointNums[ii]][1];
            }

            //cout << specValx[0] << '\t' << specValy[0] << endl;
            //cout << specValx[1] << '\t' << specValy[1] << endl;

            //cout << " lme->gausspoints.size() = " << lme->gausspoints.size() << endl;

            for(gp=0;gp<lme->gausspoints.size();gp++)
            {
              computeLagrangeBFsLine2D(nlb-1, lme->gausspoints[gp], &xx(0), &yy(0), &Nb(0), &dNb(0), detJ);

              //cout << " detJ " << detJ << '\t' << lme->gaussweights[gp] << '\t' << Nb[0] << endl;

              dvol  = lme->gaussweights[gp] * detJ;

              lagmults.setZero();
              geom.setZero();
              vel2.setZero();

              for(ii=0;ii<nlb;ii++)
              {
                geom[0]     += Nb[ii] * xx[ii];
                geom[1]     += Nb[ii] * yy[ii];
                
                vel2[0]     += Nb[ii] * specValx[ii];
                vel2[1]     += Nb[ii] * specValy[ii];

                lagmults[0] += Nb[ii] * SolnData.var3Cur(lme->pointNums[ii]*DIM);
                lagmults[1] += Nb[ii] * SolnData.var3Cur(lme->pointNums[ii]*DIM+1);
              }

              if(axsy)
              {
                dvol *= 2.0*PI*geom[0];
              }

              //if(geom[0] == 0.0 || geom[0] == 2.0)
                //vel2[0] = 24.0*(0.75-geom[1])*(geom[1]-0.25);

              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);
              //printf("Nb = %12.6f, nlb = %5d \n", Nb[0], nlb);
              //printf("vel2[0] = %12.6f, vel2[1] = %12.6f \n", vel2[0], vel2[1]);
              //printf("lagmults[0] = %12.6f, lagmults[1] = %12.6f \n", lagmults[0], lagmults[1]);
              //printf("detJ = %12.6f, gw = %12.6f \n", detJ, lme->gaussweights[gp]);

              nd = elems[findCellNumber(geom)];

              geometryToParametric(geom, param);

              knots[0] = nd->GetKnots(0);
              knots[1] = nd->GetKnots(1);

              knotBegin = nd->GetKnotBegin();
              knotIncr  = nd->GetKnotIncrement();

              GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

              if(nd->GetParent() == NULL )
              {
                Nf = NN;
                dN_dx = dNN_dx;
                dN_dy = dNN_dy;
              }
              else
              {
                Nf    = nd->SubDivMat * NN;
                dN_dx = nd->SubDivMat * dNN_dx;
                dN_dy = nd->SubDivMat * dNN_dy;
              }

              //printVector(Nf);
              //printf("\n\n");

              ind1 = lme->pointNums.size();
              ind2 = nd->GlobalBasisFuncs.size();

              Khorz.resize(ind1*DIM, ind2*ndof);    Khorz.setZero();

              Flocal.resize(ind2*ndof);   Flocal.setZero();
              Flocal2.resize(ind1*DIM);   Flocal2.setZero();

              vel(0) = vel2(0) - nd->computeValueCur(0, Nf);
              vel(1) = vel2(1) - nd->computeValueCur(1, Nf);

              for(ii=0;ii<ind1;ii++)
              {
                r = DIM*ii;
                fact1 = Nb[ii]*dvol;

                Flocal2(r)   += fact1*vel(0);
                Flocal2(r+1) += fact1*vel(1);

                fact1 = af*fact1;
                //cout << ii << '\t' << r << '\t' << fact1 << endl;

                for(jj=0;jj<ind2;jj++)
                {
                  c = ndof*jj;

                  fact2 = fact1*Nf[jj];
                  //cout << jj << '\t' << c << '\t' << fact2 << endl;
                  Khorz(r,   c)   += fact2;
                  Khorz(r+1, c+1) += fact2;
                }
              }

              for(ii=0;ii<ind2;ii++)
              {
                r = ndof*ii;

                fact = Nf[ii]*dvol;

                Flocal(r)   -= fact*lagmults[0];
                Flocal(r+1) -= fact*lagmults[1];
              }

              //printMatrix(Khorz);
              //printf("\n\n");
              //printVector(Flocal);
              //printf("\n\n");
              //printVector(Flocal3);
              //printf("\n\n");

              ind1 = lme->posIndices.size();
              ind2 = nd->forAssyVec.size();

              //printVector(lme->posIndices);
              //printf("\n\n");
              //printVector(nd->forAssyVec);
              //printf("\n\n");

              for(ii=0;ii<ind2;ii++)
              {
                c = nd->forAssyVec[ii];
                
                solverEigen->rhsVec[c] += Flocal(ii) ;
              }

              for(ii=0;ii<ind1;ii++)
              {
                r = fluidDOF + lme->posIndices[ii];
                
                solverEigen->rhsVec[r] += Flocal2(ii);

                for(jj=0;jj<ind2;jj++)
                {
                  c = nd->forAssyVec[jj];
                  
                  fact = Khorz(ii, jj);

                  solverEigen->mtx.coeffRef(r, c) += fact;
                  solverEigen->mtx.coeffRef(c, r) += fact;
                }
              }

            }//for(gp=0...
          }//for(aa=0...
        }//if(
        else
        {
          myDataIntegrateCutFEM  myData;
          
          PENALTY = ImmersedBodyObjects[bb]->GetPenaltyParameter();

          for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
          {
            lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

            for(ii=0;ii<nlb;ii++)
            {
              xx[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][0];
              yy[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][1];
            }

            for(gp=0;gp<lme->gausspoints.size();gp++)
            {
              computeLagrangeBFsLine2D(nlb-1, lme->gausspoints[gp], &xx(0), &yy(0), &Nb(0), &dNb(0), detJ);

              //cout << " detJ " << detJ << '\t' << gaussweights[gp] << endl;

              dvol  = lme->gaussweights[gp] * detJ;

              vel.setZero();

              geom.setZero();

              for(ii=0;ii<nlb;ii++)
              {
                geom[0] += Nb[ii] * xx[ii];
                geom[1] += Nb[ii] * yy[ii];
              }

              if(axsy)
              {
                dvol *= 2.0*PI*geom[0];
              }

              //cout << " uuuuuuuuuuu " << endl;

              nd = elems[findCellNumber(geom)];

              geometryToParametric(geom, param);

              //cout << " uuuuuuuuuuu " << endl;
              //cout << " dvol " << dvol << endl;

              nr = nd->forAssyVec.size();

              myData.K1 = MatrixXd::Zero(nr, nr);
              myData.F1 = VectorXd::Zero(nr);

              //if(geom[0] == 0.0 || geom[0] == 2.0)
                //vel[0] = 24.0*(0.75-geom[1])*(geom[1]-0.25);

              //if(geom[1] == 1.0 )
                //vel[0] = 1.0;

              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);
              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", param[0], param[1], param[2], dvol);
              //printf("Nb = %12.6f, nlb = %5d, ee = %5d \n", Nb[0], nlb, lme->elemNums[gp]);

              //cout << " uuuuuuuuuuu " << endl;

              myData.geom  = geom;
              myData.param = param;
              myData.dvol  = PENALTY * dvol;

              for(jj=0;jj<DIM;jj++)
              {
                myData.dir = jj;
                myData.specVal[0] = vel[jj];

                //nd->applyBoundaryConditionsAtApoint(jj, param, vel[jj], fact, myData.K1, myData.F1);
                nd->applyBoundaryConditionsAtApoint(myData);
              }

              //solverEigen->AssembleMatrixAndVector(velDOF, 0, nd->forAssyVec, nd->forAssyVec2, Klocal, Flocal);
              solverEigen->AssembleMatrixAndVector(velDOF, 0, nd->forAssyVec, myData.K1, myData.F1);
              //cout << " uuuuuuuuuuu " << endl;
            }//for(gp=0...
          }//for(aa=0...
        }//else
      }//for(bb=0;...
      //

  return;
}




void  HBSplineFEM::applyInterfaceTerms3D()
{
    int ii, jj, aa, bb, c, r, gp, nr;

    MatrixXd  Klocal, Klocal2;
    VectorXd  Flocal;
    node *nd;

    double PENALTY;
    ImmersedIntegrationElement  *lme;
    
    int nlf=(degree[0]+1)*(degree[0]+1)*(degree[0]+1);
    int nlb, ind1, ind2, nGauss;

      VectorXd  NN(nlf), dNN_dx(nlf), dNN_dy(nlf), dNN_dz(nlf), dN_dx, dN_dy, dN_dz, Nf;
      VectorXd  Flocal2, vel(DIM), vel2(DIM), lagmults(DIM), Nb, dNb, res(DIM);
      VectorXd  xx, yy, zz, specValx, specValy, specValz;
      MatrixXd  Khorz;
      myPoint  knotIncr, knotBegin;

      double  detJ, *knots[3], af, dvol, fact, fact1, fact2;

      af = SolnData.td(2);
      double  rho = GeomData.FluidProps[3];
      double  mu = GeomData.FluidProps[4];

      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      {
        nlb = ImmersedBodyObjects[0]->ImmIntgElems[0]->pointNums.size();

        Nb.resize(nlb);
        dNb.resize(nlb);
        xx.resize(nlb);
        yy.resize(nlb);
        zz.resize(nlb);
        specValx.resize(nlb);
        specValy.resize(nlb);
        specValz.resize(nlb);

        if(ImmersedBodyObjects[bb]->IsBoundaryConditionTypeLagrange())
        {
          for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
          {
            lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

            for(ii=0;ii<nlb;ii++)
            {
              xx[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][0];
              yy[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][1];
              zz[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][2];
              
              specValx[ii] = lme->GeomDataLag->specValCur[lme->pointNums[ii]][0];
              specValy[ii] = lme->GeomDataLag->specValCur[lme->pointNums[ii]][1];
              specValz[ii] = lme->GeomDataLag->specValCur[lme->pointNums[ii]][2];
            }

            //cout << " uuuuuuuuuuu " << endl;
            //cout << specValx[0] << '\t' << specValy[0] << endl;
            //cout << specValy[0] << '\t' << specValy[1] << endl;

            for(gp=0;gp<lme->gausspoints.size();gp++)
            {
              computeLagrangeBFsLine3D(nlb-1, lme->gausspoints[gp], &xx(0), &yy(0), &zz(0), &Nb(0), &dNb(0), detJ);

              //cout << " detJ " << detJ << '\t' << lme->gaussweights[gp] << endl;

              dvol  = lme->gaussweights[gp] * detJ;

              lagmults.setZero();
              geom.setZero();
              vel2.setZero();

              for(ii=0;ii<nlb;ii++)
              {
                geom[0]     += Nb[ii] * xx[ii];
                geom[1]     += Nb[ii] * yy[ii];
                geom[2]     += Nb[ii] * zz[ii];
                
                vel2[0]     += Nb[ii] * specValx[ii];
                vel2[1]     += Nb[ii] * specValy[ii];
                vel2[2]     += Nb[ii] * specValz[ii];

                lagmults[0] += Nb[ii] * SolnData.var3Cur(lme->pointNums[ii]*DIM);
                lagmults[1] += Nb[ii] * SolnData.var3Cur(lme->pointNums[ii]*DIM+1);
                lagmults[2] += Nb[ii] * SolnData.var3Cur(lme->pointNums[ii]*DIM+2);
              }

              //if(geom[0] == 0.0 || geom[0] == 2.0)
                //vel2[0] = 24.0*(0.75-geom[1])*(geom[1]-0.25);

              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);
              //printf("Nb = %12.6f, nlb = %5d \n", Nb[0], nlb);
              //printf("lagmults[0] = %12.6f, lagmults[1] = %12.6f \n", lagmults[0], lagmults[1]);
              //printf("detJ = %12.6f, gw = %5d \n", detJ, lme->gaussweights[gp]);

              //lme->computePointAtGP(gp, geom);

              nd = elems[findCellNumber(geom)];

              geometryToParametric(geom, param);

              knots[0] = nd->GetKnots(0);
              knots[1] = nd->GetKnots(1);
              knots[2] = nd->GetKnots(2);

              knotBegin = nd->GetKnotBegin();
              knotIncr  = nd->GetKnotIncrement();

              GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

              if(nd->GetParent() == NULL )
              {
                Nf = NN;
                dN_dx = dNN_dx;
                dN_dy = dNN_dy;
                dN_dz = dNN_dz;
              }
              else
              {
                Nf    = nd->SubDivMat * NN;
                dN_dx = nd->SubDivMat * dNN_dx;
                dN_dy = nd->SubDivMat * dNN_dy;
                dN_dz = nd->SubDivMat * dNN_dz;
              }

              //printVector(Nf);
              //printf("\n\n");

              ind1 = lme->pointNums.size();
              ind2 = nd->GlobalBasisFuncs.size();

              //cout << " ind1 = " << ind1 << endl;
              //cout << " ind2 = " << ind2 << endl;

              Khorz.resize(ind1*DIM, ind2*ndof);    Khorz.setZero();

              Flocal.resize(ind2*ndof);   Flocal.setZero();
              Flocal2.resize(ind1*DIM);   Flocal2.setZero();

              vel(0) = vel2(0) - nd->computeValueCur(0, Nf);
              vel(1) = vel2(1) - nd->computeValueCur(1, Nf);
              vel(2) = vel2(2) - nd->computeValueCur(2, Nf);

              for(ii=0;ii<ind1;ii++)
              {
                r = DIM*ii;
                fact1 = Nb[ii]*dvol;

                Flocal2(r)   += fact1*vel(0);
                Flocal2(r+1) += fact1*vel(1);
                Flocal2(r+2) += fact1*vel(2);

                fact1 = af*fact1;
                //cout << ii << '\t' << r << '\t' << fact1 << endl;

                for(jj=0;jj<ind2;jj++)
                {
                  c = ndof*jj;

                  fact2 = fact1*Nf[jj];
                  //cout << jj << '\t' << c << '\t' << fact2 << endl;
                  Khorz(r,   c)   += fact2;
                  Khorz(r+1, c+1) += fact2;
                  Khorz(r+2, c+2) += fact2;
                }
              }

              for(ii=0;ii<ind2;ii++)
              {
                r = ndof*ii;

                fact = Nf[ii]*dvol;

                Flocal(r)   -= fact*lagmults[0];
                Flocal(r+1) -= fact*lagmults[1];
                Flocal(r+2) -= fact*lagmults[2];
              }

              //cout << " assembling " << endl;
              //printMatrix(Khorz);
              //printf("\n\n");
              //printVector(Flocal);
              //printf("\n\n");
              //printVector(Flocal3);
              //printf("\n\n");

              ind1 = lme->posIndices.size();
              ind2 = nd->forAssyVec.size();

              //printVector(lme->posIndices);
              //printf("\n\n");
              //printVector(nd->forAssyVec);
              //printf("\n\n");

              //cout << " ind1 = " << ind1 << endl;
              //cout << " ind2 = " << ind2 << endl;

              for(ii=0;ii<ind2;ii++)
              {
                c = nd->forAssyVec[ii];
                
                solverEigen->rhsVec[c] += Flocal(ii) ;
              }

              for(ii=0;ii<ind1;ii++)
              {
                r = fluidDOF + lme->posIndices[ii];
                
                solverEigen->rhsVec[r] += Flocal2(ii);

                for(jj=0;jj<ind2;jj++)
                {
                  c = nd->forAssyVec[jj];
                  
                  fact = Khorz(ii, jj);
                  
                  //cout << ii << '\t' << r << '\t' << jj << '\t' << c << endl;

                  solverEigen->mtx.coeffRef(r, c) += fact;
                  solverEigen->mtx.coeffRef(c, r) += fact;
                }
              }

            }//for(gp=0...
          }//for(aa=0...
        }//if(
        else
        {
          myDataIntegrateCutFEM  myData;
          
          PENALTY = ImmersedBodyObjects[bb]->GetPenaltyParameter();

          for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
          {
            lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

            for(ii=0;ii<nlb;ii++)
            {
              xx[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][0];
              yy[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][1];
              zz[ii] = lme->GeomDataLag->NodePosCur[lme->pointNums[ii]][2];
            }

            for(gp=0;gp<lme->gausspoints.size();gp++)
            {
              computeLagrangeBFsLine3D(nlb-1, lme->gausspoints[gp], &xx(0), &yy(0), &zz(0), &Nb(0), &dNb(0), detJ);

              //cout << " detJ " << detJ << '\t' << gaussweights[gp] << endl;

              dvol  = lme->gaussweights[gp] * detJ;

              vel.setZero();

              geom.setZero();

              for(ii=0;ii<nlb;ii++)
              {
                geom[0] += Nb[ii] * xx[ii];
                geom[1] += Nb[ii] * yy[ii];
                geom[2] += Nb[ii] * zz[ii];
              }

              //cout << " uuuuuuuuuuu " << endl;

              nd = elems[findCellNumber(geom)];

              geometryToParametric(geom, param);

              //cout << " uuuuuuuuuuu " << endl;
              //cout << " dvol " << dvol << endl;

              nr = nd->forAssyVec.size();

              myData.K1 = MatrixXd::Zero(nr, nr);
              myData.F1 = VectorXd::Zero(nr);

              //if(geom[0] == 0.0 || geom[0] == 2.0)
                //vel[0] = 24.0*(0.75-geom[1])*(geom[1]-0.25);

              //if(geom[1] == 1.0 )
                //vel[0] = 1.0;

              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);
              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", param[0], param[1], param[2], dvol);
              //printf("Nb = %12.6f, nlb = %5d, ee = %5d \n", Nb[0], nlb, lme->elemNums[gp]);

              //cout << " uuuuuuuuuuu " << endl;

              myData.geom  = geom;
              myData.param = param;
              myData.dvol  = PENALTY * dvol;

              for(jj=0;jj<DIM;jj++)
              {
                myData.dir = jj;
                myData.specVal[0] = vel[jj];

                //nd->applyBoundaryConditionsAtApoint(jj, param, vel[jj], fact, myData.K1, myData.F1);
                nd->applyBoundaryConditionsAtApoint(myData);
              }

              //solverEigen->AssembleMatrixAndVector(velDOF, 0, nd->forAssyVec, nd->forAssyVec2, Klocal, Flocal);
              solverEigen->AssembleMatrixAndVector(velDOF, 0, nd->forAssyVec, myData.K1, myData.F1);
              //cout << " uuuuuuuuuuu " << endl;
            }//for(gp=0...
          }//for(aa=0...
        }//else
      }//for(bb=0;...


  return;
}















