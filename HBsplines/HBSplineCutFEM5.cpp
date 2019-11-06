
#include "HBSplineCutFEM.h"
#include "SolverEigen.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "BasisFunctionsLagrange.h"
#include "myDataIntegrateCutFEM.h"
#include "ImmersedIntegrationElement.h"
#include "QuadratureUtil.h"
#include "typedefs.h"


extern  MpapTime  mpapTime;



void  HBSplineCutFEM::applyInterfaceTerms2D()
{
    //cout << " HBSplineCutFEM::applyInterfaceTerms2D() ... STARTED  " << endl;

    if( ImmersedBodyObjects.size() == 0 )
      return;

    ////////////////////////////////////////////////////////
    //
    // stiffness and residual contributions 
    // from the interface terms
    //
    ////////////////////////////////////////////////////////


      node *nd, *nd2, *nd3, *nd4;

      ImmersedIntegrationElement  *lme;
      myPoly *poly;

      int  aa=0, bb=0, ii=0, jj=0, nlb=0, nlf=0, ind1=0, ind2=0;
      int  nr=0, nc=0, start=0, r=0, c=0, kk=0, gp=0, nGauss=0;
      int  TI=0, TIp1=0, TIp2=0, TJ=0, TJp1=0, TJp2=0, totnlbf2=0, elnum=0;
      int  val1=0, val2=0, nr1=0, nr2=0, start1=0, rp1=0, cp1=0, cp2=0;
      int  *tt1, *tt2;

      nlf = 1;
      for(ii=0; ii<DIM; ii++)
        nlf *= (degree[ii]+1);

      double  fact=0.0, fact1=0.0, fact2=0.0, dvol=0.0, detJ=0.0;
      double  PENALTY1=0.0, PENALTY2=0.0;
      double  pres=0.0, presPrev=0.0, x0=0.0, y0=0.0, x1=0.0, y1=0.0, rad=0.0;
      double  bb1=0.0, bb2=0.0, NitscheFact=0.0, c1=0.0, hx=0.0, hy=0.0;
      double  Ta[6], Tb[6], FlowRate=6.545*0.5;
      bool  isNitsche=true;

      AABB  bbTemp;

      double  af = SolnData.td(2);
      bool axsy = (fluidProps[2] == 1);
      double  rho = fluidProps[3];
      double  mu = fluidProps[4];

      MatrixXd  K1, grad(2,2), stress(2,2), Kc(2,3);
      VectorXd  F1, Nb;
      myPoint   vel, trac, velSpec;
      VectorXd  NN(nlf), dNN_dx(nlf), dNN_dy(nlf), N, dN_dx, dN_dy;

      myPoint  knotBegin, knotIncr;

      vector<double>  gausspoints, gaussweights;

      param.setZero();
      geom.setZero();
      normal.setZero();

      start1 = fluidDOF;
      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      {
        PENALTY1    = ImmersedBodyObjects[bb]->getPenaltyParameter();
        isNitsche   = ImmersedBodyObjects[bb]->getNitscheFlag();
        NitscheFact = ImmersedBodyObjects[bb]->getNitscheFact();

        nlb = ImmersedBodyObjects[bb]->ImmIntgElems[0]->pointNums.size();

        Nb.resize(nlb);

        nGauss = (int) cutFEMparams[1];

        getGaussPoints1D(nGauss, gausspoints, gaussweights);

        for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
        {
          lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];
          poly = ImmersedBodyObjects[bb]->ImmersedFaces[aa];

          val1 = lme->pointNums.size();
          tt1  =  &(lme->pointNums[0]);

          if( lme->isActive() )
          {
            for(gp=0;gp<nGauss;gp++)
            {
              param[0] = gausspoints[gp];
              poly->computeBasisFunctions(param, geom, Nb, detJ);

              dvol  = gaussweights[gp] * detJ;

              if(axsy)
              {
                dvol *= 2.0*PI*geom[0];
              }

              lme->computeVelocityCur(Nb, velSpec);

              // the normal from the surface of the immersed body should be away from the fluid
              // the default behaviour is assumed to be that 
              // the immersed polygon is oriented clockwise

              poly->computeNormal(geom, normal);

              //printf("normal     ... %12.6f \t %12.6f \t %12.6f \n", normal[0], normal[1], normal[2]);
              //printf("velSpec    ... %12.6f \t %12.6f \t %12.6f \n", velSpec[0], velSpec[1], velSpec[2]);
              //printf("geom = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);

              elnum = findCellNumber(geom);

              nd2 = elems[elnum];

              if( (nd2->getSubdomainId() == this_mpi_proc) && (nd2->isCutElement() || nd2->domNums[0] == 0) )
              {
                bbTemp = nd2->getAABB();
                hx = bbTemp.maxBB[0]-bbTemp.minBB[0];
                hy = bbTemp.maxBB[1]-bbTemp.minBB[1];

                // only if the element is fluid or cutcell
                if( (nd2->domNums.size() == 1) && (nd2->domNums[0] != 0) )
                {
                  if( abs(geom[1]-bbTemp.minBB[1]) < 1.0e-8 ) // bottom edge
                  {
                    nd3 = nd2->getNeighbour(EAST);

                    if( (nd3->domNums.size() == 1) && (nd3->domNums[0] != 0) )
                    {
                      nd4 = nd2->getNeighbour(SOUTH);
                      if( (nd4->domNums.size() == 1) && (nd4->domNums[0] != 0) )
                        nd = nd2->getNeighbour(WEST);
                      else
                        nd = nd4;
                    }
                    else
                      nd = nd2->getNeighbour(SOUTH);
                  }
                  else if( abs(geom[1]-bbTemp.maxBB[1]) < 1.0e-8 ) // top edge
                  {
                    nd3 = nd2->getNeighbour(EAST);

                    if( (nd3->domNums.size() == 1) && (nd3->domNums[0] != 0) )
                    {
                      nd4 = nd2->getNeighbour(NORTH);
                      if( (nd4->domNums.size() == 1) && (nd4->domNums[0] != 0) )
                        nd = nd2->getNeighbour(WEST);
                      else
                        nd = nd4;
                    }
                    else
                      nd = nd2->getNeighbour(NORTH);
                  }
                  else
                  {
                    if( abs(geom[0]-bbTemp.minBB[0]) < 1.0e-8 ) // left edge
                      nd = nd2->getNeighbour(WEST);
                    else if( abs(geom[0]-bbTemp.maxBB[0]) < 1.0e-8 ) // right edge
                      nd = nd2->getNeighbour(EAST);
                    else
                    {
                      cerr << " Boundary integration point lies in absurd element \n\n " << endl;
                    }
                  }
                }
                else
                  nd = nd2;
                //

                geometryToParametric(geom, param);

                PENALTY2 = PENALTY1;
                // to avoid singularities at the corners
                //if( nd->isBackBoundary() )
                  //PENALTY2 = PENALTY1*10000.0;

                nr = nd->forAssyVec.size();

                K1 = MatrixXd::Zero(nr, nr);
                F1 = VectorXd::Zero(nr);

                //printf("param = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", param[0], param[1], param[2], dvol);
                //printf("vx = %12.6f, vy = %12.6f, vz = %12.6f \n", velSpec[0], velSpec[1], velSpec[2]);

                // stiffness and residual terms due to Nitsche method

                knotBegin = nd->getKnotBegin();
                knotIncr  = nd->getKnotIncrement();

                GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

                if(nd->getParent() == NULL)
                {
                  N = NN;
                  dN_dx = dNN_dx;
                  dN_dy = dNN_dy;
                }
                else
                {
                  N = nd->SubDivMat*NN;
                  dN_dx = nd->SubDivMat*dNN_dx;
                  dN_dy = nd->SubDivMat*dNN_dy;
                }

                vel[0] = velSpec[0] - nd->computeValueCur(0, N);
                vel[1] = velSpec[1] - nd->computeValueCur(1, N);

                //printf("vel        ... %12.6f \t %12.6f \t %12.6f \n", vel[0], vel[1], vel[2]);

                /////////////////////////////////////////
                // compute the stresses and tractions
                /////////////////////////////////////////

                grad(0,0) = nd->computeValueCur(0, dN_dx);
                grad(0,1) = nd->computeValueCur(0, dN_dy);
                grad(1,0) = nd->computeValueCur(1, dN_dx);
                grad(1,1) = nd->computeValueCur(1, dN_dy);
                pres      = nd->computeValue(2, N);
                presPrev  = nd->computeValuePrev(2, N);

                //stress = mu*(grad+grad.transpose());
                stress = mu*grad;

                stress(0,0) -= pres;
                stress(1,1) -= pres;

                trac[0] = stress(0,0)*normal[0] + stress(0,1)*normal[1] ;
                trac[1] = stress(1,0)*normal[0] + stress(1,1)*normal[1] ;

                totnlbf2 = nr/ndof;

                for(ii=0;ii<totnlbf2;ii++)
                {
                  TI   = 3*ii;
                  TIp1 = TI+1;
                  TIp2 = TI+2;

                  bb1 = N[ii] * dvol;
                  bb2 = bb1 * PENALTY2;

                  Ta[0] = (dvol*mu)*( dN_dx[ii]*normal[0] + dN_dy[ii]*normal[1] );
                  Ta[1] = 0.0;
                  Ta[2] = -normal[0]*bb1;

                  Ta[3] = 0.0;
                  Ta[4] = Ta[0];
                  Ta[5] = -normal[1]*bb1;

                  for(jj=0; jj<6; jj++)
                    Ta[jj] *= NitscheFact;

                  for(jj=0;jj<totnlbf2;jj++)
                  {
                    TJ   = 3*jj;
                    TJp1 = TJ+1;
                    TJp2 = TJ+2;

                    fact = bb2 * af *  N[jj];
                    // stabilisation term
                    K1(TI,   TJ)   += fact;
                    K1(TIp1, TJp1) += fact;

                    // Nitsche terms

                    Tb[0] = (af*mu)*( dN_dx[jj]*normal[0] + dN_dy[jj]*normal[1] );
                    Tb[1] = 0.0;
                    Tb[2] = -normal[0]*N[jj];

                    Tb[3] = 0.0;
                    Tb[4] = Tb[0];
                    Tb[5] = -normal[1]*N[jj];

                    K1(TI, TJ)      -= (bb1*Tb[0]);
                    K1(TI, TJp1)    -= (bb1*Tb[1]);
                    K1(TI, TJp2)    -= (bb1*Tb[2]);

                    K1(TIp1, TJ)    -= (bb1*Tb[3]);
                    K1(TIp1, TJp1)  -= (bb1*Tb[4]);
                    K1(TIp1, TJp2)  -= (bb1*Tb[5]);

                    c1 = af*N[jj];

                    K1(TI,   TJ)    -= (Ta[0]*c1);
                    K1(TIp1, TJ)    -= (Ta[1]*c1);
                    K1(TIp2, TJ)    -= (Ta[2]*c1);

                    K1(TI,   TJp1)  -= (Ta[3]*c1);
                    K1(TIp1, TJp1)  -= (Ta[4]*c1);
                    K1(TIp2, TJp1)  -= (Ta[5]*c1);
                  }

                  // stabilisation terms
                  F1(TI)   += (bb2*vel[0]);
                  F1(TIp1) += (bb2*vel[1]);

                  // Nitsche terms

                  F1(TI)   += (bb1*trac[0]);
                  F1(TIp1) += (bb1*trac[1]);

                  F1(TI)   -= (Ta[0]*vel[0]);
                  F1(TIp1) -= (Ta[1]*vel[0]);
                  F1(TIp2) -= (Ta[2]*vel[0]);

                  F1(TI)   -= (Ta[3]*vel[1]);
                  F1(TIp1) -= (Ta[4]*vel[1]);
                  F1(TIp2) -= (Ta[5]*vel[1]);

                } // for(ii=0;ii<totnlbf2;ii++)

                solverPetsc->assembleMatrixAndVectorCutFEM(0, 0, nd->forAssyVec, grid_to_proc_DOF, K1, F1);
              } //if( nd2->getSubdomainId() == this_mpi_proc )
            }//for(gp=0...
          } // if( lme->isActive() )
        }//for(aa=0...
        start1 += ImmersedBodyObjects[bb]->getTotalDOF();
      }//for(bb=0;...

  //cout << " HBSplineCutFEM::applyInterfaceTerms2D() ... FINISHED  " << endl;

  if(cutFEMparams[6] > 1.0e-8)
    applyGhostPenalty2D();

  return;
}
//



/*
              if(!STAGGERED)
              {
                  //trac[0] =  stress(0,0)*normal[0] + stress(0,1)*normal[1] ;
                  //trac[1] =  stress(1,0)*normal[0] + stress(1,1)*normal[1] ;

                  //trac[0] += PENALTY*(vel[0]-velSpec[0]);
                  //trac[1] += PENALTY*(vel[1]-velSpec[1]);

                  //trac[0] -= PENALTY*(vel[0]);
                  //trac[1] -= PENALTY*(vel[1]);

                  nr1  =  lme->pointNums.size();
                  tt1  =  &(lme->pointNums[0]);

                  nr2  =  nd->GlobalBasisFuncs.size();
                  tt2  =  &(nd->GlobalBasisFuncs[0]);
                  //nr1 = nd1->forAssyVec.size();

                  for(ii=0; ii<nr1; ii++)
                  {
                    r   = start1 + tt1[ii]*2;
                    rp1 = r+1;
                    
                    bb1 = Nb[ii]*dvol;
                    bb2 = bb1 * PENALTY;

                    fact1 = -bb1*trac[0];
                    fact2 = -bb1*trac[1];

                    //VecSetValue(solverPetsc->rhsVec, r,    fact1, ADD_VALUES);
                    //VecSetValue(solverPetsc->rhsVec, rp1,  fact2, ADD_VALUES);

                    for(jj=0; jj<nr2; jj++)
                    {
                      c  = grid_to_proc_BF[tt2[jj]]*ndof;
                      cp1 = c+1;
                      cp2 = c+2;

                      fact = -af*bb2*N[jj];

                      //Ksf
                      MatSetValue(solverPetsc->mtx, r,   c,   fact, ADD_VALUES);
                      MatSetValue(solverPetsc->mtx, rp1, cp1, fact, ADD_VALUES);

                      //Kfs
                      MatSetValue(solverPetsc->mtx, c,   r,   fact, ADD_VALUES);
                      MatSetValue(solverPetsc->mtx, cp1, rp1, fact, ADD_VALUES);

                      Tb[0] = (af*mu)*( dN_dx[jj]*normal[0] + dN_dy[jj]*normal[1] );
                      Tb[1] = 0.0;
                      Tb[2] = -normal[0]*N[jj];

                      Tb[3] = 0.0;
                      Tb[4] = Tb[0];
                      Tb[5] = -normal[1]*N[jj];

                      Kc(0, 0)  = bb1*Tb[0];
                      Kc(0, 1)  = bb1*Tb[1];
                      Kc(0, 2)  = bb1*Tb[2];

                      Kc(1, 0)  = bb1*Tb[3];
                      Kc(1, 1)  = bb1*Tb[4];
                      Kc(1, 2)  = bb1*Tb[5];

                      //Ksf
                      MatSetValue(solverPetsc->mtx, r,     c,    Kc(0,0), ADD_VALUES);
                      MatSetValue(solverPetsc->mtx, r,     cp1,  Kc(0,1), ADD_VALUES);
                      MatSetValue(solverPetsc->mtx, r,     cp2,  Kc(0,2), ADD_VALUES);

                      MatSetValue(solverPetsc->mtx, rp1,   c,    Kc(1,0), ADD_VALUES);
                      MatSetValue(solverPetsc->mtx, rp1,   cp1,  Kc(1,1), ADD_VALUES);
                      MatSetValue(solverPetsc->mtx, rp1,   cp2,  Kc(1,2), ADD_VALUES);

                      Kc = Kc*NitscheFact;
                      //Kfs
                      MatSetValue(solverPetsc->mtx, c,     r,    Kc(0,0), ADD_VALUES);
                      MatSetValue(solverPetsc->mtx, cp1,   r,    Kc(0,1), ADD_VALUES);
                      MatSetValue(solverPetsc->mtx, cp2,   r,    Kc(0,2), ADD_VALUES);

                      MatSetValue(solverPetsc->mtx, c,     rp1,  Kc(1,0), ADD_VALUES);
                      MatSetValue(solverPetsc->mtx, cp1,   rp1,  Kc(1,1), ADD_VALUES);
                      MatSetValue(solverPetsc->mtx, cp2,   rp1,  Kc(1,2), ADD_VALUES);

                    } // for(jj=0; jj<nr2; jj++)
                  } // for(ii=0; ii<nr1; ii++)
              } //if(!STAGGERED)

              if(bb==-222)
              {
                if(CompareDoubles(geom[0],-1.1))
                {
                  y0 = 0.9; y1 = 1.5;
                  //if(geom[1]>=y0 && geom[1]<=y1)
                  //{
                    velSpec[0] = 20*(y1-geom[1])*(geom[1]-y0);
                    velSpec[1] = 0.0;
                  //}
                }
              }

              if(bb==-31)
              {
                x0 = 36.0;  y0 = 28.5;
                
                x1= geom[0]-x0;
                y1= geom[1]-y0;
                
                rad = sqrt(x1*x1 + y1*y1);

                velSpec[0] = - FlowRate*x1/(2.0*PI*rad*rad);
                velSpec[1] = - FlowRate*y1/(2.0*PI*rad*rad);
              }
              if(bb== -11)
              {
                x0 = -84.0;  y0 = 28.5;
                
                x1= geom[0]-x0;
                y1= geom[1]-y0;
                
                rad = sqrt(x1*x1 + y1*y1);

                velSpec[0] = - FlowRate*x1/(2.0*PI*rad*rad);
                velSpec[1] = - FlowRate*y1/(2.0*PI*rad*rad);
              }
*/


void  HBSplineCutFEM::applyInterfaceTerms3D()
{
    if( ImmersedBodyObjects.size() == 0 )
      return;

    ////////////////////////////////////////////////////////
    //
    // stiffness and residual contributions 
    // from the interface terms
    //
    ////////////////////////////////////////////////////////

    node *ndTemp;

    int  aa=0, bb=0, ii=0, jj=0, nlb=0, ind1=0, ind2=0, nr=0, nc=0;
    int  start=0, r=0, c=0, kk=0, gp=0, nGauss=0;
    int  TI=0, TIp1=0, TIp2=0, TIp3=0, TJ=0, TJp1=0, TJp2=0, TJp3=0, totnlbf2=0;

    int  nlf = (degree[0]+1)*(degree[1]+1)*(degree[2]+1);

    double  fact=0.0, fact1=0.0, fact2=0.0, fact3=0.0, dvol=0.0, PENALTY=0.0;
    double  detJ=0.0, pres=0.0, bb1=0.0, bb2=0.0, NitscheFact=0.0, c1=0.0;
    double  Ta[12], Tb[12], surfArea=0.0;

    double  af  = SolnData.td(2);
    double  rho = fluidProps[3];
    double  mu  = fluidProps[4];
    bool  isNitsche=true;

    MatrixXd  K1;
    VectorXd  F1, Nb;
    VectorXd  NN(nlf), dNN_dx(nlf), dNN_dy(nlf), dNN_dz(nlf), N, dN_dx, dN_dy, dN_dz;
    myPoint   vel, velSpec, trac;
    MatrixXd  stress(3,3);

    myPoint   knotBegin, knotIncr;

    vector<double>  gausspoints1, gausspoints2, gaussweights;

    param.setZero();
    geom.setZero();
    normal.setZero();

    ImmersedIntegrationElement  *lme;
    myPoly*  poly;

    for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
    {
        PENALTY     = ImmersedBodyObjects[bb]->getPenaltyParameter();
        isNitsche   = ImmersedBodyObjects[bb]->getNitscheFlag();
        NitscheFact = ImmersedBodyObjects[bb]->getNitscheFact();

        nlb = ImmersedBodyObjects[bb]->ImmIntgElems[0]->pointNums.size();

        Nb.resize(nlb);

        nGauss = (int) cutFEMparams[1];

        if(nlb == 3)
          getGaussPointsTriangle(nGauss, gausspoints1, gausspoints2, gaussweights);
        else
          getGaussPointsQuad(nGauss, gausspoints1, gausspoints2, gaussweights);

        surfArea = 0.0;

        for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
        {
          lme  = ImmersedBodyObjects[bb]->ImmIntgElems[aa];
          poly = ImmersedBodyObjects[bb]->ImmersedFaces[aa];

          if( lme->isActive() )
          {
            for(gp=0;gp<nGauss;gp++)
            {
              param[0] = gausspoints1[gp];
              param[1] = gausspoints2[gp];

              poly->computeBasisFunctions(param, geom, Nb, detJ);

              dvol  = gaussweights[gp] * detJ;

              poly->computeNormal();
              poly->computeNormal(param, normal);

              lme->computeVelocityCur(Nb, velSpec );

              ndTemp = elems[findCellNumber(geom)];

              if( (ndTemp->getSubdomainId() == this_mpi_proc) && (ndTemp->isCutElement() || ndTemp->domNums[0] == 0) )
              {
                geometryToParametric(geom, param);

                nr = ndTemp->forAssyVec.size();

                K1 = MatrixXd::Zero(nr, nr);
                F1 = VectorXd::Zero(nr);

                //printf("geometry   ... %12.6f \t %12.6f \t %12.6f \n", geom[0], geom[1], geom[2]);
                //printf("parameters ... %12.6f \t %12.6f \t %12.6f \n", param[0], param[1], param[2]);
                //printf("velocity   ... %12.6f \t %12.6f \t %12.6f \n", velSpec[0], velSpec[1], velSpec[2]);
                //printf("normal     ... %12.6f \t %12.6f \t %12.6f \n", normal[0], normal[1], normal[2]);
                //printf("Jacobian   ... %12.6f \t %12.6f \n", detJ, dvol);

                // stiffness and residual terms due to Nitsche method

                knotBegin = ndTemp->getKnotBegin();
                knotIncr  = ndTemp->getKnotIncrement();

                GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

                if(ndTemp->getParent() == NULL)
                {
                  N = NN;
                  dN_dx = dNN_dx;
                  dN_dy = dNN_dy;
                  dN_dz = dNN_dz;
                }
                else
                {
                  N = ndTemp->SubDivMat*NN;
                  dN_dx = ndTemp->SubDivMat*dNN_dx;
                  dN_dy = ndTemp->SubDivMat*dNN_dy;
                  dN_dz = ndTemp->SubDivMat*dNN_dz;
                }

                //cout << " uuuuuuuuuuu " << endl;

                vel[0] = velSpec[0] - ndTemp->computeValueCur(0, N);
                vel[1] = velSpec[1] - ndTemp->computeValueCur(1, N);
                vel[2] = velSpec[2] - ndTemp->computeValueCur(2, N);

                ///////////////////////////////////
                // compute the stresses and tractions
                ///////////////////////////////////

                //nd->computeVelocityAndStressCur(param, vel1, stress1);

                stress(0,0) = ndTemp->computeValueCur(0, dN_dx);
                stress(0,1) = ndTemp->computeValueCur(0, dN_dy);
                stress(0,2) = ndTemp->computeValueCur(0, dN_dz);

                stress(1,0) = ndTemp->computeValueCur(1, dN_dx);
                stress(1,1) = ndTemp->computeValueCur(1, dN_dy);
                stress(1,2) = ndTemp->computeValueCur(1, dN_dz);

                stress(2,0) = ndTemp->computeValueCur(2, dN_dx);
                stress(2,1) = ndTemp->computeValueCur(2, dN_dy);
                stress(2,2) = ndTemp->computeValueCur(2, dN_dz);

                pres        = ndTemp->computeValue(3, N);

                stress  = mu*stress;

                stress(0,0) -= pres;
                stress(1,1) -= pres;
                stress(2,2) -= pres;

                trac[0] = stress(0,0)*normal[0] + stress(0,1)*normal[1] + stress(0,2)*normal[2] ;
                trac[1] = stress(1,0)*normal[0] + stress(1,1)*normal[1] + stress(1,2)*normal[2] ;
                trac[2] = stress(2,0)*normal[0] + stress(2,1)*normal[1] + stress(2,2)*normal[2] ;

                totnlbf2 = nr/ndof;

                for(ii=0;ii<totnlbf2;ii++)
                {
                  TI   = 4*ii;
                  TIp1 = TI+1;
                  TIp2 = TI+2;
                  TIp3 = TI+3;

                  bb1 = N[ii] * dvol;
                  bb2 = bb1 * PENALTY ;

                  Ta[0]  =  (dvol*mu)*( dN_dx[ii]*normal[0] + dN_dy[ii]*normal[1] + dN_dz[ii]*normal[2] );
                  Ta[1]  =  0.0;
                  Ta[2]  =  0.0;
                  Ta[3]  =  -normal[0]*bb1;

                  Ta[4]  =  0.0;
                  Ta[5]  =  Ta[0];
                  Ta[6]  =  0.0;
                  Ta[7]  =  -normal[1]*bb1;

                  Ta[8]  =  0.0;
                  Ta[9]  =  0.0;
                  Ta[10] =  Ta[0];
                  Ta[11] =  -normal[2]*bb1;

                  for(jj=0; jj<12; jj++)
                    Ta[jj] *= NitscheFact;

                  for(jj=0;jj<totnlbf2;jj++)
                  {
                    TJ   = 4*jj;
                    TJp1 = TJ+1;
                    TJp2 = TJ+2;
                    TJp3 = TJ+3;

                    fact = af* bb2 * N[jj];
                    // stabilisation term
                    K1(TI,   TJ)   += fact;
                    K1(TIp1, TJp1) += fact;
                    K1(TIp2, TJp2) += fact;

                    // Nitsche terms
                    Tb[0]  =  af*mu*( dN_dx[jj]*normal[0] + dN_dy[jj]*normal[1] + dN_dz[jj]*normal[2] );
                    Tb[1]  =  0.0;
                    Tb[2]  =  0.0;
                    Tb[3]  =  -normal[0]*N[jj];

                    Tb[4]  =  0.0;
                    Tb[5]  =  Tb[0];
                    Tb[6]  =  0.0;
                    Tb[7]  =  -normal[1]*N[jj];

                    Tb[8]  =  0.0;
                    Tb[9]  =  0.0;
                    Tb[10] =  Tb[0];
                    Tb[11] =  -normal[2]*N[jj];


                    K1(TI,   TJ)    -= (bb1*Tb[0]);
                    K1(TI,   TJp1)  -= (bb1*Tb[1]);
                    K1(TI,   TJp2)  -= (bb1*Tb[2]);
                    K1(TI,   TJp3)  -= (bb1*Tb[3]);

                    K1(TIp1, TJ)    -= (bb1*Tb[4]);
                    K1(TIp1, TJp1)  -= (bb1*Tb[5]);
                    K1(TIp1, TJp2)  -= (bb1*Tb[6]);
                    K1(TIp1, TJp3)  -= (bb1*Tb[7]);

                    K1(TIp2, TJ)    -= (bb1*Tb[8]);
                    K1(TIp2, TJp1)  -= (bb1*Tb[9]);
                    K1(TIp2, TJp2)  -= (bb1*Tb[10]);
                    K1(TIp2, TJp3)  -= (bb1*Tb[11]);

                    c1 = af*N[jj];

                    K1(TI,   TJ)    -= (Ta[0]*c1);
                    K1(TIp1, TJ)    -= (Ta[1]*c1);
                    K1(TIp2, TJ)    -= (Ta[2]*c1);
                    K1(TIp3, TJ)    -= (Ta[3]*c1);

                    K1(TI,   TJp1)  -= (Ta[4]*c1);
                    K1(TIp1, TJp1)  -= (Ta[5]*c1);
                    K1(TIp2, TJp1)  -= (Ta[6]*c1);
                    K1(TIp3, TJp1)  -= (Ta[7]*c1);

                    K1(TI,   TJp2)  -= (Ta[8]*c1);
                    K1(TIp1, TJp2)  -= (Ta[9]*c1);
                    K1(TIp2, TJp2)  -= (Ta[10]*c1);
                    K1(TIp3, TJp2)  -= (Ta[11]*c1);
                  }

                  // stabilisation terms

                  F1(TI)   += (bb2*vel[0]);
                  F1(TIp1) += (bb2*vel[1]);
                  F1(TIp2) += (bb2*vel[2]);

                  // Nitsche terms

                  F1(TI)   += (bb1*trac[0]);
                  F1(TIp1) += (bb1*trac[1]);
                  F1(TIp2) += (bb1*trac[2]);

                  F1(TI)   -= (Ta[0]*vel[0]);
                  F1(TIp1) -= (Ta[1]*vel[0]);
                  F1(TIp2) -= (Ta[2]*vel[0]);
                  F1(TIp3) -= (Ta[3]*vel[0]);

                  F1(TI)   -= (Ta[4]*vel[1]);
                  F1(TIp1) -= (Ta[5]*vel[1]);
                  F1(TIp2) -= (Ta[6]*vel[1]);
                  F1(TIp3) -= (Ta[7]*vel[1]);

                  F1(TI)   -= (Ta[8]*vel[2]);
                  F1(TIp1) -= (Ta[9]*vel[2]);
                  F1(TIp2) -= (Ta[10]*vel[2]);
                  F1(TIp3) -= (Ta[11]*vel[2]);

                } // for(ii=0;ii<totnlbf2;ii++)

                //cout << " uuuuuuuuuuu \n\n " << endl;

                solverPetsc->assembleMatrixAndVectorCutFEM(0, 0, ndTemp->forAssyVec, grid_to_proc_DOF, K1, F1);
              } //if( ndTemp->getSubdomainId() == this_mpi_proc )
            }//for(gp=0...
          } // if( lme->isActive() )
        }//for(aa=0...
      }//for(bb=0;...
      //
      //printf("Surface area = %12.6f, \n", surfArea);

  if(abs(cutFEMparams[6]) > 1.0e-8)
    applyGhostPenalty3D();

  return;
}



void  HBSplineCutFEM::applyGhostPenalty2D()
{
    //cout << " HBSplineCutFEM::applyGhostPenalty2D() ... STARTED  " << endl;

    ////////////////////////////////////////////////////////
    //
    // stiffness and residual due to ghost-penalty terms
    // for 2D problem
    ////////////////////////////////////////////////////////

    int  ii=0, jj=0, ee=0, TI=0, TJ=0, gp1=0, gp2=0, side=0, dir=0, NUM_NEIGHBOURS=4;
    int  nlbf1=0, nlbf2=0, nr1=0, nr2=0, domTemp=0;
    int  r=0, c=0, start1=0, start2=0, size1=0, size2=0;

    double  Ta1=0.0, Ta2=0.0, Tb1=0.0, Tb2=0.0, JacMultLoc=0.0, rad=0.0, dvol=0.0, fact=0.0;
    double  bb1=0.0, bb2=0.0, PENALTY=0.0, h1=0.0, h2=0.0;
    vector<double>  boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2;
    vector<double>  gammGP(ndof+1), temp(ndof+1);

    std::vector<int>::iterator itInt;

    int nlf = (degree[0]+1) * (degree[1]+1);

    VectorXd   NN1(nlf), dNN1_dx(nlf), dNN1_dy(nlf), N1, dN1_dx, dN1_dy;
    VectorXd   NN2(nlf), dNN2_dx(nlf), dNN2_dy(nlf), N2, dN2_dx, dN2_dy;
    myPoint   normal1, normal2, hElem;
    myPoint   knotIncr1, knotBegin1, knotEnd1, knotSum1, knotIncr2, knotBegin2, knotEnd2, knotSum2;

    node *nd1, *nd2;

    VectorXd   trac1(ndof), trac2(ndof), res(ndof);
    MatrixXd   grad1(ndof, ndof), grad2(ndof, ndof);

    double  af  = SolnData.td(2);
    double  rho = fluidProps[3];
    double  mu  = fluidProps[4];
    bool   axsy = ((int)fluidProps[2] == 1);

    // find the maximum h
    // maximum h is the maximum of diameters of all the cut elements
    // where, diameter of a cut element is taken to be the diagonal of that element

    ii = elems[cutCellIds[0]]->getLevel();

    for(ee=0; ee<cutCellIds.size(); ee++)
    {
      ii = max( ii, elems[cutCellIds[ee]]->getLevel() );
    }

    // cout << " Maximum cut-cell level = " << ii << endl;
    // lengths of element sides at level '0'

    hElem[0] = GeomData.getGridLength(0)/nelem[0];
    hElem[1] = GeomData.getGridLength(1)/nelem[1];

    // lengths of element sides at level 'ii'
    fact = pow(2.0, double (ii) );

    hElem[0] /= fact;
    hElem[1] /= fact;

    h1 = sqrt( hElem[0]*hElem[0] + hElem[1]*hElem[1] );


    for(ee=0; ee<cutCellIds.size(); ee++)
    {
      nd1 = elems[cutCellIds[ee]];

      //if( nd1->isCutElement() && !(nd1->isBoundary()) )
      if( nd1->getSubdomainId() == this_mpi_proc )
      {
        if( nd1->isCutElement() )
        {
          //cout << " nd1->getID() " <<  nd1->getID() << '\t' <<  nd1->getLevel() << endl;

          knotBegin1 = nd1->getKnotBegin();
          knotEnd1   = nd1->getKnotEnd();
          knotIncr1  = nd1->getKnotIncrement();
          knotSum1   = nd1->getKnotSum();

          //cout << " nd1->isRightBoundary() = " << nd1->isRightBoundary() << endl;

          if( nd1->isBoundary() )
          {
            //gammGP[0]  = 0.0;
            //gammGP[0]  = cutFEMparams[6] *mu* h1*h1*h1;
            gammGP[0]  = cutFEMparams[6] *mu*h1;
            gammGP[1]  = gammGP[0];
           }
          else
          {
            //gammGP[0]  = cutFEMparams[6] *mu* h1*h1*h1;
            gammGP[0]  = cutFEMparams[6] * mu*h1;
            gammGP[1]  = gammGP[0];
          }

          gammGP[2]  = cutFEMparams[7] * h1*h1*h1/mu;

          bb1 = 1.0;
          bb2 = 1.0;
          for(ii=1; ii<degree[0]; ii++)
          {
            bb1 *= h1*h1;
          }
          gammGP[0] *= bb1;
          gammGP[1] *= bb1;
          gammGP[2] *= bb1;

          for(side=0; side<NUM_NEIGHBOURS; side++)
          {
            nd2 = nd1->getNeighbour(side);

            //cout << " side = " << side << '\t' << nd2->getDomainNumber() << endl;
            //if( !(nd2 == NULL) && !(nd2->isGhost()) && nd2->isLeaf() && (nd2->getDomainNumber() == -1) && (nd2->getDomainNumber() >= 1) )
            //if( (nd2 != NULL) && !(nd2->isGhost()) && !(nd2->isBoundary()) && nd2->isLeaf() && (nd2->isCutElement() || nd2->domNums[0] == 0) )
            if( (nd2 != NULL) && !(nd2->isGhost()) && nd2->isLeaf() && (nd2->isCutElement() || nd2->domNums[0] == 0) )
            {
              nr1 = nd1->forAssyVec.size();
              nr2 = nd2->forAssyVec.size();

              nlbf1 = nr1/ndof;
              nlbf2 = nr2/ndof;

              MatrixXd  K1, K2, Kc;
              VectorXd  F1, F2;

              K1 = MatrixXd::Zero(nr1, nr1);
              K2 = MatrixXd::Zero(nr2, nr2);
              Kc = MatrixXd::Zero(nr1, nr2);
              F1 = VectorXd::Zero(nr1);
              F2 = VectorXd::Zero(nr2);

              knotBegin2 = nd2->getKnotBegin();
              knotEnd2   = nd2->getKnotEnd();
              knotIncr2  = nd2->getKnotIncrement();
              knotSum2   = nd2->getKnotSum();

              GeomData.getBoundaryNormal2D(side, normal1);
              GeomData.setBoundaryGPs2D(side, boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2);

              normal2 = -normal1;

              JacMultLoc = nd1->getJacBoundary(side);

              for(gp2=0;gp2<boundaryGPs2.size();gp2++)
              {
                  param[1] = 0.5 * (knotIncr1[1] * boundaryGPs2[gp2] + knotSum1[1]);
              for(gp1=0;gp1<boundaryGPs1.size();gp1++)
              {
                  param[0] = 0.5 * (knotIncr1[0] * boundaryGPs1[gp1] + knotSum1[0]);
 
                  dvol = JacMultLoc * boundaryGWs2[gp2] * boundaryGWs1[gp1] ;

                  // multiply dvol with 0.5 as the ghost-penalty operation is applied
                  // twice on the face shared by the two elements both of which are cut elements

                  if( nd2->getDomainNumber() == -1 )
                    dvol *= 0.5;

                  //printf(" %4d \t %4d \t %12.6f \t %12.6f \n", side, dir, param[0], param[1]);

                  GeomData.computeBasisFunctionsGhostPenalty2D(knotBegin1, knotIncr1, param, NN1, dNN1_dx, dNN1_dy );

                  if(nd1->getParent() == NULL)
                  {
                    N1 = NN1;
                    dN1_dx = dNN1_dx;
                    dN1_dy = dNN1_dy;
                  }
                  else
                  {
                    N1 = nd1->SubDivMat*NN1;
                    dN1_dx = nd1->SubDivMat*dNN1_dx;
                    dN1_dy = nd1->SubDivMat*dNN1_dy;
                  }

                  GeomData.computeBasisFunctionsGhostPenalty2D(knotBegin2, knotIncr2, param, NN2, dNN2_dx, dNN2_dy );

                  if(nd2->getParent() == NULL)
                  {
                    N2 = NN2;
                    dN2_dx = dNN2_dx;
                    dN2_dy = dNN2_dy;
                  }
                  else
                  {
                    N2 = nd2->SubDivMat*NN2;
                    dN2_dx = nd2->SubDivMat*dNN2_dx;
                    dN2_dy = nd2->SubDivMat*dNN2_dy;
                  }

                  geom[0] = GeomData.computeCoord(0, param[0]);
                  geom[1] = GeomData.computeCoord(1, param[1]);

                  if(axsy)
                    dvol *= (2.0*PI*geom[0]);

                  for(ii=0; ii<ndof; ii++)
                  {
                    grad1(ii,0) = nd1->computeValueCur(ii, dN1_dx);
                    grad1(ii,1) = nd1->computeValueCur(ii, dN1_dy);

                    grad2(ii,0) = nd2->computeValueCur(ii, dN2_dx);
                    grad2(ii,1) = nd2->computeValueCur(ii, dN2_dy);

                    trac1[ii] = grad1(ii,0)*normal1[0] + grad1(ii,1)*normal1[1];
                    trac2[ii] = grad2(ii,0)*normal2[0] + grad2(ii,1)*normal2[1];

                    res[ii] = trac1[ii] + trac2[ii];
                  }

                  for(ii=0;ii<nlbf1;ii++)
                  {
                      Ta1 = dN1_dx(ii)*normal1[0] + dN1_dy(ii)*normal1[1] ;

                      fact = dvol*Ta1;

                      for(dir=0; dir<ndof; dir++)
                        temp[dir] = gammGP[dir]*fact;

                      TI = ii*ndof ;

                      for(jj=0;jj<nlbf1;jj++)
                      {
                        Tb1 = af*(dN1_dx(jj)*normal1[0] + dN1_dy(jj)*normal1[1]) ;

                        TJ = jj*ndof;
                        for(dir=0; dir<ndof; dir++)
                          K1(TI+dir, TJ+dir)  += temp[dir]*Tb1;
                      }

                      for(jj=0;jj<nlbf2;jj++)
                      {
                        Tb2 = af*(dN2_dx(jj)*normal2[0] + dN2_dy(jj)*normal2[1] );

                        TJ = jj*ndof;
                        for(dir=0; dir<ndof; dir++)
                          Kc(TI+dir, TJ+dir)  += temp[dir]*Tb2;
                      }

                      for(dir=0; dir<ndof; dir++)
                        F1(TI+dir) -= (temp[dir]*res(dir));
                  } // for(ii=0;ii<nlbf1;ii++)

                  for(ii=0;ii<nlbf2;ii++)
                  {
                      Ta2 = dN2_dx(ii)*normal2[0] + dN2_dy(ii)*normal2[1] ;
                      
                      fact = dvol*Ta2;

                      for(dir=0; dir<ndof; dir++)
                        temp[dir] = gammGP[dir]*fact;

                      TI = ii*ndof;

                      for(jj=0;jj<nlbf2;jj++)
                      {
                        Tb2 = af*( dN2_dx(jj)*normal2[0] + dN2_dy(jj)*normal2[1] );

                        TJ = jj*ndof;
                        for(dir=0; dir<ndof; dir++)
                          K2(TI+dir, TJ+dir)  += temp[dir]*Tb2;
                      }

                      for(dir=0; dir<ndof; dir++)
                        F2(TI+dir) -= (temp[dir]*res(dir));
                  } // for(ii=0;ii<nlbf2;ii++)

              }// for(gp1=0...
              }// for(gp2=0...

              //cout << " lllllllllllllll " << endl;

              for(ii=0; ii<nr1; ii++)
              {
                r = grid_to_proc_DOF[nd1->forAssyVec[ii]];

                //solverEigen->rhsVec[r] += F1(ii);
                VecSetValue(solverPetsc->rhsVec, r, F1(ii), ADD_VALUES);

                for(jj=0; jj<nr1; jj++)
                {
                  c = grid_to_proc_DOF[nd1->forAssyVec[jj]];

                  //solverEigen->mtx.coeffRef(r, c) += K1(ii, jj);
                  MatSetValue(solverPetsc->mtx, r, c, K1(ii,jj), ADD_VALUES);
                }

                for(jj=0; jj<nr2; jj++)
                {
                  //cout << ii << '\t' << jj << endl;
                  c = grid_to_proc_DOF[nd2->forAssyVec[jj]];

                  //solverEigen->mtx.coeffRef(r, c) += Kc(ii, jj);
                  //solverEigen->mtx.coeffRef(c, r) += Kc(ii, jj);
                  MatSetValue(solverPetsc->mtx, r, c, Kc(ii,jj), ADD_VALUES);
                  MatSetValue(solverPetsc->mtx, c, r, Kc(ii,jj), ADD_VALUES);
                }
              }
              //cout << " lllllllllllllll " << endl;

              for(ii=0; ii<nr2; ii++)
              {
                r = grid_to_proc_DOF[nd2->forAssyVec[ii]];

                //solverEigen->rhsVec[r] += F2(ii);
                VecSetValue(solverPetsc->rhsVec, r, F2(ii), ADD_VALUES);

                for(jj=0; jj<nr2; jj++)
                {
                  c = grid_to_proc_DOF[nd2->forAssyVec[jj]];

                  //solverEigen->mtx.coeffRef(r, c) += K2(ii, jj);
                  MatSetValue(solverPetsc->mtx, r, c, K2(ii,jj), ADD_VALUES);
                }
              }
              //cout << " lllllllllllllll " << endl;
            }//  if( neighbours[side] != NULL && !neighbours[side]->isGhost() )
          } // for(side=0; side<NUM_NEIGHBOURS
        } //if( nd1->isCutElement() )
      } //if( nd1->getSubdomainId() == this_mpi_proc )
    } //for(ee=0; ee<activeElements.size(); ee++)

    //cout << " HBSplineCutFEM::applyGhostPenalty2D() ... FINISHED  " << endl;

  return;
}




void  HBSplineCutFEM::applyGhostPenalty3D()
{
    //cout << " need to modify HBSplineCutFEM::applyGhostPenalty3D() ... " << endl;

    ////////////////////////////////////////////////////////
    //
    // stiffness and residual due to ghost-penalty terms
    // for 3D problems
    ////////////////////////////////////////////////////////

    int  ii=0, jj=0, ee=0, TI=0, TJ=0, gp=0, side=0, dir=0, NUM_NEIGHBOURS=6, nGauss=0;
    int  nlbf1=0, nlbf2=0, nr1=0, nr2=0, domTemp=0, levTemp=0;
    int  r=0, c=0, start1=0, start2=0, size1=0, size2=0;

    double  Ta1=0.0, Ta2=0.0, Tb1=0.0, Tb2=0.0, JacTemp=0.0, dvol=0.0, fact=0.0;
    double  bb1=0.0, bb2=0.0, PENALTY=0.0, h1=0.0, h2=0.0, h3=0.0;
    vector<double>  gammGP(ndof+1), temp(ndof+1);
    vector<double>  boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2, boundaryGPs3, boundaryGWs3;

    std::vector<int>::iterator itInt;

    int nlf = (degree[0]+1) * (degree[1]+1) * (degree[2]+1);

    VectorXd   NN1(nlf), dNN1_dx(nlf), dNN1_dy(nlf), dNN1_dz(nlf), N1, dN1_dx, dN1_dy, dN1_dz;
    VectorXd   NN2(nlf), dNN2_dx(nlf), dNN2_dy(nlf), dNN2_dz(nlf), N2, dN2_dx, dN2_dy, dN2_dz;
    myPoint   normal1, normal2, hElem;
    myPoint   knotIncr1, knotBegin1, knotEnd1, knotSum1, knotIncr2, knotBegin2, knotEnd2, knotSum2;
    double *gws;
    myPoint *gps;

    node *nd1, *nd2;

    VectorXd   trac1(ndof), trac2(ndof), res(ndof);
    MatrixXd   grad1(ndof, ndof), grad2(ndof, ndof);

    double  af  = SolnData.td(2);
    double  rho = fluidProps[3];
    double  mu  = fluidProps[4];

    // find the maximum h
    // maximum h is the maximum of diameters of all the cut elements
    // where, diameter of a cut element is taken to be the diagonal of that element
    
    ii = elems[cutCellIds[0]]->getLevel();
    
    for(ee=0; ee<cutCellIds.size(); ee++)
    {
      ii = max( ii, elems[cutCellIds[ee]]->getLevel() );
    }

    //cout << " Maximum cut-cell level = " << ii << endl;
    // lengths of element sides at level '0'

    hElem[0] = GeomData.getGridLength(0)/nelem[0];
    hElem[1] = GeomData.getGridLength(1)/nelem[1];
    hElem[2] = GeomData.getGridLength(2)/nelem[2];

    // lengths of element sides at level 'ii'
    fact = pow(2.0, double (ii) );

    hElem[0] /= fact;
    hElem[1] /= fact;
    hElem[2] /= fact;

    h1 = sqrt( hElem[0]*hElem[0] + hElem[1]*hElem[1] + hElem[2]*hElem[2] );

    for(ee=0; ee<cutCellIds.size(); ee++)
    {
      nd1 = elems[cutCellIds[ee]];

      if( nd1->getSubdomainId() == this_mpi_proc )
      {
        if( nd1->isCutElement() )
        {
          knotBegin1 = nd1->getKnotBegin();
          knotEnd1   = nd1->getKnotEnd();
          knotIncr1  = nd1->getKnotIncrement();
          knotSum1   = nd1->getKnotSum();

          //
          //if( nd1->isLeftBoundary() || nd1->isRightBoundary() )
          if( nd1->isBoundary() )
          {
            gammGP[0]  = cutFEMparams[6] * mu*h1;
            gammGP[1]  = gammGP[0];
            gammGP[2]  = gammGP[0];
            gammGP[3]  = cutFEMparams[7] * h1*h1*h1/mu;
          }
          else
          {
            gammGP[0]  = cutFEMparams[6] * mu*h1;
            gammGP[1]  = gammGP[0];
            gammGP[2]  = gammGP[0];
            gammGP[3]  = cutFEMparams[7] * h1*h1*h1/mu;
          }

          bb1 = 1.0;
          for(ii=1; ii<degree[0]; ii++)
          {
            bb1 *= h1*h1;
          }
          gammGP[0] *= bb1;
          gammGP[1] *= bb1;
          gammGP[2] *= bb1;
          gammGP[3] *= bb1;
          //

          for(side=0; side<NUM_NEIGHBOURS; side++)
          {
            //cout << " side = " << side << endl;
            nd2 = nd1->getNeighbour(side);

            //if( nd2 != NULL && !(nd2->isGhost()) )
            if( (nd2 != NULL) && !(nd2->isGhost()) && nd2->isLeaf() && (nd2->isCutElement() || nd2->domNums[0] == 0) )
            {
              nr1 = nd1->forAssyVec.size();
              nr2 = nd2->forAssyVec.size();

              nlbf1 = nr1/ndof;
              nlbf2 = nr2/ndof;

              MatrixXd  K1, K2, Kc;
              VectorXd  F1, F2;

              K1 = MatrixXd::Zero(nr1, nr1);
              K2 = MatrixXd::Zero(nr2, nr2);
              Kc = MatrixXd::Zero(nr1, nr2);
              F1 = VectorXd::Zero(nr1);
              F2 = VectorXd::Zero(nr2);

              knotBegin2 = nd2->getKnotBegin();
              knotEnd2   = nd2->getKnotEnd();
              knotIncr2  = nd2->getKnotIncrement();
              knotSum2   = nd2->getKnotSum();


              normal1 = GeomData.boundaryNormals[side];
              normal2 = -normal1;

              nGauss = GeomData.boundaryQuadrature3D[side].gausspoints.size();

              gps = &(GeomData.boundaryQuadrature3D[side].gausspoints[0]);
              gws = &(GeomData.boundaryQuadrature3D[side].gaussweights[0]);

              JacTemp = GeomData.boundaryJacobians[side][nd2->getLevel()];

              /*
              h1 = abs(hElem.dot(normal1));

              //if( nd1->IsBoundary() )
              //{
                //gammGP[0]  = cutFEMparams[6] * mu*h1;
                //gammGP[1]  = gammGP[0];
                //gammGP[2]  = gammGP[0];
                //gammGP[3]  = cutFEMparams[7] * h1/mu;
              //}
              //else
              //{
                gammGP[0]  = cutFEMparams[6] * mu*h1;
                gammGP[1]  = gammGP[0];
                gammGP[2]  = gammGP[0];
                gammGP[3]  = cutFEMparams[7] * h1*h1*h1/mu;
              //}
              */

              for(gp=0; gp<nGauss; gp++)
              {
                  param[2] = 0.5 * (knotIncr1[2] * gps[gp][2] + knotSum1[2]);
                  param[1] = 0.5 * (knotIncr1[1] * gps[gp][1] + knotSum1[1]);
                  param[0] = 0.5 * (knotIncr1[0] * gps[gp][0] + knotSum1[0]);

                  dvol = JacTemp * gws[gp] ;

                  if( nd2->getDomainNumber() == -1 )
                    dvol *= 0.5;

                  //printf(" %4d \t %4d \t %12.6f \t %12.6f \n", side, dir, vv, uu);

                  GeomData.computeBasisFunctionsGhostPenalty3D(knotBegin1, knotIncr1, param, NN1, dNN1_dx, dNN1_dy, dNN1_dz );

                  if(nd1->getParent() == NULL)
                  {
                    N1 = NN1;
                    dN1_dx = dNN1_dx;
                    dN1_dy = dNN1_dy;
                    dN1_dz = dNN1_dz;
                  }
                  else
                  {
                    N1 = nd1->SubDivMat*NN1;
                    dN1_dx = nd1->SubDivMat*dNN1_dx;
                    dN1_dy = nd1->SubDivMat*dNN1_dy;
                    dN1_dz = nd1->SubDivMat*dNN1_dz;
                  }

                  GeomData.computeBasisFunctionsGhostPenalty3D(knotBegin2, knotIncr2, param, NN2, dNN2_dx, dNN2_dy, dNN2_dz );

                  if(nd2->getParent() == NULL)
                  {
                    N2 = NN2;
                    dN2_dx = dNN2_dx;
                    dN2_dy = dNN2_dy;
                    dN2_dz = dNN2_dz;
                  }
                  else
                  {
                    N2 = nd2->SubDivMat*NN2;
                    dN2_dx = nd2->SubDivMat*dNN2_dx;
                    dN2_dy = nd2->SubDivMat*dNN2_dy;
                    dN2_dz = nd2->SubDivMat*dNN2_dz;
                  }

                  geom[0] = GeomData.computeCoord(0, param[0]);
                  geom[1] = GeomData.computeCoord(1, param[1]);
                  geom[2] = GeomData.computeCoord(2, param[2]);

                  for(ii=0; ii<ndof; ii++)
                  {
                    grad1(ii,0) = nd1->computeValueCur(ii, dN1_dx);
                    grad1(ii,1) = nd1->computeValueCur(ii, dN1_dy);
                    grad1(ii,2) = nd1->computeValueCur(ii, dN1_dz);

                    grad2(ii,0) = nd2->computeValueCur(ii, dN2_dx);
                    grad2(ii,1) = nd2->computeValueCur(ii, dN2_dy);
                    grad2(ii,2) = nd2->computeValueCur(ii, dN2_dz);

                    trac1[ii] = grad1(ii,0)*normal1[0] + grad1(ii,1)*normal1[1] + grad1(ii,2)*normal1[2];
                    trac2[ii] = grad2(ii,0)*normal2[0] + grad2(ii,1)*normal2[1] + grad2(ii,2)*normal2[2];

                    res[ii] = trac1[ii] + trac2[ii];
                  }

                  for(ii=0;ii<nlbf1;ii++)
                  {
                      Ta1 = dN1_dx(ii)*normal1[0] + dN1_dy(ii)*normal1[1] + dN1_dz(ii)*normal1[2] ;

                      fact = dvol*Ta1;

                      for(dir=0; dir<ndof; dir++)
                        temp[dir] = gammGP[dir]*fact;

                      TI = ii*ndof ;

                      for(jj=0;jj<nlbf1;jj++)
                      {
                        Tb1 = af*(dN1_dx(jj)*normal1[0] + dN1_dy(jj)*normal1[1] + dN1_dz(jj)*normal1[2]) ;

                        TJ = jj*ndof;
                        for(dir=0; dir<ndof; dir++)
                          K1(TI+dir, TJ+dir)  += temp[dir]*Tb1;
                      }

                      for(jj=0;jj<nlbf2;jj++)
                      {
                        Tb2 = af*(dN2_dx(jj)*normal2[0] + dN2_dy(jj)*normal2[1] + dN2_dz(jj)*normal2[2] );

                        TJ = jj*ndof;
                        for(dir=0; dir<ndof; dir++)
                          Kc(TI+dir, TJ+dir)  += temp[dir]*Tb2;
                      }

                      for(dir=0; dir<ndof; dir++)
                        F1(TI+dir) -= (temp[dir]*res(dir));
                  } // for(ii=0;ii<nlbf1;ii++)

                  for(ii=0;ii<nlbf2;ii++)
                  {
                      Ta2 = dN2_dx(ii)*normal2[0] + dN2_dy(ii)*normal2[1] + dN2_dz(ii)*normal2[2] ;

                      fact = dvol*Ta2;

                      for(dir=0; dir<ndof; dir++)
                        temp[dir] = gammGP[dir]*fact;

                      TI = ii*ndof;

                      for(jj=0;jj<nlbf2;jj++)
                      {
                        Tb2 = af*( dN2_dx(jj)*normal2[0] + dN2_dy(jj)*normal2[1] + dN2_dz(jj)*normal2[2] );

                        TJ = jj*ndof;
                        for(dir=0; dir<ndof; dir++)
                          K2(TI+dir, TJ+dir)  += temp[dir]*Tb2;
                      }

                      for(dir=0; dir<ndof; dir++)
                        F2(TI+dir) -= (temp[dir]*res(dir));
                  } // for(ii=0;ii<nlbf2;ii++)

              }// for(gp=0...

              //cout << " lllllllllllllll " << endl;

              for(ii=0; ii<nr1; ii++)
              {
                r = grid_to_proc_DOF[nd1->forAssyVec[ii]];

                //solverEigen->rhsVec[r] += F1(ii);
                VecSetValue(solverPetsc->rhsVec, r, F1(ii), ADD_VALUES);

                for(jj=0; jj<nr1; jj++)
                {
                  c = grid_to_proc_DOF[nd1->forAssyVec[jj]];

                  //solverEigen->mtx.coeffRef(r, c) += K1(ii, jj);
                  MatSetValue(solverPetsc->mtx, r, c, K1(ii,jj), ADD_VALUES);
                }

                for(jj=0; jj<nr2; jj++)
                {
                  c = grid_to_proc_DOF[nd2->forAssyVec[jj]];

                  //solverEigen->mtx.coeffRef(r, c) += Kc(ii, jj);
                  //solverEigen->mtx.coeffRef(c, r) += Kc(ii, jj);

                  MatSetValue(solverPetsc->mtx, r, c, Kc(ii,jj), ADD_VALUES);
                  MatSetValue(solverPetsc->mtx, c, r, Kc(ii,jj), ADD_VALUES);
                }
              }
              //cout << " lllllllllllllll " << endl;

              for(ii=0; ii<nr2; ii++)
              {
                r = grid_to_proc_DOF[nd2->forAssyVec[ii]];

                //solverEigen->rhsVec[r] += F2(ii);
                VecSetValue(solverPetsc->rhsVec, r, F2(ii), ADD_VALUES);

                for(jj=0; jj<nr2; jj++)
                {
                  c = grid_to_proc_DOF[nd2->forAssyVec[jj]];

                  //solverEigen->mtx.coeffRef(r, c) += K2(ii, jj);
                  MatSetValue(solverPetsc->mtx, r, c, K2(ii,jj), ADD_VALUES);
                }
              }
              //cout << " lllllllllllllll " << endl;
            }//  if( neighbours[side] != NULL && !neighbours[side]->isGhost() )
          } // for(side=0; side<NUM_NEIGHBOURS
        } //if( nd1->isCutElement() )
      } //if( nd1->getSubdomainId() == this_mpi_proc )
    } //for(ee=0; ee<activeElements.size(); ee++)

  return;
}







