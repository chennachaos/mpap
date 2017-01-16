
#include "HBSplineCutFEM.h"

#include "DataBlockTemplate.h"
#include "SolverEigen.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "NurbsUtilities.h"

#include "BasisFunctionsLagrange.h"
#include "myDataIntegrateCutFEM.h"
#include "ImmersedIntegrationElement.h"
#include "QuadratureUtil.h"
#include "typedefs.h"


extern MpapTime           mpapTime;

using namespace std;





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

      int  aa, bb, ii, jj, nlb, nlf, ind1, ind2, nr, nc, start, r, c, kk, gp, nGauss;
      int TI, TIp1, TIp2, TJ, TJp1, TJp2, totnlbf2, elnum;

      nlf = 1;
      for(ii=0; ii<DIM; ii++)
        nlf *= (degree[ii]+1);

      double  fact, fact1, fact2, dvol, PENALTY, detJ, pres, presPrev, x0, y0, FlowRate=6.545*0.5, x1, y1, rad;
      double  Ta[6], Tb[6], bb1, bb2, NitscheFact, c1, hx, hy;
      bool  isNitsche;
      
      AABB  bbTemp;

      double  af = SolnData.td(2);
      bool axsy = (fluidProps[2] == 1);
      double  rho = fluidProps[3];
      double  mu = fluidProps[4];

      MatrixXd  K1, grad(2,2), stress(2,2);
      VectorXd  F1;
      myPoint   vel, trac, velSpec;
      VectorXd  NN(nlf), dNN_dx(nlf), dNN_dy(nlf), N, dN_dx, dN_dy;

      myPoint  knotBegin, knotIncr;

      VectorXd  Nb, dNb, xx, yy;

      //cout << " axsy " << axsy << '\t' << PI << endl;

      vector<double>  gausspoints, gaussweights;

      normal.setZero();

      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      {
        PENALTY     = ImmersedBodyObjects[bb]->GetPenaltyParameter();
        isNitsche   = ImmersedBodyObjects[bb]->GetNitscheFlag();
        NitscheFact = ImmersedBodyObjects[bb]->GetNitscheFact();

        nlb = ImmersedBodyObjects[bb]->ImmIntgElems[0]->pointNumsGlobal.size();

        nGauss = (int) cutFEMparams[1];

        Nb.resize(nlb);
        dNb.resize(nlb);
        xx.resize(nlb);
        yy.resize(nlb);

        getGaussPoints1D(nGauss, gausspoints, gaussweights);
        //printVector(gausspoints);          printf("\n\n\n\n");
        //printVector(gaussweights);          printf("\n\n\n\n");

        for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
        {
          lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];
          poly = ImmersedBodyObjects[bb]->ImmersedFaces[aa];

          //cout << bb << '\t' << aa << '\t' << lme->IsActive() << endl;

          if( lme->IsActive() )
          {
            for(gp=0;gp<nGauss;gp++)
            {
              //computeLagrangeBFs1D2(nlb-1, gausspoints[gp], &xx(0), &yy(0), &Nb(0), &dNb(0), detJ);

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

              //if(bb==2)
	      //{
                //printf("normal     ... %12.6f \t %12.6f \t %12.6f \n", normal[0], normal[1], normal[2]);
                //printf("velSpec    ... %12.6f \t %12.6f \t %12.6f \n", velSpec[0], velSpec[1], velSpec[2]);
	      //}
              //cout << " dvol = " << dvol << endl;

              //printf("geom = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);

              elnum = findCellNumber(geom);

              nd2 = elems[elnum];

              bbTemp = nd2->GetAABB();
              hx = bbTemp.maxBB[0]-bbTemp.minBB[0];
              hy = bbTemp.maxBB[1]-bbTemp.minBB[1];

              //bbTemp.printSelf();

              //
              if( (nd2->domNums.size() == 1) && (nd2->domNums[0] != 0) )
              {
                //cout << " elnum = " << elnum << endl;
                //printVector(nd2->domNums);
                
                if( abs(geom[1]-bbTemp.minBB[1]) < 1.0e-8 ) // bottom edge
                {
                  nd3 = nd2->GetNeighbour(EAST);
                  
                  if( (nd3->domNums.size() == 1) && (nd3->domNums[0] != 0) )
                  {
                    nd4 = nd2->GetNeighbour(SOUTH);
                    if( (nd4->domNums.size() == 1) && (nd4->domNums[0] != 0) )
                      nd = nd2->GetNeighbour(WEST);
                    else
                      nd = nd4;
                  }
                  else
                    nd = nd2->GetNeighbour(SOUTH);
                }
                else if( abs(geom[1]-bbTemp.maxBB[1]) < 1.0e-8 ) // top edge
                {
                  nd3 = nd2->GetNeighbour(EAST);
                  
                  if( (nd3->domNums.size() == 1) && (nd3->domNums[0] != 0) )
                  {
                    nd4 = nd2->GetNeighbour(NORTH);
                    if( (nd4->domNums.size() == 1) && (nd4->domNums[0] != 0) )
                      nd = nd2->GetNeighbour(WEST);
                    else
                      nd = nd4;
                  }
                  else
                    nd = nd2->GetNeighbour(NORTH);
                }
                else
                {
                  if( abs(geom[0]-bbTemp.minBB[0]) < 1.0e-8 ) // left edge
                    nd = nd2->GetNeighbour(WEST);
                  else if( abs(geom[0]-bbTemp.maxBB[0]) < 1.0e-8 ) // right edge
                    nd = nd2->GetNeighbour(EAST);
                  else
                  {
                    //cerr << " Boundary integration point lies in absurd element \n\n " << endl;
                  }
                  //cout << " llllllllll " << endl;
                  //if( abs(geom[1]-bbTemp.minBB[1]) < 1.0e-10 )
                    //nd = nd2->GetNeighbour(SOUTH);
                  //else if( abs(geom[1]-bbTemp.maxBB[1]) < 1.0e-10 )
                    //nd = nd2->GetNeighbour(NORTH);
                }
              }
              else
                nd = nd2;
              //

              geometryToParametric(geom, param);

              nr = nd->forAssyVec.size();

              K1 = MatrixXd::Zero(nr, nr);
              F1 = VectorXd::Zero(nr);

              //printf("param = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", param[0], param[1], param[2], dvol);
              //printf("vx = %12.6f, vy = %12.6f, vz = %12.6f \n", velSpec[0], velSpec[1], velSpec[2]);

              // stiffness and residual terms due to Nitsche method
              //cout << " stiffness and residual terms due to Nitsche method " << endl;

              knotBegin = nd->GetKnotBegin();
              knotIncr  = nd->GetKnotIncrement();

              GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

              if(nd->GetParent() == NULL)
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

              //cout << " uuuuuuuuuuu " << endl;

              vel[0] = velSpec[0] - nd->computeValueCur(0, N);
              vel[1] = velSpec[1] - nd->computeValueCur(1, N);

              //printf("vel        ... %12.6f \t %12.6f \t %12.6f \n", vel[0], vel[1], vel[2]);

              ///////////////////////////////////
              // compute the stresses and tractions
              ///////////////////////////////////

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

              //cout << " vvvvvvvvvvvvv " << endl;

              totnlbf2 = nr/ndof;

              for(ii=0;ii<totnlbf2;ii++)
              {
                  TI   = ndof*ii;
                  TIp1 = TI+1;
                  TIp2 = TI+2;

                  bb1 = N[ii] * dvol;
                  bb2 = bb1 * PENALTY ;

                  Ta[0] = (dvol*mu)*( normal[0]*dN_dx(ii) + normal[1]*dN_dy(ii) );
                  Ta[1] = 0.0;
                  Ta[2] = -normal[0]*bb1;

                  Ta[3] = 0.0;
                  Ta[4] = (dvol*mu)*( normal[0]*dN_dx(ii) + normal[1]*dN_dy(ii) );
                  Ta[5] = -normal[1]*bb1;

                  for(jj=0; jj<6; jj++)
                    Ta[jj] *= NitscheFact;

                  for(jj=0;jj<totnlbf2;jj++)
                  {
                    TJ   = ndof*jj;
                    TJp1 = TJ+1;
                    TJp2 = TJ+2;

                    fact = bb2 * af *  N[jj];
                    // stabilisation term
                    K1(TI,   TJ)   += fact;
                    K1(TIp1, TJp1) += fact;

                    //K1(TIp2, TJp2) -= (1.0*N[ii]*N[jj]*dvol);

                    // Nitsche terms

                    Tb[0] = (af*mu)*( normal[0]*dN_dx(jj) + normal[1]*dN_dy(jj) );
                    Tb[1] = 0.0;
                    Tb[2] = -normal[0]*N(jj);

                    Tb[3] = 0.0;
                    Tb[4] = (af*mu)*( normal[0]*dN_dx(jj) + normal[1]*dN_dy(jj) );
                    Tb[5] = -normal[1]*N(jj);

                    K1(TI, TJ)      -= (bb1*Tb[0]);
                    K1(TI, TJp1)    -= (bb1*Tb[1]);
                    K1(TI, TJp2)    -= (bb1*Tb[2]);

                    K1(TIp1, TJ)    -= (bb1*Tb[3]);
                    K1(TIp1, TJp1)  -= (bb1*Tb[4]);
                    K1(TIp1, TJp2)  -= (bb1*Tb[5]);

                    c1 = af*N(jj);

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
                  
                  //F1(TIp2) -= (1.0*N[ii]*dvol*(presPrev-pres));

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

              //cout << " uuuuuuuuuuu " << endl;

              solverEigen->AssembleMatrixAndVectorCutFEM(0, 0, nd->forAssyVec, forAssyCutFEM, K1, F1);

              //cout << " fffffffffff " << endl;
            }//for(gp=0...
          } // if( lme->IsActive() )
        }//for(aa=0...
      }//for(bb=0;...
      //

  //cout << " HBSplineCutFEM::applyInterfaceTerms2D() ... FINISHED  " << endl;

  if(cutFEMparams[6] > 1.0e-8)
    applyGhostPenalty2D();

  return;
}
//



/*
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

      int  aa, bb, ii, jj, nlb, nlf, ind1, ind2, nr, nc, start, r, c, kk, gp, nGauss;
      int TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3, totnlbf2;

      nlf = (degree[0]+1)*(degree[1]+1)*(degree[2]+1);

      double  fact, fact1, fact2, fact3, dvol, PENALTY;
      double  detJ, pres, bb1, bb2, NitscheFact, c1;
      double  Ta[4], Tb[4], surfArea=0.0;

      double  af  = SolnData.td(2);
      double  rho = fluidProps[3];
      double  mu  = fluidProps[4];
      bool  isNitsche;

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

      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      {
        //cout << " uuuuuuuuuuu " << endl;

        PENALTY     = ImmersedBodyObjects[bb]->GetPenaltyParameter();
        isNitsche   = ImmersedBodyObjects[bb]->GetNitscheFlag();
        NitscheFact = ImmersedBodyObjects[bb]->GetNitscheFact();

        nlb = ImmersedBodyObjects[bb]->ImmIntgElems[0]->pointNumsGlobal.size();

        nGauss = (int) cutFEMparams[1];

        Nb.resize(nlb);
        //dNb.resize(nlb);          dNb_dx.resize(nlb);
        //dNb_dy.resize(nlb);          dNb_dz.resize(nlb);

        //xx.resize(nlb);          yy.resize(nlb);          zz.resize(nlb);

        getGaussPointsTriangle(nGauss, gausspoints1, gausspoints2, gaussweights);
        //printVector(gausspoints);          printf("\n\n\n\n");
        //printVector(gaussweights);          printf("\n\n\n\n");

        surfArea = 0.0;

        for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
        {
          lme  = ImmersedBodyObjects[bb]->ImmIntgElems[aa];
          poly = ImmersedBodyObjects[bb]->ImmersedFaces[aa];
          
          if( lme->IsActive() )
          {
            for(gp=0;gp<nGauss;gp++)
            {
              param[0] = gausspoints1[gp];
              param[1] = gausspoints2[gp];
              
              //cout << " Need to implement computing basis functions for boundary triangles in 3D " << endl;
              //computeLagrangeBFs1D2(nlb-1, gausspoints[gp], &xx(0), &yy(0), &Nb(0), &dNb(0), detJ);

              poly->computeBasisFunctions(param, geom, Nb, detJ);

              dvol  = gaussweights[gp] * detJ;

              //normal = poly->GetNormal();
              poly->computeNormal();
              poly->computeNormal(param, normal);

              normal *= -1.0;

              lme->computeVelocityCur(Nb, velSpec );

              ndTemp = elems[findCellNumber(geom)];

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

              //cout << " stiffness and residual terms due to Nitsche method " << endl;

              knotBegin = ndTemp->GetKnotBegin();
              knotIncr  = ndTemp->GetKnotIncrement();

              GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

              if(ndTemp->GetParent() == NULL)
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

              stress(0,0) = mu*ndTemp->computeValueCur(0, dN_dx);
              stress(0,1) = mu*ndTemp->computeValueCur(0, dN_dy);
              stress(0,2) = mu*ndTemp->computeValueCur(0, dN_dz);

              stress(1,0) = mu*ndTemp->computeValueCur(1, dN_dx);
              stress(1,1) = mu*ndTemp->computeValueCur(1, dN_dy);
              stress(1,2) = mu*ndTemp->computeValueCur(1, dN_dz);

              stress(2,0) = mu*ndTemp->computeValueCur(2, dN_dx);
              stress(2,1) = mu*ndTemp->computeValueCur(2, dN_dy);
              stress(2,2) = mu*ndTemp->computeValueCur(2, dN_dz);

              pres        =    ndTemp->computeValue(3, N);

              stress(0,0) -= pres;
              stress(1,1) -= pres;
              stress(2,2) -= pres;

              trac[0] = stress(0,0)*normal[0] + stress(0,1)*normal[1] + stress(0,2)*normal[2] ;
              trac[1] = stress(1,0)*normal[0] + stress(1,1)*normal[1] + stress(1,2)*normal[2] ;
              trac[2] = stress(2,0)*normal[0] + stress(2,1)*normal[1] + stress(2,2)*normal[2] ;

              //cout << " vvvvvvvvvvvvv " << endl;

              totnlbf2 = nr/ndof;

              for(ii=0;ii<totnlbf2;ii++)
              {
                  TI   = 4*ii;
                  TIp1 = TI+1;
                  TIp2 = TI+2;
                  TIp3 = TI+3;

                  bb1 = N[ii] * dvol;
                  bb2 = bb1 * PENALTY ;

                  Ta[0]  = (dvol*mu)*( dN_dx(ii)*normal[0] + dN_dy(ii)*normal[1] + dN_dz(ii)*normal[2] );

                  Ta[1]  = -normal[0]*bb1;
                  Ta[2]  = -normal[1]*bb1;
                  Ta[3]  = -normal[2]*bb1;
                  
                  for(jj=0; jj<4; jj++)
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
                    Tb[0]  = af*mu*( dN_dx(jj)*normal[0] + dN_dy(jj)*normal[1] + dN_dz(jj)*normal[2] );

                    Tb[1]  = -normal[0]*N[jj];
                    Tb[2]  = -normal[1]*N[jj];
                    Tb[3]  = -normal[2]*N[jj];

                    c1 = af*N[jj];

                    K1(TI, TJ)      -= (bb1*Tb[0]);
                    K1(TI, TJp3)    -= (bb1*Tb[1]);

                    K1(TI, TJ)      -= (Ta[0]*c1);
                    K1(TIp3, TJ)    -= (Ta[1]*c1);

                    K1(TIp1, TJp1)  -= (bb1*Tb[0]);
                    K1(TIp1, TJp3)  -= (bb1*Tb[2]);

                    K1(TIp1, TJp1)  -= (Ta[0]*c1);
                    K1(TIp3, TJp1)  -= (Ta[2]*c1);

                    K1(TIp2, TJp2)  -= (bb1*Tb[0]);
                    K1(TIp2, TJp3)  -= (bb1*Tb[3]);

                    K1(TIp2, TJp2)  -= (Ta[0]*c1);
                    K1(TIp3, TJp2)  -= (Ta[3]*c1);
                  }

                  // stabilisation terms

                  F1(TI)   += (bb2*vel[0]);
                  F1(TIp1) += (bb2*vel[1]);
                  F1(TIp2) += (bb2*vel[2]);

                  // Nitsche terms

                  F1(TI)   += (bb1*trac[0]);
                  F1(TI)   -= (Ta[0]*vel[0]);
                  F1(TIp3) -= (Ta[1]*vel[0]);

                  F1(TIp1) += (bb1*trac[1]);
                  F1(TIp1) -= (Ta[0]*vel[1]);
                  F1(TIp3) -= (Ta[2]*vel[1]);

                  F1(TIp2) += (bb1*trac[2]);
                  F1(TIp2) -= (Ta[0]*vel[2]);
                  F1(TIp3) -= (Ta[3]*vel[2]);

              } // for(ii=0;ii<totnlbf2;ii++)

              //cout << " uuuuuuuuuuu \n\n " << endl;

              solverEigen->AssembleMatrixAndVectorCutFEM(0, 0, ndTemp->forAssyVec, forAssyCutFEM, K1, F1);

            }//for(gp=0...
          } // if( lme->IsActive() )
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

    int  ii, jj, ee, TI, TJ, gp1, gp2, side, dir, NUM_NEIGHBOURS=4;
    int  nlbf1, nlbf2, nr1, nr2, nlf, domTemp;
    int  r, c, start1, start2, size1, size2;

    double  Ta1, Ta2, Tb1, Tb2, JacMultLoc, rad, dvol, fact;
    double  bb1, bb2, PENALTY, h1, h2;
    double  *knots1[2], *knots2[2];
    vector<double>  boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2;
    vector<double>  gammGP(ndof+1), temp(ndof+1);

    std::vector<int>::iterator itInt;

    nlf = 1;
    for(ii=0; ii<DIM; ii++)
      nlf *= (degree[ii]+1);

    VectorXd   NN1(nlf), dNN1_dx(nlf), dNN1_dy(nlf), N1, dN1_dx, dN1_dy;
    VectorXd   NN2(nlf), dNN2_dx(nlf), dNN2_dy(nlf), N2, dN2_dx, dN2_dy;
    myPoint   normal1, normal2, hElem;
    myPoint   knotIncr1, knotBegin1, knotIncr2, knotBegin2;

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

    
    ii = elems[cutCellIds[0]]->GetLevel();
    
    for(ee=0; ee<cutCellIds.size(); ee++)
    {
      ii = max( ii, elems[cutCellIds[ee]]->GetLevel() );
    }

    // cout << " Maximum cut-cell level = " << ii << endl;
    // lengths of element sides at level '0'

    hElem[0] = GeomData.GetGridLength(0)/nelem[0];
    hElem[1] = GeomData.GetGridLength(1)/nelem[1];

    // lengths of element sides at level 'ii'
    fact = pow(2.0, double (ii) );

    hElem[0] /= fact;
    hElem[1] /= fact;

    h1 = sqrt( hElem[0]*hElem[0] + hElem[1]*hElem[1] );


    for(ee=0; ee<cutCellIds.size(); ee++)
    {
      nd1 = elems[cutCellIds[ee]];

      //if( nd1->IsCutElement() && !( nd1->IsLeftBoundary() || nd1->IsRightBoundary() ) )
      if( nd1->IsCutElement() )
      {
        //cout << " nd1->GetID() " <<  nd1->GetID() << '\t' <<  nd1->GetLevel() << endl;

        knots1[0] = nd1->GetKnots(Dir1);
        knots1[1] = nd1->GetKnots(Dir2);

        knotBegin1 = nd1->GetKnotBegin();
        knotIncr1  = nd1->GetKnotIncrement();

        //cout << " nd1->IsRightBoundary() = " << nd1->IsRightBoundary() << endl;

        //
        //if( nd1->IsLeftBoundary() || nd1->IsRightBoundary() )
        if( nd1->IsBoundary() )
        {
          gammGP[0]  = 0.0;
          //gammGP[0]  = cutFEMparams[6] *mu* h1*h1*h1;
          //gammGP[0]  = cutFEMparams[6] *mu*h1/rho;
          gammGP[0]  = cutFEMparams[6] *mu*h1;
          gammGP[2]  = cutFEMparams[7] * h1*h1*h1/mu;
          gammGP[1]  = gammGP[0];
          //gammGP[2]  = gammGP[0];
        }
        else
        {
          //gammGP[0]  = cutFEMparams[6] * mu*h1/rho;
          gammGP[0]  = cutFEMparams[6] * mu*h1;
          gammGP[1]  = gammGP[0];
          gammGP[2]  = cutFEMparams[7] * h1*h1*h1/mu;
          //gammGP[2]  = gammGP[0];
        }

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
          nd2 = nd1->GetNeighbour(side);

          //cout << " side = " << side << '\t' << nd2->GetDomainNumber() << endl;
          //if( !(nd2 == NULL) && !(nd2->IsGhost()) && nd2->IsLeaf() && (nd2->GetDomainNumber() == -1) && (nd2->GetDomainNumber() >= 1) )
          //if( !(nd2 == NULL) && !(nd2->IsGhost()) && nd2->IsLeaf() )
          if( (nd2 != NULL) && !(nd2->IsGhost()) && nd2->IsLeaf() && (nd2->IsCutElement() || nd2->domNums[0] == 0) )
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

              knots2[0] = nd2->GetKnots(Dir1);
              knots2[1] = nd2->GetKnots(Dir2);

              knotBegin2 = nd2->GetKnotBegin();
              knotIncr2  = nd2->GetKnotIncrement();


              GeomData.getBoundaryNormal2D(side, normal1);
              GeomData.setBoundaryGPs2D(side, boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2);

              normal2 = -normal1;

              JacMultLoc = nd1->getJacBoundary(side);
        
              for(gp2=0;gp2<boundaryGPs2.size();gp2++)
              {
                  param[1] = 0.5 * (knots1[1][2] * boundaryGPs2[gp2] + knots1[1][3]);
              for(gp1=0;gp1<boundaryGPs1.size();gp1++)
              {
                  param[0] = 0.5 * (knots1[0][2] * boundaryGPs1[gp1] + knots1[0][3]);
 
                  dvol = JacMultLoc * boundaryGWs2[gp2] * boundaryGWs1[gp1] ;

                  // multiply dvol with 0.5 as the ghost-penalty operation is applied
                  // twice on the face shared by the two elements both of which are cut elements

                  //if( nd2->GetDomainNumber() == -1 )
                    //dvol *= 0.5;

                  //printf(" %4d \t %4d \t %12.6f \t %12.6f \n", side, dir, vv, uu);

                  GeomData.computeBasisFunctionsGhostPenalty2D(knotBegin1, knotIncr1, param, NN1, dNN1_dx, dNN1_dy );

                  if(nd1->GetParent() == NULL)
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

                  if(nd2->GetParent() == NULL)
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

                  geom[0] = GeomData.ComputeCoord(0, param[0]);
                  geom[1] = GeomData.ComputeCoord(1, param[1]);
                  //rad = sqrt(xx*xx+yy*yy);
                  //specVal = 1.0+log(2.0*rad);
                  //cout << xx << '\t' << yy << endl;
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
                r = forAssyCutFEM[nd1->forAssyVec[ii]];

                solverEigen->rhsVec[r] += F1(ii);

                for(jj=0; jj<nr1; jj++)
                {
                  c = forAssyCutFEM[nd1->forAssyVec[jj]];

                  solverEigen->mtx.coeffRef(r, c) += K1(ii, jj);
                }

                for(jj=0; jj<nr2; jj++)
                {
                  //cout << ii << '\t' << jj << endl;
                  c = forAssyCutFEM[nd2->forAssyVec[jj]];

                  solverEigen->mtx.coeffRef(r, c) += Kc(ii, jj);
                  solverEigen->mtx.coeffRef(c, r) += Kc(ii, jj);
                }
              }
              //cout << " lllllllllllllll " << endl;

              for(ii=0; ii<nr2; ii++)
              {
                r = forAssyCutFEM[nd2->forAssyVec[ii]];

                solverEigen->rhsVec[r] += F2(ii);

                for(jj=0; jj<nr2; jj++)
                {
                  c = forAssyCutFEM[nd2->forAssyVec[jj]];

                  solverEigen->mtx.coeffRef(r, c) += K2(ii, jj);
                }
              }
              //cout << " lllllllllllllll " << endl;
          }//  if( neighbours[side] != NULL && !neighbours[side]->IsGhost() )
        } // for(side=0; side<NUM_NEIGHBOURS

      } //if( nd1->IsCutElement() )
    } //for(ee=0; ee<activeElements.size(); ee++)

    //cout << " HBSplineCutFEM::applyGhostPenalty2D() ... FINISHED  " << endl;

  return;
}




void  HBSplineCutFEM::applyGhostPenalty3D()
{
    ////////////////////////////////////////////////////////
    //
    // stiffness and residual due to ghost-penalty terms
    // for 3D problems
    ////////////////////////////////////////////////////////

    int  ii, jj, ee, TI, TJ, gp, side, dir, NUM_NEIGHBOURS=6, nGauss, levTemp;
    int  nlbf1, nlbf2, nr1, nr2, nlf, domTemp;
    int  r, c, start1, start2, size1, size2;

    double  Ta1, Ta2, Tb1, Tb2, JacTemp, dvol, fact;
    double  bb1, bb2, PENALTY, h1, h2, h3;
    double  *knots1[3], *knots2[3];
    vector<double>  gammGP(ndof+1), temp(ndof+1);
    vector<double>  boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2, boundaryGPs3, boundaryGWs3;

    std::vector<int>::iterator itInt;

    nlf = 1;
    for(ii=0; ii<DIM; ii++)
      nlf *= (degree[ii]+1);

    VectorXd   NN1(nlf), dNN1_dx(nlf), dNN1_dy(nlf), dNN1_dz(nlf), N1, dN1_dx, dN1_dy, dN1_dz;
    VectorXd   NN2(nlf), dNN2_dx(nlf), dNN2_dy(nlf), dNN2_dz(nlf), N2, dN2_dx, dN2_dy, dN2_dz;
    myPoint   normal1, normal2, hElem;
    myPoint  knotIncr1, knotBegin1, knotIncr2, knotBegin2;
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

    
    ii = elems[cutCellIds[0]]->GetLevel();
    
    for(ee=0; ee<cutCellIds.size(); ee++)
    {
      ii = max( ii, elems[cutCellIds[ee]]->GetLevel() );
    }

    //cout << " Maximum cut-cell level = " << ii << endl;
    // lengths of element sides at level '0'

    hElem[0] = GeomData.GetGridLength(0)/nelem[0];
    hElem[1] = GeomData.GetGridLength(1)/nelem[1];
    hElem[2] = GeomData.GetGridLength(2)/nelem[2];

    // lengths of element sides at level 'ii'
    fact = pow(2.0, double (ii) );

    hElem[0] /= fact;
    hElem[1] /= fact;
    hElem[2] /= fact;

    h1 = sqrt( hElem[0]*hElem[0] + hElem[1]*hElem[1] + hElem[2]*hElem[2] );

    for(ee=0; ee<cutCellIds.size(); ee++)
    {
      nd1 = elems[cutCellIds[ee]];

      if( nd1->IsCutElement() )
      {
        //cout << " nd1->GetID() " <<  nd1->GetID() << '\t' <<  nd1->GetLevel() << endl;

        knots1[0] = nd1->GetKnots(Dir1);
        knots1[1] = nd1->GetKnots(Dir2);
        knots1[2] = nd1->GetKnots(Dir3);

        knotBegin1 = nd1->GetKnotBegin();
        knotIncr1  = nd1->GetKnotIncrement();

        /*
        //if( nd1->IsLeftBoundary() || nd1->IsRightBoundary() )
        if( nd1->IsBoundary() )
        {
          gammGP[0]  = cutFEMparams[6] * mu*h1;
          gammGP[1]  = gammGP[0];
          gammGP[2]  = gammGP[0];
          gammGP[3]  = cutFEMparams[7] * h1*h1*h1/mu;
          //gammGP[3]  = gammGP[0];
        }
        else
        {
          gammGP[0]  = cutFEMparams[6] * mu*h1;
          gammGP[1]  = gammGP[0];
          gammGP[2]  = gammGP[0];
          gammGP[3]  = cutFEMparams[7] * h1*h1*h1/mu;
          //gammGP[3]  = gammGP[0];
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
        */

        for(side=0; side<NUM_NEIGHBOURS; side++)
        {
          //cout << " side = " << side << endl;
          nd2 = nd1->GetNeighbour(side);

          //if( nd2 != NULL && !(nd2->IsGhost()) )
          if( (nd2 != NULL) && !(nd2->IsGhost()) && nd2->IsLeaf() && (nd2->IsCutElement() || nd2->domNums[0] == 0) )
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

              knots2[0] = nd2->GetKnots(Dir1);
              knots2[1] = nd2->GetKnots(Dir2);
              knots2[2] = nd2->GetKnots(Dir3);

              knotBegin2 = nd2->GetKnotBegin();
              knotIncr2  = nd2->GetKnotIncrement();

              normal1 = GeomData.boundaryNormals[side];
              normal2 = -normal1;

              nGauss = GeomData.boundaryQuadrature3D[side].gausspoints.size();

              gps = &(GeomData.boundaryQuadrature3D[side].gausspoints[0]);
              gws = &(GeomData.boundaryQuadrature3D[side].gaussweights[0]);

              JacTemp = GeomData.boundaryJacobians[side][nd2->GetLevel()];

              //
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

              for(gp=0; gp<nGauss; gp++)
              {
                  param[2] = 0.5 * (knots1[2][2] * gps[gp][2] + knots1[2][3]);
                  param[1] = 0.5 * (knots1[1][2] * gps[gp][1] + knots1[1][3]);
                  param[0] = 0.5 * (knots1[0][2] * gps[gp][0] + knots1[0][3]);

                  dvol = JacTemp * gws[gp] ;

                  if( nd2->GetDomainNumber() == -1 )
                    dvol *= 0.5;

                  //printf(" %4d \t %4d \t %12.6f \t %12.6f \n", side, dir, vv, uu);

                  GeomData.computeBasisFunctionsGhostPenalty3D(knotBegin1, knotIncr1, param, NN1, dNN1_dx, dNN1_dy, dNN1_dz );

                  if(nd1->GetParent() == NULL)
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

                  if(nd2->GetParent() == NULL)
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

                  geom[0] = GeomData.ComputeCoord(0, param[0]);
                  geom[1] = GeomData.ComputeCoord(1, param[1]);
                  geom[2] = GeomData.ComputeCoord(2, param[2]);
                  //cout << xx << '\t' << yy << endl;

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

                  //fact  = dvol * PENALTY;

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
                r = forAssyCutFEM[nd1->forAssyVec[ii]];

                solverEigen->rhsVec[r] += F1(ii);

                for(jj=0; jj<nr1; jj++)
                {
                  c = forAssyCutFEM[nd1->forAssyVec[jj]];

                  solverEigen->mtx.coeffRef(r, c) += K1(ii, jj);
                }

                for(jj=0; jj<nr2; jj++)
                {
                  //cout << ii << '\t' << jj << endl;
                  c = forAssyCutFEM[nd2->forAssyVec[jj]];

                  solverEigen->mtx.coeffRef(r, c) += Kc(ii, jj);
                  solverEigen->mtx.coeffRef(c, r) += Kc(ii, jj);
                }
              }
              //cout << " lllllllllllllll " << endl;

              for(ii=0; ii<nr2; ii++)
              {
                r = forAssyCutFEM[nd2->forAssyVec[ii]];

                solverEigen->rhsVec[r] += F2(ii);

                for(jj=0; jj<nr2; jj++)
                {
                  c = forAssyCutFEM[nd2->forAssyVec[jj]];

                  solverEigen->mtx.coeffRef(r, c) += K2(ii, jj);
                }
              }
              //cout << " lllllllllllllll " << endl;
          }//  if( neighbours[side] != NULL && !neighbours[side]->IsGhost() )
        } // for(side=0; side<NUM_NEIGHBOURS

      } //if( nd1->IsCutElement() )
    } //for(ee=0; ee<activeElements.size(); ee++)


  return;
}







