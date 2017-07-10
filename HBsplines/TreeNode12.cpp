
#include "TreeNode.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "Functions.h"
#include "SolutionData.h"
#include "BasisFunctionsBSpline.h"
#include "myDataIntegrateCutFEM.h"
#include "myPoly.h"

#include "ExactSolutionsElasticity.h"


extern  MpapTime  mpapTime;
extern List<TimeFunction> timeFunction;



template<>
void TreeNode<2>::calcStiffnessAndResidualCutFEMPoisson(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Poisson equation
    //
    ///////////////////////////////////////

    int ii, jj, kk, gp, nGauss, tempId, gpDomainId;

    double  JacTemp, Jac, dvol, force;
    double  fact, fact2, bb1, bb2, bb3, mu;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  N, dN_dx, dN_dy, grad(2);
    myPoint  param, geom;

    //ElasticityEx3  analy(BULK, mu);
    PoissonInterfaceEx4  analy(0.5);

    double *gws;
    myPoint *gps;
    
    if(domNums.size() > 1)
    {
      nGauss = Quadrature.gausspoints.size();
      
      gps = &(Quadrature.gausspoints[0]);
      gws = &(Quadrature.gaussweights[0]);

      tempId = domainCur;

      JacTemp = 1.0;
    }
    else
    {
      nGauss = GeomData->gausspoints.size();

      gps = &(GeomData->gausspoints[0]);
      gws = &(GeomData->gaussweights[0]);

      tempId = domNums[0];

      JacTemp = JacMultElem;
    }

    //cout << " nGauss " << nGauss << '\t' << tempId << endl;

    for(gp=0; gp<nGauss; gp++)
    {
        param[0]  = 0.5*(knots[0][2] * gps[gp][0] + knots[0][3]);
        param[1]  = 0.5*(knots[1][2] * gps[gp][1] + knots[1][3]);

        dvol = gws[gp] * JacTemp;

        //cout << gp << '\t' << gps[gp][0] << '\t' << gps[gp][1] << '\t' << gws[gp] << '\t' << dvol << endl;

        GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

        if(parent == NULL)
        {
          N = NN;
          dN_dx = dNN_dx;
          dN_dy = dNN_dy;
        }
        else
        {
          N = SubDivMat*NN;
          dN_dx = SubDivMat*dNN_dx;
          dN_dy = SubDivMat*dNN_dy;
        }

        geom[0] = GeomData->computeCoord(0, param[0]);
        geom[1] = GeomData->computeCoord(1, param[1]);

        //cout << uu << '\t' << vv << endl;
        //cout << geom[0] << '\t' << geom[1] << endl;

        //cout <<  " gpDomainId = " << gpDomainId << '\t' << totnlbf2 << endl;
        //printVector(GlobalBasisFuncs);

        //if( GeomData->polyImm.within(xx, yy) )
        if(domNums.size() > 1)
        {
          gpDomainId = QuadratureDomNums[gp];
          
        }
        else
        {
          gpDomainId = domNums[0];
        }

      if(gpDomainId == tempId)
      {
          mu  = elmDat[4+gpDomainId];

          //cout << E << '\t' << nu << endl;

          grad.setZero();
          for(ii=0;ii<totnlbf2;ii++)
          {
            jj  =  GlobalBasisFuncs[ii];

            bb1 =  SolnData->solnCFM(gpDomainId, jj);

            grad(0) += ( bb1 * dN_dx(ii) );
            grad(1) += ( bb1 * dN_dy(ii) );
          }

          force = 0.0;
          //force = analy.computeForce(0, geom[0], geom[1]);

        for(ii=0;ii<totnlbf2;ii++)
        {
          bb1 = dN_dx[ii]*dvol;
          bb2 = dN_dy[ii]*dvol;
          bb3 = N[ii]*dvol;

          Flocal[ii]   += ( bb3*force - bb1*grad[0] - bb2*grad[1]) ;

          for(jj=0;jj<totnlbf2;jj++)
          {
              Klocal(ii, jj)  +=  ( bb1*dN_dx[jj] + bb2*dN_dy[jj] ) ;
          }
        }
      }
    }//gp

    return;
}



template<>
void TreeNode<2>::applyBoundaryConditionsAtApointCutFEMPoisson(myDataIntegrateCutFEM& myData)
{
  // at this moment no jumps in velocity or tractions is considered

  // compute stiffness and force vectors corresponding 
  // to Nitsche method of applying interface conditions
  // 
  // diagonal terms
  
    PoissonEx3  analy;
    //PoissonInterfaceEx4  analy(0.5);

    int ii, jj, TI, TIp1, TIp2, TJ, TJp1, TJp2, k1, k2;

    double  fact, fact1, fact2, res, trac;
    double  Ta1, Ta2, Tb1, Tb2, mu1, mu2, specVal;
    double  bb1, bb2, u1, u2, t1, t2, pres;
    double  val1, val2, trac1, trac2;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), N, dN_dx, dN_dy, grad1(2), grad2(2);
    myPoint   normal1, normal2;

    double  af = SolnData->td(2);

    double  g1, g2, jumpi, jumpj;

    GeomData->computeBasisFunctions2D(knotBegin, knotIncr, myData.param, NN, dNN_dx, dNN_dy);

    if(parent == NULL)
    {
      N = NN;
      dN_dx = dNN_dx;
      dN_dy = dNN_dy;
    }
    else
    {
      N = SubDivMat*NN;
      dN_dx = SubDivMat*dNN_dx;
      dN_dy = SubDivMat*dNN_dy;
    }
    
    // at this moment no jumps in velocity or tractions is considered

    val1 = myData.specVal[0];
    val1 = analy.computeValue(0, myData.geom[0], myData.geom[1]);
    val2 = val1;
    
    //printf("\t %12.6f \t %12.6f \n", myData.geom[0], myData.geom[1] );
    //printf("\t %12.6f \t %12.6f \n", val1, val2 );

    //val1 -= computeValueCur(0, N);
    //val2 -= computeValue2Cur(0, N);

    grad1.setZero();
    grad2.setZero();
    
    k1 = domNums[0];
    k2 = domNums[1];
    
    for(ii=0;ii<totnlbf2;ii++)
    {
      jj   =  GlobalBasisFuncs[ii];

      bb1  =  SolnData->solnCFM(k1, jj);
      bb2  =  SolnData->solnCFM(k2, jj);
      
      val1 -= bb1*N(ii);

      grad1(0) += ( bb1 * dN_dx(ii) );
      grad1(1) += ( bb1 * dN_dy(ii) );

      val2 -= bb2*N(ii);

      grad2(0) += ( bb2 * dN_dx(ii) );
      grad2(1) += ( bb2 * dN_dy(ii) );
    }

    ///////////////////////////////////
    // compute the stresses and tractions
    ///////////////////////////////////

    // domain 1

    mu1 = elmDat[4];
    normal1 = -myData.normal;

    //grad1(0) = mu1*computeValueCur(0, dN_dx);
    //grad1(1) = mu1*computeValueCur(0, dN_dy);
    grad1 = mu1*grad1;

    trac1 = 0.0 + ( grad1(0)*normal1[0] + grad1(1)*normal1[1] );

    // domain 2

    mu2 = elmDat[5];
    normal2 = -normal1;

    //grad2(0) = mu2*computeValue2Cur(0, dN_dx);
    //grad2(1) = mu2*computeValue2Cur(0, dN_dy);
    grad2 = mu2*grad2;

    trac2 = 0.0 + ( grad2(0)*normal2[0] + grad2(1)*normal2[1] );

    //for(ii=0;ii<totnlbf;ii++)
      //printf(" \t %14.8f \n", N[ii]);

    //specVal = GeomData->analyDBC->computeValue(0, myData.geom[0], myData.geom[1]);

    g1 = g2 = 1.0;
    //g1 = g2 = 0.5;

    //cout << " myData.PENALTY = " << myData.PENALTY << endl;

    for(ii=0;ii<totnlbf2;ii++)
    {
        bb1 = N[ii] * myData.dvol;
        bb2 = bb1 * myData.PENALTY ;

        Ta1 = (myData.dvol*mu1)*( dN_dx(ii)*normal1[0] + dN_dy(ii)*normal1[1] );
        Ta2 = (myData.dvol*mu2)*( dN_dx(ii)*normal2[0] + dN_dy(ii)*normal2[1] );

        for(jj=0;jj<totnlbf2;jj++)
        {
          fact = af* bb2 * N[jj];
          // stabilisation term
          myData.K1(ii, jj)   += fact;

          myData.K2(ii, jj)   += fact;

          // Nitsche terms
          Tb1 = af*mu1*( dN_dx(jj)*normal1[0] + dN_dy(jj)*normal1[1] );

          myData.K1(ii, jj)      += (g1*bb1*Tb1);

          myData.K1(ii, jj)      += (g1*Ta1*af*N(jj))*myData.NitscheFact;

          Tb2 = af*mu2*( dN_dx(jj)*normal2[0] + dN_dy(jj)*normal2[1] );

          myData.K2(ii, jj)      += (g2*bb1*Tb2);

          myData.K2(ii, jj)      += (g2*Ta2*af*N(jj))*myData.NitscheFact;
        }

        // stabilisation terms due to jumpi
        myData.F1(ii)   += (bb2*val1);
        myData.F2(ii)   += (bb2*val2);

        // Nitsche terms

        myData.F1(ii)   += (g1*-bb1*trac1);
        myData.F1(ii)   += (g1*Ta1*val1)*myData.NitscheFact;

        myData.F2(ii)   += (g2*-bb1*trac2);
        myData.F2(ii)   += (g2*Ta2*val2)*myData.NitscheFact;
    } // for(ii=0;ii<totnlbf2;ii++)
    //

   return;
}



template<>
void TreeNode<2>::applyBoundaryConditionsAtApointCutFEMPoisson2(myDataIntegrateCutFEM& myData)
{
  // compute stiffness and force vectors corresponding 
  // to Nitsche method of applying interface conditions
  // 
  // coupling terms


   return;
}





template<>
void TreeNode<2>::applyDirichletBCsCutFEMPoisson(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  if(DirichletData.size() > 0)
  {
    if( isCutElement() )
    {
      vector<int>  domTemp;
      std::vector<int>::iterator  itTemp;

      for(vector<myPoly*>::iterator poly = subTrias.begin() ; poly != subTrias.end(); ++poly)
      {
        domTemp.push_back( (*poly)->getDomainNumber() );
      }
      //printVector(domTemp);

      //if( ! any_of( domTemp.begin(), domTemp.end(), [](int aa){return aa==0;} ) )
      if( my_any_of(domTemp, 0) )
      {
        //cout << " true " << endl;
        return;
      }
    }


    int ii, jj, aa, gp, nGauss, index, dir, side;
    double  theta, y0, y1, Ta, Tb, res, JacTemp, NitscheFact;
    double  dvol, specVal, PENALTY, Jac, fact, rad, R, bb1, bb2;
    bool  isNitsche;

    double  mu, trac;

  //ThickCylinder  analy(sss, matDat[1], matDat[2]);
  //PlateWithHole  analy(sss, matDat[1], matDat[2]);
  //ElasticityEx2  analy(matDat[0], matDat[1]);
    //ElasticityEx3  analy(BULK, mu);

    PoissonEx3  analy;

    double  af = SolnData->td(2);
    af = 1.0;

    VectorXd  N, dN_dx, dN_dy;
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), grad(2);
    myPoint  param, normal, geom;
    double *gws;
    myPoint *gps;
    
    y0 = 0.0;
    y1 = 0.41;

    for(aa=0;aa<DirichletData.size();aa++)
    {
        // printVector(DirichletData[aa]);

        isNitsche   = false;
        side        = (int) (DirichletData[aa][0] - 1);
        dir         = (int) (DirichletData[aa][1] - 1);
        specVal     = DirichletData[aa][2];
        PENALTY     = DirichletData[aa][3];
        isNitsche   = ( (int) DirichletData[aa][4] == 1 );
        NitscheFact = DirichletData[aa][5];

        //for symmetric Nitsche method -> NitscheFact = 1.0
        //for unsymmetric Nitsche method -> NitscheFact = -1.0

        normal = GeomData->boundaryNormals[side];

        if(domNums.size() > 1)
        {
          nGauss = BoundaryQuadrature[side].gausspoints.size() ;

          gps = &(BoundaryQuadrature[side].gausspoints[0]);
          gws = &(BoundaryQuadrature[side].gaussweights[0]);

          JacTemp = 1.0;
        }
        else
        {
          nGauss = GeomData->boundaryQuadrature2D[side].gausspoints.size();

          gps = &(GeomData->boundaryQuadrature2D[side].gausspoints[0]);
          gws = &(GeomData->boundaryQuadrature2D[side].gaussweights[0]);

          JacTemp = GeomData->boundaryJacobians[side][level];
        }

        for(gp=0; gp<nGauss; gp++)
        {
            param[1] = 0.5 * (knots[1][2] * gps[gp][1] + knots[1][3]);
            param[0] = 0.5 * (knots[0][2] * gps[gp][0] + knots[0][3]);

            dvol = JacTemp * gws[gp] ;

            //printf(" %4d \t %4d \t %12.6f \t %12.6f \n", side, dir, vv, uu);

            geom[0] = GeomData->computeCoord(0, param[0]);
            geom[1] = GeomData->computeCoord(1, param[1]);
            //yy = GeomData->computeCoord(1, param[1]);
            //rad = sqrt(xx*xx+yy*yy);
            //cout << xx << '\t' << yy << endl;

            //if(! GeomData->polyImm.within(xx, yy) )
            //{
              GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy );

              if(parent == NULL)
              {
                N = NN;
                dN_dx = dNN_dx;
                dN_dy = dNN_dy;
              }
              else
              {
                N = SubDivMat*NN;
                dN_dx = SubDivMat*dNN_dx;
                dN_dy = SubDivMat*dNN_dy;
              }

              //for(ii=0;ii<totnlbf;ii++)
              //printf(" \t %14.8f \n", NN[ii]);

              //printf("\n\n tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n\n\n", xx, yy, val, val, Jac, JacMult, dvol);
              //if(pout) printf(" knotsAtGPs = %12.6f \t xx = %12.6f \t yy = %12.6f \n", knotsAtGPs, xx, yy);
              //printf(" uu = %12.6f \t vv = %12.6f \t dvol = %12.6f \t volume = %12.6f \n", uu, vv, dvol, volume);

              //specVal = DirichletData[aa][2];

              specVal = analy.computeValue(0, geom[0], geom[1]);

              //force[0] = analy.forceX(geom[0], geom[1]);
              //force[1] = analy.forceY(geom[0], geom[1]);

              specVal *= timeFunction[0].prop;

              //if( GeomData->distFuncs[0]->ComputeDistance(xx, yy) > 0.0 )
              //if( xx <= 0.5 )
              //if( rad >= 0.5 )
  
              res = specVal - computeValue(dir, N);

              grad(0) = computeValue(0, dN_dx);
              grad(1) = computeValue(0, dN_dy);

              trac = grad(0)*normal[0] + grad(1)*normal[1] ;

              for(ii=0;ii<totnlbf2;ii++)
              {
                fact = (dvol*PENALTY)* N[ii] ;

                Flocal(ii) += fact*res;

                for(jj=0;jj<totnlbf2;jj++)
                {
                  Klocal(ii, jj) += fact * N[jj];
                }
              }

              // terms corresponding to Nitsche method for Stokes and Navier-Stokes
              if(isNitsche)
              {
                if(dir == 0)
                {
                  for(ii=0;ii<totnlbf2;ii++)
                  {
                    bb1 = N[ii]*dvol;

                    Ta = (mu*dvol)*( dN_dx[ii]*normal[0] + dN_dy[ii]*normal[1] );

                    for(jj=0;jj<totnlbf2;jj++)
                    {
                      Tb  =  mu*(dN_dx[jj]*normal[0] + dN_dy[jj]*normal[1] );

                      Klocal(ii, jj)   -= (bb1*Tb);

                      // Nitsche terms
                      Klocal(ii, jj) -= (Ta*af*N(jj))*NitscheFact;
                    }

                    Flocal(ii)   += (bb1*trac);

                    Flocal(ii)   -= (Ta*res)*NitscheFact;
                  }
                }
              } //if(isNitsche)
           //} // if(! GeomData->polyImm.within(xx, yy) )
        }// for(gp=0...
      } // for(aa=0;aa<DirichletData.size();aa++)
  } // if(DirichletData.size() > 0)

  return;
}




template<>
void TreeNode<2>::applyNeumannBCsCutFEMPoisson(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  //cout << " TreeNode<2>::applyBoundaryConditionsCutFEM() ...... NeumannFlag to be implemented "  << endl;
  if( NeumannData.size() > 0 )
  {
      int ii, jj, aa, gp, nGauss, index, dir, side;
      double  res, JacTemp, dvol, specVal, Jac, rad, freq;

      VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf);
      VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
      myPoint  param, geom;
      double *gws;
      myPoint *gps;

      bool   axsy = ((int)elmDat[2] == 1);


      for(aa=0;aa<NeumannData.size();aa++)
      {
        side    = (int) (NeumannData[aa][0] - 1);
        dir     = (int) (NeumannData[aa][1] - 1);

        //normal = GeomData->boundaryNormals[side];

        if(domNums.size() > 1)
        {
          //TreeNode<2>::setBoundaryGPsCutFEM(side, 0);

          nGauss = BoundaryQuadrature[side].gausspoints.size() ;

          gps = &(BoundaryQuadrature[side].gausspoints[0]);
          gws = &(BoundaryQuadrature[side].gaussweights[0]);

          JacTemp = 1.0;
        }
        else
        {
          nGauss = GeomData->boundaryQuadrature2D[side].gausspoints.size();

          gps = &(GeomData->boundaryQuadrature2D[side].gausspoints[0]);
          gws = &(GeomData->boundaryQuadrature2D[side].gaussweights[0]);

          JacTemp = GeomData->boundaryJacobians[side][level];
        }

        for(gp=0; gp<nGauss; gp++)
        {
            param[1] = 0.5 * (knots[1][2] * gps[gp][1] + knots[1][3]);
            param[0] = 0.5 * (knots[0][2] * gps[gp][0] + knots[0][3]);

            dvol = JacTemp * gws[gp] ;

            geom[0] = GeomData->computeCoord(0, param[0]);
            geom[1] = GeomData->computeCoord(1, param[1]);

            if(axsy)
            {
              if(side != 0)
                dvol *= 2.0*PI*geom[0];
            }

            GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy );

            if(parent == NULL)
            {
              N = NN;
            }
            else
            {
              N = SubDivMat*NN;
            }

            specVal = NeumannData[aa][2];

            res = dvol * specVal * timeFunction[0].prop;

            for(ii=0;ii<totnlbf2;ii++)
              Flocal(ndof*ii+dir) += (res*N(ii));

        }// for(gp=0...
      } // for(aa=0;aa<NeumannData.size();aa++)
  } // if(NeumannData.size() > 0)

  return;
}




template<>
int TreeNode<2>::calcErrorPoisson(int index, int domainCur)
{
    // computing error for Poisson (or Laplace) problem
    // 
    // CUTFEM approach
    ///////////////////////////////////////////////////////////


    int ii, jj, kk, gp, nGauss, tempId, gpDomainId;

    double  JacTemp, Jac, dvol, val;
    double  fact, fact2, bb1, bb2, bb3, mu;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  N, dN_dx, dN_dy, grad(2);
    myPoint  param, geom;

    PoissonEx3  analy;
    //PoissonInterfaceEx4  analy(0.5);

    double *gws;
    myPoint *gps;
    
    if(domNums.size() > 1)
    {
      nGauss = Quadrature.gausspoints.size();
      
      gps = &(Quadrature.gausspoints[0]);
      gws = &(Quadrature.gaussweights[0]);

      tempId = domainCur;

      JacTemp = 1.0;
    }
    else
    {
      nGauss = GeomData->gausspoints.size();

      gps = &(GeomData->gausspoints[0]);
      gws = &(GeomData->gaussweights[0]);

      tempId = domNums[0];

      JacTemp = JacMultElem;
    }

    elemError = 0.0;

  if(index == 1) // L2 norm in displacement
  {
    for(gp=0; gp<nGauss; gp++)
    {
        param[0]  = 0.5*(knots[0][2] * gps[gp][0] + knots[0][3]);
        param[1]  = 0.5*(knots[1][2] * gps[gp][1] + knots[1][3]);

        dvol = gws[gp] * JacTemp;

        //cout << gp << '\t' << gps[gp][0] << '\t' << gps[gp][1] << '\t' << gws[gp] << '\t' << dvol << endl;

        GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

        if(parent == NULL)
        {
          N = NN;
          dN_dx = dNN_dx;
          dN_dy = dNN_dy;
        }
        else
        {
          N = SubDivMat*NN;
          dN_dx = SubDivMat*dNN_dx;
          dN_dy = SubDivMat*dNN_dy;
        }

        geom[0] = GeomData->computeCoord(0, param[0]);
        geom[1] = GeomData->computeCoord(1, param[1]);

        //cout << uu << '\t' << vv << endl;
        //cout << geom[0] << '\t' << geom[1] << endl;

        //cout <<  " gpDomainId = " << gpDomainId << '\t' << totnlbf2 << endl;
        //printVector(GlobalBasisFuncs);

        //if( GeomData->polyImm.within(xx, yy) )
        if(domNums.size() > 1)
        {
          gpDomainId = QuadratureDomNums[gp];
          
        }
        else
        {
          gpDomainId = domNums[0];
        }

        if(gpDomainId == tempId)
        {
          val = analy.computeValue(0, geom[0], geom[1]);

          val -= computeValueCFM(gpDomainId, 0, N);

          fact = val*val;

          elemError += ( fact * dvol );
        }
    }//gp
  }


  if(index == 3) // L2 norm in displacement
  {
    for(gp=0; gp<nGauss; gp++)
    {
        param[0]  = 0.5*(knots[0][2] * gps[gp][0] + knots[0][3]);
        param[1]  = 0.5*(knots[1][2] * gps[gp][1] + knots[1][3]);

        dvol = gws[gp] * JacTemp;

        //cout << gp << '\t' << gps[gp][0] << '\t' << gps[gp][1] << '\t' << gws[gp] << '\t' << dvol << endl;

        GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

        if(parent == NULL)
        {
          N = NN;
          dN_dx = dNN_dx;
          dN_dy = dNN_dy;
        }
        else
        {
          N = SubDivMat*NN;
          dN_dx = SubDivMat*dNN_dx;
          dN_dy = SubDivMat*dNN_dy;
        }

        geom[0] = GeomData->computeCoord(0, param[0]);
        geom[1] = GeomData->computeCoord(1, param[1]);

        //cout << uu << '\t' << vv << endl;
        //cout << geom[0] << '\t' << geom[1] << endl;

        //cout <<  " gpDomainId = " << gpDomainId << '\t' << totnlbf2 << endl;
        //printVector(GlobalBasisFuncs);

        //if( GeomData->polyImm.within(xx, yy) )
        if(domNums.size() > 1)
        {
          gpDomainId = QuadratureDomNums[gp];
          
        }
        else
        {
          gpDomainId = domNums[0];
        }

        if(gpDomainId == tempId)
        {
          val = analy.computeValue(0, geom[0], geom[1]);

          val -= computeValueCFM(gpDomainId, 0, N);

          analy.computeDerivatives(geom[0], geom[1], &(grad(0)));

          grad(0) -= computeValueCFM(gpDomainId, 0, dN_dx);
          grad(1) -= computeValueCFM(gpDomainId, 0, dN_dy);

          fact = val*val + grad[0]*grad[0] + grad[1]*grad[1];

          elemError += ( fact * dvol );
        }
    }//gp
  }

  //if(index == 2)
    //printf(" \t element = %5d ... \t ... error   =   %12.6E \n " , id, elemError);

    return 0;
}










template<>
void TreeNode<3>::calcStiffnessAndResidualCutFEMPoisson(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{

    return;
}



template<>
void TreeNode<3>::applyBoundaryConditionsAtApointCutFEMPoisson(myDataIntegrateCutFEM& myData)
{

   return;
}



template<>
void TreeNode<3>::applyBoundaryConditionsAtApointCutFEMPoisson2(myDataIntegrateCutFEM& myData)
{
  // compute stiffness and force vectors corresponding 
  // to Nitsche method of applying interface conditions
  // 
  // coupling terms


   return;
}





template<>
void TreeNode<3>::applyDirichletBCsCutFEMPoisson(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{

  return;
}




template<>
void TreeNode<3>::applyNeumannBCsCutFEMPoisson(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{

  return;
}




template<>
int TreeNode<3>::calcErrorPoisson(int index, int domainCur)
{

    return 0;
}













