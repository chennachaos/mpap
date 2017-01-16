
#include "TreeNode.h"
#include "MpapTime.h"
#include "Functions.h"
#include "SolutionData.h"
#include "BasisFunctionsBSpline.h"
#include "myDataIntegrateCutFEM.h"

extern  MpapTime  mpapTime;



template<>
void TreeNode<1>::calcStiffnessAndResidualCutFEMPoisson(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Poisson problem
    //
    ///////////////////////////////////////

    PoissonInterfaceEx3  analy(1.0, 6.0, 0.0);
    
    int ii, jj, gp, nGauss, tempId, gpDomainId;

    double  Jac, dvol, xx, force, val, mu;

    VectorXd  N(totnlbf), NN(totnlbf), dNN_dx(totnlbf), dN_dx(totnlbf);
    
    double *gws;
    myPoint *gps, param;
    
    if(domNums.size() > 1)
    {
      nGauss = Quadrature.gausspoints.size();
      
      gps = &(Quadrature.gausspoints[0]);
      gws = &(Quadrature.gaussweights[0]);
      
      tempId = domainCur;
    }
    else
    {
      nGauss = GeomData->gausspoints.size();

      gps = &(GeomData->gausspoints[0]);
      gws = &(GeomData->gaussweights[0]);
      
      tempId = domNums[0];
    }
    
    cout << " nGauss " << nGauss << endl;

    for(gp=0; gp<nGauss; gp++)
    {
        param[0]   = 0.5*(knots[0][2] * gps[gp][0] + knots[0][3]);
        dvol = gws[gp] * JacMultElem;
        
        cout << gp << '\t' << gps[gp][0] << '\t' << gws[gp] << '\t' << dvol << endl;

        GeomData->computeBasisFunctions1D(knotBegin, knotIncr, param, NN, dNN_dx);

        if(parent == NULL)
        {
          N = NN;
          dN_dx = dNN_dx;
        }
        else
        {
          N = SubDivMat*NN;
          dN_dx = SubDivMat*dNN_dx;
        }
       
        xx = GeomData->ComputeCoord(0, param[0]);
        
        if( xx <= 0.5)
        {
          gpDomainId = 0;
          mu = elmDat[4];
        }
        else
        {
          gpDomainId = 1;
          mu = elmDat[5];
        }
        
        //if(gpDomainId == tempId)
          //dvol *= 1.0;
        //else
          //dvol *= 0.0;

      if(gpDomainId == tempId)
      {
        force = analy.computeForce(0, xx);
        //force = 0.0;
       
        if(tempId == 0)
          val = computeValue(0, dN_dx);
        else
          val = computeValue2(0, dN_dx);

        Klocal += ( ((dvol*mu)*dN_dx) * dN_dx.transpose());

        val *= mu;
        Flocal += ( dvol*(N*force - dN_dx*val) );
      }
    }
    
    return;
}



template<>
void TreeNode<1>::applyBoundaryConditionsAtApointCutFEMPoisson(myDataIntegrateCutFEM& myData)
{
  // compute stiffness and force vectors corresponding 
  // to Nitsche method of applying interface conditions
  // 
  // diagonal terms

    int ii, jj, TI;

    double  fact, fact1, fact2, res, trac, specVal;
    double  Ta1, Ta2, Tb1, Tb2, mu1, mu2, normal1, normal2;
    double  bb1, bb2, u1, u2, t1, t2, trac1, trac2;
    VectorXd   NN(totnlbf), dNN_dx(totnlbf), N, dN_dx;
    
    double  g1, g2, jumpi, jumpj;

    GeomData->computeBasisFunctions1D(knotBegin, knotIncr, myData.param, NN, dNN_dx );

    if(parent == NULL)
    {
      N = NN;
      dN_dx = dNN_dx;
    }
    else
    {
      N = SubDivMat*NN;
      dN_dx = SubDivMat*dNN_dx;
    }

    mu1 = elmDat[4];
    normal1  = 1.0;

    mu2 = elmDat[5];
    normal2 = -normal1;

    //for(ii=0;ii<totnlbf;ii++)
      //printf(" \t %14.8f \n", N[ii]);

    specVal = 0.25/(mu1+mu2);
    //specVal = 0.75;
    
    g1 = g2 = 1.0;
    //g1 = g2 = 0.5;

    jumpi = 0.0;
    jumpj = 0.0;

    u1 = specVal;
    u2 = u1 + jumpi;
    t1 = 0.0;
    t2 = t1 + jumpj;

    u1 -= computeValue(0, N);
    u2 -= computeValue2(0, N);
    
    t1 -= (mu1*computeValue(0, dN_dx) )*normal1;
    t2 -= (mu2*computeValue2(0, dN_dx) )*normal2;

    for(ii=0;ii<totnlbf2;ii++)
    {
        fact  = N[ii] * myData.dvol;
        fact1 = fact * myData.PENALTY ;

        Ta1 = (myData.dvol*myData.NitscheFact)*(mu1*dN_dx(ii))*normal1;
        Ta2 = (myData.dvol*myData.NitscheFact)*(mu2*dN_dx(ii))*normal2;

        fact2 = fact*myData.NitscheFact;

        for(jj=0;jj<totnlbf2;jj++)
        {
          // stabilisation term
          myData.K1(ii, jj) += fact1 * N[jj];
          myData.K2(ii, jj) += fact1 * N[jj];

          // Nitsche terms
          Tb1 = (mu1*dN_dx(jj))*normal1;
          Tb2 = (mu2*dN_dx(jj))*normal2;

          myData.K1(ii, jj)  -= g1*fact2*Tb1;
          myData.K1(ii, jj)  -= g1*Ta1*N(jj);

          myData.K2(ii, jj)  -= g2*fact2*Tb2;
          myData.K2(ii, jj)  -= g2*Ta2*N(jj);

          // coupling terms
          //Kc(ii, jj)  -= fact1 * N[jj];
          //Kc(ii, jj)  -= g1*fact2*Tb1;
          //Kc(ii, jj)  -= g2*Ta2*N(jj);
        }

        // stabilisation terms due to jumpi
        myData.F1(ii) += (fact1*u1);
        myData.F2(ii) += (fact1*u2);

        // Nitsche terms due to jumpi
        myData.F1(ii) -= (g1*Ta1*u1);
        myData.F2(ii) -= (g2*Ta2*u2);

        // Nitsche terms due to jumpj
        myData.F1(ii) -= (g2*(fact2)*t1);
        myData.F2(ii) -= (g1*(fact2)*t2);

        /*
        // Nitsche terms due to jumpi
        myData.F1(ii) += (g1*Ta1*jumpi);
        myData.F2(ii) -= (g2*Ta2*jumpi);

        // stabilisation terms due to jumpi
        myData.F1(ii) -= (fact1*jumpi);
        myData.F2(ii) += (fact1*jumpi);

        // Nitsche terms due to jumpj
        myData.F1(ii) += (g2*(fact2)*jumpj);
        myData.F2(ii) += (g1*(fact2)*jumpj);
        */

    } // for(ii=0;ii<totnlbf2;ii++)

   return;
}



template<>
void TreeNode<1>::applyBoundaryConditionsAtApointCutFEMPoisson2(myDataIntegrateCutFEM& myData)
{
  // compute stiffness and force vectors corresponding 
  // to Nitsche method of applying interface conditions
  // 
  // coupling terms
/*

    int ii, jj, TI;

    double  fact1, fact2, uu, dvol, xx, NitscheFact, Ta1, Ta2, Tb1, Tb2;
    double  mu1, mu2, normal1, normal2;

    VectorXd   NN(totnlbf), dNN_dx(totnlbf), N, dN_dx;
    
    double  g1=0.5, g2=0.5;

    NitscheFact = 1.0;

    uu = param[0];

    GeomData->computeBasisFunctions1D(knots[0][0], knots[0][2], uu, &NN(0), &dNN_dx(0) );

    dvol = 1.0;

    if(parent == NULL)
    {
      N = NN;
      dN_dx = dNN_dx;
    }
    else
    {
      N = SubDivMat*NN;
      dN_dx = SubDivMat*dNN_dx;
    }

    xx = GeomData->ComputeCoord(0, uu);
    
    mu1 = 1.0;
    normal1  = 1.0;

    mu2 = 1.0;
    normal2 = -1.0;


    //for(ii=0;ii<totnlbf;ii++)
      //printf(" \t %14.8f \n", N[ii]);

    for(ii=0;ii<totnlbf2;ii++)
    {
      fact1 = N[ii] * PENALTY * dvol;

      Ta1 = dvol*NitscheFact*(mu1*dN_dx(ii))*normal1;
      Ta2 = dvol*NitscheFact*(mu2*dN_dx(ii))*normal2;

      fact2 = N(ii)*dvol*NitscheFact;

      for(jj=0;jj<totnlbf2;jj++)
      {
        Klocal(ii, jj) -= fact1 * N[jj];

        Tb1 = (mu1*dN_dx(jj))*normal1;
        Tb2 = (mu2*dN_dx(jj))*normal2;

        Klocal(ii, jj)   -= g1*fact2*Tb1;

        Klocal(ii, jj)   += g2*Ta2*N(jj);
      }
    }
*/
   return;
}



template<>
void TreeNode<1>::applyDirichletBCsCutFEMPoisson(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  if(DirichletData.size() > 0)
  {
    int ii, jj, side, TI, dir, aa;
    double  specVal, PENALTY, fact, NitscheFact, res, dvol, normal, trac;
    double  xx, Ta, Tb, a, mu;
    bool isNitsche;
    myPoint  param;
    
    a = elmDat[3];
    mu = elmDat[4];
    
    AdvDiffExact1D  analy;
    //AdvDiffExact1D  analy(1.0, a, mu, 1.0, 0.0);

    VectorXd   NN(totnlbf), dNN_dx(totnlbf), N, dN_dx;

    for(aa=0;aa<DirichletData.size();aa++)
    {
        isNitsche = false;
        side        = (int) (DirichletData[aa][0] - 1);
        dir         = (int) (DirichletData[aa][1] - 1);
        specVal     = DirichletData[aa][2];
        PENALTY     = DirichletData[aa][3];
        isNitsche   = ( (int) DirichletData[aa][4] == 1 );
        NitscheFact = DirichletData[aa][5];
        
        //for symmetric Nitsche method -> NitscheFact = 1.0
        //for unsymmetric Nitsche method -> NitscheFact = -1.0

        if(side == 0)
        {
          param[0] = 0.0;
          normal = -1.0;
        }
        else
        {
          param[0] = 1.0;
          normal = 1.0;
        }

        GeomData->computeBasisFunctions1D(knotBegin, knotIncr, param, NN, dNN_dx );
        
        dvol = 1.0;

        if(parent == NULL)
        {
          N = NN;
          dN_dx = dNN_dx;
        }
        else
        {
          N = SubDivMat*NN;
          dN_dx = SubDivMat*dNN_dx;
        }

        xx = GeomData->ComputeCoord(0, param[0]);
        
        //specVal = analy.computeValue(dir, xx, 0.0);

        if( xx <= 0.5)
        {
          mu = elmDat[4];
          res = specVal - computeValue(0, N);

          //trac = (mu*computeValue(0, dN_dx) - a*computeValue(0, N) )*normal;
          trac = 0.0 - (mu*computeValue(0, dN_dx))*normal;
        }
        else
        {
          mu = elmDat[5];
          res = specVal - computeValue2(0, N);

          //trac = (mu*computeValue2(0, dN_dx) - a*computeValue2(0, N) )*normal;
          trac = 0.0 - (mu*computeValue2(0, dN_dx))*normal;
        }

            if(dir < 4)
            {
              for(ii=0;ii<totnlbf2;ii++)
              {
                fact = N[ii] * dvol * PENALTY;

                TI = ndof*ii+dir;

                Flocal(TI) += fact*res;
                //fact *= af;

                for(jj=0;jj<totnlbf2;jj++)
                {
                  Klocal(TI, ndof*jj+dir) += fact * N[jj];
                }
              }
            }

          if(isNitsche)
          {
            //terms corresponding to Nitsche method for Poisson equation
            if(dir == 0)
            {
              for(ii=0;ii<totnlbf2;ii++)
              {
                //Ta = dvol*(mu*dN_dx(ii)-a*N[ii])*normal;
                Ta = dvol*(mu*dN_dx(ii))*normal;
                fact = N(ii)*dvol;

                for(jj=0;jj<totnlbf2;jj++)
                {
                  //Tb = (mu*dN_dx(jj)-a*N[jj])*normal;
                  Tb = (mu*dN_dx(jj))*normal;

                  Klocal(ii, jj)   -= fact*Tb;

                  Klocal(ii, jj)   -= Ta*N(jj)*NitscheFact;
                }

                Flocal(ii)   -= fact*trac;

                Flocal(ii)   -= Ta*res*NitscheFact;
              }
            } // if(dir == 0)
	  }

    }
  }

  return;
}




template<>
void TreeNode<1>::applyNeumannBCsCutFEMPoisson(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  //cout << " TreeNode<1>::applyBoundaryConditions() ...... NeumannFlag to be implemented "  << endl;

  return;
}



