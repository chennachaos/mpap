
#include "TreeNode.h"
#include <string.h>
#include "MpapTime.h"
#include "Functions.h"
#include "FunctionsBiology.h"
#include "BasisFunctionsBSpline.h"
#include "myDataIntegrateCutFEM.h"

extern  MpapTime  mpapTime;


template<>
double TreeNode<1>::getVolume()
{
  volume = 1.0;
  for(int ii=0; ii<1; ii++)
    volume *=  knots[ii][2]* GeomData->getGridLength(ii);

  return  volume;
}


template<>
double TreeNode<2>::getVolume()
{
  volume = 1.0;
  for(int ii=0; ii<2; ii++)
    volume *=  knots[ii][2]* GeomData->getGridLength(ii);

  return  volume;
}


template<>
double TreeNode<3>::getVolume()
{
  volume = 1.0;
  for(int ii=0; ii<3; ii++)
    volume *=  knots[ii][2]* GeomData->getGridLength(ii);

  return  volume;
}


template<>
double TreeNode<1>::getVolumeGaussPoints(int domTemp)
{
  volume = 1.0;
  for(int ii=0; ii<1; ii++)
    volume *=  knots[ii][2]* GeomData->getGridLength(ii);

  return  volume;
}


template<>
double TreeNode<2>::getVolumeGaussPoints(int domTemp)
{
  if(domNums.size() > 1)
  {
    volume = 0.0;
    for(int gp=0; gp<Quadrature.gausspoints.size(); gp++)
    {
      volume += Quadrature.gaussweights[gp];
    }
  }
  else
  {
    volume = 1.0;
    for(int ii=0; ii<2; ii++)
      volume *=  knots[ii][2]* GeomData->getGridLength(ii);
  }

  return  volume;
}


template<>
double TreeNode<3>::getVolumeGaussPoints(int domTemp)
{
  volume = 1.0;
  for(int ii=0; ii<2; ii++)
    volume *=  knots[ii][2]* GeomData->getGridLength(ii);

  return  volume;
}



template<>
void TreeNode<1>::checkPartitionOfUnity()
{
    int ii, jj, gp1, gp2, gp3, count;
    double  JacY, JacZ, Jac, dvol;
    myPoint  param;

    VectorXd  NN(totnlbf), N;

    printf(" elemnum = %6d \t level = %2d  ...  \n ", id, level);
    printVector(GlobalBasisFuncs);
    printMatrix(SubDivMat);

    count = 0;
       
    for(gp1=0;gp1<GeomData->getNGP(0);gp1++)
    {
       param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);

       GeomData->computeBasisFunctions1D(knotBegin, knotIncr, param, NN);

       //NN = GeomData->shpfns[level][count].N;
       //count++;

       if(parent == NULL)
         N  = NN;
       else
         N  = SubDivMat*NN;
       
       //for(ii=0;ii<totnlbf2;ii++)
         //printf("\t %6.4f ", N(ii));

       //printf("\t sum %6.4f \n", N.sum());
       if( !CompareDoubles(N.sum(), 1.0) )
         printf("\t Partition of unity failed for element # = %5d \t sum = %8.4f \n", id, N.sum());
    }//gp1

    //printf("\n\n");

   return;
}



template<>
void TreeNode<2>::checkPartitionOfUnity()
{
    int ii, jj, gp1, gp2, gp3, count;
    double  JacY, JacZ, Jac, dvol;
    myPoint  param;

    VectorXd  NN(totnlbf), N;

    //printf(" elemnum = %6d ...  ", id);

    count = 0;
    for(gp2=0;gp2<GeomData->getNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       
    for(gp1=0;gp1<GeomData->getNGP(0);gp1++)
    {
       param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);

       GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

       //NN = GeomData->shpfns[level][count].N;
       //count++;

       if(parent == NULL)
         N  = NN;
       else
         N  = SubDivMat*NN;

       //for(ii=0;ii<totnlbf2;ii++)
         //printf("\t %6.4f ", N(ii));

       //printf("\t sum %6.4f \n", N.sum());
       if( !CompareDoubles(N.sum(), 1.0) )
         printf("\t Partition of unity failed for element # = %5d \t sum = %8.4f \n", id, N.sum());
    }//gp1
    }//gp2

    //printf("\n\n");

   return;
}




template<>
void TreeNode<3>::checkPartitionOfUnity()
{
    if(parent == NULL)
      return;
    
    int ii, jj, gp1, gp2, gp3, count;
    double  JacY, JacZ, Jac, dvol;
    myPoint  param;

    VectorXd  NN(totnlbf), N;

    //printf(" elemnum = %6d ...  ", id);

    count = 0;
    for(gp3=0;gp3<GeomData->getNGP(2);gp3++)
    {
       param[2]  = 0.5*(knots[2][2] * GeomData->gausspoints3[gp3] + knots[2][3]);
       JacZ = GeomData->gaussweights3[gp3] * JacMultElem;

    for(gp2=0;gp2<GeomData->getNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       JacY = GeomData->gaussweights2[gp2] * JacZ;
       
    for(gp1=0;gp1<GeomData->getNGP(0);gp1++)
    {
       param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
       dvol = GeomData->gaussweights1[gp1] * JacY;

       GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN);

       //NN = GeomData->shpfns[level][count].N;
       //count++;

       if(parent == NULL)
       {
         N  = NN;
       }
       else
       {
         N  = SubDivMat*NN;
       }
       
       //for(ii=0;ii<totnlbf2;ii++)
         //printf("\t %6.4f ", N(ii));

       //printf("\t sum %6.4f \n", N.sum());
       //if( !CompareDoubles(N.sum(), 1.0) )
         //printf("\t Partition of unity failed for element # = %5d \t sum = %8.4f \n", id, N.sum());

    }//gp1
    }//gp2
    }//gp3
    //printf("\n\n", N.sum());

   return;
}



template<>
int TreeNode<1>::calcError(int index, int domainCur)
{
    //HemkerExact  advexact;

    int ii, jj, gp;
   
    double  dvol, xi, val, grad, a, mu, xx, rho1, rho2;
    myPoint  param;

    a    = elmDat[3];
    mu   = elmDat[4];
    rho1 = elmDat[5];
    rho2 = elmDat[6];


    //AdvDiffExact1D  analy(1.0, a, mu, 1.0, 0.0);
    //AdvDiffExact1D  analy;
    //Poisson1DEx1  analy;
    //Fisher1DEx2  analy;
    FK_unsteady_Ex1 analy(mu, rho1);
    //FK_unsteady_Ex2 analy(mu, rho1, rho2);
    //Biharmonic1DEx1  analy;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), N, dN_dx;
    
    elemError = 0.0;
    volume = knots[0][2] * GeomData->getGridLength(0);

    if(index == 0) // L2 error norm
    {
      for(gp=0;gp<GeomData->getNGP(0);gp++)   // loop over Gauss points
      {
        param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp] + knots[0][3]);

        GeomData->computeBasisFunctions1D(knotBegin, knotIncr, param, NN, dNN_dx);

        dvol = GeomData->gaussweights1[gp] * JacMultElem;

        //
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
        //

        xx = GeomData->computeCoord(0, param[0]);

        val = analy.computeValue(0, xx, 0.0);
        val -= computeValue(0, N);

        //printf(" uu, xx, val and soln \t %12.8f \t %12.8f \t %12.8f \t %12.8f \n", uu, xx, analy.computeValue(0, xx, 0.0), computeValue(0, N));
        //printf(" uu, xx, val and soln \t %12.8f \t %12.8f \t %12.8f \t %12.8f \n", uu, xx, dvol, val);

        elemError += ( (val * val) * dvol );
      }
    }
    else if(index == 1) // H1 error norm
    {
      for(gp=0;gp<GeomData->getNGP(0);gp++)   // loop over Gauss points
      {
        param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp] + knots[0][3]);

        GeomData->computeBasisFunctions1D(knotBegin, knotIncr, param, NN, dNN_dx);

        dvol = GeomData->gaussweights1[gp] * JacMultElem;

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

        xx = GeomData->computeCoord(0, param[0]);

        val = analy.computeValue(0, xx, 0.0);
        val -= computeValue(0, N);

        analy.computeDerivatives(xx, 0.0, &grad);
        grad -= computeValue(0, dN_dx);

        //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \t %12.8f \n", computed, exact, diff, xx);

        elemError += ( (val*val + grad*grad) * dvol );
      }
    }
    else if(index == 2) // H2 norm
    {
      VectorXd  d2NN_dx2(totnlbf), d2N_dx2;
      double  Du;

      for(gp=0;gp<GeomData->getNGP(0);gp++)   // loop over Gauss points
      {
        param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp] + knots[0][3]);

        GeomData->computeBasisFunctions1D(knotBegin, knotIncr, param, NN, dNN_dx, d2NN_dx2);

        dvol = GeomData->gaussweights1[gp] * JacMultElem;

          if(parent == NULL)
          {
            N = NN;
            dN_dx = dNN_dx;
            d2N_dx2 = d2NN_dx2;
          }
          else
          {
            N = SubDivMat*NN;
            dN_dx = SubDivMat*dNN_dx;
            d2N_dx2 = SubDivMat*d2NN_dx2;
          }

          xx = GeomData->computeCoord(0, param[0]);

          val = analy.computeValue(0, xx, 0.0);
          val -= computeValue(0, N);

          analy.computeDerivatives(xx, 0.0, &grad);
          grad -= computeValue(0, dN_dx);

          //analy.computeDerivatives2(xx, 0.0, &Du);

          Du -= computeValue(0, d2N_dx2);

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", dx, dy, diff);
        
          elemError += ( (val*val + grad*grad + Du*Du) * dvol );
      }//gp
    }
    else // gradient based error
    {
        for(gp=0;gp<GeomData->getNGP(0);gp++)   // loop over Gauss points
        {
            param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp] + knots[0][3]);

            GeomData->computeBasisFunctions1D(knotBegin, knotIncr, param, NN, dNN_dx);
       
            dvol = GeomData->gaussweights1[gp] * JacMultElem;

            val = computeValue(0, dN_dx);

            elemError += ( val * val * dvol );
        }

       //elemError = sqrt(elemError)/volume;
    }

    //printf(" \t element = %5d ... \t ... error   =   %12.6E \n " , id, elemError);
    //printf("\n\n");

   return 0;
}




/*
template<>
void TreeNode<1>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    // LSFEM for Hemker problem
    //
    ///////////////////////////////////////
    
    int ii, jj, gp;
    double  uu, dvol, xcoord, res, fact;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), D(totnlbf),  NN(GlobalBasisFuncs.size());

    double  val1 = 0.5*knots[0][2];
    double  val2 = 0.5*knots[0][3];

    JacMult = GeomData->getJacobianFull() * val1 ;

    Klocal.setZero();
    Flocal.setZero();
    
    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);
    double  *Data = &(GeomData->soln(0));
    
    for(gp=0;gp<GeomData->getTotalNGP();gp++)   // loop over Gauss points
    {
       uu   = val1*gausspoints[gp] + val2;
       dvol = gaussweights[gp] * JacMult;

       computeBasisFunctions1D(knots[0][0], knots[0][2], uu, &N(0), &dN_dx(0), &d2N_dx2(0));

       xcoord = GeomData->computeCoord(0, uu);
       a = xcoord;

       D = a*dN_dx - mu*d2N_dx2;

       res = 0.0;
       res = mu*PI*PI*cos(PI*xcoord) - PI*xcoord*sin(PI*xcoord) ;
       if(parent == NULL)
       {
          for(ii=0;ii<totnlbf;ii++)
            res -= Data[GlobalBasisFuncs[ii]] * D(ii);
       }
       else
       {
          NN = SubDivMat * D;
          for(ii=0;ii<GlobalBasisFuncs.size();ii++)
            res -= Data[GlobalBasisFuncs[ii]] * NN(ii);
       }

       Klocal += ( (D*dvol)*D.transpose() );
       Flocal += ( (res*dvol) * D);
    }
    //printf("\n\n");
    
    return;
}
*/



/*
template<>
void TreeNode<1>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    // LSFEM for advection-diffusion
    //
    ///////////////////////////////////////

    AdvDiffExact  advexact(3.0, a, mu, s);
    advexact.setBoundaryConditions(0.0, 1.0);

    int ii, jj, gp;
   
    double  knotsAtGPs, Jac, dvol, xi, dummy, x0, xcoord, grad, res, L;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), temp(totnlbf),  NN(GlobalBasisFuncs.size());

    double  incr = knots[0][2];
    double  val1 = 0.5*incr;
    double  val2 = 0.5*knots[0][3];
    
    L = GeomData->getGridLength(0);

    JacMult = GeomData->getJacobianFull() * val1 ;

    Klocal.setZero();
    Mlocal.setZero();
    Flocal.setZero();
    
    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);
    
    for(gp=0;gp<GeomData->getTotalNGP();gp++)   // loop over Gauss points
    {
       knotsAtGPs = val1*gausspoints[gp] + val2;

       computeBasisFunctions1D(knots[0][0], incr, knotsAtGPs, &N(0), &dN_dx(0), &d2N_dx2(0));
       
       dvol = gaussweights[gp] * JacMult;

       temp = a*dN_dx - mu*d2N_dx2 ;

       s = 0.0;
       res = -s ;
       
       Klocal += ( temp*temp.transpose())*dvol;
       Flocal += ( (-res*dvol) * (temp ));
    }
    //printf("\n\n");
    
    return;
}
*/



/*
template<>
void TreeNode<1>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    // GFEM for Hemker problem
    //
    ///////////////////////////////////////

    int ii, jj, gp;
   
    double  knotsAtGPs, tau, Pe, Jac, dvol, xi, x0, xcoord, res, L, h;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), temp(totnlbf);
    
    double  incr = knots[0][2];
    double  val1 = 0.5*incr;
    double  val2 = 0.5*knots[0][3];

    JacMult = GeomData->getJacobianFull() * val1 ;

    Klocal.setZero();
    Flocal.setZero();
    
    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);

    x0 = -1.0;
    L = GeomData->getGridLength(0);
    h = incr*L;

    for(gp=0;gp<GeomData->getTotalNGP();gp++)
    {
       knotsAtGPs = val1 * gausspoints[gp] + val2;

       computeBasisFunctions1D(knots[0][0], incr, knotsAtGPs, &N(0), &dN_dx(0), &d2N_dx2(0));
       
       dvol = gaussweights[gp] * JacMult;

       xcoord = x0 + L * knotsAtGPs;
       a = xcoord;

       if( CompareDoubles(elmDat[4], 0.0) )
         tau = 0.0;
       else
       {
          Pe  =  a*h/(2.0*mu);
          xi  =  1.0/tanh(Pe) - 1.0/Pe;
          tau =  xi*h/(2.0*a);
       }

       // Least-Squares Stabilization term for GLS stabilisation
       temp = a * dN_dx - mu * d2N_dx2 ;

       s = mu*PI*PI*cos(PI*xcoord) - PI*xcoord*sin(PI*xcoord) ;
       res = - s;

       Klocal += ( (a*N + mu*dN_dx) * dN_dx.transpose()  + tau*temp*temp.transpose())*dvol;

       Flocal += ( (-res*dvol) * (N + tau * temp ));
    }
    
    return;
}
*/



template<>
void TreeNode<1>::setInitialProfile()
{
    // to compute the level set function
    // or
    // to set the initial profile for time dependent problems
    ///////////////////////////////////////

    int ii, jj, gp;
    double  dvol, xx, res, mu, rho1, rho2;

    mu   = elmDat[4];
    rho1 = elmDat[5];
    rho2 = elmDat[6];
    
    FK_unsteady_Ex1  analy(mu, rho1);
    //FK_unsteady_Ex2 analy(mu, rho1, rho2);

    VectorXd  NN(totnlbf), N;
    myPoint  param;
    
    ii = forAssyVec.size();

    //Klocal.resize(ii, ii);
    //Flocal.resize(ii);

    //Klocal.setZero();
    //Flocal.setZero();

    for(gp=0; gp<GeomData->gausspoints.size(); gp++)
    {
        param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
        dvol = GeomData->gaussweights[gp] * JacMultElem;

        GeomData->computeBasisFunctions1D(knotBegin, knotIncr, param, NN);

        if(parent == NULL)
        {
          N = NN;
        }
        else
        {
          N = SubDivMat*NN;
        }
       
        xx = GeomData->computeCoord(0, param[0]);

        //u  = computeValue(0, N);

        res = analy.computeValue(0, 0.0, xx);
        //cout << " res = " << res << endl;

        //Klocal += ( (N*dvol)*N.transpose() );
        //Flocal += ( (res*dvol) * N);
    }

    //printMatrix(Klocal);
    //printf("\n\n");
    //printVector(Flocal);

   return;
}

/*
template<>
void TreeNode<1>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    // LSFEM for lower-order Poisson problem
    //
    ///////////////////////////////////////

    int ii, jj, gp, TI, TIp1;
   
    double  knotsAtGPs, Jac, dvol, delta, beta, xx, val, Ddelta=0.0, rr=0.0, fact, sqdvol;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), temp(totnlbf), res(2);
    MatrixXd  D(nsize, ndof);
    
    double  val1 = 0.5*knots[0][2];
    double  val2 = 0.5*knots[0][3];

    JacMult = GeomData->getJacobianFull() * val1 ;
    
    beta = elmDat[3] * GeomData->getGridLength(0) * knots[0][2];
    
    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);

    temp.setZero();
    D.setZero();

    for(gp=0;gp<GeomData->getTotalNGP();gp++)
    {
       knotsAtGPs = val1 * gausspoints[gp] + val2;

       computeBasisFunctions1D(knots[0][0], knots[0][2], knotsAtGPs, &N(0), &dN_dx(0), &d2N_dx2(0));
       
       dvol = gaussweights[gp] * JacMult;

       xx = GeomData->computeCoord(0, knotsAtGPs) - 0.5*GeomData->getGridLength(0) ;

       delta = DiracDelta1(xx, beta);
       
       //delta = computeForce(0,N);
       
       //Klocal += ( (-dvol*dN_dx) * dN_dx.transpose());
       //Flocal += ( dvol*(dN_dx*computeValue(0, dN_dx) + N*(1.0-delta*elmDat[4])));
       //Flocal += ( dvol*(dN_dx*computeValue(0, dN_dx) + N*1.0));

       //delta = computeForce(0,N);

       //if(delta > 0.0)
         //cout << xx << '\t' << delta << endl;
       
       for(ii=0;ii<totnlbf;ii++)
       {
         TI = ndof*ii;
         TIp1 = TI+1;

         D(TI, 0)   = dN_dx(ii);
         D(TIp1, 0)   = -N(ii);
         D(TIp1, 1) = dN_dx(ii);
       }

       res(0) = computeValue(1, N) - computeValue(0, dN_dx);
       res(1) = 1.0 - delta*elmDat[2] - computeValue(1, dN_dx);

       Klocal += ( (dvol*D) * D.transpose());
       Flocal += ( dvol * (D*res) );
    }
    
    return;
}
*/



/*
template<>
void TreeNode<1>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    // LSFEM for interface problem
    //
    ///////////////////////////////////////

    int ii, jj, gp;
   
    double  knotsAtGPs, Jac, dvol, res, delta, beta, xx, val, Ddelta=0.0, rr=0.0, fact, sqdvol;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), temp(totnlbf);

    double  val1 = 0.5*knots[0][2];
    double  val2 = 0.5*knots[0][3];

    JacMult = GeomData->getJacobianFull() * val1 ;
    
    beta = elmDat[3] * GeomData->getGridLength(0) * knots[0][2];
    
    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);

    temp.setZero();

    for(gp=0;gp<GeomData->getTotalNGP();gp++)
    {
       knotsAtGPs = val1 * gausspoints[gp] + val2;

       computeBasisFunctions1D(knots[0][0], knots[0][2], knotsAtGPs, &N(0), &dN_dx(0), &d2N_dx2(0));
       
       dvol = gaussweights[gp] * JacMult;

       xx = GeomData->computeCoord(0, knotsAtGPs) - 0.5*GeomData->getGridLength(0) ;

       delta = DiracDelta1(xx, beta);
       
       //delta = computeForce(0,N);
       
       //Klocal += ( (-dvol*dN_dx) * dN_dx.transpose());
       //Flocal += ( dvol*(dN_dx*computeValue(0, dN_dx) + N*(1.0-delta*elmDat[4])));
       //Flocal += ( dvol*(dN_dx*computeValue(0, dN_dx) + N*1.0));

       //delta = computeForce(0,N);

       //if(delta > 0.0)
         //cout << xx << '\t' << delta << endl;

       res = 1.0;
       //res -= delta*elmDat[4];
       //res += 1.0;
       res -= computeValue(0, d2N_dx2);

       Klocal += ( (dvol*d2N_dx2) * d2N_dx2.transpose());
       Flocal += ( (res*dvol) * d2N_dx2 );
    }
    
    return;
}
*/


/*
template<>
void TreeNode<1>::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
    // LSFEM for unsteady advection-diffusion
    //
    ///////////////////////////////////////

    int ii, jj, gp;
   
    double  knotsAtGPs, Jac, dvol, res, delta, beta, xx, val, fact;
    double  a, mu, dt, am, af, gamma, d1;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), D(totnlbf);

    double  val1 = 0.5*knots[0][2];
    double  val2 = 0.5*knots[0][3];

    JacMult = GeomData->getJacobianFull() * val1 ;
    
    beta = elmDat[3] * GeomData->getGridLength(0) * knots[0][2];
    
    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);

    a = elmDat[1];
    mu  = elmDat[2];

    dt = SolnData->td(0);
    af = SolnData->td(2);
    am = SolnData->td(3);
    gamma = SolnData->td(4);
    d1 = am/af/gamma/dt;


    for(gp=0;gp<GeomData->getTotalNGP();gp++)
    {
       knotsAtGPs = val1 * gausspoints[gp] + val2;

       //computeBasisFunctions1D(knots[0][0], knots[0][2], knotsAtGPs, &N(0), &dN_dx(0), &d2N_dx2(0));
       computeBasisFunctions1D(knots[0][0], knots[0][2], knotsAtGPs, &N(0), &dN_dx(0));
       
       dvol = gaussweights[gp] * JacMult;
       
       D = d1*N + a*dN_dx ;

       res = computeValueDotCur(0,N) + a*computeValueCur(0, dN_dx) ;

       Klocal += ( ((af*dvol)*D) * D.transpose());
       Flocal -= ( (res*dvol) * D );
    }
    
    return;
}
*/


template<>
void TreeNode<1>::computeAndReturnJacobian(int index, double* position, double* normal, double* specVal, double arclen, double* data, MatrixXd& D, VectorXd& vec1, VectorXd& vec2)
{
    // LSFEM for interface problem
    //
    ///////////////////////////////////////

    int ii, jj, gp;
   
    double  Jac, dvol, res, beta, xx, Ddelta=0.0, fact=0.0, delta, fact1, delta2;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), D2(totnlbf);
    myPoint  param;
    
    double  incr = knots[0][2];
    double  val1 = 0.5*knots[0][2];
    double  val2 = 0.5*knots[0][3];

    //JacMult = GeomData->getJacobianFull() * val1 ;
    
    beta = elmDat[3] * GeomData->getGridLength(0) * incr;
    
    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);

    D2.setZero();
    Ddelta = 0.0;
    fact1 = 0.0;

    for(gp=0;gp<GeomData->getNGP(0);gp++)
    {
       param[0] = val1 * gausspoints[gp] + val2;

       GeomData->computeBasisFunctions1D(knotBegin, knotIncr, param, N, dN_dx, d2N_dx2);
       
       dvol = gaussweights[gp] * JacMultElem;

       xx = GeomData->computeCoord(0, param[0]) - position[0];
       delta  = DiracDelta1(xx, beta);
       //delta2 = DiracDelta1(0.0, beta);

       //xx = (GeomData->computeCoord(0, knotsAtGPs) - position[0])/(GeomData->getGridLength(0) * incr);
       //delta  = DiracDelta2(xx, 1.0)/(GeomData->getGridLength(0) * incr);
       //delta2 = DiracDelta2(0.0, 1.0);
       
       //Klocal += ( (dvol*dN_dx) * dN_dx.transpose());
       //Flocal += ( dvol*(-dN_dx*computeValue(0, dN_dx) - N*1.0*computeForce(0,N)));
       //Flocal += ( dvol*(-dN_dx*computeValue(0, dN_dx) - N*1.0));

       delta = computeForce(0, N);

       //cout << xx << '\t' << delta << endl;

       fact = delta * dvol;
       
       Ddelta += (delta*fact);

       D2 += (fact * d2N_dx2);
       fact1 += ((1.0-computeValue(0,d2N_dx2)-delta*data[0])*fact); // 1.0 is the body force. need to generalise this

       //fact = delta*delta*dvol*1.0;
       //Klocal += ( (fact*N)*N.transpose());
       //Flocal += ( ((0.0-computeValue(0, N))*fact) * N);
    }

    //printVector(D2);
    //printf("\n\n");

    /*
    printf("\nAAAAAAAAAAAAAAAAAAAAA\n");
    printVector(temp);
    printf("\n\n");
    printf("\t Ddelta = %12.6f \n", Ddelta);
    printf("\n\n");
    printf("\t fact = %12.6f \n", fact);
    printf("\n\n");
    */

    //printf("\n  elem %5d \n", id);
    //printMatrix(Klocal);
    //printf("\n\n");

    if(parent == NULL)
    {
      //printf("\nAAAAAAAAAAAAAAAAAAAAA\n");
      D.col(0) = D2;
    }
    else
    {
      //printf("\nBBBBBBBBBBBBBBBBBBBB  %5d \t %5d\n", nsize,nsize2);
      D.resize(nsize2,1);
      D.col(0) = SubDivMat*D2;
    }

    vec1(0) = Ddelta;
    vec2(0) = fact1;

    /*
    printf("\nAAAAAAAAAAAAAAAAAAAAA\n");
    printVector(Flocal);
    printf("\n\n");
    printf("\nAAAAAAAAAAAAAAAAAAAAA\n");
    printf("\n\n");
    printVector(temp);
    printf("\n\n");
    printf("\t Ddelta = %12.6f \n", Ddelta);
    printf("\n\n");
    printf("\t rr = %12.6f \n", rr);
    printf("\n\n");
    printf("\nAAAAAAAAAAAAAAAAAAAAA\n");
    printf("\n\n");
    */
    
    return;
}



/*
template<>
void TreeNode<1>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal)
{
    // GFEM for advection-diffusion
    //
    ///////////////////////////////////////

    int ii, jj, gp;
   
    double  uu, tau, Pe, dvol, xi, he, a, force, val;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), tempVec(totnlbf), tempVec2(totnlbf);

    he = GeomData->getGridLength(0) * knots[0][2];
//    cout << volume << '\t' << volume << endl;

    a   = elmDat[3];
    mu  = elmDat[4];

    Pe = a*he/mu/2.0;
    xi = 1.0/sqrt(1.0 + 9.0/Pe/Pe);
    //xi  =  1.0/tanh(Pe) - 1.0/Pe;

    //xi = min(1.0, Pe/3.0/degree[0]/degree[0]);

    tau = he*xi/a/2.0;

    //tau = he*he/(12.0*mu)/degree[0]/degree[0];

    tau *= elmDat[5];

    //cout << " tau " << '\t' << tau << endl;
    
    //Klocal.setZero();
    //Mlocal.setZero();
    //Flocal.setZero();

    for(gp=0;gp<GeomData->getNGP(0);gp++)
    {
       uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp] + knots[0][3]);

       GeomData->computeBasisFunctions1D(knots[0][0], knots[0][2], uu, &N(0), &dN_dx(0), &d2N_dx2(0));
       
       dvol = GeomData->gaussweights1[gp] * JacMult;
       
      
       force = 0.0;
       //s = 10.0*exp(-5.0*knotsAtGPs) - 4.0*exp(-knotsAtGPs) ;


       //tempVec = a*N + mu*dN_dx;

       //val = computeValue(0, dN_dx);

       //Klocal += ( tempVec * dN_dx.transpose() )*dvol;

       //Flocal += ( dvol*( N*force - tempVec*val ) );


       tempVec = mu*dN_dx - a*N;

       val = computeValue(0, tempVec);

       Klocal += ( dN_dx * tempVec.transpose() )*dvol;

       Flocal += ( dvol*( N*force - dN_dx*val ) );

       // stabilisation

       tempVec = a*dN_dx - mu*d2N_dx2 ;

       // Least-Squares Stabilization term for GLS stabilisation
       //tempVec2 = tempVec;

       // SUPG stabilisation
       tempVec2 = a * dN_dx ;

       Klocal += ( (tau*dvol*tempVec2)*tempVec.transpose());

       Flocal += ( (tau*dvol*tempVec2)*( force - computeValue(0, tempVec) ));

    }
    
    return;
}
*/



/*
template<>
void TreeNode<1>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Poisson problem
    //
    ///////////////////////////////////////

    //PoissonInterfaceEx3  analy(1.0, 6.0, 0.0);
    Poisson1DEx1  analy;
    
    int ii, jj, gp, nGauss, tempId, gpDomainId;

    double  uu, Jac, dvol, xx, force, val, mu;
    myPoint  param;

    mu = elmDat[4];

    VectorXd  N(totnlbf), NN(totnlbf), dNN_dx(totnlbf), dN_dx(totnlbf);

    for(gp=0; gp<GeomData->gausspoints.size(); gp++)
    {
        param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
        dvol = GeomData->gaussweights[gp] * JacMultElem;

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
       
        xx = GeomData->computeCoord(0, param[0]);
        
        force = - analy.computeForce(0, xx);
        //force = 0.0;
        //cout << " force = " << force << '\t' << dvol << endl;

        val = mu*computeValue(0, dN_dx);

        Klocal += ( ((dvol*mu)*dN_dx) * dN_dx.transpose());

        Flocal += ( dvol*(N*force - dN_dx*val) );
    }
    
    return;
}
*/


//
template<>
void TreeNode<1>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Fisher-Kolmogorov equation
    //
    ///////////////////////////////////////

    int ii, jj, gp, nGauss, tempId, gpDomainId;

    double  uu, Jac, dvol, xx, force, val, mu, rho1=1.0, rho2=1.0, u, du, udot, af, am, acceFact;
    myPoint  param;

    mu   = elmDat[4];
    rho1 = elmDat[5];
    rho2 = elmDat[6];

    //PoissonInterfaceEx3  analy(1.0, 6.0, 0.0);
    //Poisson1DEx1  analy;
    //Fisher1DEx2  analy;
    FK_unsteady_Ex1 analy(mu, rho1);
    //FK_unsteady_Ex2 analy(mu, rho1, rho2);

    VectorXd  N(totnlbf), NN(totnlbf), dNN_dx(totnlbf), dN_dx(totnlbf);

    af = SolnData->td(2);
    am = SolnData->td(1);
    acceFact = am*SolnData->td(9);

    for(gp=0; gp<GeomData->gausspoints.size(); gp++)
    {
        param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
        dvol = GeomData->gaussweights[gp] * JacMultElem;

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
       
        xx = GeomData->computeCoord(0, param[0]);
        
        force = analy.computeForce(0, xx);
        //force = 0.0;

        u  = computeValueCur(0, N);
        du = computeValueCur(0, dN_dx);
        
        udot = computeValueDotCur(0, N);

        //cout << " force = " << force << '\t' << u << '\t' << du << endl;

        //Klocal += ( dvol*( (mu*af*dN_dx) * dN_dx.transpose() + (acceFact*N)*N.transpose() ) );

        //Flocal += ( dvol*(N*(force-udot) - dN_dx*(mu*du)) );

        Klocal += ( dvol*( (mu*af*dN_dx) * dN_dx.transpose() + (acceFact - (rho1-2.0*rho2*u)*af)*N*N.transpose() ) );

        Flocal += ( dvol*(N*(force-udot+rho1*u-rho2*u*u) - dN_dx*(mu*du)) );
    }

    return;
}
//






template<>
void TreeNode<1>::applyBoundaryConditionsAtApoint(myDataIntegrateCutFEM& myData)
{
   int ii, jj, TI;
   
   double  fact;
   VectorXd  NN(totnlbf), N;

   GeomData->computeBasisFunctions1D(knotBegin, knotIncr, myData.param, NN);

   //for(ii=0;ii<totnlbf;ii++)
     //printf(" \t %14.8f \n", N[ii]);
   
   if(parent == NULL)
     N = NN;
   else
     N = SubDivMat*NN;

  myData.specVal[0] -= computeValue(myData.dir, N);

  for(ii=0;ii<totnlbf2;ii++)
  {
    fact = N[ii] * myData.PENALTY;

    TI = ndof*ii+myData.dir;

    myData.F1(TI) += fact*myData.specVal[0];
    for(jj=0;jj<totnlbf2;jj++)
      myData.K1(TI, ndof*jj+myData.dir) += fact * N[jj];
  }

  return;
}




//
template<>
void TreeNode<1>::applyDirichletBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  // Poisson equation

  if(DirichletData.size() > 0)
  {
    int ii, jj, side, TI, dir, aa;
    double  uu, specVal, PENALTY, fact, NitscheFact, res, dvol, normal, trac;
    double  xx, Ta, Tb, a, mu, rho1, rho2;
    bool isNitsche;
    myPoint  param;

    a    = elmDat[3];
    mu   = elmDat[4];
    rho1 = elmDat[5];
    rho2 = elmDat[6];

    //Fisher1DEx2  analy;
    FK_unsteady_Ex1 analy(mu, rho1);
    //FK_unsteady_Ex2 analy(mu, rho1, rho2);

    //AdvDiffExact1D  analy;
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

        xx = GeomData->computeCoord(0, param[0]);
        
        specVal = analy.computeValue(0, mpapTime.cur, xx);

        res = specVal - computeValue(0, N);

        trac = 0.0 - (mu*computeValue(0, dN_dx))*normal;
        
        //specVal = analy.computeValue(dir, xx, 0.0);

          for(ii=0;ii<totnlbf2;ii++)
          {
            fact = N[ii] * dvol * PENALTY;

            TI = ndof*ii+dir;

            Flocal(TI) += fact*res;

            for(jj=0;jj<totnlbf2;jj++)
            {
              Klocal(TI, ndof*jj+dir) += fact * N[jj];
            }
          }

          if(isNitsche)
          {
            //terms corresponding to Nitsche method for Poisson equation

              for(ii=0;ii<totnlbf2;ii++)
              {
                Ta = dvol*(mu*dN_dx(ii))*normal;
                fact = N(ii)*dvol;

                for(jj=0;jj<totnlbf2;jj++)
                {
                  Tb = (mu*dN_dx(jj))*normal;

                  Klocal(ii, jj)   -= fact*Tb;

                  Klocal(ii, jj)   -= Ta*N(jj)*NitscheFact;
                }

                Flocal(ii)   -= fact*trac;

                Flocal(ii)   -= Ta*res*NitscheFact;
              }
          }

    }
  }

  return;
}
//


/*
template<>
void TreeNode<1>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for extended Fisher-Kolmogorov equation
    //
    ///////////////////////////////////////

    int ii, jj, gp, nGauss;

    double  uu, Jac, dvol, xx, force, val, u, du, d2u, udot, af, am, acceFact;
    myPoint  param;

    double gamm = elmDat[4];
    double mu   = elmDat[4];
    double rho1 = elmDat[5];
    double rho2 = elmDat[6];

    Biharmonic1DEx1  analy;
    //PoissonInterfaceEx3  analy(1.0, 6.0, 0.0);
    //Poisson1DEx1  analy;
    //Fisher1DEx2  analy;
    //FK_unsteady_Ex1 analy(mu, rho1);
    //FK_unsteady_Ex2 analy(mu, rho1, rho2);

    VectorXd  N(totnlbf), NN(totnlbf), dNN_dx(totnlbf), dN_dx(totnlbf), d2NN_dx2(totnlbf), d2N_dx2(totnlbf);

    af = SolnData->td(2);
    am = SolnData->td(1);
    acceFact = am*SolnData->td(9);

    for(gp=0; gp<GeomData->gausspoints.size(); gp++)
    {
        param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
        dvol = GeomData->gaussweights[gp] * JacMultElem;

        GeomData->computeBasisFunctions1D(knotBegin, knotIncr, param, NN, dNN_dx, d2NN_dx2);

        if(parent == NULL)
        {
          N = NN;
          dN_dx = dNN_dx;
          d2N_dx2 = d2NN_dx2;
        }
        else
        {
          N = SubDivMat*NN;
          dN_dx = SubDivMat*dNN_dx;
          d2N_dx2 = SubDivMat*d2NN_dx2;
        }
       
        xx = GeomData->computeCoord(0, param[0]);
        
        force = analy.computeForce(0, xx);
        //force = 0.0;

        u  = computeValueCur(0, N);
        du = computeValueCur(0, dN_dx);
        d2u = computeValueCur(0, d2N_dx2);

        udot = computeValueDotCur(0, N);

        //cout << " force = " << force << '\t' << u << '\t' << du << endl;

        Klocal += ( dvol*( (gamm*af*d2N_dx2) * d2N_dx2.transpose() ) );

        Flocal += ( dvol*(N*force - d2N_dx2*(gamm*d2u)) );

        //Klocal += ( dvol*( (mu*af*d2N_dx2) * d2N_dx2.transpose() + (acceFact*N)*N.transpose() ) );

        //Flocal += ( dvol*(N*(force-udot) - d2N_dx2*(mu*d2u)) );

        //Klocal += ( dvol*( (mu*af*dN_dx) * dN_dx.transpose() + (acceFact - (rho1-2.0*rho2*u)*af)*N*N.transpose() ) );

        //Flocal += ( dvol*(N*(force-udot+rho1*u-rho2*u*u) - dN_dx*(mu*du)) );
    }

    return;
}
*/



/*
template<>
void TreeNode<1>::applyDirichletBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  // Biharmonic equation or
  // extended Fisher-Kolmogorov equation

  if(DirichletData.size() > 0)
  {
    int ii, jj, side, TI, dir, aa;
    double  uu, specVal, PENALTY, fact, NitscheFact, res, dvol, normal, trac;
    double  xx, Ta, Tb;
    bool isNitsche;
    myPoint  param;

    double  gamm = elmDat[3];
    double  mu   = elmDat[4];
    double  rho1 = elmDat[5];
    double  rho2 = elmDat[6];

    double  hx = bbox.maxBB[0]-bbox.minBB[0];
    double  hy = bbox.maxBB[1]-bbox.minBB[1];

    Biharmonic1DEx1  analy;
    //Fisher1DEx2  analy;
    //FK_unsteady_Ex1 analy(mu, rho1);
    //FK_unsteady_Ex2 analy(mu, rho1, rho2);

    //AdvDiffExact1D  analy;
    //AdvDiffExact1D  analy(1.0, a, mu, 1.0, 0.0);

    VectorXd   NN(totnlbf), dNN_dx(totnlbf), N, dN_dx;
    VectorXd  d2NN_dx2(totnlbf), d2N_dx2(totnlbf);
    VectorXd  d3NN_dx3(totnlbf), d3N_dx3(totnlbf);


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

        GeomData->computeBasisFunctions1D(knotBegin, knotIncr, param, NN, dNN_dx, d2NN_dx2, d3NN_dx3 );

        dvol = 1.0;

        if(parent == NULL)
        {
          N = NN;
          dN_dx = dNN_dx;
          d2N_dx2 = d2NN_dx2;
          d3N_dx3 = d3NN_dx3;
        }
        else
        {
          N = SubDivMat*NN;
          dN_dx = SubDivMat*dNN_dx;
          d2N_dx2 = SubDivMat*d2NN_dx2;
          d3N_dx3 = SubDivMat*d3NN_dx3;
        }

        xx = GeomData->computeCoord(0, param[0]);

        //specVal = analy.computeValue(0, mpapTime.cur, xx);
        specVal = analy.computeValue(dir, xx, 0.0);

        res = specVal - computeValue(0, N);

        trac = 0.0 - gamm*(normal*computeValue(0, d3N_dx3));

          for(ii=0;ii<totnlbf2;ii++)
          {
            fact = N[ii] * dvol * (PENALTY*gamm/hx);

            TI = ndof*ii+dir;

            Flocal(TI) += fact*res;

            for(jj=0;jj<totnlbf2;jj++)
            {
              Klocal(TI, ndof*jj+dir) += fact * N[jj];
            }
          }

          if(isNitsche)
          {
            //terms corresponding to Nitsche method for Biharmonic equation

              for(ii=0;ii<totnlbf2;ii++)
              {
                Ta = dvol*(gamm*d3N_dx3(ii))*normal;
                fact = N(ii)*dvol;

                for(jj=0;jj<totnlbf2;jj++)
                {
                  Tb = (gamm*d3N_dx3(jj))*normal;

                  Klocal(ii, jj)   += fact*Tb;

                  Klocal(ii, jj)   += Ta*N(jj)*NitscheFact;
                }

                Flocal(ii)   += fact*trac;

                Flocal(ii)   += Ta*res*NitscheFact;
              }
          }
    }
  }

  return;
}
*/


template<>
void TreeNode<1>::applyNeumannBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  //cout << " TreeNode<1>::applyBoundaryConditions() ...... NeumannFlag to be implemented "  << endl;

  return;
}


template<>
void TreeNode<1>::applyDirichletBCsLSFEM(bool flag, MatrixXd& Klocal, VectorXd& Flocal)
{
  return;
}

template<>
void TreeNode<1>::applyNeumannBCsLSFEM(bool flag, MatrixXd& Klocal, VectorXd& Flocal)
{
  return;
}


template< >
double TreeNode<1>::computeTotalBodyForce(int index, int dir)
{
    return  0.0;
}


/*
template<>
void TreeNode<1>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
//    HemkerExact  advexact(2.0, a, mu, s);
    
//    advexact.setBoundaryConditions(0.0, 1.0);

    int ii, jj, gp;
   
    double  knotsAtGPs, tau, Pe, Jac, dvol, xi, dummy, x0, xcoord, grad;
    double  computed, exact, diff, res, tau1, ahat, L, Wt, h;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), temp(totnlbf),  NN(GlobalBasisFuncs.size());
    
    double  incr = (knots[0][1] - knots[0][0]);
    double  val1 = 0.5*(knots[0][1] - knots[0][0]);
    double  val2 = 0.5*(knots[0][1] + knots[0][0]);

    JacMult = GeomData->getJacobianFull() * val1 ;

    Klocal.setZero();
    Mlocal.setZero();
    Flocal.setZero();
    
    x0 = -1.0;
    //L  = 2.0;
    L = 1.0;
    
    h = incr * L;

    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);

    for(gp=0;gp<GeomData->getTotalNGP();gp++)   // loop over Gauss points
    {
       knotsAtGPs = val1 * gausspoints[gp] + val2;

       computeBasisFunctions1D(knots[0][0], incr, knotsAtGPs, &N(0), &dN_dx(0), &d2N_dx2(0));
       
       dvol = gaussweights[gp] * JacMult;

       //xcoord = x0 + 2.0 * knotsAtGPs;
       //a = xcoord;

       temp = a * dN_dx - mu * d2N_dx2 ;

       if( CompareDoubles(elmDat[4], 0.0) )
          tau = 0.0;
       else
       {
          Pe  =  a*h/(2.0*mu);
          xi  =  1.0/tanh(Pe) - 1.0/Pe;
          tau =  xi*h/(2.0*a);
          //printf(" %14.8f\t %14.8f\t%14.8E\n", xi, tau, Pe);
       }

       //
       //printf("\n");
       res = 0.0;
       if(parent == NULL)
       {
          //printf(" %5d\t %5d\t %14.8f \t %14.8f \t %14.8f \n", ii, LocalBasisFuncs[ii], SolnData->soln(LocalBasisFuncs[ii]), dN_dx(ii), temp(ii));
          for(ii=0;ii<totnlbf;ii++)
            res += SolnData->soln(LocalBasisFuncs[ii]) * temp(ii);
       }
       else
       {
          NN = SubDivMat * temp;
          for(ii=0;ii<GlobalBasisFuncs.size();ii++)
            res += SolnData->soln(GlobalBasisFuncs[ii]) * NN(ii);
       }
       //printf("\n");
       //
       //

       computed = 0.0;
       if(parent == NULL)
       {
          for(ii=0;ii<totnlbf;ii++)
            computed += SolnData->soln(LocalBasisFuncs[ii]) * dN_dx(ii);
       }
       else
       {
          NN = SubDivMat * dN_dx;
          for(ii=0;ii<GlobalBasisFuncs.size();ii++)
            computed += SolnData->soln(GlobalBasisFuncs[ii]) * NN(ii);
       }
       //
       computed = 0.0;
       if(parent == NULL)
       {
          for(ii=0;ii<totnlbf;ii++)
            computed += SolnData->soln(LocalBasisFuncs[ii]) * N(ii);
       }
       else
       {
          NN = SubDivMat * N;
          for(ii=0;ii<GlobalBasisFuncs.size();ii++)
            computed += SolnData->soln(GlobalBasisFuncs[ii]) * NN(ii);
       }
       
       //
       exact = advexact.computeSolutionAt(knotsAtGPs);
       diff  = computed - exact;
       
       //diff *= diff;
       //

       if(counter == 0)
         Wt = 1.0;
       else       
         Wt = 1.0/pow(abs(computed),3.0);
       //

       //if( abs(res) < 1e-6 )         Wt = 0.001;
       //printf(" Wt ... %14.8f \t \t %14.8f \t \t %14.8f \t \t %14.8f \n", res, s, computed, Wt);
       
       //dvol *= Wt;

       //s = mu*PI*PI*cos(PI*xcoord) - PI*xcoord*sin(PI*xcoord) ;
       s = 0.0;
       res = -s ;

       //Klocal += ( (a * N + mu * dN_dx) * dN_dx.transpose()  + tau*temp*temp.transpose())*dvol;
       //Flocal += ( (-res*dvol) * (N + tau * temp ));

       Klocal += ( temp*temp.transpose())*dvol;
       Flocal += ( (-res*dvol) * (temp ));
    }
    
    //printf("\n\n");
    
    counter++;
    
   return;
}
*/




/*
template<>
void TreeNode<1>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    int  ii, jj, gp;
   
    double  knotsAtGPs, tau, Jac, dvol, dt, computed;
    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), temp(totnlbf);

    dt = mpapTime.dt;

    double  incr = (knots[0][1] - knots[0][0]);
    double  val1 = 0.5*(knots[0][1] - knots[0][0]);
    double  val2 = 0.5*(knots[0][1] + knots[0][0]);

    JacMult = GeomData->getJacobianFull() * val1 ;

    Klocal.setZero();
    Mlocal.setZero();
    Flocal.setZero();

    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);

    for(gp=0;gp<GeomData->getTotalNGP();gp++)   // loop over Gauss points
    {
       knotsAtGPs = val1 * gausspoints[gp] + val2;

       computeBasisFunctions1D(knots[0][0], incr, knotsAtGPs, &N(0), &dN_dx(0), &d2N_dx2(0));
       
       dvol = gaussweights[gp] * JacMult;

       computed = 0.0;
       for(ii=0;ii<totnlbf;ii++)
         computed += SolnData->soln(LocalBasisFuncs[ii]) * N(ii);

       temp = N + (dt*a) * dN_dx;

       //Klocal += ( a * temp * dN_dx.transpose() );
       
       Mlocal += (temp * temp.transpose() * dvol);
       Flocal   +=  (temp * (computed*dvol));
    }

  return;
}
*/


/*
template<>
void TreeNode<1>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    int  ii, jj, gp;
   
    double  knotsAtGPs, tau, Jac, dvol, dt, tmp1, tmp2;
    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), temp1(totnlbf), temp2(totnlbf);

    double  incr = (knots[0][1] - knots[0][0]);
    double  val1 = 0.5*(knots[0][1] - knots[0][0]);
    double  val2 = 0.5*(knots[0][1] + knots[0][0]);

    JacMult = GeomData->getJacobianFull() * val1 ;

    if( CompareDoubles(elmDat[4], 0.0) )
       tau = 0.0;
    else
    {
       tau = incr/a/2.0;
    }
    
    Klocal.setZero();
    Mlocal.setZero();
    Flocal.setZero();
    
    dt = mpapTime.dt;
    
    tmp1 = 0.5 * dt * a;

    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);

    for(gp=0;gp<GeomData->getTotalNGP();gp++)   // loop over Gauss points
    {
       knotsAtGPs = val1 * gausspoints[gp] + val2;

       computeBasisFunctions1D(knots[0][0], incr, knotsAtGPs, &N(0), &dN_dx(0), &d2N_dx2(0));
       
       dvol = gaussweights[gp] * JacMult;
       
       temp1 = N + tmp1 * dN_dx;
       temp2 = N - tmp1 * dN_dx;

       Mlocal += ( temp1 * temp1.transpose())*dvol;
       
       Klocal += ( temp1 * temp2.transpose())*dvol;
    }

  return;
}
*/

/*
template<>
void TreeNode<1>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    int  ii, jj, gp;
   
    double  knotsAtGPs, tau, Jac, dvol, dt, tmp1, tmp2;
    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf);

    double  incr = (knots[0][1] - knots[0][0]);
    double  val1 = 0.5*(knots[0][1] - knots[0][0]);
    double  val2 = 0.5*(knots[0][1] + knots[0][0]);

    JacMult = GeomData->getJacobianFull() * val1 ;

    if( CompareDoubles(elmDat[4], 0.0) )
       tau = 0.0;
    else
    {
       tau = incr/a/2.0;
    }
    
    Klocal.setZero();
    Mlocal.setZero();
    Flocal.setZero();
    
    dt = mpapTime.dt;
    
    tmp1 = 0.5 * dt * a * a;
    tmp2 = a*a*a*dt*dt/6.0;
    
    //tmp2 = 0.0;

    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);

    for(gp=0;gp<GeomData->getTotalNGP();gp++)   // loop over Gauss points
    {
       knotsAtGPs = val1 * gausspoints[gp] + val2;

       computeBasisFunctions1D(knots[0][0], incr, knotsAtGPs, &N(0), &dN_dx(0), &d2N_dx2(0));
       
       dvol = gaussweights[gp] * JacMult;

       //Klocal += ( a * N * dN_dx.transpose() + tmp1 * dN_dx * dN_dx.transpose() - tmp2 * dN_dx * d2N_dx2.transpose())*dvol;
       
       Klocal += ( a * N * dN_dx.transpose() - tmp1 * N * d2N_dx2.transpose() - tmp2 * dN_dx * d2N_dx2.transpose())*dvol;
       
       //Klocal += ( a * N * dN_dx.transpose() + tmp1 * dN_dx * dN_dx.transpose() + tmp2 * d2N_dx2 * dN_dx.transpose())*dvol;
       
       Mlocal += (N * N.transpose())*dvol;

       //for(ii=0;ii<totnlbf;ii++)
         //Mlocal(ii,ii) +=  (N(ii) + (tau*a) * dN_dx(ii) )*dvol;
    }

  return;
}
*/



/*
template<>
void TreeNode<1>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    int  ii, jj, gp;
   
    double  knotsAtGPs, tau, Jac, dvol, dt, c;
    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf);
    
    MatrixXd  A(totnlbf, totnlbf), B(totnlbf, totnlbf), C(totnlbf, totnlbf);

    double  incr = (knots[0][1] - knots[0][0]);
    double  val1 = 0.5*(knots[0][1] - knots[0][0]);
    double  val2 = 0.5*(knots[0][1] + knots[0][0]);

    JacMult = GeomData->getJacobianFull() * val1 ;

    if( CompareDoubles(elmDat[4], 0.0) )
       tau = 0.0;
    else
    {
       tau = incr/a/2.0;
    }
    
    Klocal.setZero();
    Mlocal.setZero();
    Flocal.setZero();
    
    dt = mpapTime.dt;
    
    c = 0.5 * dt * a/incr;

    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);

    for(gp=0;gp<GeomData->getTotalNGP();gp++)   // loop over Gauss points
    {
       knotsAtGPs = val1 * gausspoints[gp] + val2;

       computeBasisFunctions1D(knots[0][0], incr, knotsAtGPs, &N(0), &dN_dx(0), &d2N_dx2(0));
       
       dvol = gaussweights[gp] * JacMult;
       
       A = N * N.transpose()*dvol;
       
       B = dN_dx * dN_dx.transpose()*dvol;
       
       C = N * dN_dx.transpose()*dvol;

       Klocal += ( A + 0.5*c*(C.transpose() - C) - (c*c/6.0)*B);
       
       Mlocal += ( A + 0.5*c*(C.transpose() + C) + (c*c/3.0)*B);
    }

  return;
}
*/



template<>
int TreeNode<1>::calcLoadVector(int ind1, int ind2, double inp1, double inp2)
{
    int  ii, jj, gp;
   
    double  dvol, x0, x1, val, mid, computed, diff, Wt, Jac;
    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf);
    myPoint  param;

    double  incr = (knots[0][1] - knots[0][0]);
    double  val1 = 0.5*(knots[0][1] - knots[0][0]);
    double  val2 = 0.5*(knots[0][1] + knots[0][0]);
    
    //x0 = 0.2;    x1 = 0.12;
    //x0 = 0.3; x1 = 0.5;
    x0 = 0.1; x1 = 0.3;
    
    mid = x1-x0;

    //JacMult = GeomData->getJacobianFull() * val1 ;

    //Mlocal.setZero();
    //Flocal.setZero();

    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);

    //cout << " CCCCCCCCCCCC " << GeomData->getTotalNGP() << endl;

    for(gp=0;gp<GeomData->getTotalNGP();gp++)   // loop over Gauss points
    {
       param[0] = val1 * gausspoints[gp] + val2;

       //HB_BasisFuns(degree[0], knots[0][0], incr, knotsAtGPs, &N(0));

       GeomData->computeBasisFunctions1D(knotBegin, knotIncr, param, N, dN_dx, d2N_dx2);

       dvol = gaussweights[gp] * JacMultElem;

       //computed = 0.0;
       //for(ii=0;ii<totnlbf;ii++)
         //computed += SolnData->soln(LocalBasisFuncs[ii]) * dN_dx(ii);

       //printf("\t %12.6f\t %12.6f\t %12.6f \n", x0, x1, knotsAtGPs);
       /*
       if( abs(knotsAtGPs - x0) <= x1 )
         val = 0.5*(1.0 + cos(PI*(knotsAtGPs - x0)/x1));
       else
         val = 0.0;
       /
       if( knotsAtGPs >= x0 && knotsAtGPs <= x1 )
         val = 1.0;
       else
         val = 0.0;
       */
       if( param[0] >= x0 && param[0] <= mid )
         val = 10.0*(param[0]-x0);
       else if( param[0] >= mid && param[0] <= x1 )
         val = 10.0*(x1-param[0]);
       else
         val = 0.0;
       /*

       //diff = computed - val;
          
       if(counter == 0)
         Wt = 1.0;
       else
       {
          Wt = 1.0/pow(abs(computed),2.0);

          //if( abs(diff) < 1e-6 )         Wt = 10000.0;
          //if( CompareDoubles(computed, 0.0) )  Wt = 10000.0;
       }
       dvol *= Wt;
       */

       //printf(" Wt ... %14.8f \t \t %14.8f \t \t %14.8f \t \t %14.8f \n", computed, val, diff, Wt);
       //printf(" Wt ... %14.8f \t \t %14.8f \n", computed, Wt);

       //Mlocal += (N * N.transpose())*dvol;
       //Flocal += (val*dvol) * N ;
    }
    //printf("\n\n");
    //cout << " CCCCCCCCCCCC " << endl;
    
    //counter++;

  return 0;
}




/*
template<>
void TreeNode<1>::calcSubdivisionMatrix()
{
    nlbf2[0] = nlbf[0];
    totnlbf2 = totnlbf;
    nsize2   = nsize;

    if(parent == NULL)
    {
      GlobalBasisFuncs = LocalBasisFuncs;
    }
    else
    {
       int ii, jj, row, col, count=0;
       double *tmp[2];

       TreeNode_PTR  nd;
       
       // find the non-zero basis function numbers for all the elements in the branch starting from present node
    
       GlobalBasisFuncs.clear();
       for(ii=0;ii<LocalBasisFuncs.size();ii++)
       {
          if(LocalBasisFuncs[ii] != -1)
            GlobalBasisFuncs.push_back(LocalBasisFuncs[ii]);
       }
       
       nd = parent;

       while(nd != NULL)
       {
          for(ii=0;ii<nd->LocalBasisFuncs.size();ii++)
          {
             if(nd->LocalBasisFuncs[ii] != -1)
               GlobalBasisFuncs.push_back(nd->LocalBasisFuncs[ii]);
          }
          nd = nd->getParent();
       }
       
       // new number of basis functions

       nlbf2[0] = GlobalBasisFuncs.size();
       totnlbf2 = nlbf2[0];
       nsize2   = totnlbf2 * ndof;

       // size of the matrix 'mat' defined below is totnlbd1xtotnlbf2 and it is correct 
       // because it is take for a single element without considered the refinement

       MatrixXd  mat(totnlbf, totnlbf);

       SubDivMat = MatrixXd::Zero(totnlbf2, totnlbf);
       
       //cout << this->getID() << '\t' << this->orientation << endl;

       if( orientation == LEFT)
          mat = GeomData->coeffLeft;
       else
          mat = GeomData->coeffRight;

       count = 0;
       for(ii=0;ii<LocalBasisFuncs.size();ii++)
       {
          if(LocalBasisFuncs[ii] != -1)
          {
             SubDivMat(count++,ii) = 1.0;
          }
       }

       nd = parent;

       while(nd != NULL)
       {
          for(ii=0;ii<LocalBasisFuncs.size();ii++)
          {
             if(nd->LocalBasisFuncs[ii] != -1)
             {
               SubDivMat.row(count++) = mat.row(ii);
             }
          }
          
          //cout << nd->getID() << '\t' << nd->orientation << endl;
          if( nd->orientation == LEFT)
            mat = GeomData->coeffLeft * mat;
          else
            mat = GeomData->coeffRight * mat;

          nd = nd->getParent();
       }
       
       //printf("\n\n");
       //printMatrix(SubDivMat);       printf("\n\n");
    }

    PROCESSED = true;

    return;
}
*/

    
/*
template<>
void TreeNode<1>::calcSubdivisionMatrix()
{
    nlbf2[0] = nlbf[0];
    totnlbf2 = totnlbf;
    nsize2   = nsize;

    if(parent == NULL)
    {
      GlobalBasisFuncs = LocalBasisFuncs;
    }
    else
    {
      int ii, jj, row, col, count=0, ind1, ind2;

      TreeNode_PTR  nd, nd2;
       
      // find the non-zero basis function numbers for all the elements in the branch starting from present node
    
      GlobalBasisFuncs.clear();
      ind1=0;
      for(ii=0;ii<LocalBasisFuncs.size();ii++)
      {
        if(LocalBasisFuncs[ii] != -1)
        {
          GlobalBasisFuncs.push_back(LocalBasisFuncs[ii]);
          ind1++;
        }
      }

      ind2=0;
      for(ii=0;ii<parent->LocalBasisFuncs.size();ii++)
      {
        if(parent->LocalBasisFuncs[ii] != -1)
        {
          GlobalBasisFuncs.push_back(parent->LocalBasisFuncs[ii]);
          ind2++;
        }
      }
       
      // new number of basis functions

      nlbf2[0] = GlobalBasisFuncs.size();
      totnlbf2 = nlbf2[0];
      nsize2   = totnlbf2 * ndof;

      // size of the matrix 'mat' defined below is totnlbd1xtotnlbf2 and it is correct 
      // because it is take for a single element without considered the refinement

      MatrixXd  mat(totnlbf, totnlbf);
      VectorXd  tmpVec;

      SubDivMat = MatrixXd::Zero(totnlbf2, totnlbf);
       
      //cout << this->getID() << '\t' << this->orientation << endl;

      if( orientation == LEFT)
        mat = GeomData->coeffLeft;
      else
        mat = GeomData->coeffRight;
       
      cout << " id " << id << endl;
      printf("\n\n");
      printMatrix(mat);       printf("\n\n");
      cout << " ind1 " << ind1 << endl;

      //nd = parent;
      for(ii=0;ii<parent->LocalBasisFuncs.size();ii++)
      {
        tmpVec = mat.row(ii);
        //cout << " ii " << ii << endl;
        if(parent->LocalBasisFuncs[ii] == -1)
        {
          ind2 = 0;
          for(jj=0;jj<LocalBasisFuncs.size();jj++)
          {
            //cout << " jj " << jj << endl;
            if(LocalBasisFuncs[jj] != -1)
            {
              SubDivMat(ind2++, jj) += tmpVec(jj);
            }
          }
        }
        else
        {
          //SubDivMat.row(ind1++) = mat.row(ii);
          SubDivMat.row(ind1++) = tmpVec;
        }
      }

       printf("\n\n");
       printMatrix(SubDivMat);       printf("\n\n");
    }

    PROCESSED = true;

    return;
}
*/

    

template<>
void TreeNode<1>::calcSubdivisionMatrix()
{
    nlbf2[0] = nlbf[0];
    totnlbf2 = totnlbf;
    nsize2   = nsize;
    
    LocalBasisFuncsPrev = LocalBasisFuncs;

    if(parent == NULL)
    {
      GlobalBasisFuncs = LocalBasisFuncs;
    }
    else
    {
      int ii, jj, row, col, count=0, ind1, ind2, ind3, size2;

      TreeNode_PTR  nd, nd2;
      
      vector<int>  bfstmp;
       
      // find the non-zero basis function numbers for all the elements in the branch starting from present node
    
      size2 = LocalBasisFuncs.size();
      GlobalBasisFuncs.clear();
      ind1=0;
      for(ii=0;ii<size2;ii++)
      {
        if(LocalBasisFuncs[ii] != -1)
        {
          GlobalBasisFuncs.push_back(LocalBasisFuncs[ii]);
          ind1++;
        }
      }

      nd = parent;
      while(nd != NULL)
      {
        for(ii=0;ii<size2;ii++)
        {
          if(nd->LocalBasisFuncs[ii] != -1)
            GlobalBasisFuncs.push_back(nd->LocalBasisFuncs[ii]);
        }
        nd = nd->getParent();
      }
      cout << " id " << id << endl;
      cout << " GlobalBasisFuncs.size() " << GlobalBasisFuncs.size() << endl;
      printVector(GlobalBasisFuncs);
      printf("\n\n");

      nlbf2[0] = GlobalBasisFuncs.size();
      totnlbf2 = nlbf2[0];

      // size of the matrix 'mat' defined below is totnlbd1xtotnlbf2 and it is correct 
      // because it is take for a single element without considered the refinement

      MatrixXd  mat(totnlbf, totnlbf);
      VectorXd  tmpVec;

      SubDivMat = MatrixXd::Zero(totnlbf2, totnlbf);


      ///////////////////////////////////////////////////////////////
      //
      // Check if the all the local basis functions of the parent element are -1.
      // This happens when enough number of elements adjacent to the parent element are refined.
      // When all the elements at level 'k+1' are refined then the bfs at level 'k' are represented as
      // the linear combination of the bfs at level 'k+2'
      // 
      ///////////////////////////////////////////////////////////////

      for(ii=0;ii<size2;ii++)
      {
        if(parent->LocalBasisFuncsPrev[ii] != -1)
          bfstmp.push_back(parent->LocalBasisFuncs[ii]);
      }

      nd = parent->getParent();
      while(nd != NULL)
      {
        for(ii=0;ii<size2;ii++)
        {
          if(nd->LocalBasisFuncs[ii] != -1)
            bfstmp.push_back(nd->LocalBasisFuncs[ii]);
        }
        nd = nd->getParent();
      }

      cout << " bfstmp.size() " << bfstmp.size() << endl;
      printVector(bfstmp);
      printf("\n\n");

      //cout << this->getID() << '\t' << this->orientation << endl;
      //cout << " flag " << flag << endl;
      cout << " ind1 " << ind1 << endl;
      cout << " totnlbf2 " << totnlbf2 << endl;

      if( orientation == LEFT )
        mat = GeomData->coeffLeft;
      else
        mat = GeomData->coeffRight;

      printMatrix(mat);       printf("\n\n");

      if( parent->getLevel() > 0 )
        mat = parent->SubDivMat*mat;

      printMatrix(mat);       printf("\n\n");

      for(ii=0;ii<bfstmp.size();ii++)
      {
        tmpVec = mat.row(ii);
        cout << " ii " << ii << endl;
        if(bfstmp[ii] == -1)
        {
          ind2 = 0;
          for(jj=0;jj<size2;jj++)
          {
            cout << " jj " << jj << endl;
            if(LocalBasisFuncs[jj] != -1)
            {
              SubDivMat(ind2++, jj) += tmpVec(jj);
            }
          }
        }
        else
        {
          //SubDivMat.row(ind1++) = mat.row(ii);
          SubDivMat.row(ind1++) = tmpVec;
          cout << " bfstmp[ii] != -1 " << endl;
        }
      }

      for(ii=0;ii<SubDivMat.rows();ii++)
      {
        //cout << ii << '\t' << SubDivMat.row(ii).sum() << endl;
        if( CompareDoubles(SubDivMat.row(ii).sum(), 0.0) )
        {
          printf("\t This refinement is not compatible for element # = %5d \n", id);
          printMatrix(SubDivMat);
          exit(-1);
        }
      }

       printf("\n SubDivMat \n");
       printMatrix(SubDivMat);       printf("\n\n");
    }

    PROCESSED = true;
    
    /*
    cout << " elem # " << getID() << endl;
    cout << " nsize2   and  GlobalBasisFuncs " << nsize2 << endl;
    printVector(GlobalBasisFuncs);
    printf("\n\n");
    cout << " forAssyVec " << endl;
    printVector(forAssyVec);
    printf("\n\n");
    */

    return;
}


//
template<>
void TreeNode<2>::calcSubdivisionMatrix()
{
    nlbf2[0] = nlbf[0];
    nlbf2[1] = nlbf[1];
    totnlbf2 = totnlbf;
    nsize2   = nsize;
    
    //cout << this->getID() << '\t' << this->orientation << endl;
    //cout << " level  = " << level << endl;

    if(parent == NULL)
      GlobalBasisFuncs = LocalBasisFuncs;
    else
    {
       int ii, jj, row, col, count=0;

       TreeNode_PTR  nd, nd2;
       
       // find the non-zero basis function numbers for all the elements in the branch starting from present node
       
       GlobalBasisFuncs.clear();
       for(ii=0;ii<LocalBasisFuncs.size();ii++)
       {
          if(LocalBasisFuncs[ii] != -1)
            GlobalBasisFuncs.push_back(LocalBasisFuncs[ii]);
       }
       
       nd = parent;

       while(nd != NULL)
       {
          for(ii=0;ii<nd->LocalBasisFuncs.size();ii++)
          {
             jj = nd->LocalBasisFuncs[ii];
             if(jj != -1)
               GlobalBasisFuncs.push_back(jj);
          }
          nd = nd->getParent();
       }
       
       // new number of basis functions

       totnlbf2 = GlobalBasisFuncs.size();
       nsize2   = totnlbf2 * ndof;

       // size of the matrix 'mat' defined below is totnlbfxtotnlbf and it is correct 
       // because it is take for a single element without considering the refinement

       MatrixXd  mat(totnlbf, totnlbf);

       SubDivMat = MatrixXd::Zero(totnlbf2, totnlbf);

       //cout << this->getID() << '\t' << this->orientation << endl;
       //cout << totnlbf2 << '\t' << nsize2 << endl;

       switch(this->orientation)
       {
          case SW:
            mat = GeomData->coeffSW;
          break;

          case SE:
            mat = GeomData->coeffSE;
          break;

          case NW:
            mat = GeomData->coeffNW;
          break;

          case NE:
            mat = GeomData->coeffNE;
          break;

          default:
            cerr << " TreeNode<2>::calcSubdivisionMatrix() .... WRONG orientation " << this->orientation << " 222 " << endl;
          break;
       }

       count = 0;
       for(ii=0;ii<LocalBasisFuncs.size();ii++)
       {
          if(LocalBasisFuncs[ii] != -1)
             SubDivMat(count++,ii) = 1.0;
       }

       //cout << " aaaaaaaaaaaaaaa " << endl;
       //nd = NULL;
       nd = parent;

       while(nd != NULL)
       {
          for(ii=0;ii<LocalBasisFuncs.size();ii++)
          {
             if(nd->LocalBasisFuncs[ii] != -1)
               SubDivMat.row(count++) = mat.row(ii);
          }
          
          //cout << nd->getID() << '\t' << nd->orientation << endl;

          switch(nd->orientation)
          {
             case ENUM_PARENT:
                  // do nothing as this node is a base node
             break;

             case SW:
                  mat = GeomData->coeffSW * mat;
             break;

             case SE:
                  mat = GeomData->coeffSE * mat;
             break;

             case NW:
                  mat = GeomData->coeffNW * mat;
             break;

             case NE:
                  mat = GeomData->coeffNE * mat;
             break;

             default:
                  cerr << " TreeNode<2>::calcSubdivisionMatrix() .... WRONG orientation " << nd->orientation << endl;
             break;
          }
          nd = nd->getParent();
       }

       
       //printf("\n\n");
       //printMatrix(SubDivMat);       printf("\n\n");
       //printMatrix(SubDivMat2);       printf("\n\n");
    }

    PROCESSED = true;

  return;
}
//


       /*
       if(ndof > 1)
       {
          SubDivMat2.resize(SubDivMat.rows()*ndof, SubDivMat.cols()*ndof);
          MatrixXd  Iden = MatrixXd::Identity(ndof,ndof);
          for(ii=0;ii<SubDivMat.rows();ii++)
          {
             for(jj=0;jj<SubDivMat.cols();jj++)
             {
               SubDivMat2.block(ndof*ii,ndof*jj,ndof,ndof) = SubDivMat(ii,jj)*Iden;
             }
          }
       }
       */

/*
template<>
void TreeNode<2>::calcSubdivisionMatrix()
{
    nlbf2[0] = nlbf[0];
    nlbf2[1] = nlbf[1];
    totnlbf2 = totnlbf;
    nsize2   = nsize;

    LocalBasisFuncsPrev = LocalBasisFuncs;

    if(parent == NULL)
      GlobalBasisFuncs = LocalBasisFuncs;
    else
    {
      int ii, jj, row, col, count=0, ind1, ind2, ind3, size2;

      TreeNode_PTR  nd, nd2;
      
      vector<int>  bfstmp;
       
      // find the non-zero basis function numbers for all the elements in the branch starting from present node
    
      size2 = LocalBasisFuncs.size();
      GlobalBasisFuncs.clear();
      ind1=0;
      for(ii=0;ii<size2;ii++)
      {
        if(LocalBasisFuncs[ii] != -1)
        {
          GlobalBasisFuncs.push_back(LocalBasisFuncs[ii]);
          ind1++;
        }
      }

      nd = parent;
      while(nd != NULL)
      {
        for(ii=0;ii<size2;ii++)
        {
          if(nd->LocalBasisFuncs[ii] != -1)
            GlobalBasisFuncs.push_back(nd->LocalBasisFuncs[ii]);
        }
        nd = nd->getParent();
      }
      //cout << " id " << id << endl;
      //cout << " GlobalBasisFuncs.size() " << GlobalBasisFuncs.size() << endl;
      //printVector(GlobalBasisFuncs);
      //printf("\n\n");

      totnlbf2 = GlobalBasisFuncs.size();

      // size of the matrix 'mat' defined below is totnlbd1xtotnlbf2 and it is correct 
      // because it is take for a single element without considered the refinement

      MatrixXd  mat(totnlbf, totnlbf);
      VectorXd  tmpVec;

      SubDivMat = MatrixXd::Zero(totnlbf2, totnlbf);


      ///////////////////////////////////////////////////////////////
      //
      // Check if the all the local basis functions of the parent element are -1.
      // This happens when enough number of elements adjacent to the parent element are refined.
      // When all the elements at level 'k+1' are refined then the bfs at level 'k' are represented as
      // the linear combination of the bfs at level 'k+2'
      // 
      ///////////////////////////////////////////////////////////////

      for(ii=0;ii<size2;ii++)
      {
        if(parent->LocalBasisFuncsPrev[ii] != -1)
          bfstmp.push_back(parent->LocalBasisFuncs[ii]);
      }

      nd = parent->getParent();
      while(nd != NULL)
      {
        for(ii=0;ii<size2;ii++)
        {
          if(nd->LocalBasisFuncs[ii] != -1)
            bfstmp.push_back(nd->LocalBasisFuncs[ii]);
        }
        nd = nd->getParent();
      }

      //cout << " bfstmp.size() " << bfstmp.size() << endl;
      //printVector(bfstmp);
      //printf("\n\n");

      //cout << this->getID() << '\t' << this->orientation << endl;
      //cout << " flag " << flag << endl;
      //cout << " ind1 " << ind1 << endl;
      //cout << " totnlbf2 " << totnlbf2 << endl;

      switch(orientation)
      {
        case SW:
          mat = GeomData->coeffSW;
        break;

        case SE:
          mat = GeomData->coeffSE;
        break;

        case NW:
          mat = GeomData->coeffNW;
        break;

        case NE:
          mat = GeomData->coeffNE;
        break;

        default:
          cerr << " TreeNode<2>::calcSubdivisionMatrix() .... WRONG orientation " << this->orientation << " 222 " << endl;
        break;
      }

      //printMatrix(mat);       printf("\n\n");

      if( parent->getLevel() > 0 )
        mat = parent->SubDivMat*mat;

      //printMatrix(mat);       printf("\n\n");

      for(ii=0;ii<bfstmp.size();ii++)
      {
        tmpVec = mat.row(ii);
        //cout << " ii " << ii << endl;
        if(bfstmp[ii] == -1)
        {
          ind2 = 0;
          for(jj=0;jj<size2;jj++)
          {
            //cout << " jj " << jj << endl;
            if(LocalBasisFuncs[jj] != -1)
            {
              SubDivMat(ind2++, jj) += tmpVec(jj);
            }
          }
        }
        else
        {
          SubDivMat.row(ind1++) = tmpVec;
          //cout << " bfstmp[ii] != -1 " << endl;
        }
      }

      for(ii=0;ii<SubDivMat.rows();ii++)
      {
        //cout << ii << '\t' << SubDivMat.row(ii).sum() << endl;
        if( CompareDoubles(SubDivMat.row(ii).sum(), 0.0) )
        {
          printf("\t This refinement is not compatible for element # = %5d \n", id);
          printMatrix(SubDivMat);
          exit(-1);
        }
      }

       //printf("\n\n");
       //printMatrix(SubDivMat);       printf("\n\n");
       //printMatrix(SubDivMat2);       printf("\n\n");
    }

    PROCESSED = true;

  return;
}
*/



       /*
       if(ndof > 1)
       {
          //SubDivMat2.resize(SubDivMat.rows()*ndof, SubDivMat.cols()*ndof);
          SubDivMat2 = MatrixXd::Zero(SubDivMat.rows()*ndof, SubDivMat.cols()*ndof);
          //MatrixXd  Iden = MatrixXd::Identity(ndof,ndof);
          for(ii=0;ii<SubDivMat.rows();ii++)
          {
             for(jj=0;jj<SubDivMat.cols();jj++)
             {
               SubDivMat2.block(ndof*ii,ndof*jj,ndof,ndof) = SubDivMat(ii,jj)*Iden;
             }
          }
       }
       //
       if(ndof > 1)
       {
          SubDivMat2.resize(SubDivMat.rows()*ndof, SubDivMat.cols()*ndof);
          SubDivMat2.setZero();
          //SubDivMat2 = MatrixXd::Zero(SubDivMat.rows()*ndof, SubDivMat.cols()*ndof);
          //MatrixXd  Iden = MatrixXd::Identity(ndof,ndof);
          for(ii=0;ii<SubDivMat.rows();ii++)
          {
             for(jj=0;jj<SubDivMat.cols();jj++)
             {
               for(kk=0;kk<ndof;kk++)
               {
                 SubDivMat2(ii*ndof+kk, jj*ndof+kk) = SubDivMat(ii,jj);
               }
             }
          }
       }
       */


//
template<>
void TreeNode<3>::calcSubdivisionMatrix()
{
    nlbf2[0] = nlbf[0];
    nlbf2[1] = nlbf[1];
    nlbf2[2] = nlbf[2];
    totnlbf2 = totnlbf;
    nsize2   = nsize;

    if(parent == NULL)
      GlobalBasisFuncs = LocalBasisFuncs;
    else
    {
       int ii, jj, kk, ll, row, col, count=0;

       TreeNode_PTR  nd;
       
       // find the non-zero basis function numbers for all the elements in the branch starting from present node
       
       GlobalBasisFuncs.clear();
       for(ii=0;ii<LocalBasisFuncs.size();ii++)
       {
          if(LocalBasisFuncs[ii] != -1)
            GlobalBasisFuncs.push_back(LocalBasisFuncs[ii]);
       }
       
       nd = parent;

       while(nd != NULL)
       {
          for(ii=0;ii<nd->LocalBasisFuncs.size();ii++)
          {
             jj = nd->LocalBasisFuncs[ii];
             if(jj != -1)
               GlobalBasisFuncs.push_back(jj);
          }
          nd = nd->getParent();
       }
       
       // new number of basis functions

       totnlbf2 = GlobalBasisFuncs.size();
       nsize2   = totnlbf2 * ndof;

       // size of the matrix 'mat' defined below is totnlbfxtotnlbf and it is correct 
       // because it is take for a single element without considering the refinement

       //cout << totnlbf << '\t' << totnlbf2 << endl;

       MatrixXd  mat(totnlbf, totnlbf);

       SubDivMat = MatrixXd::Zero(totnlbf2, totnlbf);
       
       //cout << this->getID() << '\t' << this->orientation << endl;

       switch(this->orientation)
       {
          case SW_FRONT:
            mat = GeomData->coeffSW_Front;
          break;

          case SE_FRONT:
            mat = GeomData->coeffSE_Front;
          break;

          case NW_FRONT:
            mat = GeomData->coeffNW_Front;
          break;

          case NE_FRONT:
            mat = GeomData->coeffNE_Front;
          break;

          case SW_BACK:
            mat = GeomData->coeffSW_Back;
          break;

          case SE_BACK:
            mat = GeomData->coeffSE_Back;
          break;

          case NW_BACK:
            mat = GeomData->coeffNW_Back;
          break;

          case NE_BACK:
            mat = GeomData->coeffNE_Back;
          break;

          default:
            cerr << " TreeNode<3>::calcSubdivisionMatrix() .... WRONG orientation " << this->orientation << endl;
          break;
       }

       count = 0;
       for(ii=0;ii<LocalBasisFuncs.size();ii++)
       {
          if(LocalBasisFuncs[ii] != -1)
             SubDivMat(count++,ii) = 1.0;
       }

       nd = parent;

       while(nd != NULL)
       {
          for(ii=0;ii<LocalBasisFuncs.size();ii++)
          {
             if(nd->LocalBasisFuncs[ii] != -1)
               SubDivMat.row(count++) = mat.row(ii);
          }
          
          //cout << nd->getID() << '\t' << nd->orientation << endl;

          switch(nd->orientation)
          {
             case ENUM_PARENT:
                  // do nothing as this node is a base node
             break;

             case SW_FRONT:
                  mat = GeomData->coeffSW_Front * mat;
             break;

             case SE_FRONT:
                  mat = GeomData->coeffSE_Front * mat;
             break;

             case NW_FRONT:
                  mat = GeomData->coeffNW_Front * mat;
             break;

             case NE_FRONT:
                  mat = GeomData->coeffNE_Front * mat;
             break;

             case SW_BACK:
                  mat = GeomData->coeffSW_Back * mat;
             break;

             case SE_BACK:
                  mat = GeomData->coeffSE_Back * mat;
             break;

             case NW_BACK:
                  mat = GeomData->coeffNW_Back * mat;
             break;

             case NE_BACK:
                  mat = GeomData->coeffNE_Back * mat;
             break;

             default:
                  cerr << " TreeNode<2>::calcSubdivisionMatrix() .... WRONG orientation " << nd->orientation << endl;
             break;
          }
          nd = nd->getParent();
       }

       //printf("\n\n");
       //printMatrix(SubDivMat);       printf("\n\n");
       //printMatrix(SubDivMat2);       printf("\n\n");
    }

    PROCESSED = true;

    return;
}
//


/*
template<>
void TreeNode<3>::calcSubdivisionMatrix()
{
    nlbf2[0] = nlbf[0];
    nlbf2[1] = nlbf[1];
    nlbf2[2] = nlbf[2];
    totnlbf2 = totnlbf;
    nsize2   = nsize;

    LocalBasisFuncsPrev = LocalBasisFuncs;

    if(parent == NULL)
      GlobalBasisFuncs = LocalBasisFuncs;
    else
    {
      int ii, jj, row, col, count=0, ind1, ind2, ind3, size2;

      TreeNode_PTR  nd, nd2;
      
      vector<int>  bfstmp;
       
      // find the non-zero basis function numbers for all the elements in the branch starting from present node
    
      size2 = LocalBasisFuncs.size();
      GlobalBasisFuncs.clear();
      ind1=0;
      for(ii=0;ii<size2;ii++)
      {
        if(LocalBasisFuncs[ii] != -1)
        {
          GlobalBasisFuncs.push_back(LocalBasisFuncs[ii]);
          ind1++;
        }
      }

      nd = parent;
      while(nd != NULL)
      {
        for(ii=0;ii<size2;ii++)
        {
          if(nd->LocalBasisFuncs[ii] != -1)
            GlobalBasisFuncs.push_back(nd->LocalBasisFuncs[ii]);
        }
        nd = nd->getParent();
      }
      //cout << " id " << id << endl;
      //cout << " GlobalBasisFuncs.size() " << GlobalBasisFuncs.size() << endl;
      //printVector(GlobalBasisFuncs);
      //printf("\n\n");

      totnlbf2 = GlobalBasisFuncs.size();

      // size of the matrix 'mat' defined below is totnlbd1xtotnlbf2 and it is correct 
      // because it is take for a single element without considered the refinement

      MatrixXd  mat(totnlbf, totnlbf);
      VectorXd  tmpVec;

      SubDivMat = MatrixXd::Zero(totnlbf2, totnlbf);


      ///////////////////////////////////////////////////////////////
      //
      // Check if the all the local basis functions of the parent element are -1.
      // This happens when enough number of elements adjacent to the parent element are refined.
      // When all the elements at level 'k+1' are refined then the bfs at level 'k' are represented as
      // the linear combination of the bfs at level 'k+2'
      // 
      ///////////////////////////////////////////////////////////////

      for(ii=0;ii<size2;ii++)
      {
        if(parent->LocalBasisFuncsPrev[ii] != -1)
          bfstmp.push_back(parent->LocalBasisFuncs[ii]);
      }

      nd = parent->getParent();
      while(nd != NULL)
      {
        for(ii=0;ii<size2;ii++)
        {
          if(nd->LocalBasisFuncs[ii] != -1)
            bfstmp.push_back(nd->LocalBasisFuncs[ii]);
        }
        nd = nd->getParent();
      }

      //cout << " bfstmp.size() " << bfstmp.size() << endl;
      //printVector(bfstmp);
      //printf("\n\n");

      //cout << this->getID() << '\t' << this->orientation << endl;
      //cout << " flag " << flag << endl;
      //cout << " ind1 " << ind1 << endl;
      //cout << " totnlbf2 " << totnlbf2 << endl;

      switch(orientation)
      {
          case SW_BACK:
            mat = GeomData->coeffSW_Back;
          break;

          case SE_BACK:
            mat = GeomData->coeffSE_Back;
          break;

          case NW_BACK:
            mat = GeomData->coeffNW_Back;
          break;

          case NE_BACK:
            mat = GeomData->coeffNE_Back;
          break;

          case SW_FRONT:
            mat = GeomData->coeffSW_Front;
          break;

          case SE_FRONT:
            mat = GeomData->coeffSE_Front;
          break;

          case NW_FRONT:
            mat = GeomData->coeffNW_Front;
          break;

          case NE_FRONT:
            mat = GeomData->coeffNE_Front;
          break;

          default:
            cerr << " TreeNode<3>::calcSubdivisionMatrix() .... WRONG orientation " << this->orientation << endl;
          break;
      }

      //printMatrix(mat);       printf("\n\n");

      if( parent->getLevel() > 0 )
        mat = parent->SubDivMat*mat;

      //printMatrix(mat);       printf("\n\n");

      for(ii=0;ii<bfstmp.size();ii++)
      {
        tmpVec = mat.row(ii);
        //cout << " ii " << ii << endl;
        if(bfstmp[ii] == -1)
        {
          ind2 = 0;
          for(jj=0;jj<size2;jj++)
          {
            //cout << " jj " << jj << endl;
            if(LocalBasisFuncs[jj] != -1)
            {
              SubDivMat(ind2++, jj) += tmpVec(jj);
            }
          }
        }
        else
        {
          SubDivMat.row(ind1++) = tmpVec;
          //cout << " bfstmp[ii] != -1 " << endl;
        }
      }

      for(ii=0;ii<SubDivMat.rows();ii++)
      {
        //cout << ii << '\t' << SubDivMat.row(ii).sum() << endl;
        if( CompareDoubles(SubDivMat.row(ii).sum(), 0.0) )
        {
          printf("\t This refinement is not compatible for element # = %5d \n", id);
          cout << " bfstmp.size() " << bfstmp.size() << '\t' << ii << '\t' << mat(ii, 0) << '\t' << SubDivMat(ii,0) << endl;
          printVector(bfstmp);
          printf("\n\n");
          printVector(LocalBasisFuncs);
          printf("\n\n");
          cout << this->getID() << '\t' << this->orientation << endl;
          cout << " ind1 " << ind1 << endl;
          cout << " totnlbf2 " << totnlbf2 << endl;
          printf("\n\n");
          //printMatrix(mat.row(ii));
          //printf("\n\n");
          //printMatrix(SubDivMat.row(ii));
          //printf("\n\n");
          printMatrix(mat);
          printf("\n\n");
          printMatrix(SubDivMat);
          printf("\n\n");
          printf("\n\n");

          for(ii=0;ii<15;ii++)
          {
            tmpVec = mat.row(ii);
            //cout << " ii " << ii << endl;
            if(bfstmp[ii] == -1)
            {
              ind2 = 0;
              for(jj=0;jj<size2;jj++)
              {
                //cout << " jj " << jj << endl;
                if(LocalBasisFuncs[jj] != -1)
                {
                  cout << ii << '\t' << jj << '\t' << ind2++ << '\t' << tmpVec(jj) << endl;
                }
              }
            }
          }

          exit(-1);
        }
      }

       //printf("\n\n");
       //printMatrix(SubDivMat);       printf("\n\n");
       //printMatrix(SubDivMat2);       printf("\n\n");
    }

    PROCESSED = true;

  return;
}
*/




template<>
void TreeNode<1>::mapBoundaryPointDataToGlobalBodyForceVector(double* position, double* normal, double arclen, double* useThisData)
{
    //printVector(Flocal);
    //printf("\n\n");

  return;
}


template<>
void TreeNode<1>::RhsToMapResult(int a1, int b1, double* val)
{
  return;
}


template<>
void TreeNode<1>::MatrixToMapResult(int a1, int b1, SparseMatrixXd& myMat)
{
  
  return;
}






