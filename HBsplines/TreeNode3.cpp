
#include "TreeNode.h"
#include "MpapTime.h"
#include "Functions.h"
#include "FunctionsBiology.h"
#include "TimeFunction.h"

#include "BasisFunctionsBSpline.h"

extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;


template<>
void TreeNode<2>::MatrixToMapResult(int ind1, int ind2, SparseMatrixXd& globalK)
{
    int ii, jj, gp;
   
    double  uu, vv, Jac, dvol, fact;

    MatrixXd  Kl(totnlbf2, totnlbf2);
    VectorXd  NN(totnlbf), N(totnlbf2);
    myPoint  param;

    Kl.setZero();

    for(gp=0; gp<GeomData->gausspoints.size(); gp++)
    {
        param[0] = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
        param[1] = 0.5*(knots[1][2] * GeomData->gausspoints[gp][1] + knots[1][3]);

        dvol = GeomData->gaussweights[gp] * JacMultElem;

        GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

        if(parent == NULL)
          N = NN;
        else
          N = SubDivMat*NN;

        Kl += ((N*dvol)*N.transpose());
    }
    
    //printMatrix(Kl);
    //cout << " AAAAAAAAAA " << endl;

    for(ii=0;ii<totnlbf2;ii++)
    {
      for(jj=0;jj<totnlbf2;jj++)
        globalK.coeffRef(GlobalBasisFuncs[ii], GlobalBasisFuncs[jj]) += Kl(ii,jj);
    }

  return;
}



template<>
void TreeNode<2>::RhsToMapResult(int ind1, int ind2, double* rhs)
{
    int ii, jj, gp;

    double  Jac, dvol, fact, val;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  N(totnlbf2), dN_dx(totnlbf2), dN_dy(totnlbf2), Fl(totnlbf2);
    myPoint  param;

    Fl.setZero();

    for(gp=0; gp<GeomData->gausspoints.size(); gp++)
    {
        param[0]  = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
        param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints[gp][1] + knots[1][3]);

        dvol = GeomData->gaussweights[gp] * JacMultElem;

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
          
          val  = computeValue(1, dN_dx);
          val -= computeValue(0, dN_dy);
          
          val *= dvol;
          
          Fl += ((N*val));
    }

    for(ii=0;ii<totnlbf2;ii++)
      rhs[GlobalBasisFuncs[ii]] += Fl(ii);

  return;
}


/*
template<>
int TreeNode<2>::calcError(int index, int domainCur)
{
    // computing error for Poisson (or Laplace) problem
    ///////////////////////////////////////////////////////////
    
    int      ii, jj, gp1, gp2, count;
    double   uu, vv, Jac, dvol, diff, fact, grad[2], val, xx, yy, r;
    VectorXd   NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), N, dN_dx, dN_dy;
    myPoint  param, geom;

    double  mu = elmDat[4];
    double  rho = elmDat[5];

    //PoissonEx1  analy;
    FK2DsteadyEx1 analy(mu, rho);

    elemError = 0.0;

    //cout << index << '\t' << index << endl;
    //cout << af << endl;

  if(index == 0) // L2 norm
  {
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMultElem;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

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

          geom[0] = GeomData->ComputeCoord(0, param[0]);
          geom[1] = GeomData->ComputeCoord(1, param[1]);

          //val = GeomData->analyDBC->computeValue(0, geom[0], geom[1]);
          val = analy.computeValue(0, geom[0], geom[1]);

          val -= computeValue(0, N);

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f  \t %12.8f   \t %12.8f \n", vx, vx2, vy, vy2, diff);
          
          elemError += ( (val*val) * dvol );
      }//gp1
      }//gp2
 }
 else if(index == 1) // H1 norm
 {
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMultElem;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

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

          geom[0] = GeomData->ComputeCoord(0, param[0]);
          geom[1] = GeomData->ComputeCoord(1, param[1]);

          val = analy.computeValue(0, geom[0], geom[1]);
          analy.computeDerivatives(geom[0], geom[1], grad);

          //val = GeomData->analyDBC->computeValue(0, geom[0], geom[1]);
          //GeomData->analyDBC->computeDerivatives(geom[0], geom[1], grad);

          val -= computeValue(0, N);
          grad[0] -= computeValue(0, dN_dx);
          grad[1] -= computeValue(0, dN_dy);

          fact = val*val + grad[0]*grad[0] + grad[1]*grad[1];

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", dx, dy, diff);
        
          elemError += ( fact * dvol );
      }//gp1
      }//gp2
 }
 else //if(index == 2) // H2 norm
 {
    VectorXd  d2NN_dx2(totnlbf), d2NN_dy2(totnlbf), d2N(totnlbf), d2N_dx2, d2N_dy2;
    double  D2u[3], Du;

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMultElem;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

          //GeomData->computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &NN(0), &dNN_dx(0), &dNN_dy(0), &d2NN_dx2(0), &d2NN_dy2(0));

          if(parent == NULL)
          {
            N = NN;
            dN_dx = dNN_dx;
            dN_dy = dNN_dy;
            d2N_dx2 = d2NN_dx2;
            d2N_dy2 = d2NN_dy2;
          }
          else
          {
            N = SubDivMat*NN;
            dN_dx = SubDivMat*dNN_dx;
            dN_dy = SubDivMat*dNN_dy;
            d2N_dx2 = SubDivMat*d2NN_dx2;
            d2N_dy2 = SubDivMat*d2NN_dy2;
          }
          
          d2N = d2N_dx2 + d2N_dy2;

          geom[0] = GeomData->ComputeCoord(0, param[0]);
          geom[1] = GeomData->ComputeCoord(1, param[1]);

          val = analy.computeValue(0, geom[0], geom[1]);
          analy.computeDerivatives(geom[0], geom[1], grad);
          analy.computeDerivatives2(geom[0], geom[1], D2u);

          //val = GeomData->analyDBC->computeValue(0, geom[0], geom[1]);
          //GeomData->analyDBC->computeDerivatives(geom[0], geom[1], grad);
          //GeomData->analyDBC->computeDerivative2(geom[0], geom[1], D2u);

          val -= computeValue(0, N);
          grad[0] -= computeValue(0, dN_dx);
          grad[1] -= computeValue(0, dN_dy);
          Du = D2u[0] + D2u[1] - computeValue(0, d2N);

          fact = val*val + grad[0]*grad[0] + grad[1]*grad[1] + Du*Du;

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", dx, dy, diff);
        
          elemError += ( fact * dvol );
      }//gp1
      }//gp2
 }

    //printf(" \t element = %5d ... \t ... error   =   %12.6E \n " , id, elemError);

    return 0;
}
*/



/*
template<>
int TreeNode<2>::calcError(int index)
{
    // computing error for Poisson (or Laplace) problem for CUTFEM
    ///////////////////////////////////////////////////////////
    
    int      ii, jj, gp, count;
    double   uu, vv, Jac, dvol, diff, fact, grad[2], val, xx, yy, r;
    VectorXd   NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), N, dN_dx, dN_dy;

    elemError = 0.0;

    //cout << index << '\t' << index << endl;
    //cout << af << endl;

 if(index == 0) // L2 norm
 {
    for(gp=0; gp<GeomData->gausspoints.size(); gp++)
    {
        uu   = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
        vv   = 0.5*(knots[1][2] * GeomData->gausspoints[gp][1] + knots[1][3]);

        dvol = GeomData->gaussweights[gp] * JacMult;

          GeomData->computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &NN(0));

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

          xx = GeomData->ComputeCoord(0, uu);
          yy = GeomData->ComputeCoord(1, vv);
          //r  = sqrt(xx*xx+yy*yy);

          //if(r <= 0.5)
            //val = 1.0;
          //else
            //val = 1.0 + log(2.0*r);

          val = GeomData->analyDBC->computeValue(0, xx, yy);
          
          //val -= computeValue(0, N);
          if( GeomData->polyImm.distance(xx, yy) )
            val -= computeValue2(0, N);

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f  \t %12.8f   \t %12.8f \n", vx, vx2, vy, vy2, diff);
          
          elemError += ( (val*val) * dvol );
    }
 }
 else if(index == 1) // H1 norm
 {
    for(gp=0; gp<GeomData->gausspoints.size(); gp++)
    {
        uu   = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
        vv   = 0.5*(knots[1][2] * GeomData->gausspoints[gp][1] + knots[1][3]);

        dvol = GeomData->gaussweights[gp] * JacMult;

          GeomData->computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &NN(0), &dNN_dx(0), &dNN_dy(0));

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

          xx = GeomData->ComputeCoord(0, uu);
          yy = GeomData->ComputeCoord(1, vv);

          val = GeomData->analyDBC->computeValue(0, xx, yy);
          GeomData->analyDBC->computeDerivatives(xx, yy, grad);

          if( GeomData->polyImm.distance(xx, yy) )
          {
            val -= computeValue(0, N);
            grad[0] -= computeValue(0, dN_dx);
            grad[1] -= computeValue(0, dN_dy);
          }

          fact = val*val + grad[0]*grad[0] + grad[1]*grad[1];

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", dx, dy, diff);
        
          elemError += ( fact * dvol );
    }
 }
 else //if(index == 2) // H2 norm
 {
    VectorXd  d2NN_dx2(totnlbf), d2NN_dy2(totnlbf), d2N(totnlbf), d2N_dx2, d2N_dy2;
    double  D2u[3], Du;

    for(gp=0; gp<GeomData->gausspoints.size(); gp++)
    {
        uu   = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
        vv   = 0.5*(knots[1][2] * GeomData->gausspoints[gp][1] + knots[1][3]);

        dvol = GeomData->gaussweights[gp] * JacMult;

          GeomData->computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &NN(0), &dNN_dx(0), &dNN_dy(0), &d2NN_dx2(0), &d2NN_dy2(0));

          if(parent == NULL)
          {
            N = NN;
            dN_dx = dNN_dx;
            dN_dy = dNN_dy;
            d2N_dx2 = d2NN_dx2;
            d2N_dy2 = d2NN_dy2;
          }
          else
          {
            N = SubDivMat*NN;
            dN_dx = SubDivMat*dNN_dx;
            dN_dy = SubDivMat*dNN_dy;
            d2N_dx2 = SubDivMat*d2NN_dx2;
            d2N_dy2 = SubDivMat*d2NN_dy2;
          }
          
          d2N = d2N_dx2 + d2N_dy2;

          xx = GeomData->ComputeCoord(0, uu);
          yy = GeomData->ComputeCoord(1, vv);

          //val = analy.computeValue(0, xx, yy);
          //analy.computeDerivatives(xx, yy, grad);
          //analy.computeDerivatives2(xx, yy, D2u);

          val -= computeValue(0, N);
          grad[0] -= computeValue(0, dN_dx);
          grad[1] -= computeValue(0, dN_dy);
          Du = D2u[0] + D2u[1] - computeValue(0, d2N);

          fact = val*val + grad[0]*grad[0] + grad[1]*grad[1] + Du*Du;

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", dx, dy, diff);
        
          elemError += ( fact * dvol );
    }
 }

    //printf(" \t element = %5d ... \t ... error   =   %12.6E \n " , id, elemError);

    return 0;
}
*/





template< >
double TreeNode<2>::computeTotalBodyForce(int index, int dir)
{
    int      ii, jj, gp, count;
    double   uu, vv, Jac, dvol, totalVal, fact, val, xx, yy;
    VectorXd   NN(totnlbf), N;
    myPoint  param;

    totalVal = 0.0;
    count = 0;
    if(index == 0)
    {
       for(gp=0; gp<GeomData->gausspoints.size(); gp++)
       {
         param[0] = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
         param[1] = 0.5*(knots[1][2] * GeomData->gausspoints[gp][1] + knots[1][3]);

         dvol = GeomData->gaussweights[gp] * JacMultElem;

             GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

             if(parent == NULL)
               N = NN;
             else
               N = SubDivMat*NN;

             val = computeForce(dir, N);

             //printf(" val \t %12.8f \n", val);
        
             totalVal += ( val * dvol );
       }
    }
    if(index == 1)
    {
       for(gp=0; gp<GeomData->gausspoints.size(); gp++)
       {
         param[0]  = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
         param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints[gp][1] + knots[1][3]);

         dvol = GeomData->gaussweights[gp] * JacMultElem;

             //if(parent == NULL)
               //N = GeomData->shpfns[level][count].N;
             //else
               //N = SubDivMat*GeomData->shpfns[level][count].N;

             count++;

             val = computeForcePrev(dir, N);

             totalVal += ( val * dvol );
       }
    }
    if(index == 2)
    {
       for(gp=0; gp<GeomData->gausspoints.size(); gp++)
       {
         param[0]  = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
         param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints[gp][1] + knots[1][3]);

         dvol = GeomData->gaussweights[gp] * JacMultElem;

             GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

             if(parent == NULL)
               N = NN;
             else
               N = SubDivMat*NN;

             val = computeForceCur(dir, N);
        
             totalVal += ( val * dvol );
       }
    }
    //cout << totalVal << endl;

    return  totalVal;
}


//
template<>
int TreeNode<2>::calcError(int index, int domainCur)
{
    // computing error for Navier-Stokes
    ///////////////////////////////////////////////////////////
    
    //Stokes2DEx1  analy;
    //Stokes2DEx2  analy;

    //Stokes2DEx3  analy;
    Kovasznay  analy;
    //analy.SetPressure(0.0);

    //PearsonVortex analy;

    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    //KimMoinFlow  analy(rho, mu);
    //Kovasznay  analy;

    int      ii, jj, gp1, gp2, TI;
    double   Jac, dvol, diff, fact, val;
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  N, dN_dx, dN_dy;
    myPoint  param, geom;
    
    elemError = 0.0;

 if(index < 3) // L2 norm in x-velocity (index=0), y-velocity (index=1) and pressure (index=2)
 {
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMultElem;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

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

          geom[0] = GeomData->ComputeCoord(0, param[0]);
          geom[1] = GeomData->ComputeCoord(1, param[1]);

          val = analy.computeValue(index, geom[0], geom[1]);
          //val *= exp(-2.0*mu*timeFunction[0].prop);
          val -= computeValue(index, N);

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f  \t %12.8f   \t %12.8f \n", vx, vx2, vy, vy2, diff);
          
          elemError += ( (val*val) * dvol );
    }//gp1
    }//gp2
 }
 else // H1 semi-norm in velocity
 {
    MatrixXd  F(2,2);
    double  v[2];

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMultElem;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

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

          geom[0] = GeomData->ComputeCoord(0, param[0]);
          geom[1] = GeomData->ComputeCoord(1, param[1]);

          v[0] = analy.computeValue(0, geom[0], geom[1]);
          v[1] = analy.computeValue(1, geom[0], geom[1]);

          v[0] -= computeValue(0, N);
          v[1] -= computeValue(1, N);

          F.setZero();
          analy.computeDerivatives(geom[0], geom[1], &(F(0,0)));

          F(0,0) -= computeValue(0, dN_dx);
          F(0,1) -= computeValue(0, dN_dy);
          F(1,0) -= computeValue(1, dN_dx);
          F(1,1) -= computeValue(1, dN_dy);

          val = v[0]*v[0]+v[1]*v[1] + F(0,0)*F(0,0)+F(0,1)*F(0,1)+F(1,0)*F(1,0)+F(1,1)*F(1,1);

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", dx, dy, diff);
        
          elemError += ( val * dvol );
    }//gp1
    }//gp2
 }
  //if(index == 2)
    //printf(" \t element = %5d ... \t ... error   =   %12.6E \n " , id, elemError);

    return 0;
}
//




template<>
void TreeNode<2>::setInitialProfile()
{
    // to compute the level set function
    // or
    // to set the initial profile for time dependent problems
    ///////////////////////////////////////

    //Circle  circ1(0.5,0.5,0.25);
    //Circle  circ2(0.5,0.5,0.1);
    //Circle  circ3(0.4,0.6,0.15);
    //Circle  circ4(0.6,0.6,0.15);

    PearsonVortex analy(elmDat[4]);

    int ii, jj, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;
    double  Jac, dvol, fact, fact1, fact2, fact3, y0, y1, xx, yy, y2;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), dN_dy(totnlbf), d2N_dy2(totnlbf);
    //double  *N, *dN_dx, *d2N_dx2, *dN_dy, *d2N_dy2;

    y0 = 0.25;
    y1 = 0.75;
    y2 = 24.0;

    //y0 = 0.0;
    //y1 = 1.0;
    //y2 = 6.0;

    VectorXd  res(3);
    MatrixXd  D(nsize,3); D.setZero();
    myPoint  param;

    Klocal.setZero();
    Flocal.setZero();

    for(gp=0; gp<GeomData->gausspoints.size(); gp++)
    {
        param[0] = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
        param[1] = 0.5*(knots[1][2] * GeomData->gausspoints[gp][1] + knots[1][3]);

        dvol = GeomData->gaussweights[gp] * JacMultElem;

          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, N, dN_dx, dN_dy, d2N_dx2, d2N_dy2);
          //computeBasisFunctions2D(uu, vv, N, dN_dx, dN_dy, d2N_dx2, d2N_dy2);
          
          //if(circ1.checkPointLocation(uu, vv) & !circ2.checkPointLocation(uu, vv)) // || circ3.checkPointLocation(uu, vv) || circ4.checkPointLocation(uu, vv))
            //fact = 1.0;
          //else
            //fact = 0.0;

          //fact = (uu-0.5)*(uu-0.5) + (vv-0.5)*(vv-0.5);
          //printf(" \t %14.8f \t %14.8f \t %14.8f\n", f(0), f(1), f(2));

          //res(0) = analy.computeXVelocity(uu, vv);
          //res(1) = analy.computeYVelocity(uu, vv);
          //res(2) = analy.computePressure(uu, vv);

          xx = GeomData->ComputeCoord(0, param[0]);
          yy = GeomData->ComputeCoord(1, param[1]);

          if(yy <= y1 && yy >= y0)
            res(0) = y2*(yy-y0)*(y1-yy);

          for(ii=0;ii<totnlbf;ii++)
          {
             TI = 3*ii;

             D(TI,0)   = N[ii];
             D(TI+1,1) = N[ii];
             D(TI+2,2) = N[ii];
          }

          Klocal += ((D*dvol) * D.transpose());
          Flocal += (D*(dvol*res));
    } //gp

    //printf("\n\n");
    //printVector(Flocal);
    
   return;
}


/*
template<>
void TreeNode<2>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    // LSFEM for Advection-Diffusion
    //
    //////////////////////////////////////////////////
  
    int ii, jj, gp1, gp2;

    double  uu, vv, Jac, dvol, fact, res;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), dN_dy(totnlbf), d2N_dy2(totnlbf), D(totnlbf), NN;

    Klocal.setZero();
    Flocal.setZero();
    
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMult;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          //computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &N(0), &dN_dx(0), &dN_dy(0), &d2N_dx2(0), &d2N_dy2(0));
          computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &N(0), &dN_dx(0), &dN_dy(0));

          //printf(" uu and vv \t\t %12.8f \t\t %12.8f \t\t %12.8f \n ", uu, vv, dvol);

          //for(ii=0;ii<nsize;ii++)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);

          D = dN_dx + dN_dy;

          //res = sin(uu + vv)*dvol;
          res = 0.0;
          res -= computeValue(0, D);

          Klocal += ( (D*dvol) * D.transpose());
          Flocal += ( (res*dvol) * D );
       }
    }
    //printMatrix(Klocal);
    //printf("\n\n");
    //printVector(Flocal);
    //printf("\n\n");

    return;
}
*/



/*
template<>
void TreeNode<2>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    // LSFEM for Poisson's (or Laplace) problem
    ///////////////////////////////////////

    int      ii, jj, gp1, gp2;
    double   uu, vv, Jac, dvol, fact, res, xx, yy, r, beta, betax, betay, b;
    //double   N[totnlbf], dN_dx[totnlbf], d2N_dx2[totnlbf], dN_dy[totnlbf], d2N_dy2[totnlbf];
    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), dN_dy(totnlbf), d2N_dy2(totnlbf), D(totnlbf), NN;

    a = elmDat[3];

    Klocal.setZero();
    Flocal.setZero();
    
    b = -3.0;
    //b = 10.0;

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMult;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &N(0), &dN_dx(0), &dN_dy(0), &d2N_dx2(0), &d2N_dy2(0));

          //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", uu, vv, fact, dvol0);
          //printf("BasisFuns \n");
          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %12.6f  \t %12.6f  \t %12.6f \n ", N(ii), dN_dx(ii), dN_dy(ii));

          xx = GeomData->ComputeCoord(0, uu);
          yy = GeomData->ComputeCoord(1, vv);
          r  = sqrt(xx*xx+yy*yy);

          if(r <= 0.5)
          {
            beta = r*r + 1.0;
            betax = 2.0*xx;
            betay = 2.0*yy;
            //res = 8.0*r*r + 4.0;
          }
          else
          {
            beta = b;
            betax = 0.0;
            betay = 0.0;
            //res = 0.0;
          }
          //res = 8.0*r*r + 4.0;

          beta = 1.0;
          betax = betay = res = 0.0;

          D = betax*dN_dx + betay*dN_dy + beta*(d2N_dx2 + d2N_dy2);

          res += computeForce(0, N);
          //cout << " res " << res << endl;
          res -= computeValue(0, D);
          //res += 8.0*r*r + 4.0;

          Flocal += ((dvol*res)*D);
          Klocal += ((dvol*D)*D.transpose());
    }//gp1
    }//gp2

    //printMatrix(Klocal);
    //printf("\n\n");
    //printVector(Flocal);
    
   return;
}
*/


/*
template<>
void TreeNode<2>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Poisson's (or Laplace) problem
    ///////////////////////////////////////

    //PoissonEx1 analy;
    PoissonEx3 analy;

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;
   
    double  Jac, dvol, force;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  N, dN_dx, dN_dy;

    myPoint  param, geom;

    double  mu  = elmDat[4];

    //cout << " uuuuuuuuuuuuuuu " << endl;
    //cout << GeomData->GetNGP(1) << '\t' << GeomData->GetNGP(0) << endl;
    //printVector(GeomData->gausspoints1);
    //printVector(GeomData->gausspoints2);
    //printVector(GeomData->gaussweights1);
    //printVector(GeomData->gaussweights2);
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
      param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
      Jac = GeomData->gaussweights2[gp2] * JacMultElem;

      for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
      {
          param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          //cout << " uuuuuuuuuuuuuuu " << gp2 << '\t' << gp1 << endl;

          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

          //cout << " AAAAAAAAAAA " << endl;

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

          geom[0] = GeomData->ComputeCoord(0, param[0]);
          geom[1] = GeomData->ComputeCoord(1, param[1]);

          force = analy.computeForce(0, geom[0], geom[1]);
          //fact = 0.0;

        Klocal += ((dvol*mu)*(dN_dx*dN_dx.transpose()+dN_dy*dN_dy.transpose()));

        //cout << " AAAAAAAAAAA " << endl;
        Flocal += (dvol*(N*force - dN_dx*(mu*computeValue(0,dN_dx)) - dN_dy*(mu*computeValue(0,dN_dy))) );
        //cout << " AAAAAAAAAAA " << endl;
    }//gp1
    }//gp2

  //printMatrix(Klocal);
  //printVector(Flocal);

  return;
}
*/


/*
template<>
void TreeNode<2>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Fisher's equation
    ///////////////////////////////////////

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;
   
    double  Jac, dvol, force, u, du[2];

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  N, dN_dx, dN_dy;

    myPoint  param, geom;


    double  mu   = elmDat[4];
    double  rho  = elmDat[5];

    //PoissonEx1 analy;
    //PoissonEx3 analy;
    FK2DsteadyEx1 analy(mu, rho);


    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMultElem;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

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

          geom[0] = GeomData->ComputeCoord(0, param[0]);
          geom[1] = GeomData->ComputeCoord(1, param[1]);
          
          u     = computeValue(0,N);
          du[0] = computeValue(0,dN_dx);
          du[1] = computeValue(0,dN_dy);

          force = analy.computeForce(0, geom[0], geom[1]);
          force *= timeFunction[0].prop;

          //fact = 0.0;

        Klocal += ( dvol*((mu*dN_dx)*dN_dx.transpose() + (mu*dN_dy)*dN_dy.transpose() ) );
        Klocal += ( dvol*(rho*(3.0*u*u-1.0)*N)*N.transpose() ) ;

        //cout << " AAAAAAAAAAA " << endl;
        Flocal += (dvol*(N*force - dN_dx*(mu*du[0]) - dN_dy*(mu*du[1]) ) );
        Flocal -= (dvol*(N*rho*(u*u*u-u)) );
        //cout << " AAAAAAAAAAA " << endl;

    }//gp1
    }//gp2

   return;
}
*/


/*
template<>
void TreeNode<2>::applyDirichletBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  // Nitsche method for Poisson equation

  if( DirichletData.size() > 0 )
  {
    int ii, jj, aa, gp1, gp2, TI, TIp1, TIp2, index, TJ, TJp1, TJp2, dir, side;
    double  Ta, Tb, res, JacMultLoc, af, NitscheFact;
    double  dvol, specVal, PENALTY, Jac, fact, trac;
    bool  isNitsche;

    myPoint  param, geom, normal;
    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf);
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);

    vector<double>  boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2;

    bool   axsy = ((int)elmDat[2] == 1);

    double  mu  = elmDat[4];
    double  rho = elmDat[5];

    double  hx = bbox.maxBB[0]-bbox.minBB[0];
    double  hy = bbox.maxBB[1]-bbox.minBB[1];


    af = SolnData->td(2);

    //TwoDim_Ex1  analy;
    //PoissonEx2 analy;
    //PoissonEx1 analy;
    PoissonEx3 analy;
    //FK2DsteadyEx1 analy(mu, rho);

      for(aa=0;aa<DirichletData.size();aa++)
      {
       // printVector(DirichletData[aa]);

        isNitsche = false;
        side        = (int) (DirichletData[aa][0] - 1);
        dir         = (int) (DirichletData[aa][1] - 1);
        specVal     = DirichletData[aa][2];
        PENALTY     = DirichletData[aa][3];
        isNitsche   = ( (int) DirichletData[aa][4] == 1 );
        NitscheFact = DirichletData[aa][5];

        //for symmetric Nitsche method -> NitscheFact = 1.0
        //for unsymmetric Nitsche method -> NitscheFact = -1.0

        //cout << " PENALTY " << PENALTY << endl;

        GeomData->getBoundaryNormal2D(side, normal);
        GeomData->setBoundaryGPs2D(side, boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2);

        JacMultLoc = TreeNode<2>::getJacBoundary(side);
        
        for(gp2=0;gp2<boundaryGPs2.size();gp2++)
        {
            param[1] = 0.5 * (knots[1][2] * boundaryGPs2[gp2] + knots[1][3]);
        for(gp1=0;gp1<boundaryGPs1.size();gp1++)
        {
            param[0] = 0.5 * (knots[0][2] * boundaryGPs1[gp1] + knots[0][3]);
 
            dvol = JacMultLoc * boundaryGWs2[gp2] * boundaryGWs1[gp1] ;
            
            //printf(" %4d \t %4d \t %12.6f \t %12.6f \n", side, dir, vv, uu);

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

            geom[0] = GeomData->ComputeCoord(0, param[0]);
            geom[1] = GeomData->ComputeCoord(1, param[1]);
            //cout << xx << '\t' << yy << '\t' << DirichletData[aa][2] << endl;
            //cout << aa << '\t' << side << '\t' << dir << '\t' << val << endl;

            specVal = DirichletData[aa][2];

            specVal = analy.computeValue(dir, geom[0], geom[1]);
            //cout << side << '\t' << dir << '\t' << specVal << endl;

            //cout << " lllllllllll " << endl;
            specVal *= timeFunction[0].prop;
            //cout << " lllllllllll " << endl;

            //cout << aa << '\t' << side << '\t' << dir << '\t' << res << endl;
            //cout << specVal << '\t' << computeValue(dir, N) << endl;
            //cout << " PENALTY = " << PENALTY << endl;

            res = specVal - computeValue(0, N);

            trac = mu*(normal[0]*computeValue(0, dN_dx) + normal[1]*computeValue(0, dN_dy));

            for(ii=0;ii<totnlbf2;ii++)
            {
              fact = N[ii] * dvol * PENALTY;
              //fact = N[ii] * dvol * (PENALTY*mu/hx);

              TI = ndof*ii+dir;

              Flocal(TI) += fact*res;

              for(jj=0;jj<totnlbf2;jj++)
              {
                Klocal(TI, ndof*jj+dir) += fact * N[jj];
              }
            }

            if(isNitsche)
            {
              for(ii=0;ii<totnlbf2;ii++)
              {
                Ta = mu*dvol*(normal[0]*dN_dx(ii)+normal[1]*dN_dy(ii));
                fact = N(ii)*dvol;

                for(jj=0;jj<totnlbf2;jj++)
                {
                  Tb = mu*(normal[0]*dN_dx(jj)+normal[1]*dN_dy(jj));

                  Klocal(ii, jj)   -= fact*Tb;

                  Klocal(ii, jj)   -= Ta*N(jj)*NitscheFact;
                }

                Flocal(ii)   += fact*trac;
                Flocal(ii)   -= Ta*res*NitscheFact;
              }
            }
        }// for(gp1=0...
        }// for(gp2=0...
      } // for(aa=0;aa<DirichletData.size();aa++)
  } // if(DirichletData.size() > 0)

  return;
}
*/



/*
template<>
void TreeNode<2>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur=0)
{
    // GFEM for Biharmonic equation
    ///////////////////////////////////////
    //PoissonEx1 analy;
    //PoissonEx3 analy;
    BiharmonicEx1 analy;

    int      ii, jj, gp;
    double   uu, vv, Jac, dvol, fact, res, xx, yy, r, beta, betax, betay, b;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), d2NN_dx2(totnlbf), dNN_dy(totnlbf), d2NN_dy2(totnlbf), d2N(totnlbf);
    VectorXd  N, dN_dx, d2N_dx2, dN_dy, d2N_dy2;

    for(gp=0; gp<GeomData->gausspoints.size(); gp++)
    {
        uu   = 0.5*(knots[0][2] * GeomData->gausspoints[gp][0] + knots[0][3]);
        vv   = 0.5*(knots[1][2] * GeomData->gausspoints[gp][1] + knots[1][3]);

        dvol = GeomData->gausspoints[gp] * JacMult;

          GeomData->computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &NN(0), &dNN_dx(0), &dNN_dy(0), &d2NN_dx2(0), &d2NN_dy2(0));

          if(parent == NULL)
          {
            N = NN;
            dN_dx = dNN_dx;
            dN_dy = dNN_dy;
            d2N_dx2 = d2NN_dx2;
            d2N_dy2 = d2NN_dy2;
          }
          else
          {
            N = SubDivMat*NN;
            dN_dx = SubDivMat*dNN_dx;
            dN_dy = SubDivMat*dNN_dy;
            d2N_dx2 = SubDivMat*d2NN_dx2;
            d2N_dy2 = SubDivMat*d2NN_dy2;
          }

          //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", uu, vv, fact, dvol);
          //printf("BasisFuns \n");
          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %12.6f  \t %12.6f  \t %12.6f \n ", N(ii), dN_dx(ii), dN_dy(ii));

          xx = GeomData->ComputeCoord(0, uu);
          yy = GeomData->ComputeCoord(1, vv);
          //r  = sqrt(xx*xx+yy*yy);
          fact = analy.computeForce(0, xx, yy);
          //fact = 0.0;
          
          d2N = d2N_dx2 + d2N_dy2;

          //cout << " AAAAAAAAAAA " << endl;
          Klocal += ( (dvol*d2N) * d2N.transpose() );

          //cout << " AAAAAAAAAAA " << endl;
          Flocal += ( dvol*(N*fact - d2N*computeValue(0, d2N)) );

    }//gp

    //printMatrix(Klocal);
    //printf("\n\n");
    //printVector(Flocal);
    
   return;
}
*/


/*
template<>
void TreeNode<2>::applyDirichletBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  // Nitsche method for Biharmonic equation

  if( DirichletData.size() > 0 )
  {
    int ii, jj, aa, gp1, gp2, TI, TIp1, TIp2, index, TJ, TJp1, TJp2, dir, side, nU, nP;
    double theta, y0, y1, Ta, Tb, res, JacMultLoc, af, temp, NitscheFact;
    double  uu, vv, dvol, specVal, PENALTY, Jac, fact, rad, R, bb1, bb2;
    bool  isNitsche;

    myPoint  param, geom, normal, trac;
    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf);
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf);
    VectorXd  d3N_dx3(totnlbf), d3N_dy3(totnlbf), d3N_dxdy2(totnlbf), d3N_dx2dy(totnlbf);
    VectorXd  dDN_dx(totnlbf), dDN_dy(totnlbf);
    vector<double>  boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2;


    bool   axsy = ((int)elmDat[2] == 1);
    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    af = SolnData->td(2);

    BiharmonicEx1 analy;

      for(aa=0;aa<DirichletData.size();aa++)
      {
       // printVector(DirichletData[aa]);

        isNitsche = false;
        side        = (int) (DirichletData[aa][0] - 1);
        dir         = (int) (DirichletData[aa][1] - 1);
        specVal     = DirichletData[aa][2];
        PENALTY     = DirichletData[aa][3];
        isNitsche   = ( (int) DirichletData[aa][4] == 1 );
        NitscheFact = DirichletData[aa][5];

        //for symmetric Nitsche method -> NitscheFact = 1.0
        //for unsymmetric Nitsche method -> NitscheFact = -1.0

        //cout << " PENALTY " << PENALTY << endl;

        GeomData->getBoundaryNormal2D(side, normal);
        GeomData->setBoundaryGPs2D(side, boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2);

        JacMultLoc = TreeNode<2>::getJacBoundary(side);
        
        for(gp2=0;gp2<boundaryGPs2.size();gp2++)
        {
            param[1] = 0.5 * (knots[1][2] * boundaryGPs2[gp2] + knots[1][3]);
        for(gp1=0;gp1<boundaryGPs1.size();gp1++)
        {
            param[0] = 0.5 * (knots[0][2] * boundaryGPs1[gp1] + knots[0][3]);
 
            dvol = JacMultLoc * boundaryGWs2[gp2] * boundaryGWs1[gp1] ;
            
            //printf(" %4d \t %4d \t %12.6f \t %12.6f \n", side, dir, vv, uu);

            GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy );

            //GeomData->computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv,
              //                &d3N_dx3(0), &d3N_dy3(0), &d3N_dxdy2(0), &d3N_dx2dy(0));

            //dDN_dx = d3N_dx3 + d3N_dxdy2;
            //dDN_dy = d3N_dy3 + d3N_dx2dy;
            
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

            //R  = GeomData->ComputeCoord(0, 1.0);
            R = 0.5;
            geom[0] = GeomData->ComputeCoord(0, param[0]);
            geom[1] = GeomData->ComputeCoord(1, param[1]);
            //r = sqrt(xx*xx+yy*yy);
            //val = 1.0+log(2.0*r);
            //cout << xx << '\t' << yy << '\t' << DirichletData[aa][2] << endl;

            //cout << aa << '\t' << side << '\t' << dir << '\t' << val << endl;
            
            //res = specVal;

            if(axsy)
            {
                dvol *= (2.0*PI*geom[0]);
            }

            specVal = DirichletData[aa][2];

            //specVal = analy.computeValue(dir, xx, yy);
            //cout << side << '\t' << dir << '\t' << specVal << endl;

            //cout << " lllllllllll " << endl;
            specVal *= timeFunction[0].prop;
            //cout << " lllllllllll " << endl;

            //cout << aa << '\t' << side << '\t' << dir << '\t' << res << endl;
            //cout << specVal << '\t' << computeValue(dir, N) << endl;
            //cout << " PENALTY = " << PENALTY << endl;

              res = specVal - computeValue(dir, N);

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

          if(isNitsche)
          {
              res = specVal - computeValue(dir, N);

              trac[0] = mu*(normal[0]*computeValue(0, dDN_dx) + normal[1]*computeValue(0, dDN_dy));
              for(ii=0;ii<totnlbf2;ii++)
              {
                Ta = mu*dvol*(normal[0]*dDN_dx(ii)+normal[1]*dDN_dy(ii));
                fact = N(ii)*dvol;

                for(jj=0;jj<totnlbf2;jj++)
                {
                  Tb = mu*(normal[0]*dDN_dx(jj)+normal[1]*dDN_dy(jj));

                  Klocal(ii, jj)   += fact*Tb;

                  Klocal(ii, jj)   += Ta*N(jj)*NitscheFact;
                }

                Flocal(ii)   -= fact*trac[0];

                Flocal(ii)   += Ta*res*NitscheFact;
              }
          } // if(isNitsche)

        }// for(gp1=0...
        }// for(gp2=0...
      } // for(aa=0;aa<DirichletData.size();aa++)
  } // if(DirichletData.size() > 0)

  return;
}
*/


/*
template<>
void TreeNode<2>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    // LSFEM for Poisson's (or Laplace) problem
    ///////////////////////////////////////

    int      ii, jj, gp1, gp2;
    double   uu, vv, Jac, dvol, fact, res, xx, yy, r, beta, betax, betay, b;
    double   beta1, beta2, delta1, delta2, delta, position[2];

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), dN_dy(totnlbf), d2N_dy2(totnlbf), D(totnlbf), NN;

    beta1 = elmDat[5] * GeomData->GetGridLength(0) * knots[0][2];
    beta2 = elmDat[5] * GeomData->GetGridLength(1) * knots[1][2];

    //cout << " beta " << beta1 << '\t' << beta2 << endl;

    position[0] = 0.5 * GeomData->GetGridLength(0);
    position[1] = 0.5 * GeomData->GetGridLength(0);

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMult;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &N(0), &dN_dx(0), &dN_dy(0), &d2N_dx2(0), &d2N_dy2(0));

          //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", uu, vv, fact, dvol0);
          //printf("BasisFuns \n");
          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %12.6f  \t %12.6f  \t %12.6f \n ", N(ii), dN_dx(ii), dN_dy(ii));

          xx = GeomData->ComputeCoord(0, uu) - position[0];
          yy = GeomData->ComputeCoord(1, vv) - position[1];

          delta1 = DiracDelta1(xx, beta1);
          delta2 = DiracDelta1(yy, beta2);

          delta = delta1 * delta2;

          //if(delta > 0.0)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", xx, yy, delta1, delta2, delta);

          //D = 1.0*(d2N_dx2 + d2N_dy2);

          res = 5.0;
          //res -= delta*elmDat[4];
          //res -= computeForce(0, N)*elmDat[4];
          res -= computeValue(0, D);

          //Flocal += ((dvol*res)*D);
          //Klocal += ((dvol*D)*D.transpose());

          // GFEM
          Klocal += (dvol*(dN_dx*dN_dx.transpose()+dN_dy*dN_dy.transpose()));
          Flocal += ( dvol*(-5.0*N - dN_dx*computeValue(0, dN_dx)  - dN_dy*computeValue(0, dN_dy) ));
    }//gp1
    }//gp2
    
   return;
}
//


//
template<>
void TreeNode<2>::computeAndReturnJacobian(int index, double* position, double* normal, double* specVal, double arclen, double* data, MatrixXd& D, VectorXd& vec1, VectorXd& vec2)
{
    // LSFEM for Poisson problem
    //
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ;
    double  uu, vv, Jac, dvol, fact, fact1, fact2, fact3;
    double  pres, nu, b1, b2, b3, xx, yy, rr;
    double  delta1, delta2, delta, beta1, beta2, res, Ddelta;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), dN_dy(totnlbf), d2N_dy2(totnlbf), vectmp;
    MatrixXd  D1(nsize, 1), D2(nsize, 1);

    beta1 = elmDat[5] * GeomData->GetGridLength(0) * knots[0][2];
    beta2 = elmDat[5] * GeomData->GetGridLength(1) * knots[1][2];

    fact1 = 0.0;
    D2.setZero();
    Ddelta = 0.0;

    count = 0;
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMult;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          //computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &N(0), &dN_dx(0), &dN_dy(0), &d2N_dx2(0), &d2N_dy2(0));

          N = GeomData->shpfns[level][count].N;
          dN_dx = GeomData->shpfns[level][count].dN_dx;
          dN_dy = GeomData->shpfns[level][count].dN_dy;
          d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          count++;

          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);

          xx = GeomData->ComputeCoord(0, uu) - position[0];
          yy = GeomData->ComputeCoord(1, vv) - position[1];

          delta1 = DiracDelta1(xx, beta1);
          delta2 = DiracDelta1(yy, beta2);

          delta = delta1 * delta2;

          //delta = computeForce(0, N);

          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", uu, vv, delta1, delta2, delta);

          vectmp = d2N_dx2 + d2N_dy2;

          fact = delta*dvol;

          Ddelta += (delta*fact);

          D2 += (fact*vectmp);
          fact1 += (5.0 - computeValue(0, vectmp) - delta*data[0])*fact;
          //fact1 += (5.0 - computeValue(0, vectmp))*fact; // 5.0 is the body force. need to generalise this

          fact = delta*delta*dvol*1.0/knots[0][2]/knots[1][2];
          //fact = delta*delta*dvol*1.0e2;
          //Klocal += ( (fact*N)*N.transpose());
          //Flocal += ( N*((0.2-computeValue(0,N))*fact));
       }
    }

    //printf("\n\n");
    //printVector(Flocal);
    //printMatrix(D2);
    //printf("\n\n");

    if(parent == NULL)
    {
      //printf("\nAAAAAAAAAAAAAAAAAAAAA\n");
      D = D2;
    }
    else
    {
      //printf("\nBBBBBBBBBBBBBBBBBBBB  %5d \t %5d\n", nsize,nsize2);
      D.resize(nsize2,1);
      D = SubDivMat*D2;
    }

    vec1(0) = Ddelta;
    vec2(0) = fact1;

   return;
}
*/


          /*
          N = &(GeomData->shpfns[0][count].N(0));
          dN_dx = &(GeomData->shpfns[0][count].dN_dx(0));
          dN_dy = &(GeomData->shpfns[0][count].dN_dy(0));
          d2N_dx2 = &(GeomData->shpfns[0][count].d2N_dx2(0));
          d2N_dy2 = &(GeomData->shpfns[0][count].d2N_dy2(0));
          count++;
          //
          N = GeomData->shpfns[level][count].N;
          dN_dx = GeomData->shpfns[level][count].dN_dx;
          dN_dy = GeomData->shpfns[level][count].dN_dy;
          d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          count++;

          if(mpapTime.cur > 40.0 && mpapTime.cur < 60)
          {
            if(xx > 9.7 && xx < 10.3)
            { 
              if(yy < 9.4 && yy > 9.1)
              { 
                //cout << " AAAAAAAAA  " << xx << '\t' << yy << endl;
                res(0) = 0.005;
              }
              if(yy < 10.9 && yy > 10.6)
              { 
                //cout << " BBBBBBBBB  " << xx << '\t' << yy << endl;
                res(0) = -0.005;
              }
            }
          }
          */

/*
template<>
void TreeNode<2>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    // LSFEM for Navier-Stokes (steady-state)
    //
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ;
    double  uu, vv, Jac, dvol, fact, fact1, fact2, fact3;
    double  pres, nu, b1, b2, b3, xx, yy, trgradu, gamma=1.0, HH;
    double  delta, delta1, delta2, position[2], beta1, beta2;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), dN_dy(totnlbf), d2N_dy2(totnlbf);
    VectorXd  res(3), dp(2), Du(2), vel(2), vectmp;
    MatrixXd  D(nsize, 3), F(2,2);

    beta1 = elmDat[5] * GeomData->GetGridLength(0) * knots[0][2];
    beta2 = elmDat[5] * GeomData->GetGridLength(1) * knots[1][2];

    position[0] = 0.5 * GeomData->GetGridLength(0);
    position[1] = 0.5 * GeomData->GetGridLength(1);

    gamma = 1.0/knots[0][2]/knots[1][2];
    //gamma = 1.0e6;
    gamma = 1.0;

    count = 0;
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMult;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &N(0), &dN_dx(0), &dN_dy(0), &d2N_dx2(0), &d2N_dy2(0));
          //computeBasisFunctions2D(uu, vv, N, dN_dx, dN_dy, d2N_dx2, d2N_dy2);

          //N = GeomData->shpfns[level][count].N;
          //dN_dx = GeomData->shpfns[level][count].dN_dx;
          //dN_dy = GeomData->shpfns[level][count].dN_dy;
          //d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          //d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          count++;

          xx = GeomData->ComputeCoord(0, uu) - position[0];
          yy = GeomData->ComputeCoord(1, vv) - position[1];

          delta1 = DiracDelta1(xx, beta1);
          delta2 = DiracDelta1(yy, beta2);

          delta = delta1 * delta2;

          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);

          vel(0) = computeValue(0, N);
          vel(1) = computeValue(1, N);

          vectmp = d2N_dx2 + d2N_dy2;
          Du(0) = computeValue(0, vectmp);
          Du(1) = computeValue(1, vectmp);

          dp(0) = computeValue(2, dN_dx);
          dp(1) = computeValue(2, dN_dy);

          F(0,0) = computeValue(0, dN_dx);
          F(0,1) = computeValue(0, dN_dy);
          F(1,0) = computeValue(1, dN_dx);
          F(1,1) = computeValue(1, dN_dy);

          for(ii=0;ii<totnlbf;ii++)
          {
             TI   =  3*ii;
             TIp1 =  TI+1;
             TIp2 =  TI+2;

             b1 = dN_dx[ii];
             b2 = dN_dy[ii];
             b3 = N[ii];

             fact = rho*(vel(0)*b1 + vel(1)*b2) - mu*vectmp[ii];
             
             D(TI,0)   = fact + F(0,0)*b3;
             D(TIp1,0) = F(0,1)*b3;
             D(TIp2,0) = b1;

             D(TI,1)   = F(1,0)*b3;
             D(TIp1,1) = fact + F(1,1)*b3;
             D(TIp2,1) = b2;
           
             D(TI,2)   = gamma*b1;
             D(TIp1,2) = gamma*b2;
             D(TIp2,2) = 0.0;
          }

          res.setZero();
          //res(0) = computeForce(0, N) ;
          //res(1) = computeForce(1, N) ;

          //res(0) += delta*0.3;
          //res(1) -= delta*0.3;

          res(0) += (F(0,0)*vel(0) + F(0,1)*vel(1) - mu*Du(0) + dp(0));
          res(1) += (F(1,0)*vel(0) + F(1,1)*vel(1) - mu*Du(1) + dp(1));
          res(2) += (gamma*F.trace());

          //printf("\n\n \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\t%14.8f\t%14.8f\t%14.8f\t%14.8f\n\n\n", res(0), res(1), F(0,0), F(0,1), F(1,0), F(1,1), JacMult, Jac, dvol);
          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\t%14.8f\t%14.8f\n", res(0), res(1), res(2), F(0,0), F(0,1), F(1,0), F(1,1));

          if(ind1) Klocal += ((dvol*D)*D.transpose());
          Flocal -= (D*(dvol*res));
       }
    }

   return;
}
*/



/*
template<>
void TreeNode<2>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    // LSFEM for Navier-Stokes (steady-state) 4 DOF
    //
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, TIp3, count, TJ;
    double  uu, vv, Jac, dvol, fact, fact1, fact2, fact3;
    double  pres, nu, b1, b2, b3, xx, yy, trgradu, gamma=1.0, HH;
    double  delta, delta1, delta2, position[2], beta1, beta2, vort;

    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf);
    VectorXd  res(4), dp(2), vel(2), dw(2);
    MatrixXd  D(nsize, 4), F(2,2);

    beta1 = elmDat[5] * GeomData->GetGridLength(0) * knots[0][2];
    beta2 = elmDat[5] * GeomData->GetGridLength(1) * knots[1][2];

    position[0] = 0.5 * GeomData->GetGridLength(0);
    position[1] = 0.5 * GeomData->GetGridLength(1);

    gamma = 1.0/knots[0][2]/knots[1][2];
    //gamma = 1.0e6;
    gamma = 1.0;

    count = 0;
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMult;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          //computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &N(0), &dN_dx(0), &dN_dy(0), &d2N_dx2(0), &d2N_dy2(0));
          //computeBasisFunctions2D(uu, vv, N, dN_dx, dN_dy, d2N_dx2, d2N_dy2);

          N = GeomData->shpfns[level][count].N;
          dN_dx = GeomData->shpfns[level][count].dN_dx;
          dN_dy = GeomData->shpfns[level][count].dN_dy;
          //d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          //d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          count++;

          xx = GeomData->ComputeCoord(0, uu) - position[0];
          yy = GeomData->ComputeCoord(1, vv) - position[1];

          delta1 = DiracDelta1(xx, beta1);
          delta2 = DiracDelta1(yy, beta2);

          delta = delta1 * delta2;

          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);

          vel(0) = computeValue(0, N);
          vel(1) = computeValue(1, N);

          dp(0) = computeValue(3, dN_dx);
          dp(1) = computeValue(3, dN_dy);

          dw(0) = computeValue(2, dN_dx);
          dw(1) = computeValue(2, dN_dy);
          vort  = computeValue(2, N);

          F(0,0) = computeValue(0, dN_dx);
          F(0,1) = computeValue(0, dN_dy);
          F(1,0) = computeValue(1, dN_dx);
          F(1,1) = computeValue(1, dN_dy);

          for(ii=0;ii<totnlbf;ii++)
          {
             TI   =  ndof*ii;
             TIp1 =  TI+1;
             TIp2 =  TI+2;
             TIp3 =  TI+3;

             b1 = dN_dx[ii];
             b2 = dN_dy[ii];
             b3 = N[ii];

             fact = rho*(vel(0)*b1 + vel(1)*b2);
             
             D(TI,0)   = fact + F(0,0)*b3;
             D(TIp1,0) = F(0,1)*b3;
             D(TIp2,0) = mu*b2;
             D(TIp3,0) = b1;

             D(TI,1)   = F(1,0)*b3;
             D(TIp1,1) = fact + F(1,1)*b3;
             D(TIp2,1) = -mu*b1;
             D(TIp3,1) = b2;
           
             D(TI,2)   = b2;
             D(TIp1,2) = -b1;
             D(TIp2,2) = b3;
             D(TIp3,2) = 0.0;

             D(TI,3)   = b1;
             D(TIp1,3) = b2;
             D(TIp2,3) = 0.0;
             D(TIp3,3) = 0.0;
          }

          res.setZero();
          //res(0) = computeForce(0, N) ;
          //res(1) = computeForce(1, N) ;

          //res(0) += delta*0.3;
          //res(1) -= delta*0.3;

          res(0) += (F(0,0)*vel(0) + F(0,1)*vel(1) + mu*dw(1) + dp(0));
          res(1) += (F(1,0)*vel(0) + F(1,1)*vel(1) - mu*dw(0) + dp(1));
          res(2) += (vort+F(0,1)-F(1,0));
          res(3) += (F.trace());

          //printf("\n\n \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\t%14.8f\t%14.8f\t%14.8f\t%14.8f\n\n\n", res(0), res(1), F(0,0), F(0,1), F(1,0), F(1,1), JacMult, Jac, dvol);
          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\t%14.8f\t%14.8f\n", res(0), res(1), res(2), F(0,0), F(0,1), F(1,0), F(1,1));

          if(ind1) Klocal += ((dvol*D)*D.transpose());
          Flocal -= (D*(dvol*res));

       }
    }

   return;
}
*/




/*
template<>
void TreeNode<2>::computeAndReturnJacobian(int index, double* position, double* normal, double* specVal, double arclen, double* data, MatrixXd& D, VectorXd& vec1, VectorXd& vec2)
{
    // LSFEM for Navier-Stokes (steady-state)
    //
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ;
    double  uu, vv, Jac, dvol, fact, fact1, fact2, fact3;
    double  pres, nu, b1, b2, b3, xx, yy;
    double  delta1, delta2, delta, Ddelta, beta1, beta2;

    Circle  circle2(10.0,10.0,0.5);
    
    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), dN_dy(totnlbf), d2N_dy2(totnlbf);
    VectorXd  res(2), dp(2), Du(2), vel(3), vectmp, vel2(3);
    MatrixXd  D1(nsize, 2), D2(nsize, 2), F(2,2), NN(nsize,3);

    beta1 = elmDat[5] * GeomData->GetGridLength(0) * knots[0][2];
    beta2 = elmDat[5] * GeomData->GetGridLength(1) * knots[1][2];

    Ddelta = 0.0;
    fact1 = 0.0;
    fact2 = 0.0;
    D2.setZero();
    NN.setZero();
    vel.setZero();

    count = 0;
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMult;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          //computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &N(0), &dN_dx(0), &dN_dy(0), &d2N_dx2(0), &d2N_dy2(0));
          //computeBasisFunctions2D(uu, vv, N, dN_dx, dN_dy, d2N_dx2, d2N_dy2);

          N = GeomData->shpfns[level][count].N;
          dN_dx = GeomData->shpfns[level][count].dN_dx;
          dN_dy = GeomData->shpfns[level][count].dN_dy;
          d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          count++;

          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);

          xx = GeomData->ComputeCoord(0, uu) - position[0];
          yy = GeomData->ComputeCoord(1, vv) - position[1];

          delta1 = DiracDelta1(xx, beta1);
          delta2 = DiracDelta1(yy, beta2);

          delta = delta1 * delta2;

          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", xx, yy, delta1, delta2, delta);

          vel(0) = computeValue(0, N);
          vel(1) = computeValue(1, N);

          vectmp = d2N_dx2 + d2N_dy2;
          Du(0) = computeValue(0, vectmp);
          Du(1) = computeValue(1, vectmp);

          dp(0) = computeValue(2, dN_dx);
          dp(1) = computeValue(2, dN_dy);

          F(0,0) = computeValue(0, dN_dx);
          F(0,1) = computeValue(0, dN_dy);
          F(1,0) = computeValue(1, dN_dx);
          F(1,1) = computeValue(1, dN_dy);

          for(ii=0;ii<totnlbf;ii++)
          {
             TI   =  3*ii;
             TIp1 =  TI+1;
             TIp2 =  TI+2;

             b1 = dN_dx[ii];
             b2 = dN_dy[ii];
             b3 = N[ii];

             fact = rho*(vel(0)*b1 + vel(1)*b2) - mu*vectmp[ii];

             D1(TI,0)   = fact + F(0,0)*b3;
             D1(TIp1,0) = F(0,1)*b3;
             D1(TIp2,0) = b1;

             D1(TI,1)   = F(1,0)*b3;
             D1(TIp1,1) = fact + F(1,1)*b3;
             D1(TIp2,1) = b2;

             NN(TI,0)   = b3;
             NN(TIp1,1) = b3;
          }

          res.setZero();
          //res(0) = computeForce(0, N) ;
          //res(1) = computeForce(1, N) ;

          delta *= arclen;

          res(0) += (F(0,0)*vel(0) + F(0,1)*vel(1) - mu*Du(0) + dp(0) + delta*data[0]) ;
          res(1) += (F(1,0)*vel(0) + F(1,1)*vel(1) - mu*Du(1) + dp(1) + delta*data[1]);

          fact = delta*dvol;

          D2 += (fact*D1);
          Ddelta += delta*fact;

          fact1 += (0.0-res(0))*fact; // 0.0 is the body force. need to generalise this
          fact2 += (0.0-res(1))*fact;

          vel2(0) = 0.0 - vel(0);
          vel2(1) = 0.0 - vel(1);

          delta = delta1 * delta2;

          fact = delta*delta*dvol*1.0;
          //Klocal += ( (fact*NN)*NN.transpose());
          //Flocal += ( NN*(vel2*fact));
       }
    }

    //printf("\n\n");
    //printVector(Flocal);
    //printMatrix(D2);
    //printf("\n\n");

    if(parent == NULL)
    {
      D = D2;
    }
    else
    {
      D.resize(nsize2,1);
      D = SubDivMat*D2;
    }

    vec1(0) = Ddelta;
    vec1(1) = Ddelta;
    vec2(0) = fact1;
    vec2(1) = fact2;

   return;
}
*/




/*
template<>
void TreeNode<2>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    // GFEM for Navier-Stokes flow (steady-state)
    //
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;
   
    double  uu, vv, Jac, dvol, fact, nu, b1, b2, b3, b4, b5, fact1, fact2, fact3, ci, trgradu, dt, gamma;
    double  beta1, beta2, delta1, delta2, delta, position[2], xx, yy;
    double  dist, HH, pres, h1, h2, Da, Db, rad, urdr, urdr2, eps;

    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), dN_dy(totnlbf), d2N_dy2(totnlbf);
    VectorXd  res(2), dp(2), Du(2), vel(3), vel2(3), vectmp(totnlbf), force(3), stres(3);
    MatrixXd  D(nsize,2), F(2,2), FN(2,2);
    D.setZero();

    count = 0;
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMult;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &N(0), &dN_dx(0), &dN_dy(0), &d2N_dx2(0), &d2N_dy2(0));

          //N = GeomData->shpfns[level][count].N;
          //dN_dx = GeomData->shpfns[level][count].dN_dx;
          //dN_dy = GeomData->shpfns[level][count].dN_dy;
          //d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          //d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          count++;

          //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", res(0), res(1), JacMult, Jac, fact, dvol);

          //xx = GeomData->ComputeCoord(0, uu) - position[0];
          //yy = GeomData->ComputeCoord(1, vv) - position[1];

          //delta1 = DiracDelta1(xx, beta1);
          //delta2 = DiracDelta1(yy, beta2);
          //delta = delta1 * delta2;

          vel(0) = computeValue(0, N);
          vel(1) = computeValue(1, N);

          vectmp = d2N_dx2 + d2N_dy2;
          Du(0) = computeValue(0, vectmp);
          Du(1) = computeValue(1, vectmp);

          pres  = computeValue(2, N);
          dp(0) = computeValue(2, dN_dx);
          dp(1) = computeValue(2, dN_dy);

          F(0,0) = computeValue(0, dN_dx);
          F(0,1) = computeValue(0, dN_dy);
          F(1,0) = computeValue(1, dN_dx);
          F(1,1) = computeValue(1, dN_dy);

          trgradu = F.trace();

          force.setZero();
          //force(0) = analy.computeXForce(uu, vv);
          //force(1) = analy.computeYForce(uu, vv);
          //force(0) = -0.3*delta;

          force(0) -= rho*(F(0,0)*vel(0)+F(0,1)*vel(1));
          force(1) -= rho*(F(1,0)*vel(0)+F(1,1)*vel(1));

          xx = GeomData->ComputeCoord(0, uu);
          yy = GeomData->ComputeCoord(1, vv);
          if(axsy)
          {
            rad = xx;

            urdr  = vel(0)/rad;
            urdr2 = urdr/rad;
            dvol *= (2*PI*rad);
          }

          for(ii=0;ii<totnlbf;ii++)
          {
             TI = ndof*ii;
             TIp1 = TI+1;
             TIp2 = TI+2;

             // GLS stabilisation term

             b1 = dN_dx[ii];
             b2 = dN_dy[ii];
             b3 = N[ii];

             fact = rho*(vel(0)*b1 + vel(1)*b2) - mu*vectmp[ii];

             //D(TI,0)   = F(0,0)*b3 + fact;
             //D(TIp1,0) = F(0,1)*b3;
             D(TIp2,0) = b1;

             //D(TI,1)   = F(1,0)*b3;
             //D(TIp1,1) = F(1,1)*b3 + fact;
             D(TIp2,1) = b2;

             ////////////////////////////////////

             b1 = dN_dx(ii)*dvol;
             b2 = dN_dy(ii)*dvol;
             b3 = N(ii)*dvol;

             b4 = mu*b1;
             b5 = mu*b2;

             FN = (rho*b3)*F;

             for(jj=0;jj<totnlbf;jj++)
             {
               TJ = ndof*jj;
               TJp1 = TJ+1;
               TJp2 = TJ+2;

               // Diffusion term

               fact = b4*dN_dx(jj) + b5*dN_dy(jj);

               Klocal(TI,   TJ)   += fact;
               Klocal(TIp1, TJp1) += fact;

               // Convection term
               Db = rho*(vel(0)*dN_dx(jj) + vel(1)*dN_dy(jj));

               Klocal(TI,   TJ)   += (FN(0,0)*N(jj) + b3*Db);
               Klocal(TI,   TJp1) += (FN(0,1)*N(jj));
               Klocal(TIp1, TJ)   += (FN(1,0)*N(jj));
               Klocal(TIp1, TJp1) += (FN(1,1)*N(jj) + b3*Db);

               // pressure term
               Klocal(TI,   TJp2) -= (b1*N(jj));
               Klocal(TIp1, TJp2) -= (b2*N(jj));

               // continuity term
               Klocal(TIp2, TJ)   -= (b3*dN_dx(jj));
               Klocal(TIp2, TJp1) -= (b3*dN_dy(jj));
               Klocal(TIp2, TJp2) += 0.0;

               if(axsy)
               {
                  // diffusion term
                  Klocal(TI, TJ)     += (mu*b3*N(jj)/rad/rad);
                  Klocal(TI, TJp2)   -= (b3*N(jj)/rad);

                  // continuity equation
                  Klocal(TIp2, TJ)   -= (b3*N(jj)/rad);
               }
             }

             Flocal(TI)   += (b3*force(0) - b4*F(0,0) - b5*F(0,1) + b1*pres);
             Flocal(TIp1) += (b3*force(1) - b4*F(1,0) - b5*F(1,1) + b2*pres);
             Flocal(TIp2) += (b3*trgradu);

             if(axsy)
             {
                Flocal(TI)   += (-b3*(mu*vel(0)/rad/rad));
                Flocal(TI)   += (b3*pres/rad);
                Flocal(TIp2) += (b3*vel(0)/rad);
             }
          }

          //res = force;
          res.setZero();

          //res(0) -= (F(0,0)*vel(0)+F(0,1)*vel(1)-mu*Du(0) + dp(0));
          //res(1) -= (F(1,0)*vel(0)+F(1,1)*vel(1)-mu*Du(1) + dp(1));

          //res(0) -= (-mu*Du(0) + dp(0));
          //res(1) -= (-mu*Du(1) + dp(1));

          res(0) -=  dp(0);
          res(1) -=  dp(1);

          dvol *= stabParam;

          Klocal -= (dvol*(D*D.transpose() ));
          Flocal -= (D*(dvol*res));
    }//gp1
    }//gp2
    
    //if(id == 555)
      //printMatrix(Klocal);
    
    //printVector(Flocal);
    //printf("\n\n");
    
    return;
}
*/


//
template<>
void TreeNode<2>::computeAndReturnJacobian(int index, double* position, double* normal, double* specVal, double arclen, double* data, MatrixXd& D, VectorXd& vec1, VectorXd& vec2)
{
    // LSFEM for Stokes (steady-state)
    //
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, TI, TIp1, TIp2, count, TJ;
    double  uu, vv, Jac, dvol, fact, fact1, fact2, fact3;
    double  pres, nu, b1, b2, b3, xx, yy;
    double  delta1, delta2, delta, Ddelta, beta1, beta2, h1, h2, stabParam, mu;
    
    VectorXd  N(totnlbf), dN_dx(totnlbf), d2N_dx2(totnlbf), dN_dy(totnlbf), d2N_dy2(totnlbf);
    VectorXd  res(2), dp(2), Du(2), vel(3), vectmp, vel2(3);
    MatrixXd  D1(nsize, 2), D2(nsize, 2), F(2,2), NN(nsize,2);
    myPoint  param;

    beta1 = elmDat[5] * GeomData->GetGridLength(0) * knots[0][2];
    beta2 = elmDat[5] * GeomData->GetGridLength(1) * knots[1][2];
    
    mu = elmDat[4];

    fact1 = 0.0;
    fact2 = 0.0;
    Ddelta = 0.0;
    D2.setZero();
    NN.setZero();
    vel.setZero();

    count = 0;
    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       Jac = GeomData->gaussweights2[gp2] * JacMultElem;
       
       for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
       {
          param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
          dvol = GeomData->gaussweights1[gp1] * Jac;

          //GeomData->computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], uu, vv, &N(0), &dN_dx(0), &dN_dy(0), &d2N_dx2(0), &d2N_dy2(0));

          //N = GeomData->shpfns[level][count].N;
          //dN_dx = GeomData->shpfns[level][count].dN_dx;
          //dN_dy = GeomData->shpfns[level][count].dN_dy;
          //d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          //d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          count++;

          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", N[ii], dN_dx[ii], dN_dy[ii], d2N_dx2[ii], d2N_dy2[ii]);

          xx = GeomData->ComputeCoord(0, uu) - position[0];
          yy = GeomData->ComputeCoord(1, vv) - position[1];

          delta1 = DiracDelta1(xx, beta1);
          delta2 = DiracDelta1(yy, beta2);

          delta = delta1 * delta2;

          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", xx, yy, delta1, delta2, delta);

          vectmp = d2N_dx2 + d2N_dy2;

          Du(0) = computeValue(0, vectmp);
          Du(1) = computeValue(1, vectmp);

          dp(0) = computeValue(2, dN_dx);
          dp(1) = computeValue(2, dN_dy);

          for(ii=0;ii<totnlbf;ii++)
          {
             TI   =  3*ii;
             TIp1 =  TI+1;
             TIp2 =  TI+2;

             fact = -mu * vectmp[ii];

             D1(TI,0)   = fact;
             D1(TIp1,0) = 0.0;
             D1(TIp2,0) = dN_dx[ii];

             D1(TI,1)   = 0.0;
             D1(TIp1,1) = fact;
             D1(TIp2,1) = dN_dy[ii];

             NN(TI,0)   = N[ii];
             NN(TIp1,1) = N[ii];
          }

          //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", res(0), res(1), delta1, delta2, delta);

          //delta *= arclen;
          
          delta *= stabParam;

          fact = delta*dvol;

          D2 += (fact*D1+NN*dvol*delta*1.0e3);

          //Ddelta += (delta*fact);

          res.setZero();
          //res(0) = computeForce(0, N) ;
          //res(1) = computeForce(1, N) ;

          //res(0) = -mu*Du(0) + dp(0) + delta*data[0];
          //res(1) = -mu*Du(1) + dp(1) + delta*data[1];

          res(0) = -mu*Du(0) + dp(0) ;
          res(1) = -mu*Du(1) + dp(1) ;
          
          fact1  += (0.0-res(0))*fact; // 0.0 is the body force. need to generalise this
          fact2  += (0.0-res(1))*fact;

          //vel(0) = computeValue(0, N);
          //vel(1) = computeValue(1, N);

          //vel2.setZero();

          //vel2(0) = 0.0 - vel(0);
          //vel2(1) = 0.0 - vel(1);

          //delta = delta1 * delta2;

          //fact = delta*delta*dvol*1.0;
          //Klocal += ( (fact*NN)*NN.transpose());
          //Flocal += ( NN*(vel2*fact));
       }
    }

    //printf("\n\n");
    //printVector(Flocal);
    //printMatrix(D2);
    //printf("\n\n");

    /*
    if(parent == NULL)
    {
      D = D2;
    }
    else
    {
      D.resize(nsize2,2);
      D = SubDivMat2*D2;
    }
    */

    vec2(0) = fact1;
    vec2(1) = fact2;

   return;
}
//





template<>
void TreeNode<2>::unRefine()
{
   if(child == NULL)
   {
      cout << " This element has no children. So can't unrefine it ... " << endl;
      return;
   }
   
   //int  ii, jj;

//   for(ii=0;ii<NUM_CHILDREN;ii++)
//      child[ii] = NULL;

//   child = NULL;

//cout << " AAAAAAAAAA " << endl;

   return;
}


