
#include "TreeNode.h"
#include "MpapTime.h"
#include "DistFunctions.h"
#include "Functions.h"
#include "TimeFunction.h"
#include "SolutionData.h"
#include "myDataIntegrateCutFEM.h"

#include "BasisFunctionsBSpline.h"

extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;

MatrixXd  Iden = MatrixXd::Identity(4, 4);

template<>
void TreeNode<3>::unRefine()
{
  return;
}

template<>
int TreeNode<3>::calcError(int ind, int domainCur)
{
  return 0;
}

template<>
double TreeNode<3>::computeTotalBodyForce(int, int)
{
  return 0.0;
}


template<>
void TreeNode<3>::setInitialProfile()
{
  return;
}



template<>
void TreeNode<3>::MatrixToMapResult(int ind1, int ind2, SparseMatrixXd& globalK)
{
    int ii, jj, gp1, gp2, gp3;
   
    double  JacZ, JacY, dvol;

    MatrixXd  Kl(totnlbf2, totnlbf2);
    VectorXd  NN(totnlbf), N(totnlbf2);
    myPoint  param;

    Kl.setZero();

    for(gp3=0;gp3<GeomData->GetNGP(2);gp3++)
    {
       param[2]  = 0.5*(knots[2][2] * GeomData->gausspoints3[gp3] + knots[2][3]);
       JacZ = GeomData->gaussweights3[gp3] * JacMultElem;
       //cout << " ww " << ww << endl;

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       JacY = GeomData->gaussweights2[gp2] * JacZ;
       //cout << " vv " << vv << endl;
       
    for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
    {
       param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
       dvol = GeomData->gaussweights1[gp1] * JacY;

       GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, N);

       if(parent == NULL)
         N = NN;
       else
         N = SubDivMat*NN;

       Kl += ((N*dvol)*N.transpose());
    }
    }
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
void TreeNode<3>::RhsToMapResult(int ind1, int ind2, double* rhs)
{
    int ii, jj, gp1, gp2, gp3;
   
    double  JacY, JacZ, dvol, val;
    myPoint param;

    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), N(totnlbf2), dN_dx(totnlbf2), dN_dy(totnlbf2), Fl(totnlbf2);

    Fl.setZero();

    for(gp3=0;gp3<GeomData->GetNGP(2);gp3++)
    {
       param[2]  = 0.5*(knots[2][2] * GeomData->gausspoints3[gp3] + knots[2][3]);
       JacZ = GeomData->gaussweights3[gp3] * JacMultElem;
       //cout << " ww " << ww << endl;

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       JacY = GeomData->gaussweights2[gp2] * JacZ;
       //cout << " vv " << vv << endl;
       
    for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
    {
       param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
       dvol = GeomData->gaussweights1[gp1] * JacY;

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
    }
    }

    for(ii=0;ii<totnlbf2;ii++)
      rhs[GlobalBasisFuncs[ii]] += Fl(ii);

  return;
}


/*
template<>
void TreeNode<3>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Poisson's (or Laplace) problem
    ///////////////////////////////////////

    PoissonEx2 analy;

    int      ii, jj, gp1, gp2, gp3;
    double   JacZ, JacY, dvol, fact, res, xx, yy, zz, r, beta, betax, betay, b;
    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf), dN_dz(totnlbf);
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), dNN_dz(totnlbf);
    myPoint  param;

    for(gp3=0;gp3<GeomData->GetNGP(2);gp3++)
    {
       param[2]  = 0.5*(knots[2][2] * GeomData->gausspoints3[gp3] + knots[2][3]);
       JacZ = GeomData->gaussweights3[gp3] * JacMultElem;
       //cout << " ww " << ww << endl;

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       JacY = GeomData->gaussweights2[gp2] * JacZ;
       //cout << " vv " << vv << endl;
       
    for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
    {
       param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
       dvol = GeomData->gaussweights1[gp1] * JacY;

       GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

       //printf(" %5d \t %5d \t %5d n", gp3, gp2, gp1);
       //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", uu, vv, fact, dvol0);
       //printf("BasisFuns \n");
       //for(ii=0;ii<totnlbf;ii++)
       //printf(" \t %12.6f  \t %12.6f  \t %12.6f \n ", N(ii), dN_dx(ii), dN_dy(ii));

          if(parent == NULL)
          {
            N     = NN;
            dN_dx = dNN_dx;
            dN_dy = dNN_dy;
            dN_dz = dNN_dz;
          }
          else
          {
            N     = SubDivMat*NN;
            dN_dx = SubDivMat*dNN_dx;
            dN_dy = SubDivMat*dNN_dy;
            dN_dz = SubDivMat*dNN_dz;
          }


       xx = GeomData->ComputeCoord(0, param[0]);
       yy = GeomData->ComputeCoord(1, param[1]);
       zz = GeomData->ComputeCoord(2, param[2]);
       //r  = sqrt(xx*xx+yy*yy);
       fact = analy.computeForce(0, xx, yy, zz);

       Klocal += (dvol*(dN_dx*dN_dx.transpose()+dN_dy*dN_dy.transpose()+dN_dz*dN_dz.transpose()));
       Flocal += (dvol*(N*fact - dN_dx*computeValue(0,dN_dx) - dN_dy*computeValue(0,dN_dy) - dN_dz*computeValue(0,dN_dz)));
       //cout << " AAAAAAAAA " << endl;
    }//gp1
    }//gp2
    }//gp3
  
   return;
}
*/


/*
template<>
void TreeNode<3>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    // LSFEM for Poisson's (or Laplace) problem
    ///////////////////////////////////////

    Laplace_Ex1 analy;

    int      ii, jj, gp1, gp2, gp3;
    double   uu, vv, ww, JacZ, JacY, dvol, fact, res, xx, yy, zz, r, beta, betax, betay, b;
    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf), dN_dz(totnlbf);
    VectorXd  d2N_dx2(totnlbf), d2N_dy2(totnlbf), d2N_dz2(totnlbf), D(totnlbf);
    
    for(gp3=0;gp3<GeomData->GetNGP(2);gp3++)
    {
       ww  = 0.5*(knots[2][2] * GeomData->gausspoints3[gp3] + knots[2][3]);
       JacZ = GeomData->gaussweights3[gp3] * JacMult;
       //cout << " ww " << ww << endl;

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       JacY = GeomData->gaussweights2[gp2] * JacZ;
       //cout << " vv " << vv << endl;
       
    for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
    {
       uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
       dvol = GeomData->gaussweights1[gp1] * JacY;
       //cout << " uu " << uu << endl;

       computeBasisFunctions3D(knots[0][0], knots[1][0], knots[2][0], knots[0][2], knots[1][2], knots[2][2], 
					 uu, vv, ww, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0),
					 &d2N_dx2(0), &d2N_dy2(0), &d2N_dz2(0));

       //printf(" %5d \t %5d \t %5d n", gp3, gp2, gp1);
       //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", uu, vv, fact, dvol0);
       //printf("BasisFuns \n");
       //for(ii=0;ii<totnlbf;ii++)
       //printf(" \t %12.6f  \t %12.6f  \t %12.6f \n ", N(ii), dN_dx(ii), dN_dy(ii));

       xx = GeomData->ComputeCoord(0, uu);
       yy = GeomData->ComputeCoord(1, vv);
       zz = GeomData->ComputeCoord(2, ww);
       //r  = sqrt(xx*xx+yy*yy);
       res = analy.computeForceAt(xx, yy, zz);

       //Klocal += (dvol*(dN_dx*dN_dx.transpose()+dN_dy*dN_dy.transpose()+dN_dz*dN_dz.transpose()));
       //Flocal += (dvol*(N*res - dN_dx*computeValue(0,dN_dx) - dN_dy*computeValue(0,dN_dy) - dN_dz*computeValue(0,dN_dz)));
       //cout << " AAAAAAAAA " << endl;

       D = d2N_dx2 + d2N_dy2 + d2N_dz2;

       //res += computeForce(0, N);
       res -= computeValue(0, D);
       //res += 8.0*r*r + 4.0;

       Flocal += ((dvol*res)*D);
       Klocal += ((dvol*D)*D.transpose());
    }//gp1
    }//gp2
    }//gp3
  
   return;
}
*/



/*
template<>
void TreeNode<3>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Stokes flow
    //
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, gp3, TI, TIp1, TIp2, TIp3, count, TJ, TJp1, TJp2, TJp3, nU, nP;
 
    nU  = forAssyVec.size();
    nP  = forAssyVec2.size();

    double  uu, vv, ww, JacY, JacZ, Jac, dvol, fact, b1, b2, b3, b4, b5, b6, b7, b8;
    double  xx, yy, zz, pres, afdt, af, am, acceFact, muTaf, d1, c1;

    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf), dN_dz(totnlbf), tmpvec;
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), dNN_dz(totnlbf);
    VectorXd  res(3), dp(3), Du(3), vel(3), force(3), velDot(3);
    MatrixXd  F(3,3), stress(3,3);

    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = rho*am*SolnData->td(9);
    muTaf = mu*af;

    count = 0;
    for(gp3=0;gp3<GeomData->GetNGP(2);gp3++)
    {
       ww  = 0.5*(knots[2][2] * GeomData->gausspoints3[gp3] + knots[2][3]);
       JacZ = GeomData->gaussweights3[gp3] * JacMult;

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       JacY = GeomData->gaussweights2[gp2] * JacZ;
       
    for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
    {
       uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
       dvol = GeomData->gaussweights1[gp1] * JacY;

       GeomData->computeBasisFunctions3D(knots[0][0], knots[1][0], knots[2][0], knots[0][2], knots[1][2], knots[2][2], 
          uu, vv, ww, &NN(0), &dNN_dx(0), &dNN_dy(0), &(dNN_dz(0)));

          //N = GeomData->shpfns[level][count].N;
          //dN_dx = GeomData->shpfns[level][count].dN_dx;
          //dN_dy = GeomData->shpfns[level][count].dN_dy;
          //dN_dz = GeomData->shpfns[level][count].dN_dz;
          //d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          //d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          //count++;

          //xx = GeomData->ComputeCoord(0, uu);
          //yy = GeomData->ComputeCoord(1, vv);
          //zz = GeomData->ComputeCoord(2, ww);

          if(parent == NULL)
          {
            N     = NN;
            dN_dx = dNN_dx;
            dN_dy = dNN_dy;
            dN_dz = dNN_dz;
          }
          else
          {
            N     = SubDivMat*NN;
            dN_dx = SubDivMat*dNN_dx;
            dN_dy = SubDivMat*dNN_dy;
            dN_dz = SubDivMat*dNN_dz;
          }

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);
          vel(2) = computeValueCur(2, N);

          F(0,0) = computeValueCur(0, dN_dx);
          F(0,1) = computeValueCur(0, dN_dy);
          F(0,2) = computeValueCur(0, dN_dz);
          F(1,0) = computeValueCur(1, dN_dx);
          F(1,1) = computeValueCur(1, dN_dy);
          F(1,2) = computeValueCur(1, dN_dz);
          F(2,0) = computeValueCur(2, dN_dx);
          F(2,1) = computeValueCur(2, dN_dy);
          F(2,2) = computeValueCur(2, dN_dz);

          pres   = computeValue2(0, N);
          dp(0)  = computeValue2(0, dN_dx);
          dp(1)  = computeValue2(0, dN_dy);
          dp(2)  = computeValue2(0, dN_dz);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          velDot(2) = computeValueDotCur(2, N);

          velDot *= rho;
          //velDot.setZero();

          //vectmp = d2N_dx2 + d2N_dy2;
          //Du(0) = computeValueCur(0, vectmp);
          //Du(1) = computeValueCur(1, vectmp);
          //Du(2) = computeValueCur(2, vectmp);

          // this is pseudo-stress
          stress = mu*F;
          stress(0,0) -= pres;
          stress(1,1) -= pres;
          stress(2,2) -= pres;

          force.setZero();

          for(ii=0;ii<totnlbf2;ii++)
          {
             TI   = ndof*ii;
             TIp1 = TI+1;
             TIp2 = TI+2;
             TIp3 = nU+ii;

             b1 = dN_dx[ii]*dvol;
             b2 = dN_dy[ii]*dvol;
             b3 = dN_dz[ii]*dvol;
             b4 = N[ii]*dvol;

             b5 = muTaf*b1;
             b6 = muTaf*b2;
             b7 = muTaf*b3;
             b8 = af*b4;

             c1 = acceFact*b4;

             // GLS stabilisation term

             //D(TIp3,0) = dN_dx[ii];
             //D(TIp3,1) = dN_dy[ii];
             //D(TIp3,2) = dN_dz[ii];

             for(jj=0;jj<totnlbf2;jj++)
             {
               TJ   = ndof*jj;
               TJp1 = TJ+1;
               TJp2 = TJ+2;
               TJp3 = nU+jj;

               // time acceleration term
               fact = c1*N(jj) ;

               // diffusion term
               fact += b5*dN_dx(jj)+b6*dN_dy(jj)+b7*dN_dz(jj);

               Klocal(TI,   TJ)    += fact;
               Klocal(TIp1, TJp1)  += fact;
               Klocal(TIp2, TJp2)  += fact;

               // pressure term
               Klocal(TI,   TJp3)  -= (b1*N(jj));
               Klocal(TIp1, TJp3)  -= (b2*N(jj));
               Klocal(TIp2, TJp3)  -= (b3*N(jj));

               // continuity equation
               Klocal(TIp3, TJ)    -= (b8*dN_dx(jj));
               Klocal(TIp3, TJp1)  -= (b8*dN_dy(jj));
               Klocal(TIp3, TJp2)  -= (b8*dN_dz(jj));

               // stabilisation terms
               Klocal(TIp3, TJp3) -= (b1*dN_dx(jj)+b2*dN_dy(jj)+b3*dN_dz(jj))*stabParam;
             }

             Flocal(TI)   += (b4*(force(0)-velDot(0)) - b1*stress(0,0) - b2*stress(0,1) - b3*stress(0,2) );
             Flocal(TIp1) += (b4*(force(1)-velDot(1)) - b1*stress(1,0) - b2*stress(1,1) - b3*stress(1,2) );
             Flocal(TIp2) += (b4*(force(2)-velDot(2)) - b1*stress(2,0) - b2*stress(2,1) - b3*stress(2,2) );
             Flocal(TIp3) += (b4*F.trace());

             // stabilisation terms
             Flocal(TIp3) += (stabParam*(b1*dp(0)+b2*dp(1)+b3*dp(2)));
          }

          //res = force;

          //res(0) -= (-mu*Du(0) + dp(0));
          //res(1) -= (-mu*Du(1) + dp(1));
          //res(2) -= (-mu*Du(2) + dp(2));

          //res(0) -= dp(0);
          //res(1) -= dp(1);
          //res(2) -= dp(2);

          //dvol *= stabParam;

          //Klocal -= (dvol*(D*D.transpose()));
          //Flocal -= (D*(dvol*res));
    }//gp1
    }//gp2
    }//gp3

    return;
}
*/


/*
template<>
void TreeNode<3>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Stokes flow
    // full subroutine
    //
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, gp3, TI, TIp1, TIp2, TIp3, count, TJ, TJp1, TJp2, TJp3;

    double  uu, vv, ww, JacY, JacZ, Jac, dvol, fact, b1, b2, b3, b4, b5, b6, b7, b8;
    double  xx, yy, zz, pres, afdt, af, am, acceFact, muTaf, c1, Db;

    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf), dN_dz(totnlbf);
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), dNN_dz(totnlbf);
    VectorXd  res(3), res2(3), dp(3), Du(3), vel(3), force(3), velDot(3), conv(3), gradTvel(3);
    MatrixXd  grad(3,3), gradN(3,3), stress(3,3);
    myPoint  param, geom;

    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = am*SolnData->td(9);
    muTaf = mu*af;

    for(gp3=0;gp3<GeomData->GetNGP(2);gp3++)
    {
       param[2]  = 0.5*(knots[2][2] * GeomData->gausspoints3[gp3] + knots[2][3]);
       JacZ = GeomData->gaussweights3[gp3] * JacMultElem;

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       JacY = GeomData->gaussweights2[gp2] * JacZ;

    for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
    {
       param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
       dvol = GeomData->gaussweights1[gp1] * JacY;

       GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

          //N = GeomData->shpfns[level][count].N;
          //dN_dx = GeomData->shpfns[level][count].dN_dx;
          //dN_dy = GeomData->shpfns[level][count].dN_dy;
          //d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          //d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          //count++;

          //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", res(0), res(1), JacMult, Jac, fact, dvol);

          //xx = GeomData->ComputeCoord(0, uu);
          //yy = GeomData->ComputeCoord(1, vv);
          //zz = GeomData->ComputeCoord(2, ww);

          if(parent == NULL)
          {
            N     = NN;
            dN_dx = dNN_dx;
            dN_dy = dNN_dy;
            dN_dz = dNN_dz;
          }
          else
          {
            N     = SubDivMat*NN;
            dN_dx = SubDivMat*dNN_dx;
            dN_dy = SubDivMat*dNN_dy;
            dN_dz = SubDivMat*dNN_dz;
          }

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);
          vel(2) = computeValueCur(2, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(0,2) = computeValueCur(0, dN_dz);

          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);
          grad(1,2) = computeValueCur(1, dN_dz);

          grad(2,0) = computeValueCur(2, dN_dx);
          grad(2,1) = computeValueCur(2, dN_dy);
          grad(2,2) = computeValueCur(2, dN_dz);

          pres   = computeValueCur(3, N);
          //pres   = computeValue2(0, N);
          //dp(0)  = computeValue2(0, dN_dx);
          //dp(1)  = computeValue2(0, dN_dy);
          //dp(2)  = computeValue2(0, dN_dz);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          velDot(2) = computeValueDotCur(2, N);

          //vectmp = d2N_dx2 + d2N_dy2;
          //Du(0) = computeValueCur(0, vectmp);
          //Du(1) = computeValueCur(1, vectmp);
          //Du(2) = computeValueCur(2, vectmp);

          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;
          stress(2,2) -= pres;

          force.setZero();

          res2 = rho*(velDot - force) ;


          for(ii=0;ii<totnlbf2;ii++)
          {
             TI   = ndof*ii;
             TIp1 = TI+1;
             TIp2 = TI+2;
             TIp3 = TI+3;

             b1 = dN_dx[ii]*dvol;
             b2 = dN_dy[ii]*dvol;
             b3 = dN_dz[ii]*dvol;
             b4 = N[ii]*dvol;

             b5 = muTaf*b1;
             b6 = muTaf*b2;
             b7 = muTaf*b3;

             b8 = af*b4;

             c1 = (rho*acceFact)*b4;

             for(jj=0;jj<totnlbf2;jj++)
             {
               TJ   = ndof*jj;
               TJp1 = TJ+1;
               TJp2 = TJ+2;
               TJp3 = TJ+3;

               // time acceleration term
               fact = c1*N(jj) ;

               // diffusion term
               fact += (b5*dN_dx(jj)+b6*dN_dy(jj)+b7*dN_dz(jj));

               Klocal(TI,   TJ)    += fact;
               Klocal(TIp1, TJp1)  += fact;
               Klocal(TIp2, TJp2)  += fact;

               // pressure term
               Klocal(TI,   TJp3)  -= (b1*N(jj));
               Klocal(TIp1, TJp3)  -= (b2*N(jj));
               Klocal(TIp2, TJp3)  -= (b3*N(jj));

               // continuity equation
               Klocal(TIp3, TJ)    -= (b8*dN_dx(jj));
               Klocal(TIp3, TJp1)  -= (b8*dN_dy(jj));
               Klocal(TIp3, TJp2)  -= (b8*dN_dz(jj));
             }

             Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) + b3*stress(0,2) );
             Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) + b3*stress(1,2) );
             Flocal(TIp2) -= (b4*res2(2) + b1*stress(2,0) + b2*stress(2,1) + b3*stress(2,2) );
             Flocal(TIp3) += (b4*grad.trace());

             // stabilisation terms
             //Flocal(TIp3) -= (stabParam*(b1*dp(0)+b2*dp(1)+b3*dp(2)));
          }
    }//gp1
    }//gp2
    }//gp3

    return;
}
*/



//
template<>
void TreeNode<3>::calcStiffnessAndResidualGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
    // GFEM for Naiver-Stokes flow
    // full subroutine
    //
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, gp3, TI, TIp1, TIp2, TIp3, count, TJ, TJp1, TJp2, TJp3;

    double  uu, vv, ww, JacY, JacZ, Jac, dvol, fact, b1, b2, b3, b4, b5, b6, b7, b8;
    double  xx, yy, zz, pres, afdt, af, am, acceFact, muTaf, c1, Db;

    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf), dN_dz(totnlbf);
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), dNN_dz(totnlbf);
    VectorXd  res(3), res2(3), dp(3), Du(3), vel(3), force(3), velDot(3), conv(3), gradTvel(3);
    MatrixXd  grad(3,3), gradN(3,3), stress(3,3);
    myPoint  param, geom;

    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    af = SolnData->td(2);
    am = SolnData->td(3);
    acceFact = am*SolnData->td(9);
    muTaf = mu*af;

    for(gp3=0;gp3<GeomData->GetNGP(2);gp3++)
    {
       param[2]  = 0.5*(knots[2][2] * GeomData->gausspoints3[gp3] + knots[2][3]);
       JacZ = GeomData->gaussweights3[gp3] * JacMultElem;

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       JacY = GeomData->gaussweights2[gp2] * JacZ;

    for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
    {
       param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
       dvol = GeomData->gaussweights1[gp1] * JacY;

       GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

          //N = GeomData->shpfns[level][count].N;
          //dN_dx = GeomData->shpfns[level][count].dN_dx;
          //dN_dy = GeomData->shpfns[level][count].dN_dy;
          //d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          //d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          //count++;

          //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", res(0), res(1), JacMult, Jac, fact, dvol);

          //xx = GeomData->ComputeCoord(0, uu);
          //yy = GeomData->ComputeCoord(1, vv);
          //zz = GeomData->ComputeCoord(2, ww);

          if(parent == NULL)
          {
            N     = NN;
            dN_dx = dNN_dx;
            dN_dy = dNN_dy;
            dN_dz = dNN_dz;
          }
          else
          {
            N     = SubDivMat*NN;
            dN_dx = SubDivMat*dNN_dx;
            dN_dy = SubDivMat*dNN_dy;
            dN_dz = SubDivMat*dNN_dz;
          }

          vel(0) = computeValueCur(0, N);
          vel(1) = computeValueCur(1, N);
          vel(2) = computeValueCur(2, N);

          grad(0,0) = computeValueCur(0, dN_dx);
          grad(0,1) = computeValueCur(0, dN_dy);
          grad(0,2) = computeValueCur(0, dN_dz);

          grad(1,0) = computeValueCur(1, dN_dx);
          grad(1,1) = computeValueCur(1, dN_dy);
          grad(1,2) = computeValueCur(1, dN_dz);

          grad(2,0) = computeValueCur(2, dN_dx);
          grad(2,1) = computeValueCur(2, dN_dy);
          grad(2,2) = computeValueCur(2, dN_dz);

          pres   = computeValueCur(3, N);
          //pres   = computeValue2(0, N);
          //dp(0)  = computeValue2(0, dN_dx);
          //dp(1)  = computeValue2(0, dN_dy);
          //dp(2)  = computeValue2(0, dN_dz);

          velDot(0) = computeValueDotCur(0, N);
          velDot(1) = computeValueDotCur(1, N);
          velDot(2) = computeValueDotCur(2, N);

          //vectmp = d2N_dx2 + d2N_dy2;
          //Du(0) = computeValueCur(0, vectmp);
          //Du(1) = computeValueCur(1, vectmp);
          //Du(2) = computeValueCur(2, vectmp);

          // this is pseudo-stress
          stress = mu*grad;
          stress(0,0) -= pres;
          stress(1,1) -= pres;
          stress(2,2) -= pres;

          force.setZero();

          gradTvel = grad*vel ;

          res2 = rho*(velDot + gradTvel - force) ;


          for(ii=0;ii<totnlbf2;ii++)
          {
             TI   = ndof*ii;
             TIp1 = TI+1;
             TIp2 = TI+2;
             TIp3 = TI+3;

             b1 = dN_dx[ii]*dvol;
             b2 = dN_dy[ii]*dvol;
             b3 = dN_dz[ii]*dvol;
             b4 = N[ii]*dvol;

             b5 = muTaf*b1;
             b6 = muTaf*b2;
             b7 = muTaf*b3;

             b8 = af*b4;

             c1 = (rho*acceFact)*b4;

             for(jj=0;jj<totnlbf2;jj++)
             {
               TJ   = ndof*jj;
               TJp1 = TJ+1;
               TJp2 = TJ+2;
               TJp3 = TJ+3;

               // time acceleration term
               fact = c1*N(jj) ;

               // diffusion term
               fact += (b5*dN_dx(jj)+b6*dN_dy(jj)+b7*dN_dz(jj));

               Klocal(TI,   TJ)    += fact;
               Klocal(TIp1, TJp1)  += fact;
               Klocal(TIp2, TJp2)  += fact;

               // Convection term

               Db = rho*(vel(0)*dN_dx(jj) + vel(1)*dN_dy(jj) + vel(2)*dN_dz(jj));

               gradN = grad*(rho*N(jj));

               gradN(0,0) += Db;
               gradN(1,1) += Db;
               gradN(2,2) += Db;
               gradN *= af;

               Klocal(TI,   TJ)   += (b4*gradN(0,0));
               Klocal(TI,   TJp1) += (b4*gradN(0,1));
               Klocal(TI,   TJp2) += (b4*gradN(0,2));

               Klocal(TIp1, TJ)   += (b4*gradN(1,0));
               Klocal(TIp1, TJp1) += (b4*gradN(1,1));
               Klocal(TIp1, TJp2) += (b4*gradN(1,2));

               Klocal(TIp2, TJ)   += (b4*gradN(2,0));
               Klocal(TIp2, TJp1) += (b4*gradN(2,1));
               Klocal(TIp2, TJp2) += (b4*gradN(2,2));

               // pressure term
               Klocal(TI,   TJp3)  -= (b1*N(jj));
               Klocal(TIp1, TJp3)  -= (b2*N(jj));
               Klocal(TIp2, TJp3)  -= (b3*N(jj));

               // continuity equation
               Klocal(TIp3, TJ)    -= (b8*dN_dx(jj));
               Klocal(TIp3, TJp1)  -= (b8*dN_dy(jj));
               Klocal(TIp3, TJp2)  -= (b8*dN_dz(jj));
             }

             Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) + b3*stress(0,2) );
             Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) + b3*stress(1,2) );
             Flocal(TIp2) -= (b4*res2(2) + b1*stress(2,0) + b2*stress(2,1) + b3*stress(2,2) );
             Flocal(TIp3) += (b4*grad.trace());

             // stabilisation terms
             //Flocal(TIp3) -= (stabParam*(b1*dp(0)+b2*dp(1)+b3*dp(2)));
          }
    }//gp1
    }//gp2
    }//gp3

    return;
}
//



template<>
void TreeNode<3>::calcStiffnessAndResidualLSFEM(bool flag, MatrixXd& Klocal, VectorXd& Flocal)
{
    // LSFEM for Stokes flow (steady-state)
    //
    //////////////////////////////////////////////////

    //TwoDim_Stokes1  analy(1.0);
    //TwoDim_Stokes2  analy(1.0);

    int ii, jj, gp1, gp2, gp3, TI, TIp1, TIp2, TIp3, count, TJ, TJp1, TJp2, TJp3;
   
    double  JacY, JacZ, Jac, dvol, fact, b1, b2, b3, b4, b5, b6, b7, trgradu;
    double  gamma, xx, yy, zz, pres;
    myPoint  param;

    double  rho = elmDat[3];
    double  mu  = elmDat[4];


    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf), dN_dz(totnlbf);
    VectorXd  res(4), dp(3), Du(3), vel(3), force(4);
    VectorXd  d2N_dx2(totnlbf), d2N_dy2(totnlbf), d2N_dz2(totnlbf), vectmp(totnlbf);
    MatrixXd  D(nsize,4);
    D.setZero();
    
    for(gp3=0;gp3<GeomData->GetNGP(2);gp3++)
    {
       param[2]  = 0.5*(knots[2][2] * GeomData->gausspoints3[gp3] + knots[2][3]);
       JacZ = GeomData->gaussweights3[gp3] * JacMultElem;
       //cout << " ww " << ww << endl;

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       param[1]  = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       JacY = GeomData->gaussweights2[gp2] * JacZ;
       //cout << " vv " << vv << endl;
       
    for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
    {
       param[0]   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
       dvol = GeomData->gaussweights1[gp1] * JacY;
       //cout << " uu " << uu << endl;

       GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, 
					N, dN_dx, dN_dy, dN_dz, d2N_dx2, d2N_dy2, d2N_dz2);

          //N = GeomData->shpfns[level][count].N;
          //dN_dx = GeomData->shpfns[level][count].dN_dx;
          //dN_dy = GeomData->shpfns[level][count].dN_dy;
          //d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          //d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          //count++;

          //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", res(0), res(1), JacMult, Jac, fact, dvol);

          //xx = GeomData->ComputeCoord(0, param[0]);
          //yy = GeomData->ComputeCoord(1, param[1]);
          //zz = GeomData->ComputeCoord(2, param[2]);

          vel(0) = computeValue(0, N);
          vel(1) = computeValue(1, N);
          vel(2) = computeValue(2, N);
          pres   = computeValue(3, N);

          vectmp = d2N_dx2 + d2N_dy2 + d2N_dz2;
          Du(0) = computeValue(0, vectmp);
          Du(1) = computeValue(1, vectmp);
          Du(2) = computeValue(2, vectmp);

          dp(0)  = computeValue(3, dN_dx);
          dp(1)  = computeValue(3, dN_dy);
          dp(2)  = computeValue(3, dN_dz);

          trgradu = computeValue(0, dN_dx);
          trgradu += computeValue(1, dN_dy);
          trgradu += computeValue(2, dN_dz);


          force.setZero();
          //force(0) = analy.computeXForce(uu, vv);
          //force(1) = analy.computeYForce(uu, vv);
          //force(0) = -0.3*delta;

          //force(0) = computeForce(0, N);
          //force(1) = computeForce(1, N);

          for(ii=0;ii<totnlbf;ii++)
          {
             TI = ndof*ii;
             TIp1 = TI+1;
             TIp2 = TI+2;
             TIp3 = TI+3;

             fact = -mu * vectmp[ii];

             D(TI,  0) = fact;       D(TI,  1) = 0.0;        D(TI,  2) = 0.0;        D(TI,  3) = dN_dx[ii];
             D(TIp1,0) = 0.0;        D(TIp1,1) = fact ;      D(TIp1,2) = 0.0;        D(TIp1,3) = dN_dy[ii];
             D(TIp2,0) = 0.0;        D(TIp2,1) = 0.0;        D(TIp2,2) = fact ;      D(TIp2,3) = dN_dz[ii];
             D(TIp3,0) = dN_dx[ii];  D(TIp3,1) = dN_dy[ii];  D(TIp3,2) = dN_dz[ii];  D(TIp3,3) = 0.0 ;
          }

          //res.setZero();
          res = force;

          //res(0) = analy.computeXForce(uu, vv);
          //res(1) = analy.computeYForce(uu, vv);

          //res(0) = computeForce(0, N) ;
          //res(1) = computeForce(1, N) ;

          res(0) -= (-mu*Du(0) + dp(0));
          res(1) -= (-mu*Du(1) + dp(1));
          res(2) -= (-mu*Du(2) + dp(2));
          res(3) -= trgradu;

          Klocal += (dvol*(D*D.transpose()));
          Flocal += (D*(dvol*res));
    }//gp1
    }//gp2
    }//gp3
    
    //if(id == 555)
//    printMatrix(Klocal);
//    printf("\n\n");
//    printVector(Flocal);
//    printf("\n\n");

    return;
}


/*
template<>
void TreeNode<3>::calcStiffnessAndResidual(int ind1, int ind2, double inp1, double inp2)
{
    // LSFEM for Navier-Stokes flow (steady-state)
    //
    //////////////////////////////////////////////////

    int ii, jj, gp1, gp2, gp3, TI, TIp1, TIp2, TIp3, count, TJ, TJp1, TJp2, TJp3;
   
    double  uu, vv, ww, JacY, JacZ, Jac, dvol, fact, b1, b2, b3, b4, b5, b6, b7;
    double  xx, yy, zz, pres, trgradu;

    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf), dN_dz(totnlbf);
    VectorXd  res(4), dp(3), Du(3), vel(3), force(4);
    VectorXd  d2N_dx2(totnlbf), d2N_dy2(totnlbf), d2N_dz2(totnlbf), vectmp(totnlbf);
    MatrixXd  D(nsize,4), F(3,3), FN(3,3);
    D.setZero();

    for(gp3=0;gp3<GeomData->GetNGP(2);gp3++)
    {
       ww   = 0.5*(knots[2][2] * GeomData->gausspoints3[gp3] + knots[2][3]);
       JacZ = GeomData->gaussweights3[gp3] * JacMult;
       //cout << " ww " << ww << endl;

    for(gp2=0;gp2<GeomData->GetNGP(1);gp2++)
    {
       vv   = 0.5*(knots[1][2] * GeomData->gausspoints2[gp2] + knots[1][3]);
       JacY = GeomData->gaussweights2[gp2] * JacZ;
       //cout << " vv " << vv << endl;

    for(gp1=0;gp1<GeomData->GetNGP(0);gp1++)
    {
       uu   = 0.5*(knots[0][2] * GeomData->gausspoints1[gp1] + knots[0][3]);
       dvol = GeomData->gaussweights1[gp1] * JacY;
       //cout << " uu " << uu << endl;

       //computeBasisFunctions3D(knots[0][0], knots[1][0], knots[2][0], knots[0][2], knots[1][2], knots[2][2], uu, vv, ww, &N(0), &dN_dx(0), &dN_dy(0), &(dN_dz(0)));
       computeBasisFunctions3D(knots[0][0], knots[1][0], knots[2][0], knots[0][2], knots[1][2], knots[2][2], 
					 uu, vv, ww, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0),
					 &d2N_dx2(0), &d2N_dy2(0), &d2N_dz2(0));

          //N = GeomData->shpfns[level][count].N;
          //dN_dx = GeomData->shpfns[level][count].dN_dx;
          //dN_dy = GeomData->shpfns[level][count].dN_dy;
          //d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
          //d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          //count++;

          //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", res(0), res(1), JacMult, Jac, fact, dvol);

          //xx = GeomData->ComputeCoord(0, uu);
          //yy = GeomData->ComputeCoord(1, vv);
          //zz = GeomData->ComputeCoord(2, ww);

          vel(0) = computeValue(0, N);
          vel(1) = computeValue(1, N);
          vel(2) = computeValue(2, N);
          pres   = computeValue(3, N);

          vectmp = d2N_dx2 + d2N_dy2 + d2N_dz2;
          Du(0) = computeValue(0, vectmp);
          Du(1) = computeValue(1, vectmp);
          Du(2) = computeValue(2, vectmp);

          dp(0)  = computeValue(3, dN_dx);
          dp(1)  = computeValue(3, dN_dy);
          dp(2)  = computeValue(3, dN_dz);

          F(0,0) = computeValue(0, dN_dx);
          F(0,1) = computeValue(0, dN_dy);
          F(0,2) = computeValue(0, dN_dz);
          F(1,0) = computeValue(1, dN_dx);
          F(1,1) = computeValue(1, dN_dy);
          F(1,2) = computeValue(1, dN_dz);
          F(2,0) = computeValue(2, dN_dx);
          F(2,1) = computeValue(2, dN_dy);
          F(2,2) = computeValue(2, dN_dz);

          trgradu = F.trace();

          for(ii=0;ii<totnlbf;ii++)
          {
             TI = ndof*ii;
             TIp1 = TI+1;
             TIp2 = TI+2;
             TIp3 = TI+3;

             b1 = dN_dx[ii];
             b2 = dN_dy[ii];
             b3 = dN_dz[ii];

             fact = rho*(vel(0)*b1 + vel(1)*b2 + vel(2)*b3) - mu*vectmp[ii];
             
             FN = rho * N[ii] * F;
             
             D(TI,  0) = FN(0,0) + fact;  D(TI,  1) = FN(1,0);         D(TI,  2) = FN(2,0);         D(TI,  3) = b1;
             D(TIp1,0) = FN(0,1);         D(TIp1,1) = FN(1,1) + fact;  D(TIp1,2) = FN(2,1);         D(TIp1,3) = b2;
             D(TIp2,0) = FN(0,2);         D(TIp2,1) = FN(1,2);         D(TIp2,2) = FN(2,2) + fact;  D(TIp2,3) = b3;
             D(TIp3,0) = b1;              D(TIp3,1) = b2;              D(TIp3,2) = b3;              D(TIp3,3) = 0.0 ;
          }

          force.setZero();
          res = force;

          //res(0) += computeForce(0, N) ;
          //res(1) += computeForce(1, N) ;
          //res(2) += computeForce(2, N) ;

          res(0) -= (rho*(F(0,0)*vel(0) + F(0,1)*vel(1) + F(0,2)*vel(2)) - mu*Du(0) + dp(0));
          res(1) -= (rho*(F(1,0)*vel(0) + F(1,1)*vel(1) + F(1,2)*vel(2)) - mu*Du(1) + dp(1));
          res(2) -= (rho*(F(2,0)*vel(0) + F(2,1)*vel(1) + F(2,2)*vel(2)) - mu*Du(2) + dp(2));
          res(3) -= trgradu;

          Klocal += (dvol*(D*D.transpose()));
          Flocal += (D*(dvol*res));
    }//gp1
    }//gp2
    }//gp3
    
    //if(id == 555)
//    printMatrix(Klocal);
//    printf("\n\n");
//    printVector(Flocal);
//    printf("\n\n");

    return;
}
*/

template<>
void TreeNode<3>::applyBoundaryConditionsAtApoint(myDataIntegrateCutFEM& myData)
{
    int ii, jj, TI, nU;
   
    double  fact, af;
    VectorXd  NN(totnlbf), N;

    af = SolnData->td(2);

    GeomData->computeBasisFunctions3D(knotBegin, knotIncr, myData.param, NN);
    
    if(parent == NULL)
      N = NN;
    else
      N = SubDivMat*NN;

    //printf("\t param %6d \t %6d \t %12.8f \t %12.8f \t %12.8f \n", id, id, param[0], param[1], param[2]);
    //printf("\n\n");
    //printVector(N);
    //printf("\n\n");

  //for(int dd=0;dd<3;dd++)
  //{
    //dir = dd;
    //specVal = 0.0;

    if(myData.dir < 4)
    {
      myData.specVal[0] -= computeValueCur(myData.dir, N);

      for(ii=0;ii<totnlbf2;ii++)
      {
        fact  = N[ii] * myData.PENALTY;

        TI = ndof*ii+myData.dir;

        myData.F1(TI) += (fact*myData.specVal[0]);
        fact *= af;

        for(jj=0;jj<totnlbf2;jj++)
          myData.K1(TI, ndof*jj+myData.dir) += fact * N[jj];
      }
    }
    if(myData.dir == 4)
    {
      nU  = forAssyVec.size();

      myData.specVal[0] -= computeValue2(0, N);

      for(ii=0;ii<totnlbf2;ii++)
      {
        fact = N[ii] * myData.PENALTY;

        TI = nU+ii;

        myData.F1(TI) += (fact*myData.specVal[0]);
        //fact *= af;

        for(jj=0;jj<totnlbf2;jj++)
        {
          myData.K1(TI, nU+jj) += fact * N[jj];
        }
      }
    }
  //}

    return;
}



template<>
double TreeNode<3>::getJacBoundary(int side)
{
  double val;

  switch(side)
  {
      case 0:
          val = 0.25 * knots[1][2] * knots[2][2] * GeomData->GetJacobian(1) * GeomData->GetJacobian(2);
      break;

      case 1:
          val = 0.25 * knots[1][2] * knots[2][2] * GeomData->GetJacobian(1) * GeomData->GetJacobian(2);
      break;

      case 2:
          val = 0.25 * knots[0][2] * knots[2][2] * GeomData->GetJacobian(0) * GeomData->GetJacobian(2);
      break;

      case 3:
          val = 0.25 * knots[0][2] * knots[2][2] * GeomData->GetJacobian(0) * GeomData->GetJacobian(2);
      break;

      case 4:
          val = 0.25 * knots[0][2] * knots[1][2] * GeomData->GetJacobian(0) * GeomData->GetJacobian(1);
      break;

      case 5:
          val = 0.25 * knots[0][2] * knots[1][2] * GeomData->GetJacobian(0) * GeomData->GetJacobian(1);
      break;

      default :
          cout << " Invalid 'side' value in TreeNode<3>::getJacBoundary " << endl;
      break;
  } //switch(side)

  return val;
}


/*

            if(side == 0)
            {
              rad = sqrt((yy-yc)*(yy-yc)+(zz-zc)*(zz-zc));
              if(dir == 10)
              {
                if(rad > 0.25)
                  res = 0.0;
                else
                  res = specVal*(y1-yy)*(yy-y0)*(z1-zz)*(zz-z0);
                //res = 0.25 - rad;

                //cout << xx << '\t' << yy << '\t' << y0 << '\t' << y1 << '\t' << res << endl;
              }
              else if(dir == 0)
              {
                if(yy < y0 || yy > y1)
                  res = 0.0;
                else
                  res = specVal*(y1-yy)*(yy-y0);

                //cout << xx << '\t' << yy << '\t' << y0 << '\t' << y1 << '\t' << res << endl;
              }
              else
                res = specVal;
            }
            else if(side == 33)
            {
              if(dir == 0)
              {
                //res = sin(PI*xx);
                if(uu<0.5)
                  res = tanh(50.0*uu);
                else
                  res = -tanh(50.0*(uu-1.0));
              }
              else
                res = specVal;
            }
            else
              res = specVal;
*/

template<>
void TreeNode<3>::applyDirichletBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  if( DirichletData.size() > 0 )
  {
    int ii, jj, kk, gp1, gp2, gp3, TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3, nGP1, nGP2, nGP3, dir, side, index, nU, nP, aa;
    double  y0, y1, z0, z1, Ta, Tb, JacMultLoc, xc, yc, zc, rad, JacZ, JacY, bb1, bb2, bb3, af;
    double  xx, yy, zz, res, dvol, specVal, PENALTY, Jac, fact, fact1, fact2, R, NitscheFact, pres;
    bool isNitsche;

    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf), dN_dz(totnlbf);
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), dNN_dz(totnlbf);
    MatrixXd  grad(3,3), stress(3,3);
    myPoint  param, normal, trac;
    vector<double>  boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2, boundaryGPs3, boundaryGWs3;

    y0 = 0.0;
    y1 = 2.0;
    z0 = 0.25;
    z1 = 0.75;

    //theta = 0.0;
    //theta = 0.463647609;
    //y0 = 0.0;
    //y1 = 0.5;

    //y0 = 0.1;
    //y1 = 0.6;

    //y0 = 0.0;
    //y1 = 1.0;

    xc = 0.0;
    yc = 0.5;
    zc = 0.5;
    R = 0.25;
    
    Circle  circ(0.5,0.5,0.25);

    PoissonEx2 analy;
    
    af = 1.0;

      for(aa=0;aa<DirichletData.size();aa++)
      {
        //printVector(DirichletData[aa]);

        side    = (int) (DirichletData[aa][0] - 1);
        dir     = (int) (DirichletData[aa][1] - 1);
        specVal = DirichletData[aa][2];
        PENALTY = DirichletData[aa][3];
        isNitsche   = ( (int) DirichletData[aa][4] == 1 );
        NitscheFact = DirichletData[aa][5];


        normal = GeomData->boundaryNormals[side];

        //GeomData->getBoundaryNormal3D(side, normal);
        GeomData->setBoundaryGPs3D(side, boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2, boundaryGPs3, boundaryGWs3);

        JacMultLoc = TreeNode<3>::getJacBoundary(side);

        for(gp3=0;gp3<boundaryGPs3.size();gp3++)
        {
            param[2] = 0.5*(knots[2][2] * boundaryGPs3[gp3] + knots[2][3]);
            JacZ = boundaryGWs3[gp3] * JacMultLoc;

        for(gp2=0;gp2<boundaryGPs2.size();gp2++)
        {
            param[1] = 0.5*(knots[1][2] * boundaryGPs2[gp2] + knots[1][3]);
            JacY = boundaryGWs2[gp2] * JacZ;

        for(gp1=0;gp1<boundaryGPs1.size();gp1++)
        {
            param[0] = 0.5*(knots[0][2] * boundaryGPs1[gp1] + knots[0][3]);
            dvol = boundaryGWs1[gp1] * JacY;

            GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

            if(parent == NULL)
            {
              N = NN;
              dN_dx = dNN_dx;
              dN_dy = dNN_dy;
              dN_dz = dNN_dz;
            }
            else
            {
              N = SubDivMat*NN;
              dN_dx = SubDivMat*dNN_dx;
              dN_dy = SubDivMat*dNN_dy;
              dN_dz = SubDivMat*dNN_dz;
            }

            //printf("\n\n tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n\n\n", xx, yy, val, val, Jac, JacMult, dvol);
            //if(pout) printf(" knotsAtGPs = %12.6f \t xx = %12.6f \t yy = %12.6f \n", knotsAtGPs, xx, yy);
            //printf(" uu = %12.6f \t vv = %12.6f \t dvol = %12.6f \t volume = %12.6f \n", uu, vv, dvol, volume);

            xx = GeomData->ComputeCoord(0, param[0]);
            yy = GeomData->ComputeCoord(1, param[1]);
            zz = GeomData->ComputeCoord(2, param[2]);

            //rad = sqrt((yy-yc)*(yy-yc)+(zz-zc)*(zz-zc));
            //r = sqrt(xx*xx+yy*yy);
            //val = 1.0+log(2.0*r);
            //cout << xx << '\t' << yy << endl;

            res = specVal;

            //
            if(side == 0)
            {
              if(dir == 0)
              {
                //if(yy <= y0 || yy >= y1)
                  //res = 0.0;
                //else
                  res = specVal*(y1-yy)*(yy-y0);
              }
            }
            //
            /*
            if(side == 0)
            {
              if(dir == 0)
              {
                if(circ.checkPointLocation(yy, zz))
                  res = specVal*(1.0-rad*rad/R/R);
                else
                  res = 0.0;                

                //if(rad >= 0.25)
                  //res = 0.0;
                //else
                  //res = specVal*(y1-yy)*(yy-y0)*(z1-zz)*(zz-z0);
                  //res = specVal*(1.0-rad*rad/R/R);
                //res = 0.25 - rad;

                //cout << xx << '\t' << yy << '\t' << y0 << '\t' << y1 << '\t' << res << endl;
              }
              else
                res = specVal;
            }
            else
              res = specVal;
            */
            //res = analy.computeValue(dir, xx, yy, zz);


            //cout << side << '\t' << dir << '\t' << DirichletData[side][dir] << '\t' << timeFunction[0].prop << endl;

            res *= timeFunction[0].prop;
            //res *= tanh(1.0*mpapTime.cur);
            //res *= sin(2.0*PI*mpapTime.cur+PI/2.0);
            //res *= sin(2.0*PI*10.0*mpapTime.cur+PI/2.0);

            //cout << res << '\t' << computeValue(dir, N) << endl;
            //res -= computeValue(dir, N);
            //cout << res << endl;
            //cout << side << '\t' << dir << endl;
            //printf("%12.6f \t %12.6f \t %12.6f \n", timeFunction[0].prop, res, DirichletData[side][dir]);

              res -= computeValue(dir, N);
              //cout << side << '\t' << dir << endl;
              //printf("%12.6f \t %12.6f \t %12.6f \n", timeFunction[0].prop, res, DirichletData[side][dir]);


              grad(0,0) = computeValue(0, dN_dx);
              grad(0,1) = computeValue(0, dN_dy);
              grad(0,2) = computeValue(0, dN_dz);

              grad(1,0) = computeValue(1, dN_dx);
              grad(1,1) = computeValue(1, dN_dy);
              grad(1,2) = computeValue(1, dN_dz);

              grad(2,0) = computeValue(2, dN_dx);
              grad(2,1) = computeValue(2, dN_dy);
              grad(2,2) = computeValue(2, dN_dz);

              pres      = computeValue(3, N);

              stress = mu*grad;
              stress(0,0) -= pres;
              stress(1,1) -= pres;
              stress(2,2) -= pres;

              trac = stress*normal ;

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

              // terms corresponding to Nitsche method for Stokes and Navier-Stokes
              if(isNitsche)
              {
                if(dir == 0)
                {
                  for(ii=0;ii<totnlbf2;ii++)
                  {
                    TI   = ndof*ii;
                    TIp1 = TI+1;
                    TIp2 = TI+2;
                    TIp3 = TI+3;

                    Ta = mu*(normal[0]*dN_dx(ii) + normal[1]*dN_dy(ii) + normal[2]*dN_dz(ii))*dvol;
                    bb1 = N[ii]*dvol;

                    for(jj=0;jj<totnlbf2;jj++)
                    {
                      TJ   = ndof*jj;
                      TJp1 = TJ+1;
                      TJp2 = TJ+2;
                      TJp3 = TJ+3;

                      Tb = af*mu*( normal[0]*dN_dx(jj) + normal[1]*dN_dy(jj) + normal[2]*dN_dz(jj) );

                      fact = bb1*(-normal[0]*N(jj));

                      Klocal(TI, TJ)   -= (bb1*Tb);
                      Klocal(TI, TJp3) -= fact;

                      // Nitsche terms
                      Klocal(TI, TJ)   -= (Ta*af*N(jj))*NitscheFact;
                      Klocal(TIp3, TJ) -= fact*af*NitscheFact;
                    }

                    Flocal(TI)   -= (-bb1*trac[0]);
                    Flocal(TIp1) -= (0.0);
                    Flocal(TIp2) -= (0.0);

                    Flocal(TI)   -= (Ta*res)*NitscheFact;
                    Flocal(TIp3) -= (bb1*(-normal[0]*res))*NitscheFact;
                  }
                }
                if(dir == 1)
                {
                  for(ii=0;ii<totnlbf2;ii++)
                  {
                    TI   = ndof*ii;
                    TIp1 = TI+1;
                    TIp2 = TI+2;
                    TIp3 = TI+3;

                    Ta = mu*(normal[0]*dN_dx(ii) + normal[1]*dN_dy(ii) + normal[2]*dN_dz(ii))*dvol;
                    bb1 = N[ii]*dvol;

                    for(jj=0;jj<totnlbf2;jj++)
                    {
                      TJ   = ndof*jj;
                      TJp1 = TJ+1;
                      TJp2 = TJ+2;
                      TJp3 = TJ+3;

                      Tb = af*mu*( normal[0]*dN_dx(jj) + normal[1]*dN_dy(jj) + normal[2]*dN_dz(jj) );

                      fact = bb1*(-normal[1]*N(jj));

                      Klocal(TIp1, TJp1)  -= (bb1*Tb);
                      Klocal(TIp1, TJp3)  -= fact;

                      // Nitsche terms
                      Klocal(TIp1, TJp1)  -= (Ta*af*N(jj))*NitscheFact;
                      Klocal(TIp3, TJp1)  -= fact*af*NitscheFact;
                    }

                    Flocal(TI)   -= (0.0);
                    Flocal(TIp1) -= (-bb1*trac[1]);
                    Flocal(TIp2) -= (0.0);

                    Flocal(TIp1) -= (Ta*res)*NitscheFact;
                    Flocal(TIp3) -= (bb1*(-normal[1]*res))*NitscheFact;
                  }
                }
                if(dir == 2)
                {
                  for(ii=0;ii<totnlbf2;ii++)
                  {
                    TI   = ndof*ii;
                    TIp1 = TI+1;
                    TIp2 = TI+2;
                    TIp3 = TI+3;

                    Ta = mu*(normal[0]*dN_dx(ii) + normal[1]*dN_dy(ii) + normal[2]*dN_dz(ii))*dvol;
                    bb1 = N[ii]*dvol;

                    for(jj=0;jj<totnlbf2;jj++)
                    {
                      TJ   = ndof*jj;
                      TJp1 = TJ+1;
                      TJp2 = TJ+2;
                      TJp3 = TJ+3;

                      Tb = af*mu*( normal[0]*dN_dx(jj) + normal[1]*dN_dy(jj) + normal[2]*dN_dz(jj) );

                      fact = bb1*(-normal[2]*N(jj));

                      Klocal(TIp2, TJp2)  -= (bb1*Tb);
                      Klocal(TIp2, TJp3)  -= fact;

                      // Nitsche terms
                      Klocal(TIp2, TJp2)  -= (Ta*af*N(jj))*NitscheFact;
                      Klocal(TIp3, TJp2)  -= fact*af*NitscheFact;
                    }

                    Flocal(TI)   -= (0.0);
                    Flocal(TIp1) -= (0.0);
                    Flocal(TIp2) -= (-bb1*trac[2]);

                    Flocal(TIp2) -= (Ta*res)*NitscheFact;
                    Flocal(TIp3) -= (bb1*(-normal[2]*res))*NitscheFact;
                  }
                }
              } //if(isNitsche)
        }// for(gp1=0...
        }// for(gp2=0...
        }// for(gp3=0...
      } // for(aa=0;aa<DirichletData.size();aa++)
  } // if(DirichletData.size() > 0)

  return;
}

/*
template<>
void TreeNode<3>::applyNeumannBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  if( NeumannData.size() > 0 )
  {
    int ii, jj, kk, gp1, gp2, gp3, TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3, dir, side, nU, nP, aa;
    double theta, y0, y1, z0, z1, Ta, Tb, xc, yc, zc, rad;
    double  xx, yy, zz, res, dvol, specVal, PENALTY, r, fact, fact1, fact2;
    double  JacY, JacZ, JacMultLoc;

    double  rho = elmDat[3];
    double  mu  = elmDat[4];

    VectorXd  N(totnlbf), dN_dx(totnlbf), dN_dy(totnlbf), dN_dz(totnlbf);
    VectorXd  NN(totnlbf), dNN_dx(totnlbf), dNN_dy(totnlbf), dNN_dz(totnlbf);
    myPoint  param;

    vector<double>  boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2, boundaryGPs3, boundaryGWs3;

    PENALTY  = 1.0;

    mu = elmDat[4];

      for(aa=0;aa<NeumannData.size();aa++)
      {
        side    = (int) (NeumannData[aa][0] - 1);
        dir     = (int) (NeumannData[aa][1] - 1);
        specVal = NeumannData[aa][2];

        GeomData->setBoundaryGPs3D(side, boundaryGPs1, boundaryGWs1, boundaryGPs2, boundaryGWs2, boundaryGPs3, boundaryGWs3);

        JacMultLoc = TreeNode<3>::getJacBoundary(side);

        for(gp3=0;gp3<boundaryGPs3.size();gp3++)
        {
            param[2] = 0.5*(knots[2][2] * boundaryGPs3[gp3] + knots[2][3]);
            JacZ = boundaryGWs3[gp3] * JacMultLoc;

        for(gp2=0;gp2<boundaryGPs2.size();gp2++)
        {
            param[1] = 0.5*(knots[1][2] * boundaryGPs2[gp2] + knots[1][3]);
            JacY = boundaryGWs2[gp2] * JacZ;

        for(gp1=0;gp1<boundaryGPs1.size();gp1++)
        {
            param[0] = 0.5*(knots[0][2] * boundaryGPs1[gp1] + knots[0][3]);
            dvol = boundaryGWs1[gp1] * JacY;

        GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

        if(parent == NULL)
        {
          N = NN;
          dN_dx = dNN_dx;
          dN_dy = dNN_dy;
          dN_dz = dNN_dz;
        }
        else
        {
          N = SubDivMat*NN;
          dN_dx = SubDivMat*dNN_dx;
          dN_dy = SubDivMat*dNN_dy;
          dN_dz = SubDivMat*dNN_dz;
        }

        //printf("\n\n tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n\n\n", xx, yy, val, val, Jac, JacMult, dvol);
        //if(pout) printf(" knotsAtGPs = %12.6f \t xx = %12.6f \t yy = %12.6f \n", knotsAtGPs, xx, yy);
        //printf(" uu = %12.6f \t vv = %12.6f \t dvol = %12.6f \t volume = %12.6f \n", uu, vv, dvol, volume);

        xx = GeomData->ComputeCoord(0, param[0]);
        yy = GeomData->ComputeCoord(1, param[1]);
        zz = GeomData->ComputeCoord(2, param[2]);
        //r = sqrt(xx*xx+yy*yy);
        //val = 1.0+log(2.0*r);
        //cout << xx << '\t' << yy << endl;
        specVal = NeumannData[aa][2];

        specVal *= timeFunction[0].prop;
        //res *= (0.5*(1.0-cos(2.0*PI*SolnData->ElemProp.data[6]*mpapTime.cur)));
        //w1 = tanh(SolnData->ElemProp.data[5]*mpapTime.cur);
        //res *= w1;
        //w1 = 1.0;
        //w2 = 2.0*PI*SolnData->ElemProp.data[6];
        //res *= sin(w2*mpapTime.cur);
        //res *= tanh(50.0*timeFunction[0].prop);

        if(side == 10)
        {
          if(dir < 2)
          {
            if(yy <= -y0)
            {
              res = specVal*dvol;
              //cout << " yy " << yy << '\t' << res << endl;
              for(ii=0;ii<totnlbf2;ii++)
                Flocal(ndof*ii+dir) += (res*N(ii));
            }
            if(yy >= y0)
            {
              res = 0.0*dvol;
              //cout << " yy " << yy << '\t' << res << endl;
              for(ii=0;ii<totnlbf2;ii++)
                Flocal(ndof*ii+dir) += (res*N(ii));
            }
          }
        }
        else
        {
          res *= dvol;
          //cout << " yy " << yy << '\t' << res << endl;
          for(ii=0;ii<totnlbf2;ii++)
            Flocal(ndof*ii+dir) += (res*N(ii));
        }
    }// for(gp1=0...
    }// for(gp2=0...
    }// for(gp3=0...
      } // for(aa=0;aa<NeumannData.size();aa++)
  } // if(NeumannData.size() > 0)

  return;
}
*/


template<>
void TreeNode<3>::applyNeumannBCsGFEM(MatrixXd& Klocal, VectorXd& Flocal, int domainCur)
{
  if( NeumannData.size() > 0 )
  {
      int  aa, ii, jj, gp, nGauss, TI, TIp1, TIp2, index, dir, side;
      double  y0, y1, res, JacTemp;
      double  dvol, specVal, freq;

      VectorXd  N(totnlbf), NN(totnlbf);
      myPoint  param, normal, geom;
      double *gws;
      myPoint *gps;

      for(aa=0;aa<NeumannData.size();aa++)
      {
        side    = (int) (NeumannData[aa][0] - 1);
        dir     = (int) (NeumannData[aa][1] - 1);
        specVal = NeumannData[aa][2];

        normal = GeomData->boundaryNormals[side];

          nGauss = GeomData->boundaryQuadrature3D[side].gausspoints.size();

          gps = &(GeomData->boundaryQuadrature3D[side].gausspoints[0]);
          gws = &(GeomData->boundaryQuadrature3D[side].gaussweights[0]);

          JacTemp = GeomData->boundaryJacobians[side][level];

        for(gp=0; gp<nGauss; gp++)
        {
            param[2] = 0.5 * (knots[2][2] * gps[gp][2] + knots[2][3]);
            param[1] = 0.5 * (knots[1][2] * gps[gp][1] + knots[1][3]);
            param[0] = 0.5 * (knots[0][2] * gps[gp][0] + knots[0][3]);

            dvol = JacTemp * gws[gp] ;

            //GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

            GeomData->computeBasisFunctions3D(knotBegin, knotIncr, param, NN);

            if(parent == NULL)
            {
              N = NN;
            }
            else
            {
              N = SubDivMat*NN;
            }

            //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %14.8f \n", NN[ii]);

            //printf("\n\n tracX and tracY ... %12.6f \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f  \t %12.6f \n\n\n", xx, yy, val, val, Jac, JacMult, dvol);
            //if(pout) printf(" knotsAtGPs = %12.6f \t xx = %12.6f \t yy = %12.6f \n", knotsAtGPs, xx, yy);
            //printf(" uu = %12.6f \t vv = %12.6f \t dvol = %12.6f \t volume = %12.6f \n", uu, vv, dvol, volume);

            geom[0] = GeomData->ComputeCoord(0, param[0]);
            geom[1] = GeomData->ComputeCoord(1, param[1]);
            geom[2] = GeomData->ComputeCoord(2, param[2]);

            //r = sqrt(xx*xx+yy*yy);
            //val = 1.0+log(2.0*r);
            //cout << xx << '\t' << yy << endl;

            specVal = NeumannData[aa][2];

            res = dvol * specVal * timeFunction[0].prop;
            //res = dvol * specVal * tanh(20.0*mpapTime.cur);
            //cout << specVal << '\t' << axsy << '\t' << dir << endl;

            //res *= (0.5*(1.0-cos(2.0*PI*SolnData->ElemProp.data[6]*mpapTime.cur)));
            //w1 = tanh(SolnData->ElemProp.data[5]*mpapTime.cur);

            freq = 2.0*PI*10.0;
            //res = dvol* specVal * (0.5*(1.0-cos(freq*mpapTime.cur)));
            //res = dvol * specVal * sin(freq*mpapTime.cur);
            //res *= tanh(50.0*timeFunction[0].prop);

            //cout << " yy " << yy << '\t' << res << endl;
            for(ii=0;ii<totnlbf2;ii++)
              Flocal(ndof*ii+dir) += (res*N(ii));

        }// for(gp=0...
      } // for(aa=0;aa<NeumannData.size();aa++)
  } // if(NeumannData.size() > 0)

  return;
}





template<>
void TreeNode<3>::applyDirichletBCsLSFEM(bool flag, MatrixXd& Klocal, VectorXd& Flocal)
{
  return;
}
template<>
void TreeNode<3>::applyNeumannBCsLSFEM(bool flag, MatrixXd& Klocal, VectorXd& Flocal)
{
  return;
}




