
#include "TreeNode.h"
#include "MpapTime.h"
#include "Functions.h"
#include "TimeFunction.h"
#include "SolutionData.h"

#include "BasisFunctionsBSpline.h"

extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;


template<>
int TreeNode<2>::calcLoadVector(int ind1, int ind2, double inp1, double inp2)
{
   cout << " TreeNode<2>::calcLoadVector() ...... yet to be implemented " << endl;

  return 0;
}




          /*
          if(parent == NULL)
          {
            N       = GeomData->shpfns[level][count].N;
            dN_dx   = GeomData->shpfns[level][count].dN_dx;
            dN_dy   = GeomData->shpfns[level][count].dN_dy;
            d2N_dx2 = GeomData->shpfns[level][count].d2N_dx2;
            d2N_dy2 = GeomData->shpfns[level][count].d2N_dy2;
          }
          else
          {
            N       = SubDivMat*GeomData->shpfns[level][count].N;
            dN_dx   = SubDivMat*GeomData->shpfns[level][count].dN_dx;
            dN_dy   = SubDivMat*GeomData->shpfns[level][count].dN_dy;
            d2N_dx2 = SubDivMat*GeomData->shpfns[level][count].d2N_dx2;
            d2N_dy2 = SubDivMat*GeomData->shpfns[level][count].d2N_dy2;
          }
          count++;
          */




template<>
void TreeNode<2>::mapBoundaryPointDataToGlobalBodyForceVector(double* position, double* normal, double arclen, double* useThisData)
{
    int      ii, jj, gp1, gp2, TI;
    double   Jac, dvol, fact, xx, yy;
    double   beta1, beta2, delta1, delta2, delta, h1, h2;

    VectorXd  N(totnlbf), res(ndof);
    MatrixXd  D(forAssyVec.size(), ndof);
    myPoint  param;

    double  val1  = 0.5*knotIncr[0];
    double  val2  = 0.5*knotSum[0];
    double  val3  = 0.5*knotIncr[1];
    double  val4  = 0.5*knotSum[1];
    double  *gausspoints1  = &(GeomData->gausspoints1[0]);
    double  *gaussweights1 = &(GeomData->gaussweights1[0]);
    double  *gausspoints2  = &(GeomData->gausspoints2[0]);
    double  *gaussweights2 = &(GeomData->gaussweights2[0]);

    //JacMult = GeomData->getJacobianFull() * val1 * val3;

    beta1 = elmDat[5] * GeomData->getGridLength(0) * knotIncr[0];
    beta2 = elmDat[5] * GeomData->getGridLength(1) * knotIncr[1];

    //h1 = GeomData->ElemProp.data[5] * GeomData->getGridLength(0)*knotIncr[0];
    //h2 = GeomData->ElemProp.data[5] * GeomData->getGridLength(1)*knotIncr[1];

    h1 = GeomData->getGridLength(0) * knotIncr[0];
    h2 = GeomData->getGridLength(1) * knotIncr[1];

    //printf("\t data \t %12.6f \t %12.6f \n", h1, h2);
    
    fact = arclen/h1/h2;
    fact = 1.0/h1/h2;
    //fact = 1.0;
    
    res.setZero();
    
    if(ndof == 1)
      res(0) = useThisData[0];
    else
    {
      res(0) = useThisData[0];
      res(1) = useThisData[1];
    }

    for(gp2=0;gp2<GeomData->getNGP(1);gp2++)
    {
       param[1]  = val3 * gausspoints2[gp2] + val4;
       Jac = gaussweights2[gp2] * JacMultElem;

       for(gp1=0;gp1<GeomData->getNGP(0);gp1++)
       {
          param[0]  = val1 * gausspoints1[gp1] + val2;
          dvol = gaussweights1[gp1] * Jac;

          GeomData->computeBasisFunctions2D(knotBegin, knotIncr, param, N);

          //printf(" \t %14.8f \t %14.8f \t %14.8f \t %14.8f\n", uu, vv, fact, dvol0);
          //printf("BasisFuns \n");
          //for(ii=0;ii<totnlbf;ii++)
            //printf(" \t %12.6f  \t %12.6f  \t %12.6f \n ", N(ii), dN_dx(ii), dN_dy(ii));

          xx = GeomData->computeCoord(0, param[0]) - position[0];
          yy = GeomData->computeCoord(1, param[1]) - position[1];

          delta1 = DiracDelta1(xx, beta1);
          delta2 = DiracDelta1(yy, beta2);

          delta = delta1 * delta2;

          dvol *= (delta*fact);

          for(ii=0;ii<totnlbf;ii++)
          {
            TI = ndof*ii;

            for(jj=0;jj<ndof;jj++)
              D(TI+jj,jj) = N[ii];
          }

          //if(delta > 0.0)
            //printf(" \t %14.8f\t %14.8f\t%14.8f\t %14.8f\t%14.8f\n", xx, yy, delta1, delta2, delta);

          //val1 = LevelSetFunc[0][7]*normal[0]*fact;
          //val2 = LevelSetFunc[0][7]*normal[1]*fact;

          //Flocal += (dvol*(D*res));

    }//gp1
    }//gp2
    
    //printVector(Flocal);
    //printf("\n\n");

  return;
}



