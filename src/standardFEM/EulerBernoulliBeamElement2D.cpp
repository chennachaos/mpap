
#include "EulerBernoulliBeamElement2D.h"

#include "Debug.h"
#include "MpapTime.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"

using namespace std;

extern MpapTime mpapTime;


EulerBernoulliBeamElement2D::EulerBernoulliBeamElement2D()
{
  ndof   = 3;
  ndim   = 2;
  npElem = 2;
  nlbf = npElem;
  nsize  = npElem*ndof;

  if (debug) cout << " constructor EulerBernoulliBeamElement2D\n\n";
}



EulerBernoulliBeamElement2D::~EulerBernoulliBeamElement2D()
{
  if (debug) cout << " destructor EulerBernoulliBeamElement2D\n\n";
}



void EulerBernoulliBeamElement2D::prepareElemData()
{
  LagrangeElement::prepareElemData();

  //Klocal.resize(nsize, nsize);
  //Flocal.resize(nsize);
  //Flocal2.resize(4);

  return;
}




void EulerBernoulliBeamElement2D::prepareElemData2()
{
  int ii, jj, ind, kk, ind2;
  
  if(! (SolnData->STAGGERED) )
  {
    //cout << " aaaaaaaaaaaaaa " << npElem*ndim << endl;

    ind = npElem*ndim;

    forAssyVec2.resize(ind);
    kk = 0;
    for(ii=0;ii<npElem;ii++)
    {
      ind2 = nodeNums[ii]*ndim;
      for(jj=0;jj<ndim;jj++)
      {
        //cout << ii << '\t' << jj << '\t' << kk << '\t' << ind2 << endl;
        forAssyVec2[kk++] = ind2 + jj;
      }
    }
  }
  
  return;
}



int EulerBernoulliBeamElement2D::calcLoadVector()
{
  return 0;
}


/*
int EulerBernoulliBeamElement2D::calcStiffnessAndResidual()
{
    int  ii, jj, gp, kk, ind1, ind2;
    
    double  sth0, cth0, af, am, dt, d1, rho, b[2];
    double  A, E, I, aa, EA, EI, fact, fact1, fact2, fact3, fact4;
    double  x2, y2, x1, y1, u1, u2, w1, w2, h;

    VectorXd  dispC(6), accC(6), velC(6);
    MatrixXd  Mlocal(6,6), RotMat(6,6), RotMatTrans(6,6);

    elmDat = &(SolnData->ElemProp.data[0]);
    matDat = &(SolnData->MatlProp.data[0]);

    nGP1 = 1;
    b[0] = elmDat[0];
    b[1] = elmDat[1];
    rho  = elmDat[2];
    A    = elmDat[3];
    I    = elmDat[4];
    E    = elmDat[5];

    EA = E*A;
    EI = E*I;
    
    //cout << " material constants... "  << EA << '\t' << EI << '\t' << GA << endl;
    
    dt = SolnData->td(0);
    am = SolnData->td(1);
    af = SolnData->td(2);
    d1 = SolnData->td(5);
    aa = SolnData->td(10);

    //  rotate nodal displacements and compute nodal positions on element axis

    //cout << " ccccccccccc " << endl;

    x1 = GeomData->NodePosOrig[nodeNums[0]][0];
    y1 = GeomData->NodePosOrig[nodeNums[0]][1];
    x2 = GeomData->NodePosOrig[nodeNums[1]][0];
    y2 = GeomData->NodePosOrig[nodeNums[1]][1];

    for(ii=0;ii<npElem;ii++)
    {
      ind1 = ndof*ii;
      ind2 = nodeNums[ii]*ndof;

      for(kk=0;kk<ndof;kk++)
      {
        dispC(ind1+kk)  =  SolnData->var1Cur[ind2+kk];
        velC(ind1+kk)   =  SolnData->var1DotCur[ind2+kk];
        accC(ind1+kk)   =  SolnData->var1DotDotCur[ind2+kk];
      }
    }

    // compute the orientation of the element

    h    = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    sth0 = (y2-y1)/h;
    cth0 = (x2-x1)/h;

    RotMat.setZero();
    RotMat(0,0) = cth0; RotMat(0,1) = -sth0;
    RotMat(1,0) = sth0; RotMat(1,1) =  cth0;
    RotMat(2,2) = 1.0;

    RotMat(3,3) = cth0; RotMat(3,4) = -sth0;
    RotMat(4,3) = sth0; RotMat(4,4) =  cth0;
    RotMat(5,5) = 1.0;
    
    RotMatTrans = RotMat.transpose();

    //cout << x1 << '\t' << y1 << '\t' << x2 << '\t' << y2 << endl;
    //cout << sth0 << '\t' << cth0 << endl;
    //cout << " ccccccccccc " << d1 << endl;

    Klocal.setZero();
    
    fact = EA/h;

    Klocal(0,0) =  fact;
    Klocal(0,3) = -fact;
    Klocal(3,0) = -fact;
    Klocal(3,3) =  fact;

    fact = 12.0*EI/h/h/h;

    Klocal(1,1) =  fact;
    Klocal(1,4) = -fact;
    Klocal(4,1) = -fact;
    Klocal(4,4) =  fact;

    fact = 6.0*EI/h/h;

    Klocal(1,2) =  fact;
    Klocal(2,1) =  fact;
    Klocal(1,5) =  fact;
    Klocal(5,1) =  fact;

    Klocal(2,4) = -fact;
    Klocal(4,2) = -fact;
    Klocal(4,5) = -fact;
    Klocal(5,4) = -fact;

    fact = EI/h;

    Klocal(2,2) =  4.0*fact;
    Klocal(2,5) =  2.0*fact;
    Klocal(5,2) =  2.0*fact;
    Klocal(5,5) =  4.0*fact;

    //printMatrix(Klocal);	printf("\n\n\n");

    //body forces

    //inertia

    Mlocal.setZero();

    Mlocal(0,0) = 140.0;
    Mlocal(0,3) =  70.0;

    Mlocal(1,1) = 156.0;
    Mlocal(1,2) =  22.0*h;
    Mlocal(1,4) =  54.0;
    Mlocal(1,5) = -13.0*h;

    Mlocal(2,1) =  22.0*h;
    Mlocal(2,2) =   4.0*h*h;
    Mlocal(2,4) =  13.0*h;
    Mlocal(2,5) =  -3.0*h*h;

    Mlocal(3,0) = 70.0;
    Mlocal(3,3) = 140.0;

    Mlocal(4,1) = 54.0;
    Mlocal(4,2) = 13.0*h;
    Mlocal(4,4) = 156.0;
    Mlocal(4,5) = -22.0*h;

    Mlocal(5,1) = -13.0*h;
    Mlocal(5,2) = -3.0*h*h;
    Mlocal(5,4) = -22.0*h;
    Mlocal(5,5) = 4.0*h*h;

    fact = rho*A*h/420.0;

    Mlocal = fact*Mlocal;

    //printVector(dispC);
    //printMatrix(Klocal); printf("\n\n");
    //printMatrix(Mlocal); printf("\n\n");

    dispC = RotMatTrans*dispC;
    accC  = RotMatTrans*accC;
    velC  = RotMatTrans*velC;

    Flocal = -Klocal*dispC - Mlocal*accC;
    Klocal = af*Klocal + d1*Mlocal;

    //Klocal /= aa;

    Klocal = (RotMat*Klocal)*RotMatTrans;
    Flocal = RotMat*Flocal;

    //printMatrix(Klocal); printf("\n\n");
    //printVector(Flocal); printf("\n\n");

   return 0;
}
*/



//
int EulerBernoulliBeamElement2D::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
    int  ii, jj, gp, kk, ind1, ind2, TI;
    
    double  sth0, cth0, af, am, d1, rho, b[2];
    double  A, E, I, EA, EI, fact;
    double  x2, y2, x1, y1, h;

    VectorXd  vecF(6), qVec(12), qdotVec(12);
    MatrixXd  RotMat(12,12), RotMatTrans(12,12);
    MatrixXd  matM(6, 6), matC(6, 6), matK(6, 6);

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    matDat = &(SolnData->MatlProp[matType].data[0]);

    nGP  = 1;
    b[0] = elmDat[0];
    b[1] = elmDat[1];
    rho  = elmDat[2];
    A    = elmDat[3];
    I    = elmDat[4];
    E    = elmDat[5];

    EA = E*A;
    EI = E*I;
    
    //cout << " material constants... "  << EA << '\t' << EI << '\t' << GA << endl;

    am = SolnData->td(1);
    af = SolnData->td(2);
    d1 = SolnData->td(8);
    
    //cout << d1 << '\t' << af << endl;

    //  rotate nodal displacements and compute nodal positions on element axis
    
    //cout << " ccccccccccc " << endl;
    
    x1 = GeomData->NodePosOrig[nodeNums[0]][0];
    y1 = GeomData->NodePosOrig[nodeNums[0]][1];
    x2 = GeomData->NodePosOrig[nodeNums[1]][0];
    y2 = GeomData->NodePosOrig[nodeNums[1]][1];

    for(ii=0;ii<npElem;ii++)
    {
      ind1 = ndof*ii;
      ind2 = nodeNums[ii]*ndof;

      for(kk=0;kk<ndof;kk++)
      {
        qVec(ind1+kk)  =  SolnData->var1Cur[ind2+kk];
        qdotVec(ind1+kk)   =  SolnData->var1DotCur[ind2+kk];
      }
    }

    // compute the orientation of the element

    h    = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    sth0 = (y2-y1)/h;
    cth0 = (x2-x1)/h;
    
    RotMat.setZero();

    RotMat(0,0) = cth0; RotMat(0,1) = -sth0;
    RotMat(1,0) = sth0; RotMat(1,1) =  cth0;
    RotMat(2,2) = 1.0;

    RotMat(3,3) = cth0; RotMat(3,4) = -sth0;
    RotMat(4,3) = sth0; RotMat(4,4) =  cth0;
    RotMat(5,5) = 1.0;

    RotMat(6,6) = cth0; RotMat(6,7) = -sth0;
    RotMat(7,6) = sth0; RotMat(7,7) =  cth0;
    RotMat(8,8) = 1.0;

    RotMat(9,9) = cth0; RotMat(9,10) = -sth0;
    RotMat(10,9) = sth0; RotMat(10,10) =  cth0;
    RotMat(11,11) = 1.0;
    
    RotMatTrans = RotMat.transpose();

    //cout << x1 << '\t' << y1 << '\t' << x2 << '\t' << y2 << endl;
    //cout << sth0 << '\t' << cth0 << endl;
    //cout << " ccccccccccc " << d1 << endl;

    matK.setZero();
    vecF.setZero();
    
    fact = EA/h;

    matK(0,0) =  fact;
    matK(0,3) = -fact;
    matK(3,0) = -fact;
    matK(3,3) =  fact;

    fact = 12.0*EI/h/h/h;

    matK(1,1) =  fact;
    matK(1,4) = -fact;
    matK(4,1) = -fact;
    matK(4,4) =  fact;

    fact = 6.0*EI/h/h;

    matK(1,2) =  fact;
    matK(2,1) =  fact;
    matK(1,5) =  fact;
    matK(5,1) =  fact;

    matK(2,4) = -fact;
    matK(4,2) = -fact;
    matK(4,5) = -fact;
    matK(5,4) = -fact;

    fact = EI/h;

    matK(2,2) =  4.0*fact;
    matK(2,5) =  2.0*fact;
    matK(5,2) =  2.0*fact;
    matK(5,5) =  4.0*fact;


    //printMatrix(matK);	printf("\n\n\n");

    //body forces

    //inertia

    matM.setZero();
    matC.setZero();

    matM(0,0) = 140.0;
    matM(0,3) =  70.0;

    matM(1,1) = 156.0;
    matM(1,2) =  22.0*h;
    matM(1,4) =  54.0;
    matM(1,5) = -13.0*h;

    matM(2,1) =  22.0*h;
    matM(2,2) =   4.0*h*h;
    matM(2,4) =  13.0*h;
    matM(2,5) =  -3.0*h*h;

    matM(3,0) = 70.0;
    matM(3,3) = 140.0;

    matM(4,1) = 54.0;
    matM(4,2) = 13.0*h;
    matM(4,4) = 156.0;
    matM(4,5) = -22.0*h;

    matM(5,1) = -13.0*h;
    matM(5,2) = -3.0*h*h;
    matM(5,4) = -22.0*h;
    matM(5,5) = 4.0*h*h;

    fact = rho*A*h/420.0;

    matM = fact*matM;
    
    //printMatrix(matM);	printf("\n\n\n");

    qVec = RotMatTrans*qVec;
    qdotVec  = RotMatTrans*qdotVec;


    // create local stiffness matrix in state space form

    vector<int> perm1(6, 1), perm2(6, 1);

    perm1[0] = 0;  perm1[1] = 1;  perm1[2] = 2;  perm1[3] = 6;  perm1[4] = 7;  perm1[5] = 8;
    perm2[0] = 3;  perm2[1] = 4;  perm2[2] = 5;  perm2[3] = 9;  perm2[4] =10;  perm2[5] =11;

    vecF.setZero();
    for(ii=0; ii<6; ii++)
    {
      for(jj=0; jj<6; jj++)
        vecF(ii) -= matK(ii,jj)*qVec(perm1[jj]);
    }

    Klocal.setZero();
    Flocal.setZero();

    /*
    // contribution due stiffness and stress residual
    for(ii=0; ii<6; ii++)
    {
      TI = perm2[ii];

      Flocal(TI) += vecF(ii);

      for(jj=0; jj<6; jj++)
      {
        Klocal(TI, perm1[jj])  += ( af*matK(ii,jj) );
      }
    }

    // contribution from the qdot part
    A1.setZero();
    for(ii=0; ii<6; ii++)
    {
      for(jj=0; jj<6; jj++)
      {
        A1(perm1[ii], perm1[jj]) += matM(ii,jj) ;
        A1(perm2[ii], perm2[jj]) += matM(ii,jj) ;
      }
    }

    Klocal  +=  d1*A1;
    Flocal  -=  A1*qdotVec;

    // contribution from the q part
    // without stiffness contribution

    A1.setZero();
    for(ii=0; ii<6; ii++)
    {
      for(jj=0; jj<6; jj++)
      {
        A1(perm1[ii], perm2[jj])  -= matM(ii,jj) ;
        //A1(perm2[ii], perm1[jj])  += matK(ii,jj) ;
        A1(perm2[ii], perm2[jj])  += matC(ii,jj) ;
      }
    }

    Klocal += af*A1;
    Flocal -= A1*qVec;
    */

    double  fact1, fact2;

    for(ii=0; ii<6; ii++)
    {
      TI = perm2[ii];

      Flocal(TI) += vecF(ii);

      fact1=0.0;
      fact2=0.0;
      for(jj=0; jj<6; jj++)
      {
        Klocal(perm1[ii], perm1[jj])  += ( d1*matM(ii,jj) );
        Klocal(perm1[ii], perm2[jj])  -= ( af*matM(ii,jj) );

        fact1 += matM(ii,jj) * qdotVec(perm1[jj]);
        fact1 -= matM(ii,jj) * qVec(perm2[jj]);

        Klocal(perm2[ii], perm2[jj])  += ( d1*matM(ii,jj) );
        Klocal(perm2[ii], perm2[jj])  += ( af*matC(ii,jj) );
        Klocal(perm2[ii], perm1[jj])  += ( af*matK(ii,jj) );

        fact2 += matM(ii,jj) * qdotVec(perm2[jj]);
        fact2 += matC(ii,jj) * qVec(perm2[jj]);
      }
      Flocal(perm1[ii]) -= fact1;
      Flocal(perm2[ii]) -= fact2;
    }

    Klocal = (RotMat*Klocal)*RotMatTrans;

    Flocal = RotMat*Flocal;
    
    //printMatrix(Klocal);
    //printVector(Flocal);

   return 0;
}
//




int EulerBernoulliBeamElement2D::calcInternalForces()
{
  return 0;
}



void EulerBernoulliBeamElement2D::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void EulerBernoulliBeamElement2D::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void EulerBernoulliBeamElement2D::projectStress(int varindex, double* outval)
{
  return;
}


void EulerBernoulliBeamElement2D::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void EulerBernoulliBeamElement2D::projectIntVar(int index, double* outval)
{
   return;
}


int EulerBernoulliBeamElement2D::calcOutput(double u1, double v1)
{
  return 0;
}


void EulerBernoulliBeamElement2D::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}

