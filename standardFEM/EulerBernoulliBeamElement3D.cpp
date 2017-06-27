
#include "EulerBernoulliBeamElement3D.h"

#include "Debug.h"
#include "MpapTime.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "util.h"

using namespace std;

extern MpapTime mpapTime;


EulerBernoulliBeamElement3D::EulerBernoulliBeamElement3D()
{
  ndof   = 6;
  ndim   = 3;
  npElem = 2;
  nlbf   = npElem;
  nsize  = npElem*ndof;

  if (debug) cout << " constructor EulerBernoulliBeamElement3D\n\n";
}



EulerBernoulliBeamElement3D::~EulerBernoulliBeamElement3D()
{
  if (debug) cout << " destructor EulerBernoulliBeamElement3D\n\n";
}



void EulerBernoulliBeamElement3D::prepareElemData()
{
  LagrangeElement::prepareElemData();

  return;
}




void EulerBernoulliBeamElement3D::prepareElemData2()
{
  int ii, jj, ind, kk, ind2;
  
  if(!(SolnData->STAGGERED) )
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


int EulerBernoulliBeamElement3D::calcLoadVector()
{
  return 0;
}


int EulerBernoulliBeamElement3D::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
    int  ii, jj, gp, kk, ll, TI, TIp1, TIp2, count, TJ, TJp1, TJp2, ind1, ind2;
    
    double  af, am, dt, d1, rho, b[3], fact, aa;
    double  A, E, G, Ix, Iy, Iz, EA, EIy, EIz, GIx;
    double  x1, y1, z1, x2, y2, z2, h;

    VectorXd  dispC(nsize), accC(nsize), velC(nsize);
    MatrixXd  Mlocal(nsize, nsize), RotMat(nsize, nsize), RotMatTrans(nsize, nsize);

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    matDat = &(SolnData->MatlProp[matType].data[0]);

    b[0] = elmDat[0];
    b[1] = elmDat[1];
    b[2] = elmDat[2];
    rho  = elmDat[3];
    A    = elmDat[4];
    Ix   = elmDat[5];
    Iy   = elmDat[6];
    Iz   = elmDat[7];
    E    = elmDat[8];
    G    = elmDat[9];

    EA  = E*A;
    EIy = E*Iy;
    EIz = E*Iz;
    GIx = G*Ix;

    //cout << " material constants... "  << EA << '\t' << EI << '\t' << GA << endl;
    
    dt = SolnData->td(0);
    af = SolnData->td(2);
    am = SolnData->td(1);
    d1 = SolnData->td(5);
    aa = SolnData->td(10);
    
    //  rotate nodal displacements and compute nodal positions on element axis
    
    //cout << " ccccccccccc " << endl;
    
    x1 = GeomData->NodePosOrig[nodeNums[0]][0];
    y1 = GeomData->NodePosOrig[nodeNums[0]][1];
    z1 = GeomData->NodePosOrig[nodeNums[0]][2];
    x2 = GeomData->NodePosOrig[nodeNums[1]][0];
    y2 = GeomData->NodePosOrig[nodeNums[1]][1];
    z2 = GeomData->NodePosOrig[nodeNums[1]][2];

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

    double  CXx, CXy, CXz, CYx, CYy, CYz, CZx, CZy, CZz, h2;
    // compute the orientation of the element

    h    = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1));
    
    if( CompareDoubles(x1, x2) && CompareDoubles(y1, y2) )
    {
       if(z2 < z1) // because here R^T is computed as RotMat
       {
         CXx =  0.0; CXy = 0.0; CXz = 1.0;
         CYx =  0.0; CYy = 1.0; CYz = 0.0;
         CZx = -1.0; CZy = 0.0; CZz = 0.0;
       }
       else
       {
         CXx = 0.0; CXy = 0.0; CXz = -1.0;
         CYx = 0.0; CYy = 1.0; CYz = 0.0;
         CZx = 1.0; CZy = 0.0; CZz = 0.0;
       }
    }
    else
    {
      CXx = (x2-x1)/h;
      CYx = (y2-y1)/h;
      CZx = (z2-z1)/h;

      h2 = sqrt(CXx*CXx + CYx*CYx);

      CXy = -CYx/h2;
      CYy =  CXx/h2;
      CZy =  0.0;

      CXz = -CXx*CZx/h2;
      CYz = -CYx*CZx/h2;
      CZz =  h2;
    }

    RotMat.setZero();

    //RotMat(0,0) = CXx; RotMat(0,1) = CXy; RotMat(0,2) = CXz;
    //RotMat(1,0) = CYx; RotMat(1,1) = CYy; RotMat(1,2) = CYz;
    //RotMat(2,0) = CZx; RotMat(2,1) = CZy; RotMat(2,2) = CZz;

    for(ii=0;ii<4;ii++)
    {
      jj = 3*ii;
      RotMat(jj+0,jj+0) = CXx; RotMat(jj+0,jj+1) = CXy; RotMat(jj+0,jj+2) = CXz;
      RotMat(jj+1,jj+0) = CYx; RotMat(jj+1,jj+1) = CYy; RotMat(jj+1,jj+2) = CYz;
      RotMat(jj+2,jj+0) = CZx; RotMat(jj+2,jj+1) = CZy; RotMat(jj+2,jj+2) = CZz;
    }


    RotMatTrans = RotMat.transpose();

    cout << x1 << '\t' << y1 << '\t' << x2 << '\t' << y2 <<  '\t' << z1 << '\t' << z2 <<endl;
    //cout << sth0 << '\t' << cth0 << endl;
    //cout << " ccccccccccc " << endl;

    Klocal.setZero();
    
    // axial component
    fact = EA/h;

    Klocal(0,0) =  fact;
    Klocal(0,6) = -fact;
    Klocal(6,0) = -fact;
    Klocal(6,6) =  fact;

    // torsion component
    fact = GIx/h;

    Klocal(3,3) =  fact;
    Klocal(3,9) = -fact;
    Klocal(9,3) = -fact;
    Klocal(9,9) =  fact;

    // bending component x-y plane
    fact = 12.0*EIz/h/h/h;

    Klocal(1,1) =  fact;
    Klocal(1,7) = -fact;
    Klocal(7,1) = -fact;
    Klocal(7,7) =  fact;

    fact = 6.0*EIz/h/h;

    Klocal(1,5)  = fact;
    Klocal(5,1)  = fact;
    Klocal(1,11) = fact;
    Klocal(11,1) = fact;

    Klocal(5,7)  = -fact;
    Klocal(7,5)  = -fact;
    Klocal(7,11) = -fact;
    Klocal(11,7) = -fact;

    fact = EIz/h;

    Klocal(5,5)   =  4.0*fact;
    Klocal(5,11)  =  2.0*fact;
    Klocal(11,5)  =  2.0*fact;
    Klocal(11,11) =  4.0*fact;

    // bending component x-z plane

    fact = 12.0*EIy/h/h/h;

    Klocal(2,2) =  fact;
    Klocal(2,8) = -fact;
    Klocal(8,2) = -fact;
    Klocal(8,8) =  fact;

    fact = 6.0*EIy/h/h;

    Klocal(2,4)  = -fact;
    Klocal(4,2)  = -fact;
    Klocal(2,10) = -fact;
    Klocal(10,2) = -fact;

    Klocal(4,8)  = fact;
    Klocal(8,4)  = fact;
    Klocal(8,10) = fact;
    Klocal(10,8) = fact;

    fact = EIy/h;

    Klocal(4,4)   =  4.0*fact;
    Klocal(4,10)  =  2.0*fact;
    Klocal(10,4)  =  2.0*fact;
    Klocal(10,10) =  4.0*fact;

    //printMatrix(Klocal);	printf("\n\n\n");

    //body forces

    //inertia

    Mlocal.setZero();

    // axial component
    Mlocal(0,0) = 140.0;
    Mlocal(0,6) =  70.0;
    Mlocal(6,0) =  70.0;
    Mlocal(6,6) = 140.0;

    // torsion component

    Mlocal(3,3) = 140.0;
    Mlocal(3,9) =  70.0;
    Mlocal(9,3) =  70.0;
    Mlocal(9,9) = 140.0;

    // bending component x-y plane

    Mlocal(1,1) = 156.0;
    Mlocal(1,7) =  54.0 ;
    Mlocal(7,1) =  54.0;
    Mlocal(7,7) = 156.0;

    fact = 22.0*h;

    Mlocal(1,5)  =  fact;
    Mlocal(5,1)  =  fact;
    Mlocal(7,11) = -fact;
    Mlocal(11,7) = -fact;

    fact = 13.0*h;

    Mlocal(5,7)  =  fact;
    Mlocal(7,5)  =  fact;
    Mlocal(1,11) = -fact;
    Mlocal(11,1) = -fact;

    fact = h*h;

    Mlocal(5,5)   =  4.0*fact;
    Mlocal(5,11)  = -3.0*fact;
    Mlocal(11,5)  = -3.0*fact;
    Mlocal(11,11) =  4.0*fact;

    // bending component x-z plane

    Mlocal(2,2) = 156.0;
    Mlocal(2,8) =  54.0 ;
    Mlocal(8,2) =  54.0;
    Mlocal(8,8) = 156.0;

    fact = 22.0*h;

    Mlocal(2,4)  =  fact;
    Mlocal(4,2)  =  fact;
    Mlocal(8,10) = -fact;
    Mlocal(10,8) = -fact;

    fact = 13.0*h;

    Mlocal(4,8)  =  fact;
    Mlocal(8,4)  =  fact;
    Mlocal(2,10) = -fact;
    Mlocal(10,2) = -fact;

    fact = h*h;

    Mlocal(4,4)   =  4.0*fact;
    Mlocal(4,10)  = -3.0*fact;
    Mlocal(10,4)  = -3.0*fact;
    Mlocal(10,10) =  4.0*fact;


    fact = rho*A*h/420.0;

    Mlocal = fact*Mlocal;

    printMatrix(Klocal); printf("\n\n");
    //printMatrix(Mlocal); printf("\n\n");

    dispC = RotMatTrans*dispC;
    accC  = RotMatTrans*accC;
    velC  = RotMatTrans*velC;

    Flocal = -Klocal*dispC - Mlocal*accC;
    Klocal = af*Klocal + d1*Mlocal;

    //Klocal /= aa;

    Klocal = (RotMat*Klocal)*RotMatTrans;
    Flocal = RotMat*Flocal;
    
    printMatrix(RotMatTrans);	printf("\n\n\n");

    printMatrix(Klocal); printf("\n\n");
    //printVector(Flocal); printf("\n\n");

   return 0;
}
//



int EulerBernoulliBeamElement3D::calcInternalForces()
{
  return 0;
}



void EulerBernoulliBeamElement3D::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void EulerBernoulliBeamElement3D::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void EulerBernoulliBeamElement3D::projectStress(int varindex, double* outval)
{
  return;
}


void EulerBernoulliBeamElement3D::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void EulerBernoulliBeamElement3D::projectIntVar(int index, double* outval)
{
   return;
}


int EulerBernoulliBeamElement3D::calcOutput(double u1, double v1)
{
  return 0;
}


void EulerBernoulliBeamElement3D::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}

