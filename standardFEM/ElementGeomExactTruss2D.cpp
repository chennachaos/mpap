
#include "ElementGeomExactTruss2D.h"

#include "Debug.h"
#include "MpapTime.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"


using namespace std;

extern MpapTime mpapTime;


ElementGeomExactTruss2D::ElementGeomExactTruss2D(void)
{
  ndof = 2;
  ndim = 2;
  npElem = 2;
  nlbf  = npElem;
  nsize = npElem*ndof;
  
  if (debug) cout << " constructor ElementGeomExactTruss2D\n\n";
}



ElementGeomExactTruss2D::~ElementGeomExactTruss2D()
{
  if (debug) cout << " destructor ElementGeomExactTruss2D\n\n";
}



void ElementGeomExactTruss2D::prepareElemData()
{
  LagrangeElement::prepareElemData();

  int ii, jj, ind, kk;
  
  //Klocal.resize(nsize, nsize);
  //Flocal.resize(nsize);
  
  return;
}




void ElementGeomExactTruss2D::prepareElemData2()
{
  int ii, jj, ind, kk;
  
  if(!(SolnData->STAGGERED) )
  {
    //cout << " aaaaaaaaaaaaaa " << npElem*ndim << endl;
    forAssyVec2.resize(npElem*ndim);
    
    kk = 0;
    for(ii=0;ii<npElem;ii++)
    {
      ind = nodeNums[ii]*ndim;
      for(jj=0;jj<ndim;jj++)
      {
        //cout << ii << '\t' << jj << '\t' << kk << '\t' << ind << endl;
        forAssyVec2[kk++] = ind + jj;
      }
    }
  }
  //cout << " PPPPPPPPPPP " << endl;
  
  return;
}


int ElementGeomExactTruss2D::calcLoadVector()
{
  return 0;
}


int ElementGeomExactTruss2D::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
  // From NFEM course material
  // http://www.colorado.edu/engineering/cas/courses.d/NFEM.d/Home.html

    int  ii, jj, count;

    double  uxn[2], uzn[2], btn[2], b[2], af, d1, xx[2];

    double  aa, rho0, A0, A, E, F, strain, stress, stress0, fact, fact1, fact2;
    double  x0[2], y0[2], x1[2], y1[2], trac[2];
    double  h, h0, c0x, c0y, cx, cy, dx, dy, py;

    VectorXd  temp(2), accC(4), forceC(4), velC(4), dispC(4), B0(4), B(4);
    MatrixXd  Mlocal(4,4), H(4,4);
    H.setZero();
    
    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    
    finite = true;
    followerLoadFlag = false;

    finite = int(elmDat[0] == 1) ;
    b[0] = elmDat[1];
    b[1] = elmDat[2];
    rho0 = elmDat[3];
    A0   = elmDat[4];
    E    = elmDat[5];

    //dt = SolidSolnData->td(0);
    af = SolnData->td(2);
    d1 = SolnData->td(5);
    aa = SolnData->td(10);

    //  rotate nodal displacements and compute nodal positions on element axis
    
    //cout << " ccccccccccc " << endl;
    
    x0[0] = GeomData->NodePosOrig[nodeNums[0]][0];
    y0[0] = GeomData->NodePosOrig[nodeNums[0]][1];
    x0[1] = GeomData->NodePosOrig[nodeNums[1]][0];
    y0[1] = GeomData->NodePosOrig[nodeNums[1]][1];

    dispC(0) = SolnData->var1Cur[nodeNums[0]*ndof+0];
    dispC(1) = SolnData->var1Cur[nodeNums[0]*ndof+1];
    dispC(2) = SolnData->var1Cur[nodeNums[1]*ndof+0];
    dispC(3) = SolnData->var1Cur[nodeNums[1]*ndof+1];

    x1[0] = x0[0];    y1[0] = y0[0];    x1[1] = x0[1];    y1[1] = y0[1];
    
    x1[0] += dispC(0);
    y1[0] += dispC(1);
    x1[1] += dispC(2);
    y1[1] += dispC(3);


    // compute the orientation of the element
    
    dx = x0[1] - x0[0];
    dy = y0[1] - y0[0];
    
    h0 = sqrt(dx*dx+dy*dy);

    c0x = dx/h0;
    c0y = dy/h0;

    dx = x1[1] - x1[0];
    dy = y1[1] - y1[0];
   
    h = sqrt(dx*dx+dy*dy);

    cx = dx/h0;
    cy = dy/h0;
    
    //cout << c0x << '\t' << c0y << '\t' << cx << '\t' << cy << endl;

    B0(0) = -c0x;   B0(1) = -c0y;   B0(2) =  c0x;   B0(3) =  c0y;
    B0 /= h0;
    
    B = B0;
    
    if(finite)
    {
      B(0)  = -cx;    B(1)  = -cy;    B(2)  =  cx;    B(3)  =  cy;
      B  /= h0;

      H(0,0) =  1.0;    H(0,2) = -1.0;
      H(1,1) =  1.0;    H(1,3) = -1.0;
      H(2,0) = -1.0;    H(2,2) =  1.0;
      H(3,1) = -1.0;    H(3,3) =  1.0;
      
      H /= (h0*h0);
    }
    //printVector(B0);
    //printVector(B);

    //cout << rho0 << '\t' << A0 << '\t' << E << endl;
    //cout << x1 << '\t' << y1 << '\t' << x2 << '\t' << y2 << '\t' << lth0 << endl;
    //cout << " element = " << nodeNums[0] << '\t' << nodeNums[1] << endl;
    
    Klocal.setZero();
    Mlocal.setZero();
    Flocal.setZero();


    stress0 = 0.0;
    //fact = h/h0;
    //strain = 0.5*(fact*fact - 1.0);

    strain = B0.transpose()*dispC;
    
    if(finite)
      strain += (0.5*dispC.transpose()*(H*dispC));
    
    stress = stress0 + E * strain;
    F = A0*stress;
    //cout << " strain and stress " << strain << '\t' << stress << endl;

    // residual
    
    Flocal -= (F*h0)*B;

    // material stiffness
    
    fact = E*A0*h0;
    Klocal += ((fact*B)*B.transpose()) ;

    // geometric stiffness
    
    if(finite)
    {
      Klocal += ((F*h0)*H);
    }

    Klocal *= af;
    
    
    ////////////////////////////////
    //
    // pressure loading
    //
    ////////////////////////////////
    
    py = 0.0;
    fact = 0.5*py;
    Flocal(0) += fact*(y1[1]-y1[0]);
    Flocal(1) += fact*(x1[1]-x1[0]);
    Flocal(2) += fact*(y1[1]-y1[0]);
    Flocal(3) += fact*(x1[1]-x1[0]);
    
    if(followerLoadFlag)
    {
       fact *= af;
       Klocal(0,1) += -fact;
       Klocal(0,3) +=  fact;
       Klocal(1,0) +=  fact;
       Klocal(1,2) += -fact;
       Klocal(2,1) += -fact;
       Klocal(2,3) +=  fact;
       Klocal(3,0) +=  fact;
       Klocal(3,2) += -fact;
    }
    
    
    //body forces
    
    // need to be implemented

    //inertia
    
    fact1 = rho0*A0*h0/6.0;
    fact2 = 2.0*fact1;

    Mlocal(0,0) = fact2;    Mlocal(0,2) = fact1;
    Mlocal(1,1) = fact2;    Mlocal(1,3) = fact1;
    Mlocal(2,0) = fact1;    Mlocal(2,2) = fact2;
    Mlocal(3,1) = fact1;    Mlocal(3,3) = fact2;

    accC(0) =  SolnData->var1DotDotCur(nodeNums[0]*ndof+0);
    accC(1) =  SolnData->var1DotDotCur(nodeNums[0]*ndof+1);
    accC(2) =  SolnData->var1DotDotCur(nodeNums[1]*ndof+0);
    accC(3) =  SolnData->var1DotDotCur(nodeNums[1]*ndof+1);

    Klocal += (d1*Mlocal);
    Flocal -= (Mlocal*accC);

    //Klocal /= aa;

    //printMatrix(Klocal);	printf("\n\n\n");
    //printMatrix(Mlocal);     printf("\n\n\n");
    //printVector(Flocal);	printf("\n\n\n");


    forceC(0) =  SolnData->forceCur(nodeNums[0]*ndof+0);
    forceC(1) =  SolnData->forceCur(nodeNums[0]*ndof+1);
    forceC(2) =  SolnData->forceCur(nodeNums[1]*ndof+0);
    forceC(3) =  SolnData->forceCur(nodeNums[1]*ndof+1);

    //printMatrix(Klocal); printf("\n\n");

   return 0;
}


int ElementGeomExactTruss2D::calcInternalForces()
{
  return 0;
}



void ElementGeomExactTruss2D::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void ElementGeomExactTruss2D::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void ElementGeomExactTruss2D::projectStress(int varindex, double* outval)
{
  return;
}

void ElementGeomExactTruss2D::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void ElementGeomExactTruss2D::projectIntVar(int index, double* outval)
{
   return;
}


int ElementGeomExactTruss2D::calcOutput(double u1, double v1)
{
  return 0;
}



void ElementGeomExactTruss2D::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}


