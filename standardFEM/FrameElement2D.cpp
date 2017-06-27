
#include "FrameElement2D.h"

#include "Debug.h"
#include "MpapTime.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "FunctionsMaterial.h"


using namespace std;

extern MpapTime mpapTime;


FrameElement2D::FrameElement2D(void)
{
  ndof = 3;
  ndim = 2;
  npElem = 2;
  nlbf = npElem;
  nsize = npElem*ndof;
  
  if (debug) cout << " constructor FrameElement2D\n\n";
}



FrameElement2D::~FrameElement2D()
{
  if (debug) cout << " destructor FrameElement2D\n\n";
}



void FrameElement2D::prepareElemData()
{
  LagrangeElement::prepareElemData();
  
  return;
}


void FrameElement2D::prepareElemData2()
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


int FrameElement2D::calcLoadVector()
{
  return 0;
}


//
int FrameElement2D::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
    int  ii, jj, gp, kk, ll, TI, TIp1, TIp2, count, TJ, TJp1, TJp2;
    
    double  lth0, sth0, cth0, uxn[2], uzn[2], b[2], af, am, dt, d1, rho;
    
    double  ux, uz, bt, sbt, cbt, dux, duz, dbt, detJ, r1d6=1.0/6.0;
    double  A, E, G, I, nu, EA, GA, EI, kappa, dvol, fact, fact1, fact2, fact3, fact4, EAdv, GAdv, EIdv;
    double  x2, y2, x1, y1, u1, u2, w1, w2, btn[2], ddu[2], ddw[2], ddbt[2], h;
    
    nlbf = npElem;
    nsize  = npElem*ndof;
    vector<double>  N(nlbf), dN_dx(nlbf);

    VectorXd  res(6), temp(2);
    MatrixXd  Mlocal(6,6), K2(6,6), Bi(3,nsize), Bj(3,3), D(3,3), hlp(2,2), R(2,2), Rt(2,2);
    D.setZero();
    Mlocal.setZero();

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    matDat = &(SolnData->MatlProp[matType].data[0]);

    nGP  = 1;
    b[0] = elmDat[0];
    b[1] = elmDat[1];
    rho  = elmDat[2];
    A    = elmDat[3];
    I    = elmDat[4];
    E    = elmDat[5];
    nu   = elmDat[6];
    kappa= elmDat[7];
    
    G  = E/2.0/(1.0+nu);
    EA = E*A;
    EI = E*I;
    GA = G*A*kappa;
    
    //cout << " material constants... "  << EA << '\t' << EI << '\t' << GA << endl;
    
    dt = SolnData->td(0);
    af = SolnData->td(2);
    am = SolnData->td(1);
    d1 = SolnData->td(5);
    
    //af = 1.0;
    //d1 = 0.0;
    
    //  rotate nodal displacements and compute nodal positions on element axis
    
    //cout << " ccccccccccc " << endl;
    
    x1 = GeomData->NodePosOrig[nodeNums[0]][0];
    y1 = GeomData->NodePosOrig[nodeNums[0]][1];
    x2 = GeomData->NodePosOrig[nodeNums[1]][0];
    y2 = GeomData->NodePosOrig[nodeNums[1]][1];


    res(0) = SolnData->var1Cur[nodeNums[0]*ndof+0];
    res(1) = SolnData->var1Cur[nodeNums[0]*ndof+1];
    res(2) = SolnData->var1Cur[nodeNums[0]*ndof+2];
    res(3) = SolnData->var1Cur[nodeNums[1]*ndof+0];
    res(4) = SolnData->var1Cur[nodeNums[1]*ndof+1];
    res(5) = SolnData->var1Cur[nodeNums[1]*ndof+2];

    // compute the orientation of the element
    
    lth0 = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    h = lth0;
    sth0   = (y2-y1)/lth0;
    cth0   = (x2-x1)/lth0;
    
    R(0,0) =  cth0;    R(0,1) =  sth0;
    R(1,0) = -sth0;    R(1,1) =  cth0;
    
    Rt = R.transpose();
    
    uxn[0] = + u1*cth0 + w1*sth0;
    uzn[0] = - u1*sth0 + w1*cth0;
    uxn[1] = + u2*cth0 + w2*sth0;
    uzn[1] = - u2*sth0 + w2*cth0;

    //loop over Gauss points (length)

    //cout << x0 << '\t' << y0 << '\t' << x1 << '\t' << y1 << endl;
    //cout << sth0 << '\t' << cth0 << endl;
    //cout << " ccccccccccc " << endl;
    
    Klocal.setZero();
    Flocal.setZero();
    
    fact = EA/h;

    Klocal(0,0) =  fact;
    Klocal(0,3) = -fact;
    Klocal(3,0) = -fact;
    Klocal(3,3) =  fact;

    fact = EI/h;

    Klocal(2,2) =  fact;
    Klocal(2,5) = -fact;
    Klocal(5,2) = -fact;
    Klocal(5,5) =  fact;

    K2.setZero();

    K2(1,1) =  4.0;
    K2(1,2) = -2.0*h;
    K2(1,4) = -4.0;
    K2(1,5) = -2.0*h;

    K2(2,1) = -2.0*h;
    K2(2,2) =  h*h;
    K2(2,4) =  2.0*h;
    K2(2,5) =  h*h;

    K2(4,1) = -4.0;
    K2(4,2) =  2.0*h;
    K2(4,4) =  4.0;
    K2(4,5) =  2.0*h;

    K2(5,1) = -2.0*h;
    K2(5,2) =  h*h;
    K2(5,4) =  2.0*h;
    K2(5,5) =  h*h;

    fact = GA/h/4.0;

    Klocal += (fact*K2);

    Flocal -= (Klocal*res);

    Klocal *= af;

    //printMatrix(Klocal);	printf("\n\n\n");

    //cout << sth0 << '\t' << cth0 << endl;
    //rotate residual and stiffness matrix

      /*
      for(ii=0;ii<npElem;ii++)
      {
        TI   =  3*ii;
        TIp1 =  TI+1;
        TIp2 =  TI+2;
        
        hlp(1,1) = Flocal(TI);
        Flocal(TI)   = cth0 * Flocal(TI)  - sth0 * Flocal(TIp1);
        Flocal(TIp1) = sth0 * hlp(1,1)    + cth0 * Flocal(TIp1);

        for(jj=0;jj<npElem;jj++)
        {
          TJ   = 3*jj;
          TJp1 = TJ+1;
          TJp2 = TJ+2;

          hlp(0,0) = cth0 * Klocal(TI,TJ)   - sth0 * Klocal(TIp1,TJ);
          hlp(1,0) = sth0 * Klocal(TI,TJ)   + cth0 * Klocal(TIp1,TJ);
          hlp(0,1) = cth0 * Klocal(TI,TJp1) - sth0 * Klocal(TIp1,TJp1);
          hlp(1,1) = sth0 * Klocal(TI,TJp1) + cth0 * Klocal(TIp1,TJp1);

          Klocal(TI,TJ)     = hlp(0,0) * cth0 - hlp(0,1) * sth0;
          Klocal(TI,TJp1)   = hlp(0,0) * sth0 + hlp(0,1) * cth0;
          Klocal(TIp1,TJ)   = hlp(1,0) * cth0 - hlp(1,1) * sth0;
          Klocal(TIp1,TJp1) = hlp(1,0) * sth0 + hlp(1,1) * cth0;

          hlp(1,1) = Klocal(TIp2,TJ);
          Klocal(TIp2,TJ) = Klocal(TIp2,TJ) * cth0 - Klocal(TIp2,TJp1) * sth0;
          Klocal(TIp2,TJp1) = hlp(1,1)     * sth0 + Klocal(TIp2,TJp1) * cth0;

          hlp(1,1) = Klocal(TI,TJp2);
          Klocal(TI,TJp2) = cth0 * Klocal(TI,TJp2) - sth0 * Klocal(TIp1,TJp2);
          Klocal(TIp1,TJp2) = sth0 * hlp(1,1)     + cth0 * Klocal(TIp1,TJp2);
        }
      }
      */
    
    //printMatrix(Klocal);	printf("\n\n\n");

    //body forces

    /*
    fact = rho * lth0 * 0.5;
      
    b[0] *= fact;
    b[1] *= fact;

    Flocal(0) += b[0];
    Flocal(1) += b[1];

    Flocal(3) += b[0];
    Flocal(4) += b[1];
    */

    //inertia
    fact1 = rho*A*h/6.0;
    fact2 = rho*I*h/6.0;

    fact3 = 2.0*fact1;
    fact4 = 2.0*fact2;

    Mlocal(0,0) = fact3;
    Mlocal(1,1) = fact3;
    Mlocal(3,3) = fact3;
    Mlocal(4,4) = fact3;

    Mlocal(0,3) = fact1;
    Mlocal(3,0) = fact1;
    Mlocal(1,4) = fact1;
    Mlocal(4,1) = fact1;
    
    Mlocal(2,2) = fact4;
    Mlocal(5,5) = fact4;
    Mlocal(2,5) = fact2;
    Mlocal(5,2) = fact2;

    res(0) =  SolnData->var1DotDotCur(nodeNums[0]*3+0);
    res(1) =  SolnData->var1DotDotCur(nodeNums[0]*3+1);
    res(2) =  SolnData->var1DotDotCur(nodeNums[0]*3+2);
    res(3) =  SolnData->var1DotDotCur(nodeNums[1]*3+0);
    res(4) =  SolnData->var1DotDotCur(nodeNums[1]*3+1);
    res(5) =  SolnData->var1DotDotCur(nodeNums[1]*3+2);

    //printMatrix(Klocal); printf("\n\n");
    //printMatrix(Mlocal); printf("\n\n");
    
    Klocal += (d1*Mlocal);
    Flocal -= (Mlocal*res);
    
    //cout << " d1 = " << d1 << endl;
    //printMatrix(Klocal); printf("\n\n");
    //printVector(Flocal); printf("\n\n");

   return 0;
}
//




int FrameElement2D::calcInternalForces()
{
  return 0;
}



void FrameElement2D::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void FrameElement2D::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void FrameElement2D::projectStress(int varindex, double* outval)
{
  return;
}

void FrameElement2D::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void FrameElement2D::projectIntVar(int index, double* outval)
{
   return;
}


int FrameElement2D::calcOutput(double u1, double v1)
{
  return 0;
}

void FrameElement2D::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}


