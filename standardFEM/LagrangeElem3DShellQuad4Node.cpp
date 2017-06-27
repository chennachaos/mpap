
#include "LagrangeElem3DShellQuad4Node.h"

#include "Debug.h"
#include "MpapTime.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"

using namespace std;

extern MpapTime mpapTime;


LagrangeElem3DShellQuad4Node::LagrangeElem3DShellQuad4Node()
{
  degree = 1;
  ndof   = 5;
  ndim   = 2;
  npElem = 4;
  nlbf = npElem;
  nsize  = npElem*ndof;

  if (debug) cout << " constructor LagrangeElem3DShellQuad4Node\n\n";
}



LagrangeElem3DShellQuad4Node::~LagrangeElem3DShellQuad4Node()
{
  if (debug) cout << " destructor LagrangeElem3DShellQuad4Node\n\n";
}



void LagrangeElem3DShellQuad4Node::prepareElemData()
{
  LagrangeElement::prepareElemData();

  return;
}



void LagrangeElem3DShellQuad4Node::prepareElemData2()
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



int LagrangeElem3DShellQuad4Node::calcLoadVector()
{
  return 0;
}




int LagrangeElem3DShellQuad4Node::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
{
  cout << " need to be corrected " << endl;

    int  ii, jj, gp1, gp2, kk, ind1, ind2, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    double  rho, fact, fact1, fact2, dvol, bforce, param[2];
    double  Jac, dt, bb1, bb2, bb3, cc1, cc2, cc3, JacMult;
    double  af, am, d1, acceFact, h, E, nu, G, I, kappa;
    double  exx, eyy, gxy, gxz, gyz, sxx, syy, txy, txz, tyz;

    VectorXd  N(nlbf), dN_dx(nlbf), dN_dy(nlbf);
    VectorXd  dispC(nsize), velC(nsize), accC(nsize);
    MatrixXd  Mlocal(nsize, nsize), cc(3,3), bc(3,3);

    //cout << " aaaaaaaaaa " << endl;

    elmDat = &(SolnData->ElemProp[elmType].data[0]);
    matDat = &(SolnData->MatlProp[matType].data[0]);

    bforce = elmDat[1];
    rho    = elmDat[2];
    h      = elmDat[3];
    E      = elmDat[4];
    nu     = elmDat[5];
    kappa  = elmDat[6];
    
    fact = E/(1.0-nu*nu);
    G    = E/2.0/(1.0+nu);
    G    = G*h*kappa;

    I = h*h*h/12.0;

    cc.setZero();
    cc(0,0) = 1.0;   cc(0,1) = nu;
    cc(1,0) = nu;    cc(1,1) = 1.0;
    cc(2,2) = (1.0-nu)*0.5;

    fact *= I;
    cc *= fact;

    bc.setZero();

    //printf(" %14.12f \t %14.12f \t %14.12f \n ", cc[0][0], cc[0][1], cc[0][2]);
    //printf(" %14.12f \t %14.12f \t %14.12f \n ", cc[1][0], cc[1][1], cc[1][2]);
    //printf(" %14.12f \t %14.12f \t %14.12f \n\n\n ", cc[2][0], cc[2][1], cc[2][2]);

    //cout << " material constants... "  << EA << '\t' << EI << '\t' << GA << endl;

    af = SolnData->td(2);
    d1 = SolnData->td(5);
    acceFact = SolnData->td(10);

    //cout << d1 << '\t' << af << endl;

    //cout << " ccccccccccc " << endl;

    double xNode[4], yNode[4], xx, yy;

    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }
    //cout << " ccccccccccc " << endl;
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

    //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", xNode[0], xNode[1], xNode[2], xNode[3]);
    //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n\n\n ", yNode[0], yNode[1], yNode[2], yNode[3]);

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();
    Mlocal.setZero();


    //   part 1. -- contribution due to bending

  int nGP = 2 ;

  vector<double>  gausspoints1, gaussweights1;

  getGaussPoints1D(nGP, gausspoints1, gaussweights1);


  for(gp2=0;gp2<nGP;gp2++)
  {
    JacMult = gaussweights1[gp2];

    param[1] = gausspoints1[gp2];

  for(gp1=0;gp1<nGP;gp1++)
  {
        param[0] = gausspoints1[gp1];

        GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          //xx = yy= 0.0;
          //for(ii=0;ii<nlbf;ii++)
          //{
            //cout << ii << '\t' << N[ii] << '\t' << dN_dx[ii] << '\t' << dN_dy[ii] << endl;
            //xx += N[ii]*xNode[ii];
            //yy += N[ii]*yNode[ii];
          //}

        dvol = gaussweights1[gp1] * JacMult * Jac;

        exx = computeValueCur(1, dN_dx);
        eyy = computeValueCur(2, dN_dy);
        gxy = computeValueCur(1, dN_dy) + computeValueCur(2, dN_dx) ;

        sxx = cc(0,0)*exx + cc(0,1)*eyy + cc(0,2)*gxy ;
        syy = cc(1,0)*exx + cc(1,1)*eyy + cc(1,2)*gxy ;
        txy = cc(2,0)*exx + cc(2,1)*eyy + cc(2,2)*gxy ;

        for(ii=0;ii<nlbf;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*N[ii];

          //bc[0][0] = 0.0;          bc[0][1] = 0.0;          bc[0][2] = 0.0;

          //bc[1][0] = bb1 * cc[0][0] + bb2 * cc[2][0];
          //bc[1][1] = bb1 * cc[0][1] + bb2 * cc[2][1];
          //bc[1][2] = bb1 * cc[0][2] + bb2 * cc[2][2];

          //bc[2][0] = bb2 * cc[1][0] + bb1 * cc[2][0];
          //bc[2][1] = bb2 * cc[1][1] + bb1 * cc[2][1];
          //bc[2][2] = bb2 * cc[1][2] + bb1 * cc[2][2];
          
          bc(0,0) = 0.0;     bc(0,1) = 0.0;    bc(0,2) = 0.0;
          bc(1,0) = bb1;     bc(1,1) = 0.0;    bc(1,2) = bb2;
          bc(2,0) = 0.0;     bc(2,1) = bb2;    bc(2,2) = bb1;

          bc = bc*cc;

          TI   = 3*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;

          Flocal(TI)   += (bb3*bforce) ;
          Flocal(TIp1) += (bb3*bforce - bb1*sxx - bb2*txy) ;
          Flocal(TIp2) += (bb3*bforce - bb1*txy - bb2*syy) ;

          for(jj=0; jj<nlbf; jj++)
          {
              cc1 = af*dN_dx[jj];
              cc2 = af*dN_dy[jj];

              TJ   = 3*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;

              Klocal(TI, TJ)      +=  0.0;
              Klocal(TI, TJp1)    +=  (bc(0,0) * cc1 + bc(0,2) * cc2) ;
              Klocal(TI, TJp2)    +=  (bc(0,1) * cc2 + bc(0,2) * cc1) ;

              Klocal(TIp1, TJ)    +=  0.0;
              Klocal(TIp1, TJp1)  +=  (bc(1,0) * cc1 + bc(1,2) * cc2) ;
              Klocal(TIp1, TJp2)  +=  (bc(1,1) * cc2 + bc(1,2) * cc1) ;

              Klocal(TIp2, TJ)    +=  0.0;
              Klocal(TIp2, TJp1)  +=  (bc(2,0) * cc1 + bc(2,2) * cc2) ;
              Klocal(TIp2, TJp2)  +=  (bc(2,1) * cc2 + bc(2,2) * cc1) ;

              bb3 = bb3*rho;
              Mlocal(TI,   TJ)    +=  bb3*N[jj];
              Mlocal(TIp1, TJp1)  +=  bb3*N[jj];
              Mlocal(TIp2, TJp2)  +=  bb3*N[jj];
          }
        }
  }//gp1
  }//gp2

  //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);

  //   part 2. -- contribution due to shear

  nGP = 1;

  getGaussPoints1D(nGP, gausspoints1, gaussweights1);


  for(gp2=0;gp2<nGP;gp2++)
  {
    JacMult = gaussweights1[gp2];

    param[1] = gausspoints1[gp2];

  for(gp1=0;gp1<nGP;gp1++)
  {
        param[0] = gausspoints1[gp1];

        GeomData->computeBasisFunctions2D(0, 2, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          //xx = yy= 0.0;
          //for(ii=0;ii<nlbf;ii++)
          //{
            //xx += N[ii]*xNode[ii];
            //yy += N[ii]*yNode[ii];
          //}

        dvol = gaussweights1[gp1] * JacMult;

        gxz = computeValueCur(0, dN_dx) + computeValueCur(1, N);
        gyz = computeValueCur(0, dN_dy) + computeValueCur(2, N);

        txz = G*gxz;
        tyz = G*gyz;

        for(ii=0;ii<nlbf;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*N[ii];

          TI   = 3*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;

          Flocal(TI)   += (0.0 - bb1*txz - bb2*tyz) ;
          Flocal(TIp1) += (0.0 - bb3*txz) ;
          Flocal(TIp2) += (0.0 - bb3*tyz) ;

          for(jj=0; jj<nlbf; jj++)
          {
              cc1 = (G*af)*dN_dx[jj];
              cc2 = (G*af)*dN_dy[jj];
              cc3 = (G*af)*N[jj];

              TJ   = 3*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;

              Klocal(TI, TJ)      +=  (bb1*cc1 + bb2*cc2);
              Klocal(TI, TJp1)    +=  (bb1*cc3) ;
              Klocal(TI, TJp2)    +=  (bb2*cc3) ;

              Klocal(TIp1, TJ)    +=  (bb3*cc1);
              Klocal(TIp1, TJp1)  +=  (bb3*cc3) ;
              Klocal(TIp1, TJp2)  +=  0.0 ;

              Klocal(TIp2, TJ)    +=  (bb3*cc2);
              Klocal(TIp2, TJp1)  +=  0.0 ;
              Klocal(TIp2, TJp2)  +=  (bb3*cc3);
          }
        }
  }//gp1
  }//gp2

  Klocal +=  d1*Mlocal;
  Flocal -=  Mlocal*accC;

  return 0;
}
//




int LagrangeElem3DShellQuad4Node::calcInternalForces()
{
  return 0;
}



void LagrangeElem3DShellQuad4Node::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void LagrangeElem3DShellQuad4Node::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void LagrangeElem3DShellQuad4Node::projectStress(int varindex, double* outval)
{
  return;
}


void LagrangeElem3DShellQuad4Node::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void LagrangeElem3DShellQuad4Node::projectIntVar(int index, double* outval)
{
   return;
}


int LagrangeElem3DShellQuad4Node::calcOutput(double u1, double v1)
{
  return 0;
}


void LagrangeElem3DShellQuad4Node::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}

