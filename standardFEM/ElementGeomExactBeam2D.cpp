// 

#include "Debug.h"
#include "MpapTime.h"

#include "ElementGeomExactBeam2D.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"


using namespace std;

extern MpapTime mpapTime;


ElementGeomExactBeam2D::ElementGeomExactBeam2D()
{
  ndof   = 3;
  ndim   = 2;
  npElem = 2;
  nlbf   = npElem;
  nsize  = npElem*ndof;

  if (debug) cout << " constructor ElementGeomExactBeam2D\n\n";
}



ElementGeomExactBeam2D::~ElementGeomExactBeam2D()
{
  if (debug) cout << " destructor ElementGeomExactBeam2D\n\n";
}



void ElementGeomExactBeam2D::prepareElemData()
{
  LagrangeElement::prepareElemData();

  int ii, jj, ind, kk;
  
  return;
}




void ElementGeomExactBeam2D::prepareElemData2()
{
  int ii, jj, ind, kk, ind2;
  
  if(!SolnData->STAGGERED)
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



int ElementGeomExactBeam2D::calcLoadVector()
{
  return 0;
}


//
int ElementGeomExactBeam2D::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();

    int  ii, jj, gp, kk, ll, TI, TIp1, TIp2, count, TJ, TJp1, TJp2, ind1, ind2;

    double  sth0, cth0, uxn[2], uzn[2], btn[2], b[2], af, d1, x0[2], y0[2], x1[2], y1[2], xx[2];

    double  ux, uz, bt, sbt, cbt, dux, duz, dbt, detJ, aa, h;
    double  rho, A, E, G, I, nu, EA, GA, EI, K, SF, NF, BM, kappa, dvol;
    double  fact, fact1, fact2, fact3, fact4, EAdv, GAdv, EIdv;
    double  u1, u2, w1, w2, ddu[2], ddw[2], dx, dy, Nf[2];

    vector<double>  N(nlbf), dN_dx(nlbf);

    VectorXd  res(3), accC(6), velC(6), dispC(6);
    MatrixXd  Mlocal(6,6), Bi(3,6), Bj(3,3), D(3,3), RotMat(6,6), RotMatTrans(6,6);
    D.setZero();
    Bi.setZero();
    Bj.setZero();

    double *elmDat = &(SolnData->ElemProp[elmType].data[0]);

    nGP = 1;
    //b[0] = elmDat[0];
    //b[1] = elmDat[1];
    b[0] = 0.0;
    b[1] = 0.0;
    rho  = elmDat[2];
    A    = elmDat[3];
    I    = elmDat[4];
    E    = elmDat[5];
    nu   = elmDat[6];
    kappa= elmDat[7];
    
    G  = E/(2.0*(1.0+nu));
    EA = E*A;
    EI = E*I;
    GA = G*A*kappa;

    af = SolnData->td(2);
    d1 = SolnData->td(5);
    aa = SolnData->td(10);
    
    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);

    //  rotate nodal displacements and compute nodal positions on element axis
    
    
    x0[0] = GeomData->NodePosOrig[nodeNums[0]][0];
    y0[0] = GeomData->NodePosOrig[nodeNums[0]][1];
    x0[1] = GeomData->NodePosOrig[nodeNums[1]][0];
    y0[1] = GeomData->NodePosOrig[nodeNums[1]][1];

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

    //cout << " bbbbbbbbbbb " << endl;

    x1[0] = x0[0];
    y1[0] = y0[0];
    x1[1] = x0[1];
    y1[1] = y0[1];

    // compute the orientation of the element

    dx = x1[1] - x1[0];
    dy = y1[1] - y1[0];

    h  = sqrt(dx*dx+dy*dy);
    sth0  = dy/h;
    cth0  = dx/h;
    xx[0] = 0.0;
    xx[1] = h;
    
    RotMat.setZero();
    RotMat(0,0) = cth0; RotMat(0,1) = -sth0;
    RotMat(1,0) = sth0; RotMat(1,1) =  cth0;
    RotMat(2,2) = 1.0;

    RotMat(3,3) = cth0; RotMat(3,4) = -sth0;
    RotMat(4,3) = sth0; RotMat(4,4) =  cth0;
    RotMat(5,5) = 1.0;
    
    RotMatTrans = RotMat.transpose();
    
    dispC = RotMatTrans*dispC;
    accC  = RotMatTrans*accC;
    velC  = RotMatTrans*velC;

    uxn[0] = dispC(0);
    uzn[0] = dispC(1);
    btn[0] = dispC(2);
    uxn[1] = dispC(3);
    uzn[1] = dispC(4);
    btn[1] = dispC(5);

    //cout << rho << '\t' << A << '\t' << G << '\t' << EA << '\t' << EI << '\t' << GA << endl;
    //cout << x1[0] << '\t' << y1[0] << '\t' << x1[1] << '\t' << y1[1] << '\t' << h << endl;
    //cout << xx[0] << '\t' << xx[1] << "\n\n" << endl;
    //cout << sth0 << '\t' << cth0 << "\n\n" << endl;
    //cout << " ccccccccccc " << endl;
    
    Klocal.setZero();
    Flocal.setZero();

    //loop over Gauss points (length)
    for(gp=0; gp<nGP; gp++)
    {
      //compute shape functions

      computeLagrangeBFsLine1D(npElem-1, gausspoints[gp], xx, &N[0], &dN_dx[0], detJ);

      //cout << N[0] << '\t' << N[1] << endl;
      //cout << dN_dx[0] << '\t' << dN_dx[1] << '\t' << detJ << "\n" << endl;

      //compute ux, uz, beta in current Gauss point

      ux = uz = bt = dux = duz = dbt = 0.0;
      for(ii=0;ii<npElem;ii++)
      {
        ux  += uxn[ii] * N[ii];
        uz  += uzn[ii] * N[ii];
        bt  += btn[ii] * N[ii];
        dux += uxn[ii] * dN_dx[ii];
        duz += uzn[ii] * dN_dx[ii];
        dbt += btn[ii] * dN_dx[ii];
      }

      sbt = sin(bt);
      cbt = cos(bt);

      //compute average normal strain, shear strain and curvature

      fact = (1.0+dux)*cbt - duz*sbt;

      E = dux + 0.5*(dux*dux + duz*duz);
      G = (1.0+dux)*sbt + duz*cbt;
      K = dbt * fact;

      //compute material response (elastic)

        NF = EA * E;// normal force
        SF = GA * G;// shear force
        BM = EI * K;// bending moment

        //multiply with volume element

        dvol  = gaussweights[gp] * detJ;
        fact1 = dvol * af;
        NF    = NF * dvol;
        SF    = SF * dvol;
        BM    = BM * dvol;
        EAdv  = EA * fact1;
        GAdv  = GA * fact1;
        EIdv  = EI * fact1;

        fact1 = (+ SF * cbt - BM * dbt * sbt) * af;
        fact2 = (- SF * sbt - BM * dbt * cbt) * af;
        
        D(0,0) = EAdv;
        D(1,1) = GAdv;
        D(2,2) = EIdv;

        for(ii=0;ii<nlbf;ii++)
        {
          TI   =  3*ii;
          TIp1 =  TI+1;
          TIp2 =  TI+2;

          Bi(0,TI)   = (1.0+dux) * dN_dx[ii];
          Bi(0,TIp1) = duz * dN_dx[ii];
          Bi(0,TIp2) = 0.0;

          Bi(1,TI)   = sbt * dN_dx[ii];
          Bi(1,TIp1) = cbt * dN_dx[ii];
          Bi(1,TIp2) = fact * N[ii];

          Bi(2,TI)   = dbt*cbt * dN_dx[ii];
          Bi(2,TIp1) = - dbt*sbt * dN_dx[ii];
          Bi(2,TIp2) = fact * dN_dx[ii] - G*dbt * N[ii];
        }

        //printMatrix(Bi);	printf("\n\n\n");
        //printMatrix(D);	printf("\n\n\n");

        res(0) = NF;
        res(1) = SF;
        res(2) = BM;

        Klocal += (Bi.transpose()*(D*Bi));
        Flocal -= (Bi.transpose()*res);

        //printMatrix(Klocal);	printf("\n\n\n");

        //compute geometrical part of stiffness matrix

        for(ii=0;ii<nlbf;ii++)
        {
          TI   =  3*ii;
          TIp1 =  TI+1;
          TIp2 =  TI+2;

          for(jj=0;jj<nlbf;jj++)
          {
            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            Klocal(TI,TJ)     += ( dN_dx[ii]*NF*dN_dx[jj] * af);
            Klocal(TIp1,TJp1) += ( dN_dx[ii]*NF*dN_dx[jj] * af);

            fact3 =  dN_dx[ii]*BM*cbt*dN_dx[jj] * af;
            fact4 = -dN_dx[ii]*BM*sbt*dN_dx[jj] * af;

            Klocal(TI,TJp2)   += (fact3 + dN_dx[ii]*fact1*N[jj] );
            Klocal(TIp1,TJp2) += (fact4 + dN_dx[ii]*fact2*N[jj] );
            Klocal(TIp2,TJ)   += (fact3 + N[ii]*fact1*dN_dx[jj] );
            Klocal(TIp2,TJp1) += (fact4 + N[ii]*fact2*dN_dx[jj] );

            Klocal(TIp2,TJp2) += (N[ii]*(-SF*G-BM*dbt*fact)*N[jj] - dN_dx[ii]*BM*G*N[jj] - N[ii]*BM*G*dN_dx[jj]) * af;
          }//for(jj
        }//for(ii
    }//for(gp

    //cout << " element = " << nodeNums[0] << '\t' << nodeNums[1] << endl;

    //printMatrix(Klocal); printf("\n\n");
    //printVector(Flocal); printf("\n\n");

    //body forces

    //inertia
    fact1 = rho*A*h/6.0;
    fact2 = rho*I*h/6.0;

    fact3 = 2.0*fact1;
    fact4 = 2.0*fact2;
    
    Mlocal.setZero();

    Mlocal(0,0) = fact3;    Mlocal(3,3) = fact3;    Mlocal(0,3) = fact1;    Mlocal(3,0) = fact1;
    Mlocal(1,1) = fact3;    Mlocal(4,4) = fact3;    Mlocal(1,4) = fact1;    Mlocal(4,1) = fact1;
    Mlocal(2,2) = fact4;    Mlocal(5,5) = fact4;    Mlocal(2,5) = fact2;    Mlocal(5,2) = fact2;

    Klocal += (d1*Mlocal);
    Flocal -= (Mlocal*accC);

    //Klocal /= aa;

    Klocal = (RotMat*Klocal)*RotMatTrans;
    Flocal = RotMat*Flocal;

    //printMatrix(Klocal); printf("\n\n");
    //printVector(Flocal); printf("\n\n");

    return 0;
}
//




/*
int ElementGeomExactBeam2D::calcStiffnessAndResidual()
{
    int  ii, jj, gp, kk, ll, TI, TIp1, TIp2, count, TJ, TJp1, TJp2, ind1, ind2;

    double  sth0, cth0, uxn[2], uzn[2], btn[2], b[2], af, am, d1, x0[2], y0[2], x1[2], y1[2], xx[2];

    double  ux, uz, bt, sbt, cbt, dux, duz, dbt, detJ, h;
    double  rho, A, E, G, I, nu, EA, GA, EI, K, SF, NF, BM, kappa, dvol;
    double  fact, fact1, fact2, fact3, fact4, EAdv, GAdv, EIdv;
    double  u1, u2, w1, w2, ddu[2], ddw[2], dx, dy;

    double  N[nlbf], dN_dx[nlbf], Nf[2];

    VectorXd  res(3), vecF(6), dispC(6), qVec(12), qdotVec(12);
    MatrixXd  Bi(3,6), Bj(3,3), D(3,3), RotMat(12,12), RotMatTrans(12,12);
    MatrixXd  matM(6, 6), matC(6, 6), matK(6, 6);

    elmDat = &(SolnData->ElemProp.data[0]);

    nGP1 = 1;
    b[0] = elmDat[0];
    b[1] = elmDat[1];
    rho  = elmDat[2];
    A    = elmDat[3];
    I    = elmDat[4];
    E    = elmDat[5];
    nu   = elmDat[6];
    kappa= elmDat[7];
    
    G  = E/(2.0*(1.0+nu));
    EA = E*A;
    EI = E*I;
    GA = G*A*kappa;

    am = SolnData->td(1);
    af = SolnData->td(2);
    af = 1.0;
    d1 = SolnData->td(8);
    
    double  *gausspoints  = &(GeomData->gausspoints1[0]);
    double  *gaussweights = &(GeomData->gaussweights1[0]);

    //  rotate nodal displacements and compute nodal positions on element axis
    
    //cout << " ccccccccccc " << endl;
    
    x0[0] = GeomData->NodePosOrig[nodeNums[0]][0];
    y0[0] = GeomData->NodePosOrig[nodeNums[0]][1];
    x0[1] = GeomData->NodePosOrig[nodeNums[1]][0];
    y0[1] = GeomData->NodePosOrig[nodeNums[1]][1];

    for(ii=0;ii<npElem;ii++)
    {
      ind1 = ndof*ii;
      ind2 = nodeNums[ii]*ndof;

      for(kk=0;kk<ndof;kk++)
      {
        qVec(ind1+kk)      =  SolnData->var1Cur[ind2+kk];
        qdotVec(ind1+kk)   =  SolnData->var1DotCur[ind2+kk];
      }
    }

    x1[0] = x0[0];
    y1[0] = y0[0];
    x1[1] = x0[1];
    y1[1] = y0[1];

    // compute the orientation of the element

    dx = x1[1] - x1[0];
    dy = y1[1] - y1[0];

    h  = sqrt(dx*dx+dy*dy);
    sth0  = dy/h;
    cth0  = dx/h;
    xx[0] = 0.0;
    xx[1] = h;
    
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

    //qVec.segment(0,6)  = RotMatTrans*qVec.segment(0,6);
    //qVec.segment(6,6)  = RotMatTrans*qVec.segment(6,6);

    //qdotVec.segment(0,6)  = RotMatTrans*qdotVec.segment(0,6);
    //qdotVec.segment(6,6)  = RotMatTrans*qdotVec.segment(6,6);

    qVec  = RotMatTrans*qVec;
    qdotVec  = RotMatTrans*qdotVec;

    uxn[0] = qVec(0);
    uzn[0] = qVec(1);
    btn[0] = qVec(2);
    uxn[1] = qVec(6);
    uzn[1] = qVec(7);
    btn[1] = qVec(8);

    //cout << rho << '\t' << A << '\t' << G << '\t' << EA << '\t' << EI << '\t' << GA << endl;
    //cout << x1[0] << '\t' << y1[0] << '\t' << x1[1] << '\t' << y1[1] << '\t' << h << endl;
    //cout << xx[0] << '\t' << xx[1] << "\n\n" << endl;
    //cout << sth0 << '\t' << cth0 << "\n\n" << endl;
    //cout << " ccccccccccc " << endl;

    D.setZero();
    matK.setZero();
    vecF.setZero();

    //loop over Gauss points (length)
    for(gp=0;gp<nGP1;gp++)
    {
      //compute shape functions

      computeLagrangeBFs1D(npElem-1, gausspoints[gp], xx, N, dN_dx, detJ);

      //cout << N[0] << '\t' << N[1] << endl;
      //cout << dN_dx[0] << '\t' << dN_dx[1] << '\t' << detJ << "\n" << endl;

      //compute ux, uz, beta in current Gauss point

      ux = uz = bt = dux = duz = dbt = 0.0;
      for(ii=0;ii<npElem;ii++)
      {
        ux  += uxn[ii] * N[ii];
        uz  += uzn[ii] * N[ii];
        bt  += btn[ii] * N[ii];
        dux += uxn[ii] * dN_dx[ii];
        duz += uzn[ii] * dN_dx[ii];
        dbt += btn[ii] * dN_dx[ii];
      }

      sbt = sin(bt);
      cbt = cos(bt);

      //compute average normal strain, shear strain and curvature

      fact = (1.0+dux)*cbt - duz*sbt;

      E = dux + 0.5*(dux*dux + duz*duz);
      G = (1.0+dux)*sbt + duz*cbt;
      K = dbt * fact;

      //compute material response (elastic)

        NF = EA * E;// normal force
        SF = GA * G;// shear force
        BM = EI * K;// bending moment

        //multiply with volume element

        dvol  = gaussweights[gp] * detJ;
        fact1 = dvol * af;
        NF    = NF * dvol;
        SF    = SF * dvol;
        BM    = BM * dvol;
        EAdv  = EA * fact1;
        GAdv  = GA * fact1;
        EIdv  = EI * fact1;

        fact1 = (+ SF * cbt - BM * dbt * sbt) * af;
        fact2 = (- SF * sbt - BM * dbt * cbt) * af;
        
        D(0,0) = EAdv;
        D(1,1) = GAdv;
        D(2,2) = EIdv;

        for(ii=0;ii<nlbf;ii++)
        {
          TI   =  3*ii;
          TIp1 =  TI+1;
          TIp2 =  TI+2;

          Bi(0,TI)   = (1.0+dux) * dN_dx[ii];
          Bi(0,TIp1) = duz * dN_dx[ii];
          Bi(0,TIp2) = 0.0;
          Bi(1,TI)   = sbt * dN_dx[ii];
          Bi(1,TIp1) = cbt * dN_dx[ii];
          Bi(1,TIp2) = fact * N[ii];
          Bi(2,TI)   = dbt*cbt * dN_dx[ii];
          Bi(2,TIp1) = - dbt*sbt * dN_dx[ii];
          Bi(2,TIp2) = fact * dN_dx[ii] - G*dbt * N[ii];
        }

        //printMatrix(Bi);	printf("\n\n\n");
        //printMatrix(D);	printf("\n\n\n");

        res(0) = NF;
        res(1) = SF;
        res(2) = BM;

        matK += (Bi.transpose()*(D*Bi));
        vecF -= (Bi.transpose()*res);

        //printMatrix(Klocal);	printf("\n\n\n");

        //compute geometrical part of stiffness matrix

        for(ii=0;ii<nlbf;ii++)
        {
          TI   =  3*ii;
          TIp1 =  TI+1;
          TIp2 =  TI+2;

          for(jj=0;jj<nlbf;jj++)
          {
            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            matK(TI,TJ)     += dN_dx[ii]*NF*dN_dx[jj] * af;
            matK(TIp1,TJp1) += dN_dx[ii]*NF*dN_dx[jj] * af;

            fact3 =  dN_dx[ii]*BM*cbt*dN_dx[jj] * af;
            fact4 = -dN_dx[ii]*BM*sbt*dN_dx[jj] * af;

            matK(TI,TJp2)   += (fact3 + dN_dx[ii]*fact1*N[jj] );
            matK(TIp1,TJp2) += (fact4 + dN_dx[ii]*fact2*N[jj] );
            matK(TIp2,TJ)   += (fact3 + N[ii]*fact1*dN_dx[jj] );
            matK(TIp2,TJp1) += (fact4 + N[ii]*fact2*dN_dx[jj] );

            matK(TIp2,TJp2) += (N[ii]*(-SF*G-BM*dbt*fact)*N[jj] - dN_dx[ii]*BM*G*N[jj] - N[ii]*BM*G*dN_dx[jj]) * af;
          }//for(jj
        }//for(ii
    }//for(gp

    //cout << " element = " << nodeNums[0] << '\t' << nodeNums[1] << endl;

    //printMatrix(Klocal); printf("\n\n");
    //printVector(Flocal); printf("\n\n");

    //body forces

    //inertia
    fact1 = rho*A*h/6.0;
    fact2 = rho*I*h/6.0;

    fact3 = 2.0*fact1;
    fact4 = 2.0*fact2;
    
    matC.setZero();
    matM.setZero();

    matM(0,0) = fact3;    matM(3,3) = fact3;    matM(0,3) = fact1;    matM(3,0) = fact1;
    matM(1,1) = fact3;    matM(4,4) = fact3;    matM(1,4) = fact1;    matM(4,1) = fact1;
    matM(2,2) = fact4;    matM(5,5) = fact4;    matM(2,5) = fact2;    matM(5,2) = fact2;
    

    // create local stiffness matrix in state space form

    vector<int> perm1(6, 1), perm2(6, 1);

    perm1[0] = 0;  perm1[1] = 1;  perm1[2] = 2;  perm1[3] = 6;  perm1[4] = 7;  perm1[5] = 8;
    perm2[0] = 3;  perm2[1] = 4;  perm2[2] = 5;  perm2[3] = 9;  perm2[4] =10;  perm2[5] =11;

    Klocal.setZero();
    Flocal.setZero();
    af = SolnData->td(2);

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

        Klocal(perm2[ii], perm1[jj])  += ( af*matK(ii,jj) );
        Klocal(perm2[ii], perm2[jj])  += ( d1*matM(ii,jj) );
        Klocal(perm2[ii], perm2[jj])  += ( af*matC(ii,jj) );

        fact2 += matM(ii,jj) * qdotVec(perm2[jj]);
        fact2 += matC(ii,jj) * qVec(perm2[jj]);
      }
      Flocal(perm1[ii]) -= fact1;
      Flocal(perm2[ii]) -= fact2;
    }

    Klocal = (RotMat*Klocal)*RotMatTrans;

    Flocal = RotMat*Flocal;


    //printMatrix(Klocal); printf("\n\n");
    //printVector(Flocal); printf("\n\n");

   return 0;
}
*/




int ElementGeomExactBeam2D::calcInternalForces()
{
  return 0;
}



void ElementGeomExactBeam2D::discreteContourplot(int vartype, int varindex, int index, int nCol, double umin, double umax)
{
  return;
}


void ElementGeomExactBeam2D::projectToKnots(bool extrapolateFlag, int vartype, int varindex, int index)
{
  return;
}


void ElementGeomExactBeam2D::projectStress(int varindex, double* outval)
{
  return;
}



void ElementGeomExactBeam2D::projectStrain(int vartype, int varindex, double* outval)
{
  return;
}



void ElementGeomExactBeam2D::projectIntVar(int index, double* outval)
{
  return;
}


int ElementGeomExactBeam2D::calcOutput(double u1, double v1)
{
  return 0;
}



void ElementGeomExactBeam2D::toPostprocess(int vartype, int varindex, int type, SparseMatrixXd&  coeffMat, VectorXd& rhsVec)
{
  return;
}





