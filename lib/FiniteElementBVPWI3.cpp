
#include <iostream>

#include "FiniteElementBVPWI.h"
#include "DomainTree.h"
#include "DataBlockTemplate.h"
#include "Plot.h"
#include "ComputerTime.h"
#include "MathGeom.h"


extern ComputerTime computerTime;
extern DomainTree   domain;
extern Plot         plot;



using namespace std;






void FiniteElementBVPWI::calculateDerivatives1(void)
{
  // calculate partial derivatives  d gf / d uf,  d gf / d ul  and  d gl / d uf,
  // do this by looping once over the layer[0] elements

  char fct[] = "FiniteElementBVPWI::calculateDerivatives1";

  computerTime.go(fct);

  //cout << fct << "\n";

  int i, j, e, e2, nst, nenEl, ndfEl, 
      cff, cfl, clf,
      ir, nr, jr, 
      ic, nc, jc;

  BndNode  *brPtr, *bcPtr;

  for (e2=0; e2<nLayElem2; e2++)
  {
    e = layElem[0][e2] - 1;

    if (elem[e]->calcStiffnessAndResidual() != 0) prgError(1,fct,"local stiffness error!");

    nenEl = elem[e]->nen();
    ndfEl = elem[e]->ndf();
    nst   = ndfEl * nenEl;

    cff = 0;
    cfl = 0;
    clf = 0;

    for (ir=0; ir<nenEl; ir++)
    {
      nr = ixbu[e2][ir];

      if (nr > -1) // if row node is a boundary node
      {
        brPtr = &(bndNode[nr]);

        for (jr=0; jr<brPtr->freeU.n; jr++)
        {
          for (ic=0; ic<nenEl; ic++)
          {
            nc = ixbu[e2][ic];
        
            if (nc > -1) // if column node is a boundary node
            {
              bcPtr = &(bndNode[nc]);

              for (jc=0; jc<bcPtr->freeU.n; jc++)
              {
                calcF2Fu(jr,jc, *brPtr,*bcPtr, ir*ndfEl,ic*ndfEl, nst, pffu[e2][cff++]->dgfduf.x);
              }
            }
            nc = ixlu[e2][ic];
        
            if (nc > -1) // if column node is a layer node
            {
              calcF2Lu(jr, *brPtr, ir*ndfEl, ic*ndfEl, ndfEl, nst, pflu[e2][cfl++]->dgfdul.x);
            }
          }
        }
      }
    }
    for (ic=0; ic<nenEl; ic++)
    {
      nc = ixbu[e2][ic];

      if (nc > -1) // if column node is a boundary node
      {
        bcPtr = &(bndNode[nc]);

        for (jc=0; jc<bcPtr->freeU.n; jc++)
        {
          for (ir=0; ir<nenEl; ir++)
          {
            nr = ixlu[e2][ir];
       
            if (nr > -1) // if row node is a layer node
            {
              calcL2Fu(jc, *bcPtr, ir*ndfEl, ic*ndfEl, ndfEl, nst, pflu[e2][clf++]->dglduf.x);
            }
          }
        }
      }
    }
  }
  ctimDeriv1 += computerTime.stop(fct);

  doneDeriv1 = true;

  /*int dp;

  for (i=0; i<freeNode.n; i++)
    for (dp=0; dp<freeNode[i].toFree.n; dp++)
      for (j=0; j<freeNode[i].toFree[dp].n; j++)

        cout << freeNode[i].toFree[dp][j].dgfduf << " ives1\n";*/

  return;
}












void FiniteElementBVPWI::calculateDerivatives2(void)
{
  // calculate partial derivatives  d gf / d xf,  d gf / d xl  and  d gl / d xf,
  // do this by looping once over the layer[0] elements

  char fct[] = "FiniteElementBVPWI::calculateDerivatives2";

  if (!isALE()) return;

  if (maxMeshDepth < 0) return;

  //cout << fct << "\n";

  computerTime.go(fct);

  int i, j, e, e2, nstu, nenEl, ndfEl, 
      cff, cfl, clf,
      ir, nr, jr,
      ic, nc, jc;

  double fact;

  BndNode *brPtr, *bcPtr;

  for (e2=0; e2<layElem[0].n; e2++)
  {
    e = layElem[0][e2] - 1;

    if (elem[e]->calcMeshDerivatives() != 0) prgError(1,fct,"local stiffness error!");

    nenEl = elem[e]->nen();
    ndfEl = elem[e]->ndf();
    nstu  = ndfEl * nenEl;

    cff = 0;
    cfl = 0;
    clf = 0;

    if (e2 < nLayElem2)
    {
      for (ir=0; ir<nenEl; ir++)
      {
        nr = ixbu[e2][ir];

        if (nr > -1) // if row node is a boundary node
        {
          brPtr = &(bndNode[nr]);

          for (jr=0; jr<brPtr->freeU.n; jr++)
          {
            for (ic=0; ic<nenEl; ic++)
            {
              nc = ixbx[e2][ic];

              if (nc > -1) // if column node is a boundary node
              {
                bcPtr = &(bndNode[nc]);

                for (jc=0; jc<bcPtr->freeX.n; jc++)
                {
                  calcF2Fx(jr,jc,*brPtr,*bcPtr,ir*ndfEl,ic*ndm,nstu,pffx[e2][cff++]->dgfdxf.x);
                }
              }
              nc = ixlx[e2][ic];

              if (nc > -1) // if column node is a layer node
              {
                calcF2Lx(jr, *brPtr, ir*ndfEl, ic*ndm, nstu, pflx[e2][cfl++]->dgfdxl.x);
              }
            }
          }
        }
      }
    }

    for (ir=0; ir<nenEl; ir++)
    {
      nr = ixlu[e2][ir];
    
      if (nr > -1) // if row node is a layer node
      {
        for (ic=0; ic<nenEl; ic++)
        {
          nc = ixbx[e2][ic];
       
          if (nc > -1) // if column node is a boundary node
          {
            bcPtr = &(bndNode[nc]);
       
            for (jc=0; jc<bcPtr->freeX.n; jc++)
            {
              calcL2Fx(jc, *bcPtr, ir*ndfEl, ic*ndm, ndfEl, nstu, plfx[e2][clf++]->dgldxf.x);
            }
          }
        }
      }
    }
  }
  ctimDeriv2 += computerTime.stop(fct);

  doneDeriv2 = true;

  /*int dp;

  for (i=0; i<freeNode.n; i++)
    for (dp=0; dp<freeNode[i].toFree.n; dp++)
      for (j=0; j<freeNode[i].toFree[dp].n; j++)
      {
        cout << freeNode[i].toFree[dp][j].dgfduf << " u ives2\n";
        cout << freeNode[i].toFree[dp][j].dgfdxf << " x ives2\n";
      }*/

  return;
}







void FiniteElementBVPWI::calcF2Fu(int jr, int jc, BndNode &br, BndNode &bc,
                                  int ir0, int ic0, int nst, double *su)
{
  char fct[] = "FiniteElementBVPWI::calcF2Fu";

  int nlr, *DOFr, nlc, *DOFc, lr, lc;

  double ar, ac, fact;

  if (br.giveDerivN(jr) == 1)
  {
    ar = br.giveDeriv(jr)[0];

    nlr  = br.freeU[jr]->dofU->n;
    DOFr = br.freeU[jr]->dofU->x;

    if (bc.getDerivN(jc) == 1)
    {
      ac = bc.getDeriv(jc)[0];

      fact = ar * ac;

      nlc  = bc.freeU[jc]->dofU->n;
      DOFc = bc.freeU[jc]->dofU->x;

      for (lr=0; lr<nlr; lr++)
      {
        for (lc=0; lc<nlc; lc++)
        {
          su[lc*nlr+lr] += s[(ic0+DOFc[lc])*nst + ir0+DOFr[lr]] * fact;
        }
      }
    }
    else
    {
      prgError(1,fct,"not yet implemented!");
    }
  }
  else
  {
    prgError(2,fct,"not yet implemented!");
  }

  return;
}




void FiniteElementBVPWI::calcL2Fu(int jc, BndNode &bc,
                                  int ir0, int ic0, int ndfEl, int nst, double *dglduf)
{
  char fct[] = "FiniteElementBVPWI::calcL2Fu";

  int nlc, *DOFc, lr, lc;

  double ac;

  if (bc.giveDerivN(jc) == 1)
  {
    nlc  = bc.freeU[jc]->dofU->n;
    DOFc = bc.freeU[jc]->dofU->x;

    ac = bc.getDeriv(jc)[0];

    for (lr=0; lr<ndfEl; lr++)
    {
      for (lc=0; lc<nlc; lc++)
      {
        dglduf[lc*ndf+lr] += s[(ic0+DOFc[lc])*nst + ir0+lr] * ac;
      }
    }
  }
  else
  {
    prgError(1,fct,"not yet implemented!");
  }

  return;
}




void FiniteElementBVPWI::calcF2Lu(int jr, BndNode &br,
                                  int ir0, int ic0, int ndfEl, int nst, double *dgfdul)
{
  char fct[] = "FiniteElementBVPWI::calcF2Lu";

  int nlr, *DOFr, lr, lc;

  double ar;

  if (br.giveDerivN(jr) == 1)
  {
    nlr  = br.freeU[jr]->dofU->n;
    DOFr = br.freeU[jr]->dofU->x;

    ar = br.getDeriv(jr)[0];

    for (lr=0; lr<nlr; lr++)
    {
      for (lc=0; lc<ndfEl; lc++)
      {
        dgfdul[lc*nlr+lr] += s[(ic0+lc)*nst + ir0+DOFr[lr]] * ar;
      }
    }
  }
  else
  {
    prgError(1,fct,"not yet implemented!");
  }

  return;
}







void FiniteElementBVPWI::calcF2Fx(int jr, int jc, BndNode &br, BndNode &bc,
                                  int ir0, int ic0, int nst, double *sx)
{
  char fct[] = "FiniteElementBVPWI::calcF2Fx";

  int nlr, *DOFr, nlc, *DOFc, lr, lc;

  double ar, ac, fact;

  if (br.giveDerivN(jr) == 1)
  {
    ar = br.giveDeriv(jr)[0];

    nlr  = br.freeU[jr]->dofU->n;
    DOFr = br.freeU[jr]->dofU->x;

    if (bc.getDerivN(jc) == 1)
    {
      ac = bc.getDeriv(jc)[0];

      if (bc.freeX[jc]->isLagr) ac *= dxdu;

      fact = ar * ac;

      nlc  = bc.freeX[jc]->dofX->n;
      DOFc = bc.freeX[jc]->dofX->x;

      for (lr=0; lr<nlr; lr++)
      {
        for (lc=0; lc<nlc; lc++)
        {
          sx[lc*nlr+lr] += s[(ic0+DOFc[lc])*nst + ir0+DOFr[lr]] * fact;
        }
      }
    }
    else
    {
      prgError(1,fct,"not yet implemented!");
    }
  }
  else
  {
    prgError(2,fct,"not yet implemented!");
  }

  return;
}




void FiniteElementBVPWI::calcL2Fx(int jc, BndNode &bc,
                                  int ir0, int ic0, int ndfEl, int nst, double *dgldxf)
{
  char fct[] = "FiniteElementBVPWI::calcL2Fx";

  int nlc, *DOFc, lr, lc;

  double ac;

  if (bc.giveDerivN(jc) == 1)
  {
    nlc  = bc.freeX[jc]->dofX->n;
    DOFc = bc.freeX[jc]->dofX->x;

    ac = bc.getDeriv(jc)[0];

    if (bc.freeX[jc]->isLagr) ac *= dxdu;

    for (lr=0; lr<ndfEl; lr++)
    {
      for (lc=0; lc<nlc; lc++)
      {
        dgldxf[lc*ndf+lr] += s[(ic0+DOFc[lc])*nst + ir0+lr] * ac;
      }
    }
  }
  else
  {
    prgError(1,fct,"not yet implemented!");
  }

  return;
}




void FiniteElementBVPWI::calcF2Lx(int jr, BndNode &br,
                                  int ir0, int ic0, int nst, double *dgfdxl)
{
  char fct[] = "FiniteElementBVPWI::calcF2Lx";

  int nlr, *DOFr, lr, lc;

  double ar;

  if (br.giveDerivN(jr) == 1)
  {
    nlr  = br.freeU[jr]->dofU->n;
    DOFr = br.freeU[jr]->dofU->x;

    ar = br.getDeriv(jr)[0];

    for (lr=0; lr<nlr; lr++)
    {
      for (lc=0; lc<ndm; lc++)
      {
        dgfdxl[lc*nlr+lr] += s[(ic0+lc)*nst + ir0+DOFr[lr]] * ar;
      }
    }
  }
  else
  {
    prgError(1,fct,"not yet implemented!");
  }

  return;
}







void FiniteElementBVPWI::calculateDerivatives3(void)
{
  // loop over columns of  d g_all_non_free / d uf:
  //
  //   extract non-zero coefficients from dglduf,
  //   calculate corresponding column of total derivative  d ul / d uf,
  //   use that column of  d ul / d uf  in  d gf / d ul  x  d ul / d uf,
  //   add result of multiplication to  dgfduf

  char fct[] = "FiniteElementBVPWI::calculateDerivatives3";

  //cout << fct << "\n";

  if (nequ == 0) return;

  computerTime.go(fct);

  int nlc, nlr, d, i, j, ic, l, lc, ir, lr, rr, nr, n, *IDU = idu.x,
      iRHS, jRHS, nRHS, ic0, lc0;

  double *su, *dgfdul, *COL;

  DataToFree *tmpF;
  DataToBLay *tmpBL;

  if (solver->currentStatus == ASSEMBLY_OK) solver->factorise();

  if (solver->currentStatus != FACTORISE_OK) prgError(1,fct,"Solver status != FACTORISE_OK !");

  nRHS = 1;

  col.setDim(nRHS * nequ);

  ic = 0;
  lc = 0;

  do
  {
    ic0 = ic;
    lc0 = lc;

    iRHS = 0;

    //col.zero();
 
    COL = col.x;

    while (ic < freeNode.n && iRHS < nRHS)
    {
      if (!freeNode[ic].isLagr)
      {
        nlc = freeNode[ic].dofU->n;
 
        while (lc < nlc && iRHS < nRHS)
        {
          //cout << ic << "." << lc << " -> build column\n";
        
          for (i=0; i<nequ; i++) COL[i] = 0.;
        
          for (i=0; i<freeNode[ic].toBLay.n; i++)
          {
            nr = freeNode[ic].toBLay[i].nd - 1;
        
            for (lr=0; lr<ndf; lr++)
            {
              rr = IDU[nr*ndf+lr] - 1;
        
              if (rr > -1)  COL[rr] = freeNode[ic].toBLay[i].dglduf[lc*ndf+lr];
            }
          }

          lc++;

          iRHS++;

          COL += nequ;
        }
      }
      lc = 0; 
      ic++;
    }

    if (iRHS > 0)
    {
      //cout << iRHS << " -> " << col << " in\n";
    
      computerTime.go("backsubstitution3");
     
      //cout << "back 3\n";
    
      solver->solve(col.x,iRHS);  // -> column now holds  d u / d uf
     
      //cout << col << " out\n";

      ctimBack3 += computerTime.stop("backsubstitution3");
      
      jRHS = 0;
    
      COL = col.x;
 
      ic = ic0;
    
      lc = lc0;
    
      while (jRHS < iRHS)
      {
        if (!freeNode[ic].isLagr)
        {
          nlc = freeNode[ic].dofU->n;
    
          while (jRHS < iRHS)
          {
            // extract  d ul / d uf  from col and use in  d gf / d ul  x  d ul / d uf
           
            //cout << ic << "." << lc << " -> use column\n";
    
            for (d=0; d<freeNode[ic].toFree.n; d++)
            {
              for (j=0; j<freeNode[ic].toFree[d].n; j++)
              {
                tmpF = &(freeNode[ic].toFree[d][j]);
                su   = tmpF->invPtr->dgfduf.x;
                ir   = tmpF->nd;
                nlr  = freeNode[ir].dofU->n;
           
                for (i=0; i<freeNode[ir].toBLay.n; i++)
                {
                  tmpBL  = &(freeNode[ir].toBLay[i]);
                  n      = tmpBL->nd - 1;
                  dgfdul = tmpBL->dgfdul.x;
           
                  for (l=0; l<ndf; l++)
                  {
                    rr = IDU[n*ndf+l] - 1;
           
                    if (rr > -1)
                    {
                      for (lr=0; lr<nlr; lr++) su[lc*nlr+lr] -= dgfdul[l*nlr+lr] * COL[rr];
                    }
                  }
                }
              }
            }
            COL += nequ;
            jRHS++;
            if (++lc == nlc) { lc = 0; break; }
          }
        }
        if (lc == 0) ic++;
      }
    }
  } while (ic < freeNode.n);

  ctimDeriv3 += computerTime.stop(fct);

  doneDeriv3 = true;

  //cout << fct << ": quit it!\n\n";

  return;
}










void FiniteElementBVPWI::calculateDerivatives4(void)
{
  // calculate partial derivative  d ml / d xf,
  // do this by looping once over the layer elements

  char fct[] = "FiniteElementBVPWI::calculateDerivatives4";

  if (!isALE()) return;

  if (neqx == 0) return;

  if (maxMeshDepth < 1) return;

  computerTime.go(fct);

  //cout << fct << "\n";

  int clf, e, e2, e3, nenEl, nst,
      ir, nr, 
      ic, jc, nc;

  double ac;

  BndNode *bcPtr;

  for (e3=0; e3<nLayElem3-nLayElem1; e3++)
  {
    e2 = e3 + nLayElem1;

    e = layElem[0][e2] - 1;

    if (elem[e]->calcStiffnessAndResidualMesh() != 0) prgError(1,fct,"local stiffness error!");

    nenEl = elem[e]->nen();
    nst   = nenEl * ndm;

    clf = 0;

    for (ic=0; ic<nenEl; ic++)
    {
      nc = ixbx[e2][ic];

      if (nc > -1) // if column node is a boundary node (x)
      {
        bcPtr = &(bndNode[nc]);

        for (jc=0; jc<bcPtr->freeX.n; jc++)
        {
          for (ir=0; ir<nenEl; ir++)
          {
            nr = ixlx[e2][ir];
       
            if (nr > -1) // if row node is a layer node (m)
            {
              calcL2Fm(jc, *bcPtr, ir*ndm, ic*ndm, nst, plfm[e3][clf++]->dmldxf.x);
            }
          }
        }
      }
    }
  }
  ctimDeriv4 += computerTime.stop(fct);

  doneDeriv4 = true;

  return;
}









void FiniteElementBVPWI::calcL2Fm(int jc, BndNode &bc, int ir0, int ic0, int nst, double *dmldxf)
{
  char fct[] = "FiniteElementBVPWI::calcL2Fm";

  int nlc, *DOFc, lr, lc;

  double ac;

  if (bc.giveDerivN(jc) == 1)
  {
    nlc  = bc.freeX[jc]->dofX->n;
    DOFc = bc.freeX[jc]->dofX->x;

    ac = bc.getDeriv(jc)[0];

    if (bc.freeX[jc]->isLagr) ac *= dxdu;

    for (lc=0; lc<nlc; lc++)
    {
      for (lr=0; lr<ndm; lr++)
      {
        dmldxf[lc*ndm+lr] += s[(ic0+DOFc[lc])*nst + ir0+lr] * ac;
      }
    }
  }
  else
  {
    prgError(1,fct,"not yet implemented!");
  }
 
  return;
}









void FiniteElementBVPWI::calculateDerivatives5(void)
{
  // calculate  d xm / d xf

  char fct[] = "FiniteElementBVPWI::calculateDerivatives5";

  if (!isALE()) return;

  if (neqx == 0) return;

  if (maxMeshDepth < 1) return;

  computerTime.go(fct);

  //cout << fct << "\n";

  int d, i, j, ic, l, lc, ir, lr, nr, rr, nlc, nlr, *IDX = idx.x;

  double *xx;

  DataToMesh *tmpMesh;

  if (solverMesh->currentStatus == ASSEMBLY_OK) solverMesh->factorise();

  if (solverMesh->currentStatus != FACTORISE_OK) prgError(1,fct,"Solver status != FACTORISE_OK !");

  for (ic=0; ic<freeNode.n; ic++)
  {
    if (freeNode[ic].toBLayMesh.n > 0)
    {
      nlc = freeNode[ic].dofX->n;

      for (lc=0; lc<nlc; lc++)
      {
        for (i=0; i<neqx; i++) col[i] = 0.;
    
        for (i=0; i<freeNode[ic].toBLayMesh.n; i++)
        {
          xx = freeNode[ic].toBLayMesh[i].dmldxf.x;

          nr = freeNode[ic].toBLayMesh[i].nd - 1;
   
          for (lr=0; lr<ndm; lr++)
          {
            rr = IDX[nr*ndm+lr] - 1;
    
            if (rr > -1) col[rr] = xx[lc*ndm+lr];
          }
        }
        //for (i=0; i<neqx; i++) cout << col[i] << " in\n"; cout << "\n";
   
        computerTime.go("backsubstitution5");
    
        solverMesh->solve(col.x);  // -> column now holds  d xm / d xf
    
        ctimBack5 += computerTime.stop("backsubstitution5");
    
        //for (i=0; i<neqx; i++) cout << col[i] << " out\n"; cout << "\n\n\n";
    
        // extract  d xm / d xf  from col and store in freeToMesh.xx
   
        for (d=0; d<freeNode[ic].toMesh.n; d++)
        {
          for (j=0; j<freeNode[ic].toMesh[d].n; j++)
          {
            tmpMesh = &(freeNode[ic].toMesh[d][j]);
            nr      = tmpMesh->nd;
            xx      = tmpMesh->xx.x;

            for (lr=0; lr<ndm; lr++)
            {
              rr = IDX[nr*ndm+lr] - 1;
    
              if (rr > -1) xx[lc*ndm+lr] = col[rr]; //else xx[lc*ndm+lr] = 0.;
            }
          }
        }
      }
    }
  }

  ctimDeriv5 += computerTime.stop(fct);

  doneDeriv5 = true;

  return;
}












void FiniteElementBVPWI::calculateDerivatives6(void)
{
  // perform multiplication    d gf / d xl * d xl / d xf

  char fct[] = "FiniteElementBVPWI::calculateDerivatives6";

  if (!isALE()) return;

  if (neqx == 0) return;

  if (maxMeshDepth < 1) return;

  //cout << fct << "\n";

  computerTime.go(fct);

  int i0, i, ic, ir, j, d, jr, lr, lc, nlc, nlr, *IDX = idx.x;

  double *dgfdxf, *dxldxf, *dgfdxl, *dxmdxf;

  DataToFree *tmpF;
  DataToMesh *tmpM;
  
  SomeData   *tmpN;

  for (ic=0; ic<freeNode.n; ic++)
  {
    if (freeNode[ic].toMesh.n > 0)
    {

      // set f2f for current freeNode
    
      for (d=0; d<freeNode[ic].toFree.n; d++)
      {
        for (j=0; j<freeNode[ic].toFree[d].n; j++)
        {
          f2f[freeNode[ic].toFree[d][j].nd] = freeNode[ic].toFree[d][j].invPtr;
    
          if (freeNode[ic].toFree[d][j].invPtr == NULL) prgError(1,fct,"fatal error!!!");
        }
      }
    
      // perform multiplication
    
      nlc = freeNode[ic].dofX->n;
 
      for (d=0; d<freeNode[ic].toMesh.n; d++)
      {
        j = 0;

        while (j < freeNode[ic].toMesh[d].n)
        {
          tmpM = &(freeNode[ic].toMesh[d][j]); 
    
          tmpN = n2bl[tmpM->nd];
    
          if (tmpN == NULL) break;

          i0 = tmpM->nd * ndm;

          dxmdxf = tmpM->xx.x;
          
          for (jr=0; jr<tmpN->f.n; jr++)
          {
            ir = tmpN->f[jr];
          
            tmpF = f2f[ir];
          
            if (tmpF != NULL)
            {
              nlr = freeNode[ir].dofU->n;
          
              dgfdxf = tmpF->dgfdxf.x;
          
              dgfdxl = tmpN->bl[jr]->dgfdxl.x;
          
              for (i=0; i<ndm; i++)

                if (IDX[i0+i] > 0)

                  for (lr=0; lr<nlr; lr++)

                    for (lc=0; lc<nlc; lc++)

                      dgfdxf[lc*nlr+lr] -= dgfdxl[i*nlr+lr] * dxmdxf[lc*ndm+i];
            }
          }
          j++;
        }
      }
    
      // unset f2f
    
      for (d=0; d<freeNode[ic].toFree.n; d++)
      {
        for (j=0; j<freeNode[ic].toFree[d].n; j++)
        {
          f2f[freeNode[ic].toFree[d][j].nd] = NULL;
        }
      }
    }
  }
  ctimDeriv6 += computerTime.stop(fct);

  doneDeriv6 = true;

  return;
}













void FiniteElementBVPWI::calculateDerivatives7(void)
{
  // calculate  d gm / d xm

  char fct[] = "FiniteElementBVPWI::calculateDerivatives7";

  if (!isALE()) return;

  if (maxMeshDepth < 1) return;

  //cout << fct << "\n";

  computerTime.go(fct);

  int i, j, c = 0, d, e, e2, nstu, nenEl, ndfEl, *IX,
      ir, nr,
      ic, nc;

  for (i=0; i<nodeToMesh.n; i++)

    for (j=0; j<nodeToMesh[i].n; j++)

      nodeToMesh[i][j].zeroMtx();

  for (d=0; d<layElem.n; d++)
  {
    for (e2=0; e2<layElem[d].n; e2++)
    {
      e = layElem[d][e2] - 1;

      if (elem[e]->calcMeshDerivatives() != 0) prgError(1,fct,"local stiffness error!");

      nenEl = elem[e]->nen();
      ndfEl = elem[e]->ndf();
      nstu  = ndfEl * nenEl;

      IX    = elem[e]->ix;

      for (ir=0; ir<nenEl; ir++)
      {
        nr = IX[ir] - 1;

        if (nodeHasU[nr])
        {
          for (ic=0; ic<nenEl; ic++)
          {
            nc = IX[ic] - 1;

            if (nodeHasX[nc])
            {     
              calcM2M(ir*ndfEl, ic*ndm, ndfEl, nstu, pmm[c++]->xx.x);
            }
          }
        }
      }
    }
  }
  ctimDeriv7 += computerTime.stop(fct);

  doneDeriv7 = true;

  return;
}









void FiniteElementBVPWI::calcM2M(int ir0, int ic0, int ndfEl, int nstu, double *dgmdxm)
{
  char fct[] = "FiniteElementBVPWI::calcM2M";

  int lr, lc;

  for (lr=0; lr<ndfEl; lr++)
  {
    for (lc=0; lc<ndm; lc++)
    {
      dgmdxm[lc*ndfEl+lr] += s[(ic0+lc)*nstu + ir0+lr];

      //cout << s[(ic0+lc)*nstu + ir0+lr] << "\n";
    }
  }
  return;
}








void FiniteElementBVPWI::calculateDerivatives8(void)
{
  // calculate columns of  d ul / d xf
  // and perform multiplication   d gf / d ul  *  d ul / d xf

  char fct[] = "FiniteElementBVPWI::calculateDerivatives8";

  if (!isALE()) return;

  if (nequ == 0) return;

  //cout << fct << "\n";

  computerTime.go(fct);

  int nlc, nlr, d, i, j, ic, l, lc, ir, lr, rr, nr, n, nXndm, *IDU = idu.x, *IDX;

  double *dgmdxm, *dxmdxf, *dgfdul, *dgfdxf, *dgldxf;

  DataToMesh *tmpM;
  DataToFree *tmpF;
  DataToBLay *tmpBL;

  //if (solver->currentStatus == ASSEMBLY_OK) solver->factorise();

  if (solver->currentStatus != FACTORISE_OK) prgError(1,fct,"Solver status != FACTORISE_OK !");

  for (ic=0; ic<freeNode.n; ic++)
  {
    nlc = freeNode[ic].dofX->n;

    for (lc=0; lc<nlc; lc++)
    {
      // zero the column

      for (i=0; i<nequ; i++) col[i] = 0.;

      // add   d gl / d xf

      for (i=0; i<freeNode[ic].toBLay.n; i++)
      {
        nr = freeNode[ic].toBLay[i].nd - 1;

        dgldxf = freeNode[ic].toBLay[i].dgldxf.x;

        for (lr=0; lr<ndf; lr++)
        {
          rr = IDU[nr*ndf+lr] - 1;
 
          if (rr > -1)  col[rr] = dgldxf[lc*ndf+lr];
        }
      }

      // calculate and add product  d gm / d xm  *  d xm / d xf

      for (d=0; d<freeNode[ic].toMesh.n; d++)
      {
        for (j=0; j<freeNode[ic].toMesh[d].n; j++)
        {
          tmpM   = &(freeNode[ic].toMesh[d][j]);
          n      = tmpM->nd;
          dxmdxf = tmpM->xx.x;

          IDX    = idx.x + n * ndm;

          for (i=0; i<nodeToMesh[n].n; i++)
          {
            nr     = nodeToMesh[n][i].nd;

            dgmdxm = nodeToMesh[n][i].xx.x;

            for (lr=0; lr<ndf; lr++)
            {
              rr = IDU[nr*ndf+lr] - 1;

              if (rr > -1)

                for (l=0; l<ndm; l++)

                  if (IDX[l] > 0) col[rr] -= dgmdxm[l*ndf+lr] * dxmdxf[lc*ndm+l];
            }
          }
        }
      }

      //for (i=0; i<nequ; i++) cout << col[i] << " in 2\n"; cout << "\n";

      computerTime.go("backsubstitution8");

      //cout << "back 8\n";

      solver->solve(col.x);  // -> column now holds  d u / d uf

      ctimBack8 += computerTime.stop("backsubstitution8");

      //for (i=0; i<nequ; i++) cout << col[i] << " out\n"; cout << "\n\n\n";

      // extract  d ul / d xf  from col and use in  d gf / d ul  x  d ul / d xf

      for (d=0; d<freeNode[ic].toFree.n; d++)
      {
        for (j=0; j<freeNode[ic].toFree[d].n; j++)
        {
          tmpF   = &(freeNode[ic].toFree[d][j]);
          dgfdxf = tmpF->invPtr->dgfdxf.x;
          ir     = tmpF->nd;
          nlr    = freeNode[ir].dofU->n;

          for (i=0; i<freeNode[ir].toBLay.n; i++)
          {
            tmpBL  = &(freeNode[ir].toBLay[i]);
            n      = tmpBL->nd - 1;
            dgfdul = tmpBL->dgfdul.x;

            for (l=0; l<ndf; l++)
            {
              rr = IDU[n*ndf+l] - 1;

              if (rr > -1)
              {
                for (lr=0; lr<nlr; lr++) dgfdxf[lc*nlr+lr] -= dgfdul[l*nlr+lr] * col[rr];
              }
            }
          }
        }
      }
    }
  }
  ctimDeriv8 += computerTime.stop(fct);

  doneDeriv8 = true;

  /*//int dp;

  for (i=0; i<freeNode.n; i++)
    for (dp=0; dp<freeNode[i].toFree.n; dp++)
      for (j=0; j<freeNode[i].toFree[dp].n; j++)
      {
        //cout << freeNode[i].toFree[dp][j].dgfduf << " u ives8\n";
        cout << freeNode[i].toFree[dp][j].dgfdxf << " x ives8\n";
      }*/

  return;
}


