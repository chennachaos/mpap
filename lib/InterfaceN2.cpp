
#include <iostream>

#include "FunctionsProgram.h"
#include "InterfaceN.h"
#include "DomainTree.h"
#include "SolverMA41.h"
//#include "SolverPARDISO.h"
#include "Solid.h"
#include "Fluid.h"
#include "ComputerTime.h"
#include "WhatToSolveForEnum.h"
#include "MpapTime.h"
#include "FreeSurface.h"


extern DomainTree domain;
extern ComputerTime computerTime;
extern MpapTime mpapTime;






void InterfaceN::prepareMatrixPatternSub1(int dm, MatrixSparse<double> &mtx,
                                          int *NIDU, int *PIDU, int *NIDX, int *PIDX,
                                          ListArray< Vector<int> > &i2I,
                                          ListArray< Vector<int> > &i2IPos)
{
  char fct[] = "InterfaceN::prepareMatrixPatternSub1";

  int dp, ir, nc, nr, ndp, j, m, k, o, q;

  for (ir=0; ir<nFN; ir++)
  {
    nr = F2I[ir];

    if (domMode[dm] != RESOLVE) ndp = FN[ir].toFree.n; else ndp = 1;

    for (dp=0; dp<ndp; dp++)
    {
      for (j=0; j<FN[ir].toFree[dp].n; j++)
      {
        nc = F2I[FN[ir].toFree[dp][j].nd];

        m = i2I[nr].n;

        k = 0; while (k < m) if (i2I[nr][k] == nc) break; else k++;

        if (k == m) 
        {
          i2I[nr].append(nc);

          i2IPos[nr].append(mtx.x.n);

          for (o=0; o<NIDU[nr]; o++)
          {
            for (q=0; q<NIDU[nc]; q++) mtx.append(PIDU[nr]+o,PIDU[nc]+q);

            for (q=0; q<NIDX[nc]; q++) mtx.append(PIDU[nr]+o,PIDX[nc]+q);
          }
        }
        FN[ir].toFree[dp][j].pos = i2IPos[nr][k];
      }
    }
  }
  return;
}









void InterfaceN::calcStiffnessAndResidualSub1(int dm, double *fact)
{
  // d gi / d ui

  char fct[] = "InterfaceN::calcStiffnessAndResidualSub1";

  computerTime.go(fct);

  int cc, dp, ic, ir, j, k, ndp, ndfr, pp, rr, *IDR, *IDUC, *IDXC;

  double *dgfduf, *dgfdxf;

  for (ir=0; ir<nFN; ir++)
  {
    IDR = &(idu.x[F2I[ir]*ndf]);

    ndfr = FN[ir].dofU->n;

    for (j=0; j<ndf; j++)
    {
      k = IDR[j] - 1;

      if (k > -1)

        r[k] += fact[0] * FN[ir].reac[j] + fact[1] * FN[ir].reacn[j];
    }

    if (domMode[dm] != RESOLVE) ndp = FN[ir].toFree.n; else ndp = 1;
    
    for (dp=0; dp<ndp; dp++)
    {
      for (j=0; j<FN[ir].toFree[dp].n; j++)
      {
        ic     = FN[ir].toFree[dp][j].nd;
        pp     = FN[ir].toFree[dp][j].pos;
        dgfduf = FN[ir].toFree[dp][j].dgfduf.x;
        dgfdxf = FN[ir].toFree[dp][j].dgfdxf.x;

        IDUC = &(idu.x[F2I[ic]*ndf]);
        IDXC = &(idx.x[F2I[ic]*ndm]);

        for (rr=0; rr<ndf; rr++)
        {
          if (IDR[rr] > 0)
          {        
            for (cc=0; cc<ndf; cc++)

              if (IDUC[cc] > 0) s[pp++] += fact[2] * dgfduf[cc*ndfr+rr];

            for (cc=0; cc<ndm; cc++) 

              if (IDXC[cc] > 0) s[pp++] += fact[2] * dgfdxf[cc*ndfr+rr] * mpapTime.dt;
          }
        }
      }
    }
  }
  ctimCalcSub1 += computerTime.stop(fct);

  doneCalcSub1 = true;

  return;
}












void InterfaceN::prepareMatrixPatternSub2(int dm, MatrixSparse<double> &mtx,
                                          int *NIDU, int *PIDU, int *NIDX, int *PIDX)
{
  char fct[] = "InterfaceN::prepareMatrixPatternSub2";

  if (domMode[dm] == ELIMINATE) return;

  if (domMode[dm] == DUMMY) return;

  if (!domPtr[dm]->solverOK) prgWarning(1,fct,"domain solver not initialised!");

  int i, j, k, l, m, n, o, rr, ndfDom, *IDU;

  ListArray< Vector<int> > tmp;

  tmp.setDim(nequ+neqx);

  ndfDom = domPtr[dm]->ndf;

  IDU = domPtr[dm]->idu.x;

  o = nOff[dm] - 1;

  // d gi / d um,    d gm / d ui

  for (i=0; i<nFN; i++)
  {
    n = F2I[i];

    for (j=0; j<FN[i].toBLay.n; j++)
    {
      m = (FN[i].toBLay[j].nd - 1) * ndfDom;

      for (k=0; k<ndfDom; k++)
      {
        rr = IDU[m+k];

        if (rr > 0)
        {
          for (l=0; l<NIDU[n]; l++)
          {
            tmp[PIDU[n]+l-1].append(rr-1); tmp[PIDU[n]+l-1].append(mtx.x.n);

            mtx.append(rr + o, PIDU[n] + l);

            mtx.append(PIDU[n] + l, rr + o);
          }
          for (l=0; l<NIDX[n]; l++)
          {
            tmp[PIDX[n]+l-1].append(rr-1); tmp[PIDX[n]+l-1].append(mtx.x.n);

            mtx.append(rr + o, PIDX[n] + l);
          }
        }
      }
    }
  }

  colDblPos[dm] = tmp;

  //for (i=0; i<tmp.n; i++) cout << tmp[i] << "\n";

  if (domPtr[dm]->nequ == 0) return;

  // d gm / d um

  domPtr[dm]->adjustElementAssembly(mtx.x.n,DOF);

  mtx.append(nOff[dm],nOff[dm],((SolverSparse*)(domPtr[dm]->solver))->mtx);

  domPtr[dm]->solver->free();

  domPtr[dm]->solver = NULL;

  return;
}










void InterfaceN::calcStiffnessAndResidualSub2(int dm, double *fact)
{
  //   d gm / d um  d gl / d uf,   d gf / d ul

  char fct[] = "InterfaceN::calcStiffnessAndResidualSub2";

  if (domMode[dm] == ELIMINATE) return;

  if (domMode[dm] == DUMMY) return;

  computerTime.go(fct);

  int i, j, k, l, m, n, pos, ndfDom, nDoFU, nDoFX, *IDU, *IDX, *ID;

  double *R = domPtr[dm]->r.x, *dglduf, *dgfdul, *dgldxf;

  n = domPtr[dm]->nequ;

  for (i=0; i<n; i++) r[i+nOff[dm]-1] = R[i];

  ndfDom = domPtr[dm]->ndf;

  ID = domPtr[dm]->idu.x;

  IDU  = idu.x;
  IDX  = idx.x;

  pos = posMtxSub2[dm];

  for (i=0; i<nFN; i++)
  {
    n = F2I[i] * ndf;

    nDoFU = FN[i].dofU->n;
    nDoFX = FN[i].dofX->n;

    for (j=0; j<FN[i].toBLay.n; j++)
    {
      m = (FN[i].toBLay[j].nd - 1) * ndfDom;

      dglduf = FN[i].toBLay[j].dglduf.x;
      dgfdul = FN[i].toBLay[j].dgfdul.x;
      dgldxf = FN[i].toBLay[j].dgldxf.x;

      for (k=0; k<ndfDom; k++)
      {
        if (ID[m+k] > 0)
        {
          for (l=0; l<nDoFU; l++)
          {
            if (IDU[n+l] > 0)
            {
              s[pos++] += fact[3] * dglduf[l*ndfDom+k];

              s[pos++] += fact[0] * dgfdul[k*nDoFU+l];
            }
          }
          for (l=0; l<nDoFX; l++)
          {
            if (IDX[n+l] > 0)

              s[pos++] += dgldxf[l*ndfDom+k] * mpapTime.dt;
          }
        }
      }
    }
  }
  ctimCalcSub2 += computerTime.stop(fct);

  doneCalcSub2 = true;

  return;
}











void InterfaceN::prepareMatrixPatternSub3(int dm, MatrixSparse<double> &mtx,
                                          int *NIDU, int *PIDU, int *NIDX, int *PIDX)
{
  char fct[] = "InterfaceN::prepareMatrixPatternSub3";

  if (domMode[dm] == ELIMINATE) return;

  if (domMode[dm] == DUMMY) return;

  if (domPtr[dm]->neqx < 1 || !domPtr[dm]->isALE()) { colDblPos[dm].free(); return; }

  if (!domPtr[dm]->solverOK) prgWarning(1,fct,"domain solver not initialised!");

  int cc, dp, i, ic, j, k, l, l0, lc, lr, m, n, nlc, ndfDom, ndmDom, nr, o, q, rr, 
      *CDP, *IDU, *IDX, *IDUDom, *IDXDom;

  VectorArray<int> colFlag;

  ndfDom = domPtr[dm]->ndf;
  ndmDom = domPtr[dm]->ndm;

  IDUDom = domPtr[dm]->idu.x;
  IDXDom = domPtr[dm]->idx.x;

  IDU = idu.x;
  IDX = idx.x;

  o = nOff[dm] - 1;

  ListArray< ListArray<DataToMesh> > &N2M = domPtr[dm]->nodeToMesh;

  if (domMode[dm] == ELIMINATE_MESH) // mesh movement eliminated
  {
    colFlag.setDim(domPtr[dm]->nequ);

    for (ic=0; ic<nFN; ic++)
    {
      // (d gm / d xm) * (d xm / d xf)

      n = F2I[ic];
 
      nlc = FN[ic].dofX->n;

      //if (FN[ic].isLagr)
      { 
        for (lc=0; lc<nlc; lc++)
        {
          if (FN[ic].isLagr) cc = IDU[n*ndf+lc];
          else               cc = IDX[n*ndm+lc] + nequ;

          if (cc > 0)
          {
            colFlag.zero();

            CDP = colDblPos[dm][cc-1].x;
            m   = colDblPos[dm][cc-1].n;

            for (j=0; j<m; j+=2) colFlag[CDP[j]] = CDP[j+1];

            for (dp=0; dp<FN[ic].toMesh.n; dp++)
            {
              for (j=0; j<FN[ic].toMesh[dp].n; j++)
              {
                m = FN[ic].toMesh[dp][j].nd;
          
                l0 = m * ndmDom;
          
                for (k=0; k<N2M[m].n; k++)
                {
                  nr = N2M[m][k].nd * ndfDom;
          
                  for (lr=0; lr<ndfDom; lr++)
                  {
                    rr = IDUDom[nr+lr];
          
                    if (rr > 0)

                      if (colFlag[rr-1] < 1)
                      {
                        colFlag[rr-1] = mtx.x.n;
                        mtx.append(rr+o,cc);
                      }
                  }
                }
              }
            }
          }
        }
      }
      //else
      //{
      //}
    }
  }
  else // mesh movement resolved (fully monolithic procedure)
  {
    colDblPos[dm].free();

    q = o + domPtr[dm]->nequ;

    // d gm / d xm

    for (i=0; i<N2M.n; i++)
    {
      m = i * ndmDom;

      for (j=0; j<N2M[i].n; j++)
      {
        n = N2M[i][j].nd * ndfDom;

        for (k=0; k<ndmDom; k++)
        {
          cc = IDXDom[m+k];                

          if (cc > 0)
          {
            for (l=0; l<ndfDom; l++)
            {
              rr = IDUDom[n+l];

              if (rr > 0) mtx.append(rr + o, cc + q);
            }
          }
        }
      }
    }

    // d m / d uf  and   d m / d xf

    for (i=0; i<nFN; i++)
    {
      n = F2I[i];

      if (FN[i].isLagr)
      { 
        for (j=0; j<FN[i].toBLayMesh.n; j++)
        {
          m = (FN[i].toBLayMesh[j].nd - 1) * ndmDom;
 
          for (k=0; k<ndmDom; k++)
          {
            rr = IDXDom[m+k];

            if (rr > 0)  for (l=0; l<NIDU[n]; l++) mtx.append(rr + q, PIDU[n] + l);
          }
        }
      }
      else
      {
        for (j=0; j<FN[i].toBLayMesh.n; j++)
        {
          m = (FN[i].toBLayMesh[j].nd - 1) * ndmDom;
 
          for (k=0; k<ndmDom; k++)
          {
            rr = IDXDom[m+k];

            if (rr > 0)  for (l=0; l<NIDX[n]; l++) mtx.append(rr + q, PIDX[n] + l);
          }
        }
      }
    }

    // d gf / d xm

    for (i=0; i<nFN; i++)
    {
      n = F2I[i];

      for (j=0; j<FN[i].toBLay.n; j++)
      {
        m = (FN[i].toBLay[j].nd - 1) * ndmDom;

        for (k=0; k<ndmDom; k++)
        {
          cc = IDXDom[m+k];

          if (cc > 0)  for (l=0; l<NIDU[n]; l++)  mtx.append(PIDU[n]+l,cc+q);
        }
      }
    }

    // d m / d xm

    domPtr[dm]->adjustElementAssembly(mtx.x.n,MSH);

    mtx.append(q+1,q+1,((SolverSparse*)(domPtr[dm]->solverMesh))->mtx);

    domPtr[dm]->solverMesh->free();
  }
  return;
}











void InterfaceN::calcStiffnessAndResidualSub3(int dm, double *fact)
{
  //   (d gm / d xm) * (d xm / d uf),   (d gf / d xm) * (d xm / d uf)

  char fct[] = "InterfaceN::calcStiffnessAndResidualSub3";

  if (domMode[dm] == ELIMINATE) return;

  if (domMode[dm] == DUMMY) return;

  if (domPtr[dm]->neqx < 1 || !domPtr[dm]->isALE()) return;

  computerTime.go(fct);

  int cc, dp, i, ic, ir, j, k, l, l0, lc, lr, m, n, 
      nDoFU, nDoFX, nlc, nr, pos, ndfDom, ndmDom, rr,
      *CDP, *IDUDom, *IDXDom, *IDU, *IDX;

  double *dxmdxf, *dgmdxm, *dgfdxl, *dmldxf;

  DataToMesh *tmpM;

  ListArray< ListArray<DataToMesh> > &N2M = domPtr[dm]->nodeToMesh;

  VectorArray<int> colFlag;
  VectorArray<double> col;

  ndfDom = domPtr[dm]->ndf;
  ndmDom = domPtr[dm]->ndm;

  IDUDom = domPtr[dm]->idu.x;
  IDXDom = domPtr[dm]->idx.x;

  IDU = idu.x;
  IDX = idx.x;

  pos = posMtxSub3[dm];

  if (domMode[dm] == ELIMINATE_MESH) // mesh movement eliminated
  {
    colFlag.setDim(domPtr[dm]->nequ);
    col.setDim(colFlag.n);

    for (ic=0; ic<nFN; ic++)
    {
      n = F2I[ic];// * ndm;

      nlc = FN[ic].dofX->n;

      //if (FN[ic].isLagr)
      {
        for (lc=0; lc<nlc; lc++)
        {
          if (FN[ic].isLagr) cc = IDU[n*ndf+lc];
          else               cc = IDX[n*ndm+lc] + nequ;
          //cc = IDU[n+lc];

          if (cc > 0)
          {
            colFlag.zero(); 
            col.zero();

            CDP = colDblPos[dm][cc-1].x;
            m   = colDblPos[dm][cc-1].n;

            for (j=0; j<m; j+=2) colFlag[CDP[j]] = CDP[j+1];

            for (dp=0; dp<FN[ic].toMesh.n; dp++)
            {
              for (j=0; j<FN[ic].toMesh[dp].n; j++)
              {
                tmpM   = &(FN[ic].toMesh[dp][j]);
                m      = tmpM->nd;
                dxmdxf = tmpM->xx.x;
                l0     = m * ndmDom;            
      
                for (k=0; k<N2M[m].n; k++)
                {
                  nr = N2M[m][k].nd * ndfDom;
              
                  dgmdxm = N2M[m][k].xx.x;
              
                  for (lr=0; lr<ndfDom; lr++)
                  {
                    rr = IDUDom[nr+lr];

                    if (rr > 0)
                    { 
                      if (colFlag[rr-1] < 1) colFlag[rr-1] = pos++;

                      for (l=0; l<ndm; l++)
      
                        if (IDXDom[l0+l] > 0)
      
                          col[rr-1] -= dgmdxm[l*ndfDom+lr] * dxmdxf[lc*ndm+l];
                    }
                  }
                }
              }
            }
            if (FN[ic].isLagr)
              for (rr=0; rr<colFlag.n; rr++)
              {
                j = colFlag[rr];

                if (j > 0) s[j] += col[rr];
              }
            else
              for (rr=0; rr<colFlag.n; rr++)
              {
                j = colFlag[rr];

                if (j > 0) s[j] += col[rr] * mpapTime.dt;
              }
          }
        }
      }
    }
  }
  else // mesh movement resolved (fully monolithic procedure)
  {
    // d gm / d xm

    for (ir=0; ir<N2M.n; ir++)
    {
      m = ir * ndmDom;

      for (j=0; j<N2M[ir].n; j++)
      {
        nr = N2M[ir][j].nd * ndfDom;

        dgmdxm = N2M[ir][j].xx.x;

        for (k=0; k<ndmDom; k++)

          if (IDXDom[m+k] > 0)

            for (l=0; l<ndfDom; l++)

              if (IDUDom[nr+l] > 0) s[pos++] += dgmdxm[k*ndfDom+l] * fact[4];
      }
    }

    // d m / d ui

    for (i=0; i<nFN; i++)
    {
      n = F2I[i] * ndf;
  
      nDoFX = FN[i].dofX->n;

      if (FN[i].isLagr)
      {
        for (j=0; j<FN[i].toBLayMesh.n; j++)
        {
          m = (FN[i].toBLayMesh[j].nd - 1) * ndmDom;
      
          dmldxf = FN[i].toBLayMesh[j].dmldxf.x;
      
          for (k=0; k<ndmDom; k++)
          {
            if (IDXDom[m+k] > 0)
      
              for (l=0; l<nDoFX; l++)
      
                if (IDU[n+l] > 0) s[pos++] += fact[3] * dmldxf[l*ndmDom+k];
          }
        }
      }
      else
      {
        for (j=0; j<FN[i].toBLayMesh.n; j++)
        {
          m = (FN[i].toBLayMesh[j].nd - 1) * ndmDom;
      
          dmldxf = FN[i].toBLayMesh[j].dmldxf.x;
      
          for (k=0; k<ndmDom; k++)
          {
            if (IDXDom[m+k] > 0)
      
              for (l=0; l<nDoFX; l++)
      
                if (IDX[n+l] > 0) s[pos++] += dmldxf[l*ndmDom+k] * mpapTime.dt;
          }
        }
      }
    }

    // d gi / d xm

    for (i=0; i<nFN; i++)
    {
      n = F2I[i] * ndf;

      nDoFU = FN[i].dofU->n;

      for (j=0; j<FN[i].toBLay.n; j++)
      {
        m = (FN[i].toBLay[j].nd - 1) * ndmDom;
    
        dgfdxl = FN[i].toBLay[j].dgfdxl.x;
    
        for (k=0; k<ndmDom; k++)
        {
          if (IDXDom[m+k] > 0)
          {
            for (l=0; l<nDoFU; l++)
            {
              if (IDU[n+l] > 0)
              {
                s[pos++] += fact[5] * dgfdxl[k*nDoFU+l];
              }
            }
          }
        }
      }
    }

    // d m / d xm

    while (pos < posMtxSub4[dm]) s[pos++] *= fact[4];
  }
  ctimCalcSub3 += computerTime.stop(fct);

  doneCalcSub3 = true;

  return;
}











void InterfaceN::prepareMatrixPatternFreeSurface(int dm, MatrixSparse<double> &mtx,
                                                 int *NIDU, int *PIDU, int *NIDX, int *PIDX)
{
  char fct[] = "InterfaceN::prepareMatrixPatternFreeSurface";

  if (!isFreeSurface(*domPtr[dm])) return;

  int cc, e, i, ir, ic, jr, jc, kc, kr, nr, nc, numelDom, numnpDom, nenDom, nstDom, pos,
      *IX, *IDU, *IDX, *SPRS;

  VectorArray<int> *nodeNodeDom = domPtr[dm]->nodeNode;

  Element **elem = domPtr[dm]->elem;

  ListArray< Vector<int> > tmp;

  numnpDom = domPtr[dm]->numnp;
  numelDom = domPtr[dm]->numel;
  nenDom   = domPtr[dm]->nen;
  nstDom   = nenDom * ndm;

  tmp.setDim(numnpDom);

  IDU = idu.x;
  IDX = idx.x;

  // allocate nonzero coefficients in global system matrix

  for (ir=0; ir<numnpDom; ir++)
  {
    nr = F2I[ir];

    for (ic=0; ic<nodeNodeDom[ir].n+1; ic++)
    {
      if (ic < nodeNodeDom[ir].n) nc = F2I[nodeNodeDom[ir][ic]-1]; else nc = nr;

      tmp[ir].append(mtx.x.n);

      if (NIDX[nr] > 0) if (NIDX[nr] != ndm) prgError(1,fct,"NIDX[nr] != ndm");

      for (jr=0; jr<NIDX[nr]; jr++)
      {
        for (jc=0; jc<ndm; jc++)
        {
          cc = IDU[nc*ndf+jc];

          if (cc > 0)

          mtx.append(PIDX[nr]+jr,cc);
        }
        for (jc=0; jc<NIDX[nc]; jc++)

          mtx.append(PIDX[nr]+jr,PIDX[nc]+jc);
      }
    }
  }

  // set element pointers for stiffness assembly

  for (e=0; e<numelDom; e++)
  {
    IX = elem[e]->ix;

    if (elem[e]->sparse[0] != NULL) prgError(2,fct,"why is this the case?!");

    elem[e]->sparse[0] = new int [nstDom * nstDom * 2];

    SPRS = elem[e]->sparse[0];

    for (i=0; i<nstDom*nstDom*2; i++) SPRS[i] = -1;

    for (kr=0; kr<nenDom; kr++)
    {
      ir = IX[kr] - 1;
      nr = F2I[ir];

      for (kc=0; kc<nenDom; kc++)
      {
        ic = nodeNodeDom[ir].n;

        if (kr != kc) while (--ic > -1) if (nodeNodeDom[ir][ic] == IX[kc]) break;

        if (ic == -1) prgError(1,fct,"fatal error!!!");
        
        pos = tmp[ir][ic];

        nc = F2I[IX[kc]-1];

        for (jr=0; jr<ndm; jr++)
        {
          if (IDX[nr*ndm+jr] > 0)
          {
            for (jc=0; jc<ndm; jc++)

              if (IDU[nc*ndf+jc] > 0)  

                SPRS[(kc*ndm+jc)*nstDom+(kr*ndm+jr)] = pos++;

            for (jc=0; jc<ndm; jc++)

              if (IDX[nc*ndm+jc] > 0)

                SPRS[(nstDom+kc*ndm+jc)*nstDom+(kr*ndm+jr)] = pos++;
          }
        }
      }
    }
  }

  return;
}









void InterfaceN::calcStiffnessAndResidualFreeSurface(int dm, double *fact)
{
  char fct[] = "InterfaceN::calcStiffnessAndResidualFreeSurface";

  computerTime.go(fct);

  int e, i, j, jc, ic, rr, numelDom, nstDom, nenDom, *SPRS, *IX, *IDX, n1, n2;

  double *sDom, *pDom;

  Element **elem = domPtr[dm]->elem;

  ListArray< Vector<int> > tmp;

  numelDom = domPtr[dm]->numel;
  nenDom   = domPtr[dm]->nen;
  nstDom   = nenDom * ndm;
  n1       = nstDom * nstDom;
  n2       = n1 + n1;
  IDX      = idx.x;
  sDom     = domPtr[dm]->s;
  pDom     = domPtr[dm]->p;

  for (e=0; e<numelDom; e++)
  {
    elem[e]->calcStiffnessAndResidual();

    SPRS = elem[e]->sparse[0];

    IX = elem[e]->ix;

    for (i=0; i<nenDom; i++)
    {
      for (j=0; j<ndm; j++)
      {
        rr = IDX[(IX[i]-1)*ndm+j] - 1;

        if (rr > -1) r[nequ+rr] += pDom[i*ndm+j]; 
      }

      if (isLagr[F2I[IX[i]-1]])
      {
        for (jc=0; jc<ndm; jc++)

          for (rr=0; rr<nstDom; rr++)

            sDom[(i*ndm+jc)*nstDom+rr] += sDom[(i*ndm+jc+nstDom)*nstDom+rr] * dxdu;
      }
    }
    for (i= 0; i<n1; i++) if (SPRS[i] > -1) s[SPRS[i]] += sDom[i];
    for (i=n1; i<n2; i++) if (SPRS[i] > -1) s[SPRS[i]] += sDom[i] * mpapTime.dt;
  }

  ctimCalcFreeSurf += computerTime.stop(fct);

  doneCalcFreeSurf = true;

  return;
}





