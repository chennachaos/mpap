

#include "BndNode.h"
#include "FunctionsProgram.h"


using namespace std;



BndNode::BndNode(void)
{

  return;
}





BndNode::BndNode(int n, double *X, Vector<void*> &fU, Vector<void*> &fX, Vector<double> &d)
{
  int i;

  nd = n;

  ndat = 1;

  x = X;

  freeU.setDim(fU.n);
  for (i=0; i<fU.n; i++) freeU[i] = (FreeNode*)(fU[i]);

  if (fX.n == 0) xFlag = NOX;

  else if (fX == fU)
  {
    xFlag   = XEQU; 
    freeX.n = freeU.n;
    freeX.x = freeU.x;
  }
  else
  {
    xFlag = INDX;
    freeX.setDim(fX.n);
    for (i=0; i<fX.n; i++) freeX[i] = (FreeNode*)(fX[i]);
  }

  ndat = 1;
  dat.setDim(d.n);
  for (i=0; i<d.n; i++) dat[i] = d[i];

  return;
}






BndNode::~BndNode()
{
  if (xFlag == XEQU)
  {
    freeX.n = 0;
    freeX.x = NULL;
  }
  return;
}






BndNode &BndNode::operator=(BndNode &bN)
{
  x     = bN.x;
  nd    = bN.nd;
  ndat  = bN.ndat;
  xFlag = bN.xFlag; 

  dat   = bN.dat;

  freeU = bN.freeU;
  freeX = bN.freeX;

  uDep  = bN.uDep;
  xDep  = bN.xDep;

  return *this;
}






void BndNode::getUX(int ndf, int ndm, double *U, double *X, double *Xn, double *X0)
{
  char fct[] = "BndNode::getUX";

  int i, j, k, n;

  if (uDep.n > 0) n = uDep[0]->nd-1; else n = xDep[0]->nd-1;

  double uNew[MAX_DOF], xNew[MAX_DOF];

  if (xFlag == XEQU)
  {
    for (j=0; j<xDep.n; j++) { xDep[j]->uc = 0.; xDep[j]->duc = 0.; }
    for (j=0; j<uDep.n; j++) { uDep[j]->uc = 0.; uDep[j]->duc = 0.; }

    for (i=0; i<freeU.n; i++)
    {
      freeU[i]->giveUX(uNew, xNew, dat.x+i*ndat, U+n*ndf, X+n*ndm, X0+n*ndm);

      for (j=0; j<xDep.n; j++)
      {
        k = xDep[j]->dof - 1;

        xDep[j]->uc  += xNew[k];

        xDep[j]->duc += xNew[k];
      }

      for (j=0; j<uDep.n; j++)
      {
        k = uDep[j]->dof - 1;

        uDep[j]->uc  += uNew[k];

        uDep[j]->duc += uNew[k];
      }
    }

    for (j=0; j<xDep.n; j++)
    {
      k = xDep[j]->dof - 1;

      xDep[j]->uc  -= Xn[n*ndm+k];

      xDep[j]->duc -= X[n*ndm+k];
    }

    for (j=0; j<uDep.n; j++)
    {
      k = uDep[j]->dof - 1;

      uDep[j]->duc -= U[n*ndf+k];
    }
  }
  else
  {
    getX(ndf,ndm,U,X,Xn,X0);

    getU(ndf,ndm,U,X,Xn,X0);
  }

  return;
}







void BndNode::getX(int ndf, int ndm, double *U, double *X, double *Xn, double *X0)
{
  if (xDep.n == 0) return;

  int i, j, k, n = xDep[0]->nd - 1;

  double xNew[MAX_DOF];

  for (j=0; j<xDep.n; j++) { xDep[j]->uc = 0.; xDep[j]->duc = 0.; }

  for (i=0; i<freeX.n; i++)
  {
    freeX[i]->giveX(xNew, dat.x+i*ndat, U+n*ndf, X+n*ndm, X0+n*ndm);

    for (j=0; j<xDep.n; j++)
    {
      k = xDep[j]->dof - 1;

      xDep[j]->uc  += xNew[k];

      xDep[j]->duc += xNew[k];
    }
  }

  for (j=0; j<xDep.n; j++)
  {
    k = xDep[j]->dof - 1;

    xDep[j]->uc  -= Xn[n*ndm+k];

    xDep[j]->duc -= X[n*ndm+k];
  }

  return;
}







void BndNode::getU(int ndf, int ndm, double *U, double *X, double *Xn, double *X0)
{
  if (uDep.n == 0) return;

  int i, j, k, n = uDep[0]->nd - 1;

  double uNew[MAX_DOF];

  for (j=0; j<uDep.n; j++) { uDep[j]->uc = 0.; uDep[j]->duc = 0.; }

  for (i=0; i<freeU.n; i++)
  {
    freeU[i]->giveU(uNew, dat.x+i*ndat, U+n*ndf, X+n*ndm, X0+n*ndm);

    for (j=0; j<uDep.n; j++)
    {
      k = uDep[j]->dof - 1;

      uDep[j]->uc  += uNew[k];

      uDep[j]->duc += uNew[k];
    }
  }

  for (j=0; j<uDep.n; j++)
  {
    k = uDep[j]->dof - 1;

    uDep[j]->duc -= U[n*ndf+k];
  }

  return;
}







void BndNode::giveReactions(int ndf, int ndm, double *REAC, double *U, double *X, double *X0)
{
  if (uDep.n == 0) return;

  int i, n = uDep[0]->nd - 1;

  for (i=0; i<freeU.n; i++)
  {
    freeU[i]->getReactions(REAC+n*ndf, dat.x+i*ndat, U+n*ndf, X+n*ndm, X0+n*ndm);
  }
  return;
}







/*
void BndNode::calcDerivatives()
{
  int i, n = uDep[0]->nd;

  for (i=0; i<freeU.n; i++)

    freeU->dRdR(dRdR, REAC+n*ndf, dat.x+i*ndat, U+n*ndf, X+n*ndm, X0+n*ndm);

    freeU->dUdU(dUdU, );

    freeU->

  return;
}





void BndNode::getUXDeriv(int ndf, int ndm, double *U, double *X, double *Xn, double *X0)
{
  return;
}


void BndNode::getUDeriv(int ndf, int ndm, double *U, double *X, double *Xn, double *X0)
{



  return;
}


void BndNode::getXDeriv(int ndf, int ndm, double *U, double *X, double *Xn, double *X0)
{
  return;
}


void giveReactionsDeriv(int ndf, int ndm, double *REAC, double *U, double *X, double *X0)
{
  return;
}
*/





bool BndNode::checkOK(void)
{
  char fct[] = "BndNode::checkOK";

  int i, j, m;

  FreeNode **freeUTmp;
  double   *datTmp;

  // check that all free nodes in freeU are of the same type and share the same dofU

  for (i=1; i<freeU.n; i++)
  {
    if (freeU[i]->typ != freeU[0]->typ) prgError(1,fct,"inconsistent freeNode types!");
    
    if (freeU[i]->dofU != freeU[0]->dofU) prgError(1,fct,"inconsistent freeNode dofU!");
  }

  // check that all free nodes in freeX are of the same type and share the same dofX

  for (i=1; i<freeX.n; i++)
  {
    if (freeX[i]->typ != freeX[0]->typ) prgError(2,fct,"inconsistent freeNode types!");
    
    if (freeX[i]->dofX != freeX[0]->dofX) prgError(2,fct,"inconsistent freeNode dofX!");
  }

  // check that free nodes in freeX and freeU are of the same type

  if (freeU.n > 0 && freeX.n > 0) 

    if (freeU[0]->typ != freeX[0]->typ) prgError(3,fct,"different types of freeX and freeU!");

  // check that there is only ONE RIGIDBODY_FN in freeX and freeU

  if (freeU[0]->typ == RIGIDBODY_FN)

    if (freeU.n > 1 || freeX.n > 1) prgError(4,fct,"more than one RIGIDBODY_FN!");

  // check interpolation data structure and remove zero interpolation nodes

  if (freeU[0]->typ == INTERPOLATION_FN)
  {
    //if (xFlag != XEQU) prgError(10,fct,"INTERPOLATION_FN with freeU != freeX !?");

    if (ndat != 1) prgError(11,fct,"INTERPOLATION_FN with ndat != 1 !?");

    if (dat.n != freeU.n) prgError(12,fct,"INTERPOLATION_FN with dat.n != freeU.n !?");

    m = freeU.n;
    i = 0;

    while (i < m)
    {
      if (abs(dat[i]) < 1.e-8)
        { m--; for (j=i; j<m; j++) { freeU[j] = freeU[j+1]; dat[j] = dat[j+1]; } }
      else i++;
    }
    if (m < freeU.n)
    {
      freeUTmp = new FreeNode* [m];
      for (i=0; i<m; i++) freeUTmp[i] = freeU[i];
      delete [] freeU.x;
      freeU.x = freeUTmp;
      freeU.n = m;

      datTmp = new double [m];
      for (i=0; i<m; i++) datTmp[i] = dat[i];
      delete [] dat.x;
      dat.x = datTmp;
      dat.n = m;
    }
    if (xFlag == XEQU)
    {
      freeX.n = freeU.n;
      freeX.x = freeU.x;
    }
  }

  return true;
}









void BndNode::print(int ndm)
{
  int i;

  printf(" %5d (%6.3f",nd,x[0]); for (i=1; i<ndm; i++) printf(",%6.3f",x[i]); cout << ") : ";

  for (i=0; i<freeU.n; i++) { printf("%6.3f",dat[i]); freeU[i]->print(); }

  if (xFlag == INDX) { cout << "; "; for (i=0; i<freeX.n; i++) freeX[i]->print(); } cout << "\n";

  return;
}







void BndNode::generateDeps(int ndm, bool isALE)
{
  int j, n, *DOF, typ = UNDEF_FN;

  if      (freeU.n > 0) typ = freeU[0]->typ;

  else if (freeX.n > 0) typ = freeX[0]->typ;

  else prgError(1,"BndNode::generateDeps","fatal error!");

  if (typ == UNDEF_FN) prgError(2,"BndNode::generateDeps","fatal error!");

  if (typ == INTERPOLATION_FN)
  {
    uDep.setDim(freeU[0]->dofU->n);

    DOF = freeU[0]->dofU->x;

    for (j=0; j<uDep.n; j++) 
    {
      uDep[j]        = new DependentDoF();
      uDep[j]->nd    = nd;
      uDep[j]->dof   = DOF[j]+1;
      uDep[j]->uc    = 0.;
      uDep[j]->tmFct = 0;
    }

    if (isALE && xFlag != NOX) 
    {
      if (xFlag == INDX) { xDep.setDim(freeX[0]->dofX->n); DOF = freeX[0]->dofX->x; }

      else               { xDep.setDim(freeU[0]->dofX->n); DOF = freeU[0]->dofX->x; }

      for (j=0; j<xDep.n; j++)
      {
        xDep[j]        = new DependentDoF();
        xDep[j]->nd    = nd;
        xDep[j]->dof   = DOF[j]+1;
        xDep[j]->uc    = 0.;
        xDep[j]->tmFct = 0; 
      }
    }

    return;
  }

  if (typ == RIGIDBODY_FN)
  {
    if (xFlag == INDX) { n = freeX[0]->dofX->n; DOF = freeX[0]->dofX->x; }

    else               { n = freeU[0]->dofX->n; DOF = freeU[0]->dofX->x; }

    if (DOF[n-1] < ndm)
    {
      xDep.setDim(n);

      for (j=0; j<n; j++) 
      {
        xDep[j]        = new DependentDoF();
        xDep[j]->nd    = nd;
        xDep[j]->dof   = DOF[j]+1;
        xDep[j]->uc    = 0.;
        xDep[j]->tmFct = 0;
      }

      if (freeU.n > 0)
      {
        uDep.setDim(n);

        for (j=0; j<n; j++) 
        {
          uDep[j]        = new DependentDoF();
          uDep[j]->nd    = nd;
          uDep[j]->dof   = DOF[j]+1;
          uDep[j]->uc    = 0.;
          uDep[j]->tmFct = 0;
        }
      }
    }
    else 
    {
      xDep.setDim(ndm);

      for (j=0; j<ndm; j++)
      {
        xDep[j]        = new DependentDoF();
        xDep[j]->nd    = nd;
        xDep[j]->dof   = j + 1;
        xDep[j]->uc    = 0.;
        xDep[j]->tmFct = 0;
      }

      if (freeU.n > 0) 
      {
        uDep.setDim(ndm);

        for (j=0; j<ndm; j++)
        {
          uDep[j]        = new DependentDoF();
          uDep[j]->nd    = nd;
          uDep[j]->dof   = j + 1;
          uDep[j]->uc    = 0.;
          uDep[j]->tmFct = 0;
        }
      }
    }
    return;
  }

  return;
}





int BndNode::giveDerivN(int j)
{
  char fct[] = "BndNode::giveDerivN";

  if (freeU[j]->typ == INTERPOLATION_FN) return 1;
 
  else prgError(1,fct,"not yet implemented!");

  return 0;
}





double *BndNode::giveDeriv(int j)
{
  char fct[] = "BndNode::giveDeriv";

  if (freeU[j]->typ == INTERPOLATION_FN) return dat.x+j;
 
  else prgError(1,fct,"not yet implemented!");

  return 0;
}




int BndNode::getDerivN(int j)
{
  char fct[] = "BndNode::getDerivN";

  if (freeU[j]->typ == INTERPOLATION_FN) return 1;
 
  else prgError(1,fct,"not yet implemented!");

  return 0;
}





double *BndNode::getDeriv(int j)
{
  char fct[] = "BndNode::getDeriv";

  if (freeU[j]->typ == INTERPOLATION_FN) return dat.x+j;
 
  else prgError(1,fct,"not yet implemented!");

  return 0;
}



