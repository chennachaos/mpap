

#include "FreeNode.h"
#include "FunctionsProgram.h"


using namespace std;



FreeNode::FreeNode(void)
{
  idu = NULL;

  return;
}




FreeNode::FreeNode(int t, double *X, int fD, int mD, VectorArray<int> *dofUPtr, 
                                                     VectorArray<int> *dofXPtr)
{ 
  typ = t;

  isLagr = false;

  idu = NULL;

  dofU = dofUPtr;
  dofX = dofXPtr;
 
  int i, ndm = dofX->n, ndf = dofU->n;
  
   reac.setDim(ndf);  reac.zero();
  reacn.setDim(ndf); reacn.zero();
      u.setDim(ndf);     u.zero();
      x.setDim(ndm);

  for (i=0; i<ndm; i++) x.x[i] = X[(*dofX)[i]];

  freeDepth = fD;
  meshDepth = mD;

  return;
}




FreeNode::~FreeNode()
{
  if (idu != NULL) delete [] idu;

  return;
}





FreeNode &FreeNode::operator=(FreeNode &fN)
{
  dofU      = fN.dofU;
  dofX      = fN.dofX;

  nd        = fN.nd;
  typ       = fN.typ;
  isLagr    = fN.isLagr;
  freeDepth = fN.freeDepth;
  meshDepth = fN.meshDepth;

  u         = fN.u;
  x         = fN.x;
  reac      = fN.reac;
  reacn     = fN.reacn;

  toFree    = fN.toFree;
  toMesh    = fN.toMesh;
  toBLay    = fN.toBLay;

  return *this;
}







void FreeNode::giveUX(double *uNew, double *xNew, double *DAT, double *U, double *X, double *X0)
{
  char fct[] = "FreeNode::giveUX";

  int i;

  switch (typ)
  {
    case INTERPOLATION_FN:   for (i=0; i<dofU->n; i++) uNew[dofU->x[i]] = u[i] * (*DAT);

                             for (i=0; i<dofX->n; i++) xNew[dofX->x[i]] = x[i] * (*DAT);

                             break;

    default:                 prgError(1,fct,"not yet implemented for this type of free node!");
  }
  return;
}







void FreeNode::giveU(double *uNew, double *DAT, double *U, double*X, double *X0)
{
  char fct[] = "FreeNode::giveU";

  int i;

  switch (typ)
  {
    case INTERPOLATION_FN:   for (i=0; i<dofU->n; i++) uNew[dofU->x[i]] = u[i] * (*DAT);

                             break;

    default:                 prgError(1,fct,"not yet implemented for this type of free node!");
  }

  return;
}






void FreeNode::giveX(double *xNew, double *DAT, double *U, double *X, double *X0)
{
  char fct[] = "FreeNode::giveX";

  int i;

  switch (typ)
  {
    case INTERPOLATION_FN:   for (i=0; i<dofX->n; i++) xNew[dofX->x[i]] = x[i] * (*DAT);

                             break;

    default:                 prgError(1,fct,"not yet implemented for this type of free node!");
  }

  return;
}






void FreeNode::getReactions(double *REAC, double *DAT, double *U, double *X, double *X0)
{
  char fct[] = "FreeNode::getReactions";

  int i;

  switch (typ)
  {
    case INTERPOLATION_FN:   for (i=0; i<dofU->n; i++) reac[i] += REAC[dofU->x[i]] * (*DAT);

                             break;

    default:                 prgError(1,fct,"not yet implemented for this type of free node!");
  }

  return;
}






void FreeNode::print(void)
{
  printf(" (%6.3f",x.x[0]);

  for (int i=1; i<x.n; i++) printf(",%6.3f",x.x[i]);

  cout << ")  ";

  return;
}





void FreeNode::zeroDerivativeMtx(void)
{
  int d, j;

  for (d=0; d<toFree.n; d++) for (j=0; j<toFree[d].n; j++) toFree[d][j].zeroMtx();

  for (j=0; j<toBLay.n; j++) toBLay[j].zeroMtx();

  for (j=0; j<toBLayMesh.n; j++) toBLayMesh[j].zeroMtx();

  return;
}

