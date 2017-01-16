
#include <iostream>

#include "SolverBroyden.h"
#include "Debug.h"
#include "FunctionsSolver.h"
#include "FunctionsProgram.h"
#include "SolverWorkSpace.h"
#include "SolverTime.h"
#include "ComputerTime.h"
#include "MathBasic.h"


extern SolverWorkSpace solverWorkSpace;
extern SolverTime      solverTime;
extern ComputerTime    computerTime;


using namespace std;




SolverBroyden::SolverBroyden(void)
{
  if (debug) cout << " SolverBroyden constructor\n\n";

  mtx = NULL;
  dxn = NULL;
  dFn = NULL;
  Fn  = NULL;

  return;
}




SolverBroyden::~SolverBroyden()
{
  free();

  return;
}





int SolverBroyden::initialise(int n, int, int)
{
  neq = n;

  mtx = new double [neq*neq];

  dxn = new double [neq];

  dFn = new double [neq];

  Fn  = new double [neq];

  int i;

  for (i=0; i<neq*neq; i++) mtx[i] = 0.;

  for (i=0; i<neq; i++)
  {
    mtx[i+i*neq] = 1.;
    dFn[i] = 0.;
    dxn[i] = 0.;
    Fn[i]  = 0.;
  }

  currentStatus = INIT_OK;

  return 0;
}








double *SolverBroyden::factoriseAndSolve(double *Fn1, int)
{
  int i, j;

  double fact;

  // dxn - invJ(n-1) x dFn

  for (i=0; i<neq; i++)
    for (j=0; j<neq; j++)
      dxn[i] -= mtx[j*neq+i] * dFn[j];

  // divide by dFnT x dFn

  fact = 1. / dot(dFn,dFn,neq);

  for (i=0; i<neq; i++) dxn[i] *= fact;

  // multiply with dFnT, add result to invJ

  for (i=0; i<neq; i++)
    for (j=0; j<neq; j++)
      mtx[j*neq+i] += dxn[i] * dFn[j];

  // dx = - invJ x Fn

  for (i=0; i<neq; i++)
  {
    dxn[i] = 0.;
    for (j=0; j<neq; j++)
      dxn[i] -= mtx[j*neq+i] * Fn1[j];
  }
  //for (i=0; i<neq; i++) dxn[i] *= -.005;



  // update dFn

  for (i=0; i<neq; i++)
  {
    dFn[i] = Fn1[i] - Fn[i];
    Fn[i]  = Fn1[i];
  }

  return dxn;
}







void SolverBroyden::free(void)
{
  if (mtx != NULL) delete [] mtx;
  if (dxn != NULL) delete [] dxn;
  if (dFn != NULL) delete [] dFn;
  if (Fn  != NULL) delete [] Fn;

  return;
}
    



