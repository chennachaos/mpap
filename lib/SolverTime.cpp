
#include <iostream>

#include "SolverTime.h"
#include "Definitions.h"
#include "ComputerTime.h"


extern ComputerTime computerTime;



using namespace std;



SolverTime::SolverTime(void)
{
  computerTime.go("SolverTime"); 

  reset();

  return;
}



SolverTime::~SolverTime()
{
  return;
}




void SolverTime::print(void)
{
  double ctimSinceLastCall = computerTime.stop("SolverTime");

  COUT << "linear solvers: computer time statistics\n";
  COUT << "----------------------------------------------------\n";

  COUT; printf("time since last call:        %9.3f sec\n",ctimSinceLastCall);
  COUT << "----------------------------------------------------\n";

  COUT; printf("factorise:                     %7.3f sec ->%5.1f %\n",
                factorise, factorise/ctimSinceLastCall*100.);

  COUT; printf("solve:                         %7.3f sec ->%5.1f %\n",
                solve, solve/ctimSinceLastCall*100.);

  COUT; printf("factoriseAndSolve:             %7.3f sec ->%5.1f %\n",
                factoriseAndSolve, factoriseAndSolve/ctimSinceLastCall*100.);

  COUT; printf("total:                         %7.3f sec ->%5.1f %\n\n",
                total, total/ctimSinceLastCall*100.);

  reset();

  computerTime.go("SolverTime"); 

  return;
}







void SolverTime::reset(void)
{
  factorise         = 0.;
  solve             = 0.;
  factoriseAndSolve = 0.;
  total             = 0.;

  return;
}







