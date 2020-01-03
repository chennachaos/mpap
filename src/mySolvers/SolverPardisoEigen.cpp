
#include <iostream>

#include "SolverPardisoEigen.h"
#include "FunctionsSolver.h"
#include "FunctionsProgram.h"
#include "SolverTime.h"
#include "ComputerTime.h"
#include "util.h"

/*
// boost library 
#include <config.hpp>
#include <graph/adjacency_list.hpp>
#include <graph/cuthill_mckee_ordering.hpp>
#include <graph/properties.hpp>
#include <graph/bandwidth.hpp>

using namespace boost;

//typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;

typedef adjacency_list<vecS, vecS, bidirectionalS, 
     property<vertex_color_t, default_color_type,
       property<vertex_degree_t,int> > > Graph;

typedef graph_traits<Graph>::vertex_descriptor Vertex;

typedef graph_traits<Graph>::vertices_size_type size_type;

*/

extern SolverTime      solverTime;
extern ComputerTime    computerTime;
extern int             countThis;

using namespace std;


SolverPardisoEigen::SolverPardisoEigen()
{
  return;
}




SolverPardisoEigen::~SolverPardisoEigen()
{
  phase = -1; error = 0;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

}





int SolverPardisoEigen::initialise(int numProc, int matrixType, int nr)
{
  char fct[] = "SolverPARDISO::initialise";

  nRow = nCol = nr;

  soln.resize(nRow);
  soln.setZero();

  rhsVec   = soln;
  solnPrev = soln;

  //if (currentStatus != PATTERN_OK)
    //{ prgWarning(1,fct,"prepare matrix pattern first!"); return 1; }

  phase = 11; error = 0;
  int  mtxType = 11, idum;


  char *tmp;

  switch (matrixType)
  {
    case PARDISO_STRUCT_SYM: mtxType =  1; break; // real and structurally symmetric

    case PARDISO_UNSYM:      mtxType = 11; break; // real and unsymmetric

    default:                 prgWarning(1,fct,"matrix assumed to be unsymmetric!");
  }
/*
  IPARM[3] = 31;
//  IPARM[4] = 0;  // user input permutation
//  IPARM[4] = 2;  // return the permutation
//  IPARM[5] = 1;  // overwrite RHS with solution
//  IPARM[7] = 1;  // max number of iterative refinement steps
  IPARM[9] = 4;
  IPARM[18] = 1;
  //IPARM[27] = 0; // Parallel METIS ordering
  IPARM[31] = 0; // IPARM(32) = 0 -> sparse direct solver (default) 
  //IPARM[31] = 1; // IPARM(32)=1 -> iterative method  
  //IPARM[9] = 1;
*/

  SOLVER = 0;       // sparse direct solver
  //SOLVER = 1;       // multi-recursive iterative solver
  MTYPE  = mtxType; // matrix type
  MAXFCT = 1;       // maximum number of factorisations of same sparsity pattern
                    //      to be kept in memory
  MNUM   = 1;       // which factorisation to use
  NRHS   = 1;       // number of right hand sides
  MSGLVL = 0;       // output message level (1 -> print statistical information)

  IPARM[0] = 0;     // PARADISO will set IPARM to default values
  //IPARM[0] = 1;     // user input values to IPARM
  //IPARM[1] = 2;

  IPARM[2] = max(1,numProc);  // number of processors (no default value available)
/*
  tmp = getenv("OMP_NUM_THREADS");

  if (tmp != NULL) 
  {
    sscanf(tmp,"%d", &idum);
    if (idum != IPARM[2]) prgError(1,fct,"set environment variable OMP_NUM_THREADS to numProc!");
  }
  else prgError(2,fct,"set environment variable OMP_NUM_THREADS!");
*/

  pardisoinit_(PT, &MTYPE, &SOLVER, IPARM, DPARM, &error);

  if (error != 0)
  {
    if (error == -10) prgError(1,fct,"no license file found.");
    if (error == -11) prgError(2,fct,"license is expired.");
    if (error == -12) prgError(3,fct,"wrong username or hostname.");
  }

  int  *c1, *c2, ii;

  csr.resize(nRow+1);
  col.resize(mtx.nonZeros());

  c1 = mtx.outerIndexPtr();
  c2 = mtx.innerIndexPtr();


  for(ii=0;ii<=nRow;ii++)
    csr[ii] = c1[ii] + 1;

  for(ii=0;ii<mtx.nonZeros();ii++)
    col[ii] = c2[ii] + 1;

  array = mtx.valuePtr();

  perm.resize(nRow);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  if (error != 0)
  {
    COUT << "PARDISO ERROR = " << error << "\n\n";
    prgError(4,fct,"symbolic factorisation failed.");
  }

  //IPARM[4] = 0;  // user input permutation
  //IPARM[4] = 2;  // return the permutation
  IPARM[5] = 0; // do not overwrite RHS with solution
  IPARM[7] = 1; // max number of iterative refinement steps

  currentStatus = INIT_OK;

  return 0;
}





int SolverPardisoEigen::factorise()
{
  char fct[] = "SolverPARDISO::factorise";

  if (currentStatus != ASSEMBLY_OK) 
    { prgWarning(1,fct,"assemble matrix first!"); return 1; }

  phase = 22; error = 0;

  computerTime.go(fct);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  solverTime.total     -= solverTime.factorise;
  solverTime.factorise += computerTime.stop(fct);
  solverTime.total     += solverTime.factorise;

  currentStatus = FACTORISE_OK;
  
  return 0;
}






int SolverPardisoEigen::solve()
{ 
  char fct[] = "SolverPARDISO::solve";

  if (currentStatus != FACTORISE_OK)
    { prgWarning(1,fct,"factorise matrix first!"); return 1; }

  phase = 33; error = 0;

  computerTime.go(fct);
  soln.setZero();

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &rhsVec[0], &soln[0], &error, DPARM);

  solverTime.total -= solverTime.solve;
  solverTime.solve += computerTime.stop(fct);
  solverTime.total += solverTime.solve;

  return 0;
}




int SolverPardisoEigen::factoriseAndSolve()
{
  char fct[] = "SolverPARDISO::factoriseAndSolve";

  if (currentStatus != ASSEMBLY_OK) 
    { prgWarning(1,fct,"assemble matrix first!"); return 1; }

  phase = 23; error = 0;

  computerTime.go(fct);
  soln.setZero();

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &rhsVec[0], &soln[0], &error, DPARM);

  //printf("Peak memory [kB] phase 1       = %d \n", IPARM[14]);
  //printf("Permanent integer memory [kb]. = %d \n", IPARM[15]);
  //printf("Peak real memory [kB]          = %d \n", IPARM[16]);
  //printf("Number of nonzeros in LU.      = %d \n", IPARM[17]);
  //printf("Gflops for LU factorization.   = %d \n", IPARM[18]);


  solverTime.total             -= solverTime.factoriseAndSolve;
  solverTime.factoriseAndSolve += computerTime.stop(fct);
  solverTime.total             += solverTime.factoriseAndSolve;

  currentStatus = FACTORISE_OK;

  return 0;
}



int SolverPardisoEigen::free()
{
  phase = -1; error = 0;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  array = NULL;

  currentStatus = EMPTY;

  return 0;
}




/*
   Graph G(nRow);

   for(int k=0; k<mtx.outerSize(); ++k)
   for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
     add_edge(it.row(), it.col(), G);

    cout << " hhhhhhhhhhhhhhhhh " << endl;
    property_map<Graph, vertex_index_t>::type     index_map = get(vertex_index, G);

    //cout << "original bandwidth: " << bandwidth(G) << endl;

    vector<Vertex> inv_perm(num_vertices(G));
    //vector<size_type> perm(num_vertices(G));

    //reverse cuthill_mckee_ordering
    cout << " hhhhhhhhhhhhhhhhh " << endl;
    cuthill_mckee_ordering(G, inv_perm.rbegin(), get(vertex_color, G), make_degree_map(G));
    cout << " hhhhhhhhhhhhhhhhh " << endl;
    perm.resize(nrow);
    for(size_type c = 0; c != inv_perm.size(); ++c)
      perm[index_map[inv_perm[c]]] = c+1;

    //cout << "  bandwidth: " << bandwidth(G, make_iterator_property_map(&perm[0], index_map, perm[0])) << endl;

  //IPARM[4] = 0;  // user input permutation
  //IPARM[4] = 2;  // return the permutation

  perm.resize(nRow);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  //cout << " ppppppppppp " << IPARM[14] << '\t' << IPARM[15] <<  endl;
  //cout << " llllllllll " << endl;


   pFile2 = fopen("matrix-pattern-with-rcm.dat","w");

   fprintf(pFile2, "%9d \t %9d \n", nRow, nCol );

   //for(int k=0; k<10; ++k)
     //cout << k << '\t' << perm[k] << endl;

   for(int k=0; k<mtx.outerSize(); ++k)
   for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
     fprintf(pFile2, "%9d \t %9d \n", perm[it.row()], perm[it.col()]);

   fclose (pFile2);
*/





