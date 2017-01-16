
#include <iostream>
#include <assert.h>

#include "MathVector.h"
#include "Solver.h"
#include "SolverMA41.h"
//#include "SolverPARDISO.h"
#include "DataBlockTemplate.h"
#include "List.h"
#include "MyString.h"
#include "MyStringList.h"
#include "SolverWorkSpace.h"
#include "SolverTime.h"
#include "ComputerTime.h"


extern SolverWorkSpace solverWorkSpace;
extern SolverTime      solverTime;
extern ComputerTime    computerTime;


using namespace std;

/*
int main(int argc,char *argv[])
{
	// create the solver and other required objects
	/////////////////////////////////////

    char fct[] = "pardiso test";

	Solver *solver;

    MyStringList key;

    MatrixSparse<double> matxtemp;

    VectorArray<int>  CSRcol, CSRrow, compr1;

    VectorArray<double>  CSRval, rhs;

    MyString line, tmpl, *word;

    ifstream Ifile;

    char tmp[30];

    int nnodes, nnz, ii, jj;

    List<Vector<int> > lviTmp;
    List<Vector<double> > lvdTmp;


    // read the input file
    ///////////////////////////////////////

    key.addNew("nnodes=",
    	       "nnz=",
    	       "CSRcol=",
    	       "CSRval=",
    	       "CSRrow=",
    	       "RHS=",
               "END");

    //key.print();

    cout << "  Reading the input file ... \n\n " << endl;

    //Ifile->open("pardiso1_matinCSR_RHS.dat");

    Ifile.open("pardiso_ex1.dat");

    line.read(Ifile);

    //
    //MyString keyStrg;

    char *ss;
    ss = line.asCharArray();

    cout << " >" << line << "<\n\n";
    cout << " >" << line.stripToMin() << "<\n\n";

    tmpl = "nnodes=";

    printf("%s\n", ss);
    printf("%d\n", line.length());
    printf("%d\n", tmpl.length());

    printf("\t %d\n",key.whichBegins(line));

    printf("\t %d\n",(line=="nnodes"));

    printf("\t %d\n",key.whichBegins("nnodes"));
    printf("\t %d\n",key.whichBegins("nnz"));
    //
	
    do
    {
    	//line.stripToMin();

        cout << " >" << line.stripToMin() << "<\n\n";

    	lviTmp.free();
    	lvdTmp.free();

    switch(key.whichBegins(line))
    {
        case  0: //nnodes

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid node number or invalid keyword!");

            assert(lviTmp.n == 1);

            nnodes = lviTmp[0][0];

            printf("\t nnodes = %d \n\n", nnodes);

            break;

        case  1: //nnz

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid node number or invalid keyword!");

            assert(lviTmp.n == 1);

            nnz = lviTmp[0][0];

            printf("\t nnz = %d \n\n", nnz);

            break;

        case  2: //CSRcol

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid node number or invalid keyword!");

            printf("\t CSRcol \n");
            printf("\t lviTmp.n = %d \n\n", lviTmp.n);

            assert(lviTmp.n == nnz);

            CSRcol.setDim(nnz);

            for(ii=0;ii<nnz;ii++)
            	CSRcol[ii] = lviTmp[ii][0];

            break;

        case  3: //CSRval

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
               prgError(10,fct,"invalid knot vector or invalid keyword!");

            printf("\t CSRval \n");
            printf("\t lvdTmp.n = %d \n\n", lvdTmp.n);

            assert(lvdTmp.n == nnz);

            CSRval.setDim(nnz);

            for(ii=0;ii<nnz;ii++)
            	CSRval[ii] = lvdTmp[ii][0];

            break;

        case  4: //CSRrow

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid node number or invalid keyword!");

            printf("\t CSRrow \n");
            printf("\t lviTmp.n = %d \n\n", lviTmp.n);

            assert(lviTmp.n == (nnodes+1));

            CSRrow.setDim(nnodes+1);

            //for(ii=0;ii<=nnodes;ii++)
              //printf("\t %d \n", lviTmp[ii][0]);

            for(ii=0;ii<=nnodes;ii++)
            	CSRrow[ii] = lviTmp[ii][0];

            break;

        case  5: //RHS

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
               prgError(1,fct,"invalid knot vector or invalid keyword!");

            printf("\t RHS \n");
            printf("\t lvdTmp.n = %d \n\n", lvdTmp.n);

            assert(lvdTmp.n == nnodes);

            rhs.setDim(nnodes);

            for(ii=0;ii<nnodes;ii++)
            	rhs[ii] = lvdTmp[ii][0];

            break;
        case  6: //END

            break;
    }
    }
    while(prgNextDataKey(Ifile, line));

    Ifile.close();

    //
    // Assign data to matrix object in the solver
    /////////////////////////////////////////////

    int numProc, MAX_PROCESSORS;

    bool solverOK, cIO = true, IS_SYM = false;

    if(argc > 1)
        numProc = atoi(*argv[1]);
    else
    	numProc = 1;

    numProc = min(MAX_PROCESSORS,numProc);

    cout << " numProc " <<  numProc << endl;

    solver = (Solver*) new SolverPARDISO;

    solver->printInfo();

    if(IS_SYM)
    	if(solver->initialise(numProc,PARDISO_STRUCT_SYM) != 0)
    		return;
    else
    	if(solver->initialise(numProc,PARDISO_UNSYM) != 0)
    		return;

    solverOK = true;

    if(solver != NULL)
        solver->checkIO = cIO;



    // Assign data to matrix object in the solver
    ///////////////////////////////////////

    ((SolverSparse*)solver)->compr = CSRrow;

    ((SolverSparse*)solver)->mtx.row.x = CSRrow;

    ((SolverSparse*)solver)->mtx.col.x = CSRcol;

    ((SolverSparse*)solver)->mtx.x.x = CSRval;

    ((SolverSparse*)solver)->currentStatus = PATTERN_OK;

    double  *soln;

    soln = solver->factoriseAndSolve(rhs.x);

    printf("\t Solution ... \n\n");
    for(ii=0;ii<nnodes;ii++)
    	printf("\t %12.6f \n", soln[ii]);
    printf("\n");
    //


	return 0;
}
*/