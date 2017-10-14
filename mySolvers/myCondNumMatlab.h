
#ifndef MY_COND_NUMBER_MATLAB_H
#define MY_COND_NUMBER_MATLAB_H


#include "headersBasic.h"


/* $Revision: 0.0.1 $ */
/*
 *	myCondNumMatlab.h
 *
 *	This subroutine computes condition number of a sparse matrix
 *	using MATLAB subroutine 'condest'
 *	MATLAB functions are called from C++
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "engine.h"


#include "headersBasic.h"

typedef  SparseMatrix<double, RowMajor>  SparseMatrixXd;


int myCondNumMatlab(SparseMatrixXd& mtx)
{
    double  nRow[1]={mtx.rows()}, nnz[1]={mtx.nonZeros()};

    vector<double>  rows(nnz[0]);
    vector<double>  cols(nnz[0]);
    vector<double>  data(nnz[0]);

    int ii=0, k=0;

    for(k=0; k<mtx.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
      {
        //fprintf(pFile, "%9d \t %9d \t %20.16f \n", it.row(), it.col(), it.value());
        rows[ii] = it.row() + 1.0;
        cols[ii] = it.col() + 1.0;
        data[ii] = it.value();
        ii++;
      }
    }


    //double rows[nnz] = {1, 1, 2, 2, 2, 3, 3, 3, 3, 4, 5, 5 };
    //double cols[nnz] = {1, 2, 1, 3, 5, 2, 3, 4, 5, 3, 2, 5 };
    //double data[nnz] = {2.0,  3.0, 3.0, -1.0, 4.0, 4.0, -3.0, 1.0, 2.0, 2.0, 6.0, 1.0};


    /*
     * Start the MATLAB engine locally by executing the string
     * "matlab"
     *
     * For more complicated cases, use any string with whitespace,
     * and that string will be executed literally to start MATLAB
    */

    Engine *ep;

    //if(!(ep = engOpen("\0")))
    if(!(ep = engOpen("matlab -nodesktop")))
    {
        fprintf(stderr, "\nCan't start MATLAB engine\n");
        return EXIT_FAILURE;
    }

    // Create MATLAB variables for our data
    // Place the variable T into the MATLAB workspace

    mxArray *RR = mxCreateDoubleMatrix(1, nnz[0], mxREAL);
    memcpy((void *)mxGetPr(RR), (void *)&rows[0], sizeof(double)*nnz[0]);
    engPutVariable(ep, "RR", RR);

    mxArray *CC = mxCreateDoubleMatrix(1, nnz[0], mxREAL);
    memcpy((void *)mxGetPr(CC), (void *)&cols[0], sizeof(double)*nnz[0]);
    engPutVariable(ep, "CC", CC);

    mxArray *VV = mxCreateDoubleMatrix(1, nnz[0], mxREAL);
    memcpy((void *)mxGetPr(VV), (void *)&data[0], sizeof(double)*nnz[0]);
    engPutVariable(ep, "VV", VV);


    //engEvalString(ep, "nrow = 5;" );
    //engEvalString(ep, "nz   = 12;" );

    mxArray *nrow = mxCreateDoubleMatrix(1, 1, mxREAL);
    memcpy((void *)mxGetPr(nrow), (void *)nRow, sizeof(nRow));
    engPutVariable(ep, "nrow", nrow);

    mxArray *nz = mxCreateDoubleMatrix(1, 1, mxREAL);
    memcpy((void *)mxGetPr(nz), (void *)nnz, sizeof(nnz));
    engPutVariable(ep, "nz", nz);


    //engEvalString(ep, "plot(RR, CC);");

    engEvalString(ep, "K = sparse(RR, CC, VV, nrow, nrow, nz);");

    // Plot the pattern

    //engEvalString(ep, "spy(K);");

    engEvalString(ep, "condNum = condest(K,1);");

    // use fgetc() to make sure that we pause long enough to be able to see the plot

    //printf("Hit return to continue\n\n");
    //fgetc(stdin);

    // Get result of computation
    mxArray *result = engGetVariable(ep,"condNum") ;

    double  cnum = mxGetScalar(result);

    printf("\n\n Matrix condition number %12.6E \n\n\n", cnum);


    // We're done! Free memory, close MATLAB engine and exit.

    //printf("Done...\n");

    engEvalString(ep, "close;");

    mxDestroyArray(RR);
    mxDestroyArray(CC);
    mxDestroyArray(VV);
    mxDestroyArray(result);

    engClose(ep);

    return 1;
}





#endif

