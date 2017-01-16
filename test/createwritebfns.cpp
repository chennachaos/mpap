


#include "headersBasic.h"

#include "HBSplineFEM.h"
#include "MathBasic.h"
#include "DataBlockTemplate.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "PlotVTK.h"
#include "PropertyTypeEnum.h"

#include "Global.h"
#include "MyString.h"
#include "FunctionsEssGrp.h"

#include "UnixGlobal.h"

#include <omp.h>

#include <iostream>
#include <fstream>
#include <string.h>
#include <iomanip>
//#include "TMeshConnectivity.h"

#include "NURBS_1D.h"

using namespace std;

//void createAndPrintBasisFunctions1D()
int main(int argc, char **argv)
{
	VectorArray<double> U, uu, Utmp;
	
	int m, p, n, ii, jj;
	m=10;p=2;
	
	double incr = 1.0;

	/*
	if(argc > 0)
	{
	   m = atoi(argv[1]);
	   p = atoi(argv[2]);
	   incr = atof(argv[3]);
	}
	//
	
	n=m-p-1;

	U.setDim(m+1);

	for(int ii=0;ii<U.n;ii++)
           U[ii] = ((double) ii) * incr;
           
        if(m-2*p > 0)
        {
           Utmp.setDim(m-2*p+1);
           for(int ii=0;ii<Utmp.n;ii++)
              Utmp[ii] = U[p+ii];
        }

	cout << U << endl;
	cout << Utmp << endl;
        */

        vector<double> vectmp;

        p = 3;

        for(ii=0; ii<=p; ii++)
          vectmp.push_back(0.0);

        vectmp.push_back(2.0);

        for(ii=0; ii<=p; ii++)
          vectmp.push_back(4.0);

	U.setDim(vectmp.size());

        for(ii=0; ii<vectmp.size(); ii++)
          U[ii] = vectmp[ii];


	//create_vector2(Utmp, 20, uu);
	
	//cout << Utmp << endl;

	//create_vector2(U, 20, uu);
        create_vector(0.0, U[U.n-1], 0.01, uu);

        m = vectmp.size() - 1;

	n = m - p - 1;

	ListArray<VectorArray<double> > NN, NN_der;
	NN.setDim(n+1);
	NN_der.setDim(n+1);

	for(ii=0;ii<=n;ii++)
	{
            NN[ii].setDim(uu.n);
            NN_der[ii].setDim(uu.n);
            NN[ii].zero();
            NN_der[ii].zero();
	}


	int deriv = 1;

	unsigned int ROWS = deriv+1;
	double** ders1 = new double*[ROWS];

        for(ii=0;ii<ROWS;ii++)
            ders1[ii] = new double[p+1];


	int span, index;
	
	//
	for(int jj=0;jj<uu.n;jj++)
	{
		span = FindSpan(&U[0], U.n, p, uu[jj]);

		DersBasisFuns(&U[0], U.n, p, uu[jj], deriv, ders1);

		index = span - p;
		for(int ii=0;ii<=p;ii++)
		{
			NN[index+ii][jj] = ders1[0][ii];
			NN_der[index+ii][jj] = ders1[1][ii];
		}
	}
	//
	
	/*
	for(ii=0;ii<=n;ii++)
	{
		//span = FindSpan(&U[0], U.n, p, uu[jj]);

		//DersBasisFuns(&U[0], U.n, p, uu[jj], deriv, ders1);

		//index = span - p;
		for(jj=0;jj<uu.n;jj++)
		{
			NN[ii][jj] = oneBasisAlg1(U, p, ii, uu[jj]);
			
			NN_der[ii][jj] = DersOneBasisFun_recurs(U, p, ii, uu[jj], 2);
			
			//cout << ii << '\t' << jj << '\t' << uu[jj] << endl;
			
			//cout << " span =" <<  FindSpan(&U[0], U.n, p, uu[jj]) << endl;
			
			//NN[ii][jj] = OneBasisFun_recurs(U, p, ii, uu[jj]);
		}
	}
        */


	ofstream fout("NURBS_Basis_Functions.dat");

	if(fout.fail())
	{
		cout << " Could not open the Output file" << endl;
		exit(1);
	}

	fout.setf(ios::fixed);
	fout.setf(ios::showpoint);
	fout.precision(6);

	for(ii=0;ii<uu.n;ii++)
	{
		fout << setw(10) << uu[ii];
		for(jj=0;jj<=n;jj++)
			fout << setw(14) << NN[jj][ii] ;
		fout << endl;
	}

	fout.close();

	ofstream fout2("NURBS_Basis_Functions_Derivatives.dat");

	if(fout2.fail())
	{
		cout << " Could not open the Output file" << endl;
		exit(1);
	}

	fout2.setf(ios::fixed);
	fout2.setf(ios::showpoint);
	fout2.precision(6);

	for(ii=0;ii<uu.n;ii++)
	{
		fout2 << setw(10) << uu[ii];
		for(jj=0;jj<=n;jj++)
			fout2 << setw(14) << NN_der[jj][ii] ;
		fout2 << endl;
	}

	fout2.close();

	for(ii=0;ii<ROWS;ii++)
		delete [] ders1[ii];

        delete [] ders1;

	return 0;
}
