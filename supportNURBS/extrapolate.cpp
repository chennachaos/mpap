#include <iostream>
#include <math.h>
#include "NurbsShapeFunctions.h"

using namespace std;



double extrapolate(int ngp, int index, double* val)
{
	int ii, jj, ind;
	double xi=0.0, eta=0.0, u=0.0, v=0.0, output, value;

	vector<double>  L(ngp), R(ngp), temp(ngp);

	switch(ngp)
	{
		case 2:	value = 1.73205081;	break;

		case 3:	value = 1.29099445;	break;

		case 4:	value = 1.16125634;	break;

		case 5:	value = 1.1035337;	break;

		case 6:	value = 1.07242112;	break;

		case 7:	value = 1.05362097;	break;

		default:
			cerr << '\t' << " Not implemented yet " << endl;
			break;
	}

	switch(index)
	{
		case 1:	xi  = -value;		eta = -value;		break;

		case 2:	xi  =  value;		eta = -value;		break;

		case 3:	xi  = -value;		eta =  value;		break;

		case 4:	xi  =  value;		eta =  value;		break;
	}


	u = 0.5*(1+xi);
	v = 0.5*(1+eta);

	//cout << '\t' << " xi and eta " << xi << '\t' << eta << endl;
	//cout << '\t' << " u and v " << u << '\t' << v << endl;

	int n = ngp-1;
	for(ii=0;ii<ngp;ii++)
	{
		ind = n-ii;
		L[ii] = Bin(n,ii)*pow(u,ii)*pow((1-u),ind);
		R[ii] = Bin(n,ii)*pow(v,ii)*pow((1-v),ind);

		//cout << '\t' << L[ii] << '\t' << R[ii] << endl;
	}

	output = 0.0;
       ind = 0;
	for(jj=0;jj<ngp;jj++)
	{
           temp[jj] = 0.0;
           for(ii=0;ii<ngp;ii++)
           {
               temp[jj] += val[ind]*L[ii];
               ind++;
           }
           output += temp[jj]*R[jj];
	}

	return output;
}
