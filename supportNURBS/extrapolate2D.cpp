#include <iostream>
#include <math.h>
#include "NurbsShapeFunctions.h"

using namespace std;



double extrapolate2D(int ngp1, int ngp2, int index, double* val)
{
	int  ii, jj, ind, n;
	double  xi=0.0, eta=0.0, u=0.0, v=0.0, output, value1, value2;

	vector<double>  Nu(ngp1), Nv(ngp2);

	switch(ngp1)
	{
		case 2:	value1 = sqrt(3.0);		break;

		case 3:	value1 = sqrt(5.0/3.0);	break;

		case 4:	value1 = 1.0/sqrt( (3+2.0*sqrt(6.0/5.0))/7.0);	break;

		case 5:	value1 = 3.0/sqrt(5+2.0*sqrt(10.0/7.0));		break;

		default:
			cerr << '\t' << " Not implemented yet " << endl;
			break;
	}
	switch(ngp2)
	{
		case 2:	value2 = sqrt(3.0);		break;

		case 3:	value2 = sqrt(5.0/3.0);	break;

		case 4:	value2 = 1.0/sqrt( (3+2.0*sqrt(6.0/5.0))/7.0);	break;

		case 5:	value2 = 3.0/sqrt(5+2.0*sqrt(10.0/7.0));		break;

		default:
			cerr << '\t' << " Not implemented yet " << endl;
			break;
	}

	switch(index)
	{
		case 1:	xi  = -value1;		eta = -value2;		break;

		case 2:	xi  =  value1;		eta = -value2;		break;

		case 3:	xi  = -value1;		eta =  value2;		break;

		case 4:	xi  =  value1;		eta =  value2;		break;
	}


	u = 0.5*(1+xi);
	v = 0.5*(1+eta);

	//cout << '\t' << " xi and eta " << xi << '\t' << eta << endl;
	//cout << '\t' << " u and v " << u << '\t' << v << endl;

	n = ngp1-1;
	for(ii=0;ii<ngp1;ii++)
	{
		ind = n-ii;
		Nu[ii] = Bin(n,ii)*pow(u,ii)*pow((1-u),ind);
	}
	n = ngp2-1;
	for(ii=0;ii<ngp2;ii++)
	{
		ind = n-ii;
		Nv[ii] = Bin(n,ii)*pow(v,ii)*pow((1-v),ind);
		//cout << '\t' << L[ii] << '\t' << R[ii] << endl;
	}


	output = 0.0;
        ind = 0;
	for(jj=0;jj<ngp2;jj++)
	{
           for(ii=0;ii<ngp1;ii++)
           {
               output += val[ind++]*Nv[jj]*Nu[ii];
           }
	}

	return output;
}
