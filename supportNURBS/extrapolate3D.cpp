#include <iostream>
#include <math.h>
#include "NurbsShapeFunctions.h"

using namespace std;


inline double getGPforExtrapolation(int ngp1)
{
    double value;

    switch(ngp1)
    {
	case 2:	value = sqrt(3.0);		break;

	case 3:	value = sqrt(5.0/3.0);	break;

	case 4:	value = 1.0/sqrt( (3+2.0*sqrt(6.0/5.0))/7.0);	break;

	case 5:	value = 3.0/sqrt(5+2.0*sqrt(10.0/7.0));		break;

	default:
		cerr << '\t' << " Not implemented yet " << endl;
		break;
	}
	
	return value;
}


double extrapolate3D(int ngp1, int ngp2, int ngp3, int index, double* val)
{
	int  ii, jj, kk, ind, n;
	double  u, v, w, output, value1, value2, value3, temp;

	vector<double>  Nu(ngp1), Nv(ngp2), Nw(ngp3);

	value1 = getGPforExtrapolation(ngp1);
	value2 = getGPforExtrapolation(ngp2);
	value3 = getGPforExtrapolation(ngp3);

	switch(index)
	{
		case 1:	u = 0.5*(1-value1);       v = 0.5*(1-value2);       w = 0.5*(1-value3);		break;

		case 2:	u = 0.5*(1+value1);       v = 0.5*(1-value2);       w = 0.5*(1-value3);		break;

		case 3:	u = 0.5*(1-value1);       v = 0.5*(1+value2);       w = 0.5*(1-value3);		break;

		case 4:	u = 0.5*(1+value1);       v = 0.5*(1+value2);       w = 0.5*(1-value3);		break;

		case 5:	u = 0.5*(1-value1);       v = 0.5*(1-value2);       w = 0.5*(1+value3);		break;

		case 6:	u = 0.5*(1+value1);       v = 0.5*(1-value2);       w = 0.5*(1+value3);		break;

		case 7:	u = 0.5*(1-value1);       v = 0.5*(1+value2);       w = 0.5*(1+value3);		break;

		case 8:	u = 0.5*(1+value1);       v = 0.5*(1+value2);       w = 0.5*(1+value3);		break;

	}

	//cout << '\t' << " u, v and w " << u << '\t' << u << '\t' << w << endl;

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
	}
	n = ngp3-1;
	for(ii=0;ii<ngp3;ii++)
	{
            ind = n-ii;
	    Nw[ii] = Bin(n,ii)*pow(w,ii)*pow((1-w),ind);
	}

        ind = 0;
        output = 0.0;
        for(kk=0;kk<ngp3;kk++)
        {
            for(jj=0;jj<ngp2;jj++)
            {
               temp = Nw[kk] * Nv[jj];
               for(ii=0;ii<ngp1;ii++)
               {
                   output += val[ind++]*temp*Nu[ii];
               }
            }
	}

	return output;
}
