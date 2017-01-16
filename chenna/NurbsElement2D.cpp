
#include "Debug.h"
#include "NurbsElement2D.h"

using namespace std;



NurbsElement2D::NurbsElement2D(void)
{
  if (debug) cout << " constructor NurbsElement2D\n\n";

  // cout << "     NurbsElement2D: constructor ...\n\n";

  nlbf = ndof = nsize = nivGP = nGP = elenum = patchnum = 0;

  surf0  = NULL;
  surf1  = NULL;
}





NurbsElement2D::~NurbsElement2D()
{
  if(debug)   cout << " destructor NurbsElement2D\n\n ";

  if(surf0 != NULL)  // dont delete this pointer as it is just a reference
     surf0 = NULL;

  if(surf1 != NULL)  // dont delete this pointer as it is just a reference
     surf1 = NULL;

}




void NurbsElement2D::initialiseDOFvalues()
{
    int index1, index2;
    for(int ii=0;ii<nlbf;ii++)
    {
       index1 = ndof * ii;
       index2 = ndof * (surf0->IEN[elenum][ii]);
       for(int jj=0;jj<ndof;jj++)
          primvar[index1+jj] = surf0->Uinit[index2+jj];
    }


  return;
}







/*

void NurbsElement2D::prepareElemData()
{
    if(curve1 != NULL)
      nlbf = curve1->p + 1;  // no. of local basis functions
    if(surf1 != NULL)
      nlbf = surf0->nlbf;  // no. of local basis functions
      
    nsize = nlbf*ndof;

    forassembly.setDim(nsize);
    for(int ii=0;ii<nsize;ii++)
      forassembly[ii].setDim(nsize);

    tracdata.setDim(4);
    for(int ii=0;ii<4;ii++)
    {
      tracdata[ii].setDim(2);
      tracdata[ii].zero();
    }

    primvar.setDim(nsize);

    stiffness_local.setDim(nsize);
    for(int ii=0;ii<nsize;ii++)
      stiffness_local[ii].setDim(nsize);

    resi.setDim(nsize);

    // get the Gauss Points and Weights
    //GaussPoints2D(surf0->p, surf0->q, gausspoints, gaussweights);

  return;
}
*/



void NurbsElement2D::diffStiffTest(double ddd, int dig, int dig2, bool gfrmt)
{
  int    i, jnd, jj, j, k;

  double dd[6]  = {-3.*ddd, -2.*ddd, -ddd, +ddd, +2.*ddd, +3.*ddd },
         *r     = new double[6*nsize];

  ListArray<VectorArray<double> > sdiff, sdiff2;

  sdiff.setDim(nsize);
  sdiff2.setDim(nsize);
    for(int ii=0;ii<nsize;ii++)
    {
      stiffness_local[ii].zero();
      sdiff[ii].setDim(nsize);
      sdiff2[ii].setDim(nsize);
      sdiff[ii].zero();
      sdiff2[ii].zero();
    }

  int loc_num=0;
  vector<int>  uindex(nlbf);

  for(int jj=0; jj<=surf0->q; jj++)
  {
    for(int ii=0; ii<=surf0->p; ii++)
    {
      uindex[loc_num] = (startindex[0]+ii) + (surf0->ngbf1 * (startindex[1]+jj));
      loc_num++;
    }
  }
    for(int ii=0; ii<nlbf; ii++)
       cout << uindex[ii] << '\t';
  cout << endl;


  // loop over columns of s

  for (jnd=0; jnd<nlbf; jnd++) // nodes
  {
    for (jj=0; jj<ndof; jj++)  // dof
    {
      j = jnd*ndof + jj;

      // loop over perturbations
	    
      for (k=0; k<6; k++)
      {
        // apply pertubation
        surf1->updateCoordsSingleCP(uindex[jnd], jj, dd[k]);

        // calculate residual
        calcInternalForces();

        // remove pertubation

        surf1->updateCoordsSingleCP(uindex[jnd], jj, -dd[k]);

        // loop over rows of s, store residual

        for (i=0; i<nsize; i++) r[i*6+k] = resi[i];
      }

      // loop over rows of s
      for (i=0; i<nsize; i++)
        sdiff[j][i] = ( +       r[i*6+0]
                           -  9. * r[i*6+1]
                           + 45. * r[i*6+2]
                           - 45. * r[i*6+3]
                           +  9. * r[i*6+4]
                           -       r[i*6+5] ) / (60. * ddd);
      
    }
  }


  // calculate stiffness
  calcStiffnessMatrix(1.0);


  cout.setf(ios::fixed);
  cout.setf(ios::showpoint);
  cout.precision(5);

  cout << "    Analytical Stiffness Matrix   " << endl;
  cout << endl;
  cout << endl;
  for(int ii=0;ii<nsize;ii++)
  {
    cout << '\t' ;
    for(int jj=0;jj<nsize;jj++)
    {
      cout  <<  fixed << stiffness_local[ii][jj] << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;


  cout << "    diff Stiffness Matrix   " << endl;
  cout << endl;
  cout << endl;
  for(int ii=0;ii<nsize;ii++)
  {
    cout << '\t' ;
    for(int jj=0;jj<nsize;jj++)
    {
      sdiff2[ii][jj] = stiffness_local[ii][jj] - sdiff[ii][jj];
      cout  <<  fixed <<  sdiff[ii][jj] << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;

  cout << "    difference of Stiffness Matrices   " << endl;
  cout << endl;
  cout << endl;
   for(int ii=0;ii<nsize;ii++)
  {
    cout << '\t' ;
    for(int jj=0;jj<nsize;jj++)
    {
      cout  <<  fixed << sdiff2[ii][jj] << "   " ;
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;


/*
  prgCompareTwoSimpleMatrices(sdiff,                         // matrix 1
		              s,                             // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      nsize, nsize,                      // matrix dimension
			      dig,dig2,gfrmt,                // format
			      0,                             // indentation
			      false,                         // interactive
			      false);                        // row/column numbers
*/
  delete [] r;
  //delete [] sdiff;

  return;
}


