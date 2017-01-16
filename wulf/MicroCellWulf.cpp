
#include <iostream>

#include "FunctionsProgram.h"
#include "DomainTypeEnum.h"
#include "MicroCellWulf.h"
#include "DomainTree.h"
#include "DomainType.h"
#include "Debug.h"
#include "DependentDoF.h"
#include "TimeFunction.h"
#include "MathGeom.h"


extern DomainTree domain;
extern List<TimeFunction> timeFunction;


MicroCellWulf::MicroCellWulf(void)                       
{                                                  
  if (debug) std::cout << " MicroCellWulf constructor\n\n";

  // add new type
  
  DomainType *microCellWulf = domain.newType(MICROCELLWULF,SOLID);
  
  if (microCellWulf == NULL) return;  // domain type already exists

  microCellWulf->key.addNew("boundary conditions", "other stuff");
}




	                                          
MicroCellWulf::~MicroCellWulf(void)                     
{                                                 
  if (debug) std::cout << " MicroCellWulf destructor\n\n";
}






void MicroCellWulf::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString *word;
 
  char fct[] = "MicroCellWulf::readInputData";

  int nw, i;
 
  switch (domain[MICROCELLWULF].key.whichBegins(line))
  {
    case  0: cout << "     MICROCELLWULF: reading boundary conditions ...\n\n";
	     
             this->Solid::readInputData(Ifile, line); break; // this is temporary!
	     
             prgError(1,fct,"no boundary conditions for the MicroCellWulf!");

	     break;
	     
    case  1: cout << "     MICROCELLWULF: reading other stuff ...\n\n"; 
	     
    case -1: // go and inherit from SOLID
	     
	     this->Solid::readInputData(Ifile, line); 
	     
	     break;
  }
 
  return;
}






void MicroCellWulf::prepareInputData(void)
{
  // go and inherit from ancestors

  Solid::prepareInputData();
 
  
  std::cout << "     MICROCELLWULF: prepare input data ...\n\n";

 
  preparePeriodicBoundaryConditions(300);
	 
  return;
}















void MicroCellWulf::preparePeriodicBoundaryConditions(int maxNd)
{
  // This subroutine is closely associated with Mesh::prepareInputData.
  // It adds equally distributed boundary nodes to a 2D square
  // MicroCellWulf such that the 'periodic boundary conditions' can be used.
  // 
  // Depending on the details, most changes of Mesh::prepareInputData 
  // require to be accounted for in this subroutine.
  // I am happy to accept this "ugly" dependency on Mesh::prepareInputData.
  // The following is extremely specific, whereas Mesh::prepareInputData
  // is extremely generic and should not be touched to make life easier here.
  //
  // Note that all the data calculated here, could be input as data in the
  // input file (tedious!).
      
  char fct[] = "MicroCellWulf::preparePeriodicBoundaryConditions";
	
  if (ndm == 3) prgError(1,fct,"the following is not entirely appropriate for 3D!");
	
  int c, i, j, k, m, nNd, i1, i2, n1, n2, n3, n4;
	
  double dd;
  
  // find maximum edge length, max/min coordinates

  n1 = bndNd[nBndNd-1], n2 = bndNd[0]; 
	
  double xmin[3], xmax[3], tol,
	 d, dmax = sqrt((x(n2,1)-x(n1,1))*(x(n2,1)-x(n1,1))
	               +(x(n2,2)-x(n1,2))*(x(n2,2)-x(n1,2)));

  for (j=0; j<ndm; j++)
  {
    xmin[j] = x(1,1+j);
    xmax[j] = x(1,1+j);
  }
  
  for (i=1; i<nBndNd; i++)
  {
    n1 = n2;
    n2 = bndNd[i];
    d  = sqrt((x(n2,1)-x(n1,1))*(x(n2,1)-x(n1,1))
             +(x(n2,2)-x(n1,2))*(x(n2,2)-x(n1,2)));
    if (d > dmax) dmax = d;
    for (j=0; j<ndm; j++)
    {
      if (xmin[j] > x(n2,1+j)) xmin[j] = x(n2,1+j);
      if (xmax[j] < x(n2,1+j)) xmax[j] = x(n2,1+j);
    }
  }
  printf("         xmin         xmax\n");
  for (j=0; j<ndm; j++) printf("%12.5g %12.5g\n",xmin[j],xmax[j]);
  cout << "\n       max d = " << dmax << "\n\n";

  tol = dmax * 1.e-4;
  
  d = 0.; for (j=0; j<ndm; j++)  if (d < xmax[j] - xmin[j]) d = xmax[j] - xmin[j];

 
  // get corner nodes (the followin is strictly 2D)

  n1 = 0; n2 = 0; n3 = 0; n4 = 0;
  for (i=1; i<nBndNd; i++)
  {
    if (x(bndNd[i],1) + x(bndNd[i],2) < x(bndNd[n1],1) + x(bndNd[n1],2)) n1 = i;
    if (x(bndNd[i],1) - x(bndNd[i],2) > x(bndNd[n2],1) - x(bndNd[n2],2)) n2 = i;
    if (x(bndNd[i],1) + x(bndNd[i],2) > x(bndNd[n3],1) + x(bndNd[n3],2)) n3 = i;
    if (x(bndNd[i],1) - x(bndNd[i],2) < x(bndNd[n4],1) - x(bndNd[n4],2)) n4 = i;
  }
  n1 = bndNd[n1];
  n2 = bndNd[n2];
  n3 = bndNd[n3];
  n4 = bndNd[n4];

  if (abs(x(n1,1) - x(n4,1)) > tol 
   || abs(x(n1,2) - x(n2,2)) > tol
   || abs(x(n3,1) - x(n2,1)) > tol
   || abs(x(n3,2) - x(n4,2)) > tol) prgError(1,fct,"invalid MicroCellWulf geometry!");
  
  cout << "       corner nodes: " << n1 << ", " << n2 << ", " << n3 << ", " << n4 << "\n\n";
  
  
  // determine number of new nodes

  nNd = (int) floor(d/dmax+1.e-10);

  if (nNd > maxNd) nNd = maxNd;

  cout << "       divisions per edge = " << nNd << "\n\n";
  
  
  // generate new nodes	(the following is strictly 2D)

  int  oldNumnpXndm = numnp * ndm;
  
  oldNumnp = numnp; 

  numnp += 2 * (nNd-1);

  double *xx = new double[numnp*ndm];
  
  for (i=0; i<oldNumnpXndm; i++) xx[i] = x.x[i];
  delete [] x.x; 
  x.x = xx; 
  x.nRow = numnp; 

  d /= (double) nNd;
  nNd--;
  dd = d;
  oldNumnp++;
  for (i=0; i<nNd; i++)
  {
    x(oldNumnp+i,1) = xmin[0];
    x(oldNumnp+i,2) = xmax[1] - dd;

    x(oldNumnp+nNd+i,1) = xmin[0] + dd;
    x(oldNumnp+nNd+i,2) = xmin[1];

    dd += d;
  }
  oldNumnp--;
  
  // adjust lots of other data fields to increased node number

  VectorArray<int> *ndEl = new VectorArray<int> [numnp];

  double *xx0 = new double[numnp*ndm],
	 *xxn = new double[numnp*ndm],
	 *uu  = new double[numnp*ndf],
	 *uun = new double[numnp*ndf],
	 *uu3 = new double[numnp*ndf],
	 *uu4 = new double[numnp*ndf],
	 *uu5 = new double[numnp*ndf],
	 *uu6 = new double[numnp*ndf];

  int	*iidx = new int[numnp*ndm], 
	*iidu = new int[numnp*ndf], iindm = 0, iindf = 0;

  for (i=0; i<oldNumnp; i++)
  {
    for (j=nodeElem[i].n-1; j>-1; j--)
      ndEl[i][j] = nodeElem[i][j];

    nodeElem[i].free();
    
    for (j=0; j<ndm; j++)
    {
      xx0[iindm+j] = x.x[iindm+j];
      xxn[iindm+j] = x.x[iindm+j];
    }
    
    for (j=0; j<ndf; j++)
    {
      uu [iindf+j] =  u.x[iindf+j];
      uun[iindf+j] = un.x[iindf+j];
      uu3[iindf+j] = u3.x[iindf+j];
      uu4[iindf+j] = u4.x[iindf+j];
      uu5[iindf+j] = u5.x[iindf+j];
      uu6[iindf+j] = u6.x[iindf+j];
    }
    iindm += ndm;
    iindf += ndf;
  }
  
  for (i=oldNumnp; i<numnp; i++)
  {
    for (j=0; j<ndm; j++) 
    {
      xx0[iindm+j] = x.x[iindm+j];
      xxn[iindm+j] = x.x[iindm+j];
    }

    for (j=0; j<ndf; j++)
    {
       uu[iindf+j] = 0.;
      uun[iindf+j] = 0.;
      uu3[iindf+j] = 0.;
      uu4[iindf+j] = 0.;
      uu5[iindf+j] = 0.;
      uu6[iindf+j] = 0.;
    }
    iindm += ndm;
    iindf += ndf;
  }

  delete [] nodeElem; nodeElem = ndEl;
  
  delete [] idx.x; idx.x = iidx; idx.nRow = numnp;
  delete [] idu.x; idu.x = iidu; idu.nRow = numnp;
  
  delete []  x0.x;  x0.x =  xx0;  x0.nRow = numnp;
  delete []  xn.x;  xn.x =  xxn;  xn.nRow = numnp;
  delete []   u.x;   u.x =   uu;   u.nRow = numnp;
  delete []  un.x;  un.x =  uun;  un.nRow = numnp;
  delete []  u3.x;  u3.x =  uu3;  u3.nRow = numnp;
  delete []  u4.x;  u4.x =  uu4;  u4.nRow = numnp;
  delete []  u5.x;  u5.x =  uu5;  u5.nRow = numnp;
  delete []  u6.x;  u6.x =  uu6;  u6.nRow = numnp;

  // calculate uDep data

  uDep.free();
  
    // corner nodes

  m = 0;
  
  for (i=0; i<4; i++)
  {
    for (j=0; j<ndf; j++)
    {
      uDep.add(new DependentDoF);
      
      uDep[m].tmFct = 0;
      uDep[m].dof   = j+1;
      uDep[m++].ucBase  = 1.;
    }
  }
  for (j=0; j<ndf; j++)
  {
    uDep[m-8+j].nd = n1;
    uDep[m-6+j].nd = n2;
    uDep[m-4+j].nd = n3;
    uDep[m-2+j].nd = n4;
  }
 
    // original boundary nodes (positive boundary)

  for (i=0; i<nBndNd; i++)  // n4 -> n3
  {
    k = bndNd[i];
    if (k != n1 && k != n2 && k != n3 && k != n4)
    {    
       if (abs(x(k,2) - x(n3,2)) < tol)
      {
        if (x(k,1) - x(n4,1) < tol || x(n3,1) - x(k,1) < tol) 
          prgError(4,fct,"invalid MicroCellWulf geometry!");

        i1 = (int) floor((x(k,1) - x(n4,1)) / d) + oldNumnp + nNd;
	i2 = i1 + 1;
	
        if (i1 == oldNumnp + nNd)           i1 = n1;
	if (i2 == oldNumnp + nNd + nNd + 1) i2 = n2;
	
	for (j=0; j<ndf; j++)
	{
          uDep.add(new DependentDoF); 
	
          uDep[m].nd           = k;
          uDep[m].dof          = j+1;
          uDep[m].tmFct        = 0;
	  uDep[m].ucBase           = 0.;
          uDep[m].masterNd[0]  = i1;
          uDep[m].masterNd[1]  = i2;
          uDep[m].masterDoF[0] = j+1;
          uDep[m].masterDoF[1] = j+1;
          uDep[m].alpha[0]     = (x(i2,1)-x(k,1)) / d;
          uDep[m].alpha[1]     = (x(k,1)-x(i1,1)) / d;
          uDep[m].beta[0]      = uDep[m].alpha[0];
          uDep[m].beta[1]      = uDep[m].alpha[1]; m++;
	}
      }
    }
  }

  ndPosBnd2D[0] = (m - 8) / 2;
  
  for (i=0; i<nBndNd; i++)  // n3 -> n2
  {
    k = bndNd[i];
    if (k != n1 && k != n2 && k != n3 && k != n4)
    {    
      if (abs(x(k,1) - x(n3,1)) < tol)
      {
        if (x(k,2) - x(n2,2) < tol || x(n3,2) - x(k,2) < tol) 
          prgError(5,fct,"invalid MicroCellWulf geometry!");

        i1 = (int) floor((x(n3,2) - x(k,2)) / d) + oldNumnp;
	i2 = i1 + 1;

        if (i1 == oldNumnp)           i1 = n4;
	if (i2 == oldNumnp + nNd + 1) i2 = n1;
	
	for (j=0; j<ndf; j++)
	{
          uDep.add(new DependentDoF); 
	
          uDep[m].nd           = k;
          uDep[m].dof          = j+1;
          uDep[m].tmFct        = 0;
	  uDep[m].ucBase           = 0.;
          uDep[m].masterNd[0]  = i1;
          uDep[m].masterNd[1]  = i2;
          uDep[m].masterDoF[0] = j+1;
          uDep[m].masterDoF[1] = j+1;
          uDep[m].alpha[0]     = (x(k,2)-x(i2,2)) / d;
          uDep[m].alpha[1]     = (x(i1,2)-x(k,2)) / d;
          uDep[m].beta[0]      = uDep[m].alpha[0];
          uDep[m].beta[1]      = uDep[m].alpha[1]; m++;
	}
      }
    }
  }
 
  ndPosBnd2D[1] = (m - 8) / 2 - ndPosBnd2D[0];
  
    // original boundary nodes (negative boundary)

  for (i=0; i<nBndNd; i++)
  {
    k = bndNd[i];
   
    if (k != n1 && k != n2 && k != n3 && k != n4)
    {    
      // n4 -> n1
      if (abs(x(k,1) - x(n1,1)) < tol)
      {
        if (x(k,2) - x(n1,2) < tol || x(n4,2) - x(k,2) < tol) 
          prgError(2,fct,"invalid MicroCellWulf geometry!");

        i1 = (int) floor((x(n4,2) - x(k,2)) / d) + oldNumnp;
	i2 = i1 + 1;

        if (i1 == oldNumnp)           i1 = n4;
	if (i2 == oldNumnp + nNd + 1) i2 = n1;
	
	for (j=0; j<ndf; j++)
	{
          uDep.add(new DependentDoF); 
	
          uDep[m].nd           = k;
          uDep[m].dof          = j+1;
          uDep[m].tmFct        = 0;
	  uDep[m].ucBase           = 0.;
          uDep[m].masterNd[0]  = i1;
          uDep[m].masterNd[1]  = i2;
          uDep[m].masterDoF[0] = j+1;
          uDep[m].masterDoF[1] = j+1;
          uDep[m].alpha[0]     = (x(k,2)-x(i2,2)) / d;
          uDep[m].alpha[1]     = (x(i1,2)-x(k,2)) / d;
          uDep[m].beta[0]      = uDep[m].alpha[0];
          uDep[m].beta[1]      = uDep[m].alpha[1]; m++;
	}
      }

      // n1 -> n2
      if (abs(x(k,2) - x(n1,2)) < tol)
      {
        if (x(k,1) - x(n1,1) < tol || x(n2,1) - x(k,1) < tol) 
          prgError(3,fct,"invalid MicroCellWulf geometry!");

        i1 = (int) floor((x(k,1) - x(n1,1)) / d) + oldNumnp + nNd;
	i2 = i1 + 1;

        if (i1 == oldNumnp + nNd)           i1 = n1;
	if (i2 == oldNumnp + nNd + nNd + 1) i2 = n2;
	
	for (j=0; j<ndf; j++)
	{
          uDep.add(new DependentDoF); 
	
          uDep[m].nd           = k;
          uDep[m].dof          = j+1;
          uDep[m].tmFct        = 0;
	  uDep[m].ucBase           = 0.;
          uDep[m].masterNd[0]  = i1;
          uDep[m].masterNd[1]  = i2;
          uDep[m].masterDoF[0] = j+1;
          uDep[m].masterDoF[1] = j+1;
          uDep[m].alpha[0]     = (x(i2,1)-x(k,1)) / d;
          uDep[m].alpha[1]     = (x(k,1)-x(i1,1)) / d;
          uDep[m].beta[0]      = uDep[m].alpha[0];
          uDep[m].beta[1]      = uDep[m].alpha[1]; m++;
	}
      }
    }
  }

/*  for (j=0; j<uDep.n; j++)
  {
    cout << uDep[j].nd << " " << uDep[j].masterNd << " " << uDep[j].alpha << "\n";
  }*/
 
  // calculate idu
 
  for (i=0; i<numnp*ndf; i++) idu.x[i] = 0;
  
  for (i=0; i<uDep.n; i++)  idu(uDep[i].nd,uDep[i].dof) = - 1 - i;
  
  nequ = 0;
  
  for (i=1; i<numnp+1; i++)
    for (j=1; j<ndf+1; j++)
      if (idu(i,j) == 0) idu(i,j) = ++nequ; else if (idu(i,j) > 0) idu(i,j) = 0;

  // set master index in uDep
  
  for (i=0; i<uDep.n; i++)  
    for (j=0; j<uDep[i].masterNd.n; j++)
      uDep[i].master[j] = idu(uDep[i].masterNd[j],uDep[i].masterDoF[j]);

  // initialise for time zero
  
  for (i=0; i<uDep.n; i++)  uDep[i].init();

  /*cout << idu << "\n\n";

  for (i=0; i<uDep.n; i++)
  {
    cout << i+1 << ". " << uDep[i].nd << ": " << uDep[i].masterNd 
                << ", " << uDep[i].alpha << "\n";
  }*/

  return;
}


  
  
 



void MicroCellWulf::strainToBoundaryDisplacement(double *strain)
{
  //
  //  so far, 2D and small strains only  (ndf = ndm = 2)
  //

  char fct[] = "MicroCellWulf::strainToBoundaryDisplacement";
	
  int i, j = 0, k;
      
  double w, h,
	 eps11 = strain[0],
         eps22 = strain[1],
         eps12 = strain[2];
  
  // set time Function identifiers to 0

  for (i=0; i<uDep.n; i++)  uDep[i].tmFct = 0;

  if (timeFunction.n < 1) prgError(1,fct,"no time functions defined!");
  
  // corners
  
  for (i=0; i<4; i++)
  {
    k = uDep[j].nd;   uDep[j++].ucBase = x0(k,1) * eps11 + x0(k,2) * eps12;
                      uDep[j++].ucBase = x0(k,1) * eps12 + x0(k,2) * eps22;
  }

  // positive boundary

  w = x(uDep[2].nd,1) - x(uDep[0].nd,1);
  h = x(uDep[4].nd,2) - x(uDep[2].nd,2);
  
  //cout << ndPosBnd2D[0] << " " << ndPosBnd2D[1] << "\n";
  //cout << w << " = w, h = " << h << "\n\n";
  
  for (i=0; i<ndPosBnd2D[0]; i++)
  {
    k = uDep[j].nd;   
  
    //cout << k << "\n";
    
    uDep[j++].ucBase = h * eps12;
    uDep[j++].ucBase = h * eps22;
  }

  for (i=0; i<ndPosBnd2D[0]; i++)
  {
    k = uDep[j].nd;   
  
    //cout << k << "\n";
    
    uDep[j++].ucBase = w * eps11;
    uDep[j++].ucBase = w * eps12;
  }

  // negative boundary 

    // nothing to do here


  return;
}





void MicroCellWulf::getStressFromReactions(double *stress)
{
  // so far as small strains

  int i, j, n;
	
  for (i=0; i<4; i++)  stress[i] = 0.;

  n = nBndNd;
  
  for (i=0; i<n; i++)
  {
    j = bndNd[i];
    stress[0] -= reac(j,1) * x0(j,1);
    stress[1] -= reac(j,2) * x0(j,2);
    stress[2] -= reac(j,1) * x0(j,2);
    stress[3] -= reac(j,2) * x0(j,1);
  }
    
  int n2 = uDep[2].nd,
      n3 = uDep[4].nd,
      n2x = (n2-1)*ndm, 
      n2y = n2x+1,
      n3x = (n3-1)*ndm,
      n3y = n3x+1;
  
  double xc[2] = {0., 0.}, stressNorm,
	 V0 = triangleArea2D(xc,&(x0(bndNd[nBndNd-1],1)),&(x0(bndNd[0],1))),
         *X0 = x0.x,
         *U  = u.x,
	 fact  = 1. / (X0[n3x] * X0[n2y] - X0[n2x] * X0[n3y]),
         eps11 = (U[n3x] * X0[n2y] - U[n2x] * X0[n3y]) * fact,
         eps22 = (U[n2y] * X0[n3x] - U[n3y] * X0[n2x]) * fact,
	 eps12 = (U[n3y] * X0[n2y] - U[n2y] * X0[n3y]) * fact,
	 strainNorm = sqrt(eps11*eps11+eps22*eps22+4.*eps12*eps12);

  for (i=1; i<nBndNd; i++) 

    V0 += triangleArea2D(xc,&(x0(bndNd[i-1],1)),&(x0(bndNd[i],1)));
 
  //cout << " volume = " << V0 << "\n\n";
  
  for (i=0; i<4; i++) stress[i] /= V0;
	  
  
  // so far, just print stress
  
  COUT << "sig_11 = " << stress[0] << "\n";
  COUT << "sig_22 = " << stress[1] << "\n";
  COUT << "sig_12 = " << stress[2] << "\n";
  COUT << "sig_21 = " << stress[3] << "\n";
  stressNorm = sqrt(stress[0]*stress[0]+stress[1]*stress[1]+4.*stress[2]*stress[2]);
  COUT << "stressNorm = " << stressNorm << "\n";
  COUT << "strainNorm = " << strainNorm << "\n\n";

  char tmp[30];
  
  sprintf(tmp,"%12.5g %12.5g\n",strainNorm,stressNorm);
  
  ofstream file;
  
  file.open("Maziar.dat",ios_base::app);

  file << tmp;    

  file.close();
  
  return;
}

