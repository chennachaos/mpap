
#include <iostream>

#include "DomainTypeEnum.h"
#include "VortexSheet.h"
#include "DomainTree.h"
#include "DomainType.h"
#include "FunctionsProgram.h"
#include "FunctionsEssGrp.h"
#include "Plot.h"
#include "DataBlockTemplate.h"
#include "MathBasic.h"
#include "MathGeom.h"
#include "MyStringList.h"
#include "FunctionsSupport.h"
#include "FunctionsProgram.h"


extern DomainTree domain;
extern Plot plot;



using namespace std;



VortexSheet::VortexSheet(void)                       
{                                                  
  // add new type
  
  DomainType *vort = domain.newType(VORTEXSHEET,ROOTDOMAIN);

  if (vort == NULL) return;  // domain type exists already

  vort->key.addNew("coordinates","control");

  ndm = 2;

  v = -2.;

  return;
}




	                                          
VortexSheet::~VortexSheet(void)                     
{         
  return;
}





void VortexSheet::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString tmpl, *word;
 
  char tmp[30], fct[] = "VortexSheet::readInputData";

  int nw, i;

  Vector<double> dTmp;
  Vector<int>    iTmp, lTmp;
  MyStringList   sTmp;

  DataBlockTemplate t1, t2;

  switch (domain[VORTEXSHEET].key.whichBegins(line))
  {
    case  0: cout << "     VORTEXSHEET: reading coordinates ...\n\n";
	   
	     sprintf(tmp,"123 %df",ndm);  
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	     
             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
		     prgError(2,fct,"data error in 'coordinates'!");
	     
             nSheet = dTmp.dim() / 2;

             x.setDim(nSheet,ndm,true);

             for (i=0; i<nSheet+nSheet; i++) x.x[i] = dTmp[i];

             break;

    case  1: cout << "     VORTESHEET: reading control ...\n\n";
	   
             if (v > -1) prgError(1,fct,"'control' has already been read!");
	     
	     line.getNextLine(Ifile);
	     
	     nw = line.split(&word);
            
	     if (nw < 4)                    prgError(1,fct,"input error in 'control'!");
		                           
             if (!word[0].toInt(&iSep))     prgError(5,fct,"input error in 'control'!");

             if (!word[1].toDbl(&v))        prgError(2,fct,"input error in 'control'!");
		                           
             if (!word[2].toDbl(&alpha))    prgError(4,fct,"input error in 'control'!");

             if (!word[3].toDbl(&rho))      prgError(3,fct,"input error in 'control'!");
                                           
             v     = abs(v);
             rho   = abs(rho);
             alpha *= 0.017453293;

             if (v*rho < 1.-10)             prgError(6,fct,"input error in 'control'!");

             if (iSep < 1 || iSep > nSheet) prgError(7,fct,"input error in 'control'!");

             for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     
       	     line.getNextLine(Ifile);

	     break;

    case -1: // go and inherit from DOMAIN
	     
	     this->Domain::readInputData(Ifile, line); 
	     
	     break;
  }
 
  return;
}








void VortexSheet::prepareInputData(void)
{
  // call ancestor function

  Domain::prepareInputData();

  
  cout << "     VORTEXSHEET: prepare input data ...\n\n";
  
  char fct[] = "VortexSheet::prepareInputData"; 
	  
  // ........

  return;
}








void VortexSheet::prepareInteractions(void)
{
  // go and inherit from ancestors

  Domain::prepareInteractions();

  cout << "     VORTEXSHEET: preparing interactions ...\n\n"; 
      
  return;
}










void VortexSheet::doForVortexSheet(void)
{
  char fct[] = "VortexSheet::doForVortexSheet";

  int i0, i1, i, j, j0, j1, k, n = nSheet, isw = 0;

  double *X = x.x, xmx[2], xmn[2], xc, dx;

  // find min/max of geometry

  xmn[0] = X[0];  xmn[1] = X[1];
  xmx[0] = X[0];  xmx[1] = X[1];

  for (i=1; i<nSheet; i++) 
  {
    xmn[0] = min(xmn[0],X[i+i  ]);
    xmn[1] = min(xmn[1],X[i+i+1]);
    xmx[0] = max(xmx[0],X[i+i  ]);
    xmx[1] = max(xmx[1],X[i+i+1]);
  }

  xc = (xmx[0]+xmn[0]) *.5; dx = xmx[0]-xmn[0]; xmn[0] = xc-.65*dx; xmx[0] = xc+.65*dx;
  xc = (xmx[1]+xmn[1]) *.5; dx = xmx[1]-xmn[1]; xmn[1] = xc-.65*dx; xmx[1] = xc+.65*dx;

  // plot geometry

  plot.setColour(1);
  
  plot.fit(xmn,xmx,2);
  
  for (i0=0; i0<n; i0++)
  {
    i1 = i0 + 1; if (i1 == n) i1 = 0;

    plot.line(X+i0+i0,X+i1+i1);
  }

  // calculate vortex distribution

  double *rhs  = new double [n],
         *mtx  = new double [n*n],
         *gam  = new double [n],
         *Cp   = new double [n],
         li, lj, dxi, dyi, dxj, dyj, 
         xi, xj, yi, yj, fact, dy, sgn,
         sna = sin(alpha), csa = cos(alpha),
         c = 0.5 / 3.14159265358979323846,
         //wgp[3] = {5./18.,4./9.,5./18.}, 
         //N[3]   = {.5-sqrt(3./20.),.5,.5+sqrt(3./20.)};
         wgp[5] = {.5*0.23692689,.5*0.47862867,.5*0.56888889,.5*0.47862867,.5*0.23692689}, 
         N[5]   = {.5-.5*0.90617985,.5-.5*0.53846931,.5,.5+.5*0.53846931,.5+.5*0.90617985};

  int    *pos = new int[n], ngp = 5;

  for (i=0; i<n; i++)   rhs[i] = 0.;
  for (i=0; i<n*n; i++) mtx[i] = 0.; 

  for (i0=0; i0<n; i0++)
  {
    i1 = i0 + 1; if (i1 == n) i1 = 0;

    dxi = X[i1+i1]   - X[i0+i0  ];
    dyi = X[i1+i1+1] - X[i0+i0+1];

    li = sqrt(dxi*dxi+dyi*dyi);

    rhs[i0] = v * (dyi*csa-dxi*sna)/li;

    xi = .5 * (X[i1+i1  ]+X[i0+i0  ]);
    yi = .5 * (X[i1+i1+1]+X[i0+i0+1]);
  
    for (j0=0; j0<n; j0++)
    {
      j1 = j0 + 1; if (j1 == n) j1 = 0;
           
      if (j0 + 1 == iSep) sgn = -1.; else sgn = 1.;

      if (j0 == i0)
      {
        mtx[j0*n+i0] -= c * sgn;
        mtx[j1*n+i0] += c;
      }
      else
      {
        dxj = X[j1+j1  ] - X[j0+j0  ];
        dyj = X[j1+j1+1] - X[j0+j0+1];
      
        lj = sqrt(dxj*dxj+dyj*dyj);
      
        for (k=0; k<ngp; k++)
        {
          xj = (1.-N[k]) * X[j0+j0  ] + N[k] * X[j1+j1  ];
          yj = (1.-N[k]) * X[j0+j0+1] + N[k] * X[j1+j1+1];
      
          dx = xi - xj;
          dy = yi - yj;
          
          fact = (dx * dxi + dy * dyi) * c * lj / (li * (dx*dx + dy*dy)) * wgp[k];
      
          mtx[j0*n+i0] -= (1.-N[k]) * fact * sgn;
          mtx[j1*n+i0] -=   N[k]    * fact;
        }
      }
    }
  }

  //prgPrintSimpleMatrix(mtx, n,n,7,4,false,0,false);
  //for (i=0; i<n; i++) cout << rhs[i] << "\n";

  decomplr_matrix_(mtx,pos,&n,&isw);

  solve_matrix_(mtx,pos,rhs,gam,&n);

  // calculate tangential velocity, 
  // pressure coefficient,
  // total circulation,
  // lift force (Kutta-Joukowski)
  // lift and drag forces (from pressure)

  cout << "                x            y          theta         vt          -Cp\n";
  cout << "   ----+------------+------------+------------+------------+--------------\n";

  double vt, circ = 0., lift = 0., drag = 0.,
         b = 1.e+10, 
         pxmx = plot.x0Act[0] + .95*plot.dAct[0],
         pxmn = plot.x0Act[0] + .05*plot.dAct[0],
         pymx = plot.x0Act[1] + .95*plot.dAct[1],
         pymn = plot.x0Act[1] + .05*plot.dAct[1];

  MyString pathAndFile;
 
  char tmp[200];
 
  pathAndFile.append(files.projDir).append(SLASH).append(files.Ofile);

  files.Oout.open(pathAndFile.asCharArray());

  if (!files.Oout.is_open()) prgWarning(1,fct,"failed to open output file!");
  
  for (i0=0; i0<n; i0++)
  {
    i1 = i0 + 1; if (i1 == n) i1 = 0;

    xi = .5 * (X[i0+i0  ] + X[i1+i1]  );
    yi = .5 * (X[i0+i0+1] + X[i1+i1+1]);

    c  = atan2(yi,xi); if (c < 0.) c += 2*3.14159265358979323846;

    dxi = X[i1+i1]   - X[i0+i0  ];
    dyi = X[i1+i1+1] - X[i0+i0+1];

    li = sqrt(dxi*dxi+dyi*dyi);

    if (i0 + 1 == iSep) vt = -.5 * (gam[i1]-gam[i0]); else vt = -.5 * (gam[i1]+gam[i0]);

    circ -= vt * li;

    Cp[i0] = 1. - vt*vt/(v*v);

    sprintf(tmp,"    %3d %12g %12g %12g %12g %12g\n",i0+1,xi,yi,c,vt,-Cp[i0]);

    if (files.Oout.is_open()) files.Oout << tmp;

    cout << tmp;

    fact = .5 * Cp[i0] * rho * v * v;
    lift += fact * (dxi*csa+dyi*sna);
    drag += fact * (dxi*sna-dyi*csa);
 
    // calc max scaling factor for drawing Cp

    xj = xi - dyi/li * b * Cp[i0];
    if      (xj > pxmx) b = -(pxmx - xi) / dyi * li / Cp[i0];
    else if (xj < pxmn) b = -(pxmn - xi) / dyi * li / Cp[i0];

    yj = yi + dxi/li * b * Cp[i0];
    if      (yj > pymx) b = (pymx - yi) / dxi * li / Cp[i0];
    else if (yj < pymn) b = (pymn - yi) / dxi * li / Cp[i0];
  }
 
  cout << "\n";
  cout << "   circulation     = " << circ << "\n";
  cout << "   lift force (KJ) = " << circ * rho * v << "\n";
  cout << "   lift force (Cp) = " << lift << "\n";
  cout << "   drag force (Cp) = " << drag << "\n\n";

  files.Oout.close();

  // draw Cp distribution

  plot.setColour(0);

  double x0[2], x1[2];

  for (i0=0; i0<n; i0++)
  {
    i1 = i0 + 1; if (i1 == n) i1 = 0;

    x0[0] = .5 * (X[i0+i0  ] + X[i1+i1]  );
    x0[1] = .5 * (X[i0+i0+1] + X[i1+i1+1]);

    dxi = X[i1+i1]   - X[i0+i0  ];
    dyi = X[i1+i1+1] - X[i0+i0+1];

    li = sqrt(dxi*dxi+dyi*dyi);

    x1[0] = x0[0] - dyi/li * b * Cp[i0];
    x1[1] = x0[1] + dxi/li * b * Cp[i0];
   
    plot.line(x0,x1);
  }
 
  // release memory

  delete [] rhs;
  delete [] mtx;
  delete [] pos;
  delete [] gam;
  delete [] Cp;

  //double dth = 2.*3.1415927 / 16.;

  //for (i=0; i<16; i++) { cout << i+1 << " " << cos(i*dth) << " " << sin(i*dth) << "\n"; }

  return;
}


