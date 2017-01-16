
#include <iostream>

#include "DomainTypeEnum.h"
#include "LiftingLine.h"
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
#include "Aerofoil.h"


extern DomainTree domain;
extern Plot plot;



using namespace std;



LiftingLine::LiftingLine(void)                       
{                                                  
  // add new type
  
  DomainType *lift = domain.newType(LIFTINGLINE,ROOTDOMAIN);

  if (lift == NULL) return;  // domain type exists already

  lift->key.addNew("wing data","elliptic wing","control");

  ndm = 2;

  numsc = 0;

  a0 = 0.;

  return;
}




	                                          
LiftingLine::~LiftingLine(void)                     
{         
  return;
}





void LiftingLine::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString tmpl, *word;
 
  char tmp[30], fct[] = "LiftingLine::readInputData",
       *aerofoilName[] = AEROFOIL;

  int nw, i, j, n;

  Vector<double> dTmp;
  Vector<int>    iTmp, lTmp;
  MyStringList   sTmp;

  DataBlockTemplate t1, t2;

  switch (domain[LIFTINGLINE].key.whichBegins(line))
  {
    case  0: cout << "     LIFTINGLINE: reading wing data ...\n\n";

             if (numsc > 0) prgError(2,fct,"wing already defined?!");

             sprintf(tmp,"123 3f 1l");
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	     
             t1.initialise(tmpl,aerofoilName);
             t2.initialise(tmp,aerofoilName);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
               prgError(2,fct,"data error in 'wing data'!");
	     
             numsc = dTmp.dim() / 3;

             y.setDim(numsc);
             l.setDim(numsc);
             c.setDim(numsc);
             da.setDim(numsc);
             aerofoil.setDim(numsc);

             for (i=0; i<numsc; i++)
             {
               y[i] = dTmp[i*3+0];
               c[i] = dTmp[i*3+1];
               da[i] = dTmp[i*3+2];
               aerofoil[i] = lTmp[i];
               if (aerofoil[i] < 0) prgError(1,fct,"unknown aerofoil in 'wing data'!");
             }
             l[0] = y[0]; 
             for (i=1; i<numsc; i++)
             {
               l[i] = y[i] - y[i-1];
               if (l[i] <= 0.) prgError(1,fct,"invalid span coordinates in 'wing data'!");
             }

             break;

    case  1: cout << "     LIFTINGLINE: reading elliptic wing ...\n\n";

             if (numsc > 0) prgError(1,fct,"wing already defined?!");

	     line.getNextLine(Ifile);
	     
	     nw = line.split(&word);
            
	     if (nw != 4) prgError(1,fct,"input error in 'elliptic wing'!");
          
             dTmp.free();

             if (!word[0].toDbl(dTmp.append())) prgError(1,fct,"input error in 'elliptic wing'!");
             if (!word[1].toDbl(dTmp.append())) prgError(2,fct,"input error in 'elliptic wing'!");
             if (!word[3].toInt(&numsc))        prgError(5,fct,"input error in 'elliptic wing'!");

             j = word[2].which(aerofoilName);

             for (i=0; i<nw; i++) word[i].free(); delete [] word;

             if (j < 0) prgError(6,fct,"input error in 'elliptic wing'!");

       	     line.getNextLine(Ifile);

               y.setDim(numsc);
               l.setDim(numsc);
               c.setDim(numsc);
              da.setDim(numsc); da.zero();

             aerofoil.setDim(numsc);

             dTmp.append(0.);

             for (i=0; i<numsc; i++)
             {
                dTmp[2] += .5 * pi / ((double)numsc);

                y[i] = dTmp[0] * sin(dTmp[2]) * .5;
                c[i] = dTmp[1] * cos(dTmp[2] - .25 * pi / ((double)numsc));

                aerofoil[i] = j;
             }
             l[0] = y[0]; for (i=1; i<numsc; i++) l[i] = y[i] - y[i-1];

             break;

    case  2: cout << "     LIFTINGLINE: reading control ...\n\n";
	   
	     line.getNextLine(Ifile);
	     
	     nw = line.split(&word);
            
	     if (nw != 1) prgError(1,fct,"input error in 'control'!");
		                           
             if (!word[0].toDbl(&a0)) prgError(2,fct,"input error in 'control'!");

             for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     
       	     line.getNextLine(Ifile);

	     break;

    case -1: // go and inherit from DOMAIN
	     
	     this->Domain::readInputData(Ifile, line); 
	     
	     break;
  }
 
  return;
}








void LiftingLine::prepareInputData(void)
{
  // call ancestor function

  Domain::prepareInputData();

  
  cout << "     LIFTINGLINE: prepare input data ...\n\n";
  
  char fct[] = "LiftingLine::prepareInputData"; 

  if (numsc == 0) prgError(2,fct,"input coordinates and wing data!");

  // ........

  return;
}








void LiftingLine::prepareInteractions(void)
{
  // go and inherit from ancestors

  Domain::prepareInteractions();

  cout << "     LIFTINGLINE: preparing interactions ...\n\n"; 
      
  return;
}










void LiftingLine::doForLiftingLine(void)
{
  char fct[] = "LiftingLine::doForLiftingLine";

  int i;

  VectorArray<double> mtx, invMtx, rhs, dccl, ccl, ccdi, ccdp;

  VectorArray<int> pos;

  double CL, CDi, CD;

  mtx.setDim(numsc*numsc);
  invMtx.setDim(numsc*numsc);
  dccl.setDim(numsc);
  rhs.setDim(numsc);
  pos.setDim(numsc);

  ccl.setDim(numsc);
  ccdi.setDim(numsc);
  ccdp.setDim(numsc);

  S = 0.; for (i=0; i<numsc; i++) S += c[i] * l[i]; S += S;

  analyseWing(mtx.x,rhs.x,dccl.x,pos.x,ccl.x,ccdi.x,ccdp.x,CL,CDi,CD);

  printResults(ccl.x,ccdi.x,ccdp.x,CL,CDi,CD);

  return;

  if (!noGUI) plotResults(ccl.x,ccdi.x,ccdp.x,CL,CDi,CD);

  diffTest(mtx.x,rhs.x,dccl.x,pos.x,ccl.x,ccdi.x,ccdp.x);

  optimiseTwist(mtx.x,rhs.x,dccl.x,pos.x,ccl.x,ccdi.x,ccdp.x);

  return;
}











void LiftingLine::analyseWing(double *mtx, double *rhs, double *dccl, int *pos,
                              double *ccl, double *ccdi, double *ccdp,
                              double &CL, double &CDi, double &CD)
{
  int isw = 0;

  calculateMatrixAndRHS(mtx,rhs);

  decomplr_matrix_(mtx,pos,&numsc,&isw);

  solve_matrix_(mtx,pos,rhs,dccl,&numsc);

  calculateLiftAndDrag(dccl,ccl,ccdi,ccdp,CL,CDi,CD);

  return;
}











void LiftingLine::calculateMatrixAndRHS(double *mtx, double *rhs)
{
  char fct[] = "LiftingLine::calculateMatrixAndRHS";

  int i, j, n = numsc;

  double fact, fact2, yi2,
         cla, aL0,
         liftSlope[] = LIFT_SLOPE,
         zeroLiftAngle[] = ZERO_LIFT_ANGLE;

  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++) mtx[j*n+i] = 0.;

    for (j=i+1; j<n; j++) 
    {
      mtx[(j-1)*n+i] -= .5 * l[j];

      mtx[j*n+i] -= .5 * l[j];
    }

    if (i>0) mtx[(i-1)*n+i] -= .125 * l[i];

    mtx[i*n+i] -= .375 * l[i];
  }

  for (i=0; i<n; i++)
  {
    cla = liftSlope[aerofoil[i]];
    aL0 = zeroLiftAngle[aerofoil[i]];

    fact = .125 * c[i] * cla / pi * 180. / pi;

    yi2 = y[i] - .5 * l[i];

    for (j=0; j<n; j++)
    {
      if (i != j)
      {
        fact2 = log((y[j] - yi2) / (y[j] - l[j] - yi2)) / l[j];

        if (j>0)

          mtx[(j-1)*n+i] += fact * (1. + fact2 * (yi2 - y[j]));

        mtx[j*n+i] -= fact * (1. + fact2 * (yi2 - y[j] + l[j]));
      }
      fact2 = log((y[j] - l[j] + yi2) / (y[j] + yi2)) / l[j];

      if (j>0)

        mtx[(j-1)*n+i] += fact * (1. + fact2 * (yi2 + y[j]));

      mtx[j*n+i] -= fact * (1. + fact2 * (yi2 + y[j] - l[j]));
    }

    if (i>0) mtx[(i-1)*n+i] += fact;

    mtx[i*n+i] -= fact;

    rhs[i] = c[i] * cla * (a0 + da[i] - aL0);
  }
  return;
}











void LiftingLine::calculateLiftAndDrag(double *dccl, double *ccl, double *ccdi, double *ccdp,
                                       double &CL, double &CDi, double &CD)
{
  int nDP, dragPolarJpt[] = DRAG_POLAR_JPT, i, j, n = numsc;

  double *DP, cla, aL0,
         dragPolar[]     = DRAG_POLAR,
         liftSlope[]     = LIFT_SLOPE,
         zeroLiftAngle[] = ZERO_LIFT_ANGLE;

  ccl[n-1] = - (.375 * dccl[n-1] + .125 * dccl[n-2]) * l[n-1];

  for (i=n-2; i>0; i--) ccl[i] = ccl[i+1] - (.125 * dccl[i+1] + .375 * dccl[ i ]) * l[i+1]
                                          - (.375 * dccl[ i ] + .125 * dccl[i-1]) * l[ i ];

  ccl[0] = ccl[1] - (.125 * dccl[1] + .375 * dccl[ 0 ]) * l[1]
                  -  .375 * dccl[0]                     * l[0];

  for (i=0; i<n; i++)
  {
    cla =     liftSlope[aerofoil[i]];
    aL0 = zeroLiftAngle[aerofoil[i]];

    ccdi[i] = (a0 + da[i] - aL0 - ccl[i] / (cla * c[i])) * pi / 180. * ccl[i];
  }

  for (i=0; i<n; i++)
  {
    nDP = dragPolarJpt[aerofoil[i]];

    DP = dragPolar + nDP + nDP;

    nDP = dragPolarJpt[aerofoil[i]+1] - nDP;

    j = 0; while (j < nDP) if (DP[j+j+2]*c[i] > ccl[i]) break; else j++;

    if (j == nDP) j--;

    ccdp[i] = DP[j+j+1]*c[i] + (DP[j+j+3]-DP[j+j+1]) / (DP[j+j+2]-DP[j+j]) * (ccl[i]-DP[j+j]*c[i]);
  }
  
  CL  = 0.;
  CDi = 0.;
  CD  = 0.;

  for (i=0; i<n; i++)
  {
    CL  += ccl [i] * l[i];
    CDi += ccdi[i] * l[i];
    CD  += ccdp[i] * l[i];
  }

  CL += (dccl[0]) * l[0] * l[0] / 24.;

  for (i=1; i<n; i++) CL += (dccl[i]-dccl[i-1]) * l[i] * l[i] / 24.;

  CL  /= S * .5;
  CDi /= S * .5;
  CD  /= S * .5;
  CD  += CDi;

  return;
}












void LiftingLine::printResults(double *ccl, double *ccdi, double *ccdp,
                               double &CL, double &CDi, double &CD)
{
  int i, n = numsc;

  COUT << "total wing area:                 S  = " << S << "\n";
  COUT << "wing aspect ratio:              AR  = " << 4. * y[n-1]*y[n-1] / S << "\n";
  COUT << "root anlge of attack:        alpha0 = " << a0 << "\n";
  COUT << "-------------------------------------------------\n";
  COUT << "FROM LIFTING LINE THEORY:\n";
  COUT << "total lift coefficient:         CL  = " << CL << "\n";
  COUT << "induced drag coefficient:       CDi = " << CDi << "\n";
  COUT << "-------------------------------------------------\n";
  COUT << "WITH SECTION DRAG POLARS:\n";
  COUT << "skin friction drag coefficient: CDp = " << CD - CDi << "\n";
  COUT << "total drag coefficient:         CD  = " << CD << "\n";
  COUT << "lift/drag ratio:            CL / CD = " << CL / CD << "\n\n";

  COUT << "lift and drag distributions:\n";
  COUT << "        y         c(y) x cl(y)   c(y) x cdi(y)   c(y) x cd(y)\n";
  COUT << "---------------------------------------------------------------\n";
  for (i=0; i<numsc; i++)
  printf("           %12f   %12f   %12f   %12f\n",y[i]-.5*l[i],ccl[i],ccdi[i],ccdi[i]+ccdp[i]);
  cout << "\n";

  return;
}










void LiftingLine::plotResults(double *ccl, double *ccdi, double *ccdp,
                              double &CL, double &CDp, double &CDi)
{
/*  int i, n = numsc;

  // find maximum values

  double mxt = 0., mxf = 0., mxl = 0., mxd = 0., t0, f0, l0, d0, x0[2], x1[2], vsp, dgh, dgw;

  for (i=0; i<n; i++)
  {
    mxt = max(mxt,c[i]);
    mxf = max(mxf,c[i]*abs(sin((a0+da[i])*pi/180.)));
    mxl = max(mxl,ccl[i]);
    mxd = max(mxd,ccdp[i]+ccdi[i]);
  }

  vsp = 1.;
  dgh = 2.;
  dgw = 8.;

  d0  = 0.;
  l0  = dgh + vsp;
  f0  = l0 + dgh + vsp + .5 * dgw * mxf / y[n];
  t0  = f0 + vsp + .5 * dgw * (mxf + mxt) / y[n];
 
  x0[0] = 0.;
  x0[1] = 0.;

  x1[0] = dgw;
  x1[1] = t0 + .5 * dgw * mxt / y[n]; 

  plot.fit(x0,x1,2,0.1);

  plot.setColour(0); 

  // plot lift

  x0[0] = 0.;
  x0[1] = l0;
  x1[0] = dgw;
  x1[1] = l0;

  plot.line(x0,x1);

  x1[0] = 0.;
  x1[1] = l0 + ccl[0] / mxl * dgh;

  plot.line(x0,x1);

  for (i=0; i<numnp; i++)
  {
    x0[0] = y[i] / y[n] * dgw;
    x0[1] = l0 + ccl[i] / mxl * dgh;

    plot.line(x0,x1);

    x1[0] = x0[0];
    x1[1] = x0[1];
  }

  // plot drag

  plot.setColour(2); 

  x0[0] = 0.;
  x0[1] = d0;
  x1[0] = dgw;
  x1[1] = d0;

  plot.line(x0,x1);

  x1[0] = 0.;
  x1[1] = d0 + (ccdp[0] + ccdi[0]) / mxd * dgh;

  plot.line(x0,x1);

  for (i=0; i<numnp; i++)
  {
    x0[0] = y[i] / y[n] * dgw;
    x0[1] = d0 + (ccdp[i] + ccdi[i]) / mxd * dgh;

    plot.line(x0,x1);

    x1[0] = x0[0];
    x1[1] = x0[1];
  }

  x0[1] = d0;

  plot.line(x0,x1);

  plot.setColour(3);

  x0[0] = 0.;
  x0[1] = d0;
  x1[0] = 0.;
  x1[1] = d0 + ccdi[0] / mxd * dgh;

  plot.line(x0,x1);

  for (i=0; i<numnp; i++)
  {
    x0[0] = y[i] / y[n] * dgw;
    x0[1] = d0 + ccdi[i] / mxd * dgh;

    plot.line(x0,x1);

    x1[0] = x0[0];
    x1[1] = x0[1];
  }

  // plot top view

  plot.setColour(1);

  x0[0] = 0.;
  x0[1] = t0 - .5 * c[0] * dgw / y[n];
  x1[0] = 0.;
  x1[1] = t0 + .5 * c[0] * dgw / y[n];

  plot.line(x0,x1);

  for (i=1; i<numnp; i++)
  {
    x0[0] = y[i] / y[n] * dgw;
    x0[1] = t0 + .5 * c[i] * dgw / y[n];

    plot.line(x0,x1);

    x1[0] = x0[0];
    x1[1] = x0[1];
  }

  x1[0] = 0.;
  x1[1] = t0 - .5 * c[0] * dgw / y[n];

  for (i=1; i<numnp; i++)
  {
    x0[0] = y[i] / y[n] * dgw;
    x0[1] = t0 - .5 * c[i] * dgw / y[n];

    plot.line(x0,x1);

    x1[0] = x0[0];
    x1[1] = x0[1];
  }

  x0[1] = t0 + .5 * c[n] * dgw / y[n];

  plot.line(x0,x1);

  // plot front view

  plot.setColour(1);

  x0[0] = 0.;
  x0[1] = f0 - .5 * c[0] * sin((a0+da[0])*pi/180.) * dgw / y[n];
  x1[0] = 0.;
  x1[1] = f0 + .5 * c[0] * sin((a0+da[0])*pi/180.) * dgw / y[n];

  plot.line(x0,x1);

  for (i=1; i<numnp; i++)
  {
    x0[0] = y[i] / y[n] * dgw;
    x0[1] = f0 + .5 * c[i] * sin((a0+da[i])*pi/180.) * dgw / y[n];

    plot.line(x0,x1);

    x1[0] = x0[0];
    x1[1] = x0[1];
  }

  x1[0] = 0.;
  x1[1] = f0 - .5 * c[0] * sin((a0+da[0])*pi/180.) * dgw / y[n];

  for (i=1; i<numnp; i++)
  {
    x0[0] = y[i] / y[n] * dgw;
    x0[1] = f0 - .5 * c[i] * sin((a0+da[i])*pi/180.) * dgw / y[n];

    plot.line(x0,x1);

    x1[0] = x0[0];
    x1[1] = x0[1];
  }

  x0[1] = f0 + .5 * c[n] * sin((a0+da[n])*pi/180.) * dgw / y[n];

  plot.line(x1,x0);
*/
  return;
}










void LiftingLine::calculateDerivatives(double *mtx, double *rhs, double *dccl, int *pos,
                                       double *ccl, double *ccdi, double *ccdp,
                                       double *dCL, double *dCDi, double *dCD, double *ddCDi)
{
  VectorArray<double> invMtx, tmp2;

  int nDP, dragPolarJpt[] = DRAG_POLAR_JPT,
      n = numsc, n2 = n * n, i, j, k;

  double *DP, dragPolar[] = DRAG_POLAR,
         cla, liftSlope[] = LIFT_SLOPE,
         aL0, zeroLiftAngle[] = ZERO_LIFT_ANGLE, fact;

  invMtx.setDim(n2);

  tmp2.setDim(n2);

  // calculate   invMtx = d dccl / d da

  inverse_matrix_(mtx,invMtx.x,pos,&n);

  for (j=0; j<n; j++)
  {
    for (i=0; i<n; i++)
    {
      cla = liftSlope[aerofoil[j]];

      invMtx[j*n+i] *= cla * c[j];
    }
  }

  // calculate   tmp2 = d ccl / d da

  for (j=0; j<n; j++) tmp2[j*n + n-1] = - (.375*invMtx[j*n+n-1] + .125*invMtx[j*n+n-2]) * l[n-1];

  for (i=n-2; i>0; i--)
 
    for (j=0; j<n; j++)

      tmp2[j*n+i] = tmp2[j*n+i+1] - (.125 * invMtx[j*n+i+1] + .375 * invMtx[j*n+ i ]) * l[i+1]
                                  - (.375 * invMtx[j*n+ i ] + .125 * invMtx[j*n+i-1]) * l[ i ];
  for (j=0; j<n; j++)

    tmp2[j*n+0] = tmp2[j*n+1] - (.125 * invMtx[j*n+1] + .375 * invMtx[j*n+ 0 ]) * l[1]
                              -  .375 * invMtx[j*n+0]                           * l[0];

  // calculate dCL, dCDi

  for (j=0; j<n; j++)
  {
    dCL[j] = 0.;

    for (i=0; i<n; i++) dCL[j] += tmp2[j*n+i] * l[i];

    dCL[j] += (invMtx[j*n]) * l[0] * l[0] / 24.;

    for (i=1; i<n; i++) dCL[j] += (invMtx[j*n+i]-invMtx[j*n+i-1]) * l[i] * l[i] / 24.;

    dCL[j]  /= S * .5;

    dCDi[j] = 0.;

    for (i=0; i<n; i++)
    {
      cla =     liftSlope[aerofoil[i]];
      aL0 = zeroLiftAngle[aerofoil[i]];

      dCDi[j] += (a0 + da[i] - aL0 - 2. * ccl[i] / (cla * c[i])) * pi / 180. * tmp2[j*n+i] * l[i];
    }
    dCDi[j] += pi / 180. * ccl[j] * l[j];

    dCD[j] = dCDi[j];

    dCDi[j] /= S * .5;
  }

  // calculate dCDp

  for (i=0; i<n; i++)
  {
    nDP = dragPolarJpt[aerofoil[i]];

    DP = dragPolar + nDP + nDP;

    nDP = dragPolarJpt[aerofoil[i]+1] - nDP;

    k = 0; while (k < nDP) if (DP[k+k+2]*c[i] > ccl[i]) break; else k++;

    if (k == nDP) k--;

    for (j=0; j<n; j++)

      dCD[j] += (DP[k+k+3]-DP[k+k+1]) / (DP[k+k+2]-DP[k+k]) * tmp2[j*n+i] * l[i];
  }

  for (j=0; j<n; j++) dCD[j] /= S * .5;

  if (ddCDi == NULL) return;

  // calculate ddCDi

  for (j=0; j<n; j++)
  {
    for (k=0; k<n; k++)
    {
      ddCDi[k*n+j] = 0.;

      for (i=0; i<n; i++)
      {
        cla = liftSlope[aerofoil[i]];

        fact = - tmp2[j*n+i] * tmp2[k*n+i] / (cla *c[i]) * 2.;

        if (i == j) fact += tmp2[k*n+i];
        if (i == k) fact += tmp2[j*n+i];

        ddCDi[k*n+j] += fact * pi / 180. * l[i] / S * 2.;
      }
    }
  }
  
  return;
}












void LiftingLine::diffTest(double *mtx, double *rhs, double *dccl, int *pos, 
                           double *ccl, double *ccdi, double *ccdp)
{
  int i, j, k;

  double ddd = .00001,
         dd[6] = {-3.*ddd, -2.*ddd, -ddd, +ddd, +2.*ddd, +3.*ddd },
         CDi[6], CL[6], CD[6];

  VectorArray<double> RCL, RCDi, RCD, dCL, dCDi, dCD, KCDi, KCD, ddCDi;

  // calculate numerical derivatives

  RCL.setDim(numsc);
  RCDi.setDim(numsc);
  RCD.setDim(numsc);
  KCDi.setDim(numsc*numsc);
  dCL.setDim(numsc*6);
  dCDi.setDim(numsc*6);
  dCD.setDim(numsc*6);
  ddCDi.setDim(numsc*numsc);

  // calculate numerical derivatives

  for (j=0; j<numsc; j++) // loop over columns (twist angles)
  {
    for (k=0; k<6; k++) // loop over perturbations
    {
      // apply pertubation

      da[j] += dd[k];

      // calculate CL and CD

      analyseWing(mtx,rhs,dccl,pos,ccl,ccdi,ccdp,CL[k],CDi[k],CD[k]);

      // remove pertubation
  	
      da[j] -= dd[k];

    }
    RCL[j] = -(CL[0] - 9.*CL[1] + 45.*CL[2] - 45.*CL[3] +  9.*CL[4] - CL[5]) / (60. * ddd);

    RCDi[j] = -(CDi[0] - 9.*CDi[1] + 45.*CDi[2] - 45.*CDi[3] +  9.*CDi[4] - CDi[5]) / (60. * ddd);

    RCD[j] = -(CD[0] - 9.*CD[1] + 45.*CD[2] - 45.*CD[3] +  9.*CD[4] - CD[5]) / (60. * ddd);
  }

  // calculate analytical derivatives

  calculateDerivatives(mtx,rhs,dccl,pos,ccl,ccdi,ccdp,dCL.x,dCDi.x,dCD.x,ddCDi.x);

  // compare

  prgCompareTwoSimpleMatrices(RCL.x,                         // matrix 1
		              dCL.x,                         // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      1,numsc,                       // matrix dimension
			      10,6,false,                    // format
			      0,                             // indentation
			      false,                         // interactive
			      false);                        // row/column numbers

  prgCompareTwoSimpleMatrices(RCDi.x,                        // matrix 1
		              dCDi.x,                        // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      1,numsc,                       // matrix dimension
			      10,6,false,                    // format
			      0,                             // indentation
			      false,                         // interactive
			      false);                        // row/column numbers

  prgCompareTwoSimpleMatrices(RCD.x,                         // matrix 1
		              dCD.x,                         // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      1,numsc,                       // matrix dimension
			      10,6,false,                    // format
			      0,                             // indentation
			      false,                         // interactive
			      false);                        // row/column numbers

  // calculate numerical second derivative

  for (j=0; j<numsc; j++)
  {
    for (k=0; k<6; k++) // loop over perturbations
    {
      // apply pertubation

      da[j] += dd[k];

      // calculate CL and CD

      analyseWing(mtx,rhs,dccl,pos,ccl,ccdi,ccdp,CL[k],CDi[k],CD[k]);

      calculateDerivatives(mtx,rhs,dccl,pos,ccl,ccdi,ccdp,
                           dCL.x+k*numsc,dCDi.x+k*numsc,dCD.x+k*numsc,NULL);

      // remove pertubation
  	
      da[j] -= dd[k];

    }
    for (i=0; i<numsc; i++)

      KCDi[j*numsc+i] = -(      dCDi[0*numsc+i]
                          -  9.*dCDi[1*numsc+i]
                          + 45.*dCDi[2*numsc+i]
                          - 45.*dCDi[3*numsc+i]
                          +  9.*dCDi[4*numsc+i]
                          -     dCDi[5*numsc+i]) / (60. * ddd);
  }

  // compare

  prgCompareTwoSimpleMatrices(KCDi.x,                        // matrix 1
		              ddCDi.x,                       // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      numsc,numsc,                   // matrix dimension
			      10,6,false,                    // format
			      0,                             // indentation
			      false,                         // interactive
			      false);                        // row/column numbers*/
  return;
}








void LiftingLine::optimiseTwist(double *mtx, double *rhs, double *dccl, int *pos, 
                                double *ccl, double *ccdi, double *ccdp)
{
  int i, j, n = numsc, m = n + 1, isw = 0, iter = 0;

  double fact = 1., tol = 1.e-5, lam = 0., CD, CDi, CL;

  VectorArray<double> R, K, d, dCL, dCDi, dCD, ddCDi;

  VectorArray<int> Pos;

  R.setDim(m);
  K.setDim(m*m);
  d.setDim(m);
  dCL.setDim(n); 
  dCDi.setDim(n); 
  dCD.setDim(n); 
  ddCDi.setDim(n*n); 
  Pos.setDim(m);

  while (true)
  {
    analyseWing(mtx,rhs,dccl,pos,ccl,ccdi,ccdp,CL,CDi,CD);

    calculateDerivatives(mtx,rhs,dccl,pos,ccl,ccdi,ccdp,dCL.x,dCDi.x,dCD.x,ddCDi.x);
 
    for (i=0; i<n; i++) R[i] = dCD[i] + lam * dCL[i];

    R[n] = CL - fact;

    cout << R.norm() << "\n";

    if (R.norm() < tol) break;
 
    if (++iter > 10) break;

    for (i=0; i<n; i++)
    {
      for (j=0; j<n; j++)

        K[j*m+i] = ddCDi[j*n+i];

      K[n*m+i] = dCL[i];
      K[i*m+n] = dCL[i];
    }
    K[n*m+n] = 0.;

    decomplr_matrix_(K.x,Pos.x,&m,&isw);

    solve_matrix_(K.x,Pos.x,R.x,d.x,&m);

    for (i=0; i<n; i++) da[i] -= d[i];

    lam -= d[n];
  }
  cout << CL << ", " << CDi << ", " << CD << "\n";

  cout << da << "\n";

  return;
}




