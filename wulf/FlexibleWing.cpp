
#include <iostream>

#include "DomainTypeEnum.h"
#include "FlexibleWing.h"
#include "DomainTree.h"
#include "DomainType.h"
#include "FunctionsProgram.h"
#include "FunctionsEssGrp.h"
#include "Plot.h"
#include "DataBlockTemplate.h"
#include "MathBasic.h"
#include "MathGeom.h"
#include "MathMatrix.h"
#include "MyStringList.h"
#include "FunctionsSupport.h"
#include "FunctionsProgram.h"
#include "Aerofoil.h"


extern DomainTree domain;
extern Plot plot;



using namespace std;



FlexibleWing::FlexibleWing(void)                       
{                                                  
  // add new type
  
  DomainType *flexWing = domain.newType(FLEXIBLEWING,ROOTDOMAIN);

  if (flexWing == NULL) return;  // domain type exists already

  flexWing->key.addNew("wing data","elliptic wing");

  ndm = 2;

  numsc = 0;

  solver = NULL;

  return;
}




	                                          
FlexibleWing::~FlexibleWing(void)                     
{         
  if (solver != NULL) delete solver;

  return;
}









void FlexibleWing::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString tmpl, *word;
 
  char tmp[30], fct[] = "FlexibleWing::readInputData",
       *aerofoilName[] = AEROFOIL;

  int nw, i, j, n, np;

  Vector<double> dTmp;
  Vector<int>    iTmp, lTmp;
  MyStringList   sTmp;

  DataBlockTemplate t1, t2;

  switch (domain[FLEXIBLEWING].key.whichBegins(line))
  {
    case  0: cout << "     FLEXIBLEWING: reading wing data ...\n\n";

             if (numsc > 0) prgError(2,fct,"wing already defined?!");

             np = 5;

             sprintf(tmp,"123 %df 1l",np);
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	     
             t1.initialise(tmpl,aerofoilName);
             t2.initialise(tmp,aerofoilName);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
               prgError(2,fct,"data error in 'wing data'!");
	     
             numsc = dTmp.dim() / np;

             y.setDim(numsc);
             l.setDim(numsc);
             c.setDim(numsc);
             da.setDim(numsc);
             e.setDim(numsc);
             GJ.setDim(numsc);
             aerofoil.setDim(numsc);

             for (i=0; i<numsc; i++)
             {
                y[i] = dTmp[i*np+0];
                c[i] = dTmp[i*np+1];
               da[i] = dTmp[i*np+2];
                e[i] = dTmp[i*np+3];
               GJ[i] = dTmp[i*np+4];
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

    case  1: cout << "     FLEXIBLEWING: reading elliptic wing ...\n\n";

             if (numsc > 0) prgError(1,fct,"wing already defined?!");

	     line.getNextLine(Ifile);
	     
	     nw = line.split(&word);
            
	     if (nw != np) prgError(1,fct,"input error in 'elliptic wing'!");
          
             dTmp.free();

             if (!word[0].toDbl(dTmp.append())) prgError(1,fct,"input error in 'elliptic wing'!");
             if (!word[1].toDbl(dTmp.append())) prgError(2,fct,"input error in 'elliptic wing'!");
             if (!word[1].toDbl(dTmp.append())) prgError(2,fct,"input error in 'elliptic wing'!");
             if (!word[3].toInt(&numsc))        prgError(5,fct,"input error in 'elliptic wing'!");

             j = word[np-1].which(aerofoilName);

             for (i=0; i<nw; i++) word[i].free(); delete [] word;

             if (j < 0) prgError(6,fct,"input error in 'elliptic wing'!");

       	     line.getNextLine(Ifile);

               y.setDim(numsc);
               l.setDim(numsc);
               c.setDim(numsc);
              da.setDim(numsc); da.zero();
               e.setDim(numsc);
              GJ.setDim(numsc);

             aerofoil.setDim(numsc);

             dTmp.append(0.);

             for (i=0; i<numsc; i++)
             {
                dTmp[np-1] += .5 * pi / ((double)numsc);

                 y[i] = dTmp[0] * sin(dTmp[np-1]) * .5;
                 c[i] = dTmp[1] * cos(dTmp[np-1] - .25 * pi / ((double)numsc));
 
                 e[i] = dTmp[2];
                GJ[i] = dTmp[3];

                aerofoil[i] = j;
             }
             l[0] = y[0]; for (i=1; i<numsc; i++) l[i] = y[i] - y[i-1];

             break;

    case -1: // go and inherit from DOMAIN
	     
	     this->Domain::readInputData(Ifile, line); 
	     
	     break;
  }
 
  return;
}











void FlexibleWing::prepareInputData(void)
{
  // call ancestor function

  Domain::prepareInputData();

  
  cout << "     FLEXIBLEWING: prepare input data ...\n\n";
  
  char fct[] = "FlexibleWing::prepareInputData"; 

  if (numsc == 0) prgError(2,fct,"input coordinates and wing data!");

  // ........

  return;
}








void FlexibleWing::prepareInteractions(void)
{
  // go and inherit from ancestors

  Domain::prepareInteractions();

  cout << "     FLEXIBLEWING: preparing interactions ...\n\n"; 
      
  return;
}










void FlexibleWing::doForFlexibleWing(int taskFlag, bool cplFlg, double *par, MyString &fileName)
{
  char strg[200], fct[] = "FlexibleWing::doForFlexibleWing";

  int i;

  std::ofstream DvFile; // output stream

  VectorArray<double> ccl, ccdi, ccd;

  VectorArray<int> pos;

  double a0, CL, CDi, CD,
         rho, W, vmin, vmax, dv, v, tol;

  S = 0.; for (i=0; i<numsc; i++) S += c[i] * l[i]; S += S;

  prepareForLiftingLine();

  prepareForFSI();

  ccl.setDim(numsc);
  ccdi.setDim(numsc);
  ccd.setDim(numsc);

  sol.setDim(numsc+numsc);
  r0.setDim(numsc+numsc);
  ra0.setDim(numsc+numsc);

  th = sol.x;

  ai = th + numsc;

  solver = new SolverMA41;

  cout << "\n";

  switch (taskFlag)
  {
    case 1: // analyse coupled/uncoupled system; input is (alpha0, q)

            if (cplFlg) COUT << "COUPLED SYSTEM ANALYSIS (alpha0,q):\n\n";

            else        COUT << "DECOUPLED SYSTEM ANALYSIS (a0,q):\n\n";

            a0 = par[0];
            q  = par[1];

            generateSystemMatrixAndPrepareRHS(cplFlg);

            generateFullRHS(sol.x,a0);

            solver->factoriseAndSolve(sol.x);

            calculateLiftAndDrag(CL,CDi,CD);

            printResults(a0,CL,CDi,CD);

            writeResults(fileName);

            //calculateDerivatives0();

            //diffTest(a0);

            break;

    case 2: // analyse coupled/uncoupled system; input is (CL, q)

            if (cplFlg) COUT << "COUPLED SYSTEM ANALYSIS (CL,q):\n\n";

            else        COUT << "DECOUPLED SYSTEM ANALYSIS (CL,q):\n\n";

            CL = par[0];
            q  = par[1];

            generateSystemMatrixAndPrepareRHS(cplFlg);

            solver->factorise();
 
            solveForAlpha0(a0,CL);

            generateFullRHS(sol.x,a0);

            solver->solve(sol.x);

            calculateLiftAndDrag(CL,CDi,CD);

            printResults(a0,CL,CDi,CD);

            writeResults(fileName);

            break;

    case 3: // optimise coupled/uncoupled system

            if (cplFlg) COUT << "COUPLED SYSTEM OPTIMISATION (CL,q,tol):\n\n";

            else        COUT << "DECOUPLED SYSTEM OPTIMISATION (CL,q,tol):\n\n";

            CL  = par[0];
            q   = par[1];
            tol = par[2];

            a0 = 0.;

            generateSystemMatrixAndPrepareRHS(cplFlg);

            solver->factorise();

            optimiseTwist(cplFlg,CL,a0,tol);
 
            generateFullRHS(sol.x,a0);

            solver->solve(sol.x);

            calculateLiftAndDrag(CL,CDi,CD);

            printResults(a0,CL,CDi,CD);

            writeResults(fileName);

            break;

    case 4: // generate D-v diagram

            if (cplFlg) COUT << "COUPLED SYSTEM D(v) ANALYSIS (rho,W,vmin,vmax,np):\n\n";

            else        COUT << "DECOUPLED SYSTEM D(v) ANALYSIS (rho,W,vmin,vmax,np):\n\n";

            if (!fileName)

             { COUT << "calculation aborted: specify out put file name!\n\n"; break; }

            rho  = par[0];
            W    = par[1];
            vmin = par[2];
            vmax = par[3];
            dv   = (vmax - vmin) / par[4];

            v = vmin;

            DvFile.open(fileName.asCharArray());

            DvFile << "%        v            q           CL          CDi"
                   << "          CDp           CD           Di           Dp"
                   << "            D          thT           a0\n"
                   << "%        1            2            3            4"
                   << "            5            6            7            8"
                   << "            9           10           11\n";

            while (v < vmax + .5*dv)
            {
              q  = .5 * rho * v * v;
              CL = W / (q*S);

              generateSystemMatrixAndPrepareRHS(cplFlg);

              solver->factorise();

              solveForAlpha0(a0,CL);

              generateFullRHS(sol.x,a0);

              solver->solve(sol.x);

              calculateLiftAndDrag(CL,CDi,CD);

              sprintf(strg,
                "%12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g",
                  v, q, CL, CDi, CD-CDi, CD, q*S*CDi, q*S*(CD-CDi), q*S*CD, th[numsc-1], a0);

              DvFile << strg << "\n";

              v += dv;
            }

            DvFile.close();

            COUT << "data file " << fileName << " generated.\n\n";

            break;

    case 5: // generate D-v diagram for optimised wing twist (rho,W,vmin,vmax,np)

            if (cplFlg) COUT << "COUPLED SYSTEM D(v) OPTIMISATION (rho,W,vmin,vmax,np):\n\n";

            else        COUT << "DECOUPLED SYSTEM D(v) OPTIMISATION (rho,W,vmin,vmax,np):\n\n";

            if (!fileName)

             { COUT << "calculation aborted: specify out put file name!\n\n"; break; }

            rho  = par[0];
            W    = par[1];
            vmin = par[2];
            vmax = par[3];
            dv   = (vmax - vmin) / par[4];
            tol  = par[5];

            v = vmin;

            DvFile.open(fileName.asCharArray());

            DvFile << "%        v            q           CL          CDi"
                   << "          CDp           CD           Di           Dp"
                   << "            D          thT           a0\n"
                   << "%        1            2            3            4"
                   << "            5            6            7            8"
                   << "            9           10           11\n";

            while (v < vmax + .5*dv)
            {
              q  = .5 * rho * v * v;
              CL = W / (q*S);

              generateSystemMatrixAndPrepareRHS(cplFlg);

              solver->factorise();

              optimiseTwist(cplFlg,CL,a0,tol);

              generateFullRHS(sol.x,a0);

              solver->solve(sol.x);

              calculateLiftAndDrag(CL,CDi,CD);

              sprintf(strg,
                "%12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g %12.5g",
                  v, q, CL, CDi, CD-CDi, CD, q*S*CDi, q*S*(CD-CDi), q*S*CD, th[numsc-1], a0);

              DvFile << strg << "\n";

              v += dv;
            }

            DvFile.close();

            COUT << "data file " << fileName << " generated.\n\n";

            break;
  }

  return;
}











void FlexibleWing::prepareForLiftingLine(void)
{
  char fct[] = "FlexibleWing::prepareForLiftingLine";

  int i, j, k, n = numsc, isw = 0;

  double fact = .125 / pi * 180. / pi, fact2, yi2;

  VectorArray<int> pos;

  VectorArray<double> K;

  pos.setDim(n);

  K.setDim(n*n);
  A.setDim(n*n);
  B.setDim(n*n);
  w.setDim(n);

  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++) K[j*n+i] = 0.;

    yi2 = y[i] - .5 * l[i];

    for (j=0; j<n; j++)
    {
      if (i != j)
      {
        fact2 = log((y[j] - yi2) / (y[j] - l[j] - yi2)) / l[j];

        if (j>0)

          K[(j-1)*n+i] += fact * (1. + fact2 * (yi2 - y[j]));

        K[j*n+i] -= fact * (1. + fact2 * (yi2 - y[j] + l[j]));
      }
      fact2 = log((y[j] - l[j] + yi2) / (y[j] + yi2)) / l[j];

      if (j>0)

        K[(j-1)*n+i] += fact * (1. + fact2 * (yi2 + y[j]));

      K[j*n+i] -= fact * (1. + fact2 * (yi2 + y[j] - l[j]));
    }

    if (i>0) K[(i-1)*n+i] += fact;

    K[i*n+i] -= fact;
  }

  //for (i=0; i<n; i++) { for (j=0; j<n; j++) cout << " " << K[j*n+i]; cout << "\n"; } cout << "\n";

  decomplr_matrix_(K.x,pos.x,&n,&isw);

  inverse_matrix_(K.x,B.x,pos.x,&n);

  for (i=0; i<n; i++)
  {
    for (j=max(0,i-1); j<n; j++) K[j*n+i] = 0.;

    for (j=i+1; j<n; j++) 
    {
      K[(j-1)*n+i] -= .5 * l[j];

      K[j*n+i] -= .5 * l[j];
    }

    if (i>0) K[(i-1)*n+i] -= .125 * l[i];

    K[i*n+i] -= .375 * l[i];
  }

  //for (i=0; i<n; i++) { for (j=0; j<n; j++) cout << " " << K[j*n+i]; cout << "\n"; } cout << "\n";

  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
    {
      A[j*n+i] = 0.;

      for (k=max(0,i-1); k<n; k++) A[j*n+i] += K[k*n+i] * B[j*n+k];
    }
  }

  //for (i=0; i<n; i++) { for (j=0; j<n; j++) cout << " " << A[j*n+i]; cout << "\n"; } cout << "\n";

  for (i=0; i<n; i++)
  {
    w[i] = (A[i*n] + B[i*n] * l[0] / 24.) * l[0];

    for (j=1; j<n; j++) w[i] += (A[i*n+j] + (B[i*n+j] - B[i*n+j-1]) * l[j] / 24.) * l[j];

    w[i] /= S * .5;
  }

  return;
}










void FlexibleWing::prepareForFSI(void)
{
  int i, n = numsc;

  double liftSlope[] = LIFT_SLOPE, cla; 

  f1.setDim(n);
  f2.setDim(n);

  for (i=0; i<n; i++)
  {
    cla = liftSlope[aerofoil[i]];

    f1[i] = GJ[i] * pi / (180. * l[i]);

    f2[i] = e[i] * c[i] * c[i] * l[i] * cla;
  }

  return;
}












void FlexibleWing::generateSystemMatrixAndPrepareRHS(bool coupledFlag)
{
  int i, j, n = numsc;

  double cla, liftSlope[]     = LIFT_SLOPE,
         aL0, zeroLiftAngle[] = ZERO_LIFT_ANGLE,
         fact;

  VectorArray<double> f3;

  MatrixSparse<double> mtx;

  rda.free();

  f3 = f2;

  for (i=0; i<n; i++) f3[i] *= q;

  // K_FSI_theta and R_FSI

  for (i=0; i<n; i++) { r0[i] = 0.; ra0[i] = 0.; } 

  for (i=0; i<n-1; i++)
  {
    mtx.append(i+1, i+1,  6.*(f1[i]+f1[i+1])-2.*(f3[i]+f3[i+1]) );
    mtx.append(i+1, i+2, -6.*       f1[i+1] -          f3[i+1]  );
    mtx.append(i+2, i+1, -6.*       f1[i+1] -          f3[i+1]  );

    aL0 = zeroLiftAngle[aerofoil[i+1]];

    fact = 3. * f3[i+1];

    ra0[i]   += fact;
    ra0[i+1] += fact;

    rda.append(i+1,i+2,fact);
    rda.append(i+2,i+2,fact);

    fact *= (- aL0);

    r0[i]   += fact;
    r0[i+1] += fact;
  }

  mtx.append(n, n, 6.*f1[n-1]-2.*f3[n-1] );

  aL0 = zeroLiftAngle[aerofoil[0]];

  ra0[0] += 3. * f3[0];

  rda.append(1,1, 3.*f3[0]);

  r0[0] += 3. * f3[0] * (- aL0);

  if (coupledFlag)
  {
    // K_FSI_ai

    for (i=0; i<n-1; i++)
    {
      fact = 3. * f3[i+1];

      mtx.append(i+1,i+n+2,fact);
      mtx.append(i+2,i+n+2,fact);
    }

    mtx.append(1,n+1, 3.*f3[0]);

    // K_LL_theta 

    mtx.append(n+1, 1, -.5);

    for (i=1; i<n; i++)
    {
      mtx.append(i+n+1, i  , -.5 );
      mtx.append(i+n+1, i+1, -.5 );
    }
  }

  // K_LL_ai and R_LL

  for (i=0; i<n; i++)
  {
    cla = liftSlope[aerofoil[i]];

    aL0 = zeroLiftAngle[aerofoil[i]];

    for (j=0; j<i; j++) mtx.append(n+i+1, n+j+1, A[j*n+i]/(c[i]*cla) );

    mtx.append(n+i+1, n+i+1, A[i*n+i]/(c[i]*cla) + 1.);

    for (j=i+1; j<n; j++) mtx.append(n+i+1, n+j+1, A[j*n+i]/(c[i]*cla) );

    ra0[n+i] = 1.;

    rda.append(i+n+1,i+1,1.);

    r0[n+i] = - aL0;
  }

  solver->mtx = mtx;

  //cout << mtx << "\n\n";

  //solver->printInfo();

  solver->currentStatus = PATTERN_OK;

  if (solver->initialise() != 0) return; 

  solver->currentStatus = ASSEMBLY_OK;

  return;
}











void FlexibleWing::generateFullRHS(double *rhs, double a0)
{
  int i, n = numsc;

  for (i=0; i<n+n; i++) rhs[i] = r0[i] + ra0[i] * a0;

  for (i=0; i<rda.x.n; i++) rhs[rda.row[i]-1] += rda.x[i] * da[rda.col[i]-1];

  return;
}










void FlexibleWing::solveForAlpha0(double &a0, double CL)
{
  int i, j, n = numsc;

  double fact = 0.;

  a0 = CL;

  for (i=0; i<n+n; i++) sol[i] = r0[i];

  for (i=0; i<rda.x.n; i++) sol[rda.row[i]-1] += rda.x[i] * da[rda.col[i]-1];

  if (sol.norm() > 1.e-12)
  {
    solver->solve(sol.x);

    for (i=0; i<n; i++) a0 -= w[i] * ai[i];
  }
  for (i=0; i<n+n; i++) sol[i] = ra0[i];

  solver->solve(sol.x);

  for (i=0; i<n; i++) fact += w[i] * ai[i];

  a0 /= fact;

  return;
}












void FlexibleWing::calculateLiftAndDrag(double &CL,  double &CDi,  double &CD)
{
  int i, j, nDP, dragPolarJpt[] = DRAG_POLAR_JPT, n = numsc;

  double *DP, dragPolar[] = DRAG_POLAR;

  ccl.setDim(n);

  ccdi.setDim(n);

  ccd.setDim(n);

  CL = 0.;

  for (i=0; i<n; i++) CL += w[i] * ai[i];

  CDi = 0.;

  CD  = 0.;

  for (i=0; i<n; i++) 
  {
    ccl[i] = 0.;

    for (j=0; j<n; j++) ccl[i] += A[j*n+i] * ai[j];

    ccdi[i] = ai[i] * ccl[i] * pi / 180.;

    CDi += l[i] * ccdi[i];

    nDP = dragPolarJpt[aerofoil[i]];

    DP = dragPolar + nDP + nDP;

    nDP = dragPolarJpt[aerofoil[i]+1] - nDP;

    j = 0; while (j < nDP) if (DP[j+j+2]*c[i] > ccl[i]) break; else j++;

    if (j == nDP) j--;

    ccd[i] = ccdi[i] + DP[j+j+1]*c[i]
                     + (DP[j+j+3]-DP[j+j+1]) / (DP[j+j+2]-DP[j+j]) * (ccl[i]-DP[j+j]*c[i]);

    CD += l[i] * ccd[i];
  }
  CDi *= 2. / S;
  CD  *= 2. / S;

  return;
}










void FlexibleWing::printResults(double a0, double CL, double CDi, double CD)
{
  int i, n = numsc;

  COUT << "total wing area:                 S  = " << S << "\n";
  COUT << "wing aspect ratio:              AR  = " << 4. * y[n-1]*y[n-1] / S << "\n";
  COUT << "root anlge of attack:        alpha0 = " << a0 << "\n";
  COUT << "-------------------------------------------------\n";
  COUT << "FROM FSI ANALYSIS:\n";
  COUT << "wing tip twist angle:       theta_T = " << th[n-1] << "\n";
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

  COUT << "     y        c x cl    c x cdi     c x cd     theta       ai         da\n";
  COUT << "----------------------------------------------------------------------------\n";
  printf("         %10f %10f %10f %10f %10f %10f %10f\n",
            y[0]-.5*l[0],ccl[0],ccdi[0],ccd[0],.5*th[0],ai[0],da[0]);
  for (i=1; i<numsc; i++)
    printf("         %10f %10f %10f %10f %10f %10f %10f\n",
              y[i]-.5*l[i],ccl[i],ccdi[i],ccd[i],.5*(th[i]+th[i-1]),ai[i],da[i]);

  cout << "\n";

  return;
}










void FlexibleWing::writeResults(MyString &fileName)
{
  if (!fileName) return;

  int i, n = numsc;

  char strg[200];

  std::ofstream outFile; // output stream

  outFile.open(fileName.asCharArray());

  outFile << "%        y          ccl         ccdi          ccd"
          << "        theta           ai           da\n%        1"
          << "            2            3            4            5            6            7\n";

  sprintf(strg,"%12g %12g %12g %12g %12g %12g %12g",
          0., ccl[0], ccdi[0], ccd[0], th[0], ai[0], da[0]);

  outFile << strg << "\n";

  for (i=0; i<n-1; i++)
  {
    sprintf(strg,"%12g %12g %12g %12g %12g %12g %12g",
            y[i], ccl[i+1], ccdi[i], ccd[i], th[i+1], ai[i], da[i]);

    outFile << strg << "\n";

    sprintf(strg,"%12g %12g %12g %12g %12g %12g %12g",
            y[i], ccl[i+1], ccdi[i+1], ccd[i+1], th[i+1], ai[i+1], da[i+1]);

    outFile << strg << "\n";
  }
  sprintf(strg,"%12g %12g %12g %12g %12g %12g %12g",
          y[n-1], 0., ccdi[n-1], ccd[n-1], th[n-1], ai[n-1], da[n-1]);

  outFile << strg << "\n";

  sprintf(strg,"%12g %12g %12g %12g %12g %12g %12g",
          y[n-1], 0., 0., 0., th[n-1], 0., 0.);

  outFile << strg << "\n";

  outFile.close();

  COUT << "data file " << fileName << " generated.\n\n";

  return;
}










void FlexibleWing::calculateDerivatives0(void)
{
  // calculate derivatives which depend only on GJ and q, not on ai, a0 or da

  int i, j, k, n = numsc;

  double fact;

  VectorArray<double> tmp;

  tmp.setDim(n*n);

  dsol.setDim(n*(n+n));

  dCL.setDim(n);

  ddCDi.setDim(n*n);

  // calculate  d {th,ai} / d da

  for (j=0; j<n; j++)
  {
    tmp.zero();

    for (i=0; i<rda.x.n; i++) if (rda.col[i]-1 == j) tmp[rda.row[i]-1] += rda.x[i];

    solver->solve(tmp.x);

    for (i=0; i<n+n; i++) dsol[j*(n+n)+i] = tmp[i];
  }

  // calculate  d CL / da

  for (i=0; i<n; i++)
  {
    dCL[i] = 0;

    for (j=0; j<n; j++) dCL[i] += w[j] * dsol[i*(n+n)+j+n];
  }

  // calculate  d2 CD / d da2

  fact = 2. / S / 180. * pi;

  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
    {
      tmp[j*n+i] = 0.;
    
      for (k=0; k<n; k++)
      {
        tmp[j*n+i] += (l[i] * A[k*n+i] + l[k] * A[i*n+k]) * fact * dsol[j*(n+n)+n+k];
      }
    }
  }
  for (i=0; i<n; i++)
  {
    for (j=0; j<n; j++)
    {
      ddCDi[j*n+i] = 0.;

      for (k=0; k<n; k++)
      {
        ddCDi[j*n+i] += tmp[i*n+k] * dsol[j*(n+n)+k+n];
      }
    }
  }

  return;
}










void FlexibleWing::calculateDerivatives1(void)
{
  // calculate derivatives which depend on GJ, q and on ai, a0, da

  int i, j, k, nDP, dragPolarJpt[] = DRAG_POLAR_JPT, n = numsc;

  double *DP, dragPolar[] = DRAG_POLAR, ccl, fact;

  VectorArray<double> tmp;

  tmp.setDim(n);

  dCDi.setDim(n);

  dCD.setDim(n);

  // calculate  d CDi / d da

  for (i=0; i<n; i++)
  {
    tmp[i] = 0.;

    for (j=0; j<n; j++) tmp[i] += ai[j] * (l[i] * A[j*n+i] + l[j] * A[i*n+j]);
  }
  fact = 2. / S * pi / 180.;

  for (i=0; i<n; i++) tmp[i] *= fact;

  for (i=0; i<n; i++)
  {
    dCDi[i] = 0;

    for (j=0; j<n; j++) dCDi[i] += tmp[j] * dsol[i*(n+n)+j+n];
  }

  // calculate  d CD / d da

  tmp.zero();

  for (i=0; i<n; i++) 
  {
    ccl = 0.; for (j=0; j<n; j++) ccl += A[j*n+i] * ai[j];

    nDP = dragPolarJpt[aerofoil[i]];

    DP = dragPolar + nDP + nDP;

    nDP = dragPolarJpt[aerofoil[i]+1] - nDP;

    j = 0; while (j < nDP) if (DP[j+j+2]*c[i] > ccl) break; else j++;

    if (j == nDP) j--;

    for (k=0; k<n; k++) tmp[k] += l[i] * (DP[j+j+3]-DP[j+j+1]) / (DP[j+j+2]-DP[j+j]) * A[k*n+i];
  }
  fact = 2. / S;

  for (i=0; i<n; i++) tmp[i] *= fact;

  for (i=0; i<n; i++)
  {
    dCD[i] = dCDi[i];

    for (j=0; j<n; j++) dCD[i] += tmp[j] * dsol[i*(n+n)+j+n];
  }

  return;
}












void FlexibleWing::diffTest(double a0)
{
  int i, j, k;

  double ddd = .00001,
         dd[6] = {-3.*ddd, -2.*ddd, -ddd, +ddd, +2.*ddd, +3.*ddd },
         CDi[6], CL[6], CD[6];

  VectorArray<double> RCL, RCDi, RCD, KCDi, dCDi6;

  // calculate numerical derivatives

  RCL.setDim(numsc);
  RCDi.setDim(numsc);
  RCD.setDim(numsc);
  KCDi.setDim(numsc*numsc);
  dCDi6.setDim(numsc*6);

  // calculate numerical derivatives

  for (j=0; j<numsc; j++) // loop over columns (twist angles)
  {
    for (k=0; k<6; k++) // loop over perturbations
    {
      // apply pertubation

      da[j] += dd[k];

      // calculate CL and CD

      generateFullRHS(sol.x,a0);

      solver->solve(sol.x);

      calculateLiftAndDrag(CL[k],CDi[k],CD[k]);

      // remove pertubation
  	
      da[j] -= dd[k];
    }
    RCL[j] = -(CL[0] - 9.*CL[1] + 45.*CL[2] - 45.*CL[3] +  9.*CL[4] - CL[5]) / (60. * ddd);

    RCDi[j] = -(CDi[0] - 9.*CDi[1] + 45.*CDi[2] - 45.*CDi[3] +  9.*CDi[4] - CDi[5]) / (60. * ddd);

    RCD[j] = -(CD[0] - 9.*CD[1] + 45.*CD[2] - 45.*CD[3] +  9.*CD[4] - CD[5]) / (60. * ddd);
  }

  // calculate analytical derivatives

  generateFullRHS(sol.x,a0);

  solver->solve(sol.x);

  calculateDerivatives1();

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

      // calculate dCDi

      generateFullRHS(sol.x,a0);

      solver->solve(sol.x);

      calculateDerivatives1();

      for (i=0; i<numsc; i++) dCDi6[k*numsc+i] = dCDi[i];

      // remove pertubation
  	
      da[j] -= dd[k];

    }
    for (i=0; i<numsc; i++)

      KCDi[j*numsc+i] = -(      dCDi6[0*numsc+i]
                          -  9.*dCDi6[1*numsc+i]
                          + 45.*dCDi6[2*numsc+i]
                          - 45.*dCDi6[3*numsc+i]
                          +  9.*dCDi6[4*numsc+i]
                          -     dCDi6[5*numsc+i]) / (60. * ddd);
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








void FlexibleWing::optimiseTwist(bool cplFlg, double targetCL, double &a0, double tol)
{
  int i, j, n = numsc, m = n + 1, isw = 0, iter = 0;

  double lam = 0., CL, CDi, CD;

  VectorArray<double> R, K, d;

  VectorArray<int> Pos;

  R.setDim(m);
  K.setDim(m*m);
  d.setDim(m);

  Pos.setDim(m);

  calculateDerivatives0();

  while (true)
  {
    generateFullRHS(sol.x,0.);

    solver->solve(sol.x);

    calculateLiftAndDrag(CL,CDi,CD);

    calculateDerivatives1();
 
    for (i=0; i<n; i++) R[i] = dCD[i] + lam * dCL[i];

    R[n] = CL - targetCL;

    //cout << R.norm() << "\n";

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
  if (R.norm() > tol)
  {
    //COUT << "CL = " << targetCL << "\n\n";

    prgWarning(1,"FlexibleWing::optimiseTwist","convergence problems!");
  }
  a0 = da[0];

  for (i=0; i<n; i++) da[i] -= a0;

  return;
}




