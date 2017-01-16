

#include "RigidBody.h"
#include "FunctionsProgram.h"
#include "MpapTime.h"
#include "Files.h"
#include "TimeFunction.h"


extern MpapTime mpapTime;
extern Files    files;
extern List<TimeFunction> timeFunction;


using namespace std;



RigidBody::RigidBody(void)                       
{                                                  

  
  // add new type
  
  DomainType *rigidBody = domain.newType(RIGIDBODY,ROOTDOMAIN);

  if (rigidBody == NULL) return;  // domain type exists already

  rigidBody->key.addNew("properties",
                        "prescribed displacements",
                        "forces",
                        "control",
                        "data output");
  return;
}




	                                          
RigidBody::~RigidBody(void)                     
{         
  return;
}







void RigidBody::setSolver(int slv, int *parm, bool)
{
  solverOK = true;

  return;
}







bool RigidBody::converged(void)
{
  if (rNorm < tol && localStiffnessError == 0) return true;

  return false;
}






bool RigidBody::diverging(double factor)
{
  if (rNormPrev > -0.1 && rNorm / rNormPrev > factor) return true;  

  if (localStiffnessError != 0) return true;

  if (prgNAN(rNorm)) return true;
  
  return false;
}

   





void RigidBody::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString tmpl, *word;
 
  char tmp[30], fct[] = "RigidBody::readInputData",
       *wrndTypeNames[] = { "u", "du", "ddu", "reac", "x", NULL };

  int nw, i, j, nn;
 
  Vector<double> dTmp;
  Vector<int>    iTmp;
  MyStringList   sTmp;

  switch (domain[RIGIDBODY].key.whichBegins(line))
  {
    case  0: cout << "     RIGIDBODY: reading properties ...\n\n";

             switch (ndm)
             {
               case  2: i = 3; break;
               case  3: i = 6; break;
               default: prgError(1,fct,"problem with spatial dimension!");
             }
             if (ndf != i) prgError(2,fct,"invalid space dimension!");

                x.setDim(ndf);
               xn.setDim(ndf);
               x0.setDim(ndf);
                u.setDim(ndf);
               un.setDim(ndf);
               u3.setDim(ndf);
               u4.setDim(ndf);
               u5.setDim(ndf);
               u6.setDim(ndf);
              frc.setDim(ndf);
              idu.setDim(ndf);
         frcTmFct.setDim(ndf);
	     prop.setDim(ndf*3);

             frc.zero();
             frcTmFct.zero();

             // idu

	     line.getNextLine(Ifile);
             nw = line.split(&word);
	     if (nw != ndf) prgError(1,fct,"error in 'properties'!");
	     
             j = 0;
	     for (i=0; i<ndf; i++) 
             {
	       if (!word[i].toInt(&(idu.x[i]))) prgError(2,fct,"error in 'properties'!");
               if (idu[i] == 0) idu[i] = ++j; else idu[i] = 0;
             }
             for (i=0; i<nw; i++) word[i].free(); delete [] word;

             // physical properties

             for (j=0; j<3; j++)
	     {
               line.getNextLine(Ifile);
               nw = line.split(&word);
	       if (nw != ndf) prgError(3,fct,"error in 'properties'!");

	       for (i=0; i<ndf; i++) 
	         if (!word[i].toDbl(&(prop.x[j*ndf+i]))) prgError(4,fct,"error in 'properties'!");
	
               for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     }

             // stress free configuration (optional)

             line.getNextLine(Ifile);
             nw = line.split(&word);
             if (nw != ndf) goto stop0;
	     for (i=0; i<ndf; i++) if (!word[i].toDbl(&(dTmp[i]),false)) goto stop0;
             for (i=0; i<ndf; i++) { x0[i] = dTmp[i]; x[i] = dTmp[i]; word[i].free(); }
             delete [] word; dTmp.free();

             // initial displacements (optional)

             line.getNextLine(Ifile);
             nw = line.split(&word);
             if (nw != ndf) goto stop0;
	     for (i=0; i<ndf; i++) if (!word[i].toDbl(&(dTmp[i]),false)) goto stop0;
             for (i=0; i<ndf; i++) { un[i] = dTmp[i]; u[i] = dTmp[i]; word[i].free(); }
             delete [] word; dTmp.free();

             // initial velocities (optional)

             line.getNextLine(Ifile);
             nw = line.split(&word);
             if (nw != ndf) goto stop0;
	     for (i=0; i<ndf; i++) if (!word[i].toDbl(&(dTmp[i])),false) goto stop0;
             for (i=0; i<ndf; i++) { u4[i] = dTmp[i]; u3[i] = dTmp[i]; word[i].free(); }
             delete [] word; dTmp.free();

             line.getNextLine(Ifile);

             break;

             stop0:

             for (i=0; i<nw; i++) word[i].free(); delete [] word; 

             break;

    case  1: cout << "     RIGIDBODY: reading prescribed displacements ...\n\n";

             if (x.dim()<2) prgError(3,fct,"'properties' must precede 'prescribed displacements'!");

             prgError(4,fct,"'prescribed displacements' not yet available!");

             break;

    case  2: cout << "     RIGIDBODY: reading forces ...\n\n";

             if (x.dim()<2) prgError(3,fct,"'properties' must precede 'forces'!");

             j = 0;

             while (1)
             {
       	       line.getNextLine(Ifile);
	       nw = line.split(&word);
	       if (nw!=3) break; 
              
               if (!word[0].toInt(&(iTmp[0]),false)) break;
               if (!word[1].toInt(&(iTmp[1]),false)) break;
               if (!word[2].toDbl(&(dTmp[0]),false))     break;
 
               if (iTmp[0] < 1 || iTmp[0] > ndf) prgError(4,fct,"error in 'forces'!");
               frcTmFct[iTmp[0]-1] = iTmp[1];
               frc[iTmp[0]-1] = dTmp[0];       j++;
              
               for (i=0; i<nw; i++) word[i].free(); delete [] word;
             }

             iTmp.free();
             dTmp.free();

             for (i=0; i<nw; i++) word[i].free(); delete [] word; 

             if (j<1) prgError(1,fct,"input error in 'forces'!");

             break;

    case  3: cout << "     RIGIDBODY: reading controls ...\n\n";

             if (tis > -1)         prgError(1,fct,"'control' has already been read!");
	     
	     line.getNextLine(Ifile);
	     nw = line.split(&word);
	     if (nw < 3)               prgError(1,fct,"input error in 'control'!");
		     
	     if (!word[0].toDbl(&tol)) prgError(2,fct,"input error in 'control'!");

	     if (!word[1].toDbl(&charDiameter)) prgError(3,fct,"input error in 'control'!");

	     if (!word[2].toInt(&tis)) prgError(4,fct,"input error in 'control'!");
	     if (tis < 0) tis = 0;

             for (i=0; i<10; i++) td[i] = 0.;
	    
             nn = nw; if (nn > 13) nn = 13;
	     
             for (i=3; i<nn; i++) 
	       if (!word[i].toDbl(&(td[i-1]))) prgError(5,fct,"input error in 'control'!");

             for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     
       	     line.getNextLine(Ifile);

	     break;

    case  4: cout << "     RIGIDBODY: reading data output ...\n\n";

             j = 0;

             while (1)
             { 
       	       line.getNextLine(Ifile);
	       nw = line.split(&word);
	       if (nw!=3) break; 
              
               if (!word[0].toInt(&(iTmp[j+j+0]),false)) break;
               if (!word[2].toDbl(&(dTmp[j]),false)) break;
               iTmp[j+j+1] = word[1].which(wrndTypeNames);
               if (iTmp[j+j+1] < 0) break;
 
               wrndDoF[j]  = iTmp[j+j+0]; 
               if (wrndDoF[j] > ndf) prgError(2,fct,"input error in 'data output'");
               wrndType[j] = iTmp[j+j+1];
               wrndFact[j] = dTmp[j];     j++;
               
               for (i=0; i<nw; i++) word[i].free(); delete [] word;
             }

             iTmp.free();
             dTmp.free();

             for (i=0; i<nw; i++) word[i].free(); delete [] word; 

             if (j<1) prgError(1,fct,"input error in 'data output'!");

             break;

    case -1: // go and inherit from DOMAIN
	     
	     this->Domain::readInputData(Ifile, line); 
	     
	     break;
  }
 
  return;
}








void RigidBody::prepareInputData(void)
{
  // call ancestor function

  Domain::prepareInputData();

  
  cout << "     RIGIDBODY: prepare input data ...\n\n";
  
  char fct[] = "RigidBody::prepareInputData"; 
	  

  
  return;
}








void RigidBody::prepareInteractions(void)
{
  // call ancestor function

  Domain::prepareInteractions();

  
  cout << "     RIGIDBODY: preparing interactions ...\n\n";
  
  char fct[] = "RigidBody::prepareInteractions"; 
	  
  int i, j;
  

  // replace time function ids in frcTmFct

  for (i=0; i<frcTmFct.n; i++)
  {
    if (frcTmFct[i] != 0)
    {
      j = 0; while (j<timeFunction.n && frcTmFct[i]!=timeFunction[j].id) j++;

      if (j == timeFunction.n) prgError(1,"RigidBody::prepareInteractions",
                                          "invalid time function id in 'forces'!"); 
      frcTmFct[i] = j;
    }
    else
    {
      frcTmFct[i] = -1;
    }
  }
  
  return;
}







void RigidBody::setTimeParam(void)
{
  double dt = mpapTime.dt, alpf, alpm, beta, gamm, rho;
	
  for (int i=10; i<TD_DIM; i++) td[i] = 0.;

  switch (tis)
  {
    case  1: // generalised midpoint rule

             gamm = td[0];

	     td[6-1] = gamm;
	     td[7-1] = gamm;
	     td[8-1] = gamm;
	     
             td[10] = - (1. - gamm) / gamm;         
             td[11] = 1. / (gamm * dt);

             td[9-1] = td[11] * td[11];  // coefficient of u_n+1 in ddu_n+1 
	     
	     break;

    case  2: // generalised alpha-method

             rho  = td[0];
	     
             alpf = 1. / (1. + rho);
	     alpm = (2. - rho) / (1. + rho);
	     beta = .25 * (1. + alpm - alpf)*(1. + alpm - alpf);
             gamm = .5 + alpm - alpf;
	     
	     td[6-1] = alpf;
	     td[7-1] = alpf;
	     td[8-1] = alpm;
	     
             td[10] = gamm / (beta * dt);             // U_n+1   in  dU_n+1
             td[11] = - td[10];                       // U_n
             td[12] = 1. - gamm / beta;               // dU_n    
             td[13] = dt * (1. - gamm / (2. * beta)); // ddU_n
	     
             td[14] = 1. / (beta * dt * dt);          // U_n+1   in  ddU_n+1
             td[15] = - td[14];                       // U_n
             td[16] = - 1. / (beta * dt);             // dU_n
             td[17] = 1. - 1. / (2. * beta);          // ddU_n

             td[9-1]  = td[14];  // coefficient of U_n+1 in ddU_n+1 
	     
	     break;
	     
    default: prgError(1,"RigidBody::setTimeParam","invalid value of tis!");
  }

  return;
}







void RigidBody::timeUpdate(void)
{
  int i;

  double *X = x.x, *Xn = xn.x, *X0 = x0.x,
	 *U = u.x, *Un = un.x, *dU = u3.x, 
	 *dUn = u4.x, *ddU = u5.x, *ddUn = u6.x;

  // xn <- x
 
  for (i=0; i<ndf; i++) Xn[i] =  X[i];

  // un <- u
  
  for (i=0; i<ndf; i++)   Un[i] = U[i];
  
  // velocities, accelerations

  switch (tis)
  {
    case  1: // generalised midpoint rule

    case  2: // generalised alpha-method

	     for (i=0; i<ndf; i++)
	     {
	        dUn[i] =  dU[i];
	       ddUn[i] = ddU[i];
	     }
	     
	     break;

    default: prgError(1,"RigidBody::timeUpdate","invalid value of tis!");
  }
 
  // prescribed displacements  (set increments only)

  //updateUDepIncrements();
  
  // set iteration flag
  
  firstIter = true;
 
  localStiffnessError = 0;

  // update X, dU, ddU    !important!
  
  updateIterStep();

  return;   
}






void RigidBody::updateIterStep(void)
{
  int i, j;

  double *X = x.x, *X0 = x0.x, *U = u.x, 
	 *Un = un.x, *dU = u3.x, *dUn = u4.x, *ddU = u5.x, *ddUn = u6.x;

  // update coordinates
  
  for (i=0; i<ndf; i++) X[i] = X0[i] + U[i];

  switch (tis)
  {
    case  1: // generalised midpoint rule

	     for (i=0; i<ndf; i++)
	     {
                dU[i] = td[10]* dUn[i] + td[11]*( U[i]- Un[i]);
               ddU[i] = td[10]*ddUn[i] + td[11]*(dU[i]-dUn[i]);
	     }
	     
	     break;

    case  2: // generalised alpha-method

	     for (i=0; i<ndf; i++)
	     {
                dU[i] = td[10]*(U[i]-Un[i]) + td[12]*dUn[i] + td[13]*ddUn[i];
               ddU[i] = td[14]*(U[i]-Un[i]) + td[16]*dUn[i] + td[17]*ddUn[i];
	     }
	     
	     break;
	     
    default: prgError(1,"RigidBody::updateIterStep","invalid value of tis!");
  }

  return;
}








int RigidBody::calcStiffnessAndResidual(int printRes, bool zeroMtx, bool zeroRes)
{
  int i;

  double ul, dul, ddul, force;

  if (firstIter) rNorm = -1.;

  rNormPrev = rNorm;
  rNorm     = 0.;

  localStiffnessError = 0;

  for (i=0; i<ndf; i++)
  {  
    ul   = td[5] *  u[i] + (1.-td[5]) * un[i];
    dul  = td[6] * u3[i] + (1.-td[6]) * u4[i];
    ddul = td[7] * u5[i] + (1.-td[7]) * u6[i];

    force = frc[i];

    if (frcTmFct[i] > -1) force *= timeFunction[frcTmFct[i]].prop;

    reac[i] =  prop[i+ndf*0] * ddul
             + prop[i+ndf*1] * dul
             + prop[i+ndf*2] * ul
             + force;

    if (idu[i] > 0) rNorm += reac[i]*reac[i];
  }
  rNorm = sqrt(rNorm);
  firstIter = false;

  if (printRes > 1) { COUT << domain.name(this); printf("  %11.4e\n",rNorm); }

  return 0;
}






int RigidBody::factoriseSolveAndUpdate(void)
{
  int i;

  double stiff;

  for (i=0; i<ndf; i++)
  {
    if (idu[i] > 0)
    {
      stiff = prop[i+ndf*0] * td[14] * td[7]
            + prop[i+ndf*1] * td[10] * td[6]
            + prop[i+ndf*2]          * td[5];

      u[i]  -= reac[i] / stiff;
    }
  }

  return 0;
}







void RigidBody::writeNodalData(void)
{
  int         i, idf;
  double      val, fact;
  char        tmp[20];
  MyString    tmpStr, fname;
  
  char fct[] = "RigidBody::writeNodalData";
 
  if (mpapTime.write + 1.e-12 < mpapTime.cur) 
  { 
    sprintf(tmp,"\n%12.5g",mpapTime.cur);
    tmpStr.append(tmp);
    mpapTime.write = mpapTime.cur;
  }

  for (i=0; i<wrndType.n; i++)
  {
    fact = wrndFact[i];
    idf  = wrndDoF[i] - 1;
    
    switch (wrndType[i])
    {
      case  0: // u 
	     
               val = u[idf]; break;
	       
      case  1: // du 
	      
	       val = u3[idf]; break;
	       
      case  2: // ddu 
	      
	       val = u5[idf]; break;
	       
      case  3: // reac
	      
	       val = reac[idf]; break;
	       
      case  4: // x
	      
	       val = x[idf]; break;
	       
      default: prgError(1,fct,"unknown wrndType!");	     
    }
    val *= fact;
    sprintf(tmp," %12.5g",val);
    tmpStr.append(tmp);
  }
  
  if (!files.Tout.is_open()) 
  {
    fname.append(files.projDir).append(SLASH).append(files.Tfile);
    files.Tout.open(fname.asCharArray(),ios_base::app);
  }

  if (tmpStr[0] != '\n') files.Tout << ' ';
  files.Tout << tmpStr.stripToMin();

  files.Tout.close();

  return;
}

