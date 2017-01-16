
#include <iostream>
#include <fstream>

#include "MyString.h"
#include "Definitions.h"
#include "FunctionsProgram.h"
#include "FunctionsMaterial.h"
#include "FunctionsSupport.h"


using namespace std;






bool adjustF(double *F, double *Forig, double *C, double dd, int j, bool BFlag)
{
  double H[9], R[6], U[6], K[36], upd[6];

  int P[6], n = 6, isw = 0, iter = 0;

  unit_s_p_(U);

  while (1)
  {
    if (iter++ > 6) return false;

    set_s_m_(R,C);
    R[j] -= dd;
    mult_u_ss_p_(H,C,U);
    mult_s_su_ap_(R,U,H);

    //cout << norm2_s_(R) << "\n";
    if (norm2_s_(R) < 1.e-12) 
    { 
      //cout << "\n";
      if (BFlag) mult_u_su_p_(F,U,Forig);
      else       mult_u_us_p_(F,Forig,U);
      return true; 
    }

    diff_tat_p_(K,H);

    decomplr_matrix_(K,P,&n,&isw);
    solve_matrix_   (K,P,R,upd,&n);

    upd[3] *= .5;
    upd[4] *= .5;
    upd[5] *= .5;
     
    set_s_am_(U,upd);
  }
}







void matDiffTangentTest(double ddd, char *fname, int dig, int dig2, bool gfrmt)
{
  char fct[] = "matDiffTangentTest", 
       *strainType[] = {"small", "finite"},
       *ssState[]    = {"plane stress", "plane strain", "axisymmetric"},
       *matType[]    = MATERIAL_TYPE_NAMES;

  int i, j, nw, tnsr, sss, finite, errCode, matId, ndm, nivGP, isw, matDim[] = MATERIAL_DIMENSIONS;

  double *F = NULL, *iv1 = NULL, *iv2 = NULL, *matDat = NULL, tmp, dt;
  
// read data file 

  ifstream file;	

  MyString line, *word = NULL;
  
  file.open(fname);

  if (!file) { prgWarning(1,fct,"failed to open file!"); return; }
  
  // read material type and parameters

  if (!line.getNextLine(file)) { prgWarning(2,fct,"failed to read material type!"); goto error; }
 
  nw = line.split(&word);

  matId = word[0].which(matType) + 1;

  if (matId < 1) { prgWarning(3,fct,"invalid material type!"); goto error; }

  if (matId == 7 || matId == 8)
   { prgWarning(4,fct,"'diff,mat' is not applicable to multiscale materials!"); goto error; }

  if (matId == 9) { prgWarning(5,fct,"dummy material ?!"); goto error; }

  matDat = new double [nw-1];

  for (j=1; j<nw; j++)
  {
    if (!word[j].toDbl(&tmp)) { prgWarning(6,fct,"invalid material parameter!"); goto error; }
    matDat[j-1] = tmp;
  }

  for (j=0; j<nw; j++) word[j].free(); delete [] word; word = NULL;

  ndm = matDim[matId-1];

  if (ndm != 2 && ndm != 3) { prgWarning(7,fct,"2D or 3D only!"); return; }

  // read time step size dt, finite and stress/strain state flags

  if (!line.getNextLine(file)) 
  { prgWarning(8,fct,"failed to read time step dt, finite and stress/strain flags!"); goto error; }

  nw = line.split(&word);
  
  if ((nw < 3 && ndm == 2) || nw < 2) 
    { prgWarning(10,fct,"failed to read dt, finite flag and stress/strain state!"); goto error; }

  if (!word[0].toDbl(&dt)) { prgWarning(11,fct,"failed to read time step size dt!"); goto error; }

  if (!word[1].toInt(&finite)) { prgWarning(11,fct,"failed to read finite flag!"); goto error; }

  if (ndm == 2 && !word[2].toInt(&sss))
    { prgWarning(12,fct,"failed to read stress/strain state flag!"); goto error; }
  
  for (j=0; j<nw; j++) word[j].free(); delete [] word; word = NULL;
    
  // read deformation gradient

  F = new double [ndm*ndm];
  
  for (i=0; i<ndm; i++)
  {
    if (!line.getNextLine(file)) 
      { prgWarning(13,fct,"failed to read deformation gradient!"); goto error; }
    
    nw = line.split(&word);

    if (nw != ndm) 
      { prgWarning(14,fct,"failed to read deformation gradient!"); goto error; }
    
    for (j=0; j<ndm; j++) 
      if (!word[j].toDbl(&(F[j*ndm+i])))
        { prgWarning(15,fct,"failed to read deformation gradient!"); goto error; }

    for (j=0; j<nw; j++) word[j].free(); delete [] word; word = NULL;
  }

  // read internal variables

  isw = 1; // get nivGP

  if (ndm == 2) matlib2d_(matDat,F,F,F,F,F,F,F,&matId,&nivGP,&finite,&sss,&isw,&errCode);

  else          matlib3d_(matDat,F,F,F,F,F,F,&matId,&nivGP,&finite,     &isw,&errCode);

  if (nivGP > 0)
  {
    iv1 = new double [nivGP];
    iv2 = new double [nivGP];

    if (!line.getNextLine(file)) 
      { prgWarning(16,fct,"failed to read internal variables!"); goto error; }
    
    nw = line.split(&word);

    if (nw != nivGP)
      { prgWarning(17,fct,"failed to read internal variables!"); goto error; }
	
    for (j=0; j<nivGP; j++) 
      if (!word[j].toDbl(&(iv1[j])))
        { prgWarning(18,fct,"failed to read internal variables!"); goto error; }

    for (j=0; j<nw; j++) word[j].free(); delete [] word; word = NULL;
  }
  
  file.close();

  if (finite < 0) finite = 0;
  if (finite > 1) finite = 1;

  if (sss < 1) sss = 1;
  if (sss > 3) sss = 3;
 
  // print data

  cout << "\n";
  COUT << "material diff-test\n";
  COUT << "---------------------------------------------------------\n";
  COUT << "material type: " << ndm << "D, " << matType[matId-1] << "\n\n";
  COUT << strainType[finite] << " strains\n";
  if (ndm == 2) COUT << ssState[sss-1] << "\n"; cout << "\n";
  COUT << "time step size dt = " << dt << "\n\n";
  COUT << "deformation gradient F\n\n          "; 
  for (i=0; i<ndm; i++) 
    { for (j=0; j<ndm; j++) printf(" %7.3f",F[j*ndm+i]); cout << "\n          "; }
  cout << "\n";
  if (nivGP > 0)
  {
    COUT << "internal variables\n\n          ";
    for (j=0; j<nivGP; j++) printf(" %7.3f",iv1[j]); 
    cout << "\n\n";
  }

  goto diff;
  
  error:

    if (file) file.close();

    if (word != NULL) { for (i=0; i<nw; i++) word[i].free(); delete [] word; }	

    if (matDat != NULL) delete [] matDat;
    if (F      != NULL) delete [] F;
    if (iv1    != NULL) delete [] iv1;
    if (iv2    != NULL) delete [] iv2;

    return;
 
  diff: 

// do the actual job now

  int scomp = ndm+ndm, k, i1, i2, n,
      ic2D[] = {1},
      ic3D[] = {1,5,2};
	
  double dd[6]  = {-3.*ddd, -2.*ddd, -ddd, +ddd, +2.*ddd, +3.*ddd },
	 *r     = new double[6*scomp],
	 *cdiff = new double[scomp*scomp],
	 *c     = new double[scomp*scomp],
	 *h4    = new double[scomp*scomp],
	 *stre  = new double[scomp], 
	 *stre2 = new double[scomp], 
	 *Forig = new double[ndm*ndm],
         *B     = new double[ndm+ndm],
         *C     = new double[ndm+ndm],
	 cmax, cdmax, F33, detF, r2 = 2.e0;

  isw = 3;

  if (ndm == 3)
  {
    calcb_(B,F);
    calcc_(C,F);
  }

  for (i=0; i<ndm*ndm; i++) Forig[i] = F[i];
  for (i=0; i<nivGP; i++)   iv2[i] = iv1[i];
    
  n = scomp; if (ndm == 2 && sss < 3) n--;
    
  // loop over components of strain

  for (j=0; j<n; j++)
  {
    // loop over perturbations
	    
    for (k=0; k<6; k++)
    {
      // apply pertubation

      switch (finite)      
      {
        case  0: if (j < ndm) F[j*ndm+j] += dd[k];   // sig-eps
			 
                 else if (j<(ndm-1)*3) 
                 { 
                   if (ndm == 2) i = ic2D[j-ndm]; else i = ic3D[j-ndm];

                   F[i] += (dd[k] + dd[k]);
                 } 
                 else F33 += dd[k];

		 break;

        case  1: if (ndm == 3)
                 { 
                   if (!adjustF(F,Forig,B,dd[k],j,true))
                     { prgWarning(1,fct,"no convergence in adjustF!"); return; }
                 }
                 else { prgWarning(1,fct,"no 2D finite strains yet!"); return; }

                 detF = det_u_(F);

                 break;
	 
        default: prgWarning(1,fct,"what shall we do with the finite sailor?"); return;
      }
	    
      // calculate stress

      if (ndm == 2)
      
	matlib2d_(matDat,F,&F33,stre,c,iv1,iv2,&dt,
	          &matId,&nivGP,&finite,&sss,&isw,&errCode);

      else

	matlib3d_(matDat,F,stre,c,iv1,iv2,&dt,
	          &matId,&nivGP,&finite,&isw,&errCode);
     
      for (i=0; i<nivGP; i++) iv2[i] = iv1[i];
      
      // remove pertubation

      for (i=0; i<ndm*ndm; i++) F[i] = Forig[i];
     
      // loop over rows of s, store stress

      if (finite == 1) for (i=0; i<scomp; i++) r[i*6+k] = stre[i] * detF;

      else             for (i=0; i<scomp; i++) r[i*6+k] = stre[i];
    }

    // loop over rows of s

    for (i=0; i<n; i++)
    {		
      cdiff[j*scomp+i] = ( -        r[i*6+0]
                            +  9. * r[i*6+1]
                            - 45. * r[i*6+2]
                            + 45. * r[i*6+3]
                            -  9. * r[i*6+4]
                            +       r[i*6+5] ) / (60. * ddd);

      if (j >= ndm) cdiff[j*scomp+i] *= .5;
    }
  }

  // calculate tangent

  for (i=0; i<scomp*scomp; i++) c[i] = 0.;
 
  if (ndm == 2) 
          
    matlib2d_(matDat,F,&F33,stre,c,iv1,iv2,&dt,
              &matId,&nivGP,&finite,&sss,&isw,&errCode);

  else

    matlib3d_(matDat,F,stre,c,iv1,iv2,&dt,
              &matId,&nivGP,&finite,&isw,&errCode);
  
  for (i=0; i<nivGP; i++) iv2[i] = iv1[i];

  // modify numerical tangent (required for some finite strain stuff)

  switch (finite)
  {
    case 1: if (ndm == 3)  
            { 
              detF = det_u_(Forig); 
              set_s_f_(stre,stre,&detF); 
              i = 6;  pzero_(stre2,&i);
              i = 36; pzero_(h4,&i);
              taub2sigccpart1_(stre2,h4,stre,cdiff,B,&detF);
              taub2sigccpart2_(stre2,h4); 
              set_4_p_(cdiff,h4); 
            }  
            break;
  }
 
  // compare the matrices

  prgCompareTwoSimpleMatrices(cdiff,                         // matrix 1
		              c,                             // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      scomp,scomp,                   // matrix dimension
			      dig,dig2,gfrmt,                // format
			      10,                            // indentation
			      false,                         // interactive
			      false);                        // row/column numbers
  delete [] F;
  delete [] iv1;
  delete [] iv2;
  delete [] matDat; 
  delete [] r;
  delete [] Forig;
  delete [] B;
  delete [] C;
  delete [] c;
  delete [] cdiff;
  delete [] h4;
  delete [] stre;
  delete [] stre2;
  
  return;
}






