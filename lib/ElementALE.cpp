

#include "ElementALE.h"
#include "PropertyTypeEnum.h"
#include "ElementGroup.h"
#include "FunctionsProgram.h"


void ElementALE::diffAleTest(double ddd, int dig, int dig2, bool gfrmt)
{
  int    i, jnd, jj, j, k, nst = nen()*ndm();
	
  ElementGroup *eG = (ElementGroup*) elemGrp;
	
  double *s     = eG->dom->s,
         *p     = eG->dom->p,
         dd[6]  = {-3.*ddd, -2.*ddd, -ddd, +ddd, +2.*ddd, +3.*ddd },
	 *r     = new double[6*nst],
	 *sdiff = new double[nst*nst],
	 smax, sdmax;

  COUT << "element type: " << ((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->name << "\n";
  COUT << "    ALE type: " << ((ElementGroup*)elemGrp)->elemProp[ALETYPE]->name << "\n\n";
  
  // loop over columns of s

  for (jnd=0; jnd<nen(); jnd++) // nodes
  {
    for (jj=0; jj<ndm(); jj++)  // dof
    {
      j = jnd*ndm() + jj;
	    
      // loop over perturbations
	    
      for (k=0; k<6; k++)
      {
        // apply pertubation
	      
        d(jnd+1,jj+1) += dd[k];

        //eG->dom->updateIterStep();
	
        // calculate residual

        calcStiffnessAndResidualMesh();
	
        // remove pertubation
	
        d(jnd+1,jj+1) -= dd[k];
	
        //eG->dom->updateIterStep();
	
        // loop over rows of s, store residual

        for (i=0; i<nst; i++) r[i*6+k] = p[i];
      }

      // loop over rows of s

      for (i=0; i<nst; i++)
		
        sdiff[j*nst+i] = ( +       r[i*6+0]
                           -  9. * r[i*6+1]
                           + 45. * r[i*6+2]
                           - 45. * r[i*6+3]
                           +  9. * r[i*6+4]
                           -       r[i*6+5] ) / (60. * ddd);
    }
  }
 
  // calculate stiffness

  calcStiffnessAndResidualMesh();

  // compare the matrices

  prgCompareTwoSimpleMatrices(sdiff,                         // matrix 1
		              s,                             // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      nst,nst,                       // matrix dimension
			      dig,dig2,gfrmt,                // format
			      0,                             // indentation
			      false,                         // interactive
			      false);                        // row/column numbers
  delete [] r;
  delete [] sdiff;

  return;
}





