
#include "Element.h"
#include "Debug.h"
#include "ElementGroup.h"
#include "Mesh.h"
#include "FunctionsProgram.h"
#include "PropertyTypeEnum.h"
#include "FiniteElementBVP.h"


using namespace std;


static int cnt;        // this should be removed together with 
                       // the stuff below

Element::Element(void)
{
  GpDatId.setDim(1);

  GpDatId[0]   = cnt++;  // this is dummy stuff and should be removed
  //GpDatId[1] = cnt++;  // once some macro for initialising multiscale
  //GpDatId[2] = cnt++;  // problems has been implemented
  //GpDatId[3] = cnt++;  //

	
  ix = NULL;

  intVar1 = NULL;
  intVar2 = NULL;

  sparse[DOF] = NULL;
  sparse[MSH] = NULL;
  
  if (debug) cout << " constructor Element\n\n";

  return;
}





Element::~Element()
{
  if (ix != NULL) delete [] ix; ix = NULL;

  if (intVar1!=NULL) delete [] intVar1;
  if (intVar2!=NULL) delete [] intVar2;
 
  if (sparse[DOF] != NULL) delete [] sparse[DOF];
  if (sparse[MSH] != NULL) delete [] sparse[MSH];
 
  if (debug) cout << " destructor Element\n\n";

  return;
}








int &Element::idu(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndf()) prgError(1,"Element::idu","invalid parameters i,j!");

  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->idu(ix[i-1],j);
}








int &Element::idx(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndm()) prgError(1,"Element::idx","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->idx(ix[i-1],j);
}








double &Element::x0(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndm()) prgError(1,"Element::x0","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->x0(ix[i-1],j);
}








double &Element::x(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndm()) prgError(1,"Element::x","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->x(ix[i-1],j);
}








double &Element::xn(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndm()) prgError(1,"Element::xn","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->xn(ix[i-1],j);
}








double &Element::u(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndf()) prgError(1,"Element::u","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->u(ix[i-1],j);
}





double &Element::un(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndf()) prgError(1,"Element::un","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->un(ix[i-1],j);
}





double &Element::u3(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndf()) prgError(1,"Element::u3","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->u3(ix[i-1],j);
}





double &Element::u4(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndf()) prgError(1,"Element::u4","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->u4(ix[i-1],j);
}





double &Element::u5(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndf()) prgError(1,"Element::u5","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->u5(ix[i-1],j);
}





double &Element::u6(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndf()) prgError(1,"Element::u6","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->u6(ix[i-1],j);
}





double &Element::d(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndm()) prgError(1,"Element::d","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->d(ix[i-1],j);
}





double &Element::d0(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndm()) prgError(1,"Element::d0","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->d0(ix[i-1],j);
}





double &Element::v(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndm()) prgError(1,"Element::v","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->v(ix[i-1],j);
}





double &Element::vn(int i, int j)
{
  if (i<1 || j<1 || i>nen() || j>ndm()) prgError(1,"Element::vn","invalid parameters i,j!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->vn(ix[i-1],j);
}





double &Element::outp(int i, int j)   // j will be ignored!
{
  if (i<1 || i>nen()) prgError(1,"Element::outp","invalid parameter i!");
 
  return ((Mesh*)(((ElementGroup*)elemGrp)->dom))->outp[ix[i-1]-1];
}








void Element::belongsTo(void *eG, void *msh)
{
  // associated element with element group and vice versa
	
  ElementGroup *elGr = (ElementGroup*) eG;
	
  elemGrp = eG;

  elGr->elem.add(this);

  
  // check basic consistency of element with element group
  
  if (elGr->ndm == 0) { elGr->ndm = ndm(); elGr->closed = isClosed(); }
  else if (elGr->ndm != ndm())    prgError(1,"Element::belongsTo","ndm inconsistent!"); 

  if (elGr->nen == 0) elGr->nen = nen();
  else if (elGr->nen != nen())    prgError(1,"Element::belongsTo","nen inconsistent!"); 

  if (elGr->ndf == 0) elGr->ndf = ndf();
  else if (elGr->ndf != ndf())    prgError(1,"Element::belongsTo","ndf inconsistent!"); 

  if (elGr->closed != isClosed()) prgError(1,"Element::belongsTo","closed inconsistent!");
  
  
  // check basic consistency of element with domain

  Mesh *m = (Mesh*) msh;
  
  //cout << m->ndf << "," << m->ndm << "," << m->nen << "\n";
  
  if (m->ndf < ndf()) prgError(2,"Element::belongsTo","ndf inconsistent!");
  
  if (m->ndm < ndm()) prgError(2,"Element::belongsTo","ndm inconsistent!");
  
  if (m->nen < nen()) prgError(2,"Element::belongsTo","nen inconsistent!");
  
  
  return;
}














void Element::print(void)
{
  int i; 

  int nelm, nmsh, nmat, elmType, mshType, matType;
	
  ElementGroup *eG = (ElementGroup*) elemGrp;

  PropertyItem **eP = eG->elemProp;

  for (i=0; i<eG->nElemProp; i++) 
  { 
    cout << " " << eP[i]->id <<" = " << eP[i]->name << "   " << eP[i]->data << "\n";
  }  
 
  for (i=1; i<nen()+1; i++) cout << " --> " << x(i,1) << "," << x(i,2);

  cout << "\n\n";
  
  //cout << eP[ELEMENTTYPE]->name << "\n\n";

  
  return;
}









void Element::diffStiffTest(double ddd, int dig, int dig2, bool gfrmt)
{
  int    i, jnd, jj, j, k, nst = nen()*ndf(), niv = nGaussPoints() * nivGP();

  ElementGroup *eG = (ElementGroup*) elemGrp;

  double *s     = eG->dom->s,
         *p     = eG->dom->p,
         dd[6]  = {-3.*ddd, -2.*ddd, -ddd, +ddd, +2.*ddd, +3.*ddd },
	 *r     = new double[6*nst],
	 *sdiff = new double[nst*nst],
	 smax, sdmax;

  COUT << "element type: " << ((ElementGroup*)elemGrp)->elemProp[ELEMENTTYPE]->name << "\n\n";
 
  for (i=0; i<nivGP()*nGaussPoints(); i++) intVar1[i] = intVar2[i];

  // loop over columns of s

  for (jnd=0; jnd<nen(); jnd++) // nodes
  {
    for (jj=0; jj<ndf(); jj++)  // dof
    {
      j = jnd*ndf() + jj;

      // loop over perturbations
	    
      for (k=0; k<6; k++)
      {
        // apply pertubation
	      
        u(jnd+1,jj+1) += dd[k];

        eG->dom->updateIterStep();
	
        // calculate residual

        calcStiffnessAndResidual();

       if(jnd == 0 && jj == 0)
       {
          cout.setf(ios::fixed);
          cout.setf(ios::showpoint);
          cout.precision(8);
          for(int ii=0;ii<nst;ii++)
          {
              cout  << '\t' << p[ii] << endl ;
          }
               cout << endl;
               cout << endl;
       }
	
        for (i=0; i<niv; i++) intVar2[i] = intVar1[i];
  
        // remove pertubation
	
        u(jnd+1,jj+1) -= dd[k];
	
        eG->dom->updateIterStep();
	
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

  calcStiffnessAndResidual();

  for (i=0; i<niv; i++) intVar2[i] = intVar1[i];
  
  // compare the matrices

  prgCompareTwoSimpleMatrices(sdiff,                         // matrix 1
		              s,                             // matrix 2
		              "numerical differentiation",   // title matrix 1 
			      "analytical calculation",      // title matrix 2
			      "numerical - analytical",      // title matrix 1 - 2
			      nst, nst,                      // matrix dimension
			      dig,dig2,gfrmt,                // format
			      0,                             // indentation
			      false,                         // interactive
			      false);                        // row/column numbers
  delete [] r;
  delete [] sdiff;

  return;
}










double Element::getElemSizeOpt(double *xp, double *N)
{
  // make sure N points to an array with at least nen doubles

  double hh = 0., *h = ((Mesh*)(((ElementGroup*)elemGrp)->dom))->elemSizeOpt.x;

  if (containsPoint(xp,N))
  {
    for (int i=0; i<nen(); i++) hh += N[i]*h[ix[i]-1]; 

    return hh;
  }
  prgError(1,"Element::getElemSizeOpt","point not inside element!");

  return 0.;
}








void Element::assembleReactionForces(void)
{ 
  Mesh *mesh = (Mesh*)(((ElementGroup*)elemGrp)->dom); 
 
  int    indfDom, indf = 0, i, j, ndfDom = mesh->ndf;

  double *REAC = mesh->reac.x, *p = mesh->p;

  for (i=0; i<nen(); i++)  
  { 
    indfDom = (ix[i]-1) * ndfDom;

    for (j=0; j<ndf(); j++)  REAC[indfDom+j] += p[indf++];
  } 

  return;
}










void Element::adjustAssembly(int p0, int fl)
{
  if (sparse[0] == NULL) return;

  int i, j, nstu = nen() * ndf(), nstx = nen() * ndm(), nst2;

  if (fl == DOF) nst2 = nstu * nstu; 
  else           nst2 = nstx * nstx;

  for (i=0; i<nst2; i++)
  {
    //cout << sparse[fl][i];
    if (sparse[fl][i] > -1) sparse[fl][i] += p0;
    //cout << " -> " << sparse[fl][i] << "\n";
  }

  return;
}
 
