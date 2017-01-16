
#include <iostream>

#include "FunctionsProgram.h"
#include "DomainTypeEnum.h"
#include "Mesh.h"
#include "DomainTree.h"
#include "DomainType.h"
#include "Debug.h"
#include "Plot.h"
#include "MathVector.h"
#include "MathBasic.h"
#include "MathGeom.h"
#include "Element.h"
#include "DataBlockTemplate.h"
#include "Definitions.h"
#include "PropertyTypeEnum.h"
#include "List.h"
#include "TimeFunction.h"
#include "Fluid.h"
#include "Solid.h"
#include "FreeSurface.h"
#include "MicroCell.h"
#include "MicroCellWulf.h"
#include "MpapTime.h"
#include "ComputerTime.h"
#include "SolverMA41.h"
// all elements (required for Mesh::newElement)
#include "Element2D3nodedTriangle.h"
#include "Element2D6nodedTriangle.h"
#include "Element2D4nodedQuadrilateral.h"
#include "Element2D8nodedQuadrilateral.h"
#include "Element2D9nodedQuadrilateral.h"
#include "Element3D4nodedTetrahedron.h"
#include "Element3D10nodedTetrahedron.h"
#include "Element3D8nodedBrick.h"
#include "Element3D20nodedBrick.h"
#include "Element2D3nodedStabIncompFluid.h"
#include "Element2D3nodedStabIncompHighReFluid.h"
#include "Element2D3nodedStabCompFluid.h"
#include "Element2D3nodedLinearSolid.h"
#include "Element2D6nodedQuadraticSolid.h"
#include "Element2D4nodedLinearSolid.h"
#include "Element2D4nodedFBarSolid.h"
#include "Element2D8nodedQuadraticSolid.h"
#include "Element2D9nodedQuadraticSolid.h"
#include "Element3D4nodedStabIncompFluid.h"
//#include "Element3D4nodedStabIncompHighReFluid.h"
//#include "Element3D4nodedStabCompFluid.h"
#include "Element3D4nodedLinearSolid.h"
#include "Element3D10nodedQuadraticSolid.h"
#include "Element3D8nodedLinearSolid.h"
#include "Element3D8nodedFBarSolid.h"
#include "Element3D20nodedQuadraticSolid.h"
#include "Element2D2nodedLine.h"
#include "Element2D3nodedLine.h"
#include "Element2D2nodedGeomExSmallStrainBeam.h"
#include "Element2D2nodedKirchhoffBeam.h"
#include "Element2D2nodedTruss.h"
#include "Element1D2nodedLine.h"
#include "Element1D2nodedAdvectionDiffusion.h"
#include "Element1D2nodedSchroedingerCN.h"
#include "Element1D2nodedSchroedingerST.h"
#include "Element2D2nodedPressureLoad.h"
#include "Element2D3nodedPressureLoad.h"
#include "Element1D2nodedPipeFlowST.h"
#include "Element1D2nodedFlexibleWing.h"
#include "Element2D3nodedLinearPoisson.h"
#include "Element1D2nodedLinearSolidALE.h"
#include "Element2D2nodedFreeSurface.h"


extern Plot               plot;
extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern MpapTime           mpapTime;
extern ComputerTime       computerTime;


static const unsigned int pow2[8] = {1,2,4,8,16,32,64,128};


using namespace std;



Mesh::Mesh(void)                       
{                                                  
  if (debug) std::cout << " Mesh constructor\n\n";
 
  nen        = 0;
  numnp      = 0;
  numel      = 0;
  numElemGrp = 0;
  nBnd       = 0;
  nBndNd     = 0;

  lastSearchNode = 1;
  
  elem       = NULL;
  nodeElem   = NULL;
  nodeNode   = NULL;
  solverMesh = NULL;
  bnd        = NULL;
  bndNd      = NULL;
  surf3D     = NULL;

  addElemGrpProp(ELEMENTTYPE);
  
  // add new type
  
  DomainType *mesh = domain.newType(MESH,GEOMETRY);

  if (mesh == NULL) return;  // domain type exists already

  mesh->key.addNew("nodes per element", 
		   "coordinates", 
		   "elements", 
		   "boundary conditions", 
		   "element groups", 
		   "element type",
		   "prescribed displacements",
		   "dependent displacements",
		   "identical displacements",
		   "slip boundaries",
		   "forces",
		   "nodal data output",
		   "ALE type",
		   "prescribed mesh motion");

  ctimUpdateMesh = 0.;

  return;
}




	                                          
Mesh::~Mesh(void)                     
{
  if (elem != NULL) delete [] elem; elem = NULL;
    // the elements themselves are destructed by the elemGrp elem list destructor
  
  if (nodeElem != NULL) delete [] nodeElem;

  if (nodeNode != NULL) delete [] nodeNode;
  
  if (solverMesh != NULL) delete solverMesh;

  if (surf3D != NULL) delete surf3D;

  if (bnd != NULL) delete [] bnd;

  if (bndNd != NULL) delete [] bndNd;
 
  if (debug) cout << " Mesh destructor\n\n";

  return;
}





void Mesh::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString tmpl, *word;
 
  char tmp[30], fct[] = "Mesh::readInputData",

       *wrndTypeNames[] = WRND_TYPE_NAMES;

  int nw, i, ii, j, n, e, nst, tp, indx;
 
  double fact;
  
  Vector<double> dTmp, wrndFactTmp;
  Vector<int>    iTmp, lTmp, wrndTypeTmp, wrndIndxTmp;
  MyStringList   sTmp;
  List< Vector<int> > wrndNodeTmp;

  DataBlockTemplate t1, t2;

  switch (domain[MESH].key.whichBegins(line))
  {
    case  0: cout << "     MESH: reading number of nodes per element ...\n\n";
	   
             line.getNextLine(Ifile);
	     
             nw = line.split(&word);
            
	     if (nw < 1)                prgError(1,fct,"input error in 'nodes per element'!");
		     
	     if (!word[0].toInt(&nen))  prgError(4,fct,"input error in 'nodes per element'!");

             for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     
	     // allocate memory for s and p

             nst  = ndf * nen; if (ndm > ndf) nst = ndm * nen;
                 
             p  = new double [nst]; 
             s  = new double [nst*nst]; 
	      
             xl = new double [nst*5]; 
             ul = new double [nst*6]; 
  
       	     line.getNextLine(Ifile);

	     break;

    case  1: cout << "     MESH: reading coordinates ...\n\n";
	   
             if (numnp > 0) prgError(1,fct,"duplicate definition of 'coordinates'!");

             if (ndf < 1) prgError(1,fct,"'dimensions' must precede 'coordinates'!");
            
	     sprintf(tmp,"123 %di %df",ndm+ndf,ndm+3*ndf);  
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	    
             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
		     prgError(2,fct,"data error in 'coordinates'!");
	     
             numnp = dTmp.dim() / (ndm+3*ndf);
	     
              idx.setDim(numnp,ndm,true);
              idu.setDim(numnp,ndf,true);
               x0.setDim(numnp,ndm,true);
	        x.setDim(numnp,ndm,true);
	       xn.setDim(numnp,ndm,true);
	        u.setDim(numnp,ndf,true);
	       un.setDim(numnp,ndf,true);
	       u3.setDim(numnp,ndf,true);
	       u4.setDim(numnp,ndf,true);
	       u5.setDim(numnp,ndf,true);
	       u6.setDim(numnp,ndf,true);
             reac.setDim(numnp,ndf,true);
         nodeFlag.setDim(numnp);
 
	     if (isALE())
	     {
	        d.setDim(numnp,ndm,true);
	       d0.setDim(numnp,ndm,true);
	        v.setDim(numnp,ndm,true);
	       vn.setDim(numnp,ndm,true);
         reacMesh.setDim(numnp,ndm,true);
               for (i=0; i<numnp; i++)
  	       {
	         for (j=0; j<ndm; j++) 
	         {
		     d(i+1,j+1) = 0.;
		    d0(i+1,j+1) = 0.;
		     v(i+1,j+1) = 0.;
		    vn(i+1,j+1) = 0.;
		 }
	       }
	     }
	     
             i = ndm; if (ndf > ndm) i = ndf; 
             r.setDim(numnp*i); 
	     
             j = numnp; if (numel > j) j = numel; 
	     outp.setDim(i*j); 
	     
             for (i=0; i<numnp; i++)
	     {
	       for (j=0; j<ndm; j++) 
	       {
                 idx(i+1,j+1) = iTmp[i*(ndm+ndf)  +j];
                  x0(i+1,j+1) = dTmp[i*(ndm+3*ndf)+j];
                   x(i+1,j+1) = x0(i+1,j+1);
                  xn(i+1,j+1) = x0(i+1,j+1);
	       }
	       for (j=0; j<ndf; j++) 
	       {
                  idu(i+1,j+1) = iTmp[i*(ndm+ndf)  +ndm      +j];
                    u(i+1,j+1) = dTmp[i*(ndm+3*ndf)+ndm      +j];
		   un(i+1,j+1) = u(i+1,j+1);
                   u3(i+1,j+1) = dTmp[i*(ndm+3*ndf)+ndm+ndf  +j];
		   u4(i+1,j+1) = u3(i+1,j+1);
                   u5(i+1,j+1) = dTmp[i*(ndm+3*ndf)+ndm+2*ndf+j];
		   u6(i+1,j+1) = u5(i+1,j+1);
	       }
	     }
	     
             //cout << x0 << "\n\n" << idx << "\n\n" << idu << "\n\n" << u << "\n\n";
	     
	     break;

    case  2: cout << "     MESH: reading elements ...\n\n";

             if (numel > 0) prgError(1,fct,"duplicate definition of 'elements'!");

             if (ndf < 1) prgError(1,fct,"'dimensions' must precede 'elements'!");
             
             if (nen < 1) prgError(1,fct,"'nodes per element' must precede 'elements'!");
             
	     sprintf(tmp,"123 %di",nen + 1);  
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);

             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
		     prgError(2,fct,"data error in 'elements'!");

             numel = iTmp.dim() / (nen + 1);

	     ixTmp.setDim(numel,nen+1,true);
	    
             for (i=0; i<numel; i++)  
	     {
               ixTmp(i+1,nen+1) = iTmp[i*(nen+1)]; 
               for (j=0; j<nen; j++) 
               {
                 ii = iTmp[(i+1)*(nen+1)-nen+j];
                 if (ii < 0) prgError(3,fct,"data error in 'elements'!");
                 ixTmp(i+1,j+1) = ii; 
               }
	     }

             //cout << ix << "\n\n";

	     break;

    case  3: cout << "     MESH: reading boundary conditions ...\n\n";

             if (ndf < 1) prgError(1,fct,"'dimensions' must precede 'boundary conditions'!");
	     
             if (numnp < 1) prgError(1,fct,"'coordinates' must precede 'boundary conditions'!");
            
	     sprintf(tmp,"%di",1+ndm+ndf);      
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	     
             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
	       prgError(2,fct,"data error in 'boundary conditions'!");
	    
	     nw = iTmp.dim() / (1 + ndm + ndf);

	     for (i=0; i<nw; i++)
	     {
               n = iTmp[i*(ndm+ndf+1)];
	       for (j=1; j<ndm+1; j++) idx(n,j) = iTmp[i*(ndm+ndf+1)+j];
	       for (j=1; j<ndf+1; j++) idu(n,j) = iTmp[i*(ndm+ndf+1)+ndm+j];
	     }
	     
             //cout << idu << "\n\n";
	     
	     break;
	     
    case  4: cout << "     MESH: reading element groups ...\n\n";

             if (numel < 1) prgError(1,fct,"'elements' must precede 'element groups'!");
	     
	     sprintf(tmp,"%di",elemGrpProp.n + 1);

	     if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	     
             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,elemGrpPropTmp,dTmp,sTmp,lTmp))
	       prgError(2,fct,"data error in 'element groups'!");

             numElemGrp = elemGrpPropTmp.dim() / (elemGrpProp.n+1);
	     
             for (e=0; e<numel; e++)  
	     {
	       i = 0; 
	       while (i<numElemGrp && ixTmp(e+1,nen+1)!=elemGrpPropTmp[i*elemGrpProp.n+i]) i++;
	       if (i==numElemGrp) prgError(3,fct,"data error in 'element groups'!");
	       ixTmp(e+1,nen+1) = i;
	     }
	   
             if (ndm == 2)
             {
               if (surface.n > 0)
               {
                 for (e=0; e<surface.n; e++)
                 {
                   i = 0; 
	           while (i<numElemGrp && grpType[e]!=elemGrpPropTmp[i*elemGrpProp.n+i]) i++;
	           if (i==numElemGrp) prgError(4,fct,"data error in 'element groups'!");
	           grpType[e] = i;
                 }
               }
               else
               {
                 if (grpType.n > 0) prgError(1,fct,"surprise, surprise!");
                 grpType.append(-1);
               }
             }

	     break;
	     
    case  5: cout << "     MESH: reading element type ...\n\n";

             elemProp.add(new PropertyItem(ELEMENTTYPE));

	     elemProp[elemProp.n-1].readInputData(Ifile,line,"input error in 'element type'!");

	     break;
	     
    case  6: cout << "     MESH: reading prescribed displacements ...\n\n";

	     sprintf(tmp,"%di %df",1+ndf,ndf);  
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	     
             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,ubcIntTmp,ubcDblTmp,sTmp,lTmp))
		     prgError(1,fct,"data error in 'prescribed displacements'!");
	     
	     break;
	     
    case  7: cout << "     MESH: reading dependent displacements ...\n\n";

             if (numnp < 1) prgError(1,fct,"'coor' must precede 'dependent displacements'!");
	     
             n = uDep.n;  e = 0;
	     
             line.getNextLine(Ifile);
             nw = line.split(&word);
             
             if (nw < 9 || (nw-5) % 4 != 0) 
  	       prgError(1,fct,"data error in 'dependent displacements'!");

             while (e == 0)
	     {
               iTmp.free();
               dTmp.free();

               if (!word[0].toInt(iTmp.append()) || !word[1].toInt(iTmp.append())
	        || !word[2].toDbl(dTmp.append()) || !word[3].toInt(iTmp.append())
		|| word[4] != ":") break;
  	     
               uDep.add(new DependentDoF);

               uDep[n].nd     = iTmp[0]; 
               uDep[n].dof    = iTmp[1]; 
               uDep[n].ucBase = dTmp[0]; 
               uDep[n].tmFct  = iTmp[2]; 
	      
               iTmp.free();
               dTmp.free();
 
               if (uDep[n].nd < 1 || uDep[n].nd > numnp)
                 prgError(1,fct,"invalid node number in 'dependent displacements'!");
	       
               if (uDep[n].dof < 1 || uDep[n].dof > ndf)
                 prgError(1,fct,"invalid degree of freedom in 'dependent displacements'!");
	       
               j = 5;  i = 0;
	       
               while (e == 0)
               {
  	         if (!word[j++].toInt(iTmp.append())
  	          || !word[j++].toInt(iTmp.append())
         	  || !word[j++].toDbl(dTmp.append())
  	          || !word[j++].toDbl(dTmp.append())) 
  		    prgError(2,fct,"data error in 'dependent displacements'!");
  	      
                 if (iTmp[i] < 1 || iTmp[i++] > numnp)
		   prgError(2,fct,"invalid node number in 'dependent displacements'!");
		 
                 if (iTmp[i] < 1 || iTmp[i++] > ndf)
		   prgError(2,fct,"invalid degree of freedom in 'dependent displacements'!");
		 
		 if (j == nw)
  	         {
                   for (j=0; j<nw; j++) word[j].free(); delete [] word;
                   line.getNextLine(Ifile); 
                   nw = line.split(&word);

  		   if ((nw-5) % 4 == 0 && nw > 8) break;

		   if (nw % 4 != 0) e = 1; 
		   
                   j = 0;
	         }
	       }

               j = iTmp.n / 2;

               uDep[n].masterNd.setDim(j);
               uDep[n].masterDoF.setDim(j);
               uDep[n].alpha.setDim(j);
               uDep[n].beta.setDim(j);

               for (i=0; i<j; i++)
               {
                 uDep[n].masterNd[i]  = iTmp[i+i];
                 uDep[n].masterDoF[i] = iTmp[i+i+1];
                 uDep[n].alpha[i]     = dTmp[i+i];
                 uDep[n].beta[i]      = dTmp[i+i+1];
               }
	       n++;
	     }
             for (j=0; j<nw; j++) word[j].free(); delete [] word;
	
             /* for (i=0; i<uDep.n; i++)
	       cout << uDep[i].nd << "," << uDep[i].dof << "  " 
		    << uDep[i].masterNd << " " << uDep[i].masterDoF << " "
		    << uDep[i].alpha << " " << uDep[i].beta << "\n\n"; */
	     break;
	     
    case  8: cout << "     MESH: reading identical displacements ...\n\n";

             if (numnp < 1) prgError(1,fct,"'coor' must precede 'identical displacements'!");
	     
	     n = uDep.n;

	     while (1)
	     {
               line.getNextLine(Ifile);
               nw = line.split(&word);
             
	       if (nw < ndf + 2) break;

               iTmp.free();	      
 
               for (i=0; i<nw; i++)  if (!word[i].toInt(iTmp.append()))  break;
  
               for (i=ndf; i<nw; i++)	       
                 if (iTmp[i] < 1 || iTmp[i] > numnp) 
		   prgError(1,fct,"invalid node number in 'identical displacements'!");
			   
               for (j=0; j<ndf; j++)
	       {
	         if (iTmp[j] == 1)
		 {
                   for (i=ndf+1; i<nw; i++)
	           {
                     uDep.add(new DependentDoF);

                     uDep[n].masterNd.setDim(1);
                     uDep[n].masterDoF.setDim(1);
                     uDep[n].alpha.setDim(1);
                     uDep[n].beta.setDim(1);

		     uDep[n].nd           = iTmp[i];
		     uDep[n].masterNd[0]  = iTmp[ndf];
		     uDep[n].masterDoF[0] = j + 1;
	             uDep[n].dof          = j + 1; 
		     uDep[n].tmFct        = 0;
		     uDep[n].ucBase       = 0.;
		     uDep[n].alpha[0]     = 1.;
		     uDep[n++].beta[0]    = 1.;
	           }
		 }
	       }	       
               for (j=0; j<nw; j++) word[j].free(); delete [] word;
	     }		       
             for (j=0; j<nw; j++) word[j].free(); delete [] word;
		
             /*for (i=0; i<uDep.n; i++)
	       cout << uDep[i].nd << "," << uDep[i].dof << "  " 
		    << uDep[i].masterNd << " " << uDep[i].masterDoF << " "
		    << uDep[i].alpha << " " << uDep[i].beta << "\n\n"; */
	     
	     break;
	     
    case  9: cout << "     MESH: reading slip boundaries ...\n\n";

	     // ......

	     
	     break;
	     
    case 10: cout << "     MESH: reading forces ...\n\n";

             if (numnp < 1) prgError(1,fct,"'coor' must precede 'forces'!");

	     sprintf(tmp,"%di %df",1+ndf,ndf);  
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	     
             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
		     prgError(1,fct,"data error in 'forces'!");
	     
	     n = intDiv(iTmp.n,1+ndf);
	     
	     frcNd.setDim(n);
             frc.setDim(n*ndf);
	     frcTmFct.setDim(n*ndf);

             for (i=0; i<n; i++)
	     {
	       frcNd[i] = iTmp[i*(1+ndf)];

               if (frcNd[i] < 1 || frcNd[i] > numnp)
		 prgError(2,fct,"invalid node number in 'forces'!");
	       
	       for (j=0; j<ndf; j++)
	       {
	         frcTmFct[i*ndf+j] = iTmp[i*(1+ndf)+1+j];
		 frc     [i*ndf+j] = dTmp[i*ndf+j];
	       }
	     }
	     
	     break;

    case 11: cout << "     MESH: reading nodal data output...\n\n";
	     
             while (1)
	     { 	     
               line.getNextLine(Ifile);

               nw   = line.split(&word); if (nw < 3) break;

	       tp   = word[0].which(wrndTypeNames); if (tp < 0) break;
	       
               if (nw < 4 && tp < 6) break;

	       if (!word[1].toInt(&indx,false)) break;

	       if (indx < 1) indx = 1;
               if (tp < 4 && indx > ndf) prgError(1,fct,"invalid index in 'nodal data output'!");
	       if (tp > 3 && tp < 6 && indx > ndm) 
		       prgError(2,fct,"invalid index in 'nodal data output'!");
	       if (tp == 6 && ndf < 2) prgError(3,fct,"error in 'nodal data output'!");
	      
	       if (!word[2].toDbl(&fact,false)) break;

               wrndNodeTmp.add(new Vector<int>);

	       i = 0; 
	       while (i < nw-3)

                 if (!word[i+3].toInt(wrndNodeTmp.lastItem().append(),false)) break; else i++; 

               if (i < nw-3) break;
  	      
               for (i=0; i<nw; i++) word[i].free(); delete [] word;

               wrndTypeTmp.append(tp);
	       wrndIndxTmp.append(indx);
	       wrndFactTmp.append(fact);
	     }
             for (i=0; i<nw; i++) word[i].free(); delete [] word;

             wrndNode = wrndNodeTmp;
             wrndType = wrndTypeTmp;
             wrndIndx = wrndIndxTmp;
             wrndFact = wrndFactTmp;

	     //cout << wrndType << "\n"; 
	     //cout << wrndIndx << "\n"; 
	     //cout << wrndFact << "\n"; 
	     //for (i=0; i<wrndNode.n; i++) cout << "  " << wrndNode[i] << "\n"; 

             break;
	     
    case 12: cout << "     MESH: reading ALE type ...\n\n";
	     
             if (!isALE()) prgError(4,fct,"this domain is not capable of ALE!");
	     
             elemProp.add(new PropertyItem(ALETYPE));

	     elemProp[elemProp.n-1].readInputData(Ifile,line,"input error in 'ALE type'!");

	     break;
	    
    case 13: cout << "     MESH: reading prescribed mesh motion ...\n\n";
	     
             if (!isALE()) prgError(6,fct,"this domain is not capable of ALE!");
	     
	     sprintf(tmp,"%di %df",1+ndm,ndm);  
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	     
             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,xbcIntTmp,xbcDblTmp,sTmp,lTmp))
		     prgError(1,fct,"data error in 'prescribed displacements'!");
	     
	     break;

    case -1: // go and inherit from GEOMETRY
	     
	     this->Geometry::readInputData(Ifile, line); 
	     
	     break;
  }
 
  return;
}









void Mesh::prepareInputData(void)
{
  // call ancestor function

  Geometry::prepareInputData();

  
  cout << "     MESH: prepare input data ...\n\n";
  
  char fct[] = "MESH::prepareInputData"; 
	  
  int e, i, j, jj, k, *splnNen;

  ElementGroup *eG;


// generate element groups 

  if (elemGrp.n < 1)
  {
    if (numElemGrp < 1) prgError(1,fct,"you need to define at least one element group!");
	  
    for (e=0; e<numElemGrp; e++)  elemGrp.add(new ElementGroup((Domain*)this));
  }
 
// intialise element types

  char *elmTypeNames[] = ELEMENT_TYPE_NAMES;
  
  prepareElemProp(ELEMENTTYPE,elmTypeNames);


// intialise ALE types

  if (isALE())
  {  
    char *aleTypeNames[] = ALE_TYPE_NAMES;
  
    prepareElemProp(ALETYPE,aleTypeNames);
  }

// get nen for element groups from dummy element

  Element *dmyElem;

  for (e=0; e<elemGrp.n; e++)
  {
     dmyElem = newElement(elemGrp[e].elemProp[ELEMENTTYPE]->id);

     elemGrp[e].nen = dmyElem->nen();

     delete dmyElem;
  }

// generate data for adaptivity

  if (grpType.n == 1 && grpType[0] == -1) grpType.free();

  if (ndm == 2)
  {
    // pointer from element group to surface

    if (grpType.n != surface.n) prgError(1,fct,"grpType.n != surface.n");

    for (e=0; e<grpType.n; e++) elemGrp[grpType[e]].geomObj.append(e);

    // assign nen to surfaces

    for (e=0; e<elemGrp.n; e++) 
    {
      eG = &(elemGrp[e]);
      for (i=0; i<eG->geomObj.n; i++) surface[eG->geomObj[i]].dat1 = eG->nen;
    }
  }
  else if (ndm == 3) ;  // to be done one day!

  // assign nen to splines

  for (i=0; i<spline.n; i++) spline[i].dat1 = 0;

  for (e=0; e<surface.n; e++)
  {
    for (i=0; i<surface[e].suppSpln.n; i++)
    {
      splnNen = &(((GeomSpline*)(surface[e].suppSpln[i]))->dat1);
      if (*splnNen < surface[e].dat1) *splnNen = surface[e].dat1;
    }
  }

// generate elements
 
  elem = new Element* [numel];

  for (e=0; e<numel; e++)
  {	  
    eG = &(elemGrp[ixTmp(e+1,nen+1)]);

    elem[e] = newElement(eG->elemProp[ELEMENTTYPE]->id);

    for (i=0; i<elem[e]->nen(); i++) 
    {
      if (ixTmp(e+1,i+1) < 1 || ixTmp(e+1,i+1) > numnp) 
        prgError(1,fct,"invalid node number in 'elements'!");
      elem[e]->ix[i] = ixTmp(e+1,i+1);
    }
    elem[e]->belongsTo(eG,this);

    if (!elem[e]->forDomainType(domain.domainType(this))) 
	    prgError(2,fct,"element and domain type inconsistent!");
  }

  ixTmp.free();

// generate node - element connectivity
 
  generateNodeElemConnectivity();

// generate node - node connectivity

  generateNodeNodeConnectivity();

// 2D: generate list of boundary nodes
// 3D: generate list of boundary faces and boundary nodes

  Vector<int> tmpBndNd, tmpBnd;

  if      (ndm == 2) generateBoundaryData2D(tmpBndNd,tmpBnd);
 
  else if (ndm == 3) generateBoundaryData3D(tmpBndNd,tmpBnd);

  // from Vector to arrays

  nBnd   = tmpBnd.n;
  nBndNd = tmpBndNd.n;

  bnd   = new int* [nBnd+1];
  bndNd = new int  [nBndNd+1];

  bnd[0] = bndNd; for (i=1; i<nBnd; i++) bnd[i] = bndNd + tmpBnd[i]; bnd[nBnd] = bndNd + nBndNd;

  for (i=0; i<nBndNd; i++) bndNd[i] = tmpBndNd[i];
  /*
  cout << nBndNd << " = nBndNd\n\n";
  cout << nBnd   << " = nBnd\n\n";
  for (i=0; i<nBnd; i++)
  {
    for (j=bnd[i]-bnd[0]; j<bnd[i+1]-bnd[0]; j++) cout << bndNd[j] << ",";
    cout << "\n";
  }
  */
  tmpBndNd.free();
  tmpBnd.free();

// turn prescribed displacements into the format of dependent displacements

  initialiseUDep(&uDep);
  
  if (isALE()) initialiseUDep(&xDep);

// calculate idu (take into account dependent displacements)

  for (i=0; i<uDep.n; i++)  idu(uDep[i].nd,uDep[i].dof) = - 1 - i;

  nequ = 0;
  for (i=1; i<numnp+1; i++)
    for (j=1; j<ndf+1; j++)
      if (idu(i,j) == 0) idu(i,j) = ++nequ; else if (idu(i,j) > 0) idu(i,j) = 0;

  prepareAndCheckUDep(&uDep);

 
// calculate idx (take into account dependent node motion)

  neqx = 0;

  if (isALE())
  {
    for (i=0; i<xDep.n; i++)  idx(xDep[i].nd,xDep[i].dof) = - 1 - i;

    for (i=1; i<numnp+1; i++)
      for (j=1; j<ndm+1; j++)
        if (idx(i,j) == 0) idx(i,j) = ++neqx; else if (idx(i,j) > 0) idx(i,j) = 0;

    prepareAndCheckUDep(&xDep);
  }
  else
  {
    for (i=0; i<numnp*ndm; i++) idx.x[i] = 0;
  }

  return;
}








void Mesh::prepareInteractions(void)
{
  // go and inherit from ancestors

  Geometry::prepareInteractions();


  cout << "     MESH: preparing interactions ...\n\n"; 
      
  int i, j, k;

// replace time function ids
	
  // uDep and xDep

  replaceUDepTmFct(&uDep);

  replaceUDepTmFct(&xDep);

  // replace time function ids in frcTmFct

  replaceFrcTmFct();


  // pressure load element data
 
  for (i=0; i<elemGrp.n; i++)
  {
    if (elemGrp[i].elemProp[ELEMENTTYPE]->name.contains("PressureLoad")) 
    {
      k = roundToInt(elemGrp[i].elemProp[ELEMENTTYPE]->data[0]);

      j = 0; while (j<timeFunction.n && k!=timeFunction[j].id) j++;

      if (j == timeFunction.n) prgError(2,"Mesh::prepareInteractions",
                        "invalid time function id in 'Element*D*nodedPressureLoad'!");

      elemGrp[i].elemProp[ELEMENTTYPE]->data[0] = (double)j;
    }
  }


// other stuff ....

  
  return;
}








void Mesh::updateUDepIncrements(void)
{
  updateUDepInc_hlp(uDep,idu);
  
  updateUDepInc_hlp(xDep,idx);
  
  return;	
}









bool Mesh::isALE(bool flag)
{
  int i=0;

  while (i<elemGrpProp.n && elemGrpProp[i]!=ALETYPE) i++;

  if (i<elemGrpProp.n) { if (!flag) return true; else return (neqx > 0); }

  return false;
}








void Mesh::addElemGrpProp(int propType)
{
  int i=0;

  while (i < elemGrpProp.n) if (elemGrpProp[i] != propType) i++; else break;
  if (i < elemGrpProp.n) prgError(1,"Mesh::addElemGrpProp","this property is already defined!");
  
  elemGrpProp.append(propType);
	
  return;
}















void Mesh::prepareElemProp(int propType,char **propTypeNames)
{
  // elemGrpPropTmp  - List< Vector<int> > read from input file 'element groups',
  //                   to be deleted permanently in prepareInputData
  // elemGrpProp     - Vector<int>; list of property types admissible for this domain type
  // elemProp        - List<PropertyItem>; list of all property items (possibly different types)
  // elemGrp         - List<ElementGroup>; list of element groups

  char fct[] = "Mesh::prepareElemProp";

  int i, j=0, e;
  
  if (elemGrpPropTmp.n < 1) prgError(1,fct,"'element groups' missing!");

  // ensure that this propType is required for this domain type;
  // set j to pos. of propType in elemGrpProp
 
  while (j < elemGrpProp.n) if (elemGrpProp[j] != propType) j++; else break;
  if (j >= elemGrpProp.n) prgError(2,fct,"invalid element property!");
  
  // set pointers to property items in element groups

  for (e=0; e<numElemGrp; e++)
  { 
    // identify element property with correct type and id

    i = 0; while (i<elemProp.n && !(elemProp[i].type==propType && 
		  elemProp[i].id==elemGrpPropTmp[e*elemGrpProp.n+e+1+j])) i++;
    if (i == elemProp.n) prgError(2,fct,"unknown element property type id!");

    // unset value in elemGrpPropTmp

    elemGrpPropTmp[e*elemGrpProp.n+e+1+j] = -1;

    // set pointer

    elemGrp[e].addElemProp(propType,&(elemProp[i]));
  }

  // set id values of all property items of propType to index of property type name

  for (i=0; i<elemProp.n; i++)
  {
    if (elemProp[i].type==propType) 
    { 
      elemProp[i].id = elemProp[i].name.which(propTypeNames);
      if (elemProp[i].id < 0) prgError(3,fct,"unknown element property type name!");
    }
  }

  // delete elemGrpPropTmp if all values have been unset

  e = 0;
  while (e < numElemGrp)
  {
    i = 0;
    while (i < elemGrpProp.n) if (elemGrpPropTmp[e*(elemGrpProp.n+1)+i+1] < 0) i++; else break;

    if (i < elemGrpProp.n) break; else e++;
  }
  if (e == numElemGrp) elemGrpPropTmp.free();

  return;
}







Element* Mesh::newElement(int type)
{
  switch (type)
  {
    case  0: return (Element*) new Element2D3nodedTriangle;
    case  1: return (Element*) new Element2D6nodedTriangle;
    case  2: return (Element*) new Element2D4nodedQuadrilateral;
    case  3: return (Element*) new Element2D8nodedQuadrilateral;
    case  4: return (Element*) new Element2D9nodedQuadrilateral;
    case  5: return (Element*) new Element3D4nodedTetrahedron;
    case  6: return (Element*) new Element3D10nodedTetrahedron;
    case  7: return (Element*) new Element3D8nodedBrick;
    case  8: return (Element*) new Element3D20nodedBrick;
    case  9: return (Element*) new Element2D3nodedStabIncompFluid;
    case 10: return (Element*) new Element2D3nodedStabIncompHighReFluid;
    case 11: return (Element*) new Element2D3nodedStabCompFluid;
    case 12: return (Element*) new Element2D3nodedLinearSolid;
    case 13: return (Element*) new Element2D6nodedQuadraticSolid;
    case 14: return (Element*) new Element2D4nodedLinearSolid;
    case 15: return (Element*) new Element2D4nodedFBarSolid;
    case 16: return (Element*) new Element2D8nodedQuadraticSolid;
    case 17: return (Element*) new Element2D9nodedQuadraticSolid;
    case 18: return (Element*) new Element3D4nodedStabIncompFluid;
    //case 19: return (Element*) new Element3D4nodedStabIncompHighReFluid;
    //case 20: return (Element*) new Element3D4nodedStabCompFluid;
    case 21: return (Element*) new Element3D4nodedLinearSolid;
    case 22: return (Element*) new Element3D10nodedQuadraticSolid;
    case 23: return (Element*) new Element3D8nodedLinearSolid;
    case 24: return (Element*) new Element3D8nodedFBarSolid;
    case 25: return (Element*) new Element3D20nodedQuadraticSolid;
    case 26: return (Element*) new Element2D2nodedLine;
    case 27: return (Element*) new Element2D3nodedLine;
    case 28: return (Element*) new Element2D2nodedGeomExSmallStrainBeam;
    case 29: return (Element*) new Element2D2nodedKirchhoffBeam;
    case 30: return (Element*) new Element2D2nodedTruss;
    case 31: return (Element*) new Element1D2nodedLine;
    case 32: return (Element*) new Element1D2nodedAdvectionDiffusion;
    case 33: return (Element*) new Element1D2nodedSchroedingerCN;
    case 34: return (Element*) new Element1D2nodedSchroedingerST;
    case 35: return (Element*) new Element2D2nodedPressureLoad;
    case 36: return (Element*) new Element2D3nodedPressureLoad;
    case 37: return (Element*) new Element1D2nodedPipeFlowST;
    case 38: return (Element*) new Element1D2nodedFlexibleWing;
    case 39: return (Element*) new Element2D3nodedLinearPoisson;
    case 40: return (Element*) new Element1D2nodedLinearSolidALE;
    case 41: return (Element*) new Element2D2nodedFreeSurface;

    default: prgError(1,"Mesh::newElement","unknown element type name!"); return NULL;
  }
}







void Mesh::printInfo(void)
{
  int e;
	
  COUT << "ndm ...... = " << ndm        << "\n";
  COUT << "ndf ...... = " << ndf        << "\n";
  COUT << "numnp .... = " << numnp      << "\n";
  COUT << "numel .... = " << numel      << "\n";
  COUT << "nequ       = " << nequ       << "\n";
  COUT << "neqx       = " << neqx       << "\n";
  COUT << "numElemGrp = " << numElemGrp << "\n\n";

  //for (e=0; e<numel; e++) elem[e]->print();
  
  return;
}








void Mesh::findMinMaxX(double *xmn, double *xmx, bool defFlg)
{
  int i, j, k, l;

  double *X = x.x; if (!defFlg) X = x0.x;

  if (numnp < 1) return;

  if (ndm > 1 && nen > 2)
  {
    k = (bndNd[0]-1) * ndm; 

    for (j=0; j<ndm; j++)
    { 
      xmn[j] = X[k+j];
      xmx[j] = X[k+j];
    }
    
    for (i=1; i<bnd[1]-bnd[0]; i++) 
    {
      k = (bndNd[i]-1) * ndm;
      for (j=0; j<ndm; j++) 
      {
        xmn[j] = min(xmn[j],X[k+j]);
        xmx[j] = max(xmx[j],X[k+j]);
      }
    }
  }
  else
  {
    for (j=0; j<ndm; j++) { xmn[j] = X[j]; xmx[j] = X[j]; }

    for (i=0; i<numnp; i++)
    {
      for (j=0; j<ndm; j++)
      {
        xmn[j] = min(xmn[j],X[i*ndm+j]);
        xmx[j] = max(xmx[j],X[i*ndm+j]);
      }
    }
  }

  return;
}








void Mesh::findMinMaxU(int indx, double &umn, double &umx)
{
  int i, ii = indx, j;

  if (ii < 1 || numnp < 1) return;
  
  if (ii > ndf) 
  {
    umn = outp[0];
    umx = outp[0];

    for (i=1; i<numnp; i++)
    {
      if (umn > outp[i]) umn = outp[i];
      if (umx < outp[i]) umx = outp[i];
    }			  

    return;
  }
  
  umn = u(1,ii);
  umx = u(1,ii);
  
  for (i=2; i<=numnp; i++) 
  {
    if (umn > u(i,ii)) umn = u(i,ii);
    if (umx < u(i,ii)) umx = u(i,ii);
  }

  return;
}







void Mesh::plotMesh(bool defm, bool flag)
{
  if (ndm == 3) { surf3D->draw(8,plot.currStdColour,false,flag,defm); return; }

  int e;

  if (!flag) { for (e=0; e<numel; e++) elem[e]->plotOutline(defm); return; }

  int i, j;

  double *x1, *x2, *X = x.x;

  if (defm) X = x0.x;

  for (i=0; i<nBnd; i++)
  {
    x1 = X + ndm*(bndNd[bnd[i+1]-bnd[0]-1]-1);

    for (j=bnd[i]-bnd[0]; j<bnd[i+1]-bnd[0]; j++)
    {
      x2 = X + ndm*(bndNd[j]-1);
      plot.line(x1,x2);
      x1 = x2;
    }
  }
  return;
}





void Mesh::paintElemGrp(int iEG, bool defFlg)
{
  int e;
	
  if (iEG>-1 && iEG<numElemGrp) 
  {
    List<Element> &el = elemGrp[iEG].elem;
    
    for (e=0; e<el.n; e++)  el[e].paint(defFlg);
  }
  else 
  {
    for (e=0; e<numel; e++)  elem[e]->paint(defFlg);
  }
  
  return;
}







void Mesh::plotNodes(int num, int which, bool defFlg)
{
  int i, j, k;
	
  double dpt = (plot.dAct[0] + plot.dAct[1]) * .0055;
	 
  MatrixFullArray<double> *X;
	  
  if (defFlg) X = &x; else X = &x0;

  switch (which)
  {
    case 1: if (ndm < 3) // all
            {
              for (i=1; i<=numnp; i++)
              {
                if (num == 1) plot.point(&((*X)(i,1)),dpt,i);

                else plot.point(&((*X)(i,1)),dpt);
              }
            }
            else
            {
              for (i=0; i<nBndNd; i++)
              {                
                if (num == 1) surf3D->plotNodePoint(i,dpt,bndNd[i]);

                else surf3D->plotNodePoint(i,dpt);
              }
            }
            break;
       
    case 2: for (i=0; i<uDep.n; i++)  // masters
       {
	 for (j=0; j<uDep[i].master.n; j++)
	 {
       	   k = uDep[i].masterNd[j]; 
	 
           if (num == 1) plot.point(&((*X)(k,1)),dpt,k); else plot.point(&((*X)(k,1)),dpt);
	 }
       }
       break;

    case 3: for (i=0; i<uDep.n; i++)  // dependents
       {
	 k = uDep[i].nd; 
	 
         if (num == 1) plot.point(&((*X)(k,1)),dpt,k); else plot.point(&((*X)(k,1)),dpt);
       }
       break;

    default: prgError(1,"Mesh::plotNodes","invalid option 'which'!");
  }
  
  return;
}







void Mesh::plotElemNum(bool defFlg)
{
  int e;
	
  char strg[15];
  
  for (e=0; e<numel; e++)  
  {
    sprintf(strg,"%d",e+1);
	  
    elem[e]->putLabel(strg,defFlg);
  }
  
  return;
}







void Mesh::plotBoun(bool defFlg)
{
  int i, j, *IDU = idu.x, indm, indf,
      symb[] = { XDASH, YDASH, ZDASH, SQUARE, TRIANGLE, CIRCLE };

  if (ndf < ndm) for (i=0; i<3; i++) symb[i]     = symb[i+3];
  else           for (i=0; i<3; i++) symb[ndm+i] = symb[3+i];
  
  double dpt = (plot.dAct[0] + plot.dAct[1]) * .01, *X = x0.x;
	  
  if (defFlg) X = x.x;

  if (ndm < 3)
  {
    for (i=0; i<numnp; i++)
    {
      for (j=0; j<ndf; j++)
      {
        if (*IDU < 1)  plotBoun_hlp(X,dpt,symb[j]);
        IDU++;
      }
      X+=ndm;
    }
  }
  else
  {
    for (i=0; i<nBndNd; i++)
    {
      if (surf3D->nodeIsReallyVisible(i))
      {
        indm = (bndNd[i]-1) * ndm;
        indf = (bndNd[i]-1) * ndf;

        for (j=0; j<ndf; j++) if (*(IDU+indf+j) < 1) plotBoun_hlp(X+indm,dpt,symb[j]);
      }
    }
  }	
  return;
}







void Mesh::plotFixed(bool defFlg)
{
  int i, j, *IDX = idx.x, indm,
      symb[] = { XDASH, YDASH, ZDASH };

  double dpt = (plot.dAct[0] + plot.dAct[1]) * .01, *X = x0.x;
	  
  if (defFlg) X = x.x;

  if (ndm < 3)
  {
    for (i=0; i<numnp; i++)
    {
      for (j=0; j<ndm; j++)
      {
        if (*IDX < 1)  plotBoun_hlp(X,dpt,symb[j]);
        IDX++;
      }
      X+=ndm;
    }
  }
  else
  {
    for (i=0; i<nBndNd; i++)
    {
      if (surf3D->nodeIsReallyVisible(i))
      {
        indm = (bndNd[i]-1) * ndm;

        for (j=0; j<ndm; j++) if (*(IDX+indm+j) < 1) plotBoun_hlp(X+indm,dpt,symb[j]);
      }
    }
  }	
  return;
}









void Mesh::plotLoad(double scl, bool defFlg)
{
  if (ndf < ndm) prgWarning(1,"Mesh::plotLoad","ndf < ndm  --> macro ignored!");
 
  int    i, j;
  
  double f[3], mxFrc = 0., sclFct;

  if (scl < 1.e-14)
  {	  
    for (i=0; i<frcNd.n; i++)

      for (j=0; j<ndf; j++)  
      {
	f[0] = fabs(frc[i*ndf+j]) * timeFunction[frcTmFct[i*ndf+j]].prop;
	
	if (mxFrc < f[0]) mxFrc = f[0];
      }
  }
  else mxFrc = scl;

  if (mxFrc < 1.e-14) return;
  
  sclFct = (plot.dAct[0] + plot.dAct[1]) * .15 / mxFrc;

  MatrixFullArray<double> *X;
	  
  if (defFlg) X = &x; else X = &x0;

  for (i=0; i<frcNd.n; i++)  
  {
    for (j=0; j<ndm; j++)  f[j] = frc[i*ndf+j] * timeFunction[frcTmFct[i*ndf+j]].prop;
	   
    plot.arrow(&(*X)(frcNd[i],1),f,sclFct);
  }

  return;
}









void Mesh::plotReac(double scl, bool defFlg)
{
  if (ndf < ndm) prgWarning(1,"Mesh::plotReac","ndf < ndm  --> macro ignored!");
 
  int    i, j, *idu1;
  
  double f[3], mxFrc = 0., sclFct, *reac1 = &(reac(1,1));

  if (scl < 1.e-14)
  {	  
    for (i=0; i<numnp; i++)
    {
      idu1 = &(idu(i+1,1));
	    
      for (j=0; j<ndm; j++)  
      {
	if (idu1[j] < 1)
	{
	  f[0] = fabs(reac1[i*ndf+j]);    if (mxFrc < f[0]) mxFrc = f[0];
	}
      }
    }
  }
  else mxFrc = scl;

  if (mxFrc < 1.e-14) return;
  
  sclFct = (plot.dAct[0] + plot.dAct[1]) * .15 / mxFrc;

  MatrixFullArray<double> *X;
	  
  if (defFlg) X = &x; else X = &x0;

  for (i=0; i<numnp; i++)  
  {
    idu1 = &(idu(i+1,1));
	  
    j = 0; while (j<ndm && idu1[j]>0) j++;

    if (j<ndm) 
    {
      for (j=0; j<ndm; j++)  f[j] = - reac1[i*ndf+j];

      if (i==231) cout << f[0] << " " << f[1] << "\n";
      
      plot.arrow(&(*X)(i+1,1),f,sclFct);
    }
  }

  return;
}





void Mesh::printNodalData(int node)
{
  int i;
	
  if (node < 0 || node > numnp) { COUT << "invalid node number!\n\n"; return; }
  
  if (node != 0)
  {
    COUT << domain.name(this) << ": node " << node << "\n";
    COUT << "-----------------------------------------\n";
    COUT; printf("x0 = %11.4g",x0(node,1)); 
          for(i=1; i<ndm; i++) printf(" %11.4g",x0(node,1+i));  cout << "\n";
    COUT; printf("x  = %11.4g", x(node,1)); 
          for(i=1; i<ndm; i++) printf(" %11.4g", x(node,1+i));  cout << "\n";
    COUT; printf("u  = %11.4g", u(node,1)); 
          for(i=1; i<ndf; i++) printf(" %11.4g", u(node,1+i));  cout << "\n";
    COUT; printf("un = %11.4g",un(node,1)); 
          for(i=1; i<ndf; i++) printf(" %11.4g",un(node,1+i));  cout << "\n";
    COUT; printf("u3 = %11.4g",u3(node,1)); 
          for(i=1; i<ndf; i++) printf(" %11.4g",u3(node,1+i));  cout << "\n";
    COUT; printf("u4 = %11.4g",u4(node,1)); 
          for(i=1; i<ndf; i++) printf(" %11.4g",u4(node,1+i));  cout << "\n";
    COUT; printf("u5 = %11.4g",u5(node,1)); 
          for(i=1; i<ndf; i++) printf(" %11.4g",u5(node,1+i));  cout << "\n";
    COUT; printf("u6 = %11.4g",u6(node,1)); 
          for(i=1; i<ndf; i++) printf(" %11.4g",u6(node,1+i));  cout << "\n\n";


    return;
  }
  
  COUT << "what is going on?\n\n";

  return;
}






void Mesh::contourPlot(int var, int Indx, int nCol, bool flg, double &umin, double &umax, bool legd)
{
  int i, indx, e;
	
  double umn, umx, *U;
 
  switch (var)
  {
    case 1: // degree of freedom

            if (Indx > ndf) { prgWarning(1,"Mesh::contourPlot","index out of bounds!"); return; }
	    
            indx = Indx - 1;
	    
	    U = &(u(1,1));

	    umn = U[indx]; for (i=0; i<numnp; i++) if (U[i*ndf+indx] < umn) umn = U[i*ndf+indx];
	    umx = U[indx]; for (i=0; i<numnp; i++) if (U[i*ndf+indx] > umx) umx = U[i*ndf+indx];
	    
	    break;

    case 2: // initial coordinate

            if (Indx > ndm) { prgWarning(2,"Mesh::contourPlot","index out of bounds!"); return; }
	    
            indx = Indx - 1;
	    
	    U = &(x0(1,1));

	    umn = U[indx]; for (i=0; i<numnp; i++) if (U[i*ndm+indx] < umn) umn = U[i*ndm+indx];
	    umx = U[indx]; for (i=0; i<numnp; i++) if (U[i*ndm+indx] > umx) umx = U[i*ndm+indx];
	    
	    break;
	    
    case 3: // current coordinate

            if (Indx > ndm) { prgWarning(3,"Mesh::contourPlot","index out of bounds!"); return; }
	    
            indx = Indx - 1;
	    
	    U = &(x(1,1));

	    umn = U[indx]; for (i=0; i<numnp; i++) if (U[i*ndm+indx] < umn) umn = U[i*ndm+indx];
	    umx = U[indx]; for (i=0; i<numnp; i++) if (U[i*ndm+indx] > umx) umx = U[i*ndm+indx];

	    break;
	    
    case 4: // last projection

	    U = &(outp[0]);

	    umn = U[0]; for (i=1; i<numnp; i++) if (U[i] < umn) umn = U[i];
	    umx = U[0]; for (i=1; i<numnp; i++) if (U[i] > umx) umx = U[i];

	    break;
  }

  if (!flg) { umin = umn; umax = umx; }
  
  COUT << "contour plot      min           max\n";
  printf("             actual  %12.5g  %12.5g\n", umn, umx);
  printf("             plot    %12.5g  %12.5g\n\n", umin,umax);
  
  if (umin >= umax - 1.e-25) 
  { 
    COUT << "umin >= umax!  contour plot aborted!\n\n";
    return;
  }
  
  if (ndm == 2)

    for (e=0; e<numel; e++)  elem[e]->contourPlot(var,Indx,nCol,umin,umax);

  else if (ndm == 3)
  {
    surf3D->updateU(var,Indx);

    surf3D->contourPlot(nCol);
  }
  if (legd) plot.contourPlotLegend(umin,umax,umn,umx,flg,nCol);

  return;
}








void Mesh::projectToNodes(MyString &what, int indx, int indx2)
{
  int e, ee, i, nen, *ix, *nE;
	
  double fct, V;
  
  char *datTypeList[] = { "stress",
                          "intvar",
                          "vorticity",
                          "error",
                          "gradient",
                          "normOfGradient",
                          "normOfGradientSquared", NULL };
  
  int datType = what.which(datTypeList);
 
  if ((datType == 1 && !(isSolid(*this) || isMicroCellWulf(*this) || isMicroCell(*this)))
   || (datType == 2 && !isFluid(*this))) 
    { COUT << "'" << what << "' can not be projected for this domain type!\n\n"; return; }
  
  for (i=0; i<numnp; i++)  outp[i] = 0.;

  for (e=0; e<numel; e++)
  {
    switch (datType)
    {
      case  0: elem[e]->projectStress(indx);  break;
	      
      case  1: elem[e]->projectIntVar(indx);  break;
	     
      case  2: elem[e]->projectVorticity(indx); break;

      case  3: elem[e]->projectError(indx); break;

      case  4: elem[e]->projectGradient(indx,indx2); break;

      case  5: elem[e]->projectNormOfGradient(indx); break;

      case  6: elem[e]->projectNormOfGradientSquared(indx); break;

      default: COUT << "what shall we do with the drunken sailor ?\n\n";  break;
    }

    nen = elem[e]->nen();
    
    ix   = elem[e]->ix;
   
    V    = elem[e]->volume();

    if (V > 1.e-40)
    {
      for (i=0; i<nen; i++)
      {
        // calculate volumes of elements surrounding node i

        fct = 0.;
        
        nE = &(nodeElem[ix[i]-1][0]);
        
        for (ee=0; ee<nodeElem[ix[i]-1].n; ee++) fct += elem[nE[ee]]->volume();
        
        //cout << nodeElem[ix[i]-1] << " - " << ix[i]-1 << "\n";
        
        fct = V / fct;
        
        // assemble 

        outp[ix[i]-1] += p[i] * fct;
      }
    }
  }

  return;
}






void Mesh::writeNodalData(void)
{
  char fct[] = "Mesh::writeNodalData";

  int         i, j, indx;
  double      val, fact,   a1, b1, a2, b2, a1n, b1n, a2n, b2n;
  char        tmp[20];
  MyString    tmpStr;
  VectorArray<int> *node;
  
  for (i=0; i<wrndType.n; i++)
  {
    val  = 0.;
    fact = wrndFact[i];
    indx = wrndIndx[i];
    node = &(wrndNode[i]);
    
    switch (wrndType[i])
    {
      case  0: // u 
	     
               for (j=0; j<node->n; j++) val += u((*node)[j],indx); break;
	       
      case  1: // du 
	      
	       for (j=0; j<node->n; j++) val += u3((*node)[j],indx); break;
	       
      case  2: // ddu 
	      
	       for (j=0; j<node->n; j++) val += u5((*node)[j],indx); break;
	       
      case  3: // reac
	      
	       for (j=0; j<node->n; j++) val += reac((*node)[j],indx); break;
	       
      case  4: // x
	      
	       for (j=0; j<node->n; j++) val += x((*node)[j],indx); break;
	       
      case  5: // x0 
	      
	       for (j=0; j<node->n; j++) val += x0((*node)[j],indx); break;

      case  6: // int(uu+vv) (main purpose: integration of probability for Schroedinger equation)

	       if (ndm != 1) { prgWarning(1,fct,"'int(uu+vv)' available in 1D only!"); break; }

	       fact /= 3.;
	       
	       for (j=1; j<numnp; j++) val += (x(j+1,1)-x(j,1))
		       
		    //    * (u(j+1,1)*u(j+1,1)  + u(j,1)*u(j,1)
		    //     + u(j+1,2)*u(j+1,2)  + u(j,2)*u(j,2));
	       
		        * (u(j+1,1)*u(j+1,1) + u(j+1,1)*u(j,1) + u(j,1)*u(j,1)
			 + u(j+1,2)*u(j+1,2) + u(j+1,2)*u(j,2) + u(j,2)*u(j,2));
	       break; 

      case  8: // int(uu+vv) (main purpose: integration of probability for Schroedinger equation)

	       // this is the integral over the space time slab devided by dt

	       if (ndm != 1) { prgWarning(1,fct,"'int(uu+vv)' available in 1D only!"); break; }

	       fact /= 18.;
	       
	       for (j=1; j<numnp; j++) 
	       {
		 a1  = u(j,1);
		 b1  = u(j,2);
		 a2  = u(j+1,1);
		 b2  = u(j+1,2);
		 a1n = u(j,1);
		 b1n = u(j,2);
		 a2n = u(j+1,1);
		 b2n = u(j+1,2);
		  	 
		       val += (x(j+1,1)-x(j,1)) *
	((2*power(a1,2) + 2*power(a1n,2) + a1*(2*(a1n + a2) + a2n) + a1n*(a2 + 2*a2n) + 
       2*(power(a2,2) + a2*a2n + power(a2n,2)) + 2*power(b1,2) + 2*power(b1n,2) + b1n*b2 + 
       2*power(b2,2) + 2*(b1n + b2)*b2n + 2*power(b2n,2) + b1*(2*(b1n + b2) + b2n)));
	       }
	       break; 

      case  7: // outp  (whatever it holds at the moment)

	       for (j=0; j<node->n; j++) val += outp[(*node)[j]-1]; break;

	       
      default: prgError(1,fct,"unknown wrndType!");	     
    }
    val *= fact;
    sprintf(tmp," %12.5g",val);
    tmpStr.append(tmp);
  }
  
  prgWriteToTFile(tmpStr);
 
  return;
}






void Mesh::plotGaussPoints(int num, bool defFlg)
{
  for (int e=0; e<numel; e++)  elem[e]->plotGaussPoints(num,defFlg);
  
  return;
}





void Mesh::plotU1D(int indx)
{
  if (ndm > 1) return;

  double x1[3], x2[3];
  
  int e, i, i1, i2, ii = indx;

  if (ii < 1) return;

  for (e=0; e<numel; e++)
  {
    for (i=0; i<elem[e]->nen()-1; i++)
    {
      i1 = elem[e]->ix[i];
      i2 = elem[e]->ix[i+1];
	  
      x1[0] = x(i1,1);
      x2[0] = x(i2,1);

      if (ii > ndf)
      {
        x1[1] = outp[i1-1];
	x2[1] = outp[i2-1];	
      }
      else
      {      
        x1[1] = u(i1,ii);
        x2[1] = u(i2,ii);
      }

      plot.line(x1,x2);
    }
  }

  return;
}










void Mesh::reset(void)
{
  cout << " Mesh::reset():  what do you want?\n\n"; 
	
  return;
}








int Mesh::updateMesh(int cutType, bool printRes)
{
  char fct[] = "Mesh::updateMesh";
 
  int e, eG, iter, totalIter = 0, div = 0, 
      maxIter = 20, maxTotalIter = 5000, maxDiv = 10;
  bool firstIterMesh, fail = false;
  double rN, rNPrev, cutFct, cut;
  double *D = d.x, *Dconv = d0.x, *X = x.x, *Xn = xn.x, *V = v.x, *Vn = vn.x, *dd;

  List<Element> *el;
  Vector<double> cutSeq;
  int i, j, ii, *IDX = idx.x, indm, numnpXndm = numnp * ndm;
 
  if (!isALE()) { return 0; prgWarning(1,fct,"this is only for ALE domains!"); }

  if (neqx == 0) goto noDofs;

  computerTime.go(fct);

  cutSeq.append(1.);
  cutFct = 1.;
  cut = 1.;
 
  while (cutSeq.n > 0)
  { 
    firstIterMesh = true;
    iter          = 1;
    
    while (1)
    {
      if (solverMesh->currentStatus < INIT_OK) prgError(1,fct,"solver not initialised!"); 
          
      if (firstIterMesh) rN = -1.;

      solverMesh->zeroMtx();

      r.zero();
     
      reacMesh.zero();

      for (eG=0; eG<elemGrp.n; eG++)
      {
        el = &(elemGrp[eG].elem);
  	    
        for (e=0; e<el->n; e++)
        {
          fail = (*el)[e].calcStiffnessAndResidualMesh();

          if (fail) { break; }

          solverMesh->assembleElemMtx(&((*el)[e]));

          solverMesh->assembleElemVec(&((*el)[e]),firstIterMesh);
        }
      }
     
      rNPrev = rN; 
      rN     = r.norm();
      
      if (printRes) { COUT << domain.name(this); printf(" mesh  %6.4f %11.4e\n",cut,rN); }
      
      solverMesh->currentStatus = ASSEMBLY_OK;

      // check for convergence
      
      if (rN < tolMesh) { break; }

      // check number of iterations

      if (++iter > maxIter) { fail = true; break; }
      
      // factorise and solve

      dd = solverMesh->factoriseAndSolve(r.x);
      
      if (dd == NULL) { fail = true; break; }
     
      // update solution vector
     
      indm = 0;

      for (i=0; i<numnp; i++)  
      {	  	  
        for (j=0; j<ndm; j++)  {  ii = IDX[indm+j];  if (ii > 0)  D[indm+j] += dd[ii-1];  }

        indm += ndm; 
      }

      // update dependent displacements

      if (firstIterMesh)
      {
	for (i=0; i<xDep.n; i++) xDep[i].updateWithIncrement(D,ndm); 
        firstIterMesh = false;
      }
    }
    if (firstIterMesh)  // this is required only in non-standard testing situations ...
    {
      for (i=0; i<xDep.n; i++) xDep[i].updateWithIncrement(D,ndm);
      firstIterMesh = false;
    }
 
    totalIter += iter;
    
    ii = cutSeq.n;
    
    if (fail) // if no convergence
    {
      if (totalIter > maxTotalIter || div == maxDiv) break;
      div++;
      for (i=0; i<numnpXndm; i++) D[i] = Dconv[i]; 
      cutSeq[ii-1] *= 0.5;
      cut          -= cutSeq[ii-1];
      if (cutType == 1) cutSeq.append(cutSeq[ii-1]);
      else
      {
	cutFct = cutSeq[ii-1];
        j = roundToInt((1.-cut)/cutFct) + 1;
	for (i=0; i<j; i++) cutSeq[i] = cutFct;
      }
      cutFct = 0.5;
    }
    else // if converged
    {
      for (i=0; i<numnpXndm; i++) Dconv[i] = D[i];
      if (ii>1)
      {
        cutFct = cutSeq[ii-2] / cutSeq[ii-1];
	cut   += cutSeq[ii-2];
      }
      cutSeq.trunc(ii-1);
    }
    for (i=0; i<xDep.n; i++) if (xDep[i].master.n == 0) xDep[i].duc *= cutFct;
  }
 
  //if (printRes) cout << "\n";
  
  if (fail)
  { 
    COUT << "mesh update failed!\n\n";
    ctimUpdateMesh += computerTime.stop(fct);
    return 1;
  }
  
  if (printRes && totalIter > iter) 
    COUT << "total mesh update iter. steps = " << totalIter << "\n\n";

  for (i=0; i<numnpXndm; i++) 
  {
    if (IDX[i] != 0)
    {
      X[i] = Xn[i] + D[i];  
      V[i] = td[20] * D[i] + td[22] * Vn[i]; 
    }
  }

  if (ndm == 3) surf3D->updtFlagX = true;

  ctimUpdateMesh += computerTime.stop(fct);
  
  return 0;

  noDofs:

  for (i=0; i<xDep.n; i++) xDep[i].updateWithIncrement(D,ndm); 
 
  for (i=0; i<numnpXndm; i++) 
  {
    if (IDX[i] != 0)
    {
      X[i] = Xn[i] + D[i];  
      V[i] = td[20] * D[i] + td[22] * Vn[i]; 
    }
  }

  return 0;	
}











void Mesh::elementDiffStiffTestMesh(double ddd, int el, int dig, int dig2, bool gfrmt)
{
  int e = el - 1;

  if (e < 0)      e = 0; 
  if (e >= numel) e = numel - 1;

  COUT << "\n";
  COUT << "ALE-diff-test for " << domain.name(this) << ", element " << e+1 << "\n";
  COUT << "-------------------------------------------------------------\n";
  
  elem[e]->diffAleTest(ddd,dig,dig2,gfrmt);
  
  return;
}











double Mesh::getElemSizeOptInternal(double *xp)
{
  int i, e = pointInElement(xp);

  if (e < 1) prgError(1,"Mesh::getElemSizeOptInternal","no background mesh element found!");
  
  double N[100], h, d, val = elem[e-1]->getElemSizeOpt(xp,N);

  // check effect of point sources

  for (i=0; i<hPntSrc.n; i++)
  {
    d = sqrt(dist2(hPntSrc[i].x,xp,2));
    if (d < 4. * hPntSrc[i].h) h = hPntSrc[i].h;
    else                       h = (d - 4. * hPntSrc[i].h) * hPntSrc[i].maxGrad + hPntSrc[i].h;
    if (val > h) val = h;
  }

  return val;
}









void Mesh::calcElemSizeOpt(double hmn, double hmx, double errBar, bool keepSmallFlag) 
{ 
  int i;

  double fact;

  if (elemSizeOpt.n != numnp) elemSizeOpt.setDim(numnp);

  for (i=0; i<numnp; i++)
  {
    if (outp.x[i] < 0.) outp.x[i] = 0.;

    if (outp.x[i] > errBar) elemSizeOpt.x[i] = hmn;

    else 
    {
      fact = (errBar - outp.x[i]) / errBar; 

      elemSizeOpt.x[i] = hmn + (hmx - hmn) * fact * fact;
    }
    //cout << outp.x[i] << "->" << elemSizeOpt.x[i] << "\n";
  }

  if (!keepSmallFlag) return;

  for (i=0; i<numnp; i++) elemSizeOpt.x[i] = min(elemSizeOpt.x[i],elemSizeCurr.x[i]);

  return; 
}













void Mesh::smoothElemSizeOptDist(int nCircle, int nAcross, double maxGrad, double hmin)
{
  if (ndm != 2) prgError(1,"Mesh::smoothElemSizeOptDist","so far this is 2D only!");

  int i;

  GeomPoint *pntPtr;

  // update point.dat1 values

  if (geometryDiscretisedAtLeastOnce)
  {
    for (i=0; i<numnp; i++)
    {
      pntPtr = ndGmLnk[i].whichPoint();
      if (pntPtr != NULL) pntPtr->dat1 = i;
    }
  }

  // length of splines

  if (geometryDiscretisedAtLeastOnce) smoothElemSizeDist_hlp0();

  // curvature

  smoothElemSizeDist_hlp1(nCircle,maxGrad,hmin);

  // number of elements "across"

  smoothElemSizeDist_hlp2(nAcross,maxGrad,hmin);

  // ensure max gradient

  smoothElemSizeDist_hlp3(maxGrad);

  // ensure compatibility with elements to be kept

  if (geometryDiscretisedAtLeastOnce) smoothElemSizeDist_hlp4(maxGrad);

  return;
}












void Mesh::getInputElemSizeOpt(void)
{
  int i;

  elemSizeOpt.setDim(numnp);

  for (i=0; i<numnp; i++)  
  {
    elemSizeOpt.x[i] = u.x[i*ndf];

    if (elemSizeOpt.x[i] < 1.e-8) prgError(1,"Mesh::getInputElemSize","invalid element size!");

     u.x[i*ndf] = 0.;
    un.x[i*ndf] = 0.;
  }

  return;
}











void Mesh::calcElemSizeCurr(void)
{
  if (ndm != 2) prgError(1,"Mesh::calcElemSizeCurr","so far this is for 2D only");

  int i, e;

  if (elemSizeCurr.n != numnp) elemSizeCurr.setDim(numnp);

  for (i=0; i<numnp; i++)
  {
    elemSizeCurr[i] = 0.;

    for (e=0; e<nodeElem[i].n; e++) elemSizeCurr[i] += elem[nodeElem[i][e]]->diameter();

    elemSizeCurr[i] /= (double) nodeElem[i].n;
  } 

  return;
}











void Mesh::remeshElemGroups(bool showFlg)
{
  char fct[] = "Mesh::discretiseElemGroups";

  if (ndm != 2) prgError(1,fct,"so far this works only for 2D!");

  int c, i, j, k, n, dmy, jndm, jndf;

  GeomPoint *pntPtr;

  Vector<int> surf, spln, nBndNode;

  ListInfinite< MatrixFullArray<double> > xList;
              
  ListInfinite< MatrixFullArray<int> >   ixList;
              
  ListInfinite< ListInfinite<NodeGeomLink> > ndGmLnkList;

  // prepare lists of geometry objects to be discretised

  for (i=0; i<elemGrpToBeMeshed.n; i++)
    for (j=0; j<elemGrp[elemGrpToBeMeshed[i]].geomObj.n; j++) 
      surf.append(elemGrp[elemGrpToBeMeshed[i]].geomObj[j]);

  whichSplines(surf,spln);

  // update position of geometry points

  if (geometryDiscretisedAtLeastOnce)
  {
    for (i=0; i<ndGmLnk.n; i++)
    {
      if (ndGmLnk[i].id != i + 1) prgError(1,fct,"surprise, surprise!");
      pntPtr = ndGmLnk[i].whichPoint();
      if (pntPtr != NULL)
      {
        pntPtr->x[0] = x.x[i+i];
        pntPtr->x[1] = x.x[i+i+1];
      }
    }
  }

  // update h data of geometry points

  if (geometryDiscretisedAtLeastOnce)
  {
    for (i=0; i<point.n; i++)
    {
      if (ndGmLnk[point[i].dat1].whichPoint() != &(point[i])) prgError(1,fct,"fatal error!");
 
      point[i].h = elemSizeOpt[point[i].dat1];
    }
  }
  else // get initial h data for geometry points
  {
    for (i=0; i<point.n; i++)  point[i].h = getElemSizeOptInternal(point[i].x);

    for (i=0; i<spline.n; i++) spline[i].checkElemSizes();
  }

  // discretise splines

  computerTime.go("discretise splines");

  for (i=0; i<spln.n; i++) spline[spln[i]].discretise(this,spline[spln[i]].dat1);

  computerTime.stopAndPrint("discretise splines");

  // discretise surfaces

  computerTime.go("generate mesh");

  for (i=0; i<surf.n; i++)
    surface[surf[i]].generateMesh(xList[i],ixList[i],(void*)(&(ndGmLnkList[i])),dmy,dmy,
                                  *(nBndNode.append()),
                                  this,surface[surf[i]].dat1,showFlg);

  computerTime.stopAndPrint("generate mesh");

  // reorganise mesh and transfer nodal data and internal variables

  computerTime.go("join mesh patches together, do the data transfer");

  reorganiseMesh(xList,ixList,ndGmLnkList,surf);

  computerTime.stopAndPrint("join mesh patches together, do the data transfer");

  // XXXXX the following is taken from Mesh::prepareInputData XXXXX
  // XXXXX          and from Mesh::readInputData              XXXXX

  computerTime.go("prepare mesh");

  // generate node - element connectivity

  generateNodeElemConnectivity();

  // generate node - node connectivity

  generateNodeNodeConnectivity();

  // generate list of boundary nodes

  Vector<int> tmpBndNd, tmpBnd;

  generateBoundaryData2D(tmpBndNd,tmpBnd);

  if (bnd != NULL)   delete [] bnd;
  if (bndNd != NULL) delete [] bndNd;

  nBnd   = tmpBnd.n;
  nBndNd = tmpBndNd.n;

  bnd   = new int* [nBnd+1];
  bndNd = new int  [nBndNd+1];

  bnd[0] = bndNd; for (i=1; i<nBnd; i++) bnd[i] = bndNd + tmpBnd[i]; bnd[nBnd] = bndNd + nBndNd;

  for (i=0; i<nBndNd; i++) bndNd[i] = tmpBndNd[i];

  // associate geometry objects with sets of nodes (for boundary conditions etc, temporarily)

  for (j=0; j<numnp; j++) 

    for (i=0; i<ndGmLnk[j].geomObj.n; i++) 

      ((GeomObject*)(ndGmLnk[j].geomObj[i]))->nodeListTmp.append(j+1);

  // apply boundary conditions

  idx.setDim(numnp,ndm,true);
  idu.setDim(numnp,ndf,true);

  for (i=0; i<numnp*ndf; i++) idu.x[i] = 0;
  for (i=0; i<numnp*ndm; i++) idx.x[i] = 0;

  for (i=0; i<bndCond.n; i++)
  {
    for (j=0; j<((GeomObject*)(bndCond[i].geomObj))->nodeListTmp.n; j++)
    {
      n = ((GeomObject*)(bndCond[i].geomObj))->nodeListTmp[j] - 1;

      for (k=0; k<ndf; k++) if (bndCond[i].idu[k]) idu.x[n*ndf+k] = 1;
      for (k=0; k<ndm; k++) if (bndCond[i].idx[k]) idx.x[n*ndm+k] = 1;
    }
  }

  // apply prescribed displacements

  uDep.free();

  c = 0;

  for (i=0; i<presDisp.n; i++)
  {
    for (j=0; j<((GeomObject*)(presDisp[i].geomObj))->nodeListTmp.n; j++)
    {
      n = ((GeomObject*)(presDisp[i].geomObj))->nodeListTmp[j];
      {
        for (k=1; k<ndf+1; k++)
        {
          if (idu(n,k) != 0) 
          {
            uDep.add(new DependentDoF);
  	      
            uDep[c].nd      = n;
            uDep[c].dof     = k;
            uDep[c].ucBase  = presDisp[i].uBase[k-1];
            uDep[c++].tmFct = presDisp[i].tmFct[k-1];
          }
        }
      }
    }
  }

  // turn prescribed displacements into the format of dependent displacements

  initialiseUDep(&uDep);
  
  //initialiseUDep(&xDep);

  replaceUDepTmFct(&uDep);

  //replaceUDepTmFct(&xDep);
 
  // calculate idu (take into account dependent displacements)

  for (i=0; i<uDep.n; i++)  idu(uDep[i].nd,uDep[i].dof) = - 1 - i;
  
  nequ = 0;
  for (i=1; i<numnp+1; i++)
    for (j=1; j<ndf+1; j++)
      if (idu(i,j) == 0) idu(i,j) = ++nequ; else if (idu(i,j) > 0) idu(i,j) = 0;

  prepareAndCheckUDep(&uDep);
  
  // calculate idx (take into account dependent node motion)

  for (i=0; i<xDep.n; i++)  idx(xDep[i].nd,xDep[i].dof) = - 1 - i;

  neqx = 0;
  for (i=1; i<numnp+1; i++)
    for (j=1; j<ndm+1; j++)
      if (idx(i,j) == 0) idx(i,j) = ++neqx; else if (idx(i,j) > 0) idx(i,j) = 0;

  prepareAndCheckUDep(&xDep);

  // apply distributed loads

  int e, f0;

  Vector<double> fact, frcTmp;
  Vector<int> el, frcNdTmp, frcNdTmp1, frcTmFctTmp;
  List< Vector<int> > elNd, ifrc;

  for (i=0; i<distLoad.n; i++)
  {
    calcElemNodeDataAlongGeomObj(el,elNd,ifrc,frcNdTmp1,distLoad[i].geomObj);

    f0 = frcTmp.n;

    for (j=0; j<frcNdTmp1.n; j++) frcNdTmp.append(frcNdTmp1[j]);

    for (e=0; e<el.n; e++)
    {
      elem[el[e]-1]->getDistLoadFact(fact,elNd[e]);

      for (j=0; j<elNd[e].n; j++)  
        for (k=0; k<ndf; k++)
        { 
          frcTmp[f0+(ifrc[e][j])*ndf+k]     += fact[j] * distLoad[i].fdl[k];
          frcTmFctTmp[f0+(ifrc[e][j])*ndf+k] = distLoad[i].tmFct[k];
        }
    }
  }
  frc.setDim(frcTmp.n);      for (j=0; j<frcTmp.n; j++)   frc[j]      = frcTmp[j];
  frcNd.setDim(frcNdTmp.n);  for (j=0; j<frcNdTmp.n; j++) frcNd[j]    = frcNdTmp[j];
  frcTmFct.setDim(frcTmp.n); for (j=0; j<frcTmp.n; j++)   frcTmFct[j] = frcTmFctTmp[j];

  frcTmp.free();
  frcNdTmp.free();
  frcTmFctTmp.free();
 
  replaceFrcTmFct();

  //cout << frcNd << "\n" << frc << "\n" << frcTmFct << "\n\n";

  // apply nodal data output

  wrndNode.setDim(wrndOutp.n);
  wrndFact.setDim(wrndOutp.n);
  wrndIndx.setDim(wrndOutp.n);
  wrndType.setDim(wrndOutp.n);

  Vector<int> tmp;

  for (i=0; i<wrndOutp.n; i++)
  {
    for (k=0; k<wrndOutp[i].geomObj.n; k++)
    {
      for (j=0; j<((GeomObject*)(wrndOutp[i].geomObj[k]))->nodeListTmp.n; j++)
      {
        n = ((GeomObject*)(wrndOutp[i].geomObj[k]))->nodeListTmp[j];
        {
          if (!tmp.contains(n)) tmp.append(n);
        }
      }
    }
    wrndNode[i] = tmp;
    tmp.free();

    wrndFact[i] = wrndOutp[i].fact;
    wrndIndx[i] = wrndOutp[i].indx;
    wrndType[i] = wrndOutp[i].type;
  }

  // delete geometry object node lists

  for (i=0; i<point.n; i++) point[i].nodeListTmp.free();
  for (i=0; i<spline.n; i++) spline[i].nodeListTmp.free();
  for (i=0; i<surface.n; i++) surface[i].nodeListTmp.free();

 
  // domain type specific stuff:
  // -----------------------------

  // initialise internal variables

  if (isSolid(*this)) for (int e=0; e<numel; e++)  elem[e]->initialiseIntVar();



  computerTime.stopAndPrint("prepare mesh");

  // finish

  COUT << "---------------------\n";
  COUT << "numnp = " << numnp << "\n";
  COUT << "numel = " << numel << "\n\n";

  geometryDiscretisedAtLeastOnce = true;

  hPntSrc.free();

  return;
}








bool Mesh::elemSizeRatioOK(double ratio)
{
  int i, j; 

  Vector<int> delGeomObj;

  Vector<void*> keepGeomObj;

  // generate list of top geometry objects to be kept

  for (i=0; i<elemGrpToBeMeshed.n; i++)
    for (j=0; j<elemGrp[elemGrpToBeMeshed[i]].geomObj.n; j++) 
      delGeomObj.append(elemGrp[elemGrpToBeMeshed[i]].geomObj[j]);

  for (i=0; i<surface.n; i++)
  {
    j = 0; while (j < delGeomObj.n && delGeomObj[j] != i) j++;
    if (j >= delGeomObj.n) keepGeomObj.append((void*)(&(surface[i])));
  }

  // check

  i = 0;
  while (i < numnp)
  {
    if (!ndGmLnk[i].geomObj.containsAtLeastOneOf(keepGeomObj))
    {
      if (elemSizeCurr[i] > elemSizeOpt[i] * ratio) return false;
      if (elemSizeOpt[i] > elemSizeCurr[i] * ratio) return false;
    }
    i++;
  }
  
  return true;
}












void Mesh::transferNodalDataTest(double *xp)
{
  Vector<int>    node;

  Vector<double> factor;

  int i;

  getTransferDataInternal(xp,node,factor);

  plot.setColour(3);

  for (i=0; i<node.n; i++) 
  {
    plot.point(&(x.x[node[i]+node[i]-2]),plot.dPt(),node[i]);

    COUT << "node: " << node[i] << "  -> " << factor[i] << "\n";
  }
  cout << "\n";

  prgUpdateDisplay();  

  return;
}









void Mesh::writeMeshToFile(MyString &fileName, bool defm)
{
  std::ofstream out;
	
  int i, j;

  char tmp[500];

  double *X;
  
  if (defm) X = x.x;  else X = x0.x;

  out.open(fileName.asCharArray());

  if (!out) 
  { 
    prgWarning(1,"Mesh::writeMeshToFile","could not open file for writing!");
    return;
  }
  
  out << "% Mpap2 finite element mesh\n\n";

  sprintf(tmp,"coordinates | 123");
  for (i=0; i<ndm+ndf; i++) sprintf(&(tmp[17+i+i+i])," i0");
  sprintf(&(tmp[17+i+i+i])," %df\n", ndm);
  out << tmp;

  for (i=0; i<numnp; i++)
  {
    sprintf(tmp,"%6d",i+1);
    for (j=0; j<ndm; j++) sprintf(&(tmp[strlen(tmp)])," %12.5g", X[i*ndm+j]);
    out << tmp << "\n";
  }  

  sprintf(tmp,"\n\nelements | 123 i1 %di\n", nen);
  out << tmp;

  for (i=0; i<numel; i++)
  {
    sprintf(tmp,"%6d  ",i+1);
    for (j=0; j<elem[i]->nen(); j++) sprintf(&(tmp[strlen(tmp)]),"%6d", elem[i]->ix[j]);
    out << tmp << "\n";
  }

  sprintf(tmp,"\n\nboundary conditions\n");
  out << tmp;

  for (i=0; i<numnp; i++)
  {
    j = 0; while (j < ndm && idx(i+1,j+1) > 0) j++;
    if (j == ndm) { j = 0; while (j < ndf && idu(i+1,j+1) > 0) j++; }
    if (j != ndf)
    {
      sprintf(tmp,"%6d  ",i+1);
      for (j=0; j<ndm; j++) if (idx(i+1,j+1) < 1) sprintf(&(tmp[strlen(tmp)])," 1");
                                             else sprintf(&(tmp[strlen(tmp)])," 0");
      for (j=0; j<ndf; j++) if (idu(i+1,j+1) < 1) sprintf(&(tmp[strlen(tmp)])," 1");
                                             else sprintf(&(tmp[strlen(tmp)])," 0");
      out << tmp << "\n";
    }
  }

  out << "\n\n";

  out.close();

  return;
}








void Mesh::printComputerTime(bool reset, int detailFlg)
{
  Geometry::printComputerTime(reset,detailFlg);

  if (isALE(CHECKNEQX))
  {
    COUT << "----------------------------------------------------\n";

    COUT; printf("Mesh::updateMesh:              %7.3f sec ->%5.1f %\n",
                 ctimUpdateMesh, ctimUpdateMesh/ctimSinceLastCall*100.);
  }
  if (reset) ctimUpdateMesh = 0.;

  return;
}






