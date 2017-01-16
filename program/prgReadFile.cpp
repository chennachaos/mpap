
#include <iostream>
#include <fstream>

#include "FunctionsProgram.h"
#include "Files.h"
#include "Definitions.h"
#include "DomainTree.h"
#include "MyString.h"
#include "Debug.h"
#include "RunControl.h"
#include "Plot.h"
// include all domain types
#include "Domain.h"
#include "Mesh.h"
#include "FiniteElementBVP.h"
#include "Fluid.h"
#include "Solid.h"
#include "MicroCellWulf.h"
#include "MicroCell.h"
#include "Fem1D.h"
#include "BeamSection.h"
#include "BeamBending.h"
#include "RigidBody.h"
#include "Aeroplane.h"
#include "InterfaceN.h"
#include "FreeSurface.h"
#include "InterfaceGS.h"
//#include "InterfacePB.h"
#include "VortexSheet.h"
#include "LiftingLine.h"
#include "FlexibleWing.h"
#include "IsogeometricFEM.h"
#include "HBSplineFEM.h"
#include "HBSplineCutFEM.h"
#include "StandardFEM.h"
//#include "HBScutFEMElasticity.h"

using namespace std;


extern Files      files;
extern DomainTree domain;
extern bool       keep;
extern RunControl runCtrl;
extern Plot       plot;


bool prgReadFile(void)
{
  int      key, n, i;
  
  bool     tmpIfileFlag = false, runControlFlag = false, timeFunctionsFlag = false;
 
  ifstream *Ifile = new ifstream;	
 
  ofstream tmpIfile;
  
  MyString line, pathAndFile, tmpPathAndFile;

  char *domType[] = DOMAIN_TYPE_NAMES, **words;   int nWords;
  
  
  // open Ifile ..................................................................................

  pathAndFile.append(files.projDir).append(SLASH).append(files.Ifile);

  if (debug) cout << " " << pathAndFile << "\n\n";
  
  Ifile->open(pathAndFile.asCharArray());

  if (!*Ifile) prgError(1,"prgReadFile","failed to open the input file.");

  runCtrl.newStatus(LOADING);
		
  
  // check for MPAP2 identifier ..................................................................
  
  line.read(*Ifile);

  if (debug) cout << " >" << line << "<\n\n";
  if (debug) cout << " >" << line.stripToMin() << "<\n\n";
  
  if (line.stripToMin() != "MPAP2")  prgError(2,"prgReadFile","'MPAP2' identifier not found.");

  
  // if any '#include' then generate (temporary) single input file ...............................

  while (*Ifile && !line.begins("#include") && !line.begins("#INCLUDE")) 
     line.read(*Ifile).stripToMin();

  if (*Ifile) 
    { 
       if (debug) cout << " '#INCLUDE' found in Ifile. Generate Ifile.tmp ....\n\n";
       tmpIfileFlag = true;
       Ifile->close();  // close Ifile (it will be opened and read by 'prgResolveIncludeFile' )
       tmpPathAndFile.append(pathAndFile).append(".tmp");
       tmpIfile.open(tmpPathAndFile.asCharArray());  // open Ifile.tmp as output stream
       
       prgResolveIncludeFile(tmpIfile,pathAndFile.asCharArray());  // write Ifile.tmp

       tmpIfile.close();                             // close output stream Ifile.tmp
       Ifile->open(tmpPathAndFile.asCharArray());     // open Ifile.tmp as input stream
    } 
  else 
    {  if (debug) cout << " no '#INCLUDE' found in Ifile. Proceed with original Ifile ...\n\n";
       Ifile->close();                           //
       delete Ifile; Ifile = new ifstream;  // move to beginning of Ifile stream   
       Ifile->open(pathAndFile.asCharArray()); } // 
       
  // read domains and what_to_do part ............................................................

  do {
     key = prgNextDomainKey(*Ifile);
     
     if (debug) cout << " key = " << key << "\n\n";
  
     switch (key)
       {
	  case -2   : // read TIME_FUNCTIONS part
		  
		 if (timeFunctionsFlag) 
	           { cout << "   WARNING! multiple definition of TIME_FUNCTIONS!\n\n"; break; }
		 
                 cout << "   loading TIME_FUNCTIONS ...\n\n";
		 
                 prgReadTimeFunctions(*Ifile);
		 
		 timeFunctionsFlag = true; 
		 
		 break;
		 
	  case -1   : // read RUN_CONTROL part
		  
		 if (runControlFlag) 
	           { cout << "   WARNING! multiple definition of RUN_CONTROL part!\n\n"; 
	             break; }
		 
                 cout << "   loading RUN_CONTROL ...\n\n";
		 
                 prgReadRunControl(*Ifile);
		 
		 runControlFlag = true; 
		 
		 break;
		 
          case MESH: domain.newDom(new Mesh);  
		      
		 n = domain[MESH].dom.n;
		      
                 cout << "   loading MESH " << n << " ...\n\n"; 
		      
                 domain(MESH,n-1).readFile(*Ifile); 
		 
		 break;
    
          case FINITEELEMENTBVP: domain.newDom(new FiniteElementBVP);  
		      
		 n = domain[FINITEELEMENTBVP].dom.n;
		      
                 cout << "   loading FINITEELEMENTBVP " << n << " ...\n\n"; 
		      
                 domain(FINITEELEMENTBVP,n-1).readFile(*Ifile); 
		 
		 break;
    
          case FLUID: domain.newDom(new Fluid);  

		 n = domain[FLUID].dom.n;
		      
                 cout << "   loading FLUID " << n << " ...\n\n"; 
		      
                 domain(FLUID,n-1).readFile(*Ifile); 
		 
		 break;
    
          case SOLID: domain.newDom(new Solid);
	    
		 n = domain[SOLID].dom.n;
		      
                 cout << "   loading SOLID " << n << " ...\n\n"; 

		 domain(SOLID,n-1).readFile(*Ifile);
    
                 break;
	    
          case MICROCELLWULF: domain.newDom(new MicroCellWulf);
			  
		 n = domain[MICROCELLWULF].dom.n;
	    
                 cout << "   loading MICROCELLWULF " << n << " ...\n\n"; 
	    
		 domain(MICROCELLWULF,n-1).readFile(*Ifile);

                 break;
	    
          case MICROCELL: domain.newDom(new MicroCell);
			  
		 n = domain[MICROCELL].dom.n;
	    
                 cout << "   loading MICROCELL " << n << " ...\n\n"; 
	    
		 domain(MICROCELL,n-1).readFile(*Ifile);

                 break;
	    
          case FEM1D: domain.newDom(new Fem1D);

		 n = domain[FEM1D].dom.n;
	    
                 cout << "   loading FEM1D " << n << " ...\n\n"; 
	    
		 domain(FEM1D,n-1).readFile(*Ifile);

                 break;
	    
          case BEAMSECTION: domain.newDom(new BeamSection);
		
		 n = domain[BEAMSECTION].dom.n;
	    
                 cout << "   loading BEAMSECTION " << n << " ...\n\n"; 
	    
		 domain(BEAMSECTION,n-1).readFile(*Ifile);

                 break;
	    
          case BEAMBENDING: domain.newDom(new BeamBending);
		
		 n = domain[BEAMBENDING].dom.n;
	    
                 cout << "   loading BEAMBENDING " << n << " ...\n\n"; 
	    
		 domain(BEAMBENDING,n-1).readFile(*Ifile);

                 break;
	    
          case RIGIDBODY: domain.newDom(new RigidBody);
		
		 n = domain[RIGIDBODY].dom.n;
	    
                 cout << "   loading RIGIDBODY " << n << " ...\n\n"; 
	    
		 domain(RIGIDBODY,n-1).readFile(*Ifile);

                 break;
	    
          case AEROPLANE: domain.newDom(new Aeroplane);
		
		 n = domain[AEROPLANE].dom.n;
	    
                 cout << "   loading AEROPLANE " << n << " ...\n\n"; 
	    
		 domain(AEROPLANE,n-1).readFile(*Ifile);

                 break;
	    
          case INTERFACEN: domain.newDom(new InterfaceN);
	    
		 n = domain[INTERFACEN].dom.n;
			  
                 cout << "   loading INTERFACEN " << n << " ...\n\n"; 
	    
		 domain(INTERFACEN,n-1).readFile(*Ifile);
		 
                 break;
	    
          case FREESURFACE: domain.newDom(new FreeSurface);
	    
		 n = domain[FREESURFACE].dom.n;
			  
                 cout << "   loading FREESURFACE " << n << " ...\n\n"; 
	    
		 domain(FREESURFACE,n-1).readFile(*Ifile);
		 
                 break;
	    
        case INTERFACEGS: domain.newDom(new InterfaceGS);
	    
		 n = domain[INTERFACEGS].dom.n;
			  
                 cout << "   loading INTERFACEGS " << n << " ...\n\n"; 
	    
		 domain(INTERFACEGS,n-1).readFile(*Ifile);
		 
                 break;
	    
/*          case INTERFACEPB: domain.newDom(new InterfacePB);
	    
		 n = domain[INTERFACEPB].dom.n;
			  
                 cout << "   loading INTERFACEPB " << n << " ...\n\n"; 
	    
		 domain(INTERFACEPB,n-1).readFile(*Ifile);
		 
                 break;
	    */
          case VORTEXSHEET: domain.newDom(new VortexSheet);
	    
		 n = domain[VORTEXSHEET].dom.n;
			  
                 cout << "   loading VORTEXSHEET " << n << " ...\n\n"; 
	    
		 domain(VORTEXSHEET,n-1).readFile(*Ifile);
		 
                 break;
	    
          case LIFTINGLINE: domain.newDom(new LiftingLine);
	    
		 n = domain[LIFTINGLINE].dom.n;
			  
                 cout << "   loading LIFTINGLINE " << n << " ...\n\n"; 
	    
		 domain(LIFTINGLINE,n-1).readFile(*Ifile);
		 
                 break;
	    
          case FLEXIBLEWING: domain.newDom(new FlexibleWing);
	    
		 n = domain[FLEXIBLEWING].dom.n;
			  
                 cout << "   loading FLEXIBLEWING " << n << " ...\n\n"; 
	    
		 domain(FLEXIBLEWING,n-1).readFile(*Ifile);
		 
                 break;
	    
          case ISOGEOMETRICFEM:

                 domain.newDom(new IsogeometricFEM);

                n = domain[ISOGEOMETRICFEM].dom.n;

                cout << "   loading ISOGEOMETRICFEM " << n << " ...\n\n"; 

                domain(ISOGEOMETRICFEM,n-1).readFile(*Ifile);

                break;

          case  HBSPLINEFEM:

                domain.newDom(new HBSplineFEM);

                //cout << " AAAAAAAAAAAAAAA  " << HBSPLINEFEM << endl;

                n = domain[HBSPLINEFEM].dom.n;

                cout << "   loading HBSPLINEFEM " << n << " ...\n\n"; 

                domain(HBSPLINEFEM,n-1).readFile(*Ifile);

                break;

          case  HBSPLINECUTFEM:

                domain.newDom(new HBSplineCutFEM);

                //cout << " AAAAAAAAAAAAAAA  " << HBSPLINEFEM << endl;

                n = domain[HBSPLINECUTFEM].dom.n;

                cout << "   loading HBSPLINECUTFEM " << n << " ...\n\n"; 

                domain(HBSPLINECUTFEM,n-1).readFile(*Ifile);

                break;

          case  STANDARDFEM:

                domain.newDom(new StandardFEM);

                //cout << " AAAAAAAAAAAAAAA  " << StandardFEM << endl;

                n = domain[STANDARDFEM].dom.n;

                cout << "   loading STANDARDFEM " << n << " ...\n\n"; 

                domain(STANDARDFEM,n-1).readFile(*Ifile);

                break;

          /*

          case  HBSCUTFEMELASTICITY:

                domain.newDom(new HBScutFEMElasticity);

                //cout << " AAAAAAAAAAAAAAA  " << HBSPLINEFEM << endl;

                n = domain[HBSCUTFEMELASTICITY].dom.n;

                cout << "   loading HBSCUTFEMELASTICITY " << n << " ...\n\n"; 

                domain(HBSCUTFEMELASTICITY,n-1).readFile(*Ifile);

                break;
          */

	  // case OTHERS: .....

      }
  
  } while (key > -3 && Ifile); 

  
  // close file ...................................................................................

  Ifile->close(); delete Ifile;
  
  if (tmpIfileFlag && !keep) 
    {  
       tmpPathAndFile.insert(0," ").insert(0,DELETE_FILE_COMMAND);
       if (system(tmpPathAndFile.asCharArray()) != 0)
	 cout << "   WARNING! I could not remove the temporary input File!\n\n"; 
    }
  
  if (!runControlFlag) prgError(1,"prgReadFile","You must define the RUN_CONTROL section!");

  
  // prepare eventual domain interactions ........................................................

  for (i=0; i<domain.ndom; i++)  
  {
    cout << "   preparing interactions for " << domain.key(&(domain(i))) << " ...\n\n"; 

    domain(i).prepareInteractions();

    if (domain(i).ndm > plot.dim) plot.dim = domain(i).ndm;
  }

  for (i=0; i<domain.ndom; i++) domain(i).startComputerTime(); 

  // print info and finish ........................................................................
  
  if (domain.ndom > 0)
  {
    cout << "\n   inheritance of domain types loaded:\n\n";  domain.print(); 
  
    cout << "\n   number of domains:\n\n";
    for (int j, i=0; i<prgListCharArrayLength(domType);i++)
      {  j = domain.nDomainOfType(i);
         if (j > 0) printf("    %-11s : %i\n",domType[i],j); } cout << "\n\n";
  }
  cout << "   The file " << files.Ifile << " has been read successfully.\n\n";
 
  // create empty Tfile (overwriting existing one!)

  pathAndFile.free().append(files.projDir).append(SLASH).append(files.Tfile);

  tmpIfile.open(pathAndFile.asCharArray(),ios_base::out);
  tmpIfile.close();
 
  return true; 
}

