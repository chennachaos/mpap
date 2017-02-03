
#include <iostream>
#include <fstream>

#include "FunctionsProgram.h"
#include "Files.h"
#include "Definitions.h"
#include "DomainTree.h"
#include "MyString.h"
#include "Debug.h"
#include "RunControl.h"
#include "Domain.h"
//#include "IsogeometricFEM.h"
#include "HBSplineFEM.h"
#include "HBSplineCutFEM.h"
#include "StandardFEM.h"


using namespace std;

extern Files      files;
extern DomainTree domain;
extern bool       keep;
extern RunControl runCtrl;


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

    //if (domain(i).ndm > plot.dim) plot.dim = domain(i).ndm;
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

