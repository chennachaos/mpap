
#include <iostream>

#include "FunctionsProgram.h"
#include "DomainTypeEnum.h"
#include "Domain.h"
#include "DomainTree.h"
#include "DomainType.h"
#include "Debug.h"
#include "ComputerTime.h"
#include "MyString.h"



extern DomainTree domain;
extern ComputerTime computerTime;


using namespace std;


Domain::Domain(void)                       
{                                                  
  if (debug) std::cout << " Domain constructor\n\n";

  //std::cout << ++nObj << "\n\n";

  s  = NULL;
  p  = NULL;
  xl = NULL;
  ul = NULL;
  
  solverOK = false;
  
  tis = -1;

  firstIter = false;

  ndm = 1;

  // add new type
  
  DomainType *dom = domain.newType(ROOTDOMAIN,-1);

  if (dom == NULL) return;  // domain type exists already
  
  dom->key.addNew("dimensions");
 
  return;
}




	                                          
Domain::~Domain(void)                     
{
  if (s  != NULL) delete [] s;
  if (p  != NULL) delete [] p;
  if (xl != NULL) delete [] xl;
  if (ul != NULL) delete [] ul;
 
  if (debug) std::cout << " Domain destructor\n\n";

  //std::cout << --nObj << "\n\n";

  //std::cout << " destructor -- Domain " << std::endl;
}




void Domain::readFile(std::ifstream &Ifile)
{
  MyString keyStrg;

  while (prgNextDataKey(Ifile, keyStrg)) 
    readInputData(Ifile, keyStrg);

  prepareInputData();
  
  return;
}




void Domain::readInputData(std::ifstream &Ifile, MyString &line)
{
  int i, nw;
	
  MyString *word;
 
  char fct[] = "Domain::readInputData";
 
  //cout << fct << "\n";
 
  switch (domain[ROOTDOMAIN].key.whichBegins(line))
  {
    case  0: cout << "     DOMAIN: reading dimensions ...\n\n";
	   
	     line.getNextLine(Ifile);
	     
	     nw = line.split(&word);
            
	     if (nw < 2)                prgError(1,fct,"input error in 'dimensions'!");
		     
             if (!word[0].toInt(&ndf))  prgError(2,fct,"input error in 'dimensions'!");
		     
             if (!word[1].toInt(&ndm))  prgError(3,fct,"input error in 'dimensions'!");

             for (i=0; i<nw; i++) word[i].free(); delete [] word;

             if (ndm < 0 || ndf < 0 || ndm > 3) prgError(4,fct,"invalid values in 'dimensions'!");
	     
       	     line.getNextLine(Ifile);

	     break;

    case -1: std::cout << " Error in input file: \"" << line << "\"\n\n";
	
             prgError(1,"Domain::readInputData","unknown keyword!");	
  }

  return;
}




void Domain::prepareInputData(void)
{
  // nothing to be done here!
  
  //cout << " NOTHING DONE ........................................................ " << endl;
	
  return;
}




void Domain::prepareInteractions(void)
{
  // nothing to be done here!

  return;
}





void Domain::startComputerTime(void)
{
  MyString tmp = domain.name(this);

  computerTime.go(tmp.asCharArray()); 

  return;
}





void Domain::printComputerTime(bool reset, int detailFlg)
{
  MyString tmp = domain.name(this);

  ctimSinceLastCall = computerTime.stop(tmp.asCharArray());

  COUT << domain.name(this) << ": computer time statistics\n";
  COUT << "----------------------------------------------------\n";
  COUT; printf("time since last call:        %9.3f sec\n",ctimSinceLastCall);

  computerTime.go(tmp.asCharArray()); 

  return;
}

