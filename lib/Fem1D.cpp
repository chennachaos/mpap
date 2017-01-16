
#include <iostream>

#include "FunctionsProgram.h"
#include "PropertyTypeEnum.h"
#include "DomainTypeEnum.h"
#include "Fem1D.h"
#include "DomainTree.h"
#include "DomainType.h"
#include "Definitions.h"
#include "List.h"
#include "TimeFunction.h"
#include "MpapTime.h"


extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern MpapTime           mpapTime;


using namespace std;



Fem1D::Fem1D(void)                       
{                                                  
  // add new type
 
  DomainType *fem1D = domain.newType(FEM1D,FINITEELEMENTBVP);
 
  if (fem1D == NULL) return;  // domain type already exists

  fem1D->key.addNew("control");

  return;
}




	                                          
Fem1D::~Fem1D(void)                     
{                                                 
  return;
}









void Fem1D::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString *word;
 
  char fct[] = "Fem1D::readInputData";

  int nw, i;
 
  switch (domain[FEM1D].key.whichBegins(line))
  {
    case  0: cout << "     FEM1D: reading control ...\n\n";

	     line.getNextLine(Ifile);
	     
	     nw = line.split(&word);
            
	     if (nw < 1)                    prgError(1,fct,"input error in 'control'!");
		     
             if (!word[0].toDbl(&tol))      prgError(2,fct,"input error in 'control'!");
		     
             for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     
       	     line.getNextLine(Ifile);

	     break;
	     
    case -1: // go and inherit from FINITEELEMENTBVP
	     
	     this->FiniteElementBVP::readInputData(Ifile, line); 
	     
	     break;
  }
 
  return;
}





void Fem1D::prepareInputData(void)
{
  int e;
      
  // go and inherit from ancestors

  FiniteElementBVP::prepareInputData();
 
  
  std::cout << "     FEM1D: prepare input data ...\n\n";

  

  return;
}







