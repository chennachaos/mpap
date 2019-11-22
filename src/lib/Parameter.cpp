
#include <iostream>

#include "Parameter.h"
#include "Debug.h"
#include "FunctionsProgram.h"
#include "MacroQueue.h"
#include "Loop.h"
#include "List.h"
#include "DomainTree.h"
#include "Counter.h"


extern MacroQueue            macroQueue;
extern DomainTree            domain;
extern ListInfinite<Counter> counter;

using namespace std;




Parameter::Parameter(void)
{
  x = 0.;
  fctName.free();
  fctPar.free();

  if (debug) std::cout << " Parameter constructor\n\n";
  
  return;
}





Parameter::~Parameter(void)
{
  fctName.free();
  fctPar.free();

  if (debug) std::cout << " Parameter destructor\n\n";
  
  return;
}





Parameter &Parameter::operator=(Parameter &par)
{
  x       = par.x;
  fctName = par.fctName;
  fctPar  = par.fctPar;
	
  return *this;
}






void Parameter::free(void)
{
  x = 0.;
  fctName.free();
  fctPar.free(); 
  return;
}






bool Parameter::interprete(char *strg, NumberType numberType)
{
  int i, j, fct, l = strlen(strg);

  char *functions[] = { "type", NULL},
       *domains[]   = DOMAIN_TYPE_NAMES;
  
  fctName.free();
  fctPar.free();
  
  if (!scanDbl(strg,&x,false)) 
  {
    x = 0.;

    i = 1; while (strg[i]!='(' && i<l) i++; 
    
    j = i; while (strg[j]!=')' && j<l) j++;

    if (strg[0]!='$' || i == 1 || (j>i && strg[j]!=')') || j<l-1) return false;
    
    fctName.append(&(strg[1])).trunc(i-1);
    
    if (i != j) fctPar.append(&(strg[i+1])).trunc(j-i-1);

    // standard functions which can be evaluated straightaway
  
    fct = fctName.which(functions);

    switch (fct)
    {
      case 0: // type
	   
              fct = fctPar.which(domains);

	      if (fct < 0) return false;
	    
	      x = (double) fct;

              fctName.free();
              fctPar.free();
 
              break;	      
    }  
 
    return true;
  }

  switch (numberType)
  {
    case ALL_REAL    : return true;
		       break;
		       
    case NON_NEG_REAL: if (x > -1.e-15) return true;
		       break;
		       
    case ALL_INT     : i = roundToInt(x);
		       if (fabs((double) i - x) < 1.e-3) { x = (double) i; return true; }
		       break;
		       
    case NON_NEG_INT : i = roundToInt(x);
		       if (fabs((double) i - x) < 1.e-3 && i > -1) { x = (double) i; return true; }
		       break;

    default          : prgError(1,"Parameter::interprete","unknown number type!");
  }
 
  return false;  
}






double Parameter::evaluate(int act)
{
  int fct;
	
  bool outp = debug;
  
  List<Loop> *loop = &macroQueue.loop;
  
  char *functions[] = { "colour", "ntype", NULL}, 
       *colours[]   = { COLOUR_NAMES, NULL },
       *domains[]   = DOMAIN_TYPE_NAMES;
	    
  
  // constant paramter
	
  if (!fctName) return x;

  if (outp) cout << "  " << fctName << "(" << fctPar << ")  ";


  // standard functions
  
  fct = fctName.which(functions);

  switch (fct)
  {
    case 0: // colour

	    fct = fctPar.which(colours);
	    
	    if (fct < 0) goto unknownPar;
		    
            if (outp) cout << " evaluate -> " << fct+1 << "\n\n";
	    
	    return (double) fct+1;
	    
    case 1: // ntype
	    
            fct = fctPar.which(domains);

	    if (fct < 0) goto unknownPar;
	    
            if (outp) cout << " evaluate -> " << domain.nDomainOfType(fct) << "\n\n";
	    
	    return (double) domain.nDomainOfType(fct);
  }  

  
  // loop variable

  fct = 0; while (fct < loop->n && fctName != (*loop)[fct].name) fct++;
  
  if (fct < loop->n)
  { 
    if ((*loop)[fct].beg > act || (*loop)[fct].end < act) goto outOfLoop;
  
    if (outp) cout << " evaluate -> " << (*loop)[fct].cnt << "\n\n";
	    
    return (double) (*loop)[fct].cnt;
  }


  // counter variable

  fct = 0; while (fct < counter.n && fctName != counter[fct].name) fct++;

  if (fct < counter.n)
  {
    if (outp) cout << " evaluate -> " << (*loop)[fct].cnt << "\n\n";

    return (double) counter[fct].val;
  }
  

  // failures
  
  // unknownFct:
  
    prgWarning(1,"Parameter::evaluate","unknown function name! -> using  p = 0");
  
    return 0.;
    
  unknownPar:

    prgWarning(1,"Parameter::evaluate","unknown parameter name! -> using  p = 0");
    
    return 0.;

  outOfLoop:

    prgWarning(1,"Parameter::evaluate","reference to inactive loop counter! -> using  p = 0");
    
    return 0.;
}






