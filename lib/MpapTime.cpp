
#include "MpapTime.h"
#include "FunctionsProgram.h"

// ---- for macro163 ----------
#include <iostream>
#include "MyString.h"
#include "Definitions.h"
#include "DataBlockTemplate.h"
using namespace std;
// ----------------------------




void MpapTime::cut(void)
{ 
  dt *= 0.5;
  
  stack.append(dt);

  return; 
}



void MpapTime::stck(void)
{ 
  int m = stack.n;

  if (m == 1) { dt = stack[0]; return; }
 
  if (m > 1)  { dt = stack[m-1]; stack.trunc(m-1); return; }
  
  if (m < 1)  prgError(1,"MpapTime::stck","stack totally empty!");
	
  return; 
}



void MpapTime::reset(void)
{ 
  dt = -1.; 
  cur = 0.; 
  prev = -1.; 
  prev2 = -2.; 
  write = -1.;
  dtOK = false; 
  stack.free(); 
  
  return; 
}





// ---- for Deniz's macro163 -------------

void MpapTime::file(char* fname)
{ 
char fct[] = "MpapTime::file";
int i=0,j,nw;


if(Tstep.n==0)	{

  // open file 

  ifstream file;

  double dtmp;

  MyString line, *word = NULL;
  
  file.open(fname);

  if (!file) { prgError(1,fct,"failed to open file!"); }


	// read time step map
  
  while (true)
  {
    if (!line.getNextLine(file)) 
      { prgError(1,fct,"failed to read time steps!"); }
    
    nw=line.split(&word);

	 if (nw > 1) 
      { prgError(2,fct,"one step per row!"); }

	 	 if (nw == 0) {break; }

	 if (!word[0].toDbl(&dtmp)) {break; }
	 else {Tstep[i]=dtmp;  i++; }
	  

    for (j=0; j<nw; j++) word[j].free(); delete [] word; word = NULL;
  }

    nTstep=-1;
	}
  nTstep++;  	 if (nTstep>=Tstep.n) { prgWarning(3,fct,"no more time steps in file!"); nTstep--;}
  dt = Tstep[nTstep];
  dtMax=dt;
   
  
  return; 
}
