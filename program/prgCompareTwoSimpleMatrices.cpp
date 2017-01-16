
#include <iostream>
#include <cmath>

#include "Definitions.h"
#include "FunctionsProgram.h"
#include "MyString.h"
#include "FunctionsEssGrp.h"


extern bool noGUI;


using namespace std;


void prgCompareTwoSimpleMatrices(double *mtx1, double *mtx2, 
		                 char *title1, 
				 char *title2, 
				 char *titleDiff, 
				 int nrow, int ncol,
				 int dig, int dig2, bool gfrmt, 
		                 int indent, 
				 bool interactive, 
				 bool numbers)
{
  int i, action, nw, row0 = 1, col0 = 1;
	
  double mx, dmx;

  MyString input, *word;

  char *quit[] = {"q","Q","x","X",NULL}, tmp[100];
  
  if (nrow*ncol == 0) 
  {
    cout << "\n"; COUT << "no matrices! nothing to be compared!\n\n";
    return;
  }

  // calculate max_coefficient(mtx1-mtx2) / max_coefficient(mtx1,mtx2)
	
  mx = abs(mtx1[0]);
  for (i=1; i<nrow*ncol; i++) if (mx < abs(mtx1[i])) mx = abs(mtx1[i]);
  for (i=0; i<nrow*ncol; i++) if (mx < abs(mtx2[i])) mx = abs(mtx2[i]);

  dmx = abs(mtx1[0]-mtx2[0]);
  for (i=1; i<nrow*ncol; i++) if (dmx < abs(mtx1[i]-mtx2[i])) dmx = abs(mtx1[i]-mtx2[i]);
 
  cout << "\n"; 
  COUT << " max_coeff(difference)\n";
  COUT << "----------------------- = " << dmx / mx << "\n";
  COUT << "    max_coeff(both)\n\n";

  if (!interactive || (interactive && nrow < 11 && ncol < 11))
  {
	  
  // display total matrices
  
    COUT << title1 << "\n\n";
    
    prgPrintSimpleMatrix(mtx1,nrow,ncol,dig,dig2,gfrmt,indent,numbers);
    
    COUT << title2 << "\n\n";
    
    prgPrintSimpleMatrix(mtx2,nrow,ncol,dig,dig2,gfrmt,indent,numbers);
    
    COUT << titleDiff << "\n\n";
    
    prgPrintSimpleMatrixDiff(mtx1,mtx2,nrow,ncol,dig,dig2,gfrmt,indent,numbers);
  }
  else
  {

  // display matrix sections

    while (1)
    {
      COUT << title1 << "\n\n";
      
      prgPrintSimpleMatrixSec(mtx1,nrow,ncol,row0,col0,10,dig,dig2,gfrmt,indent,numbers);
      
      COUT << title2 << "\n\n";
      
      prgPrintSimpleMatrixSec(mtx2,nrow,ncol,row0,col0,10,dig,dig2,gfrmt,indent,numbers);
      
      COUT << titleDiff << "\n\n";
      
      prgPrintSimpleMatrixSecDiff(mtx1,mtx2,nrow,ncol,row0,col0,10,dig,dig2,gfrmt,indent,numbers);

      action = -1;
	
      if (noGUI) 
      {
	while (action < 0)
	{
          COUT << "input row and column numbers, 'x' for exit > ";
	  nw = input.free().input().split(&word);

	  if (nw == 1)       action = word[0].which(quit);
	  else if (nw == 2)  if (word[0].toInt(&row0,false) && word[1].toInt(&col0,false)) 
		                         action = 10;
	
          for (i=0; i<nw; i++) word[i].free(); if (nw > 0) delete [] word;
	}
	cout << "\n";

        if (action < 10) break;
      }
      else
      {
	sprintf(tmp,"%d, %d",row0,col0);
	input = tmp;
	while (action < 0) 
	{	
          if (!essGrpInquire("input row and column numbers!","OK","Exit",&input)) action = 0;
	  else
	  {
	    nw = input.split(&word);
            if (nw == 2) if (word[0].toInt(&row0,false) && word[1].toInt(&col0,false)) 
                action = 10;
	    for (i=0; i<nw; i++) word[i].free(); if (nw > 0) delete [] word;
	  }
	} 
	if (action < 10) break;
      }
    }
  }
  return;
}


