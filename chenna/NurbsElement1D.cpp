
#include "Debug.h"
#include "FunctionsProgram.h"
#include "PropertyTypeEnum.h"
#include "NurbsElement1D.h"
#include <iostream>
#include <iomanip>
#include "Plot.h"
#include "MpapTime.h"


extern Plot plot;
extern MpapTime           mpapTime;

using namespace std;



NurbsElement1D::NurbsElement1D(void)
{
  if (debug) cout << " constructor NurbsElement1D\n\n";

  // cout << "     NurbsElement1D: constructor ...\n\n";

  curve0 	= NULL;
  curve1 	= NULL;

}





NurbsElement1D::~NurbsElement1D()
{
  if (debug)   cout << " destructor NurbsElement1D\n\n";

  if(curve0 != NULL)  // dont delete this pointer as it is just a reference
     curve0 = NULL;

  if(curve1 != NULL)  // dont delete this pointer as it is just a reference
     curve1 = NULL;

}




void NurbsElement1D::initialiseDOFvalues()
{
    int index1, index2;
    for(int ii=0;ii<nlbf;ii++)
    {
       index1 = ndof * ii;
       index2 = ndof * (curve0->IEN[elenum][ii]);
       for(int jj=0;jj<ndof;jj++)
          primvar[index1+jj] = curve0->Uinit[index2+jj];
    }


  return;
}






