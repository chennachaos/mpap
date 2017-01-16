#include <iostream>
#include <math.h>
#include "NurbsMiscFunctions.h"

using namespace std;



void checkContSurfaces(NurbsSURFACE* surf1, NurbsSURFACE* surf2, int dir)
{
   int d, ii, count=0;
   ListArray<ListArray<EPOINT> > SKL1, SKL2;

   d = min(surf1->p, surf2->p);
   surf1->SurfDerPointRat(1.0, 1.0, d, SKL1);
   surf2->SurfDerPointRat(0.0, 1.0, d, SKL2);

   cout.setf(ios::fixed);
   cout.precision(8);
   for(ii=0;ii<=d;ii++)
   {
      SKL1[ii][0].print2screen();
      SKL2[ii][0].print2screen();
      cout << fixed << '\t' << (SKL2[ii][0].x - SKL1[ii][0].x)*1000 << '\t' << (SKL2[ii][0].y - SKL1[ii][0].y)*1000 << endl;
      cout << endl;
      if(SKL1[ii][0] == SKL2[ii][0])
        count = ii;
      else
        break;
   }

   cout << '\t' << " surf1->p : " << surf1->p << endl;
   cout << '\t' << " surf2->p : " << surf2->p << endl;
   cout << endl;
   cout << '\t' << " The two surfaces are 'C" << count << "' continuous in " << " direction : " << dir << endl;
   cout << endl;

  return;
}
