/*=============================================================================
        File: NurbsBASE.cpp
  Created by: Chennakesava Kadapa          (08 Jan 2011)
 Purpose    : Implementation file for the definitions of NURBS SURFACE Class

 ============================================================================*/

#include "NurbsBASE.h"
#include "DataBlockTemplate.h"
//#include "Plot.h"
#include <iomanip>
#include "MpapTime.h"


extern MpapTime           mpapTime;


using namespace std;

//extern Plot plot;




NurbsBASE::~NurbsBASE()
{
  //cout << "    NurbsBASE: destructor ...\n\n";

}


/*
//
//
//
//
//
void NurbsBASE::initializeBCdata()
{
   int ii, jj;

   nsize = nlbf*ndof;

   dispBCs.setDim(ngbf);
   gbfnums.setDim(ngbf);

  for(ii=0;ii<ngbf;ii++)
  {
     gbfnums[ii] = -5555;

     dispBCs[ii].setDim(ndof);
     for(jj=0;jj<ndof;jj++)
        dispBCs[ii][jj] = -7777; // some fancy number
  }

  forceBCs = dispBCs;

  Uinit.setDim(ndof*ngbf);
  Uinit.zero();
  
  Values = Uinit; // for the 1st constraint variable
  
  Values2 = Uinit;  // for the 2nd constraint variable
  
  toUfull.setDim(ndof*ngbf);
  toUfull.zero();

  dispBCdata.setDim(8);
  dispBCdata.zero();

  edgedata.setDim(4);
  edgedata.zero();
  intfdata.setDim(8);
  intfdata.zero();


  return;
}
*/







int  NurbsBASE::GenerateConnectivityArrays2()
{
    //LM Array
    
    int ii, e, ind, a;

    ind = ndof*nlbf;
    
    LM.setDim(nelem);
    for(ii=0;ii<nelem;ii++)
       LM[ii].setDim(ind);

    for(e=0;e<nelem;e++)
    {
        for(a=0;a<nlbf;a++)
        {
            ind = ndof*a;
            for(ii=0;ii<ndof;ii++)
	       LM[e][ind+ii] = ID[ii][IEN[e][a]];
        }
    }

    return 0;
}






void NurbsBASE::printConnectivityArrays()
{
    int ii, jj;

    cout << "     INC Array ... " << endl;
    for(ii=0;ii<2;ii++)
    {
       for(jj=0;jj<ngbf;jj++)
          cout << '\t' << INC[ii][jj] ;
       cout << endl;
    }
    cout << endl;
    cout << endl;    

    cout << "     IEN Array ... " << endl;
    for(ii=0;ii<nelem;ii++)
    {
        for(jj=0;jj<nlbf;jj++)
           cout <<  '\t' << IEN[ii][jj] ;
        cout << endl;
    }
    cout << endl;
    cout << endl;


    cout << "     ID Array ... " << endl;
    for(ii=0;ii<ndof;ii++)
    {
        for(jj=0;jj<ngbf;jj++)
           cout << '\t' << ID[ii][jj] ;
   	cout << endl;
    }
    cout << endl;
    cout << endl;

    cout << "     LM Array ... " << endl;
    for(ii=0;ii<nelem;ii++)
    {
        for(jj=0;jj<LM[ii].n;jj++)
           cout <<  '\t' << LM[ii][jj] ;
        cout << endl;
    }
    cout << endl;
    cout << endl;

  
return;
}




void  NurbsBASE::setSolidOrFluid(int ttt)
{
  PHYSICS_TYPE = true; 

  return;
}



void  NurbsBASE::setTimeIntegrationParameters(int ttt, double rho1)
{
  tis = ttt;
  rhoInfty = rho1;

  return;
}



void NurbsBASE::setTimeParam()
{
  td.resize(50);
  td.setZero();

  SetTimeParametersSolid(tis, rhoInfty, mpapTime.dt, td);

  return;
}




void NurbsBASE::updateIterStep()
{
  PHYSICS_TYPE = true;
  
  int ii, jj;

  if(PHYSICS_TYPE)
  {
    // displacement as the primary variable
    for(ii=0; ii<ndof; ii++)
    {
      for(jj=0; jj<ngbf; jj++)
      {
        ValuesDot[ii][jj]     = td[10]*Values[ii][jj] + td[11]*ValuesPrev[ii][jj] + td[12]*ValuesDotPrev[ii][jj] + td[13]*ValuesDotDotPrev[ii][jj] ;
        ValuesDotDot[ii][jj]  = td[15]*Values[ii][jj] + td[16]*ValuesPrev[ii][jj] + td[17]*ValuesDotPrev[ii][jj] + td[18]*ValuesDotDotPrev[ii][jj] ;

        // compute Current values

        ValuesCur[ii][jj]        =  td[2]*Values[ii][jj]       + (1.0-td[2])*ValuesPrev[ii][jj];
        ValuesDotCur[ii][jj]     =  td[2]*ValuesDot[ii][jj]    + (1.0-td[2])*ValuesDotPrev[ii][jj];
        ValuesDotDotCur[ii][jj]  =  td[1]*ValuesDotDot[ii][jj] + (1.0-td[1])*ValuesDotDotPrev[ii][jj];
      }
    }

    //printVector(var1Cur);
    //printVector(var1DotCur);
    //printVector(var1DotDotCur);

      //cout << " ggggggggggggg " << endl;
  }
  else
  {
    //cout << " uuuuuuuuuuuu  " << endl;
    //printVector(td);
    //var1Dot    = td[9]*var1 + td[10]*var1Prev + td[11]*var1Prev2 + td[12]*var1Prev3 + td[13]*var1Prev4 + td[15]*var1DotPrev ;
    //var2Dot    = td[9]*var2 + td[10]*var2Prev + td[15]*var2DotPrev ;

    //var1Cur    = td[2]*var1    + (1.0-td[2])*var1Prev; // velocity
    //var2Cur    = td[2]*var2    + (1.0-td[2])*var2Prev; // pressure
    //var3Cur    = td[2]*var3    + (1.0-td[2])*var3Prev; // Lagrange parameters
    //var4Cur    = td[2]*var4    + (1.0-td[2])*var4Prev; // Solid dof (for contact elements)

    //var1DotCur = td[1]*var1Dot + (1.0-td[1])*var1DotPrev;
    //var2DotCur = td[1]*var2Dot + (1.0-td[1])*var2DotPrev;
    //cout << " uuuuuuuuuuuu  " << endl;
  }

  return;
}



