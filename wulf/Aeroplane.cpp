
#include <iostream>

#include "DomainTypeEnum.h"
#include "Aeroplane.h"
#include "DomainTree.h"
#include "DomainType.h"
#include "FunctionsProgram.h"
#include "FunctionsEssGrp.h"
#include "Plot.h"
#include "MpapTime.h"


extern DomainTree domain;
extern Plot plot;
extern MpapTime mpapTime;



using namespace std;



Aeroplane::Aeroplane(void)                       
{                                                  
  // add new type
  
  DomainType *aero = domain.newType(AEROPLANE,ROOTDOMAIN);

  if (aero == NULL) return;  // domain type exists already

  aero->key.addNew("specification","screen layout");

  aDAW = NULL;

  flightProfile = NULL;

  return;
}




	                                          
Aeroplane::~Aeroplane(void)                     
{         
  if (aDAW != NULL) delete [] aDAW;

  if (flightProfile != NULL) delete flightProfile;

  return;
}





void Aeroplane::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString tmpl, *word;
 
  char tmp[30], fct[] = "Aeroplane::readInputData",
       *aeroDatType[] = {"Ixx", "Izz", "Iyy", "Ixz", "CD0", "S", "b", NULL};

  int nw, i;

  if (ndm != 1) prgError(1,fct,"ndm != 1");
 
  switch (domain[AEROPLANE].key.whichBegins(line))
  {
    case  0: cout << "     AEROPLANE: reading specification ...\n\n";
	   
             line.getNextLine(Ifile);
 
	     break;

    case  1: cout << "     AEROPLANE: reading screen layout ...\n\n";
	   
             line.getNextLine(Ifile);
 
	     break;

    case -1: // go and inherit from DOMAIN
	     
	     this->Domain::readInputData(Ifile, line); 
	     
	     break;
  }
 
  return;
}








void Aeroplane::prepareInputData(void)
{
  // call ancestor function

  Domain::prepareInputData();

  
  cout << "     AEROPLANE: prepare input data ...\n\n";
  
  char fct[] = "Aeroplane::prepareInputData"; 
	  
  // ........

  return;
}








void Aeroplane::prepareInteractions(void)
{
  // go and inherit from ancestors

  Domain::prepareInteractions();

  cout << "     AEROPLANE: preparing interactions ...\n\n"; 
      
  return;
}








void Aeroplane::initFlightSimulatorDisplay(bool*)
{
  char fct[] = "Aeroplane::initFlightSimulatorDisplay";

  essGrpGetPlotAreaGeom();

  int w = plot.wPix, h = plot.hPix;

  cout << h << "," << w << "\n\n";

  if (h < 300 || w < 400) 
  {
    COUT << "The graphical display is too small to display the flight simulator!\n\n";

    return;
  }

  aDAW = new DAWindow [4];


  aDAW[0].setup(w-202,2,200,60,1.,1.,true);

  aDAW[0].wipe(1);

  aDAW[0].frame(0,3);


  aDAW[1].setup(w-202,104,200,60,1.,1.,true);

  aDAW[1].wipe(2);

  aDAW[1].frame(3,3);


  aDAW[2].setup(w-202,206,200,60,1.,1.,true);

  aDAW[2].wipe(4);

  aDAW[2].frame(5,3);


  flightProfile = new DAFunctionPlot;

  flightProfile->setup(100,100,500,500,1,2,3,5);

 
  return;
}







int Aeroplane::factoriseSolveAndUpdate(void)
{
  flightProfile->addXY(mpapTime.cur,1.*sin(mpapTime.cur));

  flightProfile->draw();

  return 0;
}


