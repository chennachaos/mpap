
#include <iostream>
#include <ctime>


#include "ComputerTime.h"
#include "MathBasic.h"
#include "FunctionsProgram.h"
#include "Definitions.h"



ComputerTime::ComputerTime(void)
{
  clocksPerSec = (double) CLOCKS_PER_SEC;

  secPerClock  = 1. / clocksPerSec;

  reset();
  
  return;
}



ComputerTime::~ComputerTime()
{
  free();

  return;
}


    

void ComputerTime::free(void)
{
  name.free();

  t0.free();

  return;
}





void ComputerTime::reset(void)
{
  free();

  return;
}





void ComputerTime::sleep(double ms)
{
  unsigned int clck = (unsigned int) roundToInt(clocksPerSec * ms / 1000.);
	
  clock_t goal = clck + clock();
  
  while (goal > clock());
  
  return;
}


    

bool ComputerTime::go(char *strg, bool warn)
{
  //std::cout << strg << "\n";

  //std::cout << name << "\n";

  if (name.which(strg) > -1) 
  { 
    if (warn) prgWarning(1,"ComputerTime::go","time is already being taken!");
    return false;
  }

  name.addNew(strg);

  t0.append((unsigned int) clock());
      
  return true;
}




double ComputerTime::stop(char *strg)
{
  int i = name.which(strg); 
  
  if (i < 0)  { prgWarning(1,"ComputerTime::stop","time identifier unknown!");  return 0; }

  unsigned int dt = (unsigned int) clock() - t0[i];
 
  name.del(&(name[i]));  

  t0.del(i);
  
  return dt * secPerClock;
}




double ComputerTime::stopAndPrint(char *strg)
{
  int i = name.which(strg); 
  
  if (i < 0)  { prgWarning(1,"ComputerTime::stopAndPrint","time identifier unknown!");  return 0; }

  double dt = (double)((unsigned int) clock() - t0[i]) * secPerClock;

  name.del(&(name[i]));  

  t0.del(i);
 
  COUT << "computer time for '" << strg << "': " ;
  printf("%12.8E\tsec\n\n", dt);
  
  return dt;
}




