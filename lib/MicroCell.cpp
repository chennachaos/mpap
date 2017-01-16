#include "FunctionsProgram.h"
#include "MicroCell.h"
#include "TimeFunction.h"
//for copy purpose
#include "SolverMA41.h"
#include "Solver.h"
#include "MyString.h"
#include "PropertyTypeEnum.h"
#include "MathGeom.h"

extern List<TimeFunction> timeFunction;


MicroCell::MicroCell(void)                       
{                                                  

  // add new type
  
  DomainType *microCell = domain.newType(MICROCELL,SOLID);
  
  if (microCell == NULL) return;  // domain type already exists

  microCell->key.addNew(
		     "boundary conditions", 
		     "prescribed displacements",
		     "dependent displacements",
		     "identical displacements", 
	             "boundary type",
                     "orientation",
                     "clones",
                     "other stuff");
  return;
}




MicroCell::~MicroCell(void)                     
{                                                 
}
 
void MicroCell::readInputData(std::ifstream &Ifile, MyString &line)
{
  return;
}
 
void MicroCell::prepareInputData(void)
{
  return;
}
 
void MicroCell::prepareInteractions(void)
{
  return;
}

