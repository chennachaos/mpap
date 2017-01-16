
#ifndef incl_FreeSurface_h
#define incl_FreeSurface_h

#include "Fluid.h"



class FreeSurface: public Fluid
{ 
  public:

    FreeSurface(void);

    virtual ~FreeSurface();

    virtual void readInputData(std::ifstream &, MyString &);

    virtual void prepareInputData(void);

    virtual void prepareInteractions(void);

    virtual void printInfo(void);

    virtual void setSolver(int, int *parm = NULL, bool cIO = false);

    virtual void prepareForExternalSolver(void*,double*,bool);

    virtual void setLagrangianFreeNodes(VectorArray<int> &);

  private:

};



#include "DomainInlineFunctions.h"

define_reference_cast(freeSurface,FreeSurface)

define_isType(isFreeSurface,FREESURFACE)




	
#endif


