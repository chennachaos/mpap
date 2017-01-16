
#ifndef incl_PatchGroup_h
#define incl_PatchGroup_h


#include "MyString.h"
#include "List.h"
//#include "Element.h"
#include "MathVector.h"
#include "Domain.h"
#include "PropertyItem.h"



class PatchGroup: public ListItem
{
  public:

    int ndom, reparmid, eltype; //, matl; //,nPatchProp;

    VectorArray<int> polydegree, kvindex, ctrlpoints, subDiv, elevdegree ;

  //  PropertyItem patchElemProp, patchMatlProp;


    PatchGroup();

    virtual ~PatchGroup();
	  


    //ListArray<VectorArray<double> > kvectors;
};

#endif

