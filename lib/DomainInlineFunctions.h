
#ifndef incl_DomainInlineFunctions_h
#define incl_DomainInlineFunctions_h

#include "DomainTree.h"
#include "DomainTypeEnum.h"



#define define_reference_cast(type,Type)                         \
	                                                         \
inline Type &type(Domain &dom) { return *((Type*)&dom); }




#define define_isType(isType,TYPE)                               \
	                                                         \
inline bool isType(Domain &dom) { return (TYPE == type(dom)); }



/*
#define define_isDerivedFrom(isDerivedFromType,TYPE)                   \
                                                                       \
inline bool isDerivedFromType(Domain &dom) { return (dom.isDerivedFrom(TYPE)); }
*/



#endif




