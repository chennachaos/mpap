
#ifndef incl_DomainType_h
#define incl_DomainType_h

#include "Tree.h"
#include "Domain.h"
#include "List.h"
#include "MyString.h"
#include "MyStringList.h"


class DomainType: public TreeItem
{
    public: 

       DomainType(int);

       virtual ~DomainType(void);
       
       int  tid;

       MyString type;
       
       MyStringList key, var, fct;

       List<Domain> dom;
       
       bool operator==(const DomainType &);
       
       virtual char *printName(void);
};

#endif


