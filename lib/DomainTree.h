
#ifndef incl_DomainTree_h
#define incl_DomainTree_h


#include "Tree.h"
#include "Domain.h"
#include "DomainType.h"
#include "MyString.h"


class DomainTree: public Tree<DomainType>
 {
    public:

      DomainTree(void);
 
      virtual ~DomainTree(void);
	    
      int mxtid;                         // maximum id of domain types loaded

      int ndom;                          // total number of domains loaded

      void reset(void);

      int nDomainOfType(int);            // number of domains of type t loaded

      int nType();                       // number of types loaded with at least one domain
      
      DomainType *newType(int, int);     // generate new domainType t as child of t0, if it does 
                                         // not exist already; return pointer

      void delType(DomainType*);         // del domainType and all its children and all 
                                         // associated domains
     
     
      DomainType &operator[](int);       // reference to domainType t

      Domain     &operator()(int);       //

      Domain     &operator()(int, int);  // reference to domain i of type t
			        
      bool isCorrectEndString(MyString&); // compare to correct "END [DOMAINTYPE] [N]",
                                          // for reading of input file

      template<class Type> void newDom(Type*); // add new domain

      template<class Type> void delDom(Type*); // delete domain

      template<class Type> DomainType* getTypeAndId(Type*, int&, int&);

      template<class Type> char *name(Type*);  // gives:   "domainType i"

      template<class Type> char *key(Type*);   // gives:   "DOMAINTYPE i"

      template<class Type> int domainType(Type*);

    private:      

      MyString nameStr;
      
      int lastTypeAdded;

      DomainType *findType(int);

      DomainType* getTypeAndId(int, int&, int&);
};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



template<class Type> void DomainTree::newDom(Type *dom)
{
  ndom++;  (*this)[lastTypeAdded].dom.add((Domain*)dom);  return;
}





template<class Type> void DomainTree::delDom(Type *dom)
{
  int t, id;

  getTypeAndId(dom,t,id);

  ndom--;

  (*this)[t].dom.del(dom);

  return;
}





template<class Type> DomainType* DomainTree::getTypeAndId(Type *dom, int &t, int &id)
{
  this->resetTreeSearch();

  while (this->doTreeSearch())
  {
    if (this->searchItem->dom.n > 0)
    {
      id = this->searchItem->dom.getIndex(dom);
  
      if (id > -1) 
      {
        t = searchItem->tid;

        return this->searchItem;
      }
    }
  }
  return NULL;
}





template<class Type> char *DomainTree::name(Type *dom)
{
  int t, id;
	
  char tmp[10];
  
  nameStr.free().append(getTypeAndId(dom,t,id)->type);

  sprintf(tmp," %d",++id);

  nameStr.append(tmp);

  return nameStr.asCharArray();
}





template<class Type> char *DomainTree::key(Type *dom)
{
  int t, id;

  getTypeAndId(dom,t,id);

  char tmp[10], *domainKey[] = DOMAIN_KEY;

  nameStr.free().append(domainKey[t]);

  sprintf(tmp," %d",++id);

  nameStr.append(tmp);

  return nameStr.asCharArray();
}





template<class Type> int DomainTree::domainType(Type *dom)
{
  int t, id;

  getTypeAndId(dom,t,id);

  return t;  
}




//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



extern DomainTree domain;

template<class Type> int type(Type &dom)
{
  int t, id;

  domain.getTypeAndId(&dom,t,id);

  return t;
}




#endif



