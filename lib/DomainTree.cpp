
#include <iostream>

#include "FunctionsProgram.h"
#include "DomainType.h"
#include "DomainTree.h"
#include "Debug.h"





DomainTree::DomainTree()
{
  reset();
}





DomainTree::~DomainTree()
{
  free();
}





void DomainTree::reset(void)
{
  ndom          = 0;
  mxtid         = -1;
  lastTypeAdded = -1;

  return;
}





int DomainTree::nDomainOfType(int t)
{
  DomainType *domType = findType(t); 

  if (domType == NULL) return 0;

  return domType->dom.n;
}





int DomainTree::nType(void)
{
  int c = 0;

  this->resetTreeSearch();

  while (this->doTreeSearch()) if (this->searchItem->dom.n > 0) c++;

  return c;
}





DomainType *DomainTree::newType(int newtid, int partid)
{     
  lastTypeAdded = newtid;

  if (newtid > mxtid) mxtid = newtid;
  
  DomainType *domType = findType(newtid);

  if (domType == NULL)	  
  {
    domType = new DomainType(newtid);

    if (partid < 0) add(domType,NULL); else add(domType,findType(partid));

    //this->print();

    return domType;
  }  
  else return NULL;
}
  




void DomainTree::delType(DomainType *delDomType)
{
  this->del(delDomType);

  return;
}





DomainType &DomainTree::operator[](int t)
{
  DomainType *domType = findType(t);

  if (domType == NULL) prgError(1,"DomainTree::operator[]","This type does not exist!");

  return *domType; 
}





Domain &DomainTree::operator()(int i)
{
  int t, id;

  DomainType *domType = getTypeAndId(i,t,id);

  if (domType == NULL) prgError(1,"DomainTree::operator()","invalid index!");

  return domType->dom[id];
}





Domain &DomainTree::operator()(int t, int id)
{
  char fct[] = "DomainTree::operator()(type,id)";

  if (ndom < 1) prgError(1,fct,"no domains loaded!");

  DomainType *domType = findType(t);

  if (domType != NULL)
  {
    if (domType->dom.n > id && id > -1) return domType->dom[id];

    if (domType->dom.n > 0)
    {
      prgWarning(1,fct,"invalid domain id; choosing last domain of this type!");

      //id = domType->dom.n - 1;

      return domType->dom.lastItem();
    }
  }
  prgWarning(2,fct,"no domains of this type; choosing first in tree!");

  this->resetTreeSearch();

  while (this->doTreeSearch())
  {
    if (this->searchItem->dom.n > 0) break;
  }

  //t = searchItem->tid;

  //id = 0;

  return searchItem->dom[0];
}





bool DomainTree::isCorrectEndString(MyString &line)
{
  if (lastTypeAdded < 0) return false;
  
  char tmp[100], *domKey[] = DOMAIN_KEY;
  
  sprintf(tmp,"END %s %i",domKey[lastTypeAdded],domain.nDomainOfType(lastTypeAdded));

  line.stripToMin();
  
  if (line != tmp) return false;
  
  return true;
}





DomainType *DomainTree::findType(int t)
{
  this->resetTreeSearch();

  while (this->doTreeSearch()) if (this->searchItem->tid == t) return this->searchItem;

  return NULL;
}





DomainType* DomainTree::getTypeAndId(int i, int &t, int &id)
{
  if (i > ndom-1) return NULL;

  int c = 0;

  this->resetTreeSearch();

  while (this->doTreeSearch())
  {
    c += this->searchItem->dom.n; 

    if (c > i) break;
  }
  
  t = this->searchItem->tid;

  id = c - i - 1;

  return this->searchItem;
}





