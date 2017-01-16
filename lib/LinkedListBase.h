
#ifndef incl_LinkedListBase_h
#define incl_LinkedListBase_h


#include <iostream>

#include "ContainerBase.h"
#include "MathBasic.h"


using namespace std;


static int nComp;
static int nSwap;
static int nSwap2;



//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// list item base class

class ListItem
{
  public:

    ListItem *next, *prev;

    template<class Type> bool operator>(Type &b)
      { cout << "operator > not defined for this listItem!\n\n"; return false; }

    template<class Type> bool operator<(Type &b)
      { cout << "operator < not defined for this listItem!\n\n"; return false; }
};


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

//  base class for linked lists

template<class Type> class LinkedListBase: public virtual ContainerBase
{
  public:

    LinkedListBase(void);

    virtual ~LinkedListBase();

    Type &firstItem(void) { return *((Type*)first.next); }

    Type &lastItem(void) { return *((Type*)last.prev); }

    void swap(ListItem*, ListItem*);

    void swap(int, int);

    void move(int, int);

    void reverseOrder(void);

    void add(Type*);
                                 
    Type* takeOut(Type*);
                                 
    Type* takeOut(int);
                                 
    void del(Type*);
                                 
    void del(int);
                                 
    void insert(Type*,int);
                                 
    void trunc(int);
                                 
    void free(void);

    void takeOver(LinkedListBase<Type>&);

    int  getIndex(Type*);

    void quickSort(bool largestFirst = false, int *data = NULL, int m0 = -1, int m1 = -1);

  protected:

    ListItem first, last;

    Type &getItem(int);

  private:

    bool linkedListDebugMode;

    void quickSortHelp(int, int, ListItem*, ListItem*, int*,
                       bool (LinkedListBase<Type>::*)(ListItem*,ListItem*),
                       bool (LinkedListBase<Type>::*)(ListItem*,ListItem*));

    void quickSortPartition(int, int, int&, ListItem*, ListItem*, ListItem*&, int*,
                            bool (LinkedListBase<Type>::*)(ListItem*,ListItem*),
                            bool (LinkedListBase<Type>::*)(ListItem*,ListItem*));

    bool isLarger(ListItem*, ListItem*);
    bool isSmaller(ListItem*, ListItem*);

    void abort(char*);

    int c;

    ListItem *curr;
};


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


template<class Type> LinkedListBase<Type>::LinkedListBase(void)
{
  first.prev = NULL;
  first.next = &last;
  last.next  = NULL;
  last.prev  = &first;  
  curr       = &first;
  this->n    = 0;
  c          = -1;

  linkedListDebugMode = true;

  return;
}





template<class Type> LinkedListBase<Type>::~LinkedListBase()
{
  free();

  return;
}





template<class Type> void LinkedListBase<Type>::swap(ListItem *itemI, ListItem *itemJ)
{
  nSwap++;

  ListItem *tmp;

  tmp         = itemI->next;
  itemI->next = itemJ->next;
  itemJ->next = tmp;

  tmp         = itemI->prev;
  itemI->prev = itemJ->prev;
  itemJ->prev = tmp;

  if      (itemI->prev == itemI) { itemI->prev = itemJ; itemJ->next = itemI; }
  else if (itemJ->prev == itemJ) { itemJ->prev = itemI; itemI->next = itemJ; }

  itemI->prev->next = itemI;
  itemJ->prev->next = itemJ;
  itemI->next->prev = itemI;
  itemJ->next->prev = itemJ;

  if (curr == itemI) curr = itemJ; else if (curr == itemJ) curr = itemI;

  return;
}





template<class Type> void LinkedListBase<Type>::swap(int i, int j)
{
  if (linkedListDebugMode) 
  {
    if (i == j)                         abort("swap: i == j");
    if (i < 0 || j < 0)                 abort("swap: i < 0 || j < 0");
    if (i > this->n-1 || j > this->n-1) abort("swap: i > m-1 || j > m-1");
  }

  swap((ListItem*)&getItem(i),(ListItem*)&getItem(j));

  return;
}





template<class Type> void LinkedListBase<Type>::move(int i, int j)
{
  insert(takeOut(i),j);

  return;
}





template<class Type> void LinkedListBase<Type>::reverseOrder(void)
{
  ListItem *item = first.next, *dmy;

  first.next = last.prev;
  last.prev  = item;

  for (int i=0; i<this->n; i++) 
  {
    dmy        = item->next;
    item->next = item->prev;
    item->prev = dmy;
    item       = dmy;
  }

  first.next->prev = &first;
  last.prev->next  = &last;

  if (c > -1) c = this->n - 1 - c;

  return;
}





template<class Type> void LinkedListBase<Type>::add(Type *newItem)
{
  ListItem *item  = (ListItem*) newItem;

  item->prev      = last.prev;
  last.prev->next = item;
  last.prev       = item;
  item->next      = &last;

  this->n++;

  return;
}





template<class Type> Type* LinkedListBase<Type>::takeOut(Type *takeItem)
{
  if (this->n < 1) return NULL;

  ListItem *item = (ListItem*) takeItem;

  if (linkedListDebugMode) 
  {
    while (item->prev != NULL) item = item->prev;

    if (&first != item) abort("takeOut: takeItem not in List!");

    item = (ListItem*) takeItem;
  }
  
  item->prev->next = item->next;
  item->next->prev = item->prev;

  this->n--;
  curr = last.prev;
  c = this->n - 1;

  return takeItem;
}





template<class Type> Type* LinkedListBase<Type>::takeOut(int i)
{
  if (this->n < 1) return NULL;
  
  if (linkedListDebugMode)
  {
    if (i < 0)         abort("takeOut: i < 0");
    if (i > this->n-1) abort("takeOut: i > n-1");
  }

  ListItem *item = &getItem(i);

  item->prev->next = item->next;
  item->next->prev = item->prev;

  this->n--;

  if (c == this->n) { curr = last.prev; c = this->n-1; } else curr = item->next;

  return (Type*) item;
}





template<class Type> void LinkedListBase<Type>::del(Type *delItem)
{
  delete takeOut(delItem);

  return;
}





template<class Type> void LinkedListBase<Type>::del(int i)
{
  delete takeOut(i);

  return;
}





template<class Type> void LinkedListBase<Type>::insert(Type *insertItem, int i)
{
  if (linkedListDebugMode) if (i > this->n) abort("insert: i > n");

  if (i == this->n) { add(insertItem); return; }
      
  ListItem *item = (ListItem*) insertItem, *itemNext = &getItem(i);

  item->next = itemNext;
  item->prev = itemNext->prev;

  itemNext->prev->next = item;
  itemNext->prev = item;

  this->n++;
  c++;
  
  return;
}





template<class Type> void LinkedListBase<Type>::trunc(int j)
{
  while (this->n > j) del(this->n-1);
	
  return;
}





template<class Type> void LinkedListBase<Type>::free(void)
{
  ListItem *item = first.next;

  while (item->next != NULL)
  {
    item = item->next;
    delete (Type*) item->prev;
  }

  first.next = &last;
  last.prev  = &first;  
  curr       = &first;
  this->n    = 0;
  c          = -1;
	
  return;
}





template<class Type> void LinkedListBase<Type>::takeOver(LinkedListBase<Type> &list2)
{
  free();

  first.next = list2.first.next; list2.first.next = &list2.last;  first.next->prev = &first;
  last.prev  = list2.last.prev;  list2.last.prev  = &list2.first; last.prev->next = &last;
  curr       = list2.curr;       list2.curr       = &list2.first;
  this->n    = list2.n;          list2.n          = 0;
  c          = list2.c;          list2.c          = -1;

  return;
}





template<class Type> int LinkedListBase<Type>::getIndex(Type *thisItem)
{
  ListItem *item = first.next, *findItem = (ListItem*) thisItem;

  int i = 0;

  if (findItem == NULL) return -1;

  while (item != NULL && item != findItem) { item = item->next; i++; }

  if (item == findItem && i < this->n) return i;

  return -1;
}





template<class Type> Type &LinkedListBase<Type>::getItem(int i)
{	
  if (linkedListDebugMode)
  {
    if (i < 0)         abort("getItem: i < 0");
    if (i > this->n-1) abort("getItem: i > n-1");
  }
  if (i<c) 
  { 
    if (i<c-i)    { curr = first.next; for (c=0;   c<i; c++) curr = curr->next; }
    else                               for (c=c;   c>i; c--) curr = curr->prev;  
  }
  else 
  {
    if (this->n-1-i<i-c) { curr = last.prev; for (c=this->n-1; c>i; c--) curr = curr->prev; }
    else                                     for (c=c;         c<i; c++) curr = curr->next;
  }
  
  return *((Type*) curr);
}





template<class Type> void LinkedListBase<Type>::quickSort(bool largestFirst, int *data,
                                                          int m0, int m1)
{
  int n0 = 0, n1 = this->n - 1; 

  nComp  = 0;
  nSwap  = 0;
  nSwap2 = 0;

  if (m0 >  -1) n0 = max(m0,n0);
  if (m1 >  -1) n1 = min(m1,n1);
  if (n0 >= n1) return;

  ListItem *item0 = &(getItem(n0)),
           *item1 = &(getItem(n1));

  if (largestFirst) quickSortHelp(n0,n1,item0,item1,data,
                                  &LinkedListBase<Type>::isSmaller,
                                  &LinkedListBase<Type>::isLarger);

  else              quickSortHelp(n0,n1,item0,item1,data,
                                  &LinkedListBase<Type>::isLarger,
                                  &LinkedListBase<Type>::isSmaller);

  //cout << "\n  LinkedListBase<Type>::quickSort:  nComp = "
  //     << nComp << ", nSwap = " << nSwap << ", nSwap2 = " << nSwap2 << "\n\n";
 
  return;
}





template<class Type> void LinkedListBase<Type>::quickSortHelp(int a, int b, 
                                 ListItem *itemA, ListItem *itemB,int *data,
                                 bool (LinkedListBase<Type>::*isThis)(ListItem*,ListItem*),
                                 bool (LinkedListBase<Type>::*isThat)(ListItem*,ListItem*))
{
  if (a > b-1) return;
  
  if (a == b-1) { if ((this->*isThis)(itemA,itemB))

    { swap(itemA,itemB); simpleSwap(data,a,b); } return; }

  ListItem *itemC, *left = itemA->prev, *right = itemB->next;

  int p;

  quickSortPartition(a,b,p,itemA,itemB,itemC,data,isThis,isThat);

  quickSortHelp(a,p-1,left->next,itemC->prev,data,isThis,isThat);
  quickSortHelp(p+1,b,itemC->next,right->prev,data,isThis,isThat);

  return;
}





template<class Type> void LinkedListBase<Type>::quickSortPartition(int a, int b, int &p,
                                 ListItem *itemA, ListItem *itemB, ListItem *&itemP, int *data,
                                 bool (LinkedListBase<Type>::*isThis)(ListItem*,ListItem*),
                                 bool (LinkedListBase<Type>::*isThat)(ListItem*,ListItem*))
{
  p = (a + b + (a+b) % 2) / 2;

  //cout << a << " < " << p << " < " << b << "\n";

  itemP = itemA;

  ListItem *tmp,
           *top    = itemB->next,
           *bottom = itemA->prev;

  for (int i=a; i<p; i++) itemP = itemP->next;

  if (isSmaller(itemP,itemA) && isSmaller(itemP,itemB))

    { if (isSmaller(itemA,itemB)) { p = a; itemP = itemA; } else { p = b; itemP = itemB; } }

  else if (isLarger(itemP,itemA) && isLarger(itemP,itemB))

    { if (isLarger(itemA,itemB)) { p = a; itemP = itemA; } else { p = b; itemP = itemB; } }

  if (p != a) { nSwap2++; swap(itemP,itemA); simpleSwap(data,p,a); }

  top    = top->prev;
  bottom = bottom->next->next;

  int ib = a+1, it = b;

  while (1)
  {
    while ((this->*isThis)(top,itemP))    { top = top->prev; it--; }
    while (ib < it) if ((this->*isThat)(bottom,itemP)) { bottom = bottom->next; ib++; } else break;

    if (it > ib)
    {
      swap(top,bottom);
      tmp = bottom;
      bottom = top->next;
      top = tmp->prev;
      simpleSwap(data,it--,ib++);
    }
    else break;
  }
  nSwap2++; swap(itemP,top); simpleSwap(data,a,it);

  p = it;

  return;
}





template<class Type> bool LinkedListBase<Type>::isLarger(ListItem *a, ListItem *b)
{
  nComp++;

  return *((Type*)a) > *((Type*)b);
}





template<class Type> bool LinkedListBase<Type>::isSmaller(ListItem *a, ListItem *b)
{
  nComp++;

  return *((Type*)a) < *((Type*)b);
}





template<class Type> void LinkedListBase<Type>::abort(char* msg)
{
  cout << "  ERROR: LinkedListBase<Type>::" << msg << "\n\n";

  exit(0);

  return;
}



//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::




#endif

