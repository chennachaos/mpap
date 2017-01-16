
#ifndef incl_List_h
#define incl_List_h


#include "ContainerBase.h"
#include "LinkedListBase.h"
#include "MathVector.h"


// ListBase        - abstract list base class
//
// List            - linked list
// ListArray       - list based on an array
// SparseListArray - list based on an array of pointers to objects, some pointers may be NULL
// ListFixed       - list based on an array of fixed length
// ListInfinite    - infinite linked list

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

//  abstract base class for lists

template<class Type> class ListBase: public virtual ContainerBase, 
                                     public ListItem
{
  public:

    ListBase(void) { }

    virtual ~ListBase() { }

    virtual Type &operator[](int) = 0;
};


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// linked list

template<class Type> class List: public ListBase<Type>,
                                 public LinkedListBase<Type>
{
  public:

    List(void) { }

    virtual ~List() { }

    template<class Type2> List<Type> &operator=(Type2 &);

    virtual Type &operator[](int i) { return this->getItem(i); }

    virtual void reverseOrder(void) { LinkedListBase<Type>::reverseOrder(); return; }
};


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// list based on an array

template<class Type> class ListArray: public ListBase<Type>
{
  public:

    ListArray(void);

    virtual ~ListArray();

    Type *x;

    virtual void  setDim(int);

    virtual void  free(void);

    virtual Type &operator[](int);

    template<class Type2> ListArray<Type> &operator=(Type2 &);
};


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// list based on a fixed array

template<class Type, int dimension> class ListFixed: public ListBase<Type>
{
  public:

    ListFixed(void) { this->n = dimension; }

    virtual ~ListFixed() { }

    virtual Type &operator[](int i) { return x[i]; }

    template<class Type2> ListFixed<Type,dimension> &operator=(Type2&);

    Type x[dimension];
};


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// an infinite linked list

template<class Type> class ListInfinite: public List<Type>
{
  public:

    virtual Type &operator[](int i) 

      { while (i > this->n-1) this->add(new Type); return this->getItem(i); }
};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



template<class Type> template<class Type2> List<Type> &List<Type>::operator=(Type2 &b)
{
  this->free();

  for (int i=0; i<b.n; i++) { this->add(new Type); this->getItem(i) = b[i]; }

  return *this;
}



//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



template<class Type> ListArray<Type>::ListArray(void)
{
  x = NULL;
}





template<class Type> ListArray<Type>::~ListArray()
{
  //free();
}





template<class Type> void ListArray<Type>::setDim(int m)
{
  if (m != this->n)
  {
    free();

    x = new Type [m];

    this->n = m;
  }
  return;
}





template<class Type> void ListArray<Type>::free(void)
{
  if (x != NULL && this->n != 0) delete [] x;
 
  x = NULL;

  this->n = 0;

  return;
}





template<class Type> Type &ListArray<Type>::operator[](int i)
{
  if (i + 1 > this->n) cout << "ERROR in Type &ListArray<Type>::operator[]: i + 1 > n !\n\n";

  return *(x+i);
}





template<class Type> template<class Type2> ListArray<Type> &ListArray<Type>::operator=(Type2 &b)
{
  free();
  setDim(b.n);
  for (int i=0; i<b.n; i++) x[i] = b[i];
  return *this;
}



//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



template<class Type, int dimension> template<class Type2> ListFixed<Type,dimension>& 
                         ListFixed<Type,dimension>::operator=(Type2 &b)
{
  if (this->n == b.n) for (int i=0; i<this->n; i++) x[i] = b[i];

  else cout << "ERROR in Type &ListFixed<Type,dimension>::operator=:  b.n != n\n\n";

  return *this;
}


#endif

