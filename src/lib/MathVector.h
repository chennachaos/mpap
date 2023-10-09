
#ifndef incl_MathVector_h
#define incl_MathVector_h


#include <cmath>

#include "ContainerBase.h"
#include "LinkedListBase.h"
#include "MathBasic.h"



//  vector classes
//
//  VectorBase     - abstract base class
//  myVector       - vector based on a linked list
//  VectorArray    - vector based on an array
//  VectorFixed    - vector based on an array of fixed length
//  VectorInfinite - vector based on an infinite linked list


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// vector coefficient class needed for myVector

template<typename Type> class VectorCoeff: public ListItem
{
  public:

    VectorCoeff(void) { val = (Type) 0; }

    virtual ~VectorCoeff() { }

    Type val;

    template<class Type2> bool operator>(Type2 &b) { return val > b.val; }
    template<class Type2> bool operator<(Type2 &b) { return val < b.val; }
};


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// abstract vector base class

template<typename Type> class VectorBase: public virtual ContainerBase, public ListItem
{
  public:

    VectorBase(void) { }
      
    virtual ~VectorBase() { }
      
    virtual Type &operator[](int) = 0;

    virtual Type &operator()(int);

    template<class Type2> bool operator==(Type2 &);
    
    virtual Type &firstCoeff(void);

    virtual Type &lastCoeff(void);

    virtual void print(std::ostream &os = std::cout);
    
    virtual void zero(void);

    virtual Type norm(void);

    virtual bool contains(Type, int *indx = NULL);

    template<class Type2> bool containsAllOf(Type2 &);

    template<class Type2> bool containsAtLeastOneOf(Type2 &, int *indx = NULL, int *indx2 = NULL);

    template<class Type2> bool containsSomeOf(int, Type2 &);

    template<class Type2> bool intersection(Type2 &, VectorBase<int> &, VectorBase<int> &);

    virtual bool distinct(void);
};


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// vector based on a linked list

template<typename Type> class myVector: public VectorBase<Type>, 
                                      public LinkedListBase< VectorCoeff<Type> >
{
  public:

    myVector(void);
       
    virtual ~myVector();
       
    Type &operator[](int);
       
    template<class Type2> myVector<Type> &operator=(Type2 &);

    Type* append(Type val = (Type) 0);

    void insert(Type,int);

    void del(int);

    Type &firstCoeff(void);

    Type &lastCoeff(void);

    void reverseOrder(void)

      { LinkedListBase< VectorCoeff<Type> >::reverseOrder(); return; }
};


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// vector based on an array

template<typename Type> class VectorArray: public VectorBase<Type>
{
  public:

    VectorArray(void);
       
    virtual ~VectorArray();
       
    virtual Type &operator[](int);       
    
    template<class Type2> VectorArray<Type> &operator=(Type2 &);

    void setDim(int);
    
    virtual void free(void);

    Type *x;
};

//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// simple vector based on array of fixed length

template<typename Type, int dimension> class VectorFixed: public VectorBase<Type>
{
  public:

    VectorFixed(void) { this->n = dimension; }

    virtual ~VectorFixed() { }

    template<class Type2> VectorFixed<Type,dimension> &operator=(Type2 &b)

      { if (this->n == b.n) for (int i=0; i<this->n; i++) x[i] = b[i]; return *this; }

    virtual Type &operator[](int i) { return x[i]; }

    Type x[dimension];
};


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// vector based on an infinite linked list

template<typename Type> class VectorInfinite: public myVector<Type>
{
  public:

    virtual Type &operator[](int i) 

      { while (i > this->n-1) this->append((Type) 0);

        return LinkedListBase< VectorCoeff<Type> >::getItem(i).val; }
};


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


template<typename Type> Type &VectorBase<Type>::operator()(int i)
{
  return (*this)[i-1];
}



template<typename Type> Type &VectorBase<Type>::firstCoeff(void)
{
  return (*this)[0];
}




template<typename Type> Type &VectorBase<Type>::lastCoeff(void)
{
  return (*this)[this->n-1];
}




template<typename Type> void VectorBase<Type>::print(std::ostream &os)
{
  os << "{"; 
  if (this->n > 0) os << (*this)[0]; else os << " ";
  for (int i=1; i<this->n; i++) os << "," << (*this)[i]; 
  os << "}";
  return;
}




template<typename Type> void VectorBase<Type>::zero(void)
{
  for (int i=0; i<this->n; i++) (*this)[i] = (Type) 0;
  return;
}




template<typename Type> Type VectorBase<Type>::norm(void)
{
  std::cout << " Warning! norm only works for myVector<double> or myVector<int>!\n\n";
      
  return (Type) 0;
}




template<> inline double VectorBase<double>::norm(void)
{
  double y = 0.;

  for (int i=0; i<this->n; i++) y += (*this)[i] * (*this)[i];

  return sqrt(y);
}




template<> inline int VectorBase<int>::norm(void)
{
  int y = 0;

  for (int i=0; i<this->n; i++) y += (*this)[i] * (*this)[i]; 

  return roundToInt(sqrt((double)y));
}




template<> inline unsigned int VectorBase<unsigned int>::norm(void)
{
  unsigned int y = 0;

  for (int i=0; i<this->n; i++) y += (*this)[i] * (*this)[i]; 

  return (unsigned int) roundToInt(sqrt((double)y));
}




template<typename Type> bool VectorBase<Type>::contains(Type xx, int *indx)
{
  int i = 0; while (i < this->n) if ((*this)[i] == xx) break; else i++;
  if (i < this->n)
  {
    if (indx != NULL) *indx = i;
    return true;
  }
  if (indx != NULL) *indx = -1;
  return false;
}




template<typename Type> template<class Type2> bool VectorBase<Type>::containsAllOf(Type2 &vec2)
{
  for (int i=0; i<vec2.n; i++)  if (!this->contains(vec2[i])) return false;

  return true;
}




template<typename Type> template<class Type2> 
   bool VectorBase<Type>::containsAtLeastOneOf(Type2 &vec2, int *indx, int *indx2)
{
  int i = 0, j;
  while (i < this->n)
  {
    j = 0;
    while (j < vec2.n && (*this)[i] != vec2[j]) j++;
    if (j < vec2.n) 
    {
      if (indx  != NULL) *indx  = i;
      if (indx2 != NULL) *indx2 = j;
      return true;
    }
    i++;
  }
  return false;
}





template<typename Type> template<class Type2>
   bool VectorBase<Type>::containsSomeOf(int m, Type2 &vec2)
{
  int i = 0, j, c = 0;
  while (i < this->n)
  {
    j = 0;
    while (j < vec2.n && (*this)[i] != vec2[j]) j++;
    if (j < vec2.n) { c++; if (c == m) return true; }
    i++;
  }
  return false;
}





template<typename Type> bool VectorBase<Type>::distinct(void)
{
  int i, j;

  for (i=1; i<this->n; i++)
    for (j=0; j<i; j++)
      if ((*this)[i] == (*this)[j]) return false;

  return true;
}





template<typename Type> template<class Type2> bool VectorBase<Type>::operator==(Type2 &b)
{
  if (this->n != b.n) return false;

  int i = 0;

  while (i < this->n) if ((*this)[i] == b[i]) i++; else return false;

  return true;
}





template<typename Type> std::ostream &operator<<(std::ostream &os, VectorBase<Type> &vec)
{
  vec.print(os);  return os;
}




 
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


template<typename Type> myVector<Type>::myVector(void)
{
  return;
}



template<typename Type> myVector<Type>::~myVector()
{
  this->free();   

  return;
}



template<typename Type> Type &myVector<Type>::operator[](int i)
{
  return LinkedListBase< VectorCoeff<Type> >::getItem(i).val;
}



template<typename Type> template<class Type2> myVector<Type> &myVector<Type>::operator=(Type2 &b)
{
  this->free();

  for (int i=0; i<b.n; i++) append(b[i]);

  return *this;
}



template<typename Type> Type *myVector<Type>::append(Type x)
{       
  VectorCoeff<Type> *c = new VectorCoeff<Type>;

  c->val = x;

  this->add(c);

  return &(c->val);
}




template<typename Type> void myVector<Type>::insert(Type x, int pos)
{       
  VectorCoeff<Type> *c = new VectorCoeff<Type>;

  c->val = x;

  LinkedListBase< VectorCoeff<Type> >::insert(c,pos);
	
  return;
}




template<typename Type> void myVector<Type>::del(int pos)
{       
  LinkedListBase< VectorCoeff<Type> >::del(pos);
	
  return;
}




template<typename Type> Type &myVector<Type>::firstCoeff(void)
{
  return ((VectorCoeff<Type>*)(this->first.next))->val;
}




template<typename Type> Type &myVector<Type>::lastCoeff(void)
{
  return ((VectorCoeff<Type>*)(this->last.prev))->val;
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



template<typename Type> VectorArray<Type>::VectorArray(void)
{
  x = NULL;

  return;
}



template<typename Type> VectorArray<Type>::~VectorArray()
{
  free();

  return;
}



template<typename Type> Type &VectorArray<Type>::operator[](int i)
{
  if (i < 0 || i+1 > this->n) 

    cout << "  ERROR: VectorArray<Type>::operator[]: invalid index!\n\n";

  return x[i];
}






template<typename Type> template<class Type2> VectorArray<Type> &VectorArray<Type>::operator=(Type2 &b)
{
  setDim(b.n);

  for (int i=0; i<this->n; i++) x[i] = b[i];

  return *this;
}







template<typename Type> void VectorArray<Type>::setDim(int d)
{       
  if (d != this->n)
  {
    free();
    x = new Type [d];
    this->n = d;
  }
  return;
}



template<typename Type> void VectorArray<Type>::free(void)
{       
  if (x != NULL && this->n != 0) delete [] x; x = NULL; this->n = 0;
	
  return;
}




#endif



