
#ifndef incl_ProfileListArray_h
#define incl_ProfileListArray_h

#include "List.h"
#include "MathVector.h"



template<class Type> class ProfileListArray: public ListItem
{
  public:

    ProfileListArray(void);

    virtual ~ProfileListArray();

    Type *data, **ptr;

    int n;

    inline Type* operator[](int i) { return ptr[i]; };

    inline int m(int i) { return ptr[i] - ptr[i-1]; };

    template<typename Type2> void generate(List< Vector<Type2> > &);

    void print(std::ostream &os = std::cout);

  private:

};






template<class Type> ProfileListArray<Type>::ProfileListArray(void)
{
  data = NULL;
  ptr  = NULL;
  n    = 0;

  return;
}






template<class Type> ProfileListArray<Type>::~ProfileListArray()
{
  if (data != NULL) delete [] data;

  if (ptr  != NULL) delete [] ptr;

  return;
}






template<class Type> template<typename Type2> 
  void ProfileListArray<Type>::generate(List< Vector<Type2> > &tmp)
{
  int i, j;

  n = tmp.n;

  j = 0; for (i=0; i<n; i++) j += tmp[i].n;

  ptr  = new Type* [n+1];

  data = new Type [j+1];

  j = 0; for (i=0; i<n; i++) { ptr[i] = &(data[j]); j += tmp[i].n; }

  ptr[n] = &(data[j]);

  for (i=0; i<n; i++) for (j=0; j<tmp[i].n; j++) ptr[i][j].setData(tmp[i][j]);

  return;
}






template<class Type> void ProfileListArray<Type>::print(std::ostream &os)
{
  int i, j;

  for (i=0; i<n; i++)
  {
    std::cout << "{";
    for (j=0; j<m(i+1)-1; j++) std::cout << ptr[i][j].nd << ",";
    if (m(i+1) > 0) std::cout << ptr[i][j].nd << "}\n"; else std::cout << " }\n";
  }
  return;
}







template<typename Type> std::ostream &operator<<(std::ostream &os, ProfileListArray<Type> &dat)
{
    dat.print(os);  return os;
}


#endif








