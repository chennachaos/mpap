
#ifndef incl_NurbsShapeFns_h
#define incl_NurbsShapeFns_h


#include "MyString.h"
#include "List.h"
#include "NurbsUtilitiesSURFACE.h"
using namespace std;

/*
class NurbsShapeFnsBase: public ListItem
{
  public:




    NurbsShapeFnsBase() {};

    ~NurbsShapeFnsBase() {};

    virtual void zero();
    
};
*/


class NurbsShapeFns1Dtype1: public ListItem
{
  public:

    VectorArray<double> N, dN_dx;
    
    double J;

    NurbsShapeFns1Dtype1();

    NurbsShapeFns1Dtype1(int p);

    virtual ~NurbsShapeFns1Dtype1();

    void zero();
    
};






class NurbsShapeFns1Dtype2: public ListItem
{
  public:

    VectorArray<double> N, dN_dx, d2N_dx2;

    double J;
    
    NurbsShapeFns1Dtype2();

    NurbsShapeFns1Dtype2(int p);

    virtual ~NurbsShapeFns1Dtype2();

    void zero();    
};





class NurbsShapeFns2Dtype1: public ListItem
{
  public:

    VectorArray<double> N, dN_dx, dN_dy;
    
    double J;


    NurbsShapeFns2Dtype1();

    NurbsShapeFns2Dtype1(int p, int q);

    virtual ~NurbsShapeFns2Dtype1();

    void zero();

//    void calcShapeFunctions(NurbsSURFACE *surf1, int ni, int nj, double u, double v);
};



class NurbsShapeFns2Dtype2: public ListItem
{
  public:

    VectorArray<double> N, dN_dx, d2N_dx2, dN_dy, d2N_dy2;
    
    double J;


    NurbsShapeFns2Dtype2();

    NurbsShapeFns2Dtype2(int p, int q);

    virtual ~NurbsShapeFns2Dtype2();

    void zero();
    
};
#endif

