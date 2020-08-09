
#ifndef incl_MathBasic_h
#define incl_MathBasic_h


#include <cmath>
#include <cstdlib>
#include <iostream>
#include  <stdio.h>



static const double pi = 3.141592653589793238462643;




template<typename Type> inline void simpleSwap(Type *x, int a, int b)
{
  if (x == NULL) return;

  Type tmp = x[a];
  
  x[a] = x[b];
  x[b] = tmp;

  return;
}




template<typename Type> inline int roundToInt(Type x)
{
  return (int) floor(x + 0.5);
}




inline int intDiv(int a, int b)
{
  if (b == 0) return 0;
	
  return (a - a % b) / b;
}




template<typename Type> inline Type norm(Type *x, int dim)
{
  Type y = (Type) 0;

  for (int i=0; i<dim; i++) y += x[i] * x[i];
	
  return (Type) sqrt(y);
}




template<typename Type> Type inline dot3(Type *x1, Type *x2)
{
  return x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
}




template<typename Type> Type inline dot(Type *v1, Type *v2, int n)
{
  double p = 0.;

  for (int i=0; i<n; i++) p += *(++v1) * *(++v2);

  return p;
}




template<typename Type> void inline cross3(Type *c, Type *x1, Type *x2)
{
  c[0] = x1[1]*x2[2]-x1[2]*x2[1];
  c[1] = x1[2]*x2[0]-x1[0]*x2[2];
  c[2] = x1[0]*x2[1]-x1[1]*x2[0];

  return;
}




double inline myAcos(double x)
{
  double tol = 1.e-12;

  if      (x > 1. - tol) return 0.;
  else if (x < tol - 1.) return pi;
  else                   return acos(x);
}





double inline power(double x, int n)
{
  double pwr = 1.;

  int i;

  for (i=0; i<n; i++) pwr *= x;
	
  return pwr;
}





template<typename Type> inline bool contains(Type *v, Type x, int n, int *p = NULL)
{
  int i = 0;

  while (i < n && v[i] != x) i++;

  if (p != NULL) *p = i;

  if (i < n) return true;

  return false;
}





void inline generateRandomDbl(double *r, int n, double rmn = 0., double rmx = 1.)
{
  int    i, mn = RAND_MAX, mx = 0, *random = new int [n];
  double d = rmx - rmn;

  for (i=0; i<n; i++) { random[i] = rand(); if (mn>random[i]) mn = random[i]; }
  for (i=0; i<n; i++) { random[i] -= mn;    if (mx<random[i]) mx = random[i]; }
 
  d /= double(mx);

  for (i=0; i<n; i++) r[i] = double(random[i])*d+rmn;

  delete [] random;

  return;
}






void inline splitInterval(int *array, int &a, int &b, int numb)
{
  int c = a + intDiv(b-a,2);

  if (numb < array[c]) b = c; else a = c; 

  return;
}







int inline getIndex(int *array, int n, int numb)
{
  // array contains distinct, ordered integer numbers (smallest first)

  int a = 0, b = n;

  if (numb < array[a] || numb > array[b]) return -1;

  while (b - a > 4) splitInterval(array,a,b,numb); //std::cout << a << "-" << b << "\n"; }

  while (array[a] != numb) a++;

  return a;  
}









template<typename Type> inline bool isElementOf(Type x, Type *a, int n)
{
  int i = 0;

  while (i < n) if (x == a[i]) return true; else i++;

  return false;
} 






#endif

