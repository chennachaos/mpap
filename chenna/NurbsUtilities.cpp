
#include "NurbsUtilities.h"
#include "myConstants.h"

#include <iostream>

using namespace std;




void findunique(KNOTVECTOR& U, VectorArray<double>& XX)
{
  int count1 = 0;

  Vector<double> X;

  X.append(U[0]);
  for(unsigned int i=1;i<U.n;i++)
  {
    if( X[count1] < U[i] )
    {
      X.append(U[i]);
      count1++;
    }
  }

  XX = X;

}




void finduniqueInt(VectorArray<int>& U, VectorArray<int>& XX)
{
  int count1 = 0;

  Vector<int> X;

  X.append(U[0]);
  for(unsigned int i=1;i<U.n;i++)
  {
    if( X[count1] < U[i] )
    {
      X.append(U[i]);
      count1++;
    }
  }

  XX = X;

}


void create_vector(double start, double end, double incr, VectorArray<double>& uuu)
{
  int steps, ii;

  steps = static_cast<int>(((end-start)/incr) + 1);
  uuu.setDim(steps);

  uuu[0] = start;
  for(ii=1;ii<steps;ii++)
    uuu[ii] = uuu[ii-1] + incr;

  return;
}



// Refine the knot vector 'U' by subdividing every non-zero knot interval 'num' number of times

void create_vector2(KNOTVECTOR& U, int num, KNOTVECTOR& uu1)
{
  VectorArray<double> X1;
  findunique(U, X1);

  double count, incr;
  int  ii, jj;

  uu1.setDim((X1.n-1)*(num-1) + X1.n);

  uu1[0] = U[0]; count = 1;
  for(ii=0;ii<X1.n-1;ii++)
  {
     incr = (X1[ii+1] - X1[ii])/num;
     for(jj=0;jj<num;jj++)
     {
        uu1[count] = uu1[count-1] + incr;
        count++;
      }
  }

  return;
}






void GenKnotVecForRefining(KNOTVECTOR& U, int Nsub, KNOTVECTOR& X)
{
  VectorArray<double> X1;

  findunique(U, X1);

  double a, b, dummy;
  int count=0;

  dummy = (X1.n-1) * (pow(2.0,Nsub)-1);

  X.setDim((int) dummy); // typecasting

  for(unsigned int ii=0;ii<X1.n-1;ii++)
  {
    a = X1[ii];
    b = X1[ii+1];
    for(int j=1;j<=pow(2.0,Nsub)-1;j++)
    {
      X[count+j-1] = a + j*(b-a)/pow(2.0,Nsub);
    }
    count = count + (int) pow(2.0,Nsub) -1;
  }
}


void getLowerOrderKnotVector(VectorArray<double>& U, int p, VectorArray<double>& VV)
{
    Vector<double> V1;

    int m = U.n - 1, count = 0, ii;

    for(ii=0;ii<p;ii++)
	V1.append(0.0);

    for(ii=p+1;ii<m-p;ii++)
    {
	if( U[ii] < U[ii+1] ) 	// case of a single knot
	   V1.append(U[ii]);
	else 			// case of multiple knots
	{
	   V1.append(U[ii]);
	   ii++;
	}
    }

    for(ii=0;ii<p;ii++)
	V1.append(U[m]);

    VV = V1;

  return;
}



// calculates the distance between two points P1 and P2 in Euclidean domain

double CalcDist(EPOINT& P1, EPOINT& P2)
{
  double dx, dy, dz;
  dx = (P1.x - P2.x);
  dy = (P1.y - P2.y);
  dz = (P1.z - P2.z);

  return sqrt(dx * dx + dy * dy + dz * dz);
}



int EntryExists(MatrixSparse<double>& mtx1, int& r1, int& c1)
{
  for(int i=0;i<mtx1.x.n;i++)
  {
    if(mtx1.row[i] == r1)
    {
      if(mtx1.col[i] == c1)
        return i+1;
      else
      {
         for(int j=i;j<mtx1.x.n;j++)
         {
           if(mtx1.row[j] == r1)
           {
             if(mtx1.col[j] == c1)
               return j+1;
           }
         }
      }
    }
  }

  return -1;
}


void SortArray(VectorArray<double>& X1)
{
  int i, j;
  double a=0.0;

  for(j=0;j<X1.n;j++)
  {
    a = X1[j];
    i=j;
    while(i>0 && X1[i-1] > a)
    {
       X1[i] = X1[i-1];
       i--;
    }
    X1[i] = a;
  }

  return;
}


void SortArrayInt(VectorArray<int>& X1)
{
  int i, j;
  int a;

  for(j=0;j<X1.n;j++)
  {
    a = X1[j];
    i=j;
    while(i>0 && X1[i-1] > a)
    {
       X1[i] = X1[i-1];
       i--;
    }
    X1[i] = a;
  }

  return;
}




// computes the 'arctan'
double ArcTan(double x, double y)
{
  double temp = 0.0;

  if( CompareDoubles(x, 0.0) )
  {
     if( CompareDoubles(y, 0.0) )
        temp = 0.0;
     else if( y > 0.0 )
        temp = 0.5*PI;
     else // if( y < 0.0 )
        temp = 1.5*PI;
  }
  else if( x > 0.0)
  {
     if( CompareDoubles(y, 0.0) )
        temp = 0.0;
     else if( y > 0.0 )
        temp = atan(y/x);
     else
        temp = 2*PI - atan( abs(y/x) );
  }
  else // if( x < 0.0)
  {
     if( CompareDoubles(y, 0.0) )
        temp = PI;
     else if( y > 0.0 )
        temp = PI - atan( abs(y/x) );
     else
        temp = PI + atan( abs(y/x) );
  }

  return temp;
}



int binarySearch(int* sortedArray, int first, int last, int key)
{
   // function:
   //   Searches sortedArray[first]..sortedArray[last] for key.  
   // returns: index of the matching element if it finds key, 
   //         otherwise  -(index where it could be inserted)-1.
   // parameters:
   //   sortedArray in  array of sorted (ascending) values.
   //   first, last in  lower and upper subscript bounds
   //   key         in  value to search for.
   // returns:
   //   index of key, or -insertion_position -1 if key is not 
   //                 in the array. This value can easily be
   //                 transformed into the position to insert it.

   int mid;
   while(first <= last)
   {
       mid = (first + last) / 2;  // compute mid point.
       if (key > sortedArray[mid]) 
           first = mid + 1;  // repeat search in top half.
       else if (key < sortedArray[mid]) 
           last = mid - 1; // repeat search in bottom half.
       else
           return mid;     // found it. return position /////
   }
   return -1;    // failed to find key
}




// computes array XX which contains entries only in X1 but not in X2

void sub2VectorArraysInt(VectorArray<int>& X1, VectorArray<int>& X2, VectorArray<int>& XX)
{
  int i, j;
  Vector<int> temp;
  bool flag;

  for(i=0;i<X1.n;i++)
  {
     flag = false;
     for(j=0;j<X2.n;j++)
     {
        if(X1[i] == X2[j])
           flag = true;

        if(flag)
           break;
     }
     if(!flag)
       temp.append(X1[i]);
  }

  XX = temp;

  return;
}





double vonMises3D(double* stre)
{
  return  sqrt(0.5* (pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6.0*(stre[3]*stre[3]+stre[4]*stre[4]+stre[5]*stre[5])));
}




// computes the Cross product of two vectors stored in P1 and P2

void CrossProduct(EPOINT& P1, EPOINT& P2, EPOINT& P3 )
{
   P3.x  = P1.y*P2.z - P2.y*P1.z;
   P3.y  = P2.x*P1.z - P1.x*P2.z;
   P3.z  = P1.x*P2.y - P1.y*P2.x;
   return ;
}



double DotProduct(EPOINT& P1, EPOINT& P2)
{
  return ( (P1.x * P2.x) + (P1.y * P2.y) + (P1.z * P2.z) );
}







