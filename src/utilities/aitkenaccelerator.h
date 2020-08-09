

#ifndef AITKEN_ACCELERATOR_H
#define AITKEN_ACCELERATOR_H



inline  int aitken_accelerator(VectorXd&  vec)
{
  // Applies Aitken acceleration to the sequence in x,
  // returning a sequence of n-2 new approximations,
  // where n is the length of the vector x.

  int n = vec.rows();
  VectorXd  a(n);

  //for(int k=0; k<n; k++)
  //{
    //a(k) = x(k) - (x(k+1) - x(k))^2/(x(k+2) - 2*x(k+1) + x(k));
  //}

  return 1;
}