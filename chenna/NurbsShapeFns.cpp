
#include "Debug.h"
#include "NurbsShapeFns.h"


//
//
//
//
//    NurbsShapeFns1D type1           //
//
//
//
NurbsShapeFns1Dtype1::NurbsShapeFns1Dtype1()
{
  if (debug) cout << " constructor NURBS Shape Functions 1D \n\n";

  J = 0.0;

}



NurbsShapeFns1Dtype1::NurbsShapeFns1Dtype1(int p)
{
  if (debug) cout << " constructor NURBS Shape Functions 1D \n\n";

    N.setDim(p+1);          N.zero();
    dN_dx.setDim(p+1);      dN_dx.zero();
    
    J = 0.0;

}

NurbsShapeFns1Dtype1::~NurbsShapeFns1Dtype1()
{
  if (debug) cout << " destructor NURBS Shape Functions 1D \n\n";

}

void NurbsShapeFns1Dtype1::zero()
{
  J = 0.0;
  N.zero(); dN_dx.zero();
}

//
//
//
//
//    NurbsShapeFns1D type2           //
//
//
//
NurbsShapeFns1Dtype2::NurbsShapeFns1Dtype2()
{
  if (debug) cout << " constructor NURBS Shape Functions 1D \n\n";

  J = 0.0;

}



NurbsShapeFns1Dtype2::NurbsShapeFns1Dtype2(int p)
{
  if (debug) cout << " constructor NURBS Shape Functions 1D \n\n";

    N.setDim(p+1);          N.zero();
    dN_dx.setDim(p+1);      dN_dx.zero();
    d2N_dx2.setDim(p+1);    d2N_dx2.zero();
    
    J = 0.0;

}

NurbsShapeFns1Dtype2::~NurbsShapeFns1Dtype2()
{
  if (debug) cout << " destructor NURBS Shape Functions 1D \n\n";

}

void NurbsShapeFns1Dtype2::zero()
{
  J = 0.0;
  N.zero(); dN_dx.zero();   d2N_dx2.zero();
}

//
//
//
//
//    NurbsShapeFns2D type1           //
//
//
//
NurbsShapeFns2Dtype1::NurbsShapeFns2Dtype1()
{
  if (debug) cout << " constructor NURBS Shape Functions 2D \n\n";

  J = 0.0;

}



NurbsShapeFns2Dtype1::NurbsShapeFns2Dtype1(int p, int q)
{
  if (debug) cout << " constructor NURBS Shape Functions 2D \n\n";

    N.setDim((p+1)*(q+1)); 		 N.zero();
    dN_dx.setDim((p+1)*(q+1));      dN_dx.zero();
    dN_dy.setDim((p+1)*(q+1));      dN_dy.zero();
    
    J = 0.0;


}

NurbsShapeFns2Dtype1::~NurbsShapeFns2Dtype1()
{
 // if (debug) cout << " destructor NURBS Shape Functions 2D \n\n";

//cout << " destructor NURBS Shape Functions 2D \n\n";


}

void NurbsShapeFns2Dtype1::zero()
{
  J = 0.0;
  N.zero(); dN_dx.zero();   dN_dy.zero();
}


//
//
//
//
//    NurbsShapeFns2D type2           //
//
//
//
NurbsShapeFns2Dtype2::NurbsShapeFns2Dtype2()
{
  if (debug) cout << " constructor NURBS Shape Functions 2D \n\n";

  J = 0.0;

}



NurbsShapeFns2Dtype2::NurbsShapeFns2Dtype2(int p, int q)
{
  if (debug) cout << " constructor NURBS Shape Functions 2D \n\n";

    N.setDim((p+1)*(q+1));  N.zero();
    dN_dx.setDim((p+1)*(q+1));      dN_dx.zero();
    d2N_dx2.setDim((p+1)*(q+1));    d2N_dx2.zero();
    dN_dy.setDim((p+1)*(q+1));      dN_dy.zero();
    d2N_dy2.setDim((p+1)*(q+1));    d2N_dy2.zero();
    
    J = 0.0;


}

NurbsShapeFns2Dtype2::~NurbsShapeFns2Dtype2()
{
  if (debug) cout << " destructor NURBS Shape Functions 2D \n\n";

}


void NurbsShapeFns2Dtype2::zero()
{
  J = 0.0;
  N.zero(); dN_dx.zero();   dN_dy.zero();
  d2N_dx2.zero();    d2N_dy2.zero();
}
