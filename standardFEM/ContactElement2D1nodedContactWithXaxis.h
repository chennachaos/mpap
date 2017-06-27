#ifndef incl_ContactElement2D1nodedContactWithXaxis_h
#define incl_ContactElement2D1nodedContactWithXaxis_h

#include "LagrangeElement.h"


class ContactElement2D1nodedContactWithXaxis: public LagrangeElement
{
    public:

        ContactElement2D1nodedContactWithXaxis();

        ~ContactElement2D1nodedContactWithXaxis();

        virtual int getElmTypeNameNum()
        {  return 33;     }

        virtual void prepareElemData();

        virtual int  calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal);

        //void  assembleElementVector(int ind, bool flag, double* rhs);

        //void  AssembleMatrixAndVector(int, SparseMatrixXd&, double*);
};






#endif


