#ifndef incl_ContactElement2D1nodedContactWithYaxis_h
#define incl_ContactElement2D1nodedContactWithYaxis_h

#include "LagrangeElement.h"



class ContactElement2D1nodedContactWithYaxis: public LagrangeElement
{
    public:

        ContactElement2D1nodedContactWithYaxis();

        ~ContactElement2D1nodedContactWithYaxis();

        virtual int getElmTypeNameNum()
        {
           return 34;
        }

        virtual void prepareElemData();

        virtual int  calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal);

        //void  assembleElementVector(int ind, bool flag, double* rhs);

        //void  AssembleMatrixAndVector(int, SparseMatrixXd&, double*);
};






#endif


