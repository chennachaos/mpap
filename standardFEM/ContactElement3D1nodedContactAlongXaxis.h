#ifndef incl_ContactElement3D1nodedContactAlongXaxis_h
#define incl_ContactElement3D1nodedContactAlongXaxis_h

#include "LagrangeElement.h"



class ContactElement3D1nodedContactAlongXaxis: public LagrangeElement
{
    public:

        ContactElement3D1nodedContactAlongXaxis();

        ~ContactElement3D1nodedContactAlongXaxis();

        virtual int getElmTypeNameNum()
        {
           return 35;
        }

        virtual void prepareElemData();

        virtual int  calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter=false);
};






#endif


