#ifndef incl_ContactElement2D1nodedContactAlongXaxis_h
#define incl_ContactElement2D1nodedContactAlongXaxis_h

#include "LagrangeElement.h"


class ContactElement2D1nodedContactAlongXaxis: public LagrangeElement
{
    public:

        ContactElement2D1nodedContactAlongXaxis();

        ~ContactElement2D1nodedContactAlongXaxis();

        virtual int getElmTypeNameNum()
        {  return 33;     }

        virtual void prepareElemData();

        virtual int  calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter=false);
};






#endif


