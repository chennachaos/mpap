#ifndef incl_ContactElement2D1nodedContactAlongYaxis_h
#define incl_ContactElement2D1nodedContactAlongYaxis_h

#include "LagrangeElement.h"



class ContactElement2D1nodedContactAlongYaxis: public LagrangeElement
{
    public:

        ContactElement2D1nodedContactAlongYaxis();

        ~ContactElement2D1nodedContactAlongYaxis();

        virtual int getElmTypeNameNum()
        {
           return 34;
        }

        virtual void prepareElemData();

        virtual int  calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter=false);
};






#endif


