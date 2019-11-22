#ifndef incl_ContactElement3D1nodedContactAlongZaxis_h
#define incl_ContactElement3D1nodedContactAlongZaxis_h

#include "LagrangeElement.h"



class ContactElement3D1nodedContactAlongZaxis: public LagrangeElement
{
    public:

        ContactElement3D1nodedContactAlongZaxis();

        ~ContactElement3D1nodedContactAlongZaxis();

        virtual int getElmTypeNameNum()
        {
           return 37;
        }

        virtual void prepareElemData();

        virtual int  calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter=false);
};






#endif


