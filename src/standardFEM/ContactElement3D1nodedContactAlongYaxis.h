#ifndef incl_ContactElement3D1nodedContactAlongYaxis_h
#define incl_ContactElement3D1nodedContactAlongYaxis_h

#include "LagrangeElement.h"



class ContactElement3D1nodedContactAlongYaxis: public LagrangeElement
{
    public:

        ContactElement3D1nodedContactAlongYaxis();

        ~ContactElement3D1nodedContactAlongYaxis();

        virtual int getElmTypeNameNum()
        {
           return 36;
        }

        virtual void prepareElemData();

        virtual int  calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter=false);
};






#endif


