
//#include "Element.h"
#include "FunctionsProgram.h"
#include "FunctionsDeniz.h"
#include "ElementGroup.h"
#include "PropertyTypeEnum.h"

void hyplasmaterial_(double *matData, double *F, double *F33, double *stre, double *cc,
                     double *iv1, double *iv2,
                     int *matid, int *finite, int *sss, int *isw, int *error,
					 Element *elm)

{

int hyplasmatid, nmatData, iiter;
 double thkgp;

if(*sss!=1 && *sss!=2)prgError(1,"hyplasmaterial","hyplas material not available for this state!");


switch(*matid)
{
if (*isw==2) matData[0]=0.0;
/*HyplasELASTC*/case 20: hyplasmatid=1;  nmatData=3; break;
/*HyplasTRESCA*/case 21: hyplasmatid=2;  nmatData=3+2*(int) matData[0]; break;
/*HyplasVMISES*/case 22: hyplasmatid=3;  nmatData=3+2*(int) matData[0]; break;
/*HyplasMOHCOU*/case 23: hyplasmatid=4;  nmatData=6+2*(int) matData[0];
         if(*sss<2)      prgError(1,"hyplasmaterial","HyplasMOHCOU is not defined for plane stress state!"); break;
/*HyplasDRUPRA*/case 24: hyplasmatid=5;  nmatData=6+2*(int) matData[0]; break;
/*HyplasCAPDP*/case 25: hyplasmatid=6;  nmatData=7+2*(int) matData[0];
         if(*sss<2)      prgError(1,"hyplasmaterial","HyplasCAPDP is not defined for plane stress state!") ; break;
/*HyplasLEMDAM*/case 26: hyplasmatid=7;  nmatData=5+2*(int) matData[0];
         if(*sss<2)      prgError(1,"hyplasmaterial","HyplasLEMDAM is not defined for plane stress state!"); break;
/*HyplasDAMELA*/case 27: hyplasmatid=8;  nmatData=7;
         if(*sss<2)      prgError(1,"hyplasmaterial","HyplasDAMELA is not defined for plane stress state!");  break;
/*HyplasOGDEN*/case 28: hyplasmatid=20; nmatData=5+2*(int) matData[0];
         if(*sss<2)     prgWarning(1,"hyplasmaterial","HyplasOGDEN does not return consistent tangent in with plane stress state!");
         if(*finite==0)  prgError(1,"hyplasmaterial","HyplasOGDEN  is not defined for small strains!"); break;
/*HyplasPDSCRY*/case 29: hyplasmatid=40; nmatData=1+2*(int) matData[0];
         if(*sss<2)     prgError(1,"hyplasmaterial","HyplasPDSCRY is not defined for plane stress state!");
         if(*finite==0)  prgError(1,"hyplasmaterial","HyplasPDSCRY is not defined for small strains!"); break;
/*HyplasVMMIXD*/case 30: hyplasmatid=60; nmatData=3+2*(int) matData[0];
         if(*sss<2)     prgError(1,"hyplasmaterial","HyplasVMMIXD is not defined for plane stress state!");
         if(*finite!=0)  prgError(1,"hyplasmaterial","HyplasVMMIXD is not defined for finite strains!");  break;
default:                 prgError(1,"hyplasmaterial","this is not a recognized hyplas material!"); break;
}

    if (*isw==2)
{
	hyplas_(matData,F,F33,stre,cc,iv1,iv2,NULL, //8 double
              finite,sss,isw,error,&hyplasmatid, &nmatData, NULL ); //7 int

	    return;

}

/*
    cout << '\t' <<  "      AAAAA       " << endl;
    cout << '\t' <<  "      hyplasmatid       " << '\t' << hyplasmatid << endl;
    cout << '\t' <<  "      *matid       " << '\t' << *matid << endl;
    cout << '\t' <<  "      nmatData       " << '\t' << nmatData << endl;
*/

	if(elm==NULL) prgError(2,"hyplasmaterial","hyplasmaterial called without elm!!!");

	//ElementGroup *eG = (ElementGroup*) elm->elemGrp;
	//thkgp     = eG->elemProp[ELEMENTTYPE]->data[3];
	iiter        = 0;
	//if (eG->dom->firstIter) iiter = 1;

       thkgp =1.0;
	hyplas_(matData,F,F33,stre,cc,iv1,iv2,&thkgp, //8 double
              finite,sss,isw,error,&hyplasmatid, &nmatData, &iiter ); //7 int

     double cc1[4][4];
     int ii, jj, ind=0;

        for(ii=0;ii<4;ii++)
        {
           for(jj=0;jj<4;jj++)
                cc1[ii][jj] =  cc[ind++];
        }

        cc[0] = cc1[0][0];
        cc[1] = cc1[0][1];
        cc[2] = cc1[0][3];
        cc[3] = cc1[0][2];

        cc[4] = cc1[1][0];
        cc[5] = cc1[1][1];
        cc[6] = cc1[1][3];
        cc[7] = cc1[1][2];

        cc[8] = cc1[3][0];
        cc[9] = cc1[3][1];
        cc[10] = cc1[3][3];
        cc[11] = cc1[3][2];

        cc[12] = cc1[2][0];
        cc[13] = cc1[2][1];
        cc[14] = cc1[2][3];
        cc[15] = cc1[2][2];

/*
        cc[0] = cc1[0][0];
        cc[1] = cc1[0][1];
        cc[2] = cc1[0][1];
        cc[3] = cc1[0][2];

        cc[4] = cc1[1][0];
        cc[5] = cc1[1][1];
        cc[6] = cc1[1][0];
        cc[7] = cc1[1][2];

        cc[8] = cc1[1][0];
        cc[9] = cc1[1][0];
        cc[10] = cc1[1][1];
        cc[11] = cc1[3][2];

        cc[12] = cc1[2][0];
        cc[13] = cc1[2][1];
        cc[14] = cc1[2][3];
        cc[15] = cc1[2][2];


        ind=0;
        cout << " Tangent Tensor  inside HYPLAS " << endl;
        for(ii=0;ii<4;ii++)
        {
           for(jj=0;jj<4;jj++)
           {
              printf("\t%8.4f",cc[ind++]);
              //printf("\t%8.4f",cc1[ii][jj]);
           }
           printf("\n");
        }
        printf("\n");
        printf("\n");
*/


	return;
}
