
#include "ImmersedIntegrationElement.h"
#include "BasisFunctionsLagrange.h"

#include "TimeFunction.h"
#include "MpapTime.h"
#include "Functions.h"
#include "SolutionData.h"


extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;





void ImmersedIntegrationElement::computeBodyForce(bool flag, double* VelSolid)
{
    int  ii, jj, nlocal, *degree, TI, endof, *bfstmp, ind;

    double   fact, fact2, ds, h1, h2, arclen=1.0, val1, val2, PENALTY_NUM, vel[2];
    double   *knots[2], *forceTmp, temp[2], normal[2];

    PENALTY_NUM = 1.0;

    degree = elem->getDegree();
    endof  = elem->getNdof();
    //endof  = 3;

    nlocal = (degree[0] + 1)*(degree[1] + 1);
    
    VectorXd  NN(nlocal), N;

    knots[0] = elem->getKnots(0);
    knots[1] = elem->getKnots(1);
    
    //GeomData->computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], param[0], param[1], &NN(0));

    //cout << " AAAAAAAAAAA " << endl;
    if(elem->getParent() == NULL)
      N = NN;
    else
      N = elem->SubDivMat*NN;
    
    vel[0] = elem->computeValue(0, N);
    vel[1] = elem->computeValue(1, N);
    
    //cout << " AAAAAAAAAAA " << endl;

    //U[0] += (mpapTime.dt*(vel[0]-VelSolid[0]));
    U[0] += (0.5*mpapTime.dt*(dU[0] + vel[0] - VelSolid[0]));
    //cout << " AAAAAAAAAAA " << endl;
    force[0] = PENALTY_NUM*(0.0 - U[0]);
    //force[0] = PENALTY_NUM*(positionCur[0] - positionOrig[0] - U[0]);


    //U[1] += (mpapTime.dt*(vel[1]-VelSolid[1]));
    U[1] += (0.5*mpapTime.dt*(dU[1] + vel[1] - VelSolid[1]));

    force[1] = PENALTY_NUM*(0.0 - U[1]); //case1
    //force[1] = PENALTY_NUM*(positionCur[1] - positionOrig[1] - U[1]); //case2
    //force[1] = PENALTY_NUM*(0.0 - U[1] + 1.0*(VelSolid[1]-vel[1]));
    //force[1] = PENALTY_NUM*(positionCur[1] - positionOrig[1] - U[1] + 1.0*(VelSolid[1]-vel[1]) );

    //printf("\t force \t %12.6f \t %12.6f \n", positionCur[1], positionOld[1]);
    //printf("\t param \t %12.6f \t %12.6f \n\n", param[0], param[1]);
    //printf("\t vel   \t %12.6f \t %12.6f \n\n", vel[0], vel[1]);
    //printf("\t disp  \t %12.6f \t %12.6f \n\n", U[0], U[1]);
    //printf("\t force \t %12.6f \t %12.6f \n\n", force[0], force[1]);

    dU[0] = vel[0];
    dU[1] = vel[1];
    //cout << " AAAAAAAAAAA " << endl;
    
    // map the force to global bodyForce vector


    //h1 = GeomData->getGridLength(0) * knots[0][2];
    //h2 = GeomData->getGridLength(1) * knots[1][2];

    arclen = 2.0*PI*0.5/301.0;
    //arclen = 2.0*PI*3.0/101.0;
    //arclen = 2.0/400.0;

    fact = arclen/h1/h2;

    //printf("\t posIndices %5d \t %5d \t %5d \n", id, posIndices[0], posIndices[1]);
    //printf("\t force \t %12.6f \t %12.6f \t %12.6f \n", force[0], force[1], arclen);

    /*
    fact2 = sqrt((param[0]-0.5)*(param[0]-0.5) + (param[1]-0.5)*(param[1]-0.5));
    
    normal[0] = (param[0]-0.5)/fact2;
    normal[1] = (param[1]-0.5)/fact2;

    force[0] = normal[0]*1.0;
    force[1] = normal[1]*1.0;
    */

    //printf("\t force \t %12.6f \t %12.6f \t %12.6f \n", force[0], force[1], arclen);
    //printf("\t h1 \t %12.6f \t %12.6f \t %12.6f \n", h1, h2, fact);

    forceTmp = &(SolnData->bodyForce(0));

    bfstmp = &(elem->GlobalBasisFuncs[0]);
    ind = elem->GlobalBasisFuncs.size();
    
    temp[0] = force[0]*fact;
    temp[1] = force[1]*fact;

    for(ii=0;ii<ind;ii++)
    {
       TI   =  endof*bfstmp[ii];

       for(jj=0;jj<DIM;jj++)
         forceTmp[TI+jj] += (N[ii] * temp[jj]);
    }


    return;
}



void ImmersedIntegrationElement::mapDataToGlobalBodyForceVector(bool flag, double* useThisData)
{
    // if (flag == true) use the values supplied in useThisData 
    // otherwise use the values in 'force'
    //

    int  ii, jj, nlocal, *degree, TI, endof, *bfstmp, ind;

    double   fact, fact2, ds, h1, h2, arclen=1.0, val1;
    double   *knots[2], *forceTmp, temp[2], normal[2];

    degree = elem->getDegree();
    endof  = elem->getNdof();

    nlocal = (degree[0] + 1)*(degree[1] + 1);
 
    VectorXd  NN(nlocal), N;

    knots[0] = elem->getKnots(0);
    knots[1] = elem->getKnots(1);

    //cout << " param " << param[0] << '\t' << param[1] << endl;

    //GeomData->computeBasisFunctions2D(knots[0][0], knots[1][0], knots[0][2], knots[1][2], param[0], param[1], &NN(0));

    //printf("BasisFuns \n");
    //for(ii=0;ii<nlocal;ii++)
      //printf(" \t %5d \t %12.6f \t %12.6f \t %12.6f \n ", ii, N[ii], dN_dx[ii], d2N_dx2[ii]);

    //cout << " LLLLLLLLLLLLLLLLLL " << endl;

    if(elem->getParent() == NULL )
      N = NN;
    else
      N = elem->SubDivMat * NN;

    //h1 = GeomData->getGridLength(0) * knots[0][2];
    //h2 = GeomData->getGridLength(1) * knots[1][2];

    arclen = 2.0*PI*0.5/201.0;
    arclen = 2.0*PI*3.0/101.0;
    //arclen = 2.0/400.0;

    fact = arclen/h1/h2;

    //printf("\t posIndices %5d \t %5d \t %5d \n", id, posIndices[0], posIndices[1]);
    //printf("\t force \t %12.6f \t %12.6f \t %12.6f \n", force[0], force[1], arclen);

    /*
    if(flag)
    {
       for(ii=0;ii<DIM;ii++)
         temp[ii] = useThisData[ii];
    }
    else
    {
       for(ii=0;ii<DIM;ii++)
         temp[ii] = force[ii];
    }
    */

    /*
    fact2 = sqrt((param[0]-0.5)*(param[0]-0.5) + (param[1]-0.5)*(param[1]-0.5));
    
    normal[0] = (param[0]-0.5)/fact2;
    normal[1] = (param[1]-0.5)/fact2;

    force[0] = normal[0]*1.0;
    force[1] = normal[1]*1.0;
    */

    //printf("\t force \t %12.6f \t %12.6f \t %12.6f \n", force[0], force[1], arclen);
    //printf("\t h1 \t %12.6f \t %12.6f \t %12.6f \n", h1, h2, fact);

    forceTmp = &(SolnData->bodyForce(0));

    bfstmp = &(elem->GlobalBasisFuncs[0]);
    ind = elem->GlobalBasisFuncs.size();
    
    temp[0] = force[0]*fact;
    temp[1] = force[1]*fact;

    for(ii=0;ii<ind;ii++)
    {
       TI   =  endof*bfstmp[ii];

       for(jj=0;jj<DIM;jj++)
         forceTmp[TI+jj] += (N[ii] * temp[jj]);
    }

  return;
}
//


