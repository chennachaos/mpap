
#include "HBSplineCutFEM.h"
#include "MathBasic.h"
#include "DataBlockTemplate.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "PropertyTypeEnum.h"
#include "ImmersedFlexibleSolid.h"
#include "ImmersedRigidSolid.h"
#include "ImmersedSolid.h"
#include "ContactElementPointToPoint2D.h"


extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern MpapTime           mpapTime;
extern ComputerTime       computerTime;


using namespace std;





HBSplineCutFEM::HBSplineCutFEM()
{
  CUTCELL_INTEGRATION_TYPE = 1;
  
  slnTemp.resize(3);
  slnTemp.setZero();

  slnTempPrev2 = slnTemp;
  slnTempPrev  = slnTemp;
  slnTempCur   = slnTemp;

  stagParams.resize(10);
  stagParams[0] = 0;
  stagParams[1] = 1;
  stagParams[2] = 0;
  stagParams[3] = 0;
  stagParams[4] = 0;

    // add new type

  DomainType *hbsplinecutfem = domain.newType(HBSPLINECUTFEM, ROOTDOMAIN);

  if(hbsplinecutfem == NULL)
    return;  // domain type already exists

  hbsplinecutfem->key.addNew("origin", //0 
                          "grid dimensions", //1
                          "polynomial degrees", //2
                          "number of elements", //3
                          "refinement type", //4
                          "mesh refinement limits",//5
                          "elements to refine", //6
                          "elements to unrefine", //7
                          "profile to refine",//8
                          "dirichlet boundary conditions", //9
                          "neumann boundary conditions", //10
                          "point boundary conditions", //11
                          "periodic boundary conditions", //12
                          "initial conditions", //13
                          "fluid properties", //14
                          "control parameters", //15
                          "staggered or monolithic", //16
                          "Galerkin or Leastsquares", //17
                          "analytical function", //18
                          "cutfem parameters", //19
                          "immersed body data", //20
                          "immersed points", //21
                          "immersed integration elements",//22
                          "rigid body mass", //23
                          "rigid body damping",//24
                          "rigid body stiffness",//25
                          "rigid body degree of freedom", //26
                          "rigid body prescribed motion", //27
                          "rigid body preload", // 28
                          "rigid body initial force predictor", // 29
                          "rigid body motion limits", //30
                          "solid elements",//31
                          "immersed point boundary condition",//32
                          "immersed body output", //33
                          "contact elements", //34
                          "element type", //35
                          "material type");

}


HBSplineCutFEM::~HBSplineCutFEM()
{
  //std::for_each(elems.begin(), elems.end(), delete_pointed_to<elems>);
  //elems.erase(elems.begin(), elems.end());
  //cout << " elems.size() " << elems.size() << endl;

}




void  HBSplineCutFEM::readInputData(std::ifstream &Ifile, MyString &line)
{
  char fct[] = "HBSplineCutFEM::readInputData";

  MyString tmpl, *word;

  char tmp[30];

  int nw, i, j, k, n, nn, ii, bb;
  
  double fact;

  MyStringList   sTmp;
  List<Vector<int> > lviTmp;
  List<Vector<double> > lvdTmp;

  DataBlockTemplate t1, t2;

  vector<vector<double> >  vecvecDbl;

  vector<vector<int> >   vecvecInt;

  vector<double>  vecDbl;
  vector<int>  vecInt;


  switch (domain[HBSPLINECUTFEM].key.whichBegins(line))
  {
    case  0: //cout << "     HBSplineCutFEM: reading 'origin' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
                prgError(1, fct, "invalid specification of 'origin' !");

             if(lvdTmp[0].n < 1)
                prgError(2, fct, "invalid number of 'origin' coordinates !");
                
             for (i=0; i<lvdTmp[0].n; i++)
                origin[i] = lvdTmp[0][i];

             break;

    case  1: //cout << "     HBSplineCutFEM: reading 'grid dimensions' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
                prgError(1, fct, "invalid specification of 'grid dimensions' !");

             if(lvdTmp[0].n < 1)
                prgError(2, fct, "invalid number of 'grid dimensions' !");

             for (i=0; i<lvdTmp[0].n; i++)
                gridLEN[i] = lvdTmp[0][i];

             break;

    case  2: //cout << "     HBSplineCutFEM: reading 'polynomial degrees' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'polynomial degrees'!");

            if(lviTmp[0].n < 1)
              prgError(2, fct, "invalid number of 'polynomial degrees' !");

            for (i=0; i<lviTmp[0].n; i++)
              degree[i] = lviTmp[0][i];

            break;

    case  3: //cout << "     HBSplineCutFEM: reading 'number of elements' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'number of elements'!");

            if(lviTmp[0].n < 1)
               prgError(2, fct, "invalid number of 'number of elements' !");

            for (i=0; i<lviTmp[0].n; i++)
              nelem[i] = lviTmp[0][i];

            break;

    case  4: //cout << "     HBSplineCutFEM: reading 'refinement type' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'refinement type'!");

            refinementData.resize(lvdTmp.n);

            //cout << lvdTmp[0] << endl;
            refinementData.resize(lvdTmp[0].n);
            if(lvdTmp[0].n < 1)
              prgError(2, fct, "invalid number of 'refinement type' !");

            for(j=0;j<lvdTmp[0].n;j++)
              refinementData[j] = lvdTmp[0][j];

            //printVector(refinementData);

            break;

    case  5: //cout << "     HBSplineCutFEM: reading 'mesh refinement limits' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'mesh refinement limits'!");

            refineLimitVals.resize(lvdTmp.n);

            for(i=0;i<lvdTmp.n;i++)
            {
                //cout << lvdTmp[i] << endl;
                refineLimitVals[i].resize(lvdTmp[i].n);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'mesh refinement limits' !");
                
                for(j=0;j<lvdTmp[i].n;j++)
                  refineLimitVals[i][j] = lvdTmp[i][j];
            }

            //printVector(&refineLimitVals[0][0], refineLimitVals[0].size());
            //printVector(&refineLimitVals[1][0], refineLimitVals[1].size());

            break;

    case  6: //cout << "     HBSplineCutFEM: reading 'elements to refine' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'elements to refine'!");

            if(lviTmp[0].n < 1)
              prgError(2, fct, "invalid number of 'elements to refine' !");

            elemsToRefine.resize(lviTmp[0].n);

            for (i=0; i<lviTmp[0].n; i++)
              elemsToRefine[i] = lviTmp[0][i];

            break;

    case  7: //cout << "     HBSplineCutFEM: reading 'elements to unRefine' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'elements to refine'!");

            if(lviTmp[0].n < 1)
              prgError(2, fct, "invalid number of 'elements to refine' !");

            elemsToUnRefine.resize(lviTmp[0].n);

            for (i=0; i<lviTmp[0].n; i++)
              elemsToUnRefine[i] = lviTmp[0][i];

            break;

    case  8: //cout << "     HBSplineCutFEM: reading 'profile to refine' ...\n\n";

            if(!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'profile to refine'!");

            /*
            LevelSetFunc.resize(lvdTmp.n);

            for(i=0;i<lvdTmp.n;i++)
            {
                LevelSetFunc[i].resize(lvdTmp[i].n);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'profile to refine' !");
                
                for(j=0;j<lvdTmp[i].n;j++)
                  LevelSetFunc[i][j] = lvdTmp[i][j];
            }
            */

            break;

    case  9: //cout << "     HBSplineCutFEM: reading 'dirichlet boundary conditions' ...\n\n";

            if(!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'dirichlet boundary conditions'!");

            DirichletBCs.resize(lvdTmp.n);

            for(i=0;i<lvdTmp.n;i++)
            {
                DirichletBCs[i].resize(lvdTmp[i].n);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'dirichlet boundary conditions' !");
                
                for(j=0;j<lvdTmp[i].n;j++)
                  DirichletBCs[i][j] = lvdTmp[i][j];
            }

            break;

    case  10: //cout << "     HBSplineCutFEM: reading 'neumann boundary conditions' ...\n\n";

            if(!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'neumann boundary conditions'!");

            NeumannBCs.resize(lvdTmp.n);

            for(i=0;i<lvdTmp.n;i++)
            {
                NeumannBCs[i].resize(lvdTmp[i].n);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'neumann boundary conditions' !");
                
                for(j=0;j<lvdTmp[i].n;j++)
                  NeumannBCs[i][j] = lvdTmp[i][j];
            }

            break;

    case  11: //cout << "     HBSplineCutFEM: reading 'point boundary conditions' ...\n\n";

            if(!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'point boundary conditions'!");

            pointBCs.resize(lvdTmp.n);

            for(i=0;i<lvdTmp.n;i++)
            {
                pointBCs[i].resize(lvdTmp[i].n);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'point boundary conditions' !");
                
                for(j=0;j<lvdTmp[i].n;j++)
                  pointBCs[i][j] = lvdTmp[i][j];
            }

            break;

    case  12: //cout << "     HBSplineCutFEM: reading 'periodic boundary conditions' ...\n\n";

            if(!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'periodic boundary conditions'!");

            PERIODIC_BCS = true;

            break;

    case  13: //cout << "     HBSplineCutFEM: reading 'initial conditions' ...\n\n";

            if(!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'initial conditions'!");

            Iconds.resize(lvdTmp.n);

            for(i=0;i<lvdTmp.n;i++)
            {
                Iconds[i].resize(lvdTmp[i].n);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'initial conditions' !");

                for(j=0;j<lvdTmp[i].n;j++)
                  Iconds[i][j] = lvdTmp[i][j];
            }

            break;

    case  14: //cout << "     HBSplineCutFEM: reading 'fluid properties' ...\n\n";

            if(!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'fluid properties'!");

            fluidProps.resize(lvdTmp[0].n);

            for(i=0;i<lvdTmp[0].n;i++)
              fluidProps[i] = lvdTmp[0][i];

            break;

    case  15: //cout << "     HBSplineCutFEM: reading 'control parameters' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'control parameters'!");

            if( lvdTmp[0].n < 3)
              cerr <<  " Error in (( HBSplineCutFEM: reading 'control parameters' )) " << endl;

            tol      = lvdTmp[0][0];
            tis      = (int) lvdTmp[0][1];
            rhoInfty = lvdTmp[0][2];

            //cout << tol << '\t' << tis << '\t' << rhoInfty << endl;

            break;

    case  16: //cout << "     HBSplineCutFEM: reading 'staggered or monolithic' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'immersed body data'!");

            //stagParams.resize(lvdTmp.n);

            for(i=0;i<lvdTmp[0].n;i++)
              stagParams[i] =  lvdTmp[0][i] ;

            STAGGERED = (lvdTmp[0][0] == 0);
            cout << " STAGGERED " << STAGGERED << endl;

            break;

    case  17: cout << "     HBSplineCutFEM: reading 'Galerkin or Leastsquares' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'Galerkin or Leastsquares'!");

            LSFEM_FLAG = (lviTmp[0][0] == 1);
            //cout << " LSFEM_FLAG " << LSFEM_FLAG << endl;

            break;

    case  18: cout << "     HBSplineCutFEM: reading 'analytical function' ...\n\n";

            line.getNextLine(Ifile);

            anlySolnType = line;

            cout << line << endl;

            line.getNextLine(Ifile);

            break;

    case  19: //cout << "     HBSplineCutFEM: reading 'cutfem parameters' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'cutfem parameters'!");

            cutFEMparams.resize(lvdTmp[0].n);

            for(i=0;i<lvdTmp[0].n;i++)
              cutFEMparams[i] = lvdTmp[0][i];

            break;

    case  20: //cout << "     HBSplineCutFEM: reading 'immersed body data' ...\n\n";

            ImmersedSolid  *imsolid;

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'immersed body data'!");

            if(lvdTmp.n != 1)
              prgError(2, fct, "error in 'immersed body data'! Only one immersed solid definition at a time");

            if(lvdTmp[0].n < 5)
              prgError(3, fct, "invalid number of data items 'immersed body data' !");


            if( lvdTmp[0][0] == 0 )
            {
              imsolid = new ImmersedRigidSolid(DIM);
            }
            else
            {
              imsolid = new ImmersedFlexibleSolid(DIM);
            }

            imsolid->setBoundaryConditionType(1);
            imsolid->setPenaltyParameter(lvdTmp[0][2]);
            imsolid->setNitscheFlag(lvdTmp[0][3]);
            imsolid->setNitscheFact(lvdTmp[0][4]);

            ImmersedBodyObjects.push_back(imsolid);

            break;

    case  21: //cout << "     HBSplineCutFEM: reading 'immersed points' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'immersed points'!");

            vecvecDbl.resize(lvdTmp.n);

            for(i=0;i<lvdTmp.n;i++)
            {
                //cout << lvdTmp[i] << endl;
                vecvecDbl[i].resize(lvdTmp[i].n-1);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'immersed points' !");
                
                for(j=0;j<lvdTmp[i].n-1;j++)
                  vecvecDbl[i][j] = lvdTmp[i][j+1];
            }
            //for(i=0;i<vecvecDbl.size();i++)
              //printVector(&(vecvecDbl[i][0]), vecvecDbl[i].size());

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->setNodalPositions(vecvecDbl);

            break;

    case  22: //cout << "     HBSplineCutFEM: reading 'immersed integration elements' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'immersed integration elements'!");

            vecvecInt.resize(lviTmp.n);

            k = lviTmp[0].n - 1 ;

            for(i=0;i<lviTmp.n;i++)
            {
              if(lviTmp[i].n < 1)
                prgError(2, fct, "invalid number of 'immersed integration elements' !");

              vecvecInt[i].resize(k);

              vecvecInt[i][0] = lviTmp[i][1];      // whether the element is active or inactive

                for(j=2; j<lviTmp[i].n; j++)
                  vecvecInt[i][j-1] = lviTmp[i][j] - 1;
            }

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->setImmersedIntegrationElements(vecvecInt);

            break;

    case  23: //cout << "     HBSplineCutFEM: reading 'rigid body mass' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'rigid body mass'!");

            vecDbl.resize(lvdTmp[0].n);

            for(j=0; j<lvdTmp[0].n; j++)
              vecDbl[j] = lvdTmp[0][j];

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->setMass(vecDbl);

            break;

    case  24: //cout << "     HBSplineCutFEM: reading 'rigid body damping' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'rigid body damping'!");

            vecDbl.resize(lvdTmp[0].n);

            for(j=0; j<lvdTmp[0].n; j++)
              vecDbl[j] = lvdTmp[0][j];

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->setDamping(vecDbl);

            break;

    case  25: //cout << "     HBSplineCutFEM: reading 'rigid body stiffness' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'rigid body stiffness'!");

            vecDbl.resize(lvdTmp[0].n);

            for(j=0; j<lvdTmp[0].n; j++)
              vecDbl[j] = lvdTmp[0][j];

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->setStiffness(vecDbl);

            break;

    case  26: //cout << "     HBSplineCutFEM: reading 'rigid body degree of freedom' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'rigid body degree of freedom'!");

            vecInt.resize(lviTmp[0].n);

            for(j=0; j<lviTmp[0].n; j++)
              vecInt[j] = lviTmp[0][j];

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->setBoundaryConditions(vecInt);

            break;

    case  27: //cout << "     HBSplineCutFEM: reading 'rigid body prescribed motion' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'rigid body prescribed motion'!");

            //cerr <<  " Error in (( HBSplineCutFEM: reading 'rigid body prescribed motion' )) " << endl;

            vecvecInt.resize(lviTmp.n);

            for(i=0;i<lviTmp.n;i++)
            {
                vecvecInt[i].resize(lviTmp[i].n);
                if(lviTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'rigid body prescribed motion' !");
                
                for(j=0;j<lviTmp[i].n;j++)
                  vecvecInt[i][j] = lviTmp[i][j];
            }

            break;

    case  28: //cout << "     HBSplineCutFEM: reading 'rigid body preload' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'rigid body preload'!");

            vecDbl.resize(lvdTmp[0].n);

            for(j=0; j<lvdTmp[0].n; j++)
              vecDbl[j] = lvdTmp[0][j];

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->setPreload(vecDbl);

            break;

    case  29: //cout << "     HBSplineCutFEM: reading 'rigid body initial force predictor' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'rigid body initial force predictor'!");

            vecDbl.resize(lvdTmp[0].n);

            for(j=0; j<lvdTmp[0].n; j++)
              vecDbl[j] = lvdTmp[0][j];

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->setInitialForcePredictor(vecDbl);

            break;

    case  30: //cout << "     HBSplineCutFEM: reading 'rigid body motion limits' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'rigid body motion limits'!");

            vecvecDbl.resize(lvdTmp.n);

            for(i=0;i<lvdTmp.n;i++)
            {
                vecvecDbl[i].resize(lvdTmp[i].n);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'rigid body motion limits' !");
                
                for(j=0;j<lvdTmp[i].n;j++)
                  vecvecDbl[i][j] = lvdTmp[i][j];
            }

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->setRigidBodyMotionLimits(vecvecDbl);

            break;

    case  31: //cout << "     HBSplineCutFEM: reading 'solid elements' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'solid elements'!");

            vecvecInt.resize(lviTmp.n);

            for(i=0;i<lviTmp.n;i++)
            {
                //cout << lviTmp[i] << endl;
                vecvecInt[i].resize(lviTmp[i].n-1);
                if(lviTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'solid elements' !");
                
                for(j=1;j<lviTmp[i].n;j++)
                  vecvecInt[i][j-1] = lviTmp[i][j] - 1;
            }

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->setSolidElements(vecvecInt);

            break;

    case  32: //cout << "     HBSplineCutFEM: reading 'immersed point boundary condition' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'immersed point boundary condition'!");

            vecvecDbl.resize(lvdTmp.n);

            for(i=0;i<lvdTmp.n;i++)
            {
                //cout << lvdTmpTmp[i] << endl;
                vecvecDbl[i].resize(lvdTmp[i].n);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'immersed point boundary condition' !");
                
                for(j=0;j<lvdTmp[i].n;j++)
                  vecvecDbl[i][j] = lvdTmp[i][j];
            }

            //printVector(&vecvecDbl[0][0], vecvecDbl[0].size());
            //printVector(&vecvecDbl[1][0], vecvecDbl[1].size());

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->setBoundaryConditions(vecvecDbl);

            break;

    case  33: //cout << "     HBSplineCutFEM: reading 'immersed body output' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'immersed body output'!");

            vecvecInt.resize(lviTmp.n);

            for(i=0; i<lviTmp.n; i++)
            {
                vecvecInt[i].resize(lviTmp[i].n);
                if(lviTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'immersed body output' !");
                
                for(j=0;j<lviTmp[i].n;j++)
                  vecvecInt[i][j] = lviTmp[i][j];
            }

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->setDataForOutput(vecvecInt);

            break;

    case  34: //cout << "     HBSplineCutFEM: reading 'contact elements' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'contact elements'!");

            contElemData.resize(lviTmp.n);

            for(i=0;i<lviTmp.n;i++)
            {
                contElemData[i].resize(lviTmp[i].n);
                if(lviTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'contact elements' !");
                
                for(j=0;j<lviTmp[i].n;j++)
                  contElemData[i][j] = lviTmp[i][j];
            }

            break;

    case  35: //cout << "     HBSplineCutFEM: reading 'element type' ...\n\n";

            //ElemProp.add(new PropertyItem(ELEMENTTYPE));
            //ElemProp[ElemProp.n-1].readInputData(Ifile,line,"input error in 'element type'!");

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->SolnData.ElemProp.add(new PropertyItem(ELEMENTTYPE));
            ImmersedBodyObjects[bb]->SolnData.ElemProp[ImmersedBodyObjects[bb]->SolnData.ElemProp.n-1].readInputData(Ifile,line,"input error in 'element type'!");

            break;

    case  36: //cout << "     HBSplineCutFEM: reading 'material type' ...\n\n";

            //MatlProp.add(new PropertyItem(MATERIAL));
            //MatlProp[MatlProp.n-1].readInputData(Ifile,line,"input error in 'material type'!");

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->SolnData.MatlProp.add(new PropertyItem(ELEMENTTYPE));
            ImmersedBodyObjects[bb]->SolnData.MatlProp[ImmersedBodyObjects[bb]->SolnData.MatlProp.n-1].readInputData(Ifile,line,"input error in 'element type'!");

            break;

    case -1: // go and inherit from DOMAIN

            this->Domain::readInputData(Ifile,line);

            DIM = ndm;
            ndof = ndf;

            break;
  }

  return;
}




void HBSplineCutFEM::prepareInputData()
{
  //printf("\n     HBSplineCutFEM::prepareInputData()  .... STARTED ...\n");

  HBSplineBase::prepareInputData();

  grid_to_cutfem_BF.assign(gridBF1, -1);
  grid_to_cutfem_BFprev.assign(gridBF1, -1);
  grid_to_proc_BF.assign(gridBF1, -1);

  totalDOF  = gridBF1 * ndof;
  
  grid_to_cutfem_DOF.assign(totalDOF, -1);
  grid_to_cutfem_DOFprev.assign(totalDOF, -1);
  grid_to_proc_DOF.assign(totalDOF, -1);

  SolnData.initialise(totalDOF, totalDOF, 0, 0);

  SolnData.setPhysicsTypetoFluid();
  SolnData.setTimeIncrementType(tis);
  SolnData.setSpectralRadius(td[0]);
  SolnData.SetStaggeredParams(stagParams);

  int bb;

  GRID_CHANGED = true;
  IB_MOVED = false;

  numDomains = ImmersedBodyObjects.size() + 1;

  // 0 -> the domain is inactive
  // 1 -> the domain is active

  //for(int bb=0; bb<numDomains; bb++)
    //domainInclYesNo[bb] = 0;
  
  //domainInclYesNo[0] = 1;

  //GeomData.domainInclYesNo = domainInclYesNo;


  if(nImmSolids > 0)
  {
    GeomData.cutFEMparams = cutFEMparams;
    
    CUTCELL_INTEGRATION_TYPE = cutFEMparams[0];

    double minVal[]={0.0, 0.0, 0.0}, maxVal[]={0.0, 0.0, 0.0};

    for(int ii=0; ii<DIM; ii++)
    {
      minVal[ii] = GeomData.computeCoord(ii, 0.0);
      maxVal[ii] = GeomData.computeCoord(ii, 1.0);
    }
    //cout << minVal[0] << '\t' << minVal[1] << '\t' << minVal[2] << endl;
    //cout << maxVal[0] << '\t' << maxVal[1] << '\t' << maxVal[2] << endl;

    for(bb=0; bb<nImmSolids; bb++)
      ImmersedBodyObjects[bb]->adjustBoundaryPoints(minVal, maxVal);
  }


  /*
  // subroutines using VTK library
  ////////////////////////////////////////////////////////////////

  // create the points on the background grids
  // for determining whether they lie inside/outside the immersed bodies
  // Also useful for post-processing purposes.
  // Using VTK datastructure only unique points are stored to save memory and 
  // time taken for the inside/outside query

  if(DIM == 3)
  {
    int ee, ll;

    vtkIdType  ptId;

    node* nd;

    AABB  bbTemp;

    ee = nelem[0]*nelem[1]*nelem[2];
    
    if(MAX_LEVEL > 0)
      ee *= 2;
    
    double  bounds[6],  ptTemp[8][3];

    bounds[0] = computeGeometry(0, 0.0);
    bounds[1] = computeGeometry(0, 1.0);

    bounds[2] = computeGeometry(1, 0.0);
    bounds[3] = computeGeometry(1, 1.0);

    bounds[4] = computeGeometry(2, 0.0);
    bounds[5] = computeGeometry(2, 1.0);

    //cout << " cccccccccc " << endl;

    mergePoints->InitPointInsertion(pointsVTKfluidgrid, bounds);
    mergePoints->SetDivisions(nelem[0], nelem[1], nelem[2]);
    mergePoints->setTolerance(1.0e-5);

    //cout << " AAAAAAAAAAA " << endl;

    for(ee=0; ee<activeElements.size(); ee++)
    {
      nd = elems[activeElements[ee]];

      bbTemp = nd->getAABB();
      
      //cout <<  " ee = " << ee << '\t' << activeElements[ee] << endl;

      ptTemp[0][0] = bbTemp.minBB[0];    ptTemp[0][1] = bbTemp.minBB[1];    ptTemp[0][2] = bbTemp.minBB[2];
      ptTemp[1][0] = bbTemp.maxBB[0];    ptTemp[1][1] = bbTemp.minBB[1];    ptTemp[1][2] = bbTemp.minBB[2];
      ptTemp[2][0] = bbTemp.maxBB[0];    ptTemp[2][1] = bbTemp.maxBB[1];    ptTemp[2][2] = bbTemp.minBB[2];
      ptTemp[3][0] = bbTemp.minBB[0];    ptTemp[3][1] = bbTemp.maxBB[1];    ptTemp[3][2] = bbTemp.minBB[2];

      ptTemp[4][0] = bbTemp.minBB[0];    ptTemp[4][1] = bbTemp.minBB[1];    ptTemp[4][2] = bbTemp.maxBB[2];
      ptTemp[5][0] = bbTemp.maxBB[0];    ptTemp[5][1] = bbTemp.minBB[1];    ptTemp[5][2] = bbTemp.maxBB[2];
      ptTemp[6][0] = bbTemp.maxBB[0];    ptTemp[6][1] = bbTemp.maxBB[1];    ptTemp[6][2] = bbTemp.maxBB[2];
      ptTemp[7][0] = bbTemp.minBB[0];    ptTemp[7][1] = bbTemp.maxBB[1];    ptTemp[7][2] = bbTemp.maxBB[2];

      for(ll=0;ll<8;ll++)
      {
        mergePoints->InsertUniquePoint( ptTemp[ll], ptId );
        hexVTK->GetPointIds()->SetId(ll, ptId);
      }

      uGridVTKfluid->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
    }

    //cout << " AAAAAAAAAAA " << endl;

    pointsPolydata->SetPoints(pointsVTKfluidgrid);

    uGridVTKfluid->SetPoints(pointsVTKfluidgrid);
    //uGridVTKfluid->SetPoints(newPts);

    //cout << " uGridVTKfluid->GetNumberOfPoints() = " <<  uGridVTKfluid->GetNumberOfPoints() << endl;
    //cout << " uGridVTKfluid->GetNumberOfCells()  = " <<  uGridVTKfluid->GetNumberOfCells() << endl;
  }

  ////////////////////////////////////////
  // generate polygons for the immersed boundaries
  //
  ////////////////////////////////////////

  vtkSmartPointer<vtkPolyData> pointsPolydataVTK   =   vtkSmartPointer<vtkPolyData>::New();
  pointsPolydataVTK->SetPoints(pointsVTKfluidgrid);

  for(bb=0; bb<nImmSolids; bb++)
  {
    #if VTK_MAJOR_VERSION == 5
      ImmersedBodyObjects[bb]->selectEnclosedPoints->SetInput(pointsPolydataVTK);
    #else
      ImmersedBodyObjects[bb]->selectEnclosedPoints->SetInputData(pointsPolydataVTK);
    #endif

    //cout << "  HBSplineCutFEM::setImmersedFaces()  " << endl;
    ImmersedBodyObjects[bb]->setImmersedFaces();
  }

  */

  // subroutines using CGAL library
  ////////////////////////////////////////////////////////////////

  for(bb=0; bb<nImmSolids; bb++)
  {
    //cout << "  HBSplineCutFEM::setImmersedFaces()  " << endl;
    ImmersedBodyObjects[bb]->setImmersedFaces();
  }

  PetscPrintf(MPI_COMM_WORLD, "\n     HBSplineCutFEM::prepareInputData()  .... FINISHED ...\n\n");

  return;
}




void HBSplineCutFEM::printInfo()
{
  return;
}




