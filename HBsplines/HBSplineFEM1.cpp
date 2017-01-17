
#include "HBSplineFEM.h"
#include "DataBlockTemplate.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "PropertyTypeEnum.h"
#include "ImmersedIntegrationElement.h"
#include "ImmersedFlexibleSolid.h"
#include "ImmersedRigidSolid.h"
#include "ImmersedSolid.h"
#include "ContactElementPointToPoint2D.h"


extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern MpapTime           mpapTime;
extern ComputerTime       computerTime;


using namespace std;





HBSplineFEM::HBSplineFEM()
{
    stagParams.resize(3);
    stagParams[0] = 0;
    stagParams[1] = 1;
    stagParams[2] = 0;

    velDOF = presDOF = fluidDOF = IBDOF = solidDOF = 0;
    CREATE_POSTPROCESS_GRID = true;

    // add new type

    DomainType *hbsplinefem = domain.newType(HBSPLINEFEM, ROOTDOMAIN);

    if (hbsplinefem == NULL) return;  // domain type already exists

  hbsplinefem->key.addNew("origin", //0 
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
                          "solid elements",//28
                          "immersed point boundary condition",//29
                          "immersed body output", //30
                          "contact elements", //31
                          "element type", //32
                          "material type");

}


HBSplineFEM::~HBSplineFEM()
{

}




void  HBSplineFEM::readInputData(std::ifstream &Ifile, MyString &line)
{
  char fct[] = "HBSplineFEM::readInputData";

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


  switch (domain[HBSPLINEFEM].key.whichBegins(line))
  {
    case  0: cout << "     HBSplineFEM: reading 'origin' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
                prgError(1, fct, "invalid specification of 'origin' !");

             if(lvdTmp[0].n < 1)
                prgError(2, fct, "invalid number of 'origin' coordinates !");
                
             for (i=0; i<lvdTmp[0].n; i++)
                origin[i] = lvdTmp[0][i];

             break;

    case  1: cout << "     HBSplineFEM: reading 'grid dimensions' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
                prgError(1, fct, "invalid specification of 'grid dimensions' !");

             if(lvdTmp[0].n < 1)
                prgError(2, fct, "invalid number of 'grid dimensions' !");

             for (i=0; i<lvdTmp[0].n; i++)
                gridLEN[i] = lvdTmp[0][i];

             break;

    case  2: cout << "     HBSplineFEM: reading 'polynomial degrees' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'polynomial degrees'!");

            if(lviTmp[0].n < 1)
              prgError(2, fct, "invalid number of 'polynomial degrees' !");

            for (i=0; i<lviTmp[0].n; i++)
              degree[i] = lviTmp[0][i];

            break;

    case  3: cout << "     HBSplineFEM: reading 'number of elements' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'number of elements'!");

            if(lviTmp[0].n < 1)
               prgError(2, fct, "invalid number of 'number of elements' !");

            for (i=0; i<lviTmp[0].n; i++)
              nelem[i] = lviTmp[0][i];

            break;

    case  4: cout << "     HBSplineFEM: reading 'refinement type' ...\n\n";

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

    case  5: cout << "     HBSplineFEM: reading 'mesh refinement limits' ...\n\n";

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

    case  6: cout << "     HBSplineFEM: reading 'elements to refine' ...\n\n";

             if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
               prgError(1,fct,"invalid input in 'elements to refine'!");

             if(lviTmp[0].n < 1)
                prgError(2, fct, "invalid number of 'elements to refine' !");
                
             elemsToRefine.resize(lviTmp[0].n);

             for (i=0; i<lviTmp[0].n; i++)
                elemsToRefine[i] = lviTmp[0][i];

             break;

    case  7: cout << "     HBSplineFEM: reading 'elements to unRefine' ...\n\n";

             if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
               prgError(1,fct,"invalid input in 'elements to refine'!");

             if(lviTmp[0].n < 1)
                prgError(2, fct, "invalid number of 'elements to refine' !");
                
             elemsToUnRefine.resize(lviTmp[0].n);

             for (i=0; i<lviTmp[0].n; i++)
                elemsToUnRefine[i] = lviTmp[0][i];

             break;

    case  8: cout << "     HBSplineFEM: reading 'profile to refine' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
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

    case  9: cout << "     HBSplineFEM: reading 'dirichlet boundary conditions' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
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

    case  10: cout << "     HBSplineFEM: reading 'neumann boundary conditions' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
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

    case  11: cout << "     HBSplineFEM: reading 'point boundary conditions' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
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

    case  12: cout << "     HBSplineFEM: reading 'periodic boundary conditions' ...\n\n";

             if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
               prgError(1,fct,"invalid input in 'periodic boundary conditions'!");

            PERIODIC_BCS = true;
            
            break;

    case  13: cout << "     HBSplineFEM: reading 'initial conditions' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
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

    case  14: cout << "     HBSplineFEM: reading 'fluid properties' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'fluid properties'!");

            fluidProps.resize(lvdTmp[0].n);

            for(i=0;i<lvdTmp[0].n;i++)
              fluidProps[i] = lvdTmp[0][i];

            break;

    case  15: cout << "     HBSplineFEM: reading 'control parameters' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'control parameters'!");

            if( lvdTmp[0].n < 3)
              cerr <<  " Error in (( HBSplineFEM: reading 'control parameters' )) " << endl;

            tol      = lvdTmp[0][0];
            tis      = (int) lvdTmp[0][1];
            rhoInfty = lvdTmp[0][2];

            //cout << tol << '\t' << tis << '\t' << rhoInfty << endl;

            break;

    case  16: cout << "     HBSplineFEM: reading 'staggered or monolithic' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'immersed body data'!");

            //stagParams.resize(lvdTmp.n);

            for(i=0;i<lvdTmp[0].n;i++)
              stagParams[i] =  lvdTmp[0][i] ;

            STAGGERED = (lvdTmp[0][0] == 0);
            cout << " STAGGERED " << STAGGERED << endl;

            break;

    case  17: cout << "     HBSplineFEM: reading 'Galerkin or Leastsquares' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'Galerkin or Leastsquares'!");

            LSFEM_FLAG = (lviTmp[0][0] == 1);
            //cout << " LSFEM_FLAG " << LSFEM_FLAG << endl;

            break;

    case  18: cout << "     HBSplineFEM: reading 'analytical function' ...\n\n";

             line.getNextLine(Ifile);

             anlySolnType = line;
             
             cout << line << endl;

             line.getNextLine(Ifile);

             break;

    case  19: cout << "     HBSplineFEM: reading 'cutfem parameters' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'cutfem parameters'!");

            break;

    case  20: cout << "     HBSplineFEM: reading 'immersed body data' ...\n\n";

            ImmersedSolid  *imsolid;

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'immersed body data'!");

            if(lvdTmp.n != 1)
              prgError(2, fct, "error in 'immersed body data'! Only one immersed solid definition at a time");

            if(lvdTmp[0].n < 5)
              prgError(3, fct, "invalid number of data items 'immersed body data' !");


            if( (int) lvdTmp[0][0] == 0 )
            {
              imsolid = new ImmersedRigidSolid(DIM);
            }
            else
            {
              imsolid = new ImmersedFlexibleSolid(DIM);
            }

            imsolid->SetBoundaryConditionType( (int) lvdTmp[0][1] );
            imsolid->SetPenaltyParameter(lvdTmp[0][2]);
            imsolid->SetNitscheFlag(lvdTmp[0][3]);
            imsolid->SetNitscheFact(lvdTmp[0][4]);

            ImmersedBodyObjects.push_back(imsolid);

            break;

    case  21: cout << "     HBSplineFEM: reading 'immersed points' ...\n\n";

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

            ImmersedBodyObjects[bb]->SetNodalPositions(vecvecDbl);

            break;

    case  22: cout << "     HBSplineFEM: reading 'immersed integration elements' ...\n\n";

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

            ImmersedBodyObjects[bb]->SetImmersedIntegrationElements(vecvecInt);

            break;

    case  23: cout << "     HBSplineFEM: reading 'rigid body mass' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'rigid body mass'!");

            vecDbl.resize(lvdTmp[0].n);

            for(j=0; j<lvdTmp[0].n; j++)
              vecDbl[j] = lvdTmp[0][j];

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->SetMass(vecDbl);

            break;

    case  24: cout << "     HBSplineFEM: reading 'rigid body damping' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'rigid body damping'!");

            vecDbl.resize(lvdTmp[0].n);

            for(j=0; j<lvdTmp[0].n; j++)
              vecDbl[j] = lvdTmp[0][j];

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->SetDamping(vecDbl);

            break;

    case  25: cout << "     HBSplineFEM: reading 'rigid body stiffness' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'rigid body stiffness'!");

            vecDbl.resize(lvdTmp[0].n);

            for(j=0; j<lvdTmp[0].n; j++)
              vecDbl[j] = lvdTmp[0][j];

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->SetStiffness(vecDbl);

            break;

    case  26: cout << "     HBSplineFEM: reading 'rigid body degree of freedom' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'rigid body degree of freedom'!");

            vecInt.resize(lviTmp[0].n);

            for(j=0; j<lviTmp[0].n; j++)
              vecInt[j] = lviTmp[0][j];

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->SetBoundaryConditions(vecInt);

            break;

    case  27: cout << "     HBSplineFEM: reading 'rigid body prescribed motion' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'rigid body prescribed motion'!");

            //cerr <<  " Error in (( HBSplineFEM: reading 'rigid body prescribed motion' )) " << endl;

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

    case  28: cout << "     HBSplineFEM: reading 'solid elements' ...\n\n";

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

            //cout << vecvecInt[0] << endl;
            //cout << vecvecInt[1] << endl;

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;

            ImmersedBodyObjects[bb]->SetSolidElements(vecvecInt);

            break;

    case  29: cout << "     HBSplineFEM: reading 'immersed point boundary condition' ...\n\n";

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

            ImmersedBodyObjects[bb]->SetBoundaryConditions(vecvecDbl);

            break;


    case  30: cout << "     HBSplineFEM: reading 'immersed body output' ...\n\n";

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

            ImmersedBodyObjects[bb]->SetDataForOutput(vecvecInt);

            break;

    case  31: cout << "     HBSplineFEM: reading 'contact elements' ...\n\n";

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

    case  32: cout << "     HBSplineFEM: reading 'element type' ...\n\n";

            //ElemProp.add(new PropertyItem(ELEMENTTYPE));
            //ElemProp[ElemProp.n-1].readInputData(Ifile,line,"input error in 'element type'!");

            assert(ImmersedBodyObjects.size() > 0);

            bb = ImmersedBodyObjects.size() - 1;
            
            cout << " bb = " << bb << endl;

            ImmersedBodyObjects[bb]->SolnData.ElemProp.add(new PropertyItem(ELEMENTTYPE));
            ImmersedBodyObjects[bb]->SolnData.ElemProp[ImmersedBodyObjects[bb]->SolnData.ElemProp.n-1].readInputData(Ifile,line,"input error in 'element type'!");

            break;

    case  33: cout << "     HBSplineFEM: reading 'material type' ...\n\n";

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




void HBSplineFEM::prepareInputData()
{
    printf("\n     HBSplineFEM::prepareInputData()  .... STARTED ...\n");

    HBSplineBase::prepareInputData();

    for(int ii=0; ii<activeElements.size(); ii++)
    {
      elems[activeElements[ii]]->initialiseDOFvalues();
      //cout << " uuuuuuuuuuuuu " << endl;
    }


    ///////////////////////////////////////////////////
    //
    // prepare the data for immersed points
    // 
    ///////////////////////////////////////////////////

    int aa, bb, gp, r;

    ImmersedIntegrationElement *lme;

    for(bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
      for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
      {
        lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

        //cout << bb << '\t' << aa << '\t' << lme->IsActive() << '\t' << lme->gausspoints.size() << endl;

        lme->SolnData = &(SolnData);
        lme->GeomDataHBS = &(GeomData);

        if( lme->IsActive() )
        {
        //lme->elem.resize(lme->gausspoints.size());
        for(gp=0;gp<lme->gausspoints.size();gp++)
        {
          lme->computePointAtGP(gp, geom);
          r = findCellNumber(geom);

          lme->elemNums[gp] = r;
          //lme->elem[gp] = elems[r];
          //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %5d, \n", geom[0], geom[1], geom[2], r);
        }
        lme->elem = elems[lme->elemNums[0]];
        geometryToParametric(geom, lme->param);
        }
      }
    }


    velDOF  = gridBF1 * ndof;  // equal order for velocity and pressure
    presDOF = 0;

    //if( !LSFEM_FLAG && ndof > 1)
      //presDOF  = gridBF2;

    fluidDOF = velDOF + presDOF;

    IBDOF = 0;
    for(bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
      if(ImmersedBodyObjects[bb]->IsBoundaryConditionTypeLagrange())
        IBDOF += (ImmersedBodyObjects[bb]->GetNumNodes() * DIM);
    }

    solidDOF = 0;
    if(!STAGGERED)
    {
      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
        solidDOF += ImmersedBodyObjects[bb]->GetTotalDOF();

      solidDOF += contElemData.size();
    }
    //cout << " jjjjjjjjjjjjjjj  " << endl;

    totalDOF = velDOF + presDOF + IBDOF + solidDOF;


    printf("\n \t   Number of Velocity DOF    =  %5d\n\n", velDOF);
    printf("\n \t   Number of Pressure DOF    =  %5d\n\n", presDOF);
    printf("\n \t   Number of Immersed DOF    =  %5d\n\n", IBDOF);
    printf("\n \t   Number of Solid DOF       =  %5d\n\n", solidDOF);
    printf("\n \t   Total number of DOF       =  %5d\n\n", totalDOF);
    
    SolnData.initialise(velDOF, presDOF, IBDOF, solidDOF);

    SolnData.SetPhysicsTypetoFluid();
    SolnData.SetTimeIncrementType(tis);
    SolnData.SetRho(rhoInfty);
    SolnData.SetStaggeredParams(stagParams);

    for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      ImmersedBodyObjects[bb]->SolnData.stagParams = stagParams;

    soln.resize(totalDOF);
    soln.setZero();
    solnInit = soln;

    //matxGrid.free();
    //MatrixSparse<double> matxtemp;
    //prepareMatrixPatternGrid(matxtemp);
    //matxGrid = matxtemp;

    GRID_CHANGED = true;
    IB_MOVED = false;

    if(IBDOF > 0)
    {
      IB_MOVED = true;
    }

    printf("     HBSplineFEM::prepareInputData()  .... FINISHED ...\n\n");

    return;
}




void HBSplineFEM::printInfo()
{
  return;
}




