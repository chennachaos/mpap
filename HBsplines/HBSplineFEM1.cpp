
#include "HBSplineFEM.h"

//#include "headersBasic.h"

#include "MathBasic.h"
#include "DataBlockTemplate.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "PlotVTK.h"
#include "PropertyTypeEnum.h"
#include "ImmersedIntegrationElement.h"
#include "ImmersedFlexibleSolid.h"
#include "ImmersedRigidSolid.h"
#include "ImmersedSolid.h"
#include "ContactElementPointToPoint2D.h"
#include "NurbsUtilities.h"

#include <omp.h>

extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern MpapTime           mpapTime;
extern ComputerTime       computerTime;
extern PlotVTK  plotvtk;


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
                          "boundary knots", //4
                          "refinement type", //5
                          "mesh refinement limits",//6
                          "elements to refine", //7
                          "elements to unrefine", //8
                          "element type", //9
                          "material type", //10
                          "dirichlet boundary conditions", //11
                          "neumann boundary conditions", //12
                          "traction bc data", //13
                          "point boundary conditions", //14
                          "periodic boundary conditions", //15
                          "initial conditions", //16
                          "control parameters", //17
                          "profile to refine",//18
                          "staggered or monolithic", //19
                          "Galerkin or Leastsquares", //20
                          "immersed body data", //21
                          "immersed points", //22
                          "immersed integration elements",//23
                          "solid elements",//24
                          "immersed body immersed intergration elements connectivity", //25
                          "immersed body solid element connectivity", //26
                          "immersed point boundary condition",//27
                          "immersed body output", //28
                          "analytical function", //29
                          "rigid body mass", //30
                          "rigid body damping",//31
                          "rigid body stiffness",//32
                          "rigid body degree of freedom", //33
                          "fluid properties", //34
                          "contact elements", //35
                          "rigid body prescribed motion");
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

    case  4: cout << "     HBSplineFEM: reading 'boundary knots' ...\n\n";
    
             prgError(2, fct, "nothing to be read in 'boundary knots' !");

             break;


    case  5: cout << "     HBSplineFEM: reading 'refinement type' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
               prgError(1,fct,"invalid input in 'refinement type'!");
             
             refinementData.resize(lvdTmp.n);

             cout << lvdTmp[0] << endl;
             refinementData.resize(lvdTmp[0].n);
             if(lvdTmp[0].n < 1)
               prgError(2, fct, "invalid number of 'refinement type' !");
                
             for(j=0;j<lvdTmp[0].n;j++)
               refinementData[j] = lvdTmp[0][j];

             printVector(refinementData);

             break;

    case  6: cout << "     HBSplineFEM: reading 'mesh refinement limits' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
               prgError(1,fct,"invalid input in 'mesh refinement limits'!");
             
             refineLimitVals.resize(lvdTmp.n);

             for(i=0;i<lvdTmp.n;i++)
             {
                cout << lvdTmp[i] << endl;
                refineLimitVals[i].resize(lvdTmp[i].n);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'mesh refinement limits' !");
                
                for(j=0;j<lvdTmp[i].n;j++)
                  refineLimitVals[i][j] = lvdTmp[i][j];
             }
             
             printVector(&refineLimitVals[0][0], refineLimitVals[0].size());
             //printVector(&refineLimitVals[1][0], refineLimitVals[1].size());
             
             break;

    case  7: cout << "     HBSplineFEM: reading 'elements to refine' ...\n\n";

             if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
               prgError(1,fct,"invalid input in 'elements to refine'!");

             if(lviTmp[0].n < 1)
                prgError(2, fct, "invalid number of 'elements to refine' !");
                
             elemsToRefine.resize(lviTmp[0].n);

             for (i=0; i<lviTmp[0].n; i++)
                elemsToRefine[i] = lviTmp[0][i];

             break;

    case  8: cout << "     HBSplineFEM: reading 'elements to unRefine' ...\n\n";

             if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
               prgError(1,fct,"invalid input in 'elements to refine'!");

             if(lviTmp[0].n < 1)
                prgError(2, fct, "invalid number of 'elements to refine' !");
                
             elemsToUnRefine.resize(lviTmp[0].n);

             for (i=0; i<lviTmp[0].n; i++)
                elemsToUnRefine[i] = lviTmp[0][i];

             break;


    case  9: cout << "     HBSplineFEM: reading 'element type' ...\n\n";

             ElemProp.add(new PropertyItem(ELEMENTTYPE));
             ElemProp[ElemProp.n-1].readInputData(Ifile,line,"input error in 'element type'!");

             break;

    case  10: cout << "     HBSplineFEM: reading 'material type' ...\n\n";

             MatlProp.add(new PropertyItem(MATERIAL));
             MatlProp[MatlProp.n-1].readInputData(Ifile,line,"input error in 'material type'!");

             break;

    case  11: cout << "     HBSplineFEM: reading 'dirichlet boundary conditions' ...\n\n";

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

    case  12: cout << "     HBSplineFEM: reading 'neumann boundary conditions' ...\n\n";

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

    case  13: cout << "     HBSplineFEM: reading 'traction bc data' ...\n\n";
    
             prgError(2, fct, "nothing to be read in 'traction bc data' !");

             break;

    case  14: cout << "     HBSplineFEM: reading 'point boundary conditions' ...\n\n";

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

    case  15: cout << "     HBSplineFEM: reading 'periodic boundary conditions' ...\n\n";

             if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
               prgError(1,fct,"invalid input in 'periodic boundary conditions'!");

            PERIODIC_BCS = true;
            
            break;

    case  16: cout << "     HBSplineFEM: reading 'initial conditions' ...\n\n";

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

    /*
    case  17: cout << "     HBSplineFEM: reading 'control' ...\n\n";

             if (tol > -1)         prgError(1,fct,"'control' has already been read!");

             line.getNextLine(Ifile);

             nw = line.split(&word);

             //  if (nw < 3)                    prgError(1,fct,"input error in 'control'!");

             //if (!word[0].toDbl(&PENALTY_NUM))      prgError(2,fct,"input error in 'control'!");

             if (!word[1].toDbl(&tol))      prgError(2,fct,"input error in 'control'!");

             if (!word[2].toInt(&lumpType)) prgError(3,fct,"input error in 'control'!");

             if (lumpType < 0) lumpType = 0;

             if (!word[3].toInt(&tis))      prgError(4,fct,"input error in 'control'!");

             if (tis < 0) tis = 0;

             for(i=0;i<10;i++)
              td[i] = 0.0;

             nn = nw; if (nn > 13) nn = 13;

             for(i=4; i<nn; i++)
               if (!word[i].toDbl(&(td[i-4]))) prgError(5,fct,"input error in 'control'!");

             for(i=0; i<nw; i++) word[i].free(); delete [] word;

             line.getNextLine(Ifile);

             break;
    */

    case  17: cout << "     HBSplineFEM: reading 'control parameters' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'control parameters'!");

            if( lvdTmp[0].n < 3)
              cerr <<  " Error in (( HBSplineFEM: reading 'control parameters' )) " << endl;

            tol      = lvdTmp[0][0];
            tis      = (int) lvdTmp[0][1];
            rhoInfty = lvdTmp[0][2];
            
            cout << tol << '\t' << tis << '\t' << rhoInfty << endl;

            break;

    case  18: cout << "     HBSplineFEM: reading 'profile to refine' ...\n\n";

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

    case  19: cout << "     HBSplineFEM: reading 'staggered or monolithic' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'staggered or monolithic'!");

            //stagParams.resize(lvdTmp.n);

            assert(lvdTmp[0].n == 3);

            cout << lvdTmp[0].n << endl;
            cout << lvdTmp[0] << endl;
            for(i=0;i<lvdTmp[0].n;i++)
              stagParams[i] = lvdTmp[0][i];
            
            printVector(stagParams);

            STAGGERED = (lvdTmp[0][0] == 0);
            cout << " STAGGERED " << STAGGERED << endl;

            break;

    case  20: cout << "     HBSplineFEM: reading 'Galerkin or Leastsquares' ...\n\n";

             if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
               prgError(1,fct,"invalid input in 'Galerkin or Leastsquares'!");

             LSFEM_FLAG = (lviTmp[0][0] == 1);
             cout << " LSFEM_FLAG " << LSFEM_FLAG << endl;

             break;

    case  21: cout << "     HBSplineFEM: reading 'immersed body data' ...\n\n";

             break;

    case  22: cout << "     HBSplineFEM: reading 'immersed points' ...\n\n";

             break;

    case  23: cout << "     HBSplineFEM: reading 'immersed integration elements' ...\n\n";

             break;

    case  24: cout << "     HBSplineFEM: reading 'solid elements' ...\n\n";

             break;


    case  27: cout << "     HBSplineFEM: reading 'immersed point boundary condition' ...\n\n";

             break;

    case  28: cout << "     HBSplineFEM: reading 'immersed body output' ...\n\n";

             break;

    case  29: cout << "     HBSplineFEM: reading 'analytical function' ...\n\n";

             line.getNextLine(Ifile);

             anlySolnType = line;
             
             cout << line << endl;

             line.getNextLine(Ifile);

             break;

    case  30: cout << "     HBSplineFEM: reading 'rigid body mass' ...\n\n";

             break;

    case  31: cout << "     HBSplineFEM: reading 'rigid body damping' ...\n\n";

             break;

    case  32: cout << "     HBSplineFEM: reading 'rigid body stiffness' ...\n\n";

             break;

    case  33: cout << "     HBSplineFEM: reading 'rigid body degree of freedom' ...\n\n";

             break;

    case  34: cout << "     HBSplineFEM: reading 'fluid properties' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'fluid properties'!");

            fluidProps.resize(lvdTmp[0].n);

            for(i=0;i<lvdTmp[0].n;i++)
              fluidProps[i] = lvdTmp[0][i];

             break;

    case  35: cout << "     HBSplineFEM: reading 'contact elements' ...\n\n";

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

    case  36: cout << "     HBSplineFEM: reading 'rigid body prescribed motion' ...\n\n";

            break;

    case -1: // go and inherit from DOMAIN

             this->Domain::readInputData(Ifile,line);

        break;
  }

  ndof = ndf;

  return;
}




void HBSplineFEM::prepareInputData()
{
    printf("\n     HBSplineFEM::prepareInputData()  .... STARTED ...\n");

    HBSplineBase::prepareInputData();

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

        cout << bb << '\t' << aa << '\t' << lme->IsActive() << '\t' << lme->gausspoints.size() << endl;

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
          printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %5d, \n", geom[0], geom[1], geom[2], r);
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
    SolnData.SetRho(td[0]);
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

