
#include "StandardFEM.h"
#include "DataBlockTemplate.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "PropertyTypeEnum.h"
#include "NewLagrangeElement.h"
#include "List.h"
#include "Functions.h"
#include "KimMoinFlow.h"
#include <algorithm>    // std::any_of


extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern MpapTime           mpapTime;
extern ComputerTime       computerTime;


using namespace std;


StandardFEM::StandardFEM()
{
    MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);

    //cout << " this_mpi_proc " << this_mpi_proc << endl;
    //cout << " n_mpi_procs " << n_mpi_procs << endl;

    //PetscErrorCode ierr;
    //PetscInitialize(NULL,NULL,(char *)0,NULL);
    //PetscInitialize(NULL, NULL, "petsc_options.dat", NULL);

    geom.resize(3);
    param.resize(3);
    normal.resize(3);

    PENALTY_NUM = 1.0;

    ndof = nElem = 0;

    velDOF = presDOF = fluidDOF = totalDOF = 0;

    totalDOF  =  nNode = npElem = ndof = 0;

    firstIter = STAGGERED = true;

    tol = -2.0;

    solverEigen  = NULL;
    solverPetsc  = NULL;
    elems  = NULL;

    localStiffnessError = 0;
    filecount = 0;
    IterNum = 1;

    uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New(); 
    pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    lineVTK      =  vtkSmartPointer<vtkLine>::New();
    quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    hexVTK       =  vtkSmartPointer<vtkHexahedron>::New();
    vertexVTK    =  vtkSmartPointer<vtkVertex>::New();

    triaVTK      =  vtkSmartPointer<vtkTriangle>::New();
    polygonVTK   =  vtkSmartPointer<vtkPolygon>::New();
    tetraVTK     =  vtkSmartPointer<vtkTetra>::New();
    pyramidVTK   =  vtkSmartPointer<vtkPyramid>::New();
    wedgeVTK     =  vtkSmartPointer<vtkWedge>::New();

    procIdVTK    =  vtkSmartPointer<vtkIntArray>::New();
    vecVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    scaVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    scaVTK2      =  vtkSmartPointer<vtkFloatArray>::New();
    vecVTK2      =  vtkSmartPointer<vtkFloatArray>::New();
    cellDataVTK  =  vtkSmartPointer<vtkFloatArray>::New();
    cellDataVTK2 =  vtkSmartPointer<vtkFloatArray>::New();
  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    // add new type

    DomainType *standardfem = domain.newType(STANDARDFEM, ROOTDOMAIN);

    if(standardfem == NULL) return;  // domain type already exists

  standardfem->key.addNew("element type", //0
                          "material type", //1
                          "nodes", //2
                          "elements", //3
                          "prescribed boundary conditions", //4
                          "nodal forces", //5
                          "nodal data output", //6
                          "control parameters"); //7

}


StandardFEM::~StandardFEM()
{
  if(elems != NULL)
  {
    for(int ii=0;ii<nElem;ii++)
      delete elems[ii];

    delete [] elems;
    elems = NULL;
  }

  if(solverEigen != NULL)
    delete solverEigen;
  solverEigen = NULL;

  if(solverPetsc != NULL)
    delete solverPetsc;
  solverPetsc = NULL;
}




void  StandardFEM::readInputData(std::ifstream &Ifile, MyString &line)
{
  char fct[] = "StandardFEM::readInputData";

  MyString tmpl, *word;

  char tmp[30];

  int nw, i, j, k, n, nn, ii, bb;
  
  double fact;

  MyStringList   sTmp;
  List<myVector<int> > lviTmp;
  List<myVector<double> > lvdTmp;

  DataBlockTemplate t1, t2;

  switch (domain[STANDARDFEM].key.whichBegins(line))
  {
    case  0: //cout << "     StandardFEM: reading 'element type' ...\n\n";

            SolnData.ElemProp.add(new PropertyItem(ELEMENTTYPE));
            SolnData.ElemProp[SolnData.ElemProp.n-1].readInputData(Ifile,line,"input error in 'element type'!");

            //ElemProp.push_back(PropertyItem(ELEMENTTYPE));
            //ElemProp[ElemProp.size()-1].readInputData(Ifile,line,"input error in 'element type'!");

            break;

    case  1: //cout << "     StandardFEM: reading 'material type' ...\n\n";

            //cout << " llllllllllll " << SolnData.MatlProp.n << endl;

            SolnData.MatlProp.add(new PropertyItem(MATERIAL));
            //cout << " llllllllllll " << SolnData.MatlProp.n << endl;
            SolnData.MatlProp[SolnData.MatlProp.n-1].readInputData(Ifile,line,"input error in 'material type'!");

            //MatlProp.push_back(PropertyItem(MATERIAL));
            //MatlProp[MatlProp.size()-1].readInputData(Ifile,line,"input error in 'material type'!");

            //cout << " llllllllllll " << endl;

            break;

    case  2: //cout << "     StandardFEM: reading 'nodes' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1, fct, "invalid specification of 'nodes' !");

            nodePosData.resize(lvdTmp.n);

            for(i=0; i<lvdTmp.n; i++)
            {
              for(j=0;j<(lvdTmp[i].n-1);j++)
                nodePosData[i][j] = lvdTmp[i][j+1];
            }

            break;

    case  3: //cout << "     StandardFEM: reading 'elements' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'elements'!");

            elemConn.resize(lviTmp.n);

            for(i=0;i<lviTmp.n;i++)
            {
              //cout << lviTmp[i] << endl;
              elemConn[i].resize(lviTmp[i].n-1 );
              if(lviTmp[i].n < 1)
                prgError(2, fct, "invalid number of 'elements' !");

              for(j=1;j<lviTmp[i].n;j++)
                elemConn[i][j-1] = lviTmp[i][j] - 1;
             }

            break;

    case  4: //cout << "     StandardFEM: reading 'prescribed boundary conditions' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'prescribed boundary conditions'!");

            DirichletBCs.resize(lvdTmp.n);

            for(i=0;i<lvdTmp.n;i++)
            {
              if(lvdTmp[i].n < 1)
                prgError(2, fct, "invalid number of 'prescribed boundary conditions' !");

              DirichletBCs[i].resize(lvdTmp[i].n);

              DirichletBCs[i][0] = lvdTmp[i][0]-1;
              DirichletBCs[i][1] = lvdTmp[i][1]-1;
              DirichletBCs[i][2] = lvdTmp[i][2];
            }

            break;

    case  5: //cout << "     StandardFEM: reading 'nodal forces' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
               prgError(1,fct,"invalid input in 'nodal forces'!");

             nodeForcesData.resize(lvdTmp.n);

             for(i=0;i<lvdTmp.n;i++)
             {
                nodeForcesData[i].resize(lvdTmp[i].n);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'nodal forces' !");

                for(j=0;j<lvdTmp[i].n;j++)
                  nodeForcesData[i][j] = lvdTmp[i][j];
             }

             break;

    case  6: //cout << "     StandardFEM: reading 'nodal data output' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
               prgError(1,fct,"invalid input in 'nodal data output'!");

             OutputData.resize(lvdTmp.n);

             for(i=0;i<lvdTmp.n;i++)
             {
                OutputData[i].resize(lvdTmp[i].n);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'nodal data output' !");

                for(j=0;j<lvdTmp[i].n;j++)
                  OutputData[i][j] = lvdTmp[i][j];
             }

             //printVector(OutputData);

             break;

    case   7: //cout << "     StandardFEM: reading 'control parameters' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'control parameters'!");

            if( lvdTmp[0].n < 3)
              cerr <<  " Error in (( StandardFEM: reading 'control parameters' )) " << endl;

            tol      = lvdTmp[0][0];
            tis      = (int) lvdTmp[0][1];
            rhoInfty = lvdTmp[0][2];

            //cout << tol << '\t' << tis << '\t' << rhoInfty << endl;

            break;

    case -1: // go and inherit from DOMAIN

            this->Domain::readInputData(Ifile,line);

            DIM  = ndm;
            ndof = ndf;

            break;
  }

  return;
}




void StandardFEM::readInput(ifstream& fname)
{
  MyString  line;

    line = "nodes";
    cout << " nodes " << endl;

    readInputData(fname, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////


    line = "elements";

    readInputData(fname, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "prescribed boundary conditions";

    readInputData(fname, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "nodal forces";

    readInputData(fname, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "nodal data output";

    readInputData(fname, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "element type";

    readInputData(fname, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "material type";

    readInputData(fname, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "control";

    readInputData(fname, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////


  return;
}




void StandardFEM::readNodes(ifstream& fname)
{
    MyString  line;

    line = "nodes";

    readInputData(fname, line);

    for(int ii=0; ii<nodePosData.size(); ii++)
      cout << ii << '\t' << nodePosData[ii][0] << '\t' << nodePosData[ii][1] << endl;

    return;
}


void  StandardFEM::readElementConnectivity(ifstream& fname)
{
    MyString  line;
    line = "elements";

    readInputData(fname, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    return;
}


void  StandardFEM::readElementProps(ifstream& fname)
{
    MyString  line;

    line = "element type";

    readInputData(fname, line);

    return;
}



void  StandardFEM::readMaterialProps(ifstream& fname)
{
    MyString  line;
    line = "material type";

    readInputData(fname, line);

    return;
}


void StandardFEM::prepareElemProp()
{
    char fct[] = "StandardFEM::prepareElemProp()";

    if(SolnData.ElemProp.n < 1)
      prgError(1,fct,"'element property data ' missing!");

    char *elmTypeNames[] = STANDARDFEM_ELEMENT_TYPE_NAMES;

    for(int ii=0; ii<SolnData.ElemProp.n; ii++)
    {
      // assign correct id of the element type (as stored in the database) based on the element name
      SolnData.ElemProp[ii].id = SolnData.ElemProp[ii].name.which(elmTypeNames);

      if(SolnData.ElemProp[ii].id < 0)
        prgError(2,fct,"unknown element type name!");
    }

    return;
}





void StandardFEM::prepareMatlProp()
{
    char fct[] = "StandardFEM::prepareMatlProp()";

    if(SolnData.MatlProp.n < 1)
      prgError(1,fct,"'patch material property data ' missing!");

    char *matlTypeNames[] = MATERIAL_TYPE_NAMES;

    for(int ii=0; ii<SolnData.MatlProp.n; ii++)
    {
      // assign correct id of the material type (as stored in the database) based on the material name
      SolnData.MatlProp[ii].id = SolnData.MatlProp[ii].name.which(matlTypeNames);

      if(SolnData.MatlProp[ii].id < 0)
        prgError(2,fct,"unknown element type name!");
    }

    return;
}






void StandardFEM::prepareInputData()
{
    //printf("\n     StandardFEM::prepareInputData()  .... STARTED ...\n");

    int ii, jj, kk, ee, aa, bb, cc, ind, count, gp, r;
    bool  nContX=false, nContY=false, nContZ=false;

    // go and inherit from ancestors
    Domain::prepareInputData();

    assert(DIM > 0 && DIM < 4);

    //cout << " DIM  = " << DIM << endl;
    //cout << " ndof = " << ndof << endl;

    // ==================================================
    //
    // Check the  consistency of input data
    //
    // ==================================================

    //checkInputData();

    // Prepare patchElemProp data. Assign suitable id(in the database) based on the element name

    if(SolnData.ElemProp.n > 0)
      prepareElemProp();

    // Prepare patchMatlProp data. Assign suitable id(in the database) based on the material name

    if(SolnData.MatlProp.n > 0)
      prepareMatlProp();

    ///////////////////////////////////////////////////////////////////
    //
    ///////////////////////////////////////////////////////////////////

    nNode = nodePosData.size();
    nElem = elemConn.size();
    npElem = elemConn[0].size()-3;

    cout << " nNode and nElem " << nNode << '\t' << nElem << '\t' << SolnData.ElemProp[elemConn[0][0]].id << endl;

    IEN.resize(nElem);
    LM.resize(nElem);

    elems = new LagrangeElement* [nElem];

    nElem_Constraint = 0;
    for(ee=0;ee<nElem;ee++)
    {
      // contact element along X axis, 2D
      if(SolnData.ElemProp[elemConn[ee][0]].id == 33)
      {
        nContX = true;
        nElem_Constraint++;
      }
      // contact element along Y axis, 2D
      if(SolnData.ElemProp[elemConn[ee][0]].id == 34)
      {
        nContY = true;
        nElem_Constraint++;
      }
      // contact element along X axis, 3D
      if(SolnData.ElemProp[elemConn[ee][0]].id == 35)
      {
        nContX = true;
        nElem_Constraint++;
      }
      // contact element along Y axis, 3D
      if(SolnData.ElemProp[elemConn[ee][0]].id == 36)
      {
        nContY = true;
        nElem_Constraint++;
      }
      // contact element along Z axis, 3D
      if(SolnData.ElemProp[elemConn[ee][0]].id == 37)
      {
        nContZ = true;
        nElem_Constraint++;
      }

      elems[ee] = NewLagrangeElement(SolnData.ElemProp[elemConn[ee][0]].id);

      elems[ee]->elmType  =  elemConn[ee][0];
      elems[ee]->matType  =  elemConn[ee][1];
      elems[ee]->secType  =  elemConn[ee][2];

      vector<int>  vecTemp;
      for(ii=0;ii<elemConn[ee].size()-3;ii++)
        vecTemp.push_back(elemConn[ee][3+ii]);

      elems[ee]->nodeNums = vecTemp;

      IEN[ee] = vecTemp;

      elems[ee]->SolnData = &(SolnData);
      elems[ee]->GeomData = &(GeomData);
      //elems[ee]->prepareElemData();
    }

    ndof = elems[0]->getNdofPerNode();

    int ndof_temp1 = ndof;
    int ndof_temp2 = ndof_temp1;

    PetscPrintf(MPI_COMM_WORLD, "\n    nContX = %5d ...\n", nContX);
    PetscPrintf(MPI_COMM_WORLD, "\n    nContY = %5d ...\n", nContY);
    PetscPrintf(MPI_COMM_WORLD, "\n    nContZ = %5d ...\n", nContZ);

    if(nContX || nContY || nContZ)
    {
      ndof_temp2 += DIM;
    }

    NodeType.resize(nNode);
    ID.resize(nNode);

    for(ii=0;ii<nNode;ii++)
    {
      NodeType[ii].resize(ndof_temp1, false);
      ID[ii].resize(ndof_temp2, -1);
    }


    //npElem = elems[0]->nodeNums.size() ;
    //cout << " AAAAAAAAAAAAAAAAA " << npElem << endl;


    ///////////////////////////////////////////////////////////////////
    //
    // set GeomData details
    //
    ///////////////////////////////////////////////////////////////////

    GeomData.setDimension(DIM);
    GeomData.setNdof(ndof);

    if(SolnData.ElemProp[elemConn[0][1]].id == 4 || SolnData.ElemProp[elemConn[0][1]].id ==  6 || 
       SolnData.ElemProp[elemConn[0][1]].id == 9 || SolnData.ElemProp[elemConn[0][1]].id == 12 || 
       SolnData.ElemProp[elemConn[0][1]].id == 13)
      GeomData.setNGP(SolnData.ElemProp[0].data[0]);
    else
      GeomData.setNGP(1);

    GeomData.build();
    GeomData.setNodalPositions(nodePosData);

    //cout << " AAAAAAAAAAAAAAAAA " << endl;

    ///////////////////////////////////////////////////////////////////
    //
    // set SolnData details
    //
    ///////////////////////////////////////////////////////////////////

    SolnData.setTimeIncrementType(tis);
    SolnData.setSpectralRadius(rhoInfty);

    //cout << " totalDOF =  " << totalDOF << endl;
    cout << " nElem_Constraint =  " << nElem_Constraint << endl;

    totalDOF = nNode*ndof;
    totalDOF += nElem_Constraint;

    SolnData.initialise(totalDOF, 0, 0, 0);

    soln.resize(totalDOF);
    soln.setZero();

    solnInit = soln;
    reac  =  soln;

    vector<int>  solidElems, fluidElems;
    //1 2 3 4 7 8 9 10  15 18 19 20 21 
    solidElems.push_back(1);
    solidElems.push_back(2);
    solidElems.push_back(3);
    solidElems.push_back(4);
    solidElems.push_back(7);
    solidElems.push_back(8);
    solidElems.push_back(9);
    solidElems.push_back(10);
    solidElems.push_back(15);
    solidElems.push_back(18);
    solidElems.push_back(19);
    solidElems.push_back(20);
    solidElems.push_back(21);
    solidElems.push_back(26);
    solidElems.push_back(27);
    solidElems.push_back(28);
    solidElems.push_back(29);
    solidElems.push_back(30);
    solidElems.push_back(31);

    fluidElems.push_back(11);
    fluidElems.push_back(12);
    fluidElems.push_back(13);
    fluidElems.push_back(14);
    fluidElems.push_back(22);
    fluidElems.push_back(23);
    fluidElems.push_back(24);
    fluidElems.push_back(25);



    int ttt = SolnData.ElemProp[elemConn[0][0]].id;

    std::vector<int>::iterator it;
    it = find (solidElems.begin(), solidElems.end(), ttt);

    if( it != solidElems.end() )
    {
      SolnData.setPhysicsTypetoSolid();
      PHYSICS_TYPE = PHYSICS_TYPE_SOLID;
    }
    else
    {
      SolnData.setPhysicsTypetoFluid();
      PHYSICS_TYPE = PHYSICS_TYPE_FLUID;
    }

    //cout << " AAAAAAAAAAAAAAAAA " << endl;

    //printf("     StandardFEM::prepareInputData()  .... FINISHED ...\n\n");

    return;
}



void StandardFEM::prepareInteractions(void)
{
    // go and inherit from ancestors

    Domain::prepareInteractions();

    printf("     StandardFEM::prepareInteractions()  .... STARTED ...\n");

    char fct[] = "StandardFEM::prepareInteractions";

    printf("     NOTHING to be done here ...\n");

    printf("     StandardFEM::prepareInteractions()  .... FINISHED ...\n\n");

    return;
}



void StandardFEM::printInfo()
{
/*
  printf("\n Background Fluid Grid information \n");
    printf("\n ----------------------------------------------------------------------------- \n");
    printf("\t   Problem   Dimension        =  %5d\n\n", DIM);
    printf("\t   Polynomial Degree          =  %5d\t%5d\t%5d\n\n", degree[0], degree[1], degree[2]);
    printf("\t   Number of Elements         =  %5d\t%5d\t%5d\n\n", nelem[0], nelem[1], nelem[2]);
    printf("\t   Origin of the grid         =  %12.6f\t%12.6f\t%12.6f\n\n", origin[0], origin[1], origin[2]);
    printf("\t   Grid Length                =  %12.6f\t%12.6f\t%12.6f\n\n", gridLEN[0], gridLEN[1], gridLEN[2]);
    printf("\t   MAXIMUM REFINEMENT LEVEL   =  %5d\n\n", MAX_LEVEL);
    printf("\t   DOF at each Control Point  =  %5d\n\n", ndof);
    printf("\t   Total number of DOF        =  %5d\n\n", totalDOF);
    printf("\n ----------------------------------------------------------------------------- \n");
*/
    return;
}


void StandardFEM::findMinMaxX(double *xmn, double *xmx, bool defFlg)
{
    //xmn[0] = origin[0];
    //xmn[1] = origin[1];

    //xmx[0] = origin[0] + gridLEN[0];
    //xmx[1] = origin[1] + gridLEN[1];

    return;
}


double  StandardFEM::computeGeometry(const int dir, double param)
{
  return  0.0;
}


void  StandardFEM::computeGeometry(const myPoint& param, myPoint& geom)
{
  return;
}


void StandardFEM::assignBoundaryConditions()
{

}



bool StandardFEM::converged()
{
  char fct[] = "StandardFEM::converged";

  if (rNorm < tol && localStiffnessError == 0)
    return true;

  return false;
}




bool StandardFEM::diverging(double factor)
{
  if (rNormPrev > -0.1 && (rNorm / rNormPrev) > factor) return true;

  if (localStiffnessError != 0) return true;

  if (prgNAN(rNorm)) return true;

  return false;
}


void StandardFEM::printComputerTime(bool reset, int detailFlg)
{
  COUT << "----------------------------------------------------\n";

  COUT; printf("StandardFEM::calcStiffnessRes:  %7.3f sec ->%5.1f %\n",
               ctimCalcStiffRes, ctimCalcStiffRes/ctimSinceLastCall*100.);

  COUT; printf("StandardFEM::factoriseSolvUpdt: %7.3f sec ->%5.1f %\n",
               ctimFactSolvUpdt, ctimFactSolvUpdt/ctimSinceLastCall*100.);

  if (reset)
  {
    ctimFactSolvUpdt = 0.;
    ctimCalcStiffRes = 0.;
  }

  return;
}



void StandardFEM::setInitialConditions()
{
    int ii, jj, nn, dof, aa, ind;

    KimMoinFlow  analy(SolnData.ElemProp[0].data[4], SolnData.ElemProp[0].data[5], 1.0);
    double  xx, yy, zz, fact;
    double  specVal, dt=mpapTime.dt, mtwodt=-2.0*dt, mthreedt=-3.0*dt, omega=105.0;

    for(aa=0;aa<nNode;aa++)
    {
      xx = GeomData.NodePosOrig[aa][0];
      yy = GeomData.NodePosOrig[aa][1];
      zz = GeomData.NodePosOrig[aa][2];

      ind = aa*ndof;

      /*
      SolnData.var1[ind]    = analy.computeValue(0, xx, yy, 0.0);
      SolnData.var1[ind+1]  = analy.computeValue(1, xx, yy, 0.0);
      SolnData.var1[ind+2]  = analy.computeValue(2, xx, yy, 0.0);

      SolnData.var1Prev[ind]   = analy.computeValue(0, xx, yy, -dt);
      SolnData.var1Prev[ind+1] = analy.computeValue(1, xx, yy, -dt);
      SolnData.var1Prev[ind+2] = analy.computeValue(2, xx, yy, -dt);

      SolnData.var1Prev2[ind]   = analy.computeValue(0, xx, yy, mtwodt);
      SolnData.var1Prev2[ind+1] = analy.computeValue(1, xx, yy, mtwodt);
      SolnData.var1Prev2[ind+2] = analy.computeValue(2, xx, yy, mtwodt);

      SolnData.var1Prev3[ind]   = analy.computeValue(0, xx, yy, mthreedt);
      SolnData.var1Prev3[ind+1] = analy.computeValue(1, xx, yy, mthreedt);
      SolnData.var1Prev3[ind+2] = analy.computeValue(2, xx, yy, mthreedt);
      */
      //SolnData.var1Dot[ind]    = -omega*yy;
      //SolnData.var1Dot[ind+1]  =  omega*xx;

      fact = 100.0*sin(PI*yy/12.0);

      //SolnData.var1Dot[ind]    =  fact*zz;
      //SolnData.var1Dot[ind+1]  =  0.0;
      //SolnData.var1Dot[ind+2]  = -fact*xx;

      //SolnData.var1Dot[ind]    =  5.0*yy/3.0;
      //SolnData.var1Dot[ind+1]  =  0.0;
      //SolnData.var1Dot[ind+2]  =  0.0;

      //SolnData.var1[ind]    =  -0.25*yy;
      //SolnData.var1[ind+1]  =  0.0;
      //SolnData.var1[ind+2]  =  0.0;
    }

    //SolnData.updateIterStep();

  return;
}




void StandardFEM::setTimeParam()
{
  SolnData.setTimeParam();
  // solve for the initial profile
  //cout << " aaaaaaaaaaaaaaa " << endl;

  return;
}



void StandardFEM::timeUpdate()
{
  firstIter = true;
  localStiffnessError = 0;
  iterCount = 1;
  filecount++;

  //cout << " aaaaaaaaaaaaaaa " << endl;

  SolnData.timeUpdate();

  //
  int ii, jj, nn, dof, aa, ind;

  double fact1=1.0, specVal;
         //if(mpapTime.cur <= 5.0e-3)
         //{
           //fact1 = 0.5*( 1.0-cos(628.3185*mpapTime.cur));
           //fact1 = 0.01*0.5*( 1.0-cos(20.0*mpapTime.cur)) ;
           //fact = sin(628.3185*mpapTime.cur);
           //fact1 = 1.0;
           //fact = sin(20.0*mpapTime.cur);
         //}
         //else
           //fact1 = 0.0;
    //fact1 = mpapTime.dt;
    //fact1 = mpapTime.cur - mpapTime.dt;

    fact1 = timeFunction[0].prop;
    fact1 = 1.0;

  if(firstIter)
  {
    int ii, jj, nn, dof, aa, ind;

    KimMoinFlow  analy(SolnData.ElemProp[0].data[4], SolnData.ElemProp[0].data[5], 1.0);
    double  xx, yy;

    for(aa=0;aa<DirichletBCs.size();aa++)
    {
      nn  = (int) (DirichletBCs[aa][0]);
      dof = (int) (DirichletBCs[aa][1]);
      specVal = DirichletBCs[aa][2];

      //itint = find(GlobalPointNumbers.begin(), GlobalPointNumbers.end(), nn);
      //nn   = distance(GlobalPointNumbers.begin(), itint);

      xx = GeomData.NodePosOrig[nn][0];
      yy = GeomData.NodePosOrig[nn][1];

      ind = nn*ndof+dof;

      //SolnData.var1[ind] = analy.computeValue(dof, xx, yy, mpapTime.cur);

      //SolnData.var1[ind] = fact1*SolnData.var1applied[ind];

      //SolnData.var1[ind] = fact1*specVal ;
    }
  }
  //

  //cout << " aaaaaaaaaaaaaaa " << endl;

  updateIterStep();

  //cout << " aaaaaaaaaaaaaaa " << endl;

  return;
}



void StandardFEM::updateIterStep()
{
  SolnData.updateIterStep();

  //cout << " PHYSICS_TYPE = " << PHYSICS_TYPE << endl;

  if( PHYSICS_TYPE == PHYSICS_TYPE_SOLID )
  {
    int  ii, ind, bb;

    for(bb=0;bb<nNode;bb++)
    {
      for(ii=0;ii<DIM;ii++)
      {
        ind = bb*ndof+ii;
        //cout << bb << '\t' << ii << '\t' << ind << '\t' << endl;
        //cout << SolnData.var1[ind] << '\t' << SolnData.var1Cur[ind] << endl;

        GeomData.NodePosNew[bb][ii] = GeomData.NodePosOrig[bb][ii] + SolnData.var1[ind];
        GeomData.NodePosCur[bb][ii] = GeomData.NodePosOrig[bb][ii] + SolnData.var1Cur[ind];

        // second-order based formulation
        GeomData.specValNew[bb][ii] = SolnData.var1Dot[ind];
        GeomData.specValCur[bb][ii] = SolnData.var1DotCur[ind];
      }
    }
  }

  return;
}




void StandardFEM::reset()
{
  SolnData.reset();

  return;
}



void StandardFEM::printData(int, int)
{
  return;
}



void StandardFEM::addExternalForces()
{
  return;
}



void StandardFEM::writeNodalData()
{
    int ii, jj, bb, ee, type, nn1, nn2=0, dof;
    double val, fact=1.0;

    for(ee=0; ee<OutputData.size(); ee++)
    {
      type  = (int) (OutputData[ee][0]);
      nn2   = (int) (OutputData[ee][1] - 1);
      dof   = (int) (OutputData[ee][2] - 1);

      //cout << ee << '\t' << type << '\t' << dof << '\t' << fact << '\t' << nn2 << endl;

      switch(type)
      {
        case  1 : // total force on all the requested nodes

              val = 0.0;
              for(ii=0; ii<(OutputData.size()-3); ii++)
                val += SolnData.force[OutputData[ee][3+ii]*ndof+dof];

        break;

        case  2 : // total reaction on all the requested nodes

              val = 0.0;
              for(ii=0; ii<(OutputData[ee].size()-3); ii++)
              {
                nn2 =  (int) (OutputData[ee][3+ii]-1) ;
                val += SolnData.reac[nn2*ndof+dof];
              }

        break;

        case  3 : // displacement

              val = SolnData.var1[nn2*ndof+dof];

        break;

        case  4 : // velocity

              val = SolnData.var1Dot[nn2*ndof+dof];

        break;

        case  5 : // acceleration

               val = SolnData.var1DotDot[nn2*ndof+dof];

        break;

        default : 

               prgError(1,"StandardFEM::writeNodalData()","invalid value of 'type'!");

        break;
      }

      val *= fact;

      char        tmp[200];
      MyString    tmpStr;

      //sprintf(tmp," \t %12.8f \t %12.6E \n", mpapTime.cur, val);
      sprintf(tmp," \t %12.6E", val);

      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);
    }

  return;
}


  
  
  
