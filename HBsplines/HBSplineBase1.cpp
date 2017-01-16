
#include "HBSplineBase.h"

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
#include "ImmersedIntegrationElement.h"
#include "util.h"
#include "DistFunctions.h"
#include "myConstants.h"
#include "NurbsUtilities.h"
#include "StandardFEM.h"



extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern MpapTime           mpapTime;
extern ComputerTime       computerTime;


using namespace std;


typedef TreeNode<1> node1D;
typedef TreeNode<2> node2D;
typedef TreeNode<3> node3D;

template<>
int node1D::nodecount = 0;

template<>
int node2D::nodecount = 0;

template<>
int node3D::nodecount = 0;


template<>
int node1D::MAX_LEVEL = 0;

template<>
int node2D::MAX_LEVEL = 0;

template<>
int node3D::MAX_LEVEL = 0;

int ImmersedSolid::count = 0;


int ImmersedIntegrationElement::pointcount = 0;



HBSplineBase::HBSplineBase()
{
    //cout << "  HBSplineBase::HBSplineBase() ... constructor " << endl;

    root = NULL;

    for(int ii=0;ii<3;ii++)
    {
      nelem[ii]  = 0;
      degree[ii] = 0;
      nlbf[ii]   = 0;
      origin[ii] = 0.0;
      knotsLeft[ii] = 0.0;
      knotsRight[ii] = 1.0;
      gridLEN[ii]    = 1.0;
    }

    geom.setZero();
    param = geom;
    normal = geom;

    CURRENT_LEVEL = MAX_LEVEL = numActiveElem = ndof = numActiveBasis = totnumel = 0;

    gridBF1 = gridBF2 = gridBFtotal = totalDOF = 0;

    STAGGERED = true;

    LSFEM_FLAG = false;

    STATUS_BCS = PERIODIC_BCS  = false;

    firstIter = CREATE_POSTPROCESS_GRID = true;
    
    GRID_CHANGED = IB_MOVED = STABILISED = true;

    NodeNumsAtLevel.resize(10);

    tol = -2.0;

    solverEigen = NULL;
    solverPetsc = NULL;
    localStiffnessError = 0;
    filecount = 0;
    IterNum = 1;

    mapperVTK    =  vtkSmartPointer<vtkDataSetMapper>::New();
    actorVTK     =  vtkSmartPointer<vtkActor>::New();
    uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New(); 
    pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    lineVTK      =  vtkSmartPointer<vtkLine>::New();
    quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    tetVTK       =  vtkSmartPointer<vtkTetra>::New();
    hexVTK       =  vtkSmartPointer<vtkHexahedron>::New();
    vertexVTK    =  vtkSmartPointer<vtkVertex>::New();

    vecVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    scaVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    scaVTK2      =  vtkSmartPointer<vtkFloatArray>::New();
    vecVTK2      =  vtkSmartPointer<vtkFloatArray>::New();
    cellDataVTK  =  vtkSmartPointer<vtkFloatArray>::New();
    cellDataVTK2 =  vtkSmartPointer<vtkFloatArray>::New();
  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    pointsVTKfluidgrid   =   vtkSmartPointer<vtkPoints>::New();
    pointsPolydata       =   vtkSmartPointer<vtkPolyData>::New();
    mergePoints          =   vtkSmartPointer<vtkMergePoints>::New();
    uGridVTKfluid        =  vtkSmartPointer<vtkUnstructuredGrid>::New(); 

    //cout << "  HBSplineBase::HBSplineBase() ... constructor " << endl;
}


HBSplineBase::~HBSplineBase()
{
  if(root != NULL)
    root = NULL;

  if(solverEigen != NULL)
    delete solverEigen;
  solverEigen = NULL;


  if(solverPetsc != NULL)
    delete solverPetsc;
  solverPetsc = NULL;


  for(vector<node*>::iterator pObj = elems.begin(); pObj != elems.end(); ++pObj)
  {
    delete *pObj; // Note that this is deleting what pObj points to, which is a pointer
  }
  elems.clear(); // Purge the contents so no one tries to delete them again



  for(vector<ImmersedSolid*>::iterator pObj = ImmersedBodyObjects.begin(); 
      pObj != ImmersedBodyObjects.end(); ++pObj)
  {
    delete *pObj; // Note that this is deleting what pObj points to, which is a pointer
  }
  ImmersedBodyObjects.clear(); // Purge the contents so no one tries to delete them again


  for(vector<ContactElementPointToPoint2D*>::iterator pObj = contactElementObjects.begin(); 
      pObj != contactElementObjects.end(); ++pObj)
  {
    delete *pObj; // Note that this is deleting what pObj points to, which is a pointer
  }
  contactElementObjects.clear(); // Purge the contents so no one tries to delete them again

}






void HBSplineBase::prepareInputData()
{
    //printf("\n     HBSplineBase::prepareInputData()  .... STARTED ...\n");

    int ii, jj, kk, ee, aa, bb, cc, gp, r;

    // go and inherit from ancestors
    Domain::prepareInputData();

    //cout << " DIM  = " << DIM << endl;
    //cout << " ndof = " << ndof << endl;

    /*
    if(LSFEM_FLAG)
    {
      if(ndof > DIM)
      {
        printf("\n HBSplineBase : Preparing element data ...... 'ndof' can't be more than 'DIM' \n \n ");
        exit(1);
      }
    }
    */

    // ==================================================
    //
    // Check the  consistency of input data
    //
    // ==================================================

    //checkInputData();

    ///////////////////////////////////////////////////////////////////
    //
    ///////////////////////////////////////////////////////////////////

    assert(DIM > 0 && DIM < 4);

    for(ii=0;ii<DIM;ii++)
      nlbf[ii] = degree[ii] + 1;

    //cout << " AAAAAAAAAAAAAAAAA " << endl;

    ///////////////////////////////////////////////////////////////////
    //
    // set FluidSolnData details
    //
    ///////////////////////////////////////////////////////////////////
    
    GeomData.SetDimension(DIM);

    for(ii=0;ii<DIM;ii++)
    {
      GeomData.SetDegree(ii, degree[ii]);
      GeomData.SetNelem(ii, nelem[ii]);
      GeomData.SetGridLength(ii, gridLEN[ii]);
      GeomData.SetOrigin(ii, origin[ii]);
    }

    GeomData.FluidProps = fluidProps;
    //printVector(GeomData.FluidProps);

    GeomData.SetNDOF(ndf);
    GeomData.build();
    
    //cout << " BBBBBBBBBBBBBBB " << endl;

    ///////////////////////////////////////////////////////////////////
    //
    //  Build the base cartesian grid
    //
    ///////////////////////////////////////////////////////////////////

    boundaryNodes.resize(2*DIM);

    BuildBase();

    //printf("\n HBSplineBase : BuildBase        DONE \n \n ");

    /////////////////////////////////////////////////////////////
    //
    // perform local refinement if any
    //
    /////////////////////////////////////////////////////////////

    // compute the basis functions at level ZERO using algorithm #3  

    gridBF1 = 0;
    Algorithm3(0);

    /*
    for(ee=0;ee<elems.size();ee++)
    {
      if( !(elems[ee]->IsGhost()) )
      {
        //cout << " uuuuuuuuuuuuu " << endl;
        elems[ee]->calcSubdivisionMatrix();
        //cout << " uuuuuuuuuuuuu " << endl;
      }
    }
    */

    //printf("\n \t  Creating immersed boundary points \n ");

    ///////////////////////////////////////////////////
    //
    // setting up immersed bodies
    // 
    ///////////////////////////////////////////////////

    nImmSolids = ImmersedBodyObjects.size();

    //cout << " nImmSolids = " << nImmSolids << endl;

    //ImmersedSolid  *imsolid;

    if(nImmSolids > 0)
    {
      cc = 0;
      for(bb=0; bb<nImmSolids; bb++)
      {
        ImmersedBodyObjects[bb]->SetTimeIncrementType(tis);
        ImmersedBodyObjects[bb]->SetRho(rhoInfty);
        ImmersedBodyObjects[bb]->SetTolerance(tol);

        ImmersedBodyObjects[bb]->STAGGERED = STAGGERED;

        ImmersedBodyObjects[bb]->SolnData.stagParams = stagParams;

        //cout << " zzzzzzzzzzzzzzzzz " << endl;
        if( ImmersedBodyObjects[bb]->IsFlexibleBody() )
        {
          //if( ImmersedBodyData[bb][4] > 0 )
            //ImmersedBodyObjects[bb]->SolnData.ElemProp = ElemProp[ImmersedBodyData[bb][4]-1];

          //if( ImmersedBodyData[bb][5] > 0 )
            //ImmersedBodyObjects[bb]->SolnData.MatlProp = MatlProp[ImmersedBodyData[bb][5]-1];

          //for(bb=0; bb<ElemProp.n; bb++)
            //SolnData.ElemProp.add(&ElemProp[bb]);

          //for(bb=0; bb<MatlProp.n; bb++)
            //SolnData.MatlProp.add(&MatlProp[bb]);

          //cout << " 11111111111111 " << endl;
          //cout << " VVVVVVVVVVV " << endl;
        }

        //cout << " PPPPPPPPPPP " << endl;
        ImmersedBodyObjects[bb]->initialise();

        //cout << " PPPPPPPPPPP " << endl;
        GeomData.immSolidPtrs.push_back( ImmersedBodyObjects[bb] );

        GeomData.domainFixedYesNo.push_back( ImmersedBodyObjects[bb]->GetTotalDOF() == 0 ) ;
      }
    }

    cout << " Immersed solids preparation DONE ... " << endl;

    geom.setZero();
    param.setZero();

    // Refine the underlying grid if specified
    
    if(refinementData.size() > 0)
    {
      //printf("\n HBSplineBase : Refinement process ...... STARTED \n \n ");
      
      for(kk=0;kk<refinementData[1];kk++)
      {
        elemsToRefine.clear();

        switch( int (refinementData[0]) )
        {
           default :
           case 0:
              Refine(refinementData[2]);
              break;
           case 1:
              refinementforHemkerProblem();
              break;
           case 2:
              refinementforAdvDiff1D();
              break;
           case 3:
              refinementforAdvDiff2D();
              break;
           case 4:
              PointBasedRefinement(refinementData[2]);
              break;
           case 5:
              LimitBasedRefinement(kk);
              break;
        }

        findUnique(elemsToRefine);

        /*
        printf(" elemsToRefine  \n");
        for(ii=0;ii<elemsToRefine.size();ii++)
          printf("%5d \t", elemsToRefine[ii] );
        printf("\n\n\n");
        */

        MAX_LEVEL += 1;
        ApplyRefinementProcess();
        CURRENT_LEVEL += 1;
        ProcessBoundaryConditionsRefinedLevels();

        //for(ii=0;ii<elemsToRefine.size();ii++)
        //elems[elemsToRefine[ii]]->printSelf();
      }
      //printf("\n HBSplineBase : Refinement process ...... FINISHED \n \n ");
    }
    
    /////////////////////////////////
    // assign boundary conditions to the element
    //

    AssignBoundaryConditions();

    /////////////////////////////////

    if(DIM == 1)
    {
      // account for periodic boundary conditions if any
      /////////////////////////////////////////////////////
    
      cout << " PERIODIC_BCS " << PERIODIC_BCS << '\t' << gridBF1 << endl;

      if(PERIODIC_BCS)
      {
         node  *nodetmp2;
         int val;

         for(ii=0;ii<degree[0];ii++)
         {
            nodetmp2 = elems[boundaryNodes[1][0]];
            val = elems[boundaryNodes[0][0]]->LocalBasisFuncs[degree[0]-1-ii];
            for(jj=0;jj<=ii;jj++)
            {
              nodetmp2->LocalBasisFuncs[degree[0]+jj-ii] = val;
              nodetmp2 = nodetmp2->GetNeighbour(LEFT);
            }
           gridBF1 -= 1;
         }
      }
    }


    ///////////////////////////////////////////////////
    //
    // prepare element data
    // 
    ///////////////////////////////////////////////////
    
    time_t tstart, tend;

    //printf("\n HBSplineBase : Preparing element data ...... STARTED \n \n ");

    tstart = time(0);
    for(ee=0;ee<elems.size();ee++)
    {
      if( !(elems[ee]->IsGhost()) &&  elems[ee]->IsLeaf() )
        activeElements.push_back(elems[ee]->GetID());
    }

    /*
    omp_set_num_threads(4);
    printf("Max number of threads: %i \n",omp_get_max_threads());
    printf("Number of threads: %i \n",omp_get_num_threads());
    #pragma omp parallel
      printf("Number of threads: %i \n",omp_get_num_threads());

    #pragma omp parallel num_threads(8)
      printf("Hello from thread %d, nthreads %d, nprocs %d\n", omp_get_thread_num(), omp_get_num_threads(), omp_get_num_procs());
    */

    //#pragma omp parallel for private(ii)

    //printf("\n number of total  elements = %6d \n", elems.size());
    //printf("\n number of active elements = %6d \n", activeElements.size());

    for(ii=0;ii<activeElements.size();ii++)
    {
      ee = activeElements[ii];

      //cout << " elems[ee]->IsProcessed() " << '\t' << elems[ee]->GetID() << '\t' << elems[ee]->IsProcessed() << endl;
      //cout << " elems[ee]->GetID() " << '\t' << elems[ee]->GetID() << '\t' << elems[ee]->GetLevel() << endl;
      elems[ee]->prepareElemData();
      //cout << " uuuuuuuuuuuuu " << endl;
      elems[ee]->calcSubdivisionMatrix();
      //cout << " uuuuuuuuuuuuu " << endl;
      //elems[ee]->initialiseDOFvalues();
      //cout << " uuuuuuuuuuuuu " << endl;
      //elems[ee]->checkPartitionOfUnity();
    }

    tend = time(0);

    //cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;

    //printf("\n HBSplineBase : Preparing element data ...... FINISHED \n \n ");

    // total number of DOF for the background grid

    gridBF2 = gridBF1;

    //gridBFtotal = gridBF1 + gridBF2;

    printf("\n \t   Number of basis functions in the background grid  =  %5d\n\n", gridBF1);

    /////////////////////////////////////////
    // create contact elements
    ////////////////////////////////////////
    
    //cout << " jjjjjjjjjjjjjjj  " << endl;
    
    ContactElementPointToPoint2D  *contElm;

    for(aa=0; aa<contElemData.size(); aa++)
    {
      contElm = new ContactElementPointToPoint2D();

      //contElm->rigidbodypointers.push_back(ImmersedBodyObjects[0]);
      //contElm->FluidSolnData = &(FluidSolnData);
      //contElm->SolnData = &(ImmersedBodyObjects[0]->SolnData);

      contactElementObjects.push_back(contElm);

      contactElementObjects[aa]->prepareElemData();
    }

    //cout << " jjjjjjjjjjjjjjj  " << endl;

    printf("     HBSplineBase::prepareInputData()  .... FINISHED ...\n\n");

    return;
}



void HBSplineBase::refinementforHemkerProblem()
{
    //
    // local refinement for Hemker problem
    /////////////////////////////////////
 
    node *nd2;
    int  ii, kk, start, N=degree[2];

    //for(kk=0;kk<degree[1];kk++)
    //{
       elemsToRefine.clear();

       if(CURRENT_LEVEL == 0)
       {
         start = nelem[0]/2 - N + degree[0];
         //cout << N << '\t' << start << endl;
         for(ii=0;ii<2*N;ii++)
           elemsToRefine.push_back(start+ii);
       }
       else
       {
         start = nelem[0] + 2*degree[0] + 4*N*(CURRENT_LEVEL-1)+N;
         for(ii=0;ii<2*N;ii++)
           elemsToRefine.push_back(start+ii);
       }
    //}

    return;
}



void HBSplineBase::refinementforAdvDiff1D()
{
    // local refinement for advection-diffusion problem in 1D
    ///////////////////////////////////////////////////////////

    node *nd, *nd1, *nd2;
    int  ii, kk, start, N=degree[2];

    for(kk=0;kk<degree[1];kk++)
    {
      elemsToRefine.clear();

      nd2=elems[boundaryNodes[1][0]];
      elemsToRefine.push_back(nd2->GetID());
    
      for(ii=0;ii<degree[2];ii++)
      {
        nd2 = nd2->GetNeighbour(LEFT);
        elemsToRefine.push_back(nd2->GetID());
      }
    }

    return;
}



void HBSplineBase::refinementforAdvDiff2D()
{
    // refinement for advection-diffusion in 2D

    node *nd, *nd1, *nd2;
    int  ii, kk, ee;

    for(kk=0;kk<degree[2];kk++)
    {
      //cout << " PPPPPPPPPPPPPPP " << endl;

      elemsToRefine.clear();

      for(ee=0;ee<NodeNumsAtLevel[kk].size();ee++)
      {
          nd = elems[NodeNumsAtLevel[kk][ee]];
          
          //cout << " nd->GetID() " << nd->GetID() << '\t' << nd->IsGhost() << endl;
          //nd->printSelf();
          if(nd->IsLeftBoundary() && nd->IsBottomBoundary())
          {
             while(!nd->IsGhost())
             {
                //cout << " nd->GetID() " << nd->GetID() << endl;
                elemsToRefine.push_back(nd->GetID());
                nd = nd->GetNeighbour(EAST);
                nd = nd->GetNeighbour(NORTH);
                //cout << " AAAAAAAAAAAAAAAAAAAAA " << nd->GetID() << '\t' << nd->IsGhost() << endl;
             }
          }
       }

       //cout << " elemsToRefine  " << endl;
       //for(ii=0;ii<elemsToRefine.size();ii++)
         //cout << '\t' << elemsToRefine[ii] ;
       //printf("\n\n\n");

       nodes2divide.clear();
       for(ii=0;ii<elemsToRefine.size();ii++)
       {
         nodes2divide.push_back(elemsToRefine[ii]);
         nd = elems[elemsToRefine[ii]];
   
         //cout << nd->GetID() << endl;
   
         nd1 = nd->GetNeighbour(WEST);
         if(nd1 != NULL)
         {
            nodes2divide.push_back(nd1->GetID());

            nd2 = nd1->GetNeighbour(NORTH);
            if(nd2 != NULL)
              nodes2divide.push_back(nd2->GetID());

            nd2 = nd1->GetNeighbour(SOUTH);
            if(nd2 != NULL)
              nodes2divide.push_back(nd2->GetID());
         }

         nd1 = nd->GetNeighbour(EAST);
         if(nd1 != NULL)
         {
           nodes2divide.push_back(nd1->GetID());

           nd2 = nd1->GetNeighbour(NORTH);
           if(nd2 != NULL)
             nodes2divide.push_back(nd2->GetID());

           nd2 = nd1->GetNeighbour(SOUTH);
           if(nd2 != NULL)
             nodes2divide.push_back(nd2->GetID());
         }

         nd1 = nd->GetNeighbour(NORTH);
         if(nd1 != NULL)
           nodes2divide.push_back(nd1->GetID());
   
         nd1 = nd->GetNeighbour(SOUTH);
         if(nd1 != NULL)
           nodes2divide.push_back(nd1->GetID());
       }

      findUnique(nodes2divide);

       elemsToRefine = nodes2divide;
    }
    return;
}



void HBSplineBase::LimitBasedRefinement(int kk)
{
    node *nd, *nd1, *nd2;
    int  ii, ee, mm;
    double *tmp1, *tmp2, *tmp3, val[6];

    assert(refinementData[1] <= refineLimitVals.size());

    //cout << " PPPPPPPPPP " << endl;
      if(DIM == 1)
      {
        //printVector(NodeNumsAtLevel[kk]);
        for(ee=0;ee<NodeNumsAtLevel[kk].size();ee++)
        {
          nd = elems[NodeNumsAtLevel[kk][ee]];

          tmp1 = nd->GetKnots(Dir1);

          param[0] = tmp1[0];
          ComputeGeometry(param, geom);

          //cout << " ee " << ee << endl;
          for(mm=0;mm<refineLimitVals.size();mm++)
          {
            if(refineLimitVals[mm][0] == (kk+1))
            {
              if( (geom[0] >= refineLimitVals[mm][1] && geom[0] <= refineLimitVals[mm][2]) )
                elemsToRefine.push_back(nd->GetID());
            }
          }
        }
      }
      else if(DIM == 2)
      {
        for(ee=0;ee<NodeNumsAtLevel[kk].size();ee++)
        {
          nd = elems[NodeNumsAtLevel[kk][ee]];

          tmp1 = nd->GetKnots(Dir1);
          tmp2 = nd->GetKnots(Dir2);

          param[0] = tmp1[0];  param[1] = tmp2[0];
          ComputeGeometry(param, geom);

          for(mm=0;mm<refineLimitVals.size();mm++)
          {
            if(refineLimitVals[mm][0] == (kk+1))
            {
              if( (geom[0] >= refineLimitVals[mm][1] && geom[0] <= refineLimitVals[mm][2]) && (geom[1] >= refineLimitVals[mm][3] && geom[1] < refineLimitVals[mm][4]) )
              {
               // printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, cell = %5d, \n", geom[0], geom[1], geom[2], nd->GetID());
                elemsToRefine.push_back(nd->GetID());
              }
            }
          }
        }
      }
      else
      {
        for(ee=0;ee<NodeNumsAtLevel[kk].size();ee++)
        {
          nd = elems[NodeNumsAtLevel[kk][ee]];

          if(!nd->IsGhost())
          {
            tmp1 = nd->GetKnots(Dir1);
            tmp2 = nd->GetKnots(Dir2);
            tmp3 = nd->GetKnots(Dir3);

            param[0] = 0.5*(tmp1[0]+tmp1[1]);  param[1] = 0.5*(tmp2[0]+tmp2[1]);  param[2] = 0.5*(tmp3[0]+tmp3[1]);
            ComputeGeometry(param, geom);
            //val[0] = ComputeGeometry(0, tmp1[0]);
            //val[1] = ComputeGeometry(0, tmp1[1]);
            //val[2] = ComputeGeometry(1, tmp2[0]);
            //val[3] = ComputeGeometry(1, tmp2[1]);
            //val[4] = ComputeGeometry(2, tmp3[0]);
            //val[5] = ComputeGeometry(2, tmp3[1]);

            //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \t %12.6f \t %12.6f \n", val[0], val[1], val[2], val[3], val[4], val[5]);
            for(mm=0;mm<refineLimitVals.size();mm++)
            {
              //printVector(refineLimitVals[mm]);
              if(refineLimitVals[mm][0] == (kk+1))
              {
                if( (geom[0] >= refineLimitVals[mm][1] && geom[0] <= refineLimitVals[mm][2]) && (geom[1] >= refineLimitVals[mm][3] && geom[1] < refineLimitVals[mm][4])  && (geom[2] >= refineLimitVals[mm][5] && geom[2] < refineLimitVals[mm][6]) )
                {
                  //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \t %12.6f \t %12.6f \n", tmp1[0], tmp1[1], tmp2[0], tmp2[1], tmp3[0], tmp3[1]);
                  //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \t %12.6f \t %12.6f \n", val[0], val[1], val[2], val[3], val[4], val[5]);
                  //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, cell = %5d, \n", geom[0], geom[1], geom[2], nd->GetID());
                  elemsToRefine.push_back(nd->GetID());
                }
              }
            }
          }
        }
      }

    return;
}



void  HBSplineBase::addNeighbourElements(int depth)
{
  if(DIM == 1)
    addNeighbourElements1D(depth);
  else if(DIM == 2)
    addNeighbourElements2D(depth);
  else
    addNeighbourElements3D(depth);

  return;
}



void  HBSplineBase::addNeighbourElements1D(int depth)
{
  return;
}


//
void  HBSplineBase::addNeighbourElements2D(int depth)
{
    int ii, aa, bb, cc;
    node *nd, *nd1, *nd2, *nd3, *nd4;

    nodes2divide.clear();

    for(ii=0;ii<elemsToRefine.size();ii++)
    {
       nodes2divide.push_back(elemsToRefine[ii]);

       nd = elems[elemsToRefine[ii]];
   
       //cout << nd->GetID() << endl;

       nd4=nd;
       for(bb=0;bb<depth;bb++)
       {
          nd4 = nd4->GetNeighbour(NORTH);
          if(nd4 != NULL)
            nodes2divide.push_back(nd4->GetID());
          else
            break;
       }
       nd4=nd;
       for(bb=0;bb<depth;bb++)
       {
          nd4 = nd4->GetNeighbour(SOUTH);
          if(nd4 != NULL)
            nodes2divide.push_back(nd4->GetID());
          else
            break;
       }

       nd1 = nd;
       nd2 = nd;
       for(aa=0;aa<depth;aa++)
       {
          for(cc=0;cc<2;cc++)
          {
            if(cc==0)
            {
              nd1 = nd1->GetNeighbour(EAST);
              nd3 = nd1;
            }
            if(cc==1)
            {
              nd2 = nd2->GetNeighbour(WEST);
              nd3 = nd2;
            }

            nd4=nd3;
            for(bb=0;bb<depth;bb++)
            {
              nd4 = nd4->GetNeighbour(NORTH);
              if(nd4 != NULL)
                nodes2divide.push_back(nd4->GetID());
              else
                break;
            }
            nd4=nd3;
            for(bb=0;bb<depth;bb++)
            {
              nd4 = nd4->GetNeighbour(SOUTH);
              if(nd4 != NULL)
                nodes2divide.push_back(nd4->GetID());
              else
                break;
            }
          }
       }
    }

    //findUnique(nodes2divide);

    elemsToRefine = nodes2divide;

  return;
}




void  HBSplineBase::addNeighbourElements3D(int depth)
{
    int ii, jj, kk, aa, bb, cc, ee;
    node *nd, *nd1, *nd2, *nd3, *nd4;

    nodes2divide.clear();

    for(ee=0;ee<elemsToRefine.size();ee++)
    {
       nodes2divide.push_back(elemsToRefine[ee]);

       nd = elems[elemsToRefine[ee]];
   
       //cout << nd->GetID() << endl;

       nd1=nd;
       aa=0;bb=0;cc=0;
       for(jj=0;jj<depth;jj++)
       {
         nd1 = nd1->GetNeighbour(BACK);
         if(nd1 != NULL)
         {
           cc++;
           nd1 = nd1->GetNeighbour(SOUTH);
           if(nd1 != NULL)
           {
             bb++;
             nd1 = nd1->GetNeighbour(WEST);
             if(nd1 != NULL)
               aa++;
             else
               break;
           }
           else
             break;
         }
         else
           break;
       }
       
       nd2 = nd1;
       for(kk=0;kk<=(depth+cc);kk++)
       {
         if(nd2 != NULL)
         {
           nd3 = nd2;
           for(jj=0;jj<=(depth+bb);jj++)
           {
             if(nd3 != NULL)
             {
               nd4 = nd3;
               for(ii=0;ii<=(depth+aa);ii++)
               {
                 nd4 = nd4->GetNeighbour(EAST);
                 if(nd4 != NULL)
                   nodes2divide.push_back(nd4->GetID());
                 else
                   break;
               }
             }
             else
               break;

             nd3 = nd3->GetNeighbour(NORTH);
           }
         }
         else
           break;

         nd2 = nd2->GetNeighbour(FRONT);
       }
    }

    //findUnique(nodes2divide);

    elemsToRefine = nodes2divide;

  return;
}



void HBSplineBase::PointBasedRefinement(int kk)
{
  node *nd, *nd1, *nd2;
  int  ii, ee, bb, aa, cc, ll;

  ImmersedIntegrationElement *lme;

  //cout << " hhhhhhhhhhhh " << kk << endl;

  for(bb=0;bb<ImmersedBodyObjects.size();bb++)
  {
    for(aa=0; aa<ImmersedBodyObjects[bb]->GetNumNodes(); aa++)
    {
      //lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

      for(ii=0;ii<DIM;ii++)
        geom[ii] = ImmersedBodyObjects[bb]->GeomData.NodePosCur[aa][ii];
      
      //lme->computePointAtGP(ee, geom);

        cc = findCellNumber(geom);
        nd = elems[cc];

        //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, cell = %5d, \n", geom[0], geom[1], geom[2], cc);

        if(!nd->IsGhost())
        {
          elemsToRefine.push_back(nd->GetID());

          /*
          if(bb == 1)
          {
            nd1 = nd->GetNeighbour(NORTH);
            nd2 = nd->GetNeighbour(SOUTH);
            for(ll=0;ll<4;ll++)
            {
              elemsToRefine.push_back(nd1->GetID());
              elemsToRefine.push_back(nd2->GetID());
              nd1 = nd1->GetNeighbour(NORTH);
              nd2 = nd2->GetNeighbour(SOUTH);
            }
          }
          */
        }
    }
  }
    
  addNeighbourElements(kk);

  return;
}



void  HBSplineBase::Refine(int kk)
{
    //if((int)LevelSetFunc[0][0] == 1)
      //Circle   profile(LevelSetFunc[0][3], LevelSetFunc[0][4], LevelSetFunc[0][5]);
    //else
      //Ellipse  profile(LevelSetFunc[0][3], LevelSetFunc[0][4], LevelSetFunc[0][5], 0.7);
    
    //Circle   profile(LevelSetFunc[0][3], 14.8, LevelSetFunc[0][5]);
    //Circle   profile2(LevelSetFunc[0][3], 15.2, LevelSetFunc[0][5]);
    //Rectangle profile3(9.5, 14.55, 10.5, 15.45);
    //Circle   profile4(LevelSetFunc[0][3], LevelSetFunc[0][4], LevelSetFunc[0][5]);

    //Circle   profile5(LevelSetFunc[0][3], 14.9, LevelSetFunc[0][5]);
    //Circle   profile6(LevelSetFunc[0][3], 15.1, LevelSetFunc[0][5]);
    Circle   profile(10.0,15.0,0.5);


    //Polygon  profile(points);

    int ii, jj, ee, aa, bb, cc;
    double  *tmp1, *tmp2, val[2];
    bool  f1, f2, f3, f4, ff;
    node *nd, *nd1, *nd2, *nd3, *nd4;

    //for(kk=LevelSetFunc[0][2];kk<(2*LevelSetFunc[0][2]+2);kk++)
    //for(kk=CURRENT_LEVEL;kk<=LevelSetFunc[0][2];kk++)
    //for(kk=0;kk<LevelSetFunc[0][2];kk++)

       //cout << " PPPPPPPPPPPPPPP " << endl;
       elemsToRefine.clear();
       for(ee=0;ee<NodeNumsAtLevel[CURRENT_LEVEL].size();ee++)
       {
          nd = elems[NodeNumsAtLevel[CURRENT_LEVEL][ee]];

          if(!nd->IsGhost())
          {

          tmp1 = nd->GetKnots(Dir1);
          tmp2 = nd->GetKnots(Dir2);
          
          //cout << tmp1[0] << '\t' << tmp1[1] << '\t' << tmp2[0] << '\t' << tmp2[1] << endl;

          /*
          ComputeGeometry2D(tmp1[0], tmp2[0], val);
          //f1 = profile.checkPointLocation(val[0], val[1]);
          pt[0] = val[0];
          pt[1] = val[1];
          f1 = (profile.ComputeDistance(pt) > 0.0);
          
          ComputeGeometry2D(tmp1[1], tmp2[0], val);
          //f2 = profile.checkPointLocation(val[0], val[1]);
          pt[0] = val[0];
          pt[1] = val[1];
          f2 = (profile.ComputeDistance(pt) > 0.0);

          ComputeGeometry2D(tmp1[1], tmp2[1], val);
          //f3 = profile.checkPointLocation(val[0], val[1]);
          pt[0] = val[0];
          pt[1] = val[1];
          f3 = (profile.ComputeDistance(pt) > 0.0);

          ComputeGeometry2D(tmp1[0], tmp2[1], val);
          //f4 = profile.checkPointLocation(val[0], val[1]);
          pt[0] = val[0];
          pt[1] = val[1];
          f4 = (profile.ComputeDistance(pt) > 0.0);
          */

          //cout << ee << '\t' << f1 << '\t' << f2 << '\t' << f3 << '\t' << f4 << endl;

          //if( (f1 || f2 || f3 || f4) && !(f1 && f2 && f3 && f4))
            //elemsToRefine.push_back(nd->GetID());

          /*
          ComputeGeometry2D(tmp1[0], tmp2[0], val);
          f1 = (profile.checkPointLocation(val[0], val[1]) || profile2.checkPointLocation(val[0], val[1]) || profile3.checkPointLocation(val[0], val[1]) );//&& profile4.checkPointLocation(val[0], val[1]));
          //f1 = (f1 || profile4.checkPointLocation(val[0], val[1]));

          ComputeGeometry2D(tmp1[1], tmp2[0], val);
          f2 = (profile.checkPointLocation(val[0], val[1]) || profile2.checkPointLocation(val[0], val[1]) || profile3.checkPointLocation(val[0], val[1]) );//&& profile4.checkPointLocation(val[0], val[1]));
          //f2 = (f2 || profile4.checkPointLocation(val[0], val[1]));

          ComputeGeometry2D(tmp1[1], tmp2[1], val);
          f3 = (profile.checkPointLocation(val[0], val[1]) || profile2.checkPointLocation(val[0], val[1]) || profile3.checkPointLocation(val[0], val[1]) );//&& profile4.checkPointLocation(val[0], val[1]));
          //f3 = (f3 || profile4.checkPointLocation(val[0], val[1]));

          ComputeGeometry2D(tmp1[0], tmp2[1], val);
          f4 = (profile.checkPointLocation(val[0], val[1]) || profile2.checkPointLocation(val[0], val[1]) || profile3.checkPointLocation(val[0], val[1]) );//&& profile4.checkPointLocation(val[0], val[1]));
          //f4 = (f4 || profile4.checkPointLocation(val[0], val[1]));

          if( (f1 || f2 || f3 || f4) && !(f1 && f2 && f3 && f4))
            elemsToRefine.push_back(nd->GetID());
          */

          //
          param[0] = tmp1[0];
          param[1] = tmp2[0];
          ComputeGeometry(param, geom);
          f1 = profile.checkPointLocation(geom[0], geom[1]);

          param[0] = tmp1[1];
          param[1] = tmp2[0];
          ComputeGeometry(param, geom);
          f2 = profile.checkPointLocation(geom[0], geom[1]);

          param[0] = tmp1[1];
          param[1] = tmp2[1];
          ComputeGeometry(param, geom);
          f3 = profile.checkPointLocation(geom[0], geom[1]);

          param[0] = tmp1[0];
          param[1] = tmp2[1];
          ComputeGeometry(param, geom);
          f4 = profile.checkPointLocation(geom[0], geom[1]);

          if( (f1 || f2 || f3 || f4) && !(f1 && f2 && f3 && f4))
            elemsToRefine.push_back(nd->GetID());

          /*
          ComputeGeometry2D(tmp1[0], tmp2[0], val);
          f1 = profile5.checkPointLocation(val[0], val[1]);

          ComputeGeometry2D(tmp1[1], tmp2[0], val);
          f2 = profile5.checkPointLocation(val[0], val[1]);

          ComputeGeometry2D(tmp1[1], tmp2[1], val);
          f3 = profile5.checkPointLocation(val[0], val[1]);

          ComputeGeometry2D(tmp1[0], tmp2[1], val);
          f4 = profile5.checkPointLocation(val[0], val[1]);

          if( (f1 || f2 || f3 || f4) && !(f1 && f2 && f3 && f4))
            elemsToRefine.push_back(nd->GetID());
          */

          /*
          //ComputeGeometry2D(tmp1[0], tmp2[0], val);
          f1 = profile6.checkPointLocation(val[0], val[1]);

          //ComputeGeometry2D(tmp1[1], tmp2[0], val);
          f2 = profile6.checkPointLocation(val[0], val[1]);

          //ComputeGeometry2D(tmp1[1], tmp2[1], val);
          f3 = profile6.checkPointLocation(val[0], val[1]);

          //ComputeGeometry2D(tmp1[0], tmp2[1], val);
          f4 = profile6.checkPointLocation(val[0], val[1]);

          if( (f1 || f2 || f3 || f4) && !(f1 && f2 && f3 && f4))
            elemsToRefine.push_back(nd->GetID());
          */
        }
       }

       //cout << " elemsToRefine  " << endl;
       //for(ii=0;ii<elemsToRefine.size();ii++)
         //cout << '\t' << elemsToRefine[ii] ;
       //printf("\n\n\n");
      
       //for(ii=0;ii<elemsToRefine.size();ii++)
       //elems[elemsToRefine[ii]]->printSelf();
       //

  addNeighbourElements(kk);

  return;
}



void HBSplineBase::prepareInteractions()
{
    // go and inherit from ancestors

    Domain::prepareInteractions();

    //printf("     HBSplineBase::prepareInteractions()  .... STARTED ...\n");
    //printf("     NOTHING to be done here ...\n");
    //printf("     HBSplineBase::prepareInteractions()  .... FINISHED ...\n\n");

    return;
}



void HBSplineBase::printInfo()
{
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

    return;
}


void HBSplineBase::findMinMaxX(double *xmn, double *xmx, bool defFlg)
{
    xmn[0] = origin[0];
    xmn[1] = origin[1];

    xmx[0] = origin[0] + gridLEN[0];
    xmx[1] = origin[1] + gridLEN[1];

    return;
}



void  HBSplineBase::createImmersedBoundaryPoints()
{
/*
    int  IBPsize, aa, bb, ee, ind, ii, numPoints, count;

    double  ds = 1.0, fact;

    cout << " AAAAAAAAAAAA " << endl;

    geom.setZero();
    param.setZero();

    Bpoint  *bpoint;

    for(bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
      for(aa=0;aa<ImmersedBodyObjects[bb]->GetNNodes();aa++)
      {
        bpoint = new Bpoint();

        count = aa;

        for(ii=0;ii<DIM;ii++)
        {
          fact = ImmersedPointData[count][ii];
          geom[ii] = fact;
          bpoint->SetPositionOrig(ii, fact);
          bpoint->SetPositionCur(ii, fact);
        }

        //cout << IBpoints[bb]->GetPositionOrig(0) << '\t' << IBpoints[bb]->GetPositionOrig(1) << endl;

        //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, cell = %5d, \n", geom[0], geom[1], geom[2], ee);

        //ee = findCellNumber(geom);

        //bpoint->SetElementNum(ee);

        //bpoint->elem = elems[ee];

        geometryToParametric(geom, param);

        //printf("uu = %12.6f, vv = %12.6f, ww = %12.6f, cell = %5d, \n", param[0], param[1], param[2], ee);

        for(ii=0;ii<DIM;ii++)
          bpoint->SetParam(ii, param[ii]);

        bpoint->SetArcLength(ds);

        //bpoint->SolnData = &(SolnData);

        //bpoint->Setndof2(IBDOF);

        bpoint->findElements();

        bpoint->initialiseDOFvalues();

        bpoint->swapPositions();

        for(ii=0;ii<DIM;ii++)
          bpoint->SetSpecVal(ii, 0.0);

        //ImmersedBodyObjects[bb]->IBpoints.push_back(bpoint);
      }
    }
*/
  return;
}


void HBSplineBase::ProcessBoundaryConditionsRefinedLevels()
{
    assert(CURRENT_LEVEL > 0);

    //cout << "     HBSplineBase: processing boundary conditions for refined levels in 2D ...\n\n";

    int  ii, aa, bb, ind, tmpvec[6][4];
    node*  nd;
    vector<int>  bnodes;

    if(DIM == 1)
    {
      ind = 1;
      tmpvec[0][0] = 0;  // left
      tmpvec[1][0] = 1;  // right
    }
    else if(DIM == 2)
    {
      ind = 2;
      tmpvec[0][0] = 0;    tmpvec[0][1] = 2; // left
      tmpvec[1][0] = 1;    tmpvec[1][1] = 3; // right
      tmpvec[2][0] = 0;    tmpvec[2][1] = 1; // bottom
      tmpvec[3][0] = 2;    tmpvec[3][1] = 3; // top
    }
    else
    {
      ind = 4;
      tmpvec[0][0] = 0; tmpvec[0][1] = 2; tmpvec[0][2] = 4; tmpvec[0][3] = 6; // left face
      tmpvec[1][0] = 1; tmpvec[1][1] = 3; tmpvec[1][2] = 5; tmpvec[1][3] = 7; // right face
      tmpvec[2][0] = 0; tmpvec[2][1] = 1; tmpvec[2][2] = 4; tmpvec[2][3] = 5; // bottom face
      tmpvec[3][0] = 2; tmpvec[3][1] = 3; tmpvec[3][2] = 6; tmpvec[3][3] = 7; // top face
      tmpvec[4][0] = 0; tmpvec[4][1] = 1; tmpvec[4][2] = 2; tmpvec[4][3] = 3; // back face
      tmpvec[5][0] = 4; tmpvec[5][1] = 5; tmpvec[5][2] = 6; tmpvec[5][3] = 7; // front face
    }


    for(bb=0;bb<2*DIM;bb++)
    {
      bnodes.clear();
      for(aa=0;aa<boundaryNodes[bb].size();aa++)
      {
        nd = elems[boundaryNodes[bb][aa]];

        if(nd->IsLeaf())
          bnodes.push_back(nd->GetID());
        else
        {
          for(ii=0;ii<ind;ii++)
            bnodes.push_back(nd->GetChild(tmpvec[bb][ii])->GetID());
        }
      }
      findUnique(bnodes);
      boundaryNodes[bb] = bnodes;
    }

    return;
}


void HBSplineBase::AssignBoundaryConditions()
{
    //cout << "     HBSplineBase: assigning boundary conditions to the boundary elements ...\n\n";

    int  aa, bb, side;

    for(bb=0;bb<DirichletBCs.size();bb++)
    {
      side = (int) (DirichletBCs[bb][0] - 1);

      for(aa=0;aa<boundaryNodes[side].size();aa++)
        elems[boundaryNodes[side][aa]]->DirichletData.push_back(DirichletBCs[bb]);
    }

    for(bb=0;bb<NeumannBCs.size();bb++)
    {
      side = (int) (NeumannBCs[bb][0] - 1);

      for(aa=0;aa<boundaryNodes[side].size();aa++)
        elems[boundaryNodes[side][aa]]->NeumannData.push_back(NeumannBCs[bb]);
    }

    for(bb=0;bb<DerivativeBCs.size();bb++)
    {
      side = (int) (DerivativeBCs[bb][0] - 1);

      for(aa=0;aa<boundaryNodes[side].size();aa++)
        elems[boundaryNodes[side][aa]]->DerivativeBCData.push_back(DerivativeBCs[bb]);
    }

    return;
}




