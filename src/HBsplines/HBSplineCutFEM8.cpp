
#include "HBSplineCutFEM.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "Functions.h"
#include "Files.h"
#include "MyString.h"
#include "myGeomUtilities.h"
#include "myTria.h"
#include "DistFunctions.h"
#include "AABB.h"

#include "insertuniquepoint.h"


extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern Files files;


using namespace myGeom;



void HBSplineCutFEM::plotGeomAdapIntegration1D(int val1, bool flag2, int col, bool PLOT_KNOT_LINES, int* resln)
{
  return;
}


void HBSplineCutFEM::plotGeomAdapIntegration2D(int val1, bool flag2, int col, bool PLOT_KNOT_LINES, int* resln)
{
    int  bb=0, ee=0, ll=0, kk=0, ii=0, totalNGP=0;
    double  volTotal=0.0;

    vtkIdType  pt[50], cellId;

    node *ndTemp;
    myPoint  ptTemp;
    NodeOrientation  neighbour_map1[] = {EAST, NORTH, WEST, WEST};
    NodeOrientation  neighbour_map2[] = {RIGHT, LEFT};

    AABB  bbTemp;

    for(ee=0; ee<activeElements.size(); ee++)
    {
      ndTemp = elems[activeElements[ee]];

      //volTotal += ndTemp->getVolumeGaussPoints(1);

      if( ndTemp->getSubdomainId() == this_mpi_proc )
      {
        //cout << " Node # " << elems[ee]->getID() << '\t' << elems[ee]->isLeaf() << endl;
        //if( !nd->isCutElement() )
        //if( !(ndTemp->getDomainNumber() == -1) )
        if( !ndTemp->isCutElement() )
        {
          if( ndTemp->getDomainNumber() == 0 )
          {
            bbTemp = ndTemp->getAABB();

            pt[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], 0.0);
            pt[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], 0.0);
            pt[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], 0.0);
            pt[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], 0.0);

            for(ll=0;ll<4;ll++)
              quadVTK->GetPointIds()->SetId(ll, pt[ll]);

            cellDataVTK->InsertNextValue( ndTemp->getDomainNumber() );
            cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());
            //cellDataVTK2->InsertNextValue( ndTemp->getLevel() );

            uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
          }
        }
        else                                                // the element is cutCell
        //
        //if( ndTemp->getDomainNumber() == -1 )
        {
          totalNGP += ndTemp->Quadrature.gausspoints.size();

          //
          /////////////////////////////////////////////////////////////
          //  method 1 - adaptive octree
          /////////////////////////////////////////////////////////////

          typedef  typename  node::adapTreePTR  adapTreePTR;

          adapTreePTR  adapNd1, adapNd2, adapNd3;

          adapNd2 = ndTemp->adapIntegNode;

          while( adapNd2 != NULL )
          {
            if( adapNd2->isLeaf() )
            {
              if( (adapNd2->getDomainNumber() == 0) || (adapNd2->getDomainNumber() == -1) )
              //if( adapNd2->getDomainNumber() != 1 ) 
              {
                bbTemp = adapNd2->getAABB();

                pt[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], 0.0);
                pt[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], 0.0);
                pt[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], 0.0);
                pt[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], 0.0);

                for(ll=0;ll<4;ll++)
                  quadVTK->GetPointIds()->SetId(ll, pt[ll]);

                cellDataVTK->InsertNextValue( adapNd2->getDomainNumber() );
                cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());
                uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
              }

              if(adapNd2->getOrientation() == -1)
                adapNd3 = NULL;
              else
                adapNd3 = adapNd2->getNeighbour(neighbour_map1[adapNd2->getOrientation()]);

              while( adapNd2->getOrientation() == NW )
              {
                adapNd2 = adapNd2->getParent();

                if(adapNd2->getOrientation() == -1)
                  adapNd3 = NULL;
                else
      	        {
                  adapNd3 = adapNd2->getNeighbour(neighbour_map1[adapNd2->getOrientation()]);
	            }
              }

              adapNd2 = adapNd3;
            }
            else
            {
              adapNd2 = adapNd2->getChild(0);
            }
          }
          //

        /*
        /////////////////////////////////////////////////////////////
        //  method 2 - adaptive bindarytree
        /////////////////////////////////////////////////////////////

        typedef  typename  node::adapTreePTR  adapTreePTR;

        adapTreePTR  adapNd1, adapNd2, adapNd3;

        adapNd2 = ndTemp->adapIntegNode;

        while( adapNd2 != NULL )
        {
          if( adapNd2->isLeaf() )
          {
            if( (adapNd2->getDomainNumber() == 0) || (adapNd2->getDomainNumber() == -1) )
            {
              bbTemp = adapNd2->getAABB();

              pt[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], 0.0);
              pt[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], 0.0);
              pt[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], 0.0);
              pt[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], 0.0);

              for(ll=0;ll<4;ll++)
                quadVTK->GetPointIds()->SetId(ll, pt[ll]);

              cellDataVTK->InsertNextValue( adapNd2->getDomainNumber() );
              cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());
              //cellDataVTK2->InsertNextValue( adapNd2->getLevel() );
              uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
            }

            if(adapNd2->getOrientation() == -1)
              adapNd3 = NULL;
            else
              adapNd3 = adapNd2->getNeighbour(neighbour_map2[adapNd2->getOrientation()]);

            while( adapNd2->getOrientation() == RIGHT )
            //while(adapNd3 == NULL)
            {
              adapNd2 = adapNd2->getParent();

              if(adapNd2->getOrientation() == -1)
                adapNd3 = NULL;
              else
      	      {
                adapNd3 = adapNd2->getNeighbour(neighbour_map2[adapNd2->getOrientation()]);
	            }
            }

            adapNd2 = adapNd3;
          }
          else
          {
            adapNd2 = adapNd2->getChild(0);
          }
        }
        */
        }                                                   //else
      }
      //
      //cout << " aaaaaaaaaa " << endl;
    } // for(ee=0;ee<elems.size();ee++)

    PetscPrintf(MPI_COMM_WORLD, "\n Total number of Gauss points in cut cells = %d \n\n", totalNGP);
    PetscPrintf(MPI_COMM_WORLD, "\n Total area = %14.12f \n\n", volTotal);

    uGridVTK->SetPoints(pointsVTK);

    //cellDataVTK->SetName("ElemType");
    //cellDataVTK2->SetName("ElemLevel");

    uGridVTK->GetCellData()->SetScalars(cellDataVTK);
    uGridVTK->GetCellData()->AddArray(cellDataVTK2);

    return;
}


//
void HBSplineCutFEM::plotGeomAdapIntegration3D(int val1, bool flag2, int col, bool PLOT_KNOT_LINES, int* resln)
{
    double tstart = MPI_Wtime();

    //vtkSmartPointer<vtkMergePoints>   mergePoints  =  vtkSmartPointer<vtkMergePoints>::New();

    //int ii = nelem[0]*MAX_LEVEL;
    //int jj = nelem[1]*MAX_LEVEL;
    //int kk = nelem[2]*MAX_LEVEL;

    //double  tolTemp = max( max(gridLEN[0], gridLEN[1]), gridLEN[2] )*1.0e-6;

    //mergePoints->SetDivisions(ii, jj, kk);
    //mergePoints->SetTolerance(tolTemp);

    //pointsVTK->InsertNextPoint(origin[0], origin[1], origin[2]);
    //mergePoints->InitPointInsertion(pointsVTK, pointsVTK->GetBounds());


    int ii=0, jj=0, kk=0;
    int ee=0, ll=0, totalNGP=0, ind=0;

    vtkIdType  ptIds[8], id;
    double  vec[]={0.0, 0.0, 0.0};

    node* ndTemp;

    AABB  bbTemp;

    NodeOrientation  neighbour_map1[] = {EAST, NORTH, FRONT, WEST, WEST, WEST, EAST, SOUTH};
    NodeOrientation  neighbour_map2[] = {RIGHT, LEFT};

    for(ee=0; ee<activeElements.size(); ee++)
    {
      ndTemp = elems[activeElements[ee]];

      if( ndTemp->getSubdomainId() == this_mpi_proc )
      {
        if( ndTemp->domNums.size() == 1 )
        {
          if( ndTemp->domNums[0] == 0 )
          {
            bbTemp = ndTemp->getAABB();

            //InsertUniquePointsForAABB3D(mergePoints, ndTemp->getAABB(), ptIds);

            ptIds[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.minBB[2]);
            ptIds[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.minBB[2]);
            ptIds[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]);
            ptIds[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]);

            ptIds[4] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]);
            ptIds[5] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]);
            ptIds[6] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]);
            ptIds[7] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]);

            pointsVTK->InsertNextPoint(origin[0], origin[1], origin[2]);

            for(ll=0;ll<8;ll++)
              hexVTK->GetPointIds()->SetId(ll, ptIds[ll]);

            cellDataVTK->InsertNextValue( ndTemp->getDomainNumber() );
            cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());

            uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
          }
        }
        else if( ndTemp->domNums.size() > 1) // the element is cutCell
        {
          totalNGP += ndTemp->Quadrature.gausspoints.size();

          //
          /////////////////////////////////////////////////////////////
          //  method 1 - adaptive octree
          /////////////////////////////////////////////////////////////

          typedef  typename  node::adapTreePTR  adapTreePTR;

          adapTreePTR  adapNd1, adapNd2, adapNd3;

          adapNd2 = ndTemp->adapIntegNode;

          while( adapNd2 != NULL )
          {
            if( adapNd2->isLeaf() )
            {
              if( (adapNd2->getDomainNumber() == 0) || (adapNd2->getDomainNumber() == -1) )
              {
                //InsertUniquePointsForAABB3D(mergePoints, adapNd2->getAABB(), ptIds);
 
                bbTemp = adapNd2->getAABB();

                ptIds[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.minBB[2]);
                ptIds[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.minBB[2]);
                ptIds[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]);
                ptIds[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]);

                ptIds[4] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]);
                ptIds[5] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]);
                ptIds[6] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]);
                ptIds[7] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]);

                for(ll=0;ll<8;ll++)
                  hexVTK->GetPointIds()->SetId(ll, ptIds[ll]);

                cellDataVTK->InsertNextValue( adapNd2->getDomainNumber() );
                cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());

                uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
              }

              if(adapNd2->getOrientation() == -1)
                adapNd3 = NULL;
              else
                adapNd3 = adapNd2->getNeighbour(neighbour_map1[adapNd2->getOrientation()]);

              while( adapNd2->getOrientation() == SW_FRONT )
              {
                adapNd2 = adapNd2->getParent();

                if(adapNd2->getOrientation() == -1)
                  adapNd3 = NULL;
                else
                {
                  adapNd3 = adapNd2->getNeighbour(neighbour_map1[adapNd2->getOrientation()]);
                }
              }

              adapNd2 = adapNd3;
            }
            else
            {
              adapNd2 = adapNd2->getChild(0);
            }
          }
          //

          /*
          //if(2<1)//////////
          //{
          /////////////////////////////////////////////////////////////
          //  method 2 - adaptive bindarytree
          /////////////////////////////////////////////////////////////
        
          typedef  typename  node::adapTreePTR  adapTreePTR;

          adapTreePTR  adapNd1, adapNd2, adapNd3;
        
          adapNd2 = ndTemp->adapIntegNode;
        
          while( adapNd2 != NULL )
          {
            if( adapNd2->isLeaf() )
            {
              if( (adapNd2->getDomainNumber() == 0) || (adapNd2->getDomainNumber() == -1) )
              {
                //InsertUniquePointsForAABB3D(mergePoints, adapNd2->getAABB(), ptIds);

                bbTemp = adapNd2->getAABB();

                ptIds[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.minBB[2]);
                ptIds[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.minBB[2]);
                ptIds[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]);
                ptIds[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]);

                ptIds[4] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]);
                ptIds[5] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]);
                ptIds[6] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]);
                ptIds[7] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]);

                for(ll=0;ll<8;ll++)
                  hexVTK->GetPointIds()->SetId(ll, ptIds[ll]);

                cellDataVTK->InsertNextValue( adapNd2->getDomainNumber() );
                cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());

                uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
              }

              if(adapNd2->getOrientation() == -1)
                adapNd3 = NULL;
              else
                adapNd3 = adapNd2->getNeighbour(neighbour_map2[adapNd2->getOrientation()]);

              while( adapNd2->getOrientation() == RIGHT )
              {
                adapNd2 = adapNd2->getParent();

                if(adapNd2->getOrientation() == -1)
                  adapNd3 = NULL;
                else
	               {
                  adapNd3 = adapNd2->getNeighbour(neighbour_map2[adapNd2->getOrientation()]);
	               }
              }
              adapNd2 = adapNd3;
            }
            else
            {
              adapNd2 = adapNd2->getChild(0);
            }
          }
          //}/////////////
          */
        } //else
      }
    }

    PetscPrintf(MPI_COMM_WORLD, "\n Total number of Gauss points in cut cells = %d \n\n", totalNGP);

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetCellData()->SetScalars(cellDataVTK);
    uGridVTK->GetCellData()->AddArray(cellDataVTK2); 

    double tend = MPI_Wtime(); 
    PetscPrintf(MPI_COMM_WORLD, "HBSplineCutFEM::plotGeomAdapIntegration3D() took %f millisecond(s) \n ", (tend-tstart)*1000);

    return;
}
//



void  HBSplineCutFEM::postProcessAdapIntegration1D(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    int  ee=0, ii=0, jj=0, kk=0, nlocal=0, index=0;

    nlocal = degree[0]+1;

    VectorXd  N(nlocal), dN_dx(nlocal), d2N_dx2(nlocal), tempVec;
    VectorXd  NN(nlocal), dNN_dx(nlocal), d2NN_dx2(nlocal);
    myPoint  knotIncr, knotBegin, knotEnd;

    double   fact=0, xx[2], val[3], xIntf=0.5;

    vector<double>  uu, uu1, uu2;

    node* nd1;

    AdvDiffExact1D  analy;

    vector<vtkIdType>   pt(resln[0]*3);

    index = 0;
    for(ee=0; ee<activeElements.size(); ee++)
    {
        nd1 = elems[activeElements[ee]];

        knotBegin = nd1->getKnotBegin();
        knotEnd   = nd1->getKnotEnd();
        knotIncr  = nd1->getKnotIncrement();

        if( nd1->isCutElement() )
        {
          fact = (xIntf - knotBegin[0])/resln[0];
          create_vector(knotBegin[0], xIntf, fact, uu1);

          fact = (knotEnd[0] - xIntf)/resln[0];
          create_vector(xIntf, knotEnd[0], fact, uu2);

          kk = uu1.size();
          ii = uu1.size() + uu2.size();
          uu.resize(ii);
          for(ii=0; ii<=resln[0]; ii++)
          {
            uu[ii]    = uu1[ii];
            uu[kk+ii] = uu2[ii];
          }
        }
        else
        {
          fact = (knotEnd[0] - knotBegin[0])/resln[0];
          create_vector(knotBegin[0], knotEnd[0], fact, uu);
        }
          //cout << uu << endl;

          //create the coordinates of the pointsVTK (nodes in FEM)

          for(ii=0;ii<uu.size();ii++)
          {
              param[0] = uu[ii];
              GeomData.computeBasisFunctions1D(knotBegin, knotIncr, param, NN, dNN_dx, d2N_dx2);

              if(nd1->getParent() == NULL)
              {
                N = NN;
                dN_dx = dNN_dx;
                d2N_dx2 = d2NN_dx2;
              }
              else
              {
                N = nd1->SubDivMat*NN;
                dN_dx = nd1->SubDivMat*dNN_dx;
                d2N_dx2 = nd1->SubDivMat*d2NN_dx2;
              }
              
              xx[0] = computeGeometry(0, uu[ii]);

              pt[ii] = pointsVTK->InsertNextPoint(xx[0], 0.0, 0.0);
              
              //if(nd1->getDomainNumber() == 0)
              if(xx[0] <= 0.5)
              {
                val[0] = nd1->computeValue(0, N);
                scaVTK->InsertNextTuple(val);

                val[0] = nd1->computeValue(0, dN_dx);
                scaVTK2->InsertNextTuple(val);
              }
              else
              {
                val[0] = nd1->computeValue2(0, N);
                scaVTK->InsertNextTuple(val);

                val[0] = nd1->computeValue2(0, dN_dx);
                scaVTK2->InsertNextTuple(val);
              }

          } //for(ii=0;ii<uu.n;ii++)

          for(ii=0;ii<uu.size()-1;ii++)
          {
            lineVTK->GetPointIds()->SetId(0, pt[ii]);
            lineVTK->GetPointIds()->SetId(1, pt[ii+1]);

            uGridVTK->InsertNextCell(lineVTK->GetCellType(), lineVTK->GetPointIds());
          } //for(ii=0;ii<uu.n;ii++)
    }

    uGridVTK->SetPoints(pointsVTK);

    scaVTK->SetName("value");
    scaVTK2->SetName("grad");

    //assign nodal coordinates and field data to uGridVTK
    // no need to create lookup table here. All this stuff can be done in Paraview

    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->AddArray(scaVTK2);

//    cout << "   Number of pointsVTK in uGridVTK " << uGridVTK->GetNumberOfPoints()  << endl;

    return;
}




void  HBSplineCutFEM::postProcessAdapIntegration2D(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    int  dd=0, ii=0, jj=0, kk=0, ll=0, count=0, index=0;
    int  ind1=0, ind2=0, e=0, ee=0, gcount=0, ind=0, domTemp=0;
    int  nlocal = (degree[0]+1) * (degree[1] + 1);

    VectorXd  NN(nlocal), N(nlocal), dN_dx(nlocal), dN_dy(nlocal), dNN_dx(nlocal), dNN_dy(nlocal), tempVec, tempVec2, d2N_dx2(nlocal), d2N_dy2(nlocal);
    VectorXd  vectmp(nlocal), rhsTemp;
    myPoint  knotIncr, knotBegin, knotEnd;

    double   fact, uleft, uright, val1;
    NodeOrientation  neighbour_map1[] = {EAST, NORTH, WEST, WEST};
    NodeOrientation  neighbour_map2[] = {RIGHT, LEFT};
    vector<double>  uu, vv;

    vtkIdType pt[50], cellId;

    AABB  bbTemp;

    node* ndTemp;

    if(ndf == 1)
    {
      index = 0;
      for(ee=0; ee<activeElements.size(); ee++)
      {
        ndTemp = elems[activeElements[ee]];

        knotBegin = ndTemp->getKnotBegin();
        knotEnd   = ndTemp->getKnotEnd();
        knotIncr  = ndTemp->getKnotIncrement();

        if( !ndTemp->isCutElement() )
        {
          if( ndTemp->getDomainNumber() == 0 )
          {
            fact = knotIncr[0]/resln[0];
            create_vector(knotBegin[0], knotEnd[0], fact, uu);

            fact = knotIncr[1]/resln[1];
            create_vector(knotBegin[1], knotEnd[1], fact, vv);

            //create the coordinates of the pointsVTK (nodes in FEM)

            count = 0;
            for(jj=0;jj<vv.size();jj++)
            {
              param[1] = vv[jj];
              geom[1] = computeGeometry(1, vv[jj]);

              for(ii=0;ii<uu.size();ii++)
              {
                param[0] = uu[ii];
                geom[0] = computeGeometry(0, uu[ii]);

                pt[count] = pointsVTK->InsertNextPoint(geom[0], geom[1], 0.0);

                GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

                if(ndTemp->getParent() == NULL)
                  N = NN;
                else
                  N = ndTemp->SubDivMat*NN;

                if(domTemp == 0)
                  fact = ndTemp->computeValue(0, N);
                else
                  fact = ndTemp->computeValue2(0, N);

                scaVTK->InsertNextValue(fact);
                scaVTK2->InsertNextValue(fact - GeomData.analyDBC->computeValue(0, geom[0], geom[1]) );

                count++;
              } //for(ii=0;ii<uu.size();ii++)
            } //for(jj=0;jj<vv.size();jj++)

            quadVTK->GetPointIds()->SetId(0, pt[0]);
            quadVTK->GetPointIds()->SetId(1, pt[1]);
            quadVTK->GetPointIds()->SetId(2, pt[3]);
            quadVTK->GetPointIds()->SetId(3, pt[2]);

            uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());

            cellDataVTK->InsertNextValue(ndTemp->getDomainNumber() );
            cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());
	        }
        } //if( !nd->isCutElement() )
        else // the element is cutCell
        {
          vtkSmartPointer<vtkTriangle> triaVTK =  vtkSmartPointer<vtkTriangle>::New();

          myPoly *poly;

          for(ii=0; ii<ndTemp->subTrias.size(); ii++)
          {
            poly = ndTemp->subTrias[ii];

            if( poly->getDomainNumber() == 0 )
            {
              for(kk=0; kk<3; kk++)
              {
                geom = poly->GetPoint(kk);
                //cout << kk << '\t' << ptTemp[0] << '\t' << ptTemp[1] << endl;

                pt[kk] = pointsVTK->InsertNextPoint(geom[0], geom[1], 0.0);

                geometryToParametric(geom, param);
                GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

                if(ndTemp->getParent() == NULL)
                  N = NN;
                else
                  N = ndTemp->SubDivMat*NN;

                if( domTemp )
                  fact = ndTemp->computeValue2(0, N);
                else
                  fact = ndTemp->computeValue(0, N);

                scaVTK->InsertNextValue(fact);
                scaVTK2->InsertNextValue(fact - GeomData.analyDBC->computeValue(0, geom[0], geom[1]) );

                triaVTK->GetPointIds()->SetId(kk, pt[kk] );
              }

              cellDataVTK->InsertNextValue(poly->getDomainNumber());
              cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());

              uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
            }
          } //  for(ii=0; ii<nd->subTrias.size(); ii++)
          //cout << " AAAAAAAAAA " << endl;
        } // else
    }

    scaVTK->SetName("value");
    scaVTK2->SetName("force");

    uGridVTK->SetPoints(pointsVTK);

    //assign nodal coordinates and field data to uGridVTK
    // no need to create lookup table here. All this stuff can be done in Paraview

    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->AddArray(scaVTK2);

    uGridVTK->GetCellData()->SetScalars(cellDataVTK);
    uGridVTK->GetCellData()->AddArray(cellDataVTK2);

    // create a write object and write uGridVTK to it
  }
  else // for Stokes and Navier-Stokes
  {
    double vec[3]={0.0, 0.0, 0.0}; //    vec[0] = vec[1] = vec[2] = 0.0;

    vecVTK->SetNumberOfComponents(3);
    //vecVTK->SetNumberOfTuples(count);
    vecVTK2->SetNumberOfComponents(3);
    //vecVTK2->SetNumberOfTuples(count);
    //scaVTK->SetNumberOfTuples(count);
    //scaVTK2->SetNumberOfTuples(count);

    //cout << " Node aaaaaaaaaaa " << endl;

    for(ee=0; ee<activeElements.size(); ee++)
    {
      ndTemp = elems[activeElements[ee]];
      //cout << " Node # " << nd->getID() << endl;

      if( ndTemp->getSubdomainId() == this_mpi_proc )
      {
        knotBegin = ndTemp->getKnotBegin();
        knotEnd   = ndTemp->getKnotEnd();
        knotIncr  = ndTemp->getKnotIncrement();

        //if( !nd->isCutElement() )
        if( ndTemp->getDomainNumber() == 0 )
        {
            fact = knotIncr[0]/resln[0];
            create_vector(knotBegin[0], knotEnd[0], fact, uu);

            fact = knotIncr[1]/resln[1];
            create_vector(knotBegin[1], knotEnd[1], fact, vv);

            //create the coordinates of the pointsVTK (nodes in FEM)

            count = 0;
            for(jj=0;jj<vv.size();jj++)
            {
              param[1] = vv[jj];
              geom[1] = computeGeometry(1, vv[jj]);
              for(ii=0;ii<uu.size();ii++)
              {
                param[0] = uu[ii];
                geom[0] = computeGeometry(0, uu[ii]);

                pt[count++] = pointsVTK->InsertNextPoint(geom[0], geom[1], 0.0);

                GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

                if(ndTemp->getParent() == NULL)
                {
                  N = NN;
                  dN_dx = dNN_dx;
                  dN_dy = dNN_dy;
                }
                else
                {
                  N = ndTemp->SubDivMat*NN;
                  dN_dx = ndTemp->SubDivMat*dNN_dx;
                  dN_dy = ndTemp->SubDivMat*dNN_dy;
                }

                vec[0] = ndTemp->computeValue(0, N);
                vec[1] = ndTemp->computeValue(1, N);
                fact   = ndTemp->computeValue(2, N);

                vecVTK->InsertNextTuple(vec);
                //vecVTK2->InsertNextTuple(vec);
                scaVTK->InsertNextValue(fact);
                fact = ndTemp->computeValue(1, dN_dx) - ndTemp->computeValue(0, dN_dy);
                scaVTK2->InsertNextValue(fact);

              } // for(ii=0;ii<uu.size();ii++)
            } // for(jj=0;jj<vv.size();jj++)

            quadVTK->GetPointIds()->SetId(0, pt[0]);
            quadVTK->GetPointIds()->SetId(1, pt[1]);
            quadVTK->GetPointIds()->SetId(2, pt[3]);
            quadVTK->GetPointIds()->SetId(3, pt[2]);

            uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());

            //cellDataVTK->InsertNextValue(0);
            //cellDataVTK2->InsertNextValue(0);

            cellDataVTK->InsertNextValue(ndTemp->getDomainNumber() );
            cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());

            //cout << " ooooooooooooo " << endl;
        } //if( !nd->isCutElement() )

        if( ndTemp->getDomainNumber() == -1 )
        {
          //
          /////////////////////////////////////////////////////////////
          //  method 1 - octree
          /////////////////////////////////////////////////////////////

          typedef  typename  node::adapTreePTR  adapTreePTR;

          adapTreePTR  adapNd1, adapNd2, adapNd3;

          adapNd2 = ndTemp->adapIntegNode;

          while( adapNd2 != NULL )
          {
            if( adapNd2->isLeaf() )
            {
              if( (adapNd2->getDomainNumber() == 0) || (adapNd2->getDomainNumber() == -1) )
              //if( (adapNd2->getDomainNumber() <= 0) )
              {
                bbTemp = adapNd2->getAABB();

                pt[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], 0.0);
                pt[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], 0.0);
                pt[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], 0.0);
                pt[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], 0.0);

                for(ll=0;ll<4;ll++)
                  quadVTK->GetPointIds()->SetId(ll, pt[ll]);

                cellDataVTK->InsertNextValue( adapNd2->getDomainNumber() );
                //cellDataVTK2->InsertNextValue( adapNd2->getLevel() );
                cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());

                uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());

                knotBegin = ndTemp->getKnotBegin();
                knotEnd   = ndTemp->getKnotEnd();
                knotIncr  = ndTemp->getKnotIncrement();

                fact = knotIncr[0]/resln[0];
                create_vector(knotBegin[0], knotEnd[0], fact, uu);

                fact = knotIncr[1]/resln[1];
                create_vector(knotBegin[1], knotEnd[1], fact, vv);

                for(jj=0;jj<vv.size();jj++)
                {
                  param[1] = vv[jj];
                  for(ii=0;ii<uu.size();ii++)
                  {
                    param[0] = uu[ii];

                    GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

                    if(ndTemp->getParent() == NULL)
                    {
                      N = NN;
                      dN_dx = dNN_dx;
                      dN_dy = dNN_dy;
                    }
                    else
                    {
                      N = ndTemp->SubDivMat*NN;
                      dN_dx = ndTemp->SubDivMat*dNN_dx;
                      dN_dy = ndTemp->SubDivMat*dNN_dy;
                    }

                    vec[0] = ndTemp->computeValue(0, N);
                    vec[1] = ndTemp->computeValue(1, N);
                    fact   = ndTemp->computeValue(2, N);

                    vecVTK->InsertNextTuple(vec);
                    //vecVTK2->InsertNextTuple(vec);
                    scaVTK->InsertNextValue(fact);
                    fact = ndTemp->computeValue(1, dN_dx) - ndTemp->computeValue(0, dN_dy);
                    scaVTK2->InsertNextValue(fact);
                  } // for(ii=0;ii<uu.size();ii++)
                } // for(jj=0;jj<vv.size();jj++)
              } // if( (adapNd2->getDomainNumber() <= 0) )

              if(adapNd2->getOrientation() == -1)
                adapNd3 = NULL;
              else
                adapNd3 = adapNd2->getNeighbour(neighbour_map1[adapNd2->getOrientation()]);

              //cout <<  adapNd2->getID() << '\t' << adapNd2->getOrientation() << endl;
              while( adapNd2->getOrientation() == NW )
              {
                adapNd2 = adapNd2->getParent();

                if(adapNd2->getOrientation() == -1)
                  adapNd3 = NULL;
                else
                {
                  adapNd3 = adapNd2->getNeighbour(neighbour_map1[adapNd2->getOrientation()]);
                }
              } //while( adapNd2->getOrientation() == NW )

              adapNd2 = adapNd3;
            } // if( adapNd2->isLeaf() )
            else
            {
              adapNd2 = adapNd2->getChild(0);
            }
          } // while( adapNd2 != NULL )
          //
          //cout << " AAAAAAAAAA " << endl;

          /*
          /////////////////////////////////////////////////////////////
          //  method 2 - bindarytree
          /////////////////////////////////////////////////////////////

          typedef  typename  node::adapTreePTR  adapTreePTR;

          adapTreePTR  adapNd1, adapNd2, adapNd3;

          adapNd2 = ndTemp->adapIntegNode;

          while( adapNd2 != NULL )
          {
            if( adapNd2->isLeaf() )
            {
              if( (adapNd2->getDomainNumber() <= 0) )
              {
                bbTemp = adapNd2->getAABB();

                pt[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], 0.0);
                pt[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], 0.0);
                pt[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], 0.0);
                pt[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], 0.0);

                for(ll=0;ll<4;ll++)
                  quadVTK->GetPointIds()->SetId(ll, pt[ll]);

                cellDataVTK->InsertNextValue( adapNd2->getDomainNumber() );
                //cellDataVTK2->InsertNextValue( adapNd2->getLevel() );
                cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());

                uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());

                knotBegin = ndTemp->getKnotBegin();
                knotEnd   = ndTemp->getKnotEnd();
                knotIncr  = ndTemp->getKnotIncrement();

                fact = knotIncr[0]/resln[0];
                create_vector(knotBegin[0], knotEnd[0], fact, uu);

                fact = knotIncr[1]/resln[1];
                create_vector(knotBegin[1], knotEnd[1], fact, vv);

                for(jj=0;jj<vv.size();jj++)
                {
                  param[1] = vv[jj];
                  for(ii=0;ii<uu.size();ii++)
                  {
                    param[0] = uu[ii];

                    GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

                    if(ndTemp->getParent() == NULL)
                    {
                      N = NN;
                      dN_dx = dNN_dx;
                      dN_dy = dNN_dy;
                    }
                    else
                    {
                      N = ndTemp->SubDivMat*NN;
                      dN_dx = ndTemp->SubDivMat*dNN_dx;
                      dN_dy = ndTemp->SubDivMat*dNN_dy;
                    }

                    vec[0] = ndTemp->computeValue(0, N);
                    vec[1] = ndTemp->computeValue(1, N);
                    fact   = ndTemp->computeValue(2, N);

                    vecVTK->InsertNextTuple(vec);
                    //vecVTK2->InsertNextTuple(vec);
                    scaVTK->InsertNextValue(fact);
                    fact = ndTemp->computeValue(1, dN_dx) - ndTemp->computeValue(0, dN_dy);
                    scaVTK2->InsertNextValue(fact);
                  } // for(ii=0;ii<uu.size();ii++)
                } // for(jj=0;jj<vv.size();jj++)
              } // if( (adapNd2->getDomainNumber() <= 0) )

              if(adapNd2->getOrientation() == -1)
                adapNd3 = NULL;
              else
                adapNd3 = adapNd2->getNeighbour(neighbour_map2[adapNd2->getOrientation()]);

              //cout <<  adapNd2->getID() << '\t' << adapNd2->getOrientation() << endl;
              while( adapNd2->getOrientation() == RIGHT )
              {
                adapNd2 = adapNd2->getParent();

                if(adapNd2->getOrientation() == -1)
                  adapNd3 = NULL;
                else
      	        {
                  adapNd3 = adapNd2->getNeighbour(neighbour_map2[adapNd2->getOrientation()]);
	              }
              } //while( adapNd2->getOrientation() == RIGHT )

              adapNd2 = adapNd3;
            } // if( adapNd2->isLeaf() )
            else
            {
              adapNd2 = adapNd2->getChild(0);
            }
          } // while( adapNd2 != NULL )
          */          //cout << " AAAAAAAAAA " << endl;
        } // else
      }
    }

      //cout << " jjjjjjjjjjjjjjjjjj " << endl;
      vecVTK->SetName("vel");
      //vecVTK2->SetName("force");
      scaVTK->SetName("pres");
      scaVTK2->SetName("vortz");

      //assign nodal coordinates and field data to uGridVTK
      // no need to create lookup table here. All this stuff can be done in Paraview

      uGridVTK->SetPoints(pointsVTK);

      uGridVTK->GetPointData()->SetScalars(scaVTK);
      uGridVTK->GetPointData()->SetVectors(vecVTK);
      //uGridVTK->GetPointData()->AddArray(vecVTK2);
      uGridVTK->GetPointData()->AddArray(scaVTK2);

      uGridVTK->GetCellData()->SetScalars(cellDataVTK);
      uGridVTK->GetCellData()->AddArray(cellDataVTK2);

      //cout << " jjjjjjjjjjjjjjjjjj " << endl;
  }


    return;
}
//




void  HBSplineCutFEM::postProcessAdapIntegration3D(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    //cout << " jjjjjjjjjjjjjjjjjj " << endl;

    //vtkSmartPointer<vtkMergePoints>   mergePoints  =  vtkSmartPointer<vtkMergePoints>::New();
    //vtkSmartPointer<vtkPoints>        pointsVTK    = vtkSmartPointer<vtkPoints>::New();
    // using mergePoints algorithm takes 45 seconds for the Level-1 mesh for the 3D cylinder problem
    // where as without this algorithm, it only takes 0.5 seconds.

    //int ii = nelem[0]*MAX_LEVEL;
    //int jj = nelem[1]*MAX_LEVEL;
    //int kk = nelem[2]*MAX_LEVEL;

    //double  tolTemp = max( max(gridLEN[0], gridLEN[1]), gridLEN[2] )*1.0e-6;

    //mergePoints->SetDivisions(ii, jj, kk);
    //mergePoints->SetTolerance(tolTemp);

    //pointsVTK->Reset();

    //pointsVTK->InsertNextPoint(origin[0], origin[1], origin[2]);
    //mergePoints->InitPointInsertion(pointsVTK, pointsVTK->GetBounds());

    int  ii=0, jj=0, kk=0;
    int  count=0, ee=0, domTemp=0;
    vtkIdType ptIds[50], ptExists, id;

    int  nlocal = (degree[0]+1) * (degree[1] + 1) * (degree[2] + 1);

    VectorXd  NN(nlocal), dNN_dx(nlocal), dNN_dy(nlocal), dNN_dz(nlocal);
    VectorXd  N(nlocal), dN_dx(nlocal), dN_dy(nlocal), dN_dz(nlocal);
    VectorXd  vectmp(nlocal), rhsTemp;
    myPoint   knotBegin, knotEnd, knotIncr;
    myPoint   knotBeginTemp, knotEndTemp, knotIncrTemp;

    double   pres, uu[2], vv[2], ww[2];

    NodeOrientation  neighbour_map1[] = {EAST, NORTH, FRONT, WEST, WEST, WEST, EAST, SOUTH};
    NodeOrientation  neighbour_map2[] = {RIGHT, LEFT};

    node* ndTemp;


  //if(ndf == 1)
  //{
  //}
  //else // for Stokes and Navier-Stokes
  //{
      double vec[]={0.0, 0.0, 0.0};

      vecVTK->SetName("vel");
      vecVTK2->SetName("vort");
      scaVTK->SetName("pres");

      vecVTK->SetNumberOfComponents(3);
      //vecVTK->SetNumberOfTuples(count);
      vecVTK2->SetNumberOfComponents(3);
      //vecVTK2->SetNumberOfTuples(count);
      //scaVTK->SetNumberOfTuples(count);
      //scaVTK2->SetNumberOfTuples(count);

      for(ee=0; ee<activeElements.size(); ee++)
      {
        ndTemp = elems[activeElements[ee]];

        if( ndTemp->getSubdomainId() == this_mpi_proc )
        {
          knotBegin = ndTemp->getKnotBegin();
          knotEnd   = ndTemp->getKnotEnd();
          knotIncr  = ndTemp->getKnotIncrement();

          if( (ndTemp->domNums.size() == 1) && (ndTemp->domNums[0] == 0) )
          {
            uu[0] = knotBegin[0];  uu[1] = knotEnd[0];
            vv[0] = knotBegin[1];  vv[1] = knotEnd[1];
            ww[0] = knotBegin[2];  ww[1] = knotEnd[2];

            count = 0;
            for(kk=0; kk<2; kk++)
            {
                param[2] = ww[kk];
                geom[2] = computeGeometry(2, ww[kk]);

            for(jj=0; jj<2; jj++)
            {
                param[1] = vv[jj];
                geom[1] = computeGeometry(1, vv[jj]);

            for(ii=0; ii<2; ii++)
            {
                param[0] = uu[ii];
                geom[0] = computeGeometry(0, uu[ii]);

                ptIds[count++] = pointsVTK->InsertNextPoint(geom[0], geom[1], geom[2]);

                //ptExists = mergePoints->InsertUniquePoint(&geom[0], id);
                //ptIds[count++] = id;

                //if(ptExists == 1)
                //{
                  if(ndTemp->getParent() == NULL)
                  {
                    GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, N, dN_dx, dN_dy, dN_dz);
                  }
                  else
                  {
                    GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

                    N = ndTemp->SubDivMat*NN;
                    dN_dx = ndTemp->SubDivMat*dNN_dx;
                    dN_dy = ndTemp->SubDivMat*dNN_dy;
                    dN_dz = ndTemp->SubDivMat*dNN_dz;
                  }

                  vec[0] = ndTemp->computeValue(0, N);
                  vec[1] = ndTemp->computeValue(1, N);
                  vec[2] = ndTemp->computeValue(2, N);

                  pres   = ndTemp->computeValue(3, N);

                  vecVTK->InsertNextTuple(vec);
                  vecVTK2->InsertNextTuple(vec);
                  scaVTK->InsertNextValue(pres);
                //} //if(ptExists == 1)
            } // for(ii=0;ii<uu.size();ii++)
            } // for(jj=0;jj<vv.size();jj++)
            } // for(kk=0; kk<ww.size(); kk++)

            hexVTK->GetPointIds()->SetId(0, ptIds[0]);
            hexVTK->GetPointIds()->SetId(1, ptIds[1]);
            hexVTK->GetPointIds()->SetId(2, ptIds[3]);
            hexVTK->GetPointIds()->SetId(3, ptIds[2]);

            hexVTK->GetPointIds()->SetId(4, ptIds[4]);
            hexVTK->GetPointIds()->SetId(5, ptIds[5]);
            hexVTK->GetPointIds()->SetId(6, ptIds[7]);
            hexVTK->GetPointIds()->SetId(7, ptIds[6]);

            uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());

            cellDataVTK->InsertNextValue(ndTemp->getDomainNumber() );
            cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());

            //cout << " ooooooooooooo " << endl;
        } // if( (ndTemp->domNums.size() == 0) && (ndTemp->domNums[0] == 0) )
        else if( ndTemp->domNums.size() > 1 )
        {
          /////////////////////////////////////////////////////////////
          //  method 2
          /////////////////////////////////////////////////////////////
        
          typedef  typename  node::adapTreePTR  adapTreePTR;

          adapTreePTR  adapNd1, adapNd2, adapNd3;
        
          adapNd2 = ndTemp->adapIntegNode;

          while( adapNd2 != NULL )
          {
            if( adapNd2->isLeaf()  )
            {
              if( adapNd2->getDomainNumber() <= 0 && (adapNd2->getLevel() < 3) )
              {
                knotBeginTemp = adapNd2->getKnotBegin();
                knotEndTemp   = adapNd2->getKnotEnd();
                knotIncrTemp  = adapNd2->getKnotIncrement();

                uu[0] = knotBeginTemp[0];  uu[1] = knotEndTemp[0];
                vv[0] = knotBeginTemp[1];  vv[1] = knotEndTemp[1];
                ww[0] = knotBeginTemp[2];  ww[1] = knotEndTemp[2];

                count = 0;
                for(kk=0; kk<2; kk++)
                {
                  param[2] = ww[kk];
                  geom[2] = computeGeometry(2, ww[kk]);

                for(jj=0;jj<2;jj++)
                {
                  param[1] = vv[jj];
                  geom[1] = computeGeometry(1, vv[jj]);

                for(ii=0;ii<2;ii++)
                {
                  param[0] = uu[ii];
                  geom[0] = computeGeometry(0, uu[ii]);

                  ptIds[count++] = pointsVTK->InsertNextPoint(geom[0], geom[1], geom[2]);

                  //ptExists = mergePoints->InsertUniquePoint(&geom[0], id);
                  //ptIds[count++] = id;

                  //if(ptExists == 1)
                  //{
                    if(ndTemp->getParent() == NULL)
                    {
                      GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, N, dN_dx, dN_dy, dN_dz);
                    }
                    else
                    {
                      GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

                      N = ndTemp->SubDivMat*NN;
                      dN_dx = ndTemp->SubDivMat*dNN_dx;
                      dN_dy = ndTemp->SubDivMat*dNN_dy;
                      dN_dz = ndTemp->SubDivMat*dNN_dz;
                    }

                    vec[0] = ndTemp->computeValue(0, N);
                    vec[1] = ndTemp->computeValue(1, N);
                    vec[2] = ndTemp->computeValue(2, N);

                    pres   = ndTemp->computeValue(3, N);

                    vecVTK->InsertNextTuple(vec);
                    vecVTK2->InsertNextTuple(vec);
                    scaVTK->InsertNextValue(pres);
                  //} //if(ptExists == 1)
                } // for(ii=0;ii<uu.size();ii++)
                } // for(jj=0;jj<vv.size();jj++)
                } // for(kk=0; kk<ww.size(); kk++)

                hexVTK->GetPointIds()->SetId(0, ptIds[0]);
                hexVTK->GetPointIds()->SetId(1, ptIds[1]);
                hexVTK->GetPointIds()->SetId(2, ptIds[3]);
                hexVTK->GetPointIds()->SetId(3, ptIds[2]);

                hexVTK->GetPointIds()->SetId(4, ptIds[4]);
                hexVTK->GetPointIds()->SetId(5, ptIds[5]);
                hexVTK->GetPointIds()->SetId(6, ptIds[7]);
                hexVTK->GetPointIds()->SetId(7, ptIds[6]);

                uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());

                cellDataVTK->InsertNextValue(ndTemp->getDomainNumber());
                cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());

            } // if( (adapNd2->getDomainNumber() <= 0) )

            if(adapNd2->getOrientation() == -1)
              adapNd3 = NULL;
            else
              adapNd3 = adapNd2->getNeighbour(neighbour_map1[adapNd2->getOrientation()]);

            while( adapNd2->getOrientation() == SW_FRONT )
            {
              adapNd2 = adapNd2->getParent();

              if(adapNd2->getOrientation() == -1)
                adapNd3 = NULL;
              else
              {
                adapNd3 = adapNd2->getNeighbour(neighbour_map1[adapNd2->getOrientation()]);
              }
            }
            adapNd2 = adapNd3;
          }
          else
          {
            adapNd2 = adapNd2->getChild(0);
          }
          //cout << " AAAAAAAAAA " << endl;
	      }// else
        }
      } // for(ee=0; ee<activeElements.size(); ee++)

      uGridVTK->SetPoints(pointsVTK);

      uGridVTK->GetPointData()->SetScalars(scaVTK);
      uGridVTK->GetPointData()->SetVectors(vecVTK);
      //uGridVTK->GetPointData()->AddArray(vecVTK2);
      //uGridVTK->GetPointData()->AddArray(scaVTK2);

      //uGridVTK->GetCellData()->SetScalars(cellDataVTK);
      //uGridVTK->GetCellData()->AddArray(cellDataVTK2);

      //cout << " jjjjjjjjjjjjjjjjjj " << endl;
  }

  return;
}



