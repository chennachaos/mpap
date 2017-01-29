
#include "HBSplineCutFEM.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "PlotVTK.h"
#include "Functions.h"
#include "Files.h"
#include "MyString.h"
#include "headersVTK.h"
#include "myGeomUtilities.h"
#include "NurbsUtilities.h"
#include "myTria.h"
#include "DistFunctions.h"
#include "AABB.h"


extern PlotVTK  plotvtk;
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
    int  bb, ee, ll, kk, ii, totalNGP=0;
    double  volTotal=0.0;

    vtkIdType  pt[50], cellId;

    node *ndTemp;
    myPoint  ptTemp;
    NodeOrientation  neigbour_map1[] = {EAST, NORTH, WEST, WEST};
    NodeOrientation  neigbour_map2[] = {RIGHT, LEFT};

    AABB  bbTemp;


    for(ee=0; ee<activeElements.size(); ee++)
    {
      ndTemp = elems[activeElements[ee]];

      //volTotal += ndTemp->GetVolumeGaussPoints(1);

      //cout << " Node # " << elems[ee]->GetID() << '\t' << elems[ee]->IsLeaf() << endl;
      //if( !nd->IsCutElement() )
      if( !(ndTemp->GetDomainNumber() == -1) )
      {
        if( ndTemp->GetDomainNumber() == 0 )
        {
          bbTemp = ndTemp->GetAABB();

          pt[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], 0.0);
          pt[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], 0.0);
          pt[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], 0.0);
          pt[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], 0.0);

          for(ll=0;ll<4;ll++)
            quadVTK->GetPointIds()->SetId(ll, pt[ll]);

          cellDataVTK->InsertNextValue( ndTemp->GetDomainNumber() );
          cellDataVTK2->InsertNextValue(ndTemp->get_subdomain_id());
          //cellDataVTK2->InsertNextValue( ndTemp->GetLevel() );

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        }
      }
      else // the element is cutCell
      //
      //if( ndTemp->GetDomainNumber() == -1 )
      {
        totalNGP += ndTemp->Quadrature.gausspoints.size();
        
        /////////////////////////////////////////////////////////////
        //  method 1 - adaptive octree
        /////////////////////////////////////////////////////////////
        
        /*
        //AdaptiveOctree<2>  *adapNd1, *adapNd2, *adapNd3;
        //adapNd2 = ndTemp->adapIntegNode1;

        typedef  typename  node::adapTreePTR  adapTreePTR;

        adapTreePTR  adapNd1, adapNd2, adapNd3;

        adapNd2 = ndTemp->adapIntegNode;

        while( adapNd2 != NULL )
        {
          if( adapNd2->IsLeaf() )
          {
            //if( (adapNd2->GetDomainNumber() == 0) || (adapNd2->GetDomainNumber() == -1) )
            if( adapNd2->GetDomainNumber() != 1 ) 
            {
              bbTemp = adapNd2->GetAABB();

              pt[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], 0.0);
              pt[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], 0.0);
              pt[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], 0.0);
              pt[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], 0.0);

              for(ll=0;ll<4;ll++)
                quadVTK->GetPointIds()->SetId(ll, pt[ll]);

              cellDataVTK->InsertNextValue( adapNd2->GetDomainNumber() );
              uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
            }

            if(adapNd2->GetOrientation() == -1)
              adapNd3 = NULL;
            else
              adapNd3 = adapNd2->GetNeighbour(neigbour_map1[adapNd2->GetOrientation()]);

            //cout <<  adapNd2->GetID() << '\t' << adapNd2->GetOrientation() << endl;
            while( adapNd2->GetOrientation() == NW )
            //while(adapNd3 == NULL)
            {
              //cout << " inside IF " << endl;
              adapNd2 = adapNd2->GetParent();
              //cout << " inside IF " << endl;

              if(adapNd2->GetOrientation() == -1)
                adapNd3 = NULL;
              else
	      {
                //cout <<  adapNd2->GetID() << '\t' << adapNd2->GetOrientation() << endl;
                adapNd3 = adapNd2->GetNeighbour(neigbour_map1[adapNd2->GetOrientation()]);
	      }

              //cout << adapNd2->GetID() << endl;
            }

            adapNd2 = adapNd3;
          }
          else
          {
            adapNd2 = adapNd2->GetChild(0);
          }
        }
        */

        //
        /////////////////////////////////////////////////////////////
        //  method 2 - adaptive bindarytree
        /////////////////////////////////////////////////////////////
        
        typedef  typename  node::adapTreePTR  adapTreePTR;
        
        adapTreePTR  adapNd1, adapNd2, adapNd3;
        
        adapNd2 = ndTemp->adapIntegNode;

        while( adapNd2 != NULL )
        {
          if( adapNd2->IsLeaf() )
          {
            //if( (adapNd2->GetDomainNumber() == 0) || (adapNd2->GetDomainNumber() == -1) )
            //{
              bbTemp = adapNd2->GetAABB();

              pt[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], 0.0);
              pt[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], 0.0);
              pt[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], 0.0);
              pt[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], 0.0);

              for(ll=0;ll<4;ll++)
                quadVTK->GetPointIds()->SetId(ll, pt[ll]);

              cellDataVTK->InsertNextValue( adapNd2->GetDomainNumber() );
              cellDataVTK2->InsertNextValue(ndTemp->get_subdomain_id());
              //cellDataVTK2->InsertNextValue( adapNd2->GetLevel() );
              uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
            //}

            if(adapNd2->GetOrientation() == -1)
              adapNd3 = NULL;
            else
              adapNd3 = adapNd2->GetNeighbour(neigbour_map2[adapNd2->GetOrientation()]);

            //cout <<  adapNd2->GetID() << '\t' << adapNd2->GetOrientation() << endl;
            while( adapNd2->GetOrientation() == RIGHT )
            //while(adapNd3 == NULL)
            {
              //cout << " inside IF " << endl;
              adapNd2 = adapNd2->GetParent();
              //cout << " inside IF " << endl;

              if(adapNd2->GetOrientation() == -1)
                adapNd3 = NULL;
              else
	      {
                //cout <<  adapNd2->GetID() << '\t' << adapNd2->GetOrientation() << endl;
                adapNd3 = adapNd2->GetNeighbour(neigbour_map2[adapNd2->GetOrientation()]);
	      }

              //cout << adapNd2->GetID() << endl;
            }

            adapNd2 = adapNd3;
          }
          else
          {
            adapNd2 = adapNd2->GetChild(0);
          }
        }
        //

      } //else
      //
      //cout << " aaaaaaaaaa " << endl;
    } // for(ee=0;ee<elems.size();ee++)

    PetscPrintf(MPI_COMM_WORLD, "\n Total number of Gauss points in cut cells = %d \n\n", totalNGP);
    PetscPrintf(MPI_COMM_WORLD, "\n Total area = %12.8f \n\n", volTotal);

    uGridVTK->SetPoints(pointsVTK);

    //cellDataVTK->SetName("ElemType");
    //cellDataVTK2->SetName("ElemLevel");

    uGridVTK->GetCellData()->SetScalars(cellDataVTK);
    uGridVTK->GetCellData()->AddArray(cellDataVTK2);

    return;
}




void HBSplineCutFEM::plotGeomAdapIntegration3D(int val1, bool flag2, int col, bool PLOT_KNOT_LINES, int* resln)
{
    int  ii, ee, ll, kk, totalNGP=0;

    vtkIdType  ptIds[8];

    time_t tstart, tend;
    //tstart = time(0);
    
    node* ndTemp;

    AABB  bbTemp;

    NodeOrientation  neigbour_map1[] = {EAST, NORTH, FRONT, WEST, WEST, WEST, EAST, SOUTH};
    NodeOrientation  neigbour_map2[] = {RIGHT, LEFT};

    for(ee=0; ee<activeElements.size(); ee++)
    {
      ndTemp = elems[activeElements[ee]];

      if( ndTemp->domNums.size() == 1 )
      {
        bbTemp = ndTemp->GetAABB();

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

        cellDataVTK->InsertNextValue( ndTemp->GetDomainNumber() );
        cellDataVTK2->InsertNextValue(ndTemp->get_subdomain_id());

        uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
      }
      else if( ndTemp->domNums.size() > 1) // the element is cutCell
      {
        totalNGP += ndTemp->Quadrature.gausspoints.size();

        /*
        /////////////////////////////////////////////////////////////
        //  method 1 - adaptive octree
        /////////////////////////////////////////////////////////////
        
        AdaptiveOctree<3>  *adapNd1, *adapNd2, *adapNd3;
        
        adapNd2 = ndTemp->adapIntegNode1;
        
        //adapNd2 = adapNd1->GetChild(0);
        //cout << adapNd2->GetID() << endl;

        while( adapNd2 != NULL )
        {
          if( adapNd2->IsLeaf() )
          {
            //if( (adapNd2->GetDomainNumber() == 0) || (adapNd2->GetDomainNumber() == -1) )
            //{
              bbTemp = adapNd2->GetAABB();

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

              cellDataVTK->InsertNextValue( adapNd2->GetDomainNumber() );
              uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
            //}

            if(adapNd2->GetOrientation() == -1)
              adapNd3 = NULL;
            else
              adapNd3 = adapNd2->GetNeighbour(neigbour_map1[adapNd2->GetOrientation()]);

            //cout <<  adapNd2->GetID() << '\t' << adapNd2->GetOrientation() << endl;
            while( adapNd2->GetOrientation() == SW_FRONT )
            //while(adapNd3 == NULL)
            {
              //cout << " inside IF " << endl;
              adapNd2 = adapNd2->GetParent();
              //cout << " inside IF " << endl;

              if(adapNd2->GetOrientation() == -1)
                adapNd3 = NULL;
              else
	      {
                //cout <<  adapNd2->GetID() << '\t' << adapNd2->GetOrientation() << endl;
                adapNd3 = adapNd2->GetNeighbour(neigbour_map1[adapNd2->GetOrientation()]);
	      }

              //cout << adapNd2->GetID() << endl;
            }

            adapNd2 = adapNd3;
          }
          else
          {
            adapNd2 = adapNd2->GetChild(0);
          }
        }
        */

        /////////////////////////////////////////////////////////////
        //  method 2 - adaptive bindarytree
        /////////////////////////////////////////////////////////////
        
        typedef  typename  node::adapTreePTR  adapTreePTR;

        adapTreePTR  adapNd1, adapNd2, adapNd3;
        
        adapNd2 = ndTemp->adapIntegNode;
        
        //adapNd2 = adapNd1->GetChild(0);
        //cout << adapNd2->GetID() << endl;

        while( adapNd2 != NULL )
        {
          if( adapNd2->IsLeaf() )
          {
            //if( (adapNd2->GetDomainNumber() == 0) || (adapNd2->GetDomainNumber() == -1) )
            //{
              bbTemp = adapNd2->GetAABB();

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

              cellDataVTK->InsertNextValue( adapNd2->GetDomainNumber() );
              cellDataVTK2->InsertNextValue(ndTemp->get_subdomain_id());

              uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
            //}

            if(adapNd2->GetOrientation() == -1)
              adapNd3 = NULL;
            else
              adapNd3 = adapNd2->GetNeighbour(neigbour_map2[adapNd2->GetOrientation()]);

            //cout <<  adapNd2->GetID() << '\t' << adapNd2->GetOrientation() << endl;
            while( adapNd2->GetOrientation() == RIGHT )
            {
              //cout << " inside IF " << endl;
              adapNd2 = adapNd2->GetParent();
              //cout << " inside IF " << endl;

              if(adapNd2->GetOrientation() == -1)
                adapNd3 = NULL;
              else
	      {
                //cout <<  adapNd2->GetID() << '\t' << adapNd2->GetOrientation() << endl;
                adapNd3 = adapNd2->GetNeighbour(neigbour_map2[adapNd2->GetOrientation()]);
	      }

              //cout << adapNd2->GetID() << endl;
            }

            adapNd2 = adapNd3;
          }
          else
          {
            adapNd2 = adapNd2->GetChild(0);
          }
        }

      } //else
    }

    PetscPrintf(MPI_COMM_WORLD, "\n Total number of Gauss points in cut cells = %d \n\n", totalNGP);

    //cellDataVTK->SetName("ElemType");

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetCellData()->SetScalars(cellDataVTK);
    uGridVTK->GetCellData()->AddArray(cellDataVTK2);

    return;
}
//





void  HBSplineCutFEM::postProcessAdapIntegration1D(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    int  ee, ii, jj, kk, nlocal, index;

    nlocal = degree[0]+1;

    VectorXd  N(nlocal), dN_dx(nlocal), d2N_dx2(nlocal), tempVec;
    VectorXd  NN(nlocal), dNN_dx(nlocal), d2NN_dx2(nlocal);
    myPoint  knotIncr, knotBegin, knotEnd;

    double   fact, incr1, *tmp1, xx[2], val[3], xIntf=0.5;

    vector<double>  uu, uu1, uu2;

    node* nd1;

    AdvDiffExact1D  analy;

    vector<vtkIdType>   pt(resln[0]*3);

    index = 0;
    for(ee=0; ee<activeElements.size(); ee++)
    {
        nd1 = elems[activeElements[ee]];

        tmp1 = nd1->GetKnots(0);
        knotBegin = nd1->GetKnotBegin();
        knotIncr  = nd1->GetKnotIncrement();

        //printf("\t tmp[0] and tmp[1]  ... : %12.8f\t%12.8f\n", tmp[0], tmp[1] );

        if( nd1->IsCutElement() )
        {
          fact = (xIntf - tmp1[0])/resln[0];
          create_vector(tmp1[0], xIntf, fact, uu1);

          fact = (tmp1[1] - xIntf)/resln[0];
          create_vector(xIntf, tmp1[1], fact, uu2);
          
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
          fact = (tmp1[1] - tmp1[0])/resln[0];
          create_vector(tmp1[0], tmp1[1], fact, uu);
        }
          //cout << uu << endl;

          incr1 = tmp1[1] - tmp1[0];

          //create the coordinates of the pointsVTK (nodes in FEM)

          for(ii=0;ii<uu.size();ii++)
          {
              param[0] = uu[ii];
              GeomData.computeBasisFunctions1D(knotBegin, knotIncr, param, NN, dNN_dx, d2N_dx2);

              if(nd1->GetParent() == NULL)
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
              
              xx[0] = ComputeGeometry(0, uu[ii]);

              pt[ii] = pointsVTK->InsertNextPoint(xx[0], 0.0, 0.0);
              
              //if(nd1->GetDomainNumber() == 0)
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

              //outp.push_back(nd1->computeValue(0, N));
              //outp2.push_back(nd1->computeForce(0, N));
              //outp2.push_back(nd1->computeValue(0, dN_dx));
              //outp3.push_back(nd1->computeValue(0, d2N_dx2));
              //outp3.push_back(analy.computeValue(0, uu[kk], 0.0));
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
    //cout << " jjjjjjjjjjjjjjjjjj " << endl;

    int  dd, ii, jj, kk, ll, count, nlocal, index, ind1, ind2, e, ee, gcount, ind, domTemp;
    vtkIdType pt[50], cellId;

    nlocal = (degree[0]+1) * (degree[1] + 1);

    VectorXd  NN(nlocal), N(nlocal), dN_dx(nlocal), dN_dy(nlocal), dNN_dx(nlocal), dNN_dy(nlocal), tempVec, tempVec2, d2N_dx2(nlocal), d2N_dy2(nlocal);
    VectorXd  vectmp(nlocal), rhsTemp;
    myPoint  knotIncr, knotBegin;

    double   fact, uleft, uright, *tmp0, *tmp1, incr1, incr2, val1;
    NodeOrientation  neigbour_map2[] = {RIGHT, LEFT};
    vector<double>  uu, vv;

    time_t tstart, tend;
    
    AABB  bbTemp;

    node* ndTemp;

  if(ndf == 1)
  {
    index = 0;
    for(ee=0; ee<activeElements.size(); ee++)
    {
      ndTemp = elems[activeElements[ee]];
      //cout << " Node # " << nd->GetID() << endl;

        tmp0 = ndTemp->GetKnots(Dir1);
        tmp1 = ndTemp->GetKnots(Dir2);

        knotBegin = ndTemp->GetKnotBegin();
        knotIncr  = ndTemp->GetKnotIncrement();

        //printf("\t tmp[0] and tmp[1]  ... : %12.8f\t%12.8f\n", tmp[0], tmp[1] );

        incr1 = tmp0[2] ;
        incr2 = tmp1[2] ;

        if( !ndTemp->IsCutElement() )
        {
          if( ndTemp->GetDomainNumber() == 0 )
          {

          fact = incr1/resln[0];
          create_vector(tmp0[0], tmp0[1], fact, uu);

          fact = incr2/resln[1];
          create_vector(tmp1[0], tmp1[1], fact, vv);

          //create the coordinates of the pointsVTK (nodes in FEM)

          count = 0;
          for(jj=0;jj<vv.size();jj++)
          {
              param[1] = vv[jj];
              geom[1] = ComputeGeometry(1, vv[jj]);

              for(ii=0;ii<uu.size();ii++)
              {
                param[0] = uu[ii];
                geom[0] = ComputeGeometry(0, uu[ii]);

                pt[count] = pointsVTK->InsertNextPoint(geom[0], geom[1], 0.0);

                GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

                if(ndTemp->GetParent() == NULL)
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

          cellDataVTK->InsertNextValue(ndTemp->GetDomainNumber() );
          cellDataVTK2->InsertNextValue(ndTemp->get_subdomain_id());
	  }
        } //if( !nd->IsCutElement() )
        else // the element is cutCell
        {
          vtkSmartPointer<vtkTriangle> triaVTK =  vtkSmartPointer<vtkTriangle>::New();

          myPoly *poly;

          for(ii=0; ii<ndTemp->subTrias.size(); ii++)
          {
            poly = ndTemp->subTrias[ii];

            if( poly->GetDomainNumber() == 0 )
            {
              for(kk=0; kk<3; kk++)
              {
                geom = poly->GetPoint(kk);
                //cout << kk << '\t' << ptTemp[0] << '\t' << ptTemp[1] << endl;

                pt[kk] = pointsVTK->InsertNextPoint(geom[0], geom[1], 0.0);

                geometryToParametric(geom, param);
                GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

                if(ndTemp->GetParent() == NULL)
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

              cellDataVTK->InsertNextValue(poly->GetDomainNumber());
              cellDataVTK2->InsertNextValue(ndTemp->get_subdomain_id());

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
    double vec[3];
    vec[0] = vec[1] = vec[2] = 0.0;

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
        //cout << " Node # " << nd->GetID() << endl;

        knotBegin = ndTemp->GetKnotBegin();
        knotIncr  = ndTemp->GetKnotIncrement();

        //printf("\t tmp[0] and tmp[1]  ... : %12.8f\t%12.8f\n", tmp[0], tmp[1] );

        //if( !nd->IsCutElement() )
        if( ndTemp->GetDomainNumber() == 0 )
        {
            tmp0 = ndTemp->GetKnots(Dir1);
            tmp1 = ndTemp->GetKnots(Dir2);

            incr1 = tmp0[2] ;
            incr2 = tmp1[2] ;

            fact = incr1/resln[0];
            create_vector(tmp0[0], tmp0[1], fact, uu);

            fact = incr2/resln[1];
            create_vector(tmp1[0], tmp1[1], fact, vv);

            //create the coordinates of the pointsVTK (nodes in FEM)
            //cout << " ooooooooooooo " << endl;

            count = 0;
            for(jj=0;jj<vv.size();jj++)
            {
              param[1] = vv[jj];
              geom[1] = ComputeGeometry(1, vv[jj]);
              for(ii=0;ii<uu.size();ii++)
              {
                param[0] = uu[ii];
                geom[0] = ComputeGeometry(0, uu[ii]);

                pt[count++] = pointsVTK->InsertNextPoint(geom[0], geom[1], 0.0);

                GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

                if(ndTemp->GetParent() == NULL)
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

            cellDataVTK->InsertNextValue(ndTemp->GetDomainNumber() );
            cellDataVTK2->InsertNextValue(ndTemp->get_subdomain_id());

            //cout << " ooooooooooooo " << endl;
        } //if( !nd->IsCutElement() )

        if( ndTemp->GetDomainNumber() == -1 )
        {
        //
        /////////////////////////////////////////////////////////////
        //  method 2
        /////////////////////////////////////////////////////////////
        
        typedef  typename  node::adapTreePTR  adapTreePTR;
        
        adapTreePTR  adapNd1, adapNd2, adapNd3;
        
        adapNd2 = ndTemp->adapIntegNode;
        
        //adapNd2 = adapNd1->GetChild(0);
        //cout << adapNd2->GetID() << endl;

        while( adapNd2 != NULL )
        {
          if( adapNd2->IsLeaf() )
          {
            if( (adapNd2->GetDomainNumber() <= 0) )
            {
              bbTemp = adapNd2->GetAABB();

              pt[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], 0.0);
              pt[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], 0.0);
              pt[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], 0.0);
              pt[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], 0.0);

              for(ll=0;ll<4;ll++)
                quadVTK->GetPointIds()->SetId(ll, pt[ll]);

              cellDataVTK->InsertNextValue( adapNd2->GetDomainNumber() );
              //cellDataVTK2->InsertNextValue( adapNd2->GetLevel() );
              cellDataVTK2->InsertNextValue(ndTemp->get_subdomain_id());

              uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());

              tmp0 = adapNd2->GetKnots(Dir1);
              tmp1 = adapNd2->GetKnots(Dir2);

              incr1 = tmp0[2] ;
              incr2 = tmp1[2] ;

              fact = incr1/resln[0];
              create_vector(tmp0[0], tmp0[1], fact, uu);

              fact = incr2/resln[1];
              create_vector(tmp1[0], tmp1[1], fact, vv);

              for(jj=0;jj<vv.size();jj++)
              {
                param[1] = vv[jj];
                for(ii=0;ii<uu.size();ii++)
                {
                  param[0] = uu[ii];

                  GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

                  if(ndTemp->GetParent() == NULL)
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
            } // if( (adapNd2->GetDomainNumber() <= 0) )

            if(adapNd2->GetOrientation() == -1)
              adapNd3 = NULL;
            else
              adapNd3 = adapNd2->GetNeighbour(neigbour_map2[adapNd2->GetOrientation()]);

            //cout <<  adapNd2->GetID() << '\t' << adapNd2->GetOrientation() << endl;
            while( adapNd2->GetOrientation() == RIGHT )
            //while(adapNd3 == NULL)
            {
              //cout << " inside IF " << endl;
              adapNd2 = adapNd2->GetParent();
              //cout << " inside IF " << endl;

              if(adapNd2->GetOrientation() == -1)
                adapNd3 = NULL;
              else
	      {
                //cout <<  adapNd2->GetID() << '\t' << adapNd2->GetOrientation() << endl;
                adapNd3 = adapNd2->GetNeighbour(neigbour_map2[adapNd2->GetOrientation()]);
	      }

              //cout << adapNd2->GetID() << endl;
            } //while( adapNd2->GetOrientation() == RIGHT )

            adapNd2 = adapNd3;
          } // if( adapNd2->IsLeaf() )
          else
          {
            adapNd2 = adapNd2->GetChild(0);
          }
        } // while( adapNd2 != NULL )
        //          //cout << " AAAAAAAAAA " << endl;
        } // else
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

    int  dd, ii, jj, kk, ll, count, nlocal, index, ind1, ind2, e, ee, gcount, ind, domTemp;
    vtkIdType ptIds[50], cellId;

    nlocal = (degree[0]+1) * (degree[1] + 1) * (degree[2] + 1);

    VectorXd  NN(nlocal), dNN_dx(nlocal), dNN_dy(nlocal), dNN_dz(nlocal),  tempVec, tempVec2;
    VectorXd  N(nlocal), dN_dx(nlocal), dN_dy(nlocal), dN_dz(nlocal);
    VectorXd  vectmp(nlocal), rhsTemp;
    myPoint   knotIncr, knotBegin;

    double   fact, *tmp0, *tmp1, *tmp2, incr1, incr2, incr3, val1;

    NodeOrientation  neigbour_map1[] = {EAST, NORTH, FRONT, WEST, WEST, WEST, EAST, SOUTH};
    NodeOrientation  neigbour_map2[] = {RIGHT, LEFT};

    vector<int> bfs;
    vector<double>  uu, vv, ww;

    time_t tstart, tend;

    node* ndTemp;

  //if(ndf == 1)
  //{
  //}
  //else // for Stokes and Navier-Stokes
  //{
    double vec[3];
    vec[0] = vec[1] = vec[2] = 0.0;

    vecVTK->SetNumberOfComponents(3);
    //vecVTK->SetNumberOfTuples(count);
    vecVTK2->SetNumberOfComponents(3);
    //vecVTK2->SetNumberOfTuples(count);
    //scaVTK->SetNumberOfTuples(count);
    //scaVTK2->SetNumberOfTuples(count);

    cout << " Node aaaaaaaaaaa " << endl;
    
    for(ee=0; ee<activeElements.size(); ee++)
    {
        ndTemp = elems[activeElements[ee]];
        //cout << " Node # " << nd->GetID() << endl;

        knotBegin = ndTemp->GetKnotBegin();
        knotIncr  = ndTemp->GetKnotIncrement();

          if( (ndTemp->domNums.size() == 1) && (ndTemp->domNums[0] == 0) )
          {
             tmp0 = ndTemp->GetKnots(Dir1);
             tmp1 = ndTemp->GetKnots(Dir2);
             tmp2 = ndTemp->GetKnots(Dir3);

             incr1 = tmp0[2] ;
             incr2 = tmp1[2] ;
             incr3 = tmp2[2] ;

             fact = incr1/resln[0];
             create_vector(tmp0[0], tmp0[1], fact, uu);

             fact = incr2/resln[1];
             create_vector(tmp1[0], tmp1[1], fact, vv);

             fact = incr3/resln[2];
             create_vector(tmp2[0], tmp2[1], fact, ww);

             //create the coordinates of the pointsVTK (nodes in FEM)
             //cout << " ooooooooooooo " << endl;

            count = 0;
            for(kk=0; kk<ww.size(); kk++)
            {
              param[2] = ww[kk];
              geom[2] = ComputeGeometry(2, ww[kk]);

            for(jj=0;jj<vv.size();jj++)
            {
              param[1] = vv[jj];
              geom[1] = ComputeGeometry(1, vv[jj]);

              for(ii=0;ii<uu.size();ii++)
              {
                param[0] = uu[ii];
                geom[0] = ComputeGeometry(0, uu[ii]);

                ptIds[count++] = pointsVTK->InsertNextPoint(geom[0], geom[1], geom[2]);

                GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

                if(ndTemp->GetParent() == NULL)
                {
                  N = NN;
                  dN_dx = dNN_dx;
                  dN_dy = dNN_dy;
                  dN_dz = dNN_dz;
                }
                else
                {
                  N = ndTemp->SubDivMat*NN;
                  dN_dx = ndTemp->SubDivMat*dNN_dx;
                  dN_dy = ndTemp->SubDivMat*dNN_dy;
                  dN_dz = ndTemp->SubDivMat*dNN_dz;
                }

                vec[0] = ndTemp->computeValue(0, N);
                vec[1] = ndTemp->computeValue(1, N);
                vec[2] = ndTemp->computeValue(2, N);

                fact   = ndTemp->computeValue(3, N);

                vecVTK->InsertNextTuple(vec);
                vecVTK2->InsertNextTuple(vec);
                scaVTK->InsertNextValue(fact);
                //scaVTK2->InsertNextValue(fact);

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
            //cellDataVTK->InsertNextValue(0);
            //cellDataVTK2->InsertNextValue(0);

            cellDataVTK->InsertNextValue(ndTemp->GetDomainNumber() );
            cellDataVTK2->InsertNextValue(ndTemp->get_subdomain_id());

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
            if( adapNd2->IsLeaf() )
            {
              if( adapNd2->GetDomainNumber() <= 0 )
              {
                //bbTemp = adapNd2->GetAABB();

                tmp0 = adapNd2->GetKnots(Dir1);
                tmp1 = adapNd2->GetKnots(Dir2);
                tmp2 = adapNd2->GetKnots(Dir3);

                incr1 = tmp0[2] ;
                incr2 = tmp1[2] ;
                incr3 = tmp2[2] ;

                fact = incr1/resln[0];
                create_vector(tmp0[0], tmp0[1], fact, uu);

                fact = incr2/resln[1];
                create_vector(tmp1[0], tmp1[1], fact, vv);

                fact = incr3/resln[2];
                create_vector(tmp2[0], tmp2[1], fact, ww);

                //cout << " ooooooooooooo " << endl;

                count = 0;
                for(kk=0; kk<ww.size(); kk++)
                {
                  param[2] = ww[kk];
                  geom[2] = ComputeGeometry(2, ww[kk]);

                for(jj=0;jj<vv.size();jj++)
                {
                  param[1] = vv[jj];
                  geom[1] = ComputeGeometry(1, vv[jj]);

                for(ii=0;ii<uu.size();ii++)
                {
                  param[0] = uu[ii];
                  geom[0] = ComputeGeometry(0, uu[ii]);

                  ptIds[count++] = pointsVTK->InsertNextPoint(geom[0], geom[1], geom[2]);

                  GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

                  if(ndTemp->GetParent() == NULL)
                  {
                    N = NN;
                    dN_dx = dNN_dx;
                    dN_dy = dNN_dy;
                    dN_dz = dNN_dz;
                  }
                  else
                  {
                    N = ndTemp->SubDivMat*NN;
                    dN_dx = ndTemp->SubDivMat*dNN_dx;
                    dN_dy = ndTemp->SubDivMat*dNN_dy;
                    dN_dz = ndTemp->SubDivMat*dNN_dz;
                  }

                  vec[0] = ndTemp->computeValue(0, N);
                  vec[1] = ndTemp->computeValue(1, N);
                  vec[2] = ndTemp->computeValue(2, N);
                  fact   = ndTemp->computeValue(3, N);

                  vecVTK->InsertNextTuple(vec);
                  vecVTK2->InsertNextTuple(vec);
                  scaVTK->InsertNextValue(fact);
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

                cellDataVTK->InsertNextValue(ndTemp->GetDomainNumber() );
                cellDataVTK2->InsertNextValue(ndTemp->get_subdomain_id());

            } // if( (adapNd2->GetDomainNumber() <= 0) )

            if(adapNd2->GetOrientation() == -1)
              adapNd3 = NULL;
            else
              adapNd3 = adapNd2->GetNeighbour(neigbour_map2[adapNd2->GetOrientation()]);

            //cout <<  adapNd2->GetID() << '\t' << adapNd2->GetOrientation() << endl;
            while( adapNd2->GetOrientation() == RIGHT )
            //while(adapNd3 == NULL)
            {
              //cout << " inside IF " << endl;
              adapNd2 = adapNd2->GetParent();
              //cout << " inside IF " << endl;

              if(adapNd2->GetOrientation() == -1)
                adapNd3 = NULL;
              else
	      {
                //cout <<  adapNd2->GetID() << '\t' << adapNd2->GetOrientation() << endl;
                adapNd3 = adapNd2->GetNeighbour(neigbour_map2[adapNd2->GetOrientation()]);
	      }
              //cout << adapNd2->GetID() << endl;
            }
            adapNd2 = adapNd3;
          }
          else
          {
            adapNd2 = adapNd2->GetChild(0);
          }
          //cout << " AAAAAAAAAA " << endl;
	}// else
      } // for(

    //cout << " jjjjjjjjjjjjjjjjjj " << endl;
    vecVTK->SetName("vel");
    vecVTK2->SetName("vort");
    scaVTK->SetName("pres");

    //assign nodal coordinates and field data to uGridVTK
    // no need to create lookup table here. All this stuff can be done in Paraview

    uGridVTK->SetPoints(pointsVTK);

    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->SetVectors(vecVTK);
    uGridVTK->GetPointData()->AddArray(vecVTK2);
    //uGridVTK->GetPointData()->AddArray(scaVTK2);

    uGridVTK->GetCellData()->SetScalars(cellDataVTK);
    uGridVTK->GetCellData()->AddArray(cellDataVTK2);

    //cout << " jjjjjjjjjjjjjjjjjj " << endl;
  }

    return;
}



