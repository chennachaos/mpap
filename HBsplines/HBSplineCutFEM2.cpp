
#include "HBSplineCutFEM.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "Files.h"
#include "MyString.h"
#include "headersVTK.h"
#include "myGeomUtilities.h"
#include "myTria.h"
#include "myTet.h"


//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Triangulation_3.h>
//#include <CGAL/Delaunay_triangulation_3.h>
//#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <omp.h>

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern Files files;

using namespace myGeom;


void  HBSplineCutFEM::prepareCutElements()
{
  int  ee, bb, domTemp;
  bool  flag;

  node  *nd1;

  if( ImmersedBodyObjects.size() == 0 )
  {
    for(ee=0; ee<activeElements.size(); ee++)
    {
      nd1 = elems[activeElements[ee]];
      nd1->domNums.clear();
      nd1->domNums.push_back(0);
    }
    return;
  }
  
  cutCellIds.clear();

  for(int bb=0; bb<ImmersedBodyObjects.size(); bb++)
  {
    ImmersedBodyObjects[bb]->UpdateImmersedFaces();
    ImmersedBodyObjects[bb]->computeAABB(2);
  }

/*
  if(CUTCELL_INTEGRATION_TYPE == 1)
  {
    if(ndm == 1)
      prepareCutElementsSubTrias1D();
    else if(ndm == 2)
      prepareCutElementsSubTrias2D();
    else
      prepareCutElementsSubTrias3D();
  }
  else
  {
    if(ndm == 1)
      prepareCutElementsAdapIntegration1D();
    else if(ndm == 2)
      prepareCutElementsAdapIntegration2D();
    else
      prepareCutElementsAdapIntegration3D();
  }
*/

  cout << " HBSplineCutFEM::prepareCutElements() " << endl;

/*
//#pragma omp parallel for default(shared)  private(ee, nd1)
//{
  for(ee=0; ee<activeElements.size(); ee++)
  {
    nd1 = elems[activeElements[ee]];

    //cout << " gggggggggggg " << ee << endl;

    //domTemp = *std::max_element(nd1->domNums.begin(), nd1->domNums.end() ) - 1;

    //if( !GeomData.domainFixedYesNo[domTemp] )
    //{
      if( nd1->prepareCutCell(cutFEMparams) != 1)
      {
	cerr << " Error in prepareCutCell() for TreeNode<2>::prepareCutCell() " << endl;
      }
      
      //if( nd1->IsRightBoundary() && nd1->IsTopBoundary() && nd1->IsBackBoundary() )
        //printVector(nd1->domNums);

      if( nd1->IsCutElement() )
      {
        cutCellIds.push_back(activeElements[ee]);

        //cout << " hhhhhhhhhhhh " << endl;

        if(cutFEMparams[0] == 1)
          nd1->computeGaussPointsSubTrias(cutFEMparams[2], false);
        else
          nd1->computeGaussPointsAdapIntegration(cutFEMparams[3], cutFEMparams[4], false, true);
      }
      //cout << " ffffffffffff " << endl;
    //}
  }
//}
*/

  //GeomData.domainFixedYesNo[0] = 0;
  //GeomData.domainFixedYesNo[1] = 0;

//#pragma omp parallel for default(shared)  private(ee, nd1)
//{
  for(ee=0; ee<activeElements.size(); ee++)
  {
    nd1 = elems[activeElements[ee]];

    //cout << " gggggggggggg " << ee << endl;
    
    flag = false;

    if( nd1->domNums.size() == 0 )
    {
      flag = true;
    }
    else if( nd1->domNums.size() == 1 )
    {
      domTemp = nd1->domNums[0] ;

      if( domTemp == 0 )
      {
        flag = true;
      }
      else
      {
        if( !GeomData.domainFixedYesNo[domTemp-1] )
          flag = true;
      }
    }
    else
    {
      domTemp = *std::max_element(nd1->domNums.begin(), nd1->domNums.end() ) - 1;

      if( !GeomData.domainFixedYesNo[domTemp] )
        flag = true;
    } 


    if( flag )
    {
      if( nd1->prepareCutCell(cutFEMparams) != 1)
      {
        cerr << " Error in prepareCutCell() for TreeNode<2>::prepareCutCell() " << endl;
      }
      
      //if( nd1->IsRightBoundary() && nd1->IsTopBoundary() && nd1->IsBackBoundary() )
        //printVector(nd1->domNums);

      if( nd1->IsCutElement() )
      {
        cutCellIds.push_back(activeElements[ee]);

        //cout << " hhhhhhhhhhhh " << endl;

        if(cutFEMparams[0] == 1)
          nd1->computeGaussPointsSubTrias(cutFEMparams[2], false);
        else
          nd1->computeGaussPointsAdapIntegration(cutFEMparams[3], cutFEMparams[4], false, true);
      }
      //cout << " ffffffffffff " << endl;
    }
  }
//}

  //int resln[]={1,1,1};

  //plotGeom(1, 1, 10, 0, resln);

  cout << " HBSplineCutFEM::prepareCutElements() " << endl;

  return;
}

  
void  HBSplineCutFEM::prepareCutElementsAdapIntegration1D()
{
  return;
}



void  HBSplineCutFEM::prepareCutElementsAdapIntegration2D()
{
  return;
}


void  HBSplineCutFEM::prepareCutElementsAdapIntegration3D()
{
  return;
}





void  HBSplineCutFEM::prepareCutElementsSubTrias1D()
{
  return;
}




void  HBSplineCutFEM::prepareCutElementsSubTrias2D()
{
  int  ee, bb, ii, jj, domTemp;

  node  *nd1;

  AABB  bbTemp;
  vector<myPoint>  ptOut;
  vector<int>  cornerInOut(4);


  for(ee=0; ee<activeElements.size(); ee++)
  {
    nd1 = elems[activeElements[ee]];

    bbTemp = nd1->GetAABB();
    
    ptOut.clear();

    domTemp = 0;
    for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
    {
      if( ImmersedBodyObjects[bb]->doAABBintersect(bbTemp) )
      {
        domTemp = ImmersedBodyObjects[bb]->doIntersect2D(bbTemp, true, cornerInOut, ptOut) ;
        break;
      }
    }

    //if( !std::any_of(domTemp.begin(), domTemp.end(), [](int i){return i==-1;}) )
    //{
      //nd1->SetDomainNumber( *std::max_element(domTemp.begin(), domTemp.end()) );
    
    //nd1->SetDomainNumber(domTemp);

    if(domTemp == -1)
    {
      cutCellIds.push_back(activeElements[ee]);

      // remove the old subtriangulation

      nd1->clearSubtriangulation();

      // find all the intersection points

      int dd;
      vtkIdType pt[20], ii, id, kk;
      double  ptTemp[3];
      myPoint  ptNew;
      std::vector<int>::iterator itInt;
      vector<int>  vecTemp;

      vtkSmartPointer<vtkPoints>       pointsLoc    =  vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkPolyData>     polyDataLoc  =  vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkPolyData>     polyDataLoc2 =  vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkDelaunay2D>   delaunay     =  vtkSmartPointer<vtkDelaunay2D>::New();

      //cout << " AAAAAAAAA " << endl;

      vector<myPoint>  ptVec;
      
      //cout << " Number of intersection points = " << ptOut.size() << endl;

      for(ii=0; ii<ptOut.size(); ii++)
      {
        ptNew = ptOut[ii];
        //cout << ii << '\t' << ptNew[0] << '\t' << ptNew[1] << '\t' << ptNew[2] << endl;
        if( !( pointExists(ptVec, ptNew) ) )
        {
          ptVec.push_back(ptNew);
          pointsLoc->InsertNextPoint( ptNew[0], ptNew[1], ptNew[2] );
        }
      }

      //cout << " domTemp = " << activeElements[ee] << '\t' << domTemp << '\t' << ptOut.size() << '\t' << ptVec.size() <<  endl;

      //cout << " AAAAAAAAA " << endl;
      
      // generate subtriangulation using Delauny procedure in VTK library

      polyDataLoc->SetPoints(pointsLoc);

      //delaunay->BoundingTriangulationOn();
      //delaunay->SetTolerance(0.0010);
      //delaunay->SetAlpha(0.0);
      delaunay->SetOffset(1000.0);

      //delaunay->SetInputData(polyDataLoc);
      delaunay->SetInput(polyDataLoc);
      delaunay->Update();

      polyDataLoc2 = delaunay->GetOutput();

      vtkIdType  cellId2;
      vtkCell  *cellVTK2;
      myPoint  pt1, pt2, pt3, pt4;

      for(cellId2=0; cellId2<polyDataLoc2->GetNumberOfCells(); cellId2++)
      {
        cellVTK2 = polyDataLoc2->GetCell(cellId2);

        polyDataLoc2->GetPoint(cellVTK2->GetPointId(0), &pt1[0]);
        polyDataLoc2->GetPoint(cellVTK2->GetPointId(1), &pt2[0]);
        polyDataLoc2->GetPoint(cellVTK2->GetPointId(2), &pt3[0]);

        myPoly  *poly = new myTria(pt1, pt2, pt3);

        poly->centroid(pt4);
        //pt4[1] += 1.0e-6;

        dd = 0;
        for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
        {
          if( ImmersedBodyObjects[bb]->within(pt4) )
            dd = bb+1;
        }
        //cout << cellId2 << '\t' << dd << endl;

        vecTemp.push_back(dd);
        
        poly->SetDomainNumber(dd);

        nd1->subTrias.push_back(poly);
      } //  for(cellId2=0; cellId2<polyDataLoc2->GetNumberOfCells(); cellId2++)
      
      if( std::equal(vecTemp.begin()+1, vecTemp.end(), vecTemp.begin()) )
      {
        //nd1->SetDomainNumber(vecTemp[0]);
        nd1->clearSubtriangulation();
      }
      
      // once the subtriangulation is generated 
      // find the Gausspoints for the triangles which lie 
      // in the domain #0 (fluid domain)

      nd1->computeGaussPointsSubTrias(cutFEMparams[2]);
    } //else
  }

  return;
}





/*
void  HBSplineCutFEM::prepareCutElements3D()
{
  int  ee, bb, ii, jj, pp;

  node  *nd1, *nd2;

  std::vector<int>  domTemp(ImmersedBodyObjects.size());

  AABB  bbTemp;

  vtkIdType  cellId;
  vtkCell   *cellVTK;

  std::vector<int>  vecTemp(8);
  std::vector<int>::iterator  minVal, maxVal;


  for(int bb=0; bb<ImmersedBodyObjects.size(); bb++)
  {
    ImmersedBodyObjects[bb]->selectEnclosedPoints->SetInputData(pointsPolydata);
    ImmersedBodyObjects[bb]->selectEnclosedPoints->Update();
  }

  //for(ee=0; ee<activeElements.size(); ee++)
  //{
      //nd1 = elems[activeElements[ee]];

  for(cellId=0; cellId<uGridVTKfluid->GetNumberOfCells(); cellId++)
  {
    //cout << " cellId = " << cellId << endl;

    nd1 = elems[activeElements[cellId]];

    bbTemp = nd1->GetAABB();

    cellVTK = uGridVTKfluid->GetCell(cellId);

    for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
    {
      if( ImmersedBodyObjects[bb]->doAABBintersect(bbTemp) )
      {
        //cout << " aaaaaaa " << cellId << '\t' << bb << endl;

        for(pp=0; pp<cellVTK->GetNumberOfPoints(); pp++)
        {
          //cout << " pp = " << pp << endl;
          vecTemp[pp] = ImmersedBodyObjects[bb]->selectEnclosedPoints->IsInside(cellVTK->GetPointId(pp));
        }

        if( std::equal(vecTemp.begin()+1, vecTemp.end(), vecTemp.begin()) )
        {
          if( vecTemp[0] == 1 )
            domTemp[bb] = bb+1;
          else
            domTemp[bb] = 0;
        }
        else
          domTemp[bb] = -1;
      }
      else
        domTemp[bb] = 0;
    }

    if( !std::any_of(domTemp.begin(), domTemp.end(), [](int i){return i==-1;}) )
    {
      nd1->SetDomainNumber( *std::max_element(domTemp.begin(), domTemp.end()) );
    } 
    else
    {
      nd1->SetDomainNumber(-1);

      //if(2 < 1)
      //{
      //cout << " subtriangulation " << endl;
      // remove the old subtriangulation
      nd1->clearSubtriangulation();

      // find all the intersection points

      int dd, ii, jj, inserted, kk;
      vtkIdType pt[100], id;
      double  ptTemp[3];
      myPoint  ptNew;
      vector<myPoint>  ptList;
      vector<int>  vecTemp;

      vtkSmartPointer<vtkPoints>       pointsLoc    =  vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkPoints>       pointsLoc2   =  vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkPolyData>     polyDataCell =  vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkPolyData>     polyDataLoc2 =  vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkPolyData>     polyDataLoc3 =  vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkDelaunay3D>   delaunay     =  vtkSmartPointer<vtkDelaunay3D>::New();

      vtkSmartPointer<vtkCellArray>    polyList     =  vtkSmartPointer<vtkCellArray>::New();
      vtkSmartPointer<vtkPolygon>      polygonVTK   =  vtkSmartPointer<vtkPolygon>::New();

      vtkSmartPointer<vtkUnstructuredGrid> uGridLoc =   vtkSmartPointer<vtkUnstructuredGrid>::New();


      pointsLoc->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.minBB[2]);
      pointsLoc->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.minBB[2]);
      pointsLoc->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]);
      pointsLoc->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]);

      pointsLoc->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]);
      pointsLoc->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]);
      pointsLoc->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]);
      pointsLoc->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]);

      for(ii=0; ii<pointsLoc->GetNumberOfPoints(); ii++)
      {
        pointsLoc->GetPoint(ii, &(ptNew[0]) );
        //cout << ptNew[0] << '\t' << ptNew[1] << endl;

        ptList.push_back(ptNew);
      }

      polygonVTK->GetPointIds()->SetNumberOfIds(3);

      //int verts[12][3] ={ {0,1,2}, {2,1,3},{2,3,7}, {7,3,1}, {7,1,5}, {5,1,4}, {5,4,7}, {7,4,6}, {7,6,2}, {2,6,4}, {2,4,0}, {0,4,1} };
      int verts[12][3] ={ {0,2,1}, {2,3,1},{2,7,3}, {7,1,3}, {7,5,1}, {5,4,1}, {5,7,4}, {7,6,4}, {7,2,6}, {2,4,6}, {2,0,4}, {0,1,4} };

      for(ii=0; ii<12; ii++)
      {
        for(jj=0; jj<3; jj++)
          polygonVTK->GetPointIds()->SetId(jj, verts[ii][jj]);

        polyList->InsertNextCell(polygonVTK);
      }

      polyDataCell->SetPoints(pointsLoc);
      polyDataCell->SetPolys(polyList);

      vtkSmartPointer<vtkSelectEnclosedPoints>  selectEnclosedPointsCell  =   vtkSmartPointer<vtkSelectEnclosedPoints>::New();

      selectEnclosedPointsCell->CheckSurfaceOn();

      selectEnclosedPointsCell->SetTolerance(0.000001);
      selectEnclosedPointsCell->Initialize(polyDataCell);

      //cout << " BBBBBBBBB  " << endl;

      if(1 > 2)
      {
      vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOperation  =  vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();

      booleanOperation->SetTolerance(1.0e-10);

      //booleanOperation->SetOperationToUnion();
      //booleanOperation->SetOperationToIntersection();
      booleanOperation->SetOperationToDifference();

      //booleanOperation->ReorientDifferenceCellsOn();

      booleanOperation->SetInputData( 0, polyDataCell );
      booleanOperation->SetInputData( 1, ImmersedBodyObjects[0]->polyData );
      //cout << " CCCCCCCCC " << endl;
      booleanOperation->Update();

      polyDataLoc2 = booleanOperation->GetOutput();
      }//

      //
      vtkSmartPointer<vtkIntersectionPolyDataFilter>  intersectionOperation  =  vtkSmartPointer<vtkIntersectionPolyDataFilter>::New();

      intersectionOperation->SplitFirstOutputOn();
      intersectionOperation->SplitSecondOutputOn();

      //intersectionOperation->SplitFirstOutputOff();
      //intersectionOperation->SplitSecondOutputOff();

      intersectionOperation->SetInputData( 0, polyDataCell );
      intersectionOperation->SetInputData( 1, ImmersedBodyObjects[0]->polyData );
      intersectionOperation->Update();

      polyDataLoc2 = intersectionOperation->GetOutput(1);
      //

  vtkSmartPointer<vtkXMLPolyDataWriter>  writerPolyData2   =     vtkSmartPointer<vtkXMLPolyDataWriter>::New();

  char fname1[50];
  sprintf(fname1,"%s%d%s", "boolean3D-", cellId,".vtp");

  //Write the file.
  writerPolyData2->SetFileName(fname1);
  writerPolyData2->SetInputData(polyDataLoc2);
  writerPolyData2->Write();

      //cout << " AAAAAAAAA " << endl;

      //
      for(ii=0; ii<polyDataLoc2->GetNumberOfPoints(); ii++)
      {
        polyDataLoc2->GetPoint(ii, &(ptNew[0]) );
        //cout << ptNew[0] << '\t' << ptNew[1] << '\t' << ptNew[2] << endl;
        inserted = pointExists(ptList, ptNew);

        if(!inserted)
        {
          ptList.push_back(ptNew);
          pointsLoc2->InsertNextPoint( ptNew[0], ptNew[1], ptNew[2] );
        }
      }
      //cout << " AAAAAAAAA " << endl;
      //
      //
      for(ii=0; ii<ImmersedBodyObjects[0]->polyData->GetNumberOfPoints(); ii++)
      {
        ImmersedBodyObjects[0]->polyData->GetPoint(ii, &(ptNew[0]) );
        
        if( selectEnclosedPointsCell->IsInsideSurface( ptNew[0], ptNew[1], ptNew[2] ) )
        {
          //cout << ptNew[0] << '\t' << ptNew[1] << '\t' << ptNew[2] << endl;
          pointsLoc2->InsertNextPoint( ptNew[0], ptNew[1], ptNew[2] );
        }
      }
      //
      //
      if(2 > 3)
      {
      cout << " lllllllllllll " << endl;

      for(ii=0; ii<vecTemp.size(); ii++)
      {
        if(vecTemp[ii] == 0)
        {
          pointsVTKfluidgrid->GetPoint(cellVTK->GetPointId(ii), &(ptNew[0]));
          cout << ptNew[0] << '\t' << ptNew[1] << '\t' << ptNew[2] << endl;
          pointsLoc2->InsertNextPoint( ptNew[0], ptNew[1], ptNew[2] );
        }
      }
      }//


      ////////////////////////////////////////
      // 
      // generate subtriangulation using Delauny procedure in VTK library
      // 
      ////////////////////////////////////////

      polyDataLoc3->SetPoints(pointsLoc2);

      //delaunay->BoundingTriangulationOn();
      //delaunay->SetTolerance(0.0010);
      //delaunay->SetAlpha(0.0);
      delaunay->SetOffset(1000.0);

      delaunay->SetInputData(polyDataLoc3);
      //delaunay->SetInputData(polyDataLoc2);
      delaunay->Update();

      uGridLoc = delaunay->GetOutput();

      vtkIdType  cellId2;
      vtkCell  *cellVTK2;
      myPoint  pt1, pt2, pt3, pt4;

      //cout << " Number of tetras = " << uGridLoc->GetNumberOfCells() << endl;
      
      for(cellId2=0; cellId2<uGridLoc->GetNumberOfCells(); cellId2++)
      {
        cellVTK2 = uGridLoc->GetCell(cellId2);

        uGridLoc->GetPoint(cellVTK2->GetPointId(0), &pt1[0]);
        uGridLoc->GetPoint(cellVTK2->GetPointId(1), &pt2[0]);
        uGridLoc->GetPoint(cellVTK2->GetPointId(2), &pt3[0]);
        uGridLoc->GetPoint(cellVTK2->GetPointId(3), &pt4[0]);

        myPoly  *poly = new myTet(pt1, pt2, pt3, pt4);

        poly->centroid(pt1);

        //for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
          //domTemp[bb] = ImmersedBodyObjects[bb]->within(pt1) ;

        //dd = myGeom::myFindInt(domTemp, 1) + 1;

        dd =  ImmersedBodyObjects[0]->selectEnclosedPoints2->IsInsideSurface( pt1[0], pt1[1], pt1[2] ) ;

        //cout << pt1.norm() << '\t' << dd << endl;

        vecTemp.push_back(dd);
        
        poly->SetDomainNumber(dd);

        nd1->subTrias.push_back(poly);
      } //  for(cellId2=0; cellId2<polyDataLoc2->GetNumberOfCells(); cellId2++)
      
      //if( std::equal(vecTemp.begin()+1, vecTemp.end(), vecTemp.begin()) )
      //{
        //nd1->SetDomainNumber(vecTemp[0]);
        //nd1->clearSubtriangulation();
      //}
      //

      //}
    } //else
  }
  return;
}
*/


/*
void  HBSplineCutFEM::prepareCutElements3D()
{
  int  ee, bb, ii, jj, pp;

  node  *nd1, *nd2;

  std::vector<int>  domTemp(ImmersedBodyObjects.size());

  AABB  bbTemp;

  vtkIdType  cellId;
  vtkCell   *cellVTK;

  std::vector<int>  vecTemp(8);
  std::vector<int>::iterator  minVal, maxVal;


  for(int bb=0; bb<ImmersedBodyObjects.size(); bb++)
  {
    ImmersedBodyObjects[bb]->selectEnclosedPoints->SetInputData(pointsPolydata);
    ImmersedBodyObjects[bb]->selectEnclosedPoints->Update();
  }

  //for(ee=0; ee<activeElements.size(); ee++)
  //{
      //nd1 = elems[activeElements[ee]];

  for(cellId=0; cellId<uGridVTKfluid->GetNumberOfCells(); cellId++)
  {
    //cout << " cellId = " << cellId << endl;

    nd1 = elems[activeElements[cellId]];

    bbTemp = nd1->GetAABB();

    cellVTK = uGridVTKfluid->GetCell(cellId);

    for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
    {
      if( ImmersedBodyObjects[bb]->doAABBintersect(bbTemp) )
      {
        //cout << " aaaaaaa " << cellId << '\t' << bb << endl;

        for(pp=0; pp<cellVTK->GetNumberOfPoints(); pp++)
        {
          //cout << " pp = " << pp << endl;
          vecTemp[pp] = ImmersedBodyObjects[bb]->selectEnclosedPoints->IsInside(cellVTK->GetPointId(pp));
        }

        if( std::equal(vecTemp.begin()+1, vecTemp.end(), vecTemp.begin()) )
        {
          if( vecTemp[0] == 1 )
            domTemp[bb] = bb+1;
          else
            domTemp[bb] = 0;
        }
        else
          domTemp[bb] = -1;
      }
      else
        domTemp[bb] = 0;
    }

    if( !std::any_of(domTemp.begin(), domTemp.end(), [](int i){return i==-1;}) )
    {
      nd1->SetDomainNumber( *std::max_element(domTemp.begin(), domTemp.end()) );
    } 
    else
    {
      nd1->SetDomainNumber(-1);
    }
  }
  return;
}
*/


/*
void  HBSplineCutFEM::prepareCutElements3D()
{
  int  ee, bb, ii, jj;

  node  *nd1;

  int  domTemp;

  AABB  bbTemp;

  for(int bb=0; bb<ImmersedBodyObjects.size(); bb++)
  {
    //ImmersedBodyObjects[bb]->computeFaceNormals();
    for(ii=0; ii<ImmersedBodyObjects[bb]->ImmersedFaces.size(); ii++)
      ImmersedBodyObjects[bb]->ImmersedFaces[ii]->computeNormal();
  }

  for(ee=0; ee<activeElements.size(); ee++)
  {
    nd1 = elems[activeElements[ee]];

    bbTemp = nd1->GetAABB();
    
    //cout << " ee = " << ee << endl;

    //for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
      //domTemp[bb] = ImmersedBodyObjects[bb]->doIntersect3D(bbTemp) ;

    domTemp = ImmersedBodyObjects[0]->doIntersect3D(bbTemp) ;
    //cout << " domTemp = " << domTemp[0] << endl;

     nd1->SetDomainNumber(domTemp);
  }
  return;
}
*/



void  HBSplineCutFEM::prepareCutElementsSubTrias3D()
{
  int  ee, bb, ii, jj, dd, domTemp, zz=0;

  node  *ndTemp;

  AABB  bbTemp;

  vector<myPoint>  ptOut;
  vector<int>  vecTemp(8);
  myPoint  ptNew;


  for(ee=0; ee<activeElements.size(); ee++)
  {
    ndTemp = elems[activeElements[ee]];

    bbTemp = ndTemp->GetAABB();

    ptOut.clear();

    domTemp = 0;
    for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
    {
      if( ImmersedBodyObjects[bb]->doAABBintersect(bbTemp) )
        domTemp = ImmersedBodyObjects[bb]->doIntersect3D(bbTemp, true, vecTemp, ptOut) ;
    }

    //ndTemp->SetDomainNumber(domTemp);

    //if( (domTemp == -1) && flag && (ee==56) )
    //if( (domTemp == -1) && flag )

    if(domTemp == -1)
    {
      cutCellIds.push_back(activeElements[ee]);

      vtkSmartPointer<vtkPoints>       pointsLoc    =  vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkPolyData>     polyDataLoc  =  vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkPolyData>     polyDataLoc2 =  vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkDelaunay2D>   delaunay2d   =  vtkSmartPointer<vtkDelaunay2D>::New();
      vtkSmartPointer<vtkDelaunay3D>   delaunay3d   =  vtkSmartPointer<vtkDelaunay3D>::New();

      vtkSmartPointer<vtkUnstructuredGrid> uGridLoc =  vtkSmartPointer<vtkUnstructuredGrid>::New();

      vector<myPoint>  ptVec;
      
      //cout << " Number of intersection points = " << ptOut.size() << endl;

      for(ii=0; ii<ptOut.size(); ii++)
      {
        ptNew = ptOut[ii];
        //cout << ii << '\t' << ptNew[0] << '\t' << ptNew[1] << '\t' << ptNew[2] << endl;
        if( !( pointExists(ptVec, ptNew) ) )
        {
          ptVec.push_back(ptNew);
          pointsLoc->InsertNextPoint( ptNew[0], ptNew[1], ptNew[2] );
        }
      }
      
      //cout << " Number of intersection points = " << ptVec.size() << endl;
      //cout << " performing Delaunay in 3D " << endl;

      ////////////////////////////////////////
      // 
      // generate subtriangulation using Delauny procedure in VTK library
      // 
      ////////////////////////////////////////

      polyDataLoc->SetPoints(pointsLoc);

      vtkIdType  cellId2;
      vtkCell  *cellVTK2;
      myPoint  pt1, pt2, pt3, pt4;

      vector<int>  vecTemp;

      //delaunay3d->BoundingTriangulationOn();
      //delaunay3d->SetAlpha(0.0);

      /*
      //delaunay2d->SetTolerance(0.0010);
      delaunay2d->SetOffset(10.0);
      delaunay2d->SetInputData(polyDataLoc);
      delaunay2d->Update();

      polyDataLoc2 = delaunay2d->GetOutput();

      //cout << " kkkkkkkkkkkkkkkk " << endl;

      for(cellId2=0; cellId2<polyDataLoc2->GetNumberOfCells(); cellId2++)
      {
        cellVTK2 = polyDataLoc2->GetCell(cellId2);

        polyDataLoc2->GetPoint(cellVTK2->GetPointId(0), &pt1[0]);
        polyDataLoc2->GetPoint(cellVTK2->GetPointId(1), &pt2[0]);
        polyDataLoc2->GetPoint(cellVTK2->GetPointId(2), &pt3[0]);
        //cout << " lllllllllllll " << endl;
        myPoly  *poly = new myTria(pt1, pt2, pt3);

        poly->centroid(pt4);

        domTemp = 0;
        for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
        {
          if( ImmersedBodyObjects[bb]->within(pt4) )
            domTemp = bb+1;
        }

        vecTemp.push_back(domTemp);
        
        poly->SetDomainNumber(domTemp);

        nd1->subTrias.push_back(poly);

        //cout << " lllllllllllll " << endl;
      } //  for(cellId2=0; cellId2<polyDataLoc2->GetNumberOfCells(); cellId2++)

      vtkSmartPointer<vtkXMLPolyDataWriter>  writerPolyData2 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

      char fname1[50];
      sprintf(fname1,"%s%d%s", "delaunay3d-surf-", ee,".vtp");

      //Write the file.
      writerPolyData2->SetFileName(fname1);
      writerPolyData2->SetInputData(polyDataLoc2);
      writerPolyData2->Write();
      */


      // Construct the surface and create isosurface.	
      /*
      vtkSmartPointer<vtkSurfaceReconstructionFilter> surf = vtkSmartPointer<vtkSurfaceReconstructionFilter>::New();
      surf->SetInputData(polyDataLoc);

      vtkSmartPointer<vtkXMLImageDataWriter> writer =  vtkSmartPointer<vtkXMLImageDataWriter>::New();
      writer->SetFileName("surftemp.vti");
      writer->SetInputConnection(surf->GetOutputPort());
      writer->Write();
      */

      //
      delaunay3d->SetTolerance(0.0010);
      delaunay3d->SetOffset(1000.0);
      //delaunay3d->SetInputData(polyDataLoc);
      delaunay3d->SetInput(polyDataLoc);
      delaunay3d->Update();

      uGridLoc = delaunay3d->GetOutput();

      //cout << " Number of tetras = " << uGridLoc->GetNumberOfCells() << endl;

      for(cellId2=0; cellId2<uGridLoc->GetNumberOfCells(); cellId2++)
      {
        cellVTK2 = uGridLoc->GetCell(cellId2);

        uGridLoc->GetPoint(cellVTK2->GetPointId(0), &pt1[0]);
        uGridLoc->GetPoint(cellVTK2->GetPointId(1), &pt2[0]);
        uGridLoc->GetPoint(cellVTK2->GetPointId(2), &pt3[0]);
        uGridLoc->GetPoint(cellVTK2->GetPointId(3), &pt4[0]);

        myPoly  *poly = new myTet(pt1, pt2, pt3, pt4);

        poly->centroid(pt1);

        domTemp = 0;
        for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
        {
          if( ImmersedBodyObjects[bb]->within(pt4) )
            domTemp = bb+1;
        }

        vecTemp.push_back(domTemp);
        
        poly->SetDomainNumber(domTemp);

        ndTemp->subTrias.push_back(poly);
      } //  for(cellId2=0; cellId2<uGridLoc->GetNumberOfCells(); cellId2++)
      //

      /*
      typedef CGAL::Exact_predicates_inexact_constructions_kernel         K;
      typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned, K>    Vb;
      typedef CGAL::Triangulation_data_structure_3<Vb>                    Tds;
      //Use the Fast_location tag. Default or Compact_location works too.
      //typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location> Delaunay;
      typedef CGAL::Delaunay_triangulation_3<K, Tds>                      Delaunay;
      typedef Delaunay::Point  Point;

      std::vector<std::pair<Point, unsigned> > points;
      for(ii=0; ii<ptVec.size(); ii++)
        points.push_back( std::make_pair(Point(ptVec[ii][0], ptVec[ii][1], ptVec[ii][2]), ii) );

      //std::vector<Point>  points;
      //for(ii=0; ii<ptVec.size(); ii++)
        //points.push_back(Point(ptVec[ii][0], ptVec[ii][1], ptVec[ii][2]));

      Delaunay T( points.begin(),points.end() );

      std::ofstream oFileT("output",std::ios::out);
      // writing file output;
      oFileT << T;
      
      //cout << " Number of points = " << ptVec.size() << '\t' << T.number_of_vertices() << endl;
      //cout << " Number of tets = " << T.number_of_cells() << endl;

      ii=0;
      Delaunay::Finite_cells_iterator  cit;
      for(cit = T.finite_cells_begin(); cit != T.finite_cells_end(); ++cit)
      {
        //cout << " ii = " << ii++ << endl;

        //cout << T.tetrahedron(cit) << endl;
        //cout << T.tetrahedron(cit).vertex(0) << endl;
        
        //for(jj=0; jj<4; jj++) 
          //cout << jj << '\t' << cit->vertex(jj)->info() << '\t' << endl;

        pt1 = ptVec[cit->vertex(0)->info()];
        pt2 = ptVec[cit->vertex(1)->info()];
        pt3 = ptVec[cit->vertex(2)->info()];
        pt4 = ptVec[cit->vertex(3)->info()];

        myPoly  *poly = new myTet(pt1, pt2, pt3, pt4);

        poly->centroid(pt1);

        domTemp = 0;
        for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
        {
          if( ImmersedBodyObjects[bb]->within(pt4) )
            domTemp = bb+1;
        }

        vecTemp.push_back(domTemp);
        
        poly->SetDomainNumber(domTemp);

        ndTemp->subTrias.push_back(poly);

      } //  for(cit = T.finite_cells_begin(); cit != T.finite_cells_end(); ++cit)
      */

      //if( std::equal(vecTemp.begin()+1, vecTemp.end(), vecTemp.begin()) )
      //{
        //nd1->SetDomainNumber(vecTemp[0]);
        //nd1->clearSubtriangulation();
      //}
      //
      ndTemp->computeGaussPointsSubTrias(cutFEMparams[2]);
    }
    zz++;
  }  //  for(ee=0; ee<activeElements.size(); ee++)

  return;
}



void  HBSplineCutFEM::computeGaussPoints()
{
  if(CUTCELL_INTEGRATION_TYPE == 1)
  {
    if(ndm == 1)
      computeGaussPointsSubTrias1D();
    else if(ndm == 2)
      computeGaussPointsSubTrias2D();
    else
      computeGaussPointsSubTrias3D();
  }
  else
  {
    if(ndm == 1)
      computeGaussPointsAdapIntegration1D();
    else if(ndm == 2)
      computeGaussPointsAdapIntegration2D();
    else
      computeGaussPointsAdapIntegration3D();
  }

  return;
}


void  HBSplineCutFEM::computeGaussPointsSubTrias1D()
{
  return;
}

void  HBSplineCutFEM::computeGaussPointsAdapIntegration1D()
{
  return;
}


void  HBSplineCutFEM::computeGaussPointsSubTrias2D()
{
  return;
}


void  HBSplineCutFEM::computeGaussPointsAdapIntegration2D()
{
  return;
}



void  HBSplineCutFEM::computeGaussPointsSubTrias3D()
{
  return;
}



void  HBSplineCutFEM::computeGaussPointsAdapIntegration3D()
{
  return;
}






