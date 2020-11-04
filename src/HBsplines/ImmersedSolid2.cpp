#include "ImmersedSolid.h"
#include "SolutionData.h"
#include "GeomDataLagrange.h"
#include "ImmersedIntegrationElement.h"
#include "myLine.h"
#include "myTria.h"
#include "myQuad.h"
#include "Files.h"
#include "util.h"


extern Files files;

using namespace myGeom;

//#define CUTCELL3D_ALGO_VTK


void ImmersedSolid::adjustBoundaryPoints(double* minVal, double* maxVal)
{
    //cout << " ImmersedSolid::adjustBoundaryPoints .... needs to be modified " << endl;
//#ifdef CUTCELL3D_ALGO_VTK
    int ii=0, jj=0;
    double  tol1=1.0e-8, tol2=1.0e-4;

    myPoint  pt1, pt2;

    for(ii=0; ii<GeomData.NodePosOrig.size(); ii++)
    {
      pt1 = GeomData.NodePosOrig[ii];

      for(jj=0; jj<DIM; jj++)
      {
        if( abs(pt1[jj]-minVal[jj]) < tol1 )
        {
          GeomData.NodePosOrig[ii][jj] -= tol2;
          GeomData.NodePosCur[ii][jj]  -= tol2;
        }

        if( abs(pt1[jj]-maxVal[jj]) < tol1 )
        {
          GeomData.NodePosOrig[ii][jj] += tol2;
          GeomData.NodePosCur[ii][jj]  += tol2;
        }
      }
    }
//#endif
    return;
}


/*
void  ImmersedSolid::setImmersedFaces()
{
  // using VTK library
  //cout << "  ImmersedSolid::setImmersedFaces()  " << endl;

  // set the immersed faces for cutFEM purposes

  myPoly* poly;
  myPoint  pt1, pt2, pt3, pt4, normal;
  int ii, jj, kk;

  vtkSmartPointer<vtkPoints>     pointsVTK   =  vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkLine>       lineVTK     =  vtkSmartPointer<vtkLine>::New();
  vtkSmartPointer<vtkTriangle>   triaVTK     =  vtkSmartPointer<vtkTriangle>::New();
  vtkSmartPointer<vtkQuad>       quadVTK     =  vtkSmartPointer<vtkQuad>::New();
  vtkSmartPointer<vtkCellArray>  polyList    =  vtkSmartPointer<vtkCellArray>::New();
  vtkSmartPointer<vtkVertex>     vertexVTK2  =  vtkSmartPointer<vtkVertex>::New();
  vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK2  = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkFloatArray>   triaStatus = vtkSmartPointer<vtkFloatArray>::New();
  vtkSmartPointer<vtkFloatArray>   triaNormal = vtkSmartPointer<vtkFloatArray>::New();

  vtkIdType  pts[10];

  //cout << " nImmInt = " << nImmInt << '\t' << GeomData.NodePosOrig.size() << '\t' << DIM << endl;

  if(DIM == 2)
  {
    for(ii=0; ii<GeomData.NodePosOrig.size(); ii++)
    {
      pt1 = GeomData.NodePosOrig[ii];
      pts[0] = pointsVTK->InsertNextPoint(pt1[0], pt1[1], 0.0);

      vertexVTK2->GetPointIds()->SetId(0, pts[0]);

      uGridVTK2->InsertNextCell(vertexVTK2->GetCellType(), vertexVTK2->GetPointIds());
    }

    for(ii=0; ii<nImmInt; ii++)
    {
      pt1 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[0]];
      pt2 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[1]];

      poly = new myLine(pt1, pt2);

      ImmersedFaces.push_back(poly);

      lineVTK->GetPointIds()->SetId(0, ImmIntgElems[ii]->pointNums[0]);
      lineVTK->GetPointIds()->SetId(1, ImmIntgElems[ii]->pointNums[1]);

      polyList->InsertNextCell(lineVTK);
      
      uGridVTK2->InsertNextCell(lineVTK->GetCellType(), lineVTK->GetPointIds());

      //vertexVTK2->GetPointIds()->SetId(0, ii);

      //polyList->InsertNextCell(vertexVTK2);
    }

    char fname[200];
    sprintf(fname,"%s%s%d%s", files.Ofile.asCharArray(), "-immersedsolid-", id,".vtu");

    uGridVTK2->SetPoints(pointsVTK);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    writerUGridVTK->SetFileName(fname);
    writerUGridVTK->SetInputData(uGridVTK2);
    writerUGridVTK->Write();

    // The area of a convex polygon is 
    // positive if the points are arranged in a counterclockwise order, and
    // negative if they are in clockwise order

    double area = 0.0;
    for(kk=0; kk<nImmInt; kk++)
    {
      ii = ImmIntgElems[kk]->pointNums[0];
      jj = ImmIntgElems[kk]->pointNums[1];
      area += GeomData.NodePosOrig[ii][0] * GeomData.NodePosOrig[jj][1];
      area -= GeomData.NodePosOrig[jj][0] * GeomData.NodePosOrig[ii][1];
    }

    area = 0.5*area;

    printf("\n Polygon # %d area = %16.12f \n\n", id, area );
  }

  if(DIM == 3)
  {
    double  vec[3];
    triaNormal->SetNumberOfComponents(3);

    for(ii=0; ii<GeomData.NodePosOrig.size(); ii++)
    {
      pt1 = GeomData.NodePosOrig[ii];

      pts[0] = pointsVTK->InsertNextPoint(pt1[0], pt1[1], pt1[2]);
    }

    for(ii=0; ii<nImmInt; ii++)
    {
      if(ImmIntgElems[ii]->pointNums.size() == 3) // triangle
      {
        pt1 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[0]];
        pt2 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[1]];
        pt3 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[2]];

        poly = new myTria(pt1, pt2, pt3);

        ImmersedFaces.push_back(poly);

        triaVTK->GetPointIds()->SetId(0, ImmIntgElems[ii]->pointNums[0]);
        triaVTK->GetPointIds()->SetId(1, ImmIntgElems[ii]->pointNums[1]);
        triaVTK->GetPointIds()->SetId(2, ImmIntgElems[ii]->pointNums[2]);

        polyList->InsertNextCell(triaVTK);

        triaStatus->InsertNextValue(ImmIntgElems[ii]->isActive());

        poly->computeNormal();

        poly->computeNormal(pt1, normal);

        vec[0] = normal[0];  vec[1] = normal[1];  vec[2] = normal[2];

        triaNormal->InsertNextTuple(vec);
      }
      else if(ImmIntgElems[ii]->pointNums.size() == 4) // quad
      {
        pt1 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[0]];
        pt2 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[1]];
        pt3 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[2]];
        pt4 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[3]];

        poly = new myQuad(pt1, pt2, pt3, pt4);

        ImmersedFaces.push_back(poly);

        quadVTK->GetPointIds()->SetId(0, ImmIntgElems[ii]->pointNums[0]);
        quadVTK->GetPointIds()->SetId(1, ImmIntgElems[ii]->pointNums[1]);
        quadVTK->GetPointIds()->SetId(2, ImmIntgElems[ii]->pointNums[3]);
        quadVTK->GetPointIds()->SetId(3, ImmIntgElems[ii]->pointNums[2]);

        polyList->InsertNextCell(quadVTK);

        triaStatus->InsertNextValue(ImmIntgElems[ii]->isActive());

        poly->computeNormal();

        poly->computeNormal(pt1, normal);

        vec[0] = normal[0];  vec[1] = normal[1];  vec[2] = normal[2];

        triaNormal->InsertNextTuple(vec);
      }
      else
      {
        cout << " Error in ImmersedSolid::setImmersedFaces()  " << endl;
      }
    }

    polyDataVTK->SetPoints(pointsVTK);
    polyDataVTK->SetPolys(polyList);

    triaStatus->SetName("Active");
    triaNormal->SetName("normal");

    polyDataVTK->GetCellData()->SetScalars(triaStatus);
    polyDataVTK->GetCellData()->SetVectors(triaNormal);

    vtkSmartPointer<vtkXMLPolyDataWriter>  writerPolyData =     vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    //vtkSmartPointer<vtkMassProperties>  massProp =     vtkSmartPointer<vtkMassProperties>::New();

    char fname1[200];
    sprintf(fname1,"%s%s%d%s", files.Ofile.asCharArray(), "-immersedsolid-", id,".vtp");

    writerPolyData->SetFileName(fname1);
    writerPolyData->SetInputData(polyDataVTK);
    selectEnclosedPoints->SetSurfaceData(polyDataVTK);
    //massProp->SetInputData(polyDataVTK);
    writerPolyData->Write();

    selectEnclosedPoints->SetTolerance(0.000001);
    //selectEnclosedPoints->Update();
    selectEnclosedPoints->Initialize(polyDataVTK);

    //massProp->Update();

    //printf("\n Polyhedron # %d volume = %12.8f \n\n", id, massProp->getVolume() );
    //printf("\n Polyhedron # %d Surface area = %12.8f \n\n", id, massProp->GetSurfaceArea() );
  }

  return;
}


*/


//
// set the immersed faces for cutFEM purposes
void  ImmersedSolid::setImmersedFaces()
{
    //cout << "  ImmersedSolid::setImmersedFaces()  " << endl;

    //////////////////
    // CGAL
    //////////////////

    myPoly* poly;
    myPoint  pt1, pt2, pt3, pt4, normal;
    int ii=0, jj=0, kk=0;

    vtkSmartPointer<vtkPoints>     pointsVTK   =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkLine>       lineVTK     =  vtkSmartPointer<vtkLine>::New();
    vtkSmartPointer<vtkTriangle>   triaVTK     =  vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkQuad>       quadVTK     =  vtkSmartPointer<vtkQuad>::New();
    vtkSmartPointer<vtkCellArray>  polyList    =  vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkVertex>     vertexVTK2  =  vtkSmartPointer<vtkVertex>::New();
    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK2  = vtkSmartPointer<vtkUnstructuredGrid>::New();

    vtkSmartPointer<vtkFloatArray>   triaStatus = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>   triaNormal = vtkSmartPointer<vtkFloatArray>::New();

    vtkIdType  pts[10];

    //cout << " nImmInt = " << nImmInt << '\t' << GeomData.NodePosOrig.size() << '\t' << DIM << endl;

    if(DIM == 2)
    {
      for(ii=0; ii<GeomData.NodePosOrig.size(); ii++)
      {
        pt1 = GeomData.NodePosOrig[ii];

        pts[0] = pointsVTK->InsertNextPoint(pt1[0], pt1[1], 0.0);

        vertexVTK2->GetPointIds()->SetId(0, pts[0]);

        uGridVTK2->InsertNextCell(vertexVTK2->GetCellType(), vertexVTK2->GetPointIds());
      }

      for(ii=0; ii<nImmInt; ii++)
      {
        pt1 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[0]];
        pt2 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[1]];

        poly = new myLine(pt1, pt2);

        ImmersedFaces.push_back(poly);

        lineVTK->GetPointIds()->SetId(0, ImmIntgElems[ii]->pointNums[0]);
        lineVTK->GetPointIds()->SetId(1, ImmIntgElems[ii]->pointNums[1]);

        polyList->InsertNextCell(lineVTK);

        uGridVTK2->InsertNextCell(lineVTK->GetCellType(), lineVTK->GetPointIds());

        //vertexVTK2->GetPointIds()->SetId(0, ii);

        //polyList->InsertNextCell(vertexVTK2);
      }

      char fname[500];
      sprintf(fname,"%s%s%s%s%d%s", files.projDir.asCharArray(), "/", files.Ofile.asCharArray(), "-immersedsolid-", id,".vtu");

      uGridVTK2->SetPoints(pointsVTK);

      vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

      writerUGridVTK->SetFileName(fname);
      writerUGridVTK->SetInputData(uGridVTK2);
      writerUGridVTK->Write();

      // The area of a convex polygon is 
      // positive if the points are arranged in a counterclockwise order, and
      // negative if they are in clockwise order

      double area = 0.0;
      for(kk=0; kk<nImmInt; kk++)
      {
        ii = ImmIntgElems[kk]->pointNums[0];
        jj = ImmIntgElems[kk]->pointNums[1];
        area += GeomData.NodePosOrig[ii][0] * GeomData.NodePosOrig[jj][1];
        area -= GeomData.NodePosOrig[jj][0] * GeomData.NodePosOrig[ii][1];
      }
      area = 0.5*area;
 
      //printf("\n Polygon # %d area = %16.12f \n\n", id, area );
    }

    if(DIM == 3)
    {
        int  nodesPerFace = ImmIntgElems[0]->pointNums.size(); // assumed that all the faces are of same type - triangle or quadrilateral
        double  vec[3];
        triaNormal->SetNumberOfComponents(3);

        pointsCGAL.clear();

        for(ii=0; ii<GeomData.NodePosOrig.size(); ii++)
        {
          pt1 = GeomData.NodePosOrig[ii];

          pts[0] = pointsVTK->InsertNextPoint(pt1[0], pt1[1], pt1[2]);

          pointsCGAL.push_back(CGAL_Point(pt1[0], pt1[1], pt1[2]));
        }

        facesCGAL.clear();

        for(ii=0; ii<nImmInt; ii++)
        {
          for(jj=0; jj<nodesPerFace; jj++)
            facesCGAL.push_back(ImmIntgElems[ii]->pointNums[jj]);

          if(nodesPerFace == 3) // triangle
          {
            pt1 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[0]];
            pt2 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[1]];
            pt3 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[2]];

            poly = new myTria(pt1, pt2, pt3);

            ImmersedFaces.push_back(poly);

            triaVTK->GetPointIds()->SetId(0, ImmIntgElems[ii]->pointNums[0]);
            triaVTK->GetPointIds()->SetId(1, ImmIntgElems[ii]->pointNums[1]);
            triaVTK->GetPointIds()->SetId(2, ImmIntgElems[ii]->pointNums[2]);

            polyList->InsertNextCell(triaVTK);

            triaStatus->InsertNextValue(ImmIntgElems[ii]->isActive());

            poly->computeNormal();

            poly->computeNormal(pt1, normal);

            vec[0] = normal[0];  vec[1] = normal[1];  vec[2] = normal[2];

            triaNormal->InsertNextTuple(vec);
          }
          else if(nodesPerFace == 4) // quad
          {
            pt1 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[0]];
            pt2 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[1]];
            pt3 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[2]];
            pt4 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[3]];

            poly = new myQuad(pt1, pt2, pt3, pt4);

            ImmersedFaces.push_back(poly);

            quadVTK->GetPointIds()->SetId(0, ImmIntgElems[ii]->pointNums[0]);
            quadVTK->GetPointIds()->SetId(1, ImmIntgElems[ii]->pointNums[1]);
            quadVTK->GetPointIds()->SetId(2, ImmIntgElems[ii]->pointNums[3]);
            quadVTK->GetPointIds()->SetId(3, ImmIntgElems[ii]->pointNums[2]);

            polyList->InsertNextCell(quadVTK);

            triaStatus->InsertNextValue(ImmIntgElems[ii]->isActive());

            poly->computeNormal();

            poly->computeNormal(pt1, normal);

            vec[0] = normal[0];  vec[1] = normal[1];  vec[2] = normal[2];

            triaNormal->InsertNextTuple(vec);
          }
          else
          {
            cout << " Error in ImmersedSolid::setImmersedFaces()  " << endl;
          }
        }

        treeCGAL.clear();

        // build a polyhedron from the loaded arrays
        //CGAL_Polyhedron  polyhedronTemp;
        polyhedron_builder<HalfedgeDS>  poly_builder( pointsCGAL, facesCGAL, nodesPerFace );
        polyhedronTemp.delegate( poly_builder );

        treeCGAL.insert(polyhedronTemp.facets_begin(), polyhedronTemp.facets_end(), polyhedronTemp);
        treeCGAL.accelerate_distance_queries();

        //CGAL_Point_inside point_inside_tester(treeCGAL);

        point_inside_tester = new CGAL_Point_inside(treeCGAL);

        //PetscPrintf(MPI_COMM_WORLD, " treeCGAL size = %d, \n", treeCGAL.size() );

        //if( checkBoundedSideCGAL( (*point_inside_tester)(CGAL_Point(0.5, 0.2, 0.2) ) )
        //cout << " 11111 ---- 1111111111 " << endl; }
        //res = (*point_inside_tester)(CGAL_Point(0.458333,0.164,0.0205));

        polyDataVTK->SetPoints(pointsVTK);
        polyDataVTK->SetPolys(polyList);

        triaStatus->SetName("Active");
        triaNormal->SetName("normal");

        polyDataVTK->GetCellData()->SetScalars(triaStatus);
        polyDataVTK->GetCellData()->SetVectors(triaNormal);

        vtkSmartPointer<vtkXMLPolyDataWriter>  writerPolyData =     vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        //vtkSmartPointer<vtkMassProperties>  massProp =     vtkSmartPointer<vtkMassProperties>::New();

        char fname1[200];
        sprintf(fname1,"%s%s%s%s%d%s", files.projDir.asCharArray(), "/", files.Ofile.asCharArray(), "-immersedsolid-", id,".vtp");

        writerPolyData->SetFileName(fname1);
        writerPolyData->SetInputData(polyDataVTK);
        //massProp->SetInputData(polyDataVTK);
        writerPolyData->Write();
    }

    return;
}


void ImmersedSolid::updateImmersedFaces()
{
    //cout << " ImmersedSolid::updateImmersedFaces() ... " << endl;

    vtkSmartPointer<vtkPoints>     pointsVTK   =  vtkSmartPointer<vtkPoints>::New();

    int ii=0, jj=0;

    myPoint  pt1;
    vtkIdType  ptId;

    pointsCGAL.clear();

    for(ii=0; ii<GeomData.NodePosCur.size(); ii++)
    {
      pt1 = GeomData.NodePosCur[ii];

      ptId = pointsVTK->InsertNextPoint(pt1[0], pt1[1], pt1[2]);

      pointsCGAL.push_back(CGAL_Point(pt1[0], pt1[1], pt1[2]));
    }

    int  nodesPerFace = ImmIntgElems[0]->pointNums.size(); // assumed that all the faces are of same type - triangle or quadrilateral

    if(DIM == 3)
    {
        polyDataVTK->SetPoints(pointsVTK);

#ifdef CUTCELL3D_ALGO_VTK
        selectEnclosedPoints->Initialize(polyDataVTK);
#else

        // build a polyhedron from the loaded arrays
        //CGAL_Polyhedron  polyhedronTemp;
        polyhedronTemp.clear();
        polyhedron_builder<HalfedgeDS>  poly_builder( pointsCGAL, facesCGAL, nodesPerFace );
        polyhedronTemp.delegate( poly_builder );

        treeCGAL.clear();
        //PetscPrintf(MPI_COMM_WORLD, " treeCGAL size = %d \n", treeCGAL.size() );
        treeCGAL.insert(polyhedronTemp.facets_begin(), polyhedronTemp.facets_end(), polyhedronTemp);
        treeCGAL.accelerate_distance_queries();

        //PetscPrintf(MPI_COMM_WORLD, " treeCGAL size = %d \n", treeCGAL.size() );

        if(point_inside_tester != NULL)
          delete  point_inside_tester;

        point_inside_tester = NULL;

        point_inside_tester = new CGAL_Point_inside(treeCGAL);

        //cout << "checking = " << checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(0.0, 0.0, 0.0))) << endl;
#endif
    } // if(DIM == 3)

    for(ii=0; ii<nImmInt; ii++)
    {
      for(jj=0; jj<nodesPerFace; jj++)
        ImmersedFaces[ii]->updatePoint(jj, GeomData.NodePosCur[ImmIntgElems[ii]->pointNums[jj]] );

      ImmersedFaces[ii]->computeNormal();
      ImmersedFaces[ii]->computeAABB();
    }

    return;
}
//



int ImmersedSolid::within(myPoint& ptTemp)
{
  // check if the point is outside/inside the boundingbox of the polygon
  // if the point is outside the boundingbox then
  // the point is outside the polygon
  // otherwise, the point may be inside the polygon and needs further checking

  if( !(bbox.within(ptTemp)) )
    return 0;

  int i=0, c=0;

  if(DIM == 2)
  {
    ray1.updateOrigin(ptTemp);

    vector<myPoint>  ptOut;

    c=0;
    for(i=0; i<ImmersedFaces.size(); i++)
    {
      if( ImmersedFaces[i]->IntersectWithRay(ray1, ptOut ) )
        c = !c;
    }
  }
  else
  {
#ifdef CUTCELL3D_ALGO_VTK
    c =  selectEnclosedPoints->IsInsideSurface( ptTemp[0], ptTemp[1], ptTemp[2] );
#else
    c = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(ptTemp[0], ptTemp[1], ptTemp[2])));
#endif
  }

  return c;
}


int  ImmersedSolid::doAABBintersect(AABB&  bb2)
{
  if( bbox.doIntersect(bb2) )
    return 1;
  else
    return 0;
}



int  ImmersedSolid::doIntersect2D(AABB&  bbTemp, bool flag, vector<int>& vecCorners, vector<myPoint>&  ptOut)
{
    // code where both corners and edges are checked

    // return value = 0       ---> the cell is the domain '0' which is the default background grid
    // return value = -1      ---> the cell is cut by the current immersed body
    // return value = (bb+1)  ---> the cell is completely inside the immersed body 'bb'

    // this method is used to test whether 
    // a background grid cell is cut by the immersed body (IB)
    // 
    // bbTemp is the boundingbox of the cell to be tested

    // check if the boundingbox 'bb2' intersects the boundingbox of the IB
    // if the boundingbox 'bb2' does not intersect the boundingbox of the IB then
    // the IB does not intersect the cell

    // if the boundingbox 'bb2' intersects the boundingbox of the IB then
    // the IB may intersect the cell and needs further testing
  
    if( !(bbox.doIntersect(bbTemp)) )
      return 0;

    int  ii=0, jj=0, c=0, pp=0;
    myPoint  ptTemp;
    vector<myPoint>  ptVec(4);

    vector<int>  vecEdges(4);

    ptVec[0][0] = bbTemp.minBB[0];    ptVec[0][1] = bbTemp.minBB[1];    ptVec[0][2] = 0.0;
    ptVec[1][0] = bbTemp.maxBB[0];    ptVec[1][1] = bbTemp.minBB[1];    ptVec[1][2] = 0.0;
    ptVec[2][0] = bbTemp.minBB[0];    ptVec[2][1] = bbTemp.maxBB[1];    ptVec[2][2] = 0.0;
    ptVec[3][0] = bbTemp.maxBB[0];    ptVec[3][1] = bbTemp.maxBB[1];    ptVec[3][2] = 0.0;

    for(pp=0; pp<4; pp++)
    {
      ptTemp = ptVec[pp];

      ray1.updateOrigin(ptTemp);

      c=0;
      for(ii=0; ii<ImmersedFaces.size(); ii++)
      {
        if( ImmersedFaces[ii]->IntersectWithRay(ray1, ptOut ) )
          c = !c;
      }
      vecCorners[pp] = c;
    } // for(pp=0; pp<4; pp++)

    //printVector(vecTemp);

    for(ii=0; ii<4; ii++)
    {
      //if(vecTemp[ii] == 0)
        ptOut.push_back(ptVec[ii]);
    }

    //int verts[4][2] ={ {0,1}, {1,3}, {3,2}, {2,0} };
    ////int verts[4][2] ={ {0,1}, {1,2}, {2,3}, {3,0} };

    vector<myLine>  cellLines;

    for(ii=0; ii<4; ii++)
    {
      cellLines.push_back( myLine(ptVec[verts[ii][0]], ptVec[verts[ii][1]]) );
      cellLines[ii].computeNormal();
      cellLines[ii].computeAABB();
    }

    for(ii=0; ii<ImmersedFaces.size(); ii++)
    {
      // check if the boundingbox intersects any of the line segments
      // if it intersects then find the intersection points

      if( ImmersedFaces[ii]->doIntersect(bbTemp) )
      {
        for(jj=0; jj<2; jj++)
        {
          ptTemp = ImmersedFaces[ii]->GetPoint(jj);
          if( bbTemp.within(ptTemp) )
            ptOut.push_back(ptTemp);
        }

        for(jj=0; jj<4; jj++)
          vecEdges[jj] = ImmersedFaces[ii]->IntersectWithLine(cellLines[jj], ptOut);
      }
    } // for(ii=0; ii<ImmersedFaces.size(); ii++)
    cellLines.clear();

    int val = -1;
    if( std::equal(vecCorners.begin()+1, vecCorners.end(), vecCorners.begin()) )
    {
      //if ( std::any_of(vecEdges.begin(), vecEdges.end(), [](int i){return i==1;}) )
      if( my_any_of(vecEdges, 1) )
      {
        val = -1;
      }
      else
      {
        if( vecCorners[0] == 1 )
          val = id+1;
        else
          val = 0;
      }
    }

    return val;
}




int  ImmersedSolid::doIntersect2Dfor3D(int sideTemp, double coord3, AABB&  bbTemp, bool flag, vector<int>& vecTemp, vector<myPoint>&  ptOut)
{
    // return value = 0       ---> the cell is the domain '0' which is the default background grid
    // return value = -1      ---> the cell is cut by the current immersed body
    // return value = (bb+1)  ---> the cell is completely inside the immersed body 'bb'
  
    // this method is used to test whether 
    // a background grid cell is cut by the immersed body (IB)
    // 
    // bbTemp is the boundingbox of the cell to be tested

    // check if the boundingbox 'bbTemp' intersects the boundingbox of the IB
    // if the boundingbox 'bbTemp' does not intersect the boundingbox of the IB then
    // the IB does not intersect the cell

    // if the boundingbox 'bbTemp' intersects the boundingbox of the IB then
    // the IB may intersect the cell and needs further testing

    // this subroutine is for faces of boundary elements in 3D

    if( !(bbox.doIntersect(bbTemp)) )
      return 0;

    //bbTemp.printSelf();
    //cout << " coord3 = " << coord3 << endl;
    //cout << " sideTemp = " << sideTemp << endl;

    int  val=0;

#ifdef CUTCELL3D_ALGO_VTK
    // using VTK subroutines

    if(sideTemp == 0 || sideTemp == 1)
    {
      vecTemp[0] = selectEnclosedPoints->IsInsideSurface( coord3, bbTemp.minBB[1], bbTemp.minBB[2] );
      vecTemp[1] = selectEnclosedPoints->IsInsideSurface( coord3, bbTemp.maxBB[1], bbTemp.minBB[2] );
      vecTemp[2] = selectEnclosedPoints->IsInsideSurface( coord3, bbTemp.minBB[1], bbTemp.maxBB[2] );
      vecTemp[3] = selectEnclosedPoints->IsInsideSurface( coord3, bbTemp.maxBB[1], bbTemp.maxBB[2] );
    }
    else if(sideTemp == 2 || sideTemp == 3)
    {
      vecTemp[0] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], coord3, bbTemp.minBB[2] );
      vecTemp[1] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], coord3, bbTemp.minBB[2] );
      vecTemp[2] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], coord3, bbTemp.maxBB[2] );
      vecTemp[3] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], coord3, bbTemp.maxBB[2] );
    }
    else if(sideTemp == 4 || sideTemp == 5)
    {
      vecTemp[0] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], bbTemp.minBB[1], coord3 );
      vecTemp[1] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], bbTemp.minBB[1], coord3 );
      vecTemp[2] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], bbTemp.maxBB[1], coord3 );
      vecTemp[3] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], bbTemp.maxBB[1], coord3 );
    }
    //
#else

    // using CGAL subroutines

    if(sideTemp == 0 || sideTemp == 1)
    {
      vecTemp[0] = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(coord3, bbTemp.minBB[1], bbTemp.minBB[2])));
      vecTemp[1] = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(coord3, bbTemp.maxBB[1], bbTemp.minBB[2])));
      vecTemp[2] = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(coord3, bbTemp.minBB[1], bbTemp.maxBB[2])));
      vecTemp[3] = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(coord3, bbTemp.maxBB[1], bbTemp.maxBB[2])));
    }
    else if(sideTemp == 2 || sideTemp == 3)
    {
      vecTemp[0] = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(bbTemp.minBB[0], coord3, bbTemp.minBB[2])));
      vecTemp[1] = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(bbTemp.maxBB[0], coord3, bbTemp.minBB[2])));
      vecTemp[2] = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(bbTemp.minBB[0], coord3, bbTemp.maxBB[2])));
      vecTemp[3] = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(bbTemp.maxBB[0], coord3, bbTemp.maxBB[2])));
    }
    else if(sideTemp == 4 || sideTemp == 5)
    {
      vecTemp[0] = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(bbTemp.minBB[0], bbTemp.minBB[1], coord3)));
      vecTemp[1] = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(bbTemp.maxBB[0], bbTemp.minBB[1], coord3)));
      vecTemp[2] = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(bbTemp.minBB[0], bbTemp.maxBB[1], coord3)));
      vecTemp[3] = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(bbTemp.maxBB[0], bbTemp.maxBB[1], coord3)));
    }

#endif
    //printVector(vecTemp);

    val = -1;
    if( std::equal(vecTemp.begin()+1, vecTemp.end(), vecTemp.begin()) )
    {
      if( vecTemp[0] == 1 )
        val = id+1;
      else
        val = 0;
    }

    return val;
}




int  ImmersedSolid::doIntersect3D(AABB&  bbTemp, bool flag, vector<int>& vecTemp, vector<myPoint>&  ptOut)
{
  // this method is used to test whether 
  // a background grid cell is cut by the immersed body (IB)
  // 
  // bbTemp is the boundingbox of the cell to be tested

  // return value = 0       ---> the cell is the domain '0' which is the default background grid
  // return value = -1      ---> the cell is cut by the current immersed body
  // return value = (bb+1)  ---> the cell is completely inside the immersed body 'bb'

  // check if the boundingbox 'bbTemp' intersects the boundingbox of the IB
  // If the boundingbox 'bbTemp' does not intersect the boundingbox of the IB then
  // the IB does not intersect the cell

    if( !(bbox.doIntersect(bbTemp)) )
      return 0;

    // If the boundingbox of the cell intersects the boundingbox of the IB then
    // the IB may intersect the cell and needs further testing

#ifdef CUTCELL3D_ALGO_VTK

    // using VTK subroutines
    
    vecTemp[0] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.minBB[2] );
    vecTemp[1] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.minBB[2] );
    vecTemp[2] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.minBB[2] );
    vecTemp[3] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.minBB[2] );

    vecTemp[4] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.maxBB[2] );
    vecTemp[5] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.maxBB[2] );
    vecTemp[6] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2] );
    vecTemp[7] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2] );

#else

    // using CGAL subroutines

    CGAL::Bounded_side  res;

    vecTemp[0] = checkBoundedSideCGAL( (*point_inside_tester)
                  (CGAL_Point(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.minBB[2])) );

    vecTemp[1] = checkBoundedSideCGAL( (*point_inside_tester)
                  (CGAL_Point(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.minBB[2])) );

    vecTemp[2] = checkBoundedSideCGAL( (*point_inside_tester)
                  (CGAL_Point(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.minBB[2])) );

    vecTemp[3] = checkBoundedSideCGAL( (*point_inside_tester)
                  (CGAL_Point(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.minBB[2])) );

    vecTemp[4] = checkBoundedSideCGAL( (*point_inside_tester)
                  (CGAL_Point(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.maxBB[2])) );
    
    vecTemp[5] = checkBoundedSideCGAL( (*point_inside_tester)
                  (CGAL_Point(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.maxBB[2])) );

    vecTemp[6] = checkBoundedSideCGAL( (*point_inside_tester)
                  (CGAL_Point(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2])) );

    vecTemp[7] = checkBoundedSideCGAL( (*point_inside_tester)
                  (CGAL_Point(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2])) );
#endif

    //printVector(vecTemp);
    int val = -1;
    if( std::equal(vecTemp.begin()+1, vecTemp.end(), vecTemp.begin()) )
    {
      if( vecTemp[0] == 1 )
        val = id+1;
      else
        val = 0;
    }

    return val;
}



double  ImmersedSolid::distanceFromPoint(myPoint&  pt)
{
  return 0.0;
}


double  ImmersedSolid::distanceFromPoint(double xx, double yy, double zz)
{
  return 0.0;
}



