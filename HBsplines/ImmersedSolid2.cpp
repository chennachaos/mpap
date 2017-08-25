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



void ImmersedSolid::adjustBoundaryPoints(double* minVal, double* maxVal)
{
    int ii=0, jj=0;
    double  tol1=1.0e-8, tol2=1.0e-6;

    myPoint  pt1, pt2;

    for(ii=0; ii<GeomData.NodePosOrig.size(); ii++)
    {
      pt1 = GeomData.NodePosOrig[ii];
      //cout << pt1[0] << '\t' << pt1[1] << '\t' << pt1[2] << endl;
      
      for(jj=0; jj<DIM; jj++)
      {
        if( abs(pt1[jj]-minVal[jj]) < tol1 )
        //if( CompareDoubles(pt1[jj], minVal[jj]) )
        {
          //cout << pt1[0] << '\t' << pt1[1] << '\t' << pt1[2] << endl;
          GeomData.NodePosOrig[ii][jj] -= tol2;
          GeomData.NodePosCur[ii][jj]  -= tol2;
        }

        if( abs(pt1[jj]-maxVal[jj]) < tol1 )
        //if( CompareDoubles(pt1[jj], maxVal[jj]) )
        {
          GeomData.NodePosOrig[ii][jj] += tol2;
          GeomData.NodePosCur[ii][jj]  += tol2;
        }
      }
    }

    return;
}


/*
void  ImmersedSolid::setImmersedFaces()
{
  cout << "  ImmersedSolid::setImmersedFaces()  " << endl;

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

  cout << " nImmInt = " << nImmInt << '\t' << GeomData.NodePosOrig.size() << '\t' << DIM << endl;

  if(DIM == 2)
  {
    for(ii=0; ii<GeomData.NodePosOrig.size(); ii++)
    {
      pt1 = GeomData.NodePosOrig[ii];
      //cout << pt1[0] << '\t' << pt1[1] << '\t' << pt1[2] << endl;
      pts[0] = pointsVTK->InsertNextPoint(pt1[0], pt1[1], 0.0);

      vertexVTK2->GetPointIds()->SetId(0, pts[0]);

      uGridVTK2->InsertNextCell(vertexVTK2->GetCellType(), vertexVTK2->GetPointIds());
    }

    for(ii=0; ii<nImmInt; ii++)
    {
      pt1 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[0]];
      pt2 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[1]];

      //cout << " ii = " << ii << '\t' << ImmIntgElems[ii]->pointNums[0] << '\t' << ImmIntgElems[ii]->pointNums[1] << endl;

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

#if VTK_MAJOR_VERSION == 5
    writerUGridVTK->SetInput(uGridVTK2);
#else
    writerUGridVTK->SetInputData(uGridVTK2);
#endif

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
      //cout << pt1[0] << '\t' << pt1[1] << '\t' << pt1[2] << endl;
      pts[0] = pointsVTK->InsertNextPoint(pt1[0], pt1[1], pt1[2]);
    }

    for(ii=0; ii<nImmInt; ii++)
    {
      //cout << ii << '\t' << ImmIntgElems[ii]->pointNums.size() << endl;

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

#if VTK_MAJOR_VERSION == 5
    writerPolyData->SetInput(polyDataVTK);
    selectEnclosedPoints->SetSurface(polyDataVTK);
    //massProp->SetInput(polyDataVTK);
#else
    writerPolyData->SetInputData(polyDataVTK);
    selectEnclosedPoints->SetSurfaceData(polyDataVTK);
    //massProp->SetInputData(polyDataVTK);
#endif

    writerPolyData->Write();

    //selectEnclosedPoints->Update();
    selectEnclosedPoints->Initialize(polyDataVTK);

    //massProp->Update();

    //printf("\n Polyhedron # %d volume = %12.8f \n\n", id, massProp->getVolume() );
    //printf("\n Polyhedron # %d Surface area = %12.8f \n\n", id, massProp->GetSurfaceArea() );
  }

  return;
}


void ImmersedSolid::updateImmersedFaces()
{
  //if(totalDOF == 0)
    //return;

  //cout << " totalDOF = " << totalDOF << endl;
  //cout << " nImmInt = " << nImmInt << endl;

  vtkSmartPointer<vtkPoints>     pointsVTK   =  vtkSmartPointer<vtkPoints>::New();

  int ii, jj;

  myPoint  pt1;
  vtkIdType  ptId;

  for(ii=0; ii<GeomData.NodePosCur.size(); ii++)
  {
    pt1 = GeomData.NodePosCur[ii];
    //cout << pt1[0] << '\t' << pt1[1] << '\t' << pt1[2] << endl;
    ptId = pointsVTK->InsertNextPoint(pt1[0], pt1[1], pt1[2]);
  }

  if(DIM == 3)
  {
    polyDataVTK->SetPoints(pointsVTK);

    selectEnclosedPoints->Initialize(polyDataVTK);
  }

  for(ii=0; ii<nImmInt; ii++)
  {
    for(jj=0; jj<ImmIntgElems[ii]->pointNums.size(); jj++)
      ImmersedFaces[ii]->updatePoint(jj, GeomData.NodePosCur[ImmIntgElems[ii]->pointNums[jj]] );

    ImmersedFaces[ii]->computeNormal();
    ImmersedFaces[ii]->computeAABB();
  }

  return;
}

*/



void  ImmersedSolid::setImmersedFaces()
{
  //cout << "  ImmersedSolid::setImmersedFaces()  " << endl;
  // set the immersed faces for cutFEM purposes

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
      //cout << pt1[0] << '\t' << pt1[1] << '\t' << pt1[2] << endl;
      pts[0] = pointsVTK->InsertNextPoint(pt1[0], pt1[1], 0.0);

      vertexVTK2->GetPointIds()->SetId(0, pts[0]);

      uGridVTK2->InsertNextCell(vertexVTK2->GetCellType(), vertexVTK2->GetPointIds());
    }

    for(ii=0; ii<nImmInt; ii++)
    {
      pt1 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[0]];
      pt2 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[1]];

      //cout << " ii = " << ii << '\t' << ImmIntgElems[ii]->pointNums[0] << '\t' << ImmIntgElems[ii]->pointNums[1] << endl;

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

#if VTK_MAJOR_VERSION == 5
    writerUGridVTK->SetInput(uGridVTK2);
#else
    writerUGridVTK->SetInputData(uGridVTK2);
#endif

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
    double  vec[3];
    triaNormal->SetNumberOfComponents(3);

    for(ii=0; ii<GeomData.NodePosOrig.size(); ii++)
    {
      pt1 = GeomData.NodePosOrig[ii];
      //cout << pt1[0] << '\t' << pt1[1] << '\t' << pt1[2] << endl;
      pts[0] = pointsVTK->InsertNextPoint(pt1[0], pt1[1], pt1[2]);
    }

    treeCGAL.clear();
    trianglesCGAL.clear();

    for(ii=0; ii<nImmInt; ii++)
    {
      //cout << ii << '\t' << ImmIntgElems[ii]->pointNums.size() << endl;

      if(ImmIntgElems[ii]->pointNums.size() == 3) // triangle
      {
        pt1 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[0]];
        pt2 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[1]];
        pt3 = GeomData.NodePosOrig[ImmIntgElems[ii]->pointNums[2]];

        poly = new myTria(pt1, pt2, pt3);

        ImmersedFaces.push_back(poly);

        CGAL_Point  pta(pt1[0], pt1[1], pt1[2]);
        CGAL_Point  ptb(pt2[0], pt2[1], pt2[2]);
        CGAL_Point  ptc(pt3[0], pt3[1], pt3[2]);

        trianglesCGAL.push_back(CGAL_Triangle(pta, ptb, ptc));
        //treeCGAL.insert(CGAL_Triangle(pta, ptb, ptc));

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

        CGAL_Point  pta(pt1[0], pt1[1], pt1[2]);
        CGAL_Point  ptb(pt2[0], pt2[1], pt2[2]);
        CGAL_Point  ptc(pt3[0], pt3[1], pt3[2]);
        CGAL_Point  ptd(pt4[0], pt4[1], pt4[2]);

        trianglesCGAL.push_back(CGAL_Triangle(pta, ptb, ptd));
        trianglesCGAL.push_back(CGAL_Triangle(pta, ptd, ptc));

        //treeCGAL.insert(CGAL_Triangle(pta, ptb, ptc));
        //treeCGAL.insert(CGAL_Triangle(pta, ptc, ptd));

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

    // constructs AABB tree
    treeCGAL.insert(trianglesCGAL.begin(), trianglesCGAL.end());
    treeCGAL.accelerate_distance_queries();

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

#if VTK_MAJOR_VERSION == 5
    writerPolyData->SetInput(polyDataVTK);
    //massProp->SetInput(polyDataVTK);
#else
    writerPolyData->SetInputData(polyDataVTK);
    //massProp->SetInputData(polyDataVTK);
#endif

    writerPolyData->Write();
  }

  return;
}


void ImmersedSolid::updateImmersedFaces()
{
  //if(totalDOF == 0)
    //return;

  vtkSmartPointer<vtkPoints>     pointsVTK   =  vtkSmartPointer<vtkPoints>::New();

  int ii=0, jj=0;

  myPoint  pt1, pt2, pt3, pt4, normal;
  vtkIdType  ptId;

  for(ii=0; ii<GeomData.NodePosCur.size(); ii++)
  {
    pt1 = GeomData.NodePosCur[ii];
    //cout << pt1[0] << '\t' << pt1[1] << '\t' << pt1[2] << endl;
    ptId = pointsVTK->InsertNextPoint(pt1[0], pt1[1], pt1[2]);
  }

  if(DIM == 3)
  {
    //cout << " bbbbbbbbbbbbbb " << endl;

    polyDataVTK->SetPoints(pointsVTK);

    //////////////////
    // VTK
    //////////////////
    //selectEnclosedPoints->Initialize(polyDataVTK);

    //////////////////
    // CGAL
    //////////////////
    
    treeCGAL.clear();
    trianglesCGAL.clear();
    
    for(ii=0; ii<nImmInt; ii++)
    {
      //cout << ii << '\t' << ImmIntgElems[ii]->pointNums.size() << endl;

      if(ImmIntgElems[ii]->pointNums.size() == 3) // triangle
      {
        pt1 = GeomData.NodePosCur[ImmIntgElems[ii]->pointNums[0]];
        pt2 = GeomData.NodePosCur[ImmIntgElems[ii]->pointNums[1]];
        pt3 = GeomData.NodePosCur[ImmIntgElems[ii]->pointNums[2]];

        CGAL_Point  pta(pt1[0], pt1[1], pt1[2]);
        CGAL_Point  ptb(pt2[0], pt2[1], pt2[2]);
        CGAL_Point  ptc(pt3[0], pt3[1], pt3[2]);

        trianglesCGAL.push_back(CGAL_Triangle(pta, ptb, ptc));
        //treeCGAL.insert(CGAL_Triangle(pta, ptb, ptc));
      }
      else if(ImmIntgElems[ii]->pointNums.size() == 4) // quad
      {
        pt1 = GeomData.NodePosCur[ImmIntgElems[ii]->pointNums[0]];
        pt2 = GeomData.NodePosCur[ImmIntgElems[ii]->pointNums[1]];
        pt3 = GeomData.NodePosCur[ImmIntgElems[ii]->pointNums[2]];
        pt4 = GeomData.NodePosCur[ImmIntgElems[ii]->pointNums[3]];

        CGAL_Point  pta(pt1[0], pt1[1], pt1[2]);
        CGAL_Point  ptb(pt2[0], pt2[1], pt2[2]);
        CGAL_Point  ptc(pt3[0], pt3[1], pt3[2]);
        CGAL_Point  ptd(pt4[0], pt4[1], pt4[2]);

        trianglesCGAL.push_back(CGAL_Triangle(pta, ptb, ptd));
        trianglesCGAL.push_back(CGAL_Triangle(pta, ptd, ptc));
      }
      else
      {
        cout << " Error in ImmersedSolid::setImmersedFaces()  " << endl;
      }
    }

    treeCGAL.insert(trianglesCGAL.begin(), trianglesCGAL.end());
    treeCGAL.accelerate_distance_queries();
  } // if(DIM == 3)

  //cout << " treeCGAL size = " << treeCGAL.size() << endl;

  for(ii=0; ii<nImmInt; ii++)
  {
    for(jj=0; jj<ImmIntgElems[ii]->pointNums.size(); jj++)
      ImmersedFaces[ii]->updatePoint(jj, GeomData.NodePosCur[ImmIntgElems[ii]->pointNums[jj]] );

    ImmersedFaces[ii]->computeNormal();
    ImmersedFaces[ii]->computeAABB();
  }

  return;
}



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
    //return  selectEnclosedPoints->IsInsideSurface( ptTemp[0], ptTemp[1], ptTemp[2] );

    //CGAL_Point  ray_begin(ptTemp[0], ptTemp[1], ptTemp[2]);
    //CGAL_Point  ray_end(1000.0, ptTemp[1], ptTemp[2]);
    // counts #intersections
    //CGAL_Ray ray_query(ray_begin, ray_end);

    c = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                 CGAL_Point(ptTemp[0], ptTemp[1], ptTemp[2]),
                 CGAL_Point(1000.0,    ptTemp[1], ptTemp[2]))) % 2;
    //std::cout << nints << " intersections with ray query" << std::endl;
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
  //bbox.printSelf();
  
  if( !(bbox.doIntersect(bbTemp)) )
    return 0;

    //cout << " ImmersedFaces.size() = " << ImmersedFaces.size() << endl;

    int  val=0, ii=0, jj=0, c=0, pp=0;
    myPoint  ptTemp;
    vector<myPoint>  ptVec(4);
    
    //vector<int>  vecCorners(4);
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

          //cout << " ii = " << ii << endl;

          if( ImmersedFaces[ii]->doIntersect(bbTemp) )
          {
            //cout << " face " << ii << " intersects with cell AABB " << endl;
            for(jj=0; jj<2; jj++)
            {
              ptTemp = ImmersedFaces[ii]->GetPoint(jj);
              if( bbTemp.within(ptTemp) )
                ptOut.push_back(ptTemp);
            }

            for(jj=0; jj<4; jj++)
            {
              vecEdges[jj] = ImmersedFaces[ii]->IntersectWithLine(cellLines[jj], ptOut);
            }
          }
        } // for(ii=0; ii<ImmersedFaces.size(); ii++)
        cellLines.clear();

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
    else
    {
      val = -1;
    } // else

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

    if(sideTemp == 0 || sideTemp == 1)
    {
      //vecTemp[0] = selectEnclosedPoints->IsInsideSurface( coord3, bbTemp.minBB[1], bbTemp.minBB[2] );
      //vecTemp[1] = selectEnclosedPoints->IsInsideSurface( coord3, bbTemp.maxBB[1], bbTemp.minBB[2] );
      //vecTemp[2] = selectEnclosedPoints->IsInsideSurface( coord3, bbTemp.minBB[1], bbTemp.maxBB[2] );
      //vecTemp[3] = selectEnclosedPoints->IsInsideSurface( coord3, bbTemp.maxBB[1], bbTemp.maxBB[2] );

      vecTemp[0] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(coord3, bbTemp.minBB[1], bbTemp.minBB[2]), 
                    CGAL_Point(1000.0, bbTemp.minBB[1], bbTemp.minBB[2]))) % 2;

      vecTemp[1] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(coord3, bbTemp.maxBB[1], bbTemp.minBB[2]), 
                    CGAL_Point(1000.0, bbTemp.maxBB[1], bbTemp.minBB[2]))) % 2;

      vecTemp[2] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(coord3, bbTemp.minBB[1], bbTemp.maxBB[2]), 
                    CGAL_Point(1000.0, bbTemp.minBB[1], bbTemp.maxBB[2]))) % 2;

      vecTemp[3] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(coord3, bbTemp.maxBB[1], bbTemp.maxBB[2]), 
                    CGAL_Point(1000.0, bbTemp.maxBB[1], bbTemp.maxBB[2]))) % 2;

    }
    else if(sideTemp == 2 || sideTemp == 3)
    {
      //vecTemp[0] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], coord3, bbTemp.minBB[2] );
      //vecTemp[1] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], coord3, bbTemp.minBB[2] );
      //vecTemp[2] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], coord3, bbTemp.maxBB[2] );
      //vecTemp[3] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], coord3, bbTemp.maxBB[2] );

      vecTemp[0] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.minBB[0], coord3, bbTemp.minBB[2]), 
                    CGAL_Point(1000.0,          coord3, bbTemp.minBB[2]))) % 2;

      vecTemp[1] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.maxBB[0], coord3, bbTemp.minBB[2]), 
                    CGAL_Point(1000.0,          coord3, bbTemp.minBB[2]))) % 2;

      vecTemp[2] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.minBB[0], coord3, bbTemp.maxBB[2]), 
                    CGAL_Point(1000.0,          coord3, bbTemp.maxBB[2]))) % 2;

      vecTemp[3] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.maxBB[0], coord3, bbTemp.maxBB[2]), 
                    CGAL_Point(1000.0,          coord3, bbTemp.maxBB[2]))) % 2;
    }
    else if(sideTemp == 4 || sideTemp == 5)
    {
      //vecTemp[0] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], bbTemp.minBB[1], coord3 );
      //vecTemp[1] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], bbTemp.minBB[1], coord3 );
      //vecTemp[2] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], bbTemp.maxBB[1], coord3 );
      //vecTemp[3] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], bbTemp.maxBB[1], coord3 );

      vecTemp[0] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.minBB[0], bbTemp.minBB[1], coord3), 
                    CGAL_Point(1000.0,          bbTemp.minBB[1], coord3))) % 2;

      vecTemp[1] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.maxBB[0], bbTemp.minBB[1], coord3), 
                    CGAL_Point(1000.0,          bbTemp.minBB[1], coord3))) % 2;

      vecTemp[2] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.minBB[0], bbTemp.maxBB[1], coord3), 
                    CGAL_Point(1000.0,          bbTemp.maxBB[1], coord3))) % 2;

      vecTemp[0] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.maxBB[0], bbTemp.maxBB[1], coord3), 
                    CGAL_Point(1000.0,          bbTemp.maxBB[1], coord3))) % 2;
    }
    //printVector(vecTemp);

    if( std::equal(vecTemp.begin()+1, vecTemp.end(), vecTemp.begin()) )
    {
      if( vecTemp[0] == 1 )
        val = id+1;
      else
        val = 0;
    }
    else
    {
      val = -1;
    } // else

  return val;
}




int  ImmersedSolid::doIntersect3D(AABB&  bbTemp, bool flag, vector<int>& vecTemp, vector<myPoint>&  ptOut)
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
  
    //bb2.printSelf();

    if( !(bbox.doIntersect(bbTemp)) )
      return 0;

    int val=0, ii=0, jj=0, c=0, pp=0;
    myPoint  ptTemp;
    vector<myPoint>  ptVec(8);

    //bbTemp.printSelf();
    
    //vecTemp[0] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.minBB[2] );
    //vecTemp[1] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.minBB[2] );
    //vecTemp[2] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.minBB[2] );
    //vecTemp[3] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.minBB[2] );

    //vecTemp[4] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.maxBB[2] );
    //vecTemp[5] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.maxBB[2] );
    //vecTemp[6] = selectEnclosedPoints->IsInsideSurface( bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2] );
    //vecTemp[7] = selectEnclosedPoints->IsInsideSurface( bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2] );

    vecTemp[0] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.minBB[2]), 
                    CGAL_Point(100.0,           bbTemp.minBB[1], bbTemp.minBB[2]))) % 2;

    vecTemp[1] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.minBB[2]), 
                    CGAL_Point(100.0,           bbTemp.minBB[1], bbTemp.minBB[2]))) % 2;

    vecTemp[2] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]), 
                    CGAL_Point(100.0,           bbTemp.maxBB[1], bbTemp.minBB[2]))) % 2;

    vecTemp[3] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]), 
                    CGAL_Point(100.0,           bbTemp.maxBB[1], bbTemp.minBB[2]))) % 2;

    vecTemp[4] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]), 
                    CGAL_Point(100.0,           bbTemp.minBB[1], bbTemp.maxBB[2]))) % 2;

    vecTemp[5] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]), 
                    CGAL_Point(100.0,           bbTemp.minBB[1], bbTemp.maxBB[2]))) % 2;

    vecTemp[6] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]), 
                    CGAL_Point(100.0,           bbTemp.maxBB[1], bbTemp.maxBB[2]))) % 2;

    vecTemp[7] = treeCGAL.number_of_intersected_primitives(CGAL_Ray(
                    CGAL_Point(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]), 
                    CGAL_Point(100.0,           bbTemp.maxBB[1], bbTemp.maxBB[2]))) % 2;


    //printVector(vecTemp);

    if( std::equal(vecTemp.begin()+1, vecTemp.end(), vecTemp.begin()) )
    {
      if( vecTemp[0] == 1 )
        val = id+1;
      else
        val = 0;
    }
    else
    {
      val = -1;
      
      /*
      if(flag)
      {
        ///////////////////////////////////
        // find the intersection points between the current background grid cell 
        // and the immersed polyhedron

        //int verts[12][3] ={ {0,1,2}, {2,1,3},{2,3,7}, {7,3,1}, {7,1,5}, {5,1,4}, {5,4,7}, {7,4,6}, {7,6,2}, {2,6,4}, {2,4,0}, {0,4,1} };
        int verts[12][3] ={ {0,2,1}, {2,3,1},{2,7,3}, {7,1,3}, {7,5,1}, {5,4,1}, {5,7,4}, {7,6,4}, {7,2,6}, {2,4,6}, {2,0,4}, {0,1,4} };

        vector<myTria>  cellTrias;
        
        for(ii=0; ii<8; ii++)
        {
          if(vecTemp[ii] == 0)
            ptOut.push_back(ptVec[ii]);
        }

        for(ii=0; ii<12; ii++)
        {
          cellTrias.push_back(myTria(ptVec[verts[ii][0]], ptVec[verts[ii][1]], ptVec[verts[ii][2]]));
          cellTrias[ii].computeNormal();
          cellTrias[ii].computeAABB();
        }
        //cout << " cell triangles constructed " << endl;

        for(ii=0; ii<ImmersedFaces.size(); ii++)
        {
          if( ImmersedFaces[ii]->doIntersect(bbTemp) )
          {
            //cout << " face " << ii << " intersects with cell AABB " << endl;
            for(jj=0; jj<3; jj++)
            {
              ptTemp = ImmersedFaces[ii]->GetPoint(jj);
              if( bbTemp.within(ptTemp) )
                ptOut.push_back(ptTemp);
            }

            for(jj=0; jj<12; jj++)
            {
              pp = ImmersedFaces[ii]->IntersectWithTriangle(cellTrias[jj], ptOut);
            }
          }
        } // for(ii=0; ii<ImmersedFaces.size(); ii++)
        cellTrias.clear();
      } // if(flag)
      */
    } // else

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



