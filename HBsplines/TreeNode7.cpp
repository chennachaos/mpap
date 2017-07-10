
#include "TreeNode.h"
#include "headersVTK.h"
#include "DistFunctions.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "myTria.h"



template<>
int TreeNode<1>::prepareCutCell(vector<double>& cutFEMparams)
{
  double  x0, x1, xif=0.5;
  
  x0 = GeomData->computeCoord(0, knots[0][0]);
  x1 = GeomData->computeCoord(0, knots[0][1]);
  
  cout << x0 << '\t' << x1 << endl;
  
  bool f0, f1;
  
  f0 = ( x0 <= xif );
  f1 = ( x1 <= xif );
  
  cout << f0 << '\t' << f1 << endl;
  
  if( f0 == f1 )
  {
    if(f0)
      domNums[0] = 0;
    else
      domNums[0] = 1;
  }
  else
  {
    domNums[0] = -1;
  }
  
  //cout << " domainNum " << domainNums << endl;

  return 1;
}





template<>
int TreeNode<2>::prepareCutCell(vector<double>& cutFEMparams)
{
  int  ee, bb, ii, jj, domTemp;

  vector<myPoint>  ptOut;
  vector<int>  cornerInOut(4);

  //cout << " cell id = " << id << endl;
  
  if( (int) cutFEMparams[0] == 2 ) // adaptive integration
  {
    GeomData->doIntersect2D(bbox, false, cornerInOut, ptOut, domNums) ;
    ptOut.clear();
  }
  else  // subtriangulation
  {
    GeomData->doIntersect2D(bbox, true, cornerInOut, ptOut, domNums) ;

    //printVector(domNums);

    if(domNums.size() > 1)
    {
      if( ptOut.size() == 0 )
        return -1;

      // remove the old subtriangulation

      clearSubtriangulation();

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
        ptNew[2] = 0.0;
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
      //delaunay->setTolerance(0.0010);
      //delaunay->SetAlpha(0.0);
      delaunay->SetOffset(1000.0);

#if VTK_MAJOR_VERSION == 5
      delaunay->SetInput(polyDataLoc);
#else
      delaunay->SetInputData(polyDataLoc);
#endif

      delaunay->Update();

      polyDataLoc2 = delaunay->GetOutput();

      vtkIdType  cellId2;
      vtkCell  *cellVTK2;
      myPoint  pt1, pt2, pt3, pt4;
      
      if( polyDataLoc2->GetNumberOfCells() == 0 )
        return -2;

      for(cellId2=0; cellId2<polyDataLoc2->GetNumberOfCells(); cellId2++)
      {
        cellVTK2 = polyDataLoc2->GetCell(cellId2);

        polyDataLoc2->GetPoint(cellVTK2->GetPointId(0), &pt1[0]);
        polyDataLoc2->GetPoint(cellVTK2->GetPointId(1), &pt2[0]);
        polyDataLoc2->GetPoint(cellVTK2->GetPointId(2), &pt3[0]);

        myPoly  *poly = new myTria(pt1, pt2, pt3);

        poly->centroid(pt4);

        dd = GeomData->within(pt4) ;

        //cout << cellId2 << '\t' << dd << endl;

        vecTemp.push_back(dd);
        
        poly->SetDomainNumber(dd);

        subTrias.push_back(poly);
      } //  for(cellId2=0; cellId2<polyDataLoc2->GetNumberOfCells(); cellId2++)
      
      if( std::equal(vecTemp.begin()+1, vecTemp.end(), vecTemp.begin()) )
      {
        domNums.clear();
        domNums.push_back(vecTemp[0]);

        clearSubtriangulation();
      }
      
      // once the subtriangulation is generated 
      // find the Gausspoints for the triangles which lie 
      // in the domain #0 (fluid domain)
    }
  } //else

  return 1;
}



template<>
int TreeNode<3>::prepareCutCell(vector<double>& cutFEMparams)
{
  int  ee, bb, ii, jj, domTemp;

  vector<myPoint>  ptOut;
  vector<int>  cornerInOut(8);
  
  if( (int) cutFEMparams[0] == 2 ) // adaptive integration
  {
    return GeomData->doIntersect3D(bbox, false, cornerInOut, ptOut, domNums) ;

    //if( isRightBoundary() && isTopBoundary() && isBackBoundary() )
    //{
      //bbox.printSelf();
      //printVector(cornerInOut);
      //printVector(domNums);
    //}
  }
  else // subtriangulation
  {
    cerr << " TreeNode<3>::prepareCutCell() ... Subtriangulation is not implemented for 3D problems ... " << endl;
    //GeomData->doIntersect3D(bbox, true, cornerInOut, ptOut, domNums) ;
  }

  return 1;
}


template<>
int TreeNode<1>::computeGaussPointsSubTrias(int nGP, int refLev2, int inclFlag, int flag2)
{
  return 1;
}





template<>
int TreeNode<2>::computeGaussPointsSubTrias(int nGP, int inclFlag, int flag1, int flag2)
{
  if(domNums.size() == 1)
    return 1;

  if( subTrias.size() == 0 )
  {
    cout << " Error in TreeNode<2>::computeGaussPointsSubTrias  " << endl;
    return -1;
  }

  int  ii, nc, gp, domTemp;
  double area11, area12, area21, area22, area, area1, area2, temp;
  area11 = area12 = area21 = area22 = 0.0;
  myPoint  param, geom, ptTemp;

  area = (bbox.maxBB[0]-bbox.minBB[0])*(bbox.maxBB[1]-bbox.minBB[1]);

  vector<myPoint>  gpsLoc;
  vector<double>  gwsLoc;

  Quadrature.reset();
  
  if(inclFlag)
  {
    QuadratureDomNums.clear();

    for(vector<myPoly*>::iterator poly = subTrias.begin() ; poly != subTrias.end(); ++poly)
    {
      (*poly)->getGaussPointsCUTFEM(nGP, gpsLoc, gwsLoc );
      
      domTemp = (*poly)->getDomainNumber() ; 
      //area2 = (*poly)->volume();
      
      //cout << area << '\t' << area2 << '\t' << area/area2 << endl;

        //if( area2 > 0.001*area )
        //{
        for(gp=0; gp<gpsLoc.size(); gp++)
        {
          for(ii=0; ii<2; ii++)
          {
            // point coordinates in physcial domain
            ptTemp[ii] = gpsLoc[gp][ii];

            // physcial domain to parametric domain
            param[ii] = GeomData->computeParam(ii, ptTemp[ii]);

            // parametric domain to integration master-quadrilateral domain
            ptTemp[ii] = (2.0*param[ii] - knots[ii][3])/knots[ii][2];
          }

          Quadrature.gausspoints.push_back(ptTemp);
          Quadrature.gaussweights.push_back(gwsLoc[gp]);

          //param[0]  = 0.5*(knots[0][2] * Quadrature.gausspoints[gp][0] + knots[0][3]);
          //param[1]  = 0.5*(knots[1][2] * Quadrature.gausspoints[gp][1] + knots[1][3]);

          //geom[0] = GeomData->computeCoord(0, param[0]);
          //geom[1] = GeomData->computeCoord(1, param[1]);

          //QuadratureDomNums.push_back( GeomData->within(geom) ) ;
          QuadratureDomNums.push_back( domTemp ) ;
        }
    }
  }
  else
  {
    for(vector<myPoly*>::iterator poly = subTrias.begin() ; poly != subTrias.end(); ++poly)
    {
      if( (*poly)->getDomainNumber() == 0 )
      {
        (*poly)->getGaussPointsCUTFEM(nGP, gpsLoc, gwsLoc );
      
        //area2 = (*poly)->volume();
      
        //cout << area << '\t' << area2 << '\t' << area/area2 << endl;

        //if( area2 > 0.001*area )
        //{
        for(gp=0; gp<gpsLoc.size(); gp++)
        {
          for(ii=0; ii<2; ii++)
          {
            // point coordinates in physcial domain
            geom[ii] = gpsLoc[gp][ii];

            // physcial domain to parametric domain
            param[ii] = GeomData->computeParam(ii, geom[ii]);

            // parametric domain to integration master-quadrilateral domain
            ptTemp[ii] = (2.0*param[ii] - knots[ii][3])/knots[ii][2];
          }

          //cout << " gp = " << gp << endl;

          Quadrature.gausspoints.push_back(ptTemp);
          Quadrature.gaussweights.push_back(gwsLoc[gp]);
        }
      } // if( (*poly)->getDomainNumber() == 0 )
    }
  }

  //temp = area11 + area12 + area21 + area22;
  //cout << " area1 = " << area11 << '\t' << area12 << endl;
  //cout << " area2 = " << area21 << '\t' << area22 << endl;
  //cout << " total area = " << area << '\t' << (area1+area2) << '\t' << temp << endl;

  if( isBoundary() )
  {
    //cout << " boundary element " << endl;

    BoundaryQuadrature.resize(4);
    int aa, side;

    for(side=0; side<4; side++)
      BoundaryQuadrature[side].reset();

    if(DirichletData.size() > 0)
    {
      for(aa=0;aa<DirichletData.size();aa++)
      {
        side  = (int) (DirichletData[aa][0] - 1);

        if( BoundaryQuadrature[side].gausspoints.size() == 0 )
          TreeNode<2>::computeGaussPointsAdapIntegrationBoundary(side, 0, 0, inclFlag, flag2);
      }
    }

    if(NeumannData.size() > 0)
    {
      for(aa=0;aa<NeumannData.size();aa++)
      {
        side  = (int) (NeumannData[aa][0] - 1);
        
        if( BoundaryQuadrature[side].gausspoints.size() == 0 )
          TreeNode<2>::computeGaussPointsAdapIntegrationBoundary(side, 0, 0, inclFlag, flag2);
      }
    }
    //cout << " boundary element " << endl;
  }

  return 1;
}



template<>
int TreeNode<3>::computeGaussPointsSubTrias(int nGP, int refLev2, int inclFlag, int flag2)
{
  if(domNums.size() == 1)
    return 1;
  
  if( subTrias.size() == 0 )
  {
    cout << " Error in TreeNode<3>::computeGaussPointsSubTrias  " << endl;
    return -1;
  }

  //cout << " AAAAAAAAAA " << endl;
  
  int  ii, nc, gp;
  double area11, area12, area21, area22, area, area1, area2, temp;
  myPoint  ptTemp, param;

  area = (bbox.maxBB[0]-bbox.minBB[0])*(bbox.maxBB[1]-bbox.minBB[1])*(bbox.maxBB[2]-bbox.minBB[2]);

  vector<myPoint>  gpsLoc;
  vector<double>  gwsLoc;

  //myPoly *poly;

  Quadrature.reset();

  area11 = area12 = area21 = area22 = 0.0;
  for(vector<myPoly*>::iterator poly = subTrias.begin() ; poly != subTrias.end(); ++poly)
  {
    if( (*poly)->getDomainNumber() == 0 )
    {
      (*poly)->getGaussPointsCUTFEM(nGP, gpsLoc, gwsLoc );
      
      area2 = (*poly)->volume();
      
      //cout << area << '\t' << area2 << '\t' << area/area2 << endl;

        //if( area2 > 0.001*area )
        //{
        for(gp=0; gp<gpsLoc.size(); gp++)
        {
          for(ii=0; ii<3; ii++)
          {
            // point coordinates in physcial domain
            ptTemp[ii] = gpsLoc[gp][ii];

            // physcial domain to parametric domain
            param[ii] = GeomData->computeParam(ii, ptTemp[ii]);

            // parametric domain to integration master-quadrilateral domain
            ptTemp[ii] = (2.0*param[ii] - knots[ii][3])/knots[ii][2];
          }

          Quadrature.gausspoints.push_back(ptTemp);
          Quadrature.gaussweights.push_back(gwsLoc[gp]);
        }
      //}
    }
  }


  if( isBoundary() )
  {
    BoundaryQuadrature.resize(6);
    int aa, side;

    for(side=0; side<6; side++)
      BoundaryQuadrature[side].reset();

    if(DirichletData.size() > 0)
    {
      for(aa=0;aa<DirichletData.size();aa++)
      {
        side  = (int) (DirichletData[aa][0] - 1);

        if( BoundaryQuadrature[side].gausspoints.size() == 0 )
          TreeNode<3>::computeGaussPointsAdapIntegrationBoundary(side, 0, 0, inclFlag, flag2);
      }
    }

    if(NeumannData.size() > 0)
    {
      for(aa=0;aa<NeumannData.size();aa++)
      {
        side  = (int) (NeumannData[aa][0] - 1);
        
        if( BoundaryQuadrature[side].gausspoints.size() == 0 )
          TreeNode<3>::computeGaussPointsAdapIntegrationBoundary(side, 0, 0, inclFlag, flag2);
      }
    }
  }


  return 1;
}






template<>
int TreeNode<1>::computeGaussPointsAdapIntegration(int refLev1, int refLev2, int inclFlag, int flag2)
{
  return 1;
}



template<>
int TreeNode<2>::computeGaussPointsAdapIntegration(int refLev1, int refLev2, int inclFlag, int flag2)
{
  // refLev1   --->  depth of subdivison for adaptive integration
  // refLev2   --->  depth of subdivison for Gauss point merging, refLev2 levels from refLev1
  // inclFlag  --->  flag to specify whether to include all the Gauss points or only those belonging 
  //                 to fluid domain, domain '0'

  //if(refLev1 == 0)
    //return;

  if(domNums.size() == 1)
    return 1;

  if(adapIntegNode != NULL)
    delete  adapIntegNode;

  adapIntegNode = NULL;

  //adapIntegNode = new AdaptiveBinarytree<2>(0);
  adapIntegNode = new AdaptiveOctree<2>(0);

  adapIntegNode->setKnots(knots[0][0], knots[0][1], knots[1][0], knots[1][1]);

  adapIntegNode->GeomData = GeomData;
  adapIntegNode->domNums = domNums;
  //adapIntegNode->setSplitDirection(-1);

  adapIntegNode->prepareData();
  adapIntegNode->subDivide(refLev1);

  Quadrature.reset();

  //for(int dd=0; dd<domNums.size(); dd++)
  //{
    //if( GeomData->domainInclYesNo[dd] )
      //adapIntegNode->computeGaussPoints(refLev2, inclFlag, 1, 1, Quadrature);
      adapIntegNode->computeGaussPoints(refLev2, 0, 1, 1, Quadrature);
  //}

  // parametric domain to integration master-quadrilateral domain
  int  ii, gp;
  for(gp=0; gp<Quadrature.gausspoints.size(); gp++)
  {
    for(ii=0; ii<2; ii++)
      Quadrature.gausspoints[gp][ii] = (2.0*Quadrature.gausspoints[gp][ii] - knots[ii][3])/knots[ii][2];
  }

  if(inclFlag)
  {
    myPoint  param, geom;
    QuadratureDomNums.clear();
    for(gp=0; gp<Quadrature.gausspoints.size(); gp++)
    {
      param[0]  = 0.5*(knots[0][2] * Quadrature.gausspoints[gp][0] + knots[0][3]);
      param[1]  = 0.5*(knots[1][2] * Quadrature.gausspoints[gp][1] + knots[1][3]);

      geom[0] = GeomData->computeCoord(0, param[0]);
      geom[1] = GeomData->computeCoord(1, param[1]);

      QuadratureDomNums.push_back( GeomData->within(geom) ) ;
    }
  }

  if( isBoundary() )
  {
    BoundaryQuadrature.resize(4);
    int aa, side;

    for(side=0; side<4; side++)
      BoundaryQuadrature[side].reset();

    if(DirichletData.size() > 0)
    {
      for(aa=0;aa<DirichletData.size();aa++)
      {
        side  = (int) (DirichletData[aa][0] - 1);

        if( BoundaryQuadrature[side].gausspoints.size() == 0 )
          TreeNode<2>::computeGaussPointsAdapIntegrationBoundary(side, refLev1, refLev2, inclFlag, flag2);
      }
    }

    if(NeumannData.size() > 0)
    {
      for(aa=0;aa<NeumannData.size();aa++)
      {
        side  = (int) (NeumannData[aa][0] - 1);
        
        if( BoundaryQuadrature[side].gausspoints.size() == 0 )
          TreeNode<2>::computeGaussPointsAdapIntegrationBoundary(side, refLev1, refLev2, inclFlag, flag2);
      }
    }
  }

  //cout << id << '\t' << Quadrature.gausspoints.size() << endl;

  return 1;
}
//



template<>
int TreeNode<3>::computeGaussPointsAdapIntegration(int refLev1, int refLev2, int inclFlag, int flag2)
{
  // refLev1   --->  depth of subdivison for adaptive integration
  // refLev2   --->  depth of subdivison for Gauss point merging, refLev2 levels from refLev1
  // inclFlag  --->  flag to specify whether to include all the Gauss points or only those belonging 
  //                 to fluid domain, domain '0'

  if(domNums.size() == 1)
    return 1;


  if(adapIntegNode != NULL)
    delete  adapIntegNode;

  adapIntegNode = NULL;

  //adapIntegNode = new AdaptiveBinarytree<3>(0);
  adapIntegNode = new AdaptiveOctree<3>(0);
  
  adapIntegNode->setKnots(knots[0][0], knots[0][1], knots[1][0], knots[1][1], knots[2][0], knots[2][1]);

  adapIntegNode->GeomData = GeomData;
  adapIntegNode->domNums = domNums;
  //adapIntegNode->setSplitDirection(-1);

  adapIntegNode->prepareData();
  adapIntegNode->subDivide(refLev1);
 
  Quadrature.reset();

  //for(int dd=0; dd<domNums.size(); dd++)
  //{
    //if( GeomData->domainInclYesNo[dd] )
      //adapIntegNode->computeGaussPoints(refLev2, inclFlag, 1, 1, Quadrature);
      adapIntegNode->computeGaussPoints(refLev2, 0, 1, 1, Quadrature);
  //}


  // parametric domain to integration master-quadrilateral domain
  int  ii, gp;
  for(gp=0; gp<Quadrature.gausspoints.size(); gp++)
  {
    for(ii=0; ii<3; ii++)
      Quadrature.gausspoints[gp][ii] = (2.0*Quadrature.gausspoints[gp][ii] - knots[ii][3])/knots[ii][2];
  }

  if(inclFlag)
  {
    myPoint  param, geom;
    QuadratureDomNums.clear();
    for(gp=0; gp<Quadrature.gausspoints.size(); gp++)
    {
      param[0]  = 0.5*(knots[0][2] * Quadrature.gausspoints[gp][0] + knots[0][3]);
      param[1]  = 0.5*(knots[1][2] * Quadrature.gausspoints[gp][1] + knots[1][3]);
      param[2]  = 0.5*(knots[2][2] * Quadrature.gausspoints[gp][2] + knots[2][3]);

      geom[0] = GeomData->computeCoord(0, param[0]);
      geom[1] = GeomData->computeCoord(1, param[1]);
      geom[2] = GeomData->computeCoord(1, param[2]);

      QuadratureDomNums.push_back( GeomData->within(geom) ) ;
    }
  }


  if( isBoundary() )
  {
    BoundaryQuadrature.resize(6);
    int aa, side;

    for(side=0; side<6; side++)
      BoundaryQuadrature[side].reset();

    if(DirichletData.size() > 0)
    {
      for(aa=0;aa<DirichletData.size();aa++)
      {
        side  = (int) (DirichletData[aa][0] - 1);

        //cout << " DirichletData ... " << aa << '\t' << side << endl;
        if( BoundaryQuadrature[side].gausspoints.size() == 0 )
          TreeNode<3>::computeGaussPointsAdapIntegrationBoundary(side, refLev1, refLev2, inclFlag, flag2);
        //printVector(domNums);
      }
    }

    if(NeumannData.size() > 0)
    {
      for(aa=0;aa<NeumannData.size();aa++)
      {
        side  = (int) (NeumannData[aa][0] - 1);

        //cout << " NeumannData ... " << aa << '\t' << side << endl;

        if( BoundaryQuadrature[side].gausspoints.size() == 0 )
          TreeNode<3>::computeGaussPointsAdapIntegrationBoundary(side, refLev1, refLev2, inclFlag, flag2);
      }
    }
  }

  //cout << id << '\t' << Quadrature.gausspoints.size() << endl;

  //cout << " AAAAAAAAAA " << endl;

  return 1;
}


template<>
int TreeNode<1>::checkCutCellValidityAdapIntegration()
{
  return 1;
}



template<>
int TreeNode<2>::checkCutCellValidityAdapIntegration()
{
  if(domNums.size() == 1)
    return 1;
  
  vector<int>  vectmp(4);
  double  vv[2], uu[2];
  
  uu[0] = knots[0][0] + 0.25*( knots[0][1] - knots[0][0]);
  uu[1] = knots[0][0] + 0.75*( knots[0][1] - knots[0][0]);

  vv[0] = knots[1][0] + 0.25*( knots[1][1] - knots[1][0]);
  vv[1] = knots[1][0] + 0.75*( knots[1][1] - knots[1][0]);

  //cout << " AAAAAAAAAA " << endl;

  //cout << " AAAAAAAAAA .... " << val << endl;

  myPoint  param, geom;
  
  param[1] = vv[0];  param[0] = uu[0];

  GeomData->computeCoord(param, geom);
  vectmp[0] = GeomData->within(geom) ;

  param[1] = vv[0];  param[0] = uu[1];

  GeomData->computeCoord(param, geom);
  vectmp[1] = GeomData->within(geom) ;

  param[1] = vv[1];  param[0] = uu[0];

  GeomData->computeCoord(param, geom);
  vectmp[2] = GeomData->within(geom) ;

  param[1] = vv[1];  param[0] = uu[1];

  GeomData->computeCoord(param, geom);
  vectmp[3] = GeomData->within(geom) ;

  if( std::equal(vectmp.begin()+1, vectmp.end(), vectmp.begin()) )
  {
    domNums.push_back(vectmp[0]);
    return 0;
  }

  return 1;
}





template<>
int TreeNode<3>::checkCutCellValidityAdapIntegration()
{
  return 1;
}


