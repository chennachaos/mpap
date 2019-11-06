
#include "TreeNode.h"
#include "headersVTK.h"
#include "DistFunctions.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "myTria.h"



template<>
int TreeNode<1>::prepareCutCell(vector<double>& cutFEMparams)
{
  // clear the previous subtriangulation/adaptive integration
  clearSubtriangulation();


  double  x0 = GeomData->computeCoord(0, knotBegin[0]);
  double  x1 = GeomData->computeCoord(0, knotEnd[0]);
  double  xif=0.5;

  bool f0 = ( x0 <= xif );
  bool f1 = ( x1 <= xif );

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

  return 1;
}





template<>
int TreeNode<2>::prepareCutCell(vector<double>& cutFEMparams)
{
  // clear the previous subtriangulation/adaptive integration
  clearSubtriangulation();


  int  ee=0, bb=0, ii=0, jj=0, domTemp=0;
  int  id, kk, dd;
  vtkIdType pt[20];
  myPoint  ptNew;
  vector<int>  cornerInOut(4), vecTemp;
  vector<myPoint>  ptVec, ptOut;


  if( (int) cutFEMparams[0] == 2 ) // adaptive integration
  {
    GeomData->doIntersect2D(bbox, false, cornerInOut, ptOut, domNums) ;
  }
  else  // subtriangulation
  {
    GeomData->doIntersect2D(bbox, true, cornerInOut, ptOut, domNums) ;

    if(domNums.size() > 1)
    {
      if( ptOut.size() == 0 )
      {
        cerr << " TreeNode<2>::prepareCutCell() ... Subtriangulation error! " << endl;
        exit(1);
        return -1;
      }

      // find all the intersection points

      vtkSmartPointer<vtkPoints>       pointsLoc    =  vtkSmartPointer<vtkPoints>::New();
      vtkSmartPointer<vtkPolyData>     polyDataLoc  =  vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkPolyData>     polyDataLoc2 =  vtkSmartPointer<vtkPolyData>::New();
      vtkSmartPointer<vtkDelaunay2D>   delaunay     =  vtkSmartPointer<vtkDelaunay2D>::New();

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

        vecTemp.push_back(dd);

        poly->SetDomainNumber(dd);

        subTrias.push_back(poly);
      } //  for(cellId2=0; cellId2<polyDataLoc2->GetNumberOfCells(); cellId2++)

      // if all the subtriangles are inside the same solid, 
      // then there is no need for subtriangulation
      // Hence, the subtriangulation is deleted
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

  //ptOut.clear();

  return 1;
}



template<>
int TreeNode<3>::prepareCutCell(vector<double>& cutFEMparams)
{
  // clear the previous subtriangulation/adaptive integration
  clearSubtriangulation();


  vector<myPoint>  ptOut;
  vector<int>  cornerInOut(8);

  if( (int) cutFEMparams[0] == 2 ) // adaptive integration
  {
    return GeomData->doIntersect3D(bbox, false, cornerInOut, ptOut, domNums) ;
  }
  else // subtriangulation
  {
    cerr << " TreeNode<3>::prepareCutCell() ... Subtriangulation is not implemented for 3D problems ... " << endl;
  }

  return 0;
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

  int  ii=0, nc=0, gp=0, domTemp=0, aa=0, side=0;
  double  area11=0.0, area12=0.0, area21=0.0, area22=0.0, area1=0.0, area2=0.0, temp=0.0;
  myPoint  param, geom, ptTemp;

  double  area = (bbox.maxBB[0]-bbox.minBB[0])*(bbox.maxBB[1]-bbox.minBB[1]);

  vector<myPoint>  gpsLoc;
  vector<double>  gwsLoc;

  Quadrature.reset();

  if(inclFlag)
  {
    //QuadratureDomNums.clear();

    for(vector<myPoly*>::iterator poly = subTrias.begin() ; poly != subTrias.end(); ++poly)
    {
      (*poly)->getGaussPointsCUTFEM(nGP, gpsLoc, gwsLoc );

      domTemp = (*poly)->getDomainNumber() ; 
      //area2 = (*poly)->volume();

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
            ptTemp[ii] = (2.0*param[ii] - knotSum[ii])/knotIncr[ii];
          }

          Quadrature.gausspoints.push_back(ptTemp);
          Quadrature.gaussweights.push_back(gwsLoc[gp]);

          //QuadratureDomNums.push_back( domTemp ) ;
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
        //cout << " area2 = " << area2 << endl;

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
            ptTemp[ii] = (2.0*param[ii] - knotSum[ii])/knotIncr[ii];
          }

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
    BoundaryQuadrature.resize(4);

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
  }

  return 1;
}



template<>
int TreeNode<3>::computeGaussPointsSubTrias(int nGP, int refLev2, int inclFlag, int flag2)
{
  cerr << " TreeNode<3>::computeGaussPointsSubTrias ... is not implemented " << endl;
  return 0;


  if(domNums.size() == 1)
    return 1;

  if( subTrias.size() == 0 )
  {
    cout << " Error in TreeNode<3>::computeGaussPointsSubTrias  " << endl;
    return -1;
  }

  //cout << " AAAAAAAAAA " << endl;

  int  ii=0, nc=0, gp=0, aa=0, side=0;
  double area11=0.0, area12=0.0, area21=0.0, area22=0.0, area1=0.0, area2=0.0, temp=0.0;
  myPoint  ptTemp, param;

  double  area = (bbox.maxBB[0]-bbox.minBB[0])*(bbox.maxBB[1]-bbox.minBB[1])*(bbox.maxBB[2]-bbox.minBB[2]);

  vector<myPoint>  gpsLoc;
  vector<double>  gwsLoc;

  Quadrature.reset();

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
            ptTemp[ii] = (2.0*param[ii] - knotSum[ii])/knotIncr[ii];
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

  // if the element is inside only one domain then 
  // there is no need to compute quadrature points using cutfem approach
  if(domNums.size() == 1)
    return 1;

  if(adapIntegNode != NULL)
    delete  adapIntegNode;

  adapIntegNode = NULL;

  //adapIntegNode = new AdaptiveBinarytree<2>(0);
  adapIntegNode = new AdaptiveOctree<2>(0);

  adapIntegNode->setKnots(knotBegin[0], knotEnd[0], knotBegin[1], knotEnd[1]);

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
  int  ii=0, gp=0, aa=0, side=0;
  for(gp=0; gp<Quadrature.gausspoints.size(); gp++)
  {
    for(ii=0; ii<2; ii++)
      Quadrature.gausspoints[gp][ii] = (2.0*Quadrature.gausspoints[gp][ii] - knotSum[ii])/knotIncr[ii];
  }

  myPoint  param, geom;
  if(inclFlag)
  {
    //QuadratureDomNums.clear();
    for(gp=0; gp<Quadrature.gausspoints.size(); gp++)
    {
      param[0]  = 0.5*(knotIncr[0] * Quadrature.gausspoints[gp][0] + knotSum[0]);
      param[1]  = 0.5*(knotIncr[1] * Quadrature.gausspoints[gp][1] + knotSum[1]);

      geom[0] = GeomData->computeCoord(0, param[0]);
      geom[1] = GeomData->computeCoord(1, param[1]);

      //QuadratureDomNums.push_back( GeomData->within(geom) ) ;
    }
  }

  if( isBoundary() )
  {
    BoundaryQuadrature.resize(4);

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
  
  adapIntegNode->setKnots(knotBegin[0], knotEnd[0], knotBegin[1], knotEnd[1], knotBegin[2], knotEnd[2]);

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
  int  ii=0, gp=0, aa=0, side=0;
  for(gp=0; gp<Quadrature.gausspoints.size(); gp++)
  {
    for(ii=0; ii<3; ii++)
      Quadrature.gausspoints[gp][ii] = (2.0*Quadrature.gausspoints[gp][ii] - knotSum[ii])/knotIncr[ii];
  }

  myPoint  param, geom;
  if(inclFlag)
  {
    //QuadratureDomNums.clear();
    for(gp=0; gp<Quadrature.gausspoints.size(); gp++)
    {
      param[0]  = 0.5*(knotIncr[0] * Quadrature.gausspoints[gp][0] + knotSum[0]);
      param[1]  = 0.5*(knotIncr[1] * Quadrature.gausspoints[gp][1] + knotSum[1]);
      param[2]  = 0.5*(knotIncr[2] * Quadrature.gausspoints[gp][2] + knotSum[2]);

      geom[0] = GeomData->computeCoord(0, param[0]);
      geom[1] = GeomData->computeCoord(1, param[1]);
      geom[2] = GeomData->computeCoord(1, param[2]);

      //QuadratureDomNums.push_back( GeomData->within(geom) ) ;
    }
  }


  if( isBoundary() )
  {
    BoundaryQuadrature.resize(6);

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

  uu[0] = knotBegin[0] + 0.25*( knotEnd[0] - knotBegin[0]);
  uu[1] = knotBegin[0] + 0.75*( knotEnd[0] - knotBegin[0]);

  vv[0] = knotBegin[1] + 0.25*( knotEnd[1] - knotBegin[1]);
  vv[1] = knotBegin[1] + 0.75*( knotEnd[1] - knotBegin[1]);

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


