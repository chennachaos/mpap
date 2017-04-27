
#include "AdaptiveOctree.h"
#include "util.h"


/////////////////////////////////////////////////////////////////
//
// compute Quadrature points for cut cells using adaptive integration
//
// rule #1: quadrature points are computed only for those cells which either lie
//          completely inside fluid domain or are cut by the interface
//          
//  IF the current cell is a LEAF 
//     THEN compute the quadrature points of that cell
//  ELSE
//     traverse each CHILD of the cell and compute quadrature points of each CHILD



template<>
void AdaptiveOctree<2>::computeGaussPoints(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
  if( IsLeaf() )
  {
    myPoint  param, geom;
    int  gp, ii;

    //vector<myPoint>  ptOut;
    //vector<int>  vecTemp(4);

    //if(level == MAX_LEVEL)
      //GeomData->doIntersect2D(bbox, false, vecTemp, ptOut, domNums) ;


      // the cell is not a cut cell
      if( (domNums.size() == 1) && (domNums[0] == inclDom) )
      {
        for(gp=0; gp<GeomData->gausspoints.size(); gp++)
        {
          // parametric domain to integration master-quadrilateral domain

          for(ii=0; ii<2; ii++)
            param[ii] = 0.5*(knots[ii][2] * GeomData->gausspoints[gp][ii] + knots[ii][3]);

          GeomData->ComputeCoord(param, geom);

          //if( GeomData->within(geom) == inclDom )
          //{
            quadTemp.gausspoints.push_back(param);
            quadTemp.gaussweights.push_back(GeomData->gaussweights[gp]*JacMultElem);
          //}
        }
      }
      else // the cell is a cut cell
      {
          // merge Gauss points from all the levels below LEVEL #

          myPoint  ptTemp;
          int count=0;

          // find the Gauss points of the current cell that lie inside the fluid domain (domain '0')
          //
          for(gp=0; gp<GeomData->gausspoints.size(); gp++)
          {
            // parametric domain to integration master-quadrilateral domain

            for(ii=0; ii<2; ii++)
              param[ii] = 0.5*(knots[ii][2] * GeomData->gausspoints[gp][ii] + knots[ii][3]);

            //cout << gp << '\t' << GeomData->gausspoints[gp][0] << '\t' << GeomData->gausspoints[gp][1] << endl;
            //cout << gp << '\t' << param(0) << '\t' << param(1) << endl;

            GeomData->ComputeCoord(param, geom);

            if( GeomData->within(geom) == inclDom )
            {
              quadTemp.gausspoints.push_back(param);
              quadTemp.gaussweights.push_back(GeomData->gaussweights[gp]*JacMultElem);
              count++;
            }
          } // for(gp=0; gp<GeomData->gausspoints.size(); gp++)
          
          //cout << " fact " << endl;
          if(mergeFlag && refLev2)
          {
            //cout << wt << '\t' << mergeFlag << '\t' << refLev2 << '\t' << vecPts.size() << endl;

            double  wt=0.0, fact=0.0;
            int  nGPsMerge;

            this->mergeGaussPoints(refLev2, inclDom, 1, nGPsMerge, ptTemp, wt);

            // adjust weights from the lower levels to the weights of the GPs from the current cell

            if( count == 0 )
            {
              if( nGPsMerge != 0 )
              {
                quadTemp.gausspoints.push_back(ptTemp);
                quadTemp.gaussweights.push_back(wt);
              }
            }
            else
            {
              if( nGPsMerge != 0 )
              {
                //cout << " wt = " << wt << '\t' << count << '\t' << fact << " \n\n " << endl;
                fact = wt/count;

                ii = quadTemp.gaussweights.size()-1;

                for(gp=0; gp<count; gp++)
                  quadTemp.gaussweights[ii-gp] = fact;
              }
            }
          } // if(mergeFlag && refLev2)
      }
    //}
  }
  else
  {
    for(int gp=0; gp<NUM_CHILDREN; gp++)
    {
      child[gp]->computeGaussPoints(refLev2, inclDom, chkFlag, mergeFlag, quadTemp);
    }
  }

  return;
}



template<>
void AdaptiveOctree<3>::computeGaussPoints2Dfor3D(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
cout << " IsLeaf() ... 3D " << IsLeaf() << endl;

  return;
}



template<>
void AdaptiveOctree<2>::computeGaussPoints2Dfor3D(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
  //cout << " IsLeaf() " << IsLeaf() << endl;
  if( IsLeaf() )
  {
    myPoint  param, geom;
    int  gp, ii;

    //vector<myPoint>  ptOut;
    //vector<int>  vecTemp(4);

    //if(level == MAX_LEVEL)
      //GeomData->doIntersect2D(bbox, false, vecTemp, ptOut, domNums) ;

     //cout << " AdaptiveOctree<2>::computeGaussPoints2Dfor3D " << endl;

      if( (domNums.size() == 1) && (domNums[0] == inclDom) )
      {
        for(gp=0; gp<GeomData->gausspoints2D.size(); gp++)
        {
          // parametric domain to integration master-quadrilateral domain

          for(ii=0; ii<2; ii++)
            param[ii] = 0.5*(knots[ii][2] * GeomData->gausspoints2D[gp][ii] + knots[ii][3]);
          
          map2DPointTo3DPoint(sideTemp, param, param3);

           //GeomData->ComputeCoord(param, geom);
           //cout << " Geom " << geom[0] << '\t' << geom[1] << '\t' << geom[2] << endl;

            //if( GeomData->within(geom) == inclDom )
            //{
              quadTemp.gausspoints.push_back(param);
              quadTemp.gaussweights.push_back(GeomData->gaussweights2D[gp]*JacMultElem);
            //}
        }
      }
      else 
      {
        // merge Gauss points from all the levels below LEVEL 5
        
          myPoint  ptTemp;
          int count=0;

          // find the Gauss points of the current cell that lie inside the fluid domain (domain '0')
          //
          for(gp=0; gp<GeomData->gausspoints2D.size(); gp++)
          {
            // parametric domain to integration master-quadrilateral domain

            for(ii=0; ii<2; ii++)
              param[ii] = 0.5*(knots[ii][2] * GeomData->gausspoints2D[gp][ii] + knots[ii][3]);

            //cout << gp << '\t' << GeomData->gausspoints[gp][0] << '\t' << GeomData->gausspoints[gp][1] << endl;
            //cout << gp << '\t' << param(0) << '\t' << param(1) << endl;

            map2DPointTo3DPoint(sideTemp, param, param3);

            GeomData->ComputeCoord(param, geom);

              if( GeomData->within(geom) == inclDom )
              {
                quadTemp.gausspoints.push_back(param);
                quadTemp.gaussweights.push_back(GeomData->gaussweights2D[gp]*JacMultElem);
                count++;
              }
          } // for(gp=0; gp<GeomData->gausspoints.size(); gp++)
          
          //cout << " fact " << endl;
          if(mergeFlag && refLev2)
          {
            //cout << wt << '\t' << mergeFlag << '\t' << refLev2 << '\t' << vecPts.size() << endl;

            double  wt=0.0, fact=0.0;
            int nGPsMerge;

            this->mergeGaussPoints(refLev2, inclDom, 1, nGPsMerge, ptTemp, wt);

            // adjust weights from the lower levels to the weights of the GPs from the current cell

            if( count == 0 )
            {
              if( nGPsMerge != 0 )
              {
                quadTemp.gausspoints.push_back(ptTemp);
                quadTemp.gaussweights.push_back(wt);
              }
            }
            else
            {
              if( nGPsMerge != 0 )
              {
                //cout << " wt = " << wt << '\t' << count << '\t' << fact << " \n\n " << endl;
                fact = wt/count;

                ii = quadTemp.gaussweights.size()-1;

                for(gp=0; gp<count; gp++)
                  quadTemp.gaussweights[ii-gp] = fact;
              }
            }
          } // if(mergeFlag && refLev2)
      }
    //}
  }
  else
  {
    for(int gp=0; gp<NUM_CHILDREN; gp++)
    {
      child[gp]->computeGaussPoints2Dfor3D(refLev2, inclDom, chkFlag, mergeFlag, quadTemp);
    }
  }

  return;
}



template<>
void AdaptiveOctree<2>::computeGaussPointsForMerging(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
  if( IsLeaf() )
  {
    myPoint  param, geom;
    int  gp, ii;
    double  wt=0.0, fact=0.0;

        //
        for(gp=0; gp<GeomData->gausspoints.size(); gp++)
        {
          // parametric domain to integration master-quadrilateral domain

          for(ii=0; ii<2; ii++)
            param[ii] = 0.5*(knots[ii][2] * GeomData->gausspoints[gp][ii] + knots[ii][3]);

          GeomData->ComputeCoord(param, geom);

          if( GeomData->within(geom) == inclDom )
          {
            quadTemp.gausspoints.push_back(param);
            quadTemp.gaussweights.push_back(GeomData->gaussweights[gp]*JacMultElem);
          }
        }
        //

        /*
          // parametric domain to integration master-quadrilateral domain

          for(ii=0; ii<2; ii++)
            param[ii] = 0.5*(knots[ii][2] * 0.0 + knots[ii][3]);

          GeomData->ComputeCoord(param, geom);

          if( GeomData->within(geom) == inclDom )
          {
            quadTemp.gausspoints.push_back(param);
            quadTemp.gaussweights.push_back(GeomData->gaussweights[gp]*JacMultElem);
          }
        */
  }
  else
  {
    for(int gp=0; gp<NUM_CHILDREN; gp++)
    {
      child[gp]->computeGaussPointsForMerging(refLev2, inclDom, chkFlag, mergeFlag, quadTemp);
    }
  }

  return;
}




template<>
void AdaptiveOctree<2>::computeGaussPointsForMerging2Dfor3D(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
  if( IsLeaf() )
  {
    myPoint  param, geom;
    int  gp, ii;
    double  wt=0.0, fact=0.0;

        //
        for(gp=0; gp<GeomData->gausspoints2D.size(); gp++)
        {
          // parametric domain to integration master-quadrilateral domain

          for(ii=0; ii<2; ii++)
            param[ii] = 0.5*(knots[ii][2] * GeomData->gausspoints2D[gp][ii] + knots[ii][3]);

          map2DPointTo3DPoint(sideTemp, param, param3);

          GeomData->ComputeCoord(param, geom);

          if( GeomData->within(geom) == inclDom )
          {
            quadTemp.gausspoints.push_back(param);
            quadTemp.gaussweights.push_back(GeomData->gaussweights2D[gp]*JacMultElem);
          }
        }
        //

        /*
          // parametric domain to integration master-quadrilateral domain

          for(ii=0; ii<2; ii++)
            param[ii] = 0.5*(knots[ii][2] * 0.0 + knots[ii][3]);

          map2DPointTo3DPoint(sideTemp, param, param3);

          GeomData->ComputeCoord(param, geom);

          if( GeomData->within(geom) == inclDom )
          {
            quadTemp.gausspoints.push_back(param);
            quadTemp.gaussweights.push_back(GeomData->gaussweights[gp]*JacMultElem);
          }
        */
  }
  else
  {
    for(int gp=0; gp<NUM_CHILDREN; gp++)
    {
      child[gp]->computeGaussPointsForMerging2Dfor3D(refLev2, inclDom, chkFlag, mergeFlag, quadTemp);
    }
  }

  return;
}


template<>
void AdaptiveOctree<3>::computeGaussPoints(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
  if( IsLeaf() )
  {
    myPoint  param, geom;
    int  gp, ii;
    //vector<myPoint>  ptOut;
    //vector<int>  vecTemp(8);

    //if(level == MAX_LEVEL)
      //GeomData->doIntersect2D(bbox, false, vecTemp, ptOut, domNums) ;

      //cout << " cut cell " << endl;
      //printVector(domNums);

      if( (domNums.size() == 1) && (domNums[0] == inclDom) )
      {
        for(gp=0; gp<GeomData->gausspoints.size(); gp++)
        {
          // parametric domain to integration master-quadrilateral domain

          for(ii=0; ii<3; ii++)
            param[ii] = 0.5*(knots[ii][2] * GeomData->gausspoints[gp][ii] + knots[ii][3]);

          GeomData->ComputeCoord(param, geom);

            //if( GeomData->within(geom) == inclDom )
            //{
              quadTemp.gausspoints.push_back(param);
              quadTemp.gaussweights.push_back(GeomData->gaussweights[gp]*JacMultElem);
            //}
        }
      }
      else 
      {
        //cout << " cut cell " << endl;
          // merge Gauss points from all the levels below refLev1
        
          int count=0;

          // find the Gauss points of the current cell that lie inside the domain inclDom
          //
          for(gp=0; gp<GeomData->gausspoints.size(); gp++)
          {
            // parametric domain to integration master-quadrilateral domain

            for(ii=0; ii<3; ii++)
              param[ii] = 0.5*(knots[ii][2] * GeomData->gausspoints[gp][ii] + knots[ii][3]);

            //cout << gp << '\t' << GeomData->gausspoints[gp][0] << '\t' << GeomData->gausspoints[gp][1] << endl;
            //cout << gp << '\t' << param(0) << '\t' << param(1) << endl;

            GeomData->ComputeCoord(param, geom);

              if( GeomData->within(geom) == inclDom )
              {
                quadTemp.gausspoints.push_back(param);

                //vecWts.push_back(GeomData->gaussweights[gp]*JacMultElem);
                quadTemp.gaussweights.push_back(GeomData->gaussweights[gp]*JacMultElem);
                count++;
              }
          } // for(gp=0; gp<GeomData->gausspoints.size(); gp++)

          //cout << " fact " << endl;
          if(mergeFlag && refLev2)
          {
            //cout << wt << '\t' << mergeFlag << '\t' << refLev2 << '\t' << vecPts.size() << endl;

            double  wt=0.0, fact=0.0;
            myPoint  ptTemp;
            int nGPsMerge;

            this->mergeGaussPoints(refLev2, inclDom, 1, nGPsMerge, ptTemp, wt);

            // adjust weights from the lower levels to the weights of the GPs from the current cell

            if( count == 0 )
            {
              if( nGPsMerge != 0 )
              {
                quadTemp.gausspoints.push_back(ptTemp);
                quadTemp.gaussweights.push_back(wt);
              }
            }
            else
            {
              if( nGPsMerge != 0 )
              {
                //cout << " wt = " << wt << '\t' << count << '\t' << fact << " \n\n " << endl;
                fact = wt/count;

                ii = quadTemp.gaussweights.size()-1;

                for(gp=0; gp<count; gp++)
                  quadTemp.gaussweights[ii-gp] = fact;
              }
            }
          } // if(mergeFlag && refLev2)
      }
    //}
  }
  else
  {
    for(int gp=0; gp<NUM_CHILDREN; gp++)
    {
      child[gp]->computeGaussPoints(refLev2, inclDom, chkFlag, mergeFlag, quadTemp);
    }
  }

  return;
}



template<>
void AdaptiveOctree<3>::computeGaussPointsForMerging(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
  if( IsLeaf() )
  {
    myPoint  param, geom;
    int  gp, ii;

        //
        for(gp=0; gp<GeomData->gausspoints.size(); gp++)
        {
          for(ii=0; ii<3; ii++)
            param[ii] = 0.5*(knots[ii][2] * GeomData->gausspoints[gp][ii] + knots[ii][3]);

          GeomData->ComputeCoord(param, geom);

          if( GeomData->within(geom) == inclDom )
          {
            quadTemp.gausspoints.push_back(param);
            quadTemp.gaussweights.push_back(GeomData->gaussweights[gp]*JacMultElem);
          }
        }
        //

        /*
          for(ii=0; ii<3; ii++)
            param[ii] = 0.5*(knots[ii][2] * 0.0 + knots[ii][3]);

          GeomData->ComputeCoord(param, geom);

          if( GeomData->within(geom) == inclDom )
          {
            quadTemp.gausspoints.push_back(param);
            quadTemp.gaussweights.push_back(8.0*JacMultElem);
          }
        */
  }
  else
  {
    for(int gp=0; gp<NUM_CHILDREN; gp++)
    {
      child[gp]->computeGaussPointsForMerging(refLev2, inclDom, chkFlag, mergeFlag, quadTemp);
    }
  }

  return;
}





