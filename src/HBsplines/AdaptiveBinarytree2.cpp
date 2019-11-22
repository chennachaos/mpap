
#include "AdaptiveBinarytree.h"
#include "util.h"


/*
template<>
void AdaptiveBinarytree<2>::computeGaussPoints(int flag, GaussQuadrature&  quadTemp)
{
  if( isLeaf() )
  {
    myPoint  param, geom;
    int  gp, ii, bb;
    vector<myPoint>  ptOut;
    vector<int>  vecTemp(4);

    //cout << " domainNum = " << domainNum << endl;

    GeomData->doIntersect2D(bbox, false, vecTemp, ptOut, domNums) ;

    //if(domNums.size() > 1)
    //{
      //cout << " gps " << GeomData->gausspoints.size() << endl;
      //cout << knotBegin[0] << '\t' << knotEnd[0] << '\t' << knotIncr[0] << '\t' << knotSum[0] << endl;
      //cout << knotBegin[1] << '\t' << knotEnd[1] << '\t' << knotIncr[1] << '\t' << knotSum[1] << endl;

      if( (level < 5) || (level == MAX_LEVEL) )
      {
        for(gp=0; gp<GeomData->gausspoints.size(); gp++)
        {
          // parametric domain to integration master-quadrilateral domain

          for(ii=0; ii<2; ii++)
            param[ii] = 0.5*(knotIncr[ii] * GeomData->gausspoints[gp][ii] + knotSum[ii]);

          //cout << gp << '\t' << GeomData->gausspoints[gp][0] << '\t' << GeomData->gausspoints[gp][1] << endl;
          //cout << gp << '\t' << param(0) << '\t' << param(1) << endl;

          GeomData->computeCoord(param, geom);

          if(flag) // add only those Gauss points which are in domain '0' // for fluid problems
          {
            if( !(GeomData->within(geom)) )
            {
              quadTemp.gausspoints.push_back(param);
              quadTemp.gaussweights.push_back(GeomData->gaussweights[gp]*JacMultElem);
            }
          }
          else
          {
              quadTemp.gausspoints.push_back(param);
              quadTemp.gaussweights.push_back(GeomData->gaussweights[gp]*JacMultElem);
          }
        }
      }
      else 
      {
        // add only one Gauss point for the cells which are below level 5

          for(ii=0; ii<2; ii++)
            param[ii] = 0.5*(knotIncr[ii] * 0.0 + knotSum[ii]);

          GeomData->computeCoord(param, geom);

          if(flag)
          {
            if( !(GeomData->within(geom)) )
            {
              quadTemp.gausspoints.push_back(param);
              quadTemp.gaussweights.push_back(4.0*JacMultElem);
            }
          }
          else
          {
              quadTemp.gausspoints.push_back(param);
              quadTemp.gaussweights.push_back(4.0*JacMultElem);
          }
      }
    //}
  }
  else
  {
    for(int gp=0; gp<NUM_CHILDREN; gp++)
    {
      //if(child[gp] != NULL)
      //if( child[gp]->isLeaf() )
        child[gp]->computeGaussPoints(flag, quadTemp);
    }
  }

  return;
}
*/



template<>
void AdaptiveBinarytree<2>::computeGaussPoints(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
    if( isLeaf() )
    {
      myPoint  param, geom;
      int  gp=0, ii=0;

      //vector<myPoint>  ptOut;
      //vector<int>  vecTemp(4);

      //if(level == MAX_LEVEL)
        //GeomData->doIntersect2D(bbox, false, vecTemp, ptOut, domNums) ;

      if( (domNums.size() == 1) && (domNums[0] == inclDom) )
      {
        for(gp=0; gp<GeomData->gausspoints.size(); gp++)
        {
          // parametric domain to integration master-quadrilateral domain

          for(ii=0; ii<2; ii++)
            param[ii] = 0.5*(knotIncr[ii] * GeomData->gausspoints[gp][ii] + knotSum[ii]);

          GeomData->computeCoord(param, geom);

            //if( GeomData->within(geom) == inclDom )
            //{
              quadTemp.gausspoints.push_back(param);
              quadTemp.gaussweights.push_back(GeomData->gaussweights[gp]*JacMultElem);
            //}
        }
      }
      else 
      {
        // merge Gauss points from all the levels below 
        
          myPoint  ptTemp;
          int count=0;

          // find the Gauss points of the current cell that lie inside the fluid domain (domain '0')
          //
          for(gp=0; gp<GeomData->gausspoints.size(); gp++)
          {
            // parametric domain to integration master-quadrilateral domain

            for(ii=0; ii<2; ii++)
              param[ii] = 0.5*(knotIncr[ii] * GeomData->gausspoints[gp][ii] + knotSum[ii]);

            //cout << gp << '\t' << GeomData->gausspoints[gp][0] << '\t' << GeomData->gausspoints[gp][1] << endl;
            //cout << gp << '\t' << param(0) << '\t' << param(1) << endl;

            GeomData->computeCoord(param, geom);

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
void AdaptiveBinarytree<3>::computeGaussPoints2Dfor3D(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
    cout << " isLeaf() ... 3D " << isLeaf() << endl;
    return;
}



template<>
void AdaptiveBinarytree<2>::computeGaussPoints2Dfor3D(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
    if( isLeaf() )
    {
      myPoint  param, geom;
      int  gp=0, ii=0;

      //vector<myPoint>  ptOut;
      //vector<int>  vecTemp(4);

      //if(level == MAX_LEVEL)
        //GeomData->doIntersect2D(bbox, false, vecTemp, ptOut, domNums) ;

      //cout << " AdaptiveBinarytree<2>::computeGaussPoints2Dfor3D " << endl;

      if( (domNums.size() == 1) && (domNums[0] == inclDom) )
      {
        for(gp=0; gp<GeomData->gausspoints2D.size(); gp++)
        {
          // parametric domain to integration master-quadrilateral domain

          for(ii=0; ii<2; ii++)
            param[ii] = 0.5*(knotIncr[ii] * GeomData->gausspoints2D[gp][ii] + knotSum[ii]);
          
          map2DPointTo3DPoint(sideTemp, param, param3);

           //GeomData->computeCoord(param, geom);
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
              param[ii] = 0.5*(knotIncr[ii] * GeomData->gausspoints2D[gp][ii] + knotSum[ii]);

            //cout << gp << '\t' << GeomData->gausspoints[gp][0] << '\t' << GeomData->gausspoints[gp][1] << endl;
            //cout << gp << '\t' << param(0) << '\t' << param(1) << endl;

            map2DPointTo3DPoint(sideTemp, param, param3);

            GeomData->computeCoord(param, geom);

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
void AdaptiveBinarytree<2>::computeGaussPointsForMerging(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
    if( isLeaf() )
    {
        myPoint  param, geom;
        int  gp=0, ii=0;

        //
        for(gp=0; gp<GeomData->gausspoints.size(); gp++)
        {
          // parametric domain to integration master-quadrilateral domain

          for(ii=0; ii<2; ii++)
            param[ii] = 0.5*(knotIncr[ii] * GeomData->gausspoints[gp][ii] + knotSum[ii]);

          GeomData->computeCoord(param, geom);

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
            param[ii] = 0.5*(knotIncr[ii] * 0.0 + knotSum[ii]);

          GeomData->computeCoord(param, geom);

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
void AdaptiveBinarytree<2>::computeGaussPointsForMerging2Dfor3D(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
    if( isLeaf() )
    {
        myPoint  param, geom;
        int  gp=0, ii=0;
        double  wt=0.0, fact=0.0;

        //
        for(gp=0; gp<GeomData->gausspoints2D.size(); gp++)
        {
          // parametric domain to integration master-quadrilateral domain

          for(ii=0; ii<2; ii++)
            param[ii] = 0.5*(knotIncr[ii] * GeomData->gausspoints2D[gp][ii] + knotSum[ii]);

          map2DPointTo3DPoint(sideTemp, param, param3);

          GeomData->computeCoord(param, geom);

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
            param[ii] = 0.5*(knotIncr[ii] * 0.0 + knotSum[ii]);

          map2DPointTo3DPoint(sideTemp, param, param3);

          GeomData->computeCoord(param, geom);

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
void AdaptiveBinarytree<3>::computeGaussPoints(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
    if( isLeaf() )
    {
        myPoint  param, geom;
        int  gp=0, ii=0;

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
              param[ii] = 0.5*(knotIncr[ii] * GeomData->gausspoints[gp][ii] + knotSum[ii]);

            GeomData->computeCoord(param, geom);

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
              param[ii] = 0.5*(knotIncr[ii] * GeomData->gausspoints[gp][ii] + knotSum[ii]);

            //cout << gp << '\t' << GeomData->gausspoints[gp][0] << '\t' << GeomData->gausspoints[gp][1] << endl;
            //cout << gp << '\t' << param(0) << '\t' << param(1) << endl;

            GeomData->computeCoord(param, geom);

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
void AdaptiveBinarytree<3>::computeGaussPointsForMerging(int refLev2, int inclDom, int chkFlag, int mergeFlag, GaussQuadrature&  quadTemp)
{
    if( isLeaf() )
    {
        myPoint  param, geom;
        int  gp=0, ii=0;

        //
        for(gp=0; gp<GeomData->gausspoints.size(); gp++)
        {
          for(ii=0; ii<3; ii++)
            param[ii] = 0.5*(knotIncr[ii] * GeomData->gausspoints[gp][ii] + knotSum[ii]);

          GeomData->computeCoord(param, geom);

          if( GeomData->within(geom) == inclDom )
          {
            quadTemp.gausspoints.push_back(param);
            quadTemp.gaussweights.push_back(GeomData->gaussweights[gp]*JacMultElem);
          }
        }
        //

        /*
          for(ii=0; ii<3; ii++)
            param[ii] = 0.5*(knotIncr[ii] * 0.0 + knotSum[ii]);

          GeomData->computeCoord(param, geom);

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





