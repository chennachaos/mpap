
#include "AdaptiveBinarytree.h"


template<> int AdaptiveBinarytree<1>::nodecount = 0;
template<> int AdaptiveBinarytree<2>::nodecount = 0;
template<> int AdaptiveBinarytree<3>::nodecount = 0;

template<> int AdaptiveBinarytree<1>::MAX_LEVEL = 0;
template<> int AdaptiveBinarytree<2>::MAX_LEVEL = 0;
template<> int AdaptiveBinarytree<3>::MAX_LEVEL = 0;



template<>
void AdaptiveBinarytree<2>::subDivide(int refLev)
{
    //assert("Current Cell is a LEAF and Subdivision starts" && isLeaf());

    if( !isLeaf() || (refLev < 1) )
      return;


    vector<myPoint>  ptOut;
    vector<int>  cornerInOut(4);

    if(sideTemp == -1)
      GeomData->doIntersect2D(bbox, false, cornerInOut, ptOut, domNums) ;
    else
      GeomData->doIntersect2Dfor3D(sideTemp, coord3, bbox, false, cornerInOut, ptOut, domNums) ;

    if( domNums.size() == 1 )
      return;


    if(level == (2*refLev) )
      return;

    int ii=0, jj=0, bb=0, splitDirTemp=0;

    // find the split direction based on the Inside/Outside from of the current cell
    // and the split direction of the parent

    int  levTemp = level+1;

    if(splitDir == -1)
    {
      if( (cornerInOut[0] == cornerInOut[1]) && (cornerInOut[2] == cornerInOut[3]) )
        splitDirTemp = 1;
      else //if( (cornerInOut[0] == cornerInOut[2]) && (cornerInOut[1] == cornerInOut[3]) )
        splitDirTemp = 0;
    }
    else
    {
      if(splitDir == 0)
        splitDirTemp = 1;
      else
        splitDirTemp = 0;
    }
    
    NUM_CHILDREN = 2;
    child      = new AdaptiveBinarytree_PTR[NUM_CHILDREN];

    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii] = new AdaptiveBinarytree<2>(levTemp);

      child[ii]->setParent(this);
    }

    child[FIRST]->orientation = FIRST;
    child[SECOND]->orientation = SECOND;


    double  mid = 0.5*(knots[splitDirTemp][1] + knots[splitDirTemp][0]);

    if(splitDirTemp == 0)
    {
      child[FIRST]->setKnots(Dir1, knots[0][0], mid);
      child[FIRST]->setKnots(Dir2, knots[1][0], knots[1][1]);

      child[SECOND]->setKnots(Dir1, mid,        knots[0][1]);
      child[SECOND]->setKnots(Dir2, knots[1][0], knots[1][1]);
    }

    if(splitDirTemp == 1)
    {
      child[FIRST]->setKnots(Dir1, knots[0][0], knots[0][1]);
      child[FIRST]->setKnots(Dir2, knots[1][0], mid);

      child[SECOND]->setKnots(Dir1, knots[0][0], knots[0][1]);
      child[SECOND]->setKnots(Dir2, mid,         knots[1][1]);
    }


    child[FIRST]->setNeighbour(SECOND,  child[SECOND]);
    child[SECOND]->setNeighbour(FIRST, child[FIRST]);


    /*
    AdaptiveBinarytree_PTR  tmpnode;

    if(neighbours[FIRST] != NULL)
    {
      if(neighbours[FIRST]->isLeaf())
        child[FIRST]->setNeighbour(FIRST, NULL);
      else
      {
        tmpnode = neighbours[FIRST]->getChild(SECOND);

        child[FIRST]->setNeighbour(FIRST, tmpnode);
        tmpnode->setNeighbour(SECOND,child[FIRST]);
      }
    }
    
    if(neighbours[SECOND] != NULL)
    {
      if(neighbours[SECOND]->isLeaf())
        child[SECOND]->setNeighbour(SECOND, NULL);
      else
      {
        tmpnode = neighbours[SECOND]->getChild(FIRST);

        child[SECOND]->setNeighbour(SECOND, tmpnode);
        tmpnode->setNeighbour(FIRST,child[SECOND]);
      }
    }
    */

    //cout << " bbbbbbbbbbbb " << endl;

    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii]->GeomData = GeomData;
      child[ii]->setLevel(levTemp);
      child[ii]->setSideTemp(sideTemp);
      child[ii]->setParam3(param3);
      child[ii]->setCoord3(coord3);
      child[ii]->prepareData();
      child[ii]->setSplitDirection(splitDirTemp);
      //cout << " eeeeeeeeee " << endl;
      child[ii]->subDivide(refLev);
    }

    return;
}





template<>
void AdaptiveBinarytree<3>::subDivide(int refLev)
{
    //assert("Current Cell is a LEAF and Subdivision starts" && isLeaf());

    if( !isLeaf() || (refLev < 1) )
      return;


    vector<myPoint>  ptOut;
    vector<int>  cornerInOut(8);

    GeomData->doIntersect3D(bbox, false, cornerInOut, ptOut, domNums) ;

    if( domNums.size() == 1 )
      return;

    if(level == (3*refLev) )
      return;

    int ii=0, jj=0, bb=0, splitDirTemp=0;

    // find the split direction based on the Inside/Outside from of the current cell
    // and the split direction of the parent

    int  levTemp = level+1;

    //if( (cornerInOut[0] == cornerInOut[1]) && (cornerInOut[2] == cornerInOut[3]) )
      //splitDirTemp = 1;
    //else if( (cornerInOut[0] == cornerInOut[2]) && (cornerInOut[1] == cornerInOut[3]) )
      //splitDirTemp = 0;
    //else
    //{
      if(splitDir == 0)
        splitDirTemp = 1;
      else if(splitDir == 1)
        splitDirTemp = 2;
      else
        splitDirTemp = 0;

      //levTemp = level;
    //}
    
    NUM_CHILDREN = 2;
    child      = new AdaptiveBinarytree_PTR[NUM_CHILDREN];

    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii] = new AdaptiveBinarytree<3>(levTemp);

      child[ii]->setParent(this);
    }

    child[FIRST]->orientation = FIRST;
    child[SECOND]->orientation = SECOND;


    double  mid = 0.5*(knots[splitDirTemp][1] + knots[splitDirTemp][0]);

    if(splitDirTemp == 0)
    {
      child[FIRST]->setKnots(Dir1, knots[0][0], mid);
      child[FIRST]->setKnots(Dir2, knots[1][0], knots[1][1]);
      child[FIRST]->setKnots(Dir3, knots[2][0], knots[2][1]);

      child[SECOND]->setKnots(Dir1, mid,        knots[0][1]);
      child[SECOND]->setKnots(Dir2, knots[1][0], knots[1][1]);
      child[SECOND]->setKnots(Dir3, knots[2][0], knots[2][1]);
    }

    if(splitDirTemp == 1)
    {
      child[FIRST]->setKnots(Dir1, knots[0][0], knots[0][1]);
      child[FIRST]->setKnots(Dir2, knots[1][0], mid);
      child[FIRST]->setKnots(Dir3, knots[2][0], knots[2][1]);

      child[SECOND]->setKnots(Dir1, knots[0][0], knots[0][1]);
      child[SECOND]->setKnots(Dir2, mid,         knots[1][1]);
      child[SECOND]->setKnots(Dir3, knots[2][0], knots[2][1]);
    }

    if(splitDirTemp == 2)
    {
      child[FIRST]->setKnots(Dir1, knots[0][0], knots[0][1]);
      child[FIRST]->setKnots(Dir2, knots[1][0], knots[1][1]);
      child[FIRST]->setKnots(Dir3, knots[2][0], mid);

      child[SECOND]->setKnots(Dir1, knots[0][0], knots[0][1]);
      child[SECOND]->setKnots(Dir2, knots[1][0], knots[1][1]);
      child[SECOND]->setKnots(Dir3, mid,         knots[2][1]);
    }


    child[FIRST]->setNeighbour(SECOND,  child[SECOND]);
    child[SECOND]->setNeighbour(FIRST, child[FIRST]);


    //
    AdaptiveBinarytree_PTR  tmpnode;

    if(neighbours[FIRST] != NULL)
    {
      if(neighbours[FIRST]->isLeaf())
        child[FIRST]->setNeighbour(FIRST, NULL);
      else
      {
        tmpnode = neighbours[FIRST]->getChild(SECOND);

        child[FIRST]->setNeighbour(FIRST, tmpnode);
        tmpnode->setNeighbour(SECOND,child[FIRST]);
      }
    }
    
    if(neighbours[SECOND] != NULL)
    {
      if(neighbours[SECOND]->isLeaf())
        child[SECOND]->setNeighbour(SECOND, NULL);
      else
      {
        tmpnode = neighbours[SECOND]->getChild(FIRST);

        child[SECOND]->setNeighbour(SECOND, tmpnode);
        tmpnode->setNeighbour(FIRST,child[SECOND]);
      }
    }
    //

    //cout << " ddddddddddd " << endl;

    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii]->GeomData = GeomData;
      //child[ii]->setLevel(levTemp);
      child[ii]->prepareData();
      child[ii]->setSplitDirection(splitDirTemp);
      child[ii]->subDivide(refLev);
    }

    return;
}



template<>
void AdaptiveBinarytree<2>::printSelf()
{
    printf("\t   ID          = %5d\n", id);
    printf("\t   Level       = %5d\n", level);

    if(parent == NULL)
      printf("\t   Parent      = %5d\n", -1);
    else
      printf("\t   Parent      = %5d\n", parent->getID());

    printf("\t   Parameters ... \n");
    for(int ii=0;ii<2;ii++)
    {
      printf("\t\t direction #%5d ---> \t%12.6f \t %12.6f\n", (ii+1), knots[ii][0], knots[ii][1]);
    }

    printf("\n\t   Neighbours ... \n");
    if(neighbours != NULL)
    {
      if(neighbours[FIRST] != NULL)
        printf("\t\t FIRST   neighbour ID   = %5d \n", neighbours[FIRST]->getID());

      if(neighbours[SECOND] != NULL)
        printf("\t\t SECOND   neighbour ID   = %5d \n", neighbours[SECOND]->getID());
    }
    else
      printf("\t\t No Neighbours \n\n");

    printf("\n\t   children ... \n");
    if(child != NULL)
    {
      if(child[FIRST] != NULL)
        printf("\t\t FIRST   child ID   = %5d \n", child[FIRST]->getID());

      if(child[SECOND] != NULL)
        printf("\t\t SECOND   child ID   = %5d \n", child[SECOND]->getID());
    }
    else
      printf("\t\t No children \n");

    printf("\n\n");

   return;
}



template<>
void AdaptiveBinarytree<3>::printSelf()
{
    printf("\t   ID          = %5d\n", id);
    printf("\t   Level       = %5d\n", level);

    if(parent == NULL)
      printf("\t   Parent      = %5d\n", -1);
    else
      printf("\t   Parent      = %5d\n", parent->getID());

    printf("\t   Basis Functions --->  \n");

    printf("\t   Parameters ... \n");
    for(int ii=0;ii<3;ii++)
    {
      printf("\t\t direction #%5d ---> \t%12.6f \t %12.6f\n", (ii+1), knots[ii][0], knots[ii][1]);
    }

    printf("\n\t   Neighbours ... \n");
    if(neighbours != NULL)
    {
      if(neighbours[FIRST] != NULL)
        printf("\t\t FIRST   neighbour ID   = %5d \n", neighbours[FIRST]->getID());

        if(neighbours[SECOND] != NULL)
          printf("\t\t SECOND   neighbour ID   = %5d \n", neighbours[SECOND]->getID());
    }
    else
      printf("\t\t No Neighbours \n\n");

    printf("\n\t   children ... \n");
    if(child != NULL)
    {
      if(child[FIRST] != NULL)
        printf("\t\t FIRST   child ID   = %5d \n", child[FIRST]->getID());

      if(child[SECOND] != NULL)
        printf("\t\t SECOND   child ID   = %5d \n", child[SECOND]->getID());
    }
    else
      printf("\t\t No children \n");

    printf("\n\n");

   return;
}




template<>
void AdaptiveBinarytree<2>::mergeGaussPoints(int refLev2, int inclDom, int dummy1, int& nGPsMerge, myPoint& ptTemp, double& wt)
{
    AdaptiveBinarytree<2>  *adapIntegNodeLocal = new AdaptiveBinarytree<2>(0);

    adapIntegNodeLocal->setKnots(knots[0][0], knots[0][1], knots[1][0], knots[1][1]);

    adapIntegNodeLocal->GeomData = GeomData;
    adapIntegNodeLocal->domNums = domNums;
    adapIntegNodeLocal->setSplitDirection(splitDir);

    adapIntegNodeLocal->setSideTemp(sideTemp);
    adapIntegNodeLocal->setParam3(param3);
    adapIntegNodeLocal->setCoord3(coord3);

    adapIntegNodeLocal->prepareData();
    adapIntegNodeLocal->subDivide(refLev2);
  
    GaussQuadrature  QuadratureLocal;

    int mergeFlag=0;

    if(sideTemp == -1)
      adapIntegNodeLocal->computeGaussPointsForMerging(0, inclDom, 1, mergeFlag, QuadratureLocal);
    else
      adapIntegNodeLocal->computeGaussPointsForMerging2Dfor3D(0, inclDom, 1, mergeFlag, QuadratureLocal);

    // parametric domain to integration master-quadrilateral domain
    int  ii=0, gp=0;
    myPoint  param, geom;

    ptTemp.setZero();
    param.setZero();
    wt = 0.0;
    for(gp=0; gp<QuadratureLocal.gausspoints.size(); gp++)
    {
      for(ii=0; ii<2; ii++)
        param[ii] = QuadratureLocal.gausspoints[gp][ii] ;

      //cout << gp << '\t' << geom[0] << '\t' << geom[1] << '\t' << QuadratureLocal.gaussweights[gp] << endl;
      wt +=  QuadratureLocal.gaussweights[gp];
    
      ptTemp += param;
    }

    nGPsMerge = QuadratureLocal.gausspoints.size();

    ptTemp /= nGPsMerge;

    delete  adapIntegNodeLocal;

    return;
}





template<>
void AdaptiveBinarytree<3>::mergeGaussPoints(int refLev2, int inclDom, int dummy1, int& nGPsMerge, myPoint& ptTemp, double& wt)
{
    AdaptiveBinarytree<3>  *adapIntegNodeLocal = new AdaptiveBinarytree<3>(0);

    adapIntegNodeLocal->setKnots(knots[0][0], knots[0][1], knots[1][0], knots[1][1], knots[2][0], knots[2][1]);

    adapIntegNodeLocal->GeomData = GeomData;
    adapIntegNodeLocal->domNums = domNums;
    adapIntegNodeLocal->setSplitDirection(-1);

    adapIntegNodeLocal->prepareData();
    adapIntegNodeLocal->subDivide(refLev2);
  
    GaussQuadrature  QuadratureLocal;

    int mergeFlag=0;

    adapIntegNodeLocal->computeGaussPointsForMerging(0, inclDom, 1, mergeFlag, QuadratureLocal);

    // parametric domain to integration master-quadrilateral domain
    int  ii=0, gp=0;
    myPoint  param, geom;

    ptTemp.setZero();
    param.setZero();

    wt = 0.0;
    for(gp=0; gp<QuadratureLocal.gausspoints.size(); gp++)
    {
      for(ii=0; ii<3; ii++)
        param[ii] = QuadratureLocal.gausspoints[gp][ii] ;

      wt +=  QuadratureLocal.gaussweights[gp];
    
      ptTemp += param;
    }

    nGPsMerge = QuadratureLocal.gausspoints.size();

    ptTemp /= nGPsMerge;

    delete  adapIntegNodeLocal;

    return;
}








