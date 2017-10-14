
#include "AdaptiveOctree.h"

template<> int AdaptiveOctree<1>::nodecount = 0;
template<> int AdaptiveOctree<2>::nodecount = 0;
template<> int AdaptiveOctree<3>::nodecount = 0;

template<> int AdaptiveOctree<1>::MAX_LEVEL = 0;
template<> int AdaptiveOctree<2>::MAX_LEVEL = 0;
template<> int AdaptiveOctree<3>::MAX_LEVEL = 0;



template<>
void AdaptiveOctree<1>::subDivide(int lev)
{
    //assert("Current Cell is a LEAF and Subdivision starts" && isLeaf());
    
    if( !isLeaf() )
      return;

    if(level == lev)
      return;

    vector<myPoint>  ptOut;
    vector<int>  vecTemp(4);

    GeomData->doIntersect2D(bbox, false, vecTemp, ptOut, domNums) ;

    if( domNums.size() == 1 )
      return;


    NUM_CHILDREN = 2;
    child  = new AdaptiveOctree_PTR[NUM_CHILDREN];

    int  ii=0, temp = level+1;
    
    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii] = new AdaptiveOctree<1>(temp);

      child[ii]->setParent(this);
    }

    double  mid = 0.5*(knotBegin[0] + knotEnd[0]);

    AdaptiveOctree_PTR  tmpnode;
    
    child[LEFT]->orientation = LEFT;
    child[RIGHT]->orientation = RIGHT;
    
    //cout << " AAAAAAAAAAAAAA " << endl;

    child[LEFT]->setKnots(Dir1,  knotBegin[0], mid);
    child[RIGHT]->setKnots(Dir1, mid, knotEnd[0]);

    // link the children

    child[LEFT]->setNeighbour(RIGHT, child[RIGHT]);
    child[RIGHT]->setNeighbour(LEFT, child[LEFT]);
    
    if(neighbours[LEFT] != NULL)
    {
      if(neighbours[LEFT]->isLeaf())
        child[LEFT]->setNeighbour(LEFT, NULL);
      else
      {
        tmpnode = neighbours[LEFT]->getChild(RIGHT);

        child[LEFT]->setNeighbour(LEFT, tmpnode);
        tmpnode->setNeighbour(RIGHT,child[LEFT]);
      }
    }
    
    if(neighbours[RIGHT] != NULL)
    {
      if(neighbours[RIGHT]->isLeaf())
        child[RIGHT]->setNeighbour(RIGHT, NULL);
      else
      {
        tmpnode = neighbours[RIGHT]->getChild(LEFT);

        child[RIGHT]->setNeighbour(RIGHT, tmpnode);
        tmpnode->setNeighbour(LEFT,child[RIGHT]);
      }
    }

    for(ii=0;ii<NUM_CHILDREN;ii++)
      child[ii]->prepareData();

    if(level < lev)
    {
      for(ii=0;ii<NUM_CHILDREN;ii++)
        child[ii]->subDivide(lev);
    }
    
    return;
}




template<>
void AdaptiveOctree<2>::subDivide(int refLev)
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

    if(level == refLev )
      return;


    int  ii=0, jj=0, bb=0;

    NUM_CHILDREN = 4;

    child  = new AdaptiveOctree_PTR[NUM_CHILDREN];

    int  levTemp = level+1;

    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii] = new AdaptiveOctree<2>(levTemp);

      child[ii]->setParent(this);
    }

    child[SW]->orientation = SW;
    child[SE]->orientation = SE;
    child[NW]->orientation = NW;
    child[NE]->orientation = NE;


    double  midu=0.0, midv=0.0;
    
    midu = 0.5*(knotEnd[0] + knotBegin[0]);
    midv = 0.5*(knotEnd[1] + knotBegin[1]);

    child[SW]->setKnots(Dir1, knotBegin[0], midu);
    child[SW]->setKnots(Dir2, knotBegin[1], midv);

    child[SE]->setKnots(Dir1, midu,        knotEnd[0]);
    child[SE]->setKnots(Dir2, knotBegin[1], midv);

    child[NW]->setKnots(Dir1, knotBegin[0], midu);
    child[NW]->setKnots(Dir2, midv,        knotEnd[1]);

    child[NE]->setKnots(Dir1, midu,        knotEnd[0]);
    child[NE]->setKnots(Dir2, midv,        knotEnd[1]);


    child[SW]->setNeighbour(EAST,  child[SE]);
    child[SW]->setNeighbour(NORTH, child[NW]);
    child[SE]->setNeighbour(WEST,  child[SW]);
    child[SE]->setNeighbour(NORTH, child[NE]);
    child[NW]->setNeighbour(EAST,  child[NE]);
    child[NW]->setNeighbour(SOUTH, child[SW]);
    child[NE]->setNeighbour(WEST,  child[NW]);
    child[NE]->setNeighbour(SOUTH, child[SE]);

    //
    AdaptiveOctree_PTR  tmpnode;

    if(neighbours[WEST] != NULL)
    {
      if(neighbours[WEST]->isLeaf())
      {
        child[SW]->setNeighbour(WEST, NULL);
        child[NW]->setNeighbour(WEST, NULL);
      }
      else
      {
        tmpnode = neighbours[WEST]->getChild(SE);

        child[SW]->setNeighbour(WEST, tmpnode);
        tmpnode->setNeighbour(EAST, child[SW]);

        tmpnode = neighbours[WEST]->getChild(NE);

        child[NW]->setNeighbour(WEST, tmpnode);
        tmpnode->setNeighbour(EAST, child[NW]);
      }
    }

    if(neighbours[EAST] != NULL)
    {
      if(neighbours[EAST]->isLeaf())
      {
        child[SE]->setNeighbour(EAST, NULL);
        child[NE]->setNeighbour(EAST, NULL);
      }
      else
      {
        tmpnode = neighbours[EAST]->getChild(SW);

        child[SE]->setNeighbour(EAST, tmpnode);
        tmpnode->setNeighbour(WEST, child[SE]);

        tmpnode = neighbours[EAST]->getChild(NW);

        child[NE]->setNeighbour(EAST, tmpnode);
        tmpnode->setNeighbour(WEST, child[NE]);
      }
    }

    if(neighbours[NORTH] != NULL)
    {
      if(neighbours[NORTH]->isLeaf())
      {
        child[NW]->setNeighbour(NORTH, NULL);
        child[NE]->setNeighbour(NORTH, NULL);
      }
      else
      {
        tmpnode = neighbours[NORTH]->getChild(SW);

        child[NW]->setNeighbour(NORTH, tmpnode);
        tmpnode->setNeighbour(SOUTH, child[NW]);

        tmpnode = neighbours[NORTH]->getChild(SE);

        child[NE]->setNeighbour(NORTH, tmpnode);
        tmpnode->setNeighbour(SOUTH, child[NE]);
      }
    }

    if(neighbours[SOUTH] != NULL)
    {
      if(neighbours[SOUTH]->isLeaf())
      {
        child[SW]->setNeighbour(SOUTH, NULL);
        child[SE]->setNeighbour(SOUTH, NULL);
      }
      else
      {
        tmpnode = neighbours[SOUTH]->getChild(NW);

        child[SW]->setNeighbour(SOUTH, tmpnode);
        tmpnode->setNeighbour(NORTH, child[SW]);

        tmpnode = neighbours[SOUTH]->getChild(NE);

        child[SE]->setNeighbour(SOUTH, tmpnode);
        tmpnode->setNeighbour(NORTH, child[SE]);
      }
    }
    //

    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii]->GeomData = GeomData;
      child[ii]->setLevel(levTemp);
      child[ii]->setSideTemp(sideTemp);
      child[ii]->setParam3(param3);
      child[ii]->setCoord3(coord3);
      child[ii]->prepareData();

      child[ii]->subDivide(refLev);
    }

    return;
}





template<>
void AdaptiveOctree<3>::subDivide(int refLev)
{
    //assert("Current Cell is a LEAF and Subdivision starts" && isLeaf());

    if( !isLeaf() || (refLev < 1) )
      return;

    vector<myPoint>  ptOut;
    vector<int>  cornerInOut(8);

    GeomData->doIntersect3D(bbox, false, cornerInOut, ptOut, domNums) ;

    if( domNums.size() == 1 )
      return;

    if(level == refLev )
      return;

    int  ii=0, jj=0, bb=0;

    NUM_CHILDREN = 8;
    child = new AdaptiveOctree_PTR[NUM_CHILDREN];

    int  levTemp = level+1;
    
    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii] = new AdaptiveOctree<3>(levTemp);

      child[ii]->setParent(this);
    }

    child[SW_FRONT]->orientation = SW_FRONT;
    child[SE_FRONT]->orientation = SE_FRONT;
    child[NW_FRONT]->orientation = NW_FRONT;
    child[NE_FRONT]->orientation = NE_FRONT;

    child[SW_BACK]->orientation = SW_BACK;
    child[SE_BACK]->orientation = SE_BACK;
    child[NW_BACK]->orientation = NW_BACK;
    child[NE_BACK]->orientation = NE_BACK;
    
    
    ///////////////////////////////////////////////
    //
    // set the knot values

    double  midu = 0.5*(knotBegin[0] + knotEnd[0]);
    double  midv = 0.5*(knotBegin[1] + knotEnd[1]);
    double  midw = 0.5*(knotBegin[2] + knotEnd[2]);

    child[SW_BACK ]->setKnots(Dir1, knotBegin[0], midu);
    child[SW_FRONT]->setKnots(Dir1, knotBegin[0], midu);
    child[NW_BACK ]->setKnots(Dir1, knotBegin[0], midu);
    child[NW_FRONT]->setKnots(Dir1, knotBegin[0], midu);

    child[SE_BACK ]->setKnots(Dir1, midu,        knotEnd[0]);
    child[SE_FRONT]->setKnots(Dir1, midu,        knotEnd[0]);
    child[NE_BACK ]->setKnots(Dir1, midu,        knotEnd[0]);
    child[NE_FRONT]->setKnots(Dir1, midu,        knotEnd[0]);

    child[SW_BACK ]->setKnots(Dir2, knotBegin[1], midv);
    child[SW_FRONT]->setKnots(Dir2, knotBegin[1], midv);
    child[SE_BACK ]->setKnots(Dir2, knotBegin[1], midv);
    child[SE_FRONT]->setKnots(Dir2, knotBegin[1], midv);

    child[NW_BACK ]->setKnots(Dir2, midv,        knotEnd[1]);
    child[NW_FRONT]->setKnots(Dir2, midv,        knotEnd[1]);
    child[NE_BACK ]->setKnots(Dir2, midv,        knotEnd[1]);
    child[NE_FRONT]->setKnots(Dir2, midv,        knotEnd[1]);

    child[SW_BACK ]->setKnots(Dir3, knotBegin[2], midw);
    child[SE_BACK ]->setKnots(Dir3, knotBegin[2], midw);
    child[NW_BACK ]->setKnots(Dir3, knotBegin[2], midw);
    child[NE_BACK ]->setKnots(Dir3, knotBegin[2], midw);

    child[SW_FRONT]->setKnots(Dir3, midw,        knotEnd[2]);
    child[SE_FRONT]->setKnots(Dir3, midw,        knotEnd[2]);
    child[NW_FRONT]->setKnots(Dir3, midw,        knotEnd[2]);
    child[NE_FRONT]->setKnots(Dir3, midw,        knotEnd[2]);


    ///////////////////////////////////////////////
    //
    // set the neighbours among the children


    child[SW_FRONT]->setNeighbour(EAST,  child[SE_FRONT]);
    child[SW_FRONT]->setNeighbour(NORTH, child[NW_FRONT]);
    child[SW_FRONT]->setNeighbour(BACK,  child[SW_BACK]);

    child[SE_FRONT]->setNeighbour(WEST,  child[SW_FRONT]);
    child[SE_FRONT]->setNeighbour(NORTH, child[NE_FRONT]);
    child[SE_FRONT]->setNeighbour(BACK,  child[SE_BACK]);

    child[NW_FRONT]->setNeighbour(EAST,  child[NE_FRONT]);
    child[NW_FRONT]->setNeighbour(SOUTH, child[SW_FRONT]);
    child[NW_FRONT]->setNeighbour(BACK,  child[NW_BACK]);

    child[NE_FRONT]->setNeighbour(WEST,  child[NW_FRONT]);
    child[NE_FRONT]->setNeighbour(SOUTH, child[SE_FRONT]);
    child[NE_FRONT]->setNeighbour(BACK,  child[NE_BACK]);

    child[SW_BACK]->setNeighbour(EAST,  child[SE_BACK]);
    child[SW_BACK]->setNeighbour(NORTH, child[NW_BACK]);
    child[SW_BACK]->setNeighbour(FRONT, child[SW_FRONT]);

    child[SE_BACK]->setNeighbour(WEST,  child[SW_BACK]);
    child[SE_BACK]->setNeighbour(NORTH, child[NE_BACK]);
    child[SE_BACK]->setNeighbour(FRONT, child[SE_FRONT]);

    child[NW_BACK]->setNeighbour(EAST,  child[NE_BACK]);
    child[NW_BACK]->setNeighbour(SOUTH, child[SW_BACK]);
    child[NW_BACK]->setNeighbour(FRONT, child[NW_FRONT]);

    child[NE_BACK]->setNeighbour(WEST,  child[NW_BACK]);
    child[NE_BACK]->setNeighbour(SOUTH, child[SE_BACK]);
    child[NE_BACK]->setNeighbour(FRONT, child[NE_FRONT]);

    ///////////////////////////////////////////////
    //
    // set the neighbours with the children of neighbours

    AdaptiveOctree_PTR  tmpnode, nd;

    nd = neighbours[WEST];
    if(nd != NULL)
    {
      if(nd->isLeaf())
      {
        child[SW_FRONT]->setNeighbour(WEST, NULL);
        child[NW_FRONT]->setNeighbour(WEST, NULL);
        child[SW_BACK ]->setNeighbour(WEST, NULL);
        child[NW_BACK ]->setNeighbour(WEST, NULL);
      }
      else
      {
        tmpnode = nd->getChild(SE_FRONT);

        child[SW_FRONT]->setNeighbour(WEST, tmpnode);
        tmpnode->setNeighbour(EAST, child[SW_FRONT]);

        tmpnode = nd->getChild(NE_FRONT);

        child[NW_FRONT]->setNeighbour(WEST, tmpnode);
        tmpnode->setNeighbour(EAST, child[NW_FRONT]);

        tmpnode = nd->getChild(SE_BACK);

        child[SW_BACK]->setNeighbour(WEST, tmpnode);
        tmpnode->setNeighbour(EAST, child[SW_BACK]);

        tmpnode = nd->getChild(NE_BACK);

        child[NW_BACK]->setNeighbour(WEST, tmpnode);
        tmpnode->setNeighbour(EAST, child[NW_BACK]);
      }
    }

    nd = neighbours[EAST];
    if(nd != NULL)
    {
      if(nd->isLeaf())
      {
        child[SE_FRONT]->setNeighbour(EAST, NULL);
        child[NE_FRONT]->setNeighbour(EAST, NULL);
        child[SE_BACK ]->setNeighbour(EAST, NULL);
        child[NE_BACK ]->setNeighbour(EAST, NULL);
      }
      else
      {
        tmpnode = nd->getChild(SW_FRONT);

        child[SE_FRONT]->setNeighbour(EAST, tmpnode);
        tmpnode->setNeighbour(WEST, child[SE_FRONT]);

        tmpnode = nd->getChild(NW_FRONT);

        child[NE_FRONT]->setNeighbour(EAST, tmpnode);
        tmpnode->setNeighbour(WEST, child[NE_FRONT]);

        tmpnode = nd->getChild(SW_BACK);

        child[SE_BACK]->setNeighbour(EAST, tmpnode);
        tmpnode->setNeighbour(WEST, child[SE_BACK]);

         tmpnode = nd->getChild(NW_BACK);

        child[NE_BACK]->setNeighbour(EAST, tmpnode);
        tmpnode->setNeighbour(WEST, child[NE_BACK]);
      }
    }

    nd = neighbours[SOUTH];
    if(nd != NULL)
    {
      if(nd->isLeaf())
      {
        child[SW_FRONT]->setNeighbour(SOUTH, NULL);
        child[SE_FRONT]->setNeighbour(SOUTH, NULL);
        child[SW_BACK ]->setNeighbour(SOUTH, NULL);
        child[SE_BACK ]->setNeighbour(SOUTH, NULL);
      }
      else
      {
        tmpnode = nd->getChild(NW_FRONT);

        child[SW_FRONT]->setNeighbour(SOUTH, tmpnode);
        tmpnode->setNeighbour(NORTH, child[SW_FRONT]);

        tmpnode = nd->getChild(NE_FRONT);

        child[SE_FRONT]->setNeighbour(SOUTH, tmpnode);
        tmpnode->setNeighbour(NORTH, child[SE_FRONT]);

        tmpnode = nd->getChild(NW_BACK);

        child[SW_BACK]->setNeighbour(SOUTH, tmpnode);
        tmpnode->setNeighbour(NORTH, child[SW_BACK]);

        tmpnode = nd->getChild(NE_BACK);

        child[SE_BACK]->setNeighbour(SOUTH, tmpnode);
        tmpnode->setNeighbour(NORTH, child[SE_BACK]);
      }
    }

    nd = neighbours[NORTH];
    if(nd != NULL)
    {
      if(nd->isLeaf())
      {
        child[NW_FRONT]->setNeighbour(NORTH, NULL);
        child[NE_FRONT]->setNeighbour(NORTH, NULL);
        child[NW_BACK ]->setNeighbour(NORTH, NULL);
        child[NE_BACK ]->setNeighbour(NORTH, NULL);
      }
      else
      {
        tmpnode = nd->getChild(SW_FRONT);

        child[NW_FRONT]->setNeighbour(NORTH, tmpnode);
        tmpnode->setNeighbour(SOUTH, child[NW_FRONT]);

        tmpnode = nd->getChild(SE_FRONT);

        child[NE_FRONT]->setNeighbour(NORTH, tmpnode);
        tmpnode->setNeighbour(SOUTH, child[NE_FRONT]);

        tmpnode = nd->getChild(SW_BACK);

        child[NW_BACK]->setNeighbour(NORTH, tmpnode);
        tmpnode->setNeighbour(SOUTH, child[NW_BACK]);

        tmpnode = nd->getChild(SE_BACK);

        child[NE_BACK]->setNeighbour(NORTH, tmpnode);
        tmpnode->setNeighbour(SOUTH, child[NE_BACK]);
      }
    }

    nd = neighbours[BACK];
    if(nd != NULL)
    {
      if(nd->isLeaf())
      {
        child[SW_BACK]->setNeighbour(BACK, NULL);
        child[SE_BACK]->setNeighbour(BACK, NULL);
        child[NW_BACK]->setNeighbour(BACK, NULL);
        child[NE_BACK]->setNeighbour(BACK, NULL);
      }
      else
      {
        tmpnode = nd->getChild(SW_FRONT);

        child[SW_BACK]->setNeighbour(BACK, tmpnode);
        tmpnode->setNeighbour(FRONT, child[SW_BACK]);

        tmpnode = nd->getChild(SE_FRONT);

        child[SE_BACK]->setNeighbour(BACK, tmpnode);
        tmpnode->setNeighbour(FRONT, child[SE_BACK]);

        tmpnode = nd->getChild(NW_FRONT);

        child[NW_BACK]->setNeighbour(BACK, tmpnode);
        tmpnode->setNeighbour(FRONT, child[NW_BACK]);

        tmpnode = nd->getChild(NE_FRONT);

        child[NE_BACK]->setNeighbour(BACK, tmpnode);
        tmpnode->setNeighbour(FRONT, child[NE_BACK]);
      }
    }

    nd = neighbours[FRONT];
    if(nd != NULL)
    {
      if(nd->isLeaf())
      {
        child[SW_FRONT]->setNeighbour(FRONT, NULL);
        child[SE_FRONT]->setNeighbour(FRONT, NULL);
        child[NW_FRONT]->setNeighbour(FRONT, NULL);
        child[NE_FRONT]->setNeighbour(FRONT, NULL);
      }
      else
      {
        tmpnode = nd->getChild(SW_BACK);

        child[SW_FRONT]->setNeighbour(FRONT, tmpnode);
        tmpnode->setNeighbour(BACK, child[SW_FRONT]);

        tmpnode = nd->getChild(SE_BACK);

        child[SE_FRONT]->setNeighbour(FRONT, tmpnode);
        tmpnode->setNeighbour(BACK, child[SE_FRONT]);

        tmpnode = nd->getChild(NW_BACK);

        child[NW_FRONT]->setNeighbour(FRONT, tmpnode);
        tmpnode->setNeighbour(BACK, child[NW_FRONT]);

        tmpnode = nd->getChild(NE_BACK);

        child[NE_FRONT]->setNeighbour(FRONT, tmpnode);
        tmpnode->setNeighbour(BACK, child[NE_FRONT]);
      }
    }

    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii]->GeomData = GeomData;
      //child[ii]->setLevel(levTemp);
      child[ii]->prepareData();

      child[ii]->subDivide(refLev);
    }

    return;
}




template<>
void AdaptiveOctree<1>::printSelf()
{
    printf("\t   ID          = %5d\n", id);
    printf("\t   Level       = %5d\n", level);

    if(parent == NULL)
      printf("\t   Parent      = %5d\n", -1);
    else
      printf("\t   Parent      = %5d\n", parent->getID());

    printf("\t   children    = %5d\n", NUM_CHILDREN);

    printf("\t   Parameters ... \n");
    printf("\t\t direction #%5d ---> \t%12.6f \t %12.6f\n", 1, knotBegin[0], knotEnd[0]);

    printf("\n\t   Neighbours ... \n");
    if(neighbours != NULL)
    {
      if(neighbours[LEFT] != NULL)
        printf("\t\t LEFT  neighbour ID   = %5d \n", neighbours[LEFT]->getID());

      if(neighbours[RIGHT] != NULL)
        printf("\t\t RIGHT neighbour ID   = %5d \n", neighbours[RIGHT]->getID());
    }
    else
      printf("\t\t No Neighbours \n");

    printf("\n\t   Children ... \n");
    if(child != NULL)
    {
      if(child[LEFT] != NULL)
        printf("\t\t LEFT  neighbour ID   = %5d \n", child[LEFT]->getID());

      if(child[RIGHT] != NULL)
        printf("\t\t RIGHT neighbour ID   = %5d \n", child[RIGHT]->getID());
    }
    else
      printf("\t\t No Children \n");

    printf("\n\n");
    return;
}




template<>
void AdaptiveOctree<2>::printSelf()
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
      printf("\t\t direction #%5d ---> \t%12.6f \t %12.6f\n", (ii+1), knotBegin[ii], knotEnd[ii]);
    }

    printf("\n\t   Neighbours ... \n");
    if(neighbours != NULL)
    {
      if(neighbours[EAST] != NULL)
        printf("\t\t EAST   neighbour ID   = %5d \n", neighbours[EAST]->getID());

      if(neighbours[WEST] != NULL)
        printf("\t\t WEST   neighbour ID   = %5d \n", neighbours[WEST]->getID());

      if(neighbours[NORTH] != NULL)
        printf("\t\t NORTH  neighbour ID   = %5d \n", neighbours[NORTH]->getID());

      if(neighbours[SOUTH] != NULL)
        printf("\t\t SOUTH  neighbour ID   = %5d \n", neighbours[SOUTH]->getID());

    }
    else
      printf("\t\t No Neighbours \n\n");

    printf("\n\t   children ... \n");
    if(child != NULL)
    {
      //printf("\t  # of children    = %5d\n", NUM_CHILDREN);
               
      if(child[SW] != NULL)
        printf("\t\t SW   child ID   = %5d \n", child[SW]->getID());

      if(child[SE] != NULL)
        printf("\t\t SE   child ID   = %5d \n", child[SE]->getID());

      if(child[NW] != NULL)
        printf("\t\t NW   child ID   = %5d \n", child[NW]->getID());

      if(child[NE] != NULL)
        printf("\t\t NE   child ID   = %5d \n", child[NE]->getID());

    }
    else
      printf("\t\t No children \n");
    
    printf("\n\n");

    return;
}



template<>
void AdaptiveOctree<3>::printSelf()
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
      printf("\t\t direction #%5d ---> \t%12.6f \t %12.6f\n", (ii+1), knotBegin[ii], knotEnd[ii]);
    }

    printf("\n\t   Neighbours ... \n");
    if(neighbours != NULL)
    {
      if(neighbours[WEST]  != NULL)
        printf("\t\t WEST   neighbour ID   = %5d \n", neighbours[WEST]->getID());
      else
        printf("\t\t WEST   neighbour ID   = %5d \n", -1);

      if(neighbours[EAST]  != NULL)
        printf("\t\t EAST   neighbour ID   = %5d \n", neighbours[EAST]->getID());
      else
        printf("\t\t EAST   neighbour ID   = %5d \n", -1);

      if(neighbours[SOUTH] != NULL)
        printf("\t\t SOUTH  neighbour ID   = %5d \n", neighbours[SOUTH]->getID());
      else
        printf("\t\t SOUTH  neighbour ID   = %5d \n", -1);

      if(neighbours[NORTH] != NULL)
        printf("\t\t NORTH  neighbour ID   = %5d \n", neighbours[NORTH]->getID());
      else
        printf("\t\t NORTH  neighbour ID   = %5d \n", -1);

      if(neighbours[FRONT] != NULL)
        printf("\t\t FRONT  neighbour ID   = %5d \n", neighbours[FRONT]->getID());
      else
        printf("\t\t FRONT  neighbour ID   = %5d \n", -1);

      if(neighbours[BACK]  != NULL)
        printf("\t\t BACK   neighbour ID   = %5d \n", neighbours[BACK]->getID());
      else
        printf("\t\t BACK   neighbour ID   = %5d \n", -1);
    }
    else
      printf("\t\t No Neighbours \n\n");

    printf("\n\t   children ... \n");
    if(child != NULL)
    {
      //printf("\t  # of children    = %5d\n", NUM_CHILDREN);
               
      if(child[SW] != NULL)
        printf("\t\t SW   child ID   = %5d \n", child[SW]->getID());

      if(child[SE] != NULL)
        printf("\t\t SE   child ID   = %5d \n", child[SE]->getID());

      if(child[NW] != NULL)
        printf("\t\t NW   child ID   = %5d \n", child[NW]->getID());

      if(child[NE] != NULL)
        printf("\t\t NE   child ID   = %5d \n", child[NE]->getID());

    }
    else
      printf("\t\t No children \n");
    
    printf("\n\n");

    return;
}




template<>
void AdaptiveOctree<2>::mergeGaussPoints(int refLev2, int inclDom, int dummy1, int& nGPsMerge, myPoint& ptTemp, double& wt)
{
    AdaptiveOctree<2>  *adapIntegNodeLocal = new AdaptiveOctree<2>(0);

    adapIntegNodeLocal->setKnots(knotBegin[0], knotEnd[0], knotBegin[1], knotEnd[1]);

    //cout << " AAAAAAAAAA " << endl;

    adapIntegNodeLocal->GeomData = GeomData;
    adapIntegNodeLocal->domNums = domNums;

    //cout << " AAAAAAAAAA .... " << refLev << endl;

    adapIntegNodeLocal->setSideTemp(sideTemp);
    adapIntegNodeLocal->setParam3(param3);
    adapIntegNodeLocal->setCoord3(coord3);
    adapIntegNodeLocal->prepareData();
    adapIntegNodeLocal->subDivide(refLev2);
  
    //cout << " AAAAAAAAAA .... " << refLev << endl;

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

      //GeomData->computeCoord(param, geom);
      //ptTemp += geom;

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
void AdaptiveOctree<3>::mergeGaussPoints(int refLev2, int inclDom, int dummy1, int& nGPsMerge, myPoint& ptTemp, double& wt)
{
    AdaptiveOctree<3>  *adapIntegNodeLocal = new AdaptiveOctree<3>(0);

    adapIntegNodeLocal->setKnots(knotBegin[0], knotEnd[0], knotBegin[1], knotEnd[1], knotBegin[2], knotEnd[2]);

    adapIntegNodeLocal->GeomData = GeomData;
    adapIntegNodeLocal->domNums = domNums;

    adapIntegNodeLocal->prepareData();
    adapIntegNodeLocal->subDivide(refLev2);
  
    //cout << " AAAAAAAAAA .... " << refLev << endl;

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

      //GeomData->computeCoord(param, geom);
      //ptTemp += geom;

      wt +=  QuadratureLocal.gaussweights[gp];
    
      ptTemp += param;
    }

    nGPsMerge = QuadratureLocal.gausspoints.size();

    ptTemp /= nGPsMerge;

    delete  adapIntegNodeLocal;

   return;
}











