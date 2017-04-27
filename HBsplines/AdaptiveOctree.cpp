
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
    //assert("Current Cell is a LEAF and Subdivision starts" && IsLeaf());
    
    if( !IsLeaf() )
      return;

    if(level == lev)
      return;


    vector<myPoint>  ptOut;
    vector<int>  vecTemp(4);

    GeomData->doIntersect2D(bbox, false, vecTemp, ptOut, domNums) ;

    if( domNums.size() == 1 )
      return;

    //cout << " bbbbbbbbbbb " << endl;


    NUM_CHILDREN = 2;
    child  = new AdaptiveOctree_PTR[NUM_CHILDREN];

    int ii, temp;
    
    temp = level+1;
    
    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii] = new AdaptiveOctree<1>(temp);

      child[ii]->SetParent(this);
    }

    double  mid(0.0);

    mid = 0.5*(knots[0][0] + knots[0][1]);

    AdaptiveOctree_PTR  tmpnode;
    
    child[LEFT]->orientation = LEFT;
    child[RIGHT]->orientation = RIGHT;
    
    //cout << " AAAAAAAAAAAAAA " << endl;

    child[LEFT]->SetKnots(Dir1, knots[0][0], mid);
    child[RIGHT]->SetKnots(Dir1, mid, knots[0][1]);

    // link the children

    child[LEFT]->SetNeighbour(RIGHT, child[RIGHT]);
    child[RIGHT]->SetNeighbour(LEFT, child[LEFT]);
    
    if(neighbours[LEFT] != NULL)
    {
      if(neighbours[LEFT]->IsLeaf())
        child[LEFT]->SetNeighbour(LEFT, NULL);
      else
      {
        tmpnode = neighbours[LEFT]->GetChild(RIGHT);

        child[LEFT]->SetNeighbour(LEFT, tmpnode);
        tmpnode->SetNeighbour(RIGHT,child[LEFT]);
      }
    }
    
    if(neighbours[RIGHT] != NULL)
    {
      if(neighbours[RIGHT]->IsLeaf())
        child[RIGHT]->SetNeighbour(RIGHT, NULL);
      else
      {
        tmpnode = neighbours[RIGHT]->GetChild(LEFT);

        child[RIGHT]->SetNeighbour(RIGHT, tmpnode);
        tmpnode->SetNeighbour(LEFT,child[RIGHT]);
      }
    }

    //cout << " bbbbbbbbbbb " << endl;
    for(ii=0;ii<NUM_CHILDREN;ii++)
      child[ii]->prepareData();

    //cout << " bbbbbbbbbbb " << endl;

    if(level < lev)
    {
      for(ii=0;ii<NUM_CHILDREN;ii++)
        child[ii]->subDivide(lev);
    }
    //cout << " bbbbbbbbbbb " << endl;
    
    return;
}




template<>
void AdaptiveOctree<2>::subDivide(int refLev)
{
    //assert("Current Cell is a LEAF and Subdivision starts" && IsLeaf());

    if( !IsLeaf() || (refLev < 1) )
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


    int ii, levTemp, jj, bb;

    NUM_CHILDREN = 4;

    child      = new AdaptiveOctree_PTR[NUM_CHILDREN];

    levTemp = level+1;

    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii] = new AdaptiveOctree<2>(levTemp);

      child[ii]->SetParent(this);
    }

    child[SW]->orientation = SW;
    child[SE]->orientation = SE;
    child[NW]->orientation = NW;
    child[NE]->orientation = NE;


    double  midu(0.0), midv(0.0);
    
    midu = 0.5*(knots[0][1] + knots[0][0]);
    midv = 0.5*(knots[1][1] + knots[1][0]);

    child[SW]->SetKnots(Dir1, knots[0][0], midu);
    child[SW]->SetKnots(Dir2, knots[1][0], midv);

    child[SE]->SetKnots(Dir1, midu,        knots[0][1]);
    child[SE]->SetKnots(Dir2, knots[1][0], midv);

    child[NW]->SetKnots(Dir1, knots[0][0], midu);
    child[NW]->SetKnots(Dir2, midv,        knots[1][1]);

    child[NE]->SetKnots(Dir1, midu,        knots[0][1]);
    child[NE]->SetKnots(Dir2, midv,        knots[1][1]);


    child[SW]->SetNeighbour(EAST,  child[SE]);
    child[SW]->SetNeighbour(NORTH, child[NW]);
    child[SE]->SetNeighbour(WEST,  child[SW]);
    child[SE]->SetNeighbour(NORTH, child[NE]);
    child[NW]->SetNeighbour(EAST,  child[NE]);
    child[NW]->SetNeighbour(SOUTH, child[SW]);
    child[NE]->SetNeighbour(WEST,  child[NW]);
    child[NE]->SetNeighbour(SOUTH, child[SE]);

    //
    AdaptiveOctree_PTR  tmpnode;

    if(neighbours[WEST] != NULL)
    {
      if(neighbours[WEST]->IsLeaf())
      {
        child[SW]->SetNeighbour(WEST, NULL);
        child[NW]->SetNeighbour(WEST, NULL);
      }
      else
      {
        tmpnode = neighbours[WEST]->GetChild(SE);

        child[SW]->SetNeighbour(WEST, tmpnode);
        tmpnode->SetNeighbour(EAST, child[SW]);

        tmpnode = neighbours[WEST]->GetChild(NE);

        child[NW]->SetNeighbour(WEST, tmpnode);
        tmpnode->SetNeighbour(EAST, child[NW]);
      }
    }

    if(neighbours[EAST] != NULL)
    {
      if(neighbours[EAST]->IsLeaf())
      {
        child[SE]->SetNeighbour(EAST, NULL);
        child[NE]->SetNeighbour(EAST, NULL);
      }
      else
      {
        tmpnode = neighbours[EAST]->GetChild(SW);

        child[SE]->SetNeighbour(EAST, tmpnode);
        tmpnode->SetNeighbour(WEST, child[SE]);

        tmpnode = neighbours[EAST]->GetChild(NW);

        child[NE]->SetNeighbour(EAST, tmpnode);
        tmpnode->SetNeighbour(WEST, child[NE]);
      }
    }

    if(neighbours[NORTH] != NULL)
    {
      if(neighbours[NORTH]->IsLeaf())
      {
        child[NW]->SetNeighbour(NORTH, NULL);
        child[NE]->SetNeighbour(NORTH, NULL);
      }
      else
      {
        tmpnode = neighbours[NORTH]->GetChild(SW);

        child[NW]->SetNeighbour(NORTH, tmpnode);
        tmpnode->SetNeighbour(SOUTH, child[NW]);

        tmpnode = neighbours[NORTH]->GetChild(SE);

        child[NE]->SetNeighbour(NORTH, tmpnode);
        tmpnode->SetNeighbour(SOUTH, child[NE]);
      }
    }

    if(neighbours[SOUTH] != NULL)
    {
      if(neighbours[SOUTH]->IsLeaf())
      {
        child[SW]->SetNeighbour(SOUTH, NULL);
        child[SE]->SetNeighbour(SOUTH, NULL);
      }
      else
      {
        tmpnode = neighbours[SOUTH]->GetChild(NW);

        child[SW]->SetNeighbour(SOUTH, tmpnode);
        tmpnode->SetNeighbour(NORTH, child[SW]);

        tmpnode = neighbours[SOUTH]->GetChild(NE);

        child[SE]->SetNeighbour(SOUTH, tmpnode);
        tmpnode->SetNeighbour(NORTH, child[SE]);
      }
    }
    //

    //cout << " bbbbbbbbbbbb " << endl;

    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii]->GeomData = GeomData;
      child[ii]->SetLevel(levTemp);
      child[ii]->SetSideTemp(sideTemp);
      child[ii]->SetParam3(param3);
      child[ii]->SetCoord3(coord3);
      child[ii]->prepareData();

      //cout << " eeeeeeeeee " << endl;
      child[ii]->subDivide(refLev);
    }

    //cout << " ddddddddddd " << endl;

    return;
}





template<>
void AdaptiveOctree<3>::subDivide(int refLev)
{
    //assert("Current Cell is a LEAF and Subdivision starts" && IsLeaf());

    if( !IsLeaf() || (refLev < 1) )
      return;


    vector<myPoint>  ptOut;
    vector<int>  cornerInOut(8);

    GeomData->doIntersect3D(bbox, false, cornerInOut, ptOut, domNums) ;

    if( domNums.size() == 1 )
      return;

    if(level == refLev )
      return;

    int ii, levTemp, jj, bb;

    NUM_CHILDREN = 8;
    child = new AdaptiveOctree_PTR[NUM_CHILDREN];

    levTemp = level+1;
    
    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii] = new AdaptiveOctree<3>(levTemp);

      child[ii]->SetParent(this);
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
    // set the knots


    double  midu, midv, midw;
    
    midu = 0.5*(knots[0][0] + knots[0][1]);
    midv = 0.5*(knots[1][0] + knots[1][1]);
    midw = 0.5*(knots[2][0] + knots[2][1]);

    child[SW_BACK ]->SetKnots(Dir1, knots[0][0], midu);
    child[SW_FRONT]->SetKnots(Dir1, knots[0][0], midu);
    child[NW_BACK ]->SetKnots(Dir1, knots[0][0], midu);
    child[NW_FRONT]->SetKnots(Dir1, knots[0][0], midu);

    child[SE_BACK ]->SetKnots(Dir1, midu,        knots[0][1]);
    child[SE_FRONT]->SetKnots(Dir1, midu,        knots[0][1]);
    child[NE_BACK ]->SetKnots(Dir1, midu,        knots[0][1]);
    child[NE_FRONT]->SetKnots(Dir1, midu,        knots[0][1]);


    child[SW_BACK ]->SetKnots(Dir2, knots[1][0], midv);
    child[SW_FRONT]->SetKnots(Dir2, knots[1][0], midv);
    child[SE_BACK ]->SetKnots(Dir2, knots[1][0], midv);
    child[SE_FRONT]->SetKnots(Dir2, knots[1][0], midv);

    child[NW_BACK ]->SetKnots(Dir2, midv,        knots[1][1]);
    child[NW_FRONT]->SetKnots(Dir2, midv,        knots[1][1]);
    child[NE_BACK ]->SetKnots(Dir2, midv,        knots[1][1]);
    child[NE_FRONT]->SetKnots(Dir2, midv,        knots[1][1]);


    child[SW_BACK ]->SetKnots(Dir3, knots[2][0], midw);
    child[SE_BACK ]->SetKnots(Dir3, knots[2][0], midw);
    child[NW_BACK ]->SetKnots(Dir3, knots[2][0], midw);
    child[NE_BACK ]->SetKnots(Dir3, knots[2][0], midw);

    child[SW_FRONT]->SetKnots(Dir3, midw,        knots[2][1]);
    child[SE_FRONT]->SetKnots(Dir3, midw,        knots[2][1]);
    child[NW_FRONT]->SetKnots(Dir3, midw,        knots[2][1]);
    child[NE_FRONT]->SetKnots(Dir3, midw,        knots[2][1]);


    ///////////////////////////////////////////////
    //
    // set the neighbours among the children


    child[SW_FRONT]->SetNeighbour(EAST,  child[SE_FRONT]);
    child[SW_FRONT]->SetNeighbour(NORTH, child[NW_FRONT]);
    child[SW_FRONT]->SetNeighbour(BACK,  child[SW_BACK]);

    child[SE_FRONT]->SetNeighbour(WEST,  child[SW_FRONT]);
    child[SE_FRONT]->SetNeighbour(NORTH, child[NE_FRONT]);
    child[SE_FRONT]->SetNeighbour(BACK,  child[SE_BACK]);

    child[NW_FRONT]->SetNeighbour(EAST,  child[NE_FRONT]);
    child[NW_FRONT]->SetNeighbour(SOUTH, child[SW_FRONT]);
    child[NW_FRONT]->SetNeighbour(BACK,  child[NW_BACK]);

    child[NE_FRONT]->SetNeighbour(WEST,  child[NW_FRONT]);
    child[NE_FRONT]->SetNeighbour(SOUTH, child[SE_FRONT]);
    child[NE_FRONT]->SetNeighbour(BACK,  child[NE_BACK]);

    child[SW_BACK]->SetNeighbour(EAST,  child[SE_BACK]);
    child[SW_BACK]->SetNeighbour(NORTH, child[NW_BACK]);
    child[SW_BACK]->SetNeighbour(FRONT, child[SW_FRONT]);

    child[SE_BACK]->SetNeighbour(WEST,  child[SW_BACK]);
    child[SE_BACK]->SetNeighbour(NORTH, child[NE_BACK]);
    child[SE_BACK]->SetNeighbour(FRONT, child[SE_FRONT]);

    child[NW_BACK]->SetNeighbour(EAST,  child[NE_BACK]);
    child[NW_BACK]->SetNeighbour(SOUTH, child[SW_BACK]);
    child[NW_BACK]->SetNeighbour(FRONT, child[NW_FRONT]);

    child[NE_BACK]->SetNeighbour(WEST,  child[NW_BACK]);
    child[NE_BACK]->SetNeighbour(SOUTH, child[SE_BACK]);
    child[NE_BACK]->SetNeighbour(FRONT, child[NE_FRONT]);

    ///////////////////////////////////////////////
    //
    // set the neighbours with the children of neighbours

    AdaptiveOctree_PTR  tmpnode, nd;

    nd = neighbours[WEST];
    if(nd != NULL)
    {
      if(nd->IsLeaf())
      {
        child[SW_FRONT]->SetNeighbour(WEST, NULL);
        child[NW_FRONT]->SetNeighbour(WEST, NULL);
        child[SW_BACK ]->SetNeighbour(WEST, NULL);
        child[NW_BACK ]->SetNeighbour(WEST, NULL);
      }
      else
      {
        tmpnode = nd->GetChild(SE_FRONT);

        child[SW_FRONT]->SetNeighbour(WEST, tmpnode);
        tmpnode->SetNeighbour(EAST, child[SW_FRONT]);

        tmpnode = nd->GetChild(NE_FRONT);

        child[NW_FRONT]->SetNeighbour(WEST, tmpnode);
        tmpnode->SetNeighbour(EAST, child[NW_FRONT]);

        tmpnode = nd->GetChild(SE_BACK);

        child[SW_BACK]->SetNeighbour(WEST, tmpnode);
        tmpnode->SetNeighbour(EAST, child[SW_BACK]);

        tmpnode = nd->GetChild(NE_BACK);

        child[NW_BACK]->SetNeighbour(WEST, tmpnode);
        tmpnode->SetNeighbour(EAST, child[NW_BACK]);
      }
    }

    nd = neighbours[EAST];
    if(nd != NULL)
    {
      if(nd->IsLeaf())
      {
        child[SE_FRONT]->SetNeighbour(EAST, NULL);
        child[NE_FRONT]->SetNeighbour(EAST, NULL);
        child[SE_BACK ]->SetNeighbour(EAST, NULL);
        child[NE_BACK ]->SetNeighbour(EAST, NULL);
      }
      else
      {
        tmpnode = nd->GetChild(SW_FRONT);

        child[SE_FRONT]->SetNeighbour(EAST, tmpnode);
        tmpnode->SetNeighbour(WEST, child[SE_FRONT]);

        tmpnode = nd->GetChild(NW_FRONT);

        child[NE_FRONT]->SetNeighbour(EAST, tmpnode);
        tmpnode->SetNeighbour(WEST, child[NE_FRONT]);

        tmpnode = nd->GetChild(SW_BACK);

        child[SE_BACK]->SetNeighbour(EAST, tmpnode);
        tmpnode->SetNeighbour(WEST, child[SE_BACK]);

         tmpnode = nd->GetChild(NW_BACK);

        child[NE_BACK]->SetNeighbour(EAST, tmpnode);
        tmpnode->SetNeighbour(WEST, child[NE_BACK]);
      }
    }

    nd = neighbours[SOUTH];
    if(nd != NULL)
    {
      if(nd->IsLeaf())
      {
        child[SW_FRONT]->SetNeighbour(SOUTH, NULL);
        child[SE_FRONT]->SetNeighbour(SOUTH, NULL);
        child[SW_BACK ]->SetNeighbour(SOUTH, NULL);
        child[SE_BACK ]->SetNeighbour(SOUTH, NULL);
      }
      else
      {
        tmpnode = nd->GetChild(NW_FRONT);

        child[SW_FRONT]->SetNeighbour(SOUTH, tmpnode);
        tmpnode->SetNeighbour(NORTH, child[SW_FRONT]);

        tmpnode = nd->GetChild(NE_FRONT);

        child[SE_FRONT]->SetNeighbour(SOUTH, tmpnode);
        tmpnode->SetNeighbour(NORTH, child[SE_FRONT]);

        tmpnode = nd->GetChild(NW_BACK);

        child[SW_BACK]->SetNeighbour(SOUTH, tmpnode);
        tmpnode->SetNeighbour(NORTH, child[SW_BACK]);

        tmpnode = nd->GetChild(NE_BACK);

        child[SE_BACK]->SetNeighbour(SOUTH, tmpnode);
        tmpnode->SetNeighbour(NORTH, child[SE_BACK]);
      }
    }

    nd = neighbours[NORTH];
    if(nd != NULL)
    {
      if(nd->IsLeaf())
      {
        child[NW_FRONT]->SetNeighbour(NORTH, NULL);
        child[NE_FRONT]->SetNeighbour(NORTH, NULL);
        child[NW_BACK ]->SetNeighbour(NORTH, NULL);
        child[NE_BACK ]->SetNeighbour(NORTH, NULL);
      }
      else
      {
        tmpnode = nd->GetChild(SW_FRONT);

        child[NW_FRONT]->SetNeighbour(NORTH, tmpnode);
        tmpnode->SetNeighbour(SOUTH, child[NW_FRONT]);

        tmpnode = nd->GetChild(SE_FRONT);

        child[NE_FRONT]->SetNeighbour(NORTH, tmpnode);
        tmpnode->SetNeighbour(SOUTH, child[NE_FRONT]);

        tmpnode = nd->GetChild(SW_BACK);

        child[NW_BACK]->SetNeighbour(NORTH, tmpnode);
        tmpnode->SetNeighbour(SOUTH, child[NW_BACK]);

        tmpnode = nd->GetChild(SE_BACK);

        child[NE_BACK]->SetNeighbour(NORTH, tmpnode);
        tmpnode->SetNeighbour(SOUTH, child[NE_BACK]);
      }
    }

    nd = neighbours[BACK];
    if(nd != NULL)
    {
      if(nd->IsLeaf())
      {
        child[SW_BACK]->SetNeighbour(BACK, NULL);
        child[SE_BACK]->SetNeighbour(BACK, NULL);
        child[NW_BACK]->SetNeighbour(BACK, NULL);
        child[NE_BACK]->SetNeighbour(BACK, NULL);
      }
      else
      {
        tmpnode = nd->GetChild(SW_FRONT);

        child[SW_BACK]->SetNeighbour(BACK, tmpnode);
        tmpnode->SetNeighbour(FRONT, child[SW_BACK]);

        tmpnode = nd->GetChild(SE_FRONT);

        child[SE_BACK]->SetNeighbour(BACK, tmpnode);
        tmpnode->SetNeighbour(FRONT, child[SE_BACK]);

        tmpnode = nd->GetChild(NW_FRONT);

        child[NW_BACK]->SetNeighbour(BACK, tmpnode);
        tmpnode->SetNeighbour(FRONT, child[NW_BACK]);

        tmpnode = nd->GetChild(NE_FRONT);

        child[NE_BACK]->SetNeighbour(BACK, tmpnode);
        tmpnode->SetNeighbour(FRONT, child[NE_BACK]);
      }
    }

    nd = neighbours[FRONT];
    if(nd != NULL)
    {
      if(nd->IsLeaf())
      {
        child[SW_FRONT]->SetNeighbour(FRONT, NULL);
        child[SE_FRONT]->SetNeighbour(FRONT, NULL);
        child[NW_FRONT]->SetNeighbour(FRONT, NULL);
        child[NE_FRONT]->SetNeighbour(FRONT, NULL);
      }
      else
      {
        tmpnode = nd->GetChild(SW_BACK);

        child[SW_FRONT]->SetNeighbour(FRONT, tmpnode);
        tmpnode->SetNeighbour(BACK, child[SW_FRONT]);

        tmpnode = nd->GetChild(SE_BACK);

        child[SE_FRONT]->SetNeighbour(FRONT, tmpnode);
        tmpnode->SetNeighbour(BACK, child[SE_FRONT]);

        tmpnode = nd->GetChild(NW_BACK);

        child[NW_FRONT]->SetNeighbour(FRONT, tmpnode);
        tmpnode->SetNeighbour(BACK, child[NW_FRONT]);

        tmpnode = nd->GetChild(NE_BACK);

        child[NE_FRONT]->SetNeighbour(FRONT, tmpnode);
        tmpnode->SetNeighbour(BACK, child[NE_FRONT]);
      }
    }

    //cout << " ddddddddddd " << endl;

    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
      child[ii]->GeomData = GeomData;
      //child[ii]->SetLevel(levTemp);
      child[ii]->prepareData();

      child[ii]->subDivide(refLev);
    }

    //cout << " bbbbbbbbbbb " << endl;

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
              printf("\t   Parent      = %5d\n", parent->GetID());

            printf("\t   children    = %5d\n", NUM_CHILDREN);

            printf("\t   Parameters ... \n");
            printf("\t\t direction #%5d ---> \t%12.6f \t %12.6f\n", 1, knots[0][0], knots[0][1]);

            printf("\n\t   Neighbours ... \n");
            if(neighbours != NULL)
            {
               if(neighbours[LEFT] != NULL)
                 printf("\t\t LEFT  neighbour ID   = %5d \n", neighbours[LEFT]->GetID());

               if(neighbours[RIGHT] != NULL)
                 printf("\t\t RIGHT neighbour ID   = %5d \n", neighbours[RIGHT]->GetID());
            }
            else
              printf("\t\t No Neighbours \n");

            printf("\n\t   Children ... \n");
            if(child != NULL)
            {
               if(child[LEFT] != NULL)
                 printf("\t\t LEFT  neighbour ID   = %5d \n", child[LEFT]->GetID());

               if(child[RIGHT] != NULL)
                 printf("\t\t RIGHT neighbour ID   = %5d \n", child[RIGHT]->GetID());
            }
            else
              printf("\t\t No Children \n");
            printf("\n\n");
}






template<>
void AdaptiveOctree<2>::printSelf()
{
            printf("\t   ID          = %5d\n", id);
            printf("\t   Level       = %5d\n", level);

            if(parent == NULL)
              printf("\t   Parent      = %5d\n", -1);
            else
              printf("\t   Parent      = %5d\n", parent->GetID());

            printf("\t   Parameters ... \n");
            for(int ii=0;ii<2;ii++)
            {
               printf("\t\t direction #%5d ---> \t%12.6f \t %12.6f\n", (ii+1), knots[ii][0], knots[ii][1]);
            }

            printf("\n\t   Neighbours ... \n");
            if(neighbours != NULL)
            {
               if(neighbours[EAST] != NULL)
                 printf("\t\t EAST   neighbour ID   = %5d \n", neighbours[EAST]->GetID());

               if(neighbours[WEST] != NULL)
                 printf("\t\t WEST   neighbour ID   = %5d \n", neighbours[WEST]->GetID());

               if(neighbours[NORTH] != NULL)
                 printf("\t\t NORTH  neighbour ID   = %5d \n", neighbours[NORTH]->GetID());

               if(neighbours[SOUTH] != NULL)
                 printf("\t\t SOUTH  neighbour ID   = %5d \n", neighbours[SOUTH]->GetID());

            }
            else
              printf("\t\t No Neighbours \n\n");

            printf("\n\t   children ... \n");
            if(child != NULL)
            {
               //printf("\t  # of children    = %5d\n", NUM_CHILDREN);
               
               if(child[SW] != NULL)
                 printf("\t\t SW   child ID   = %5d \n", child[SW]->GetID());

               if(child[SE] != NULL)
                 printf("\t\t SE   child ID   = %5d \n", child[SE]->GetID());

               if(child[NW] != NULL)
                 printf("\t\t NW   child ID   = %5d \n", child[NW]->GetID());

               if(child[NE] != NULL)
                 printf("\t\t NE   child ID   = %5d \n", child[NE]->GetID());

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
              printf("\t   Parent      = %5d\n", parent->GetID());

            printf("\t   Basis Functions --->  \n");

            printf("\t   Parameters ... \n");
            for(int ii=0;ii<3;ii++)
            {
               printf("\t\t direction #%5d ---> \t%12.6f \t %12.6f\n", (ii+1), knots[ii][0], knots[ii][1]);
            }

            printf("\n\t   Neighbours ... \n");
            if(neighbours != NULL)
            {
               if(neighbours[WEST]  != NULL)
                 printf("\t\t WEST   neighbour ID   = %5d \n", neighbours[WEST]->GetID());
               else
                 printf("\t\t WEST   neighbour ID   = %5d \n", -1);

               if(neighbours[EAST]  != NULL)
                 printf("\t\t EAST   neighbour ID   = %5d \n", neighbours[EAST]->GetID());
               else
                 printf("\t\t EAST   neighbour ID   = %5d \n", -1);

               if(neighbours[SOUTH] != NULL)
                 printf("\t\t SOUTH  neighbour ID   = %5d \n", neighbours[SOUTH]->GetID());
               else
                 printf("\t\t SOUTH  neighbour ID   = %5d \n", -1);

               if(neighbours[NORTH] != NULL)
                 printf("\t\t NORTH  neighbour ID   = %5d \n", neighbours[NORTH]->GetID());
               else
                 printf("\t\t NORTH  neighbour ID   = %5d \n", -1);

               if(neighbours[FRONT] != NULL)
                 printf("\t\t FRONT  neighbour ID   = %5d \n", neighbours[FRONT]->GetID());
               else
                 printf("\t\t FRONT  neighbour ID   = %5d \n", -1);

               if(neighbours[BACK]  != NULL)
                 printf("\t\t BACK   neighbour ID   = %5d \n", neighbours[BACK]->GetID());
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
                 printf("\t\t SW   child ID   = %5d \n", child[SW]->GetID());

               if(child[SE] != NULL)
                 printf("\t\t SE   child ID   = %5d \n", child[SE]->GetID());

               if(child[NW] != NULL)
                 printf("\t\t NW   child ID   = %5d \n", child[NW]->GetID());

               if(child[NE] != NULL)
                 printf("\t\t NE   child ID   = %5d \n", child[NE]->GetID());

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

//  cout << knots[0][0] << '\t' << knots[0][1] << '\t' << knots[0][2] << '\t' << knots[0][3] << endl;
//  cout << knots[1][0] << '\t' << knots[1][1] << '\t' << knots[1][2] << '\t' << knots[1][3] << endl;

  adapIntegNodeLocal->SetKnots(knots[0][0], knots[0][1], knots[1][0], knots[1][1]);

  //cout << " AAAAAAAAAA " << endl;

  adapIntegNodeLocal->GeomData = GeomData;
  adapIntegNodeLocal->domNums = domNums;

  //cout << " AAAAAAAAAA .... " << refLev << endl;

  adapIntegNodeLocal->SetSideTemp(sideTemp);
  adapIntegNodeLocal->SetParam3(param3);
  adapIntegNodeLocal->SetCoord3(coord3);

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
  int  ii, gp;
  myPoint  param, geom;
  wt = 0.0;
  ptTemp.setZero();
  param.setZero();
  for(gp=0; gp<QuadratureLocal.gausspoints.size(); gp++)
  {
    for(ii=0; ii<2; ii++)
      param[ii] = QuadratureLocal.gausspoints[gp][ii] ;

    //GeomData->ComputeCoord(param, geom);
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

//  cout << knots[0][0] << '\t' << knots[0][1] << '\t' << knots[0][2] << '\t' << knots[0][3] << endl;
//  cout << knots[1][0] << '\t' << knots[1][1] << '\t' << knots[1][2] << '\t' << knots[1][3] << endl;

  adapIntegNodeLocal->SetKnots(knots[0][0], knots[0][1], knots[1][0], knots[1][1], knots[2][0], knots[2][1]);

  //cout << " AAAAAAAAAA " << endl;

  adapIntegNodeLocal->GeomData = GeomData;
  adapIntegNodeLocal->domNums = domNums;

  //cout << " AAAAAAAAAA .... " << refLev << endl;
  
  adapIntegNodeLocal->prepareData();
  adapIntegNodeLocal->subDivide(refLev2);
  
  //cout << " AAAAAAAAAA .... " << refLev << endl;

  GaussQuadrature  QuadratureLocal;

  int mergeFlag=0;

  adapIntegNodeLocal->computeGaussPointsForMerging(0, inclDom, 1, mergeFlag, QuadratureLocal);

  // parametric domain to integration master-quadrilateral domain
  int  ii, gp;
  myPoint  param, geom;
  wt = 0.0;
  ptTemp.setZero();
  param.setZero();
  for(gp=0; gp<QuadratureLocal.gausspoints.size(); gp++)
  {
    for(ii=0; ii<3; ii++)
      param[ii] = QuadratureLocal.gausspoints[gp][ii] ;

    //GeomData->ComputeCoord(param, geom);
    //ptTemp += geom;

    wt +=  QuadratureLocal.gaussweights[gp];
    
    ptTemp += param;
  }

  nGPsMerge = QuadratureLocal.gausspoints.size();

  ptTemp /= nGPsMerge;

  delete  adapIntegNodeLocal;

  return;
}











