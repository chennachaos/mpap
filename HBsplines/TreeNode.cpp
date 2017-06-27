
#include "TreeNode.h"

#include "GeomDataHBSplines.h"
#include "SolutionData.h"




template <>
bool  TreeNode<1>::isBoundary()
{
  return (isLeftBoundary() || isRightBoundary() );
}


template <>
bool  TreeNode<2>::isBoundary()
{
  return (isLeftBoundary() || isRightBoundary() || isBottomBoundary() || isTopBoundary());
}

template <>
bool  TreeNode<3>::isBoundary()
{
  return (isLeftBoundary() || isRightBoundary() || isBottomBoundary() || isTopBoundary() || isFrontBoundary() || isBackBoundary() );
}




template<>
void TreeNode<1>::subDivide()
{
    //assert("Current Cell is a LEAF and Subdivision starts" && isLeaf());
    
    if( !isLeaf() )
      return;
    
    NUM_CHILDREN = 2;
    child  = new TreeNode_PTR[NUM_CHILDREN];

    int ii, temp;
    
    temp = level+1;
    
    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
       child[ii] = new TreeNode<1>(temp);

       child[ii]->setDegree(degree);
       child[ii]->setParent(this);
    }
    
    if(isGhost())
    {
       for(ii=0;ii<NUM_CHILDREN;ii++)
         child[ii]->setGhostOn();
    }
    
    double  mid(0.0);

    mid = 0.5*(knots[0][0] + knots[0][1]);

    TreeNode_PTR  tmpnode;
    
    child[LEFT]->orientation = LEFT;
    child[RIGHT]->orientation = RIGHT;
    
    //cout << " AAAAAAAAAAAAAA " << endl;

    child[LEFT]->setKnots(Dir1, knots[0][0], mid);
    child[RIGHT]->setKnots(Dir1, mid, knots[0][1]);

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
}



template<>
void TreeNode<2>::subDivide()
{
    //assert("Current Cell is a LEAF and Subdivision starts" && isLeaf());

    if( !isLeaf() )
      return;

    int ii, temp, jj;

    NUM_CHILDREN = 4;
    child      = new TreeNode_PTR[NUM_CHILDREN];

    temp = level+1;
    
    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
       child[ii] = new TreeNode<2>(temp);

       child[ii]->setDegree(degree);
       
       child[ii]->setParent(this);
    }
    //
    if(isGhost())
    {
       for(ii=0;ii<NUM_CHILDREN;ii++)
         child[ii]->setGhostOn();
    }
    //

    child[SW]->orientation = SW;
    child[SE]->orientation = SE;
    child[NW]->orientation = NW;
    child[NE]->orientation = NE;

    //cout << " child[SW]->orientation " << '\t' << SW << '\t' << child[SW]->orientation << endl;
    //cout << " child[SE]->orientation " << '\t' << SE << '\t' << child[SE]->orientation << endl;
    //cout << " child[NW]->orientation " << '\t' << NW << '\t' << child[NW]->orientation << endl;
    //cout << " child[NE]->orientation " << '\t' << NE << '\t' << child[NE]->orientation << endl;

    double  midu(0.0), midv(0.0);
    
    midu = 0.5*(knots[0][1] + knots[0][0]);
    midv = 0.5*(knots[1][1] + knots[1][0]);

    child[SW]->setKnots(Dir1, knots[0][0], midu);
    child[SW]->setKnots(Dir2, knots[1][0], midv);

    child[SE]->setKnots(Dir1, midu,        knots[0][1]);
    child[SE]->setKnots(Dir2, knots[1][0], midv);

    child[NW]->setKnots(Dir1, knots[0][0], midu);
    child[NW]->setKnots(Dir2, midv,        knots[1][1]);

    child[NE]->setKnots(Dir1, midu,        knots[0][1]);
    child[NE]->setKnots(Dir2, midv,        knots[1][1]);


    child[SW]->setNeighbour(EAST,  child[SE]);
    child[SW]->setNeighbour(NORTH, child[NW]);
    child[SE]->setNeighbour(WEST,  child[SW]);
    child[SE]->setNeighbour(NORTH, child[NE]);
    child[NW]->setNeighbour(EAST,  child[NE]);
    child[NW]->setNeighbour(SOUTH, child[SW]);
    child[NE]->setNeighbour(WEST,  child[NW]);
    child[NE]->setNeighbour(SOUTH, child[SE]);

    TreeNode_PTR  tmpnode;

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
    
    return;
}





template<>
void TreeNode<3>::subDivide()
{
    //assert("Current Cell is a LEAF and Subdivision starts" && isLeaf());

    if( !isLeaf() )
      return;

    int ii, temp, jj;

    NUM_CHILDREN = 8;
    child = new TreeNode_PTR[NUM_CHILDREN];

    temp = level+1;
    
    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
       child[ii] = new TreeNode<3>(temp);

       child[ii]->setDegree(degree);
       child[ii]->setParent(this);
    }

    if(this->isGhost())
    {
       for(ii=0;ii<NUM_CHILDREN;ii++)
         child[ii]->setGhostOn();
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

    child[SW_BACK ]->setKnots(Dir1, knots[0][0], midu);
    child[SW_FRONT]->setKnots(Dir1, knots[0][0], midu);
    child[NW_BACK ]->setKnots(Dir1, knots[0][0], midu);
    child[NW_FRONT]->setKnots(Dir1, knots[0][0], midu);

    child[SE_BACK ]->setKnots(Dir1, midu,        knots[0][1]);
    child[SE_FRONT]->setKnots(Dir1, midu,        knots[0][1]);
    child[NE_BACK ]->setKnots(Dir1, midu,        knots[0][1]);
    child[NE_FRONT]->setKnots(Dir1, midu,        knots[0][1]);


    child[SW_BACK ]->setKnots(Dir2, knots[1][0], midv);
    child[SW_FRONT]->setKnots(Dir2, knots[1][0], midv);
    child[SE_BACK ]->setKnots(Dir2, knots[1][0], midv);
    child[SE_FRONT]->setKnots(Dir2, knots[1][0], midv);

    child[NW_BACK ]->setKnots(Dir2, midv,        knots[1][1]);
    child[NW_FRONT]->setKnots(Dir2, midv,        knots[1][1]);
    child[NE_BACK ]->setKnots(Dir2, midv,        knots[1][1]);
    child[NE_FRONT]->setKnots(Dir2, midv,        knots[1][1]);


    child[SW_BACK ]->setKnots(Dir3, knots[2][0], midw);
    child[SE_BACK ]->setKnots(Dir3, knots[2][0], midw);
    child[NW_BACK ]->setKnots(Dir3, knots[2][0], midw);
    child[NE_BACK ]->setKnots(Dir3, knots[2][0], midw);

    child[SW_FRONT]->setKnots(Dir3, midw,        knots[2][1]);
    child[SE_FRONT]->setKnots(Dir3, midw,        knots[2][1]);
    child[NW_FRONT]->setKnots(Dir3, midw,        knots[2][1]);
    child[NE_FRONT]->setKnots(Dir3, midw,        knots[2][1]);


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

    TreeNode_PTR  tmpnode, nd;

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



    return;
}



template<>
void TreeNode<1>::unRefine()
{
   if(child == NULL)
   {
      cout << " This element has no children. So can't unrefine it ... " << endl;
      return;
   }
   
   int  ii, jj;

    TreeNode_PTR  tmpnode, nd;

    // LEFT side children

    tmpnode = child[LEFT]->getNeighbour(LEFT);

    if(tmpnode != NULL)
       tmpnode->setNeighbour(RIGHT, NULL);

    child[LEFT]->setNeighbour(LEFT, NULL);
    child[LEFT]->setNeighbour(RIGHT, NULL);

    // RIGHT side children

    tmpnode = child[RIGHT]->getNeighbour(RIGHT);

    if(tmpnode != NULL)
       tmpnode->setNeighbour(LEFT, NULL);

    child[RIGHT]->setNeighbour(LEFT, NULL);
    child[RIGHT]->setNeighbour(RIGHT, NULL);

    child[LEFT]->deactivate();
    child[RIGHT]->deactivate();

cout << " AAAAAAAAAA " << endl;

//   for(ii=0;ii<NUM_CHILDREN;ii++)
//      child[ii] = NULL;

//   child = NULL;

   return;
}





/*
template<>
bool  TreeNode<1>::isGhost()
{
    if( (neighbours[RIGHT]) || (neighbours[LEFT] == NULL) )
      return true;
      
    TreeNode_PTR  ndtmp1, ndtmp2;

    ndtmp1 = neighbours[RIGHT];
    ndtmp2 = neighbours[LEFT];

    for(int ii=0;ii<degree[0]-1;ii++)
    {
       ndtmp1 = ndtmp1->getNeighbour(RIGHT);
       ndtmp2 = ndtmp2->getNeighbour(LEFT);
       
       if( (ndtmp1 == NULL) || (ndtmp2 == NULL) )
       {
         break;
         return true;
       }
    }

    return false;
}
*/


template<>
void TreeNode<1>::printSelf()
{
      if(!isGhost())
      {
            printf("\t   ID          = %5d\n", id);
            printf("\t   Level       = %5d\n", level);
            printf("\t   Degree      = %5d\n", degree[0]);
            printf("\t   Ghost       = %5d\n", isGhost());
            if(parent == NULL)
              printf("\t   Parent      = %5d\n", -1);
            else
              printf("\t   Parent      = %5d\n", parent->getID());

            printf("\t   children    = %5d\n", NUM_CHILDREN);
            printf("\t   Basis Functions --->  \n");
            if( !LocalBasisFuncs.empty() )
            {
               printf("\t\t");
               for(int jj=0;jj<LocalBasisFuncs.size();jj++)
                 printf("\t%5d\t", LocalBasisFuncs[jj]);
               printf("\n");
            }
            else
              printf(" NONE \n");
            printf("\t   Parameters ... \n");
            printf("\t\t direction #%5d ---> \t%12.6f \t %12.6f\n", 1, knots[0][0], knots[0][1]);

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

    }
}






template<>
void TreeNode<2>::printSelf()
{
            printf("\t   ID          = %5d\n", id);
            printf("\t   Level       = %5d\n", level);
            printf("\t   Degree      = ");
            for(int ii=0;ii<2;ii++)
              printf("%5d\t", degree[ii]);
            printf("\n");
            if(parent == NULL)
              printf("\t   Parent      = %5d\n", -1);
            else
              printf("\t   Parent      = %5d\n", parent->getID());

            printf("\t   Basis Functions --->  \n");
            if( !LocalBasisFuncs.empty() )
            {
               printf("\t\t");
               for(int jj=0;jj<LocalBasisFuncs.size();jj++)
                 printf("\t%5d\t", LocalBasisFuncs[jj]);
               printf("\n");
            }
            else
              printf(" NONE \n");
            printf("\t   Parameters ... \n");
            for(int ii=0;ii<2;ii++)
            {
               printf("\t\t direction #%5d ---> \t%12.6f \t %12.6f\n", (ii+1), knots[ii][0], knots[ii][1]);
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
void TreeNode<3>::printSelf()
{
            printf("\t   ID          = %5d\n", id);
            printf("\t   Level       = %5d\n", level);
            printf("\t   Degree      = ");
            for(int ii=0;ii<3;ii++)
              printf("%5d\t", degree[ii]);
            printf("\n");
            if(parent == NULL)
              printf("\t   Parent      = %5d\n", -1);
            else
              printf("\t   Parent      = %5d\n", parent->getID());

            printf("\t   Basis Functions --->  \n");
            if( !LocalBasisFuncs.empty() )
            {
               printf("\t\t");
               for(int jj=0;jj<LocalBasisFuncs.size();jj++)
                 printf("\t%5d\t", LocalBasisFuncs[jj]);
               printf("\n");
            }
            else
              printf(" NONE \n");
            printf("\t   Parameters ... \n");
            for(int ii=0;ii<3;ii++)
            {
               printf("\t\t direction #%5d ---> \t%12.6f \t %12.6f\n", (ii+1), knots[ii][0], knots[ii][1]);
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
bool  TreeNode<1>::pointLiesInside(const myPoint& pt)
{
  if( (pt[0] >= knots[0][0]) && (pt[0] <= knots[0][1]) )
    return true;
  else
    return false;
}


template<>
bool  TreeNode<2>::pointLiesInside(const myPoint& pt)
{
  if( (pt[0] >= knots[0][0]) && (pt[0] <= knots[0][1]) )
  {
    if( (pt[1] >= knots[1][0]) && (pt[1] <= knots[1][1]) )
      return true;
    else
      return false;
  }
  else
    return false;
}



template<>
bool  TreeNode<3>::pointLiesInside(const myPoint& pt)
{
  if( (pt[0] >= knots[0][0]) && (pt[0] <= knots[0][1]) )
  {
    if( (pt[1] >= knots[1][0]) && (pt[1] <= knots[1][1]) )
    {
      if( (pt[2] >= knots[2][0]) && (pt[2] <= knots[2][1]) )
        return true;
      else
        return false;
    }
    else
      return false;
  }
  else
    return false;
}












