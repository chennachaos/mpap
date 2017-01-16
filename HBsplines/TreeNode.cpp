
#include "TreeNode.h"

#include "GeomDataHBSplines.h"
#include "SolutionData.h"




template <>
bool  TreeNode<1>::IsBoundary()
{
  return (IsLeftBoundary() || IsRightBoundary() );
}


template <>
bool  TreeNode<2>::IsBoundary()
{
  return (IsLeftBoundary() || IsRightBoundary() || IsBottomBoundary() || IsTopBoundary());
}

template <>
bool  TreeNode<3>::IsBoundary()
{
  return (IsLeftBoundary() || IsRightBoundary() || IsBottomBoundary() || IsTopBoundary() || IsFrontBoundary() || IsBackBoundary() );
}




template<>
void TreeNode<1>::subDivide()
{
    //assert("Current Cell is a LEAF and Subdivision starts" && IsLeaf());
    
    if( !IsLeaf() )
      return;
    
    NUM_CHILDREN = 2;
    child  = new TreeNode_PTR[NUM_CHILDREN];

    int ii, temp;
    
    temp = level+1;
    
    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
       child[ii] = new TreeNode<1>(temp);

       child[ii]->SetDegree(degree);
       child[ii]->SetParent(this);
    }
    
    if(IsGhost())
    {
       for(ii=0;ii<NUM_CHILDREN;ii++)
         child[ii]->SetGhostOn();
    }
    
    double  mid(0.0);

    mid = 0.5*(knots[0][0] + knots[0][1]);

    TreeNode_PTR  tmpnode;
    
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
}



template<>
void TreeNode<2>::subDivide()
{
    //assert("Current Cell is a LEAF and Subdivision starts" && IsLeaf());

    if( !IsLeaf() )
      return;

    int ii, temp, jj;

    NUM_CHILDREN = 4;
    child      = new TreeNode_PTR[NUM_CHILDREN];

    temp = level+1;
    
    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
       child[ii] = new TreeNode<2>(temp);

       child[ii]->SetDegree(degree);
       
       child[ii]->SetParent(this);
    }
    //
    if(IsGhost())
    {
       for(ii=0;ii<NUM_CHILDREN;ii++)
         child[ii]->SetGhostOn();
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

    TreeNode_PTR  tmpnode;

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
    
    return;
}





template<>
void TreeNode<3>::subDivide()
{
    //assert("Current Cell is a LEAF and Subdivision starts" && IsLeaf());

    if( !IsLeaf() )
      return;

    int ii, temp, jj;

    NUM_CHILDREN = 8;
    child = new TreeNode_PTR[NUM_CHILDREN];

    temp = level+1;
    
    for(ii=0;ii<NUM_CHILDREN;ii++)
    {
       child[ii] = new TreeNode<3>(temp);

       child[ii]->SetDegree(degree);
       child[ii]->SetParent(this);
    }

    if(this->IsGhost())
    {
       for(ii=0;ii<NUM_CHILDREN;ii++)
         child[ii]->SetGhostOn();
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

    TreeNode_PTR  tmpnode, nd;

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

    tmpnode = child[LEFT]->GetNeighbour(LEFT);

    if(tmpnode != NULL)
       tmpnode->SetNeighbour(RIGHT, NULL);

    child[LEFT]->SetNeighbour(LEFT, NULL);
    child[LEFT]->SetNeighbour(RIGHT, NULL);

    // RIGHT side children

    tmpnode = child[RIGHT]->GetNeighbour(RIGHT);

    if(tmpnode != NULL)
       tmpnode->SetNeighbour(LEFT, NULL);

    child[RIGHT]->SetNeighbour(LEFT, NULL);
    child[RIGHT]->SetNeighbour(RIGHT, NULL);

    child[LEFT]->Deactivate();
    child[RIGHT]->Deactivate();

cout << " AAAAAAAAAA " << endl;

//   for(ii=0;ii<NUM_CHILDREN;ii++)
//      child[ii] = NULL;

//   child = NULL;

   return;
}





/*
template<>
bool  TreeNode<1>::IsGhost()
{
    if( (neighbours[RIGHT]) || (neighbours[LEFT] == NULL) )
      return true;
      
    TreeNode_PTR  ndtmp1, ndtmp2;

    ndtmp1 = neighbours[RIGHT];
    ndtmp2 = neighbours[LEFT];

    for(int ii=0;ii<degree[0]-1;ii++)
    {
       ndtmp1 = ndtmp1->GetNeighbour(RIGHT);
       ndtmp2 = ndtmp2->GetNeighbour(LEFT);
       
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
      if(!IsGhost())
      {
            printf("\t   ID          = %5d\n", id);
            printf("\t   Level       = %5d\n", level);
            printf("\t   Degree      = %5d\n", degree[0]);
            printf("\t   Ghost       = %5d\n", IsGhost());
            if(parent == NULL)
              printf("\t   Parent      = %5d\n", -1);
            else
              printf("\t   Parent      = %5d\n", parent->GetID());

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
              printf("\t   Parent      = %5d\n", parent->GetID());

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
              printf("\t   Parent      = %5d\n", parent->GetID());

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












