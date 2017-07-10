
#include "HBSplineBase.h"



void  HBSplineBase::subDivide(int nodenum)
{
    int ii, jj, ll;

    node  *nd, *nd2;

    nd = elems[nodenum];
    if( nd->isLeaf() )
    {
      nd->subDivide();

      ll = nd->getLevel() + 1;

      for(jj=0;jj<nd->getNumberOfChildren();jj++)
      {
         nd2 = nd->getChild(jj);

         nd2->SolnData = &(SolnData);
         nd2->GeomData = &(GeomData);

         elems.push_back(nd2);
         NodeNumsAtLevel[ll].push_back(nd2->getID());
      }
    }

    return;
}



void  HBSplineBase::applyRefinementProcess()
{
    //printf(" Current level = %2d \n", CURRENT_LEVEL);
    //printf(" Maximum level = %2d \n\n", MAX_LEVEL);

    algorithm1(CURRENT_LEVEL);
    algorithm2(CURRENT_LEVEL);
    algorithm3(CURRENT_LEVEL+1);

    return;
}



void  HBSplineBase::algorithm1(int lev)
{
    // nucleus operation for each element: create new leaves (children) of level k+1
    //
    //////////////////////////////////////////////////////////////////////////////////////
   
    //printf(" algorithm1 for level = %2d STARTED \n ",lev);

    if(DIM == 1)
      algorithm1_1D(lev);
    else if(DIM == 2)
      algorithm1_2D(lev);
    else
      algorithm1_3D(lev);


    findUnique(nodes2divide);

    int ii;
    /*
    printf(" nodes2divide  \n");
    for(ii=0;ii<nodes2divide.size();ii++)
       printf("%5d \t", nodes2divide[ii] );
    printf("\n\n\n");
    */

    for(ii=0;ii<nodes2divide.size();ii++)
      subDivide(nodes2divide[ii]);
    
    //findUnique(boundaryNodes);

    /*
    printf(" boundaryNodes .. : \n");
    for(int jj=0;jj<boundaryNodes.size();jj++)
       cout << '\t' << boundaryNodes[jj] ;
    printf("\n\n");
    */

    //printf(" algorithm1 for level = %2d FINISHED \n ",lev);

    return;
}



void  HBSplineBase::algorithm2(int lev)
{
    // remove linear dependence between basis functions of level k and level k+1.
    
    //printf("\n\n\n  algorithm2 for level = %2d STARTED \n ",lev);

    /*
    printf("\n\n\n");
    printf(" Nodes at level = %5d \n", lev);
    for(int ii=0;ii<NodeNumsAtLevel[lev].size();ii++)
       cout << '\t' << NodeNumsAtLevel[lev][ii] ;
    printf("\n\n\n");
    */

    if(DIM == 1)
      algorithm2_1D(lev);
    else if(DIM == 2)
      algorithm2_2D(lev);
    else
      algorithm2_3D(lev);

    findUnique(VacantBFs);

    /*
    printf("\n\n\n");
    printf(" Vacant Basis Functions \n");
    for(int ii=0;ii<VacantBFs.size();ii++)
       cout << '\t' << VacantBFs[ii] ;
    printf("\n\n\n");
    */

    //printf(" algorithm2 for level = %2d FINISHED \n", lev);

    return;
}



void  HBSplineBase::algorithm3(int lev)
{
    // assign numbers to the new basis functions of level k+1
    
    //printf("\n\n\n algorithm3 for level = %2d STARTED \n ",lev);

    if(DIM == 1)
      algorithm3_1D(lev);
    else if(DIM == 2)
      algorithm3_2D(lev);
    else
      algorithm3_3D(lev);

    //printf("\n \t   Total number of basis functions    =  %5d\n\n", gridBF1);

    //printf(" algorithm3 for level = %2d FINISHED \n ",lev);

    return;
}




void  HBSplineBase::addGhostNodes(node* nd, int dir)
{
    nd = nd->getNeighbour(dir);

    if( nd != NULL )
    {
      if(nd->isGhost())
      {
          nodes2divide.push_back(nd->getID());
          nd = nd->getNeighbour(dir);
         
          while( nd != NULL)
          {
             nodes2divide.push_back(nd->getID());
             nd = nd->getNeighbour(dir);
          }
      }
    }

    return;
}


void  HBSplineBase::algorithm1_1D(int lev)
{
    int ee, nodenum;
    node  *nd;

    printVector(elemsToRefine);
    nodes2divide.clear();
    for(ee=0;ee<elemsToRefine.size();ee++)
    {
      //cout << " ee " << ee << '\t' << elemsToRefine[ee] << endl;
      nodenum = elemsToRefine[ee];

      nodes2divide.push_back(nodenum);

      // if the current node is a ghost node then take all of its ghost neighbours
      ////////////////////////////////////////////////////////

      nd = elems[nodenum];
      //nd->printSelf();
      //cout << " AAAAAAAAA " << endl;
      if(nd->isLeftBoundary())
        addGhostNodes(nd, LEFT);

      //cout << " AAAAAAAAA " << endl;
      nd = elems[nodenum];
      if(nd->isRightBoundary())
        addGhostNodes(nd, RIGHT);
      //cout << " ee " << ee << '\t' << elemsToRefine[ee] << endl;
    }

    return;
}


void  HBSplineBase::algorithm1_2D(int lev)
{
    int  ee, nodenum;
    node  *nd;

    nodes2divide.clear();
    for(ee=0;ee<elemsToRefine.size();ee++)
    {
       nodenum = elemsToRefine[ee];

       nodes2divide.push_back(nodenum);
       
       //
       // if the current node is a boundary/ghost node then take all the ghost neighbours of it
       ////////////////////////////////////////////////////////
       
       nd = elems[nodenum];
       if(nd->isLeftBoundary())
         addGhostNodes(nd, WEST);
       
       nd = elems[nodenum];
       if(nd->isRightBoundary())
         addGhostNodes(nd, EAST);
       
       nd = elems[nodenum];
       if(nd->isTopBoundary())
         addGhostNodes(nd, NORTH);

       nd = elems[nodenum];
       if(nd->isBottomBoundary())
         addGhostNodes(nd, SOUTH);


       // bottom-left corner
       nd = elems[nodenum];
       if(nd->isLeftBoundary() && nd->isBottomBoundary())
       {
          while(nd != NULL)
          {
            addGhostNodes(nd, WEST);
            nd = nd->getNeighbour(SOUTH);
          }
       }
       
       // bottom-right corner
       nd = elems[nodenum];
       if(nd->isRightBoundary() && nd->isBottomBoundary())
       {
          while(nd != NULL)
          {
            addGhostNodes(nd, EAST);
            nd = nd->getNeighbour(SOUTH);
          }
       }

       // top-right corner
       nd = elems[nodenum];
       if(nd->isRightBoundary() && nd->isTopBoundary())
       {
          while(nd != NULL)
          {
            addGhostNodes(nd, EAST);
            nd = nd->getNeighbour(NORTH);
          }
       }

       // top-left corner
       nd = elems[nodenum];
       if(nd->isLeftBoundary() && nd->isTopBoundary())
       {
          while(nd != NULL)
          {
            addGhostNodes(nd, WEST);
            nd = nd->getNeighbour(NORTH);
          }
       }
    }

    return;
}


void  HBSplineBase::algorithm1_3D(int lev)
{
  // still pending
    
    int  ee, nodenum;
    node  *nd, *nd1;

    nodes2divide.clear();
    for(ee=0;ee<elemsToRefine.size();ee++)
    {
       nodenum = elemsToRefine[ee];

      nodes2divide.push_back(nodenum);

      //
      // if the current node is a boundary/ghost node then take all the ghost neighbours of it
      ////////////////////////////////////////////////////////

      nd = elems[nodenum];
      if(nd->isLeftBoundary())
        addGhostNodes(nd, WEST);

      nd = elems[nodenum];
      if(nd->isRightBoundary())
        addGhostNodes(nd, EAST);

      nd = elems[nodenum];
      if(nd->isTopBoundary())
        addGhostNodes(nd, NORTH);

      nd = elems[nodenum];
      if(nd->isBottomBoundary())
        addGhostNodes(nd, SOUTH);

      nd = elems[nodenum];
      if(nd->isBackBoundary())
        addGhostNodes(nd, BACK);

      nd = elems[nodenum];
      if(nd->isFrontBoundary())
        addGhostNodes(nd, FRONT);

      // front-left edge
      nd = elems[nodenum];
      if(nd->isFrontBoundary() && nd->isLeftBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, FRONT);
          nd = nd->getNeighbour(WEST);
        }
      }
       
      // front-right edge
      nd = elems[nodenum];
      if(nd->isFrontBoundary() && nd->isRightBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, FRONT);
          nd = nd->getNeighbour(EAST);
        }
      }
       
      // front-bottom edge
      nd = elems[nodenum];
      if(nd->isFrontBoundary() && nd->isBottomBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, FRONT);
          nd = nd->getNeighbour(SOUTH);
        }
      }

      // front-top edge
      nd = elems[nodenum];
      if(nd->isFrontBoundary() && nd->isTopBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, FRONT);
          nd = nd->getNeighbour(NORTH);
        }
      }

      // back-left edge
      nd = elems[nodenum];
      if(nd->isBackBoundary() && nd->isLeftBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, BACK);
          nd = nd->getNeighbour(WEST);
        }
      }
       
      // back-right edge
      nd = elems[nodenum];
      if(nd->isBackBoundary() && nd->isRightBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, BACK);
          nd = nd->getNeighbour(EAST);
        }
      }

      // back-bottom edge
      nd = elems[nodenum];
      if(nd->isBackBoundary() && nd->isBottomBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, BACK);
          nd = nd->getNeighbour(SOUTH);
        }
      }

      // back-top edge
      nd = elems[nodenum];
      if(nd->isBackBoundary() && nd->isTopBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, BACK);
          nd = nd->getNeighbour(NORTH);
        }
      }

      //cout << " AAAAAAAAA " << endl;
      // front-left-bottom corner
      nd = elems[nodenum];
      if(nd->isFrontBoundary() && nd->isLeftBoundary() && nd->isBottomBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, SOUTH);
            nd1 = nd1->getNeighbour(WEST);
          }
          nd = nd->getNeighbour(FRONT);
        }
      }

      // front-left-top corner
      nd = elems[nodenum];
      if(nd->isFrontBoundary() && nd->isLeftBoundary() && nd->isTopBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, NORTH);
            nd1 = nd1->getNeighbour(WEST);
          }
          nd = nd->getNeighbour(FRONT);
        }
      }

      // front-right-bottom corner
      nd = elems[nodenum];
      if(nd->isFrontBoundary() && nd->isRightBoundary() && nd->isBottomBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, SOUTH);
            nd1 = nd1->getNeighbour(EAST);
          }
          nd = nd->getNeighbour(FRONT);
        }
      }

      // front-right-top corner
      nd = elems[nodenum];
      if(nd->isFrontBoundary() && nd->isRightBoundary() && nd->isTopBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, NORTH);
            nd1 = nd1->getNeighbour(EAST);
          }
          nd = nd->getNeighbour(FRONT);
        }
      }

      // back-left-bottom corner
      nd = elems[nodenum];
      if(nd->isBackBoundary() && nd->isLeftBoundary() && nd->isBottomBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, SOUTH);
            nd1 = nd1->getNeighbour(WEST);
          }
          nd = nd->getNeighbour(BACK);
        }
      }

      // back-left-top corner
      nd = elems[nodenum];
      if(nd->isBackBoundary() && nd->isLeftBoundary() && nd->isTopBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, NORTH);
            nd1 = nd1->getNeighbour(WEST);
          }
          nd = nd->getNeighbour(BACK);
        }
      }

      // back-right-bottom corner
      nd = elems[nodenum];
      if(nd->isBackBoundary() && nd->isRightBoundary() && nd->isBottomBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, SOUTH);
            nd1 = nd1->getNeighbour(EAST);
          }
          nd = nd->getNeighbour(BACK);
        }
      }

      // back-right-top corner
      nd = elems[nodenum];
      if(nd->isBackBoundary() && nd->isRightBoundary() && nd->isTopBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, NORTH);
            nd1 = nd1->getNeighbour(EAST);
          }
          nd = nd->getNeighbour(BACK);
        }
      }
    }

    return;
}



void  HBSplineBase::algorithm2_1D(int lev)
{
    int ee, ii, jj, val, temp1, temp2;
    bool flag1, flag2;

    node  *nd, *nd1;

    for(ee=0;ee<NodeNumsAtLevel[lev].size();ee++)
    {
       nd = elems[NodeNumsAtLevel[lev][ee]];
      
       // if all the (degree[0]+1) elements are 
       // ghost elements (or NULL) then do not remove any number
      
       flag1 = true;
       
       if(nd->isGhost())
         flag2 = false;
       else
         flag2 = true;
      
       nd1 = nd;
       for(ii=0;ii<=degree[0];ii++)
       {
          if( nd1 == NULL || nd1->isLeaf() )
          {
             flag1 = false;
             break;
          }
          
          if(!nd1->isGhost())
            flag2 = true;

          nd1 = nd1->getNeighbour(RIGHT);
       }

       //cout << nd->getID() << '\t' << flag1 << '\t' << flag2 << endl;

       if( flag1 && flag2 )
       {
          temp1 = (degree[0]+1) - 1;
          temp2 = 0;

          nd1 = nd;
          val = nd1->LocalBasisFuncs[temp1];

          if(val != -1)
          {
             VacantBFs.push_back(val);
             for(ii=0;ii<=degree[0];ii++)
             {
                nd1->LocalBasisFuncs[temp1-temp2++] = -1;
                nd1 = nd1->getNeighbour(RIGHT);
             }
          }
       }
    }

    return;
}


void  HBSplineBase::algorithm2_2D(int lev)
{
    int ii, jj, val, ee, temp1, temp2;
    bool flag1, flag2;

    node  *nd, *nd1, *nd2;

    //cout << " Node # " << nd1->getID() << '\t' << nd1->isGhost() << endl;

    for(ee=0;ee<NodeNumsAtLevel[lev].size();ee++)
    {
       nd = elems[NodeNumsAtLevel[lev][ee]];
      
       flag1 = true;
      
       // if the current node is a ghost node then check if all the neighbours are ghost nodes
       // if all the neighbours are ghost nodes then don't remove any basis function

       if(nd->isGhost())
         flag2 = false;
       else
         flag2 = true;

       nd2 = nd;
       for(jj=0;jj<=degree[1];jj++)
       {
         nd1 = nd2;
         for(ii=0;ii<=degree[0];ii++)
         {
            if( nd1 == NULL || nd1->isLeaf() )
            {
              flag1 = false;
              break;
            }

            if( !nd1->isGhost() )
              flag2 = true;

            nd1 = nd1->getNeighbour(EAST);
         }
         if(!flag1)
           break;
         
         nd2 = nd2->getNeighbour(NORTH);
      }

      //cout << " nd->getID() ... : "  << '\t'<< nd->getID() << '\t' << nd->isGhost() << '\t' << flag1 << '\t' << flag2 << endl;

      if( flag1 && flag2 )
      {
         temp1 = (degree[0]+1)*(degree[1]+1) - 1;
         temp2 = 0;
         nd2 = nd;

         val = nd2->LocalBasisFuncs[temp1];

         if(val != -1)
         {
            VacantBFs.push_back(val);

            for(jj=0;jj<=degree[1];jj++)
            {
               nd1 = nd2;
               for(ii=0;ii<=degree[0];ii++)
               {
                  nd1->LocalBasisFuncs[temp1-temp2++] = -1;
                  nd1 = nd1->getNeighbour(EAST);
               }
               nd2 = nd2->getNeighbour(NORTH);
            }
         }
      }
   }

   return;
}


void  HBSplineBase::algorithm2_3D(int lev)
{
    int ii, jj, kk, val, ee, temp1, temp2;
    bool flag1, flag2;

    node  *nd, *nd1, *nd2, *nd3;

    //cout << " Node # " << nd1->getID() << '\t' << nd1->isGhost() << endl;

    for(ee=0;ee<NodeNumsAtLevel[lev].size();ee++)
    {
       nd = elems[NodeNumsAtLevel[lev][ee]];
      
       flag1 = true;
      
       // if the current node is a ghost node then check if all the neighbours are ghost nodes
       // if all the neighbours are ghost nodes then don't remove any basis function

       if(nd->isGhost())
         flag2 = false;
       else
         flag2 = true;

       nd3 = nd;
       for(kk=0;kk<=degree[2];kk++)
       {
         nd2 = nd3;
         for(jj=0;jj<=degree[1];jj++)
         {
           nd1 = nd2;
           for(ii=0;ii<=degree[0];ii++)
           {
             if( nd1 == NULL || nd1->isLeaf() )
             {
               flag1 = false;
               break;
             }

             if( !nd1->isGhost() )
               flag2 = true;

             nd1 = nd1->getNeighbour(EAST);
           }
           if(!flag1)
             break;
         
           nd2 = nd2->getNeighbour(NORTH);
        }
        if(!flag1)
          break;

        nd3 = nd3->getNeighbour(FRONT);
      }

      //cout << " nd->getID() ... : "  << '\t'<< nd->getID() << '\t' << nd->isGhost() << '\t' << flag1 << '\t' << flag2 << endl;

      if( flag1 && flag2 )
      {
         temp1 = (degree[0]+1)*(degree[1]+1)*(degree[2]+1) - 1;
         temp2 = 0;
         nd3 = nd;

         val = nd3->LocalBasisFuncs[temp1];

         if(val != -1)
         {
            VacantBFs.push_back(val);

            for(kk=0;kk<=degree[2];kk++)
            {
              nd2 = nd3;
              for(jj=0;jj<=degree[1];jj++)
              {
                nd1 = nd2;
                for(ii=0;ii<=degree[0];ii++)
                {
                  nd1->LocalBasisFuncs[temp1-temp2++] = -1;
                  nd1 = nd1->getNeighbour(EAST);
                }
                nd2 = nd2->getNeighbour(NORTH);
              }
              nd3 = nd3->getNeighbour(FRONT);
            }
         }
      }
   }

    return;
}


void  HBSplineBase::algorithm3_1D(int lev)
{
    int  ee, ii, jj, val;
    bool flag1, flag2;

    node  *nd, *nd1;

    // flag1 is to check if any of the (degree[0]+1) RIGHT neighbours is NULL.
    //
    // If flag1 is 'TRUE' then NONE of the (degree[0]+1) RIGHT neighbours is NULL and a basis function can be assigned.
    // If flag1 is 'FALSE' then ONE of the (degree[0]+1) RIGHT neighbours is NULL and a basis function can not be assigned.


    // flag2 is to check if all of the (degree[0]+1) RIGHT neighbours are GHOST cells.
    //
    // If flag2 is 'TRUE' then ALL of the (degree[0]+1) RIGHT neighbours are not GHOST cells and a basis function can be assigned.
    // Otherwise, a basis function can not be assigned.
    
    // a basis function can be assigned only when both the flags are true

    for(ee=0;ee<NodeNumsAtLevel[lev].size();ee++)
    {
       nd = elems[NodeNumsAtLevel[lev][ee]];
      
       // if all the (degree[0]+1) elements are 
       // ghost elements (or NULL) then do not assign any number

       flag1 = true;
       
       if(nd->isGhost())
         flag2 = false;
       else
         flag2 = true;
      
       nd1 = nd;
       for(jj=0;jj<=degree[0];jj++)
       {
         if( nd1 == NULL )
         {
           flag1 = false;
           break;
         }

         if(!nd1->isGhost())
          flag2 = true;

         nd1 = nd1->getNeighbour(RIGHT);
       }

       //cout << nd->getID() << '\t' << flag1 << '\t' << flag2 << endl;

       if( flag1 && flag2 )
       {
         nd1 = nd;
         if(VacantBFs.empty())
           val = gridBF1++;
         else
         {
           val = VacantBFs[0];
           VacantBFs.erase(VacantBFs.begin());
         }

         for(jj=0;jj<=degree[0];jj++)
         {
           nd1->LocalBasisFuncs[degree[0]-jj] = val;
           nd1 = nd1->getNeighbour(RIGHT);
         }
       }
    }

   return;
}



void  HBSplineBase::algorithm3_2D(int lev)
{
    int ii, jj, val, ee, temp1, temp2;
    bool flag1, flag2;

    node  *nd, *nd1, *nd2;

    for(ee=0;ee<NodeNumsAtLevel[lev].size();ee++)
    {
       nd = elems[NodeNumsAtLevel[lev][ee]];
      
       flag1 = true;

       if(nd->isGhost())
         flag2 = false;
       else
         flag2 = true;

       // if the current node is a ghost node then check if all the neighbours are ghost nodes
       // if all the neighbours are ghost nodes then don't add a new basis function

       nd2 = nd;
       for(jj=0;jj<=degree[1];jj++)
       {
          nd1 = nd2;
          for(ii=0;ii<=degree[0];ii++)
          {
             if( nd1 == NULL )
             {
               flag1 = false;
               break;
             }

             if( !nd1->isGhost() )
               flag2 = true;

             nd1 = nd1->getNeighbour(EAST);
          }
          if(!flag1)
            break;
         
          nd2 = nd2->getNeighbour(NORTH);
       }
      
       //cout << " nd->getID() ... : "  << '\t'<< nd->getID() << '\t' << nd->isGhost() << '\t' << flag1 << '\t' << flag2 << endl;

       if( flag1 && flag2 )
       {
          if(VacantBFs.empty())
            val = gridBF1++;
          else
          {
            val = VacantBFs[0];
            VacantBFs.erase(VacantBFs.begin());
          }
         
          temp1 = (degree[0]+1)*(degree[1]+1)-1;
          temp2 = 0;
          nd2 = nd;
          for(jj=0;jj<=degree[1];jj++)
          {
             nd1 = nd2;
             for(ii=0;ii<=degree[0];ii++)
             {
                nd1->LocalBasisFuncs[temp1-temp2++] = val;
                nd1 = nd1->getNeighbour(EAST);
             }
             nd2 = nd2->getNeighbour(NORTH);
          }
       }
    }
   
    return;
}




void  HBSplineBase::algorithm3_3D(int lev)
{
    int ii, jj, kk, val, ee, temp1, temp2;
    
    bool flag1, flag2;

    node  *nd, *nd1, *nd2, *nd3;

    //for(ee=0;ee<elems.size();ee++)
      //elems[ee]->printSelf();
    
    //printVector(NodeNumsAtLevel[lev]);

    for(ee=0;ee<NodeNumsAtLevel[lev].size();ee++)
    {
       //cout << " ee = " << ee << endl;
       nd = elems[NodeNumsAtLevel[lev][ee]];
       //nd->printSelf();
      
       flag1 = true;

       if(nd->isGhost())
         flag2 = false;
       else
         flag2 = true;

       // if the current node is a ghost node then check if all the neighbours are ghost nodes
       // if all the neighbours are ghost nodes then don't add a new basis function

       nd3 = nd;
       for(kk=0;kk<=degree[2];kk++)
       {
         nd2 = nd3;
         for(jj=0;jj<=degree[1];jj++)
         {
           nd1 = nd2;
           for(ii=0;ii<=degree[0];ii++)
           {
             if( nd1 == NULL )
             {
               flag1 = false;
               break;
             }
             if( !nd1->isGhost() )
               flag2 = true;

             nd1 = nd1->getNeighbour(EAST);
           }
           if(!flag1)
             break;
         
           nd2 = nd2->getNeighbour(NORTH);
         }
         if(!flag1)
           break;

         nd3 = nd3->getNeighbour(FRONT);
       }
      
       //cout << " nd->getID() ... : "  << '\t'<< nd->getID() << '\t' << nd->isGhost() << '\t' << flag1 << '\t' << flag2 << endl;

       if( flag1 && flag2 )
       {
          if(VacantBFs.empty())
            val = gridBF1++;
          else
          {
            val = VacantBFs[0];
            VacantBFs.erase(VacantBFs.begin());
          }
          //cout << " val " << val << endl;
         
          temp1 = (degree[0]+1)*(degree[1]+1)*(degree[2]+1)-1;
          temp2 = 0;
          nd3 = nd;
          for(kk=0;kk<=degree[2];kk++)
          {
            nd2 = nd3;
            for(jj=0;jj<=degree[1];jj++)
            {
              nd1 = nd2;
              for(ii=0;ii<=degree[0];ii++)
              {
                nd1->LocalBasisFuncs[temp1-temp2++] = val;
                nd1 = nd1->getNeighbour(EAST);
              }
              nd2 = nd2->getNeighbour(NORTH);
            }
            nd3 = nd3->getNeighbour(FRONT);
          }
       }
    }

    return;
}








void  HBSplineBase::unRefine(int nodenum)
{

   vector<int>  nodes2unrefine;
   
   int ii, jj, ll, count = 0;

   printf(" \n \n unrefining .............. \n \n ");

   node  *nd, *nd2;
   nd = elems[nodenum];

   nodes2unrefine.push_back(nodenum);

   nd = nd->getNeighbour(LEFT);

      if(nd->isGhost())
      {
         // take all ghost nodes to the left of this node

         nd = nd->getNeighbour(LEFT);
         while( nd != NULL)
         {
            nodes2unrefine.push_back(nd->getID());
            nd = nd->getNeighbour(LEFT);
         }
      }

   nd = elems[nodenum];

   nd = nd->getNeighbour(RIGHT);

      if(nd->isGhost())
      {
         // take all ghost nodes to the left of this node

         nd = nd->getNeighbour(RIGHT);
         while( nd != NULL)
         {
            nodes2unrefine.push_back(nd->getID());
            nd = nd->getNeighbour(RIGHT);
         }
      }

   sort(nodes2unrefine.begin(), nodes2unrefine.end());
   nodes2unrefine.erase(unique(nodes2unrefine.begin(), nodes2unrefine.end()), nodes2unrefine.end());

   printf(" nodes2unrefine  ");
   for(ii=0;ii<nodes2unrefine.size();ii++)
      printf("\t %5d", nodes2unrefine[ii] );
   printf("\n\n\n");

   for(ii=0;ii<nodes2unrefine.size();ii++)
   {
      nd = elems[nodes2unrefine[ii]];

      nd->printSelf();
      nd->unRefine();
      nd->activate();
      /*
      for(jj=0;jj<2;jj++)
      {
         cout << "  nd->getChild(jj)->getID()  " << nd->getChild(jj)->getID() << endl;
         elems.erase(elems.begin()+nd->getChild(jj)->getID()+count);
         count--;
      }
      */
   }


  return;
}





void HBSplineBase::performAdaptiveRefinement(double eTol)
{
    computeElementErrors(2);
    
    totalError *= eTol;
    
    int  ii;
    node  *nd, *nd1, *nd2;
    
    elemsToRefine.clear();
    
    if(DIM == 1)
    {    
    for(ii=0;ii<elems.size();ii++)
    {
       nd = elems[ii];
       
       if( !(nd->isGhost()) &&  nd->isLeaf() ) //&& (elems[ii]->getLevel() == MAX_LEVEL) )
       {
          if( nd->getError() > totalError )
          {
            //cout << " nd->isProcessed() " << '\t' << nd->getID() << '\t' << nd->isProcessed() << endl;

            elemsToRefine.push_back( nd->getID() );

            // add all the immediate neighbours
            
            nd2 = nd->getNeighbour(LEFT);
            if( nd2 != NULL )
              elemsToRefine.push_back(nd2->getID());

            nd2 = nd->getNeighbour(RIGHT);
            if( nd2 != NULL )
              elemsToRefine.push_back(nd2->getID());
          }
       }
    }
    }
    else if(DIM == 2)
    {
    for(ii=0;ii<elems.size();ii++)
    {
       nd = elems[ii];
       
       if( !(nd->isGhost()) &&  nd->isLeaf() && (nd->getLevel() == CURRENT_LEVEL) )
       {
          if( nd->getError() > totalError )
          {
            //cout << " nd->isProcessed() " << '\t' << nd->getID() << '\t' << nd->isProcessed() << endl;

            elemsToRefine.push_back( nd->getID() );

            // add all the immediate neighbours
            
            nd1 = nd->getNeighbour(WEST);
            if(nd1 != NULL)
            {
               elemsToRefine.push_back(nd1->getID());

               nd2 = nd1->getNeighbour(NORTH);
               if(nd2 != NULL)
                 elemsToRefine.push_back(nd2->getID());

               nd2 = nd1->getNeighbour(SOUTH);
               if(nd2 != NULL)
                 elemsToRefine.push_back(nd2->getID());
            }

            nd1 = nd->getNeighbour(EAST);
            if(nd1 != NULL)
            {
               elemsToRefine.push_back(nd1->getID());

               nd2 = nd1->getNeighbour(NORTH);
               if(nd2 != NULL)
                 elemsToRefine.push_back(nd2->getID());

               nd2 = nd1->getNeighbour(SOUTH);
               if(nd2 != NULL)
                 elemsToRefine.push_back(nd2->getID());
            }

            nd1 = nd->getNeighbour(NORTH);
            if(nd1 != NULL)
              elemsToRefine.push_back(nd1->getID());

            nd1 = nd->getNeighbour(SOUTH);
            if(nd1 != NULL)
              elemsToRefine.push_back(nd1->getID());
          }
       }
    }
    }

    sort(elemsToRefine.begin(), elemsToRefine.end());
   
    elemsToRefine.erase(unique(elemsToRefine.begin(), elemsToRefine.end()), elemsToRefine.end());

    //cout << " elemsToRefine  " << endl;
    //for(ii=0;ii<elemsToRefine.size();ii++)
      //cout << '\t' << elemsToRefine[ii] ;
    //cout << endl;   cout << endl;   cout << endl;

    /*
    // check for linear dependence when middle element is not refined
    
    nd = elems[0];
    
    while(nd1 != NULL)
    {
       if(nd->isLeaf())
         nd = nd->getNeighbour(RIGHT);
       else
         nd = nd->getChild(LEFT);

       if(nd->getLevel() == CURRENT_LEVEL)
         nd1 == NULL;
    }
    
    nd1 = nd->getNeighbour(RIGHT);

    cout << " AAAAAAAAAAAAAAAAA " << nd->getID() << endl;
    
    
    while( nd1 != NULL && !(nd1->isGhost()) )
    {
       cout << " CCCCCCCCCCCCCCCCC " << endl;
       nd2 = nd1->getNeighbour(RIGHT);
       
       cout << " AAAAAAAAAAAAAAAAA " << nd1->getID() << '\t' << nd2->getID() << endl;
       
       if( nd2 != NULL && !(nd2->isGhost()) )
       {
         if( binary_search(elemsToRefine.begin(), elemsToRefine.end(), nd->getID()) && binary_search(elemsToRefine.begin(), elemsToRefine.end(), nd2->getID()) )
          elemsToRefine.push_back(nd1->getID());
       }
       cout << " AAAAAAAAAAAAAAAAA " << endl;
       nd  = nd->getNeighbour(RIGHT);
       nd1 = nd1->getNeighbour(RIGHT);
       cout << " BBBBBBBBBBBBBBBBB " << endl;
    }
    */

    
    /*
    elemsToRefine.clear();
    
    nd2=elems[boundaryNodes[1]];
    elemsToRefine.push_back(nd2->getID());
    
    for(ii=0;ii<degree[2];ii++)
    {
       nd2 = nd2->getNeighbour(LEFT);
       elemsToRefine.push_back(nd2->getID());
    }
    */
    /*
    elemsToRefine.clear();
    
    int start, N=degree[2];
    
    if(CURRENT_LEVEL == 0)
    {
      start = nelem[0]/2 - N;
      cout << N << '\t' << start << endl;
      for(ii=0;ii<2*N;ii++)
        elemsToRefine.push_back(start+ii);
    }
    else
    {
      start = nelem[0] + 2*degree[0] + 4*N*(CURRENT_LEVEL-1)+N;
      for(ii=0;ii<2*N;ii++)
        elemsToRefine.push_back(start+ii);
    }

    cout << " elemsToRefine  " << endl;
    for(ii=0;ii<elemsToRefine.size();ii++)
      cout << '\t' << elemsToRefine[ii] ;
    cout << endl;   cout << endl;   cout << endl;
    */
    
    MAX_LEVEL += 1;
    applyRefinementProcess();
    CURRENT_LEVEL += 1;

    prepareMatrixPattern();
    
    return;
}



