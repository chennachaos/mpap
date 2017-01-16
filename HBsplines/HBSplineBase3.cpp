
#include "HBSplineBase.h"



void  HBSplineBase::subDivide(int nodenum)
{
    int ii, jj, ll;

    node  *nd, *nd2;

    nd = elems[nodenum];
    if( nd->IsLeaf() )
    {
      nd->subDivide();

      ll = nd->GetLevel() + 1;

      for(jj=0;jj<nd->GetNumberOfChildren();jj++)
      {
         nd2 = nd->GetChild(jj);

         nd2->SolnData = &(SolnData);
         nd2->GeomData = &(GeomData);

         elems.push_back(nd2);
         NodeNumsAtLevel[ll].push_back(nd2->GetID());
      }
    }

    return;
}



void  HBSplineBase::ApplyRefinementProcess()
{
    //printf(" Current level = %2d \n", CURRENT_LEVEL);
    //printf(" Maximum level = %2d \n\n", MAX_LEVEL);

    Algorithm1(CURRENT_LEVEL);
    Algorithm2(CURRENT_LEVEL);
    Algorithm3(CURRENT_LEVEL+1);

    return;
}



void  HBSplineBase::Algorithm1(int lev)
{
    // nucleus operation for each element: create new leaves (children) of level k+1
    //
    //////////////////////////////////////////////////////////////////////////////////////
   
    //printf(" Algorithm1 for level = %2d STARTED \n ",lev);

    if(DIM == 1)
      Algorithm1_1D(lev);
    else if(DIM == 2)
      Algorithm1_2D(lev);
    else
      Algorithm1_3D(lev);


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

    //printf(" Algorithm1 for level = %2d FINISHED \n ",lev);

    return;
}



void  HBSplineBase::Algorithm2(int lev)
{
    // remove linear dependence between basis functions of level k and level k+1.
    
    //printf("\n\n\n  Algorithm2 for level = %2d STARTED \n ",lev);

    /*
    printf("\n\n\n");
    printf(" Nodes at level = %5d \n", lev);
    for(int ii=0;ii<NodeNumsAtLevel[lev].size();ii++)
       cout << '\t' << NodeNumsAtLevel[lev][ii] ;
    printf("\n\n\n");
    */

    if(DIM == 1)
      Algorithm2_1D(lev);
    else if(DIM == 2)
      Algorithm2_2D(lev);
    else
      Algorithm2_3D(lev);

    findUnique(VacantBFs);

    /*
    printf("\n\n\n");
    printf(" Vacant Basis Functions \n");
    for(int ii=0;ii<VacantBFs.size();ii++)
       cout << '\t' << VacantBFs[ii] ;
    printf("\n\n\n");
    */

    //printf(" Algorithm2 for level = %2d FINISHED \n", lev);

    return;
}



void  HBSplineBase::Algorithm3(int lev)
{
    // assign numbers to the new basis functions of level k+1
    
    //printf("\n\n\n Algorithm3 for level = %2d STARTED \n ",lev);

    if(DIM == 1)
      Algorithm3_1D(lev);
    else if(DIM == 2)
      Algorithm3_2D(lev);
    else
      Algorithm3_3D(lev);

    //printf("\n \t   Total number of basis functions    =  %5d\n\n", gridBF1);

    //printf(" Algorithm3 for level = %2d FINISHED \n ",lev);

    return;
}




void  HBSplineBase::addGhostNodes(node* nd, int dir)
{
    nd = nd->GetNeighbour(dir);

    if( nd != NULL )
    {
      if(nd->IsGhost())
      {
          nodes2divide.push_back(nd->GetID());
          nd = nd->GetNeighbour(dir);
         
          while( nd != NULL)
          {
             nodes2divide.push_back(nd->GetID());
             nd = nd->GetNeighbour(dir);
          }
      }
    }

    return;
}


void  HBSplineBase::Algorithm1_1D(int lev)
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
      if(nd->IsLeftBoundary())
        addGhostNodes(nd, LEFT);

      //cout << " AAAAAAAAA " << endl;
      nd = elems[nodenum];
      if(nd->IsRightBoundary())
        addGhostNodes(nd, RIGHT);
      //cout << " ee " << ee << '\t' << elemsToRefine[ee] << endl;
    }

    return;
}


void  HBSplineBase::Algorithm1_2D(int lev)
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
       if(nd->IsLeftBoundary())
         addGhostNodes(nd, WEST);
       
       nd = elems[nodenum];
       if(nd->IsRightBoundary())
         addGhostNodes(nd, EAST);
       
       nd = elems[nodenum];
       if(nd->IsTopBoundary())
         addGhostNodes(nd, NORTH);

       nd = elems[nodenum];
       if(nd->IsBottomBoundary())
         addGhostNodes(nd, SOUTH);


       // bottom-left corner
       nd = elems[nodenum];
       if(nd->IsLeftBoundary() && nd->IsBottomBoundary())
       {
          while(nd != NULL)
          {
            addGhostNodes(nd, WEST);
            nd = nd->GetNeighbour(SOUTH);
          }
       }
       
       // bottom-right corner
       nd = elems[nodenum];
       if(nd->IsRightBoundary() && nd->IsBottomBoundary())
       {
          while(nd != NULL)
          {
            addGhostNodes(nd, EAST);
            nd = nd->GetNeighbour(SOUTH);
          }
       }

       // top-right corner
       nd = elems[nodenum];
       if(nd->IsRightBoundary() && nd->IsTopBoundary())
       {
          while(nd != NULL)
          {
            addGhostNodes(nd, EAST);
            nd = nd->GetNeighbour(NORTH);
          }
       }

       // top-left corner
       nd = elems[nodenum];
       if(nd->IsLeftBoundary() && nd->IsTopBoundary())
       {
          while(nd != NULL)
          {
            addGhostNodes(nd, WEST);
            nd = nd->GetNeighbour(NORTH);
          }
       }
    }

    return;
}


void  HBSplineBase::Algorithm1_3D(int lev)
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
      if(nd->IsLeftBoundary())
        addGhostNodes(nd, WEST);

      nd = elems[nodenum];
      if(nd->IsRightBoundary())
        addGhostNodes(nd, EAST);

      nd = elems[nodenum];
      if(nd->IsTopBoundary())
        addGhostNodes(nd, NORTH);

      nd = elems[nodenum];
      if(nd->IsBottomBoundary())
        addGhostNodes(nd, SOUTH);

      nd = elems[nodenum];
      if(nd->IsBackBoundary())
        addGhostNodes(nd, BACK);

      nd = elems[nodenum];
      if(nd->IsFrontBoundary())
        addGhostNodes(nd, FRONT);

      // front-left edge
      nd = elems[nodenum];
      if(nd->IsFrontBoundary() && nd->IsLeftBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, FRONT);
          nd = nd->GetNeighbour(WEST);
        }
      }
       
      // front-right edge
      nd = elems[nodenum];
      if(nd->IsFrontBoundary() && nd->IsRightBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, FRONT);
          nd = nd->GetNeighbour(EAST);
        }
      }
       
      // front-bottom edge
      nd = elems[nodenum];
      if(nd->IsFrontBoundary() && nd->IsBottomBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, FRONT);
          nd = nd->GetNeighbour(SOUTH);
        }
      }

      // front-top edge
      nd = elems[nodenum];
      if(nd->IsFrontBoundary() && nd->IsTopBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, FRONT);
          nd = nd->GetNeighbour(NORTH);
        }
      }

      // back-left edge
      nd = elems[nodenum];
      if(nd->IsBackBoundary() && nd->IsLeftBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, BACK);
          nd = nd->GetNeighbour(WEST);
        }
      }
       
      // back-right edge
      nd = elems[nodenum];
      if(nd->IsBackBoundary() && nd->IsRightBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, BACK);
          nd = nd->GetNeighbour(EAST);
        }
      }

      // back-bottom edge
      nd = elems[nodenum];
      if(nd->IsBackBoundary() && nd->IsBottomBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, BACK);
          nd = nd->GetNeighbour(SOUTH);
        }
      }

      // back-top edge
      nd = elems[nodenum];
      if(nd->IsBackBoundary() && nd->IsTopBoundary())
      {
        while(nd != NULL)
        {
          addGhostNodes(nd, BACK);
          nd = nd->GetNeighbour(NORTH);
        }
      }

      //cout << " AAAAAAAAA " << endl;
      // front-left-bottom corner
      nd = elems[nodenum];
      if(nd->IsFrontBoundary() && nd->IsLeftBoundary() && nd->IsBottomBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, SOUTH);
            nd1 = nd1->GetNeighbour(WEST);
          }
          nd = nd->GetNeighbour(FRONT);
        }
      }

      // front-left-top corner
      nd = elems[nodenum];
      if(nd->IsFrontBoundary() && nd->IsLeftBoundary() && nd->IsTopBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, NORTH);
            nd1 = nd1->GetNeighbour(WEST);
          }
          nd = nd->GetNeighbour(FRONT);
        }
      }

      // front-right-bottom corner
      nd = elems[nodenum];
      if(nd->IsFrontBoundary() && nd->IsRightBoundary() && nd->IsBottomBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, SOUTH);
            nd1 = nd1->GetNeighbour(EAST);
          }
          nd = nd->GetNeighbour(FRONT);
        }
      }

      // front-right-top corner
      nd = elems[nodenum];
      if(nd->IsFrontBoundary() && nd->IsRightBoundary() && nd->IsTopBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, NORTH);
            nd1 = nd1->GetNeighbour(EAST);
          }
          nd = nd->GetNeighbour(FRONT);
        }
      }

      // back-left-bottom corner
      nd = elems[nodenum];
      if(nd->IsBackBoundary() && nd->IsLeftBoundary() && nd->IsBottomBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, SOUTH);
            nd1 = nd1->GetNeighbour(WEST);
          }
          nd = nd->GetNeighbour(BACK);
        }
      }

      // back-left-top corner
      nd = elems[nodenum];
      if(nd->IsBackBoundary() && nd->IsLeftBoundary() && nd->IsTopBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, NORTH);
            nd1 = nd1->GetNeighbour(WEST);
          }
          nd = nd->GetNeighbour(BACK);
        }
      }

      // back-right-bottom corner
      nd = elems[nodenum];
      if(nd->IsBackBoundary() && nd->IsRightBoundary() && nd->IsBottomBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, SOUTH);
            nd1 = nd1->GetNeighbour(EAST);
          }
          nd = nd->GetNeighbour(BACK);
        }
      }

      // back-right-top corner
      nd = elems[nodenum];
      if(nd->IsBackBoundary() && nd->IsRightBoundary() && nd->IsTopBoundary() )
      {
        while(nd != NULL)
        {
          nd1 = nd;
          while( nd1 != NULL )
          {
            addGhostNodes(nd1, NORTH);
            nd1 = nd1->GetNeighbour(EAST);
          }
          nd = nd->GetNeighbour(BACK);
        }
      }
    }

    return;
}



void  HBSplineBase::Algorithm2_1D(int lev)
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
       
       if(nd->IsGhost())
         flag2 = false;
       else
         flag2 = true;
      
       nd1 = nd;
       for(ii=0;ii<=degree[0];ii++)
       {
          if( nd1 == NULL || nd1->IsLeaf() )
          {
             flag1 = false;
             break;
          }
          
          if(!nd1->IsGhost())
            flag2 = true;

          nd1 = nd1->GetNeighbour(RIGHT);
       }

       //cout << nd->GetID() << '\t' << flag1 << '\t' << flag2 << endl;

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
                nd1 = nd1->GetNeighbour(RIGHT);
             }
          }
       }
    }

    return;
}


void  HBSplineBase::Algorithm2_2D(int lev)
{
    int ii, jj, val, ee, temp1, temp2;
    bool flag1, flag2;

    node  *nd, *nd1, *nd2;

    //cout << " Node # " << nd1->GetID() << '\t' << nd1->IsGhost() << endl;

    for(ee=0;ee<NodeNumsAtLevel[lev].size();ee++)
    {
       nd = elems[NodeNumsAtLevel[lev][ee]];
      
       flag1 = true;
      
       // if the current node is a ghost node then check if all the neighbours are ghost nodes
       // if all the neighbours are ghost nodes then don't remove any basis function

       if(nd->IsGhost())
         flag2 = false;
       else
         flag2 = true;

       nd2 = nd;
       for(jj=0;jj<=degree[1];jj++)
       {
         nd1 = nd2;
         for(ii=0;ii<=degree[0];ii++)
         {
            if( nd1 == NULL || nd1->IsLeaf() )
            {
              flag1 = false;
              break;
            }

            if( !nd1->IsGhost() )
              flag2 = true;

            nd1 = nd1->GetNeighbour(EAST);
         }
         if(!flag1)
           break;
         
         nd2 = nd2->GetNeighbour(NORTH);
      }

      //cout << " nd->GetID() ... : "  << '\t'<< nd->GetID() << '\t' << nd->IsGhost() << '\t' << flag1 << '\t' << flag2 << endl;

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
                  nd1 = nd1->GetNeighbour(EAST);
               }
               nd2 = nd2->GetNeighbour(NORTH);
            }
         }
      }
   }

   return;
}


void  HBSplineBase::Algorithm2_3D(int lev)
{
    int ii, jj, kk, val, ee, temp1, temp2;
    bool flag1, flag2;

    node  *nd, *nd1, *nd2, *nd3;

    //cout << " Node # " << nd1->GetID() << '\t' << nd1->IsGhost() << endl;

    for(ee=0;ee<NodeNumsAtLevel[lev].size();ee++)
    {
       nd = elems[NodeNumsAtLevel[lev][ee]];
      
       flag1 = true;
      
       // if the current node is a ghost node then check if all the neighbours are ghost nodes
       // if all the neighbours are ghost nodes then don't remove any basis function

       if(nd->IsGhost())
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
             if( nd1 == NULL || nd1->IsLeaf() )
             {
               flag1 = false;
               break;
             }

             if( !nd1->IsGhost() )
               flag2 = true;

             nd1 = nd1->GetNeighbour(EAST);
           }
           if(!flag1)
             break;
         
           nd2 = nd2->GetNeighbour(NORTH);
        }
        if(!flag1)
          break;

        nd3 = nd3->GetNeighbour(FRONT);
      }

      //cout << " nd->GetID() ... : "  << '\t'<< nd->GetID() << '\t' << nd->IsGhost() << '\t' << flag1 << '\t' << flag2 << endl;

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
                  nd1 = nd1->GetNeighbour(EAST);
                }
                nd2 = nd2->GetNeighbour(NORTH);
              }
              nd3 = nd3->GetNeighbour(FRONT);
            }
         }
      }
   }

    return;
}


void  HBSplineBase::Algorithm3_1D(int lev)
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
       
       if(nd->IsGhost())
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

         if(!nd1->IsGhost())
          flag2 = true;

         nd1 = nd1->GetNeighbour(RIGHT);
       }

       //cout << nd->GetID() << '\t' << flag1 << '\t' << flag2 << endl;

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
           nd1 = nd1->GetNeighbour(RIGHT);
         }
       }
    }

   return;
}



void  HBSplineBase::Algorithm3_2D(int lev)
{
    int ii, jj, val, ee, temp1, temp2;
    bool flag1, flag2;

    node  *nd, *nd1, *nd2;

    for(ee=0;ee<NodeNumsAtLevel[lev].size();ee++)
    {
       nd = elems[NodeNumsAtLevel[lev][ee]];
      
       flag1 = true;

       if(nd->IsGhost())
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

             if( !nd1->IsGhost() )
               flag2 = true;

             nd1 = nd1->GetNeighbour(EAST);
          }
          if(!flag1)
            break;
         
          nd2 = nd2->GetNeighbour(NORTH);
       }
      
       //cout << " nd->GetID() ... : "  << '\t'<< nd->GetID() << '\t' << nd->IsGhost() << '\t' << flag1 << '\t' << flag2 << endl;

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
                nd1 = nd1->GetNeighbour(EAST);
             }
             nd2 = nd2->GetNeighbour(NORTH);
          }
       }
    }
   
    return;
}




void  HBSplineBase::Algorithm3_3D(int lev)
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

       if(nd->IsGhost())
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
             if( !nd1->IsGhost() )
               flag2 = true;

             nd1 = nd1->GetNeighbour(EAST);
           }
           if(!flag1)
             break;
         
           nd2 = nd2->GetNeighbour(NORTH);
         }
         if(!flag1)
           break;

         nd3 = nd3->GetNeighbour(FRONT);
       }
      
       //cout << " nd->GetID() ... : "  << '\t'<< nd->GetID() << '\t' << nd->IsGhost() << '\t' << flag1 << '\t' << flag2 << endl;

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
                nd1 = nd1->GetNeighbour(EAST);
              }
              nd2 = nd2->GetNeighbour(NORTH);
            }
            nd3 = nd3->GetNeighbour(FRONT);
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

   nd = nd->GetNeighbour(LEFT);

      if(nd->IsGhost())
      {
         // take all ghost nodes to the left of this node

         nd = nd->GetNeighbour(LEFT);
         while( nd != NULL)
         {
            nodes2unrefine.push_back(nd->GetID());
            nd = nd->GetNeighbour(LEFT);
         }
      }

   nd = elems[nodenum];

   nd = nd->GetNeighbour(RIGHT);

      if(nd->IsGhost())
      {
         // take all ghost nodes to the left of this node

         nd = nd->GetNeighbour(RIGHT);
         while( nd != NULL)
         {
            nodes2unrefine.push_back(nd->GetID());
            nd = nd->GetNeighbour(RIGHT);
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
      nd->Activate();
      /*
      for(jj=0;jj<2;jj++)
      {
         cout << "  nd->GetChild(jj)->GetID()  " << nd->GetChild(jj)->GetID() << endl;
         elems.erase(elems.begin()+nd->GetChild(jj)->GetID()+count);
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
       
       if( !(nd->IsGhost()) &&  nd->IsLeaf() ) //&& (elems[ii]->GetLevel() == MAX_LEVEL) )
       {
          if( nd->GetError() > totalError )
          {
            //cout << " nd->IsProcessed() " << '\t' << nd->GetID() << '\t' << nd->IsProcessed() << endl;

            elemsToRefine.push_back( nd->GetID() );

            // add all the immediate neighbours
            
            nd2 = nd->GetNeighbour(LEFT);
            if( nd2 != NULL )
              elemsToRefine.push_back(nd2->GetID());

            nd2 = nd->GetNeighbour(RIGHT);
            if( nd2 != NULL )
              elemsToRefine.push_back(nd2->GetID());
          }
       }
    }
    }
    else if(DIM == 2)
    {
    for(ii=0;ii<elems.size();ii++)
    {
       nd = elems[ii];
       
       if( !(nd->IsGhost()) &&  nd->IsLeaf() && (nd->GetLevel() == CURRENT_LEVEL) )
       {
          if( nd->GetError() > totalError )
          {
            //cout << " nd->IsProcessed() " << '\t' << nd->GetID() << '\t' << nd->IsProcessed() << endl;

            elemsToRefine.push_back( nd->GetID() );

            // add all the immediate neighbours
            
            nd1 = nd->GetNeighbour(WEST);
            if(nd1 != NULL)
            {
               elemsToRefine.push_back(nd1->GetID());

               nd2 = nd1->GetNeighbour(NORTH);
               if(nd2 != NULL)
                 elemsToRefine.push_back(nd2->GetID());

               nd2 = nd1->GetNeighbour(SOUTH);
               if(nd2 != NULL)
                 elemsToRefine.push_back(nd2->GetID());
            }

            nd1 = nd->GetNeighbour(EAST);
            if(nd1 != NULL)
            {
               elemsToRefine.push_back(nd1->GetID());

               nd2 = nd1->GetNeighbour(NORTH);
               if(nd2 != NULL)
                 elemsToRefine.push_back(nd2->GetID());

               nd2 = nd1->GetNeighbour(SOUTH);
               if(nd2 != NULL)
                 elemsToRefine.push_back(nd2->GetID());
            }

            nd1 = nd->GetNeighbour(NORTH);
            if(nd1 != NULL)
              elemsToRefine.push_back(nd1->GetID());

            nd1 = nd->GetNeighbour(SOUTH);
            if(nd1 != NULL)
              elemsToRefine.push_back(nd1->GetID());
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
       if(nd->IsLeaf())
         nd = nd->GetNeighbour(RIGHT);
       else
         nd = nd->GetChild(LEFT);

       if(nd->GetLevel() == CURRENT_LEVEL)
         nd1 == NULL;
    }
    
    nd1 = nd->GetNeighbour(RIGHT);

    cout << " AAAAAAAAAAAAAAAAA " << nd->GetID() << endl;
    
    
    while( nd1 != NULL && !(nd1->IsGhost()) )
    {
       cout << " CCCCCCCCCCCCCCCCC " << endl;
       nd2 = nd1->GetNeighbour(RIGHT);
       
       cout << " AAAAAAAAAAAAAAAAA " << nd1->GetID() << '\t' << nd2->GetID() << endl;
       
       if( nd2 != NULL && !(nd2->IsGhost()) )
       {
         if( binary_search(elemsToRefine.begin(), elemsToRefine.end(), nd->GetID()) && binary_search(elemsToRefine.begin(), elemsToRefine.end(), nd2->GetID()) )
          elemsToRefine.push_back(nd1->GetID());
       }
       cout << " AAAAAAAAAAAAAAAAA " << endl;
       nd  = nd->GetNeighbour(RIGHT);
       nd1 = nd1->GetNeighbour(RIGHT);
       cout << " BBBBBBBBBBBBBBBBB " << endl;
    }
    */

    
    /*
    elemsToRefine.clear();
    
    nd2=elems[boundaryNodes[1]];
    elemsToRefine.push_back(nd2->GetID());
    
    for(ii=0;ii<degree[2];ii++)
    {
       nd2 = nd2->GetNeighbour(LEFT);
       elemsToRefine.push_back(nd2->GetID());
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
    ApplyRefinementProcess();
    CURRENT_LEVEL += 1;

    prepareMatrixPattern();
    
    return;
}



