
#include "HBSplineBase.h"

//using namespace std;



void  HBSplineBase::BuildBase()
{
  if(DIM == 1)
    BuildBase1D();
  else if(DIM == 2)
    BuildBase2D();
  else
    BuildBase3D();
  
  //cout << " HBSplineBase::BuildBase() " << endl;

  return;
}


void  HBSplineBase::BuildBase1D()
{
    int  kk, nn, ii, jj;

    double  du = (knotsRight[0] - knotsLeft[0])/nelem[0];
    double  val1(0.0), val2(0.0);
    
    node  *nodetmp;
    
    nn = nelem[0] + 2*degree[0];
    
    // create the elements/nodes
    /////////////////////////////////////

    val1 = -du*degree[0];
    for(ii=0;ii<nn;ii++)
    {
        nodetmp = new node(0);
        
        //nodetmp->SetDegree(degree);
        
        elems.push_back(nodetmp);
        
        //elems.push_back(std::tr1::unqiue_ptr<node> node(0));
        
        elems[ii]->SetDegree(degree);
        
        NodeNumsAtLevel[0].push_back(nodetmp->GetID());

        val2 = val1 + du;

        elems[ii]->SetKnots(Dir1, val1, val2);
        
        val1 = val2;
    }

    // set the neighbours
    /////////////////////////////////////
    
    nn -= 1;

    elems[0]->SetNeighbour(RIGHT, elems[1]);
    elems[0]->SolnData = &(SolnData);
    elems[0]->GeomData = &(GeomData);

    elems[nn]->SetNeighbour(LEFT, elems[nn-1]);
    elems[nn]->SolnData = &(SolnData);
    elems[nn]->GeomData = &(GeomData);

    for(ii=1;ii<nn;ii++)
    {
       elems[ii]->SetNeighbour(LEFT,  elems[ii-1]);
       elems[ii]->SetNeighbour(RIGHT, elems[ii+1]);
        
       elems[ii]->SolnData = &(SolnData);
       elems[ii]->GeomData = &(GeomData);
    }

    // turn ON the ghost elements
    /////////////////////////////////////

    for(ii=0;ii<degree[0];ii++)
    {
       elems[ii]->SetGhostOn();
       elems[nn-ii]->SetGhostOn();
    }

    // identify the boundary nodes
    /////////////////////////////////////

    boundaryNodes[0].push_back(degree[0]);
    boundaryNodes[1].push_back(nelem[0]+degree[0]-1);

    return;
}



void  HBSplineBase::BuildBase2D()
{
    int  ii, jj, kk, ll, nn0, nn1, count=0, nor, sou, ind1, ind2;

    double  du = (knotsRight[0] - knotsLeft[0])/nelem[0];
    double  dv = (knotsRight[1] - knotsLeft[1])/nelem[1];
    
    double  u1, u2, v1, v2;

    node  *nodetmp;

    // create the elements/nodes
    /////////////////////////////////////
    
    nn0 = nelem[0] + 2*degree[0];
    nn1 = nelem[1] + 2*degree[1];

    count=0;
    v1 = -dv*degree[1];
    
    for(jj=0;jj<nn1;jj++)
    {
        v2 = v1 + dv ;
        u1 = -du*degree[0];
        
        for(ii=0;ii<nn0;ii++)
        {
            nodetmp = new node(0);

            //nodetmp = std::make_unique<node>(0); // won't leak, frees itself

            nodetmp->SetDegree(degree);
        
            u2 = u1 + du;

            nodetmp->SetKnots(Dir1, u1, u2);
            nodetmp->SetKnots(Dir2, v1, v2);

            nodetmp->SolnData = &(SolnData);
            nodetmp->GeomData = &(GeomData);

            elems.push_back(nodetmp);
            
            NodeNumsAtLevel[0].push_back(nodetmp->GetID());
        
            u1 = u2;
            count++;
        }
        v1 = v2;
    }

    //cout <<  " total number of elements = " << elems.size() << endl;

    root = elems[0];

    // set the neighbours
    /////////////////////////////////////

    for(jj=0;jj<nn1;jj++)
    {
        ind1 = nn0 * jj;
        nor  = nn0 * (jj+1);
        sou  = nn0 * (jj-1);

        if(jj == 0)
        {
           for(ii=0;ii<nn0;ii++)
           {
              ind2 = ind1 + ii;
              
              elems[ind2]->SetNeighbour(NORTH, elems[nor+ii]);
              if(ii == 0)
                elems[ind2]->SetNeighbour(EAST,  elems[ind2+1]);
              else if(ii == (nn0-1))
                elems[ind2]->SetNeighbour(WEST,  elems[ind2-1]);
              else
              {
                elems[ind2]->SetNeighbour(EAST,  elems[ind2+1]);
                elems[ind2]->SetNeighbour(WEST,  elems[ind2-1]);
              }
           }
        }
        else if(jj == (nn1-1))
        {
           for(ii=0;ii<nn0;ii++)
           {
              ind2 = ind1 + ii;

              elems[ind2]->SetNeighbour(SOUTH, elems[sou+ii]);
              if(ii == 0)
                elems[ind2]->SetNeighbour(EAST,  elems[ind2+1]);
              else if(ii == (nn0-1))
                elems[ind2]->SetNeighbour(WEST,  elems[ind2-1]);
              else
              {
                elems[ind2]->SetNeighbour(EAST,  elems[ind2+1]);
                elems[ind2]->SetNeighbour(WEST,  elems[ind2-1]);
              }
           }
        }
        else
        {
           for(ii=0;ii<nn0;ii++)
           {
              ind2 = ind1 + ii;
              elems[ind2]->SetNeighbour(NORTH, elems[nor+ii]);
              elems[ind2]->SetNeighbour(SOUTH, elems[sou+ii]);
              if(ii == 0)
                elems[ind2]->SetNeighbour(EAST,  elems[ind2+1]);
              else if(ii == (nn0-1))
                elems[ind2]->SetNeighbour(WEST,  elems[ind2-1]);
              else
              {
                elems[ind2]->SetNeighbour(EAST,  elems[ind2+1]);
                elems[ind2]->SetNeighbour(WEST,  elems[ind2-1]);
              }
           }
        }
    }

    // turn ON the ghost elements
    /////////////////////////////////////
    
    // direction #1

    for(jj=0;jj<degree[1];jj++)
    {
        ind1 = nn0 * jj;
        ind2 = nn0*(nn1-1-jj);

        for(ii=0;ii<nn0;ii++)
        {
           elems[ind1+ii]->SetGhostOn();
           elems[ind2+ii]->SetGhostOn();
        }
    }
    
    // direction #2

    for(jj=0;jj<(nn1-2*degree[1]);jj++)
    {
        ind1 = nn0 * (degree[1]+jj);
        ind2 = ind1 + nn0 - 1;

        for(ii=0;ii<degree[0];ii++)
        {
          elems[ind1+ii]->SetGhostOn();
          elems[ind2-ii]->SetGhostOn();
        }
    }

    // identify the boundary nodes
    /////////////////////////////////////

    // left and right columns

    for(ii=0;ii<nelem[1];ii++)
    {
      ind1 = nn0*(degree[1]+ii) + degree[0];
      ind2 = ind1 + nelem[0] - 1;

      boundaryNodes[0].push_back(ind1);
      boundaryNodes[1].push_back(ind2);
    }
    
    // top and bottom rows

    ind1 = nn0*degree[1] + degree[0];
    ind2 = nn0*(nn1-degree[1]-1) + degree[0];

    for(ii=0;ii<nelem[0];ii++)
    {
      boundaryNodes[2].push_back(ind1+ii);
      boundaryNodes[3].push_back(ind2+ii);
    }
    
    //findUnique(boundaryNodes);

    // account for periodic boundary conditions if any
    /////////////////////////////////////////////////////

    return;
}




void  HBSplineBase::BuildBase3D()
{
    int  ii, jj, kk, ll, ee, nn0, nn1, nn2, count=0, eas, wes, nor, sou, fro, bac, ind1, ind2, nelm;

    double  du = (knotsRight[0] - knotsLeft[0])/nelem[0];
    double  dv = (knotsRight[1] - knotsLeft[1])/nelem[1];
    double  dw = (knotsRight[2] - knotsLeft[2])/nelem[2];

    double  u1, u2, v1, v2, w1, w2;

    node  *nodetmp;

    // create the elements/nodes
    /////////////////////////////////////
    
    nn0 = nelem[0] + 2*degree[0];
    nn1 = nelem[1] + 2*degree[1];
    nn2 = nelem[2] + 2*degree[2];
    
    //cout << nn0 << '\t' << nn1 << '\t' << nn2 << endl;
    //cout << degree[0] << '\t' << degree[1] << '\t' << degree[2] << endl;
    //cout << du << '\t' << dv << '\t' << dw << endl;

    w1 = -dw*degree[2];

    for(kk=0;kk<nn2;kk++)
    {
       w2 = w1 + dw;
       v1 = -dv*degree[1];

       for(jj=0;jj<nn1;jj++)
       {
          v2 = v1 + dv ;
          u1 = -du*degree[0];

          for(ii=0;ii<nn0;ii++)
          {
             nodetmp = new node(0);
        
             nodetmp->SetDegree(degree);
        
             u2 = u1 + du;

             nodetmp->SetKnots(Dir1, u1, u2);
             nodetmp->SetKnots(Dir2, v1, v2);
             nodetmp->SetKnots(Dir3, w1, w2);

             nodetmp->SolnData = &(SolnData);
             nodetmp->GeomData = &(GeomData);

             elems.push_back(nodetmp);
            
             NodeNumsAtLevel[0].push_back(nodetmp->GetID());
        
             u1 = u2;
          }
          v1 = v2;
       }
       w1 = w2;
    }
    
    //cout <<  " total number of elements = " << elems.size() << endl;
    
    // set the neighbours
    /////////////////////////////////////
    //
    // ee = nn0*nn1*kk + nn0*jj + ii

    nelm = elems.size();

    for(ee=0;ee<nelm;ee++)
    {
      if((ee ) % nn0 == 0)
        wes = - 1;
      else
      {
        wes = ee - 1;
        elems[ee]->SetNeighbour(WEST,   elems[wes]);
      }

      if( (ee+1) % nn0 == 0)
        eas = -1;
      else
      {
        eas = ee + 1;
        elems[ee]->SetNeighbour(EAST,   elems[eas]);
      }

      if( (ee%(nn0*nn1)-nn0) < 0 )
        sou = -1;
      else
      {
        sou = ee - nn0;
        elems[ee]->SetNeighbour(SOUTH,  elems[sou]);
      }

      if( (ee%(nn0*nn1)+nn0) >= nn0*nn1 )
        nor = -1;
      else
      {
        nor = ee + nn0;
        elems[ee]->SetNeighbour(NORTH,  elems[nor]);
      }

      if( (ee - nn0*nn1) < 0 )
        bac = -1;
      else
      {
        bac = ee - nn0*nn1;
        elems[ee]->SetNeighbour(BACK,   elems[bac]);
      }

      if( (ee + nn0*nn1) >= nelm )
        fro = -1;
      else
      {
        fro = ee + nn0*nn1;
        elems[ee]->SetNeighbour(FRONT,   elems[fro]);
      }

      //printf("%5d \t %5d \t %5d \t %5d \t %5d \t %5d \t %5d \n",ee,wes,eas,sou,nor,bac,fro);
    }

    // turn ON the ghost elements
    /////////////////////////////////////

    // left degree[0] blocks

    for(kk=0;kk<nn2;kk++)
    {
      ind2 = nn0*nn1*kk;
      for(jj=0;jj<nn1;jj++)
      {
        ind1 = ind2 + nn0 * jj;
        for(ii=0;ii<degree[0];ii++)
          elems[ind1+ii]->SetGhostOn();
      }
    }

    // right degree[0] blocks

    for(kk=0;kk<nn2;kk++)
    {
      ind2 = nn0*nn1*kk;
      for(jj=0;jj<nn1;jj++)
      {
        ind1 = ind2 + nn0 * jj;
        for(ii=(nn0-degree[0]);ii<nn0;ii++)
          elems[ind1+ii]->SetGhostOn();
      }
    }

    // bottom degree[1] blocks

    for(kk=0;kk<nn2;kk++)
    {
      ind2 = nn0*nn1*kk;
      for(jj=0;jj<degree[1];jj++)
      {
        ind1 = ind2 + nn0 * jj;
        for(ii=0;ii<nn0;ii++)
          elems[ind1+ii]->SetGhostOn();
      }
    }

    // top degree[1] blocks

    for(kk=0;kk<nn2;kk++)
    {
      ind2 = nn0*nn1*kk;
      for(jj=(nn1-degree[1]);jj<nn1;jj++)
      {
        ind1 = ind2 + nn0 * jj;
        for(ii=0;ii<nn0;ii++)
          elems[ind1+ii]->SetGhostOn();
      }
    }

    // back degree[2] blocks

    for(kk=0;kk<degree[2];kk++)
    {
      ind2 = nn0*nn1*kk;
      for(jj=0;jj<nn1;jj++)
      {
        ind1 = ind2 + nn0 * jj;
        for(ii=0;ii<nn0;ii++)
          elems[ind1+ii]->SetGhostOn();
      }
    }

    // front degree[2] blocks

    for(kk=(nn2-degree[2]);kk<nn2;kk++)
    {
      ind2 = nn0*nn1*kk;
      for(jj=0;jj<nn1;jj++)
      {
        ind1 = ind2 + nn0 * jj;
        for(ii=0;ii<nn0;ii++)
          elems[ind1+ii]->SetGhostOn();
      }
    }

    // identify the boundary nodes
    /////////////////////////////////////

    // left block

    ii = degree[0];

    for(kk=degree[2];kk<(nn2-degree[2]);kk++)
    {
      ind2 = nn0*nn1*kk;
      for(jj=degree[1];jj<(nn1-degree[1]);jj++)
      {
        ind1 = ind2 + nn0 * jj;
        boundaryNodes[0].push_back(ind1+ii);
      }
    }

    // right block

    ii = nn0 - degree[0] - 1;

    for(kk=degree[2];kk<(nn2-degree[2]);kk++)
    {
      ind2 = nn0*nn1*kk;
      for(jj=degree[1];jj<(nn1-degree[1]);jj++)
      {
        ind1 = ind2 + nn0 * jj;
        boundaryNodes[1].push_back(ind1+ii);
      }
    }

    // bottom block

    jj = degree[1];

    for(kk=degree[2];kk<(nn2-degree[2]);kk++)
    {
      ind2 = nn0*nn1*kk;
      ind1 = ind2 + nn0 * jj;
      for(ii=degree[0];ii<(nn0-degree[0]);ii++)
        boundaryNodes[2].push_back(ind1+ii);
    }

    // top block

    jj = nn1 - degree[1] - 1;

    for(kk=degree[2];kk<(nn2-degree[2]);kk++)
    {
      ind2 = nn0*nn1*kk;
      ind1 = ind2 + nn0 * jj;
      for(ii=degree[0];ii<(nn0-degree[0]);ii++)
        boundaryNodes[3].push_back(ind1+ii);
    }

    // back block

    kk = degree[2];

    ind2 = nn0*nn1*kk;
    for(jj=degree[1];jj<(nn1-degree[1]);jj++)
    {
      ind1 = ind2 + nn0 * jj;
      for(ii=degree[0];ii<(nn0-degree[0]);ii++)
        boundaryNodes[4].push_back(ind1+ii);
    }

    // front block

    kk = nn2 - degree[2] - 1;

    ind2 = nn0*nn1*kk;
    for(jj=degree[1];jj<(nn1-degree[1]);jj++)
    {
      ind1 = ind2 + nn0 * jj;
      for(ii=degree[0];ii<(nn0-degree[0]);ii++)
        boundaryNodes[5].push_back(ind1+ii);
    }

    //findUnique(boundaryNodes);

    // account for periodic boundary conditions if any
    /////////////////////////////////////////////////////

    return;
}









