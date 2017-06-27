
#ifndef incl_AdaptiveBinarytree_h
#define incl_AdaptiveBinarytree_h



#include "headersBasic.h"
#include "util.h"
#include "GeomDataHBSplines.h"
#include "GaussQuadrature.h"

#include "AABB.h"
#include "AdaptiveOctree.h"


template <int DIM>
class AdaptiveBinarytree
{
    typedef AdaptiveBinarytree<DIM>*  AdaptiveBinarytree_PTR;

    //typedef Matrix<double, DIM, 1>  knotPoint;

    private:
    
        static int  nodecount, MAX_LEVEL;

        int  level, id, NUM_CHILDREN, NUM_NEIGHBOURS;
        NodeOrientation  orientation;

        double  knots[3][4], volume, JacMultElem, coord3, param3;

        myPoint  knotBegin, knotEnd, knotIncr;

        AdaptiveBinarytree_PTR  *child, parent, *neighbours;

        int splitDir, splitCount, sideTemp;

        AABB  bbox;

    public:

        vector<int>  domNums;

        GeomDataHBSplines  *GeomData;

        GaussQuadrature  Quadrature;

        AdaptiveBinarytree();
        
        AdaptiveBinarytree(int lev=0);
        
        ~AdaptiveBinarytree();
        
        int  getLevel()
        {  return level; }
           
        void setLevel(int lev)
        {  level = lev; }

        bool isLeaf()
        {  return (child == NULL); }

        int getID()
        {  return id;}
        
        int getDimension()
        {  return DIM;}

        void SetSideTemp(int dd)
        {  sideTemp = dd;  return ;}

        int GetSideTemp()
        {  return sideTemp;}

        void SetParam3(double dd)
        {  param3 = dd;  return ;}

        int GetParam3()
        {  return param3; }

        void SetCoord3(double dd)
        {  coord3 = dd;  return ;}

        int GetCoord3()
        {  return coord3; }

        double  getVolume()
        {  return volume; }

        void setKnots(double* knots_)
        {  knots = knots_; }
        
        double* getKnots(int ind)
        {  return knots[ind]; }
        
        void setParent(AdaptiveBinarytree_PTR  parent1)
        {  parent = parent1; }
        
        AdaptiveBinarytree_PTR  getParent()
        {  return parent; }

        NodeOrientation  GetOrientation()
        {  return  orientation; }

        int getNumberOfChildren()
        {  return NUM_CHILDREN; }

        int getNumberOfNeighbours()
        {  return NUM_NEIGHBOURS; }

        void setNeighbour(int ind, AdaptiveBinarytree_PTR node1)
        {
           assert(ind < pow(2,DIM));
           neighbours[ind] = node1;
        }

        AdaptiveBinarytree_PTR  getNeighbour(int ind) 
        {  return neighbours[ind]; }

        AdaptiveBinarytree_PTR  getChild(int ind)
        {  return child[ind];	}
        
        static int  getCount()
        { return nodecount; }

        void setKnots(int index, double val0, double val1)
        {
           knots[index][0] = val0;
           knots[index][1] = val1;
        }

        void setKnots(double u0, double u1, double v0, double v1)
        {
           knots[0][0] = u0;
           knots[0][1] = u1;
           knots[1][0] = v0;
           knots[1][1] = v1;
        }

        void setKnots(double u0, double u1, double v0, double v1, double w0, double w1)
        {
           knots[0][0] = u0;
           knots[0][1] = u1;
           knots[1][0] = v0;
           knots[1][1] = v1;
           knots[2][0] = w0;
           knots[2][1] = w1;
        }

        AABB& getAABB()
        {  return  bbox; }

        double  getKnotspan(int dir)
        {  return (knots[dir][1] - knots[dir][0]);  }

        double  getKnotAt(int dir, int loc)
        {  return  knots[dir][loc];  }

        //void SetDomainNumber(int  dd)
        //{ domainNum = dd;  } 

        //int getDomainNumber()
        //{  return domainNum; }

        int getDomainNumber()
        {  return ( (domNums.size() == 1) ? domNums[0] : -1); }

        //bool isCutElement()
        //{  return (domainNum == -1); }

        bool isCutElement()
        {  return (domNums.size() > 1); }

        void SetSplitDirection(int  dd)
        { splitDir = dd;  } 

        int GetSplitDirection()
        {  return  splitDir; }

        void SetSplitCount(int  dd)
        { splitCount = dd;  } 

        int GetSplitCount()
        {  return  splitCount; }

        bool  within(myPoint& pt);
        
        void  subDivide(int lev);

        void  plotSelf();
        
        void  printSelf();

        void  prepareData();

        void  computeGaussPoints(int, int, int, int, GaussQuadrature&  quadTemp);

        void  computeGaussPoints2Dfor3D(int, int, int, int, GaussQuadrature&  quadTemp);

        void  computeGaussPointsForMerging(int, int, int, int, GaussQuadrature&  quadTemp);
        
        void  computeGaussPointsForMerging2Dfor3D(int, int, int, int, GaussQuadrature&  quadTemp);

        void  mergeGaussPoints(int, int, int, int&, myPoint&, double&);

};



template <int DIM>
AdaptiveBinarytree<DIM>::AdaptiveBinarytree(): level(0), NUM_CHILDREN(0)
{
    child = NULL;
    neighbours = NULL;
    parent = NULL;
    sideTemp = -1;
    coord3 = 0.0;

    orientation = ENUM_PARENT;

    NUM_NEIGHBOURS = 2;

    neighbours = new AdaptiveBinarytree_PTR[NUM_NEIGHBOURS];
    for(int ii=0;ii<NUM_NEIGHBOURS;ii++)
      neighbours[ii] = NULL;

    id = nodecount++;
}




template <int DIM>
AdaptiveBinarytree<DIM>::AdaptiveBinarytree(int lev):level(lev), NUM_CHILDREN(0)
{
    child = NULL;
    neighbours = NULL;
    parent = NULL;
    sideTemp = -1;
    coord3 = 0.0;

    orientation = ENUM_PARENT;

    if(level > MAX_LEVEL)
      MAX_LEVEL = level;

    NUM_NEIGHBOURS = 2;

    neighbours = new AdaptiveBinarytree_PTR[NUM_NEIGHBOURS];
    for(int ii=0;ii<NUM_NEIGHBOURS;ii++)
      neighbours[ii] = NULL;

    id = nodecount++;
}



template <int DIM>
AdaptiveBinarytree<DIM>::~AdaptiveBinarytree()
{
//cout << " AdaptiveBinarytree<DIM>::~AdaptiveBinarytree() " << '\t' << nodecount << endl;
  --nodecount;
//
  if(child != NULL)
  {
    for(int ii=0;ii<NUM_CHILDREN;ii++)
      delete child[ii];

    delete [] child;
    child = NULL;
  }

  //if(neighbours != NULL)
  //{
    //for(int ii=0;ii<NUM_NEIGHBOURS;ii++)
      //delete neighbours[ii];

    //delete [] neighbours;
    neighbours = NULL;
  //}

  //if(parent != NULL)
  //{
    //delete parent;
    parent = NULL;
  //}
  //Quadrature.reset();
//
}



template<int DIM>
void AdaptiveBinarytree<DIM>::prepareData()
{
    //cout << knots[0][0] << '\t' << knots[0][1] << '\t' << knots[0][2] << '\t' << knots[0][3] << endl;
    //cout << knots[1][0] << '\t' << knots[1][1] << '\t' << knots[1][2] << '\t' << knots[1][3] << endl;

    JacMultElem = GeomData->GetJacobianFull();

    for(int ii=0;ii<DIM;ii++)
    {
      //cout << " ii = " << ii << endl;
      knots[ii][2] = knots[ii][1] - knots[ii][0];
      knots[ii][3] = knots[ii][1] + knots[ii][0];

      knotBegin[ii] = knots[ii][0];
      knotEnd[ii]   = knots[ii][1];
      knotIncr[ii]  = knots[ii][2];
      //bbox.printSelf();
      //cout << " ii = " << ii << endl;

      bbox.minBB[ii] = GeomData->ComputeCoord(ii, knots[ii][0]);
      bbox.maxBB[ii] = GeomData->ComputeCoord(ii, knots[ii][1]);
      //cout << bbox.minBB[ii] << '\t' << bbox.maxBB[ii] << endl;

      JacMultElem *= (0.5*knots[ii][2]);
    }


    if(sideTemp != -1)
    {
      switch(sideTemp)
      {
         case 0:
         case 1:

                bbox.minBB[0] = coord3 ;
                bbox.minBB[1] = GeomData->ComputeCoord(1, knots[0][0]) ;
                bbox.minBB[2] = GeomData->ComputeCoord(2, knots[1][0]) ;

                bbox.maxBB[0] = coord3 ;
                bbox.maxBB[1] = GeomData->ComputeCoord(1, knots[0][1]) ;
                bbox.maxBB[2] = GeomData->ComputeCoord(2, knots[1][1]) ;
                
                JacMultElem /= GeomData->GetJacobian(0);

        break;

        case 2:
        case 3:

                bbox.minBB[0] = GeomData->ComputeCoord(0, knots[0][0]) ;
                bbox.minBB[1] = coord3 ;
                bbox.minBB[2] = GeomData->ComputeCoord(2, knots[1][0]) ;

                bbox.maxBB[0] = GeomData->ComputeCoord(0, knots[0][1]) ;
                bbox.maxBB[1] = coord3 ;
                bbox.maxBB[2] = GeomData->ComputeCoord(2, knots[1][1]) ;
                
                JacMultElem /= GeomData->GetJacobian(1);

        break;

        case 4:
        case 5:

                bbox.minBB[0] = GeomData->ComputeCoord(0, knots[0][0]) ;
                bbox.minBB[1] = GeomData->ComputeCoord(1, knots[1][0]) ;
                bbox.minBB[2] = coord3 ;

                bbox.maxBB[0] = GeomData->ComputeCoord(0, knots[0][1]) ;
                bbox.maxBB[1] = GeomData->ComputeCoord(1, knots[1][1]) ;
                bbox.maxBB[2] = coord3 ;
                
                JacMultElem /= GeomData->GetJacobian(2);

        break;

        default :

            cout << " Invalid 'sideTemp' value in AdaptiveBinarytree<DIM>::prepareData() " << endl;
        break;
      } //switch(sideTemp)
    }

    return;
}


template<int DIM>
bool AdaptiveBinarytree<DIM>::within(myPoint& pt)
{
  return  bbox.within(pt);
}


#ifdef _DOMAIN1D
typedef AdaptiveBinarytree<1>  adapBinarytree;
#endif


#ifdef _DOMAIN2D
typedef AdaptiveBinarytree<2>  adapBinarytree;
#endif


#ifdef _DOMAIN3D
typedef AdaptiveBinarytree<3>  adapBinarytree;
#endif


#endif

