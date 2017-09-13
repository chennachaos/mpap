
#ifndef incl_AdaptiveOctree_h
#define incl_AdaptiveOctree_h


#include "headersBasic.h"
#include "util.h"
#include "GeomDataHBSplines.h"
#include "GaussQuadrature.h"
#include "AABB.h"


enum NodeOrientation
{
  ENUM_PARENT = -1,
  FIRST=0, SECOND=1,
  LEFT=0, RIGHT=1,
  SW=0, SE=1, NW=2, NE=3,
  WEST=0, EAST=1, NORTH=2, SOUTH=3, BACK=4, FRONT=5,
  SW_BACK =0, SE_BACK =1, NW_BACK =2, NE_BACK=3,
  SW_FRONT=4, SE_FRONT=5, NW_FRONT=6, NE_FRONT=7
};

enum {Dir1=0, Dir2, Dir3};

enum elemTypeCutFEM { ELEM_INSIDE=0, ELEM_OUTSIDE=1, ELEM_CUT=2 };


template <int DIM>
class AdaptiveOctree
{
    typedef AdaptiveOctree<DIM>*  AdaptiveOctree_PTR;

    //typedef Matrix<double, DIM, 1>  knotPoint;

    private:
    
        static int  nodecount, MAX_LEVEL;

        int  level, id, NUM_CHILDREN, NUM_NEIGHBOURS, sideTemp;
        NodeOrientation  orientation;

        double  knots[3][4], volume, JacMultElem, coord3, param3;

        myPoint  knotBegin, knotEnd, knotIncr;

        AdaptiveOctree_PTR  *child, parent;//, *neighbours;

        AdaptiveOctree_PTR  neighbours[2*DIM];

        AABB  bbox;

    public:

        vector<int>  domNums;

        GeomDataHBSplines  *GeomData;

        GaussQuadrature  Quadrature;

        AdaptiveOctree();
        
        AdaptiveOctree(int lev=0);
        
        ~AdaptiveOctree();
        
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

        void setSideTemp(int dd)
        {  sideTemp = dd;  return ;}

        int getSideTemp()
        {  return sideTemp;}

        void setParam3(double dd)
        {  param3 = dd;  return ;}

        int getParam3()
        {  return param3; }

        void setCoord3(double dd)
        {  coord3 = dd;  return ;}

        int getCoord3()
        {  return coord3; }

        double  getVolume()
        {  return volume; }

        void setParent(AdaptiveOctree_PTR  parent1)
        {  parent = parent1; }
        
        AdaptiveOctree_PTR  getParent()
        {  return parent; }
        
        NodeOrientation  getOrientation()
        {  return  orientation; }
        
        int getNumberOfChildren()
        {  return NUM_CHILDREN; }

        int getNumberOfNeighbours()
        {  return NUM_NEIGHBOURS; }
        
        void setNeighbour(int ind, AdaptiveOctree_PTR node1)
        {
           assert(ind < pow(2,DIM));
           neighbours[ind] = node1;
        }
        
        AdaptiveOctree_PTR  getNeighbour(int ind) 
        {  return neighbours[ind]; }

        AdaptiveOctree_PTR  getChild(int ind)
        {  return child[ind];	}
        
        static int  getCount()
        { return nodecount; }

        void setKnots(double* knots_)
        {  knots = knots_; }
        
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

        double* getKnots(int ind)
        {  return knots[ind]; }

        double  getKnotspan(int dir)
        {  return  knots[dir][2];  }

        double  getKnotAt(int dir, int loc)
        {  return  knots[dir][loc];  }

        myPoint& getKnotBegin()
        {  return knotBegin; }

        myPoint& getKnotEnd()
        {  return knotEnd; }

        myPoint& getKnotIncrement()
        {  return knotIncr; }

        AABB& getAABB()
        {  return  bbox; }

        int getDomainNumber()
        {  return ( (domNums.size() == 1) ? domNums[0] : -1); }

        bool isCutElement()
        {  return (domNums.size() > 1); }

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
AdaptiveOctree<DIM>::AdaptiveOctree(): level(0), NUM_CHILDREN(0)
{
    child = NULL;
    //neighbours = NULL;
    parent = NULL;
    sideTemp = -1;
    coord3 = 0.0;

    orientation = ENUM_PARENT;

    NUM_NEIGHBOURS = 2*DIM;

    //neighbours = new AdaptiveOctree_PTR[NUM_NEIGHBOURS];
    for(int ii=0;ii<NUM_NEIGHBOURS;ii++)
      neighbours[ii] = NULL;

    id = nodecount++;
}




template <int DIM>
AdaptiveOctree<DIM>::AdaptiveOctree(int lev):level(lev), NUM_CHILDREN(0)
{
    child = NULL;
    //neighbours = NULL;
    parent = NULL;
    sideTemp = -1;
    coord3 = 0.0;

    orientation = ENUM_PARENT;

    if(level > MAX_LEVEL)
      MAX_LEVEL = level;

    NUM_NEIGHBOURS = 2*DIM;

    //neighbours = new AdaptiveOctree_PTR[NUM_NEIGHBOURS];
    for(int ii=0;ii<NUM_NEIGHBOURS;ii++)
      neighbours[ii] = NULL;

    id = nodecount++;
}



template <int DIM>
AdaptiveOctree<DIM>::~AdaptiveOctree()
{
  //cout << " AdaptiveOctree<DIM>::~AdaptiveOctree() " << '\t' << nodecount << endl;
  --nodecount;

  if(child != NULL)
  {
    for(int ii=0;ii<NUM_CHILDREN;ii++)
      delete child[ii];

    delete [] child;
    child = NULL;
  }

  //if(neighbours != NULL)
  //{
    for(int ii=0;ii<NUM_NEIGHBOURS;ii++)
      neighbours[ii] = NULL;
      //delete neighbours[ii];

    //delete [] neighbours;
    //neighbours = NULL;
  //}

  //if(parent != NULL)
  //{
    //delete parent;
    parent = NULL;
  //}

  Quadrature.reset();
}



template<int DIM>
void AdaptiveOctree<DIM>::prepareData()
{
    JacMultElem = GeomData->getJacobianFull();

    for(int ii=0;ii<DIM;ii++)
    {
      knots[ii][2] = knots[ii][1] - knots[ii][0];
      knots[ii][3] = knots[ii][1] + knots[ii][0];

      knotBegin[ii] = knots[ii][0];
      knotEnd[ii]   = knots[ii][1];
      knotIncr[ii]  = knots[ii][2];

      bbox.minBB[ii] = GeomData->computeCoord(ii, knots[ii][0]);
      bbox.maxBB[ii] = GeomData->computeCoord(ii, knots[ii][1]);

      JacMultElem *= (0.5*knots[ii][2]);
    }


    if(sideTemp != -1)
    {
      switch(sideTemp)
      {
         case 0:
         case 1:

                bbox.minBB[0] = coord3 ;
                bbox.minBB[1] = GeomData->computeCoord(1, knots[0][0]) ;
                bbox.minBB[2] = GeomData->computeCoord(2, knots[1][0]) ;

                bbox.maxBB[0] = coord3 ;
                bbox.maxBB[1] = GeomData->computeCoord(1, knots[0][1]) ;
                bbox.maxBB[2] = GeomData->computeCoord(2, knots[1][1]) ;
                
                JacMultElem /= GeomData->getJacobian(0);

        break;

        case 2:
        case 3:

                bbox.minBB[0] = GeomData->computeCoord(0, knots[0][0]) ;
                bbox.minBB[1] = coord3 ;
                bbox.minBB[2] = GeomData->computeCoord(2, knots[1][0]) ;

                bbox.maxBB[0] = GeomData->computeCoord(0, knots[0][1]) ;
                bbox.maxBB[1] = coord3 ;
                bbox.maxBB[2] = GeomData->computeCoord(2, knots[1][1]) ;
                
                JacMultElem /= GeomData->getJacobian(1);

        break;

        case 4:
        case 5:

                bbox.minBB[0] = GeomData->computeCoord(0, knots[0][0]) ;
                bbox.minBB[1] = GeomData->computeCoord(1, knots[1][0]) ;
                bbox.minBB[2] = coord3 ;

                bbox.maxBB[0] = GeomData->computeCoord(0, knots[0][1]) ;
                bbox.maxBB[1] = GeomData->computeCoord(1, knots[1][1]) ;
                bbox.maxBB[2] = coord3 ;
                
                JacMultElem /= GeomData->getJacobian(2);

        break;

        default :

            cout << " Invalid 'sideTemp' value in AdaptiveBinarytree<DIM>::prepareData() " << endl;
        break;
      } //switch(sideTemp)
    }

    return;
}


template<int DIM>
bool AdaptiveOctree<DIM>::within(myPoint& pt)
{
  return  bbox.within(pt);
}


#ifdef _DOMAIN1D
typedef AdaptiveOctree<1>  adapOctree;
#endif


#ifdef _DOMAIN2D
typedef AdaptiveOctree<2>  adapOctree;
#endif


#ifdef _DOMAIN3D
typedef AdaptiveOctree<3>  adapOctree;
#endif


#endif

