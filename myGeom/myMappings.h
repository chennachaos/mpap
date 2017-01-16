
#ifndef  MY_MAPPINGS_H
#define  MY_MAPPINGS_H

#include "headersEigen.h"


// global variables for HBSpline methods

  // go generate cell edges in 2D from corners of axis-aligned bounding box 'AABB'
static  int verts[4][2] ={ {0,1}, {1,3}, {3,2}, {2,0} };
  //int verts[4][2] ={ {0,1}, {1,2}, {2,3}, {3,0} };


  // for mapping 2D faces in 3D to global 3D coordinates
static  int coord_map_HBS_2Dto3D[6][2] ={ {1,2}, {1,2}, {0,2}, {0,2}, {0,1}, {0,1} };
  //int bbox_map[3] = {};

  
  
void map2DPointTo3DPoint(int side, myPoint& ptTemp, double val3)
{
    switch(side)
    {
        case 0:
        case 1:
                ptTemp[2] = ptTemp[1] ;
                ptTemp[1] = ptTemp[0] ;
                ptTemp[0] = val3 ;

        break;

        case 2:
        case 3:

                ptTemp[0] = ptTemp[0] ;
                ptTemp[2] = ptTemp[1] ;
                ptTemp[1] = val3 ;

        break;

        case 4:
        case 5:

                ptTemp[0] = ptTemp[0] ;
                ptTemp[1] = ptTemp[1] ;
                ptTemp[2] = val3 ;

        break;

        default :

            cout << " Invalid 'side' in map2DPointTo3DPoint in myMappings.h " << endl;
        break;
    } //switch(side)
  return;
}
  
  
  
  
  
  
  
  
  
  
  
#endif