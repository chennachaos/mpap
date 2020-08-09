
#ifndef INSERT_UNIQUE_POINT_h
#define INSERT_UNIQUE_POINT_h

#include "headersVTK.h"


int  InsertUniquePointsForAABB2D(vtkMergePoints* mergePoints, const AABB& bbTemp, vtkIdType* ptIds)
{
    vtkIdType  ind;
    double  vec[]={0.0, 0.0, 0.0};
    
    vec[0] = bbTemp.minBB[0];  vec[1] = bbTemp.minBB[1];
    ind = mergePoints->InsertUniquePoint(vec, ptIds[0]);

    vec[0] = bbTemp.maxBB[0];  vec[1] = bbTemp.minBB[1];
    ind = mergePoints->InsertUniquePoint(vec, ptIds[1]);

    vec[0] = bbTemp.maxBB[0];  vec[1] = bbTemp.maxBB[1];
    ind = mergePoints->InsertUniquePoint(vec, ptIds[2]);

    vec[0] = bbTemp.minBB[0];  vec[1] = bbTemp.maxBB[1];
    ind = mergePoints->InsertUniquePoint(vec, ptIds[3]);

    return 1;
}


int  InsertUniquePointsForAABB3D(vtkMergePoints* mergePoints, const AABB& bbTemp, vtkIdType* ptIds)
{
    vtkIdType  ind;
    double  vec[]={0.0, 0.0, 0.0};
    
    vec[0] = bbTemp.minBB[0];  vec[1] = bbTemp.minBB[1];  vec[2] = bbTemp.minBB[2];
    ind = mergePoints->InsertUniquePoint(vec, ptIds[0]);

    vec[0] = bbTemp.maxBB[0];  vec[1] = bbTemp.minBB[1];  vec[2] = bbTemp.minBB[2];
    ind = mergePoints->InsertUniquePoint(vec, ptIds[1]);

    vec[0] = bbTemp.maxBB[0];  vec[1] = bbTemp.maxBB[1];  vec[2] = bbTemp.minBB[2];
    ind = mergePoints->InsertUniquePoint(vec, ptIds[2]);

    vec[0] = bbTemp.minBB[0];  vec[1] = bbTemp.maxBB[1];  vec[2] = bbTemp.minBB[2];
    ind = mergePoints->InsertUniquePoint(vec, ptIds[3]);

    vec[0] = bbTemp.minBB[0];  vec[1] = bbTemp.minBB[1];  vec[2] = bbTemp.maxBB[2];
    ind = mergePoints->InsertUniquePoint(vec, ptIds[4]);

    vec[0] = bbTemp.maxBB[0];  vec[1] = bbTemp.minBB[1];  vec[2] = bbTemp.maxBB[2];
    ind = mergePoints->InsertUniquePoint(vec, ptIds[5]);

    vec[0] = bbTemp.maxBB[0];  vec[1] = bbTemp.maxBB[1];  vec[2] = bbTemp.maxBB[2];
    ind = mergePoints->InsertUniquePoint(vec, ptIds[6]);

    vec[0] = bbTemp.minBB[0];  vec[1] = bbTemp.maxBB[1];  vec[2] = bbTemp.maxBB[2];
    ind = mergePoints->InsertUniquePoint(vec, ptIds[7]);

    return 1;
}
              /*
              ptIds[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.minBB[2]);
              ptIds[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.minBB[2]);
              ptIds[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]);
              ptIds[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]);

              ptIds[4] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]);
              ptIds[5] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]);
              ptIds[6] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]);
              ptIds[7] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]);
              */








#endif