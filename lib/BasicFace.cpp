
#include "BasicFace.h"
#include "Plot.h"
#include "FunctionsEssGrp.h"
#include "MathGeom.h"
#include "Files.h"


extern Plot plot;
extern Files files;


static const unsigned int pow2[8] = {1,2,4,8,16,32,64,128};



BasicFace::BasicFace(void)
{
  return;
}






BasicFace::~BasicFace()
{
  return;
}





void BasicFace::calcNormal(float *x)
{
  float *x1  = &(x[ix[0]+ix[0]+ix[0]]),
         *x2  = &(x[ix[1]+ix[1]+ix[1]]),
         *x3  = &(x[ix[2]+ix[2]+ix[2]]),
         a[3] = { x2[0] - x1[0], x2[1] - x1[1], x2[2] - x1[2] },
         b[3] = { x3[0] - x2[0], x3[1] - x2[1], x3[2] - x2[2] };

  normal[0] = a[1] * b[2] - a[2] * b[1];
  normal[1] = a[2] * b[0] - a[0] * b[2];
  normal[2] = a[0] * b[1] - a[1] * b[0];

  return;
}







bool BasicFace::facingOK(float *observer, float *x)
{
  int i = ix[0] + ix[0] + ix[0];

  if ( (observer[0]-x[i  ])*normal[0]
      +(observer[1]-x[i+1])*normal[1]
      +(observer[2]-x[i+2])*normal[2] < 0.)  return false;

  else return true;
}






void BasicFace::simplePaint(int colFace, int colEdge, int *scrx, float *picx, bool edgesOnly)
{
  int i, iPolyPnt[10];

  i = ix[0] + ix[0];
  iPolyPnt[0] = scrx[i++];
  iPolyPnt[1] = scrx[i];

  i = ix[1] + ix[1];
  iPolyPnt[2] = scrx[i++];
  iPolyPnt[3] = scrx[i];

  i = ix[2] + ix[2];
  iPolyPnt[4] = scrx[i++];
  iPolyPnt[5] = scrx[i];

  if (!edgesOnly) { essGrpSetColour(colFace,true); essGrpFillPoly(iPolyPnt,3); }

  essGrpSetColour(colEdge,true);

  if (edgeBits & 1) essGrpDrawLine(iPolyPnt[0],iPolyPnt[1],iPolyPnt[2],iPolyPnt[3]);
  if (edgeBits & 2) essGrpDrawLine(iPolyPnt[2],iPolyPnt[3],iPolyPnt[4],iPolyPnt[5]);
  if (edgeBits & 4) essGrpDrawLine(iPolyPnt[4],iPolyPnt[5],iPolyPnt[0],iPolyPnt[1]);
 
  if (!plot.psOpen) return;

  plot.xy2DPS(picx+ix[0]+ix[0],iPolyPnt);
  plot.xy2DPS(picx+ix[1]+ix[1],iPolyPnt+2);
  plot.xy2DPS(picx+ix[2]+ix[2],iPolyPnt+4);

  if (!edgesOnly)
  {
    files.Pout << 'C' << colFace+1 << "\n";

    files.Pout << iPolyPnt[0] << " " << iPolyPnt[1] << " m ";
    for (i=1; i<3; i++) files.Pout << iPolyPnt[i+i] - iPolyPnt[i+i-2] << " " 
                                   << iPolyPnt[i+i+1] - iPolyPnt[i+i-1] << " v ";
    files.Pout << "f\n";
  }

  files.Pout << 'C' << colEdge+1 << "\n";

  if (edgeBits & 1)
    files.Pout << iPolyPnt[0]             << " " << iPolyPnt[1]             << " m " 
               << iPolyPnt[2]-iPolyPnt[0] << " " << iPolyPnt[3]-iPolyPnt[1] << " V\n";

  if (edgeBits & 2)
    files.Pout << iPolyPnt[2]             << " " << iPolyPnt[3]             << " m " 
               << iPolyPnt[4]-iPolyPnt[2] << " " << iPolyPnt[5]-iPolyPnt[3] << " V\n";

  if (edgeBits & 4)
    files.Pout << iPolyPnt[4]             << " " << iPolyPnt[5]             << " m " 
               << iPolyPnt[0]-iPolyPnt[4] << " " << iPolyPnt[1]-iPolyPnt[5] << " V\n";
  return;
}






float BasicFace::calcZDepth(float *x)
{
  int i   = ix[0] + ix[0] + ix[0],
      ii  = ix[1] + ix[1] + ix[1],
      iii = ix[2] + ix[2] + ix[2];

  float xh[3] = { (x[i++] + x[ii++] + x[iii++])*.33333333333333,
                  (x[i++] + x[ii++] + x[iii++])*.33333333333333,
                  (x[i  ] + x[ii  ] + x[iii  ])*.33333333333333 };
 
  return plot.perspective.zDepth(xh);
}





bool BasicFace::frontContains(float *PX, float *px)
{
  if (pointInTriangle2D(PX+ix[0]+ix[0],PX+ix[1]+ix[1],PX+ix[2]+ix[2],px)) return true;  

  return false;
}





bool BasicFace::backContains(float *PX, float *px)
{
  if (pointInTriangle2D(PX+ix[2]+ix[2],PX+ix[1]+ix[1],PX+ix[0]+ix[0],px)) return true;  

  return false;
}





bool BasicFace::frontContains(int *SX, int *sx)
{
  if (pointInTriangle2D(SX+ix[2]+ix[2],SX+ix[1]+ix[1],SX+ix[0]+ix[0],sx)) return true;  

  return false;
}





bool BasicFace::backContains(int *SX, int *sx)
{
  if (pointInTriangle2D(SX+ix[0]+ix[0],SX+ix[1]+ix[1],SX+ix[2]+ix[2],sx)) return true;  

  return false;
}





void BasicFace::contourPlot(int nCol, float *X, float *U, float umn, float umx, 
                            int *SX, float *PX, bool backFaceFlg)
{
  float u1 = U[ix[0]],
        u2 = U[ix[1]],
        u3 = U[ix[2]],
        *x1 = X + ix[0]+ix[0]+ix[0],
        *x2 = X + ix[1]+ix[1]+ix[1],
        *x3 = X + ix[2]+ix[2]+ix[2];

  if (backFaceFlg) plot.triangleContourPlot(x3, x2, x1, u3, u2, u1, umn, umx, nCol);

  else             plot.triangleContourPlot(x1, x2, x3, u1, u2, u3, umn, umx, nCol);

  int ix0 = ix[0] + ix[0],
      ix1 = ix[1] + ix[1],
      ix2 = ix[2] + ix[2];

  essGrpSetColour(WHITE,true);

  if (edgeBits &  8) essGrpDrawLine(SX[ix0], SX[ix0+1], SX[ix1], SX[ix1+1]);
  if (edgeBits & 16) essGrpDrawLine(SX[ix1], SX[ix1+1], SX[ix2], SX[ix2+1]);
  if (edgeBits & 32) essGrpDrawLine(SX[ix2], SX[ix2+1], SX[ix0], SX[ix0+1]);

  if (!plot.psOpen) return;

  int iPolyPnt[10];

  plot.xy2DPS(PX+ix0,iPolyPnt);
  plot.xy2DPS(PX+ix1,iPolyPnt+2);
  plot.xy2DPS(PX+ix2,iPolyPnt+4);

  files.Pout << 'C' << WHITE+1 << "\n";

  if (edgeBits & 8)
    files.Pout << iPolyPnt[0]             << " " << iPolyPnt[1]             << " m " 
               << iPolyPnt[2]-iPolyPnt[0] << " " << iPolyPnt[3]-iPolyPnt[1] << " V\n";

  if (edgeBits & 16)
    files.Pout << iPolyPnt[2]             << " " << iPolyPnt[3]             << " m " 
               << iPolyPnt[4]-iPolyPnt[2] << " " << iPolyPnt[5]-iPolyPnt[3] << " V\n";

  if (edgeBits & 32)
    files.Pout << iPolyPnt[4]             << " " << iPolyPnt[5]             << " m " 
               << iPolyPnt[0]-iPolyPnt[4] << " " << iPolyPnt[1]-iPolyPnt[5] << " V\n";

  return;
}




/*
void BasicFace::backContourPlot(int nCol, float *X, float *U, float umn, float umx,
                                int *SX, float *PX)
{
  float u1 = U[ix[0]],
        u2 = U[ix[1]],
        u3 = U[ix[2]],
        *x1 = X + ix[0]+ix[0]+ix[0],
        *x2 = X + ix[1]+ix[1]+ix[1],
        *x3 = X + ix[2]+ix[2]+ix[2];

  plot.triangleContourPlot(x1, x2, x3, u1, u2, u3, umn, umx, nCol);

  int ix0 = ix[0] + ix[0],
      ix1 = ix[1] + ix[1],
      ix2 = ix[2] + ix[2];

  plot.setColour(WHITE);

  if (edgeBits &  8) essGrpDrawLine(SX[ix0], SX[ix0+1], SX[ix1], SX[ix1+1]);
  if (edgeBits & 16) essGrpDrawLine(SX[ix1], SX[ix1+1], SX[ix2], SX[ix2+1]);
  if (edgeBits & 32) essGrpDrawLine(SX[ix2], SX[ix2+1], SX[ix0], SX[ix0+1]);

  if (!plot.psOpen) return;

  int iPolyPnt[10];

  plot.xy2DPS(PX+ix0,iPolyPnt);
  plot.xy2DPS(PX+ix1,iPolyPnt+2);
  plot.xy2DPS(PX+ix2,iPolyPnt+4);

  files.Pout << 'C' << WHITE+1 << "\n";

  if (edgeBits & 8)
    files.Pout << iPolyPnt[0]             << " " << iPolyPnt[1]             << " m " 
               << iPolyPnt[2]-iPolyPnt[0] << " " << iPolyPnt[3]-iPolyPnt[1] << " V\n";

  if (edgeBits & 16)
    files.Pout << iPolyPnt[2]             << " " << iPolyPnt[3]             << " m " 
               << iPolyPnt[4]-iPolyPnt[2] << " " << iPolyPnt[5]-iPolyPnt[3] << " V\n";

  if (edgeBits & 32)
    files.Pout << iPolyPnt[4]             << " " << iPolyPnt[5]             << " m " 
               << iPolyPnt[0]-iPolyPnt[4] << " " << iPolyPnt[1]-iPolyPnt[5] << " V\n";

  return;
}
*/

