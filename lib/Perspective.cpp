
#include <iostream>
#include <cmath>

#include "Perspective.h"
#include "MathMatrix.h"
#include "DomainTree.h"
#include "Plot.h"


extern DomainTree domain;
extern Plot       plot;


using namespace std;


Perspective::Perspective()
{
  uX[0] = 1.;
  uX[1] = 0.;
  uX[2] = 0.;

  M1.setDim(3,3);
  M2.setDim(3,3);
   M.setDim(3,3);

  rotMode = ROTZ;

  objSurfPtr = NULL;
  nObjSurf   = 0;

  return;
}





Perspective::~Perspective(void)
{
  if (objSurfPtr != NULL) delete [] objSurfPtr;

  return;
}






void Perspective::addObjSurf(ObjectSurface *ptr)
{
  if (nObjSurf == 0)
  {
    objSurfPtr = new ObjectSurface* [1];

    objSurfPtr[0] = ptr;

    nObjSurf = 1;

    return;
  }

  ObjectSurface **tmp = new ObjectSurface* [nObjSurf+1];

  for (int i=0; i<nObjSurf; i++) tmp[i] = objSurfPtr[i];

  tmp[nObjSurf++] = ptr;

  delete [] objSurfPtr;

  objSurfPtr = tmp;

  return;
}







void Perspective::removeAllObjSurf(void)
{
  nObjSurf = 0;

  if (objSurfPtr != NULL) delete [] objSurfPtr;

  objSurfPtr = NULL;

  return;
}







void Perspective::setNewPersFlag(void)
{
  boundingBox.newPersFlag = true;

  for (int i=0; i<nObjSurf; i++) objSurfPtr[i]->newPersFlag = true;

  return;
}






float Perspective::fit(float *xmax, float *xmin)
{
  int i;

  float optViewDistanceFact = 1.5, 
        fact, xx[24], 
        picMax[2] = {-1.e12,-1.e12},
        picCoor[2];

  for (i=0; i<3; i++) 
  {
    view[i] = xmin[i] - xmax[i];
    observer[i] = xmin[i] - view[i] * (optViewDistanceFact + 1.);
  }

  reset();

  xx[ 0] = xmin[0]; xx[ 1] = xmin[1]; xx[ 2] = xmin[2];
  xx[ 3] = xmin[0]; xx[ 4] = xmin[1]; xx[ 5] = xmax[2];
  xx[ 6] = xmin[0]; xx[ 7] = xmax[1]; xx[ 8] = xmin[2];
  xx[ 9] = xmin[0]; xx[10] = xmax[1]; xx[11] = xmax[2];
  xx[12] = xmax[0]; xx[13] = xmin[1]; xx[14] = xmin[2];
  xx[15] = xmax[0]; xx[16] = xmin[1]; xx[17] = xmax[2];
  xx[18] = xmax[0]; xx[19] = xmax[1]; xx[20] = xmin[2];
  xx[21] = xmax[0]; xx[22] = xmax[1]; xx[23] = xmax[2];

  for (i=0; i<8; i++)
  {
    pictureCoor(xx+i+i+i,picCoor); 
    if (abs(picCoor[0]) > picMax[0]) { picMax[0] = abs(picCoor[0]); }
    if (abs(picCoor[1]) > picMax[1]) { picMax[1] = abs(picCoor[1]); }
  }

  fact = 0.45 / picMax[0];

  for (i=0; i<3; i++) view[i] *= fact;

  viewDotObserver = view[0] * observer[0] + view[1] * observer[1] + view[2] * observer[2];
  viewDotView     = view[0] * view[0]     + view[1] * view[1]     + view[2] * view[2];

  return picMax[1] / picMax[0];
}









void Perspective::reset(void)
{
  setNewPersFlag();

  viewDotObserver = view[0] * observer[0] + view[1] * observer[1] + view[2] * observer[2];
  viewDotView     = view[0] * view[0]     + view[1] * view[1]     + view[2] * view[2];

  float l, L = 1. / sqrt(viewDotView), g[3];

  switch (rotMode)
  {
    case ROTFREE:

    case ROTOBS:

    case ROTZ: l = sqrt(view[0] * view[0] + view[1] * view[1]);

               if (l > 1.e-8)
               {
                 l = 1. / l;

                 uX[0] = + view[1] * l;
                 uX[1] = - view[0] * l;
                 uX[2] = 0.;
               
                 uY[0] = + uX[1] * view[2] * L;
                 uY[1] = - uX[0] * view[2] * L;
                 uY[2] = + L / l;
               }
               else 
               { 
                 uX[0] = 1.;
                 uX[1] = 0.;
                 uX[2] = 0.;

                 uY[0] = 0.;
                 uY[1] = 1.;  
                 uY[2] = 0.;
               }

               break;

    case ROTX: l = sqrt(view[1] * view[1] + view[2] * view[2]);

               if (l > 1.e-8)
               {
                 l = 1. / l;

                 uX[0] = 0.;
                 uX[1] = + view[2] * l;
                 uX[2] = - view[1] * l;
               
                 uY[0] = + L / l;
                 uY[1] = + uX[2] * view[0] * L;
                 uY[2] = - uX[1] * view[0] * L;
               }
               else 
               { 
                 uX[0] = 0.;
                 uX[1] = 1.;
                 uX[2] = 0.;

                 uY[0] = 0.;
                 uY[1] = 0.;  
                 uY[2] = 1.;
               }
               break;

    case ROTY: l = sqrt(view[2] * view[2] + view[0] * view[0]);

               if (l > 1.e-8)
               {
                 l = 1. / l;

                 uX[0] = - view[2] * l;
                 uX[1] = 0.;
                 uX[2] = + view[0] * l;
               
                 uY[0] = - uX[2] * view[1] * L;
                 uY[1] = + L / l;
                 uY[2] = + uX[0] * view[1] * L;
               }
               else 
               { 
                 uX[0] = 0.;
                 uX[1] = 0.;
                 uX[2] = 1.;

                 uY[0] = 1.;
                 uY[1] = 0.;  
                 uY[2] = 0.;
               }
               break;
  }
  return;
}







void Perspective::pictureCoor(float *point, float *picCoor)
{
  float ray[3] = { observer[0] - point[0], observer[1] - point[1], observer[2] - point[2] },
        l      = view[0] * ray[0] + view[1] * ray[1] + view[2] * ray[2],
        picVec[3];

  if (abs(l) < 1.e-8) return;

  l = viewDotView / l;

  picVec[0] = ray[0] * l - view[0];
  picVec[1] = ray[1] * l - view[1];
  picVec[2] = ray[2] * l - view[2];
  
  picCoor[0] = uX[0] * picVec[0] + uX[1] * picVec[1] + uX[2] * picVec[2];
  picCoor[1] = uY[0] * picVec[0] + uY[1] * picVec[1] + uY[2] * picVec[2];

  return;
}









void Perspective::calcNew(float s1, float s2, float *dAct, int button)
{
  float xRot[3], dmy1[3], dmy2[3], dmy3[3],
        xs, xC, ys, yC, xyC, zs, zC, yzC, zxC, 
        s, c, omc, theta, S, C, fact1, fact2, dx, dy;

  int i, j;

  setNewPersFlag();

  // react to mouse input following 'pers'

  if (button == 1)   // rotate 
  {
    if (rotMode != ROTFREE)
    {
      theta = 0.78539816 * s2;

      s     = sin(theta);
      c     = cos(theta);
      omc   = 1. - c;

      theta = 0.78539816 * s1 * 1.5;

      S     = sin(theta);
      C     = cos(theta);

      switch (rotMode)
      {
        case ROTZ: fact1 = C * uX[0] + S * uX[1];
                   fact2 = S * uX[0] - C * uX[1];
              
                   M(1,1) = c * C + omc * uX[0] * fact1;
                   M(1,2) = c * S + omc * uX[1] * fact1;
                   M(1,3) = s * fact2;

                   M(2,1) = - c * S - omc * uX[0] * fact2;
                   M(2,2) = c * C - omc * uX[1] * fact2;
                   M(2,3) = s * fact1;

                   M(3,1) = s * uX[1];
                   M(3,2) = - s * uX[0];
                   M(3,3) = c;

                   break;

        case ROTX: S = -S; s = -s;

                   fact1 = C * uX[2] + S * uX[1];
                   fact2 = S * uX[2] - C * uX[1];
              
                   M(1,1) = c;
                   M(1,2) = - s * uX[2];
                   M(1,3) = s * uX[1];
                       
                   M(2,1) = s * fact1;
                   M(2,2) = c * C - omc * uX[1] * fact2;
                   M(2,3) = - c * S - omc * uX[2] * fact2;

                   M(3,1) = s * fact2;
                   M(3,2) = c * S + omc * uX[1] * fact1;
                   M(3,3) = c * C + omc * uX[2] * fact1;

                   break;

        case ROTY: s = -s; S = -S;

                   fact1 = C * uX[0] + S * uX[2];
                   fact2 = S * uX[0] - C * uX[2];
              
                   M(1,1) = c * C + omc * uX[0] * fact1;
                   M(1,2) = s * fact2;
                   M(1,3) = c * S + omc * uX[2] * fact1;

                   M(2,1) = s * uX[2];
                   M(2,2) = c;
                   M(2,3) = - s * uX[0];

                   M(3,1) = - c * S - omc * uX[0] * fact2;
                   M(3,2) = s * fact1;
                   M(3,3) = c * C - omc * uX[2] * fact2;

                   break;
      }
    }
    else
    {
      theta = 0.78539816 * s2;

      c = cos(theta); s = -sin(theta); C = 1-c;

      xs  = uX[0]*s;   ys = uX[1]*s;   zs = uX[2]*s;
      xC  = uX[0]*C;   yC = uX[1]*C;   zC = uX[2]*C;
      xyC = uX[0]*yC; yzC = uX[1]*zC; zxC = uX[2]*xC;

      M1(1,1) = uX[0]*xC + c;
      M1(1,2) = xyC - zs;
      M1(1,3) = zxC + ys;

      M1(2,1) = xyC + zs;
      M1(2,2) = uX[1]*yC + c;
      M1(2,3) = yzC - xs;

      M1(3,1) = zxC - ys;
      M1(3,2) = yzC + xs;
      M1(3,3) = uX[2]*zC + c;

      theta = 0.78539816 * s1 * 1.5;

      c = cos(theta); s = -sin(theta); C = 1-c;

      xs  = uY[0]*s;   ys = uY[1]*s;   zs = uY[2]*s;
      xC  = uY[0]*C;   yC = uY[1]*C;   zC = uY[2]*C;
      xyC = uY[0]*yC; yzC = uY[1]*zC; zxC = uY[2]*xC;

      M2(1,1) = uY[0]*xC + c;
      M2(1,2) = xyC - zs;
      M2(1,3) = zxC + ys;

      M2(2,1) = xyC + zs;
      M2(2,2) = uY[1]*yC + c;
      M2(2,3) = yzC - xs;

      M2(3,1) = zxC - ys;
      M2(3,2) = yzC + xs;
      M2(3,3) = uY[2]*zC + c;

      for (i=1; i<4; i++)
        for (j=1; j<4; j++)
          M(i,j) = M1(i,1)*M2(1,j) + M1(i,2)*M2(2,j) + M1(i,3)*M2(3,j);
    }

    fact1 = (view[0] * xCntr[0] + view[1] * xCntr[1] + view[2] * xCntr[2] - viewDotObserver) 
               / viewDotView;

    fact2 = 0.;

    for (i=0; i<3; i++)
    {
      dmy1[i] = 0.;
      dmy2[i] = 0.;
      dmy3[i] = 0.;
      for (j=0; j<3; j++) 
      {
        dmy1[i] += M(1+i,1+j) * view[j];
        dmy2[i] += M(1+i,1+j) * uX[j];
        dmy3[i] += M(1+i,1+j) * uY[j];
      }
      xRot[i] = observer[i] + view[i] * fact1;
      fact2 += (observer[i] - xRot[i]) * (observer[i] - xRot[i]);
    }
    fact2 = sqrt(fact2 / viewDotView);
    viewDotObserver = 0.;
    for (i=0; i<3; i++) 
    {
      observer[i] = xRot[i] - dmy1[i] * fact2;
      view[i] = dmy1[i];
      uX[i]   = dmy2[i];
      uY[i]   = dmy3[i];
      viewDotObserver += view[i] * observer[i];
    }
    return;
  }

  // translate

  if (button == 3)
  {
    dx = s1 * dAct[0];
    dy = s2 * dAct[1];
   
    viewDotObserver = 0.; 
    for (i=0; i<3; i++) 
    {
      observer[i] -= (dx * uX[i] - dy * uY[i]);
      viewDotObserver += view[i] * observer[i];
    }

    return;
  }

  // zoom

  if (button == 4) fact1 = 1.05; else fact1 = 0.95;

  if (zoomMode == VIEW)
  {
    for (i=0; i<3; i++) view[i] *= fact1;    
    viewDotView *= (fact1 * fact1);
    viewDotObserver *= fact1;
  }
  else
  {
    fact1 = ((xCntr[0]-observer[0])*view[0]
           + (xCntr[1]-observer[1])*view[1]
           + (xCntr[2]-observer[2])*view[2]) / viewDotView * (1.-fact1);
    for (i=0; i<3; i++) observer[i] -= fact1 * view[i];
    viewDotObserver = observer[0]*view[0] + observer[1]*view[1] + observer[2]*view[2];
  }
  return;
}







void Perspective::plotChangePersObjSurf(bool backFaceFlag, bool defmFlag)
{
  changePersObjSurfPtr->draw(8,3,false,backFaceFlag,defmFlag);

  drawCoorAxes();

  return;
}







void Perspective::print(float *dDes)
{
  // to be invoked after pers,mous  (unix: grpDrawingAreaMouseInput)

  COUT << "new perspective:\n\n";

  printf("          pers,,,,2,,,,,,%d,%9.3g,%9.3g,%9.3g &\n",
            (int)rotMode+1,observer[0],observer[1],observer[2]);
  printf("                           %9.3g,%9.3g,%9.3g &\n",view[0],view[1],view[2]);
  printf("                           %9.3g,%9.3g,%9.3g &\n",uX[0],uX[1],uX[2]);
  printf("                           %9.3g,%9.3g,%9.3g &\n",uY[0],uY[1],uY[2]);
  printf("                           %9.3g\n\n",dDes[1]/dDes[0]);

  return;
}







float Perspective::zDepth(float *x)
{
  return view[0]*x[0] + view[1]*x[1] + view[2]*x[2] - viewDotObserver;
}






void Perspective::prepare(Vector<double> &p, int i0)
{
  int i;

  for (i=0; i<3; i++)
  {
    observer[i] = (float)p[i0+  i];
    view[i]     = (float)p[i0+3+i];
    uX[i]       = (float)p[i0+6+i];
    uY[i]       = (float)p[i0+9+i];
  }

  float fact = 1. / sqrt(uX[0]*uX[0] + uX[1]*uX[1] + uX[2]*uX[2]);
  uX[0] *= fact; uX[1] *= fact; uX[2] *= fact;

  fact = 1. / sqrt(uY[0]*uY[0] + uY[1]*uY[1] + uY[2]*uY[2]);
  uY[0] *= fact; uY[1] *= fact; uY[2] *= fact;

  viewDotObserver = view[0] * observer[0] + view[1] * observer[1] + view[2] * observer[2];
  viewDotView     = view[0] * view[0]     + view[1] * view[1]     + view[2] * view[2];

  return;
}






void Perspective::searchCone(double *cone, double *axis, float *picCoor, float r)
{
  axis[0] = (double)(view[0] + uX[0] * picCoor[0] + uY[0] * picCoor[1]);
  axis[1] = (double)(view[1] + uX[1] * picCoor[0] + uY[1] * picCoor[1]);
  axis[2] = (double)(view[2] + uX[2] * picCoor[0] + uY[2] * picCoor[1]);

  *cone   = (double)(r/sqrt(axis[0]*axis[0]+axis[1]*axis[1]+axis[2]*axis[2]));

  return;
}






void Perspective::generateBoundingBox(double *xmn, double *xmx)
{
  boundingBox.generateBoundingBox(xmn,xmx);

  changePersObjSurfPtr = &boundingBox;

  return;
}







void Perspective::drawCoorAxes(void)
{
  int xscr0[2] = { 60, plot.hPix - 60 }, xscr1[2], i;

  float picCoor0[2], picCoor1[2], x0[3], x1[3], fact, dx[2];

  char *tmp[3] = {"x","y","z"};

  for (i=0; i<3; i++)
  {
    x0[i] = observer[i]+view[i];
    x1[i] = x0[i];
  }

  for (i=0; i<3; i++)
  {
    x1[i] += 1.;

    pictureCoor(x1,picCoor1);

    dx[0] = picCoor1[0];
    dx[1] = picCoor1[1];

    fact = 40. * plot.dAct[0] / (plot.w * sqrt(dx[0]*dx[0]+dx[1]*dx[1]));

    picCoor1[0] = dx[0] * fact;
    picCoor1[1] = dx[1] * fact;

    plot.xy2D(picCoor1,xscr1);

    xscr1[0] += 60 - (int)(plot.w * .5);
    xscr1[1] -= 60 - (int)(plot.h * .5);

    essGrpDrawLine(xscr0[0],xscr0[1],xscr1[0],xscr1[1]);

    xscr1[0] += (int)((float)(xscr1[0] - xscr0[0]) * .3);
    xscr1[1] += (int)((float)(xscr1[1] - xscr0[1]) * .3);

    essGrpPutText(xscr1,tmp[i],5);

    x1[i] -= 1.;
  }

  return;
}


