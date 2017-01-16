
#ifndef incl_Plot_h
#define incl_Plot_h


#include "MathVector.h"
#include "Perspective.h"
#include "ColoursEnum.h"
#include "FunctionsEssGrp.h"
#include "Files.h"


extern bool  noGUI;
extern Files files;


class Plot
{
  public:

    Plot(void);
    
    ~Plot();
	  
    int   dim, nNamedColours, hPix, wPix, currStdColour;
	  
    float w, h, wdh, hdw; 
	   
    Perspective perspective;

    float wFact1D, hFact1D, wFact1DPS, hFact1DPS, x1D, u1D, dx1D, du1D;
	   
    float wFactAct, hFactAct, wFactDesPS, hFactDesPS;
                           
    float psFact, psBoundingBoxArea;
    
    float dAct[2], dDes[2], x0Act[2], x0Des[2];

    float YPixPerXPix, selectSearchRadius;
    
    bool  selectBox1, suppressCopyPixmap, vert1D, psOpen;
    
    int   sBx, sBy;
    
    float *polyPnt;

    int   *iPolyPnt;

    void  calcAspectRatio(int, int, int, int); 
    
    bool  operator!(void);
    
    void  fit(double *, double *, int, double fact = 0.4);
    
    void  adjustToNewSize(void);
    
    void  resetCoor(void);	  
    
    void  reset(void);	  
    
    void  setColour(int, bool stdCol = true);

    void  poly(int);

    void point(float*, float dpt = -1., int pn = 0); 

    void point(double*, double dpt = -1., int pn = 0);

    void dash(double*, double, int);

    void dash(float*, float, int);

    void line(float*, float*);
 
    void line(double*, double*);

    void fillCircle(double*, double&);

    void fillCircle(float*, float&);

    void circle(double *, double &);

    void circle(float*, float&);

    void triangle(float*, float&);

    void triangle(double*, double&);
    
    void square(double *, double &);

    void square(float*, float&);
    
    void putText(double*, char*, int, bool BGFlag = false);

    void putText(float*, char*, int, bool BGFlag = false);

    void arrow(double*, double *, double &);

    void arrow(float*, float*, float&);

    void triangleContourPlot(double*, double*, double*, 
		             double,  double,  double,
		             double,  double,  int);

    void triangleContourPlot(float*, float*, float*, 
		             float,  float,  float,
		             float,  float,  int);

    void contourPlotLegend(double, double, double, double, bool, int);

    void adjust1DSettings(bool, double, double, bool, double, double, bool);

    bool psOpenFile(void); 

    bool psCloseFile(void); 

    void drawDesFrame(void);

    void wipe(void);

    float dPt(void) { return (dAct[0] + dAct[1]) * .0055; }

    template<typename Type> void xy3D(Type*, int*);

    template<typename Type> void xy1D(Type*, int*);

    template<typename Type> void xy2D(Type*, int*);

    template<typename Type> void xy2DInv(Type*, int*);

    template<typename Type> void xy1DPS(Type*, int*);

    template<typename Type> void xy2DPS(Type*, int*);

    template<typename Type> void xy3DPS(Type*, int*);

  private:

    bool clipLine2D(int*, int*);

    bool clipLine2DPS(int*, int*);

    bool clipPoint2D(int*);
};






template<typename Type> inline void Plot::xy1D(Type *x, int *ix)
{
  Type          x00 = x1D, x01 = u1D;
  if (vert1D) { x00 = u1D; x01 = x1D; }
  
  ix[0] =        roundToInt((x[0]-x00) * wFact1D);
  ix[1] = hPix - roundToInt((x[1]-x01) * hFact1D);
  
  return;
}

template<typename Type> inline void Plot::xy1DPS(Type *x, int *ix)
{
  Type          x00 = x1D, x01 = u1D;
  if (vert1D) { x00 = u1D; x01 = x1D; }
  
  ix[0] = roundToInt((x[0]-x00) * wFact1DPS);
  ix[1] = roundToInt((x[1]-x01) * hFact1DPS);

  return;
}

template<typename Type> inline void Plot::xy2D(Type *x, int *ix)
{
  ix[0] =        roundToInt((x[0]-x0Act[0]) * wFactAct);
  ix[1] = hPix - roundToInt((x[1]-x0Act[1]) * hFactAct);

  return;
}

template<typename Type> inline void Plot::xy2DPS(Type *x, int *ix)
{
  ix[0] = roundToInt((x[0]-x0Des[0]) * wFactDesPS);
  ix[1] = roundToInt((x[1]-x0Des[1]) * hFactDesPS);

  return;
}

template<typename Type> inline void Plot::xy2DInv(Type *x, int *ix)
{
  x[0] = (float)   ix[0]     / wFactAct + x0Act[0];
  x[1] = (float)(hPix-ix[1]) / hFactAct + x0Act[1];

  return;
}

template<typename Type> inline void Plot::xy3D(Type *x, int *ix)
{
  float picCoor[2];

  perspective.pictureCoor(x,picCoor);

  xy2D(picCoor,ix);

  return;
}

template<typename Type> inline void Plot::xy3DPS(Type *x, int *ix)
{
  float picCoor[2];

  perspective.pictureCoor(x,picCoor);

  xy2DPS(picCoor,ix);

  return;
}

inline void Plot::triangleContourPlot(double *x1, double *x2, double *x3, 
                                      double u1, double u2, double u3,
                                      double umn, double umx, int nd)
{
  float x1f[3], x2f[3], x3f[3], 
        u1f = (float)u1, u2f = (float)u2, u3f = (float)u3, 
        umnf = (float)umn, umxf = (float)umx;

  for (int i=0; i<dim; i++) 
  {
    x1f[i] = (float)x1[i];
    x2f[i] = (float)x2[i];
    x3f[i] = (float)x3[i];
  }
  triangleContourPlot(x1f,x2f,x3f,u1f,u2f,u3f,umnf,umxf,nd);

  return;
}

inline void Plot::point(double *x, double dpt, int pn)
{ 
  float xf[3] = { (float)x[0], (float)x[1], (float)x[2] }, dptf = (float) dpt;
  point(xf,dptf,pn); 
  return; 
}

inline void Plot::dash(double*x, double d, int i)
{ 
  float xf[3] = { (float)x[0], (float)x[1], (float)x[2] }, df = (float)d;
  dash(xf,df,i); 
  return; 
}

inline void Plot::line(double*x1, double*x2)
{ 
  float x1f[3] = { (float)x1[0], (float)x1[1], (float)x1[2] },
        x2f[3] = { (float)x2[0], (float)x2[1], (float)x2[2] };
  line(x1f,x2f); 
  return; 
}

inline void Plot::fillCircle(double *x, double &d)
{
  float xf[3] = { (float)x[0], (float)x[1], (float)x[2] }, df = (float)d;
  fillCircle(xf,df);
  d = (double)df;
  return;
}

inline void Plot::circle(double *x, double &d)
{
  float xf[3] = { (float)x[0], (float)x[1], (float)x[2] }, df = (float)d;
  circle(xf,df);
  d = (double)df;
  return;
}

inline void Plot::triangle(double *x, double &d)
{
  float xf[3] = { (float)x[0], (float)x[1], (float)x[2] }, df = (float)d;
  triangle(xf,df);
  d = (double)df;
  return;
}

inline void Plot::square(double *x, double &d)
{
  float xf[3] = { (float)x[0], (float)x[1], (float)x[2] }, df = (float)d;
  square(xf,df);
  d = (double)df;
  return;
}

inline void Plot::putText(double *x, char *s, int p, bool BGFlag)
{
  float xf[3] = { (float)x[0], (float)x[1], (float)x[2] };
  putText(xf,s,p,BGFlag);
  return;
}

inline void Plot::arrow(double *x1, double *x2, double &d)
{
  float x1f[3] = { (float)x1[0], (float)x1[1], (float)x1[2] },
        x2f[3] = { (float)x2[0], (float)x2[1], (float)x2[2] }, df = (float)d;
  arrow(x1f,x2f,df);
  d = (double)df;
  return;
}




#endif

