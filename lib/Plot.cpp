
#include <iostream>

#include "Plot.h"
#include "FunctionsEssGrp.h"
#include "FunctionsProgram.h"
#include "Definitions.h"
#include "MathBasic.h"
#include "MyString.h"
#include "Files.h"
#include "MpapTime.h"


extern bool     noGUI;
extern bool     aspectRatioCorr;
extern Files    files;
extern MpapTime mpapTime;


using namespace std;



Plot::Plot(void)
{
  dim = 0;

  reset();

  hPix = - 1;
  
  nNamedColours = 0;

  psOpen = false;

  psFact = 10.;
 
  psBoundingBoxArea = 500.*300.; 

  YPixPerXPix = 1.;

  w = 400;
  h = 300;
  
  polyPnt  = new float [200];
  iPolyPnt = new int   [200];

  return;
}




Plot::~Plot()
{
  delete [] polyPnt;
  delete [] iPolyPnt;

  return;
}




void Plot::calcAspectRatio(int ww, int wMM, int hh, int hMM)
{
  YPixPerXPix = (float)hh/(float)hMM / ((float)ww/(float)wMM); 
  
  return;
}




void Plot::fit(double *xMn, double *xMx, int domDim, double factIn)
{
  float dx, 
        xmn[3] = { (float)xMn[0],  (float)xMn[1], (float)xMn[2] },
        xmx[3] = { (float)xMx[0],  (float)xMx[1], (float)xMx[2] },
        fact = (float) factIn;

  dim = domDim;
  
  switch (dim)
  {
    case 1: // nothing to be done here
	    
	    break;

    case 2: dx = fact * 0.5 * (xmx[0] - xmn[0] + xmx[1] - xmn[1]);

            dDes[0] = xmx[0] - xmn[0] + dx;
	    dDes[1] = xmx[1] - xmn[1] + dx;
	    
	    x0Des[0] = 0.5 * (xmx[0] + xmn[0] - dDes[0]);
	    x0Des[1] = 0.5 * (xmx[1] + xmn[1] - dDes[1]);
	    
	    break;
	    
    case 3: dDes[0] = 1.0;
            dDes[1] = perspective.fit(xmx,xmn);

            x0Des[0] = -.5;
            x0Des[1] = -.5 * dDes[1];

	    break;
  }

  adjustToNewSize();
  
  return;
}





void Plot::adjustToNewSize(void)
{
  essGrpGetPlotAreaGeom();
    
  if (aspectRatioCorr) h = h / YPixPerXPix;
    
  wdh = w / h;
  hdw = h / w;
  
  switch (dim)
  {
    case 1: // nothing to be done
	    
	    break;

    case 3: perspective.setNewPersFlag();

    case 2: x0Act[0] = x0Des[0];
            x0Act[1] = x0Des[1];
  
            dAct[0] = dDes[0];
            dAct[1] = dDes[1];
			     
	    if (dAct[0]/dAct[1] > wdh) 
	         { x0Act[1] += (dAct[1]  - hdw * dAct[0]) * 0.5;  dAct[1] = hdw * dAct[0]; }
	    else { x0Act[0] += (dAct[0]  - wdh * dAct[1]) * 0.5;  dAct[0] = wdh * dAct[1]; }
	    
	    wFactAct = w / dAct[0];
	    hFactAct = h / dAct[1];

	    wFactDesPS = sqrt(psBoundingBoxArea/(dDes[0]*dDes[1])) * psFact;
	    hFactDesPS = wFactDesPS; 
	    
	    if (aspectRatioCorr) hFactAct = hFactAct * YPixPerXPix;
	    
	    break;
  }

  selectSearchRadius = (dAct[0] + dAct[1]) / (w + h) * 5;

  //cout << dAct[0] << "," << dAct[1] << ";" << w << ";" << h << "->" << selectSearchRadius << "\n";

  return;
}




bool Plot::operator!(void)
{
  if (dDes[0] < 0) return true;  return false;
}





void Plot::reset(void)
{
  resetCoor();
  
  selectBox1 = false;

  suppressCopyPixmap = false;
 
  dim = 0;
 
  perspective.rotMode = ROTZ;

  perspective.removeAllObjSurf();
 
  return;
}





void Plot::resetCoor(void)
{
  dDes[0] = -1.;
  
  return;
}





void Plot::setColour(int col, bool stdCol)
{
  if (stdCol) currStdColour = col;

  essGrpSetColour(col, stdCol);

  if (!psOpen) return;
    
  if (stdCol) files.Pout << 'C' << col+1 << "\n";
  else        files.Pout << col << " c\n"; 
	  
  return;
}





void Plot::poly(int npt)
{
  int i;

  float picCoor[2];
 
  if (dim == 1)      
  {       }
  else if (dim == 2) 
  { 
    for (i=0; i<npt; i++) xy2D(polyPnt+i+i+i,iPolyPnt+i+i); 
  }
  else               
  { 
    for (i=0; i<npt; i++) 
    { 
      perspective.pictureCoor(polyPnt+i+i+i,picCoor); 
      xy2D(picCoor,iPolyPnt+i+i); 
    }
  }

  essGrpFillPoly(&(iPolyPnt[0]),npt);
  
  if (!psOpen) return;
 
 
  if (dim == 1)      
  {       }
  else if (dim == 2) 
  { 
    for (i=0; i<npt; i++) xy2DPS(&(polyPnt[i*3]),&(iPolyPnt[i*2])); 
  }
  else
  {
    for (i=0; i<npt; i++) 
    { 
      perspective.pictureCoor(&(polyPnt[i*3]),picCoor);
      xy2DPS(picCoor,&(iPolyPnt[i*2]));
    }
  }
          
  files.Pout << iPolyPnt[0] << " " << iPolyPnt[1] << " m ";
  for (i=1; i<npt; i++) files.Pout << iPolyPnt[i*2] - iPolyPnt[(i-1)*2] << " " 
                                   << iPolyPnt[i*2+1] - iPolyPnt[i*2-1] << " v ";
  files.Pout << "f\n";
	    
  return;
}






void Plot::triangleContourPlot(float *x1, float *x2, float *x3, 
		               float u1, float u2, float u3,
			       float umn, float umx, int nd)
{
  int ix[10], imn, imx, imm, mxCol = 511, col;
	
  float u[4], uc, du, Du;

  // abort, if the triangle lies outside the visible window
  
  if (dim == 3)
  {
    // this check is done by ObjectSurface

    //perspective.pictureCoor(x1,u); 
    //if (u[0]>x0Act[0]+dAct[0] || u[0]<x0Act[0] || u[1]>x0Act[1]+dAct[1] || u[1]<x0Act[1])
    //{
    //  perspective.pictureCoor(x2,u); 
    //  if (u[0]>x0Act[0]+dAct[0] || u[0]<x0Act[0] || u[1]>x0Act[1]+dAct[1] || u[1]<x0Act[1])
    //  {
    //    perspective.pictureCoor(x3,u); 
    //    if (u[0]>x0Act[0]+dAct[0] || u[0]<x0Act[0] || u[1]>x0Act[1]+dAct[1] || u[1]<x0Act[1])
    //      return;
    //  }
    //}
  }
  else if (dim == 2)
  {
    u[0] = x0Act[0] - 0.1*dAct[0];
    u[1] = x0Act[1] - 0.1*dAct[1];
    u[2] = x0Act[0] + 1.1*dAct[0];
    u[3] = x0Act[1] + 1.1*dAct[1];

    if (x1[0] > u[2] || x1[0] < u[0] || x1[1] > u[3] || x1[1] < u[1])
    {
      if (x2[0] > u[2] || x2[0] < u[0] || x2[1] > u[3] || x2[1] < u[1])
      {
        if (x3[0] > u[2] || x3[0] < u[0] || x3[1] > u[3] || x3[1] < u[1])
          return;
      }
    }
  }
  else
  {
    prgWarning(1,"Plot::triangleContourPlot","dimension must be 2 or 3!"); return;
  }
 
  // rearrange u

  u[0] = u1;
  u[1] = u2;
  u[2] = u3;

  // find the nodes with the smallest, the largest and the intermedium u
	
  imn = 0; if (u[imn] > u[1]) imn = 1; if (u[imn] > u[2]) imn = 2;
  imx = 0; if (u[imx] < u[1]) imx = 1; if (u[imx] < u[2]) imx = 2;

  imm = 0; if (imm==imn || imm==imx) imm = 1; if (imm==imn || imm==imx) imm = 2;
  
  // calculate du

  du = (umx - umn) / (float) nd;
  
  Du = ((float)(nd-1)) * du;
  
  // find the largest u_contour < umn

  uc = umn + du * (float) roundToInt((u[imn] - umn) / du);  if (uc > u[imn]) uc -= du;

  if (uc > u[imn]) cout << "  uc > u[imn]  !!\n\n";
  
  // check whether only one colour is required

  if (u[imx] <= uc + du*1.0001)
  {
    if (dim == 2) { xy2D(x1,ix); xy2D(x2,&ix[2]); xy2D(x3,&ix[4]); }
    else          { xy3D(x1,ix); xy3D(x2,&ix[2]); xy3D(x3,&ix[4]); }

    col = roundToInt((uc - umn) / Du * (float)mxCol);

    essGrpSetColour(col,false);

    essGrpFillPoly(ix,3);

    if (!psOpen) return;

    if (dim == 2) { xy2DPS(x1,ix); xy2DPS(x2,&ix[2]); xy2DPS(x3,&ix[4]); }
    else          { xy3DPS(x1,ix); xy3DPS(x2,&ix[2]); xy3DPS(x3,&ix[4]); }

    files.Pout << ix[4]-ix[2] << " " << ix[5]-ix[3] << " " 
	       << ix[2]-ix[0] << " " << ix[3]-ix[1] << " " 
	       << ix[0] << " " << ix[1] << " " << col << " O\n";
    return;
  }  

  // otherwise do the complicated job ....
  
  //         * imn
  //         |\
  //         |A\
  //         | /\
  //         |/  \    \
  //         | B /\   beta 
  //    |    |  /  \    \
  //  alpha  | /   /\    \
  //    |    |/ C /  \
  //    |    |   / D /\
  //         |  /   /E \
  //         *----------* imx
  //        imm
  
  int i, nb = 1, jmn = imn*dim, jmx = imx*dim, jmm = imm*dim, dim2 = dim + dim, dim3 = dim2 + dim;
 
  float s[150], alph, beta, dalp, dbet, r1 = 1.0, xx[10], 
	du1 =     u[imx] - u[imn],
	du2 = abs(u[imx] - u[imm]),
	du3 = abs(u[imm] - u[imn]);
 
  for (i=0; i<dim; i++) { xx[0+i] = x1[i]; xx[dim+i] = x2[i]; xx[dim2+i] = x3[i]; }

  for (i=0; i<dim; i++) s[i] = xx[jmn+i];

  col = roundToInt((uc - umn) / Du * (float)mxCol);

  uc += du;
  
  beta = (uc - u[imn]) / du1;
  dbet = du / du1;
	  
  // deal with triangle A and with all quadrilaterals B

  if (du3 < 1.e-5 * du) alph = 2.0;
  else
  {
    alph = (uc - u[imn]) / du3;
    dalp = du / du3;
  }
 
  while (alph <= r1) 
  {
    for (i=0; i<dim; i++)
    {
      s[nb*dim+i]     = (r1-beta) * xx[jmn+i] + beta * xx[jmx+i];
      s[(nb+1)*dim+i] = (r1-alph) * xx[jmn+i] + alph * xx[jmm+i];
    }

    if (dim == 2) { for (i=0; i<nb+2; i++) xy2D(s+i+i  ,ix+i+i); }
    else          { for (i=0; i<nb+2; i++) xy3D(s+i+i+i,ix+i+i); }
    
    essGrpSetColour(col,false);

    essGrpFillPoly(ix,nb+2);

    if (psOpen)
    {
      if (dim == 2) { for (i=0; i<nb+2; i++) xy2DPS(s+i+i  ,ix+i+i); }
      else          { for (i=0; i<nb+2; i++) xy3DPS(s+i+i+i,ix+i+i); }

      if (nb == 2) files.Pout << ix[6]-ix[4] << " " << ix[7]-ix[5] << " ";
      
      files.Pout << ix[4]-ix[2] << " " << ix[5]-ix[3] << " " 
	         << ix[2]-ix[0] << " " << ix[3]-ix[1] << " " 
	         << ix[0] << " " << ix[1] << " " << col;
     
      if (nb == 2) files.Pout << " P\n"; else files.Pout << " O\n";
    }
    
    alph += dalp;
    beta += dbet;

    for (i=0; i<dim; i++) s[dim+i] = s[nb*dim+i];
    for (i=0; i<dim; i++) s[i]     = s[(nb+1)*dim+i];

    nb  = 2; 
    
    col = roundToInt((uc - umn) / Du * (float)mxCol);
    
    uc += du;
  }  

  // deal with pentagonal C
 
  if (beta <= r1 && du2>1.e-4*du1)
  {
    alph = (uc - u[imm]) / du2;
    dalp = du / du2;
    
    for (i=0; i<dim; i++)
    {
      s[nb*dim+i]     = (r1-beta) * xx[jmn+i] + beta * xx[jmx+i];
      s[(nb+1)*dim+i] = (r1-alph) * xx[jmm+i] + alph * xx[jmx+i];
      s[(nb+2)*dim+i] = xx[jmm+i];
    }
    
    if (dim == 2) { for (i=0; i<nb+3; i++) xy2D(s+i+i,  ix+i+i); }
    else          { for (i=0; i<nb+3; i++) xy3D(s+i+i+i,ix+i+i); }
    
    essGrpSetColour(col,false);

    essGrpFillPoly(ix,nb+3);

    if (psOpen)
    {
      if (dim == 2) { for (i=0; i<nb+3; i++) xy2DPS(s+i+i,  ix+i+i); }
      else          { for (i=0; i<nb+3; i++) xy3DPS(s+i+i+i,ix+i+i); }

      if (nb == 2) files.Pout << ix[8]-ix[6] << " " << ix[9]-ix[7] << " ";
      
      files.Pout << ix[6]-ix[4] << " " << ix[7]-ix[5] << " " 
	         << ix[4]-ix[2] << " " << ix[5]-ix[3] << " " 
	         << ix[2]-ix[0] << " " << ix[3]-ix[1] << " " 
	         << ix[0] << " " << ix[1] << " " << col;
     
      if (nb == 2) files.Pout << " Q\n"; else files.Pout << " P\n";
    }
    
    alph += dalp;
    beta += dbet;

    for (i=0; i<dim; i++)
    {
      s[dim+i] = s[nb*dim+i];
      s[i]     = s[(nb+1)*dim+i];
    }
    
    nb  = 2; 
    
    col = roundToInt((uc - umn) / Du * (float)mxCol);
    
    uc += du;
  }  
  else if (u[imx]-uc+du>1.e-4*du1) 
  {
    for (i=0; i<dim; i++)
    {
      s[nb*dim+i]     = xx[jmx+i];
      s[(nb+1)*dim+i] = xx[jmm+i];
    }

    if (dim == 2) { for (i=0; i<nb+2; i++) xy2D(s+i+i,  ix+i+i); }
    else          { for (i=0; i<nb+2; i++) xy3D(s+i+i+i,ix+i+i); }
    
    essGrpSetColour(col,false);

    essGrpFillPoly(ix,nb+2);

    if (psOpen)
    {
      if (dim == 2) { for (i=0; i<nb+2; i++) xy2DPS(s+i+i,  ix+i+i); }
      else          { for (i=0; i<nb+2; i++) xy3DPS(s+i+i+i,ix+i+i); }

      if (nb == 2) files.Pout << ix[6]-ix[4] << " " << ix[7]-ix[5] << " ";
      
      files.Pout << ix[4]-ix[2] << " " << ix[5]-ix[3] << " " 
	         << ix[2]-ix[0] << " " << ix[3]-ix[1] << " " 
	         << ix[0] << " " << ix[1] << " " << col;
     
      if (nb == 2) files.Pout << " P\n"; else files.Pout << " O\n";
    }
  } 
  
  // deal with all quadrilaterals D

  while (alph <= r1)
  {
    for (i=0; i<dim; i++)
    {
      s[dim2+i] = (r1-beta) * xx[jmn+i] + beta * xx[jmx+i];
      s[dim3+i] = (r1-alph) * xx[jmm+i] + alph * xx[jmx+i];
    }
    
    if (dim == 2) { for (i=0; i<4; i++) xy2D(s+i+i,  ix+i+i); }
    else          { for (i=0; i<4; i++) xy3D(s+i+i+i,ix+i+i); }
    
    essGrpSetColour(col,false);

    essGrpFillPoly(ix,4);

    if (psOpen)
    {
      if (dim == 2) { for (i=0; i<4; i++) xy2DPS(s+i+i,  ix+i+i); }
      else          { for (i=0; i<4; i++) xy3DPS(s+i+i+i,ix+i+i); }

      files.Pout << ix[6]-ix[4] << " " << ix[7]-ix[5] << " " 
	         << ix[4]-ix[2] << " " << ix[5]-ix[3] << " " 
	         << ix[2]-ix[0] << " " << ix[3]-ix[1] << " " 
	         << ix[0] << " " << ix[1] << " " << col << " P\n";
    }
    
    alph += dalp;
    beta += dbet;

    for (i=0; i<dim; i++) { s[dim+i] = s[dim2+i]; s[i] = s[dim3+i]; }

    col = roundToInt((uc - umn) / Du * (float)mxCol);
    
    uc += du;
  }
  
  // deal with triangle E

  if (u[imx]-uc+du > 1.e-4*du1)
  {
    for (i=0; i<dim; i++) s[dim2+i] = xx[jmx+i];

    if (dim == 2) { for (i=0; i<3; i++) xy2D(s+i+i,  ix+i+i); }
    else          { for (i=0; i<3; i++) xy3D(s+i+i+i,ix+i+i); }
    
    essGrpSetColour(col,false);

    essGrpFillPoly(ix,3);

    if (psOpen)
    {
      if (dim == 2) { for (i=0; i<3; i++) xy2DPS(s+i+i,  ix+i+i); }
      else          { for (i=0; i<3; i++) xy3DPS(s+i+i+i,ix+i+i); }

      files.Pout << ix[4]-ix[2] << " " << ix[5]-ix[3] << " " 
	         << ix[2]-ix[0] << " " << ix[3]-ix[1] << " " 
	         << ix[0] << " " << ix[1] << " " << col << " O\n";
    }
  }

  return;
}









void Plot::contourPlotLegend(double umin, double umax, double umn, double umx, bool flg, int nCol)
{
  if (noGUI) { cout << "    So far, the legend is only available on the screen!\n\n"; return; }

  int iw = (int) w, i, c, mxCol = 511, y1, y2, 
      
                                // all lengths in pixel

      xtr[2]   = { iw-8, 8 },   // top right corner of legend area
      dla[2]   = { 140, 250 },  // width and height of legend area
      
      tlo[2]   = { 14, 10 },    // top/bottom left offset of colour box
      
      wcb      = 40,            // width of colour box

      htxt     = 18,            // height of text
      
      hcb      = dla[1] - 2 * tlo[1];
     
  if (flg) { tlo[1] += htxt; hcb -= 2 * htxt; }
  
  int xcrn[10] = { xtr[0],        xtr[1], 
	           xtr[0]-dla[0], xtr[1], 
		   xtr[0]-dla[0], xtr[1]+dla[1], 
		   xtr[0],        xtr[1]+dla[1] },
  
      xtl[2]   = { xtr[0]-dla[0]+tlo[0], xtr[1]+tlo[1] };
  
  float fctY = (float) hcb / (float) nCol,
         fctC = (float) mxCol / (float) (nCol-1);
  
  // clean legend area
	  
  setColour(8);
  
  essGrpFillPoly(xcrn,4);  

  // plot colour box
  
  for (i=0; i<nCol; i++)
  {
    y1 = xtl[1] + roundToInt((float) i * fctY);
    y2 = xtl[1] + roundToInt((float) (i+1) * fctY);  
    
    xcrn[0] = xtl[0];
    xcrn[1] = y1;
    xcrn[2] = xtl[0]+wcb;
    xcrn[3] = y1;
    xcrn[4] = xtl[0]+wcb;
    xcrn[5] = y2;
    xcrn[6] = xtl[0];
    xcrn[7] = y2;
   
    c = 511 - roundToInt((float) i * fctC); //cout << c << "\n";
    
    setColour(c,false);
    
    essGrpFillPoly(xcrn,4);  
  }

  // draw a few lines to make it look nice

  setColour(7);
 
  y1 = xtl[1] + roundToInt((float) hcb * 0.5);
  
  essGrpDrawLine(xtl[0],     xtl[1],     xtl[0]+wcb+10, xtl[1]    ); 
  essGrpDrawLine(xtl[0],     y1,         xtl[0]+wcb+10, y1        );
  essGrpDrawLine(xtl[0],     xtl[1]+hcb, xtl[0]+wcb+10, xtl[1]+hcb);
  
  essGrpDrawLine(xtl[0],     xtl[1],     xtl[0],        xtl[1]+hcb); 
  essGrpDrawLine(xtl[0]+wcb, xtl[1],     xtl[0]+wcb,    xtl[1]+hcb); 
 
  // put labels

  char strg[30];

  sprintf(strg,"%-10.4g",umax);  
  xcrn[0] = xtl[0]+wcb+15;  
  xcrn[1] = xtl[1];  
  essGrpPutText(xcrn,strg,4);
  
  sprintf(strg,"%-10.4g",0.5*(umin+umax));  
  xcrn[1] = y1;  
  essGrpPutText(xcrn,strg,4);
  
  sprintf(strg,"%-10.4g",umin);  
  xcrn[1] = xtl[1]+hcb;  
  essGrpPutText(xcrn,strg,4);
 
  if (!flg) return;

  xcrn[0] = xtl[0] - roundToInt((float)tlo[0] * 0.5);
  y1      = 3      + roundToInt((float)htxt * 0.5);
  
  xcrn[1] = xtl[1] - y1;
  sprintf(strg,"max = %-10.4g",umx);
  essGrpPutText(xcrn,strg,1);
	  
  xcrn[1] = xtl[1]+hcb+y1;
  sprintf(strg,"min = %-10.4g",umn);
  essGrpPutText(xcrn,strg,7);
	  
  return;
}









void Plot::adjust1DSettings(bool xcur, double xmn, double xmx, 
		            bool ucur, double umn, double umx, bool vert)
{
  if (hPix < 0)
  {
    essGrpGetPlotAreaGeom();
    if (aspectRatioCorr) h = h / YPixPerXPix;
  }
  
  if (!xcur) { x1D  = xmn; dx1D  = xmx - xmn; }

  if (!ucur) { u1D  = umn; du1D  = umx - umn; }
  
  vert1D  = vert;

  if (!vert1D)
  {
    wFact1D = w / dx1D;
    hFact1D = h / du1D;
  }
  else
  {
    wFact1D = w / du1D;
    hFact1D = h / dx1D;
  }
  
  return;
}










bool Plot::psOpenFile(void)
{ 
  MyString pathAndFile;
  
  char tmp[20];
 
  sprintf(tmp,".%04d.eps",++files.nps);
  
  pathAndFile.append(files.projDir).append(SLASH).append(files.Pfile).append(tmp);

  while (prgFileExist(pathAndFile)) pathAndFile.trunc(pathAndFile.length()-4).append(".new.eps");
  
  COUT << &((pathAndFile.asCharArray())[files.projDir.length()+1]) << " open for plot output.\n\n";
 
  pathAndFile.append(".tmp");

  files.Pout.open(pathAndFile.asCharArray());

  if (files.Pout.is_open()) return true;
  
  prgWarning(1,"Plot::psOpenFile","could not open file!");
 
  return false; 
}







bool Plot::psCloseFile(void)
{ 
  ofstream Pout2; 
  ifstream Pout3;
  
  MyString fileName, pathAndFile;
  
  char c, tmp[200];
 
  // close temporary file

  files.Pout.close();

  // open final eps file

  sprintf(tmp,".%04d.eps",files.nps);
  
  pathAndFile.append(files.projDir).append(SLASH).append(files.Pfile).append(tmp);

  while (prgFileExist(pathAndFile)) pathAndFile.trunc(pathAndFile.length()-4).append(".new.eps");
  
  fileName.append(&((pathAndFile.asCharArray())[files.projDir.length()+1]));
	  
  Pout2.open(pathAndFile.asCharArray());

  // calculate bounding box
  
  int psXmx = roundToInt(dDes[0] * wFactDesPS / psFact),
      psYmx = roundToInt(dDes[1] * hFactDesPS / psFact);
      
//  cout << '\t' << psXmx << '\t' << psYmx << endl;
//  psXmx = 500;  psYmx = 600;
  
  sprintf(tmp,"0  0  %d  %d",psXmx,psYmx);

  // write header stuff
  
  Pout2 << "%!PS-Adobe-3.0 EPSF-3.0\n";
  Pout2 << "%%BoundingBox:    " << tmp << "\n";
  Pout2 << "%%Title:       MPAP2 (" << fileName << ")\n";
  Pout2 << "%%\n";
  Pout2 << "%%  generated with MPAP2\n";
  Pout2 << "%%\n";
  Pout2 << "%%  at mpapTime.cur = " << mpapTime.cur << "\n";
  Pout2 << "%%\n";
  Pout2 << "%%  Colour flag: 0 -> black/white\n";
  Pout2 << "%%               1 -> gray scale\n";
  Pout2 << "%%               2 -> colour\n";
  Pout2 << "%%\n";
  Pout2 << "%%=======================================\n";
  Pout2 << "/Colour            2  def\n";
  Pout2 << "/TextScale       100  def\n";
  Pout2 << "/BaseLineWidth   0.1  def\n";
  Pout2 << "/PointSize        30  def\n";
  Pout2 << "/BlackWhiteSwap  true def\n";
  Pout2 << "%%=======================================\n";
  Pout2 << "/m  {moveto}         bind def\n";
  Pout2 << "/p  {rmoveto}        bind def\n";
  Pout2 << "/v  {rlineto}        bind def\n";
  Pout2 << "/l  {lineto}         bind def\n";
  Pout2 << "/s  {stroke}         bind def\n";
  Pout2 << "/lw {setlinewidth}   bind def\n";
  Pout2 << "/sc {setrgbcolor}    bind def\n";
  Pout2 << "/g  {setgray}        bind def\n";
  Pout2 << "/n  {newpath}        bind def\n";
  Pout2 << "/cp {closepath}      bind def\n";
  Pout2 << "/f  {fill}           bind def\n";


  Pout2 << "/fct { 2 Colour eq { 0.0078125 } { 0.001953125 } ifelse } def\n";
  Pout2 << "/c { 2 Colour eq { dup 128 lt { fct mul 0 2 1 roll 1 sc} {dup 256 lt\n";
  Pout2 << "  {-128 add fct mul neg 1 add 0 1 3 2 roll sc} {dup 384 lt\n";
  Pout2 << "  {-256 add fct mul 1 0 sc}{-384 add fct mul neg 1 add 1 exch 0 sc}\n";
  Pout2 << "  ifelse}ifelse}ifelse} { fct mul g } ifelse } def\n";
  Pout2 << "/BL { BlackWhiteSwap { 0 } { 1 } ifelse} def\n";
  Pout2 << "/WH { BlackWhiteSwap { 1 } { 0 } ifelse} def\n";
  Pout2 << "/C1 { 2 Colour eq {s 1 0 0 sc}  { 1 Colour eq {s .578 g}{s BL g} ifelse} ifelse} def %% red\n";
  Pout2 << "/C2 { 2 Colour eq {s 0 0 1 sc}  { 1 Colour eq {s .578 g}{s BL g} ifelse} ifelse} def %% blue\n";
  Pout2 << "/C3 { 2 Colour eq {s 0 1 0 sc}  { 1 Colour eq {s .578 g}{s BL g} ifelse} ifelse} def %% green\n";
  Pout2 << "/C4 { 2 Colour eq {s 1 1 0 sc}  { 1 Colour eq {s .816 g}{s BL g} ifelse} ifelse} def %% yellow\n";
  Pout2 << "/C5 { 2 Colour eq {s 0 1 1 sc}  { 1 Colour eq {s .816 g}{s BL g} ifelse} ifelse} def %% cyan\n";
  Pout2 << "/C6 { 2 Colour eq {s 1 0 1 sc}  { 1 Colour eq {s .816 g}{s BL g} ifelse} ifelse} def %% magenta\n";
  Pout2 << "/C7 { 2 Colour eq {s .5 .5 1 sc}{ 1 Colour eq {s .707 g}{s BL g} ifelse} ifelse} def %% lightblue\n";
  Pout2 << "/C8 { s BL g} def %% black\n";
  Pout2 << "/C9 { s WH g} def %% white\n";
  Pout2 << "/bg { BlackWhiteSwap { } { 0 g 0 0 m " << psXmx*psFact << " 0 v 0 "
                                                  << psYmx*psFact << " v -"
                                                  << psXmx*psFact << " 0 v f } ifelse } def\n";
  Pout2 << "/V {v s} def\n";
  Pout2 << "/O {c m v v f} def\n";
  Pout2 << "/P {c m v v v f} def\n";
  Pout2 << "/Q {c m v v v v f} def\n";
  Pout2 << "/R {c m v v v v v f} def\n";
  Pout2 << "/S {c m v v v v v v f} def\n";
  Pout2 << "/T {c m v v v v v v v f} def\n";
  Pout2 << "/CH { 0 360 s arc s} def\n";
  Pout2 << "/CF { 0 360 s arc f} def\n";
  Pout2 << "/H {/Helvetica findfont TextScale scalefont setfont} def\n";
  Pout2 << "/w {stringwidth pop exch div neg 0 p} def\n";
  Pout2 << "/x {0 exch TextScale 0.7 mul exch div neg p} def\n";
  Pout2 << "/blw BaseLineWidth def\n";
  Pout2 << "/PS { PointSize } def\n";
  Pout2 << "/P0 { PS CF} def\n";
  Pout2 << "/PN { 3 1 roll 2 copy PS CF PS .8 mul add exch\n";
  Pout2 << "  PS .8 mul add exch moveto show} def\n";
  Pout2 << "%%==== set scale and line style =========\n";
  Pout2 << "0.1  0.1  scale\n";
  Pout2 << "1 setlinecap\n";
  Pout2 << "1 setlinejoin\n";
  Pout2 << "[] 0 setdash\n";
  Pout2 << "blw lw H bg\n";
  Pout2 << "%%==== plot-macros ======================\n";

  // copy content of temporary file, close and delete tmp file
 
  Pout3.open(pathAndFile.append(".tmp").asCharArray());

  while (Pout3)
  {
    Pout3.get(c);
//    cout << c ;
    Pout2 << c;
  } 
  cout << endl;

  Pout3.close();

  pathAndFile.insert(0," ").insert(0,DELETE_FILE_COMMAND);
  if (system(pathAndFile.asCharArray()) != 0)
    std::cout << "   WARNING! I could not remove the temporary eps file!\n\n"; 

  // write conclusion

  Pout2 << "%%==== close file =======================\n";
  Pout2 << "s\n";
  Pout2 << "showpage\n";
  Pout2 << "%%EOF\n";

  // close eps file
  
  Pout2.close();
  
  COUT << fileName << " closed.\n\n";
	
  return true; 
} 








void Plot::drawDesFrame(void)
{
  int ix1[2], ix2[2];

  ix1[0] =        roundToInt((x0Des[0]-x0Act[0]) * wFactAct);
  ix1[1] = hPix - roundToInt((x0Des[1]-x0Act[1]) * hFactAct);
  
  ix2[0] =        roundToInt((x0Des[0]+dDes[0]-x0Act[0]) * wFactAct);
  ix2[1] = hPix - roundToInt((x0Des[1]+dDes[1]-x0Act[1]) * hFactAct);

  setColour(4);

  essGrpDrawLine(ix1[0],ix1[1],ix2[0],ix1[1]);
  essGrpDrawLine(ix2[0],ix1[1],ix2[0],ix2[1]);
  essGrpDrawLine(ix2[0],ix2[1],ix1[0],ix2[1]);
  essGrpDrawLine(ix1[0],ix2[1],ix1[0],ix1[1]);

  return;
}








void Plot::wipe(void)
{
  essGrpWipe();

  return;
}








void Plot::point(float *x, float dpt, int pn)
{
  int i, ix[2], id[2];
  char txt[30];
  float d = dpt, picCoor[2]; if (d < 0.) d = dPt();

  id[0] = roundToInt(d * wFactAct);
  id[1] = roundToInt(d * hFactAct);

  if (dim == 1)      {    }
  else if (dim == 2) { xy2D(x,ix); }
  else               { perspective.pictureCoor(x,picCoor); xy2D(picCoor,ix); }

  //cout << ix[0] << "," << ix[1] << "\n";
  //cout << id[0] << "," << id[1] << "\n";

  if (clipPoint2D(ix)) essGrpFillCircle(ix,id);
	   
  ix[0] += roundToInt(0.6 * (float)id[0]);
  ix[1] -= roundToInt(0.6 * (float)id[1]);

  if (pn > 0) { sprintf(txt,"%d",pn); essGrpPutText(ix,txt,1); }

  if (!psOpen) return;

  if (dim == 1)      {    }
  else if (dim == 2) { xy2DPS(x,ix); }
  else               {    }
  
  files.Pout << ix[0] << " " << ix[1];

  if (pn > 0) files.Pout << " (" << txt << ") PN\n"; else files.Pout << " P0\n";

  return;
}





void Plot::dash(float *x, float d, int ii)
{
  int ix[4], id[2];

  float x2[3] = { x[0], x[1], x[2] }, xp1[2], xp2[2], vp[2], fact;

  if (dim == 1)      
  {    

  }
  else if (dim == 2) 
  {
    xp1[0] = x[0];
    xp1[1] = x[1];
    xp2[0] = x[0];
    xp2[1] = x[1];

    xp1[ii] -= .5 * d;
    xp2[ii] += .5 * d;
  }
  else               
  {
    x2[ii] += 1.e4;

    perspective.pictureCoor(x, xp1);
    perspective.pictureCoor(x2,xp2);

    vp[0] = xp2[0] - xp1[0];
    vp[1] = xp2[1] - xp1[1];

    fact = .5*d/sqrt(vp[0]*vp[0]+vp[1]*vp[1]);

    vp[0] *= fact;
    vp[1] *= fact;

    xp2[0] = xp1[0] + vp[0];
    xp2[1] = xp1[1] + vp[1];

    xp1[0] -= vp[0];
    xp1[1] -= vp[1];
  }

  xy2D(xp1,ix);
  xy2D(xp2,ix+2);

  essGrpDrawLine(ix[0],ix[1],ix[2],ix[3]);
	   
  if (!psOpen) return;

  return;
}








void Plot::line(float *x1, float *x2)
{
  int ix1[2], ix2[2];

  if (dim == 1)      { xy1D(x1,ix1); xy1D(x2,ix2); }
  else if (dim == 2) { xy2D(x1,ix1); xy2D(x2,ix2); }
  else               { xy3D(x1,ix1); xy3D(x2,ix2); }

  if (clipLine2D(ix1,ix2)) essGrpDrawLine(ix1[0],ix1[1],ix2[0],ix2[1]);
  
  if (!psOpen) return;

  if (dim == 1)      { xy1DPS(x1,ix1); xy1DPS(x2,ix2); }
  else if (dim == 2) { xy2DPS(x1,ix1); xy2DPS(x2,ix2); }
  else               { xy3DPS(x1,ix1); xy3DPS(x2,ix2); }
  
  if (clipLine2DPS(ix1,ix2))

    files.Pout << ix1[0]        << " " << ix1[1]        << " m " 
               << ix2[0]-ix1[0] << " " << ix2[1]-ix1[1] << " V\n";

  return;
}









void Plot::circle(float *x, float &d)
{
  int ix[2], id[2];

  if (dim == 1)      {    }
  else if (dim == 2) xy2D(x,ix);
  else               xy3D(x,ix);

  if (clipPoint2D(ix))
  {
    id[0] = roundToInt(d * wFactAct);
    id[1] = roundToInt(d * hFactAct);

    essGrpDrawCircle(ix,id);
  }	   
  if (!psOpen) return;
/*
  if (dim == 1)      {    }
  else if (dim == 2) { xy2DPS(x,ix); 
	               id[0] = roundToInt(d * wFactDesPS *.5); }
  else               {    }
  
  files.Pout << ix[0] << " " << ix[1] << " " << id[0] << " CH\n";
*/
return;
}






void Plot::fillCircle(float *x, float &d)
{
  int ix[2], id[2];

  if (dim == 1)      {    }
  else if (dim == 2) xy2D(x,ix);
  else               xy3D(x,ix);

  id[0] = roundToInt(d * wFactAct);
  id[1] = roundToInt(d * hFactAct);

  essGrpFillCircle(ix,id);
	   
  if (!psOpen) return;
/*
  if (dim == 1)      {    }
  else if (dim == 2) { xy2DPS(x,ix); 
	               id[0] = roundToInt(d * wFactDesPS *.5); }
  else               {    }
  
  files.Pout << ix[0] << " " << ix[1] << " " << id[0] << " CF\n";
*/
  return;
}






void Plot::triangle(float *x, float &d)
{
  int i, ix[8], id2d3[2], id1d3[2], id1d2[2];
  
  float xx[3],
        d2d3 = -d * 0.6666666,
        d1d3 = -d * 0.3333333,
        d1d2 =  d * 0.5;

  if      (dim == 1) {   }
  else if (dim == 2) xy2D(x,ix);
  else               xy3D(x,ix);

  if (clipPoint2D(ix))
  {
    id2d3[0] = roundToInt(d2d3 * wFactAct);
    id2d3[1] = roundToInt(d2d3 * hFactAct);

    id1d3[0] = roundToInt(d1d3 * wFactAct);
    id1d3[1] = roundToInt(d1d3 * hFactAct);

    id1d2[0] = roundToInt(d1d2 * wFactAct);
    id1d2[1] = roundToInt(d1d2 * hFactAct);

    ix[2] = ix[0]+id1d2[0];
    ix[3] = ix[1]-id1d3[1];
    ix[4] = ix[0];
    ix[5] = ix[1]+id2d3[1];
    ix[6] = ix[0]-id1d2[0];
    ix[7] = ix[3];
    ix[0] = ix[6];
    ix[1] = ix[7];
   
    for (i=0; i<5; i+=2) essGrpDrawLine(ix[i],ix[i+1],ix[i+2],ix[i+3]);
  }
  if (!psOpen) return;
/*
  if (dim == 1)      {   }
  else if (dim == 2) { xx[0] = x[0];      xx[1] = x[1]+d2d3; xy2DPS(xx, ix);
	               xx[0] = x[0]+d1d2; xx[1] = x[1]-d1d3; xy2DPS(xx,&ix[2]);
	               xx[0] = x[0]-d1d2;                    xy2DPS(xx,&ix[4]); }
  else if (dim == 3) {   }
  
  ix[6] = ix[0];
  ix[7] = ix[1];
	    
  files.Pout << ix[0] << " " << ix[1] << " m";
  for (i=2; i<6; i+=2) files.Pout << " " << ix[i]-ix[i-2] << " " << ix[i+1]-ix[i-1] << " v";
  files.Pout << "\n";
*/
  return;
}






void Plot::square(float *x, float &d)
{
  int i, ix[10], id[2];
  
  float xx[3], d1d2 = d * 0.5;
  
  if      (dim == 1) {   }
  else if (dim == 2) xy2D(x,ix);
  else               xy3D(x,ix);

  if (clipPoint2D(ix))
  {
    id[0] = roundToInt(d1d2 * wFactAct);
    id[1] = roundToInt(d1d2 * hFactAct);

    ix[2] = ix[0]+id[0];
    ix[3] = ix[1]-id[1];
    ix[4] = ix[2];
    ix[5] = ix[1]+id[1];
    ix[6] = ix[0]-id[0];
    ix[7] = ix[5];
    ix[8] = ix[6];
    ix[9] = ix[3];
    ix[0] = ix[8];
    ix[1] = ix[9];
	    
    for (i=0; i<7; i+=2) essGrpDrawLine(ix[i],ix[i+1],ix[i+2],ix[i+3]);
  }
  if (!psOpen) return;
/*
  if (dim == 1)      {    }
  else if (dim == 2) { xx[0] = x[0]-d1d2; xx[1] = x[1]-d1d2; xy2DPS(xx,&ix[0]);
	               xx[0] = x[0]+d1d2;                    xy2DPS(xx,&ix[2]);
		                          xx[1] = x[1]+d1d2; xy2DPS(xx,&ix[4]); }
  else               {    }

  ix[6] = ix[0];
  ix[7] = ix[5];
  ix[8] = ix[0];
  ix[9] = ix[1];

  files.Pout << ix[0] << " " << ix[1] << " m";
  for (i=2; i<8; i+=2) files.Pout << " " << ix[i]-ix[i-2] << " " << ix[i+1]-ix[i-1] << " v";
  files.Pout << "\n";
*/
  return;
}






void Plot::putText(float *x, char *txt, int posFlg, bool BGFlag)
{
  int i, ix[2], id[2], pF = posFlg;
  
  if (dim == 1)      {             }
  else if (dim == 2) { xy2D(x,ix); }
  else               {             }

  essGrpPutText(ix,txt,posFlg,BGFlag);

  if (!psOpen) return;

  if (dim == 1)      {               }
  else if (dim == 2) { xy2DPS(x,ix); }
  else               {               }

  files.Pout << ix[0] << " " << ix[1] << " m ";
 
  if (posFlg < 1) pF = 1; if (posFlg > 9) pF = 9;
 
  i = pF - intDiv(pF-1,3) * 3; if (i == 1) i = 0; else if (i == 3) i = 1;

  if (i > 0) files.Pout << i << " (" << txt << ") w ";

  i = intDiv(pF+2,3); if (i == 1) i = 0; else if (i == 3) i = 1;

  if (i > 0) files.Pout << i << " x ";
  
  files.Pout << "(" << txt << ") show\n";

  return;
}





void Plot::arrow(float *x0, float *dl, float &sclFct)
{
  if (noGUI) return;

  float d1[2], d2[2];
 
  switch (dim)
  {
    case 1: break;
	    
    case 2: d1[0] = - dl[0] * sclFct;
	    d1[1] = - dl[1] * sclFct;
	   
	    d2[0] = - d1[1];
            d2[1] = + d1[0];
	           
	    polyPnt [0] = (float) x0[0];
            polyPnt [1] = (float) x0[1];
	   
            polyPnt [3] = (float)(x0[0] + 0.2*d1[0] - 0.070*d2[0]);
            polyPnt [4] = (float)(x0[1] + 0.2*d1[1] - 0.070*d2[1]);
           
            polyPnt [6] = (float)(x0[0] + 0.2*d1[0] - 0.011*d2[0]);
            polyPnt [7] = (float)(x0[1] + 0.2*d1[1] - 0.011*d2[1]);
           
            polyPnt [9] = (float)(x0[0] +     d1[0] - 0.011*d2[0]);
            polyPnt[10] = (float)(x0[1] +     d1[1] - 0.011*d2[1]);
           
            polyPnt[12] = (float)(x0[0] +     d1[0] + 0.011*d2[0]);
            polyPnt[13] = (float)(x0[1] +     d1[1] + 0.011*d2[1]);
           
            polyPnt[15] = (float)(x0[0] + 0.2*d1[0] + 0.011*d2[0]);
            polyPnt[16] = (float)(x0[1] + 0.2*d1[1] + 0.011*d2[1]);
           
            polyPnt[18] = (float)(x0[0] + 0.2*d1[0] + 0.070*d2[0]);
            polyPnt[19] = (float)(x0[1] + 0.2*d1[1] + 0.070*d2[1]);

	    break;
	    
    case 3: break;
  }

  polyPnt[21] = polyPnt[0];
  polyPnt[22] = polyPnt[1];
  polyPnt[23] = polyPnt[2];

  polyPnt[25] = (polyPnt [9] + polyPnt[12]) * 0.5;
  polyPnt[26] = (polyPnt[10] + polyPnt[13]) * 0.5;
  polyPnt[27] = (polyPnt[11] + polyPnt[14]) * 0.5;
  
  poly(8);

  line(&(polyPnt[0]),&(polyPnt[25]));
  
  return;
}







bool Plot::clipLine2D(int *ix1, int *ix2)
{
  unsigned int clipBit1 = 0, clipBit2 = 0;

  if (ix1[0] < 0) clipBit1 = 1;            else if (ix1[0] > wPix) clipBit1 = 2;
  if (ix1[1] < 0) clipBit1 = clipBit1 | 4; else if (ix1[1] > hPix) clipBit1 = clipBit1 | 8;

  if (ix2[0] < 0) clipBit2 = 1;            else if (ix2[0] > wPix) clipBit2 = 2;
  if (ix2[1] < 0) clipBit2 = clipBit2 | 4; else if (ix2[1] > hPix) clipBit2 = clipBit2 | 8;

  if (clipBit1 == 0 && clipBit2 == 0) return true;

  //cout << wPix << "," << hPix << "\n";
  //cout << ix1[0] << "," << ix1[1] << " -> " << clipBit1 << "\n";
  //cout << ix2[0] << "," << ix2[1] << " -> " << clipBit2 << "\n\n";

  if (clipBit1 & 1 && clipBit2 & 1) return false;
  if (clipBit1 & 2 && clipBit2 & 2) return false;
  if (clipBit1 & 4 && clipBit2 & 4) return false;
  if (clipBit1 & 8 && clipBit2 & 8) return false;

  return true;
}







bool Plot::clipLine2DPS(int *ix1, int *ix2)
{


  return true;
}







bool Plot::clipPoint2D(int *ix)
{
  if (ix[0] < 0)    return false;
  if (ix[0] > wPix) return false;
  if (ix[1] < 0)    return false;
  if (ix[1] > hPix) return false;

  return true;
}  

