
#include "GeomSurface.h"
#include "Plot.h"
#include "FunctionsProgram.h"
#include "AdvancingFront2D.h"
#include "MyString.h"
#include "FunctionsEssGrp.h"


extern Plot plot;



GeomSurface::GeomSurface(void)
{
  hasBeenDiscretised = false;

  return;
}




GeomSurface::~GeomSurface()
{
  return;
}





bool GeomSurface::initialise(void)
{
  // THE FOLLOWING IS 2D STUFF ONLY !!!

  int i0, i1, i, j, n = suppSpln.n, iLoop = 0;

  Vector<void*> *suppPnt, p0, p1;

  void *p;

  double xmin, xmn, area, xp[2] = {0.,0.};

  // order support spline list, identify loops

  for (i=0; i<suppSpln.n; i++)
  {
    suppPnt = &(((GeomSpline*)(suppSpln[i]))->suppPnt);

    p0.append((*suppPnt)[0]);
    p1.append((*suppPnt)[suppPnt->n-1]);

    orientation.append(FORWARD);
  }

  iLoop = 0;

  loop.append(0);

  i = loop[iLoop];

  while (i < n-1)
  {
    while (i < n-1)
    {
      if (orientation[i] == FORWARD) p = p1[i]; else p = p0[i];

      j = i + 1;
      while (j < n && p != p0[j]) j++;
      if (j < n) 
      {
        suppSpln.move(j,++i);
        p0.move(j,i);
        p1.move(j,i);
      }
      else
      {
        j = i + 1;
        while (j < n && p != p1[j]) j++;
        if (j < n)
        {
          suppSpln.move(j,++i);
          p0.move(j,i);
          p1.move(j,i);
          orientation[i] = BACKWARD;
        }
        else
        {
        if ((p0[loop[iLoop]] != p1[i] && orientation[i] == FORWARD) || 
            (p0[loop[iLoop]] != p0[i] && orientation[i] == BACKWARD))
          prgError(1,"GeomSurface::initialise","surface not closed!");

        loop.append(++i); iLoop++;

        break;
        }
      }
    }
  }
  if ((p0[loop[iLoop]] != p1[i] && orientation[i] == FORWARD) ||  
      (p0[loop[iLoop]] != p0[i] && orientation[i] == BACKWARD))
    prgError(2,"GeomSurface::initialise","surface not closed!");

  loop.append(suppSpln.n);

  // identify outer boundary loop and move to top of list  

  xmin = ((GeomSpline*)(suppSpln[0]))->xmin();
  j = 0;

  for (iLoop=0; iLoop<loop.n-1; iLoop++)
  {
    xmn = ((GeomSpline*)(suppSpln[loop[iLoop]]))->xmin();
    for (i=loop[iLoop]; i<loop[iLoop+1]; i++)
      if (xmn > ((GeomSpline*)(suppSpln[i]))->xmin()) xmn = ((GeomSpline*)(suppSpln[i]))->xmin();

    if (xmin > xmn) { xmin = xmn; j = iLoop; }
  }

  //cout << j << "\n";

  for (i=loop[j]; i<loop[j+1]; i++) { suppSpln.move(i,i-loop[j]); orientation.move(i,i-loop[j]); }

  n = loop[j+1] - loop[j];

  for (i=j; i>0; i--) loop[i] = loop[i-1] + n;

  //cout << loop << "\n\n";

  // correct orientations if necessary 

  for (iLoop=0; iLoop<loop.n-1; iLoop++)
  {
    area = 0.;
    for (i=loop[iLoop]; i<loop[iLoop+1]; i++)
    {
      if (orientation[i] == BACKWARD) area -= ((GeomSpline*)(suppSpln[i]))->area(xp);
      else                            area += ((GeomSpline*)(suppSpln[i]))->area(xp);
    }
    if ((area < 0. && iLoop == 0) || (area > 0. && iLoop > 0))
    {
      i0 = loop[iLoop];
      i1 = loop[iLoop+1];
      for (i=i0; i<i1; i++)
      {
        if (orientation[i] == BACKWARD) orientation[i] = FORWARD; else orientation[i] = BACKWARD;
      }
      for (i=i1-1; i>i0+1; i--) { suppSpln.move(i0+1,i); orientation.move(i0+1,i); }
    }
  }

  return true;
}








void GeomSurface::draw(int surfaceId)
{
/*  int i, iLoop;

  plot.wipe();
  plot.setColour(3);

  for (iLoop=0; iLoop<loop.n-1; iLoop++)
  {
    for (i=loop[iLoop]; i<loop[iLoop+1]; i++)
    {
      ((GeomSpline*)(suppSpln[i]))->draw(0);
      if (orientation[i] == BACKWARD) cout << "BACKWARD\n"; else cout << "FORWARD\n";
      prgUpdateDisplayAndWait();
    }
    cout << "\n";
  }
 */
  return;
}











void GeomSurface::generateMesh(MatrixFullArray<double> &x, 
                               MatrixFullArray<int> &ix, 
                               void *ndGmLnk,
                               int &numel, int &numnp, int &nBndNd,
                               Domain *dom, int nd, bool showFlg)
{
  if (nd < 3 || nd > 9 || nd == 5 || nd == 7)
    prgError(1,"GeomSurface::discretise","invalid nd!");

  MyString ch;

  AdvancingFront2D front;

  front.generateMesh(x,ix,ndGmLnk,nBndNd,*this,dom,nd,showFlg);

  numel = ix.nRow;
  numnp = x.nRow;

  hasBeenDiscretised = true;

  //plot.setColour(3);
  //front.drawPolygons();

  return;
}


