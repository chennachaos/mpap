
#include <iostream>

#include "FunctionsProgram.h"
#include "DomainTypeEnum.h"
#include "Geometry.h"
#include "DomainTree.h"
#include "DomainType.h"
#include "Debug.h"
#include "MyString.h"
#include "DataBlockTemplate.h"
#include "MathBasic.h"
#include "MathGeom.h"



extern DomainTree domain;


static const unsigned int pow2[8] = {1,2,4,8,16,32,64,128};


using namespace std;



Geometry::Geometry(void)                       
{                                                  
  if (debug) std::cout << " Geometry constructor\n\n";
 
  // add new type
  
  DomainType *geom = domain.newType(GEOMETRY,ROOTDOMAIN);

  if (geom == NULL) return;  // domain type exists already

  geom->key.addNew("geometry points", 
                   "geometry splines", 
		   "geometry surfaces",
                   "geometry boundary conditions",
                   "geometry prescribed displacements",
                   "geometry distributed loads",
                   "geometry nodal data output");

  geometryDiscretisedAtLeastOnce = false;

  return;
}




	                                          
Geometry::~Geometry(void)                     
{         
  if (debug) cout << " Geometry destructor\n\n";

  return;
}





void Geometry::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString tmpl, *word;
 
  char tmp[30], fct[] = "Geometry::readInputData",
       *splineGeometry[] = {"spline", "line", "arc", NULL},
       *geomObjectType[] = {"point", "spline", "surface", "volume", NULL},
       *wrndTypeNames[] = WRND_TYPE_NAMES;

  int nw, i, ii, j, tp, n, e, nst, imin, pmin, prp, indx;
 
  double x1[2], x2[2], x3[2], xm[2], fact, alph12, alph23, alph, s, c, dmin, d;

  Vector<double> dTmp;
  Vector<int>    iTmp, lTmp;
  MyStringList   sTmp;

  DataBlockTemplate t1, t2;

  switch (domain[GEOMETRY].key.whichBegins(line))
  {
    case  0: cout << "     GEOMETRY: reading geometry points ...\n\n";
	   
             if (ndf < 1) prgError(1,fct,"'dimensions' must precede 'geometry points'!");
            
	     sprintf(tmp,"123 %df",ndm);  
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	     
             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
		     prgError(2,fct,"data error in 'geometry points'!");
	     
             n = dTmp.dim() / ndm;
	    
             for (i=0; i<n; i++)
             {
               point.add(new GeomPoint);

               point[i].x[0] = dTmp[i+i];
               point[i].x[1] = dTmp[i+i+1];
             }
 
	     break;

    case  1: cout << "     GEOMETRY: reading geometry splines ...\n\n";
	  
             if (point.n < 1) prgError(1,fct,"'geometry points' must precede 'geometry splines'!");
             if (spline.n > 0) prgError(2,fct,"multiple definition of 'geometry splines'!");

             n = 0;

             while (1)
             {
               line.getNextLine(Ifile);
               nw = line.split(&word);

               if (nw < 4)                   break;
               if (!word[0].toInt(&i,false)) break;
               if (i != n + 1)               break;
 
	       tp = word[1].which(splineGeometry); 
               if (tp < 0) prgError(3,fct,"error in 'geometry splines'!");

               iTmp.free();

               i = 2;
               while (i < nw && word[i].toInt(iTmp.append(),false)) i++;
               if (i < nw) break;

               spline.add(new GeomSpline);

               for (i=0; i<nw; i++) word[i].free(); delete [] word;

               switch (tp)
               {
                 case  0: // spline

                          for (i=0; i<nw-2; i++)
                          {
                            if (iTmp[i] < 1 || iTmp[i] > point.n) 
                              prgError(4,fct,"error in 'geometry splines'!");

                            spline[n].suppPnt.append((void*) &(point[iTmp[i]-1]));
                          } 

                          break;

                 case  1: // line

                          if (nw != 4) prgError(5,fct,"error in 'geometry splines'!");

                          for (i=0; i<nw-2; i++)
                          {
                            if (iTmp[i] < 1 || iTmp[i] > point.n) 
                              prgError(4,fct,"error in 'geometry splines'!");

                            spline[n].suppPnt.append((void*) &(point[iTmp[i]-1]));
                          } 

                          break;

                 case  2: // arc

                          if (nw != 5) prgError(6,fct,"error in 'geometry splines'!");

                          if (iTmp[0] < 1 || iTmp[0] > point.n ||
                              iTmp[1] < 1 || iTmp[1] > point.n ||
                              iTmp[2] < 1 || iTmp[2] > point.n)
                            prgError(4,fct,"error in 'geometry splines'!");

                          x1[0] = point[iTmp[0]-1].x[0];
                          x1[1] = point[iTmp[0]-1].x[1];
                          x2[0] = point[iTmp[1]-1].x[0];
                          x2[1] = point[iTmp[1]-1].x[1];
                          x3[0] = point[iTmp[2]-1].x[0];
                          x3[1] = point[iTmp[2]-1].x[1];

                          fact = 0.5 / (x1[0] * (x2[1]-x3[1])
                                      + x2[0] * (x3[1]-x1[1])
                                      + x3[0] * (x1[1]-x2[1]) );

                          xm[0] =   (x1[0] * x1[0] * (x2[1]-x3[1]) 
                                   + x2[0] * x2[0] * (x3[1]-x1[1])
                                   + x3[0] * x3[0] * (x1[1]-x2[1]) 
                                   - (x1[1]-x2[1]) * (x2[1]-x3[1]) * (x3[1]-x1[1]) ) * fact;

                          xm[1] = - (x1[1] * x1[1] * (x2[0]-x3[0]) 
                                   + x2[1] * x2[1] * (x3[0]-x1[0])
                                   + x3[1] * x3[1] * (x1[0]-x2[0]) 
                                   - (x1[0]-x2[0]) * (x2[0]-x3[0]) * (x3[0]-x1[0]) ) * fact;

                          fact = (xm[0]-x1[0])*(xm[0]-x1[0]) + (xm[1]-x1[1])*(xm[1]-x1[1]);

                          dmin = fact * 1000.;

                          fact = 1./ fact;

                          alph12 = myAcos(((x1[0]-xm[0])*(x2[0]-xm[0])
                                         + (x1[1]-xm[1])*(x2[1]-xm[1])) * fact);

                          alph23 = myAcos(((x2[0]-xm[0])*(x3[0]-xm[0])
                                         + (x2[1]-xm[1])*(x3[1]-xm[1])) * fact);

                          alph   = myAcos(((x3[0]-xm[0])*(x1[0]-xm[0])
                                         + (x3[1]-xm[1])*(x1[1]-xm[1])) * fact);

                          if (abs(alph12 + alph23 - alph) > 1.e-6) alph = 6.2831853 - alph;

                          e    = roundToInt(alph * 40. / 3.1415927);

                          if (triangleArea2D(x1,x2,x3) < 0.) alph = -alph;

                          fact = alph / (double) e;
                          alph = 0.;

                          spline[n].suppPnt.append((void*) &(point[iTmp[0]-1]));
                           
                          for (i=1; i<e; i++)
                          {
                            point.add(new GeomPoint);

                            alph += fact;
                            s = sin(alph);
                            c = cos(alph);

                            x3[0] = xm[0] + (x1[0]-xm[0]) * c - (x1[1]-xm[1]) * s;
                            x3[1] = xm[1] + (x1[1]-xm[1]) * c + (x1[0]-xm[0]) * s;

                            point[point.n-1].x[0] = x3[0];
                            point[point.n-1].x[1] = x3[1];

                            spline[n].suppPnt.append((void*) &(point[point.n-1]));

                            d = (x2[0]-x3[0]) * (x2[0]-x3[0]) + (x2[1]-x3[1]) * (x2[1]-x3[1]);

                            if (d < dmin) { dmin = d; imin = i; pmin = point.n-1; }
                          } 
                          spline[n].suppPnt.append((void*) &(point[iTmp[2]-1]));

                          point.del(pmin);
                          spline[n].suppPnt[imin] = (void*) &(point[iTmp[1]-1]);

                          break;

                 default: prgError(7,fct,"error in 'geometry splines'!");
               }
               n++;
             }

             for (i=0; i<nw; i++) word[i].free(); delete [] word;

             if (n < 1) prgError(8,fct,"error in 'geometry splines'!");

       	     break;

    case  2: cout << "     GEOMETRY: reading geometry surfaces ...\n\n";

             if (point.n < 1) prgError(1,fct,"'geometry points' must precede 'geometry surfaces'!");
             if (spline.n<1) prgError(2,fct,"'geometry splines' must precede 'geometry surfaces'!");
             if (surface.n > 0) prgError(3,fct,"multiple definition of 'geometry surfaces'!");

             n = 0;

             if (ndm == 2) prp = 1; else prp = 0;

             if (prp && grpType.n > 0) prgError(4,fct,
               "in 2D, 'geometry surfaces' must precede group type (e.g.'element groups')!");

             while (1)
             {
               line.getNextLine(Ifile);
               nw = line.split(&word);

               if (nw < 2 + prp)             break;
               if (!word[0].toInt(&i,false)) break;
               if (i != n + 1)               break;

               iTmp.free();
 
               i = 1 + prp;
               while (i < nw && word[i].toInt(iTmp.append(),false)) i++;
               if (i < nw) break;

               surface.add(new GeomSurface);

               if (prp) if (!word[1].toInt(grpType.append())) 
                 prgError(4,fct,"error in 'geometry surfaces'!");

               for (i=0; i<nw; i++) word[i].free(); delete [] word;

               for (i=0; i<nw-1-prp; i++)
               {
                 if (iTmp[i] < 1 || iTmp[i] > spline.n) 
                   prgError(5,fct,"error in 'geometry surfaces'!");

                 surface[n].suppSpln.append((void*) &(spline[iTmp[i]-1]));
               }

               n++;
             }

             for (i=0; i<nw; i++) word[i].free(); delete [] word;

             if (n < 1) prgError(6,fct,"error in 'geometry surfaces'!");

             break;

    case  3: cout << "     GEOMETRY: reading geometry boundary conditions ...\n\n";

             if (ndf<1) prgError(1,fct,"'dimensions' must precede 'geometry boundary conditions'!");
           
             if (surface.n == 0) prgError(1,fct,
                "'geometry point/spline/surface' must precede 'geometry boundary conditions'!");
 
	     sprintf(tmp,"%di 1s",ndm+ndf+1);  
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	     
             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
		     prgError(2,fct,"data error in 'geometry boundary conditions'!");
	     
             n = sTmp.dim();

	     for (i=0; i<n; i++)
             {
               ii = i*(1+ndm+ndf);

               bndCond.add(new GeomBndCond);

               for (j=0; j<ndm; j++) bndCond[i].idx.append(iTmp[ii+1+j]);
               for (j=0; j<ndf; j++) bndCond[i].idu.append(iTmp[ii+1+ndm+j]);
 
               switch (sTmp[i].which(geomObjectType))
               {
                 case  0: if (iTmp[ii] > point.n) goto error3;
                          bndCond[i].geomObj = (void*)(&(point[iTmp[ii]-1]));
                          break;
                 case  1: if (iTmp[ii] > spline.n) goto error3;
                          bndCond[i].geomObj = (void*)(&(spline[iTmp[ii]-1]));
                          break;
                 case  2: if (iTmp[ii] > surface.n) goto error3;
                          bndCond[i].geomObj = (void*)(&(surface[iTmp[ii]-1]));
                          break;
                 //case  3: if (iTmp[ii] > volume.n) goto error3;
                 //         bndCond[i].geomObj = (void*)(&(volume.[iTmp[ii]-1]));
                 //         break;
                 default: prgError(4,fct,"invalid geometry object type!");
               }
             }

             break;

             error3: prgError(4,fct,"data error in 'geometry boundary conditions'!");

    case  4: cout << "     GEOMETRY: reading geometry prescribed displacements ...\n\n";

             if (ndf<1) 
               prgError(1,fct,"'dimensions' must precede 'geometry prescribed displacements'!");

             if (surface.n == 0) prgError(1,fct,
              "'geometry point/spline/surface' must precede 'geometry prescribed displacements'!");

	     sprintf(tmp,"%di %df 1s",1+ndf,ndf);

             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);

             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);

             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
               prgError(2,fct,"data error in 'geometry prescribed displacements'!");

             n = sTmp.dim();

	     for (i=0; i<n; i++)
             {
               ii = i*(1+ndf);

               presDisp.add(new GeomPresDisp);

               for (j=0; j<ndf; j++) presDisp[i].tmFct.append(iTmp[ii+1+j]);
               for (j=0; j<ndf; j++) presDisp[i].uBase.append(dTmp[i*ndf+j]);
 
               switch (sTmp[i].which(geomObjectType))
               {
                 case  0: if (iTmp[ii] > point.n) goto error4;
                          presDisp[i].geomObj = (void*)(&(point[iTmp[ii]-1]));
                          break;
                 case  1: if (iTmp[ii] > spline.n) goto error4;
                          presDisp[i].geomObj = (void*)(&(spline[iTmp[ii]-1]));
                          break;
                 case  2: if (iTmp[ii] > surface.n) goto error4;
                          presDisp[i].geomObj = (void*)(&(surface[iTmp[ii]-1]));
                          break;
                 //case  3: if (iTmp[ii] > volume.n) goto error4;
                 //         presDisp[i].geomObj = (void*)(&(volume.[iTmp[i*(ndf+1)]-1]));
                 //         break;
                 default: prgError(5,fct,"invalid geometry object type!");
               }
             }

             break;

             error4: prgError(5,fct,"data error in 'geometry prescribed displacements'!");

    case  5: cout << "     GEOMETRY: reading geometry distributed loads ...\n\n";

             if (ndf<1) 
               prgError(1,fct,"'dimensions' must precede 'geometry distributed loads'!");

             if (surface.n == 0) prgError(1,fct,
              "'geometry point/spline/surface' must precede 'geometry distributed loads'!");

	     sprintf(tmp,"%di %df 1s",1+ndf,ndf);

             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);

             t1.initialise(tmpl);
	     t2.initialise(tmp);
             t1.expandToMatch(t2);

             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
               prgError(2,fct,"data error in 'geometry distributed loads'!");

             n = sTmp.dim();

	     for (i=0; i<n; i++)
             {
               distLoad.add(new GeomDistLoad);

               for (j=0; j<ndf; j++) distLoad[i].tmFct[j] = iTmp[i*(ndf+1)+1+j];
               for (j=0; j<ndf; j++) distLoad[i].fdl[j]   = dTmp[i*ndf+j];
 
               switch (sTmp[i].which(geomObjectType))
               {
                 case  0: if (iTmp[ii] > point.n) goto error5;
                          distLoad[i].geomObj = (void*)(&(point[iTmp[i*(ndf+1)]-1]));
                          break;
                 case  1: if (iTmp[ii] > spline.n) goto error5;
                          distLoad[i].geomObj = (void*)(&(spline[iTmp[i*(ndf+1)]-1]));
                          break;
                 case  2: if (iTmp[ii] > surface.n) goto error5;
                          distLoad[i].geomObj = (void*)(&(surface[iTmp[i*(ndf+1)]-1]));
                          break;
                 //case  3: if (iTmp[ii] > volume.n) goto error5;
                 //         distLoad[i].geomObj = (void*)(&(volume.[iTmp[i*(ndf+1)]-1]));
                 //         break;
                 default: prgError(3,fct,"invalid geometry object type!");
               }
             }

             break;

             error5: prgError(5,fct,"data error in 'geometry distributed loads'!");

    case  6: cout << "     GEOMETRY: reading geometry nodal data output ...\n\n";

             if (ndf<1) 
               prgError(1,fct,"'dimensions' must precede 'geometry nodal data output'!");

             if (surface.n == 0) prgError(1,fct,
              "'geometry point/spline/surface' must precede 'gemetry nodal data output'!");

             n = 0;

             while (1)
             { 	     
               line.getNextLine(Ifile);

               nw = line.split(&word); if (nw < 5) break;

	       tp = word[0].which(wrndTypeNames); if (tp < 0) break;
	       
	       if (!word[1].toInt(&indx,false)) break;

	       if (!word[2].toDbl(&fact,false)) break;

               wrndOutp.add(new GeomNodalDataOutput);

	       i = 0; 
	       while (i < nw-3)
               {
                 if (!word[i+3].toInt(iTmp,false)) break;
                 
                 switch (word[i+4].which(geomObjectType))
                 {
                   case  0: for (j=0; j<iTmp.n; j++)
                            {
                              if (iTmp[j] > point.n) goto error6;
                              if (!wrndOutp[n].geomObj.contains((void*)(&(point[iTmp[j]-1]))))
                                wrndOutp[n].geomObj.append((void*)(&(point[iTmp[j]-1])));
                            }
                            break;
                   case  1: for (j=0; j<iTmp.n; j++)
                            {
                              if (iTmp[j] > spline.n) goto error6;
                              if (!wrndOutp[n].geomObj.contains((void*)(&(spline[iTmp[j]-1]))))
                                wrndOutp[n].geomObj.append((void*)(&(spline[iTmp[j]-1])));
                            }
                            break;
                   case  2: for (j=0; j<iTmp.n; j++)
                            {
                              if (iTmp[j] > surface.n) goto error6;
                              if (!wrndOutp[n].geomObj.contains((void*)(&(surface[iTmp[j]-1]))))
                                wrndOutp[n].geomObj.append((void*)(&(surface[iTmp[j]-1])));
                            }
                            break;
                   default: prgError(3,fct,"invalid geometry object type!");
                 }
                 i += 2;
               }
               if (i < nw - 3) break;

               wrndOutp[n].fact = fact;
               wrndOutp[n].indx = indx;
               wrndOutp[n].type = tp;

               n++;
	     }

             delete [] word;

             break;

             error6: prgError(5,fct,"data error in 'geometry nodal data output'!");

    case -1: // go and inherit from DOMAIN
	     
	     this->Domain::readInputData(Ifile, line); 
	     
	     break;
  }
 
  return;
}








void Geometry::prepareInputData(void)
{
  // call ancestor function

  Domain::prepareInputData();

  
  cout << "     GEOMETRY: prepare input data ...\n\n";
  
  char fct[] = "Geometry::prepareInputData"; 
	  
  int i, j, n;

  for (i=0; i<surface.n; i++)  surface[i].initialise();

  // abort if some points are not used for splines

  for (i=0; i<point.n; i++) point[i].dat1 = 0;

  for (i=0; i<spline.n; i++) 
    for (j=0; j<spline[i].suppPnt.n; j++)
      if (((GeomPoint*)(spline[i].suppPnt[j]))->dat1 == 0) 
        ((GeomPoint*)(spline[i].suppPnt[j]))->dat1 = 1;

  for (i=0; i<point.n; i++) if (point[i].dat1 == 0) prgError(1,fct,"obsolete geometry points!");

  // abort if some splines are not used for surfaces

  for (i=0; i<spline.n; i++) spline[i].dat1 = 0;

  for (i=0; i<surface.n; i++) 
    for (j=0; j<surface[i].suppSpln.n; j++)
      if (((GeomSpline*)(surface[i].suppSpln[j]))->dat1 == 0) 
        ((GeomSpline*)(surface[i].suppSpln[j]))->dat1 = 1;

  for (i=0; i<spline.n; i++) if (spline[i].dat1 == 0) prgError(1,fct,"obsolete geometry splines!");

  // sort boundary conditions according to inverse priority

  n = bndCond.n; 

  i = 0; while (i < n)
  { if (bndCond[i].whichSurface() != NULL) { bndCond.move(i,bndCond.n-1); n--; } else i++; }

  i = 0; while (i < n)
  { if (bndCond[i].whichSpline() != NULL) { bndCond.move(i,bndCond.n-1); n--; } else i++; }

  i = 0; while (i < n)
  { if (bndCond[i].whichPoint() != NULL) { bndCond.move(i,bndCond.n-1); n--; } else i++; }

  return;
}








void Geometry::plotGeometry(unsigned int bitFlag)
{
  int i;

  if (bitFlag & pow2[0]) 
  {
    if (bitFlag & pow2[1])
      for (i=0; i<point.n; i++) point[i].draw(i+1);
    else
      for (i=0; i<point.n; i++) point[i].draw(0);
  }

  if (bitFlag & pow2[2]) 
  {
    if (bitFlag & pow2[3])
      for (i=0; i<spline.n; i++) spline[i].draw(i+1);
    else
      for (i=0; i<spline.n; i++) spline[i].draw(0);
  }

  if (bitFlag & pow2[4]) 
  {
    if (bitFlag & pow2[5])
      for (i=0; i<surface.n; i++) surface[i].draw(i+1);
    else
      for (i=0; i<surface.n; i++) surface[i].draw(0);
  }

  return;
}








void Geometry::findGeomMinMaxX(double *xmn, double *xmx)
{
  if (point.n < 1) 
  {
    prgWarning(1,"Geometry::findGeomMinMaxX","no geometry data available!"); 
    return;
  }

  int i, j;

  for (j=0; j<ndm; j++)
  {
    xmn[j] = point[0].x[j];
    xmx[j] = point[0].x[j];
  }

  for (i=1; i<point.n; i++)
    for (j=0; j<ndm; j++)
    {
      if (xmn[j] > point[i].x[j]) xmn[j] = point[i].x[j];
      if (xmx[j] < point[i].x[j]) xmx[j] = point[i].x[j];
    }

  return;
}










void Geometry::whichSurfaces(Vector<int> &surf, Vector<int> &spln)
{
  prgError(1,"Geometry::whichSurfaces","3D Geometry not yet implemented!");

  return;
}










void Geometry::whichSplines(Vector<int> &surf, Vector<int> &spln)
{
  int i, j, k;

  spln.free();

  Vector<void*> splnPtr;

  for (i=0; i<spline.n; i++) splnPtr.append((void*)(&(spline[i])));

  for (i=0; i<surface.n; i++) 
  {
    if (surface[i].hasBeenDiscretised)
    {
      j = 0;
      while (j < surf.n) if (surf[j] != i) j++; else break;
      if (j == surf.n)
      {
        for (j=0; j<surface[i].suppSpln.n; j++)
        {
          k = 0;
          while (k < splnPtr.n && splnPtr[k] != surface[i].suppSpln[j]) k++;
          if (k < splnPtr.n) splnPtr.del(k);
        }
      }
    }
  }

  for (i=0; i<splnPtr.n; i++) spln.append(spline.getIndex((GeomSpline*)(splnPtr[i])));

  return;
}










void Geometry::printComputerTime(bool reset, int detailFlg)
{
  Domain::printComputerTime(reset,detailFlg);

  return;
}

