

#include <cmath>

#include "ObjectSurface.h"
#include "Plot.h"
#include "Mesh.h"
#include "FunctionsProgram.h"
#include "Element.h"
#include "ComputerTime.h"


static const unsigned int pow2[8] = {1,2,4,8,16,32,64,128}, all = 255;


extern ComputerTime computerTime;
extern Plot plot;






void SurfaceEdge::setVisibility(BasicFace **face, bool vis)
{
  if (vis)
  {
    face[f1]->edgeBits = face[f1]->edgeBits | pow2[j1+3];
    face[f2]->edgeBits = face[f2]->edgeBits | pow2[j2+3];
  }
  else
  {
    face[f1]->edgeBits = face[f1]->edgeBits & (all-pow2[j1+3]);
    face[f2]->edgeBits = face[f2]->edgeBits & (all-pow2[j2+3]);
  }
  return;
}







bool SurfaceEdge::isKink(BasicFace **face, float cosAlph)
{
  if ( abs(     dot3(face[f1]->normal,face[f2]->normal)
         / sqrt(dot3(face[f1]->normal,face[f1]->normal)*
                dot3(face[f2]->normal,face[f2]->normal))) < cosAlph) return true;

  return false;
}







ObjectSurface::ObjectSurface(void)
{
  view     = plot.perspective.view;
  observer = plot.perspective.observer;

  bnd2face      = NULL;
  edge          = NULL;
  face          = NULL;
  x             = NULL;
  picx          = NULL;
  u             = NULL;
  scrx          = NULL;
  nodeIsVisible = NULL;
  faceIsVisible = NULL;
  perm          = NULL;

  updtFlagX        = true;
  updtFlagU        = true;
  newPersFlag      = true;
  showBackFaces    = false;
  showDeformed     = false;
  carefulCheckDone = false;

  // this object surface is not associated with a specific domain
  // (it is probably the bounding box belonging to Perspective)

  dom = NULL;

  return;
}






ObjectSurface::ObjectSurface(Domain *domain)
{
  view     = plot.perspective.view;
  observer = plot.perspective.observer;

  bnd2face      = NULL;
  edge          = NULL;
  face          = NULL;
  x             = NULL;
  picx          = NULL;
  u             = NULL;
  scrx          = NULL;
  nodeIsVisible = NULL;
  faceIsVisible = NULL;
  perm          = NULL;

  updtFlagX        = true;
  updtFlagU        = true;
  newPersFlag      = true;
  showBackFaces    = false;
  showDeformed     = false;
  carefulCheckDone = false;

  // associate this object surface with Domain
  // and add it to the list of domain surfaces handled by Perspective

  dom = domain;

  plot.perspective.addObjSurf(this);

  return;
}






ObjectSurface::~ObjectSurface()
{
  free();

  return;
}





void ObjectSurface::free(void)
{
  paintSequence.free();
 
  dblSortValue.free();

  faceList.free();

  if (face != NULL) delete [] face; face = NULL;

  if (bnd2face != NULL) delete [] bnd2face; bnd2face = NULL;

  if (edge != NULL) delete [] edge; edge = NULL;

  if (x != NULL) delete [] x; x = NULL;

  if (scrx != NULL) delete [] scrx; scrx = NULL;

  if (picx != NULL) delete [] picx; picx = NULL;

  if (u != NULL) delete [] u; u = NULL;

  if (nodeIsVisible != NULL) delete [] nodeIsVisible; nodeIsVisible = NULL;

  if (faceIsVisible != NULL) delete [] faceIsVisible; faceIsVisible = NULL;

  if (perm != NULL) delete [] perm;

  return;
}







void ObjectSurface::addElementFace(int e, int iFace, Vector<int> &face, VectorArray<int> &nd2bnd)
{
  Element *elem = ((Mesh*)dom)->elem[e];

  int i, j, n = elem->nBasicFacesPerFace(), m = faceList.n, *IX;

  unsigned int *edgeBits;

  Vector<int> lFace;

  for (i=0; i<face.n; i++) lFace.append(nd2bnd[face[i]-1]);

  for (i=0; i<n; i++) 
  {
    faceList.add(new BasicFace);

    edgeBits = &(faceList[m].edgeBits);
    IX       = faceList[m++].ix;

    elem->defineBasicFace(iFace,i,IX,edgeBits);

    for (j=0; j<3; j++) IX[j] = lFace[IX[j]];
  }

  return;
}







void ObjectSurface::finalise(Vector<int> &bndNd)
{
  if (face != NULL) prgError(1,"ObjectSurface::finalise","face != NULL !");

  int i, j, k, l;

  // faces

  numfc = faceList.n;

  face = new BasicFace* [numfc];

  paintSequence.free();
  dblSortValue.free();

  for (i=0; i<numfc; i++)
  {
    face[i] = &(faceList[i]);
    paintSequence.append(i);
    dblSortValue.append();
  }

  // allocate memory

  numnp = bndNd.n;

  x = new float [3*numnp];

  picx = new float [2*numnp];

  scrx = new int [2*numnp];

  u = new float [numnp];

  nodeIsVisible = new char [numnp];

  faceIsVisible = new int [numfc];

  perm = new int [numfc];

  // generate node to basic face connectivity

  bnd2face = new Vector<int> [bndNd.n];

  for (i=0; i<numfc; i++)

    for (j=0; j<3; j++)

      bnd2face[faceList[i].ix[j]].append(i);

  // order the pointers to faces in bnd2face

  for (i=0; i<numnp; i++)
  {
    for (j=0; j<bnd2face[i].n-1; j++)
    {
      k = j+1;
      while (face[bnd2face[i][j]]->n2(i) != face[bnd2face[i][k]]->n1(i)) k++;
      if (k > j+1)
      {
        l                = bnd2face[i][k];
        bnd2face[i][k]   = bnd2face[i][j+1];
        bnd2face[i][j+1] = l;
      }
    }
  }

  return;
}







void ObjectSurface::generateBoundingBox(double *Xmn, double *Xmx)
{
  int i, j;

  float *X, xmn[3] = { (float)Xmn[0], (float)Xmn[1], (float)Xmn[2] },
            xmx[3] = { (float)Xmx[0], (float)Xmx[1], (float)Xmx[2] };

  // delete all stuff that is currently there

  free();

  // allocate memory

  numnp = 8;

  numfc = 12;

  face = new BasicFace* [numfc];

  for (i=0; i<numfc; i++)
  {
    faceList.add(new BasicFace);

    face[i] = &(faceList[i]);

    paintSequence[i] = i;
  }

  x = new float [3*numnp];

  picx = new float [2*numnp];

  scrx = new int [2*numnp];

  u = new float [numnp];

  nodeIsVisible = new char [numnp];

  faceIsVisible = new int [numfc];

  // generate surface of a box

  face[ 0]->ix[0] = 0; face[ 0]->ix[1] = 1; face[ 0]->ix[2] = 4; face[ 0]->edgeBits = 5; //101
  face[ 1]->ix[0] = 1; face[ 1]->ix[1] = 2; face[ 1]->ix[2] = 5; face[ 1]->edgeBits = 5;
  face[ 2]->ix[0] = 2; face[ 2]->ix[1] = 3; face[ 2]->ix[2] = 6; face[ 2]->edgeBits = 5;
  face[ 3]->ix[0] = 3; face[ 3]->ix[1] = 0; face[ 3]->ix[2] = 7; face[ 3]->edgeBits = 5;
  face[ 4]->ix[0] = 4; face[ 4]->ix[1] = 1; face[ 4]->ix[2] = 5; face[ 4]->edgeBits = 6; //011
  face[ 5]->ix[0] = 5; face[ 5]->ix[1] = 2; face[ 5]->ix[2] = 6; face[ 5]->edgeBits = 6;
  face[ 6]->ix[0] = 6; face[ 6]->ix[1] = 3; face[ 6]->ix[2] = 7; face[ 6]->edgeBits = 6;
  face[ 7]->ix[0] = 7; face[ 7]->ix[1] = 0; face[ 7]->ix[2] = 4; face[ 7]->edgeBits = 6;
  face[ 8]->ix[0] = 3; face[ 8]->ix[1] = 1; face[ 8]->ix[2] = 0; face[ 8]->edgeBits = 6;
  face[ 9]->ix[0] = 1; face[ 9]->ix[1] = 3; face[ 9]->ix[2] = 2; face[ 9]->edgeBits = 6;
  face[10]->ix[0] = 4; face[10]->ix[1] = 5; face[10]->ix[2] = 7; face[10]->edgeBits = 5;
  face[11]->ix[0] = 7; face[11]->ix[1] = 5; face[11]->ix[2] = 6; face[11]->edgeBits = 6;
  
  x[ 0] = xmn[0]; x[ 1] = xmn[1]; x[ 2] = xmn[2];
  x[ 3] = xmx[0]; x[ 4] = xmn[1]; x[ 5] = xmn[2];
  x[ 6] = xmx[0]; x[ 7] = xmx[1]; x[ 8] = xmn[2];
  x[ 9] = xmn[0]; x[10] = xmx[1]; x[11] = xmn[2];
  x[12] = xmn[0]; x[13] = xmn[1]; x[14] = xmx[2];
  x[15] = xmx[0]; x[16] = xmn[1]; x[17] = xmx[2];
  x[18] = xmx[0]; x[19] = xmx[1]; x[20] = xmx[2];
  x[21] = xmn[0]; x[22] = xmx[1]; x[23] = xmx[2];

  for (i=0; i<numfc; i++) for(j=0; j<3; j++) face[i]->normal[j] = 0.;

  face[ 0]->normal[1] = -1.;
  face[ 1]->normal[0] = +1.;
  face[ 2]->normal[1] = +1.;
  face[ 3]->normal[0] = -1.;
  face[ 4]->normal[1] = -1.;
  face[ 5]->normal[0] = +1.;
  face[ 6]->normal[1] = +1.;
  face[ 7]->normal[0] = -1.;
  face[ 8]->normal[2] = -1.;
  face[ 9]->normal[2] = -1.;
  face[10]->normal[2] = +1.;
  face[11]->normal[2] = +1.;

  newPersFlag = true;

  return;
}






void ObjectSurface::draw(int col1, int col2, bool edgesOnly, bool backFaceFlag, bool defmFlag)
{
  updateVisible(backFaceFlag,defmFlag);

  for (int i=0; i<nPaint; i++)

    face[paintSequence[i]]->simplePaint(col1,col2,scrx,picx,edgesOnly);

  return;
}





#define CLIP_MARGIN 200

void ObjectSurface::updateVisible(bool backFaceFlag, bool defmFlag)
{
  int i, ii, j, k, *IX, *SX = scrx, ix, ix2, ix3,
      mn  = - CLIP_MARGIN,
      mxx = plot.wPix + CLIP_MARGIN,
      mxy = plot.hPix + CLIP_MARGIN;

  char *vs;

  float z, *PX;

  if (!newPersFlag) newPersFlag = (backFaceFlag != showBackFaces);

  if (!updtFlagX) updtFlagX = (defmFlag != showDeformed);

  showBackFaces = backFaceFlag;

  showDeformed = defmFlag;

  if (dom == NULL) updtFlagX = false;

  if (updtFlagX)
  {
    updateCoor();

    updtFlagX = false;
  }
  else if (!newPersFlag) return;

  newPersFlag = false;

  carefulCheckDone = false;

  // assemble all faces except those lying behind the observer
  // or facing the wrong way in paintSequence

  nPaint = numfc;

  i = 0;

  while (i < nPaint)
  {
    ii = paintSequence[i];
 
    if (backFaceFlag != face[ii]->facingOK(observer,x))
    {
      z = face[ii]->calcZDepth(x);

      if (z > 0.) dblSortValue[i++] = z;

      else
      {
        dblSortValue.del(i);   dblSortValue.append(-1.);
        paintSequence.del(i);  paintSequence.append(ii);  nPaint--;
      }
    }
    else
    {
      dblSortValue.del(i);   dblSortValue.append(-1.);
      paintSequence.del(i);  paintSequence.append(ii);  nPaint--;
    }
  }

  // calculate the screen coordinates of the nodes that define these faces

  for (i=0; i<numnp; i++) nodeIsVisible[i] = -1;  // not checked; if left at -1 then not visible

  i = 0;

  while (i < nPaint)
  {
    ii = paintSequence[i];

    IX = face[ii]->ix;

    k = 0;

    while (k < 3)
    {
      ix = *IX;
      vs = nodeIsVisible + ix;

      if (*vs == -1)
      {
        ix2 = ix  + ix;
        ix3 = ix2 + ix;
        
        SX = scrx + ix2;
        PX = picx + ix2;
      
        plot.perspective.pictureCoor(x + ix3, PX);
      
        plot.xy2D(PX,SX);
      
        if (*SX < mn || *SX > mxx || *(SX+1) < mn || *(SX+1) > mxy)
        { 
          dblSortValue.del(i);   dblSortValue.append(-1.);
          paintSequence.del(i);  paintSequence.append(ii);  nPaint--; i--;
      
          *vs = 0;    // not visible
      
          break;
        }
        else *vs = 2; // probably visible
      }
      k++; IX++;
    }
    i++;
  }

  // sort according to z depth

  //computerTime.go("quickSort");

  for (i=0; i<nPaint; i++) perm[i] = paintSequence[i];

  dblSortValue.quickSort(true,perm,0,nPaint-1);

  for (i=0; i<nPaint; i++) paintSequence[i] = perm[i];

  //computerTime.stopAndPrint("quickSort");

  // set face visibility

  for (i=0; i<numfc; i++) faceIsVisible[i] = -1;

  for (i=0; i<nPaint; i++) faceIsVisible[paintSequence[i]] = i;

  return;
}







void ObjectSurface::updateCoor(void)
{
  int *bndNd = ((Mesh*)dom)->bndNd;

  int i, j, indm = 0;

  float  *X = x;
  double *domX = ((Mesh*)dom)->x.x;
 
  if (!showDeformed) domX = ((Mesh*)dom)->x0.x;

  for (i=0; i<numnp; i++) for (j=0; j<3; j++) *(X++) = (float)domX[(bndNd[i]-1)*3+j];

  for (i=0; i<numfc; i++) face[i]->calcNormal(x);

  return;
}






bool ObjectSurface::getNeighboursOnSameGeometry(int j0, Vector<int> &tmp, float cosAlph)
{
  float dot[30], *n1, *n2, rn1, rn2;

  int c = 0, i, nd[30], n = bnd2face[j0].n;

  n1 = face[bnd2face[j0][n-1]]->normal;

  rn1 = n1[0]*n1[0] + n1[1]*n1[1] + n1[2]*n1[2];

  for (i=0; i<n; i++)
  {
    n2 = face[bnd2face[j0][i]]->normal;

    rn2 = n2[0]*n2[0] + n2[1]*n2[1] + n2[2]*n2[2];

    dot[i] = (n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]) / sqrt(rn1 * rn2);

    nd[i] = face[bnd2face[j0][i]]->n1(j0);

    if (dot[i] < cosAlph) c++;

    n1 = n2;

    rn1 = rn2;
  }
  if (c == 0) { for (i=0; i<n; i++) tmp.append(nd[i]); return true; }

  else if (c == 2) 
  { 
    for (i=0; i<n; i++) 

      if (dot[i] < cosAlph) tmp.append(nd[i]);

    return true; 
  }
 
  return false;
}







void ObjectSurface::plotNodePoint(int i, float dpt, int nb)
{
  if (!nodeIsReallyVisible(i)) return;

  float d = dpt; if (d < 0.) d = plot.dPt();

  int sx[2], sd[2];

  char txt[30];

  sd[0] = roundToInt(d * plot.wFactAct);
  sd[1] = roundToInt(d * plot.hFactAct);

  sx[0] = scrx[i+i];
  sx[1] = scrx[i+i+1];

  essGrpFillCircle(sx,sd);

  if (nb > -1) 
  {
    sx[0] += roundToInt(0.6 * (float)sd[0]);
    sx[1] -= roundToInt(0.6 * (float)sd[1]);

    sprintf(txt,"%d",nb); essGrpPutText(sx,txt,1);
  }

  if (!plot.psOpen) return;

  cout << "ObjectSurface::plotNodePoint: Implement it!!!\n\n"; 

  return;
} 








bool ObjectSurface::nodeIsReallyVisible(int i)
{
  if (nodeIsVisible[i] < 1) return false;  // not visible

  if (nodeIsVisible[i] == 1) return true;  // definitely visible

  // (nodeIsVisible == 2)                  // possibly visible -> check properly

  int k = nPaint - 1;
  float *px = picx+i+i;

  if (!showBackFaces)
  {
    while (!bnd2face[i].contains(paintSequence[k]))

      if (face[paintSequence[k]]->frontContains(picx,px)) { nodeIsVisible[i] = 0; return false; }

      else k--;
  }
  else
  {
    while (!bnd2face[i].contains(paintSequence[k]))

      if (face[paintSequence[k]]->backContains(picx,px)) { nodeIsVisible[i] = 0; return false; }

      else k--;
  }
  nodeIsVisible[i] = 1;

  return true;
}









void ObjectSurface::updateU(int dataType, int indx)
{
  char fct[] = "ObjectSurface::updateU";

  int *bndNd = ((Mesh*)dom)->bndNd;

  int i, j, inc,
      ndf = ((Mesh*)dom)->ndf,
      ndm = ((Mesh*)dom)->ndm;

  double *domU;
 
  switch (dataType)
  {
    case 1: // degree of freedom

            if (indx > ndf) prgError(1,fct,"invalid indx");

            domU = ((Mesh*)dom)->u.x + indx - 1;

            inc = ndf;

            break;

    case 2: // initial coordinate

            if (indx > ndm) prgError(2,fct,"invalid indx");

            domU = ((Mesh*)dom)->x0.x + indx - 1;

            inc = ndm;

            break;

    case 3: // current coordinate

            if (indx > ndm) prgError(3,fct,"invalid indx");

            domU = ((Mesh*)dom)->x.x + indx - 1;

            inc = ndm;

            break;

    case 4: // last projection

            domU = ((Mesh*)dom)->outp.x;

            inc = 1;

            break;
  }

  for (i=0; i<numnp; i++) u[i] = (float)domU[(bndNd[i]-1)*inc];

  umn = u[0];
  umx = u[0];

  for (i=0; i<numnp; i++) { umn = min(umn,u[i]); umx = max(umx,u[i]); }

  return;
}







void ObjectSurface::contourPlot(int nCol)
{
  updateVisible(showBackFaces,showDeformed);

  if (!carefulCheckDone) carefulCheck();

  for (int i=0; i<nPaint; i++)

    face[paintSequence[i]]->contourPlot(nCol,x,u,umn,umx,scrx,picx,showBackFaces);

  return;
}






void ObjectSurface::carefulCheck(void)
{
  char fct[] = "ObjectSurface::carefulCheck";

  int i, j, k, l, n0, n1, f1, f2;

  float cosAlph = .5;

  // remove hidden surfaces from paintSequence

  for (i=0; i<nPaint; i++)
  {
    k = paintSequence[i];

    j = 0; while (j < 3) if (nodeIsReallyVisible(face[k]->ix[j])) break; else j++;

    if (j == 3) 
    {
      dblSortValue.del(i);   dblSortValue.append(-1.);
      paintSequence.del(i);  paintSequence.append(k);  nPaint--;

      faceIsVisible[k] = -1;
    }
    else 
    {
      faceIsVisible[k] = i;
    }    
  }

  // generate edge array if it has not been done before

  if (edge == NULL)
  {
    numed = 0;

    edge = new SurfaceEdge [numfc * 3 / 2];

    for (n0=0; n0<numnp; n0++)
    {
      f2 = bnd2face[n0][bnd2face[n0].n-1];

      for (i=0; i<bnd2face[n0].n; i++)
      {
        f1 = bnd2face[n0][i];

        n1 = face[f1]->n1(n0);

        if (n1 > n0)
        {
          l = 0; while (face[f1]->ix[l] != n0) l++;
          edge[numed].j1 = l;

          l = 0; while (face[f2]->ix[l] != n1) l++; 
          edge[numed].j2 = l;

          edge[numed].f1 = f1;
          edge[numed++].f2 = f2;
        }
        f2 = f1;
      }
    }

    if (numed != numfc * 3 / 2) prgError(1,fct,"numed != numfc * 3 / 2");
  }

  // set visible edge bits

  for (i=0; i<numed; i++)
  {
    f1 = edge[i].f1;
    f2 = edge[i].f2;

    if (faceIsVisible[f1] < 0 && faceIsVisible[f2] < 0) edge[i].setVisibility(face,false);

    else if (faceIsVisible[f1] < 0 || faceIsVisible[f2] < 0) edge[i].setVisibility(face,true);

    else if (edge[i].isKink(face,cosAlph)) edge[i].setVisibility(face,true); 

    else edge[i].setVisibility(face,false);
  }

  carefulCheckDone = true;

  return;
}









void ObjectSurface::writeInterfaceInterpolations(VectorArray<int> &all2free, ofstream &out)
{
  int c = 0, i, j, *bndNd = ((Mesh*)dom)->bndNd;

  for (i=0; i<numfc; i++)
  {
    j = 0; while (j < 3) if (all2free[bndNd[face[i]->ix[j]]-1] < 1) break; else j++;

    if (j == 3) 
    {
      c++;

      for (j=0; j<3; j++) out << " " << all2free[bndNd[face[i]->ix[j]]-1];
      out << "\n";
    }
  }

  COUT << c << " linear interpolations found from " << numfc << " triangular surfaces.\n\n";

  return;
}

