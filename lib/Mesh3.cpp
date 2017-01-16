
#include "Mesh.h"
#include "Plot.h"
#include "PropertyTypeEnum.h"
#include "MathGeom.h"


extern Plot plot;








  



int Mesh::pointInElement(double *xp)
{
  int i;

  for (i=0; i<numel; i++) elem[i]->flag = false;
  for (i=0; i<numnp; i++) nodeFlag[i] = false;

  return pointInElement_hlp(guessNearestNode(lastSearchNode,xp),xp);
}













int Mesh::pointInElement_hlp(int nd, double *xp)
{
  //cout << nd << "\n";

  int i, e, ind = nd - 1;

  // check all elements attached to node nd which have not yet been visited

  for (i=0; i<nodeElem[ind].n; i++)
  {
    e = nodeElem[ind][i];
    if (!elem[e]->flag)
    {
      if (elem[e]->containsPoint(xp,NULL)) 
      {
        //cout << "\n";

        lastSearchNode = nd;
        return e + 1;
      }
      else elem[e]->flag = true;
    }
  }
  nodeFlag[ind] = true;

  // generate list of nodes attached to node nd, which have not yet been visited

  Vector<int>      nodeTmp1;
  VectorArray<int> nodeTmp;

  for (i=0; i<nodeNode[ind].n; i++) 
    if (!nodeFlag[nodeNode[ind][i]-1]) nodeTmp1.append(nodeNode[ind][i]);

  if (nodeTmp1.n == 0) return -1;

  nodeTmp = nodeTmp1; nodeTmp1.free();

  // sort this list according to direction

  Vector<double> dotProd;

  double *x0 = &(x.x[ind+ind]), *xi;

  for (i=0; i<nodeTmp.n; i++)
  {
    xi = &(x.x[nodeTmp[i]+nodeTmp[i]-2]);
    dotProd.append(((xi[0]-x0[0])*(xp[0]-x0[0])+(xi[1]-x0[1])*(xp[1]-x0[1])) 
             / sqrt((xi[0]-x0[0])*(xi[0]-x0[0])+(xi[1]-x0[1])*(xi[1]-x0[1])));
  }
  dotProd.quickSort(true, nodeTmp.x);

  // visit all nodes in the list

  for (i=0; i<nodeTmp.n; i++) 
  {
    e = pointInElement_hlp(nodeTmp[i],xp);
    if (e > -1) return e;
  }
  return -1;
}












int Mesh::guessNearestNode(int nd, double *xp)
{
  int i, k, imn0 = nd - 1, imn = imn0;

  double d = dist2(x.x+imn0*ndm,xp,ndm), dtrial;

  while (1)
  {
    for (i=0; i<nodeNode[imn0].n; i++)
    {
      k = nodeNode[imn0][i] - 1;
      dtrial = dist2(xp,x.x+k*ndm,ndm);
      if (dtrial < d * .999) { imn = k; d = dtrial; }
    }
    if (imn0 == imn) return imn + 1;

    imn0 = imn;
  }
}












void Mesh::smoothElemSizeDist_hlp0(void)
{
  int i;

  for (i=0; i<point.n; i++) point[i].h = 1.e20;

  for (i=0; i<spline.n; i++) spline[i].checkElemSizes();

  for (i=0; i<point.n; i++) elemSizeOpt.x[point[i].dat1] = point[i].h;

  return;
}










void Mesh::smoothElemSizeDist_hlp1(int nCircle, double maxGrad, double hmin)
{
  // ensure resolution of curved boundaries (->nCircle element edges on a circle)

  int e, i, j, j1, j2, k;

  double *x0, *x1, *x2, h, fact = pi / (double)nCircle, fact1;

  Vector<void*> *suppPntPtr;

  if (geometryDiscretisedAtLeastOnce)
  {
    for (i=0; i<spline.n; i++)
    {
      suppPntPtr = &(spline[i].suppPnt);

      if (suppPntPtr->n > 2)
      {
        x0 = &(x.x[((GeomPoint*)((*suppPntPtr)[0]))->dat1 * 2]);
      
        j1 = ((GeomPoint*)((*suppPntPtr)[1]))->dat1;
      
        x1 = &(x.x[j1+j1]);
      
        for (j=2; j<suppPntPtr->n; j++)
        {
          j2 = ((GeomPoint*)((*suppPntPtr)[j]))->dat1;
          x2 = &(x.x[j2+j2]);
      
          fact1 = abs(x0[0]*(x2[1]-x1[1]) + x1[0]*(x0[1]-x2[1]) + x2[0]*(x1[1]-x0[1]));
      
          if (fact1 > 1.e-12)
          {
            h = sqrt(((x0[0]-x1[0])*(x0[0]-x1[0]) + (x0[1]-x1[1])*(x0[1]-x1[1]))
                   * ((x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]))
                   * ((x2[0]-x0[0])*(x2[0]-x0[0]) + (x2[1]-x0[1])*(x2[1]-x0[1])))
                     / fact1 * fact;
      
            if (h < hmin) h = hmin;
      
            if (elemSizeOpt.x[j1] > h) elemSizeOpt.x[j1] = h;
          }
          j1 = j2;
          x0 = x1;
          x1 = x2;
        }
      }
    }
  }
  else
  {
    for (i=0; i<spline.n; i++)
    {
      suppPntPtr = &(spline[i].suppPnt);

      x0 = ((GeomPoint*)((*suppPntPtr)[0]))->x;

      x1 = ((GeomPoint*)((*suppPntPtr)[1]))->x;

      for (j=2; j<suppPntPtr->n; j++)
      {
        x2 = ((GeomPoint*)((*suppPntPtr)[j]))->x;

        fact1 = abs(x0[0]*(x2[1]-x1[1]) + x1[0]*(x0[1]-x2[1]) + x2[0]*(x1[1]-x0[1]));

        if (fact1 > 1.e-12)
        {
          h = sqrt(((x0[0]-x1[0])*(x0[0]-x1[0]) + (x0[1]-x1[1])*(x0[1]-x1[1]))
                 * ((x1[0]-x2[0])*(x1[0]-x2[0]) + (x1[1]-x2[1])*(x1[1]-x2[1]))
                 * ((x2[0]-x0[0])*(x2[0]-x0[0]) + (x2[1]-x0[1])*(x2[1]-x0[1])))
                   / fact1 * fact;

          if (h < hmin) h = hmin;

          if (getElemSizeOptInternal(x1) > 1.001*h) hPntSrc.add(new HPointSource(x1,h,maxGrad));
        }
        x0 = x1;
        x1 = x2;
      }
    }
  }
  return;
}








void Mesh::smoothElemSizeDist_hlp2(int nAcross, double maxGrad, double hmin)
{
  // ensure at least nAcross elements at narrow parts of the domain

  int e, i, j, j0, j1, j2, k, k0, l, m, mnm, iLoop;

  double *x0, *x1, *x2, *xk, *xk0, *xk2, *xx, xd[2], x0d[2], x1d[2], w[5], 
         fact = (double)nAcross, fact1, d, h, h0, h1, mnd, rad2;

  bool chDirFlg;

  Vector<void*> *suppPntPtr, allPnt, seqPnt;

  GeomSurface *surfPtr;

  fact1 = 1. / fact; 
  if      (fact1 > .50 * maxGrad) fact1 -= .25 * maxGrad;
  else if (fact1 > .25 * maxGrad) fact1  = .25 * maxGrad;

  if (geometryDiscretisedAtLeastOnce)
  {
    for (i=0; i<spline.n; i++)
    {
      suppPntPtr = &(spline[i].suppPnt);

      x0 = &(x.x[((GeomPoint*)((*suppPntPtr)[0]))->dat1 * 2]);

      j1 = ((GeomPoint*)((*suppPntPtr)[1]))->dat1;
      j0 = j1;

      x1 = &(x.x[j1+j1]);

      for (j=2; j<suppPntPtr->n; j++)
      {
        rad2 = elemSizeOpt.x[j1] * fact; rad2 *= rad2;

        j2 = ((GeomPoint*)((*suppPntPtr)[j]))->dat1;
        x2 = &(x.x[j2+j2]);

        w[0]  = x0[1] - x2[1];
        w[1]  = x2[0] - x0[0];
        w[2]  = sqrt(w[0]*w[0]+w[1]*w[1]);
        w[3]  = 1. / w[2];
        w[0] *= w[3];
        w[1] *= w[3];
        w[4]  = .4  * ((w[0]*(x0[0]-x1[0])+w[1]*(x0[1]-x1[1]))
                     /sqrt((x0[0]-x1[0])*(x0[0]-x1[0])+(x0[1]-x1[1])*(x0[1]-x1[1]))
                     +(w[0]*(x2[0]-x1[0])+w[1]*(x2[1]-x1[1]))
                     /sqrt((x2[0]-x1[0])*(x2[0]-x1[0])+(x2[1]-x1[1])*(x2[1]-x1[1]))) - .001;

        chDirFlg = true;

        shoot:

        //plot.wipe();
        //plot.setColour(1);
        //plotMesh(true,false);
        //plot.setColour(0);
        //plot.point(x1,plot.dPt(),j1+1);

        k   = j1;
        xk  = x1;
        k0  = -1;

        while (1)
        {
          mnm = -1;
          mnd = 0.;
          for (l=0; l<nodeNode[k].n; l++)
          {
            m = nodeNode[k][l] - 1;
            if (m != k0)
            {
              xx = &(x.x[m+m]);
              d  = abs(w[0]*(x1[1]-xx[1])-w[1]*(x1[0]-xx[0]));
              if ((d < mnd || mnm == -1) && w[0]*(xx[0]-xk[0])+w[1]*(xx[1]-xk[1]) > 0.2*d)
                { mnd = d; mnm = m; }
            }
          }
          if (mnm == -1) break;
          k0 = k;
          k  = mnm;
          xk = &(x.x[k+k]);

          //plot.setColour(3);
          //plot.point(xk,plot.dPt(),0);//k+1);

          if (ndGmLnk[k].whichPoint() != NULL) break;

          if (dist2(xk,x1,2) > rad2) { mnm = -1; break; }
        }
        if (mnm > -1 && k != j1 && k != j2 && k != j0)
        {
          l = 0; while (ndGmLnk[nodeNode[k][l]-1].whichPoint() == NULL) l++;
          xk0 = &(x.x[(nodeNode[k][l++]-1)*2]);
          while (ndGmLnk[nodeNode[k][l]-1].whichPoint() == NULL) l++;
          xk2 = &(x.x[(nodeNode[k][l]-1)*2]);

          //cout << (w[0]*(xk0[0]-xk[0])+w[1]*(xk0[1]-xk[1])) /
          //  (w[2]*sqrt((xk0[0]-xk[0])*(xk0[0]-xk[0])+(xk0[1]-xk[1])*(xk0[1]-xk[1]))) << ","
          //  << (w[0]*(xk2[0]-xk[0])+w[1]*(xk2[1]-xk[1])) /
          //  (w[2]*sqrt((xk2[0]-xk[0])*(xk2[0]-xk[0])+(xk2[1]-xk[1])*(xk2[1]-xk[1]))) << "\n";

          if ((w[0]*(xk0[0]-xk[0])+w[1]*(xk0[1]-xk[1])) /
            sqrt((xk0[0]-xk[0])*(xk0[0]-xk[0])+(xk0[1]-xk[1])*(xk0[1]-xk[1])) > w[4]
           && (w[0]*(xk2[0]-xk[0])+w[1]*(xk2[1]-xk[1])) /
            sqrt((xk2[0]-xk[0])*(xk2[0]-xk[0])+(xk2[1]-xk[1])*(xk2[1]-xk[1])) > w[4])
          {
            //plot.setColour(4);
            //plot.point(xk0);
            //plot.point(xk2);

            h = fact1 * sqrt(dist2(x1,xk,2));
            
            if (h < hmin) h = hmin;

            if (elemSizeOpt.x[j1] > h) 
            {
              elemSizeOpt.x[j1] = h;

              //COUT << "optimal element size reduced by 'across' criteria!\n\n";

              //prgUpdateDisplayAndWait();
            }
          }
        }
        //prgUpdateDisplayAndWait();
        
        if (chDirFlg)
        {
          w[0] = - w[0];
          w[1] = - w[1];
          w[4] = - w[4] - .002;

          chDirFlg = false;

          goto shoot;
        }
        j0 = j1;
        j1 = j2;
        x0 = x1;
        x1 = x2;
      }
    }
  }
  else
  {
    // loop over surfaces

    for (i=0; i<surface.n; i++)
    {
      surfPtr = &(surface[i]);

      // generate list of all points 

      allPnt.free();

      for (iLoop=0; iLoop<surfPtr->loop.n-1; iLoop++)
      {
        for (k=surfPtr->loop[iLoop]; k<surfPtr->loop[iLoop+1]; k++)
        {
          suppPntPtr = &(((GeomSpline*)(surfPtr->suppSpln[k]))->suppPnt);

          for (j=0; j<suppPntPtr->n; j++) allPnt.append((*suppPntPtr)[j]);
        }
      }

      for (iLoop=0; iLoop<surfPtr->loop.n-1; iLoop++)
      {
        // generate sequence of geometry points

        seqPnt.free();

        for (k=surfPtr->loop[iLoop]; k<surfPtr->loop[iLoop+1]; k++)
        {
          suppPntPtr = &(((GeomSpline*)(surfPtr->suppSpln[k]))->suppPnt);

          if (surfPtr->orientation[k] == FORWARD)
               for (j=0; j<suppPntPtr->n-1; j++) seqPnt.append((*suppPntPtr)[j]);
          else for (j=suppPntPtr->n-1; j>0; j--) seqPnt.append((*suppPntPtr)[j]);
        }
        seqPnt.append(seqPnt[0]);
        seqPnt.append(seqPnt[1]);
        seqPnt.append(seqPnt[2]);

        // test element size for each point

        x0 = ((GeomPoint*)(seqPnt[1]))->x;
        h0 = getElemSizeOptInternal(x0);

        for (j=2; j<seqPnt.n-1; j++)
        { 
          x1 = ((GeomPoint*)(seqPnt[j]))->x;
          h1 = getElemSizeOptInternal(x1);

          w[0] = x0[1] - x1[1];
          w[1] = x1[0] - x0[0];
          w[2] = sqrt(w[0]*w[0]+w[1]*w[1]);
          w[3] = 1. / w[2];
          w[0] *= w[3];
          w[1] *= w[3];

          x0d[0] = x0[0] + w[0] * h0 * fact;
          x0d[1] = x0[1] + w[1] * h0 * fact;

          x1d[0] = x1[0] + w[0] * h1 * fact;
          x1d[1] = x1[1] + w[1] * h1 * fact;

          for (k=0; k<allPnt.n; k++)
          {
            if (allPnt[k] != seqPnt[j-1] && allPnt[k] != seqPnt[j]
             && allPnt[k] != seqPnt[j-2] && allPnt[k] != seqPnt[j+1])
            {
              xk = ((GeomPoint*)(allPnt[k]))->x;
              if (pointInQuadrilateral2D(x0,x1,x1d,x0d,xk))
              {
                d =   ((x1[0]-x0[0])*(xk[1]-x0[1])-(x1[1]-x0[1])*(xk[0]-x0[0])) 
                    / ((x1[0]-x0[0])*w[1]-(x1[1]-x0[1])*w[0]);
               
                xd[0] = xk[0] - w[0] * d;
                xd[1] = xk[1] - w[1] * d; 

               /* plot.wipe();
                plot.setColour(0);
                plotGeometry(5);
                plot.setColour(3);
                plot.line(x1,x1d);
                plot.line(x1d,x0d);
                plot.line(x0d,x0);
                plot.line(x0,x1);
                plot.setColour(4);
                plot.point(xk);
                plot.setColour(5);
                plot.point(xd);
                prgUpdateDisplayAndWait(); */

                h = d * fact1;

                if (h < hmin) h = hmin;

                //if (getElemSizeOptInternal(xd) > h) hPntSrc.add(new HPointSource(xd,h,maxGrad));
              }
            /* else
              {
                plot.wipe();
                plot.setColour(0);
                plotGeometry(5);
                plot.setColour(3);
                plot.line(x1,x1d);
                plot.line(x1d,x0d);
                plot.line(x0d,x0);
                plot.line(x0,x1);
                plot.setColour(4);
                plot.point(xk);
                prgUpdateDisplayAndWait();
              } */
            }
          }
          h0 = h1;
          x0 = x1;
        }
      }
    }
  }
  return;
}









void Mesh::smoothElemSizeDist_hlp3(double maxGrad)
{
  int i, j, k;

  double h0, fact;

  bool changed = true;

  // max gradient

  while (changed)
  {
    changed = false;

    for (i=0; i<numnp; i++)
    {
      h0 = elemSizeOpt.x[i];

      for (j=0; j<nodeNode[i].n; j++)
      {
        k = nodeNode[i][j] - 1;

        fact = elemSizeOpt.x[k] + maxGrad * sqrt(dist2(&(x.x[i+i]),&(x.x[k+k]),2));

        if (h0 > fact * 1.0001) { h0 = fact; elemSizeOpt.x[i] = h0; changed = true; }
      }
    }
  }
  return;
}











void Mesh::smoothElemSizeDist_hlp4(double maxGrad)
{
  if (!geometryDiscretisedAtLeastOnce) return;

  int i, j, k;

  double h0, fact;

  bool changed = true, changed2 = false;

  Vector<int> delGeomObj;

  Vector<void*> keepGeomObj;

  // generate list of top geometry objects to be kept

  for (i=0; i<elemGrpToBeMeshed.n; i++)
    for (j=0; j<elemGrp[elemGrpToBeMeshed[i]].geomObj.n; j++)
      delGeomObj.append(elemGrp[elemGrpToBeMeshed[i]].geomObj[j]);

  for (i=0; i<surface.n; i++)
  {
    j = 0; while (j < delGeomObj.n && delGeomObj[j] != i) j++;
    if (j >= delGeomObj.n) keepGeomObj.append((void*)(&(surface[i])));
  }

  if (keepGeomObj.n == 0) return;

  // set elemSizeOpt of nodes to be kept to elemSizeCurr

  if (elemSizeCurr.n != numnp) prgError(1,"Mesh::smoothElemSizeDist_hlp4","fatal error!");

  for (i=0; i<numnp; i++)
    if (ndGmLnk[i].geomObj.containsAtLeastOneOf(keepGeomObj)) 
      { 
        elemSizeOpt[i] = elemSizeCurr[i];
        //cout << elemSizeOpt[i] << "\n";
        //plot.setColour(9);
        //plot.point(&(x.x[i+i]));
        //prgUpdateDisplayAndWait();
      }

  // max gradient

  while (changed)
  {
    changed = false;

    for (i=0; i<numnp; i++)
    {
      h0 = elemSizeOpt[i];

      for (j=0; j<nodeNode[i].n; j++)
      {
        k = nodeNode[i][j] - 1;

        fact = elemSizeOpt[k] - maxGrad * sqrt(dist2(&(x.x[i+i]),&(x.x[k+k]),2));

        if (h0 < fact * 0.9999) { h0 = fact; elemSizeOpt[i] = h0; changed = true; }
      }
    }
    if (changed) changed2 = true;
  }

  if (changed2) COUT << "Warning! Some optimal element size values had to be increased!\n\n";

  return;
}











void Mesh::reorganiseMesh(ListInfinite< MatrixFullArray<double> >  &xList,
                          ListInfinite< MatrixFullArray<int> >     &ixList,
                          ListInfinite< ListInfinite<NodeGeomLink> > &ndGmLnkList,
                          Vector<int> &delSurf)
{
  char fct[] = "Mesh::reorganiseMesh";

  if (ndm != 2) prgError(1,fct,"so far this is 2D only!");

  int e, i, ii, j, k, l, m, n, numnpKeep, numnpNewAll, numelKeep, numelNew;

  ElementGroup *eGPtr;

  Element *elPtr, **elemNew;

  GeomPoint *pntPtr;

  NodeGeomLink *ndGmLnkPtr;

  Vector<void*> keepTopGeomObj;

  Vector<int> delN, delP;

  MatrixFullArray<double> xNew;

  List<Vector<int> > transferNode;

  List<Vector<double> > transferFact;

  // generate list of top geometry objects to be kept

  for (i=0; i<surface.n; i++) 
  {
    j = 0; while (j < delSurf.n && delSurf[j] != i) j++;
    if (j == delSurf.n) keepTopGeomObj.append((void*)&(surface[i]));
  }

  // delete all nodes not associated with keepTopGeomObj

  delN.append(0);
  j = 1;
  i = 0;
  while (i < ndGmLnk.n) 
  {
    if (!keepTopGeomObj.containsAtLeastOneOf(ndGmLnk[i].geomObj))
    {
      delN.append(1);
      ndGmLnk.del(i);
      while (i < ndGmLnk.n && !keepTopGeomObj.containsAtLeastOneOf(ndGmLnk[i].geomObj))
      {
        delN[j]++;
        ndGmLnk.del(i); 
      }
      delP.append(i);
      j++;
    }
    else i++;
  }
  for (i=1; i<delN.n; i++) 
  { 
    delN[i]   += delN[i-1];
    delP[i-1] += delN[i];
  }
  delP.append(numnp + 10);

  // update geometry point dat1

  for (i=0; i<point.n; i++) point[i].dat1 = -1;

  numnpNewAll = ndGmLnk.n;

  for (i=0; i<numnpNewAll; i++)
  {
    pntPtr = ndGmLnk[i].whichPoint();
    if (pntPtr != NULL) pntPtr->dat1 = i;
  }

  // modify ndGmLnkList (new nodes added, old nodes identified, nodes deleted as required)

  for (i=0; i<ndGmLnkList.n; i++)
  {
    for (j=0; j<ndGmLnkList[i].n; j++)
    {
      ndGmLnkPtr = &(ndGmLnkList[i][j]);

      pntPtr = ndGmLnkPtr->whichPoint();

      if (pntPtr == NULL) ndGmLnkPtr->id = ++numnpNewAll;  // internal node (new node)

      else if (pntPtr->dat1 > -1) ndGmLnkPtr->id = pntPtr->dat1 + 1; // old boundary node

      else { pntPtr->dat1 = numnpNewAll; ndGmLnkPtr->id = ++numnpNewAll; } // new boundary node
    }
  }

  // retrieve old coordinates to be kept

  xNew.setDim(numnpNewAll,ndm,true);

  if (ndGmLnk.n > 0)
  {
    numnpKeep = numnp - delN.lastCoeff();

    for (i=0; i<numnpKeep; i++)
    {
      ii            = ndGmLnk[i].id - 1;
      xNew.x[i+i]   = x.x[ii+ii  ];
      xNew.x[i+i+1] = x.x[ii+ii+1];
      ndGmLnk[i].id = i + 1;
    }
  }
  else numnpKeep = 0;

  // add new nodes and get new coordinates from xList

  m = numnpKeep;

  for (i=0; i<ndGmLnkList.n; i++)
  {
    ii = 0;
    for (j=0; j<ndGmLnkList[i].n; j++)
    {
      ndGmLnkPtr = &(ndGmLnkList[i][j]);

      if (ndGmLnkPtr->id > m)
      {
        xNew.x[m+m]   = xList[i].x[ii  ];
        xNew.x[m+m+1] = xList[i].x[ii+1];
        ndGmLnk.add(new NodeGeomLink(*ndGmLnkPtr));
        m++;
      }
      else
      {
        for (k=0; k<ndGmLnkPtr->geomObj.n; k++)
          if (!ndGmLnk[ndGmLnkPtr->id-1].geomObj.contains(ndGmLnkPtr->geomObj[k]))
            ndGmLnk[ndGmLnkPtr->id-1].geomObj.append(ndGmLnkPtr->geomObj[k]);
      }
      ii += 2;
    }
  }
  if (m != numnpNewAll) prgError(1,fct,"fatal error!");

  // reorganise connectivity of elements to be kept

  numelKeep = 0;

  for (i=0; i<elemGrp.n; i++)
  {
    if (!elemGrpToBeMeshed.contains(i)) 
    {
      eGPtr = &(elemGrp[i]);
      numelKeep += eGPtr->elem.n;
      for (e=0; e<eGPtr->elem.n; e++)
      {
        for (j=0; j<eGPtr->elem[e].nen(); j++)
        {
          k = 0; while (eGPtr->elem[e].ix[j] > delP[k]) k++;
          eGPtr->elem[e].ix[j] -= delN[k];
        }
      }
    }
  }
  
  // get new elements

  numelNew = 0;

  for (i=0; i<ixList.n; i++) numelNew += ixList[i].nRow;

  elemNew = new Element* [numelNew];

  n = 0;

  for (i=0; i<elemGrpToBeMeshed.n; i++)
  {
    eGPtr = &(elemGrp[elemGrpToBeMeshed[i]]);

    for (j=0; j<eGPtr->geomObj.n; j++)
    {
      if (!delSurf.contains(eGPtr->geomObj[j],&m)) prgError(2,fct,"fatal error!");

      for (e=0; e<ixList[m].nRow; e++)
      {
        elPtr = newElement(eGPtr->elemProp[ELEMENTTYPE]->id);

        elPtr->elemGrp = (void*)eGPtr;

        elemNew[n++] = elPtr;

        for (l=0; l<elPtr->nen(); l++) 
        {
          k = ixList[m](e+1,l+1);

          elPtr->ix[l] = ndGmLnkList[m][k-1].id;
        }
      }
    }
  }
  if (n != numelNew) prgError(3,fct,"fatal error!");

  // get internal variables for new elements

  transferInternalVariables(elemNew,numelNew);

  // get node numbers and weighting factors for data transfer to new nodes

  for (i=numnpKeep; i<numnpNewAll; i++)
  {
    ii = i - numnpKeep;
    transferNode.add(new Vector<int>);
    transferFact.add(new Vector<double>);
    pntPtr = ndGmLnk[i].whichPoint();
    if (geometryDiscretisedAtLeastOnce && pntPtr != NULL)
    {
      if (pntPtr->dat1 != i) prgError(4,fct,"fatal error!");
      for (j=0; j<pntPtr->transferNode.n; j++)
      {
        transferNode[ii].append(pntPtr->transferNode[j]);
        transferFact[ii].append(pntPtr->transferFact[j]);
      }
    }
    else
      getTransferDataInternal(&(xNew.x[i+i]),transferNode[ii],transferFact[ii]);
  }

  // generate new nodal data arrays

  transfer_hlp( u,transferNode,transferFact,numnpKeep,numnpNewAll,ndf);
  transfer_hlp(un,transferNode,transferFact,numnpKeep,numnpNewAll,ndf);
  transfer_hlp(u3,transferNode,transferFact,numnpKeep,numnpNewAll,ndf);
  transfer_hlp(u4,transferNode,transferFact,numnpKeep,numnpNewAll,ndf);
  transfer_hlp(u5,transferNode,transferFact,numnpKeep,numnpNewAll,ndf);
  transfer_hlp(u6,transferNode,transferFact,numnpKeep,numnpNewAll,ndf);

  transfer_hlp(x0,transferNode,transferFact,numnpKeep,numnpNewAll,ndm);
  transfer_hlp(xn,transferNode,transferFact,numnpKeep,numnpNewAll,ndm);

  x.takeOver(xNew);

  reac.setDim(numnpNewAll,ndf,true);

  if (ndm > ndf) i = ndm; else i = ndf;
  r.setDim(numnpNewAll*i);

  //if (numelNew + numelKeep > numnpNewAll) i *= (numelNew + numelKeep); else i *= numnpNewAll;
  //outp.setDim(i);

  outp.setDim(numnpNewAll*i); outp.zero();

  nodeFlag.setDim(numnpNewAll);

  if (isALE())
  {
    transfer_hlp( d,transferNode,transferFact,numnpKeep,numnpNewAll,ndm);
    transfer_hlp(d0,transferNode,transferFact,numnpKeep,numnpNewAll,ndm);
    transfer_hlp( v,transferNode,transferFact,numnpKeep,numnpNewAll,ndm);
    transfer_hlp(vn,transferNode,transferFact,numnpKeep,numnpNewAll,ndm);

    reacMesh.setDim(numnpNewAll,ndm,true);
  }

  domainSpecificNodalDataTransfer(numnpNewAll);

  lastSearchNode = 1;

  // delete elements associated with delTopGeomObj and finish

  for (i=0; i<elemGrpToBeMeshed.n; i++) elemGrp[elemGrpToBeMeshed[i]].elem.free();

  for (e=0; e<numelNew; e++) elemNew[e]->belongsTo(elemNew[e]->elemGrp,this);

  delete [] elemNew;

  delete [] elem;

  elem = new Element* [numelKeep + numelNew];

  numel = 0;
  for (i=0; i<elemGrp.n; i++)
    { for (e=0; e<elemGrp[i].elem.n; e++) elem[numel++] = &(elemGrp[i].elem[e]); }

  if (numel != numelKeep + numelNew) prgError(5,fct,"fatal error!");

  numnp = numnpNewAll;

  return;
}












void Mesh::getTransferDataInternal(double *xp, Vector<int> &node, Vector<double> &factor)
{
  int i, e = pointInElement(xp);

  double N[100];

  if (e < 1) prgError(1,"Mesh::getTransferDataInternal","no background mesh element found!");

  elem[--e]->containsPoint(xp,N);

  node.free();
  factor.free();

  for (i=0; i<elem[e]->nen(); i++) 
  {
    node.append(elem[e]->ix[i]);
    factor.append(N[i]);
  }
  //node.  trunc(elem[e]->nen());
  //factor.trunc(elem[e]->nen());

  return;
}











void Mesh::transfer_hlp(MatrixFullArray<double> &dat,
                        List< Vector<int> >       &transferNode,
                        List< Vector<double> >    &transferFact,
                        int numnpKeep, int numnpNewAll, int ndfm)
{
  int i, j, k, ii, kk, indfm = 0;

  MatrixFullArray<double> datNew;

  datNew.setDim(numnpNewAll,ndfm,true);

  for (i=0; i<numnpKeep; i++)
  { 
    ii = (ndGmLnk[i].id - 1) * ndfm;
    for (j=0; j<ndfm; j++) datNew.x[indfm+j] = dat.x[ii+j];
    indfm += ndfm;
  }
  for (i=numnpKeep; i<numnpNewAll; i++)
  { 
    ii = i - numnpKeep;
    for (j=0; j<ndfm; j++) datNew.x[indfm+j] = 0.;
    for (k=0; k<transferNode[ii].n; k++)
    {
      kk = (transferNode[ii][k] - 1) * ndfm;
      for (j=0; j<ndfm; j++) datNew.x[indfm+j] += dat.x[kk+j] * transferFact[ii][k];
    }
    indfm += ndfm;
  }
  dat.takeOver(datNew);

  return;
}











void Mesh::transferInternalVariables(Element **elemNew, int numelNew)
{
  COUT << "Mesh::transferInternalVariables: NOT YET IMPLEMENTED!\n\n";

  return;
}









void Mesh::calcElemNodeDataAlongGeomObj(Vector<int> &el, 
                                        List< Vector<int> > &elNd,
                                        List< Vector<int> > &iNd ,
                                        Vector<int> &nd,
                                        void *geomObj)
{
  int e, i, j, k;

  for (j=0; j<numnp; j++)
  {
    if (ndGmLnk[j].geomObj.contains(geomObj))
    {
      nd.append(j+1);
      for (i=0; i<nodeElem[j].n; i++)
      {
        e = nodeElem[j][i];
        if (!el.contains(e+1,&k))
        {
          el.append(e+1);
          elNd.add(new Vector<int>);
          elNd[elNd.n-1].append(j+1);
          iNd.add(new Vector<int>);
          iNd[iNd.n-1].append(nd.n-1);
        }
        else
        {
          elNd[k].append(j+1);
          iNd[k].append(nd.n-1);
        }
      }
    }
  }

  for (i=0; i<el.n; i++)
  {
    if (elNd[i].n < 2) 
    {
      el.del(i);
      elNd.del(i);
      iNd.del(i);
    }
  }

  //cout << nd << "\n";
  //for (i=0; i<el.n; i++)
  //  cout << el[i] << ":" << elNd[i] << iNd[i] << "\n";
  //cout << "\n";

  return;
}
