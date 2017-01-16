
#include <iostream>

#include "FiniteElementBVPWI.h"
#include "MathGeom.h"
#include "Fluid.h"
#include "SolverSparse.h"



using namespace std;




void FiniteElementBVPWI::findBndNodesForInterpolations(VectorArray<bool> &nodeIsBndNode,
                                                       bool)
{
  char fct[] = "FiniteElementBVPWI::findBndNodesForInterpolations";

  Vector<int> nextNode;

  Vector<double> fact;
  
  double xc[3], factc;

  int e = 0, i, j, k = 0, k0;

  while (e < intpTmp.n)
  {
    // calculate centroid of interpolation element

    for (j=0; j<ndm; j++) xc[j] = 0.;

    factc = 1. / (double) intpTmp[e].n;

    for (i=0; i<intpTmp[e].n; i++)
      for (j=0; j<ndm; j++) xc[j] += freeNode[intpTmp[e][i]-1].x.x[j];

    for (j=0; j<ndm; j++) xc[j] *= factc;

    // find one node in interpolation element

    while (1)
    {
      if (k == 0) k = bndNd[0];

      k0 = k;
      while (!nodeFlag[k-1]) 
      {
        k++; if (k > numnp) k = 1;
        if (k == k0)
        { 
          prgWarning(1,fct,"no boundary node found in interface interpolation element!");
          goto jump;
        }
      }

      k = guessNearestMarkedNode(k,xc);

      k = findOneNodeInIntpElem(k,intpTmp[e],fact);

      if (k > 0) break;
    }

    // find all other nodes in interpolation element

    while (k > 0) 
    {
      if (!nodeIsBndNode[k-1])
      {
        connectNodeToIntf(k,intpTmp[e],fact);
        nodeIsBndNode[k-1] = true;
      }

      k = findAnotherNodeInIntpElem(k,nextNode,intpTmp[e],fact);
    }
    k = nodeFlagChanged[nodeFlagChanged.n-1] + 1;

    jump:

    if (show) cout << "  nodeFlagChanged.n = " << nodeFlagChanged.n << "\n";

    for (i=0; i<nodeFlagChanged.n; i++) nodeFlag[nodeFlagChanged[i]] = true;
    nodeFlagChanged.free();

    e++;
  }
  if (show) cout << "\n";

  return;
}











int FiniteElementBVPWI::guessNearestMarkedNode(int nd, double *xp)
{
  char fct[] = "FiniteElementBVPWI::guessNearestMarkedNode";

  if (!nodeFlag[nd-1]) prgError(1,fct,"start node not marked!");

  int i, k, imn0 = nd - 1, imn = imn0;

  double d = dist2(x.x+imn0*ndm,xp,ndm), dtrial;

  while (1)
  {
    for (i=0; i<nodeNode[imn0].n; i++)
    {
      k = nodeNode[imn0][i] - 1;
      if (nodeFlag[k])
      {
        dtrial = dist2(xp,x.x+k*ndm,ndm);
        if (dtrial < d - 1.e-12) { imn = k; d = dtrial; }
      }
    }
    if (imn0 == imn) return imn + 1;

    imn0 = imn;
  }
  prgError(2,fct,"fatal error!");
}









int FiniteElementBVPWI::findOneNodeInIntpElem(int k, Vector<int> &intpElem, Vector<double> &fact)
{
  int i, j, n, n0;

  Vector<int> nextNode;

  nextNode.append(k);

  while (nextNode.n > 0)
  {
    n0 = nextNode[0];

    nextNode.del(0);

    if (nodeInIntpElem(n0,intpElem,fact)) return n0;

    nodeFlag[n0-1] = false;

    nodeFlagChanged.append(n0-1);

    for (i=0; i<nodeNode[n0-1].n; i++)
    {
      n = nodeNode[n0-1][i];

      if (nodeFlag[n-1])
      {
        j = 0; while (j < nextNode.n) if (nextNode[j] != n) j++; else break;

        if (j == nextNode.n) nextNode.append(n);
      }
    }
  }
  return 0;
}






void FiniteElementBVPWI::connectNodeToIntf(int k, Vector<int> &intpElem, Vector<double> &fact)
{
  int i;

  Vector<void*> vTmp, vTmp2;

  for (i=0; i<intpElem.n; i++) vTmp.append((void*)(&freeNode[intpElem[i]-1]));

  if (freeNode[intpElem[0]-1].dofX->n == 0)

    bndNodeTmp.add(new BndNode(k,x.x+(k-1)*ndm,vTmp,vTmp2,fact));

  else bndNodeTmp.add(new BndNode(k,x.x+(k-1)*ndm,vTmp,vTmp,fact));

  return;
}






int FiniteElementBVPWI::findAnotherNodeInIntpElem(int n0, 
                                                  Vector<int> &nextNode,
                                                  Vector<int> &intpElem,
                                                  Vector<double> &fact)
{
  int i, j, n;

  nodeFlag[n0-1] = false;

  nodeFlagChanged.append(n0-1);

  for (i=0; i<nodeNode[n0-1].n; i++)
  {
    n = nodeNode[n0-1][i];

    if (nodeFlag[n-1])
    {
      j = 0; while (j < nextNode.n) if (nextNode[j] != n) j++; else break;

      if (j == nextNode.n) nextNode.append(n);
    }
  }

  n = nextNode[0];

  while (!nodeInIntpElem(n,intpElem,fact))
  {
    nodeFlag[n-1] = false;

    nodeFlagChanged.append(n-1);

    nextNode.del(0);
    if (nextNode.n == 0) break; else n = nextNode[0];
  }

  if (nextNode.n > 0) { nextNode.del(0); return n; } else return 0;
}









bool FiniteElementBVPWI::nodeInIntpElem(int nd, Vector<int> &intpElem, Vector<double> &fact)
{
  if (intpElem.n != ndm)
    prgError(2,"FiniteElementBVPWI::inInterpolElem","so far, only linear interpolation!");

  if (fact.n != ndm) { fact.free(); for (int i=0; i<ndm; i++) fact.append(); }

  double a[3], alpha, d, l, t = 1.e-4, area,
         *xp = x.x+(nd-1)*ndm,
         *x1 = freeNode[intpElem[0]-1].x.x, 
         *x2 = freeNode[intpElem[1]-1].x.x, *x3;

  if (ndm == 2) 
  {
    if (pointOnEdge2D(x1,x2,xp,&alpha,&d,&l)) { fact[0] = 1.-alpha; fact[1] = alpha; }

    else return false;
  }
  else 
  {
    x3 = freeNode[intpElem[2]-1].x.x;

    if (pointInTriangle3D(x1,x2,x3,xp,a,&d,&area)) { for (int i=0; i<3; i++) fact[i] = a[i]; }

    else return false;
  }
  return true;
}




