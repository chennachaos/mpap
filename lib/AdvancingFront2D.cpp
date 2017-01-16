
#include "AdvancingFront2D.h"
#include "MathGeom.h"
#include "FunctionsProgram.h"
#include "FunctionsElement.h"
#include "GeomSpline.h"
#include "Plot.h"
#include "ComputerTime.h"
#include "NodeGeomLink.h"


extern Plot plot;
extern ComputerTime computerTime;


typedef Poly2D Triangle2D;


using namespace std;





static void *notYetMeshed, *surfPtr;

static double hfact;

static int nd;

static ListInfinite<NodeGeomLink> *ndGmLnkPtr;

static Domain *dom;

static bool show;

//static int thisOne;




//==================================================================================================



Node2D::Node2D(void)
{ 
  bndFlg = false; 

  nOptCnct = 6;

  np = 0;

  return; 
}



Node2D::~Node2D()
{
  return;
}




void Node2D::draw(void) 
{ 
  double d = (plot.dAct[0] + plot.dAct[1]) * .0055; 
  plot.point(x,d); 
  return; 
}



//==================================================================================================



Edge2D::Edge2D(void)
{
  p0 = notYetMeshed;
  p1 = notYetMeshed;

  return;
}






Edge2D::Edge2D(Node2D *n0, Node2D *n1, double l, double *x1, double *x2)
{
  c0 = n0; 
  c1 = n1;

  p0 = notYetMeshed;
  p1 = notYetMeshed;

  length = l;

  normal[0] = x1[1] - x2[1];
  normal[1] = x2[0] - x1[0];

  return;
}







Edge2D::~Edge2D()
{
  return;
}






void Edge2D::draw(void) 
{ 
  plot.line(c0->x,c1->x); 

  return; 
}





//==================================================================================================



Poly2D::Poly2D(Node2D *n0, Node2D *n1, Node2D *n2, Edge2D *e0, Edge2D *e1, Edge2D *e2)
{ 
  // generate triangle

  // poly->node
  node.append((void*)n0);
  node.append((void*)n1);
  node.append((void*)n2);

  // poly->edge
  edge.append((void*)e0);
  edge.append((void*)e1);
  edge.append((void*)e2);

  // edge->poly
  if      (e0->p1 == notYetMeshed) e0->p1 = (void*) this;
  else if (e0->p0 == notYetMeshed) e0->p0 = (void*) this;
  else prgError(1,"Poly2D::Poly2D","fatal error!");

  if      (e1->p1 == notYetMeshed) e1->p1 = (void*) this;
  else if (e1->p0 == notYetMeshed) e1->p0 = (void*) this;
  else prgError(2,"Poly2D::Poly2D","fatal error!");

  if      (e2->p1 == notYetMeshed) e2->p1 = (void*) this;
  else if (e2->p0 == notYetMeshed) e2->p0 = (void*) this;
  else prgError(3,"Poly2D::Poly2D","fatal error!");

  // node->poly
  n0->poly.append((void*)this);
  n1->poly.append((void*)this);
  n2->poly.append((void*)this);

  return; 
}








Poly2D::Poly2D(Vector<void*> &patch, ListInfinite<Node2D> &nodeList,
                                     ListInfinite<Edge2D> &edgeList,
                                     ListInfinite<Poly2D> &polyList)
{
  // generate new polygon by joining those in patch,
  // delete polygons in patch,
  // update all connectivities

  Vector<void*> *tmp, *tmp2;

  Poly2D *p;

  Node2D *c, *c0, *c1;

  Edge2D *e;

  int i, j, k, l, m;

  void *p0, *p1;

  bool localDebug = false;

  node.free();
  edge.free();

  // collect all edges

  if (localDebug) cout << "A\n";

  for (i=0; i<patch.n; i++)
  {
    tmp = &(((Poly2D*)(patch[i]))->edge);

    for (j=0; j<tmp->n; j++) edge.append((*tmp)[j]);
  }

  // remove internal (double) edges from this polygon and from edgeList,
  // remove internal nodes from nodeList and
  // update node->node connectivities

  if (localDebug) cout << "B\n";

  i = 0;
  while (i < edge.n - 1)
  {
    j = i + 1;
    while (j < edge.n && edge[j] != edge[i]) j++;
    if (j < edge.n)
    {
      edge.del(j);
      while (j < edge.n && edge[j] != edge[i]) j++;
      if (j < edge.n) prgError(1,"Poly2D::Poly2D","fatal error!");

      e  = (Edge2D*)(edge[i]);

      c0 = e->c0;
      c1 = e->c1;

      k = 0;
      while (k < c0->node.n && c0->node[k] != (void*)c1) k++;
      if (k < c0->node.n) c0->node.del(k);
      else prgError(2,"Poly2D::Poly2D","fatal error!");

      k = 0;
      while (k < c1->node.n && c1->node[k] != (void*)c0) k++;
      if (k < c1->node.n) c1->node.del(k);
      else prgError(3,"Poly2D::Poly2D","fatal error!");

      if (c0->node.n == 0) nodeList.del(c0);
      if (c1->node.n == 0) nodeList.del(c1);

      edgeList.del(e);
      edge.del(i);
    }
    else i++;
  }

  // order edges and get nodes

  if (localDebug) cout << "C\n";

  node.append((void*)(((Edge2D*)(edge[0]))->c0));
  i = 0;
  c = ((Edge2D*)(edge[i]))->c1;
  while (i < edge.n - 1)
  {
    j = i + 1;
    while (j < edge.n)
    {
      e = (Edge2D*)(edge[j]);
      if (e->c0 == c || e->c1 == c) break; else j++;
    }
    if (j < edge.n)
    {
      node.append((void*) c); i++;
      if (e->c0 == c) c = e->c1; else c = e->c0; 
      edge.move(j,i);
    }
    else prgError(4,"Poly2D::Poly2D","fatal error!");
  }
  if (c != (Node2D*)(node[0])) prgError(5,"Poly2D::Poly2D","fatal error!");

  if (edge.n != node.n) prgError(6,"Poly2D::Poly2D","fatal error!");

  // ensure anticlockwise orientation

  if (localDebug) cout << "D\n";

  double xp[2] = {0.,0.}, *x0,
           *x1 = ((Node2D*)(node[0]))->x,
          area = triangleArea2D(xp,((Node2D*)(node[node.n-1]))->x,x1);

  for (i=1; i<node.n; i++)
  {
    x0 = x1;
    x1 = ((Node2D*)(node[i]))->x;
    area += triangleArea2D(xp,x0,x1);
  }
  if (area < 0.)
  {
    node.reverseOrder();
    node.move(node.n-1,0);
    edge.reverseOrder();
  }

  for (i=0; i<edge.n-1; i++)
    if (!((Edge2D*)(edge[i]))->is((Node2D*)(node[i]),(Node2D*)(node[i+1])))
      prgError(7,"Poly2D::Poly2D","fatal error!");
  if (!((Edge2D*)(edge[i]))->is((Node2D*)(node[i]),(Node2D*)(node[0])))
    prgError(8,"Poly2D::Poly2D","fatal error!");

  // reorder nodes and edges if polygon is a triangle

  if (localDebug) cout << "E\n";

  arrangeIfTriangle();

  // update node->poly connectivity

  if (localDebug) cout << "F\n";

  for (i=0; i<node.n; i++)
  {
    tmp = &(((Node2D*)(node[i]))->poly);
    for (j=0; j<patch.n; j++)
    {
      k = 0;
      while (k < tmp->n && patch[j] != (*tmp)[k]) k++;
      if (k < tmp->n) tmp->del(k);
    }
    tmp->append((void*)this);
  }

  // update edge->poly connectivity

  if (localDebug) cout << "G\n";

  for (i=0; i<edge.n; i++)
  {
    p0 = ((Edge2D*)(edge[i]))->p0;
    p1 = ((Edge2D*)(edge[i]))->p1;

    j = 0; 
    while (j < patch.n && p0 != patch[j] && p1 != patch[j]) j++;
    if (j < patch.n)
    {
      if (p0 == patch[j]) ((Edge2D*)(edge[i]))->p0 = (void*)this;
      else                ((Edge2D*)(edge[i]))->p1 = (void*)this;
    }
  }

  // delete polygons in polyList

  if (localDebug) cout << "H\n";

  for (i=0; i<patch.n; i++) polyList.del((Poly2D*)(patch[i]));

  return;
}









Poly2D::~Poly2D()
{
  return;
}









void Poly2D::triangleInfo(Edge2D *e0, Edge2D **e1, Edge2D **e2, Node2D **c2)
{
  if (edge.n != 3) prgError(1,"Poly2D::triangleInfo","fatal error!");

  if      (edge[0] == (void*)e0) 
      { *e1 = (Edge2D*)(edge[1]); *e2 = (Edge2D*)(edge[2]); *c2 = (Node2D*)(node[0]); }

  else if (edge[1] == (void*)e0) 
      { *e1 = (Edge2D*)(edge[2]); *e2 = (Edge2D*)(edge[0]); *c2 = (Node2D*)(node[1]); }

  else if (edge[2] == (void*)e0) 
      { *e1 = (Edge2D*)(edge[0]); *e2 = (Edge2D*)(edge[1]); *c2 = (Node2D*)(node[2]); }

  else prgError(2,"Poly2D::triangleInfo","fatal error!");

  return;
}








bool Poly2D::checkConsistency(ListInfinite<Node2D> &nodeList,
                              ListInfinite<Edge2D> &edgeList,
                              ListInfinite<Poly2D> &polyList)
{
  return true;

  Node2D *j0, *j1;

  Edge2D *e;

  Poly2D *p;

  int i, j, k;

  // for each node of the polygon: check consistency of node->node connectivities with edge list
  //                               check consistency with node->poly connectivities

  for (i=0; i<node.n; i++)
  {
    j0 = (Node2D*)(node[i]);

    for (j=0; j<j0->node.n; j++)
    {
      j1 = (Node2D*)(j0->node[j]);

      k = 0;
      while (k < edgeList.n && !edgeList[k].is(j0,j1)) k++;
      if (k == edgeList.n) prgError(1,"Poly2D::checkConsistency","fatal error!");
    }

    j = 0;
    while (j < j0->poly.n && j0->poly[j] != (void*)this) j++;
    if (j == j0->poly.n) prgError(2,"Poly2D::checkConsistency","fatal error!");
  }

  // for each edge of the polygon check consistency of edge->node and edge->poly connectivities

  if (node.n == 3) node.move(0,2);

  j0 = (Node2D*)(node[0]);
  for (i=0; i<edge.n-1; i++)
  {
    e  = (Edge2D*)(edge[i]);
    j1 = (Node2D*)(node[i+1]);
    if (!e->is(j0,j1)) prgError(3,"Poly2D::checkConsistency","fatal error!");
    j0 = j1;
  }
  e = (Edge2D*)(edge[edge.n-1]);
  if (!e->is(j0,(Node2D*)(node[0]))) prgError(4,"Poly2D::checkConsistency","fatal error!");

  if (node.n == 3) node.move(2,0);

  return true;

  // check area  (at some stages, e.g. before first smoothing, negative areas may be allowed!)

  double xp[2] = {0.,0.}, *x0,
           *x1 = ((Node2D*)(node[0]))->x,
          area = triangleArea2D(xp,((Node2D*)(node[node.n-1]))->x,x1);

  for (i=1; i<node.n; i++)
  {
    x0 = x1;
    x1 = ((Node2D*)(node[i]))->x;
    area += triangleArea2D(xp,x0,x1);
  }
  if (area < 0.) prgError(5,"Poly2D::checkConsistency","fatal error!");

  return true;
}









Edge2D* Poly2D::splitInTwo(ListInfinite<Node2D> &nodeList, 
                           ListInfinite<Edge2D> &edgeList, 
                           ListInfinite<Poly2D> &polyList, Node2D *cIn, int n1)
{
  Node2D *c0, *c1, *c;

  Edge2D *ei, *e;
 
  Poly2D *p0 = this, *p1;

  int i, j, i1, n;

  if (node.n < 4) prgError(1,"Poly2D::splitInTwo","fatal error!");

  if (cIn != NULL)
  { 
    j = 0;
    while (j < node.n && node[j] != (void*)cIn) j++;
    if (j == node.n) prgError(2,"Poly2D::splitInTwo","fatal error!");
  }
  else
  {
    cout << "search for node with largest internal angle!\nImplement it!!!!\n";
    exit(0);
  }
  c0 = (Node2D*)(node[j]);
  i  = node.n - 1;

  while (node[0] != (void*)c0) { node.move(0,i); edge.move(0,i); }

  if (cIn == NULL)
  {
    cout << "search for best node c1!\nImplement it!!!!\n";
    exit(0);
  }
  else 
  {
    i1 = n1;
    c1 = (Node2D*)(node[i1]);
  }

  // generate new polygon

  p1 = new Poly2D;

  polyList.add(p1);

  // generate new edge

  ei = new Edge2D;

  edgeList.add(ei);

  // update edge->node connectivity

  ei->c0 = c0;
  ei->c1 = c1;

  // update edge->poly connectivity

  ei->p0 = (void*)p1;
  ei->p1 = (void*)p0;

  for (i=i1; i<edge.n; i++)
  { 
    e = (Edge2D*)(edge[i]);
    if      (e->p0 == (void*)p0) e->p0 = (void*)p1;
    else if (e->p1 == (void*)p0) e->p1 = (void*)p1;
    else prgError(3,"Poly2D::splitInTwo","fatal error!");
  }

  // update node->node connectivity

  c0->node.append((void*)c1);
  c1->node.append((void*)c0);

  // update node->poly connectivity

  c0->poly.append((void*)p1);
  c1->poly.append((void*)p1);

  for (i=i1+1; i<node.n; i++)
  { 
    c = (Node2D*)(node[i]);
    j = 0;
    while (c->poly[j] != (void*)p0) j++;
    c->poly[j] = (void*)p1;
  }

  // update poly->edge connectivity

  p1->edge.append((void*)ei);

  n = edge.n;
  for (i=i1; i<n; i++) 
  {
    p1->edge.append(p0->edge[i1]);
    p0->edge.del(i1);
  }
  p0->edge.append((void*)ei);

  // update poly->node connectivity 

  p1->node.append((void*)c0);
  p1->node.append((void*)c1);

  n = node.n;
  for (i=i1+1; i<n; i++)
  {
    p1->node.append(p0->node[i1+1]);
    p0->node.del(i1+1);
  }

  // ensure special node / edge order for triangles

  p0->arrangeIfTriangle();
  p1->arrangeIfTriangle();

  p0->checkConsistency(nodeList,edgeList,polyList);
  p1->checkConsistency(nodeList,edgeList,polyList);

  return ei;
}







void Poly2D::arrangeIfTriangle(void)
{
  if (node.n == 3)
  {
    if (!((Edge2D*)(edge[0]))->is((Node2D*)(node[1]),(Node2D*)(node[2])))
    {
      edge.move(0,2);
      if (!((Edge2D*)(edge[0]))->is((Node2D*)(node[1]),(Node2D*)(node[2])))
      {
        edge.move(0,2);
        if (!((Edge2D*)(edge[1]))->is((Node2D*)(node[2]),(Node2D*)(node[0]))) 
        {  
          edge.move(1,2);
        }
      }
    }
    if (!((Edge2D*)(edge[0]))->is((Node2D*)(node[1]),(Node2D*)(node[2])))
      prgError(1,"Poly2D::arrangeIfTriangle","fatal error!");
    if (!((Edge2D*)(edge[1]))->is((Node2D*)(node[2]),(Node2D*)(node[0])))
      prgError(2,"Poly2D::arrangeIfTriangle","fatal error!");
    if (!((Edge2D*)(edge[2]))->is((Node2D*)(node[0]),(Node2D*)(node[1])))
      prgError(3,"Poly2D::arrangeIfTriangle","fatal error!");
  }
  return;
}









double Poly2D::maxInternalAngle(int *nn)
{
  double maxAngle = -1., v1[2], v2[2], *x0, *x1, *x2, fact;

  int i;

  node.append(node[0]);
  node.append(node[1]);

  x0 = ((Node2D*)(node[0]))->x;
  x1 = ((Node2D*)(node[1]))->x;

  for (i=2; i<node.n; i++)
  {
    x2 = ((Node2D*)(node[i]))->x;

    v1[0] = x1[0] - x0[0];
    v1[1] = x1[1] - x0[1];
    v2[0] = x2[0] - x1[0];
    v2[1] = x2[1] - x1[1];

    fact = boundaryAngleBetweenTwoEdges2D(v1,v2,2);

    if (maxAngle < fact) maxAngle = fact;

    x0 = x1;
    x1 = x2;
  }

  node.del(node.n-1);
  node.del(node.n-1);

  if (nn != NULL)
  {
    if (i >= node.n) i -= node.n;
    *nn = i;
  }

  return maxAngle;
}









double Poly2D::aspectRatio(void)
{
  if (node.n == 3) return triangleAspectRatio(((Node2D*)node[0])->x,
                                              ((Node2D*)node[1])->x,
                                              ((Node2D*)node[2])->x) * 2.;
  else if (node.n == 4)
    return (triangleAspectRatio(((Node2D*)node[0])->x,
                                ((Node2D*)node[1])->x,
                                ((Node2D*)node[2])->x) 
          + triangleAspectRatio(((Node2D*)node[1])->x,
                                ((Node2D*)node[2])->x,
                                ((Node2D*)node[3])->x) 
          + triangleAspectRatio(((Node2D*)node[2])->x,
                                ((Node2D*)node[3])->x,
                                ((Node2D*)node[0])->x) 
          + triangleAspectRatio(((Node2D*)node[3])->x,
                                ((Node2D*)node[0])->x,
                                ((Node2D*)node[1])->x)) * 0.60355339;
  else prgError(1,"Poly2D::aspectRatio","invalid node number!");

  return 0.;
}








void Poly2D::splitIntoQuads(ListInfinite<Node2D> &nodeList,
                            ListInfinite<Edge2D> &edgeList,
                            ListInfinite<Poly2D> &polyList)
{
  int i, j, k, l, me, mp, n = node.n, nd2 = roundToInt(.5*(double)n);

  Node2D *cNew;

  Edge2D *eNew, *e0, *e1, *e2, *e3;

  Poly2D *pNew;

  Vector<void*> *tmp;

  void *p;

  if (nd2 * 2 != n) prgError(1,"Poly2D::splitIntoQuads","fatal error!");
  if (n < 6)        prgError(2,"Poly2D::splitIntoQuads","fatal error!");

  // generate centre node

  cNew = new Node2D;

  nodeList.add(cNew);

  cNew->x[0] = 0.;
  cNew->x[1] = 0.;

  for (i=0; i<node.n; i++)
  {
    cNew->x[0] += ((Node2D*)(node[i]))->x[0];
    cNew->x[1] += ((Node2D*)(node[i]))->x[1];
  }
  cNew->x[0] /= (double)n;
  cNew->x[1] /= (double)n;

  // update node->node connectivity

  for (i=1; i<n; i+=2)
  {
    cNew->node.append(node[i]);
    ((Node2D*)(node[i]))->node.append((void*)cNew);
  }

  // generate edges and polygons and 
  // update edge->node, edge->poly, poly->node, poly->edge and node->poly connectivities

  me = edgeList.n;

  mp = polyList.n;

  eNew = new Edge2D;

  edgeList.add(eNew);

  eNew->c0 = (Node2D*)(node[1]);
  eNew->c1 = cNew;

  e0 = eNew;

  for (i=3; i<n; i+=2)
  {
    eNew = new Edge2D;

    edgeList.add(eNew);

    eNew->c0 = (Node2D*)(node[i]);
    eNew->c1 = cNew;

    pNew = new Poly2D;

    polyList.add(pNew);

    cNew->poly.append((void*)pNew);
    l = 0;
    tmp = &(((Node2D*)(node[i-2]))->poly);
    while ((*tmp)[l] != (void*)this) l++;
    (*tmp)[l] = (void*)pNew;
    l = 0;
    tmp = &(((Node2D*)(node[i-1]))->poly);
    while ((*tmp)[l] != (void*)this) l++;
    (*tmp)[l] = (void*)pNew;
    tmp = &(((Node2D*)(node[i]))->poly);
    tmp->append((void*)pNew);

    pNew->node.append((void*)cNew);
    pNew->node.append(node[i-2]);
    pNew->node.append(node[i-1]);
    pNew->node.append(node[i]);

    pNew->edge.append((void*)e0);
    pNew->edge.append(edge[i-2]);
    pNew->edge.append(edge[i-1]);
    pNew->edge.append((void*)eNew);

    e0 = eNew;
  }

  k = 0;

  e2 = &(edgeList[me]);
  e2->p0 = (void*)this;

  for (i=3; i<n; i+=2)
  {
    e0 = (Edge2D*)(edge[i-2]);
    e1 = (Edge2D*)(edge[i-1]);

    p  = (void*)(&(polyList[mp+k]));

    if (e0->p0 == (void*)this) e0->p0 = p;
    else                       e0->p1 = p;

    if (e1->p0 == (void*)this) e1->p0 = p;
    else                       e1->p1 = p;

    e3 = &(edgeList[me+ ++k]);

    e2->p1 = p;
    e3->p0 = p;

    e2 = e3;
  }
  e2->p1 = (void*)this;

  cNew                  ->poly.append((void*)this);
  ((Node2D*)(node[1]  ))->poly.append((void*)this);

  node[3] = node[1];
  node[2] = node[0];
  node[1] = node[n-1];
  node[0] = (void*)cNew;
  node.trunc(4);

  edge[3] = (void*)(&(edgeList[me]));
  edge[2] = edge[0];
  edge[1] = edge[n-1];
  edge[0] = (void*)eNew;
  edge.trunc(4);

  return;
}









void Poly2D::draw(bool colFlg)
{
  int i;

  plot.line(((Node2D*)(node[node.n-1]))->x,((Node2D*)(node[0]))->x);

  for (i=1; i<node.n; i++)
    plot.line(((Node2D*)(node[i]))->x,((Node2D*)(node[i-1]))->x);

  return;
}



//==================================================================================================







AdvancingFront2D::AdvancingFront2D(void)
{
  return;
}





AdvancingFront2D::~AdvancingFront2D()
{
  return;
}





#include "FunctionsEssGrp.h"

void AdvancingFront2D::generateMesh(MatrixFullArray<double> &x,
                                    MatrixFullArray<int> &ix,
                                    void *ndGmLnk,
                                    int &nBndNd,
                                    GeomSurface &surf, Domain *domain, int nen, bool showFlg)
{
  cout << "   --- AdvancingFront2D::generateMesh (begin) ---\n\n";

  if (nen < 3 || nen > 9 || nen == 5 || nen == 7) 
    prgError(1,"AdvancingFront2D::generateMesh","invalid nen!");

  // set static variables

  if (nen == 3 || nen == 6) hfact = 1.; else hfact = 2.;

  surfPtr = (void*)(&surf);

  ndGmLnkPtr = (ListInfinite<NodeGeomLink>*)ndGmLnk;

  dom = domain;

  nd = nen;

  Poly2D dmyPoly;
  notYetMeshed = (void*)(&dmyPoly);

  show = showFlg;

  // generate mesh

  //.......................................................

  computerTime.go("generateInitialFront, loop generateNextTriangle");

  generateInitialFront();

  if (showFlg) 
  {
    printf("          show mode: press <enter> repeatedly! ");
    prgUpdateDisplay(true);
    cout << "\n";
  }

  cout << "          number of triangles =       ";

  int    i, j;
  double xmn[3], xmx[3];      

  if (showFlg) 
  {
    printf("\b\b\b\b0");
    if (!plot) 
    {
      for (j=0; j<2; j++) { xmx[j] = node[0].x[j]; xmn[j] = node[0].x[j]; }
      for (i=1; i<node.n; i++) 
      { 
        for (j=0; j<2; j++)
        { 
          if (xmx[j] < node[i].x[j]) xmx[j] = node[i].x[j]; 
          if (xmn[j] > node[i].x[j]) xmn[j] = node[i].x[j]; 
        }
      }
      plot.fit(xmn,xmx,2);
    }
    plot.wipe();
    plot.setColour(3);
    draw();
    for (i=0; i<node.n; i++) if (node[i].np > 0) node[i].draw();

    //plot.setColour(2); node[thisOne].draw();

    prgUpdateDisplay(true);

    plot.wipe(); 
    plot.setColour(3); 
    while (generateNextTriangle() > 0)
    {
      plot.wipe(); 
      plot.setColour(3); 
      draw(); 
      drawPolygons();
      for (i=0; i<node.n; i++) if (node[i].np > 0) node[i].draw();

    //plot.setColour(2); node[thisOne].draw();

      prgUpdateDisplay(true); 
    }
  }
  else while (generateNextTriangle() > 0);

  cout << "\n\n";

  computerTime.stopAndPrint("generateInitialFront, loop generateNextTriangle");

  //.......................................................

  computerTime.go("remove3in1, remove8in6, remove8to1");

  remove3in1();

  remove8in6();

  remove8to1();

  computerTime.stopAndPrint("remove3in1, remove8in6, remove8to1");

  showAndWait();

  //.......................................................

  computerTime.go("loop swapEdges, remove4in2");

  while (swapEdges());

  remove4in2();

  computerTime.stopAndPrint("loop swapEdges, remove4in2");

  showAndWait();

  //.......................................................

  computerTime.go("smooth1, smooth2");

  smooth1(2);
  smooth2(2);

  computerTime.stopAndPrint("smooth1, smooth2");

  showAndWait();

  //.......................................................

  if (nd == 6)
  {
    computerTime.go("subdivideAllEdges");

    subdivideAllEdges();

    computerTime.stopAndPrint("subdivideAllEdges");

    showAndWait();
  }

  if (nd == 3 || nd == 6) goto finish;

  //.......................................................

  computerTime.go("convertToQuads"); 

  convertToQuads();

  computerTime.stopAndPrint("convertToQuads");
 
  showAndWait();

  //.......................................................

  computerTime.go("correct6to1, correct2to1, correct7to1");

  correct6to1();

  correct2to1();

  correct7to1();

  computerTime.stopAndPrint("correct6to1, correct2to1, correct7to1");

  showAndWait();

  //.......................................................

  computerTime.go("smooth1, smooth2");

  smooth1(2);
  smooth2(3);

  computerTime.stopAndPrint("smooth1, smooth2");

  showAndWait();

  if (nd == 4) goto finish;

  // ......................................................

  computerTime.go("subdivideAllEdges");

  subdivideAllEdges();

  computerTime.stopAndPrint("subdivideAllEdges");

  if (nd == 8) goto finish;

  // ......................................................

  computerTime.go("generateMidpoints");

  generateMidpoints();

  computerTime.stopAndPrint("generateMidpoints");

  // ......................................................

  finish:

  int nen05 = int(0.5*(double)nen);

  x.setDim(node.n,2,true);

  double *X = x.x;

  for (i=0; i<node.n; i++) 
  {
    X[i+i]   = node[i].x[0];
    X[i+i+1] = node[i].x[1];
    node[i].np = i+1;
  }

  ix.setDim(poly.n,nen+1,true);

  int *IX = ix.x;

  if (nen < 5)
  {
    for (i=0; i<poly.n; i++)
    {
      if (poly[i].node.n != nen) prgError(2,"AdvancingFront2D::generateMesh","fatal error!");
      IX[i*(nen+1)+nen] = 1;
      for (j=0; j<nen; j++) IX[i*(nen+1)+j] = ((Node2D*)(poly[i].node[j]))->np;
    }
  }
  else
  {
    for (i=0; i<poly.n; i++)
    {
      if (poly[i].node.n != nen) prgError(2,"AdvancingFront2D::generateMesh","fatal error!");
      IX[i*(nen+1)+nen] = 1;
      for (j=0; j<nen05; j++)
      {
        IX[i*(nen+1)+j]       = ((Node2D*)(poly[i].node[j+j]))->np;
        IX[i*(nen+1)+nen05+j] = ((Node2D*)(poly[i].node[j+j+1]))->np;
      }
      if (nen == 9) IX[i*(nen+1)+8] = ((Node2D*)(poly[i].node[8]))->np;
    }
  }

  for (i=nBndNode; i<node.n; i++)
  {
    (*ndGmLnkPtr)[i].isOn(surfPtr);
    (*ndGmLnkPtr)[i].id = i + 1;
  }

  nBndNd = nBndNode;
  
  if (nen == 3 || nen == 6)
    cout << "          " << poly.n << " triangles and " << node.n << " nodes generated.\n\n";
  else
    cout << "          " << poly.n << " quadrilaterals and " << node.n << " nodes generated.\n\n";
/*
  double ar, armn = 1., armx = 0.;

  for (i=0; i<poly.n; i++)
  { 
    ar = poly[i].aspectRatio(); 
    if (ar < armn) armn = ar;
    if (ar > armx) armx = ar;
  }

  cout << "          best element aspect ratio  = " << armx << "\n";
  cout << "          worst element aspect ratio = " << armn << "\n\n";
*/
  cout << "   ---- AdvancingFront2D::generateMesh (end) ----\n\n";

  return;
}







void AdvancingFront2D::showAndWait(void)
{
  if (!show) return;
  plot.wipe();
  plot.setColour(3);
  drawPolygons();
  prgUpdateDisplay(true);
  return;
}










void AdvancingFront2D::generateInitialFront(void)
{
  int i, j, m, m0, n = 0, n0, mult, iLoop;

  double *X, *x0, *x1, *x2, v1[2], v2[2], fact;

  Vector<void*> *suppPnt;

  GeomSurface *surf = (GeomSurface*)surfPtr;

  Vector<int> &loop = surf->loop;

  void *suppPntj, *splni;

  if      (nd == 3) mult = 1;
  else if (nd == 4) mult = 2;
  else if (nd == 6) mult = 2;
  else if (nd == 8) mult = 4;
  else if (nd == 9) mult = 4;
  else prgError(1,"AdvancingFront2D::generateInitialFront","fatal error!");

  for (iLoop=0; iLoop<loop.n-1; iLoop++)
  {
    // generate points (use geometry points)

    n0 = n;

    for (i=loop[iLoop]; i<loop[iLoop+1]; i++)
    {
      m0 = n;

      if (i > loop[iLoop]) (*ndGmLnkPtr)[n].isOn(splni);

      splni   = surf->suppSpln[i];
      suppPnt = &(((GeomSpline*)splni)->suppPnt);

      if (surf->orientation[i] == FORWARD) 
      {
        for (j=0; j<suppPnt->n-1; j++)
        {
          suppPntj = (*suppPnt)[j];
          X = node[n].x;
          node[n].bndFlg = true;
          (*ndGmLnkPtr)[n].isOn(splni,suppPntj,surfPtr);
          (*ndGmLnkPtr)[n].id = ++n;
          X[0] = ((GeomPoint*)suppPntj)->x[0];
          X[1] = ((GeomPoint*)suppPntj)->x[1];
          X[2] = hfact * ((GeomPoint*)suppPntj)->h;
        } 
      }
      else
      {
        for (j=suppPnt->n-1; j>0; j--)
        {
          suppPntj = (*suppPnt)[j];
          X = node[n].x;
          node[n].bndFlg = true;
          (*ndGmLnkPtr)[n].isOn(splni,suppPntj,surfPtr);
          (*ndGmLnkPtr)[n].id = ++n;
          X[0] = ((GeomPoint*)suppPntj)->x[0];
          X[1] = ((GeomPoint*)suppPntj)->x[1];
          X[2] = hfact * ((GeomPoint*)suppPntj)->h;
        }
      }
    }
    (*ndGmLnkPtr)[n0].isOn(splni);

    // generate edges

    m = frontEdge.n;

    frontEdge[m].c0 = &(node[node.n - mult]);

    frontEdge[m].c1 = &(node[n0]);

    for (j=1; j<mult; j++) frontEdge[m].midpoint.append((void*)(&(node[node.n - mult + j])));

    node[n0].np = 2;
    m++;

    for (i=n0+mult; i<node.n; i+=mult) 
    {
      frontEdge[m].c0 = &(node[i - mult]);
      frontEdge[m].c1 = &(node[i]);
      for (j=1; j<mult; j++) frontEdge[m].midpoint.append((void*)(&(node[i - mult + j])));
      node[i].np = 2;
      m++;
    }

  }

  // calculate edge data

  for (i=0; i<frontEdge.n; i++) frontEdge[i].c0->dx.setDim(4);

  for (i=0; i<frontEdge.n; i++)
  {
    x0 = frontEdge[i].c0->x;
    x1 = frontEdge[i].c1->x;
    frontEdge[i].length    = sqrt(dist2(x0,x1,2));
    frontEdge[i].normal[0] = x0[1] - x1[1];
    frontEdge[i].normal[1] = x1[0] - x0[0];
    frontEdge[i].c0->node.append((void*)frontEdge[i].c1);
    frontEdge[i].c1->node.append((void*)frontEdge[i].c0);
    frontEdge[i].p1 = NULL;

    frontEdge[i].c1->dx[3] = x1[1] - x0[1];
    frontEdge[i].c0->dx[0] = x1[0] - x0[0];
    frontEdge[i].c0->dx[1] = x1[1] - x0[1];
    frontEdge[i].c1->dx[2] = x1[0] - x0[0];
  }
  i = node.n - mult;
  suppPntj        = node[i].node[0];
  node[i].node[0] = node[i].node[1];
  node[i].node[1] = suppPntj;

  // calculate optimal connectivity for boundary nodes

  for (i=0; i<node.n; i+=mult)
  {
    x0 = ((Node2D*)(node[i].node[0]))->x;
    x1 =            node[i]           .x;
    x2 = ((Node2D*)(node[i].node[1]))->x;

    v1[0] = x1[0] - x0[0];
    v1[1] = x1[1] - x0[1];
    v2[0] = x2[0] - x1[0];
    v2[1] = x2[1] - x1[1];

    fact = boundaryAngleBetweenTwoEdges2D(v1,v2,2);

    if      (fact < 1.3) node[i].nOptCnct = 2;
    else if (fact < 2.618) node[i].nOptCnct = 3;
    else if (fact < 3.666) node[i].nOptCnct = 4;
    else if (fact < 4.713) node[i].nOptCnct = 5;
    else if (fact < 5.760) node[i].nOptCnct = 6;
    else                   node[i].nOptCnct = 7;
  }

  //for (i=0; i<frontEdge.n; i++)
  //  cout << frontEdge[i].c0->x[0] << "," << frontEdge[i].c0->x[1] << "\n";


  // sort edges shortest first

  frontEdge.quickSort();

  // prepare for generateTriangle

  nFail = 0;

  nBndNode = node.n;

  nEdge = frontEdge.n;

  return;
}











int AdvancingFront2D::generateNextTriangle(void)
{
  //cout << "generateNextTriangle\n";

  if (show) printf("                                   ");
  printf("\b\b\b\b\b\b\b%7d",poly.n+1);

  int i, n;

  Edge2D *e0 = &(frontEdge[0]), *e1, *e2;

  Node2D *c0 = e0->c0, *c1 = e0->c1, *c2;

  double *x0 = c0->x, *x1 = c1->x, *x2,
         *normal = e0->normal, length  = e0->length,
         xIdeal[2], v0[2], v1[2], fact, fact1;

  Vector<int>      nodeList;
  Vector<double>   nodeDist;
  VectorArray<int> perm;

  // generate ideal point

  fact = .5 * (x0[2] + x1[2]) / length;
  if      (fact <  .55) fact =  .55;
  else if (fact > 2.  ) fact = 2.  ;

  fact1 = fact * 0.71 * length; // 0.8

  fact *= 0.8660254;

  xIdeal[0] = .5 * (x0[0] + x1[0]) + normal[0] * fact;
  xIdeal[1] = .5 * (x0[1] + x1[1]) + normal[1] * fact;
 
  // generate ordered list of points on current front within critical distance

  for (i=0; i<node.n; i++)
  { 
    if (node[i].np > 0)
    {
      fact = sqrt(dist2(xIdeal,node[i].x,2));
      if (fact < fact1) 
      {
        //cout << "\n" << node[i].x[0] << "," << node[i].x[1] << "," << i << "\n";

        nodeList.append(i);
        nodeDist.append(fact);
      }
    }
  }

  //cout << nodeList << " in\n";

  if (nodeList.n > 0)
  {
    perm = nodeList;

    nodeDist.quickSort(false,perm.x);

    nodeList = perm;

    //cout << nodeList << " out\n";

    //bubbleSort(nodeList.n,nodeList,nodeDist);

    // test all nodes in the list 

    //thisOne = nodeList[0];

    i = 0;
    while (i < nodeList.n)
      if (!triangleOK(x0,x1,node[nodeList[i]].x,&node[nodeList[i]],c0,c1)) i++; else break;
    if (i < nodeList.n) 
    {
      c2 = &(node[nodeList[i]]);
      x2 = c2->x;
      goto useExistingPoint;
    }
    else goto searchForOtherExistingPoint;
  }
  else
  {
    // test ideal point 

    if (triangleOK(x0,x1,xIdeal))
    {
      //cout << node.n << "\n";
      c2 = &(node[node.n]);
      x2 = c2->x;
      x2[0] = xIdeal[0];
      x2[1] = xIdeal[1];
      x2[2] = hfact * dom->getElemSizeOptInternal(x2);
      goto useIdealPoint;
    }
    else goto searchForOtherExistingPoint;
  }
  
  // ............................................................

  searchForOtherExistingPoint:

  //cout << "searchForOtherExistingPoint\n";

  c2 = NULL;
  fact1 = 10;
  for (i=0; i<node.n; i++)
  { 
    if (node[i].np > 0)
    {
      x2 = node[i].x;
      if (x2 != x0 && x2!=x1)
      {
        if (triangleOK(x0,x1,x2,&node[i],c0,c1)) 
        {
          v0[0] = x0[0] - x2[0];
          v0[1] = x0[1] - x2[1];
          v1[0] = x1[0] - x2[0];
          v1[1] = x1[1] - x2[1];

          fact = cosAngleBetweenVectors(v0,v1,2);

          if (fact < fact1) { fact1 = fact; c2 = &(node[i]); }
        }
      }
    }
  }
  if (c2 != NULL) 
  { 
    x2 = c2->x; 
    goto useExistingPoint; 
  }

  // ............................................................

  // failure

  //cout << "  no element generated!\n";

  if (nFail + 1 == node.n) prgError(1,"AdvancingFront2D::generateElement","fatal failure!");

  frontEdge.move(0,nFail);    // undo previous failure swap
  frontEdge.move(++nFail,0);  // get first virgin edge to the front

  return frontEdge.n;

  // ............................................................

  useExistingPoint:

  //cout << "useExistingPoint\n";

  edge.add(frontEdge.takeOut(e0));

  fact = sqrt(dist2(x2,x1,2));
  i = 0;
  n = frontEdge.n;
  while (i < n && frontEdge[i].length < fact * 1.1 && !frontEdge[i].isStrictly(c1,c2)) i++;
  if (i == n || !frontEdge[i].isStrictly(c1,c2))
  {
    nEdge++;

    // increase front_node->front_node counters
    c2->np++;
    //c1->np++;

    // update node->node connectivity
    c2->node.append((void*)c1);
    c1->node.append((void*)c2);

    // generate new frontEdge
    e1 = new Edge2D(c2,c1,fact,x2,x1);
    frontEdge.add(e1); n++;

    // move it to right location in the list
    i = 0;
    while (i < n && frontEdge[i].length < fact) i++;
    if (i < n) frontEdge.move(n-1,i);
  }
  else if (i < n && frontEdge[i].isStrictly(c1,c2))
  { 
    // decrease front_node->front_node counters
    c1->np -= 2;
    c2->np--;

    // move edge from frontEdge to edge list
    e1 = &(frontEdge[i]);
    edge.add(frontEdge.takeOut(e1));
  }

  fact = sqrt(dist2(x0,x2,2));
  i = 0;
  n = frontEdge.n;
  while (i < n && frontEdge[i].length < fact * 1.1 && !frontEdge[i].isStrictly(c2,c0)) i++;
  if (i == n || !frontEdge[i].isStrictly(c2,c0))
  {
    nEdge++;

    //c0->np++;
    c2->np++;

    c0->node.append((void*)c2);
    c2->node.append((void*)c0);

    e2 = new Edge2D(c0,c2,fact,x0,x2);
    frontEdge.add(e2); n++;

    i = 0;
    while (i < n && frontEdge[i].length < fact) i++;
    if (i < n) frontEdge.move(n-1,i);
  }
  else if (i < n && frontEdge[i].isStrictly(c2,c0)) 
  { 
    c2->np--;
    c0->np -= 2;

    e2 = &(frontEdge[i]);
    edge.add(frontEdge.takeOut(e2));
  }

  // generate new triangle and update connectivities
  //   node->poly
  //   poly->node
  //   poly->edge
  //   edge->poly
  poly.add(new Triangle2D(c2,c0,c1,e0,e1,e2));

  nFail = 0;

  return frontEdge.n;

  // ............................................................

  useIdealPoint:

  //cout << "useIdealPoint\n";

  edge.add(frontEdge.takeOut(e0));

  fact = sqrt(dist2(x2,x0,2));

  n = frontEdge.n;

  nEdge += 2;
  n     += 2;

  e1 = new Edge2D(c2,c1,fact,x2,x1);
  frontEdge.add(e1);
  
  e2 = new Edge2D(c0,c2,fact,x0,x2);
  frontEdge.add(e2);

  // set front_node->front_node counter
  c2->np = 2;

  // update node->node connectivity
  c0->node.append((void*)c2);
  c2->node.append((void*)c0);
  c1->node.append((void*)c2);
  c2->node.append((void*)c1);

  // move edges to right location in the list
  i = 0;
  while (i < n && frontEdge[i].length < fact) i++;
  if (i < n) { frontEdge.move(n-1,i); frontEdge.move(n-1,i); }

  // generate new triangle and update connectivities
  //   node->poly
  //   poly->node
  //   poly->edge
  //   edge->poly
  poly.add(new Triangle2D(c2,c0,c1,e0,e1,e2));

  nFail = 0;

  return frontEdge.n;
}












bool AdvancingFront2D::triangleOK(double *x0, double *x1, double *x2, Node2D *ndPtr, Node2D *ndPtr0, Node2D *ndPtr1)
{
//  cout << "\n" << x0[0] << "," << x0[1] << " x0\n";
//  cout << x1[0] << "," << x1[1] << " x1\n";
//  cout << x2[0] << "," << x2[1] << " x2\n";

  if (triangleArea2D(x0,x1,x2) < 0) return false;

  int i;

  double *xx0, *xx1;

  // check for intersecting edges

  for (i=0; i<frontEdge.n; i++)
  {
    xx0 = frontEdge[i].c0->x;
    xx1 = frontEdge[i].c1->x;
    if (edgesIntersect2D(x0,x2,xx0,xx1)) return false;
    if (edgesIntersect2D(x1,x2,xx0,xx1)) return false;
  }

//  cout << "no intersecting lines\n";

  // check for active points inside the triangle

  double longestEdgeSquared = (x0[0]-x1[0])*(x0[0]-x1[0])+(x0[1]-x1[1])*(x0[1]-x1[1]),
         fact               = (x1[0]-x2[0])*(x1[0]-x2[0])+(x1[1]-x2[1])*(x1[1]-x2[1]);

  if (fact > longestEdgeSquared) longestEdgeSquared = fact;

  fact = (x2[0]-x0[0])*(x2[0]-x0[0])+(x2[1]-x0[1])*(x2[1]-x0[1]);
  if (fact > longestEdgeSquared) longestEdgeSquared = fact;

  longestEdgeSquared *= 1.1;

  int c = 0;

  for (i=0; i<node.n; i++)
  {
    if (node[i].np > 0)
    {
      if (dist2(x0,node[i].x,2) < longestEdgeSquared)
      {
        if (pointInTriangle2D(x0,x1,x2,node[i].x)) { c++; if (c == 4) return false; }
      }
    }
  }

//  cout << "no points in the triangle\n";
 
  if (ndPtr == NULL) return true;

  if (ndPtr->dx.n == 0) return true; 

  double *dx = ndPtr->dx.x;

//  cout << dx[0] << "," << dx[1] << " dx\n";
//  cout << "" << (x1[0]-x0[0])*dx[0]+(x1[1]-x0[1])*dx[1] << ","
//       << (x1[0]-x0[0])*dx[2]+(x1[1]-x0[1])*dx[3] << " dot\n";

  if ((x1[0]-x0[0])*dx[0]+(x1[1]-x0[1])*dx[1] > 0.)
    if ((x1[0]-x0[0])*dx[2]+(x1[1]-x0[1])*dx[3] > 0.)
    {
      i = 0; while (i < ndPtr->node.n) if (ndPtr->node[i] != (void*)ndPtr0 && ndPtr->node[i] != (void*)ndPtr1) i++; else break;
      if (i == ndPtr->node.n) return false;
    }

  return true;

}










void AdvancingFront2D::remove3in1(void)
{
  int i = nBndNode, j;
 
  Vector<void*> patch;

  Node2D *cDel;

  while (i < node.n)
  {
    if (node[i].node.n == 3)
    {
      if (node[i].poly.n != 3) prgError(1,"AdvancingFront2D::remove3in1","fatal error!");

      cDel = &(node[i]);

      patch.free(); for (j=0; j<3; j++) patch.append(cDel->poly[j]);

      poly.add(new Poly2D(patch,node,edge,poly));
    }
    else i++;
  }

  return;
}










void AdvancingFront2D::remove4in2(void)
{
  int i = nBndNode, j;
 
  Vector<void*> patch;

  Node2D *cDel;

  Poly2D* pNew;

  while (i < node.n)
  {
    if (node[i].node.n == 4)
    {
      if (node[i].poly.n != 4) prgError(1,"AdvancingFront2D::remove4in2","fatal error!");

      j = 0; while (j < 4 && !((Node2D*)(node[i].node[j]))->bndFlg) j++;

      if (j == 4)
      {
        cDel = &(node[i]);

        patch.free(); for (j=0; j<4; j++) patch.append(cDel->poly[j]);

        pNew = new Poly2D(patch,node,edge,poly);

        poly.add(pNew);

        pNew->splitInTwo(node,edge,poly,(Node2D*)(pNew->node[0]),2);
      }
      else i++;
    }
    else i++;
  }

  return;
}









void AdvancingFront2D::remove8in6(void)
{
  int i = nBndNode, j;
 
  Vector<void*> patch;

  Node2D *c0, *c1;

  Poly2D *pNew;

  for (i=0; i<edge.n; i++)
  {
    c0 = edge[i].c0;
    c1 = edge[i].c1;

    if (!c0->bndFlg && !c1->bndFlg && c0->poly.n == 5 && c1->poly.n == 5)
    {
      patch.free(); for (j=0; j<5; j++) patch.append(c1->poly[j]);

      pNew = new Poly2D(patch,node,edge,poly);

      poly.add(pNew);

      pNew->splitInTwo(node,edge,poly,c0,2);

      if (c0->poly.n != 5) prgError(1,"AdvancingFront2D::remove8in6","fatal error!");

      j = 0; while (((Poly2D*)(c0->poly[j]))->node.n != 4) j++;

      ((Poly2D*)(c0->poly[j]))->splitInTwo(node,edge,poly,c0,2);
    }
  }

  return;
}









void AdvancingFront2D::remove8to1(void)
{
  int i, j, k;

  Vector<void*> patch;

  Node2D *c0, *cNew0, *cNew1;

  Edge2D *eNew, *eNew0, *eNew1;

  Poly2D *pNew, *pNew0, *pNew1, *pNew2;

  for (i=0; i<node.n; i++)
  {
    if (node[i].node.n == 8 && !node[i].bndFlg)
    {
      patch.free(); for (j=0; j<8; j++) patch.append(node[i].poly[j]);

      pNew0 = new Poly2D(patch,node,edge,poly);

      poly.add(pNew0);

      pNew0->maxInternalAngle(&k);

      c0 = (Node2D*)(pNew0->node[k]);

      eNew0 = pNew0->splitInTwo(node,edge,poly,c0,3);

      if (eNew0->p0 == (void*)pNew0) pNew1 = (Poly2D*)eNew0->p1;
      else                           pNew1 = (Poly2D*)eNew0->p0;

      eNew1 = pNew1->splitInTwo(node,edge,poly,c0,3);

      if (eNew1->p0 == (void*)pNew1) pNew2 = (Poly2D*)eNew1->p1;
      else                           pNew2 = (Poly2D*)eNew1->p0;

      cNew0 = subdivideEdge(eNew0);

      cNew1 = subdivideEdge(eNew1);
 
      eNew = pNew0->splitInTwo(node,edge,poly,cNew0,2);

      if (eNew->p0 == (void*)pNew0) pNew = (Poly2D*)eNew->p1;
      else                          pNew = (Poly2D*)eNew->p0;

      pNew->splitInTwo(node,edge,poly,cNew0,2);

      eNew = pNew1->splitInTwo(node,edge,poly,cNew0,2);

      if (eNew->p0 == (void*)pNew0) pNew = (Poly2D*)eNew->p1;
      else                          pNew = (Poly2D*)eNew->p0;

      pNew->splitInTwo(node,edge,poly,cNew0,3);

      pNew->splitInTwo(node,edge,poly,cNew1,2);

      pNew2->splitInTwo(node,edge,poly,cNew1,3);

      pNew2->splitInTwo(node,edge,poly,cNew1,2);
    }
  }
  return;
}









bool AdvancingFront2D::swapEdges(void)
{
  Triangle2D *t0, *t1;

  Node2D *c0, *c1, *c2, *c3;

  Edge2D *e0, *e1, *e2, *e3, *ei;

  int i, j, curr, altn, n0, n1, n2, n3;

  bool swapped = false;

  for (i=0; i<edge.n; i++)
  {
    if (!edge[i].onBoundary())
    {
//      c3--e2---c1
//       |      /|
//       | t0  / |
//       |    /  |
//      e3  ei   e1
//       |  /    |
//       | / t1  |
//       |/      |
//      c0---e0--c2

      ei = &(edge[i]);

      c0 = ei->c0;
      c1 = ei->c1;

      if (!c0->bndFlg || !c1->bndFlg)
      {
        t0 = (Triangle2D*)(ei->p0);
        t1 = (Triangle2D*)(ei->p1);

        t0->triangleInfo(ei,&e2,&e3,&c3);
        t1->triangleInfo(ei,&e0,&e1,&c2);

        n0 = c0->node.n - c0->nOptCnct;
        n1 = c1->node.n - c1->nOptCnct;
        n2 = c2->node.n - c2->nOptCnct;
        n3 = c3->node.n - c3->nOptCnct;

        curr = n0 * n0 + n1 * n1 + n2 * n2 + n3 * n3;

        altn = (n0-1)*(n0-1) + (n1-1)*(n1-1) + (n2+1)*(n2+1) + (n3+1)*(n3+1);

        if (c0->bndFlg) { curr += n0 * n0; altn += (n0-1)*(n0-1); }
        if (c1->bndFlg) { curr += n1 * n1; altn += (n1-1)*(n1-1); }
        if (c2->bndFlg) { curr += n2 * n2; altn += (n2+1)*(n2+1); }
        if (c3->bndFlg) { curr += n3 * n3; altn += (n3+1)*(n3+1); }

        if      (altn < curr && c2->bndFlg && c1->bndFlg && c3->bndFlg)
        { if (triangleAspectRatio(c2->x,c1->x,c3->x) < 0.3) altn = curr + curr; }
        else if (altn < curr && c2->bndFlg && c3->bndFlg && c0->bndFlg)
        { if (triangleAspectRatio(c2->x,c3->x,c0->x) < 0.3) altn = curr + curr; }

        if (altn < curr)
        {
          //cout << c0->x[0] << "," << c0->x[1] << "," << c0->nOptCnct << "\n";
          //cout << c1->x[0] << "," << c1->x[1] << "," << c1->nOptCnct << "\n\n";

          ei->normal[0] = c2->x[1] - c3->x[1];
          ei->normal[1] = c3->x[0] - c2->x[0];

          // update edge->node connectivity

          ei->c0 = c2;
          ei->c1 = c3;

          // update edge->poly connectivity

          if      (e0->p0 == (void*)t1) e0->p0 = (void*)t0;
          else if (e0->p1 == (void*)t1) e0->p1 = (void*)t0;
          else prgError(1,"AdvancingFront2D::swapEdges","fatal error");

          if      (e2->p0 == (void*)t0) e2->p0 = (void*)t1;
          else if (e2->p1 == (void*)t0) e2->p1 = (void*)t1;
          else prgError(2,"AdvancingFront2D::swapEdges","fatal error");

          // update poly->node connectivity

          t0->node[0] = (void*)c2;
          t0->node[1] = (void*)c3;
          t0->node[2] = (void*)c0;

          t1->node[0] = (void*)c2;
          t1->node[1] = (void*)c1;
          t1->node[2] = (void*)c3;

          // update poly->edge connectivity
                               
          t0->edge[0] = (void*)e3;
          t0->edge[1] = (void*)e0;
          t0->edge[2] = (void*)ei; 

          t1->edge[0] = (void*)e2;
          t1->edge[1] = (void*)ei;
          t1->edge[2] = (void*)e1;

          // update node->node connectivity

          j = 0; while (c0->node[j] != (void*)c1) j++; c0->node.del(j);
          j = 0; while (c1->node[j] != (void*)c0) j++; c1->node.del(j);

          c2->node.append((void*)c3);
          c3->node.append((void*)c2);

          // update node->poly connectivity

          j = 0; while (c0->poly[j] != (void*)t1) j++; c0->poly.del(j);
          j = 0; while (c1->poly[j] != (void*)t0) j++; c1->poly.del(j);

          c2->poly.append((void*)t0);
          c3->poly.append((void*)t1);

          t0->checkConsistency(node,edge,poly);
          t1->checkConsistency(node,edge,poly);

          swapped = true;
        }
      }
    }
  }
  return swapped;
}










void AdvancingFront2D::smooth1(int nIter)
{
  int i, j, iter;

  double xa[2];

  for (iter=0; iter<nIter; iter++)
  {
    for (i=nBndNode; i<node.n; i++)
    {
      xa[0] = 0.; 
      xa[1] = 0.;

      for (j=0; j<node[i].node.n; j++)
      {
        xa[0] += ((Node2D*)(node[i].node[j]))->x[0];
        xa[1] += ((Node2D*)(node[i].node[j]))->x[1];
      }
      node[i].x[0] = xa[0] / (double) node[i].node.n;
      node[i].x[1] = xa[1] / (double) node[i].node.n;
    }
  }

  Node2D *c0, *c1;

  for (i=0; i<edge.n; i++)
  {
    c0 = edge[i].c0;
    c1 = edge[i].c1;
    edge[i].normal[0] = c0->x[1] - c1->x[1];
    edge[i].normal[1] = c1->x[0] - c0->x[0];
  }

  return;
}












void AdvancingFront2D::smooth2(int nIter)
{
  int i, j, k, l, iter, it, mxit = 5, nst, ndm = 2, dmyInt;

  double xc[2], x[16], p[2], pp[8], s[64], dp[3], fact, dmyDbl, tol = 1.e-5,
         fact3 = 5., fact4 = 0.42893219;

  Poly2D *polyj;

  Vector<void*> *nd;

  for (i=8; i<16; i++) x[i] = 0.;

  for (iter=0; iter<nIter; iter++)
  {
    for (i=nBndNode; i<node.n; i++)
    {
      xc[0] = node[i].x[0];
      xc[1] = node[i].x[1];

      it = 0;

      while (it < mxit)
      {
        it++;

        p[0] = 0.;
        p[1] = 0.;

        dp[0] = 0.;
        dp[1] = 0.;
        dp[2] = 0.;

        for (j=0; j<node[i].poly.n; j++)
        {
          polyj = (Poly2D*)(node[i].poly[j]);

          k = 0;

          while (polyj->node[k] != (void*)(&(node[i]))) k++;

          nd = &(polyj->node);

          x[0] = ((Node2D*)((*nd)[0]))->x[0];
          x[1] = ((Node2D*)((*nd)[0]))->x[1];
          x[2] = ((Node2D*)((*nd)[1]))->x[0];
          x[3] = ((Node2D*)((*nd)[1]))->x[1];
          x[4] = ((Node2D*)((*nd)[2]))->x[0];
          x[5] = ((Node2D*)((*nd)[2]))->x[1];

          if (polyj->node.n == 3)
          {
            x[6] = 0.; x[7] = 0.;
            for (l=0; l<6; l++)  pp[l] = 0.;
            for (l=0; l<36; l++)  s[l] = 0.;

            nst = 6; fact = fact3;

            ale2d3nodedtriangleaspectratio_(&dmyDbl,x,s,pp,&ndm,&nst,&dmyInt);
          }
          else if (polyj->node.n == 4)
          {
            x[6] = ((Node2D*)((*nd)[3]))->x[0];
            x[7] = ((Node2D*)((*nd)[3]))->x[1];

            nst = 8; fact = fact4;

            ale2d4nodedquadrilateralaspectratio_(&dmyDbl,x,s,pp,&ndm,&nst,&dmyInt);
          }
          else prgError(1,"AdvancingFront2D::smooth2","fatal error!");

          nst++;

          dp[0] += fact * s[(k+k)*nst]; 
          dp[1] += fact * s[(k+k+1)*nst];
          dp[2] += fact * s[(k+k)*nst+1];

          p[0]  += fact * pp[k+k];
          p[1]  += fact * pp[k+k+1];
        }

        fact = 1. / (dp[0] * dp[1] - dp[2] * dp[2]);

        node[i].x[0] += (p[0]*dp[1] - p[1]*dp[2]) * fact;
        node[i].x[1] += (p[1]*dp[0] - p[0]*dp[2]) * fact;

        //cout << sqrt(p[0]*p[0]+p[1]*p[1]) << "\n";

        if (sqrt(p[0]*p[0]+p[1]*p[1]) < tol) break;
      }
      //cout << nst - 1 << " --------------\n";

      if (sqrt(p[0]*p[0]+p[1]*p[1]) > tol) 
      {
        node[i].x[0] = xc[0];
        node[i].x[1] = xc[1];
      }
    }
  }

  Node2D *c0, *c1;

  for (i=0; i<edge.n; i++)
  {
    c0 = edge[i].c0;
    c1 = edge[i].c1;
    edge[i].normal[0] = c0->x[1] - c1->x[1];
    edge[i].normal[1] = c1->x[0] - c0->x[0];
  }

  return;
}











void AdvancingFront2D::convertToQuads(void)
{
  // generate macro mesh by joining triangles to form quads,
  // minimum maxInternalAngle first 
  // (mesh of triangles should be smooth at this stage)

  generateMacroMesh();   // -> triangles + quadrilaterals

  // enhance macro mesh,
  // for the following cases see "M.Sc. Advanced Finite Element Course", page 44/45

  // search for situation 1 and fix it

  enhance1_splitUglyQuadsAtBoundary();

  // search for situation 2 and fix it

  enhance2_joinPolygonsWith2CommonEdges();

  // search for situation 3 and fix it

  enhance3_joinQuadAndTriangle();   // -> triangles + quadrilaterals + pentagons

  // search for situation 4 and fix it

  enhance4_removeQuadAndJoinDiagonals();

  // search for situation 5 and fix it

  cout << "             search for situation 5 and fix it!  NOT YET CODED!\n\n";

  // search for situation 6 and fix it

  cout << "             search for situation 6 and fix it!  NOT YET CODED!\n\n";

  // subdivide all edges

  cout << "             subdivide all edges ...\n\n";

  subdivideAllEdges();

  // subdivide all polygons

  cout << "             split into Quads ...\n\n";

  int i, n = poly.n;
  for (i=0; i<n; i++) poly[i].splitIntoQuads(node,edge,poly);

  // check for consistency

  for (i=0; i<poly.n; i++) poly[i].checkConsistency(node,edge,poly);

  return;
}












void AdvancingFront2D::generateMacroMesh(void)
{
  int i, j, k;

  Poly2D *p0, *p1;

  Vector<void*> patch;

  for (i=0; i<edge.n; i++)
  {
    if (!edge[i].onBoundary())
    {
      p0 = (Poly2D*)(edge[i].p0);
      p1 = (Poly2D*)(edge[i].p1);
      j = 0; while (p0->edge[j] != (void*)(&(edge[i]))) j++;
      k = 0; while (p1->edge[k] != (void*)(&(edge[i]))) k++; k += 2; if (k > 3) k = 1;
      p1->node.append(p0->node[j]);
      p1->node.move(3,k);
      edge[i].length = p1->maxInternalAngle();
      p1->node.del(k);
    }
    else edge[i].length = -1.;
  }

  edge.quickSort();

  i = 0;
  while (i < edge.n)
  {
    p0 = (Poly2D*)(edge[i].p0);
    p1 = (Poly2D*)(edge[i].p1);

    if (p0 != NULL && p1 != NULL && p0->node.n == 3 && p1->node.n == 3)
    {
      patch.free();

      patch.append((void*)p0);
      patch.append((void*)p1);

      poly.add(new Poly2D(patch,node,edge,poly));
    }
    else i++;
  }

  return;
}












void AdvancingFront2D::enhance1_splitUglyQuadsAtBoundary(void)
{
  int i;

  double *x0, *x1, *x2, v1[2], v2[2];

  Poly2D *p;

  for (i=0; i<nBndNode; i++)
  {
    if (!node[i].bndFlg) 
      prgError(1,"AdvancingFront2D::enhance1_splitUglyQuadsAtBoundary","fatal error!");

    if (node[i].node.n == 2)
    {
      if (node[i].poly.n == 1) p = (Poly2D*)(node[i].poly[0]);
      else prgError(2,"AdvancingFront2D::enhance1_splitUglyQuadsAtBoundary","fatal error!");

      if (p->node.n == 4)
      {
        x0 = ((Node2D*)(node[i].node[0]))->x;
        x1 =            node[i]           .x;
        x2 = ((Node2D*)(node[i].node[1]))->x;

        v1[0] = x1[0] - x0[0];
        v1[1] = x1[1] - x0[1];
        v2[0] = x2[0] - x1[0];
        v2[1] = x2[1] - x1[1];

        if (boundaryAngleBetweenTwoEdges2D(v1,v2,2) > 2.6179939) 
          p->splitInTwo(node,edge,poly,&(node[i]),2);
      }
    }
  }
  return;
}













void AdvancingFront2D::enhance2_joinPolygonsWith2CommonEdges(void)
{
  int i = 0;

  Node2D *cDel;

  Vector<void*> patch;

  while (i < node.n)
  {
    if (!node[i].bndFlg && node[i].poly.n == 2)
    { 
      cout << " internal node with poly.n = 2 found ! -> join polygons\n";
      cout << " CODED, BUT NOT YET TESTED!\n\n";

      cDel = &(node[i]);

      patch.free();

      patch.append(cDel->poly[0]);
      patch.append(cDel->poly[1]);

      poly.add(new Poly2D(patch,node,edge,poly));
    }
    else i++;
  }

  return;
}









void AdvancingFront2D::enhance3_joinQuadAndTriangle(void)
{
  int i = 0, j, k;

  Edge2D *e;

  Poly2D *p0, *p;

  Vector<void*> patch;

  double fact, fact1;

  smooth1(2);
  smooth2(1);

  while (i < poly.n)
  {
    if (poly[i].node.n == 3)
    {
      p0   = NULL;
      fact = 10.;
      for (j=0; j<3; j++)
      {
        e = (Edge2D*)(poly[i].edge[j]);
        if (!e->onBoundary())
        {
          if (e->p0 == (void*)(&(poly[i]))) p = (Poly2D*)(e->p1); else p = (Poly2D*)(e->p0);
          if (p->node.n == 4)
          {
            k = 0;
            while (p->edge[k] != (void*)e) k++; 
            p->node.append(poly[i].node[j]);
            p->node.move(4,++k);
            fact1 = p->maxInternalAngle();
            if (fact1 < fact) { fact = fact1; p0 = p; }
            p->node.del(k);
          }
        }
      }
      if (p0 != NULL && fact < 2.6179939)
      {
        patch.free();
        patch.append((void*)(&(poly[i])));
        patch.append((void*) p0);
        poly.add(new Poly2D(patch,node,edge,poly));
        i -= 1;
      }
    }
    i++;
  }

  return;
}









void AdvancingFront2D::enhance4_removeQuadAndJoinDiagonals(void)
{
  int i = 0;

  Node2D *d0, *d1;

  while (i < poly.n)
  {
    if (poly[i].node.n == 4)
    {
      d0 = (Node2D*)(poly[i].node[0]);
      d1 = (Node2D*)(poly[i].node[2]);

      if (!d0->bndFlg && !d1->bndFlg && d0->poly.n == 3 && d1->poly.n == 3)
        enhance4_help(poly[i],0,2);
      else
      {
        d0 = (Node2D*)(poly[i].node[1]);
        d1 = (Node2D*)(poly[i].node[3]);

        if (!d0->bndFlg && !d1->bndFlg && d0->poly.n == 3 && d1->poly.n == 3)
          enhance4_help(poly[i],1,3);
      }
    }
    i++;
  }
  return;
}












void AdvancingFront2D::enhance4_help(Poly2D &p, int i0, int i2)
{
  int i1 = i0 + 1,
      i3 = i2 + 1, j;
 
  if (i3 > 3) i3 = 0;

  Node2D *c0 = (Node2D*)(p.node[i0]),
         *c1 = (Node2D*)(p.node[i1]),
         *c2 = (Node2D*)(p.node[i2]),
         *c3 = (Node2D*)(p.node[i3]),
         *c4;

  Edge2D *e0 = (Edge2D*)(p.edge[i0]),
         *e1 = (Edge2D*)(p.edge[i1]),
         *e2 = (Edge2D*)(p.edge[i2]),
         *e3 = (Edge2D*)(p.edge[i3]),
         *e4;

  Poly2D *p0, *p3;

  void *vp = (void*)(&p);

  if (!e0->is(c0,c1)) prgError(1,"AdvancingFront2D::removeQuadAndJoinDiagonals","fatal error!");
  if (!e1->is(c1,c2)) prgError(2,"AdvancingFront2D::removeQuadAndJoinDiagonals","fatal error!");
  if (!e2->is(c2,c3)) prgError(3,"AdvancingFront2D::removeQuadAndJoinDiagonals","fatal error!");
  if (!e3->is(c3,c0)) prgError(4,"AdvancingFront2D::removeQuadAndJoinDiagonals","fatal error!");

  // update e1

  if      (e0->p0 == vp) p0 = (Poly2D*)(e0->p1);
  else if (e0->p1 == vp) p0 = (Poly2D*)(e0->p0);
  else prgError(5,"AdvancingFront2D::removeQuadAndJoinDiagonals","fatal error!");

  if      (e1->p0 == vp) e1->p0 = (void*)p0;
  else if (e1->p1 == vp) e1->p1 = (void*)p0;
  else prgError(6,"AdvancingFront2D::removeQuadAndJoinDiagonals","fatal error!");

  // update e2

  if      (e3->p0 == vp) p3 = (Poly2D*)(e3->p1);
  else if (e3->p1 == vp) p3 = (Poly2D*)(e3->p0);
  else prgError(7,"AdvancingFront2D::removeQuadAndJoinDiagonals","fatal error!");

  if      (e2->p0 == vp) e2->p0 = (void*)p3;
  else if (e2->p1 == vp) e2->p1 = (void*)p3;
  else prgError(8,"AdvancingFront2D::removeQuadAndJoinDiagonals","fatal error!");

  // update c1

  j = 0; while (c1->node[j] != (void*)c0) j++; c1->node.del(j);

  j = 0; while (c1->poly[j] != (void*)vp) j++; c1->poly.del(j);

  // update c3

  j = 0; while (c3->node[j] != (void*)c0) j++; c3->node.del(j);

  j = 0; while (c3->poly[j] != (void*)vp) j++; c3->poly.del(j);

  // update p0

  j = 0;
  while (p0->node[j] != (void*)c0) j++; 
  p0->node[j] = (void*)c2;

  j = 0;
  while (p0->edge[j] != (void*)e0) j++; 
  p0->edge[j] = (void*)e1;

  // update p3 and find e4

  j = 0;
  while (p3->node[j] != (void*)c0) j++; 
  p3->node[j] = (void*)c2;

  j = 0;
  while (p3->edge[j] != (void*)e3) j++; 
  p3->edge[j] = (void*)e2;

  if (j == 0) j = p3->edge.n;
  e4 = (Edge2D*)(p3->edge[j-1]);

  // update e4 and find c4

  if      (e4->c0 == c0) { e4->c0 = c2; c4 = e4->c1; }
  else if (e4->c1 == c0) { e4->c1 = c2; c4 = e4->c0; }
  else prgError(10,"AdvancingFront2D::removeQuadAndJoinDiagonals","fatal error!");

  // update c4

  j = 0;
  while (c4->node[j] != (void*)c0) j++;
  c4->node[j] = (void*)c2;

  // update c2

  c2->node.append((void*)c4);

  c2->poly.append((void*)p0);
  c2->poly.append((void*)p3);
 
  j = 0; 
  while (c2->poly[j] != vp) j++;
  c2->poly.del(j);

  // delete e0, e3, c0, p

  edge.del(e0);
  edge.del(e3);

  node.del(c0);

  poly.del(&p);

  return;
}











void AdvancingFront2D::subdivideAllEdges(void)
{
  int i, k = edge.n;

  for (i=0; i<k; i++) subdivideEdge(&(edge[i]));

  return;
}









Node2D* AdvancingFront2D::subdivideEdge(Edge2D *e)
{
  int j;

  Node2D *c0, *c1, *cNew;

  Edge2D *eNew;

  Poly2D *p0, *p1;

  void *spline;

  c0 = e->c0;
  c1 = e->c1;

  p0 = (Poly2D*)(e->p0);
  p1 = (Poly2D*)(e->p1);

  if (p0 != NULL) p0->checkConsistency(node,edge,poly);
  if (p1 != NULL) p1->checkConsistency(node,edge,poly);

  // generate new edge

  eNew = new Edge2D;

  edge.add(eNew);

  // generate new node or get midpoint (boundary edges)

  if (!e->onBoundary()) 
  {
    cNew = new Node2D;

    node.add(cNew);

    cNew->x[0] = .5 * (c0->x[0] + c1->x[0]);
    cNew->x[1] = .5 * (c0->x[1] + c1->x[1]);
  }
  else
  {
    if (e->midpoint.n == 1) 
    {
      cNew = (Node2D*)(e->midpoint[0]);
      e->midpoint.free();
    }
    else if (e->midpoint.n == 3)
    {
      cNew = (Node2D*)(e->midpoint[1]);
      eNew->midpoint.append(e->midpoint[0]);
      e->midpoint.del(0);
      e->midpoint.del(0);
    }
    else prgError(1,"AdvancingFront2D::subdivideEdge","fatal error!");
  }

  // update node->poly connectivity

  if (p0 != NULL) cNew->poly.append((void*)p0);
  if (p1 != NULL) cNew->poly.append((void*)p1);

  // update node->node connectivity

  cNew->node.append((void*)c0);
  cNew->node.append((void*)c1);

  j = 0; while (c0->node[j] != (void*)c1) j++; c0->node[j] = (void*)cNew;
  j = 0; while (c1->node[j] != (void*)c0) j++; c1->node[j] = (void*)cNew;

  // update edge->node connectivity

  eNew->c0 = c0;
  eNew->c1 = cNew;

  e->c0    = cNew;

  // update edge->poly connectivity

  eNew->p0 = e->p0;
  eNew->p1 = e->p1;

  // update polygons

  if (p0 != NULL)
  {
    // rearrange nodes and edges for triangles
    if (p0->node.n == 3) p0->node.move(0,2);

    // update poly->node connectivity
    j = 0;
    while (p0->node[j] != (void*)c0) j++;
    p0->node.append((void*)cNew);
    p0->node.move(p0->node.n-1,j+1);

    // update poly->edge connectivity
    j = 0;
    while (p0->edge[j] != (void*)e) j++;
    p0->edge.append((void*)eNew);
    p0->edge.move(p0->edge.n-1,j);

    p0->checkConsistency(node,edge,poly);
  }

  if (p1 != NULL)
  {
    // rearrange nodes and edges for triangles
    if (p1->node.n == 3) p1->node.move(0,2);

    // update poly->node connectivity
    j = 0;
    while (p1->node[j] != (void*)c1) j++;
    p1->node.append((void*)cNew);
    p1->node.move(p1->node.n-1,j+1);

    // update poly->edge connectivity
    j = 0;
    while (p1->edge[j] != (void*)e) j++;
    p1->edge.append((void*)eNew);
    p1->edge.move(p1->edge.n-1,j+1);

    p1->checkConsistency(node,edge,poly);
  }

  return cNew;
}








void AdvancingFront2D::correct6to1(void)
{
  int i, j, k, l, m, n, 
      addCnctvy4by2[12] = {1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1}, mn4by2, kmn4by2,
      addCnctvy3by3[12] = {1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0}, mn3by3, kmn3by3;

  Vector<void*> patch;

  Node2D *ni, *c0, *c1;

  Edge2D *ei;

  Poly2D *pNew;

  bool changed = true;

  while (changed)
  {
    changed = false;

    for (i=0; i<node.n; i++)
    {
      if (node[i].poly.n == 6 && !node[i].bndFlg)
      {
        j = 0;
        while (j<node[i].node.n && !((Node2D*)(node[i].node[j]))->bndFlg) j++;
        if (j == node[i].node.n)
        {
          patch.free();

          for (j=0; j<node[i].poly.n; j++) patch.append(node[i].poly[j]);

          pNew = new Poly2D(patch,node,edge,poly);

          poly.add(pNew);

          for (j=0; j<6; j++) pNew->node.append(pNew->node[j]);
          mn4by2 = 1000;
          for (k=0; k<6; k++)
          {
            m = 0;
            for (j=0; j<12; j++) 
            {
              l = ((Node2D*)(pNew->node[j+k]))->node.n + addCnctvy4by2[j] - 4;
              m += l*l;
            }
            if (m < mn4by2)
            {
              mn4by2  = m;
              kmn4by2 = k;
            }
          }
          mn3by3 = 1000;
          for (k=0; k<3; k++)
          {
            m = 0;
            for (j=0; j<12; j++) 
            {
              l = ((Node2D*)(pNew->node[j+k]))->node.n + addCnctvy3by3[j] - 4;
              m += l*l;
            }
            if (m < mn3by3)
            {
              mn3by3  = m;
              kmn3by3 = k;
            }
          }
          for (j=0; j<6; j++) pNew->node.del(12);

          if (mn4by2 < mn3by3)
          {
            ei = pNew->splitInTwo(node,edge,poly,(Node2D*)(pNew->node[kmn4by2]),6);
            ni = subdivideEdge(ei);       
            if (ni->poly.n != 2) prgError(1,"AdvancingFront2D::correct6to1","fatal error!");
            ((Poly2D*)(ni->poly[0]))->splitIntoQuads(node,edge,poly);
            ((Poly2D*)(ni->poly[1]))->splitIntoQuads(node,edge,poly);
          }
          else
          {
            ei = pNew->splitInTwo(node,edge,poly,(Node2D*)(pNew->node[kmn3by3]),10);
            ni = subdivideEdge(ei);
            ei = pNew->splitInTwo(node,edge,poly,ni,10);
            ni = subdivideEdge(ei);
            ei = pNew->splitInTwo(node,edge,poly,ni,9);
            ei = pNew->splitInTwo(node,edge,poly,ni,8);
            ni = subdivideEdge(ei);
            ei = pNew->splitInTwo(node,edge,poly,ni,7);
            ei = pNew->splitInTwo(node,edge,poly,ni,6);
            ni = subdivideEdge(ei);
            ei = pNew->splitInTwo(node,edge,poly,ni,5);
            ei = pNew->splitInTwo(node,edge,poly,ni,3);
          }
          changed = true;
        }
      }
    }
  }

  return;
}










void AdvancingFront2D::correct2to1(void)
{
  int i, j;

  Vector<void*> patch;

  for (i=0; i<node.n; i++)
  {
    if (node[i].poly.n == 2 && !node[i].bndFlg)
    {
      patch.free();

      for (j=0; j<node[i].poly.n; j++) patch.append(node[i].poly[j]);

      poly.add(new Poly2D(patch,node,edge,poly));
    }
  }
  return;
}








void AdvancingFront2D::correct7to1(void)
{
  Vector<void*> patch;

  int i, j, k, l, m,  
      addCnctvy[14] = {1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0}, mn, kmn;

  Poly2D* pNew;

  Edge2D* ei;

  Node2D* ni;

  for (i=0; i<node.n; i++)
  {
    if (node[i].poly.n == 7 && !node[i].bndFlg)
    {
      j = 0;
      while (j<node[i].node.n && !((Node2D*)(node[i].node[j]))->bndFlg) j++;
      if (j == node[i].node.n)
      {
        patch.free();

        for (j=0; j<node[i].poly.n; j++) patch.append(node[i].poly[j]);

        pNew = new Poly2D(patch,node,edge,poly);

        poly.add(pNew);

        for (j=0; j<14; j++) pNew->node.append(pNew->node[j]);
        mn = 1000;
        for (k=0; k<14; k++)
        {
          m = 0;
          for (j=0; j<14; j++) 
          {
            l = ((Node2D*)(pNew->node[j+k]))->node.n + addCnctvy[j] - 4;
            m += l*l;
          }
          if (m < mn)
          {
            mn  = m;
            kmn = k;
          }
        }
        for (j=0; j<14; j++) pNew->node.del(14);

        ei = pNew->splitInTwo(node,edge,poly,(Node2D*)(pNew->node[kmn]),12);
        ni = subdivideEdge(ei);
        ei = pNew->splitInTwo(node,edge,poly,ni,12);
        ni = subdivideEdge(ei);
        ei = pNew->splitInTwo(node,edge,poly,ni,12);
        ni = subdivideEdge(ei);
        ei = pNew->splitInTwo(node,edge,poly,ni,11);
        ei = pNew->splitInTwo(node,edge,poly,ni,9);
        ei = pNew->splitInTwo(node,edge,poly,ni,8);
        ni = subdivideEdge(ei);
        ei = pNew->splitInTwo(node,edge,poly,ni,8);
        ni = subdivideEdge(ei);
        ei = pNew->splitInTwo(node,edge,poly,ni,7);
        ei = pNew->splitInTwo(node,edge,poly,ni,6);
        ni = subdivideEdge(ei);
        ei = pNew->splitInTwo(node,edge,poly,ni,5);
        ei = pNew->splitInTwo(node,edge,poly,ni,4);
        ni = subdivideEdge(ei);
        ei = pNew->splitInTwo(node,edge,poly,ni,3);
      }
    }
  }
  return;
}








void AdvancingFront2D::generateMidpoints(void)
{
  int i, j;

  double *X;

  Node2D *cNew;

  for (i=0; i<poly.n; i++)
  {
    cNew = new Node2D;

    node.add(cNew);

    X = cNew->x;

    X[0] = 0.;
    X[1] = 0.;

    for (j=0; j<poly[i].node.n; j++)
    {
      X[0] += ((Node2D*)(poly[i].node[j]))->x[0];
      X[1] += ((Node2D*)(poly[i].node[j]))->x[1];
    }

    X[0] /= (double)(poly[i].node.n);
    X[1] /= (double)(poly[i].node.n);

    poly[i].node.append((void*)cNew);
  }

  return;
}











void AdvancingFront2D::draw(void)
{
  int i;

  for (i=0; i<frontEdge.n; i++) frontEdge[i].draw();

  return;
}









void AdvancingFront2D::drawPolygons(void)
{
  int i;

  for (i=0; i<poly.n; i++) poly[i].draw();

  return;
}










void AdvancingFront2D::drawAllConnectivities(void)
{
  int i, j;

  cout << "  next: draw polygons\n";
  prgUpdateDisplay(true);
  plot.wipe();
  plot.setColour(3);
  drawPolygons();
  
  cout << "  next: draw edges\n";
  prgUpdateDisplay(true);
  plot.setColour(4);
  for (i=0; i<edge.n; i++) edge[i].draw();

  cout << "  next: draw node->node connectivities\n";
  prgUpdateDisplay(true);
  plot.setColour(5);
  for (i=0; i<node.n; i++)
    for (j=0; j<node[i].node.n; j++)
      plot.line(node[i].x,((Node2D*)(node[i].node[j]))->x);

  cout << "  next: draw edges\n";
  prgUpdateDisplay(true);
  plot.setColour(4);
  for (i=0; i<edge.n; i++) edge[i].draw();

  cout << "  next: draw polygons\n";
  prgUpdateDisplay(true);
  plot.wipe();
  plot.setColour(3);
  drawPolygons();

  return;
}

