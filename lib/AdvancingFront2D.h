
#ifndef incl_AdvancingFront2D_h
#define incl_AdvancingFront2D_h

#include "MathVector.h"
#include "MathMatrix.h"
#include "GeomSurface.h"
#include "Domain.h"
#include "List.h"



class Node2D: public ListItem
{
  public:

    Node2D(void);

    virtual ~Node2D();

    double x[3];   // 0,1 -> coordinates;  2 -> h

    int np, nOptCnct;

    bool bndFlg;

    Vector<void*> node, poly;

    VectorArray<double> dx;   // vectors of initial front edges

    void draw(void);
};









class Edge2D: public ListItem
{
  public:

    Edge2D(void);

    Edge2D(Node2D*, Node2D*, double, double*, double*);

    virtual ~Edge2D();

    Node2D *c0, *c1;

    Vector<void*> midpoint;

    void *p0, *p1;

    double length, normal[2];

    bool operator<(Edge2D &edge) { return (length < edge.length); }

    bool operator>(Edge2D &edge) { return (length > edge.length); }

    bool is(Node2D *k0, Node2D *k1) { return ((c0 == k0 && c1 == k1) || (c0 == k1 && c1 == k0)); }

    bool isStrictly(Node2D *k0, Node2D *k1) { return (c0 == k0 && c1 == k1); }

    bool onBoundary(void) { return (p0 == NULL || p1 == NULL); }

    bool isBetween(void *g0, void *g1) { return ((p0==g0 && p1==g1) || (p0==g1 && p1==g0)); }

    void draw(void);
};








class Poly2D: public ListItem
{
  public:

    Poly2D(void) { return; }

    Poly2D(Node2D*, Node2D*, Node2D*, Edge2D*, Edge2D*, Edge2D*);

    Poly2D(Vector<void*>&, ListInfinite<Node2D>&, ListInfinite<Edge2D>&, ListInfinite<Poly2D>&);

    virtual ~Poly2D();

    Vector<void*> node, edge;

    void triangleInfo(Edge2D *e0, Edge2D **e1, Edge2D **e2, Node2D **c2);

    bool checkConsistency(ListInfinite<Node2D>&, ListInfinite<Edge2D>&, ListInfinite<Poly2D>&);

    Edge2D* splitInTwo(ListInfinite<Node2D>&, ListInfinite<Edge2D>&, ListInfinite<Poly2D>&, 
                    Node2D* cIn = NULL, int n1 = 0);

    void splitIntoQuads(ListInfinite<Node2D>&, ListInfinite<Edge2D>&, ListInfinite<Poly2D>&);

    void arrangeIfTriangle(void);

    double maxInternalAngle(int *nn = NULL);

    double aspectRatio(void);

    void draw(bool colFlg = false);
};









class AdvancingFront2D
{ 
  public:

    AdvancingFront2D(void);

    virtual ~AdvancingFront2D();

    void generateMesh(MatrixFullArray<double>&, MatrixFullArray<int>&, void*, int&,
                      GeomSurface&, Domain*, int, bool showFlg = false);

    void drawPolygons(void);

  private:

    int nBndNode;
    
    ListInfinite<Node2D> node;

    ListInfinite<Edge2D> frontEdge, edge;

    ListInfinite<Poly2D> poly;

    int nFail, nEdge;

    void generateInitialFront(void);

    int  generateNextTriangle(void);

    bool triangleOK(double*, double*, double*, Node2D *ndPtr = NULL, Node2D *ndPtr0 = NULL, Node2D *ndPtr1 = NULL);

    void remove3in1(void);

    void remove4in2(void);

    void remove8in6(void);

    void remove8to1(void);

    bool swapEdges(void);

    void smooth1(int);

    void smooth2(int);

    void smoothTrianglesOnly(int);

    void smoothQuadsOnly(int);

    void convertToQuads(void);

    void generateMacroMesh(void);

    void enhance1_splitUglyQuadsAtBoundary(void);

    void enhance2_joinPolygonsWith2CommonEdges(void);

    void enhance3_joinQuadAndTriangle(void);

    void enhance4_removeQuadAndJoinDiagonals(void);

    void enhance4_help(Poly2D&, int, int);

    Node2D* subdivideEdge(Edge2D*);

    void subdivideAllEdges(void);

    void generateMidpoints(void);

    void correct6to1(void);

    void correct2to1(void);

    void correct7to1(void);

    void showAndWait(void);

    void draw(void);

    void drawAllConnectivities(void);
};

#endif



