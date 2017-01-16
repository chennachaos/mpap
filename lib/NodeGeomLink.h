
#ifndef incl_NodeGeomLink_h
#define incl_NodeGeomLink_h

#include "GeomPoint.h"
#include "GeomSpline.h"
#include "GeomSurface.h"
#include "MathVector.h"




class NodeGeomLink: public ListItem
{
  public:

    NodeGeomLink(void);

    NodeGeomLink(NodeGeomLink &);

    virtual ~NodeGeomLink() { return; }

    int id;

    Vector<void*> geomObj;

    void isOn   (void*, void *g1 = NULL, void *g2 = NULL, void *g3 = NULL);
    void isNotOn(void*, void *g1 = NULL, void *g2 = NULL, void *g3 = NULL);

    bool onThis(GeomObject*);
    bool onThis(void*);

    bool onSurface(void);
    GeomSurface* whichSurface(void);

    bool onSpline(void);
    GeomSpline* whichSpline(void);

    bool onPoint(void);
    GeomPoint* whichPoint(void);    
};



#endif



