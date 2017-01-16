
#ifndef incl_SelectNode_h
#define incl_SelectNode_h

#include "List.h"
#include "MathVector.h"


class SelectNode
{
  public:

    SelectNode(void);
    virtual ~SelectNode();
 
    int flag;               // how to select

    double cosAlpha, alpha; // angle to define surfaces / edges
   
    bool defm,              // select from deformed mesh
         prnt,              // print nodeList
         show,              // plot nodes
         numb,              // plot node numbers
         dslt;              // deselect flag

    double x[7];            // search volume data

    Vector<int> nodeList;   // list of currently selected nodes

    void free(void);
};


#endif



