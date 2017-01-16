
#include "SelectNode.h"
#include "DomainTree.h"
#include "Mesh.h"
#include "MathVector.h"
#include "FunctionsProgram.h"


extern DomainTree domain;




SelectNode::SelectNode(void)
{
  return;
}




SelectNode::~SelectNode()
{
  return;
}





void SelectNode::free(void)
{
  nodeList.free();

  return;
}

