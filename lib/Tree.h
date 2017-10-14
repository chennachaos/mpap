
#ifndef incl_Tree_h
#define incl_Tree_h


#include <iostream>
#include <cstring>
#include <string.h>


using namespace std;


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

//  tree item base class

class TreeItem
{
  public:

    TreeItem *parent, **child;

    void     printTree(int indent)
             {  
               int i;
               if (parent != NULL && parent->child[0] != this)
                 for (i=1;i<indent;i++) std::cout << " ";
               std::cout << " --> " << printName();
               i = 0;
               while (child[i] != NULL) 
               {
                 child[i]->printTree(indent+strlen(printName())+5); 
                 if (child[++i] != NULL) std::cout << "\n";
               }
               return;
             }
  private:

    virtual char *printName(void) { return NULL; }
};


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

// Tree - a tree data structure (slow; should only be used for small trees)

template<class Type> class Tree
{
  public:

    Tree(void);

    virtual ~Tree();

    int n;

    Type* searchItem;

    bool add(Type *, Type *);  // add an object to the tree as child of ...

    void del(Type *);          // remove and delete an object and all its descendents 

    void free(void);

    void print(void);

    void resetTreeSearch(void);

    bool doTreeSearch(void);

  private:

    TreeItem *root;
};

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::



template<class Type> Tree<Type>::Tree(void)
{
  root = NULL; n = 0; return;
}





template<class Type> Tree<Type>::~Tree()
{
  free(); return;
}





template<class Type> bool Tree<Type>::add(Type *newItem, Type *parentItem)
{  
  if (parentItem == NULL && root != NULL) return false;

  int j, k;

  TreeItem **tmp;

  if (root == NULL) 
  {
    root = (TreeItem*)newItem;

    newItem->parent = NULL;
  }
  else
  {
    j = 0; while (parentItem->child[j] != NULL) j++;

    tmp = new TreeItem* [j+2];

    for (k=0; k<j; k++) tmp[k] = parentItem->child[k];

    delete [] parentItem->child;

    tmp[j] = (TreeItem*)newItem;

    tmp[j+1] = NULL;

    parentItem->child = tmp;

    newItem->parent = (TreeItem*)parentItem;
  }

  newItem->child = new (TreeItem*);

  newItem->child[0] = NULL;

  n++;

  return true;
}





template<class Type> void Tree<Type>::del(Type *delItem)
{
  if (root == NULL) return;  

  while (delItem->child[0] != NULL) del((Type*)delItem->child[0]);  

  int j, k, l;

  TreeItem **tmp;

  if (delItem->parent != NULL) 
  { 
    j = 0; while (delItem->parent->child[j] != NULL) j++;

    tmp = new TreeItem* [j];

    l = 0;
    for (k=0; k<=j; k++) 
      if (delItem->parent->child[k] != delItem) tmp[l++] = delItem->parent->child[k];

    delete [] delItem->parent->child;

    delItem->parent->child = tmp;
  }

  delete delItem;

  n--;

  if (n == 0) root = NULL;

  return;
}





template<class Type> void Tree<Type>::free(void)
{
  del((Type*)root);

  return;
}





template<class Type> void Tree<Type>::print(void)
{
  if (root == NULL) return;

  root->printTree(1);

  cout << "\n\n";

  return;
}





template<class Type> void Tree<Type>::resetTreeSearch(void)
{
  searchItem = NULL;

  return;
}




template<class Type> bool Tree<Type>::doTreeSearch(void)
{
  if (searchItem == NULL)
  {
    if (root != NULL) { searchItem = (Type*)root; return true; } else return false;
  }

  if (searchItem->child[0] != NULL) { searchItem = (Type*)(searchItem->child[0]); return true; }

  TreeItem *par = searchItem->parent,
           *itm = (TreeItem*)searchItem;
  int i;
  while (par != NULL)
  {
    i = 0; while (par->child[i] != itm) i++;

    if (par->child[i+1] != NULL) { searchItem = (Type*)(par->child[i+1]); return true; }

    itm = par;
    par = par->parent;
  }
  searchItem = NULL;

  return false;
}



#endif

