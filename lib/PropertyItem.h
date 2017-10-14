
#ifndef incl_PropertyItem_h
#define incl_PropertyItem_h


#include "List.h"
#include "MathVector.h"
#include "MyString.h"


class PropertyItem: public ListItem
{
  public:

    PropertyItem(void);

    PropertyItem(int);

    ~PropertyItem();
	  
    int                 id, type;

    VectorArray<double> data;

    MyString            name;

    void readInputData(std::ifstream &, MyString &, char *);    

  private:

};



#endif


