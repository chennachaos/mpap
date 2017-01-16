
#ifndef incl_Files_h
#define incl_Files_h

#include <fstream>

#include "MyString.h"


class Files 

{ 
  public:
	  
    MyString projDir;   // project directory
  
    MyString Ifile;     // input data file

    MyString Ofile;     // output data file
 
    MyString Tfile;     // time data file

    MyString Pfile;     // base name for eps files

    int      nps;       // counter for eps files
  
    std::ofstream Tout, Pout, Oout; // output streams

    void reset(void) { projDir.free();
	               Ifile.free();
		       Ofile.free();
		       Tfile.free();
		       Pfile.free();
		       nps = 0;
                       if (Tout.is_open()) Tout.close();
                       if (Pout.is_open()) Pout.close();
                       if (Oout.is_open()) Oout.close();
		       return; }
};

#endif

