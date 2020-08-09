
#include <iostream>
#include <fstream>

#include "MyString.h"
#include "Debug.h"
#include "FunctionsProgram.h"
#include "Definitions.h"
#include "RunControl.h"
#include "MathVector.h"


extern RunControl runCtrl;

using namespace std;


void prgReadRunControl(ifstream &Ifile)
{
   MyString line;

   MyStringList *batch;
   
   int key;

   char *runCtrlCmd[] = RUN_CONTROL_COMMANDS;
   
   line.read(Ifile).stripToMin();

   while (Ifile && line != "END RUN_CONTROL")  
    { 
       key = line.which(runCtrlCmd);

       switch (key) 
        {
   	   case -1: if (line.length()>0) cout << "   ignoring  '" << line << "' ...\n\n"; break;
		    
	   case  0: runCtrl.cmd.append(- runCtrl.batch.n);
		   
                    runCtrl.batch.add(new MyStringList);
			   
                    batch = &(runCtrl.batch[runCtrl.batch.n-1]);
		    
                    batch->addNew("loop,,1");
		    
	            while (Ifile && line != "ctrl,end" && line != "end" && 
				    line != "ctrl,end,,1" && line != "end,,1")
		    {
		       batch->addNew(line.read(Ifile).stripToMin().asCharArray());
		       
		       if (!Ifile || line == "END RUN_CONTROL" || line.which(runCtrlCmd) > -1) 
			  prgError(1,"prgReadRunControl","Syntax error: 'ctrl,end' not found!");
		    }
		    batch->addNew("next"); batch->move(batch->n-1,batch->n-2);

                    //for (int i=0; i<batch->n; i++) cout << (*batch)[i] << "\n"; cout << "\n";
		    
		    break;
		    
           default: runCtrl.cmd.append(key); break;
	}
	   
       line.read(Ifile).stripToMin(); 
    }

   if (!Ifile) prgError(1,"prgReadRunControl","'END RUN_CONTROL' not found!");
   
   if (runCtrl.cmd.n < 1) prgError(1,"prgReadRunControl","control data section is empty!");
   
   if (debug) cout << " number of runControlCommands = " << runCtrl.cmd.n << "\n\n";
   if (debug) cout << " runCtrl.cmd = " << runCtrl.cmd << "\n\n";
   if (debug) for (int i=0; i<runCtrl.batch.n; i++)  { runCtrl.batch[i].print(); cout << "\n"; }
	  
   return;
}



