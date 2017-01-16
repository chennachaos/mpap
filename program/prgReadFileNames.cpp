
#include <iostream>

#include "FunctionsProgram.h"
#include "Files.h"
#include "Definitions.h"
#include "Debug.h"
#include "MyString.h"


extern Files files;
extern bool  lastProj;


using namespace std;


bool prgReadFileNames(void)
{
  MyString answ;
  char *no [] = NOO;
  char *yes [] = YESS;
  char *quit [] = QUIT;

  bool ok = false;
  
  
  if (prgLastProject())
  {  
     cout << "       project directory   : " << files.projDir << "\n";
     cout << "       input file name     : " << files.Ifile << "\n";
     cout << "       output file name    : " << files.Ofile << "\n";
     cout << "       time plot file name : " << files.Tfile << "\n";
     cout << "       eps plot file name  : " << files.Pfile << "\n\n";  

     if (!lastProj)
     {
       do { cout << "       Are file names correct ? (y/n/x) "; answ.input(); } 
       while ( answ.which(no)<0 && answ.which(yes)<0 && answ.which(quit)<0 );
       cout << "\n\n";
       if (answ.which(quit) >= 0) return false;
       if (answ.which(yes)  >= 0 && prgFileExist(files.projDir,files.Ifile)) ok = true;
       else if (answ.which(no) < 0) cout << "       No such project directory or input file!\n\n";
     }
     else
     {
       if (prgFileExist(files.projDir,files.Ifile)) return true; 
       else { cout << "       No such project directory or input file!\n"; return false; }
     }
  }

  
  while (!ok) {
	  
     cout << "       project directory   : "; files.projDir.input(); 
     do { cout << "       input file name     : "; files.Ifile.input(); } 
       while (files.Ifile.length()<1);

     files.Ofile.free().append(files.Ifile)[0] = 'O';
     files.Tfile.free().append(files.Ifile)[0] = 'T';
     files.Pfile.free().append(files.Ifile)[0] = 'P';

     cout <<"       output file name    : "<<files.Ofile<<" : ";files.Ofile.inputKeepIfReturn();
     cout <<"       time plot file name : "<<files.Tfile<<" : ";files.Tfile.inputKeepIfReturn();
     cout <<"       eps plot file name  : "<<files.Pfile<<" : ";files.Pfile.inputKeepIfReturn(); 
     cout <<"\n";
    
     do { cout << "       Are file names correct ? (y/n/x) "; answ.input(); }
     while ( answ.which(no)<0 && answ.which(yes)<0 && answ.which(quit)<0 );
     cout << "\n\n";   
     if (answ.which(quit) >= 0) return false;
     if (answ.which(yes)  >= 0 && prgFileExist(files.projDir,files.Ifile)) ok = true;
     else if (answ.which(no) < 0) cout << "       No such project directory or input file!\n\n";
  }


  if (debug) 
    {  cout << " " << files.projDir << "\n";
       cout << " " << files.Ifile << "\n";
       cout << " " << files.Ofile << "\n";
       cout << " " << files.Tfile << "\n";
       cout << " " << files.Pfile << "\n\n"; }

       
  prgWriteLastProject();
       
  return true;
}


