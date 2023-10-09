
#ifndef incl_MacroDialogData_h
#define incl_MacroDialogData_h


#include "MyStringList.h"
#include "List.h"
#include "MathVector.h"
#include "Definitions.h"


enum {TXTF, TBTN, RBOX, LABL, LIST };



class MacroDialogData
{
  public: 
	   
    MacroDialogData(void); 
    ~MacroDialogData(); 
       
    bool selDm, frmRBoxFlg, frmBBoxFlg;
   
    MyStringList  strgList, labl, tBtn, txtF;
    
    MyString      dfltStrg, strgTxtFLabl;
    
    List<MyStringList> list, rBox;
    
    myVector<double> dflt;
    myVector<int>    inputType, nBttn, txtFCol, bttnCol;
   
    int strgTxtFCol;
    
    bool operator!(void);
    
    void selectDomain(void);
   
    void stringList(char *s1, char *s2 = NULL, char *s3 = NULL, char *s4 = NULL, 
         	   char *s5 = NULL, char *s6 = NULL, char *s7 = NULL, char *s8 = NULL, 
         	   char *s9 = NULL, char *s10 = NULL, char *s11 = NULL, char *s12 = NULL); 
   
    void stringTextField(char *lbl, char *dflt = NULL, int c = 15);
    
    void addList(char *s1, char *s2 = NULL, char *s3 = NULL, char *s4 = NULL, 
         	char *s5 = NULL, char *s6 = NULL, char *s7 = NULL, char *s8 = NULL, 
         	char *s9 = NULL, char *s10 = NULL, char *s11 = NULL, char *s12 = NULL); 
   
    void addList(char **s);
   
    void addRadioBox(char *s1, char *s2 = NULL, char *s3 = NULL, char *s4 = NULL, 
         	    char *s5 = NULL, char *s6 = NULL, char *s7 = NULL, char *s8 = NULL, 
         	    char *s9 = NULL, char *s10 = NULL, char *s11 = NULL, char *s12 = NULL);
    
    void addLabel(char *s);
   
    void addToggleButton(char *s, bool dflt = false);
   
    void addTextField(char *s, double x = 0, int c = 4);
   
    void nextButtonBox(void);
    
    void setButtonColumnDim(int c1, int c2 = 4, int c3 = 4, int c4 = 4, int c5 = 4, 
         	           int c6 = 4, int c7 = 4, int c8 = 4, int c9 = 4, int c10 = 4,
         		   int c11 = 4, int c12 = 4); 
   
    void frameRadioBox(void);
    
    void frameButtonBox(void);
   
    void getDfltWords(void); 
   
  private:
   
    void autoNextButtonBox(void);
   
};


#endif

