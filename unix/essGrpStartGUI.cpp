
#include <Xm/Xm.h> 
#include <Xm/MainW.h>
#include <Xm/RowColumn.h>
#include <Xm/Frame.h>
#include <Xm/LabelG.h>
#include <Xm/DrawingA.h>
#include <Xm/PushB.h>
#include <Xm/PushBG.h>
#include <Xm/CascadeB.h>
#include <Xm/MessageB.h>
#include <Xm/SeparatoG.h>
#include <Xm/Form.h>
#include <Xm/DialogS.h>
#include <Xm/TextF.h>
#include <Xm/List.h>
#include <Xm/ScrolledW.h>

#include "MacroList.h"
#include "FunctionsUnix.h"
#include "FunctionsEssGrp.h"
#include "RunControl.h"
#include "FunctionsProgram.h"
#include "UnixGUI.h"
#include "Plot.h"
#include "Files.h"



extern UnixGUI    unixGUI;
extern MacroList  macro;
extern RunControl runCtrl;
extern bool       noGUI;
extern bool       lastProj;
extern Plot       plot;
extern Files      files;

 

void essGrpStartGUI(int argc, char **argv)
{    
  if (noGUI) prgError(1,"essGrpStartGUI","'noGUI' is set!");
      
  int          i, j, j0 = 0, colour, coloursAllocated;

  Widget       mainWindow,
	       menuBar,
	       drawingAreaFrame, 
	       cmdAndInfoForm, 
	       cmdListSW,
	       commandList,
	       commandText, 
	       infoForm,
	       infoFrame,
	       timeLabel1,
	       timeLabel2,
	       projectLabel,
	       statusLabel,
	       projNameLabel,
	       separator1, 
	       separator2,
	       projectMenuButton,     
	       projectMenu,     
	       helpMenuButton,     
	       closeButton,
	       openButton,
	       exitButton,
	       helpMenu,     
	       helpButton, 
	       aboutButton,
               *macroMenu       = new Widget [macro.ntype],
               *macroMenuButton = new Widget [macro.ntype],
               *macroButton     = new Widget [macro.n];

  Colormap     cmap;
   
  XColor       col, unused;

  char         *namedColours[] = {COLOURS,NULL};
  
  Dimension    width, height;
  
  char         tmp[50];
 
  XEvent       event;

  XmString     XmStr;

  Font         nFont;

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
// CHECK FOR LAST PROJECT IF REQUIRED 
  
  if (lastProj) 
  {
    if (!prgLastProject()) return;
    if (!prgFileExist(files.projDir,files.Ifile)) 
      { std::cout << "    Invalid project directory or input file!\n";  return; }
  }
  
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//  MAIN WINDOW
  
  XtSetLanguageProc(NULL, (XtLanguageProc)NULL, myNULL);

  unixGUI.topLevel = XtVaAppInitialize(&(unixGUI.app), "XMpap2", NULL, 0, &argc, argv, NULL,
		                       XmNtitle, "MPAP2", myNULL);

  mainWindow  = XtVaCreateManagedWidget("mainWindow", xmMainWindowWidgetClass,
		                        unixGUI.topLevel, myNULL);


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//  MENUBAR
  
  menuBar           = XmCreateMenuBar(mainWindow,"menuBar",NULL,0);
 
  // create project pulldown menu
  
  projectMenu       = XmCreatePulldownMenu(menuBar, "projectMenu", NULL, 0);
  
  projectMenuButton = XtVaCreateManagedWidget("project", 
		                              xmCascadeButtonWidgetClass, menuBar, 
		                              XmNsubMenuId, projectMenu, myNULL);

  openButton        = XtVaCreateManagedWidget("open", 
		                              xmPushButtonGadgetClass, projectMenu, myNULL); 

  closeButton       = XtVaCreateManagedWidget("close", 
		                              xmPushButtonGadgetClass, projectMenu, myNULL); 
      
  exitButton        = XtVaCreateManagedWidget("exit", 
		                              xmPushButtonGadgetClass, projectMenu, myNULL); 
      
  // create macro menus

  for (i=0; i<macro.ntype; i++)
  {
     sprintf(tmp,"macroMenu%d",i);
     
     macroMenu[i]       = XmCreatePulldownMenu(menuBar, tmp, NULL, 0);
 
     sprintf(tmp,macro[macro.jp[i]].type.asCharArray());
     
     macroMenuButton[i] = XtVaCreateManagedWidget(tmp, xmCascadeButtonWidgetClass, menuBar, 
		                                  XmNsubMenuId, macroMenu[i], myNULL);
     
     for (j=j0; j<=macro.jp[i]; j++)
     {
	j0 = macro.jp[i] + 1;
	
	macroButton[j] = XtVaCreateManagedWidget(macro[j].name.asCharArray(), 
			                         xmPushButtonWidgetClass, macroMenu[i], myNULL);
     }	
  }

  // create help pulldown menu
  
  helpMenu          = XmCreatePulldownMenu(menuBar, "helpMenu", NULL, 0);
  
  helpMenuButton    = XtVaCreateManagedWidget("help", 
		                              xmCascadeButtonWidgetClass, menuBar, 
		                              XmNsubMenuId, helpMenu, myNULL);
    
  aboutButton       = XtVaCreateManagedWidget("about", 
		                              xmPushButtonWidgetClass, helpMenu, myNULL); 
      
  helpButton        = XtVaCreateManagedWidget("help", 
		                              xmPushButtonWidgetClass, helpMenu, myNULL); 

  XtVaSetValues(menuBar, XmNmenuHelpWidget, helpMenuButton, myNULL);

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//  DRAWING AREA
  
  // compare "Motif Programming Manual" V6A, page 340
  
  drawingAreaFrame = XtVaCreateManagedWidget("drawingAreaFrame",
		                             xmFrameWidgetClass,mainWindow, myNULL);

  unixGUI.drawingArea = XtVaCreateManagedWidget("drawingArea",
		                   xmDrawingAreaWidgetClass,drawingAreaFrame,
			     XtVaTypedArg, XmNbackground, XmRString, "black", 6, 
                             XmNheight,     400,
                             XmNwidth,      600, myNULL);

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//  COMMAND AND INFO WIDGETS

  cmdAndInfoForm = XtVaCreateManagedWidget("cmdAndInfoForm", xmFormWidgetClass, mainWindow, myNULL);
  
  commandText    = XtVaCreateManagedWidget("commandText", xmTextFieldWidgetClass, cmdAndInfoForm,
       		                           XmNleftOffset,       8,
					   XmNrightOffset,      6,
					   XmNbottomOffset,     8,
		                           XmNbottomAttachment, XmATTACH_FORM,     
                                           XmNleftAttachment,   XmATTACH_FORM,
                                           XmNrightAttachment,  XmATTACH_POSITION, 
                                           XmNrightPosition, 65, myNULL);

  cmdListSW      = XtVaCreateManagedWidget("cmdListSW", xmScrolledWindowWidgetClass, cmdAndInfoForm,
		                           XmNspacing,          2,
		                           XmNleftOffset,       10,
					   XmNrightOffset,      8,
					   XmNtopOffset,        10,
       					   XmNbottomOffset,     8,
                                           XmNbottomAttachment, XmATTACH_WIDGET,
				           XmNbottomWidget,     commandText,	   
                                           XmNtopAttachment,    XmATTACH_FORM,     
                                           XmNleftAttachment,   XmATTACH_FORM,
                                           XmNrightAttachment,  XmATTACH_POSITION, 
                                           XmNrightPosition,    65, myNULL);

  commandList    = XtVaCreateManagedWidget("commandList", xmListWidgetClass, cmdListSW,
		                           XmNvisibleItemCount, 4, myNULL);
  
  infoFrame      = XtVaCreateManagedWidget("infoFrame", xmFrameWidgetClass, cmdAndInfoForm,
                   	                   XmNtopOffset,        10,
		                           XmNleftOffset,       5,
		                           XmNrightOffset,      10,
		                           XmNbottomOffset,     10,
                                           XmNtopAttachment,    XmATTACH_FORM,     
		                           XmNbottomAttachment, XmATTACH_FORM,     
                                           XmNrightAttachment,  XmATTACH_FORM,
                                           XmNleftAttachment,   XmATTACH_POSITION, 
		                           XmNleftPosition, 65, myNULL);
 
  infoForm       = XtVaCreateManagedWidget("infoForm", xmFormWidgetClass, infoFrame, myNULL);
  
  XmStr = XmStringCreateSimple("Time = ");
 
  timeLabel1     = XtVaCreateManagedWidget("timeLabel1", xmLabelGadgetClass, infoForm,
		                           XmNalignment,        XmALIGNMENT_END,
		                           XmNlabelString,      XmStr,
                                           XmNleftAttachment,   XmATTACH_FORM,     
		                           XmNrightAttachment,  XmATTACH_POSITION,     
					   XmNrightPosition,    50,
                                           XmNtopAttachment,    XmATTACH_FORM,
				           XmNbottomAttachment, XmATTACH_POSITION,
				           XmNbottomPosition,   26, myNULL);
  XmStringFree(XmStr);
  
  timeLabel2     = XtVaCreateManagedWidget("timeLabel2", xmLabelGadgetClass, infoForm,
		                           XmNalignment,        XmALIGNMENT_BEGINNING,
					   XmNtopOffset,        6,
                                           XmNleftAttachment,   XmATTACH_POSITION,
				           XmNleftPosition,     50,	   
		                           XmNrightAttachment,  XmATTACH_FORM,     
                                           XmNtopAttachment,    XmATTACH_FORM,
				           XmNbottomAttachment, XmATTACH_POSITION,
				           XmNbottomPosition,   26, myNULL);
  
  statusLabel    = XtVaCreateManagedWidget("statusLabel",xmLabelGadgetClass, infoForm,
         		                   XmNleftAttachment,   XmATTACH_FORM,     
	                                   XmNrightAttachment,  XmATTACH_FORM,     
                                           XmNbottomAttachment, XmATTACH_FORM, 
                                           XmNtopAttachment,    XmATTACH_POSITION, 
				           XmNtopPosition,      74, myNULL);

  separator1     = XtVaCreateManagedWidget("separator1",xmSeparatorGadgetClass, infoForm,
         		                   XmNleftAttachment,   XmATTACH_FORM,     
	                                   XmNrightAttachment,  XmATTACH_FORM,     
                                           XmNtopAttachment,    XmATTACH_WIDGET, 
				           XmNtopWidget,        timeLabel1, myNULL);
 
  separator2     = XtVaCreateManagedWidget("separator1",xmSeparatorGadgetClass, infoForm,
         		                   XmNleftAttachment,   XmATTACH_FORM,     
	                                   XmNrightAttachment,  XmATTACH_FORM,     
                                           XmNbottomAttachment, XmATTACH_WIDGET, 
				           XmNbottomWidget,     statusLabel, myNULL);
  
  XmStr = XmStringCreateSimple("Project:");

  projectLabel   = XtVaCreateManagedWidget("Project",xmLabelGadgetClass, infoForm,
		                           XmNlabelString,      XmStr,
         		                   XmNleftAttachment,   XmATTACH_FORM,     
	                                   XmNrightAttachment,  XmATTACH_FORM,     
                                           XmNtopAttachment,    XmATTACH_WIDGET, 
                                           XmNbottomAttachment, XmATTACH_POSITION, 
				           XmNtopWidget,        separator1,
				           XmNbottomPosition,   50, myNULL);
  XmStringFree(XmStr);
  
  projNameLabel  = XtVaCreateManagedWidget("projNameLabel",xmLabelGadgetClass, infoForm,
         		                   XmNleftAttachment,   XmATTACH_FORM,     
	                                   XmNrightAttachment,  XmATTACH_FORM,     
                                           XmNtopAttachment,    XmATTACH_POSITION, 
                                           XmNbottomAttachment, XmATTACH_WIDGET, 
				           XmNbottomWidget,     separator2,
				           XmNtopPosition,      50, myNULL);

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//  MANAGE AND DISPLAY
  
  // manage menu and main window

  XtManageChild(menuBar);

  XtVaSetValues(mainWindow,
                XmNmenuBar,       menuBar,
                XmNworkWindow,    drawingAreaFrame,
                XmNmessageWindow, cmdAndInfoForm, myNULL);

  // flush to screen
  
  XtRealizeWidget(unixGUI.topLevel);

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//  PREPARE FOR DRAWING AND ALLOCATE PIXMAP FOR DRAWINGAREA
  
  // generate graphical context
  
  unixGUI.gc = XCreateGC(XtDisplay(unixGUI.drawingArea), 
		         RootWindowOfScreen(XtScreen(unixGUI.drawingArea)), 0, NULL);
 
  // generate pixmap

//cout << width << '\t' << height << endl;
  
  XtVaGetValues(unixGUI.drawingArea, XmNwidth, &width, XmNheight, &height, myNULL);

  unixGUI.pixmap = XCreatePixmap(XtDisplay(unixGUI.drawingArea), 
		         RootWindowOfScreen(XtScreen(unixGUI.drawingArea)), width, height,
 		         DefaultDepthOfScreen(XtScreen(unixGUI.drawingArea)));
  // wipe it
	 
  XSetForeground(XtDisplay(unixGUI.drawingArea), unixGUI.gc, 
		 BlackPixelOfScreen (XtScreen(unixGUI.drawingArea)));

  XFillRectangle(XtDisplay(unixGUI.drawingArea), unixGUI.pixmap, unixGUI.gc, 0, 0, width, height);  
  
  // get screen dimension and aspect ratio

  unixGUI.width    = XDisplayWidth(XtDisplay(unixGUI.topLevel),0);
  unixGUI.widthMM  = XDisplayWidthMM(XtDisplay(unixGUI.topLevel),0);

  unixGUI.height   = XDisplayHeight(XtDisplay(unixGUI.topLevel),0);
  unixGUI.heightMM = XDisplayHeightMM(XtDisplay(unixGUI.topLevel),0);

  printf("    width of screen  = %5d pixel, %4d mm\n",   unixGUI.width,  unixGUI.widthMM);
  printf("    height of screen = %5d pixel, %4d mm\n\n", unixGUI.height, unixGUI.heightMM);
  
  plot.calcAspectRatio(unixGUI.width, unixGUI.widthMM, unixGUI.height, unixGUI.heightMM);

  cmap = DefaultColormapOfScreen(XtScreen(unixGUI.topLevel));


  // initialise text font

  nFont              = XLoadFont (XtDisplay(unixGUI.drawingArea),"*helvetica-medium-r-normal--12*");
  
  unixGUI.fontStruct = XQueryFont(XtDisplay(unixGUI.drawingArea),nFont);

  XSetFont(XtDisplay(unixGUI.drawingArea),unixGUI.gc,nFont);
  
  
  // initialise named colours
  
  while (namedColours[plot.nNamedColours] != NULL)
  {  
    if (!XAllocNamedColor(XtDisplay(unixGUI.topLevel), 
			  cmap, namedColours[plot.nNamedColours], &col, &unused)) break;
		    
    unixGUI.colourPixelValue[plot.nNamedColours++] = col.pixel;
  }
  if (namedColours[plot.nNamedColours] != NULL)  
  { 
    std::cout << "    Unable to allocate colour '" << namedColours[plot.nNamedColours] << "'\n\n"; 
    exit(0); 
  }

  // initialise contour plot colours

  coloursAllocated = plot.nNamedColours;
  
  for (colour=0; colour<128; colour++)
  {
    col.red   = 0;
    col.green = colour * 512;
    col.blue  = 65535;
    if (XAllocColor(XtDisplay(unixGUI.topLevel), cmap, &col))
      unixGUI.colourPixelValue[coloursAllocated++] = col.pixel;
    else std::cout << "    Unable to allocate colour " << coloursAllocated << "\n\n";
  }

  for (colour=0; colour<128; colour++)
  {
    col.red   = 0;
    col.green = 65535;
    col.blue  = 65535 - colour * 512;
    if (XAllocColor(XtDisplay(unixGUI.topLevel), cmap, &col))
      unixGUI.colourPixelValue[coloursAllocated++] = col.pixel;
    else std::cout << "    Unable to allocate colour " << coloursAllocated << "\n\n";
  }

  for (colour=0; colour<128; colour++)
  {
    col.red   = colour * 512;
    col.green = 65535;
    col.blue  = 0;
    if (XAllocColor(XtDisplay(unixGUI.topLevel), cmap, &col))
      unixGUI.colourPixelValue[coloursAllocated++] = col.pixel;
    else std::cout << "    Unable to allocate colour " << coloursAllocated << "\n\n";
  }

  for (colour=0; colour<128; colour++)
  {
    col.red   = 65535;
    col.green = 65535 - colour * 512;
    col.blue  = 0;
    if (XAllocColor(XtDisplay(unixGUI.topLevel), cmap, &col))
      unixGUI.colourPixelValue[coloursAllocated++] = col.pixel;
    else std::cout << "    Unable to allocate colour " << coloursAllocated << "\n\n";
  }

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//  CALLBACKS

  XtAddCallback(openButton, XmNactivateCallback, grpOpen, myNULL);

  XtAddCallback(closeButton, XmNactivateCallback, grpClose, myNULL); 

  XtAddCallback(exitButton, XmNactivateCallback, grpExit, myNULL);

  XtAddCallback(aboutButton, XmNactivateCallback, grpPopupAbout, myNULL);

  XtAddCallback(helpButton, XmNactivateCallback, grpPopupHelp, myNULL);

  XtAddCallback(unixGUI.drawingArea, XmNinputCallback, grpDrawingAreaMouseInput, myNULL);

  for (i=0; i<macro.n; i++)

    {  XtAddCallback(macroButton[i], XmNactivateCallback, grpPopupMacro, (XtPointer) i);  }

  XtAddCallback(commandText, XmNactivateCallback,        grpCommandEntered,    commandList);
  XtAddCallback(commandList, XmNdefaultActionCallback,   grpCommandEntered,    commandText);
  XtAddCallback(commandList, XmNbrowseSelectionCallback, grpCommandListToText, commandText);
	
  XtAddCallback(unixGUI.drawingArea, XmNexposeCallback, grpDrawingAreaExposed, myNULL);

  XtAddCallback(unixGUI.drawingArea, XmNresizeCallback, grpDrawingAreaResized, myNULL);
  
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//  DISPLAY LABELS

  runCtrl.newStatus(NOPROJECT);
  
  essGrpWriteTime();
 
  runCtrl.newMode(NOPROJ);
  
  essGrpWriteProject();

  std::cout << "    graphical user interface is set-up.\n\n";
  
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//
//  HAND OVER CONTROL TO GUI

  // if required open last project straightaway
  
  if (lastProj) grpOpenProject(); 
  
  // application loop
  
  while (!runCtrl.quit)
  { 	  
    XtAppNextEvent(unixGUI.app, &event);
    XtDispatchEvent(&event);
  }

  XFreePixmap(XtDisplay(unixGUI.topLevel), unixGUI.pixmap);

  // XtDestroyApplicationContext(unixGUI.app); // this sometimes causes a crash,
                                               // but according to the book (V 6a),
					       // it is not necessary anyway
  delete [] macroMenu;
  delete [] macroMenuButton;
  delete [] macroButton;

  return;
}

 


