
#ifndef incl_PlotVTK_h
#define incl_PlotVTK_h

#include "vtkSmartPointer.h"
#include "vtkActor.h"
#include "vtkConeSource.h"
#include "vtkGlyph3D.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkSphereSource.h"
#include "vtkXOpenGLRenderWindow.h"
#include "vtkXRenderWindowInteractor.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPoints.h"
#include <vtkAxesActor.h>


#include "MyString.h"
#include "Domain.h"



class PlotVTK
{
  public:
  
    bool  ActiveFlag;
/*
    vtkXOpenGLRenderWindow   *renWindow;

    vtkXRenderWindowInteractor  *renIntr;

    vtkRenderer  *rendr;

    vtkUnstructuredGrid    *uGrid ;
*/
    vtkSmartPointer<vtkXOpenGLRenderWindow>   renWindow;

    vtkSmartPointer<vtkXRenderWindowInteractor>  renIntr;

    vtkSmartPointer<vtkRenderer>  rendr;

    vtkSmartPointer<vtkUnstructuredGrid>    uGrid ;

    vtkSmartPointer<vtkPoints>             points;

   vtkSmartPointer<vtkAxesActor> axes;


    PlotVTK();

    ~PlotVTK();
    
    void  write2file(MyString &fileName);
    
    void  VtkGL2PSExporter(MyString &fileName, int);
    
    void  writeImages(MyString &fileName, int);
    
    void  reset();
    
    void  set();
    
    void  setBackgroundColor(double* );
    
    void  clearWindow();

    void  addAxes(double,double,double);
    
    void  removeAxes();

    void  VtkReflect(int);

};

#endif

