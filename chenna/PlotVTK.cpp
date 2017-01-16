//#include "Plot.h"
#include "PlotVTK.h"
#include "UnixGUI.h"
#include "Files.h"
#include "NurbsUtilities.h"


#include "vtkCamera.h"
#include "vtkQuadric.h"
#include "vtkSampleFunction.h"
#include "vtkContourFilter.h"
#include "vtkOutlineFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkImageData.h" 
#include "vtkAxesActor.h"
#include "vtkActor2DCollection.h"
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkGL2PSExporter.h>
#include <vtkWindowToImageFilter.h>
#include <vtkAxesActor.h>
#include <vtkReflectionFilter.h>
#include <vtkDataSetMapper.h>
#include <vtkLookupTable.h>
#include <vtkMapper.h>
#include <vtkObject.h>
#include "vtkSmartPointer.h"

#include <vtkImageData.h>
#include <vtkBMPWriter.h>
#include <vtkPNGWriter.h>
#include <vtkJPEGWriter.h>
#include <vtkTIFFWriter.h>
#include <vtkPNMWriter.h>
#include <vtkPostScriptWriter.h>

//#include <vtkImageMapper3D.h>
#include <vtkImageCanvasSource2D.h>
#include <vtkImageActor.h>
#include <vtkInteractorStyleImage.h>

extern Files    files;
extern UnixGUI    unixGUI;
//extern Plot plot;

using namespace std;



PlotVTK::PlotVTK()
{
   ActiveFlag = false;

   renWindow  =  vtkSmartPointer<vtkXOpenGLRenderWindow>::New();
   renIntr    =  vtkSmartPointer<vtkXRenderWindowInteractor>::New();
   uGrid      =  vtkSmartPointer<vtkUnstructuredGrid>::New();
   rendr      =  vtkSmartPointer<vtkRenderer>::New();
   points     =  vtkSmartPointer<vtkPoints>::New();
   axes       = vtkSmartPointer<vtkAxesActor>::New();

   renIntr->SetRenderWindow(renWindow);
   renWindow->AddRenderer(rendr);


}


PlotVTK::~PlotVTK()
{
  // no need for the following commands as 'vtkSmartPointer's are used 
  //renWindow->Delete();
  //renIntr->Delete();
  //rendr->Delete();
  //uGrid->Delete();

    clearWindow();

}




void PlotVTK::reset()
{
    ActiveFlag = false;
     
    clearWindow();

    renIntr->Disable();
  
    //plot.wipe();

  return;
}




void PlotVTK::set()
{
     //cout << "  PlotVTK::set() ... " << ActiveFlag << endl;

  if(ActiveFlag)  // deactivate PlotVTK and transfer the control to basic Xt
  {
     reset();
    
     cout << "  VTK is deactivated and GUI control transfered to basic Xt ... " << endl;
     cout << endl;
     return;
  }
  else
  {
     ActiveFlag = true;

     cout << "  VTK is active now ... " << endl;     cout << endl;

     Display *display;

     // get the display connection and give it to the renderer
     display = XtDisplay(unixGUI.topLevel);
     renWindow->SetDisplayId(display);

     rendr->SetBackground(0.4,0.1,0.2);

     rendr->SetBackground(.2,0.5,.1);
     
     //rendr->SetViewport(0.0, 0.0, 0.5, 1.0);

     renIntr->SetWidget(unixGUI.drawingArea);
     renIntr->Initialize(unixGUI.app);

     rendr->ResetCamera();
     rendr->Clear();

     //rendr->GetActiveCamera()->SetViewPlaneNormal(0,0,1);

     rendr->GetActiveCamera()->ParallelProjectionOn();

     renWindow->Render();
     
     renIntr->Enable();

//     XtAppMainLoop(unixGUI.app);

     return;
  }
}



void PlotVTK::setBackgroundColor(double* color)
{
    rendr->SetBackground(color[0], color[1], color[2]);

    renWindow->Render();

   return;
}



void PlotVTK::clearWindow()
{
    vtkActorCollection  *actorCollection = rendr->GetActors();
    actorCollection->InitTraversal();
  
    while(actorCollection->GetNumberOfItems() != 0)
    {
      vtkActor* nextActor = actorCollection->GetNextActor();
      rendr->RemoveActor(nextActor);
    }
    //cout << " actorCollection->GetNumberOfItems() " << actorCollection->GetNumberOfItems() << endl;

    vtkActor2DCollection  *actor2DCollection = rendr->GetActors2D();
    actor2DCollection->InitTraversal();
  
    while(actor2DCollection->GetNumberOfItems() != 0)
    {
      vtkActor2D* nextActor2D = actor2DCollection->GetNextActor2D();
      rendr->RemoveActor2D(nextActor2D);
    }
    
    points->Reset();
    uGrid->Reset();
    rendr->ResetCamera();
    rendr->Clear();

     rendr->ResetCamera();

     renWindow->Render();

//     renIntr->Enable();

   return;
}


void  PlotVTK::write2file(MyString &fileName)
{

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(fileName.asCharArray());
    //writer->SetInputData(uGrid);
    writer->SetInput(uGrid);
    writer->Write();


   return;
}



void  PlotVTK::VtkGL2PSExporter(MyString &filename, int index)
{
  //PS_FILE=0; EPS_FILE=1; PDF_FILE=2; TEX_FILE=3; SVG_FILE=4;
  
  vtkSmartPointer<vtkGL2PSExporter>  exporter = vtkSmartPointer<vtkGL2PSExporter>::New();

  MyString pathAndFile;
  
  char tmp[20];
 
  sprintf(tmp,".%04d",++files.nps);
  
  pathAndFile.append(files.projDir).append(SLASH).append(files.Pfile).append(tmp);

  while (prgFileExist(pathAndFile)) pathAndFile.trunc(pathAndFile.length()-4).append(".new");
  
//  COUT << &((pathAndFile.asCharArray())[files.projDir.length()+1]) << " open for plot output.\n\n";

  exporter->SetRenderWindow(renWindow);
  exporter->SetFilePrefix(pathAndFile.asCharArray());
  exporter->SetFileFormat(index);
  exporter->DrawBackgroundOff();
  exporter->CompressOn();
  exporter->SetSortToBSP();

  exporter->Update();
  exporter->Write();  

   return;
}




void  PlotVTK::writeImages(MyString &filename, int index)
{
   vtkSmartPointer<vtkWindowToImageFilter>  ImageFilter =     vtkSmartPointer<vtkWindowToImageFilter>::New();

   ImageFilter->SetInput(renWindow);

   vtkSmartPointer<vtkImageWriter>   imageWriter;


   if(index == 0) // PNG
   {
       imageWriter  =  vtkSmartPointer<vtkPNGWriter>::New();
       filename.append(".png");
   }
   if(index == 1) // JPEG
   {
       imageWriter  =  vtkSmartPointer<vtkJPEGWriter>::New();
       filename.append(".jpeg");
   }
   if(index == 2) // TIFF
   {
       imageWriter  =  vtkSmartPointer<vtkTIFFWriter>::New();
       filename.append(".tiff");
   }
   if(index == 3) // BMP
   {
       imageWriter  =  vtkSmartPointer<vtkBMPWriter>::New();
       filename.append(".bmp");
   }
   if(index == 4) // PNM
   {
       imageWriter  =  vtkSmartPointer<vtkPNMWriter>::New();
       filename.append(".pnm");
   }
   if(index == 5) // PS
   {
       imageWriter  =  vtkSmartPointer<vtkPostScriptWriter>::New();
       filename.append(".ps");
   }

   //imageWriter->SetInputData(ImageFilter->GetOutput());
   imageWriter->SetInput(ImageFilter->GetOutput());
   imageWriter->SetFileName(filename.asCharArray());
   imageWriter->Write();

   return;
}






void PlotVTK::addAxes(double val1, double val2, double val3)
{
   axes->SetTotalLength(val1,val2,val3);
   rendr->AddActor(axes);

   return;
}


void PlotVTK::removeAxes()
{

   rendr->RemoveActor(axes);

   return;
}




void PlotVTK::VtkReflect(int index)
{
   vtkSmartPointer<vtkReflectionFilter>  reflectionFilter =   vtkSmartPointer<vtkReflectionFilter>::New();
   vtkSmartPointer<vtkLookupTable>              hueLut    =   vtkSmartPointer<vtkLookupTable>::New();
   vtkSmartPointer<vtkDataSetMapper>     reflectionMapper =   vtkSmartPointer<vtkDataSetMapper>::New();
   vtkSmartPointer<vtkActor>              reflectionActor =   vtkSmartPointer<vtkActor>::New();
   vtkSmartPointer<vtkActorCollection>   actorCollection  =   vtkSmartPointer<vtkActorCollection>::New();
   vtkSmartPointer<vtkActor>              nextActor =   vtkSmartPointer<vtkActor>::New();

    double  *range;
   
    range = uGrid->GetScalarRange();

    hueLut->SetHueRange(0.66667, 0.0);
    hueLut->SetTableRange(range[0], range[1]);
    hueLut->SetNumberOfColors(10);
    hueLut->SetRampToLinear();
    hueLut->SetScaleToLinear();
    hueLut->Build();



    actorCollection = rendr->GetActors();
    actorCollection->InitTraversal();
    
    vtkIndent  indt;
    
    int  nn = actorCollection->GetNumberOfItems(), ii ;
    
    //cout <<  " actorCollection->GetNumberOfItems() " << nn << endl;

    reflectionMapper->SetScalarModeToUsePointData();
    reflectionMapper->InterpolateScalarsBeforeMappingOn();
    reflectionMapper->ScalarVisibilityOn();
    reflectionMapper->SetColorModeToMapScalars();
    reflectionMapper->SetScalarRange(range[0], range[1]);
    reflectionMapper->SetLookupTable(hueLut);


    for(ii=0;ii<nn;ii++)
    {
       nextActor = actorCollection->GetNextActor();

       reflectionFilter->SetInputConnection(nextActor->GetMapper()->GetInput()->GetProducerPort());
       //reflectionFilter->SetInputData(nextActor->GetMapper()->GetInput());
       
       reflectionFilter->CopyInputOn();

       reflectionFilter->SetPlane(index);

       reflectionFilter->Update();

       reflectionMapper->SetInputConnection(reflectionFilter->GetOutputPort());

    }

    reflectionActor->SetMapper(reflectionMapper);
    rendr->AddActor(reflectionActor);
    renWindow->Render();

   return;
}








