
#include "headersVTK.h"
//#include "headersBasic.h"

//#include "Global.h"
//#include "UnixGlobal.h"


#include "metis.h"



using namespace std;

//typedef int idx_t ;


typedef double REAL;



int main(int argc, char* argv[])
{
    //////////////////////////////////////////////
    //
    // declare vtk variables
    //
    //////////////////////////////////////////////


        vtkSmartPointer<vtkDataSetMapper>        mapperVTK;
        vtkSmartPointer<vtkActor>                actorVTK;
        vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK;
        vtkSmartPointer<vtkPoints>               pointsVTK;
        vtkSmartPointer<vtkVertex>               vertexVTK;
        vtkSmartPointer<vtkLine>                 lineVTK;
        vtkSmartPointer<vtkQuad>                 quadVTK;
        vtkSmartPointer<vtkHexahedron>           hexVTK;

        vtkSmartPointer<vtkTriangle>             triaVTK;
        vtkSmartPointer<vtkPolygon>              polygonVTK;
        vtkSmartPointer<vtkTetra>                tetraVTK;
        vtkSmartPointer<vtkPyramid>              pyramidVTK;
        vtkSmartPointer<vtkWedge>                wedgeVTK;


        vtkSmartPointer<vtkIntArray>          nodeInOutVTK;
        vtkSmartPointer<vtkIntArray>          procIdVTK;
        vtkSmartPointer<vtkIntArray>          cellOrientationVTK;
        //vtkSmartPointer<vtkDoubleArray>          vectors, vectors2, scalars, scalars2;
        vtkSmartPointer<vtkFloatArray>          vecVTK, vecVTK2, distVTK, scaVTK2, cellDataVTK, cellDataVTK2;
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK;

        vtkSmartPointer<vtkExtractEdges>     extractEdgesVTK;

    mapperVTK    =  vtkSmartPointer<vtkDataSetMapper>::New();
    actorVTK     =  vtkSmartPointer<vtkActor>::New();
    uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New(); 
    pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    lineVTK      =  vtkSmartPointer<vtkLine>::New();
    quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    hexVTK       =  vtkSmartPointer<vtkHexahedron>::New();
    vertexVTK    =  vtkSmartPointer<vtkVertex>::New();

    triaVTK      =  vtkSmartPointer<vtkTriangle>::New();
    polygonVTK   =  vtkSmartPointer<vtkPolygon>::New();
    tetraVTK     =  vtkSmartPointer<vtkTetra>::New();
    pyramidVTK   =  vtkSmartPointer<vtkPyramid>::New();
    wedgeVTK     =  vtkSmartPointer<vtkWedge>::New();

    vecVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    distVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    scaVTK2      =  vtkSmartPointer<vtkFloatArray>::New();
    vecVTK2      =  vtkSmartPointer<vtkFloatArray>::New();
    cellDataVTK  =  vtkSmartPointer<vtkFloatArray>::New();
    cellDataVTK2 =  vtkSmartPointer<vtkFloatArray>::New();
  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    nodeInOutVTK       =  vtkSmartPointer<vtkIntArray>::New();
    procIdVTK       =  vtkSmartPointer<vtkIntArray>::New();
  cellOrientationVTK    =  vtkSmartPointer<vtkIntArray>::New();

     extractEdgesVTK  =    vtkSmartPointer<vtkExtractEdges>::New();


    //vtkSmartPointer<vtkXMLPolyDataWriter> writer =  vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    // Add the polygon to a list of polygons
    vtkSmartPointer<vtkCellArray> polyList  =   vtkSmartPointer<vtkCellArray>::New();

    // Create a PolyData
    vtkSmartPointer<vtkPolyData> polyData =   vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPolyData> polyData2 =   vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkXMLPolyDataWriter>  writerPolyData =     vtkSmartPointer<vtkXMLPolyDataWriter>::New();

    //////////////////////////////////////////////
    //
    // generate base triangulation
    //
    //////////////////////////////////////////////

    double  x0, x1, y0, y1;

    int nEx, nEy, nNx, nNy, nNode, nElem, npElem;

    int  ee, ii, jj, kk, ind, ind1, ind2, ind3, n1, n2, n3, n4, e1;
    idx_t nParts  = 2;

    x0 = -1.6;    x1 =  1.6;    nEx = 5;
    y0 = -1.6;    y1 =  1.6;    nEy = 5;

    x0 = 0.0;    x1 =  1.0;    nEx = 5;
    y0 = 0.0;    y1 =  1.0;    nEy = 5;

//
    if(argc == 1)
    {
      cout << " Enter Input data " << endl;
      exit(0);
    }
    else
    {
       x0  = atof(argv[1]);
       x1  = atof(argv[2]);
       nEx = atoi(argv[3]);

       y0  = atof(argv[4]);
       y1  = atof(argv[5]);
       nEy = atoi(argv[6]);
       nParts    = atoi(argv[7]);
    }
//
      npElem = 3;
      nElem = nEx*nEy*2;

    nNx = nEx+1;
    nNy = nEy+1;

    nNode = nNx*nNy;


    uGridVTK->Reset();
    pointsVTK->Reset();
    cellDataVTK->Reset();
    cellDataVTK->Reset();
    procIdVTK->Reset();
    nodeInOutVTK->Reset();



    REAL dx = (x1-x0)/nEx;
    REAL dy = (y1-y0)/nEy;
    REAL  xx, yy, fact;

    cout << x0 << '\t' << x1 << endl;
    cout << y0 << '\t' << y1 << endl;

    cout << " nEx    = " << nEx << endl;
    cout << " nEy    = " << nEy << endl;
    cout << " nNx    = " << nNx << endl;
    cout << " nNy    = " << nNy << endl;
    cout << " npElem = " << npElem << endl;
    cout << " nNode  = " << nNode << endl;
    cout << " nElem  = " << nElem << endl;



    vtkIdType pt[10];

    // variables for METIS

    idx_t nWeights  = 1;

    idx_t objval;
    idx_t *xadj, *adjncy;
    idx_t numflag=0;
    idx_t ncommon=2;

    idx_t  eptr[nElem+1];
    idx_t  eind[nElem*npElem];
    idx_t  epart[nElem];
    idx_t  npart[nNode];


    procIdVTK->SetName("procId");
    procIdVTK->SetNumberOfTuples(nElem);

    //////////////////////////////////////////////
    //
    // create nodes/points
    //
    //////////////////////////////////////////////

    cout << " Creating nodes " << endl;

    yy = y0;
    for(jj=0; jj<nNy; jj++)
    {
      xx = x0;
      for(ii=0; ii<nNx; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(xx, yy, 0.0);

        xx += dx;
      }
      yy += dy;
    }

    //////////////////////////////////////////////
    //
    // create cells/elements
    //
    //////////////////////////////////////////////

    cout << " Creating elements " << endl;


    for(ii=0; ii<=nElem; ii++)
      eptr[ii] = npElem*ii;



    ind=0;
    for(jj=0; jj<nEy; jj++)
    {
      ind2 = nNx*jj;
      ind3 = nNx*(jj+1);

      for(ii=0; ii<nEx; ii++)
      {
        n1 = ind2 + ii;
        n2 = n1+1;
        n3 = ind3 + ii;
        n4 = n3+1;

        polygonVTK->GetPointIds()->SetNumberOfIds(3); //make a triangle

        polygonVTK->GetPointIds()->SetId(0, n1);
        polygonVTK->GetPointIds()->SetId(1, n2);
        polygonVTK->GetPointIds()->SetId(2, n4);

        polyList->InsertNextCell(polygonVTK);

        kk = ind*npElem;

        eind[kk+0] = n1;
        eind[kk+1] = n2;
        eind[kk+2] = n4;

        ind++;

        polygonVTK->GetPointIds()->SetNumberOfIds(3); //make a triangle

        polygonVTK->GetPointIds()->SetId(0, n1);
        polygonVTK->GetPointIds()->SetId(1, n4);
        polygonVTK->GetPointIds()->SetId(2, n3);

        polyList->InsertNextCell(polygonVTK);

        kk = ind*npElem;

        eind[kk+0] = n1;
        eind[kk+1] = n4;
        eind[kk+2] = n3;

        ind++;
      }
    }


      for(e1=0; e1<nElem; e1++)
      {
        kk = e1*npElem;

        cout << e1 << '\t' ;
        for(ii=0; ii<npElem; ii++)
          cout << eind[kk+ii] << '\t' ;

        cout << endl;

        ind++;
      }

      cout << " \n\n\n\n " << endl;
    
    cout << " partitioning the mesh " << endl;


    idx_t options[METIS_NOPTIONS];

    METIS_SetDefaultOptions(options);

    //options[METIS_OPTION_NSEPS] = 10;

    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;

    options[METIS_OPTION_NUMBERING] = 0;


    //if(nParts > 1)
    //{

        cout <<  " BBBBBBBBBBB " << nElem << '\t' << nNode << endl;

        // METIS partition routine
        int ret = METIS_PartMeshNodal(&nElem, &nNode, eptr, eind, NULL, NULL, 
     				       &nParts, NULL, options, &objval, epart, npart);

        //int ret = METIS_PartMeshDual(&nElem, &nNode, epart, eind, NULL, NULL, 
	//			      &ncommon, &nParts, NULL, options, &objval, epart, npart);

       if(ret == METIS_OK)
         std::cout << " METIS partition routine successful "  << std::endl;
       else
         std::cout << " METIS partition routine FAILED "  << std::endl;

       for(jj=0; jj<nElem; jj++)
         procIdVTK->InsertTuple1(jj, epart[jj]);
    //}
    //else
    //{
      //for(jj=0; jj<nElem; jj++)
      //   procIdVTK->InsertTuple1(jj, 0);
    //}


    //////////////////////////////////////////////
    //
    // setup and write polyData
    //
    //////////////////////////////////////////////

    cout << " Writing polyData  " << endl;

    char fname1[100];

    sprintf(fname1,"%s%d%s%d%s%d%s", "mesh-tria-",nEx,"-",nEy,"-",nParts,".vtp");


    polyData->SetPoints(pointsVTK);
    polyData->SetPolys(polyList);
    //polyData->GetPointData()->SetScalars(distVTK);
    //polyData->GetPointData()->AddArray(nodeInOutVTK);
    //polyData->GetCellData()->AddArray(cellOrientationVTK);
    polyData->GetCellData()->AddArray(procIdVTK);

    polyData->BuildLinks();
    polyData->BuildCells();

    //Write the file.
    writerPolyData->SetFileName(fname1);
    writerPolyData->SetInput(polyData);
    writerPolyData->Write();

    cout << " polyData->GetNumberOfPoints() = " <<  polyData->GetNumberOfPoints() << endl;
    cout << " polyData->GetNumberOfCells()  = " <<  polyData->GetNumberOfCells() << endl;
    cout << " polyData->GetNumberOfVerts()  = " <<  polyData->GetNumberOfVerts() << endl;
    cout << " polyData->GetNumberOfLines()  = " <<  polyData->GetNumberOfLines() << endl;
    cout << " polyData->GetNumberOfPolys()  = " <<  polyData->GetNumberOfPolys() << endl;
    cout << " polyData->GetNumberOfStrips() = " <<  polyData->GetNumberOfStrips() << endl;




  return 0;
}
















