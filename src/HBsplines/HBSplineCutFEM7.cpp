
#include "HBSplineCutFEM.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "Functions.h"
#include "Files.h"
#include "MyString.h"
#include "myGeomUtilities.h"
#include "myTria.h"
#include "DistFunctions.h"
#include "AABB.h"

#include "vtkXMLPUnstructuredGridWriter.h"

extern ComputerTime       computerTime;
extern MpapTime mpapTime;
extern Files files;

using namespace myGeom;



void HBSplineCutFEM::plotGeom(int val1, bool flag2, int col, bool PLOT_KNOT_LINES, int* resln)
{
    PetscPrintf(MPI_COMM_WORLD, "     HBSplineCutFEM: plotgeometry ...\n\n");

    double tstart = MPI_Wtime();

    uGridVTK->Reset();
    pointsVTK->Reset();

    scaVTK->Reset();
    scaVTK2->Reset();

    vecVTK->Reset();
    vecVTK2->Reset();

    cellDataVTK->Reset();
    cellDataVTK2->Reset();

    int ii=0;
    while(uGridVTK->GetPointData()->GetNumberOfArrays() != 0)
      uGridVTK->GetPointData()->RemoveArray(ii++);

    ii=0;
    while(uGridVTK->GetCellData()->GetNumberOfArrays() != 0)
      uGridVTK->GetCellData()->RemoveArray(ii++);


    cellDataVTK->SetNumberOfComponents(1);
    cellDataVTK2->SetNumberOfComponents(1);
    //vecVTK->SetNumberOfTuples(count);

    if(CUTCELL_INTEGRATION_TYPE == 1)
    {
      if(ndm == 1)
        plotGeomSubTrias1D(val1, flag2, col, PLOT_KNOT_LINES, resln);
      else if(ndm == 2)
        plotGeomSubTrias2D(val1, flag2, col, PLOT_KNOT_LINES, resln);
      else
        plotGeomSubTrias3D(val1, flag2, col, PLOT_KNOT_LINES, resln);
    }
    else
    {
      if(ndm == 1)
        plotGeomAdapIntegration1D(val1, flag2, col, PLOT_KNOT_LINES, resln);
      else if(ndm == 2)
        plotGeomAdapIntegration2D(val1, flag2, col, PLOT_KNOT_LINES, resln);
      else
        plotGeomAdapIntegration3D(val1, flag2, col, PLOT_KNOT_LINES, resln);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    uGridVTK->SetPoints(pointsVTK);

    cellDataVTK->SetName("ElemType");
    cellDataVTK2->SetName("SubdomId");

    uGridVTK->GetCellData()->SetScalars(cellDataVTK);
    uGridVTK->GetCellData()->AddArray(cellDataVTK2);

    uGridVTK->Squeeze();

    MPI_Barrier(MPI_COMM_WORLD);

    /////////////////////
    // write VTK files

    // write parallel vtu master (pvtu) file

    if(DIM == 3)
    {
      char fnameMaster[500];

      cout << "files.projDir = " << files.projDir << endl;
      //sprintf(fnameMaster,"%s%s", files.Ofile.asCharArray(),"-geom.pvtu");
      sprintf(fnameMaster,"%s%s%s%s", files.projDir.asCharArray(), "/", files.Ofile.asCharArray(),"-geom.pvtu");

      vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writeruGridP  =  vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();

      writeruGridP->SetFileName(fnameMaster);
      writeruGridP->SetInputData(uGridVTK);
      writeruGridP->SetNumberOfPieces(n_mpi_procs);
      writeruGridP->SetStartPiece(0);
      writeruGridP->SetEndPiece(n_mpi_procs-1);
      writeruGridP->Write();

      MPI_Barrier(MPI_COMM_WORLD);
    }

    // write individual vtu files
    char fnameLocal[500];

    //sprintf(fnameLocal,"%s%s%d%s", files.Ofile.asCharArray(),"-geom_",this_mpi_proc,".vtu");
    //sprintf(fnameLocal,"%s%s%s%s%d%s", files.projDir.asCharArray(), "/", files.Ofile.asCharArray(),"-geom_",this_mpi_proc,".vtu");
    sprintf(fnameLocal,"%s%s%s%s%s", files.projDir.asCharArray(), "/", files.Ofile.asCharArray(),"-geom",".vtu");

    writerUGridVTK->SetFileName(fnameLocal);
    writerUGridVTK->SetNumberOfPieces(1);
    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    //plotGaussPointsElement();
    //plotGaussPointsDirichletBoundary();
    //plotGaussPointsNeumannBoundary();

    double tend = MPI_Wtime();
    PetscPrintf(MPI_COMM_WORLD, "\n HBSplineCutFEM::plotGeom() took %f millisecond(s) \n ", (tend-tstart)*1000);

    return;
}



void HBSplineCutFEM::plotGeomSubTrias1D(int val1, bool flag2, int col, bool PLOT_KNOT_LINES, int* resln)
{ 
    int  ii=0, jj=0, n1=0, n2=0, ind1=0, ind2=0, ll=0;

    vtkIdType pt0, pt1, pt2;
    myPoint   knotBegin, knotEnd;

    for(ii=0;ii<elems.size();ii++)
    {
      //if( elems[ii]->isLeaf() && !(elems[ii]->isGhost()) &&  elems[ii]->isActive())
      if( elems[ii]->isActive() )
      {
          knotBegin = elems[ii]->getKnotBegin();
          knotEnd   = elems[ii]->getKnotEnd();

          param[0] = knotBegin[0];
          computeGeometry(param, geom);
          pt0 = pointsVTK->InsertNextPoint(geom[0], 0.0, 0.0);

          param[0] = knotEnd[0];
          computeGeometry(param, geom);
          pt1 = pointsVTK->InsertNextPoint(geom[0], 0.0, 0.0);

          vertexVTK->GetPointIds()->SetId(0, pt0);
          uGridVTK->InsertNextCell(vertexVTK->GetCellType(), vertexVTK->GetPointIds());

          vertexVTK->GetPointIds()->SetId(0, pt1);
          uGridVTK->InsertNextCell(vertexVTK->GetCellType(), vertexVTK->GetPointIds());

          lineVTK->GetPointIds()->SetId(0,pt0);
          lineVTK->GetPointIds()->SetId(1,pt1);
          uGridVTK->InsertNextCell(lineVTK->GetCellType(), lineVTK->GetPointIds());
       }
    }

    return;
}




void HBSplineCutFEM::plotGeomSubTrias2D(int val1, bool flag2, int col, bool PLOT_KNOT_LINES, int* resln)
{ 
    int  ee=0, ii=0, kk=0, ll=0, typetemp=0, totalNGP=0;

    vtkIdType  ptIds[20],  ptId, cellId;
    node *ndTemp;
    myPoint  ptTemp;

    AABB  bbTemp;

    for(ee=0;ee<activeElements.size();ee++)
    {
      ndTemp = elems[activeElements[ee]];

      if( ndTemp->getSubdomainId() == this_mpi_proc )
      {
        if( !ndTemp->isCutElement() )
        {
          bbTemp = ndTemp->getAABB();

          ptIds[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], 0.0);
          ptIds[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], 0.0);
          ptIds[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], 0.0);
          ptIds[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], 0.0);

          for(ll=0;ll<4;ll++)
            quadVTK->GetPointIds()->SetId(ll, ptIds[ll]);

          cellDataVTK->InsertNextValue(ndTemp->getDomainNumber());
          cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        }
        else // the element is cutCell
        {
          totalNGP += ndTemp->Quadrature.gausspoints.size();

          vtkSmartPointer<vtkTriangle> triaVTK =  vtkSmartPointer<vtkTriangle>::New();

          myPoly *poly;

          for(ii=0; ii<ndTemp->subTrias.size(); ii++)
          {
            poly = ndTemp->subTrias[ii];

            for(kk=0; kk<3; kk++)
            {
              ptTemp = poly->GetPoint(kk);

              ptId = pointsVTK->InsertNextPoint(ptTemp[0], ptTemp[1], 0.0);

              triaVTK->GetPointIds()->SetId(kk, ptId );
            }

            cellDataVTK->InsertNextValue(poly->getDomainNumber());
            cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());

            uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());

          } //  for(ii=0; ii<nd->subTrias.size(); ii++)
        } //else
      }
    } // for(ee=0;ee<elems.size();ee++)

    return;
}




void HBSplineCutFEM::plotGeomSubTrias3D(int val1, bool flag2, int col, bool PLOT_KNOT_LINES, int* resln)
{ 
    int  ee=0, ii=0, kk=0, ll=0, typetemp=0, totalNGP=0;

    vtkIdType  ptIds[20],  ptId, cellId;
    node *ndTemp;
    myPoint  ptTemp;

    AABB  bbTemp;

    for(ee=0; ee<activeElements.size(); ee++)
    {
      ndTemp = elems[activeElements[ee]];

      if( ndTemp->getSubdomainId() == this_mpi_proc )
      {
        if( !ndTemp->isCutElement() )
        {
          bbTemp = ndTemp->getAABB();

          ptIds[0] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.minBB[2]);
          ptIds[1] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.minBB[2]);
          ptIds[2] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]);
          ptIds[3] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.minBB[2]);

          ptIds[4] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]);
          ptIds[5] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.minBB[1], bbTemp.maxBB[2]);
          ptIds[6] = pointsVTK->InsertNextPoint(bbTemp.maxBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]);
          ptIds[7] = pointsVTK->InsertNextPoint(bbTemp.minBB[0], bbTemp.maxBB[1], bbTemp.maxBB[2]);

          for(ll=0;ll<8;ll++)
            hexVTK->GetPointIds()->SetId(ll, ptIds[ll]);

          cellDataVTK->InsertNextValue( ndTemp->getDomainNumber() );
          cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());

          uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
        }
        else // the element is cutCell
        {
          totalNGP += ndTemp->Quadrature.gausspoints.size();

          vtkSmartPointer<vtkTetra> tetVTK =  vtkSmartPointer<vtkTetra>::New();

          myPoly *poly;

          for(ii=0; ii<ndTemp->subTrias.size(); ii++)
          {
            poly = ndTemp->subTrias[ii];

            for(kk=0; kk<4; kk++)
            {
              ptTemp = poly->GetPoint(kk);

              ptId = pointsVTK->InsertNextPoint(ptTemp[0], ptTemp[1], ptTemp[2]);

              tetVTK->GetPointIds()->SetId(kk, ptId );
            }

            cellDataVTK->InsertNextValue(poly->getDomainNumber());
            cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());

            uGridVTK->InsertNextCell(tetVTK->GetCellType(), tetVTK->GetPointIds());

          } //  for(ii=0; ii<nd->subTrias.size(); ii++)
        } //else
      }
    } // for(ee=0;ee<elems.size();ee++)

    return;
}



void  HBSplineCutFEM::postProcessFlow(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    if( (filecount % nCol) !=  0)
        return;

    if(this_mpi_proc != 0)
        return;


    double tstart = MPI_Wtime();

    uGridVTK->Reset();
    pointsVTK->Reset();

    scaVTK->Reset();
    scaVTK2->Reset();

    vecVTK->Reset();
    vecVTK2->Reset();

    cellDataVTK->Reset();
    cellDataVTK2->Reset();

    for(int ii=0; ii<uGridVTK->GetPointData()->GetNumberOfArrays(); ii++)
      uGridVTK->GetPointData()->RemoveArray(ii);

    for(int ii=0; ii<uGridVTK->GetCellData()->GetNumberOfArrays(); ii++)
      uGridVTK->GetCellData()->RemoveArray(ii);

    if(CUTCELL_INTEGRATION_TYPE == 1)
    {
      if(DIM == 1)
        postProcessSubTrias1D(vartype, vardir, nCol, umnxflag, umin, umax, resln);
      else if(DIM == 2)
        postProcessSubTrias2D(vartype, vardir, nCol, umnxflag, umin, umax, resln);
      else
        postProcessSubTrias3D(vartype, vardir, nCol, umnxflag, umin, umax, resln);
    }
    else
    {
      if(DIM == 1)
        postProcessAdapIntegration1D(vartype, vardir, nCol, umnxflag, umin, umax, resln);
      else if(DIM == 2)
        postProcessAdapIntegration2D(vartype, vardir, nCol, umnxflag, umin, umax, resln);
      else
        postProcessAdapIntegration3D(vartype, vardir, nCol, umnxflag, umin, umax, resln);
    }

    cellDataVTK2->SetName("SubdomId");

    //MPI_Barrier(MPI_COMM_WORLD);

    /////////////////////
    // write VTK files

    // write parallel vtu master (pvtu) file

    if(DIM == 3)
    {
      char fnameMaster[500];

      //sprintf(fnameMaster,"%s%s%06d%s", files.Ofile.asCharArray(),"-",filecount,".pvtu");
      sprintf(fnameMaster,"%s%s%s%s%06d%s", files.projDir.asCharArray(), "/", files.Ofile.asCharArray(),"-",filecount,".pvtu");

      vtkSmartPointer<vtkXMLPUnstructuredGridWriter> writeruGridP  =  vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();

      writeruGridP->SetFileName(fnameMaster);
      writeruGridP->SetInputData(uGridVTK);
      writeruGridP->SetNumberOfPieces(n_mpi_procs);
      writeruGridP->SetStartPiece(0);
      writeruGridP->SetEndPiece(n_mpi_procs-1);
      writeruGridP->Write();

      MPI_Barrier(MPI_COMM_WORLD);
    }

    // write individual vtu files

    char fnameLocal[500];

    //sprintf(fnameLocal,"%s%s%06d%s", files.Ofile.asCharArray(),"-",filecount,".vtu");
    sprintf(fnameLocal,"%s%s%s%s%06d%s", files.projDir.asCharArray(), "/", files.Ofile.asCharArray(),"-",filecount,".vtu"); 

    //sprintf(fnameLocal,"%s%s%06d%s%d%s", files.Ofile.asCharArray(),"-",filecount,"_",this_mpi_proc,".vtu");

    writerUGridVTK->SetFileName(fnameLocal);
    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
      ImmersedBodyObjects[bb]->postProcess(filecount);
    }

    double tend = MPI_Wtime(); 
    PetscPrintf(MPI_COMM_WORLD, "HBSplineCutFEM::postProcessFlow() took %f millisecond(s) \n ", (tend-tstart)*1000);

    return;
}



void  HBSplineCutFEM::postProcessSubTrias1D(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
  return;
}





void  HBSplineCutFEM::postProcessSubTrias2D(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    int  dd=0, ii=0, jj=0, kk=0, ll=0, count=0, index=0;
    int  ind1=0, ind2=0, e=0, ee=0, gcount=0, ind=0, domTemp=0;
    int nlocal = (degree[0]+1) * (degree[1] + 1);

    VectorXd  NN(nlocal), N(nlocal), dN_dx(nlocal), dN_dy(nlocal), dNN_dx(nlocal), dNN_dy(nlocal), tempVec, tempVec2, d2N_dx2(nlocal), d2N_dy2(nlocal);
    VectorXd  vectmp(nlocal), rhsTemp;
    myPoint  knotIncr, knotBegin, knotEnd;

    double   fact=0.0, incr1=0.0, incr2=0.0, val1=0.0;

    vector<double>  uu, vv;

    node* ndTemp;

    vtkIdType pt[50], cellId;

    if(ndf == 1)
    {
      index = 0;
      for(ee=0; ee<activeElements.size(); ee++)
      {
          ndTemp = elems[activeElements[ee]];

          knotBegin = ndTemp->getKnotBegin();
          knotEnd   = ndTemp->getKnotEnd();
          knotIncr  = ndTemp->getKnotIncrement();

          if( !(ndTemp->isCutElement()) )
          //if( ndTemp->getDomainNumber() < 10 )
          {
            //if( ndTemp->getDomainNumber() == 0 )
            //{
              fact = knotIncr[0]/resln[0];
              create_vector(knotBegin[0], knotEnd[0], fact, uu);

              fact = knotIncr[1]/resln[1];
              create_vector(knotBegin[1], knotEnd[1], fact, vv);

              //create the coordinates of the pointsVTK (nodes in FEM)

              count = 0;
              for(jj=0;jj<vv.size();jj++)
              {
                param[1] = vv[jj];
                geom[1] = computeGeometry(1, vv[jj]);

                for(ii=0;ii<uu.size();ii++)
                {
                  param[0] = uu[ii];
                  geom[0] = computeGeometry(0, uu[ii]);

                  pt[count] = pointsVTK->InsertNextPoint(geom[0], geom[1], 0.0);

                  GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

                  if(ndTemp->getParent() == NULL)
                    N = NN;
                  else
                    N = ndTemp->SubDivMat*NN;

                  if(domTemp == 0)
                    fact = ndTemp->computeValue(0, N);
                  else
                    fact = ndTemp->computeValue2(0, N);

                  scaVTK->InsertNextValue(fact);
                  scaVTK2->InsertNextValue(fact - GeomData.analyDBC->computeValue(0, geom[0], geom[1]) );

                  count++;
                } //for(ii=0;ii<uu.size();ii++)
              } //for(jj=0;jj<vv.size();jj++)

              quadVTK->GetPointIds()->SetId(0, pt[0]);
              quadVTK->GetPointIds()->SetId(1, pt[1]);
              quadVTK->GetPointIds()->SetId(2, pt[3]);
              quadVTK->GetPointIds()->SetId(3, pt[2]);

              //cout << ndTemp->getID() << '\t' << ndTemp->getSubdomainId() << endl;

              uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
              cellDataVTK->InsertNextValue(0);
              cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());
	    //}
        } //if( !nd->isCutElement() )
        else // the element is cutCell
        {
          vtkSmartPointer<vtkTriangle> triaVTK =  vtkSmartPointer<vtkTriangle>::New();

          myPoly *poly;

          for(ii=0; ii<ndTemp->subTrias.size(); ii++)
          {
            poly = ndTemp->subTrias[ii];

            if( poly->getDomainNumber() == 0 )
            {
              for(kk=0; kk<3; kk++)
              {
                geom = poly->GetPoint(kk);

                pt[kk] = pointsVTK->InsertNextPoint(geom[0], geom[1], 0.0);

                geometryToParametric(geom, param);
                GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

                if(ndTemp->getParent() == NULL)
                  N = NN;
                else
                  N = ndTemp->SubDivMat*NN;

                if( domTemp )
                  fact = ndTemp->computeValue2(0, N);
                else
                  fact = ndTemp->computeValue(0, N);

                scaVTK->InsertNextValue(fact);
                scaVTK2->InsertNextValue(fact - GeomData.analyDBC->computeValue(0, geom[0], geom[1]) );

                triaVTK->GetPointIds()->SetId(kk, pt[kk] );
              }

              //cellDataVTK->InsertNextValue(poly->getDomainNumber());

              uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
            }
          } //  for(ii=0; ii<nd->subTrias.size(); ii++)
        } // else
      }

      scaVTK->SetName("value");
      scaVTK2->SetName("force");

      uGridVTK->SetPoints(pointsVTK);

      //assign nodal coordinates and field data to uGridVTK
      // no need to create lookup table here. All this stuff can be done in Paraview

      uGridVTK->GetPointData()->SetScalars(scaVTK);
      uGridVTK->GetPointData()->AddArray(scaVTK2);
      uGridVTK->GetCellData()->SetScalars(cellDataVTK2);
    }
    else // for Stokes and Navier-Stokes
    {
      double vec[3]={0.0, 0.0 ,0.0};

      vecVTK->SetNumberOfComponents(3);
      //vecVTK->SetNumberOfTuples(count);
      vecVTK2->SetNumberOfComponents(3);
      //vecVTK2->SetNumberOfTuples(count);
      //scaVTK->SetNumberOfTuples(count);
      //scaVTK2->SetNumberOfTuples(count);

      for(ee=0; ee<activeElements.size(); ee++)
      {
        ndTemp = elems[activeElements[ee]];

        //if( ndTemp->getSubdomainId() == this_mpi_proc )
        //{
          knotBegin = ndTemp->getKnotBegin();
          knotEnd   = ndTemp->getKnotEnd();
          knotIncr  = ndTemp->getKnotIncrement();

          if( !(ndTemp->isCutElement()) )
          //if( ndTemp->getDomainNumber() <= 10 )
          {
            if( ndTemp->getDomainNumber() == 0 )
            {
              fact = knotIncr[0]/resln[0];
              create_vector(knotBegin[0], knotEnd[0], fact, uu);

              fact = knotIncr[1]/resln[1];
              create_vector(knotBegin[1], knotEnd[1], fact, vv);

              //create the coordinates of the pointsVTK (nodes in FEM)

              count = 0;
              for(jj=0;jj<vv.size();jj++)
              {
                param[1] = vv[jj];
                geom[1] = computeGeometry(1, vv[jj]);
                for(ii=0;ii<uu.size();ii++)
                {
                  param[0] = uu[ii];
                  geom[0] = computeGeometry(0, uu[ii]);

                  pt[count++] = pointsVTK->InsertNextPoint(geom[0], geom[1], 0.0);

                  GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

                  if(ndTemp->getParent() == NULL)
                  {
                    N = NN;
                    dN_dx = dNN_dx;
                    dN_dy = dNN_dy;
                  }
                  else
                  {
                    N = ndTemp->SubDivMat*NN;
                    dN_dx = ndTemp->SubDivMat*dNN_dx;
                    dN_dy = ndTemp->SubDivMat*dNN_dy;
                  }

                  vec[0] = ndTemp->computeValue(0, N);
                  vec[1] = ndTemp->computeValue(1, N);

                  vecVTK->InsertNextTuple(vec);

                  vec[0] = ndTemp->computeValueDot(0, N);
                  vec[1] = ndTemp->computeValueDot(1, N);

                  vecVTK2->InsertNextTuple(vec);

                  fact   = ndTemp->computeValue(2, N);
                  scaVTK->InsertNextValue(fact);
                  fact = ndTemp->computeValue(1, dN_dx) - ndTemp->computeValue(0, dN_dy);
                  scaVTK2->InsertNextValue(fact);

                }                                           // for(ii=0;ii<uu.size();ii++)
              }                                             // for(jj=0;jj<vv.size();jj++)

              quadVTK->GetPointIds()->SetId(0, pt[0]);
              quadVTK->GetPointIds()->SetId(1, pt[1]);
              quadVTK->GetPointIds()->SetId(2, pt[3]);
              quadVTK->GetPointIds()->SetId(3, pt[2]);

              uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
              cellDataVTK->InsertNextValue(0);
              cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());
            }
          }                                                 //if( !nd->isCutElement() )
          else                                              // the element is cutCell
          //if( ndTemp->getDomainNumber() == -1 )
          {
            vtkSmartPointer<vtkTriangle> triaVTK =  vtkSmartPointer<vtkTriangle>::New();

            myPoly *poly;

            for(ii=0; ii<ndTemp->subTrias.size(); ii++)
            {
              poly = ndTemp->subTrias[ii];

              if( poly->getDomainNumber() == 0 )
              {
                for(kk=0; kk<3; kk++)
                {
                  geom = poly->GetPoint(kk);
                  //cout << kk << '\t' << ptTemp[0] << '\t' << ptTemp[1] << endl;

                  pt[kk] = pointsVTK->InsertNextPoint(geom[0], geom[1], 0.0);

                  geometryToParametric(geom, param);
                  GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

                  if(ndTemp->getParent() == NULL)
                  {
                    N = NN;
                    dN_dx = dNN_dx;
                    dN_dy = dNN_dy;
                  }
                  else
                  {
                    N = ndTemp->SubDivMat*NN;
                    dN_dx = ndTemp->SubDivMat*dNN_dx;
                    dN_dy = ndTemp->SubDivMat*dNN_dy;
                  }

                  vec[0] = ndTemp->computeValue(0, N);
                  vec[1] = ndTemp->computeValue(1, N);

                  vecVTK->InsertNextTuple(vec);

                  vec[0] = ndTemp->computeValueDot(0, N);
                  vec[1] = ndTemp->computeValueDot(1, N);

                  vecVTK2->InsertNextTuple(vec);

                  fact   = ndTemp->computeValue(2, N);
                  scaVTK->InsertNextValue(fact);
                  fact = ndTemp->computeValue(1, dN_dx) - ndTemp->computeValue(0, dN_dy);
                  scaVTK2->InsertNextValue(fact);

                  triaVTK->GetPointIds()->SetId(kk, pt[kk] );
                }

                cellDataVTK->InsertNextValue(poly->getDomainNumber());
                cellDataVTK2->InsertNextValue(0);

                uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
              } // if( domainInclYesNo[domTemp] )
            }                                               //  for(ii=0; ii<nd->subTrias.size(); ii++)
          } // else
        //}
      }

      //cout << " jjjjjjjjjjjjjjjjjj " << endl;
      vecVTK->SetName("vel");
      //vecVTK2->SetName("acce");
      scaVTK->SetName("pres");
      scaVTK2->SetName("vortz");

      //assign nodal coordinates and field data to uGridVTK
      // no need to create lookup table here. All this stuff can be done in Paraview

      uGridVTK->SetPoints(pointsVTK);

      uGridVTK->GetPointData()->SetScalars(scaVTK);
      uGridVTK->GetPointData()->SetVectors(vecVTK);
      //uGridVTK->GetPointData()->AddArray(vecVTK2);
      uGridVTK->GetPointData()->AddArray(scaVTK2);

      uGridVTK->GetCellData()->SetScalars(cellDataVTK2);
    }

    return;
}




void  HBSplineCutFEM::postProcessSubTrias3D(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    int  dd=0, ii=0, jj=0, kk=0, ll=0, count=0, index=0;
    int  ind1=0, ind2=0, e=0, ee=0, gcount=0, ind=0, domTemp=0;
    int  nlocal = (degree[0]+1) * (degree[1] + 1) * (degree[2] + 1);

    VectorXd  NN(nlocal), dNN_dx(nlocal), dNN_dy(nlocal), dNN_dz(nlocal),  tempVec, tempVec2;
    VectorXd  N(nlocal), dN_dx(nlocal), dN_dy(nlocal), dN_dz(nlocal);
    VectorXd  vectmp(nlocal), rhsTemp;
    myPoint  knotIncr, knotBegin, knotEnd;

    double   fact=0.0, val1=0.0;

    vector<double>  uu, vv, ww;

    node*  ndTemp;

    vtkIdType pt[50], cellId;

    if(ndf == 1)
    {
      index = 0;
      for(ee=0; ee<activeElements.size(); ee++)
      {
        ndTemp = elems[activeElements[ee]];

        knotBegin = ndTemp->getKnotBegin();
        knotEnd   = ndTemp->getKnotEnd();
        knotIncr  = ndTemp->getKnotIncrement();

        if( !ndTemp->isCutElement() )
        {
          fact = knotIncr[0]/resln[0];
          create_vector(knotBegin[0], knotEnd[0], fact, uu);

          fact = knotIncr[1]/resln[1];
          create_vector(knotBegin[1], knotEnd[1], fact, vv);

          fact = knotIncr[2]/resln[2];
          create_vector(knotBegin[2], knotEnd[2], fact, ww);

          //create the coordinates of the pointsVTK (nodes in FEM)

          count = 0;
          for(jj=0;jj<vv.size();jj++)
          {
              param[1] = vv[jj];
              geom[1] = computeGeometry(1, vv[jj]);

              for(ii=0;ii<uu.size();ii++)
              {
                param[0] = uu[ii];
                geom[0] = computeGeometry(0, uu[ii]);

                pt[count] = pointsVTK->InsertNextPoint(geom[0], geom[1], 0.0);

                GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

                if(ndTemp->getParent() == NULL)
                  N = NN;
                else
                  N = ndTemp->SubDivMat*NN;

                if(domTemp == 0)
                  fact = ndTemp->computeValue(0, N);
                else
                  fact = ndTemp->computeValue2(0, N);

                scaVTK->InsertNextValue(fact);
                scaVTK2->InsertNextValue(fact - GeomData.analyDBC->computeValue(0, geom[0], geom[1]) );

                count++;
              } //for(ii=0;ii<uu.size();ii++)
          } //for(jj=0;jj<vv.size();jj++)

          quadVTK->GetPointIds()->SetId(0, pt[0]);
          quadVTK->GetPointIds()->SetId(1, pt[1]);
          quadVTK->GetPointIds()->SetId(2, pt[3]);
          quadVTK->GetPointIds()->SetId(3, pt[2]);

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
          //cellDataVTK->InsertNextValue(0);
          //cellDataVTK2->InsertNextValue(0);

        } //if( !nd->isCutElement() )
        else // the element is cutCell
        {
          vtkSmartPointer<vtkTriangle> triaVTK =  vtkSmartPointer<vtkTriangle>::New();

          myPoly *poly;

          for(ii=0; ii<ndTemp->subTrias.size(); ii++)
          {
            poly = ndTemp->subTrias[ii];

            if( poly->getDomainNumber() == 0 )
            {
              for(kk=0; kk<3; kk++)
              {
                geom = poly->GetPoint(kk);

                pt[kk] = pointsVTK->InsertNextPoint(geom[0], geom[1], 0.0);

                geometryToParametric(geom, param);
                GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

                if(ndTemp->getParent() == NULL)
                  N = NN;
                else
                  N = ndTemp->SubDivMat*NN;

                if( domTemp )
                  fact = ndTemp->computeValue2(0, N);
                else
                  fact = ndTemp->computeValue(0, N);

                scaVTK->InsertNextValue(fact);
                scaVTK2->InsertNextValue(fact - GeomData.analyDBC->computeValue(0, geom[0], geom[1]) );

                triaVTK->GetPointIds()->SetId(kk, pt[kk] );
              }

              cellDataVTK->InsertNextValue(poly->getDomainNumber());

              uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
            }
          } //  for(ii=0; ii<nd->subTrias.size(); ii++)
        } // else
      }

      scaVTK->SetName("value");
      scaVTK2->SetName("force");

      uGridVTK->SetPoints(pointsVTK);

      //assign nodal coordinates and field data to uGridVTK
      // no need to create lookup table here. All this stuff can be done in Paraview

      uGridVTK->GetPointData()->SetScalars(scaVTK);
      uGridVTK->GetPointData()->AddArray(scaVTK2);
      // create a write object and write uGridVTK to it
    }
    else // for Stokes and Navier-Stokes
    {
      double vec[3]={0.0, 0.0, 0.0};
      //vec[0] = vec[1] = vec[2] = 0.0;

      vecVTK->SetNumberOfComponents(3);
      //vecVTK->SetNumberOfTuples(count);
      vecVTK2->SetNumberOfComponents(3);
      //vecVTK2->SetNumberOfTuples(count);
      //scaVTK->SetNumberOfTuples(count);
      //scaVTK2->SetNumberOfTuples(count);

      for(ee=0; ee<activeElements.size(); ee++)
      {
        ndTemp = elems[activeElements[ee]];

        if( ndTemp->getSubdomainId() == this_mpi_proc )
        {
          knotBegin = ndTemp->getKnotBegin();
          knotEnd   = ndTemp->getKnotEnd();
          knotIncr  = ndTemp->getKnotIncrement();

          if( !(ndTemp->isCutElement()) )
          {
            fact = knotIncr[0]/resln[0];
            create_vector(knotBegin[0], knotEnd[0], fact, uu);

            fact = knotIncr[1]/resln[1];
            create_vector(knotBegin[1], knotEnd[1], fact, vv);

            fact = knotIncr[2]/resln[2];
            create_vector(knotBegin[2], knotEnd[2], fact, ww);

            //create the coordinates of the pointsVTK (nodes in FEM)

            count = 0;
            for(kk=0; kk<ww.size(); kk++)
            {
                param[2] = ww[kk];
                geom[2] = computeGeometry(2, ww[kk]);

            for(jj=0;jj<vv.size();jj++)
            {
                param[1] = vv[jj];
                geom[1] = computeGeometry(1, vv[jj]);

            for(ii=0;ii<uu.size();ii++)
            {
                param[0] = uu[ii];
                geom[0] = computeGeometry(0, uu[ii]);

                pt[count++] = pointsVTK->InsertNextPoint(geom[0], geom[1], geom[2]);

                GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

                if(ndTemp->getParent() == NULL)
                {
                  N = NN;
                  dN_dx = dNN_dx;
                  dN_dy = dNN_dy;
                  dN_dz = dNN_dz;
                }
                else
                {
                  N = ndTemp->SubDivMat*NN;
                  dN_dx = ndTemp->SubDivMat*dNN_dx;
                  dN_dy = ndTemp->SubDivMat*dNN_dy;
                  dN_dz = ndTemp->SubDivMat*dNN_dz;
                }

                vec[0] = ndTemp->computeValue(0, N);
                vec[1] = ndTemp->computeValue(1, N);
                vec[2] = ndTemp->computeValue(2, N);
                fact   = ndTemp->computeValue(3, N);

                vecVTK->InsertNextTuple(vec);
                //vecVTK2->InsertNextTuple(vec);
                scaVTK->InsertNextValue(fact);
                scaVTK2->InsertNextValue(fact);

            } // for(ii=0;ii<uu.size();ii++)
            } // for(jj=0;jj<vv.size();jj++)
            } // for(kk=0; kk<ww.size(); kk++)

            hexVTK->GetPointIds()->SetId(0, pt[0]);
            hexVTK->GetPointIds()->SetId(1, pt[1]);
            hexVTK->GetPointIds()->SetId(2, pt[3]);
            hexVTK->GetPointIds()->SetId(3, pt[2]);

            hexVTK->GetPointIds()->SetId(4, pt[4]);
            hexVTK->GetPointIds()->SetId(5, pt[5]);
            hexVTK->GetPointIds()->SetId(6, pt[7]);
            hexVTK->GetPointIds()->SetId(7, pt[6]);

            uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
            //cellDataVTK->InsertNextValue(0);
            cellDataVTK2->InsertNextValue(ndTemp->getSubdomainId());
          }                                                 //if( !nd->isCutElement() )
          else                                              // the element is cutCell
          {
            vtkSmartPointer<vtkTetra> tetVTK =  vtkSmartPointer<vtkTetra>::New();

            myPoly *poly;

            for(ii=0; ii<ndTemp->subTrias.size(); ii++)
            {
              poly = ndTemp->subTrias[ii];
              domTemp = poly->getDomainNumber();
              if( domainInclYesNo[domTemp] )
	          {
                for(kk=0; kk<4; kk++)
                {
                  geom = poly->GetPoint(kk);

                  pt[kk] = pointsVTK->InsertNextPoint(geom[0], geom[1], geom[2]);

                  geometryToParametric(geom, param);
                  GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

                  if(ndTemp->getParent() == NULL)
                  {
                    N = NN;
                    dN_dx = dNN_dx;
                    dN_dy = dNN_dy;
                    dN_dz = dNN_dz;
                  }
                  else
                  {
                    N = ndTemp->SubDivMat*NN;
                    dN_dx = ndTemp->SubDivMat*dNN_dx;
                    dN_dy = ndTemp->SubDivMat*dNN_dy;
                    dN_dz = ndTemp->SubDivMat*dNN_dz;
                  }

                  vec[0] = ndTemp->computeValue(0, N);
                  vec[1] = ndTemp->computeValue(1, N);
                  vec[2] = ndTemp->computeValue(2, N);
                  fact   = ndTemp->computeValue(3, N);

                  vecVTK->InsertNextTuple(vec);
                  //vecVTK2->InsertNextTuple(vec);
                  scaVTK->InsertNextValue(fact);
                  scaVTK2->InsertNextValue(fact);

                  tetVTK->GetPointIds()->SetId(kk, pt[kk] );
                }

                cellDataVTK->InsertNextValue(poly->getDomainNumber());

                uGridVTK->InsertNextCell(tetVTK->GetCellType(), tetVTK->GetPointIds());
              }                                             // if( domainInclYesNo[domTemp] )
            }                                               //  for(ii=0; ii<nd->subTrias.size(); ii++)
          }                                                 // else
        }
      }

      vecVTK->SetName("vel");
      //vecVTK2->SetName("force");
      scaVTK->SetName("pres");
      //scaVTK2->SetName("vortz");

      //assign nodal coordinates and field data to uGridVTK
      // no need to create lookup table here. All this stuff can be done in Paraview

      uGridVTK->SetPoints(pointsVTK);

      uGridVTK->GetPointData()->SetScalars(scaVTK);
      uGridVTK->GetPointData()->SetVectors(vecVTK);
      //uGridVTK->GetPointData()->AddArray(vecVTK2);
      //uGridVTK->GetPointData()->AddArray(scaVTK2);

      uGridVTK->GetCellData()->SetScalars(cellDataVTK2);
    }

    return;
}



void HBSplineCutFEM::plotGaussPointsElement()
{
    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK2   =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK2  =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK2  =  vtkSmartPointer<vtkVertex>::New();

    int  ee=0, ll=0, ii=0, gp=0, nGauss=0;

    vtkIdType  ptId;
    double  volume=0.0, *gws;
    myPoint *gps;
    node*  nd1;

    myPoint   knotIncr, knotBegin, knotSum;

    param.setZero();

    for(ee=0;ee<activeElements.size();ee++)
    {
        nd1 = elems[activeElements[ee]];

        knotBegin = nd1->getKnotBegin();
        knotIncr  = nd1->getKnotIncrement();
        knotSum   = nd1->getKnotSum();

        nGauss=0;
        if( nd1->getDomainNumber() == -1 )
        {
          nGauss = nd1->Quadrature.gausspoints.size();

          gps = &(nd1->Quadrature.gausspoints[0]);
          gws = &(nd1->Quadrature.gaussweights[0]);
        }

        if( nd1->getDomainNumber() == 0 )
        {
          nGauss = GeomData.gausspoints.size();

          gps = &(GeomData.gausspoints[0]);
          gws = &(GeomData.gaussweights[0]);

          volume += nd1->getVolume();
        } // else

          for(gp=0; gp<nGauss; gp++)
          {
              for(ii=0; ii<DIM; ii++)
                param[ii]  = 0.5*(knotIncr[ii] * gps[gp][ii] + knotSum[ii]);

              if( nd1->getDomainNumber() == -1 )
                volume += gws[gp];

              computeGeometry(param, geom);

              ptId = pointsVTK2->InsertNextPoint(geom[0], geom[1], geom[2]);

              vertexVTK2->GetPointIds()->SetId(0, ptId);

              uGridVTK2->InsertNextCell(vertexVTK2->GetCellType(), vertexVTK2->GetPointIds());
          } // for(gp=0;
    }

    PetscPrintf(MPI_COMM_WORLD, "\n Area (or Volume) from Gausspoints = %12.8f \n\n", volume );

    uGridVTK2->SetPoints(pointsVTK2);

    //scaVTK3->SetName("dist");
    //uGridVTK2->GetPointData()->SetScalars(scaVTK3);

    char fname[200];

    sprintf(fname,"%s%s", files.Ofile.asCharArray(),"-GPs.vtu");

    writerUGridVTK->SetFileName(fname);
    writerUGridVTK->SetInputData(uGridVTK2);
    writerUGridVTK->Write();

    return;
}





void HBSplineCutFEM::plotGaussPointsDirichletBoundary()
{
    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK2   =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK2  =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK2  =  vtkSmartPointer<vtkVertex>::New();

    int  aa=0, ee=0, ll=0, ii=0, gp=0, side=0, nGauss=0, levTemp=0;

    vtkIdType  ptId;
    double  volume=0.0, *gws, JacTemp;
    myPoint *gps;
    node*  nd1;
    myPoint   knotIncr, knotBegin, knotSum;

    param.setZero();

    for(ee=0;ee<activeElements.size();ee++)
    {
      nd1 = elems[activeElements[ee]];

      if( nd1->isBoundary() && ( (nd1->domNums.size()>1) || (nd1->domNums[0] == 0)) )
      {
        knotBegin = nd1->getKnotBegin();
        knotIncr  = nd1->getKnotIncrement();
        knotSum   = nd1->getKnotSum();

        levTemp = nd1->getLevel();

        if( nd1->DirichletData.size() > 0)
        {
          for(aa=0;aa<nd1->DirichletData.size();aa++)
          {
             side  = (int) (nd1->DirichletData[aa][0] - 1);

             if(nd1->domNums.size() > 1)
             {
                nGauss = nd1->BoundaryQuadrature[side].gausspoints.size() ;

                gps = &(nd1->BoundaryQuadrature[side].gausspoints[0]);
                gws = &(nd1->BoundaryQuadrature[side].gaussweights[0]);

                JacTemp = 1.0;
              }
              else
              {
                nGauss = GeomData.boundaryQuadrature3D[side].gausspoints.size();

                gps = &(GeomData.boundaryQuadrature3D[side].gausspoints[0]);
                gws = &(GeomData.boundaryQuadrature3D[side].gaussweights[0]);

                JacTemp = GeomData.boundaryJacobians[side][levTemp];
              }

              for(gp=0; gp<nGauss; gp++)
              {
                  for(ii=0; ii<DIM; ii++)
                    param[ii]  = 0.5*(knotIncr[ii] * gps[gp][ii] + knotSum[ii]);

                  //if( nd1->getDomainNumber() == -1 )
                    volume += gws[gp] * JacTemp;

                  computeGeometry(param, geom);

                  ptId = pointsVTK2->InsertNextPoint(geom[0], geom[1], geom[2]);

                  vertexVTK2->GetPointIds()->SetId(0, ptId);

                  uGridVTK2->InsertNextCell(vertexVTK2->GetCellType(), vertexVTK2->GetPointIds());
              } // for(gp=0;
	        } //for(aa=0;aa<nd1->DirichletData.size();aa++)
        }
      }
    }

    PetscPrintf(MPI_COMM_WORLD, "\n Boundary area for Dirichlet boundaries from Gausspoints = %12.8f \n\n", volume );

    uGridVTK2->SetPoints(pointsVTK2);

    char fname[200];

    sprintf(fname,"%s%s", files.Ofile.asCharArray(),"-GPs-DirichletBoundary.vtu");

    writerUGridVTK->SetFileName(fname);
    writerUGridVTK->SetInputData(uGridVTK2);
    writerUGridVTK->Write();

    return;
}



void HBSplineCutFEM::plotGaussPointsNeumannBoundary()
{
    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK2   =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK2  =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK2  =  vtkSmartPointer<vtkVertex>::New();

    int  aa=0, ee=0, ll=0, ii=0, gp=0, side=0, nGauss=0, levTemp=0;

    vtkIdType  ptId;
    double  volume=0.0, *gws, JacTemp=0.0;
    myPoint *gps;
    node*  nd1;

    myPoint   knotIncr, knotBegin, knotSum;

    param.setZero();

    for(ee=0;ee<activeElements.size();ee++)
    {
      nd1 = elems[activeElements[ee]];

      if( nd1->isBoundary() && ( (nd1->domNums.size()>1) || (nd1->domNums[0] == 0)) )
      {
        knotBegin = nd1->getKnotBegin();
        knotIncr  = nd1->getKnotIncrement();
        knotSum   = nd1->getKnotSum();

        levTemp = nd1->getLevel();

        if( nd1->NeumannData.size() > 0)
        {
          for(aa=0;aa<nd1->NeumannData.size();aa++)
          {
             side  = (int) (nd1->NeumannData[aa][0] - 1);

             if(nd1->domNums.size() > 1)
             {
                nGauss = nd1->BoundaryQuadrature[side].gausspoints.size() ;

                gps = &(nd1->BoundaryQuadrature[side].gausspoints[0]);
                gws = &(nd1->BoundaryQuadrature[side].gaussweights[0]);

                JacTemp = 1.0;
              }
              else
              {
                nGauss = GeomData.boundaryQuadrature3D[side].gausspoints.size();

                gps = &(GeomData.boundaryQuadrature3D[side].gausspoints[0]);
                gws = &(GeomData.boundaryQuadrature3D[side].gaussweights[0]);

                JacTemp = GeomData.boundaryJacobians[side][levTemp];
              }

              for(gp=0; gp<nGauss; gp++)
              {
                  for(ii=0; ii<DIM; ii++)
                    param[ii]  = 0.5*(knotIncr[ii] * gps[gp][ii] + knotSum[ii]);

                  //if( nd1->getDomainNumber() == -1 )
                    volume += gws[gp] * JacTemp;

                  computeGeometry(param, geom);

                  ptId = pointsVTK2->InsertNextPoint(geom[0], geom[1], geom[2]);

                  vertexVTK2->GetPointIds()->SetId(0, ptId);

                  uGridVTK2->InsertNextCell(vertexVTK2->GetCellType(), vertexVTK2->GetPointIds());
              } // for(gp=0;
      	  } //for(aa=0;aa<nd1->DirichletData.size();aa++)
        }
      }
    }

    PetscPrintf(MPI_COMM_WORLD, "\n Boundary area for Neumann boundaries from Gausspoints = %12.8f \n\n", volume );

    uGridVTK2->SetPoints(pointsVTK2);

    char fname[200];

    sprintf(fname,"%s%s", files.Ofile.asCharArray(),"-GPs-NeumannBoundary.vtu");

    writerUGridVTK->SetFileName(fname);
    writerUGridVTK->SetInputData(uGridVTK2);
    writerUGridVTK->Write();

    return;
}



