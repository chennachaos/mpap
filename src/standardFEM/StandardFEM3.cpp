
#include "StandardFEM.h"
#include "MpapTime.h"
#include "Functions.h"
#include "Files.h"
#include "MyString.h"
#include "KimMoinFlow.h"
#include "LagrangeElement.h"

extern MpapTime mpapTime;
extern Files files;



void StandardFEM::plotGeom(int a1, bool b1, int c1, bool d1, int* ffff)
{
    int  dd, ii, jj, kk, ll, nlocal, index, ind1, ind2, e, ee, count, gcount, ind;
    double xx, yy, zz;

    vtkIdType pt[8];

    uGridVTK->Reset();
    pointsVTK->Reset();
    procIdVTK->Reset();

    procIdVTK->SetName("procId");
    procIdVTK->SetNumberOfTuples(nElem);

    count = GeomData.NodePosOrig.size();

    if(DIM == 2)
    {
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
      {
        xx = GeomData.NodePosOrig[ii][0];
        yy = GeomData.NodePosOrig[ii][1];

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, 0.0);
      }

      if(elems[0]->nodeNums.size() == 3) // tria
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<3;ll++)
            triaVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());

          procIdVTK->InsertTuple1(ee, elems[ee]->getSubdomainId());
        }
      }
      else // quad
      {
        for(ee=0;ee<nElem;ee++)
        {
          quadVTK->GetPointIds()->SetId(0, elems[ee]->nodeNums[0]);
          quadVTK->GetPointIds()->SetId(1, elems[ee]->nodeNums[1]);
          quadVTK->GetPointIds()->SetId(2, elems[ee]->nodeNums[3]);
          quadVTK->GetPointIds()->SetId(3, elems[ee]->nodeNums[2]);

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        
          procIdVTK->InsertTuple1(ee, elems[ee]->getSubdomainId());
        }
      }
    }
    else
    {
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
      {
        xx = GeomData.NodePosOrig[ii][0];
        yy = GeomData.NodePosOrig[ii][1];
        zz = GeomData.NodePosOrig[ii][2];

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, zz);
      }

      if(elems[0]->nodeNums.size() == 4) // tet
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<4;ll++)
            tetraVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());

          procIdVTK->InsertTuple1(ee, elems[ee]->getSubdomainId());
        }
      }
      else
      {
        for(ee=0;ee<nElem;ee++)
        {
          hexVTK->GetPointIds()->SetId(0, elems[ee]->nodeNums[0]);
          hexVTK->GetPointIds()->SetId(1, elems[ee]->nodeNums[1]);
          hexVTK->GetPointIds()->SetId(2, elems[ee]->nodeNums[3]);
          hexVTK->GetPointIds()->SetId(3, elems[ee]->nodeNums[2]);

          hexVTK->GetPointIds()->SetId(4, elems[ee]->nodeNums[4]);
          hexVTK->GetPointIds()->SetId(5, elems[ee]->nodeNums[5]);
          hexVTK->GetPointIds()->SetId(6, elems[ee]->nodeNums[7]);
          hexVTK->GetPointIds()->SetId(7, elems[ee]->nodeNums[6]);

          uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
        
          procIdVTK->InsertTuple1(ee, elems[ee]->getSubdomainId());
        }
      }
    }

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetCellData()->AddArray(procIdVTK);

    // create a write object and write uGridVTK to it

    char fname[200];

    //VTKfilename = files.Ofile.asCharArray();

    sprintf(fname,"%s%s", files.Ofile.asCharArray(),"-Geom.vtu");

    //cout << VTKfilename << endl;
    //sprintf(fname,"%s%s%06d%s", VTKfilename, "-", filecount, ".vtu");

    writerUGridVTK->SetFileName(fname);
    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    return;
}




void  StandardFEM::postProcess(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    //cout << " StandardFEM::postProcess " << endl;

    if(this_mpi_proc != 0)
      return;

    int  dd, ii, jj, kk, ll, nlocal, index, ind1, ind2, e, ee, count, gcount, ind;
    double vec[3], vec2[3], xx, yy, zz;

    vtkIdType pt[4];

    //ii=0;
    //while(uGridVTK->GetPointData()->GetNumberOfArrays() != 0)
      //uGridVTK->GetPointData()->RemoveArray(ii++);

    ii=0;
    while(uGridVTK->GetCellData()->GetNumberOfArrays() != 0)
      uGridVTK->GetCellData()->RemoveArray(ii++);

    uGridVTK->Reset();
    pointsVTK->Reset();
    //procIdVTK->Reset();

    vecVTK->Reset();
    vecVTK2->Reset();

    scaVTK->Reset();
    scaVTK2->Reset();

    cellDataVTK->Reset();
    cellDataVTK2->Reset();

    //cout << " jjjjjjjjjjjjjjjjjj " << endl;

    count = GeomData.NodePosOrig.size();

    procIdVTK->SetName("procIde");
    procIdVTK->SetNumberOfTuples(nElem);

    scaVTK->SetNumberOfTuples(count);

    //if(ndof > 1)
    //{
      vecVTK->SetNumberOfComponents(3);
      vecVTK->SetNumberOfTuples(count);
      vecVTK2->SetNumberOfComponents(3);
      vecVTK2->SetNumberOfTuples(count);
      scaVTK2->SetNumberOfTuples(count);
    //}

    //cout << " BBBBBBBBBBBB " << PHYSICS_TYPE << endl;

  if(PHYSICS_TYPE == PHYSICS_TYPE_SOLID)
  {
    // subroutine for solid problem

    vecVTK->SetName("force");
    vecVTK2->SetName("disp");
    scaVTK->SetName("pres");

    //cout << " BBBBBBBBBBBB " << endl;

    if(DIM == 2)
    {
      vec[2] = 0.0;
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
      {
        //xx = GeomData.NodePosOrig[ii][0];
        //yy = GeomData.NodePosOrig[ii][1];

        xx = GeomData.NodePosNew[ii][0];
        yy = GeomData.NodePosNew[ii][1];

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, 0.0);

        kk = node_map_old_to_new[ii]*ndof;

        vec[0] = SolnData.force[kk];
        vec[1] = SolnData.force[kk+1];

        vecVTK->InsertTuple(ii, vec);

        vec[0] = SolnData.var1[kk];
        vec[1] = SolnData.var1[kk+1];

        vecVTK2->InsertTuple(ii, vec);

        //vec[0] = analy.computeForce(0, xx, yy);
        //vec[1] = analy.computeForce(1, xx, yy);

        //vec[0] = SolnData.var1Cur[kk];
        //vec[1] = SolnData.var1Cur[kk+1];

        //vec[0] = SolnData.var1DotCur[kk];
        //vec[1] = SolnData.var1DotCur[kk+1];
      }

      for(ee=0;ee<nElem;ee++)
      {
        if(elems[ee]->nodeNums.size() == 3)
        {
          for(ll=0;ll<3;ll++)
            triaVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());

          procIdVTK->InsertTuple1(ee, elems[ee]->getSubdomainId());
        }
        else if(elems[ee]->nodeNums.size() == 4)
        {
          quadVTK->GetPointIds()->SetId(0, elems[ee]->nodeNums[0]);
          quadVTK->GetPointIds()->SetId(1, elems[ee]->nodeNums[1]);
          quadVTK->GetPointIds()->SetId(2, elems[ee]->nodeNums[3]);
          quadVTK->GetPointIds()->SetId(3, elems[ee]->nodeNums[2]);

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());

          procIdVTK->InsertTuple1(ee, elems[ee]->getSubdomainId());
        }
        else
        {
          //cerr << " Wrong element type ... " << endl;
        }
      }
    }
    else //if(DIM == 3)
    {
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
      {
        //xx = GeomData.NodePosOrig[ii][0];
        //yy = GeomData.NodePosOrig[ii][1];
        //zz = GeomData.NodePosOrig[ii][2];

        xx = GeomData.NodePosNew[ii][0];
        yy = GeomData.NodePosNew[ii][1];
        zz = GeomData.NodePosNew[ii][2];

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, zz);

        kk = node_map_old_to_new[ii]*ndof;

        //vec[0] = SolnData.force[kk];
        //vec[1] = SolnData.force[kk+1];
        //vec[2] = SolnData.force[kk+2];

        vec[0] = SolnData.var1Dot[kk];
        vec[1] = SolnData.var1Dot[kk+1];
        vec[2] = SolnData.var1Dot[kk+2];

        vecVTK->InsertTuple(ii, vec);

        vec[0] = SolnData.var1[kk];
        vec[1] = SolnData.var1[kk+1];
        vec[2] = SolnData.var1[kk+2];

        vecVTK2->InsertTuple(ii, vec);

        //vec[0] = analy.computeForce(0, xx, yy);
        //vec[1] = analy.computeForce(1, xx, yy);

        //vec[0] = SolnData.var1Cur[kk];
        //vec[1] = SolnData.var1Cur[kk+1];

        //vec[0] = SolnData.var1DotCur[kk];
        //vec[1] = SolnData.var1DotCur[kk+1];

        if(ndof == 4)
        {
          scaVTK->SetTuple1(ii, SolnData.var1[kk+3]);
        }
      }

      for(ee=0;ee<nElem;ee++)
      {
        if(elems[ee]->nodeNums.size() == 4) // tet
        {
          for(ll=0;ll<4;ll++)
            tetraVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());

          procIdVTK->InsertTuple1(ee, elems[ee]->getSubdomainId());
        }
        else if(elems[ee]->nodeNums.size() == 8) // hex
        {
          hexVTK->GetPointIds()->SetId(0, elems[ee]->nodeNums[0]);
          hexVTK->GetPointIds()->SetId(1, elems[ee]->nodeNums[1]);
          hexVTK->GetPointIds()->SetId(2, elems[ee]->nodeNums[3]);
          hexVTK->GetPointIds()->SetId(3, elems[ee]->nodeNums[2]);

          hexVTK->GetPointIds()->SetId(4, elems[ee]->nodeNums[4]);
          hexVTK->GetPointIds()->SetId(5, elems[ee]->nodeNums[5]);
          hexVTK->GetPointIds()->SetId(6, elems[ee]->nodeNums[7]);
          hexVTK->GetPointIds()->SetId(7, elems[ee]->nodeNums[6]);

          uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());

          procIdVTK->InsertTuple1(ee, elems[ee]->getSubdomainId());
        }
        else
        {
        }
      }
    }

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->SetVectors(vecVTK);
    uGridVTK->GetPointData()->AddArray(vecVTK2);
}
else
{
    // subroutine for fluid problem

    //cout << " jjjjjjjjjjjjjjjjjj " << endl;
    //cout << " BBBBBBBBBBBB " << endl;

    KimMoinFlow  analy(SolnData.ElemProp[0].data[4], SolnData.ElemProp[0].data[5], 1.0);

    if(ndof == 1)
    {
      scaVTK->SetName("soln");
    }
    else
    {
      vecVTK->SetName("acce");
      vecVTK2->SetName("vel");
      scaVTK->SetName("pres");
      scaVTK2->SetName("vortz");
    }

    //cout << " BBBBBBBBBBBB " << endl;
    if(DIM == 2)
    {
      vec[2] = 0.0;
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
      {
        xx = GeomData.NodePosCur[ii][0];
        yy = GeomData.NodePosCur[ii][1];

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, 0.0);

        if(ndof == 1)
        {
          scaVTK->SetTuple1(ii, SolnData.var1Cur[node_map_old_to_new[ii]]);
          //scaVTK2->SetTuple1(ii, 0.0);
        }
        else
        {
          kk = node_map_old_to_new[ii]*ndof;
          vec[0] = analy.computeValue(0, xx, yy, mpapTime.cur) - SolnData.var1Cur[kk];
          vec[1] = analy.computeValue(1, xx, yy, mpapTime.cur) - SolnData.var1Cur[kk+1];

          vecVTK->InsertTuple(ii, vec);

          vec[0] = SolnData.var1Cur[kk];
          vec[1] = SolnData.var1Cur[kk+1];

          vecVTK2->InsertTuple(ii, vec);

          scaVTK->SetTuple1(ii, SolnData.var1Cur[kk+2]);
        }
      }
      //cout << " AAAAAAAAAAA " << endl;
      if(elems[0]->nodeNums.size() == 3)
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<3;ll++)
            triaVTK->GetPointIds()->SetId(ll, SolnData.node_map_new_to_old[elems[ee]->nodeNums[ll]]);

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
          
          procIdVTK->InsertTuple1(ee, elems[ee]->getSubdomainId());
        }
      }
      else
      {
        for(ee=0;ee<nElem;ee++)
        {
          quadVTK->GetPointIds()->SetId(0, SolnData.node_map_new_to_old[elems[ee]->nodeNums[0]]);
          quadVTK->GetPointIds()->SetId(1, SolnData.node_map_new_to_old[elems[ee]->nodeNums[1]]);
          quadVTK->GetPointIds()->SetId(2, SolnData.node_map_new_to_old[elems[ee]->nodeNums[3]]);
          quadVTK->GetPointIds()->SetId(3, SolnData.node_map_new_to_old[elems[ee]->nodeNums[2]]);

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
          
          procIdVTK->InsertTuple1(ee, elems[ee]->getSubdomainId());
        }
      }
    }
    else //if(DIM == 3)
    {
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
      {
        xx = GeomData.NodePosOrig[ii][0];
        yy = GeomData.NodePosOrig[ii][1];
        zz = GeomData.NodePosOrig[ii][2];

        //xx = GeomData.NodePosNew[ii][0];
        //yy = GeomData.NodePosNew[ii][1];
        //zz = GeomData.NodePosNew[ii][2];

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, zz);

        if(ndof == 1)
        {
          scaVTK->SetTuple1(ii, SolnData.var1[node_map_old_to_new[ii]]);
        }
        else
        {
          kk = node_map_old_to_new[ii]*ndof;
          scaVTK->SetTuple1(ii, SolnData.var1[kk+3]);

          //vec[0] = SolnData.force[kk];
          //vec[1] = SolnData.force[kk+1];
          //vec[2] = SolnData.force[kk+2];

          vec[0] = SolnData.var1Dot[kk];
          vec[1] = SolnData.var1Dot[kk+1];
          vec[2] = SolnData.var1Dot[kk+2];

          vecVTK->InsertTuple(ii, vec);

          vec[0] = SolnData.var1[kk];
          vec[1] = SolnData.var1[kk+1];
          vec[2] = SolnData.var1[kk+2];

          vecVTK2->InsertTuple(ii, vec);

          //vec[0] = analy.computeForce(0, xx, yy);
          //vec[1] = analy.computeForce(1, xx, yy);

          //vec[0] = SolnData.var1Cur[ii*ndof+0];
          //vec[1] = SolnData.var1Cur[ii*ndof+1];

          //vec[0] = SolnData.var1DotCur[ii*ndof];
          //vec[1] = SolnData.var1DotCur[ii*ndof+1];
	}
      }

      //for(jj=0; jj<DirichletBCs.size(); jj++)
      //{
        //vec[0] = 1; vec[1] = 0.0; vec[2]=0.0;
        //vecVTK->InsertTuple(DirichletBCs[jj][0], vec);
      //}

      if(elems[0]->nodeNums.size() == 4) // tet
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<4;ll++)
            tetraVTK->GetPointIds()->SetId(ll, SolnData.node_map_new_to_old[elems[ee]->nodeNums[ll]]);

          uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());

          procIdVTK->InsertTuple1(ee, elems[ee]->getSubdomainId());
        }
      }
      else
      {
        for(ee=0;ee<nElem;ee++)
        {
          hexVTK->GetPointIds()->SetId(0, SolnData.node_map_new_to_old[elems[ee]->nodeNums[0]]);
          hexVTK->GetPointIds()->SetId(1, SolnData.node_map_new_to_old[elems[ee]->nodeNums[1]]);
          hexVTK->GetPointIds()->SetId(2, SolnData.node_map_new_to_old[elems[ee]->nodeNums[3]]);
          hexVTK->GetPointIds()->SetId(3, SolnData.node_map_new_to_old[elems[ee]->nodeNums[2]]);

          hexVTK->GetPointIds()->SetId(4, SolnData.node_map_new_to_old[elems[ee]->nodeNums[4]]);
          hexVTK->GetPointIds()->SetId(5, SolnData.node_map_new_to_old[elems[ee]->nodeNums[5]]);
          hexVTK->GetPointIds()->SetId(6, SolnData.node_map_new_to_old[elems[ee]->nodeNums[7]]);
          hexVTK->GetPointIds()->SetId(7, SolnData.node_map_new_to_old[elems[ee]->nodeNums[6]]);

          uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());

          procIdVTK->InsertTuple1(ee, elems[ee]->getSubdomainId());
        }
      }
    }

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetCellData()->SetScalars(procIdVTK);

    if(ndof > 1)
    {
      uGridVTK->GetPointData()->SetVectors(vecVTK);
      uGridVTK->GetPointData()->AddArray(vecVTK2);
      uGridVTK->GetPointData()->AddArray(scaVTK2);
    }
    //cout << " AAAAAAAAAAA " << endl;
}

    // create a write object and write uGridVTK to it

    //cout << " AAAAAAAAAAA " << endl;

    char fname[200];

    //VTKfilename = files.Ofile.asCharArray();

    sprintf(fname,"%s%s%06d%s", files.Ofile.asCharArray(),"-",filecount, ".vtu");

    //cout << VTKfilename << endl;
    //sprintf(fname,"%s%s%06d%s", VTKfilename, "-", filecount, ".vtu");

    writerUGridVTK->SetFileName(fname);
    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    return;
}





