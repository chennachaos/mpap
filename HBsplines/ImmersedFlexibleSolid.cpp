
#include "ImmersedFlexibleSolid.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "ImmersedIntegrationElement.h"
#include "SolutionData.h"
#include "NewLagrangeElement.h"
#include "Files.h"
#include "StandardFEM.h"
#include "headersVTK.h"


extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;

extern Files files;


ImmersedFlexibleSolid::ImmersedFlexibleSolid()
{
  nElem = ndof = npElem = DIM = nsize = 0;  

  solver = NULL;
  elems = NULL;
}


ImmersedFlexibleSolid::ImmersedFlexibleSolid(int dd)
{
  DIM = dd;
  nElem = ndof = npElem = nsize = 0;

  solver = NULL;
  elems = NULL;
}



ImmersedFlexibleSolid::~ImmersedFlexibleSolid()
{
  if(solver != NULL)
    delete solver;
  solver = NULL;

  if(elems != NULL)
  {
    for(int ii=0;ii<nElem;ii++)
      delete elems[ii];

    delete [] elems;
    elems = NULL;
  }
}


void ImmersedFlexibleSolid::printInfo()
{
    printf("\n Immersed Flexible Body information \n");
    printf("\n ----------------------------------------------------------------------------- \n");
    printf("\t   Problem   Dimension        =  %5d\n\n", DIM);
    printf("\t   Number of Elements         =  %5d\n\n", nElem);
    printf("\t   Number of Nodes            =  %5d\n\n", nNode);
    printf("\t   DOF at each node           =  %5d\n\n", ndof);
    printf("\t   Total number of DOF        =  %5d\n\n", totalDOF);
    printf("\n ----------------------------------------------------------------------------- \n");

    return;
}




void  ImmersedFlexibleSolid::prepareElemProp()
{
    char fct[] = "ImmersedFlexibleSolid::prepareElemProp";

    if(SolnData.ElemProp.n < 1)
      prgError(1,fct,"'element property data ' missing!");

    char *elmTypeNames[] = STANDARDFEM_ELEMENT_TYPE_NAMES;

    for(int ii=0; ii<SolnData.ElemProp.n; ii++)
    {
      // assign correct id of the element type (as stored in the database) based on the element name
      SolnData.ElemProp[ii].id = SolnData.ElemProp[ii].name.which(elmTypeNames);

      //cout << " SolnData.ElemProp[ii].id = " << SolnData.ElemProp[ii].id << endl;
      if(SolnData.ElemProp[ii].id < 0)
        prgError(2,fct,"unknown element type name!");
    }

    return;
}




void  ImmersedFlexibleSolid::prepareMatlProp()
{
    char fct[] = "ImmersedFlexibleSolid::prepareMatlProp";

    if(SolnData.MatlProp.n < 1)
      prgError(1,fct,"'patch material property data ' missing!");

    char *matlTypeNames[] = MATERIAL_TYPE_NAMES;

    for(int ii=0; ii<SolnData.MatlProp.n; ii++)
    {
      // assign correct id of the material type (as stored in the database) based on the material name
      SolnData.MatlProp[ii].id = SolnData.MatlProp[ii].name.which(matlTypeNames);

      if(SolnData.MatlProp[ii].id < 0)
        prgError(2,fct,"unknown element type name!");
    }

    return;
}




void  ImmersedFlexibleSolid::setSolidElements(vector<vector<int> >& elemConn)
{
    nElem = elemConn.size() ;

    // Prepare patchElemProp data. Assign suitable id(in the database) based on the element name

    if(SolnData.ElemProp.n > 0)
      prepareElemProp();

    // Prepare patchMatlProp data. Assign suitable id(in the database) based on the material name

    if(SolnData.MatlProp.n > 0)
      prepareMatlProp();

    int   ee=0, ii=0, jj=0, kk=0, ind=0;
    bool  nContX=false, nContY=false;

    //cout << " nelem and  ndof " << nElem << '\t' << ndof << '\t' << npElem << endl;
    //cout << " SolnData.ElemProp[0].id = " << SolnData.ElemProp[0].id << endl;
    //cout << " SolnData.ElemProp[1].id = " << SolnData.ElemProp[1].id << endl;

    IEN.resize(nElem);
    LM.resize(nElem);

    elems = new LagrangeElement* [nElem];

    nElem_Constraint = 0;
    for(ee=0;ee<nElem;ee++)
    {
      if(SolnData.ElemProp[elemConn[ee][0]].id == 33)
      {
        nContX = true;
        nElem_Constraint++;
      }

      if(SolnData.ElemProp[elemConn[ee][0]].id == 34)
      {
        nContY = true;
        nElem_Constraint++;
      }

      elems[ee] = NewLagrangeElement(SolnData.ElemProp[elemConn[ee][0]].id);

      elems[ee]->elmType = elemConn[ee][0];
      elems[ee]->matType = elemConn[ee][1];
      elems[ee]->secType = elemConn[ee][2];

      vector<int>  vecTemp;
      for(ii=0;ii<elemConn[ee].size()-3;ii++)
        vecTemp.push_back(elemConn[ee][3+ii]);

      elems[ee]->nodeNums = vecTemp;

      IEN[ee] = vecTemp;

      elems[ee]->SolnData = &(SolnData);
      elems[ee]->GeomData = &(GeomData);
      //elems[ee]->prepareElemData();
    }

    ndof = elems[0]->getNdofPerNode();

    int ndof_temp1 = ndof;
    int ndof_temp2 = ndof_temp1;

    PetscPrintf(MPI_COMM_WORLD, "\n    nContX = %5d ...\n", nContX);
    PetscPrintf(MPI_COMM_WORLD, "\n    nContY = %5d ...\n", nContY);

    if(nContX || nContY)
    {
      ndof_temp2 += DIM;
    }

    NodeType.resize(nNode);
    ID.resize(nNode);

    for(ii=0;ii<nNode;ii++)
    {
      NodeType[ii].resize(ndof_temp1, -7777);
      ID[ii].resize(ndof_temp2, -1);
    }

    ///////////////////////////////////////////////////////////////////
    //
    // set GeomData details
    //
    ///////////////////////////////////////////////////////////////////

    GeomData.setDimension(DIM);
    GeomData.setNdof(ndof);
    GeomData.setNGP(1);
    GeomData.build();

    ///////////////////////////////////////////////////////////////////
    //
    // set SolnData details
    //
    ///////////////////////////////////////////////////////////////////

    //cout << " nElem_Constraint =  " << nElem_Constraint << endl;

    totalDOF = nNode*ndof;
    totalDOF += nElem_Constraint;

    SolnData.initialise(totalDOF, 0, 0, 0);
    SolnData.setPhysicsTypetoSolid();
    
    fluidAcce.resize(totalDOF);
    fluidAcce.setZero();
    fluidAccePrev = fluidAcce;
    fluidAcceCur = fluidAcce;

    soln.resize(totalDOF);
    soln.setZero();

    return;
}




void  ImmersedFlexibleSolid::setNodalPositions(vector<vector<double> >&  datatemp)
{
    nNode = datatemp.size() ;

    int  ii=0, jj=0;
    double  val=0.0;

    GeomData.NodePosOrig.resize(nNode);
    GeomData.NodePosCur.resize(nNode);
    GeomData.NodePosNew.resize(nNode);  
    GeomData.NodePosPrev.resize(nNode);
    GeomData.NodePosPrevCur.resize(nNode);  

    GeomData.specValNew.resize(nNode);
    GeomData.specValCur.resize(nNode);

    GeomData.acceNew.resize(nNode);
    GeomData.acceCur.resize(nNode);

    for(ii=0;ii<nNode;ii++)
    {
      for(jj=0;jj<DIM;jj++)
      {
        val = datatemp[ii][jj];
        GeomData.NodePosOrig[ii][jj] = val;
        GeomData.NodePosCur[ii][jj]  = val;
        GeomData.NodePosNew[ii][jj]  = val;
        GeomData.NodePosPrev[ii][jj]  = val;
        GeomData.NodePosPrevCur[ii][jj]  = val;

        GeomData.specValNew[ii][jj]  = 0.0;
        GeomData.specValCur[ii][jj]  = 0.0;

        GeomData.acceNew[ii][jj]  = 0.0;
        GeomData.acceCur[ii][jj]  = 0.0;
      }
    }

    //cout << " ImmersedFlexibleSolid::SetNodalPositions " << endl;

    return;
}





/*
void  ImmersedFlexibleSolid::setNodeType(vector<vector<int> >& datatemp)
{  
  int ii, jj, val;
  vector<int>::iterator itint;

  for(ii=0;ii<datatemp.size();ii++)
  {
    //cout << datatemp[ii][0] << '\t' << datatemp[ii][1] << endl;
    val = datatemp[ii][0] - 1;

    cout << val << endl;

    for(jj=1;jj<datatemp[0].size();jj++)
    {
      if(datatemp[ii][jj] == 1)
        NodeType[val][jj-1] = datatemp[ii][jj];
    }
  }
  
  return;
}
*/



void  ImmersedFlexibleSolid::setBoundaryConditions(vector<vector<double> >& datatemp)
{  
    int ii=0, jj=0, val=0;
  
    DirichletBCs = datatemp;

    //cout << " DirichletBCs.size() = " << DirichletBCs.size() << endl;
  
    for(ii=0;ii<datatemp.size();ii++)
    {
      //cout << datatemp[ii][0] << '\t' << datatemp[ii][1] << endl;
      val = datatemp[ii][0] - 1;

      //cout << " val = " << val << endl;

      NodeType[val][datatemp[ii][1]-1] = datatemp[ii][2];
    }

    return;
}



void ImmersedFlexibleSolid::computeInitialAcceleration()
{
    // compute initial accelerations

    //cout << " ImmersedFlexibleSolid::computeInitialAcceleration() " << endl;

    return;
}



void ImmersedFlexibleSolid::writeOutput()
{
    if(BC_ENFORCE_TYPE == BC_ENFORCE_TYPE_PENALTY)
    {
      cout << " can't print the output data when PENALTY method is used for Immersed Bodies " << endl;
      return;
    }

    int  bb=0, type=0, nn2=0, dof=0, ii=0, kk=0;
    double  val_out=0.0;

    for(bb=0; bb<OutputData.size(); bb++)
    {
      type  =  OutputData[bb][0];
      nn2   =  OutputData[bb][1] - 1;
      dof   =  OutputData[bb][2] - 1;

      char        tmp[100];
      MyString    tmpStr;

      //cout << id << '\t' << type << '\t' << nn1 << '\t' << nn2 << '\t' << dof << '\t' << fact << endl;

      switch(type)
      {
        case  1 : // total force on the

            //if(isBoundaryConditionTypeLagrange())
              //computeTotalForce();

            sprintf(tmp," \t %12.6E", totalForce[dof]);

        break;

        case  2 : // force on individual boundary point

             sprintf(tmp," \t %12.6E", SolnData.force[nn2*ndof+dof]);

        break;

        case  3 : // displacement

             sprintf(tmp," \t %12.6E", SolnData.var1[nn2*ndof+dof]);

        break;

        case  4 : // velocity

             sprintf(tmp," \t %12.6E", SolnData.var1Dot[nn2*ndof+dof]);

        break;

        case  5 : // acceleration

             sprintf(tmp," \t %12.6E", SolnData.var1DotDot[nn2*ndof+dof]);

        break;

        case  6 : // sum of contact forces

             kk = nNode*ndof;
             val_out=0.0;
             for(ii=0; ii<nElem_Constraint; ii++)
             {
               val_out += SolnData.var1[kk+ii];
             }
             sprintf(tmp," \t %12.6E", val_out);
        break;

        default : 

             prgError(1,"HBSplineFEM::writeImmersedSolidOutput()","invalid value of 'type'!");

        break;
      }

      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);
    }

    return;
}



void ImmersedFlexibleSolid::postProcess(int index)
{
    //cout << " ImmersedFlexibleSolid::postProcess(int index) " << nNode << endl;

    int  ii=0, jj=0, kk=0, ll=0, ee=0;

    vtkSmartPointer<vtkPoints>              pointsVTK  =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnstructuredGrid>    uGridVTK   =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkVertex>              vertVTK    =  vtkSmartPointer<vtkVertex>::New();
    vtkSmartPointer<vtkLine>              lineVTK    =  vtkSmartPointer<vtkLine>::New();
    vtkSmartPointer<vtkTriangle>          triaVTK    =  vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkQuad>              quadVTK    =  vtkSmartPointer<vtkQuad>::New();
    vtkSmartPointer<vtkTetra>             tetraVTK   =  vtkSmartPointer<vtkTetra>::New();
    vtkSmartPointer<vtkHexahedron>        hexVTK     =     vtkSmartPointer<vtkHexahedron>::New();

    vtkSmartPointer<vtkFloatArray>          scaVTK   =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>          vecVTK   =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>          vecVTK2  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>          vecVTK3  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>          vecVTK4  =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    double vec[3], vec2[3], xx=0.0, yy=0.0, zz=0.0;

    scaVTK->SetName("pres");
    vecVTK->SetName("disp");
    vecVTK2->SetName("velo");
    vecVTK3->SetName("acce");
    vecVTK4->SetName("force");

    vecVTK->SetNumberOfComponents(3);
    vecVTK->SetNumberOfTuples(nNode);
    vecVTK2->SetNumberOfComponents(3);
    vecVTK2->SetNumberOfTuples(nNode);
    vecVTK3->SetNumberOfComponents(3);
    vecVTK3->SetNumberOfTuples(nNode);
    vecVTK4->SetNumberOfComponents(3);
    vecVTK4->SetNumberOfTuples(nNode);

    vtkIdType   pt[2];

    if(DIM == 2)
    {
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
      {
        xx = GeomData.NodePosNew[ii][0];
        yy = GeomData.NodePosNew[ii][1];

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, 0.0);

        //kk = node_map_old_to_new[ii]*ndof;
        kk = ii*ndof;

        vec[0] = SolnData.var1[kk];
        vec[1] = SolnData.var1[kk+1];
        vec[2] = 0.0;

        vecVTK->InsertTuple(ii, vec);

        vec[0] = SolnData.var1Dot[kk];
        vec[1] = SolnData.var1Dot[kk+1];
        vec[2] = 0.0;

        vecVTK2->InsertTuple(ii, vec);

        vec[0] = SolnData.var1DotDot[kk];
        vec[1] = SolnData.var1DotDot[kk+1];
        vec[2] = 0.0;

        vecVTK3->InsertTuple(ii, vec);

        vec[0] = SolnData.force[kk];
        vec[1] = SolnData.force[kk+1];
        vec[2] = 0.0;

        vecVTK4->InsertTuple(ii, vec);
      }

      for(ee=0;ee<nElem;ee++)
      {
        if(elems[ee]->nodeNums.size() == 2) // line - for beam and truss elements
        {
          lineVTK->GetPointIds()->SetId(0, elems[ee]->nodeNums[0]);
          lineVTK->GetPointIds()->SetId(1, elems[ee]->nodeNums[1]);

          uGridVTK->InsertNextCell(lineVTK->GetCellType(), lineVTK->GetPointIds());
        }
        else if(elems[ee]->nodeNums.size() == 3) // triangle element
        {
          triaVTK->GetPointIds()->SetId(0, elems[ee]->nodeNums[0]);
          triaVTK->GetPointIds()->SetId(1, elems[ee]->nodeNums[1]);
          triaVTK->GetPointIds()->SetId(2, elems[ee]->nodeNums[2]);

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
        }
        else if(elems[ee]->nodeNums.size() == 4) // quad element
        {
          quadVTK->GetPointIds()->SetId(0, elems[ee]->nodeNums[0]);
          quadVTK->GetPointIds()->SetId(1, elems[ee]->nodeNums[1]);
          quadVTK->GetPointIds()->SetId(2, elems[ee]->nodeNums[3]);
          quadVTK->GetPointIds()->SetId(3, elems[ee]->nodeNums[2]);

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        }
        else
        {
          //cout << " ERROR in 'ImmersedFlexibleSolid::postProcess()' ... this type of element is not valid " << endl;
        }
      }
    }
    else //if(DIM == 3)
    {
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
      {
        xx = GeomData.NodePosNew[ii][0];
        yy = GeomData.NodePosNew[ii][1];
        zz = GeomData.NodePosNew[ii][2];

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, zz);

        //kk = node_map_old_to_new[ii]*ndof;
        kk = ii*ndof;

        vec[0] = SolnData.var1[kk];
        vec[1] = SolnData.var1[kk+1];
        vec[2] = SolnData.var1[kk+2];

        vecVTK->InsertTuple(ii, vec);

        vec[0] = SolnData.var1Dot[kk];
        vec[1] = SolnData.var1Dot[kk+1];
        vec[2] = SolnData.var1Dot[kk+2];

        vecVTK2->InsertTuple(ii, vec);

        vec[0] = SolnData.var1DotDot[kk];
        vec[1] = SolnData.var1DotDot[kk+1];
        vec[2] = SolnData.var1DotDot[kk+2];

        vecVTK3->InsertTuple(ii, vec);

        vec[0] = SolnData.force[kk];
        vec[1] = SolnData.force[kk+1];
        vec[2] = SolnData.force[kk+2];

        vecVTK4->InsertTuple(ii, vec);

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
        }
        else
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
        }
      }
    }

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetPointData()->SetVectors(vecVTK);
    uGridVTK->GetPointData()->AddArray(vecVTK2);
    uGridVTK->GetPointData()->AddArray(vecVTK3);
    uGridVTK->GetPointData()->AddArray(vecVTK4);

    char fname[50];

    sprintf(fname,"%s%s%02d%s%06d%s", files.Ofile.asCharArray(),"-FB",id,"-", index, ".vtu");

    writerUGridVTK->SetFileName(fname);

#if VTK_MAJOR_VERSION == 5
    writerUGridVTK->SetInput(uGridVTK);
#else
    writerUGridVTK->SetInputData(uGridVTK);
#endif

    writerUGridVTK->Write();

  return;
}


int ImmersedFlexibleSolid::updatePointPositions()
{
  if(DIM == 1)
    updatePointPositions1D();
  else if(DIM == 2)
    updatePointPositions2D();
  else if(DIM == 3)
    updatePointPositions3D();

  return (totalDOF>0);
}



void ImmersedFlexibleSolid::updatePointPositions1D()
{
  return;
}



void ImmersedFlexibleSolid::updatePointPositions2D()
{
    int  ii=0, ind=0, bb=0;

    //cout << " nNode = " << nNode << '\t' << ndof << endl;

    for(bb=0;bb<nNode;bb++)
    {
      for(ii=0;ii<DIM;ii++)
      {
        ind = bb*ndof+ii;

        //cout << GeomData.NodePosNew[bb][ii] << '\t' << GeomData.NodePosOrig[bb][ii] << '\t' << SolnData.var1[ind] << endl;

        GeomData.NodePosNew[bb][ii] = GeomData.NodePosOrig[bb][ii] + SolnData.var1[ind];
        GeomData.NodePosCur[bb][ii] = GeomData.NodePosOrig[bb][ii] + SolnData.var1Cur[ind];

        // second-order based formulation
        GeomData.specValNew[bb][ii] = SolnData.var1Dot[ind];
        GeomData.specValCur[bb][ii] = SolnData.var1DotCur[ind];

        GeomData.acceNew[bb][ii] = SolnData.var1DotDot[ind];
        GeomData.acceCur[bb][ii] = SolnData.var1DotDotCur[ind];
      }
    }

    //cout << " nNode = " << nNode << '\t' << ndof << endl;

    return;
}


void ImmersedFlexibleSolid::updatePointPositions3D()
{
    int  ii=0, ind=0, bb=0;

    //cout << " nNode = " << nNode << '\t' << ndof << endl;

    for(bb=0; bb<nNode; bb++)
    {
      for(ii=0; ii<3; ii++)
      {
        ind = bb*ndof+ii;

        GeomData.NodePosNew[bb][ii] = GeomData.NodePosOrig[bb][ii] + SolnData.var1[ind];
        GeomData.NodePosCur[bb][ii] = GeomData.NodePosOrig[bb][ii] + SolnData.var1Cur[ind];

        // second-order based formulation
        GeomData.specValNew[bb][ii] = SolnData.var1Dot[ind];
        GeomData.specValCur[bb][ii] = SolnData.var1DotCur[ind];
      }
    }

    //cout << " nNode = " << nNode << '\t' << ndof << endl;

    return;
}




void ImmersedFlexibleSolid::updateForce()
{
    int  aa=0, ii=0, jj=0, ind1=0, ind2=0;

    SolnData.forceTemp.setZero();

    if(isBoundaryConditionTypeLagrange())
    {
      for(aa=0;aa<ImmIntgElems.size();aa++)
      {
        //cout << " aa " << aa << endl;
        ImmIntgElems[aa]->integrateForceFlexible(0, 0, SolnData.forceTemp);
        //printVector(vectemp);
      }
    }

   //printf("\n\n");    printVector(SolnData.forceTemp);

   SolnData.interpolateForce();

   return;
}



void ImmersedFlexibleSolid::updateForce(double* data)
{
    int  ii=0, jj=0, ind1=0, ind2=0;

    SolnData.forceTemp.setZero();

    for(ii=0;ii<nNode;ii++)
    {
      ind1 = ndof*ii;
      ind2 = DIM*ii;

      // second-order based formulation
      for(jj=0;jj<DIM;jj++)
        SolnData.forceTemp[ind1+jj] = -data[ind2+jj];
    }

    //printf("\n\n");    printVector(SolnData.forceTemp);

    SolnData.interpolateForce();
    //cout << " ndof = " << ndof << endl;
    fluidAcceCur = SolnData.td[2]*fluidAcce + (1.0-SolnData.td[2])*fluidAccePrev;
    //cout << " ndof = " << ndof << endl;
    return;
}




void  ImmersedFlexibleSolid::updateDisplacement(double* data)
{
    ///////////////////////////////
    // velocity for the solid element nodes
    // velocity is considered as the primary variable for the solid problem
    // in order to be consistent with the fluid problem

    //cout << " totalDOF " << totalDOF << endl;
    soln.setZero();
    // update solution vector
    for(int kk=0;kk<totalDOF;kk++)
    {
      //cout << kk << '\t' << assy4r[kk] << endl;
      //cout << kk << '\t' << data[kk] << endl;
      soln[assy4r[kk]] = data[kk];
    }

    SolnData.var1Dot = soln;

    //cout << " totalDOF " << totalDOF << endl;
    //printVector(SolidSolnData.disp);

    return;
}



