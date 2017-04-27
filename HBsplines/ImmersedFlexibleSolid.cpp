
#include "ImmersedFlexibleSolid.h"

#include "MpapTime.h"
#include "TimeFunction.h"
#include "ImmersedIntegrationElement.h"
#include "../utilities/SolutionData.h"
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




void  ImmersedFlexibleSolid::SetSolidElements(vector<vector<int> >& elemConn)
{
  nElem = elemConn.size() ;

  // Prepare patchElemProp data. Assign suitable id(in the database) based on the element name

  if(SolnData.ElemProp.n > 0)
    prepareElemProp();

  // Prepare patchMatlProp data. Assign suitable id(in the database) based on the material name

  if(SolnData.MatlProp.n > 0)
    prepareMatlProp();


  int   ee, ii, jj, kk, ind;

  cout << " nelem and  ndof " << nElem << '\t' << ndof << '\t' << npElem << endl;
  cout << " SolnData.ElemProp.id = " << SolnData.ElemProp[0].id << endl;

  elems = new LagrangeElement* [nElem];

  for(ee=0;ee<nElem;ee++)
  {
    elems[ee] = NewLagrangeElement(SolnData.ElemProp[elemConn[ee][1]].id);

    elems[ee]->elmType = elemConn[ee][0];
    elems[ee]->matType = elemConn[ee][1];
    elems[ee]->secType = elemConn[ee][2];

    vector<int>  vecTemp;
    for(ii=0;ii<elemConn[ee].size()-3;ii++)
      vecTemp.push_back(elemConn[ee][3+ii]);

    elems[ee]->nodeNums = vecTemp;

    elems[ee]->SolnData = &(SolnData);
    elems[ee]->GeomData = &(GeomData);
    //elems[ee]->prepareElemData();
  }

  ndof = elems[0]->GetDOFPerNode();
  npElem = elems[0]->get_nodes_per_element();

  nsize = ndof*npElem;

  NodeType.resize(nNode);
  ID.resize(nNode);
  for(ii=0;ii<nNode;ii++)
  {
    NodeType[ii].resize(ndof);
    ID[ii].resize(ndof);
    for(jj=0;jj<ndof;jj++)
    {
      NodeType[ii][jj] = -7777;
      ID[ii][jj] = -1;
    }
  }

    ///////////////////////////////////////////////////////////////////
    //
    // set GeomData details
    //
    ///////////////////////////////////////////////////////////////////

    GeomData.SetDimension(DIM);
    GeomData.SetNdof(ndof);
    GeomData.SetNGP(1);
    GeomData.build();

    cout << " AAAAAAAAAAAAAAAAA " << endl;

    ///////////////////////////////////////////////////////////////////
    //
    // set SolnData details
    //
    ///////////////////////////////////////////////////////////////////

    SolnData.initialise(nNode*ndof, 0, 0, 0);
    SolnData.SetPhysicsTypetoSolid();
    
    fluidAcce.resize(nNode*ndof);
    fluidAcce.setZero();
    fluidAccePrev = fluidAcce;
    fluidAcceCur = fluidAcce;

    soln.resize(nNode*ndof);
    soln.setZero();

    return;
}




void  ImmersedFlexibleSolid::SetNodalPositions(vector<vector<double> >&  datatemp)
{
  nNode = datatemp.size() ;

  int ii, jj;
  double  val;

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
void  ImmersedFlexibleSolid::SetNodeType(vector<vector<int> >& datatemp)
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



void  ImmersedFlexibleSolid::SetBoundaryConditions(vector<vector<double> >& datatemp)
{  
  int ii, jj, val;
  
  DirichletBCs = datatemp;
  
  cout << " DirichletBCs.size() = " << DirichletBCs.size() << endl;
  
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

  int  bb, type, nn2=0, dof;

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

            //if(IsBoundaryConditionTypeLagrange())
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

  int ii, jj, kk, ll, ee;

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

    double vec[3], vec2[3], xx, yy, zz;

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
      vec[2] = 0.0;
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
      {
        //xx = GeomData.NodePosOrig[ii][0];
        //yy = GeomData.NodePosOrig[ii][1];

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

      if(elems[0]->nodeNums.size() == 2) // line - for beam and truss elements
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<elems[ee]->nodeNums.size();ll++)
            lineVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(lineVTK->GetCellType(), lineVTK->GetPointIds());
        }
      }
      else if(elems[0]->nodeNums.size() == 3) // triangle element
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<3;ll++)
            triaVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
        }
      }
      else if(elems[0]->nodeNums.size() == 4) // quad element
      {
        for(ee=0;ee<nElem;ee++)
        {
          quadVTK->GetPointIds()->SetId(0, elems[ee]->nodeNums[0]);
          quadVTK->GetPointIds()->SetId(1, elems[ee]->nodeNums[1]);
          quadVTK->GetPointIds()->SetId(2, elems[ee]->nodeNums[3]);
          quadVTK->GetPointIds()->SetId(3, elems[ee]->nodeNums[2]);

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        }
      }
      else
      {
        cout << " ERROR in 'ImmersedFlexibleSolid::postProcess()' ... this type of element is not valid " << endl;
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

      if(elems[0]->nodeNums.size() == 4) // tet
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<4;ll++)
            tetraVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());
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
  int  ii, ind, bb;

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
  int  ii, ind, bb;

  cout << " nNode = " << nNode << '\t' << ndof << endl;

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
   int aa, ii, jj, ind1, ind2;

   SolnData.forceTemp.setZero();

    if(IsBoundaryConditionTypeLagrange())
    {
      for(aa=0;aa<ImmIntgElems.size();aa++)
      {
        //cout << " aa " << aa << endl;
        ImmIntgElems[aa]->IntegrateForceFlexible(0, 0, SolnData.forceTemp);
        //printVector(vectemp);
      }
    }

   //printf("\n\n");    printVector(SolnData.forceTemp);

   SolnData.interpolateForce();

   return;
}



void ImmersedFlexibleSolid::updateForce(double* data)
{
  int ii, jj, ind1, ind2;

  SolnData.forceTemp.setZero();
  cout << " ndof = " << ndof << endl;

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

  fluidAcceCur = SolnData.td[2]*fluidAcce + (1.0-SolnData.td[2])*fluidAccePrev;

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

