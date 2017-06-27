
#include "HBSplineFEM.h"
#include "MpapTime.h"
#include "Functions.h"
#include "Files.h"
#include "MyString.h"
#include "headersVTK.h"


extern MpapTime mpapTime;
extern Files files;



void HBSplineFEM::plotGeom(int val1, bool flag2, int col, bool PLOT_KNOT_LINES, int* resln)
{
    //cout << "     HBSplineFEM: plotgeometry ...\n\n";

    pointsVTK->Reset();
    //lineVTK->Reset();
    //quadVTK->Reset();
    uGridVTK->Reset();
    //actorVTK->Reset();
    //writerUGridVTK->Reset();

    if(ndm == 1)
      plotGeom1D(val1, flag2, col, PLOT_KNOT_LINES, resln);
    else if(ndm == 2)
      plotGeom2D(val1, flag2, col, PLOT_KNOT_LINES, resln);
    else
      plotGeom3D(val1, flag2, col, PLOT_KNOT_LINES, resln);

    char fname[50];

    sprintf(fname,"%s%s", files.Ofile.asCharArray(),"-geom.vtu");

    writerUGridVTK->SetFileName(fname);

#if VTK_MAJOR_VERSION == 5
    writerUGridVTK->SetInput(uGridVTK);
#else
    writerUGridVTK->SetInputData(uGridVTK);
#endif

    writerUGridVTK->Write();

    return;
}



void HBSplineFEM::plotGeom1D(int val1, bool flag2, int col, bool PLOT_KNOT_LINES, int* resln)
{ 
    int ii, jj, n1, n2, ind1, ind2, ll;

    vtkIdType pt0, pt1, pt2;
    double  *tmp;

    for(ii=0;ii<elems.size();ii++)
    {
       //if( elems[ii]->isLeaf() && !(elems[ii]->isGhost()) &&  elems[ii]->isActive())
       if( elems[ii]->isActive() )
       {
          //cout << " Node # " << elems[ii]->getID() << '\t' << elems[ii]->isLeaf() << '\t' << elems[ii]->isGhost() << endl;

          tmp = elems[ii]->getKnots(0);

          param[0] = tmp[0];
          ComputeGeometry(param, geom);
          pt0 = pointsVTK->InsertNextPoint(geom[0], 0.0, 0.0);

          param[0] = tmp[1];
          ComputeGeometry(param, geom);
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

    uGridVTK->SetPoints(pointsVTK);

    return;
}




void HBSplineFEM::plotGeom2D(int val1, bool flag2, int col, bool PLOT_KNOT_LINES, int* resln)
{
    int ee, ll;

    vtkIdType pt[4];
    double  *tmp1, *tmp2, dist;

    double  x1, y1, x2, y2;

    for(ee=0;ee<elems.size();ee++)
    {
       //cout << " Node # " << elems[ee]->getID() << '\t' << elems[ee]->isLeaf() << endl;
       if( elems[ee]->isLeaf() && !(elems[ee]->isGhost()) )
       //if( elems[ee]->isLeaf() )
       {
          tmp1 = elems[ee]->getKnots(Dir1);
          tmp2 = elems[ee]->getKnots(Dir2);
          
          //cout << tmp1[0] << '\t' << tmp1[1] << '\t' << tmp2[0] << '\t' << tmp2[1] << endl;

          x1 = ComputeGeometry(0, tmp1[0]);
          x2 = ComputeGeometry(0, tmp1[1]);

          y1 = ComputeGeometry(1, tmp2[0]);
          y2 = ComputeGeometry(1, tmp2[1]);


          pt[0] = pointsVTK->InsertNextPoint(x1, y1, 0.0);
          pt[1] = pointsVTK->InsertNextPoint(x2, y1, 0.0);
          pt[2] = pointsVTK->InsertNextPoint(x2, y2, 0.0);
          pt[3] = pointsVTK->InsertNextPoint(x1, y2, 0.0);

          for(ll=0;ll<4;ll++)
             quadVTK->GetPointIds()->SetId(ll, pt[ll]);

          //cellDataVTK->InsertNextValue(typetemp);
          //cellDataVTK2->InsertNextValue(elems[ee]->isGhost());

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
       }
    }

    uGridVTK->SetPoints(pointsVTK);

    return;
}


void HBSplineFEM::plotGeom3D(int val1, bool flag2, int col, bool PLOT_KNOT_LINES, int* resln)
{
    int ee, ll;

    vtkIdType pt[8];
    double  *tmp1, *tmp2, *tmp3, xx[2], yy[2], zz[2];

    time_t tstart, tend;
    tstart = time(0);

    for(ee=0;ee<elems.size();ee++)
    {
       //cout << " Node # " << elems[ee]->getID() << '\t' << elems[ee]->isLeaf() << endl;
       if( elems[ee]->isLeaf() && !(elems[ee]->isGhost()) )
       {
          //elems[ee]->printSelf();
          tmp1 = elems[ee]->getKnots(Dir1);
          tmp2 = elems[ee]->getKnots(Dir2);
          tmp3 = elems[ee]->getKnots(Dir3);
          
          //cout << tmp1[0] << '\t' << tmp1[1] << '\t' << tmp2[0] << '\t' << tmp2[1] << '\t' << tmp3[0] << '\t' << tmp3[1]<< endl;

          xx[0] = ComputeGeometry(0, tmp1[0]);
          xx[1] = ComputeGeometry(0, tmp1[1]);
          yy[0] = ComputeGeometry(1, tmp2[0]);
          yy[1] = ComputeGeometry(1, tmp2[1]);
          zz[0] = ComputeGeometry(2, tmp3[0]);
          zz[1] = ComputeGeometry(2, tmp3[1]);

          pt[0] = pointsVTK->InsertNextPoint(xx[0], yy[0], zz[0]);
          pt[1] = pointsVTK->InsertNextPoint(xx[1], yy[0], zz[0]);
          pt[2] = pointsVTK->InsertNextPoint(xx[1], yy[0], zz[1]);
          pt[3] = pointsVTK->InsertNextPoint(xx[0], yy[0], zz[1]);
          pt[4] = pointsVTK->InsertNextPoint(xx[0], yy[1], zz[0]);
          pt[5] = pointsVTK->InsertNextPoint(xx[1], yy[1], zz[0]);
          pt[6] = pointsVTK->InsertNextPoint(xx[1], yy[1], zz[1]);
          pt[7] = pointsVTK->InsertNextPoint(xx[0], yy[1], zz[1]);
          //

          for(ll=0;ll<8;ll++)
            hexVTK->GetPointIds()->SetId(ll, pt[ll]);
          
          cellDataVTK->InsertNextValue(elems[ee]->getLevel());
          
          uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
       }
    }

    tend = time(0);
    cout << " plotGeom3D "<< difftime(tend, tstart) <<" second(s)."<< endl;

    //cellDataVTK2->SetName("IsGhost");
    cellDataVTK->SetName("Level");

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetCellData()->SetScalars(cellDataVTK);
    //uGridVTK->GetCellData()->AddArray(cellDataVTK2);

    return;
}




void  HBSplineFEM::postProcessFlow(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    if(CREATE_POSTPROCESS_GRID)
    {
      if(DIM == 1)
        createPostProcessGrid1D(vartype, vardir, nCol, umnxflag, umin, umax, resln);
      else if(DIM == 2)
        createPostProcessGrid2D(vartype, vardir, nCol, umnxflag, umin, umax, resln);
      else
        createPostProcessGrid3D(vartype, vardir, nCol, umnxflag, umin, umax, resln);

       CREATE_POSTPROCESS_GRID = false;
    }

    if(DIM == 1)
      postProcess1D(vartype, vardir, nCol, umnxflag, umin, umax, resln);
    else if(DIM == 2)
      postProcess2D(vartype, vardir, nCol, umnxflag, umin, umax, resln);
    else
      postProcess3D(vartype, vardir, nCol, umnxflag, umin, umax, resln);

    return;
}



void  HBSplineFEM::createPostProcessGrid1D(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
  //cout << " HBSplineFEM::createPostProcessGrid1D ... yet to be implemented " << endl;
  return;
}



// for the general mesh (including Hierarchical splines)
void  HBSplineFEM::createPostProcessGrid2D(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    if(resln[0] > 100 || resln[1] > 100)
    {
      cout << " resolution is too high for this postprocess algorithm "  << endl;
      exit(1);
    }
    
    uGridVTK->Reset();
    pointsVTK->Reset();
    cellDataVTK->Reset();
    cellDataVTK2->Reset();
    scaVTK->Reset();
    scaVTK2->Reset();
    vecVTK->Reset();
    vecVTK2->Reset();

    int  dd, ii, jj, kk, ll, count, nlocal, index, ind1, ind2, ee, gcount, e;

    nlocal = (degree[0]+1) * (degree[1] + 1);

    double   fact, uleft, uright, *tmp1, *tmp2, incr1, incr2;
    vector<double>  xx, yy, uu, vv;

    node* nd1;

    vtkIdType pt[4];
    
    //cout << " Llllllllll " << endl;
//
    globalK2 *= 0.0;
    for(e=0;e<activeElements.size();e++)
    {
       //cout << elems[ee]->getID() << endl;
       elems[activeElements[e]]->MatrixToMapResult(1,1, globalK2);
    }
    
    //cout << " Llllllllll " << endl;

    solver3.compute(globalK2);
    if(solver3.info() != Success)
    {
      // decomposition failed
      cout << " decomposition failed " << endl;
      return;
    }
//
    //cout << " Llllllllll " << endl;

    index = 0;
    for(e=0;e<activeElements.size();e++)
    {
       ee = activeElements[e];

           nd1 = elems[ee];
           //cout << " Node # " << nd1->getID() << endl;

           tmp1 = nd1->getKnots(0);
           tmp2 = nd1->getKnots(1);

           //printf("\t tmp[0] and tmp[1]  ... : %12.8f\t%12.8f\n", tmp[0], tmp[1] );
               
           fact = (tmp1[1] - tmp1[0])/resln[0];
           create_vector(tmp1[0], tmp1[1], fact, uu);

           fact = (tmp2[1] - tmp2[0])/resln[1];
           create_vector(tmp2[0], tmp2[1], fact, vv);
           
           //cout << uu << endl;           cout << vv << endl;

           //cellDataVTK->InsertNextValue(0);

           for(jj=0;jj<vv.size();jj++)
           {
              for(ii=0;ii<uu.size();ii++)
              {
                 param[0] = uu[ii]; param[1] = vv[jj];
                 ComputeGeometry(param, geom);

                 xx.push_back(geom[0]);
                 yy.push_back(geom[1]);

                 pointsVTK->InsertNextPoint(geom[0], geom[1], 0.0);
              }
           }

           // create the connectivity of elements/nodes
           // use triangle/quadrilateral as required
              
           for(jj=0;jj<vv.size()-1;jj++)
           {
              ind1 = index + (resln[0]+1) * jj ;
              ind2 = index + (resln[0]+1) * (jj+1) ;
                 
              for(ii=0;ii<uu.size()-1;ii++)
              {
                 pt[0] = ind1 + ii;
                 pt[1] = pt[0] + 1;
                 pt[3] = ind2 + ii;
                 pt[2] = pt[3] + 1;

                 for(ll=0;ll<4;ll++)
                    quadVTK->GetPointIds()->SetId(ll, pt[ll]);
          
                 uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
              }
           }
           index = pointsVTK->GetNumberOfPoints();
    }

    //cout << count << '\t' << xx.size() << '\t' << yy.size() << endl;

    count = xx.size();

    if(ndf == 1)
    {
      scaVTK->SetName("value");
      scaVTK->SetNumberOfTuples(count);

      scaVTK2->SetName("force");
      scaVTK2->SetNumberOfTuples(count);

      for(ii=0;ii<count;ii++)
      {
        scaVTK->SetTuple1(ii, 0.0);
        scaVTK2->SetTuple1(ii, 0.0);
      }
    
      //assign nodal coordinates and field data to uGridVTK
      // no need to create lookup table here. All this stuff can be done in Paraview

      uGridVTK->SetPoints(pointsVTK);
      uGridVTK->GetPointData()->SetScalars(scaVTK);
    }
    if(ndf > 1) // for Stokes and Navier-Stokes
    {
       vecVTK->SetName("vel");
       //vecVTK2->SetName("grad");
       vecVTK2->SetName("force");
       scaVTK->SetName("pres");
       scaVTK2->SetName("vortz");

       vecVTK->SetNumberOfComponents(3);
       vecVTK->SetNumberOfTuples(count);
       vecVTK2->SetNumberOfComponents(3);
       vecVTK2->SetNumberOfTuples(count);
       scaVTK->SetNumberOfTuples(count);
       scaVTK2->SetNumberOfTuples(count);
    
       double vec[3];
    
       vec[0] = vec[1] = vec[2] = 0.0;
       for(ii=0;ii<count;ii++)
       {
         vecVTK->InsertTuple(ii, vec);
         vecVTK2->InsertTuple(ii, vec);
         scaVTK->SetTuple1(ii, 0.0);
         scaVTK2->SetTuple1(ii, 0.0);
       }
    
       //assign nodal coordinates and field data to uGridVTK
       // no need to create lookup table here. All this stuff can be done in Paraview

       uGridVTK->SetPoints(pointsVTK);
       uGridVTK->GetPointData()->SetScalars(scaVTK);
       uGridVTK->GetPointData()->SetVectors(vecVTK);
       uGridVTK->GetPointData()->AddArray(vecVTK2);
       uGridVTK->GetPointData()->AddArray(scaVTK2);
       // create a write object and write uGridVTK to it
    }

    char fname[50];

    sprintf(fname,"%s%s", files.Ofile.asCharArray(),"-Grid2D.vtu");

    writerUGridVTK->SetFileName(fname);

#if VTK_MAJOR_VERSION == 5
    writerUGridVTK->SetInput(uGridVTK);
#else
    writerUGridVTK->SetInputData(uGridVTK);
#endif

    writerUGridVTK->Write();
    
    return;
}
//



void  HBSplineFEM::postProcess1D(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    int ii, jj, kk, ll, span, nlocal, index, ind1, nodenum, lev;

    nlocal = degree[0]+1;

    VectorXd  N(nlocal), dN_dx(nlocal), d2N_dx2(nlocal), tempVec;
    VectorXd  NN(nlocal), dNN_dx(nlocal), d2NN_dx2(nlocal);
    myPoint  knotIncr, knotBegin;

    double   fact, uleft, uright, *tmp;

    vector<double>  outp, u1, u2, outp2, outp3, uu;
    vector<int> bfs;

    node* nd1;
    
    nd1 = elems[degree[0]];
    //nd1 = elems[0];
    AdvDiffExact1D  analy;

    while(nd1 != NULL)
    {
       //cout << " Node # " << nd1->getID() << '\t' << nd1->isGhost() << '\t' << nd1->isLeaf() << endl;
       if( !(nd1->isGhost()) &&  nd1->isLeaf() )
       {
           tmp = nd1->getKnots(0);

           knotBegin = nd1->getKnotBegin();
           knotIncr  = nd1->getKnotIncrement();
           //printf("\t tmp[0] and tmp[1]  ... : %12.8f\t%12.8f\n", tmp[0], tmp[1] );
               
           fact = (tmp[1] - tmp[0])/resln[0];
           //fact = tmp[1] - tmp[0];
           //
           if( CompareDoubles(tmp[1], 1.0) )
             create_vector(tmp[0], tmp[1], fact, uu);
           else
             create_vector(tmp[0], tmp[1] - fact, fact, uu);
           //
           /*
           if( CompareDoubles(tmp[1], 0.0) )
             create_vector(tmp[0], tmp[1] - fact, fact, uu);
           else
             create_vector(tmp[0], tmp[1], fact, uu);
           */

              bfs = nd1->GlobalBasisFuncs;
              for(kk=0;kk<uu.size();kk++)
              {
                 //HB_BasisFuns(degree[0], tmp[0], (tmp[1]-tmp[0]), uu[kk], &N(0), &dN_dx(0), &d2N_dx2(0));

                 param[0] = uu[kk];
                 GeomData.computeBasisFunctions1D(knotBegin, knotIncr, param, NN, dNN_dx, d2NN_dx2);
                 
                 if(nd1->getParent() == NULL)
                 {
                   N = NN;
                   dN_dx = dNN_dx;
                   d2N_dx2 = d2NN_dx2;
                 }
                 else
                 {
                   N = nd1->SubDivMat*NN;
                   dN_dx = nd1->SubDivMat*dNN_dx;
                   d2N_dx2 = nd1->SubDivMat*d2NN_dx2;
                 }

                 param[0] = uu[kk];
                 ComputeGeometry(param, geom);
                 u1.push_back(geom[0]);
                 outp.push_back(nd1->computeValue(0, N));
                 //outp2.push_back(nd1->computeForce(0, N));
                 outp2.push_back(nd1->computeValue(0, dN_dx));
                 //outp2.push_back(nd1->computeValue(0, N));
                 //outp3.push_back(nd1->computeValue(0, d2N_dx2));
                 outp3.push_back(analy.computeValue(0, uu[kk], 0.0));
              }

          while(nd1->getNeighbour(RIGHT) == NULL  && (nd1->getLevel() > 0) )
          {
            nd1 = nd1->getParent();
          }
          nd1 = nd1->getNeighbour(RIGHT);
       }
       else if( !(nd1->isLeaf()) )
       {
            nd1 = nd1->getChild(LEFT);
       }
       else
         nd1 = NULL;
    }

    char fname[50];

    //sprintf(fname,"%s%06d%s", "BSpline-",filecount, ".vtu");
    sprintf(fname,"%s%s%06d%s", files.Ofile.asCharArray(),"-",filecount, ".dat");

    //ofstream fout("advdiff-result-HB.dat");
    
    ofstream fout(fname);

    if(fout.fail())
    {
      cout << " Could not open the Output file" << endl;
      exit(1);
    }

    fout.setf(ios::fixed);
    fout.setf(ios::showpoint);
    fout.precision(8);
   
    for(ii=0;ii<u1.size();ii++)
       fout << u1[ii] << setw(20) << outp[ii] << setw(20) << outp2[ii] << setw(20) << outp3[ii] << endl;

    fout.close();

    //printData(5,1);

    //plotvtk.clearWindow();
    pointsVTK->Reset();
    uGridVTK->Reset();

    vtkIdType pt0, pt1;

    pt0 = pointsVTK->InsertNextPoint(u1[0], outp[0], 0.0);

    for(ii=1;ii<u1.size()-1;ii++)
    {
        pt1 = pointsVTK->InsertNextPoint(u1[ii+1], outp[ii+1], 0.0);
          
        lineVTK->GetPointIds()->SetId(0,pt0);
        lineVTK->GetPointIds()->SetId(1,pt1);
        uGridVTK->InsertNextCell(lineVTK->GetCellType(), lineVTK->GetPointIds());
          
        pt0 = pt1;
    }

    uGridVTK->SetPoints(pointsVTK);

    return;
}


//
// for the general mesh (including Hierarchical splines)
void  HBSplineFEM::postProcess2D(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    //uGridVTK->Reset();
    //pointsVTK->Reset();

    int  dd, ii, jj, kk, ll, count, nlocal, index, ind1, ind2, e, ee, gcount, ind;

    nlocal = (degree[0]+1) * (degree[1] + 1);

    VectorXd  NN(nlocal), N(nlocal), dN_dx(nlocal), dN_dy(nlocal), dNN_dx(nlocal), dNN_dy(nlocal), tempVec, tempVec2, d2N_dx2(nlocal), d2N_dy2(nlocal);
    VectorXd  Du(2), dp(2), vectmp(nlocal), R(2), vel(2), rhsTemp;
    MatrixXd  F(2,2);
    myPoint  knotIncr, knotBegin;

    double   fact, uleft, uright, *tmp1, *tmp2, geom[3], incr1, incr2, val1;

    vector<int> bfs;
    vector<double>  xx, yy, uu, vv;
    vector<vector<double> >  outp2, outp3;

    time_t tstart, tend;

    outp2.resize(ndf+2);
    outp3.resize(ndf+2);

    node* nd1;

    dd = velDOF/ndof;

    rhsTemp.resize(dd);
    rhsTemp.setZero();

    //cout << " Llllllllll " << endl;

    for(e=0;e<activeElements.size();e++)
       elems[activeElements[e]]->RhsToMapResult(1,1, &(rhsTemp(0)));

    //cout << " Llllllllll " << endl;

    //tstart = time(0);
    //SolnData.vorticity = solver2.solve(rhsTemp);
    //tend = time(0);
    //printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );

    //tstart = time(0);

    SolnData.vorticity = solver3.solve(rhsTemp);

    //tend = time(0);
    //printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );
    
    //BiharmonicEx1 analy;
    PoissonEx3  analy;

  if(ndf == 1)
  {
    index = 0;
    for(e=0;e<activeElements.size();e++)
    {
           nd1 = elems[activeElements[e]];
           //cout << " Node # " << nd1->getID() << endl;

           tmp1 = nd1->getKnots(0);
           tmp2 = nd1->getKnots(1);

           knotBegin = nd1->getKnotBegin();
           knotIncr  = nd1->getKnotIncrement();
           //printf("\t tmp[0] and tmp[1]  ... : %12.8f\t%12.8f\n", tmp[0], tmp[1] );
               
           fact = (tmp1[1] - tmp1[0])/resln[0];
           create_vector(tmp1[0], tmp1[1], fact, uu);

           fact = (tmp2[1] - tmp2[0])/resln[1];
           create_vector(tmp2[0], tmp2[1], fact, vv);
           
           //cout << uu << endl;           cout << vv << endl;

           incr1 = tmp1[1] - tmp1[0];
           incr2 = tmp2[1] - tmp2[0];
           
           //create the coordinates of the pointsVTK (nodes in FEM)

           count = 0;
           for(jj=0;jj<vv.size();jj++)
           {
              param[1] = vv[jj];
              for(ii=0;ii<uu.size();ii++)
              {
                 param[0] = uu[ii];
                 GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

                 if(nd1->getParent() == NULL)
                   N = NN;
                 else
                   N = nd1->SubDivMat*NN;

                 fact = nd1->computeValue(0, N);
                 outp2[0].push_back(fact);
                 //outp3[0].push_back(nd1->computeForce(0, N));
                 //outp3[0].push_back(fact - analy.computeValue(0, uu[ii], vv[jj]));
                 //outp3[0].push_back(analy.computeValue(0, uu[ii], vv[jj]));
                 outp3[0].push_back(0.0);

                 count++;
              }
           }
    }

    for(ii=0;ii<outp2[0].size();ii++)
      scaVTK->SetTuple1(ii, outp2[0][ii]);

    for(ii=0;ii<outp3[0].size();ii++)
      scaVTK2->SetTuple1(ii, outp3[0][ii]);

    scaVTK->SetName("value");
    scaVTK2->SetName("force");

    //assign nodal coordinates and field data to uGridVTK
    // no need to create lookup table here. All this stuff can be done in Paraview

    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->AddArray(scaVTK2);
    // create a write object and write uGridVTK to it
  }
  else // for Stokes and Navier-Stokes
  {
    index = 0;
    for(e=0;e<activeElements.size();e++)
    {
           nd1 = elems[activeElements[e]];
           //cout << " Node # " << nd1->getID() << endl;

           tmp1 = nd1->getKnots(0);
           tmp2 = nd1->getKnots(1);

           knotBegin = nd1->getKnotBegin();
           knotIncr  = nd1->getKnotIncrement();
           //printf("\t tmp[0] and tmp[1]  ... : %12.8f\t%12.8f\n", tmp1[0], tmp1[1] );
               
           fact = (tmp1[1] - tmp1[0])/resln[0];
           create_vector(tmp1[0], tmp1[1], fact, uu);

           fact = (tmp2[1] - tmp2[0])/resln[1];
           create_vector(tmp2[0], tmp2[1], fact, vv);
           
           incr1 = tmp1[1] - tmp1[0];
           incr2 = tmp2[1] - tmp2[0];
           
           //create the coordinates of the pointsVTK (nodes in FEM)
           //cout << " ooooooooooooo " << endl;

           for(jj=0;jj<vv.size();jj++)
           {
              param[1] = vv[jj];
              for(ii=0;ii<uu.size();ii++)
              {
                 param[0] = uu[ii];
                 GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy);

                 if(nd1->getParent() == NULL)
                 {
                   N = NN;
                   dN_dx = dNN_dx;
                   dN_dy = dNN_dy;
                 }
                 else
                 {
                   N = nd1->SubDivMat*NN;
                   dN_dx = nd1->SubDivMat*dNN_dx;
                   dN_dy = nd1->SubDivMat*dNN_dy;
                 }

                 for(dd=0;dd<ndof;dd++)
                   outp2[dd].push_back(nd1->computeValue(dd, N));
                 
                 //if(!LSFEM_FLAG)
                   //outp2[2].push_back(nd1->computeValue2(0, N));

                 //fact = nd1->computeValue(1, dN_dx);
                 //fact -= nd1->computeValue(0, dN_dy);
                 
                 if(ndof == 3)
                 {
                   fact = nd1->computeVorticity(N);
                   outp2[3].push_back(fact);
                 }

                 outp3[0].push_back(nd1->computeForce(0, N));
                 outp3[1].push_back(nd1->computeForce(1, N));

                 //outp3[0].push_back(nd1->computeValue(0, dN_dx));
                 //outp3[1].push_back(nd1->computeValue(0, dN_dy));
                 //outp3[2].push_back(nd1->computeValue(1, dN_dx));
                 //outp3[3].push_back(nd1->computeValue(1, dN_dy));
                 //count++;
              }
           }
           //cout << " ooooooooooooo " << endl;
    }

    //cout << " jjjjjjjjjjjjjjjjjj " << endl;
    vecVTK->SetName("vel");
    //vecVTK2->SetName("grad");
    vecVTK2->SetName("force");
    scaVTK->SetName("pres");
    scaVTK2->SetName("vortz");
    
    double vec[3], vec2[3];
    
    vec[2] = 0.0;
    for(ii=0;ii<outp2[0].size();ii++)
    {
      vec[0] = outp2[0][ii];
      vec[1] = outp2[1][ii];
      vecVTK->InsertTuple(ii, vec);
      scaVTK->SetTuple1(ii, outp2[2][ii]);
      scaVTK2->SetTuple1(ii, outp2[3][ii]);
    }

    for(ii=0;ii<outp2[0].size();ii++)
    {
      vec2[0] = outp3[0][ii];
      vec2[1] = outp3[1][ii];
      //vec2[2] = outp3[2][ii];
      //vec2[3] = outp3[3][ii];

      vecVTK2->InsertTuple(ii, vec2);
    }
    //cout << " jjjjjjjjjjjjjjjjjj " << endl;
    
    //assign nodal coordinates and field data to uGridVTK
    // no need to create lookup table here. All this stuff can be done in Paraview

    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->SetVectors(vecVTK);
    uGridVTK->GetPointData()->AddArray(vecVTK2);
    uGridVTK->GetPointData()->AddArray(scaVTK2);
  }


    char fname[50];

    //sprintf(fname,"%s%06d%s", "BSpline-",filecount, ".vtu");
    sprintf(fname,"%s%s%06d%s", files.Ofile.asCharArray(),"-",filecount, ".vtu");

    writerUGridVTK->SetFileName(fname);

#if VTK_MAJOR_VERSION == 5
    writerUGridVTK->SetInput(uGridVTK);
#else
    writerUGridVTK->SetInputData(uGridVTK);
#endif

    writerUGridVTK->Write();

    for(ee=0;ee<ImmersedBodyObjects.size();ee++)
    {
      if( ImmersedBodyObjects[ee]->IsFlexibleBody() )
        ImmersedBodyObjects[ee]->postProcess(filecount);
    }
    
    return;
}
//





void  HBSplineFEM::createPostProcessGrid3D(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    if(resln[0] > 10 || resln[1] > 10 || resln[1] > 10)
    {
      cout << " resolution is too high for this postprocess algorithm "  << endl;
      exit(1);
    }
    
    uGridVTK->Reset();
    pointsVTK->Reset();

    int  dd, ii, jj, kk, ll, count, nlocal, index, ind1, ind2, ind3, e, ee, gcount;
    int  n1, n2, n3, nn, ind4, ind5, ind6;

    nlocal = (degree[0]+1) * (degree[1] + 1) * (degree[2] + 1);

    double   fact, uleft, uright, *tmp1, *tmp2, *tmp3, incr1, incr2, incr3;

    vector<double>  xx, yy, zz, uu, vv, ww;

    node* nd1;

    vtkIdType pt[8];
    
    index = 0;
    for(e=0;e<activeElements.size();e++)
    {
           nd1 = elems[activeElements[e]];

           tmp1 = nd1->getKnots(0);
           tmp2 = nd1->getKnots(1);
           tmp3 = nd1->getKnots(2);

           //printf("\t tmp[0] and tmp[1]  ... : %12.8f\t%12.8f\n", tmp[0], tmp[1] );
               
           fact = (tmp1[1] - tmp1[0])/resln[0];
           create_vector(tmp1[0], tmp1[1], fact, uu);

           fact = (tmp2[1] - tmp2[0])/resln[1];
           create_vector(tmp2[0], tmp2[1], fact, vv);

           fact = (tmp3[1] - tmp3[0])/resln[2];
           create_vector(tmp3[0], tmp3[1], fact, ww);

           //cout << uu << endl;           cout << vv << endl;

           for(kk=0;kk<ww.size();kk++)
           {
             for(jj=0;jj<vv.size();jj++)
             {
               for(ii=0;ii<uu.size();ii++)
               {
                 param[0] = uu[ii]; param[1] = vv[jj]; param[2] = ww[kk];
                 ComputeGeometry(param, geom);

                 xx.push_back(geom[0]);
                 yy.push_back(geom[1]);
                 zz.push_back(geom[2]);

                 pointsVTK->InsertNextPoint(geom[0], geom[1], geom[2]);
               }
             }
           }

           // create the connectivity of elements/nodes
           // use triangle/quadrilateral as required

           n1 = uu.size();
           n2 = vv.size();
           n3 = ww.size();
    
           nn = n1*n2;

           for(kk=0;kk<n3-1;kk++)
           {
               ind5 = index + nn*kk;
               ind6 = index + nn*(kk+1);

               for(jj=0;jj<n2-1;jj++)
               {
                  ind1 = ind5 + n1*jj;
                  ind2 = ind5 + n1*(jj+1);
       
                  ind3 = ind6 + n1*jj;
                  ind4 = ind6 + n1*(jj+1);

                  for(ii=0;ii<n1-1;ii++)
                  {
                     pt[0] = ind1+ii;         pt[4] = ind3+ii;
                     pt[1] = pt[0]+1;         pt[5] = pt[4]+1;
                     pt[3] = ind2+ii;         pt[7] = ind4+ii;
                     pt[2] = pt[3]+1;         pt[6] = pt[7]+1;
              
                     for(ll=0;ll<8;ll++)
                        hexVTK->GetPointIds()->SetId(ll, pt[ll]);
          
                     uGridVTK->InsertNextCell(hexVTK->GetCellType(), hexVTK->GetPointIds());
                  }
               }
           }
           index = pointsVTK->GetNumberOfPoints();
    }

    //cout << count << '\t' << xx.size() << '\t' << yy.size() << endl;

    count = xx.size();

    if(ndf == 1)
    {
      scaVTK->SetName("value");
      scaVTK->SetNumberOfTuples(count);

      scaVTK2->SetName("force");
      scaVTK2->SetNumberOfTuples(count);

      for(ii=0;ii<count;ii++)
      {
        scaVTK->SetTuple1(ii, 0.0);
        scaVTK2->SetTuple1(ii, 0.0);
      }
    
      //assign nodal coordinates and field data to uGridVTK
      // no need to create lookup table here. All this stuff can be done in Paraview

      uGridVTK->SetPoints(pointsVTK);
      uGridVTK->GetPointData()->SetScalars(scaVTK);
    }
    else // for Stokes and Navier-Stokes
    {
       vecVTK->SetName("vel");
       vecVTK2->SetName("vortz");
       scaVTK->SetName("pres");

       vecVTK->SetNumberOfComponents(3);
       vecVTK->SetNumberOfTuples(count);
       vecVTK2->SetNumberOfComponents(3);
       vecVTK2->SetNumberOfTuples(count);
       scaVTK->SetNumberOfTuples(count);
       scaVTK2->SetNumberOfTuples(count);
    
       double vec[3];
    
       vec[0] = vec[1] = vec[2] = 0.0;
       for(ii=0;ii<count;ii++)
       {
         vecVTK->InsertTuple(ii, vec);
         vecVTK2->InsertTuple(ii, vec);
         scaVTK->SetTuple1(ii, 0.0);
         scaVTK2->SetTuple1(ii, 0.0);
       }
    
       //assign nodal coordinates and field data to uGridVTK
       // no need to create lookup table here. All this stuff can be done in Paraview

       uGridVTK->SetPoints(pointsVTK);
       uGridVTK->GetPointData()->SetScalars(scaVTK);
       uGridVTK->GetPointData()->SetVectors(vecVTK);
       uGridVTK->GetPointData()->AddArray(vecVTK2);
       //uGridVTK->GetPointData()->AddArray(scaVTK2);
       // create a write object and write uGridVTK to it
    }

    char fname[50];

    sprintf(fname,"%s%s", files.Ofile.asCharArray(),"-Grid2D.vtu");

    writerUGridVTK->SetFileName(fname);

#if VTK_MAJOR_VERSION == 5
    writerUGridVTK->SetInput(uGridVTK);
#else
    writerUGridVTK->SetInputData(uGridVTK);
#endif

    writerUGridVTK->Write();
    
    return;
}



void  HBSplineFEM::postProcess3D(int vartype, int vardir, int nCol, bool umnxflag, double umin, double umax, int* resln)
{
    //cout << " HBSplineFEM::postProcess3D " << endl;

    uGridVTK->Reset();
    pointsVTK->Reset();

    int  dd, ii, jj, kk, ll, count, nlocal, index, ind1, ind2, e, ee, gcount, ind;

    nlocal = (degree[0]+1) * (degree[1] + 1) * (degree[2] + 1);

    VectorXd  NN(nlocal), dNN_dx(nlocal), dNN_dy(nlocal), dNN_dz(nlocal), d2NN_dx2(nlocal), d2NN_dy2(nlocal);
    VectorXd  N(nlocal), dN_dx(nlocal), dN_dy(nlocal), dN_dz(nlocal), d2N_dx2(nlocal), d2N_dy2(nlocal);
    myPoint  knotIncr, knotBegin;

    double   fact, uleft, uright, *tmp1, *tmp2, *tmp3, geom[3], incr1, incr2, incr3, val1;

    vector<double>  uu, vv, ww;
    vector<vector<double> >  outp2, outp3;

    node* nd1;

    outp2.resize(ndf+2);
    outp3.resize(ndf+2);

  if(ndf == 1)
  {
    index = 0;
    for(e=0;e<activeElements.size();e++)
    {
           nd1 = elems[activeElements[e]];

           //cout << " Node # " << nd1->getID() << endl;

           tmp1 = nd1->getKnots(0);
           tmp2 = nd1->getKnots(1);
           tmp3 = nd1->getKnots(2);

           knotBegin = nd1->getKnotBegin();
           knotIncr  = nd1->getKnotIncrement();
           //printf("\t tmp[0] and tmp[1]  ... : %12.8f\t%12.8f\n", tmp[0], tmp[1] );
               
           fact = (tmp1[1] - tmp1[0])/resln[0];
           create_vector(tmp1[0], tmp1[1], fact, uu);

           fact = (tmp2[1] - tmp2[0])/resln[1];
           create_vector(tmp2[0], tmp2[1], fact, vv);

           fact = (tmp3[1] - tmp3[0])/resln[2];
           create_vector(tmp3[0], tmp3[1], fact, ww);

           //cout << uu << endl;
           //cout << vv << endl;
           //cout << ww << endl;

           incr1 = tmp1[1] - tmp1[0];
           incr2 = tmp2[1] - tmp2[0];
           incr3 = tmp3[1] - tmp3[0];
           
           //create the coordinates of the pointsVTK (nodes in FEM)

           for(kk=0;kk<ww.size();kk++)
           {
             param[2] = ww[kk];
             for(jj=0;jj<vv.size();jj++)
             {
               param[1] = vv[jj];
               for(ii=0;ii<uu.size();ii++)
               {
                 param[0] = uu[ii];
                 GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, NN);

                 if(nd1->getParent() == NULL)
                   N = NN;
                 else
                   N = nd1->SubDivMat*NN;

                 //printVector(N);
                 //cout << kk << '\t' << jj << '\t' <<  ii << endl;
                 //cout << nd1->computeValue(0, N) << endl;
                 outp2[0].push_back(nd1->computeValue(0, N));
                 //outp3[0].push_back(nd1->computeForce(0, N));
                 outp3[0].push_back(0.0);
                 //cout << " aaaaaaaaaaa " << endl;
               }
             }
           }
    }

    for(ii=0;ii<outp2[0].size();ii++)
      scaVTK->SetTuple1(ii, outp2[0][ii]);

    for(ii=0;ii<outp3[0].size();ii++)
      scaVTK2->SetTuple1(ii, outp3[0][ii]);

    scaVTK->SetName("value");
    scaVTK2->SetName("force");

    //assign nodal coordinates and field data to uGridVTK
    // no need to create lookup table here. All this stuff can be done in Paraview

    uGridVTK->GetPointData()->SetScalars(scaVTK);
    //uGridVTK->GetPointData()->AddArray(scaVTK2);
    // create a write object and write uGridVTK to it
  }
  else // for Stokes and Navier-Stokes
  {
    index = 0;
    for(e=0;e<activeElements.size();e++)
    {
           nd1 = elems[activeElements[e]];

           //cout << " Node # " << nd1->getID() << endl;

           tmp1 = nd1->getKnots(0);
           tmp2 = nd1->getKnots(1);
           tmp3 = nd1->getKnots(2);

           knotBegin = nd1->getKnotBegin();
           knotIncr  = nd1->getKnotIncrement();
           //printf("\t tmp[0] and tmp[1]  ... : %12.8f\t%12.8f\n", tmp[0], tmp[1] );
               
           fact = (tmp1[1] - tmp1[0])/resln[0];
           create_vector(tmp1[0], tmp1[1], fact, uu);

           fact = (tmp2[1] - tmp2[0])/resln[1];
           create_vector(tmp2[0], tmp2[1], fact, vv);

           fact = (tmp3[1] - tmp3[0])/resln[2];
           create_vector(tmp3[0], tmp3[1], fact, ww);

           //cout << uu << endl;
           //cout << vv << endl;
           //cout << ww << endl;

           incr1 = tmp1[1] - tmp1[0];
           incr2 = tmp2[1] - tmp2[0];
           incr3 = tmp3[1] - tmp3[0];
           //create the coordinates of the pointsVTK (nodes in FEM)

           for(kk=0;kk<ww.size();kk++)
           {
             param[2] = ww[kk];
             for(jj=0;jj<vv.size();jj++)
             {
               param[1] = vv[jj];
               for(ii=0;ii<uu.size();ii++)
               {
                 param[0] = uu[ii];
                 GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, NN, dNN_dx, dNN_dy, dNN_dz);

                 if(nd1->getParent() == NULL)
                 {
                   N = NN;
                   dN_dx = dNN_dx;
                   dN_dy = dNN_dy;
                   dN_dz = dNN_dz;
                 }
                 else
                 {
                   N = nd1->SubDivMat*NN;
                   dN_dx = nd1->SubDivMat*dNN_dx;
                   dN_dy = nd1->SubDivMat*dNN_dy;
                   dN_dz = nd1->SubDivMat*dNN_dz;
                 }

                 for(dd=0;dd<ndf;dd++)
                   outp2[dd].push_back(nd1->computeValue(dd, N));

                 //outp2[3].push_back(nd1->computeValue2(0, N));

                 //outp3[0].push_back(nd1->computeForce(0, N));
                 //outp3[1].push_back(nd1->computeForce(1, N));

                 outp3[0].push_back(nd1->computeValue(2, dN_dy) - nd1->computeValue(1, dN_dz));
                 outp3[1].push_back(nd1->computeValue(0, dN_dz) - nd1->computeValue(2, dN_dx));
                 outp3[2].push_back(nd1->computeValue(1, dN_dx) - nd1->computeValue(0, dN_dy));

                 //printVector(NN);
                 //cout << kk << '\t' << jj << '\t' <<  ii << endl;
                 //cout << nd1->computeValue(0, NN) << endl;
                 //outp3[0].push_back(nd1->computeForce(0, N));
                 //cout << " aaaaaaaaaaa " << endl;
               }
             }
           }
    }

    vecVTK->SetName("vel");
    vecVTK2->SetName("vortz");
    scaVTK->SetName("pres");
    
    double vec[3], vec2[3];
    
    for(ii=0;ii<outp2[0].size();ii++)
    {
      vec[0] = outp2[0][ii];
      vec[1] = outp2[1][ii];
      vec[2] = outp2[2][ii];
      vecVTK->InsertTuple(ii, vec);
      scaVTK->SetTuple1(ii, outp2[3][ii]);
      //scaVTK2->SetTuple1(ii, outp2[4][ii]);

      vec2[0] = outp3[0][ii];
      vec2[1] = outp3[1][ii];
      vec2[2] = outp3[2][ii];

      vecVTK2->InsertTuple(ii, vec2);
    }
    
    //assign nodal coordinates and field data to uGridVTK
    // no need to create lookup table here. All this stuff can be done in Paraview

    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->SetVectors(vecVTK);
    uGridVTK->GetPointData()->AddArray(vecVTK2);
    //uGridVTK->GetPointData()->AddArray(scaVTK2);
  }

    //if(tis > 10)
      //FluidSolnData.soln = rhsVec;

    char fname[50];

    //sprintf(fname,"%s%06d%s", "BSpline-",filecount, ".vtu");
    sprintf(fname,"%s%s%06d%s", files.Ofile.asCharArray(),"-",filecount, ".vtu");

    writerUGridVTK->SetFileName(fname);

#if VTK_MAJOR_VERSION == 5
    writerUGridVTK->SetInput(uGridVTK);
#else
    writerUGridVTK->SetInputData(uGridVTK);
#endif

    writerUGridVTK->Write();
    
    return;
}





void HBSplineFEM::plotGaussPoints()
{
    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK2;
    vtkSmartPointer<vtkPoints>               pointsVTK2;
    vtkSmartPointer<vtkVertex>               vertexVTK2;

    vtkSmartPointer<vtkFloatArray>           scaVTK3, scaVTK4;

    uGridVTK2     =  vtkSmartPointer<vtkUnstructuredGrid>::New(); 
    pointsVTK2    =  vtkSmartPointer<vtkPoints>::New();
    vertexVTK2    =  vtkSmartPointer<vtkVertex>::New();

    scaVTK3    =  vtkSmartPointer<vtkFloatArray>::New();
    scaVTK4    =  vtkSmartPointer<vtkFloatArray>::New();

    int ee, ll, gp;

    vtkIdType  ptId;
    double  *tmp0, *tmp1, *tmp2;
    node*  nd1;

    param.setZero();

    for(ee=0;ee<activeElements.size();ee++)
    {
      nd1 = elems[activeElements[ee]];

          tmp0 = nd1->getKnots(Dir1);
          tmp1 = nd1->getKnots(Dir2);
          tmp2 = nd1->getKnots(Dir3);

          //getGaussPoints1D(nGP, gausspoints1, gaussweights1);
          //getGaussPoints1D(nGP, gausspoints2, gaussweights2);

          //cout << tmp0[0] << '\t' << tmp0[1] << '\t' << tmp1[0] << '\t' << tmp1[1] << endl;

          for(gp=0; gp<GeomData.gausspoints.size(); gp++)
          {
              param[0]  = 0.5*(tmp0[2] * GeomData.gausspoints[gp][0] + tmp0[3]);
              param[1]  = 0.5*(tmp1[2] * GeomData.gausspoints[gp][1] + tmp1[3]);

              ComputeGeometry(param, geom);

              //cout << param[0] << '\t' << param[1] << endl;
              //cout << geom[0] << '\t' << geom[1] << endl;

              ptId = pointsVTK2->InsertNextPoint(geom[0], geom[1], geom[2]);

              vertexVTK2->GetPointIds()->SetId(0, ptId);

              uGridVTK2->InsertNextCell(vertexVTK2->GetCellType(), vertexVTK2->GetPointIds());
          } // for(gp=0;
    }

    uGridVTK2->SetPoints(pointsVTK2);
    
    //cellDataVTK2->SetName("IsGhost");
    //scaVTK3->SetName("dist");
    //scaVTK4->SetName("InOut");

    //uGridVTK2->GetPointData()->SetScalars(scaVTK3);
    //uGridVTK2->GetPointData()->AddArray(scaVTK4);

    char fname[200];

    sprintf(fname,"%s%s", files.Ofile.asCharArray(),"-Gausspoints.vtu");

    writerUGridVTK->SetFileName(fname);

#if VTK_MAJOR_VERSION == 5
    writerUGridVTK->SetInput(uGridVTK2);
#else
    writerUGridVTK->SetInputData(uGridVTK2);
#endif

    writerUGridVTK->Write();

    return;
}



