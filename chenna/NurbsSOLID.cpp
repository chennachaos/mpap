/*=============================================================================
        File: NurbsSOLID.cpp
  Created by: Chennakesava Kadapa          (08 Jan 2011)
 Purpose    : Implementation file for the definitions of NURBS SOLID Class

 ============================================================================*/


#include "NurbsSOLID.h"
#include "DataBlockTemplate.h"
#include "PlotVTK.h"
#include <iomanip>
#include <assert.h>
#include "QuadratureUtil.h"

#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkPolyDataMapper.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>  
#include <vtkPolygon.h>  
#include <vtkActor2D.h>
#include <vtkHexahedron.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkProperty.h>
#include <vtkFloatArray.h>
#include <vtkLine.h>

using namespace std;


extern PlotVTK plotvtk;
extern MpapTime mpapTime;



NurbsSOLID::NurbsSOLID(CNET3D& Pw1, KNOTVECTOR& U1, KNOTVECTOR& V1, KNOTVECTOR& W1, DEGREE p1, DEGREE q1, DEGREE r1)
{
  Pw = Pw1;
  U  = U1;  V  = V1;  W  = W1;
  p  = p1;  q  = q1;  r  = r1;
}




NurbsSOLID::NurbsSOLID(const NurbsSOLID& solid1)
{
  Pw = solid1.Pw;
  U  = solid1.U;  V  = solid1.V;  W  = solid1.W;
  p  = solid1.p;  q  = solid1.q;  r  = solid1.r;
}





NurbsSOLID::~NurbsSOLID()
{
   //cout << "    NurbsSOLID: destructor ...\n\n";
}




NurbsSOLID& NurbsSOLID::operator=(const NurbsSOLID &rhs)
{
  this->Pw = rhs.Pw;
  this->U  = rhs.U;  this->V  = rhs.V;  this->W  = rhs.W;
  this->p  = rhs.p;  this->q  = rhs.q;  this->r  = rhs.r;

  return *this;
}



void NurbsSOLID::initializeBCdata()
{
   int ii, jj;

  tracflag.setDim(18);
  for(ii=0;ii<18;ii++)
  tracflag[ii] = false;

  VectorArray<double> X1, X2, X3;
  findunique(U, X1);
  findunique(V, X2);
  findunique(W, X3);

  nelem1 = X1.n-1; //no. of elements in first direction
  nelem2 = X2.n-1; //no. of elements in second direction
  nelem3 = X3.n-1; //no. of elements in third direction

  nelem = nelem1*nelem2*nelem3;

  ngbf1 = U.n - p - 1;
  ngbf2 = V.n - q - 1;
  ngbf3 = W.n - r - 1;

  ngbf = ngbf1 * ngbf2 * ngbf3;
  
  nlbf1 = p+1;
  nlbf2 = q+1;
  nlbf3 = r+1;

  nlbf = nlbf1 * nlbf2 * nlbf3;
  
  ngbf1m2  = ngbf1*ngbf2;
  nlbf1m2  = nlbf1*nlbf2;
  nelem1m2 = nelem1*nelem2;
  
  nsize = nlbf*ndof;

  dispBCs.setDim(ngbf);
  gbfnums.setDim(ngbf);

  for(ii=0;ii<ngbf;ii++)
  {
     gbfnums[ii] = -5555;

     dispBCs[ii].setDim(ndof);
     for(jj=0;jj<ndof;jj++)
        dispBCs[ii][jj] = -7777; // some fancy number
  }

  forceBCs = dispBCs;

  // check if the surface is closed

  Uinit.setDim(ndof*ngbf);
  Uinit.zero();

  Values.setDim(ndof);
  ValuesPrev.setDim(ndof);
  ValuesPrev2.setDim(ndof);
  ValuesPrev3.setDim(ndof);
  ValuesPrev4.setDim(ndof);

  ValuesDot.setDim(ndof);
  ValuesDotDot.setDim(ndof);
  ValuesDotPrev.setDim(ndof);
  ValuesDotDotPrev.setDim(ndof);

  ValuesCur.setDim(ndof);
  ValuesDotCur.setDim(ndof);
  ValuesDotDotCur.setDim(ndof);


  for(ii=0;ii<ndof;ii++)
  {
    Values[ii].setDim(ngbf);           Values[ii].zero();
    ValuesPrev[ii].setDim(ngbf);       ValuesPrev[ii].zero();
    ValuesPrev2[ii].setDim(ngbf);      ValuesPrev2[ii].zero();
    ValuesPrev3[ii].setDim(ngbf);      ValuesPrev3[ii].zero();
    ValuesPrev4[ii].setDim(ngbf);      ValuesPrev4[ii].zero();

    ValuesDot[ii].setDim(ngbf);        ValuesDot[ii].zero();
    ValuesDotPrev[ii].setDim(ngbf);    ValuesDotPrev[ii].zero();
    ValuesDotDot[ii].setDim(ngbf);     ValuesDotDot[ii].zero();
    ValuesDotDotPrev[ii].setDim(ngbf); ValuesDotDotPrev[ii].zero();


    ValuesCur[ii].setDim(ngbf);        ValuesCur[ii].zero();
    ValuesDotCur[ii].setDim(ngbf);     ValuesDotCur[ii].zero();
    ValuesDotDotCur[ii].setDim(ngbf);  ValuesDotDotCur[ii].zero();
  }



  toUfull.setDim(ndof*ngbf);
  toUfull.zero();

  dispBCdata.setDim(8);
  dispBCdata.zero();

  edgedata.setDim(6);
  for(ii=0;ii<6;ii++)
     edgedata[ii] = -1;
  intfdata.setDim(12);
  for(ii=0;ii<12;ii++)
      intfdata[ii] = -1;

  face_rot_angle.resize(6);
  for(ii=0;ii<6;ii++)
    face_rot_angle[ii] = 0.0;

   PP.setDim(ngbf3);
   for(ii=0;ii<ngbf3;ii++)
   {
      PP[ii].setDim(ngbf1);
      for(jj=0;jj<ngbf1;jj++)
        PP[ii][jj].setDim(ngbf2);
   }

   face_rotation_flag = false;

   boundary_faces.setDim(6);

     // face 1
     boundary_faces[0].U = U;
     boundary_faces[0].V = V;
     boundary_faces[0].p = p;
     boundary_faces[0].q = q;

     // face 2
     boundary_faces[1].U = U;
     boundary_faces[1].V = V;
     boundary_faces[1].p = p;
     boundary_faces[1].q = q;
        
     // face 3
     boundary_faces[2].U = U;
     boundary_faces[2].V = W;
     boundary_faces[2].p = p;
     boundary_faces[2].q = r;
        
     // face 4
     boundary_faces[3].U = U;
     boundary_faces[3].V = W;
     boundary_faces[3].p = p;
     boundary_faces[3].q = r;
        
     // face 5
     boundary_faces[4].U = W;
     boundary_faces[4].V = V;
     boundary_faces[4].p = r;
     boundary_faces[4].q = q;
        
     // face 6
     boundary_faces[5].U = W;
     boundary_faces[5].V = V;
     boundary_faces[5].p = r;
     boundary_faces[5].q = q;
     
     for(ii=0;ii<6;ii++)
     {
        boundary_faces[ii].ndof = 2;
        boundary_faces[ii].initializeBCdata();
     }

//    computeNET();

//cout << " NurbsSOLID::initializeBCdata() .... DONE " << endl;

  return;
}



void NurbsSOLID::updateValues(int ind, double* uu)
{
  ind -= 1;

  if(ind > ndof)
    cerr << " ERROR! in  NurbsSURFACE::updateValues .... " << endl;

  for(int ii=0;ii<ngbf;ii++)
    Values[ind][ii] += uu[toUfull[ii]];

  return;
}



void NurbsSOLID::updateCoordinates(double* uu)
{
  int index, count, ii, jj, kk;
  double  temp;

  count=0;
  for(kk=0;kk<ngbf3;kk++)
  {
    for(jj=0;jj<ngbf2;jj++)
    {
      for(ii=0;ii<ngbf1;ii++)
      {
        index = gbfnums[count]*ndof;

        temp = Pw[kk][ii][jj].w;

        Pw[kk][ii][jj].x += (uu[index]   * temp);
        Pw[kk][ii][jj].y += (uu[index+1] * temp);
        Pw[kk][ii][jj].z += (uu[index+2] * temp);

        count++;
      }
    }
  }

  computeNET();
  return;
}



void NurbsSOLID::resetGeometry(double* uu)
{
  int index, count, ii, jj, kk;

  double  temp;

  count=0;
  for(kk=0;kk<ngbf3;kk++)
  {
    for(jj=0;jj<ngbf2;jj++)
    {
      for(ii=0;ii<ngbf1;ii++)
      {
        index = gbfnums[count]*ndof;

        temp = Pw[kk][ii][jj].w;

        Pw[kk][ii][jj].x = (uu[index]   * temp);
        Pw[kk][ii][jj].y = (uu[index+1] * temp);
        Pw[kk][ii][jj].z = (uu[index+2] * temp);

        count++;
      }
    }
  }

  computeNET();
  return;
}


void NurbsSOLID::addInitDOFvalues()
{
  int index, count=0, ii, jj, kk;

  double temp;
  
  face_rotation_flag = 1;
  cout << " face_rotation_flag = " << face_rotation_flag << endl;

  if(face_rotation_flag)
  {
    int     ind1, ind2, ind3, side, dir;
    double  disp_val, tot_angle2, theta, cs, sn, val1, val2;

    tot_angle2 = PI*(1.0)/180.0;

    theta = mpapTime.dt * tot_angle2;

    cs = cos(theta);
    sn = sin(theta);

    ind3 = ngbf2-1;

    EPOINT EP;

    for(kk=0;kk<ngbf3;kk++)
    {
      ind1 = ngbf1m2 * (kk+1) - ngbf1 ;
      for(jj=0;jj<ngbf1;jj++)
      {
        ind2 = ind1 + jj;
        EP = Pw[kk][jj][ind3].CalcEuclid();

        //val1 = EP.x * cs - EP.z * sn - EP.x;
        //val2 = EP.x * sn + EP.z * cs - EP.z;

        val1 =  EP.x * cs + EP.z * sn - EP.x;
        val2 = -EP.x * sn + EP.z * cs - EP.z;

        //printf("\t%12.8f\t%12.8f\n\n", val1, val2);

        //temp = Pw[kk][ii][jj].w;                              

        Pw[kk][jj][ind3].x += val1;
        Pw[kk][jj][ind3].z += val2;

        dir=0;
        Uinit[ndof*ind2 + dir] = val1;
        Values[dir][ind2]   = val1;

        dir=2;
        Uinit[ndof*ind2 + dir] = val2;
        Values[dir][ind2]   = val2;
      }
    }
  }
  else
  {
    count=0;
    for(kk=0;kk<ngbf3;kk++)
    {
      for(jj=0;jj<ngbf2;jj++)
      {
        for(ii=0;ii<ngbf1;ii++)
        {
          index = count*ndof;

          temp = mpapTime.dt * Pw[kk][ii][jj].w;
    
          Pw[kk][ii][jj].x += (temp * Uinit[index]);
          Pw[kk][ii][jj].y += (temp * Uinit[index+1]);
          Pw[kk][ii][jj].z += (temp * Uinit[index+2]);

          count++;
        }
      }
    }
  }              

  return;
}


void NurbsSOLID::geomToVector(double* outp)
{
  int  ind1, ind3, ii, jj, kk;
  EPOINT EP;

  ind1=0;
  for(kk=0;kk<ngbf3;kk++)
  {
    for(jj=0;jj<ngbf2;jj++)
    {
      for(ii=0;ii<ngbf1;ii++)
      {
        EP = Pw[kk][ii][jj].CalcEuclid();

        ind3 = 3*gbfnums[ind1];

        outp[ind3]   = EP.x;
        outp[ind3+1] = EP.y;
        outp[ind3+2] = EP.z;
        ind1++;
      }
    }
  }
  return;
}




void NurbsSOLID::computeNET()
{
  int  ii, jj, kk;

  for(kk=0;kk<ngbf3;kk++)
  {
    for(jj=0;jj<ngbf2;jj++)
    {
      for(ii=0;ii<ngbf1;ii++)
        PP[kk][ii][jj] = Pw[kk][ii][jj].CalcEuclid();
    }
  }
  return;
}





void NurbsSOLID::updateCoordsSingleCP(int num, int dir, double val)
{
  int ii, jj, kk, n1, n2;

  n1 = ngbf1*ngbf2;

  kk = num/n1;
  n2 = num%n1;

  ii = n2%ngbf1;
  jj = n2/ngbf1;
   
//   cout << num << '\t' << dir << '\t' << ii << '\t' << jj << '\t' << kk << endl;

  if(dir == 0)
    Pw[kk][ii][jj].x += (val * Pw[kk][ii][jj].w) ;
  else if(dir == 1)
    Pw[kk][ii][jj].y += (val * Pw[kk][ii][jj].w) ;
  else
    Pw[kk][ii][jj].z += (val * Pw[kk][ii][jj].w) ;

  return;
}






CPOINT NurbsSOLID::SolidPoint(double u, double v, double w)
{
  vector<double>  Nu(p+1), Nv(q+1), Nw(r+1);

  int uspan = FindSpan(&(U[0]), U.n, p, u) ;
  int vspan = FindSpan(&(V[0]), V.n, q, v) ;
  int wspan = FindSpan(&(W[0]), W.n, r, w) ;

  BasisFuns(&(U[0]), U.n, p, u, &Nu[0]) ;
  BasisFuns(&(V[0]), V.n, q, v, &Nv[0]) ;
  BasisFuns(&(W[0]), W.n, r, w, &Nw[0]) ;

  int row, col, plane, ind1, ind2, kk, ii, jj;

  row = uspan-p;
  col = vspan-q;
  plane = wspan-r;

  CPOINT  CP1, CP2, Sw;
  Sw = 0.0;

  for(kk=0;kk<=r;kk++)
  {
    CP1 = 0.0;
    ind1 = plane+kk;

    for(jj=0;jj<=q;jj++)
    {
      CP2 = 0.0;
      ind2 = col+jj;
      for(ii=0;ii<=p;ii++)
        CP2 = CP2 + Nu[ii]*Pw[ind1][row+ii][ind2];

        CP1 = CP1 + Nv[jj]*CP2;
    }
    Sw = Sw + Nw[kk]*CP1;
  }

  return Sw;
}




void NurbsSOLID::PlotControlPoints(int col)
{
    int   ii, jj, kk, ll, n1, n2, n3, ind1, ind2, ind3, ind4, ind5, ind6, nn;
    double  color[3];

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkCellArray>  vertices = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();

    EPOINT  EP;

    // Calculate the Points on the solid
    for(kk=0;kk<ngbf3;kk++)
    {
        for(ii=0;ii<ngbf1;ii++)
        {
           for(jj=0;jj<ngbf2;jj++)
           {
              vtkIdType  pid[1];
              
              EP = Pw[kk][ii][jj].CalcEuclid();
              pid[0] = points->InsertNextPoint(EP.x, EP.y, EP.z);
          
              vertices->InsertNextCell(1, pid);
           }
        }
    }

    // Add the points to a polydata
    vtkPolyData  *polydata = vtkPolyData::New();
    polydata->SetPoints( points );
    
    polydata->SetVerts(vertices);
    
    vtkSmartPointer<vtkActor>  actor1 = vtkSmartPointer<vtkActor>::New();

    vtkSmartPointer<vtkDataSetMapper>  mapper1  =  vtkSmartPointer<vtkDataSetMapper>::New();

    //mapper1->SetInputData(polydata);
    mapper1->SetInput(polydata);

    actor1->SetMapper(mapper1);

    getColorValue(col, color);
    actor1->GetProperty()->SetColor(color[0], color[1], color[2]);
    actor1->GetProperty()->SetPointSize(4);

    plotvtk.rendr->AddActor(actor1);

    plotvtk.renWindow->Render();

  return;
}



void NurbsSOLID::PlotElements(int col, bool PLOT_KNOT_LINES, int* resln)
{
   PlotElementsVTK(col, PLOT_KNOT_LINES, resln);
   return;
}



void NurbsSOLID::PlotElementsVTK(int col, bool PLOT_KNOT_LINES, int* resln)
{
    int   ii, jj, kk, ll, n1, n2, n3, ind1, ind2, ind3, ind4, ind5, ind6, nn, count;

/*
    VectorArray<double>  uu, vv, ww;

    if(nelem1 > 20)
      findunique(U, uu);
    else
      create_vector2(U, 10, uu);

    if(nelem2 > 20)
      findunique(V, vv);
    else
      create_vector2(V, 10, vv);

    if(nelem3 > 20)
      findunique(W, ww);
    else
      create_vector2(W, 10, ww);


    vtkSmartPointer<vtkPoints>            points    =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkUnstructuredGrid>  uGrid     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkHexahedron>        hex       =  vtkSmartPointer<vtkHexahedron>::New();
    vtkSmartPointer<vtkDataSetMapper>     mapper1   =  vtkSmartPointer<vtkDataSetMapper>::New();
    vtkSmartPointer<vtkActor>             actor1    =  vtkSmartPointer<vtkActor>::New();

    static vtkIdType pts[8], pt1, pt2;

    EPOINT  EP1, EP2;

    n1 = uu.n;
    n2 = vv.n;
    n3 = ww.n;
    
    nn = n1*n2;

    // Calculate the Points on the solid
    for(kk=0;kk<n3;kk++)
    {
        for(jj=0;jj<n2;jj++)
        {
           for(ii=0;ii<n1;ii++)
           {
              EP1 = SolidPoint(uu[ii], vv[jj], ww[kk]).CalcEuclid();
              points->InsertNextPoint(EP1.x, EP1.y, EP1.z);
           }
        }
    }

    uGrid->SetPoints(points);

    for(kk=0;kk<n3-1;kk++)
    {
        ind5 = nn*kk;
        ind6 = nn*(kk+1);

        for(jj=0;jj<n2-1;jj++)
        {
           ind1 = ind5 + n1*jj;
           ind2 = ind5 + n1*(jj+1);
       
           ind3 = ind6 + n1*jj;
           ind4 = ind6 + n1*(jj+1);

           for(ii=0;ii<n1-1;ii++)
           {
              pts[0] = ind1+ii;          pts[4] = ind3+ii;
              pts[1] = pts[0]+1;         pts[5] = pts[4]+1;
              pts[3] = ind2+ii;          pts[7] = ind4+ii;
              pts[2] = pts[3]+1;         pts[6] = pts[7]+1;
              
              for(ll=0;ll<8;ll++)
                 hex->GetPointIds()->SetId(ll, pts[ll]);
          
              uGrid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
           }
        }
    }


    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName("solidmesh.vtu");
    writer->SetInput(uGrid);

    writer->Write();

    mapper1->SetInputConnection(uGrid->GetProducerPort());

    actor1->SetMapper(mapper1);
    
    actor1->GetProperty()->SetColor(0.0, 1.0, 0.0); //(R,G,B)
    
    actor1->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0); //(R,G,B)
    actor2->GetProperty()->SetEdgeColor(0.0, 0.0, 0.0); //(R,G,B)
    
    actor1->GetProperty()->SetLineWidth(2.0);
    
//    actor1->GetProperty()->EdgeVisibilityOn();

//    plotvtk.rendr->AddActor(actor1);
*/

     int  reslntmp[2];

            // face 1
            CNET Pw1;
            ind1=0;
            Pw1.setDim(ngbf1);
            for(ii=0;ii<ngbf1;ii++)
            {
                Pw1[ii].setDim(ngbf2);        
                for(jj=0;jj<ngbf2;jj++)
                   Pw1[ii][jj] = Pw[ind1][ii][jj];
            }
            boundary_faces[0].Pw = Pw1;
            
            reslntmp[0] = resln[0];
            reslntmp[1] = resln[1];
            boundary_faces[0].computeNET();
            boundary_faces[0].PlotElementsVTK(col, PLOT_KNOT_LINES, reslntmp);

            // face 2
            ind1 = ngbf3-1;
            Pw1.setDim(ngbf1);
            for(ii=0;ii<ngbf1;ii++)
            {
                Pw1[ii].setDim(ngbf2);        
                for(jj=0;jj<ngbf2;jj++)
                   Pw1[ii][jj] = Pw[ind1][ii][jj];
            }
            boundary_faces[1].Pw = Pw1;

            reslntmp[0] = resln[0];
            reslntmp[1] = resln[1];
            boundary_faces[1].computeNET();
            boundary_faces[1].PlotElementsVTK(col, PLOT_KNOT_LINES, reslntmp);


            // face 3
            ind1 = 0;
            Pw1.setDim(ngbf1);
            for(ii=0;ii<ngbf1;ii++)
            {
                Pw1[ii].setDim(ngbf3);        
                for(jj=0;jj<ngbf3;jj++)
                   Pw1[ii][jj] = Pw[jj][ii][ind1];
            }
            boundary_faces[2].Pw = Pw1;

            reslntmp[0] = resln[0];
            reslntmp[1] = resln[2];
            boundary_faces[2].computeNET();
            boundary_faces[2].PlotElementsVTK(col, PLOT_KNOT_LINES, reslntmp);

            // face 4
            ind1 = ngbf2-1;
            Pw1.setDim(ngbf1);
            for(ii=0;ii<ngbf1;ii++)
            {
                Pw1[ii].setDim(ngbf3);        
                for(jj=0;jj<ngbf3;jj++)
                   Pw1[ii][jj] = Pw[jj][ii][ind1];
            }
            boundary_faces[3].Pw = Pw1;

            reslntmp[0] = resln[0];
            reslntmp[1] = resln[2];
            boundary_faces[3].computeNET();
            boundary_faces[3].PlotElementsVTK(col, PLOT_KNOT_LINES, reslntmp);

            // face 5
            ind1 = 0;
            Pw1.setDim(ngbf3);
            for(ii=0;ii<ngbf3;ii++)
            {
                Pw1[ii].setDim(ngbf2);        
                for(jj=0;jj<ngbf2;jj++)
                   Pw1[ii][jj] = Pw[ii][ind1][jj];
            }
            boundary_faces[4].Pw = Pw1;

            reslntmp[0] = resln[2];
            reslntmp[1] = resln[1];
            boundary_faces[4].computeNET();
            boundary_faces[4].PlotElementsVTK(col, PLOT_KNOT_LINES, reslntmp);

            // face 6
            ind1 = ngbf1-1;
            Pw1.setDim(ngbf3);
            for(ii=0;ii<ngbf3;ii++)
            {
                Pw1[ii].setDim(ngbf2);        
                for(jj=0;jj<ngbf2;jj++)
                  Pw1[ii][jj] = Pw[ii][ind1][jj];
            }
            boundary_faces[5].Pw = Pw1;

            reslntmp[0] = resln[2];
            reslntmp[1] = resln[1];
            boundary_faces[5].computeNET();
            boundary_faces[5].PlotElementsVTK(col, PLOT_KNOT_LINES, reslntmp);


     plotvtk.uGrid->SetPoints(plotvtk.points);


//    plotvtk.rendr->ResetCamera();

//    plotvtk.renWindow->Render();

   return;
}


/*
void NurbsSOLID::PlotElementsVTK()
{
    VectorArray<double>  uu, vv, ww;

    findunique(U, uu);
    findunique(V, vv);
    findunique(W, ww);
      
    int   ii, jj, kk, ll, n1, n2, n3, ind1, ind2, ind3, ind4, ind5, ind6, nn;

    // Create points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkCellArray>  vertices = vtkSmartPointer<vtkCellArray>::New();

    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();

    static vtkIdType pts[8];

    EPOINT  EP;

    n1 = uu.n;
    n2 = vv.n;
    n3 = ww.n;
    
    nn = n1*n2;
    
    // Calculate the Points on the solid
    for(kk=0;kk<n3;kk++)
    {
        for(jj=0;jj<n2;jj++)
        {
           for(ii=0;ii<n1;ii++)
           {
              vtkIdType  pid[1];
              
              EP = SolidPoint(uu[ii], vv[jj], ww[kk]).CalcEuclid();
              pid[0] = points->InsertNextPoint(EP.x, EP.y, EP.z);
          
              vertices->InsertNextCell(1, pid);
           }
        }
    }
    
//    cout << "   AAAAAAAAAAAA " << endl;

    vtkSmartPointer<vtkHexahedron>   hex  =  vtkSmartPointer<vtkHexahedron>::New();

    for(kk=0;kk<n3-1;kk++)
    {
        ind5 = nn*kk;
        ind6 = nn*(kk+1);

        for(jj=0;jj<n2-1;jj++)
        {
           ind1 = ind5 + n1*jj;
           ind2 = ind5 + n1*(jj+1);
       
           ind3 = ind6 + n1*jj;
           ind4 = ind6 + n1*(jj+1);

           for(ii=0;ii<n1-1;ii++)
           {
              pts[0] = ind1+ii;          pts[4] = ind3+ii;
              pts[1] = pts[0]+1;         pts[5] = pts[4]+1;
              pts[3] = ind2+ii;          pts[7] = ind4+ii;
              pts[2] = pts[3]+1;         pts[6] = pts[7]+1;
              
              //for(ll=0;ll<8;ll++)   cout << '\t' << pts[ll] ;
              //cout << endl;
              
              for(ll=0;ll<8;ll++)
                 hex->GetPointIds()->SetId(ll, pts[ll]);
              
              polys->InsertNextCell(hex);
              //uGrid->InsertNextCell(hex->GetCellType(), hex->GetPointIds());
           }
        }
    }


    // Add the points to a polydata
    vtkPolyData  *polydata = vtkPolyData::New();
    polydata->SetPoints( points );
    
    polydata->SetVerts(vertices);
    
    polydata->SetPolys(polys);
    polys->Delete();


    vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName("solidmesh.vtp");
    writer->SetInput(polydata);
    writer->Write();
 

    vtkSmartPointer<vtkActor>  actor1 = vtkSmartPointer<vtkActor>::New();

    vtkSmartPointer<vtkDataSetMapper>  mapper1  =  vtkSmartPointer<vtkDataSetMapper>::New();

    mapper1->SetInput(polydata);

     actor1->SetMapper(mapper1);

     plotvtk.rendr->AddActor(actor1);

     plotvtk.renWindow->Render();

   return;
}
*/







int NurbsSOLID::GenerateConnectivityArrays1(int& ntoteqs)
{
    nGP1  = (int) ElemProp.data[0] ;
    nGP2  = (int) ElemProp.data[1] ;
    nGP3  = (int) ElemProp.data[2] ;
    
    nGP = nGP1*nGP2*nGP3;
    nGP1m2 = nGP1*nGP2;
    
    int  e, a, ni, nj, nk, ii, jj, kk, ll, ind, col, pp, qq, rr;
    int  ind1, ind2, ind3;

    double  temp;
    getGaussPoints1D(nGP1, gausspoints1, gaussweights1);
    getGaussPoints1D(nGP2, gausspoints2, gaussweights2);
    getGaussPoints1D(nGP3, gausspoints3, gaussweights3);

    gaussweights.resize(nGP);
    
    ind=0;
    for(kk=0;kk<nGP3;kk++)
    {
       for(jj=0;jj<nGP2;jj++)
       {
          temp = gaussweights3[kk] * gaussweights2[jj];
          for(ii=0;ii<nGP1;ii++)
             gaussweights[ind++] = temp * gaussweights1[ii];
       }
    }

    // INC Array

    INC.setDim(3);
    for(ii=0;ii<3;ii++)
      INC[ii].setDim(ngbf);

    col=0;
    for(kk=0;kk<ngbf3;kk++)
    {
       for(jj=0;jj<ngbf2;jj++)
       {
          for(ii=0;ii<ngbf1;ii++)
          {
             INC[0][col] = ii;
             INC[1][col] = jj;
             INC[2][col] = kk;
             col++;
          }
       }
    }

    // IEN Array

    IEN.setDim(nelem);
    for(ii=0;ii<nelem;ii++)
      IEN[ii].setDim(nlbf);

    
    e = 0;
    for(kk=0;kk<ngbf3;kk++)
    {
       if( !CompareDoubles(W[kk], W[kk+1]) )
       {
           nk = kk - r;
           for(jj=q;jj<ngbf2;jj++)
           {
              if( !CompareDoubles(V[jj], V[jj+1]) )
              {
                 nj = jj - q;
                 for(ii=p;ii<ngbf1;ii++)
                 {
                     if( !CompareDoubles(U[ii], U[ii+1]) )
                     {
                         a=0;
                         ni = ii - p;
                         for(rr=0;rr<=r;rr++)
                         {
                             ind2 = ngbf1m2*(nk+rr);
                             for(qq=0;qq<=q;qq++)
                             {
                                 ind1 = ind2 + ngbf1*(nj+qq) + ni;
                                 for(pp=0;pp<=p;pp++)
                                 {
                                    IEN[e][a] = ind1+pp;
                                    a++;
                                 }
                             }
                         }
                         e++;
                     }
                 }
              }
           }
       }
    }

    // ID Array

    ID.setDim(ndof);
    for(ii=0;ii<ndof;ii++)
      ID[ii].setDim(ngbf);

    for(ii=0;ii<ngbf;ii++)
    {
        for(jj=0;jj<ndof;jj++)
        {
            if( (int) dispBCs[ii][jj] == -7777)
            {
               ID[jj][ii] = ntoteqs;
               ntoteqs++;
            }
            else
               ID[jj][ii] = -1;
        }
    }

  return 0;
}





int  NurbsSOLID::gbfNumFace1(int ii, int jj)
{
  assert(ii < ngbf1);
  assert(jj < ngbf2);
  
  return (jj*ngbf1+ii);
}



int  NurbsSOLID::gbfNumFace2(int ii, int jj)
{
  assert(ii < ngbf1);
  assert(jj < ngbf2);
  
  return (ngbf1m2*(ngbf3-1) + jj*ngbf1 + ii);
}


int  NurbsSOLID::gbfNumFace3(int ii, int jj)
{
  assert(ii < ngbf1);
  assert(jj < ngbf3);
  
  return (ngbf1m2 * jj + ii);
}


int  NurbsSOLID::gbfNumFace4(int ii, int jj)
{
  assert(ii < ngbf1);
  assert(jj < ngbf3);
  
  return (ngbf1m2 * (jj+1) - ngbf1 + ii);
}


int  NurbsSOLID::gbfNumFace5(int ii, int jj)
{
  assert(ii < ngbf3);
  assert(jj < ngbf2);
  
  return (ngbf1 * jj + ngbf1m2 * ii);
}


int  NurbsSOLID::gbfNumFace6(int ii, int jj)
{
  assert(ii < ngbf3);
  assert(jj < ngbf2);
  
  return (ngbf1 * (jj+1) - 1 + ngbf1m2 * ii);
}


void  NurbsSOLID::computeGBFnumbers(int index, int&  count)
{
  if(index == 1)
  {
    int  ii, jj, kk, ind = 0;

    for(kk=0;kk<ngbf3;kk++)
    {
      for(jj=0;jj<ngbf2;jj++)
      {
        for(ii=0;ii<ngbf1;ii++)
          gbfnums[ind++] = count++;
      }
    }
  }

  if(index == 2)
  {
    for(int ii=0;ii<ngbf;ii++)
    {
      if(gbfnums[ii] == -5555)
        gbfnums[ii] = count++;
    }
  }

  return;
}


void  NurbsSOLID::computeToUfull()
{
    int  ii, index1, index2, jj;

    for(ii=0;ii<ngbf;ii++)
    {
        index1 = gbfnums[ii]*ndof;
        index2 = ii*ndof;
        for(jj=0;jj<ndof;jj++)
           toUfull[index2+jj] = index1 + jj;
    }

   return;
}





 


