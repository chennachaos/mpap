
#include <iostream>

#include "FunctionsProgram.h"
#include "FunctionsEssGrp.h"
#include "Mesh.h"
#include "Plot.h"
#include "PropertyTypeEnum.h"
#include "MathGeom.h"


extern Plot plot;


using namespace std;






void StructuredQuad2D::addCentroidAndArea(double *a, Vector<double> &x)
{
  double xx[8] = { x[ix[0]+ix[0]-2], x[ix[0]+ix[0]-1],
                   x[ix[1]+ix[1]-2], x[ix[1]+ix[1]-1],
                   x[ix[2]+ix[2]-2], x[ix[2]+ix[2]-1],
                   x[ix[3]+ix[3]-2], x[ix[3]+ix[3]-1] },
         area1 = triangleArea2D(xx,xx+2,xx+4),
         area2 = triangleArea2D(xx+4,xx+6,xx);

  a[0] += .3333333 * (area1 * (xx[0]+xx[2]+xx[4]) + area2 * (xx[4]+xx[6]+xx[0]));
  a[1] += .3333333 * (area1 * (xx[1]+xx[3]+xx[5]) + area2 * (xx[5]+xx[7]+xx[1]));

  a[2] += area1 + area2; 

  //cout << a[0] << "," << a[1] << "," << a[2] << "\n";

  return;
}





bool StructuredQuad2D::split(StructuredQuad2D *q0, 
                             int n0,
                             int &numnp, 
                             StructuredQuad2D *q, 
                             StructuredEdge2D *e, 
                             Vector<double> &x, 
                             List<StructuredQuad2D> &quad,
                             List<StructuredEdge2D> &edge,
                             Vector<StructuredEdge2D*> &n2e)
{
  int i, j, k;

  quad.add(new StructuredQuad2D);

  StructuredQuad2D &newQuad = quad[quad.n-1];

  if (q != NULL) { i = 0; while (q->iq[i] != this) i++; q->iq[i] = &newQuad; }

  x.append(0.); x.append(0.); numnp++;

  i = 0; while (ix[i] != n0) i++;

  newQuad.ix[0] = numnp - 1;
  newQuad.iq[0] = this;
  newQuad.ie[0] = NULL;

  i++; if (i > 3) i = 0; 

  newQuad.ix[1] = ix[i]; 
  newQuad.iq[1] = q;
  newQuad.ie[1] = e;
  ix[i] = numnp - 1; 

  i++; if (i > 3) i = 0;

  newQuad.ix[2] = ix[i];
  newQuad.iq[2] = iq[i];
  newQuad.ie[2] = ie[i];
  if (iq[i] != NULL) { j = 0; while (iq[i]->iq[j] != this) j++; iq[i]->iq[j] = &newQuad; }
  ix[i] = numnp;
  iq[i] = &newQuad;
  ie[i] = NULL;

  i++; if (i > 3) i = 0;

  newQuad.ix[3] = numnp;
  newQuad.iq[3] = iq[i];
  newQuad.ie[3] = ie[i];

  n2e.append(ie[i]);

  /*int l;
  for (k=0; k<quad.n; k++)
  {
    cout << k+1 << ":"; for (l=0; l<4; l++) cout << " " << quad[k].ix[l]; 
    cout << ":"; for (l=0; l<4; l++)
    {
      if (quad[k].iq[l] == NULL) j = -1; 
      else { j = 0; while (&(quad[j]) != quad[k].iq[l]) j++; }
      cout << " " << j+1;
    }
    cout << ":"; for (l=0; l<4; l++)
    {
      if (quad[k].ie[l] == NULL) j = -1; 
      else { j = 0; while (&(edge[j]) != quad[k].ie[l]) j++; }
      cout << " " << j+1;
    }
    cout << "\n";
  }*/

  if (iq[i] != q0)
  {
    if (iq[i] != NULL) iq[i]->split(q0,ix[i],numnp,&newQuad,ie[i],x,quad,edge,n2e);

    return false;
  }

  StructuredEdge2D *e0;
  StructuredQuad2D *q1;

  int n1;

  x.del(x.n-1); x.del(x.n-1); n2e.del(n2e.n-1); numnp--; 

  k = 0; for (j=0; j<4; j++) k = max(k,q0->ix[j]);

  j = 0; while (k-1 != q0->ix[j]) j++;

  n1 = q0->ix[j]; 

  e0 = q0->ie[j];

  j++; if (j > 3) j = 0;
 
  q1 = q0->iq[j];
  
  i--; if (i < 0) i = 3;  ix[i] = n1;

  newQuad.ix[3] = n1;

  newQuad.ie[3] = e0;

  newQuad.iq[3] = q1;

  j = 0; while (q1->ix[j] != n1) j++; j++; if (j > 3) j = 0;

  q1->iq[j] = &newQuad;

  q1->ie[j] = e0;

  /*for (k=0; k<quad.n; k++)
  {
    cout << k+1 << ":"; for (l=0; l<4; l++) cout << " " << quad[k].ix[l]; 
    cout << ":"; for (l=0; l<4; l++)
    {
      if (quad[k].iq[l] == NULL) j = -1; 
      else { j = 0; while (&(quad[j]) != quad[k].iq[l]) j++; }
      cout << " " << j+1;
    }
    cout << ":"; for (l=0; l<4; l++)
    {
      if (quad[k].ie[l] == NULL) j = -1; 
      else { j = 0; while (&(edge[j]) != quad[k].ie[l]) j++; }
      cout << " " << j+1;
    }
    cout << "\n";
  }*/

  return true;
}






bool subdivideStructuredMesh2D(Vector<double> &x,
                               List<StructuredQuad2D> &quad,
                               List<StructuredEdge2D> &edge, int nen, bool keep)
{
  char fct[] = "subdivideStructuredMesh2D";

  int c, i, j, k, q, numnp = x.n / 2, qj, qk, jn, nj0, nj1, kn, nk0, nk1,
      m, m0, m1, n0, n1;

  Vector<StructuredEdge2D*> node2edge;

  Vector<int> tmp, *node2quad = new Vector<int> [numnp];

  bool goAgain, loopClosed;

  double fact, a[3], xmn[2] = { 1.e+20, 1.e+20 }, xmx[2] = { -1.e+20, -1.e+20 }, tol2, *DAT;

  cout << "      generate structured 2D mesh:\n\n";
  COUT << "generate connectivity ...\n\n";

  // generate node to quad connectivity

  for (q=0; q<quad.n; q++)

    for (i=0; i<4; i++) 
    {
      node2quad[quad[q].ix[i]-1].append(q);
      quad[q].iq[i] = NULL;
      quad[q].ie[i] = NULL;
    }

  // generate quad to quad connectivity

  for (i=0; i<numnp; i++)
  {
    for (j=0; j<node2quad[i].n; j++)
    {
      qj = node2quad[i][j];

      nj0 = quad[qj].ix[3];

      for (jn=0; jn<4; jn++)
      {
        nj1 = quad[qj].ix[jn];

        for (k=j+1; k<node2quad[i].n; k++)
        {
          qk = node2quad[i][k];

          nk0 = quad[qk].ix[3];

          for (kn=0; kn<4; kn++)
          {
            nk1 = quad[qk].ix[kn];

            if (nk1 == nj0 && nk0 == nj1) 
            {
              quad[qj].iq[jn] = &(quad[qk]); 
              quad[qk].iq[kn] = &(quad[qj]);
            }
            nk0 = nk1;
          }
        }
        nj0 = nj1;
      }
    }
  }

  // generate quad to edge connectivity

  for (i=0; i<edge.n; i++)
  {
    n0 = edge[i].ix[0];
    n1 = edge[i].ix[1];

    for (j=0; j<node2quad[n0-1].n; j++)
    {
      q = node2quad[n0-1][j];
      m0 = quad[q].ix[3];
      for (k=0; k<4; k++)
      {
        m1 = quad[q].ix[k];
        if ((m0 == n0 && m1 == n1) || (m0 == n1 && m1 == n0)) quad[q].ie[k] = &edge[i];
        m0 = m1;
      }
    }
  }

  delete [] node2quad;
 
  // generate missing boundary edges

  for (q=0; q<quad.n; q++)
  {
    for (i=0; i<4; i++)
    {
      if (quad[q].iq[i] == NULL && quad[q].ie[i] == NULL) 
      {
        j = i - 1; if (j < 0) j = 3;
        edge.add(new StructuredEdge2D);
        edge[edge.n-1].ix[0]  = quad[q].ix[j];
        edge[edge.n-1].ix[1]  = quad[q].ix[i];
        edge[edge.n-1].ix[2]  = 0;
        edge[edge.n-1].dat[0] = 0.;
        quad[q].ie[i] = &edge[edge.n-1];
      }
      else if (quad[q].iq[i] == NULL && quad[q].ie[i] != NULL)
      {
        if (quad[q].ie[i]->ix[2] < 0) 
          { prgWarning(2,fct,"boundary edge requires ix[2] > 0!"); return false; }
      }
    }
  }

  // generate missing inter element edges if required

  if (keep)
  {
    for (q=0; q<quad.n; q++)
    {
      for (i=0; i<4; i++)
      {
        if (quad[q].iq[i] != NULL && quad[q].ie[i] == NULL)
        {
          j = i - 1; if (j < 0) j = 3;
          edge.add(new StructuredEdge2D);
          edge[edge.n-1].ix[0]  = quad[q].ix[j];
          edge[edge.n-1].ix[1]  = quad[q].ix[i];
          edge[edge.n-1].ix[2]  = 0;
          edge[edge.n-1].dat[0] = 0.;
          quad[q].ie[i] = &edge[edge.n-1];
          m = quad[q].iq[i]->ix[3]; 
          k = 0; while (!(m==quad[q].ix[j] && quad[q].iq[i]->ix[k]==quad[q].ix[i])
                     && !(m==quad[q].ix[i] && quad[q].iq[i]->ix[k]==quad[q].ix[j]))
           { m = quad[q].iq[i]->ix[k]; k++; }
          quad[q].iq[i]->ie[k] = &edge[edge.n-1];
        }
      }
    }
  }

  // move inter element edges back

  for (q=0; q<quad.n; q++)
  {
    for (i=0; i<4; i++)
    {
      if (quad[q].iq[i] != NULL && quad[q].ie[i] != NULL) edge.add(edge.takeOut(quad[q].ie[i]));
    }
  }

  // move edges with ix[2] < 0 back

  j = edge.n; i = 0;

  while (i < j) if (edge[i].ix[2] < 0) { edge.add(edge.takeOut(i)); j--; } else i++;

  // output quad to node, quad to quad and quad to edge connectivities

  /*for (q=0; q<quad.n; q++)
  {
    cout << q+1 << " :"; for (i=0; i<4; i++) cout << " " << quad[q].ix[i]; 
    cout << " :"; for (i=0; i<4; i++)
    {
      if (quad[q].iq[i] == NULL) j = -1; else { j = 0; while (&(quad[j]) != quad[q].iq[i]) j++; }
      cout << " " << j+1;
    }
    cout << " :"; for (i=0; i<4; i++)
    {
      if (quad[q].ie[i] == NULL) j = -1; else { j = 0; while (&(edge[j]) != quad[q].ie[i]) j++; }
      cout << " " << j+1;
    }
    cout << "\n";
  }*/

  // initialise node to edge connectivity

  for (i=0; i<numnp; i++) node2edge.append(NULL);

  // split quads

  for (i=0; i<edge.n; i++)
  {
    c = 0; for (j=0; j<numnp; j++) if (node2edge[j] == &edge[i]) c++;

    if (c == 0)
    {
      q = 0; while (q < quad.n)
      {
        j = 0; while (j < 4) if (quad[q].ix[j] == edge[i].ix[0]) break; else j++;

        if (j < 4)
        {
          k = 0; while (k < 4) if (quad[q].ix[k] == edge[i].ix[1]) break; else k++;

          if (k < 4)
          {
            if (abs(k-j) == 1) k = min(k,j);
            else if (abs(k-j) == 3) k = 3;
            else { prgWarning(2,fct,"invalid edge!"); return false; }

            m = k+1; if (m > 3) m = 0;
            
            k = quad[q].ix[k];

            break;
          }
        }
        q++;
      }
      if (q == quad.n) { prgWarning(3,fct,"invalid edge!"); return false; }

      for (j=0; j<abs(edge[i].ix[2])-1; j++)
      {
        x.append(0.); x.append(0.); numnp++; node2edge.append(&edge[i]);
        loopClosed = quad[q].split(&quad[q],k,numnp,NULL,&edge[i],x,quad,edge,node2edge);
        if (j==0)
        {
          if (quad[q].iq[m] != NULL && !loopClosed) 
            { prgWarning(1,fct,"boundary edge missing to complement internal edge!"); return false; }
        }
      }
    }
    else if (c + 1 != abs(edge[i].ix[2]) && edge[i].ix[2] != 0) 
      { prgWarning(1,fct,"inconsistency in edge subdivisions!"); return false; }
  }

  //for (i=0; i<numnp; i++)
  //{
  //  cout << i + 1;
  //  if (node2edge[i] != NULL) { j = 0; while (node2edge[i] != &edge[j]) j++; cout << " " << j+1; }
  //  cout << "\n";
  //} 

  // generate node to node connectivity

  Vector<int> *node2node = new Vector<int> [numnp];

  for (q=0; q<quad.n; q++)
  {
    for (j=0; j<4; j++) 
    {
      jn = quad[q].ix[j];
      k = j+1; if (k > 3) k = 0;
      {
        kn = quad[q].ix[k];
        if (!node2node[jn-1].contains(kn))
          { node2node[jn-1].append(kn); node2node[kn-1].append(jn); }
      }
    }
  }

  // correct node to node connectivity for nodes on edges and position nodes on edges

  COUT << "position nodes on edges ...\n\n";

  for (i=edge.n-1; i>-1; i--)
  {
    if (edge[i].ix[2] > -1)
    {
      DAT = edge[i].dat;

      tmp.free();
   
      n0 = edge[i].ix[0] - 1;
      n1 = edge[i].ix[1] - 1;

      node2edge[n0] = &edge[i];
      node2edge[n1] = &edge[i];

      for (j=0; j<2; j++)
      {
        xmx[j] = max(xmx[j],x[n0+n0+j]);
        xmn[j] = min(xmn[j],x[n0+n0+j]);
      }

      tol2 = dist2(xmx,xmn,2) * 1.e-14;

      for (j=0; j<numnp; j++)
      {
        if (node2edge[j] == &edge[i] && j != n0 && j != n1)
        {
          tmp.append(j);

          k = 0; while (k < node2node[j].n)

            if (node2edge[node2node[j][k]-1] != &edge[i]) node2node[j].del(k); else k++;
        }
      }
   
      goAgain = true;

      if (DAT[0] > 1.e-10) // transform to cylindrical coordinates
      {
        circleThroughTwoPoints(DAT+1,*DAT,x[n0+n0],x[n0+n0+1],x[n1+n1],x[n1+n1+1]);
        
        x[n0+n0]   = DAT[0];
        x[n0+n0+1] = DAT[3];
        x[n1+n1]   = DAT[0];
        x[n1+n1+1] = DAT[4];

        for (j=0; j<tmp.n; j++)
        {
          k        = tmp[j];
          x[k+k]   = DAT[0];
          x[k+k+1] = DAT[3];
        }
      }
      while (goAgain)
      {
        goAgain = false;

        for (j=0; j<tmp.n; j++)
        {
          k = tmp[j];
          n0 = (node2node[k][0]-1); n0 += n0;
          n1 = (node2node[k][1]-1); n1 += n1; k += k;
          a[0] = .5 * (x[n0] + x[n1]);
          a[1] = .5 * (x[n0+1] + x[n1+1]);
          if ((a[0]-x[k])*(a[0]-x[k])+(a[1]-x[k+1])*(a[1]-x[k+1])>tol2) goAgain = true;
          x[k]   = a[0];
          x[k+1] = a[1];
        }
      }
      if (DAT[0] > 1.e-10) // transform to cartesian coordinates
      {
        tmp.append(edge[i].ix[0]-1);
        tmp.append(edge[i].ix[1]-1);

        for (j=0; j<tmp.n; j++)
        {
          k = tmp[j];
          x[k+k]   = DAT[1] + DAT[0] * cos(x[k+k+1]);
          x[k+k+1] = DAT[2] + DAT[0] * sin(x[k+k+1]);
        }
      }
    }
    else for (j=0; j<numnp; j++) if (node2edge[j] == &edge[i]) node2edge[j] = NULL;
  }

  // regenerate node to quad connectivity

  node2quad = new Vector<int> [numnp];

  for (q=0; q<quad.n; q++) for (i=0; i<4; i++) node2quad[quad[q].ix[i]-1].append(q);

  // position internal nodes

  COUT << "position internal nodes (1/2) ...\n\n";

  for (q=0; q<300; q++)
  {
    for (i=0; i<numnp; i++)
    {
      if (node2edge[i] == NULL)
      {
        a[0] = 0.;
        a[1] = 0.;
        for (j=0; j<node2node[i].n; j++)
        {
          k = node2node[i][j] - 1;
          a[0] += x[k+k];
          a[1] += x[k+k+1];
        }
        fact = 1./ (double) node2node[i].n;
        x[i+i]   = a[0] * fact;
        x[i+i+1] = a[1] * fact;
      }
    }
  }
  delete [] node2node;

  COUT << "position internal nodes (2/2) ...\n\n";

  goAgain = true; q = 0; tol2 *= 1.e+3;

  while (goAgain && q<1000)
  {
    goAgain = false; q++;

    for (i=0; i<numnp; i++)
    {
      if (node2edge[i] == NULL)
      {
        a[0] = 0.;
        a[1] = 0.;
        a[2] = 0.;
        for (j=0; j<node2quad[i].n; j++) quad[node2quad[i][j]].addCentroidAndArea(a,x);
        if (a[2] > 1.e-12) 
        {
          fact = 1./ a[2];
          a[0] *= fact;
          a[1] *= fact;
          if ((a[0]-x[i+i])*(a[0]-x[i+i])+(a[1]-x[i+i+1])*(a[1]-x[i+i+1])>tol2) goAgain = true;
          x[i+i]   += 1.2*(a[0]-x[i+i]);
          x[i+i+1] += 1.2*(a[1]-x[i+i+1]);
        }
        else goAgain = true;
      }
    }
  }
  delete [] node2quad;

  if (nen > 4)
  {
    COUT << "generate mid nodes ...\n\n";

    for (q=0; q<quad.n; q++) for (i=4; i<9; i++) quad[q].ix[i] = 0;

    for (q=0; q<quad.n; q++)
    {
      for (i=0; i<4; i++)
      {
        if (quad[q].ix[4+i] == 0)
        {
          k = i+1; if (k > 3) k = 0; 

          n0 = quad[q].ix[i] - 1;
          n1 = quad[q].ix[k] - 1;

          x.append(.5 * (x[n0+n0]  +x[n1+n1]));
          x.append(.5 * (x[n0+n0+1]+x[n1+n1+1]));

          quad[q].ix[4+i] = ++numnp;

          if (quad[q].ie[k] != NULL)

            if (quad[q].ie[k]->dat[0] > 1.e-10 && quad[q].ie[k]->ix[2] > 0) 
            {
              DAT    = quad[q].ie[k]->dat;
              j      = numnp + numnp - 2;
              fact   = DAT[0] / sqrt((x[j]-DAT[1])*(x[j]-DAT[1])+(x[j+1]-DAT[2])*(x[j+1]-DAT[2]));
              x[j]   = DAT[1] + (x[j]  -DAT[1])*fact;
              x[j+1] = DAT[2] + (x[j+1]-DAT[2])*fact;
            }

          if (quad[q].iq[k] != NULL)
          {
            j = 0; while (quad[q].iq[k]->iq[j] != &quad[q]) j++; j--; if (j < 0) j = 3;

            quad[q].iq[k]->ix[4+j] = numnp;
          }
        }
      }
    }
    if (nen > 8)
    {
      for (q=0; q<quad.n; q++)
      {
        a[0] = 0.;
        a[1] = 0.;
        for (i=0; i<4; i++)
        {
          k = quad[q].ix[i] - 1;
          a[0] += x[k+k];
          a[1] += x[k+k+1];
        }
        x.append(.25 * a[0]);
        x.append(.25 * a[1]);

        quad[q].ix[8] = ++numnp;
      }
    }
  }

  COUT << "The structured mesh has been generated.\n\n";

  return true;
}








bool Mesh::import(int fmt, char *fileName, int ndmIn, int ndfIn)
{
  MyString lineStrg, *word;

  char fct[] = "Mesh::import";

  List< Vector<int> > iTmp;
  List< Vector<double> > dTmp;

  Vector<int> tmp;
  VectorInfinite<int> tmp2;

  int i, j, nw, ii, nst, i0, nenEl, elType;

  bool gmshOld, gmshMessy, keep;

  List<StructuredQuad2D>   quad;
  List<StructuredEdge2D>   edge;

  Vector<double> xTmp;

  ndm = ndmIn;

  ndf = ndfIn;

  // open file 

  std::ifstream in;
	
  in.open(fileName);

  if (!in) { prgWarning(1,fct,"could not open mesh file!"); return false; }
 
  // read file 

  switch (fmt)
  {
    case 1: // netGen

            if (ndm == 2) 
            { 
              while (in && lineStrg.getNextLine(in) != "surfaceelementsgi") continue; 
              i0 = 4;
            }              

            if (ndm == 3) 
            { 
              while (in && lineStrg.getNextLine(in) != "volumeelements") continue;
              i0 = 1;
            }

            if (!in) { prgWarning(2,fct,"error in mesh file / wrong ndm?"); return false; }

            lineStrg.getNextLine(in);  // jump over line containing numel

            nen = 0;

            COUT << "number of elements =         ";

            while (1)
            {
              lineStrg.getNextLine(in);

              nw = lineStrg.split(&word);
 
              if (nw > 1)
              {
                iTmp.add(new Vector<int>);

                if (!word[i0].toInt(&j,false)) { iTmp.del(iTmp.n-1); break; } 

                if (nw < j + i0 + 1)           { iTmp.del(iTmp.n-1); break; }

                nen = max(nen,j);

                i = 0; while (i < j)
                  if (word[i0+i+1].toInt(iTmp[iTmp.n-1].append(),false)) i++; else break;
                if (i < j) { iTmp.del(iTmp.n-1); break; }

                for (i=0; i<nw; i++) word[i].free(); delete [] word;

                if (nen == 4 && ndm == 3) 
                { 
                  i                 = iTmp[iTmp.n-1][2];
                  iTmp[iTmp.n-1][2] = iTmp[iTmp.n-1][3];
                  iTmp[iTmp.n-1][3] = i;
                }
              }
              else break;

              printf("\b\b\b\b\b\b\b\b\b%9d",iTmp.n);
            }
            cout << "\n\n";

            for (i=0; i<nw; i++) word[i].free(); delete [] word;

            if (!(ndm == 3 && nen == 4) && !(ndm == 2 && nen == 3))
            {
              prgWarning(3,fct,"I don't know how to handle this (ndm,nen) in the netgen format");
              return false;
            }

            if (ndm == 2) 
            { 
              while (in && lineStrg.getNextLine(in) != "volumeelements") continue;
              if (!in) { prgWarning(3,fct,"error in mesh file / wrong ndm?"); return false; }
              lineStrg.getNextLine(in);
              nw = lineStrg.split(&word);
              word[0].toInt(&j);
              if (j != 0) { prgWarning(4,fct,"this is a 3D mesh!"); return false; }
              for (i=0; i<nw; i++) word[i].free(); delete [] word;
            }

            numel = iTmp.n;

            ixTmp.setDim(numel,nen+1,true);   ixTmp.zero();

            for (i=0; i<numel; i++) for (j=0; j<nen; j++) ixTmp(i+1,j+1) = iTmp[i][j];

            while (lineStrg.getNextLine(in) != "points") continue;

            lineStrg.getNextLine(in);  // jump over line containing numnp

            COUT << "number of nodes =         ";

            while (1)
            {
              lineStrg.getNextLine(in);

              nw = lineStrg.split(&word);  if (nw < ndm) break;

              if (nw > 1)
              {
                dTmp.add(new Vector<double>);

                i = 0; while (i < ndm) 
                  if (word[i].toDbl(dTmp[dTmp.n-1].append(),false)) i++; else break;
                if (i < ndm) { dTmp.del(dTmp.n-1); break; }

                for (i=0; i<nw; i++) word[i].free(); delete [] word;
              }
              else break;

              printf("\b\b\b\b\b\b\b\b\b%9d",dTmp.n);
            }
            cout << "\n\n";

            for (i=0; i<nw; i++) word[i].free(); delete [] word;

            numnp = dTmp.n;

            x.setDim(numnp,ndm,true);

            for (i=0; i<numnp; i++) for (j=0; j<ndm; j++) x(i+1,j+1) = dTmp[i][j];

            break;

    case 2: // Gmsh

            while (in)
            {
              lineStrg.getNextLine(in);
              if (lineStrg == "$Nodes" || lineStrg == "$NOD") break;
            }

            if (!in) { prgWarning(10,fct,"error in file"); return false; }

            if (lineStrg == "$NOD") gmshOld = true; else gmshOld = false;

            lineStrg.getNextLine(in);  // jump over line containing numnp

            COUT << "number of nodes =         ";

            while (1)
            {
              lineStrg.getNextLine(in);

              nw = lineStrg.split(&word);

              if (nw < 4) break;

              else
              {
                if (!word[0].toInt(tmp.append())) break;

                dTmp.add(new Vector<double>);

                i = 0; while (i < ndm) 
                  if (word[i+1].toDbl(dTmp[dTmp.n-1].append(),false)) i++; else break;
                if (i < ndm) { dTmp.del(dTmp.n-1); break; }

                for (i=0; i<nw; i++) word[i].free(); delete [] word;
              }
              printf("\b\b\b\b\b\b\b\b\b%9d",dTmp.n);
            }
            cout << "\n\n";

            for (i=0; i<nw; i++) word[i].free(); delete [] word;

            numnp = dTmp.n;

            x.setDim(numnp,ndm,true);

            for (i=0; i<numnp; i++) for (j=0; j<ndm; j++) x(i+1,j+1) = dTmp[i][j];

            i = 0; while (i < numnp) if (tmp[i] == i+1) i++; else break;

            if (i == numnp) gmshMessy = false; else gmshMessy = true;

            if (gmshMessy) for (i=0; i<tmp.n; i++) tmp2[tmp[i]] = i + 1;

            while (in) 
            { 
              lineStrg.getNextLine(in);
              if (lineStrg == "$Elements" || lineStrg == "$ELM") break;
            }

            if (!in) { prgWarning(11,fct,"error in file"); return false; }

            lineStrg.getNextLine(in);  // jump over line containing numel

            nen = 0;

            COUT << "number of elements =         ";

            while (1)
            {
              lineStrg.getNextLine(in);

              nw = lineStrg.split(&word);
 
              if (nw > 4)
              {
                if (!word[1].toInt(&elType,false)) break;

                if ((elType == 2 && ndm == 2) || (elType == 4 && ndm == 3)) // 3 noded triangle or
                                                                            // 4 noded tetrahedron
                {
                  if      (elType == 4) nenEl = 4;
                  else if (elType == 2) nenEl = 3;

                  nen = max(nen,nenEl);

                  iTmp.add(new Vector<int>);

                  if (!gmshOld)
                  {
                    if (!word[2].toInt(&i0,false)) { iTmp.del(iTmp.n-1); break; } 

                    i0 += 3;
                  }
                  else i0 = 5;

                  if (nw < i0 + nenEl) { iTmp.del(iTmp.n-1); break; }

                  i = 0; while (i < nenEl)
                    if (word[i0+i].toInt(iTmp[iTmp.n-1].append(),false)) i++; else break;
                  if (i < nenEl) { iTmp.del(iTmp.n-1); break; }

                  if (gmshMessy)
                    for (i=0; i<nenEl; i++) { iTmp[iTmp.n-1][i] = tmp2[iTmp[iTmp.n-1][i]]; }

                  for (i=0; i<nenEl; i++) 
                    if (iTmp[iTmp.n-1][i] < 1 || iTmp[iTmp.n-1][i] > numnp)
                    { 
                      cout << iTmp[iTmp.n-1][i] << "\n";
                      prgError(1,fct,"error in file");
                    }                     
                }
                //else if (!(elType ==  1 || elType ==  3
                // || elType ==  8 || elType ==  9 || elType == 10
                // || elType == 15 || elType == 16 || (elType > 19 && elType < 29)))
                //{
                //  prgWarning(1,fct,"import of this element type not yet implemented"); 
                //
                //  return false;
                //}
                //else // point, line or area element
              }
              else break;

              for (i=0; i<nw; i++) word[i].free(); delete [] word;

              printf("\b\b\b\b\b\b\b\b\b%9d",iTmp.n);
            }
            cout << "\n\n";

            for (i=0; i<nw; i++) word[i].free(); delete [] word;

            numel = iTmp.n;

            ixTmp.setDim(numel,nen+1,true);   ixTmp.zero();

            for (i=0; i<numel; i++) for (j=0; j<iTmp[i].n; j++) ixTmp(i+1,j+1) = iTmp[i][j];

            //cout << numel << "\n" << ixTmp << "\n";

            break;

    case 3: // mpap2 / mesh

            while (in)
            {
              lineStrg.getNextLine(in);
              if (lineStrg.begins("coordinates")) break;
            }

            if (!in) { prgWarning(10,fct,"error in file"); return false; }

            COUT << "number of nodes =         ";

            while (1)
            {
              lineStrg.getNextLine(in);

              nw = lineStrg.split(&word);

              if (nw != ndm + 1) break;

              dTmp.add(new Vector<double>);

              i = 0; while (i < ndm) 
                if (word[i+1].toDbl(dTmp[dTmp.n-1].append(),false)) i++; else break;
              if (i<ndm) { dTmp.del(dTmp.n-1); break; }

              for (i=0; i<nw; i++) word[i].free(); delete [] word;

              printf("\b\b\b\b\b\b\b\b\b%9d",dTmp.n);
            }
            cout << "\n\n";

            for (i=0; i<nw; i++) word[i].free(); delete [] word;

            numnp = dTmp.n;

            x.setDim(numnp,ndm,true);

            for (i=0; i<numnp; i++) for (j=0; j<ndm; j++) x(i+1,j+1) = dTmp[i][j];

            while (in) 
            { 
              if (lineStrg.begins("elements")) break;
              lineStrg.getNextLine(in);
            }

            if (!in) { prgWarning(110,fct,"error in file"); return false; }

            nen = 0;

            COUT << "number of elements =         ";

            while (1)
            {
              lineStrg.getNextLine(in);

              nw = lineStrg.split(&word);
 
              if (nw < 3) break;

              nenEl = nw - 1;

              nen = max(nen,nenEl);

              iTmp.add(new VectorInfinite<int>);

              i = 0; while (i < nenEl)
                if (word[1+i].toInt(&iTmp[iTmp.n-1][i],false)) i++; else break;
              if (i < nenEl) { iTmp.del(iTmp.n-1); break; }

              for (i=0; i<nenEl; i++) 
                if (iTmp[iTmp.n-1][i] < 1 || iTmp[iTmp.n-1][i] > numnp)
                { 
                  cout << iTmp[iTmp.n-1][i] << "\n";
                  prgError(1,fct,"error in file");
                }
                    
              for (i=0; i<nw; i++) word[i].free(); delete [] word;

              printf("\b\b\b\b\b\b\b\b\b%9d",iTmp.n);
            }
            cout << "\n\n";

            for (i=0; i<nw; i++) word[i].free(); delete [] word;

            numel = iTmp.n;

            ixTmp.setDim(numel,nen+1,true);   ixTmp.zero();

            for (i=0; i<numel; i++) for (j=0; j<iTmp[i].n; j++) ixTmp(i+1,j+1) = iTmp[i][j];

            //cout << numel << "\n" << ixTmp << "\n";

            break;

    case 4: // mpap2 / structured 2D

            if (ndm != 2) { prgWarning(10,fct,"ndm != 2"); return false; }

            while (in)
            {
              lineStrg.getNextLine(in);
              if (lineStrg.begins("nodes per element")) break;
            }

            if (!in) { prgWarning(5,fct,"error in file"); return false; }

            lineStrg.getNextLine(in);

            nw = lineStrg.split(&word);

            if (nw != 1) { prgWarning(6,fct,"error in file"); return false; }

            if (!word[0].toInt(&nen,false)) { prgWarning(7,fct,"error in file"); return false; }

            if (nen != 4 && nen != 8 && nen != 9)
              { prgWarning(8,fct,"choose 3, 4, 8 or 9 nodes per element!"); return false; }

            while (in)
            {
              lineStrg.getNextLine(in);
              if (lineStrg.begins("coordinates")) break;
            }

            if (!in) { prgWarning(10,fct,"error in file"); return false; }

            while (1)
            {
              lineStrg.getNextLine(in);

              nw = lineStrg.split(&word);

              if (nw != ndm + 1) break;

              i = 0; while (i < ndm) 
                if (word[i+1].toDbl(xTmp.append(),false)) i++; else break;
              if (i < ndm) { xTmp.del(xTmp.n-1); break; }

              for (i=0; i<nw; i++) word[i].free(); delete [] word;
            }
            for (i=0; i<nw; i++) word[i].free(); delete [] word;

            numnp = xTmp.n / 2;

            while (in) 
            { 
              if (lineStrg.begins("quads")) break;
              lineStrg.getNextLine(in);
            }

            if (!in) { prgWarning(115,fct,"error in file"); return false; }

            while (1)
            {
              lineStrg.getNextLine(in);

              nw = lineStrg.split(&word);
 
              if (nw != 5) break;

              quad.add(new StructuredQuad2D);

              i = 0; while (i < 4)
                if (word[1+i].toInt(quad[quad.n-1].ix+i,false)) i++; else break;
              if (i < 4) { quad.del(quad.n-1); break; }

              for (i=0; i<4; i++) 
                if (quad[quad.n-1].ix[i] < 1 || quad[quad.n-1].ix[i] > numnp)
                { 
                  cout << quad[quad.n-1].ix[i] << "\n";
                  prgWarning(1,fct,"invalid node number in 'quads'");
                  return false;
                }
                    
              if (word != NULL) { for (i=0; i<nw; i++) word[i].free(); delete [] word; }
              word = NULL;
            }
            if (word != NULL) { for (i=0; i<nw; i++) word[i].free(); delete [] word; } word = NULL;

            while (in) 
            { 
              if (lineStrg.begins("lines") || lineStrg.begins("circles") || 
                  lineStrg.begins("keep")) break;
              lineStrg.getNextLine(in);
            }

            if (lineStrg.begins("lines"))
            {
              while (1)
              {
                lineStrg.getNextLine(in);
           
                if (!in) break;
 
                nw = lineStrg.split(&word);
            
                if (nw != 3) break;
            
                edge.add(new StructuredEdge2D); edge[edge.n-1].dat[0] = 0.;
            
                i = 0; while (i < 3)
                  if (word[i].toInt(edge[edge.n-1].ix+i,false)) i++; else break;
                if (i < 3) { edge.del(edge.n-1); break; }
            
                for (i=0; i<2; i++) 
                  if (edge[edge.n-1].ix[i] < 1 || edge[edge.n-1].ix[i] > numnp)
                  { 
                    cout << edge[edge.n-1].ix[i] << "\n";
                    prgWarning(1000,fct,"invalid node number in 'line'");
                    return false;
                  }
                      
                if (word != NULL) { for (i=0; i<nw; i++) word[i].free(); delete [] word; }
                word = NULL;
              }
              if (word != NULL) { for (i=0; i<nw; i++) word[i].free(); delete [] word; }

              if (edge.n < 1) { prgWarning(1001,fct,"no lines specified!"); return false; }

              while (in) 
              { 
                if (lineStrg.begins("circles") || lineStrg.begins("keep")) break;
                lineStrg.getNextLine(in);
              }
            }
            if (lineStrg.begins("circles"))
            {
              while (1)
              {
                lineStrg.getNextLine(in);
           
                if (!in) break;
 
                nw = lineStrg.split(&word);
            
                if (nw != 4) break;
            
                edge.add(new StructuredEdge2D);
            
                i = 0; while (i < 2)
                  if (word[i].toInt(edge[edge.n-1].ix+i,false)) i++; else break;
                if (i < 2) { edge.del(edge.n-1); break; }
                if (!word[3].toInt(edge[edge.n-1].ix+2,false)) { edge.del(edge.n-1); break; }
            
                if (!word[2].toDbl(edge[edge.n-1].dat,false)) { edge.del(edge.n-1); break; }
            
                for (i=0; i<2; i++) 
                  if (edge[edge.n-1].ix[i] < 1 || edge[edge.n-1].ix[i] > numnp)
                  { 
                    cout << edge[edge.n-1].ix[i] << "\n";
                    prgWarning(1010,fct,"invalid node number in 'circles'");
                    return false;
                  }
                      
                if (word != NULL) { for (i=0; i<nw; i++) word[i].free(); delete [] word; }
                word = NULL;
              }
              if (word != NULL) { for (i=0; i<nw; i++) word[i].free(); delete [] word; }

              if (edge.n < 1) { prgWarning(1020,fct,"no circles specified!"); return false; }

              while (in) 
              { 
                if (lineStrg.begins("keep")) break;
                lineStrg.getNextLine(in);
              }
            }

            keep = false;

            if (lineStrg.begins("keep")) keep = true;

            if (!subdivideStructuredMesh2D(xTmp,quad,edge,nen,keep)) return false;

            numnp = xTmp.n / 2;

            x.setDim(numnp,ndm,true);

            for (i=0; i<numnp+numnp; i++) x.x[i] = xTmp[i];

            numel = quad.n;

            ixTmp.setDim(numel,nen+1,true);   ixTmp.zero();

            for (i=0; i<numel; i++) for (j=0; j<nen; j++) ixTmp(i+1,j+1) = quad[i].ix[j];

            //cout << x << ixTmp << "\n";

            break;

    default: cout << "      unknown mesh format!\n\n"; return false;
  }

  // close file

  in.close();

  // generate missing stuff

  completeData();

  // prepare data, finalise mesh

  prepareInputData();

  cout << "   The file " << fileName << " has been read and the mesh imported.\n\n";

  return true;
}












void generateSimple2DBlockOfQuads(int nen, int n1, int n2, int *ix, double *x)
{
  int i, j, k, n, nn, *IX = ix;

  double *X = x, dxi, deta, xi, eta;

  dxi  = 2./(double)n1;
  deta = 2./(double)n2;
  eta  = -1.;
  for (j=0; j<=n2; j++)
  {
    xi   = -1.;
    for (i=0; i<=n1; i++)
    {
      *X = xi;  X++;
      *X = eta; X++;
      xi += dxi;
    }
    eta += deta;
  }

  for (j=0; j<n2; j++)
    for (i=0; i<n1; i++)
    {
      *IX = j*(n1+1)+i+1;     IX++;
      *IX = j*(n1+1)+i+2;     IX++;
      *IX = (j+1)*(n1+1)+i+2; IX++;
      *IX = (j+1)*(n1+1)+i+1; IX++;
      for (k=4; k<nen; k++)   IX++;
      *IX = 0;                IX++;
    }

  if (nen == 4) return;

  eta = -1.;
  for (j=0; j<=n2; j++)
  {
    xi   = -1.+dxi*.5;
    for (i=0; i<n1; i++)
    {
      *X = xi;  X++;
      *X = eta; X++;
      xi += dxi;
    }
    eta += deta;
  }

  eta = -1.+deta*.5;
  for (j=0; j<n2; j++)
  {
    xi   = -1.;
    for (i=0; i<=n1; i++)
    {
      *X = xi;  X++;
      *X = eta; X++;
      xi += dxi;
    }
    eta += deta;
  }

  IX = ix;
  n  = (n1+1) * (n2+1);
  nn = n + n1*(n2+1);
  for (j=0; j<n2; j++)
  {
    for (i=0; i<n1; i++)
    {
      for (k=0; k<4; k++)   IX++;
      *IX = ++n;            IX++;
      *IX = nn+2;           IX++;
      *IX = n + n1;         IX++;
      *IX = ++nn;           IX++; 
      for (k=7; k<nen; k++) IX++;
    }
    nn++;
  }

  if (nen == 8) return;

  eta = -1.+deta*.5;
  for (j=0; j<n2; j++)
  {
    xi   = -1.+dxi*.5;
    for (i=0; i<n1; i++)
    {
      *X = xi;  X++;
      *X = eta; X++;
      xi += dxi;
    }
    eta += deta;
  }

  IX = ix;
  n  = (n1+1)*(n2+1) + n1*(n2+1) + (n1+1)*n2;
  for (j=0; j<n2; j++)
  {
    for (i=0; i<n1; i++)
    {
      for (k=0; k<8; k++) IX++;
      *IX = ++n;          IX++; IX++;
    }
  }

  return;
}









void generateSimple2DRingOfQuads(int nen, int n1, int n2, int *ix, double *x)
{
  int i, j, k, n, nn, *IX = ix;

  double *X = x, dxi, deta, xi, eta;

  dxi  = 1./(double)n1;
  deta = 1./(double)n2;
  eta  = 0.;
  for (j=0; j<=n2; j++)
  {
    xi   = 0.;
    for (i=0; i<n1; i++)
    {
      *X = xi;  X++;
      *X = eta; X++;
      xi += dxi;
    }
    eta += deta;
  }

  for (j=0; j<n2; j++)
  {
    for (i=0; i<n1-1; i++)
    {
      *IX = j*n1+i+1;       IX++;
      *IX = j*n1+i+2;       IX++;
      *IX = (j+1)*n1+i+2;   IX++;
      *IX = (j+1)*n1+i+1;   IX++;
      for (k=4; k<nen; k++) IX++;
      *IX = 0;              IX++;
    }
    *IX = (j+1)*n1;         IX++;
    *IX = j*n1+1;           IX++;
    *IX = (j+1)*n1+1;       IX++;
    *IX = (j+2)*n1;         IX++;
    for (k=4; k<nen; k++)   IX++;
    *IX = 0;                IX++;
  }
  if (nen == 4) return;

  eta = 0.;
  for (j=0; j<=n2; j++)
  {
    xi   = dxi*.5;
    for (i=0; i<n1; i++)
    {
      *X = xi;  X++;
      *X = eta; X++;
      xi += dxi;
    }
    eta += deta;
  }

  eta = deta*.5;
  for (j=0; j<n2; j++)
  {
    xi   = 0.;
    for (i=0; i<n1; i++)
    {
      *X = xi;  X++;
      *X = eta; X++;
      xi += dxi;
    }
    eta += deta;
  }

  IX = ix;
  n  = n1 * (n2+1);
  nn = n + n1*(n2+1);
  for (j=0; j<n2; j++)
  {
    for (i=0; i<n1-1; i++)
    {
      for (k=0; k<4; k++)   IX++;
      *IX = ++n;            IX++;
      *IX = nn+2;           IX++;
      *IX = n + n1;         IX++;
      *IX = ++nn;           IX++; 
      for (k=7; k<nen; k++) IX++;
    }
    for (k=0; k<4; k++)   IX++;
    *IX = ++n;            IX++;
    *IX = nn+2-n1;        IX++;
    *IX = n + n1;         IX++;
    *IX = ++nn;           IX++; 
    for (k=7; k<nen; k++) IX++;
  }

  if (nen == 8) return;

  eta = deta*.5;
  for (j=0; j<n2; j++)
  {
    xi   = dxi*.5;
    for (i=0; i<n1; i++)
    {
      *X = xi;  X++;
      *X = eta; X++;
      xi += dxi;
    }
    eta += deta;
  }

  IX = ix;
  n  = (n1+n1)*(n2+1) + n1*n2;
  for (j=0; j<n2; j++)
  {
    for (i=0; i<n1; i++)
    {
      for (k=0; k<8; k++) IX++;
      *IX = ++n;          IX++; IX++;
    }
  }

  return;
}









bool Mesh::generateSimple2DMesh(int geom, int elType, int ndfIn, int *iData, double *dData)
{
  int e, i, j, n1, n2, nst, *IX;

  double *X, *X0, eta, xi, N1, N2, N3, N4, rad, alp, t = 1.e-10, ri, ro;

  ndf = ndfIn;

  ndm = 2;

  Vector<double> R;

  Vector<int> tmp1, tmp2;

  n1 = roundToInt(iData[0]);
  n2 = roundToInt(iData[1]);

  if      (elType == 1)
  { 
    nen = 4; 
    if (geom < 3)       numnp = (n1+1)*(n2+1);
    else if (geom == 3) numnp = n1*(n2+1);
  }
  else if (elType == 2) 
  {
    nen = 8;
    if (geom < 3)       numnp = (n1+1)*(n2+1) + (n1+1)*n2 + n1*(n2+1);
    else if (geom == 3) numnp = (n1+n1)*(n2+1) + n1*n2;
  }
  else if (elType == 3)
  {
    nen = 9;
    if (geom < 3)       numnp = (n1+1)*(n2+1) + (n1+1)*n2 + n1*(n2+1) + n1*n2;
    else if (geom == 3) numnp = (n1+n1)*(n2+1) + n1*n2 + n1*n2;
  }

  numel = n1 * n2;

  x.setDim(numnp,ndm,true);

  ixTmp.setDim(numel,nen+1,true);

  if (geom < 3)
  {
    generateSimple2DBlockOfQuads(nen,n1,n2,ixTmp.x,x.x);

    for (i=0; i<numnp; i++)
    {
      xi  = x.x[i+i];
      eta = x.x[i+i+1];
      N1 = .25*(1.-xi)*(1.-eta);
      N2 = .25*(1.+xi)*(1.-eta);
      N3 = .25*(1.+xi)*(1.+eta);
      N4 = .25*(1.-xi)*(1.+eta);
      x.x[i+i]   = N1*dData[0] + N2*dData[2] + N3*dData[4] + N4*dData[6];
      x.x[i+i+1] = N1*dData[1] + N2*dData[3] + N3*dData[5] + N4*dData[7];
    }
    if (geom == 2)
    {
      for (i=0; i<numnp; i++)
      {
        rad = x.x[i+i]; 
        alp = x.x[i+i+1]*0.017453293;
        x.x[i+i]   = rad * cos(alp);
        x.x[i+i+1] = rad * sin(alp);
      }
    }
  }
  else if (geom == 3)
  {
    generateSimple2DRingOfQuads(nen,n1,n2,ixTmp.x,x.x);
 
    for (i=0; i<numnp; i++)
    {
      rad = x.x[i+i+1]*(dData[1]-dData[0]) + dData[0]; 
      alp = -x.x[i+i]*6.2831853;
      x.x[i+i]   = rad * cos(alp) + dData[2];
      x.x[i+i+1] = rad * sin(alp) + dData[3];
    }
  }

  // generate missing stuff

  completeData();

  // prepare data, finalise mesh

  prepareInputData();

  return true;
}









 
void Mesh::interactiveNodeSelection(int XX, int YY, bool deselectFlag)
{
  char fct[] = "Mesh::interactiveNodeSelection";

  int XXX[2] = { XX, YY }, i, imn, j, m, k, mNG = -1, mNB = -1;

  double xp[3], *X = x0.x, d2, dmn, cone, axis[3], da, dn, l, obs[3];

  float picCoor[2];

  Vector<int> &nodeList = selectNode.nodeList, tmp, thisList;
  
  if (XX < 0)
  {
    if (selectNode.prnt) cout << "    " << nodeList.n << " selected nodes: " << nodeList << "\n\n";

    return;
  }

  if (selectNode.defm) X = x.x;

  if (ndm == 2)
  {
  // find node from mouse coordinates in 2D

    plot.xy2DInv(xp,XXX);

    if (selectNode.flag == 2) 
    {
      dmn = dist2(xp,X,2);

      imn = 0;

      for (i=1; i<numnp; i++)
      {
        d2 = dist2(xp,X+i+i,2);

        if (d2 < dmn) { imn = i; dmn = d2; }
      }
      mNG = imn+1;
    }
    else
    {
      dmn = dist2(xp,X+(bndNd[0]-1)*ndm,2);
      imn = 0;

      for (i=0; i<nBndNd; i++)
      {
        d2 = dist2(xp,X+(bndNd[i]-1)*ndm,2);

        if (d2 < dmn) { imn = i; dmn = d2; }
      }
      mNB = imn;
      mNG = bndNd[imn];
    }
  }
  else if (ndm == 3)
  {
  // find node from mouse coordinates in 3D

    plot.xy2DInv(picCoor,XXX);

    plot.perspective.searchCone(&cone,axis,picCoor,plot.selectSearchRadius);

    dmn = 1.e20;

    obs[0] = (double)plot.perspective.observer[0];
    obs[1] = (double)plot.perspective.observer[1];
    obs[2] = (double)plot.perspective.observer[2];

    for (i=0; i<nBndNd; i++)
    {
      pointAndLine3D(obs,axis,X+(bndNd[i]-1)*ndm,&da,&dn,&l);

      da *= l;

      if (dn/da < cone)
      {
        if (surf3D->nodeIsReallyVisible(i))

          if (da > 0. && dmn > da) { imn  = i; dmn = da; }
      }
    }
    if (dmn < 1.e19) { mNB = imn; mNG = bndNd[imn]; }
  }

  if (mNG < 0) { COUT << "no nodes selected!\n\n"; return; }

  if (selectNode.flag == 1)
  {

  // geometry mode: find other nodes on the same surface / edge

    for (i=0; i<nBndNd; i++) nodeFlag[bndNd[i]-1] = false;

    tmp.append(mNB);

    while (tmp.n > 0)
    {
      m = bndNd[tmp[0]];

      if (!nodeFlag[m-1]) 
      {
        if (!deselectFlag) 
          { if (!nodeList.contains(m))   { nodeList.append(m); thisList.append(m); } }
        else 
          { if (nodeList.contains(m,&j)) { nodeList.del(j);    thisList.append(m); } }

        nodeFlag[m-1] = true;

        if (ndm == 3) surf3D->getNeighboursOnSameGeometry(tmp[0],tmp,(float)selectNode.cosAlpha);

        else          getNeighboursOnSameGeometry2D(tmp[0],tmp);
      }
      tmp.del(0);
    }
  }
  else 
  {

  // single nodes: add node to list

    if (!deselectFlag) 
      { if (!nodeList.contains(mNG)) { nodeList.append(mNG); thisList.append(mNG); } }

    else
     { if (nodeList.contains(mNG,&j)) { nodeList.del(j); thisList.append(mNG); } }
  }

  // show selected nodes

  if (!deselectFlag) plot.setColour(GREEN); else plot.setColour(RED);

  if (ndm == 3)
  {
    tmp.free(); for (i=0; i<numnp+1; i++) tmp.append(0);

    for (i=0; i<nBndNd; i++) tmp[bndNd[i]] = i;

    if (selectNode.numb)

      for (i=0; i<thisList.n; i++)

        surf3D->plotNodePoint(tmp[thisList[i]],plot.dPt(),thisList[i]);

    else

      for (i=0; i<thisList.n; i++)

        surf3D->plotNodePoint(tmp[thisList[i]],plot.dPt());
  }
  else if (ndm == 2)
  {
    if (selectNode.numb)

      for (i=0; i<thisList.n; i++) plot.point(X+(thisList[i]-1)*ndm,plot.dPt(),thisList[i]);

    else

      for (i=0; i<thisList.n; i++) plot.point(X+(thisList[i]-1)*ndm,plot.dPt());
  }

  essGrpUpdateDisplay();

  if (!plot.suppressCopyPixmap) essGrpCopyPixmap();

  return;
}










void Mesh::simpleNodeSelection(void)
{
  char fct[] = "Mesh::simpleNodeSelection";

  int i, j;

  double *X = x0.x, dd, alpha, *SNX = selectNode.x, t = 1.e-10;

  bool deselectFlag = selectNode.dslt;

  Vector<int> &nodeList = selectNode.nodeList, newList;

  if (selectNode.defm) X = x.x;

  if (selectNode.flag == 3) // select all
  {
    nodeList.free();

    if (!deselectFlag) for (i=1; i<=numnp; i++) nodeList.append(i);
  }
  else if (selectNode.flag == 4) // select all inside box
  {
    for (i=0; i<numnp; i++) 
    {
      j = 0; 
      while (j<ndm)
        if (X[i*ndm+j] < SNX[j] || X[i*ndm+j] > SNX[j+ndm]) break; else j++;
      if (j == ndm) 
      {
        if (!deselectFlag) { if (!nodeList.contains(i+1)) newList.append(i+1); }

        else if (nodeList.contains(i+1,&j)) nodeList.del(j);
      }
    }
  }
  else if (selectNode.flag == 5) // select all inside cylinder
  {
    for (i=0; i<numnp; i++)
    {
      if (ndm == 2) pointOnEdge2D(SNX,SNX+2,X+i+i,&alpha,&dd);
      else          pointOnEdge3D(SNX,SNX+3,X+i+i+i,&alpha,&dd);

      if (alpha > -t && alpha < 1.+t && dd < *(SNX+ndm+ndm)+t)
      {
        if (!deselectFlag) { if (!nodeList.contains(i+1)) newList.append(i+1); }

        else if (nodeList.contains(i+1,&j)) nodeList.del(j);
      }
    }
  }
  else if (selectNode.flag == 6) // select all inside sphere
  {
    *(SNX+ndm) = *(SNX+ndm) * *(SNX+ndm);

    for (i=0; i<numnp; i++) 
    {
      if (dist2(SNX,X+i*ndm,ndm)<*(SNX+ndm)+t)
      {
        if (!deselectFlag) { if (!nodeList.contains(i+1)) newList.append(i+1); }

        else if (nodeList.contains(i+1,&j)) nodeList.del(j);
      }
    }
  }
  else prgError(1,fct,"What do you want?");

  if (!deselectFlag) for (i=0; i<newList.n; i++) nodeList.append(newList[i]);

  // print current node selection

  if (selectNode.prnt) cout << "    " << nodeList.n << " selected nodes: " << nodeList << "\n\n";

  return;
}









void Mesh::setBoundaryConditions(bool *fix, bool remove)
{
  int i, j, *IDU = idu.x;

  Vector<int> &nodeList = selectNode.nodeList;

  if (!remove)
  {
    for (j=0; j<ndf; j++)
      if (fix[j]) for (i=0; i<nodeList.n; i++) IDU[(nodeList[i]-1)*ndf+j] = 0;
  }
  else
  {
    for (j=0; j<ndf; j++)
      if (fix[j]) for (i=0; i<nodeList.n; i++) IDU[(nodeList[i]-1)*ndf+j] = 1;
  }
  return;
}









void Mesh::fixNodeMotion(bool *fix, bool remove)
{
  int i, j, *IDX = idx.x;

  Vector<int> &nodeList = selectNode.nodeList;

  if (!remove)
  {
    for (j=0; j<ndm; j++)
      if (fix[j]) for (i=0; i<nodeList.n; i++) IDX[(nodeList[i]-1)*ndm+j] = 0;
  }
  else
  {
    for (j=0; j<ndm; j++)
      if (fix[j]) for (i=0; i<nodeList.n; i++) IDX[(nodeList[i]-1)*ndm+j] = 1;
  }
  return;
}










void Mesh::setFreeNodes(int flag)
{
  Vector<int> &nodeList = selectNode.nodeList;

  int i, j, bnd0;

  if (flag == 1)
  {
    for (i=0; i<nodeList.n; i++)
      if (!tmpFreeNode.contains(nodeList[i])) tmpFreeNode.append(nodeList[i]);

    return;
  }
  else if (flag == 2)
  {
    for (i=0; i<nodeList.n; i++)
      if (tmpFreeNode.contains(nodeList[i],&j)) tmpFreeNode.del(j);

    return;
  }

  char fct[] = "Mesh::setFreeNodes";

  double *X = x0.x;

  std::ofstream out;

  out.open("mpap2.interface");

  if (!out) { prgWarning(1,fct,"could not open file for writing!"); return; }

  if (tmpFreeNode.n == 0) { COUT << "tmpFreeNode.n == 0 -> nothing to be done!\n\n"; return; }

  for (i=0; i<tmpFreeNode.n; i++)
  {
    for (j=0; j<ndm; j++) out << " " << X[ndm*(tmpFreeNode[i]-1)+j];
    out << "\n";
  }
  out << "===\n";

  VectorArray<int> all2free;

  all2free.setDim(numnp);

  for (i=0; i<numnp; i++) all2free[i] = 0;

  for (i=0; i<tmpFreeNode.n; i++) all2free[tmpFreeNode[i]-1] = i+1;

  if (ndm == 3) surf3D->writeInterfaceInterpolations(all2free,out);

  else if (ndm == 2) 
  {
    prgWarning(1,fct,"ndm = 2; THIS HAS NOT BEEN TESTED!\n\n");

    for (i=0; i<nBnd; i++)
    {
      bnd0 = *(bnd[i+1]);

      for (j=bnd[i]-bnd[0]; j<bnd[i+1]-bnd[i]; j++)
      {
        if (all2free[bnd0-1] > 0 && all2free[bndNd[j]-1] > 0)

          out << " " << all2free[bnd0-1] << " " << all2free[bndNd[j]-1] << "\n";

        bnd0 = bndNd[j];
      }
    }
    out << "\n";
  }
 
  return;
}












bool Mesh::getNeighboursOnSameGeometry2D(int j0, Vector<int> &tmp)
{
  char fct[] = "Mesh::getNeighboursOnSameGeometry2D";

  double dot, n1[2], n2[2], rn1, rn2, *X = x.x;

  if (!selectNode.defm) X = x0.x;

  int i = 0, first, last, j1, j2, m, m1, m2;

  while (i < nBnd) if (j0 < bnd[i+1] - bndNd) break; else i++; 

  if (i == nBnd) prgError(1,fct,"invalid j0?!");

  first = bnd[i]   - bndNd;
  last  = bnd[i+1] - bndNd - 1;

  if (j0 > first) j1 = j0 - 1; else j1 = last;

  if (j0 < last) j2 = j0 + 1; else j2 = first;

  m1 = (bndNd[j1]-1)*ndm;
  m2 = (bndNd[j2]-1)*ndm;
  m  = (bndNd[j0]-1)*ndm;

  n1[0] = X[m+1] - X[m1+1];
  n1[1] = X[m1] - X[m];
  rn1   = n1[0]*n1[0] + n1[1]*n1[1];

  n2[0] = X[m2+1] - X[m+1];
  n2[1] = X[m] - X[m2];
  rn2   = n2[0]*n2[0] + n2[1]*n2[1];

  dot   = (n1[0]*n2[0] + n1[1]*n2[1]) / sqrt(rn1 * rn2);

  if (dot >= selectNode.cosAlpha) { tmp.append(j1); tmp.append(j2); return true; }

  if (dot < -.99999 && selectNode.alpha >= 180.) { tmp.append(j1); tmp.append(j2); return true; }

  return false;
}










bool Mesh::extrude(Domain *dom, double *pos, double *sze, int np, int ndfIn)
{
  char fct[] = "Mesh::extrude";

  if (dom->ndm != 2) { prgWarning(1,fct,"extrudes 2D meshes only!"); return false; }

  if (pos[0] > pos[np-1]) { prgWarning(2,fct,"pos[0] > pos[np-1]"); return false; }

  int nen2D   = mesh(*dom).nen,
      numnp2D = mesh(*dom).numnp,
      numel2D = mesh(*dom).numel;

  Element **elem2D = mesh(*dom).elem;

  if (nen2D != 4 && nen2D != 8 && nen2D != 9) 
    { prgWarning(3,fct,"this does not work for triangular elements!"); return false; }

  int e, i, j, k, m, n, n1, n2, nenEl, *IX, *IX2D,
      nc = 0, nf = 0, ne = 0;

  double fact, dh;

  Vector<double> h;

  VectorArray<int> cnd, end, fnd;

  ndm = 3;

  ndf = ndfIn;

  // generate list of corner nodes

  cnd.setDim(numnp2D);

  for (i=0; i<numnp2D; i++) cnd[i] = -1;

  for (e=0; e<numel2D; e++)
  {
    IX2D = elem2D[e]->ix;
    for (i=0; i<4; i++) if (cnd[IX2D[i]-1] < 0) cnd[IX2D[i]-1] = nc++;
  }

  // generate list of edge nodes

  if (nen2D > 4)
  {
    end.setDim(numnp2D);

    for (i=0; i<numnp2D; i++) end[i] = -1;

    for (e=0; e<numel2D; e++)

      if (elem2D[e]->nen() > 4) 
      {
        IX2D = elem2D[e]->ix;
        for (i=4; i<8; i++) if (end[IX2D[i]-1] < 0) end[IX2D[i]-1] = ne++;
      }
  }

  // generate list of face nodes

  if (nen2D > 8)
  {
    fnd.setDim(numnp2D);

    for (i=0; i<numnp2D; i++) fnd[i] = -1;

    for (e=0; e<numel2D; e++)

      if (elem2D[e]->nen() > 8) { IX2D = elem2D[e]->ix; fnd[IX2D[8]-1] = nf++; }
  }

  // calculate layer positions

  h.free(); h.append(pos[0]);

  n = 0;
  m = 1;

  while (h[n] < pos[np-1])
  {
    while (h[n] < pos[m])
    {
      dh = (h[n]-pos[m-1]) / (pos[m]-pos[m-1]) * (sze[m]-sze[m-1]) + sze[m-1];
      h.append(h[n] + dh);
      n++;
    }
    m++; if (m == np) break;
  }
  if (h[n]-pos[np-1] > pos[np-1]-h[n-1]) n--;
  fact = (pos[np-1]-pos[0]) / (h[n]-pos[0]); n++;
  for (i=1; i<n; i++) h[i] = h[0] + (h[i]-h[0])*fact;
  h.trunc(n);

  // allocate memory

  if      (nen2D == 4) { nen =  8; numnp = numnp2D * h.n;                  }
  else if (nen2D == 8) { nen = 20; numnp = numnp2D * h.n + nc * (h.n - 1); }
  else if (nen2D == 9) { nen = 27; numnp = numnp2D * (h.n + h.n - 1);      }

  numel = numel2D * (h.n - 1);

  x.setDim(numnp,ndm,true);

  ixTmp.setDim(numel,nen+1,true);
  
  // new corner nodes

  double *X, *X2D = mesh(*dom).x.x;

  for (j=0; j<h.n; j++)
  {
    for (i=0; i<numnp2D; i++)
    {
      if (cnd[i] > -1)
      {
         X = x.x + 3 * (cnd[i] + j*nc);
        *X = X2D[i+i];   X++;
        *X = X2D[i+i+1]; X++;
        *X = h[j];
      }
    }
  }
  m = nc * h.n;

  // new elements (corners only)

  IX = ixTmp.x;
  n = 1;
  for (j=0; j<h.n-1; j++)
  {
    for (e=0; e<numel2D; e++)
    {
      IX2D = elem2D[e]->ix;

      for (i=0; i<4; i++)     { *IX = cnd[IX2D[i]-1] + n;      IX++; }
      for (i=0; i<4; i++)     { *IX = cnd[IX2D[i]-1] + n + nc; IX++; }
      for (i=8; i<nen+1; i++) { *IX = 0;                       IX++; }
    }
    n += nc;
  }
  
  if (nen2D > 4)
  {
    // new edge nodes

    for (j=0; j<h.n; j++)
    {
      for (i=0; i<numnp2D; i++)
      {
        if (end[i] > -1)
        {
           X = x.x + 3 * (m + end[i] + j*ne);
          *X = X2D[i+i];   X++;
          *X = X2D[i+i+1]; X++;
          *X = h[j];
        }
      }
    }
    m += ne * h.n;

    for (j=0; j<h.n-1; j++)
    {
      for (i=0; i<numnp2D; i++)
      {
        if (cnd[i] > -1)
        {
           X = x.x + 3 * (m + cnd[i] + j*nc);
          *X = X2D[i+i];   X++;
          *X = X2D[i+i+1]; X++;
          *X = (h[j]+h[j+1])*.5;
        }
      }
    }

    // new elements (edge nodes)

    IX = ixTmp.x;
    n1 = nc * h.n + 1;
    n2 = ne * h.n + n1;
    for (j=0; j<h.n-1; j++)
    {
      for (e=0; e<numel2D; e++)
      {
        IX2D = elem2D[e]->ix;

        for (i=0; i<8; i++)                                         IX++;
        for (i=0; i<4; i++)     { *IX = end[IX2D[i+4]-1] + n1;      IX++; }
        for (i=0; i<4; i++)     { *IX = end[IX2D[i+4]-1] + n1 + ne; IX++; }
        for (i=0; i<4; i++)     { *IX = cnd[IX2D[i]-1]   + n2;      IX++; }
        for (i=20; i<nen+1; i++)                                    IX++;
      }
      n1 += ne;
      n2 += nc;
    }

    if (nen2D > 8)
    {
      // new face nodes


      cout << "  extrude 9-noded to 27-noded: why not implement it now?\n\n"; return false;


      // new elements (face nodes)




    }
  }

  // generate missing stuff

  completeData();

  // reset plot

  plot.reset();

  // prepare data, finalise mesh

  prepareInputData();

  cout << "   The 2D mesh has been extruded successfully.\n\n";

  return true;
}













void Mesh::completeData(void)
{
  int i, j, nst;

  nst = ndf * nen; if (ndm > ndf) nst = ndm * nen;
     
  p  = new double [nst]; 
  s  = new double [nst*nst]; 
  
  xl = new double [nst*5]; 
  ul = new double [nst*6]; 

   idx.setDim(numnp,ndm,true);
   idu.setDim(numnp,ndf,true);
    x0.setDim(numnp,ndm,true);
    xn.setDim(numnp,ndm,true);
     u.setDim(numnp,ndf,true);
    un.setDim(numnp,ndf,true);
    u3.setDim(numnp,ndf,true);
    u4.setDim(numnp,ndf,true);
    u5.setDim(numnp,ndf,true);
    u6.setDim(numnp,ndf,true);
  reac.setDim(numnp,ndf,true);

  nodeFlag.setDim(numnp);
 
  r.setDim(numnp*max(ndf,ndm)); 
  
  outp.setDim(max(ndf,ndm)*max(numel,numnp)); 
 
  idu.zero();
  idx.zero();
  u.zero(); 
  un.zero(); 
  u3.zero(); 
  u4.zero(); 
  u5.zero(); 
  u6.zero(); 

  for (i=1; i<=numnp; i++)
  {
    for (j=1; j<=ndm; j++) 
    {
       x0(i,j) = x(i,j);
       xn(i,j) = x(i,j);
    }
  }

  elemGrpPropTmp.free();
  elemGrpPropTmp.append(1);
  elemGrpPropTmp.append(1);

  numElemGrp = 1;

  elemProp.add(new PropertyItem(ELEMENTTYPE));

  elemProp[0].id = 1;

  if (ndm == 2 && nen ==  3) elemProp[0].name = "2D3nodedTriangle";
  if (ndm == 2 && nen ==  4) elemProp[0].name = "2D4nodedQuadrilateral";
  if (ndm == 2 && nen ==  8) elemProp[0].name = "2D8nodedQuadrilateral";
  if (ndm == 2 && nen ==  9) elemProp[0].name = "2D9nodedQuadrilateral";
  if (ndm == 3 && nen ==  4) elemProp[0].name = "3D4nodedTetrahedron";
  if (ndm == 3 && nen ==  8) elemProp[0].name = "3D8nodedBrick";
  if (ndm == 3 && nen == 20) elemProp[0].name = "3D20nodedBrick";
  if (ndm == 3 && nen == 27) elemProp[0].name = "3D27nodedBrick";

  return;
}








void Mesh::writeInputFile(char *fileName, int domType,
                          bool append, bool loadInterface, bool timeAndRun)
{
  char fct[] = "Mesh::writeInputFile";

  char *domKey[] = DOMAIN_KEY, tmp[200],
       *domTypeKey = domKey[domType-1];

  int c, e, i, j, nw, *IDU = idu.x, *IDX = idx.x, *IX, nenEl;

  double *X = x.x;

  Vector<char*> old;

  MyString line, *word;
 
  std::ofstream out;

  std::ifstream in;

  if (append)
  {
    in.open(fileName);

    if (in) 
    {
      while (1) { old.append(new char [5000]); in.read(old[old.n-1],4990); if (!in) break; }

      in.close();
    }
  }

  out.open(fileName);

  if (!out) 
  { 
    prgWarning(1,fct,"could not open file for writing!");
    for (i=0; i<old.n; i++) delete [] old[i];
    return;
  }

  for (i=0; i<old.n; i++) { out << old[i]; delete [] old[i]; }

  if (!append) out << "MPAP2\n\n";

  out << "\nBEGIN " << domTypeKey
      << " 1\n\n\ndimensions\n! ndf ndm\n   " << ndf << "   " << ndm 
      << "\n\n\nnodes per element\n  " << nen << "\n\n\ncoordinates\n";

  // coordinates and boundary conditions

  for (i=0; i<numnp; i++)
  {
    // node number
    out << i+1;
    // node moving or fixed
    for (j=0; j<ndm; j++) { if (*IDX > 0) out << " " << 0; else out << " " << 1; IDX++; }
    // boundary conditions
    for (j=0; j<ndf; j++) { if (*IDU > 0) out << " " << 0; else out << " " << 1; IDU++; }
    // coordinates
    for (j=0; j<ndm; j++) { out << " " << *X; X++; }
    // new line
    out << "\n";
  }
  if (loadInterface)
  {
    in.open("mpap2.interface");

    if (!in) prgWarning(1,fct,"could not open mpap2.interface");

    else
    {
      c = numnp;
      while (1)
      {
        line.getNextLine(in);
        if (line == "===") break;
        out << " " << ++c;
        for (j=0; j<ndm+ndf; j++) out << " 0"; 
        out << " " << line << "\n";
      }
    }
  }

  // elements

  out << "\n\nelements\n";

  for (e=0; e<numel; e++)
  {
    // element number
    out << e+1 << " 1";
    // element node connectivity
    nenEl = elem[e]->nen();
    IX = elem[e]->ix;
    for (i=0; i<nenEl; i++) { out << " " << *IX; IX++; }
    // new line
    out << "\n";
  }

  // prescribed displacements

  if (tmpPresDisp.n > 0)
  {
    out << "\n\nprescribed displacements\n";

    for (i=0; i<tmpPresDisp.n; i++)
    {
      out << tmpPresDisp[i].nd;
      for (j=0; j<ndf; j++) out << " 1";
      for (j=0; j<ndf; j++) out << " " << tmpPresDisp[i].u[j];
      out << "\n";
    }
  }

  // interface stuff

  if (tmpFreeNode.n > 0)
  {
    out << "\n\ninterface nodes\n! dof    freeDepth  meshDepth\n";
    for (j=0; j<ndm; j++) out << " 1";
    for (j=ndm; j<ndf; j++) out << " 0";
    out << "  : " << numnp << "  " << numnp << "  :\n! node numbers\n";
    i = 0;
    while (i < tmpFreeNode.n)
    {
      j = 0;
      while (j < 10)
      {
        out << " " << tmpFreeNode[i]; j++; i++; if (i == tmpFreeNode.n) break;
      }
      out << "\n";
    }
  }

  if (loadInterface)
  {
    if (!in) prgWarning(1,fct,"problem with mpap2.interface");

    out << "\n\ninterface nodes\n! dof    freeDepth  meshDepth\n";
    for (j=0; j<ndm; j++) out << " 1";
    for (j=ndm; j<ndf; j++) out << " 0";
    out << "  : " << numnp << "  " << numnp << "  :\n! node numbers\n";
    i = numnp+1;
    while (i < c+1)
    {
      j = 0;
      while (j < 10) { out << " " << i; j++; i++; if (i == c+1) break; }
      out << "\n";
    }
    out << "\n\ninterface interpolations\n";

    while (in)
    {
      line.getNextLine(in);
      if (!in) break;
      nw = line.split(&word);
      for (j=0; j<nw; j++) { word[j].toInt(&i); out << " " << i + numnp; }
      out << "\n";
      for (i=0; i<nw; i++) word[i].free(); delete [] word;
    } 
  }

  if (domType-1 == FLUID)
  {
    // element groups, element type, ALE type, control

    out << "\n\nelement groups\n  1 1 1\n\n\nelement type\n";

    if (ndm == 2) out << "1 2D3nodedStabIncompFluid\n";
    else          out << "1 3D4nodedStabIncompFluid\n";
  
    out << "   1. ,   0. ,   1. ,0.3333,  30. , 0.1  ,   0. ,   0. ,   0. ,   0.\n";
    out << "   1.e-2,  1.0 ,   0. ,   0. ,   0.\n\n\nALE type\n 1 fixed\n\n\ncontrol\n";
    out << "  1.e-10   1.e-10    2    0.5\n\n\nEND FLUID 1\n\n\n";
  }
  else if (domType-1 == SOLID)
  {
    // element groups, element type, material type, control

    out << "\n\nelement groups\n  1  1 1\n\n\nelement type\n";
    out << " 1  2D4nodedLinearSolid  4  1  2  1.  0. 0. 1.\n\n\nmaterial\n";
    out << " 1  NeoHookeElasticity  3000, 100, 1\n\n\ncontrol\n";
    out << "!   tol    lump_mass  tis   tis_param\n";
    out << "  1.e-8        0       2      0.9\n\n\nEND SOLID  1\n\n\n";
  }
  else prgWarning(1,fct,"no element groups and properties written for this domain type!");

  if (timeAndRun)
  {
    out << "BEGIN TIME_FUNCTIONS\n!\n! lam(t) = p1 + p2*t + p3*sin(p4*t+p5)\n";
    out << "!\n! id  t0    t1    p1   p2   p3   p4   p5\n";
    out <<    "   1  0.  10000.   1.   0.   0.   0.   0.\n\nEND TIME_FUNCTIONS\n\n\n";
    out << "BEGIN RUN_CONTROL\n\nBATCH\nsolv\ndt,,.1\nmesh\ndraw\nloop,,100\n";
    out << "  time\n  updt\n  loop,,20\n    tang\n    if,conv\n      xlop\n    ndif\n";
    out << "  next\n  prnt\n  wipe\n  mesh\n  draw\n  if,mous\n    wait,mous\n  ndif\n";
    out << "next\nend\n\nINTER\n\nEND RUN_CONTROL\n\n";
  }

  out.close();

  return;
}









void Mesh::setTmpPresDisp(int flag, double *val, double *xp, double rad)
{
  char fct[] = "Mesh::setTmpPresDisp";

  prgWarning(1,fct,"be careful: poor implementation!"); 

  Vector<int> &nodeList = selectNode.nodeList;

  int i, j, k, n = 0, *IDU = idu.x;

  double d2;

  if (flag == 1) // uniform displacement  
  {
    for (i=0; i<nodeList.n; i++)
    {
      k = nodeList[i] - 1;

      tmpPresDisp.add(new PrescribedDisplacement);

      tmpPresDisp[n].nd = k + 1;

      for (j=0; j<ndf; j++)
      {
        if (IDU[k*ndf+j] < 1) tmpPresDisp[n].u[j] = val[j];
        else                  tmpPresDisp[n].u[j] = 0.;
      }
      n++;
    }
  }
  else
  {
    for (i=0; i<nodeList.n; i++)
    {
      k = nodeList[i] - 1;

      d2 = dist2(x.x + k*ndm,xp,ndm);

      tmpPresDisp.add(new PrescribedDisplacement);

      tmpPresDisp[n].nd = k + 1;

      for (j=0; j<ndf; j++)
      {
        if (IDU[k*ndf+j] < 1) tmpPresDisp[n].u[j] = val[j] * (1. - d2/(rad*rad));
        else                  tmpPresDisp[n].u[j] = 0.;
      }
      n++;
    }
  }

  return;
}

