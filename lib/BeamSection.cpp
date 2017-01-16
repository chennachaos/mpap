
#include <iostream>

#include "BeamSection.h"
#include "FunctionsProgram.h"
#include "FunctionsSupport.h"
#include "DataBlockTemplate.h"
#include "MathGeom.h"



using namespace std;


static bool boomsOnly;






double SecondMomentOfArea::principalDirection(void)
{
  return .5 * atan(-(yz+yz)/(yy-zz));
}







void SecondMomentOfArea::rotate(double theta, SecondMomentOfArea &I)
{
  double c = cos(theta + theta),
         s = sin(theta + theta);

  yy = .5 * (I.yy+I.zz) + .5 * (I.yy-I.zz) * c - I.yz * s;

  zz = .5 * (I.yy+I.zz) - .5 * (I.yy-I.zz) * c + I.yz * s;

  yz =                  + .5 * (I.yy-I.zz) * s + I.yz * c;

  return;
}










void Node::calcCoorG(double *coorC, double &area)
{
  area += B;

  coorC[0] += B * coor[0];
  coorC[1] += B * coor[1];

  return;
}









void Node::calcI(double *coorC, SecondMomentOfArea &I)
{
  double y    = coor[0] - coorC[0],
         z    = coor[1] - coorC[1];

  I.yy += z * z * B;

  I.zz += y * y * B;

  I.yz += y * z * B;

  return;
}








void Node::draw(double mxBoom, int nBoom, int p)
{
  double d;

  int n = 50;
  
  if (nBoom > n) d = (plot.dDes[0]+plot.dDes[1]) / (double)nBoom * sqrt(B/mxBoom);
  else           d = (plot.dDes[0]+plot.dDes[1]) / (double)n     * sqrt(B/mxBoom);

  plot.fillCircle(coor,d);

  if (p < 1) return;

  double x[3];

  char tmp[10];

  sprintf(tmp,"%5d",p);

  n = 0; while (tmp[n] == ' ') n++;

  if (d < 1.e-10) d = .7 * plot.dPt();

  x[0] = coor[0] + d *.5;
  x[1] = coor[1] + d *.5;

  plot.putText(x,&(tmp[n]),1,true);
 
  return;
}








void Node::calcShearFlow(double *coorC, SecondMomentOfArea &I)
{
  double y    = coor[0] - coorC[0],
         z    = coor[1] - coorC[1],
         fact = B / (I.yy * I.zz - I.yz * I.yz);

  q[0] = fact * (I.yz * z - I.yy * y);

  q[1] = fact * (I.yz * y - I.zz * z);

  q[2] = 0.;

  return;
}









void Flange::calcGeom(List<Node> &node)
{
  double *x1 = node[n1].coor, *x2 = node[n2].coor, x0[2] = {0.,0.}, alpha;

  analysePointAndEdge2D(x1,x2,x0,&alpha,&p,&l);

  return;
}







void Flange::calcCoorG(List<Node> &node, double *coorC, double &area)
{
  double *x1 = node[n1].coor, *x2 = node[n2].coor;

  area += l * t;

  coorC[0] += l * t * .5 * (x1[0] + x2[0]);
  coorC[1] += l * t * .5 * (x1[1] + x2[1]);

  return;
}









void Flange::calcI(List<Node> &node, double *coorC, SecondMomentOfArea &I)
{
  double *x1  = node[n1].coor, 
         *x2  = node[n2].coor,
         dy   = x2[0] - x1[0],
         dz   = x2[1] - x1[1],
         area = l * t,
         fact = area / 12.,
         Iyy  = dz * dz * fact,
         Izz  = dy * dy * fact,
         Iyz  = dy * dz * fact,
         y    = .5 * (x1[0] + x2[0]) - coorC[0],
         z    = .5 * (x1[1] + x2[1]) - coorC[1];

  I.yy += Iyy + z * z * area;

  I.zz += Izz + y * y * area;

  I.yz += Iyz + y * z * area;

  return;
}









void Flange::calcShearFlow(List<Node> &node, double *coorC, SecondMomentOfArea &I)
{
  double *x1  = node[n1].coor, 
         *x2  = node[n2].coor,
         ty1  = (x1[0] - coorC[0]) * t,
         ty2  = (x2[0] - coorC[0]) * t,
         tz1  = (x1[1] - coorC[1]) * t,
         tz2  = (x2[1] - coorC[1]) * t,
         fact = l / (I.yy * I.zz - I.yz * I.yz);

  if (!boomsOnly)
  {
    q1[0] = 0.;
    qc[0] = fact * (I.yz * (3.*tz1 + tz2) - I.yy * (3.*ty1 + ty2)) * .125;
    q2[0] = fact * (I.yz * (   tz1 + tz2) - I.yy * (   ty1 + ty2)) * .5;
    Q [0] = l * (4. * qc[0] + q2[0]) / 6.;

    q1[1] = 0.;
    qc[1] = fact * (I.yz * (3.*ty1 + ty2) - I.zz * (3.*tz1 + tz2)) * .125;
    q2[1] = fact * (I.yz * (   ty1 + ty2) - I.zz * (   tz1 + tz2)) * .5;
    Q [1] = l * (4. * qc[1] + q2[1]) / 6.;
  }
  else
  {
    q1[0] = 0.;
    qc[0] = 0.;
    q2[0] = 0.;
    Q [0] = 0.;

    q1[1] = 0.;
    qc[1] = 0.;
    q2[1] = 0.;
    Q [1] = 0.;
  }
  q1[2] = 0.;
  qc[2] = 0.;
  q2[2] = 0.;
  Q [2] = 0.;

  //cout << q2[0] << "," << q2[1] << "\n";

  return;
}









void Flange::qForSandT(double *S, double T)
{
  q1[3] = q1[0] * S[0] + q1[1] * S[1] + q1[2] * T;
  qc[3] = qc[0] * S[0] + qc[1] * S[1] + qc[2] * T;
  q2[3] = q2[0] * S[0] + q2[1] * S[1] + q2[2] * T;
  Q [3] = Q [0] * S[0] + Q [1] * S[1] + Q [2] * T;

  return;
}









void Flange::tauForSandT(void)
{
  q1[4] = q1[3] / t;
  qc[4] = qc[3] / t;
  q2[4] = q2[3] / t;

  return;
}









void Flange::sigForMandN(List<Node> &node, double *M, double N, 
                         SecondMomentOfArea &I, double *coorC, double area)
{
  double *x1  = node[n1].coor, 
         *x2  = node[n2].coor,
         fact = - 1./(I.yy*I.zz - I.yz*I.yz),
         y, z; 

  y     = x1[0] - coorC[0];
  z     = x1[1] - coorC[1];
  q1[5] = fact * ((I.zz*M[0]+I.yz*M[1])*z - (I.yy*M[1]+I.yz*M[0])*y) - N / area;

  y     = .5 * (x1[0] + x2[0]) - coorC[0];
  z     = .5 * (x1[1] + x2[1]) - coorC[1];
  qc[5] = fact * ((I.zz*M[0]+I.yz*M[1])*z - (I.yy*M[1]+I.yz*M[0])*y) - N / area;

  y     = x2[0] - coorC[0];
  z     = x2[1] - coorC[1];
  q2[5] = fact * ((I.zz*M[0]+I.yz*M[1])*z - (I.yy*M[1]+I.yz*M[0])*y) - N / area;

  return;
}









void Flange::sigVonMises(void)
{
  q1[6] = - sqrt(q1[5]*q1[5] + 3.*q1[4]*q1[4]);
  qc[6] = - sqrt(qc[5]*qc[5] + 3.*qc[4]*qc[4]);
  q2[6] = - sqrt(q2[5]*q2[5] + 3.*q2[4]*q2[4]);

  return;
}









void Flange::draw(List<Node> &node, double mxThick, int p)
{
  double *x1, *x2, v[2], fact;

  x1 = node[n1].coor;
  x2 = node[n2].coor;

  plot.line(x1,x2);
  
  v[0] = x1[1] - x2[1];
  v[1] = x2[0] - x1[0];
 
  fact = 0.3 * plot.dPt() * t / (mxThick * l);

  v[0] *= fact;
  v[1] *= fact;

  plot.polyPnt[ 0] = x1[0] - v[0];
  plot.polyPnt[ 1] = x1[1] - v[1];
  plot.polyPnt[ 3] = x1[0] + v[0];
  plot.polyPnt[ 4] = x1[1] + v[1];
  plot.polyPnt[ 6] = x2[0] + v[0];
  plot.polyPnt[ 7] = x2[1] + v[1];
  plot.polyPnt[ 9] = x2[0] - v[0];
  plot.polyPnt[10] = x2[1] - v[1];
 
  plot.poly(4);

  return;
}













void Flange::plotAndPrint(List<Node> &node, int ii, char *str,
                          double mxq, double scale, double mxl, int n)
{
  int i, m = roundToInt(((double)n)*l/mxl);

  if (m < 4) m = 4;

  double *x1, *x2, x[2], q[2], q0[2], v[2],
         xi = 0, fact, dxi = 1. / (double) m, N1, Nc, N2;

  if (abs(qc[ii] - .5 * (q1[ii]+q2[ii])) < 1.e-15)
  {
    printf("(linear)    %s(xi=0.0000) = %9.5f\n",str,-q1[ii]);
  }
  else
  {
    printf("(quadratic) %s(xi=0.0000) = %9.5f\n",str,-q1[ii]);
    xi = .25 * (3.*q1[ii]+q2[ii]-4.*qc[ii]) / (q1[ii]+q2[ii]-qc[ii]-qc[ii]);
    if (xi > .5+1.e-8) 
      printf("                                    %s(xi=0.5000) = %9.5f\n",str,-qc[ii]);
    if (xi > 1.e-8 && xi < 1.-1.e-8)
    {
      N1 = 1. + 2. * xi*xi - 3. * xi;
      Nc = 4. * (xi - xi*xi);
      N2 = 2. * xi*xi - xi;
      printf("                                    %s(xi=%6.4f) = %9.5f\n",
             str,xi,-N1*q1[ii]-Nc*qc[ii]-N2*q2[ii]);
    }
    if (xi < .5-1.e-8) 
      printf("                                    %s(xi=0.5000) = %9.5f\n",str,-qc[ii]);
  }
  printf("                                    %s(xi=1.0000) = %9.5f\n",str,-q2[ii]);
  if (*str == 'q') printf("                                             Q   = %9.5f\n",-Q[ii]);

  x1 = node[n1].coor;
  x2 = node[n2].coor;
  
  plot.line(x1,x2);

  v[0] = x1[1] - x2[1];
  v[1] = x2[0] - x1[0];
 
  fact = scale / (mxq * l);

  v[0] *= fact;
  v[1] *= fact;

  q[0] = x1[0] + v[0] * q1[ii];
  q[1] = x1[1] + v[1] * q1[ii];

  plot.line(x1,q);

  xi = 0.;

  for (i=0; i<m; i++)
  {
    q0[0] = q[0];
    q0[1] = q[1];

    xi += dxi;

    N1 = 1. + 2. * xi*xi - 3. * xi;
    Nc = 4. * (xi - xi*xi);
    N2 = 2. * xi*xi - xi;

    x[0] = xi * x2[0] + (1.-xi) * x1[0];
    x[1] = xi * x2[1] + (1.-xi) * x1[1];

    fact = q1[ii] * N1 + qc[ii] * Nc + q2[ii] * N2;

    q[0] = x[0] + v[0] * fact;
    q[1] = x[1] + v[1] * fact; 

    plot.line(x,q);
    plot.line(q,q0);
  }

  return;
}











void Flange::plotLabels(List<Node> &node, int ii, double mxq, double scale)
{
  char tmp[30];

  double fact, *x1, *x2, q[2], v[2];

  int j;

  x1 = node[n1].coor;
  x2 = node[n2].coor;
  
  v[0] = x1[1] - x2[1];
  v[1] = x2[0] - x1[0];
 
  fact = scale / (mxq * l) * .75;

  v[0] *= fact;
  v[1] *= fact;

  q[0] = x1[0] + v[0] * q1[ii];
  q[1] = x1[1] + v[1] * q1[ii];

  sprintf(tmp,"%7.3f",-q1[ii]);

  j = 0; while (tmp[j] == ' ') j++;

  plot.putText(q,&(tmp[j]),5,true);

  q[0] = x2[0] + v[0] * q2[ii];
  q[1] = x2[1] + v[1] * q2[ii];

  sprintf(tmp,"%7.3f",-q2[ii]);

  j = 0; while (tmp[j] == ' ') j++;

  plot.putText(q,&(tmp[j]),5,true);

  return;
}













void BSLoop::sort(List<Node> &node, List<Flange> &flange)
{
  int i = 0, j, j0 = 1, nfirst = flange[flng[0]].n1, n = flange[flng[0]].n2;

  //flng[0]++;

  while (1)
  {
    j = j0;
    while (j < flng.n) if (flange[flng[j]].n1 != n) j++; else break;

    if (j == flng.n)
    {
      j = j0; 
      while (j < flng.n) if (flange[flng[j]].n2 != n) j++; else break;
      if (j < flng.n)
      {
        flng.move(j,j0);
        n = flange[flng[j0]].n1;
        flng[j0] = - flng[j0] - 1;
        j0++;
      }
      else prgError(1,"BSLoop::sort","fatal error!");
    }
    else 
    { 
      flng.move(j,j0);
      n = flange[flng[j0]].n2;
      //flng[j0]++;
      j0++;
    }
    if (n == nfirst) break;
  }

  double x0[2] = { 0., 0.}, *x1, *x2;

  area = 0.;

  for (i=0; i<flng.n; i++)
  {
    n = flng[i];
    if (n < 0) n = - 1 - n;
    x1 = node[flange[n].n1].coor;
    x2 = node[flange[n].n2].coor;
    if (flng[i] < 0) area -= triangleArea2D(x1,x2,x0);
    else             area += triangleArea2D(x1,x2,x0);
  }
  if (area < 0.) 
  { 
    flng.reverseOrder(); 
    area = - area; 
    for (i=0; i<flng.n; i++) flng[i] = - flng[i] - 1;
  }

  //cout << area << " = area\n";

  //cout << flng << "\n";

  return;
}












BeamSection::BeamSection(void)                       
{                                                  
  ndm   = 2;

  numnp = 0;
  numfl = 0;

  S[0]  = 0.;
  S[1]  = 0.;
  M[0]  = 0.;
  M[1]  = 0.;
  T     = 0.;
  N     = 0.;

  // add new type
  
  DomainType *beamSection = domain.newType(BEAMSECTION,ROOTDOMAIN);

  if (beamSection == NULL) return;  // domain type exists already

  beamSection->key.addNew("nodes", 
                          "flanges",
                          "internal forces");
  return;
}










BeamSection::~BeamSection(void)                     
{         
  return;
}










void BeamSection::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString tmpl, *word;
 
  char fct[] = "BeamSection::readInputData";

  int i, j;
 
  Vector<double> dTmp;
  Vector<int>    iTmp, lTmp;
  MyStringList   sTmp;

  DataBlockTemplate t1, t2;

  switch (domain[BEAMSECTION].key.whichBegins(line))
  {
    case  0: cout << "     BEAMSECTION: reading nodes ...\n\n";
	   
             if (!line.copyAfter('|',tmpl)) tmpl.free().append("123 3f");
	    
             t1.initialise(tmpl);
	     t2.initialise("123 3f");
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
		     prgError(2,fct,"data error in 'nodes'!");
	     
             numnp = dTmp.dim() / 3;
	     
	     for (i=0; i<numnp; i++) node.add(new Node(-dTmp[i+i+i],-dTmp[i+i+i+1],dTmp[i+i+i+2]));

	     break;
     	     
    case  1: cout << "     BEAMSection: reading flanges ...\n\n";
	   
             if (!line.copyAfter('|',tmpl)) tmpl.free().append("123 2i 1f");
	     
             t1.initialise(tmpl);
	     t2.initialise("123 2i 1f");
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
		     prgError(2,fct,"data error in 'flanges'!");
	     
             numfl = iTmp.dim() / 2;
	     
	     for (i=0; i<numfl; i++) flange.add(new Flange(iTmp[i+i],iTmp[i+i+1],dTmp[i]));
	     
	     break;

    case  2: cout << "     BEAMSection: reading internal forces ...\n\n";

             if (!line.copyAfter('|',tmpl)) tmpl.free().append("8f");
	     
             t1.initialise(tmpl);
	     t2.initialise("8f");
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
		     prgError(2,fct,"data error in 'internal forces'!");
	     
             if (dTmp.dim() / 8 != 1) prgError(3,fct,"data error in 'internal forces'!");;
	     
             S[0]     =   dTmp[0];
             S[1]     =   dTmp[1];
             coorA[0] = - dTmp[2];
             coorA[1] = - dTmp[3];
             T        =   dTmp[4];
             M[0]     =   dTmp[5];
             M[1]     =   dTmp[6];
             N        =   dTmp[7];

	     break;

    case -1: // go and inherit from DOMAIN
	     
	     this->Domain::readInputData(Ifile, line); 
	     
	     break;
  }
 
  return;
}









void BeamSection::prepareInputData(void)
{
  // call ancestor function

  Domain::prepareInputData();
  
  cout << "     BEAMSection: prepare input data ...\n\n";
  
  char fct[] = "BeamSection::prepareInputData"; 

  int i, j;

  // perform consistency checks
     
  for (i=0; i<numfl; i++)
    if (flange[i].n1 >= numnp || flange[i].n2 >=numnp || flange[i].n1 <0 || flange[i].n2 < 0) 
      prgError(1,fct,"invalid node numbers in 'flanges'!");

  // calculate flange lengths and distance p

  for (i=0; i<numfl; i++) flange[i].calcGeom(node);

  // get maximum thickness, maximum length and largest boom for plots
  
  double mnThick;

  mxThick  = -1;      for (i=0; i<numfl; i++) if (mxThick  < flange[i].t) mxThick  = flange[i].t;
  mnThick  = mxThick; for (i=0; i<numfl; i++) if (mnThick  > flange[i].t) mnThick  = flange[i].t;
  mxLength = -1;      for (i=0; i<numfl; i++) if (mxLength < flange[i].l) mxLength = flange[i].l;
 
  if (mnThick/mxThick > 0.7) mxThick *= 2.; 
 
  mxBoom = node[0].B;
  for (i=1; i<node.n; i++) if (mxBoom < node[i].B) mxBoom = node[i].B;

  // calculate node-flange and node-node connectivity

  for (i=0; i<numnp; i++) nodeFlange.add(new Vector<int>);

  for (i=0; i<numfl; i++)
  {  
    nodeFlange[flange[i].n1].append(i);
    nodeFlange[flange[i].n2].append(i);
  }
  for (i=0; i<numnp; i++)
    if (nodeFlange[i].n == 0) prgError(2,fct,"node without any attached flanges found!");

  for (i=0; i<numnp; i++)
  {
    nodeNode.add(new Vector<int>);
    for (j=0; j<nodeFlange[i].n; j++) 
    {
      if (flange[nodeFlange[i][j]].n1 == i) nodeNode[i].append(flange[nodeFlange[i][j]].n2);
      else                                  nodeNode[i].append(flange[nodeFlange[i][j]].n1);
    }
  }

  // calculate loop data

  calcBSLoops();

  // perform analysis

  boomsOnly = false;

  analysis();

  return;
}










void BeamSection::prepareInteractions(void)
{
  // go and inherit from ancestors

  Domain::prepareInteractions();


  cout << "   BeamSection: preparing interactions for " << domain.name(this) << " ...\n\n"; 
 
  return;
}












void BeamSection::calcBSLoops(void)
{
  int i, j, m, n, n1, n2, mxj;

  double mxAngle, angle, v1[2], v2[2], fact = 1.;

  Vector<int> tmp;

  start:

  i = 0;
  while (loop.n < numfl-numnp+1)
  {
    tmp.free();
    tmp.append(i);
    n1 = flange[i].n1;
    n2 = flange[i].n2;
    while (1)
    {
      if (nodeNode[n2].n == 1) { m = -1; break; }
      if (nodeNode[n2].n == 2) 
      { if (nodeNode[n2][0] == n1) mxj = 1; else mxj = 0; }
      else
      {
        v1[0] = node[n2].coor[0] - node[n1].coor[0];
        v1[1] = node[n2].coor[1] - node[n1].coor[1];
        mxAngle = -1.e4;
        for (j=0; j<nodeNode[n2].n; j++)
        { 
          n = nodeNode[n2][j];
          if (n != n1)
          {
            v2[0] = node[n].coor[0] - node[n2].coor[0];
            v2[1] = node[n].coor[1] - node[n2].coor[1];
            angle = fact * boundaryAngleBetweenTwoEdges2D(v1,v2,2);
            if (angle > mxAngle) { mxAngle = angle; mxj = j; }
          }
        }
      }
      j = nodeFlange[n2][mxj];
      if (tmp.contains(j,&m)) break;
      tmp.append(j);
      n1 = n2;
      n2 = nodeNode[n2][mxj]; 
    }
    if (m == 0 && isBSLoop(tmp) && isNewBSLoop(tmp))
    {
      loop.add(new BSLoop(tmp));
      for (j=0; j<tmp.n; j++) tmp[j]++;
      //cout << i+1 << ":" << m << ":" << tmp << "\n\n";
    }
    i++;
    if (i == numfl) break;
  }
  if (fact > 0.) { fact = -1.; goto start; }

  for (i=0; i<loop.n; i++) loop[i].sort(node,flange);

  return;
}













bool BeamSection::isBSLoop(Vector<int> &tmp)
{
  int i, j, n1, n2;

  Vector<int> nd, nn;

  for (i=0; i<tmp.n; i++)
  {
    n1 = flange[tmp[i]].n1;
    n2 = flange[tmp[i]].n2;

    j = 0; while (j < nd.n && nd[j] != n1) j++;
    if (j == nd.n) { nd.append(n1); nn.append(1); }
    else nn[j]++;

    j = 0; while (j < nd.n && nd[j] != n2) j++;
    if (j == nd.n) { nd.append(n2); nn.append(1); }
    else nn[j]++;
  }

  for (i=0; i<nn.n; i++) if (nn[i] != 2) return false;

  return true;
}









bool BeamSection::isNewBSLoop(Vector<int> &tmp)
{
  int i = 0, j;

  for (i=0; i<loop.n; i++)
  {
    j = 0;
    while (j < tmp.n && loop[i].flng.contains(tmp[j])) j++;
    if (j == tmp.n) return false;
  }
  return true;
}










void BeamSection::analysis(void)
{
  int i, j, n;

  double fact, torque;

  // calculate location of centre of gravity

  coorC[0] = 0.;
  coorC[1] = 0.;

  area = 0.;

  if (!boomsOnly) for (i=0; i<numfl; i++) flange[i].calcCoorG(node,coorC,area);

  for (i=0; i<numnp; i++) node[i].calcCoorG(coorC,area);

  coorC[0] /= area;
  coorC[1] /= area;

  // calculate second moments of area and principal direction

  I.zero();

  if (!boomsOnly) for (i=0; i<numfl; i++) flange[i].calcI(node,coorC,I);

  for (i=0; i<numnp; i++) node[i].calcI(coorC,I);

  alpha = I.principalDirection();

  pI.rotate(alpha,I);

  // calculate basic shear flows

  for (i=0; i<numfl; i++) flange[i].calcShearFlow(node,coorC,I);

  for (i=0; i<numnp; i++)   node[i].calcShearFlow(coorC,I);

  calcShearFlow(0);
  calcShearFlow(1);
  if (loop.n > 0) calcShearFlow(2);

  // calculate location of shear centre

  coorS[0] = 0.;
  coorS[1] = 0.;

  for (i=0; i<numfl; i++)
  {
    coorS[0] += flange[i].Q[1] * flange[i].p;
    coorS[1] -= flange[i].Q[0] * flange[i].p;
  } 

  // update T by adding torque from external shear force

  torque = T + (coorS[0] - coorA[0]) * S[1] - (coorS[1] - coorA[1]) * S[0];

  // calculate torsion constant
 
  J = 0.;
  if (loop.n > 0)
  {
    for (j=0; j<loop[0].flng.n; j++)
    {
      n = loop[0].flng[j];
      if (n < 0) { n = - 1 - n; fact = -1.; } else fact = 1.;      
      J += fact * flange[n].Q[2] / flange[n].t;
    }
    J = - 2. * loop[0].area / J;
  }
  else for (i=0; i<numfl; i++)  J += flange[i].t * flange[i].t * flange[i].t * flange[i].l / 3.;

  // calculate shear flow due to S in A and T

  for (i=0; i<numfl; i++) flange[i].qForSandT(S,torque);
  
  // calculate shear stress tau due to S in A and T

  for (i=0; i<numfl; i++) flange[i].tauForSandT();

  // calculate axial stress due to M and N

  for (i=0; i<numfl; i++) flange[i].sigForMandN(node,M,N,I,coorC,area);

  // calculate von Mises stress due to N, S in A, M, T

  for (i=0; i<numfl; i++) flange[i].sigVonMises();

  return;
}











void BeamSection::calcShearFlow(int ii)
{
  if (ii < 0 || ii > 2) prgError(1,"BeamSection::calcShearFlow","invalid value of ii!");

  int i, j, n, *P, r, isw = 1;

  double *K, *R, *q0, fact;

  K  = new double [numfl*numfl];
  R  = new double [numfl];
  q0 = new double [numfl];

  P  = new int [numfl];

  for (i=0; i<numfl; i++)
  {
    for (j=0; j<numfl; j++) K[i*numfl+j] = 0.;
    R [i] = 0.;
    q0[i] = 0.;
  } 

  for (i=0; i<numnp-1; i++)
  {
    for (j=0; j<numfl; j++)
    {
      if (flange[j].n1 == i)
      {
        K[i+j*numfl] += -1.;
      }
      else if (flange[j].n2 == i)
      {
        R[i]         -= flange[j].q2[ii];
        K[i+j*numfl] += +1.;
      }
    }
    R[i] -= node[i].q[ii];
  }
  if (ii < 2)
  {
    for (i=0; i<loop.n; i++)
    {
      r = numnp - 1 + i;
      for (j=0; j<loop[i].flng.n; j++)
      {
        n = loop[i].flng[j];
        if (n < 0) { n = - 1 - n; fact = -1.; } else fact = 1.;
        R[r]         -= fact * flange[n].Q[ii] / flange[n].t;
        K[r+n*numfl] += fact * flange[n].l / flange[n].t;
      }
    }
  }
  else
  {
    for (i=1; i<loop.n; i++)
    {
      r = numnp - 2 + i;
      for (j=0; j<loop[0].flng.n; j++)
      {
        n = loop[0].flng[j];
        if (n < 0) { n = - 1 - n; fact = -1.; } else fact = 1.;
        K[r+n*numfl] += fact * flange[n].l / (flange[n].t * loop[0].area);

        //cout << 0 << "," << loop[0].flng[j] << ":" << loop[0].area << "\n";
      }
      for (j=0; j<loop[i].flng.n; j++)
      {
        n = loop[i].flng[j];
        if (n < 0) { n = - 1 - n; fact = -1.; } else fact = 1.;
        K[r+n*numfl] -= fact * flange[n].l / (flange[n].t * loop[i].area);

        //cout << i << "," << loop[i].flng[j] << ":" << loop[i].area << "\n";
      }
    }
    R[numfl-1] = -1.;
    for (i=0; i<numfl; i++) K[numfl-1+i*numfl] = flange[i].l * flange[i].p;
  }

  //prgPrintSimpleMatrix(K, numfl, numfl, 12, 5, true, 5, false);

  decomplr_matrix_(K,P,&numfl,&isw);

  solve_matrix_(K,P,R,q0,&numfl);

  for (i=0; i<numfl; i++)
  {
    flange[i].q1[ii] += q0[i];
    flange[i].qc[ii] += q0[i];
    flange[i].q2[ii] += q0[i];
    flange[i].Q [ii]  = flange[i].l * (flange[i].q1[ii]
                                + 4. * flange[i].qc[ii]
                                     + flange[i].q2[ii]) / 6.;
  }

  delete [] K, R, q0, P;

  //cout << "\n";
  //for (i=0; i<numfl; i++)
  //  cout << flange[i].q1[ii] << " -> " << flange[i].qc[ii] << " -> " << flange[i].q2[ii] << "\n";
  //cout << "\n";

  return;
}











void BeamSection::plotAndPrint(int ii, char *str, double scale, int n, bool labels)
{
  int i;

  double mxq = 0.;

  for (i=0; i<numfl; i++)
  {
    if (abs(flange[i].q1[ii]) > mxq) mxq = abs(flange[i].q1[ii]);
    if (abs(flange[i].qc[ii]) > mxq) mxq = abs(flange[i].qc[ii]);
    if (abs(flange[i].q2[ii]) > mxq) mxq = abs(flange[i].q2[ii]);
  }

  for (i=0; i<numfl; i++) 
  {
    COUT; printf("   flange %2d: ",i+1);
    flange[i].plotAndPrint(node,ii,str,mxq,scale,mxLength,n);
  }
  cout << "\n";

  if (labels) for (i=0; i<numfl; i++) flange[i].plotLabels(node,ii,mxq,scale);

  return;
}










void BeamSection::doForSection(bool section, bool CA, bool CS, bool bOnly, bool geom, bool SinA, 
                               bool MinCA, bool labels, int distr, double scl, int n)
{
  double a, b, fact, cc, xmx[3], xmn[3], d, ss, x1[2], x2[2], coor[2];

  int i, col = plot.currStdColour;

  if (bOnly && !boomsOnly)      { boomsOnly = true;  analysis(); }

  else if (!bOnly && boomsOnly) { boomsOnly = false; analysis(); }

  if (section)
  {
    for (i=0; i<numfl; i++)
      if (labels) flange[i].draw(node, mxThick,i+1); else flange[i].draw(node, mxThick,0);

    for (i=0; i<numnp; i++) 
      if (labels) node[i].draw(mxBoom,node.n,i+1); else node[i].draw(mxBoom,node.n,0);
  }

  if (CA)
  {
    plot.point(coorC);
    if (labels)
    {
      coor[0] = coorC[0] + .4*plot.dPt();
      coor[1] = coorC[1] + .4*plot.dPt();
      plot.putText(coor,"CA",1,true);
    }

    if (!boomsOnly) COUT << "location of the centroid CA:\n\n";
    else            COUT << "location of the centroid CA (booms only):\n\n";
    COUT; printf("   y = %f\n",  coorC[0]);
    COUT; printf("   z = %f\n\n",coorC[1]);
  }

  if (CS)
  {
    plot.point(coorS);
    if (labels)
    {
      coor[0] = coorS[0] + .4*plot.dPt();
      coor[1] = coorS[1] + .4*plot.dPt();
      plot.putText(coor,"CS",1,true);
    }
    if (!boomsOnly) COUT << "location of the shear centre CS:\n\n";
    else            COUT << "location of the shear centre CS (booms only):\n\n";
    COUT; printf("   y = %f\n",  -coorS[0]);
    COUT; printf("   z = %f\n\n",-coorS[1]);
  }

  if (geom)
  {
    d = 0.7 * diameter;

    x1[0] = coorC[0] + cos(alpha) * d;
    x1[1] = coorC[1] + sin(alpha) * d;
    x2[0] = coorC[0] - cos(alpha) * d;
    x2[1] = coorC[1] - sin(alpha) * d;
    plot.line(x1,x2);    

    x1[0] = coorC[0] + sin(alpha) * d;
    x1[1] = coorC[1] - cos(alpha) * d;
    x2[0] = coorC[0] - sin(alpha) * d;
    x2[1] = coorC[1] + cos(alpha) * d;
    plot.line(x1,x2);    

    if (!boomsOnly) COUT << "second moments of area:\n\n";
    else            COUT << "second moments of area (booms only):\n\n";
    COUT; printf("    Iyy  = %f\n",  I.yy);
    COUT; printf("    Izz  = %f\n",  I.zz);
    COUT; printf("    Iyz  = %f\n\n",I.yz);

    if (!boomsOnly) COUT << "with respect to principal axes:\n\n";
    else            COUT << "with respect to principal axes (booms only):\n\n";
    COUT; printf("    I11  = %f\n",  pI.yy);
    COUT; printf("    I22  = %f\n",  pI.zz);
    COUT; printf("    I12  = %f\n\n",0.   );

    if (!boomsOnly) COUT << "orientation of strong axis:\n\n";
    else            COUT << "orientation of strong axis (booms only):\n\n";
    COUT; printf("   alpha = %f\n\n",180.*alpha/pi);

    if (!boomsOnly) COUT << "torsion constant:\n\n";
    else            COUT << "torsion constant (booms only):\n\n";
    COUT; printf("      J  = %f\n\n",J);
  }

  if (SinA)
  {
    d = sqrt((plot.dDes[0]*plot.dDes[0]+plot.dDes[1]*plot.dDes[1])/(S[0]*S[0]+S[1]*S[1])) * .2;
    coor[0] = - S[0] * d;
    coor[1] = - S[1] * d;
    d = 1.;
    plot.arrow(coorA,coor,d);
  }

  if (MinCA)
  {
    d = sqrt((plot.dDes[0]*plot.dDes[0]+plot.dDes[1]*plot.dDes[1])/(M[0]*M[0]+M[1]*M[1])) * .2;
    coor[0] = - M[0] * d;
    coor[1] = - M[1] * d;
    x1[0]   = coorC[0] + .5 * coor[0];
    x1[1]   = coorC[1] + .5 * coor[1];
    d = 1.;
    plot.arrow(x1,coor,d);
  }

  d = diameter * 0.15 * scl;

  switch (distr)
  {
    case 1: COUT << "shear flow q due to some Sy through CS  [Sy/length]:\n\n"; 

            if (boomsOnly) COUT << "  (booms only)\n\n";
 
            plotAndPrint(0,"q",d,n,labels); break;

    case 2: COUT << "shear flow q due to some Sz through CS  [Sz/length)]:\n\n"; 
 
            if (boomsOnly) COUT << "  (booms only)\n\n";
 
            plotAndPrint(1,"q",d,n,labels); break;

    case 3: if (loop.n < 1) { COUT << "q(T) == 0 for open sections!\n\n"; break; }

            COUT << "shear flow q due to torque  [T/(length^2)]:\n\n"; 
 
            if (boomsOnly) COUT << "  (booms only)\n\n";
 
            plotAndPrint(2,"q",d,n,labels); break;

    case 4: COUT << "shear flow q due to S through A and T [force/length]:\n\n";

            if (boomsOnly) COUT << "  (booms only)\n\n";
 
            if (S[0]*S[0]+S[1]*S[1]+T*T < 1.e-10) 
              { COUT << "  Define choose non-zero internal forces in the input file!\n\n"; break; }

            plotAndPrint(3,"q",d,n,labels); break;

    case 5: COUT << "shear stress tau due to S through A and T [force/(length^2)]:\n\n";

            if (boomsOnly) COUT << "  (booms only)\n\n";
 
            if (S[0]*S[0]+S[1]*S[1]+T*T < 1.e-10) 
              { COUT << "  Define non-zero internal forces in the input file!\n\n"; break; }

            plotAndPrint(4,"tau",d,n,labels); break;

    case 6: COUT << "axial stress sig due to M and N [force/(length^2)]:\n\n";

            if (boomsOnly)
              { COUT << "  The selected option 'booms only' does not make sense!\n\n"; break; }
 
            if (M[0]*M[0]+M[1]*M[1]+N*N < 1.e-10) 
              { COUT << "  Define non-zero internal forces in the input file!\n\n"; break; }

            plotAndPrint(5,"sig",d,n,labels); 

            d = 0.7 * diameter;
        
            a    = I.yy*M[1] + I.yz*M[0];
            b    = I.zz*M[0] + I.yz*M[1];
            fact = atan(a/b); 
            cc   = cos(fact);
            ss   = sin(fact);

            COUT; printf("   orientation of neutral axis:  alpha_n = %f\n\n",180.*fact/pi);

            fact = - N * (I.yy * I.zz - I.yz * I.yz) / (area * (a * ss + b * cc));

            x1[0] = coorC[0] - fact * ss + cc * d;
            x1[1] = coorC[1] + fact * cc + ss * d;
            x2[0] = coorC[0] - fact * ss - cc * d;
            x2[1] = coorC[1] + fact * cc - ss * d;
            plot.line(x1,x2);    
            
            break;

    case 7: COUT "von Mises stress sig* due to N, S, M and S [force/(length^2)]:\n\n";

            if (boomsOnly)
              { COUT << "  The selected option 'booms only' does not make sense!\n\n"; break; }

            if (S[0]*S[0]+S[1]*S[1]+T*T+M[0]*M[0]+M[1]*M[1]+N*N < 1.e-10) 
              { COUT << "  Define non-zero internal forces in the input file!\n\n"; break; }

            plotAndPrint(6,"sig*",d,n,labels); break;
  }

  return;
}










void BeamSection::findMinMaxX(double *xmn, double *xmx, bool defFlg)
{
  int i, j;

  if (numnp < 1) return;
  
  xmn[0] = node[0].coor[0];
  xmn[1] = node[0].coor[1];
  xmx[0] = node[0].coor[0];
  xmx[1] = node[0].coor[1];
  
  for (i=1; i<numnp; i++) 
    for (j=0; j<2; j++) 
    {
      if (xmn[j] > node[i].coor[j]) xmn[j] = node[i].coor[j];
      if (xmx[j] < node[i].coor[j]) xmx[j] = node[i].coor[j];
    }

  for (j=0; j<2; j++) 
  {
    if (xmn[j] > coorC[j]) xmn[j] = coorC[j];
    if (xmx[j] < coorC[j]) xmx[j] = coorC[j];
  }

  diameter = sqrt((xmn[0]-xmx[0])*(xmn[0]-xmx[0])+(xmn[1]-xmx[1])*(xmn[1]-xmx[1]));

  return;
}











