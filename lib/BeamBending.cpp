
#include <iostream>

#include "BeamBending.h"
#include "FunctionsProgram.h"
#include "FunctionsSupport.h"
#include "DataBlockTemplate.h"




using namespace std;






BeamBending::BeamBending(void)                       
{                                                  
  Iyy = -1.;

  ndm = 2;
 
  for (int i=0; i<10; i++) h0[i] = -1.;

  // add new type
  
  DomainType *beamBending = domain.newType(BEAMBENDING,ROOTDOMAIN);

  if (beamBending == NULL) return;  // domain type exists already

  beamBending->key.addNew("second moments of area","bearings and loads");

  return;
}










BeamBending::~BeamBending(void)                     
{         
  return;
}










void BeamBending::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString tmpl, *word;
 
  char fct[] = "BeamBending::readInputData",
       *SHLType[] = {"w","v","bz","by","My","Mz","Sz","Sy",
                     "wj","vj","bzj","byj","Ly","Lz","Fz","Fy",
                     "pz0","pz1","py0","py1",NULL};
  int i, j, nw;
 
  Vector<double>   dTmp;
  Vector<int>      iTmp, tmp2, lTmp;
  VectorArray<int> perm;
  MyStringList     sTmp;

  DataBlockTemplate t1, t2;

  switch (domain[BEAMBENDING].key.whichBegins(line))
  {
    case  0: cout << "     BEAMBENDING: reading second moments of area ...\n\n";

             if (Iyy > 0.) prgError(1,fct,"'second moments of area' have already been read.");

	     line.getNextLine(Ifile);
	     
	     nw = line.split(&word);
            
	     if (nw < 3) prgError(1,fct,"input error in 'second moments of area'!");
		     
             if (!word[0].toDbl(&Iyy)) prgError(2,fct,"input error in 'second moments of area'!");
             if (!word[1].toDbl(&Izz)) prgError(3,fct,"input error in 'second moments of area'!");
             if (!word[2].toDbl(&Iyz)) prgError(4,fct,"input error in 'second moments of area'!");
             if (Iyy < 0. || Izz < 0.) prgError(5,fct,"input error in 'second moments of area'!");
	     
             for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     
       	     line.getNextLine(Ifile);

             break;

    case  1: cout << "     BEAMBENDING: reading bearings and loads ...\n\n";

             if (x.n > 0) prgError(1,fct,"data error in 'bearings and loads'!");

             if (!line.copyAfter('|',tmpl)) tmpl.free().append("1f 1s 1f");

             t1.initialise(tmpl);
	     t2.initialise("1f 1s 1f");
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
	       prgError(2,fct,"data error in 'bearings and loads'!");

             j = sTmp.n;
             for (i=0; i<j; i++)
             {
               tmp2.append(i);
               dTmp.move(i+i,i);
             }

             if (perm.n < j) perm.setDim(j);
             for (i=0; i<j; i++) perm.x[i] = tmp2[i];
             dTmp.quickSort(false,perm.x,0,j-1);
             for (i=0; i<j; i++) tmp2[i] = perm.x[i];

             for (i=0; i<j; i++)
             {
                 x.append(dTmp[i]);
               typ.append(sTmp[tmp2[i]].which(SHLType));
               val.append(dTmp[j+tmp2[i]]);
               if (typ[i] < 0) prgError(3,fct,"invalid type identifier!");
             }

	     break;

    case -1: // go and inherit from DOMAIN
	     
	     this->Domain::readInputData(Ifile, line); 
	     
	     break;
  }
 
  return;
}











void BeamBending::prepareInputData(void)
{
  // call ancestor function

  Domain::prepareInputData();
  
  cout << "     BEAMBENDING: prepare input data ...\n\n";
  
  char fct[] = "BeamBending::prepareInputData"; 

  int i, j, n;

  double y, xpos;

  Vector<int> tmp;

  if (Iyy < 0.) prgError(1,fct,"second moments of area undefined!");

  // find minimum section length

  if (x.n < 1) return;

  minl = x[x.n-1] - x[0];

  for (i=1; i<x.n; i++) if (x[i]-x[i-1] > 1.e-12 && x[i]-x[i-1] < minl) minl = x[i]-x[i-1];

  //cout << minl << " = minl\n\n";

  // find number of sections, number of integration constants and find section lengths

  nSect = 0;
  for (i=1; i<x.n; i++) if (x[i]-x[i-1] > 1.e-12) { length.append(x[i]-x[i-1]); nSect++; }

  nIntC = nSect * 8;

  // check for double entry errors

  for (i=0; i<x.n; i++)
  {
    tmp.append(typ[i]);

    if (i == x.n - 1 || x[i+1]-x[i] > 1.e-12)
    {
      if (!tmp.distinct()) prgError(1,fct,"bending: double data");
      tmp.free();
    }
  }

  // get distributed loads

  for (i=0; i<nSect+nSect; i++) { py.append(0.); pz.append(0.); }

  n = 0;
  for (i=0; i<x.n; i++)
  {
    if (i > 0 && x[i]-x[i-1] > 1.e-12) n++;

    switch (typ[i])
    {
     case 16: if (n < nSect) pz[n+n]   = val[i]; else prgError(1,fct,"bending: pz error"); break;
     case 17: if (n > 0)     pz[n+n-1] = val[i]; else prgError(2,fct,"bending: pz error"); break;
     case 18: if (n < nSect) py[n+n]   = val[i]; else prgError(3,fct,"bending: py error"); break;
     case 19: if (n > 0)     py[n+n-1] = val[i]; else prgError(4,fct,"bending: py error"); break;
    }
    if (typ[i] > 15) typ[i] = -1;
  }

  // perform checks and insert missing boundary conditions
  // delete pz and py entries

  checkAndCompleteBC();

  // calculate the integration constants

  analysis();

  // determine maxima and minima of bending data  w v bz by My Mz Sz Sy

  for (i=0; i<8; i++) { mn[i] = 0.; mx[i] = 0.; }

  for (i=0; i<nSect; i++)
  {
    xpos = 0.;
    for (j=0; j<21; j++)
    {
      y = w(xpos,i);  if (y < mn[0]) mn[0] = y; if (y > mx[0]) mx[0] = y;
      y = v(xpos,i);  if (y < mn[1]) mn[1] = y; if (y > mx[1]) mx[1] = y;
      y = bz(xpos,i); if (y < mn[2]) mn[2] = y; if (y > mx[2]) mx[2] = y;
      y = by(xpos,i); if (y < mn[3]) mn[3] = y; if (y > mx[3]) mx[3] = y;
      y = My(xpos,i); if (y < mn[4]) mn[4] = y; if (y > mx[4]) mx[4] = y;
      y = Mz(xpos,i); if (y < mn[5]) mn[5] = y; if (y > mx[5]) mx[5] = y;
      y = Sz(xpos,i); if (y < mn[6]) mn[6] = y; if (y > mx[6]) mx[6] = y;
      y = Sy(xpos,i); if (y < mn[7]) mn[7] = y; if (y > mx[7]) mx[7] = y;

      xpos += 0.05 * length[i];
    }
  }
  return;
}










void BeamBending::checkAndCompleteBC(void)
{
  char fct[] = "BeamBending::checkAndComplete";

  Vector<int> tmp, tmp2;

  int i, j, n = 0, p;

  for (i=0; i<x.n; i++)
  {
    if (typ[i] > -1) tmp.append(typ[i]);

    if (i == x.n - 1 || x[i+1]-x[i] > 1.e-12)
    {
      if (n == 0 || n == nSect)
      {
        tmp2.free(); for (j=0; j<4; j++) tmp2.append(j + j);
        if (!tmp.containsSomeOf(2,tmp2))
        {
          if (n == 0) prgError(1,fct,"XZ-plane: invalid left bearing!");
          else        prgError(1,fct,"XZ-plane: invalid right bearing!");
        }
        for (j=0; j<4; j++) tmp2[j] = j + j + 1;
        if (!tmp.containsSomeOf(2,tmp2))
        {
          if (n == 0) prgError(1,fct,"XY-plane: invalid left bearing!");
          else        prgError(1,fct,"XY-plane: invalid right bearing!");
        }
        j = 0; while (j < tmp.n && tmp[j] < 8) j++;
        if (j < tmp.n) prgError(2,fct,"use w, v, bz, by, My, Mz, Sz, Sy only on beam ends!" );
      }
      else
      {
        tmp2.free(); for (j=0; j<4; j++) tmp2.append(j + j);
        if (tmp.containsSomeOf(2,tmp2)) 
          prgError(1,fct,"XZ-plane: overconstrained intersection support!");
        if (tmp.containsAtLeastOneOf(tmp2,&p))
        {
          if (tmp[p] == 0 || tmp[p] == 6)
          {
            if (!tmp.contains(10)) { x.insert(x[i],i); typ.insert(10,i); val.insert(0.,i++); }
            if (!tmp.contains(12)) { x.insert(x[i],i); typ.insert(12,i); val.insert(0.,i++); }
          }
          else
          {
            if (!tmp.contains( 8)) { x.insert(x[i],i); typ.insert( 8,i); val.insert(0.,i++); }
            if (!tmp.contains(14)) { x.insert(x[i],i); typ.insert(14,i); val.insert(0.,i++); }
          }
        }
        else
        {
          for (j=0; j<4; j++)
            if (!tmp.contains(8+j+j))
              { x.insert(x[i],i); typ.insert(8+j+j,i); val.insert(0.,i++); }
        }

        tmp2.free(); for (j=0; j<4; j++) tmp2.append(1 + j + j);
        if (tmp.containsSomeOf(2,tmp2)) 
          prgError(1,fct,"XY-plane: overconstrained intersection support!");
        if (tmp.containsAtLeastOneOf(tmp2,&p))
        {
          if (tmp[p] == 1 || tmp[p] == 7)
          {
            if (!tmp.contains(11)) { x.insert(x[i],i); typ.insert(11,i); val.insert(0.,i++); }
            if (!tmp.contains(13)) { x.insert(x[i],i); typ.insert(13,i); val.insert(0.,i++); }
          }
          else
          {
            if (!tmp.contains( 9)) { x.insert(x[i],i); typ.insert( 9,i); val.insert(0.,i++); }
            if (!tmp.contains(15)) { x.insert(x[i],i); typ.insert(15,i); val.insert(0.,i++); }
          }
        }
        else
        {
          for (j=0; j<4; j++)
            if (!tmp.contains(9+j+j))
              { x.insert(x[i],i); typ.insert(9+j+j,i); val.insert(0.,i++); }
        }
      }
      tmp.free();
      n++;
    }
  }

  // delete pz and py entries

  i = 0; while (i < x.n) { if (typ[i] < 0) { x.del(i); typ.del(i); val.del(i); } else i++; }

  //cout << " x:   " << x << "\n typ: " << typ << "\n val: " << val << "\n\n";

  if (x[4]    -x[3]     < .5*minl) prgError(11,fct,"error in data!");
  if (x[x.n-4]-x[x.n-5] < .5*minl) prgError(12,fct,"error in data!");

  return;
}











void BeamBending::analysis(void)
{
  // calculate the integration constants for bending

  char fct[] = "BeamBending::analysis";

  int i, j, n = 0, r = 1;

  bmtx.setDim(nIntC,nIntC);
  brhs.setDim(nIntC);

  for (i=0; i<nIntC; i++) { brhs[i] = 0.; for (j=0; j<nSect*8; j++) bmtx(i+1,j+1) = 0.; }

  // left bearing

  for (i=0; i<4; i++)
  {
    switch (typ[i])
    {
      case 0:  addW0(r++,n,val[i],1.); break;
      case 1:  addV0(r++,n,val[i],1.); break;
      case 2: addBz0(r++,n,val[i],1.); break;
      case 3: addBy0(r++,n,val[i],1.); break;
      case 4: addMy0(r++,n,val[i],1.); break;
      case 5: addMz0(r++,n,val[i],1.); break;
      case 6: addSz0(r++,n,val[i],1.); break;
      case 7: addSy0(r++,n,val[i],1.); break;
    }
  }

  // supports, loads and hinges

  for (i=4; i<typ.n-4; i++)
  {
    switch (typ[i])
    {
      case  0:  addW1(r++,n,val[i],1.);  addW0(r++,n+1,val[i],1.); break;
      case  1:  addV1(r++,n,val[i],1.);  addV0(r++,n+1,val[i],1.); break;
      case  2: addBz1(r++,n,val[i],1.); addBz0(r++,n+1,val[i],1.); break;
      case  3: addBy1(r++,n,val[i],1.); addBy0(r++,n+1,val[i],1.); break;
      case  4: addMy1(r++,n,val[i],1.); addMy0(r++,n+1,val[i],1.); break;
      case  5: addMz1(r++,n,val[i],1.); addMz0(r++,n+1,val[i],1.); break;
      case  6: addSz1(r++,n,val[i],1.); addSz0(r++,n+1,val[i],1.); break;
      case  7: addSy1(r++,n,val[i],1.); addSy0(r++,n+1,val[i],1.); break;
      case  8:  addW1(r,n,val[i],1.);  addW0(r++,n+1,0.,-1.); break;
      case  9:  addV1(r,n,val[i],1.);  addV0(r++,n+1,0.,-1.); break;
      case 10: addBz1(r,n,val[i],1.); addBz0(r++,n+1,0.,-1.); break;
      case 11: addBy1(r,n,val[i],1.); addBy0(r++,n+1,0.,-1.); break;
      case 12: addMy1(r,n,val[i],1.); addMy0(r++,n+1,0.,-1.); break;
      case 13: addMz1(r,n,val[i],1.); addMz0(r++,n+1,0.,-1.); break;
      case 14: addSz1(r,n,val[i],1.); addSz0(r++,n+1,0.,-1.); break;
      case 15: addSy1(r,n,val[i],1.); addSy0(r++,n+1,0.,-1.); break;
      default: prgError(2,fct,"Cannot handle this b.c. yet OR b.c. is rubbish!");
    }
    if (x[i+1]-x[i] > 1.e-12) n++; 
  }

  // right bearing

  for (i=typ.n-4; i<typ.n; i++)
  {
    switch (typ[i])
    {
      case 0:  addW1(r++,n,val[i],1.); break;
      case 1:  addV1(r++,n,val[i],1.); break;
      case 2: addBz1(r++,n,val[i],1.); break;
      case 3: addBy1(r++,n,val[i],1.); break;
      case 4: addMy1(r++,n,val[i],1.); break;
      case 5: addMz1(r++,n,val[i],1.); break;
      case 6: addSz1(r++,n,val[i],1.); break;
      case 7: addSy1(r++,n,val[i],1.); break;
    }
  }

  //cout << bmtx << "\n\n";
  //cout << brhs << "\n\n";

  i = 1;
  j = nSect * 8;

  VectorArray<int> P;

  A.setDim(j);
  P.setDim(j);

  decomplr_matrix_(bmtx.x,P.x,&j,&i);
  solve_matrix_(bmtx.x,P.x,brhs.x,A.x,&j);

  //cout << A << "\n\n";

  return;
}












void BeamBending::prepareInteractions(void)
{
  // go and inherit from ancestors

  Domain::prepareInteractions();

  cout << "   BeamBending: preparing interactions for " << domain.name(this) << " ...\n\n"; 
 
  return;
}










void BeamBending::findMinMaxXBB(bool syst, bool defl, bool rota, bool bmom, bool shfo, 
                                bool XZ, bool XY, double *xmn, double *xmx)
{
  refHdW = .5 / (x[x.n-1]-x[0]) * minl;

  if (refHdW > .40) refHdW = .40;
  if (refHdW < .15) refHdW = .15;

  double dd = 0.1 * refHdW;

  dh[0] = .6 * refHdW; dh[1] = dh[0];
  dh[2] = .6 * refHdW; dh[3] = dh[2];
  dh[4] = 1. * refHdW; dh[5] = dh[4];
  dh[6] = .9 * refHdW; dh[7] = dh[6];
  dh[8] = .9 * refHdW; dh[9] = dh[8];

  double h = dd;

  for (int i=0; i<10; i++) h0[i] = -1.;

  if (shfo) { if (XY) { h0[7] = h; h += dh[7]+dd; } if (XZ) { h0[6] = h; h += dh[6]+dd; } }
  if (bmom) { if (XY) { h0[5] = h; h += dh[5]+dd; } if (XZ) { h0[4] = h; h += dh[4]+dd; } }
  if (rota) { if (XY) { h0[3] = h; h += dh[3]+dd; } if (XZ) { h0[2] = h; h += dh[2]+dd; } }
  if (defl) { if (XY) { h0[1] = h+2.*dd; h += dh[1]+5.*dd; } 
              if (XZ) { h0[0] = h+2.*dd; h += dh[0]+5.*dd; } }
  if (syst) { if (XY) { h0[9] = h; h += dh[9]+dd; } if (XZ) { h0[8] = h; h += dh[8]+dd; } }

  xmn[0] = -.3;
  xmx[0] = 1.1;

  xmn[1] = 0.;
  xmx[1] = h;

  //cout << xmn[1] << " -> " << xmx[1] << "\n\n";
 
  return;
}








void BeamBending::doForBending(bool syst, bool defl, bool rota, bool bmom, bool shfo, 
                               bool XZ, bool XY, bool labels, int n)
{
  double xmx[3], xmn[3],
         mxp = 0., mxF = 0., mxv, mnv;

  int i, c = plot.currStdColour;
 
  Vector<int> tmp;

  if ((XZ && syst && h0[8] < 0.) || (XY && syst && h0[9] < 0.)
   || (XZ && defl && h0[0] < 0.) || (XY && defl && h0[1] < 0.)
   || (XZ && rota && h0[2] < 0.) || (XY && rota && h0[3] < 0.)
   || (XZ && bmom && h0[4] < 0.) || (XY && bmom && h0[5] < 0.)
   || (XZ && shfo && h0[6] < 0.) || (XY && shfo && h0[7] < 0.))
  {
    findMinMaxXBB(syst,defl,rota,bmom,shfo,XZ,XY,xmn,xmx);
    plot.fit(xmn,xmx,ndm,0.);
    plot.setColour(8); plot.wipe(); plot.setColour(c);
  }

  if (syst)
  {
    tmp.free(); tmp.append(6); tmp.append(14);
    if (XZ) for (i=0;i<typ.n;i++) if (tmp.contains(typ[i])) if (abs(val[i])>mxF) mxF=abs(val[i]);
    tmp.free(); tmp.append(7); tmp.append(15);
    if (XY) for (i=0;i<typ.n;i++) if (tmp.contains(typ[i])) if (abs(val[i])>mxF) mxF=abs(val[i]);
    if (XZ) for (i=0;i<nSect+nSect;i++) if (abs(pz[i]) > mxp) mxp = abs(pz[i]);
    if (XY) for (i=0;i<nSect+nSect;i++) if (abs(py[i]) > mxp) mxp = abs(py[i]);
    if (XZ) plotSystem(0,n,mxp,mxF);
    if (XY) plotSystem(1,n,mxp,mxF);
  }
  if (defl)
  { 
    mxv = 0.; mnv = 0.;
    if (XZ) { if (mxv < mx[0]) mxv = mx[0]; if (mnv > mn[0]) mnv = mn[0]; }
    if (XY) { if (mxv < mx[1]) mxv = mx[1]; if (mnv > mn[1]) mnv = mn[1]; }
    if (XZ) plotAndPrint(0,n,false,labels,mxv,mnv);
    if (XY) plotAndPrint(1,n,false,labels,mxv,mnv);
  }
  if (rota)
  {
    mxv = 0.; mnv = 0.;
    if (XZ) { if (mxv < mx[2]) mxv = mx[2]; if (mnv > mn[2]) mnv = mn[2]; }
    if (XY) { if (mxv < mx[3]) mxv = mx[3]; if (mnv > mn[3]) mnv = mn[3]; }
    if (XZ) plotAndPrint(2,n,true,labels,mxv,mnv);
    if (XY) plotAndPrint(3,n,true,labels,mxv,mnv);
  }
  if (bmom)
  {
    mxv = 0.; mnv = 0.;
    if (XZ) { if (mxv < mx[4]) mxv = mx[4]; if (mnv > mn[4]) mnv = mn[4]; }
    if (XY) { if (mxv < mx[5]) mxv = mx[5]; if (mnv > mn[5]) mnv = mn[5]; }
    if (XZ) plotAndPrint(4,n,true,labels,mxv,mnv);     
    if (XY) plotAndPrint(5,n,true,labels,mxv,mnv);     
  }
  if (shfo)
  {
    mxv = 0.; mnv = 0.;
    if (XZ) { if (mxv < mx[6]) mxv = mx[6]; if (mnv > mn[6]) mnv = mn[6]; }
    if (XY) { if (mxv < mx[7]) mxv = mx[7]; if (mnv > mn[7]) mnv = mn[7]; }
    if (XZ) plotAndPrint(6,n,true,labels,mxv,mnv); 
    if (XY) plotAndPrint(7,n,true,labels,mxv,mnv);
  }
  cout << "\n";

  return;
}







void BeamBending::plotSystem(int plane, int n, double mxp, double mxF)
{
  int i, j, m;

  double x1[2], x2[2], 
         y0 = h0[8+plane] + 0.5 * dh[8+plane], 
         x0 = 0., dx, y, yfact, xfact,xpos, scl = 1.,
         mfact = (double)(n * nSect) / (x[x.n-1] - x[0]);

  // system axis

  x1[0] =  .0; x1[1] = y0;
  x2[0] = 1.0; x2[1] = y0;
  plot.line(x1,x2);

  // distributed loads

  if (mxp > 1.e-12)
  {
    yfact = 0.4 * dh[8] / mxp;

    for (i=0; i<nSect; i++)
    {
      x1[1] = y0;
      m     = roundToInt(length[i] * mfact);
      dx    = 1. / (double) m;
      xfact = length[i] / (x[x.n-1] - x[0]);
      xpos  = 0.;
      for (j=0; j<m+1; j++)
      {
        if (plane == 0) y = pz[i+i] * (1. - xpos) + pz[i+i+1] * xpos;
        else            y = py[i+i] * (1. - xpos) + py[i+i+1] * xpos;

        x2[0] = x0 + xfact * xpos;
        x2[1] = y0 + yfact * y;
        plot.line(x1,x2);

        x1[0] = x2[0];
        x1[1] = y0;
        x2[0] = 0.;
        x2[1] = y0 - x2[1];
        plot.arrow(x1,x2,scl);
        x1[1] = y0 - x2[1];

        xpos += dx;
      }
      x0 += xfact;
    }
  }

  // bearings

  int i0 = 0;
  Vector<int> tmp, tmp2, tmp3;

  if (mxF > 1.e-12) yfact = .5 * dh[8] / mxF; else yfact = 0.; // set yfact for point loads

  x1[0] = 0.;
  x1[1] = y0;
  for (i=0; i<4; i++) tmp.append(typ[i]);

  cont1:

  tmp2.append(0); tmp2.append(4); tmp3.append(1); tmp3.append(5);
  if ((plane == 0 && tmp.containsAllOf(tmp2)) || (plane == 1 && tmp.containsAllOf(tmp3)))
  { drawSimpleSupp(1.,x1); drawMomHinge(1.,x1); }

  tmp2[0] = 0; tmp2[1] = 2; tmp3[0] = 1; tmp3[1] = 3;
  if ((plane == 0 && tmp.containsAllOf(tmp2)) || (plane == 1 && tmp.containsAllOf(tmp3)))
    drawClamp(1.,x1);
  
  tmp2[0] = 2; tmp2[1] = 4; tmp3[0] = 3; tmp3[1] = 5;
  if ((plane == 0 && tmp.containsAllOf(tmp2)) || (plane == 1 && tmp.containsAllOf(tmp3)))
    drawShearHinge(1.,x1);
  
  if ((plane == 0 && tmp.contains(6,&i)) || (plane == 1 && tmp.contains(7,&i))) 
    if (abs(val[i0+i]) > 1.e-12) drawPointLoad(1.,x1,val[i0+i]*yfact);

  if ((plane == 0 && tmp.contains(4,&i)) || (plane == 1 && tmp.contains(5,&i)))
    if (abs(val[i0+i]) > 1.e-12) COUT << "drawing of moment loads not yet implemented!\n\n";

  if (x1[0] < 0.1)
  {
    x1[0] = 1.;
    i0 = typ.n - 4;
    tmp.free();
    for (i=typ.n-4; i<typ.n; i++) tmp.append(typ[i]);
    goto cont1;
  }

  // hinges, intersection supports and point loads

  for (i=4; i<typ.n-4; i++)
  {
    x1[0] = x[i] / (x[x.n-1] - x[0]);

    if ((plane == 0 && typ[i] ==  0) || (plane == 1 && typ[i] ==  1)) drawSimpleSupp(1.,x1);

    if ((plane == 0 && typ[i] ==  4) || (plane == 1 && typ[i] ==  5)) drawMomHinge(1.,x1);

    if ((plane == 0 && typ[i] ==  6) || (plane == 1 && typ[i] ==  7)) drawShearHinge(1.,x1);

    if ((plane == 0 && typ[i] == 14) || (plane == 1 && typ[i] == 15))
      if (abs(val[i]) > 1.e-12) drawPointLoad(1.,x1,val[i]*yfact);

    if ((plane == 0 && typ[i] == 12) || (plane == 1 && typ[i] == 13))
      if (abs(val[i]) > 1.e-12) COUT << "drawing of moment loads not yet implemented!\n\n";
  }

  // coordinate system

  x1[0] = -.19; x1[1] = y0 + .04; 
  drawCoorSys(.7,x1,plane);

  return;
}











void BeamBending::plotAndPrint(int tp, int n, bool fill, bool label, double mxv, double mnv)
{
  int i, j, m;

  double x1[2], x2[2], 
         dy = mxv - mnv, 
         y0 = dh[tp],
         x0 = 0,
         xpos, 
         dx,
         y, 
         xfact = 1. / (x[x.n-1] - x[0]),
         yfact,
         mfact = (double)(n * nSect) * xfact,
         corr = .5 * ((mx[tp]+mn[tp]) - dy);

  if (dy >= 1.e-12)
  {
    yfact = - dh[tp] / dy;
    if (tp == 5) { yfact = - yfact; y0 = 0.; }
    y0 += h0[tp] - yfact * corr;

    for (i=0; i<nSect; i++)
    {
      x1[0] = x0;
      x1[1] = y0;
      xpos  = 0.;
      m     = roundToInt(length[i] * mfact);
      dx    = length[i] / (double) m;
      for (j=0; j<m+1; j++)
      {
        switch (tp)
        {
          case 0: y =  w(xpos,i); break;
          case 1: y =  v(xpos,i); break;
          case 2: y = bz(xpos,i); break;
          case 3: y = by(xpos,i); break;
          case 4: y = My(xpos,i); break;
          case 5: y = Mz(xpos,i); break;
          case 6: y = Sz(xpos,i); break;
          case 7: y = Sy(xpos,i); break;
        }
        x2[0] = x0 + xfact * xpos;
        x2[1] = y0 + yfact * y;

        if (j>0) plot.line(x1,x2);

        x1[0] = x2[0];
        if (fill) { x1[1] = y0; plot.line(x1,x2); }
        x1[1] = x2[1];

        xpos += dx;
      }
      x0 += length[i] * xfact;
    }
  }
  else y0 = h0[tp] + 0.5 * dh[tp]; 

  x1[0] =  .0; x1[1] = y0;
  x2[0] = 1.0; x2[1] = y0;
  plot.line(x1,x2);

  x1[0] = -.15;
  x1[1] = h0[tp] + dh[tp] * .5;

  double scl = .03 + (refHdW - .15) / .25 * .04;

  switch (tp)
  {
    case 0:  drawW(scl,x1); break;
    case 1:  drawV(scl,x1); break;
    case 2: drawBz(scl,x1); break;
    case 3: drawBy(scl,x1); break;
    case 4: drawMy(scl,x1); break;
    case 5: drawMz(scl,x1); break;
    case 6: drawSz(scl,x1); break;
    case 7: drawSy(scl,x1); break;
  }

  // print function

  int i0, p0, p1;

  double I2 = Iyy * Izz - Iyz * Iyz, l;

  for (i=0; i<nSect; i++)
  {
    l  = length[i];
    i0 = i * 8;
    p0 = i + i;
    p1 = p0 + 1;

    switch (tp)
    {
      case 0: printf("          w(x%d) = %+9.5f x^5\n",i+1,
                     ((py[p0]-py[p1])*Iyz-(pz[p0]-pz[p1])*Izz)/(120.*I2*l));
              printf("                  %+9.5f x^4\n",
                     -(py[p0]*Iyz-pz[p0]*Izz)/(24.*I2));
              printf("                  %+9.5f x^3\n",A[i0+0]/6.);
              printf("                  %+9.5f x^2\n",A[i0+1]*.5);
              printf("                  %+9.5f x\n",A[i0+2]);
              printf("                  %+9.5f\n",A[i0+3]);
              break;
      case 1: printf("          v(x%d) = %+9.5f x^5\n",i+1,
                     -((py[p0]-py[p1])*Iyy-(pz[p0]-pz[p1])*Iyz)/(120.*I2*l));
              printf("                  %+9.5f x^4\n",+(py[p0]*Iyy-pz[p0]*Iyz)/(24.*I2));
              printf("                  %+9.5f x^3\n",+A[i0+4]/6.);
              printf("                  %+9.5f x^2\n",+A[i0+5]*.5);
              printf("                  %+9.5f x\n",+A[i0+6]);
              printf("                  %+9.5f\n",+A[i0+7]);
              break;
      case 2: printf("         bz(x%d) = %+9.5f x^4\n",i+1,
                     -((py[p0]-py[p1])*Iyz-(pz[p0]-pz[p1])*Izz)/(24.*I2*l));
              printf("                  %+9.5f x^3\n",+(py[p0]*Iyz-pz[p0]*Izz)/(6.*I2));
              printf("                  %+9.5f x^2\n",-A[i0+0]*.5);
              printf("                  %+9.5f x\n",-A[i0+1]);
              printf("                  %+9.5f\n",-A[i0+2]);
              break;
      case 3: printf("         by(x%d) = %+9.5f x^4\n",i+1,
                     -((py[p0]-py[p1])*Iyy-(pz[p0]-pz[p1])*Iyz)/(24.*I2*l));
              printf("                  %+9.5f x^3\n",+(py[p0]*Iyy-pz[p0]*Iyz)/(6.*I2));
              printf("                  %+9.5f x^2\n",+A[i0+4]*.5);
              printf("                  %+9.5f x\n",+A[i0+5]);
              printf("                  %+9.5f\n",+A[i0+6]);
              break;
      case 4: printf("         My(x%d) = %+9.5f x^3\n",i+1,+(pz[p0]-pz[p1])/(6.*l));
              printf("                  %+9.5f x^2\n",-pz[p0]*.5);
              printf("                  %+9.5f x\n",-(A[i0]*Iyy+A[i0+4]*Iyz));
              printf("                  %+9.5f\n",-(A[i0+1]*Iyy+A[i0+5]*Iyz));
              break;
      case 5: printf("         Mz(x%d) = %+9.5f x^3\n",i+1,-(py[p0]-py[p1])/(6.*l));
              printf("                  %+9.5f x^2\n",+py[p0]*.5);
              printf("                  %+9.5f x\n",+(A[i0]*Iyz+A[i0+4]*Izz));
              printf("                  %+9.5f\n",+(A[i0+1]*Iyz+A[i0+5]*Izz));
              break;
      case 6: printf("         Sz(x%d) = %+9.5f x^2\n",i+1,+(pz[p0]-pz[p1])*.5/l);
              printf("                  %+9.5f x\n",-pz[p0]);
              printf("                  %+9.5f\n",-(A[i0]*Iyy+A[i0+4]*Iyz));
              break;
      case 7: printf("         Sy(x%d) = %+9.5f x^2\n",i+1,+(py[p0] - py[p1])*.5/l);
              printf("                  %+9.5f x\n",-py[p0]);
              printf("                  %+9.5f\n",-(A[i0]*Iyz+A[i0+4]*Izz));
              break;
    }
  }
  cout << "\n";

  if (!label) return;

  // plot labels

  char tmp[30];

  x0 = 0;

  if (dy >= 1.e-12)
  {
    for (i=0; i<nSect; i++)
    {
      switch (tp)
      {
        case 0: y =  w(0.,i); break;
        case 1: y =  v(0.,i); break;
        case 2: y = bz(0.,i); break;
        case 3: y = by(0.,i); break;
        case 4: y = My(0.,i); break;
        case 5: y = Mz(0.,i); break;
        case 6: y = Sz(0.,i); break;
        case 7: y = Sy(0.,i); break;
      }
      x1[0] = x0;
      x1[1] = y0 + yfact * y * .7;

      sprintf(tmp,"%+8.4f",y);
      j = 0; while (tmp[j] == ' ') j++;
      plot.putText(x1,&(tmp[j]),5,true);

      x0 += length[i] * xfact;

      switch (tp)
      {
        case 0: y =  w(length[i],i); break;
        case 1: y =  v(length[i],i); break;
        case 2: y = bz(length[i],i); break;
        case 3: y = by(length[i],i); break;
        case 4: y = My(length[i],i); break;
        case 5: y = Mz(length[i],i); break;
        case 6: y = Sz(length[i],i); break;
        case 7: y = Sy(length[i],i); break;
      }
      x1[0] = x0;
      x1[1] = y0 + yfact * y * .7;

      sprintf(tmp,"%+8.4f",y);
      j = 0; while (tmp[j] == ' ') j++;
      plot.putText(x1,&(tmp[j]),5,true);
    }
  }
  else 
  {
    x1[0] = 0.5;
    x1[1] = y0;
    sprintf(tmp,"%+8.4f",0.);
    j = 0; while (tmp[j] == ' ') j++;
    plot.putText(x1,&(tmp[j]),5,true);
  }

  return;
}











double BeamBending::w(double xpos, int i)
{
  double I2 = Iyy * Izz - Iyz * Iyz, 
         xpos2 = xpos * xpos, 
         xpos3 = xpos2 * xpos, 
         xpos4 = xpos3 * xpos, 
         xpos5 = xpos4 * xpos, 
         l = length[i];

  int i0 = i * 8, p0 = i+i, p1 = p0+1;

  return + ((py[p0]-py[p1]) * Iyz - (pz[p0]-pz[p1]) * Izz) * xpos5 / (120. * I2 * l)
         - (py[p0] * Iyz - pz[p0] * Izz) * xpos4 / (24. * I2)
         + A[i0+0] * xpos3 / 6.
         + A[i0+1] * xpos2 * .5
         + A[i0+2] * xpos
         + A[i0+3];
}







double BeamBending::v(double xpos, int i)
{
  double I2 = Iyy * Izz - Iyz * Iyz, 
         xpos2 = xpos * xpos, 
         xpos3 = xpos2 * xpos, 
         xpos4 = xpos3 * xpos, 
         xpos5 = xpos4 * xpos, 
         l = length[i];

  int i0 = i * 8, p0 = i+i, p1 = p0+1;

  return - ((py[p0]-py[p1]) * Iyy - (pz[p0]-pz[p1]) * Iyz) * xpos5 / (120. * I2 * l)
         + (py[p0] * Iyy - pz[p0] * Iyz) * xpos4 / (24. * I2)
         + A[i0+4] * xpos3 / 6.
         + A[i0+5] * xpos2 * .5
         + A[i0+6] * xpos
         + A[i0+7];
}







double BeamBending::bz(double xpos, int i)
{
  double I2 = Iyy * Izz - Iyz * Iyz, 
         xpos2 = xpos * xpos, 
         xpos3 = xpos2 * xpos, 
         xpos4 = xpos3 * xpos, 
         l = length[i];

  int i0 = i * 8, p0 = i+i, p1 = p0+1;

  return - ((py[p0]-py[p1]) * Iyz - (pz[p0]-pz[p1]) * Izz) * xpos4 / (24. * I2 * l)
         + (py[p0] * Iyz - pz[p0] * Izz) * xpos3 / (6. * I2)
         - A[i0+0] * xpos2 * .5
         - A[i0+1] * xpos
         - A[i0+2];
}







double BeamBending::by(double xpos, int i)
{
  double I2 = Iyy * Izz - Iyz * Iyz, 
         xpos2 = xpos * xpos, 
         xpos3 = xpos2 * xpos, 
         xpos4 = xpos3 * xpos, 
         l = length[i];

  int i0 = i * 8, p0 = i+i, p1 = p0+1;

  return - ((py[p0]-py[p1]) * Iyy - (pz[p0]-pz[p1]) * Iyz) * xpos4 / (24. * I2 * l)
         + (py[p0] * Iyy - pz[p0] * Iyz) * xpos3 / (6. * I2)
         + A[i0+4] * xpos2 * .5
         + A[i0+5] * xpos
         + A[i0+6];
}







double BeamBending::My(double xpos, int i)
{
  double I2 = Iyy * Izz - Iyz * Iyz, 
         xpos2 = xpos * xpos, xpos3 = xpos2 * xpos, l = length[i];

  int i0 = i * 8, p0 = i+i, p1 = p0+1;

  return + (pz[p0] - pz[p1]) * xpos3 / (6.*l)
         - pz[p0] * xpos2 * .5
         - (A[i0] * Iyy + A[i0+4] * Iyz) * xpos
         - (A[i0+1] * Iyy + A[i0+5] * Iyz);
}







double BeamBending::Mz(double xpos, int i)
{
  double I2 = Iyy * Izz - Iyz * Iyz, 
         xpos2 = xpos * xpos, xpos3 = xpos2 * xpos, l = length[i];

  int i0 = i * 8, p0 = i+i, p1 = p0+1;

  return - (py[p0] - py[p1]) * xpos3 / (6.*l)
         + py[p0] * xpos2 * .5
         + (A[i0] * Iyz + A[i0+4] * Izz) * xpos
         + (A[i0+1] * Iyz + A[i0+5] * Izz);
}







double BeamBending::Sz(double xpos, int i)
{
  double I2 = Iyy * Izz - Iyz * Iyz, xpos2 = xpos * xpos, l = length[i];

  int i0 = i * 8, p0 = i+i, p1 = p0+1;

  return + (pz[p0] - pz[p1]) * xpos2 * .5 / l
         - pz[p0] * xpos
         - (A[i0] * Iyy + A[i0+4] * Iyz);
}







double BeamBending::Sy(double xpos, int i)
{
  double I2 = Iyy * Izz - Iyz * Iyz, xpos2 = xpos * xpos, l = length[i];

  int i0 = i * 8, p0 = i+i, p1 = p0+1;

  return + (py[p0] - py[p1]) * xpos2 * .5 / l
         - py[p0] * xpos
         - (A[i0] * Iyz + A[i0+4] * Izz);
}






void BeamBending::drawCoorSys(double scl, double *xx, int plane)
{
  char tmp[4]; tmp[1] = '\0';

  double x1[3], x2[3], r = .028 * scl, sqrt05r = sqrt(.125) * r, d = 1.;

  int c = plot.currStdColour;

  x1[0] = xx[0] + .1 * scl; x1[1] = xx[1]; x2[0] = +.1 * scl; x2[1] = 0.; plot.arrow(x1,x2,d);
  x1[0] = xx[0]; x1[1] = xx[1] - .1 * scl; x2[0] = .0; x2[1] = -.1 * scl; plot.arrow(x1,x2,d);
  plot.setColour(8); plot.fillCircle(xx,r);
  plot.setColour(c); plot.circle(xx,r);
  if (plane == 0) { x1[0] = .3 * r; plot.fillCircle(xx,x1[0]); tmp[0] = 'y'; }
  else            
  { 
    x1[0] = xx[0] - sqrt05r; x1[1] = xx[1] - sqrt05r;
    x2[0] = xx[0] + sqrt05r; x2[1] = xx[1] + sqrt05r; plot.line(x1,x2);
    x1[0] = xx[0] - sqrt05r; x1[1] = xx[1] + sqrt05r;
    x2[0] = xx[0] + sqrt05r; x2[1] = xx[1] - sqrt05r; plot.line(x1,x2); tmp[0] = 'z';
  }
  x1[0] = xx[0] - .015 * scl; x1[1] = xx[1] + .015 * scl;
  plot.putText(x1,tmp,3,true);
  x1[0] = xx[0] + .105 * scl; x1[1] = xx[1] + .015 * scl;
  tmp[0] = 'x'; plot.putText(x1,tmp,0,true);
  x1[0] = xx[0] - .015 * scl; x1[1] = xx[1] - .105 * scl;
  if (plane == 0) tmp[0] = 'z'; else tmp[0] = 'y'; plot.putText(x1,tmp,9,true);

  return;
}






void BeamBending::drawClamp(double scl, double *xx)
{
  float *pnt = plot.polyPnt, f1 = .005 * scl, f2 = .025 * scl;

  pnt[0] = (float)xx[0] - f1; pnt[ 1] = (float)xx[1] - f2;
  pnt[3] = (float)xx[0] - f1; pnt[ 4] = (float)xx[1] + f2; 
  pnt[6] = (float)xx[0] + f1; pnt[ 7] = (float)xx[1] + f2; 
  pnt[9] = (float)xx[0] + f1; pnt[10] = (float)xx[1] - f2; plot.poly(4);

  return;
}






void BeamBending::drawSimpleSupp(double scl, double *xx)
{
  float *pnt = plot.polyPnt;

  int c = plot.currStdColour;

  pnt[0] = (float)xx[0]         ; pnt[ 1] = (float)xx[1]         ;
  pnt[3] = (float)xx[0]+.015*scl; pnt[ 4] = (float)xx[1]-.025*scl;
  pnt[6] = (float)xx[0]-.015*scl; pnt[ 7] = (float)xx[1]-.025*scl; 
  pnt[9] = (float)xx[0]         ; pnt[10] = (float)xx[1]         ;

  plot.setColour(8);
  plot.poly(4); 
  plot.setColour(c);
  for (int i=0; i<3; i++) plot.line(&(pnt[i*3]),&(pnt[i*3+3]));

  return;
}






void BeamBending::drawMomHinge(double scl, double *xx)
{
  int c = plot.currStdColour;

  double d = 0.0125 * scl;

  plot.setColour(8); plot.fillCircle(xx,d);
  plot.setColour(c); plot.circle(xx,d);

  return;
}






void BeamBending::drawShearHinge(double scl, double *xx)
{
  float *pnt = plot.polyPnt, f1 = .007*scl, f2 = .025*scl;

  int c = plot.currStdColour;

  pnt[0] = (float)xx[0] - f1; pnt[ 1] = (float)xx[1] - f2;
  pnt[3] = (float)xx[0] - f1; pnt[ 4] = (float)xx[1] + f2; 
  pnt[6] = (float)xx[0] + f1; pnt[ 7] = (float)xx[1] + f2; 
  pnt[9] = (float)xx[0] + f1; pnt[10] = (float)xx[1] - f2;

  plot.setColour(8);
  plot.poly(4);
  plot.setColour(c);
  plot.line(&(pnt[0]),&(pnt[3]));
  plot.line(&(pnt[6]),&(pnt[9]));

  return;
}






void BeamBending::drawPointLoad(double scl, double *xx, double F)
{
  double x1[3];

  x1[0] = 0.;  
  if (xx[0] < .5*minl/(x[x.n-1]-x[0])) x1[1] = F; else x1[1] = - F;  
  plot.arrow(xx,x1,scl);

  return;
}







void BeamBending::drawW(double scl, double *xx)
{
  double x1[3], x2[3]; 

  x1[0] = scl + scl;
  plot.circle(xx,x1[0]);

  x1[0] = +.50 * scl + xx[0]; x1[1] = + .4 * scl + xx[1];
  x2[0] = +.25 * scl + xx[0]; x2[1] = - .4 * scl + xx[1]; plot.line(x1,x2);
  x1[0] =            + xx[0]; x1[1] = + .1 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = -.25 * scl + xx[0]; x2[1] = - .4 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.50 * scl + xx[0]; x1[1] = + .4 * scl + xx[1]; plot.line(x1,x2);

  return;
}






void BeamBending::drawV(double scl, double *xx)
{
  double x1[3], x2[3]; 

  x1[0] = scl + scl;
  plot.circle(xx,x1[0]);

  x1[0] = +.4 * scl + xx[0]; x1[1] = + .4 * scl + xx[1];
  x2[0] =             xx[0]; x2[1] = - .4 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.4 * scl + xx[0]; x1[1] = + .4 * scl + xx[1]; plot.line(x1,x2);

  return;
}






void BeamBending::drawBz(double scl, double *xx)
{
  double x1[3], x2[3]; 

  x1[0] = scl + scl;
  plot.circle(xx,x1[0]);

  x1[0] = -.4 * scl + xx[0]; x1[1] = - .7 * scl + xx[1];
  x2[0] = -.4 * scl + xx[0]; x2[1] = + .6 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.3 * scl + xx[0]; x1[1] = + .7 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = -.1 * scl + xx[0]; x2[1] = + .7 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.0 * scl + xx[0]; x1[1] = + .6 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = -.0 * scl + xx[0]; x2[1] = + .4 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.1 * scl + xx[0]; x1[1] = + .3 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = +.1 * scl + xx[0]; x2[1] = + .1 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = +.1 * scl + xx[0]; x1[1] = - .2 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = +.0 * scl + xx[0]; x2[1] = - .3 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.2 * scl + xx[0]; x1[1] = - .3 * scl + xx[1]; plot.line(x1,x2);

  x1[0] = +.3 * scl + xx[0]; x1[1] = - .2 * scl + xx[1];
  x2[0] = +.6 * scl + xx[0]; x2[1] = - .2 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = +.3 * scl + xx[0]; x1[1] = - .5 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = +.6 * scl + xx[0]; x2[1] = - .5 * scl + xx[1]; plot.line(x1,x2);

  return;
}






void BeamBending::drawBy(double scl, double *xx)
{
  double x1[3], x2[3]; 

  x1[0] = scl + scl;
  plot.circle(xx,x1[0]);

  x1[0] = -.4 * scl + xx[0]; x1[1] = - .7 * scl + xx[1];
  x2[0] = -.4 * scl + xx[0]; x2[1] = + .6 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.3 * scl + xx[0]; x1[1] = + .7 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = -.1 * scl + xx[0]; x2[1] = + .7 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.0 * scl + xx[0]; x1[1] = + .6 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = -.0 * scl + xx[0]; x2[1] = + .4 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.1 * scl + xx[0]; x1[1] = + .3 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = +.1 * scl + xx[0]; x2[1] = + .1 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = +.1 * scl + xx[0]; x1[1] = - .2 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = +.0 * scl + xx[0]; x2[1] = - .3 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.2 * scl + xx[0]; x1[1] = - .3 * scl + xx[1]; plot.line(x1,x2);

  x1[0] = +.3 * scl + xx[0]; x1[1] = - .20 * scl + xx[1];
  x2[0] = +.5 * scl + xx[0]; x2[1] = - .50 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = +.7 * scl + xx[0]; x1[1] = - .20 * scl + xx[1];
  x2[0] = +.3 * scl + xx[0]; x2[1] = - .80 * scl + xx[1]; plot.line(x1,x2);

  return;
}






void BeamBending::drawMy(double scl, double *xx)
{
  double x1[3], x2[3]; 

  x1[0] = scl + scl;
  plot.circle(xx,x1[0]);

  x1[0] = -.60 * scl + xx[0]; x1[1] = - .40 * scl + xx[1];
  x2[0] = -.60 * scl + xx[0]; x2[1] = + .60 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.25 * scl + xx[0]; x1[1] = + .10 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = +.10 * scl + xx[0]; x2[1] = + .60 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = +.10 * scl + xx[0]; x1[1] = - .40 * scl + xx[1]; plot.line(x1,x2);

  x1[0] = +.3 * scl + xx[0]; x1[1] = - .20 * scl + xx[1];
  x2[0] = +.5 * scl + xx[0]; x2[1] = - .50 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = +.7 * scl + xx[0]; x1[1] = - .20 * scl + xx[1];
  x2[0] = +.3 * scl + xx[0]; x2[1] = - .80 * scl + xx[1]; plot.line(x1,x2);

  return;
}






void BeamBending::drawMz(double scl, double *xx)
{
  double x1[3], x2[3]; 

  x1[0] = scl + scl;
  plot.circle(xx,x1[0]);

  x1[0] = -.60 * scl + xx[0]; x1[1] = - .40 * scl + xx[1];
  x2[0] = -.60 * scl + xx[0]; x2[1] = + .60 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.25 * scl + xx[0]; x1[1] = + .10 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = +.10 * scl + xx[0]; x2[1] = + .60 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = +.10 * scl + xx[0]; x1[1] = - .40 * scl + xx[1]; plot.line(x1,x2);

  x1[0] = +.3 * scl + xx[0]; x1[1] = - .2 * scl + xx[1];
  x2[0] = +.6 * scl + xx[0]; x2[1] = - .2 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = +.3 * scl + xx[0]; x1[1] = - .5 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = +.6 * scl + xx[0]; x2[1] = - .5 * scl + xx[1]; plot.line(x1,x2);

  return;
}






void BeamBending::drawSz(double scl, double *xx)
{
  double x1[3], x2[3]; 

  x1[0] = scl + scl;
  plot.circle(xx,x1[0]);

  x1[0] = +.1 * scl + xx[0]; x1[1] = + .5 * scl + xx[1];
  x2[0] = -.0 * scl + xx[0]; x2[1] = + .6 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.5 * scl + xx[0]; x1[1] = + .6 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = -.6 * scl + xx[0]; x2[1] = + .5 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.6 * scl + xx[0]; x1[1] = + .2 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = -.5 * scl + xx[0]; x2[1] = + .1 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.0 * scl + xx[0]; x1[1] = + .1 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = +.1 * scl + xx[0]; x2[1] = + .0 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = +.1 * scl + xx[0]; x1[1] = - .3 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = -.0 * scl + xx[0]; x2[1] = - .4 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.5 * scl + xx[0]; x1[1] = - .4 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = -.6 * scl + xx[0]; x2[1] = - .3 * scl + xx[1]; plot.line(x1,x2);

  x1[0] = +.3 * scl + xx[0]; x1[1] = - .2 * scl + xx[1];
  x2[0] = +.6 * scl + xx[0]; x2[1] = - .2 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = +.3 * scl + xx[0]; x1[1] = - .5 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = +.6 * scl + xx[0]; x2[1] = - .5 * scl + xx[1]; plot.line(x1,x2);

  return;
}






void BeamBending::drawSy(double scl, double *xx)
{
  double x1[3], x2[3]; 

  x1[0] = scl + scl;
  plot.circle(xx,x1[0]);

  x1[0] = +.1 * scl + xx[0]; x1[1] = + .5 * scl + xx[1];
  x2[0] = -.0 * scl + xx[0]; x2[1] = + .6 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.5 * scl + xx[0]; x1[1] = + .6 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = -.6 * scl + xx[0]; x2[1] = + .5 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.6 * scl + xx[0]; x1[1] = + .2 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = -.5 * scl + xx[0]; x2[1] = + .1 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.0 * scl + xx[0]; x1[1] = + .1 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = +.1 * scl + xx[0]; x2[1] = + .0 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = +.1 * scl + xx[0]; x1[1] = - .3 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = -.0 * scl + xx[0]; x2[1] = - .4 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = -.5 * scl + xx[0]; x1[1] = - .4 * scl + xx[1]; plot.line(x1,x2);
  x2[0] = -.6 * scl + xx[0]; x2[1] = - .3 * scl + xx[1]; plot.line(x1,x2);

  x1[0] = +.3 * scl + xx[0]; x1[1] = - .20 * scl + xx[1];
  x2[0] = +.5 * scl + xx[0]; x2[1] = - .50 * scl + xx[1]; plot.line(x1,x2);
  x1[0] = +.7 * scl + xx[0]; x1[1] = - .20 * scl + xx[1];
  x2[0] = +.3 * scl + xx[0]; x2[1] = - .80 * scl + xx[1]; plot.line(x1,x2);

  return;
}









void BeamBending::addW0(int rw, int iSect, double c, double fct)
{
  brhs[rw-1] += c;
  bmtx(rw,iSect*8+4) += fct;
  return;
}






void BeamBending::addW1(int rw, int iSect, double c, double fct)
{
  int i0 = iSect * 8, p0 = iSect + iSect, p1 = p0 + 1;
  double l = length[iSect], l2 = l * l, l3 = l2 * l, l4 = l3 * l,
         I2 = Iyy * Izz - Iyz * Iyz;
  brhs[rw-1] += c - fct * l4*(Izz*(4.*pz[p0]+pz[p1])-Iyz*(4.*py[p0]+py[p1]))/(120.*I2);
  bmtx(rw,i0+1) += fct * l3/6.;
  bmtx(rw,i0+2) += fct * l2*.5;
  bmtx(rw,i0+3) += fct * l;
  bmtx(rw,i0+4) += fct;
  return;
}






void BeamBending::addV0(int rw, int iSect, double c, double fct)
{
  brhs[rw-1] += c;
  bmtx(rw,iSect*8+8) += fct;
  return;
}






void BeamBending::addV1(int rw, int iSect, double c, double fct)
{
  int i0 = iSect * 8, p0 = iSect + iSect, p1 = p0 + 1;
  double l = length[iSect], l2 = l * l, l3 = l2 * l, l4 = l3 * l,
         I2 = Iyy * Izz - Iyz * Iyz;
  brhs[rw-1] += c-fct*l4*(Iyy*(4.*py[p0]+py[p1])-Iyz*(4.*pz[p0]+pz[p1]))/(120.*I2);
  bmtx(rw,i0+5) += fct * l3/6.;
  bmtx(rw,i0+6) += fct * l2*.5;
  bmtx(rw,i0+7) += fct * l;
  bmtx(rw,i0+8) += fct;
  return;
}






void BeamBending::addBz0(int rw, int iSect, double c, double fct)
{
  brhs[rw-1] += c;
  bmtx(rw,iSect*8+3) -= fct;
  return;
}






void BeamBending::addBz1(int rw, int iSect, double c, double fct)
{
  int i0 = iSect * 8, p0 = iSect + iSect, p1 = p0 + 1;
  double l = length[iSect], l2 = l * l, l3 = l2 * l,
         I2 = Iyy * Izz - Iyz * Iyz;
  brhs[rw-1] += c-fct*l3*(Iyz*(3.*py[p0]+py[p1])-Izz*(3.*pz[p0]+pz[p1]))/(24.*I2);
  bmtx(rw,i0+1) -= fct * l2*.5;
  bmtx(rw,i0+2) -= fct * l;
  bmtx(rw,i0+3) -= fct;
  return;
}






void BeamBending::addBy0(int rw, int iSect, double c, double fct)
{
  brhs[rw-1] += c;
  bmtx(rw,iSect*8+7) += fct;
  return;
}






void BeamBending::addBy1(int rw, int iSect, double c, double fct)
{
  int i0 = iSect * 8, p0 = iSect + iSect, p1 = p0 + 1;
  double l = length[iSect], l2 = l * l, l3 = l2 * l,
         I2 = Iyy * Izz - Iyz * Iyz;
  brhs[rw-1] += c-fct*l3*(Iyy*(3.*py[p0]+py[p1])-Iyz*(3.*pz[p0]+pz[p1]))/(24.*I2);
  bmtx(rw,i0+5) += fct * l2*.5;
  bmtx(rw,i0+6) += fct * l;
  bmtx(rw,i0+7) += fct;
  return;
}






void BeamBending::addMy0(int rw, int iSect, double c, double fct)
{
  brhs[rw-1] += c;
  bmtx(rw,iSect*8+2) -= fct * Iyy; 
  bmtx(rw,iSect*8+6) -= fct * Iyz;
  return;
}






void BeamBending::addMy1(int rw, int iSect, double c, double fct)
{
  int i0 = iSect * 8, p0 = iSect + iSect, p1 = p0 + 1;
  double l = length[iSect], l2 = l * l;
  brhs[rw-1] += c + fct * l2 * (2.*pz[p0]+pz[p1]) / 6.;
  bmtx(rw,i0+1) -= fct * Iyy * l; 
  bmtx(rw,i0+2) -= fct * Iyy; 
  bmtx(rw,i0+5) -= fct * Iyz * l; 
  bmtx(rw,i0+6) -= fct * Iyz;
  return;
}






void BeamBending::addMz0(int rw, int iSect, double c, double fct)
{
  brhs[rw-1] += c;
  bmtx(rw,iSect*8+2) += fct * Iyz; 
  bmtx(rw,iSect*8+6) += fct * Izz;
  return;
}






void BeamBending::addMz1(int rw, int iSect, double c, double fct)
{
  int i0 = iSect * 8, p0 = iSect + iSect, p1 = p0 + 1;
  double l = length[iSect], l2 = l * l;
  brhs[rw-1] += c - fct * l2 * (2.*py[p0]+py[p1]) / 6.;
  bmtx(rw,i0+1) += fct * Iyz * l; 
  bmtx(rw,i0+2) += fct * Iyz; 
  bmtx(rw,i0+5) += fct * Izz * l; 
  bmtx(rw,i0+6) += fct * Izz;
  return;
}






void BeamBending::addSz0(int rw, int iSect, double c, double fct)
{
  brhs[rw-1] += c;
  bmtx(rw,iSect*8+1) -= fct * Iyy; 
  bmtx(rw,iSect*8+5) -= fct * Iyz;
  return;
}






void BeamBending::addSz1(int rw, int iSect, double c, double fct)
{
  int i0 = iSect * 8, p0 = iSect + iSect, p1 = p0 + 1;
  double l = length[iSect];
  brhs[rw-1] += c + fct * l * (pz[p0]+pz[p1]) * .5;
  bmtx(rw,i0+1) -= fct * Iyy; 
  bmtx(rw,i0+5) -= fct * Iyz;
  return;
}






void BeamBending::addSy0(int rw, int iSect, double c, double fct)
{
  brhs[rw-1] += c;
  bmtx(rw,iSect*8+1) -= fct * Iyz; 
  bmtx(rw,iSect*8+5) -= fct * Izz;
  return;
}






void BeamBending::addSy1(int rw, int iSect, double c, double fct)
{
  int i0 = iSect * 8, p0 = iSect + iSect, p1 = p0 + 1;
  double l = length[iSect];
  brhs[rw-1] += c + fct * l * (py[p0]+py[p1]) * .5;
  bmtx(rw,i0+1) -= fct * Iyz; 
  bmtx(rw,i0+5) -= fct * Izz;
  return;
}

