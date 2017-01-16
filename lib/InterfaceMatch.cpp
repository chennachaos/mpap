
#include <iostream>

#include "FunctionsProgram.h"
#include "InterfaceMatch.h"
#include "DomainTree.h"
#include "DomainType.h"
#include "DataBlockTemplate.h"
#include "MathGeom.h"
#include "MpapTime.h"
#include "FreeSurface.h"



extern DomainTree domain;
extern MpapTime mpapTime;






InterfaceMatch::InterfaceMatch(void)
{                                                  
  show = false;

  tol = -2.;

  ndm = -1;

  // add new type
  
  DomainType *interfaceMatch = domain.newType(INTERFACEMATCH,ROOTDOMAIN);
  
  if (interfaceMatch == NULL) return;  // domain type already exists

  interfaceMatch->key.addNew("domains","boundary conditions","show","control");

  return;
}




	                                          
InterfaceMatch::~InterfaceMatch(void)
{                                                 
  return;
}







void InterfaceMatch::readInputData(std::ifstream &Ifile, MyString &line)
{
  MyString tmpl, *word;
 
  char tmp[30], fct[] = "InterfaceMatch::readInputData";

  char *domKey[]        = DOMAIN_KEY,
       *domModeName[]   = { "RESOLVE", "ELIMINATE", "ELIMINATE_MESH", "RETURN_FORCES", 
                            "RETURN_DISPLACEMENTS", "DUMMY" },
       *uTypeName[]     = { "VELOCITY", "DISPLACEMENT", NULL },
       *initGuessName[] = { "KEEP_VELOCITY", "KEEP_DISPLACEMENT", NULL };

  int nw, i, j, nn;
	
  Vector<double> dTmp;
  Vector<int>    iTmp, lTmp;
  MyStringList   sTmp;

  DataBlockTemplate t1, t2;

  switch (domain[INTERFACEMATCH].key.whichBegins(line))
  {
    case  0: cout << "     INTERFACEMATCH: reading domains ...\n\n";

	     sprintf(tmp,"1l 1i 1s 1f");  
	     
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
	     
             t1.initialise(tmpl,domKey);
	     t2.initialise(tmp,domKey);
             t1.expandToMatch(t2);
	     
             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
		     prgError(1,fct,"error in 'domains'!");
	    
             domType.setDim(lTmp.n);
             domId.setDim(lTmp.n);
             domMode.setDim(lTmp.n);
             domPtr.setDim(lTmp.n);
             domReacTime.setDim(lTmp.n);
 
             for (i=0; i<lTmp.n; i++)
             {
               domType[i]     = lTmp[i];
               domId[i]       = iTmp[i] - 1;
               domMode[i]     = sTmp[i].which(domModeName);
               domReacTime[i] = dTmp[i];

               if (domMode[i] < 0) prgError(1,fct,"invalid domain mode!");

               for (j=0; j<i; j++) 
                 if (domType[i] == domType[j] && domId[i] == domId[j]) 
                   prgError(2,fct,"domain found twice!");
             }

	     break;

    case  1: cout << "     INTERFACEMATCH: reading boundary conditions ...\n\n";

             if (ndm < 1) prgError(1,fct,"'boundary conditions' must be preceded by 'dimensions'");

             sprintf(tmp,"1l %di", ndf+ndm+2);
             
             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);
         
             t1.initialise(tmpl,domKey);
             t2.initialise(tmp,domKey);
             t1.expandToMatch(t2);

             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
                     prgError(1,fct,"error in 'domain data'!");

             nn = ndf + ndm + 3;
             nw = nn - 1;

             bcTmp.setDim(lTmp.n * nn);

             for (i=0; i<lTmp.n; i++)
             {
               bcTmp[i*nn  ] = lTmp[i];
               bcTmp[i*nn+1] = iTmp[i*nw  ] - 1;
               bcTmp[i*nn+2] = iTmp[i*nw+1];
               for (j=3; j<nn; j++) 
               { 
                 bcTmp[i*nn+j] = iTmp[i*nw+j-1];
                 if (bcTmp[i*nn+j] == 0) bcTmp[i*nn+j] = 1; 
                                    else bcTmp[i*nn+j] = 0;
               }
             }  

	     break;

    case  2: show = true;

             line.getNextLine(Ifile);

             break;

    case  3: cout << "     INTERFACEMATCH: reading control ...\n\n";

             if (tol > -1)    prgError(1,fct,"'control' has already been read!");
	     
	     line.getNextLine(Ifile);
	     
	     nw = line.split(&word);
            
	     if (nw < 5)                    prgError(1,fct,"input error in 'control'!");
		     
             if (!word[0].toDbl(&tol))      prgError(2,fct,"input error in 'control'!");

             uType = word[1].which(uTypeName);

             if (uType < 0) prgError(3,fct,"invalid uType name in 'control'!");

             initGuess = word[2].which(initGuessName);

             if (initGuess < 0) prgError(3,fct,"invalid initial guess name in 'control'!");

             if (!word[3].toInt(&tis))      prgError(5,fct,"input error in 'control'!");

             for (i=4; i<nw; i++) 
             {
               if (!word[i].toDbl(&td[i-4])) prgError(6,fct,"input error in 'control'!");
             }

             for (i=0; i<nw; i++) word[i].free(); delete [] word;
	     
       	     line.getNextLine(Ifile);

	     break;

    case -1: // go and inherit from DOMAIN
	     
	     this->Domain::readInputData(Ifile,line); 
	     
	     break;
  }
 
  return;
}








void InterfaceMatch::prepareInputData(void)
{
  char fct[] = "InterfaceMatch::prepareInputData";

  // go and inherit from ancestors

  Domain::prepareInputData();
  
  std::cout << "     INTERFACEMATCH: prepare input data ...\n\n";
  
// intialise whatever

  return;
}









void InterfaceMatch::prepareInteractions(void)
{
  // go and inherit from ancestors

  Domain::prepareInteractions();

  cout << "     INTERFACEMATCH: preparing interactions ...\n\n"; 

  char fct[] = "InterfaceMatch::prepareInteractions";
 
  int dm, i, j;

  // set domain pointers and check and set domain isConnected flag
 
  for (i=0; i<domType.n; i++)
  {
    if (domain.nDomainOfType(domType[i]) <= domId[i]) prgError(1,fct,"domain does not exist!");

    domPtr[i] = (FiniteElementBVPWNI*)&(domain(domType[i],domId[i]));

    if (domPtr[i]->isConnected)

      prgError(2,fct,"domain is already connected to another interface!");

    else domPtr[i]->isConnected = true;
  }
  //for (i=0; i<domType.n; i++) cout << domain.key(domPtr[i]) << "\n";

  // set ndm, ndf and dimensions of intfNodeToFreeNode

  freeNodeToIntfNode.setDim(domType.n);

  for (dm=0; dm<domType.n; dm++)
  {
    if (ndm != domPtr[dm]->ndm) prgError(2,fct,"domains have different ndm!");

    for (i=0; i<freeNode(dm).n; i++)

      if (ndf < freeNode(dm)[i].dofU->n) prgError(2,fct,"domains have inconsistent ndf!");

    freeNodeToIntfNode[dm].setDim(freeNode(dm).n);
  }

  // get nodal coordinates and match free nodes, set numnp

  getCoorAndMatchNodes();

  // set boundary conditions, generate idu and idx and set Lagrangian nodes for each domain

  generateIdux();

  setLagrangianNodes();

  nequ = 0;
  neqx = 0;

  for (i=0; i<numnp*ndf; i++) if (idu.x[i] == 1) idu.x[i] = ++nequ; else idu.x[i] = 0;
  for (i=0; i<numnp*ndm; i++) if (idx.x[i] == 1) idx.x[i] = ++neqx; else idx.x[i] = 0;

  // initialise kinematic data

  initKinData();  

  // calculate traction force time instants

  for (dm=0; dm<domType.n; dm++)
  {
    if (domReacTime[dm] < 1.e-10)
    {  
      if (domPtr[dm]->tis != 2)
        prgError(3,fct,"use generalised-alpha method or set traction froce time instant!");

      domReacTime[dm] = 1. / (1. + domPtr[dm]->td[0]);
    }
  }

  return;
}









void InterfaceMatch::getCoorAndMatchNodes(void)
{
  char fct[] = "InterfaceMatch::getCoorAndMatchNodes";

  int dm, i, j, j0, n; 

  // generate nodes from first domain

  freeNodeToIntfNode[0].setDim(freeNode(0).n);

  for (n=0; n<freeNode(0).n; n++) 
  {
    freeNodeToIntfNode[0][n] = n;

    node.add(new MoreData(n,freeNode(0)[n]));
  }

  for (n=0; n<freeNode(0).n; n++) 
    for (i=0; i<freeNode(0)[n].toFree[0].n; i++)
      node[n].depPtr.append(&(node[freeNode(0)[n].toFree[0][i].nd]));

  numnp = n;

  // match or add nodes from other domains

  compCount = 0;

  for (dm=1; dm<domType.n; dm++)
  {
    searchFrom.free();

    changed.free();

    F2I = freeNodeToIntfNode[dm].x;

    for (i=0; i<freeNode(dm).n; i++) F2I[i] = -1;

    findNodeMatches(dm);

    addNonMatchingNodes(dm);

    extendConnectivity(dm);

    for (i=0; i<numnp; i++) node[i].fnd = false;
  }

  // perform match count

  Vector<int> matchCnt;

  for (i=0; i<numnp; i++) matchCnt.append(0);

  j = 0;

  for (dm=0; dm<domType.n; dm++)
  {
    F2I = freeNodeToIntfNode[dm].x;

    for (i=0; i<freeNodeToIntfNode[dm].n; i++)
    {
      if (F2I[i] < 0)       prgError(1,fct,"node mismatch!");
      if (F2I[i] > numnp-1) prgError(2,fct,"node mismatch!");
      matchCnt[F2I[i]]++;
      if (matchCnt[F2I[i]] > j) j = matchCnt[F2I[i]];
    }
  }

  for (i=0; i<j+1; i++) matchCnt.append(0);
  for (i=0; i<numnp; i++) matchCnt[numnp+matchCnt[i]]++;

  COUT << "free node match count:\n\n";
  for (i=1; i<j+1; i++) { COUT; printf("%5d nodes -> match(%d)\n",matchCnt[numnp+i],i); }
  cout << "\n";

  COUT << "  number of comparisons: " << compCount << "\n\n";

  node.free();

  // generate intfNodeToFreeNode

  ListArray< Vector<int> > intfNodeToFreeNodeTmp;

  intfNodeToFreeNodeTmp.setDim(numnp);

  for (dm=0; dm<domType.n; dm++)
  {
    for (i=0; i<freeNodeToIntfNode[dm].n; i++)
    {
      if (freeNodeToIntfNode[dm][i] > -1)
      {
        intfNodeToFreeNodeTmp[freeNodeToIntfNode[dm][i]].append(dm);
        intfNodeToFreeNodeTmp[freeNodeToIntfNode[dm][i]].append(i);
      }
    }
  }
  intfNodeToFreeNode = intfNodeToFreeNodeTmp; intfNodeToFreeNodeTmp.free();

  // output and check (should be removed eventually)

  for (i=0; i<numnp; i++)
  {
    if (intfNodeToFreeNode[i].n < 1) prgError(1,fct,"no freeNodes for this interface node?!");

    if (show) 
    {
      printf("%4i: ",i);
      for (j=0; j<intfNodeToFreeNode[i].n; j+=2)
        printf("%12s,%5i;",
               domain.name(domPtr[intfNodeToFreeNode[i][j]]),intfNodeToFreeNode[i][j+1]);
      cout << "\n";
    }

    dm = intfNodeToFreeNode[i][0];

    j0 = intfNodeToFreeNode[i][1];

    n = freeNode(dm)[j0].x.n;

    for (j=2; j<intfNodeToFreeNode[i].n; j+=2)
    {
      if (n != freeNode(intfNodeToFreeNode[i][j])[intfNodeToFreeNode[i][j+1]].x.n) 
        prgError(2,fct,"mismatch of spatial dimension!");

      if (dist2(freeNode(dm)[j0].x.x,
        freeNode(intfNodeToFreeNode[i][j])[intfNodeToFreeNode[i][j+1]].x.x, n) > 10.e-10)
          prgError(2,fct,"nodes do not match!");
    }
  }
  if (show) cout << "\n";

  return;
}













void InterfaceMatch::findNodeMatches(int d)
{
  int i, j, k;

  for (i=0; i<freeNode(d).n; i++)
  {
    if (F2I[i] == -1)
    {
      j = 0; 

      while (j < node.n) 
      {
        if (!node[j].fnd)
        {
          compCount++;

          if (freeNode(d)[i].x.n == node[j].x.n)

            if (dist2(freeNode(d)[i].x.x,node[j].x.x,node[j].x.n) < 1.e-10) break; 
        }
        j++;
      }

      if (j < node.n) 
      { 
        F2I[i] = j;

        //cout << freeNode(d)[i].x << " = " << node[j].x << " a\n";

        node[j].fnd = true;

        for (j=0; j<freeNode(d)[i].toFree[0].n; j++)
        {
          k = freeNode(d)[i].toFree[0][j].nd;

          if (F2I[k] == -1)
          {
            searchFrom.append(F2I[i]);

            findNodeForMe(d,k);
          }
        }
      }
      else F2I[i] = - 2; // no match in nd exists!
    }
  }
  return;
}














void InterfaceMatch::findNodeForMe(int d, int kd)
{
  int i = 0, j, k;//, xxx = searchFrom[0], yyy, zzz;

  Vector<MoreData*> &dP = node[searchFrom[0]].depPtr;

  searchFrom.del(0);

  while (i < dP.n)
  {
    if (!(dP[i]->fnd) && !(dP[i]->chk))
    {
      compCount++;

      if (dP[i]->x.n == freeNode(d)[kd].x.n)

        if (dist2(dP[i]->x.x,freeNode(d)[kd].x.x,dP[i]->x.n) < 1.e-10) break;

      dP[i]->chk = true;

      searchFrom.append(dP[i]->c);

      changed.append(dP[i]->c);
    }
    i++;
  }

  if (i < dP.n) 
  {
    F2I[kd] = dP[i]->c; 

    //yyy = kd;
    //zzz = node[xxx].depPtr[i]->c;

    //cout << freeNode(d)[kd].x << " = " << node[node[xxx].depPtr[i]->c].x << " b1\n";

    dP[i]->fnd = true;

    for (j=0; j<changed.n; j++) node[changed[j]].chk = false;

    changed.free();

    for (j=0; j<freeNode(d)[kd].toFree[0].n; j++) 
    {
      k = freeNode(d)[kd].toFree[0][j].nd;

      if (F2I[k] == -1) 
      {
        searchFrom.free();
        searchFrom.append(dP[i]->c);
        findNodeForMe(d,k);
      }
    }
    //cout << freeNode(d)[yyy].x << " = " << node[zzz].x << " b2\n";

  }
  else if (searchFrom.n > 0) findNodeForMe(d,kd);

  return;
}












void InterfaceMatch::addNonMatchingNodes(int d)
{
  int i;

  for (i=0; i<freeNode(d).n; i++)
  {
    if (F2I[i] == -2)
    {
      node.add(new MoreData(numnp,freeNode(d)[i]));

      F2I[i] = numnp++;
    }
    else if (F2I[i] == -1) 

      prgError(1,"InterfaceMatch::addNonMatchingNodes","this should never have happened!");
  }
  return;
}












void InterfaceMatch::extendConnectivity(int d)
{
  int i, j;

  MoreData *ptr;

  for (i=0; i<freeNode(d).n; i++)
  {
    for (j=0; j<freeNode(d)[i].toFree[0].n; j++)
    {
      ptr = &(node[F2I[freeNode(d)[i].toFree[0][j].nd]]);

      if (!(node[F2I[i]].depPtr.contains(ptr))) node[F2I[i]].depPtr.append(ptr);
    }
  }
  return;
}






void InterfaceMatch::printInfo(void)
{ 
  COUT << "ndm ...... = " << ndm   << "\n";
  COUT << "ndf ...... = " << ndf   << "\n";
  COUT << "numnp .... = " << numnp << "\n";
  COUT << "nequ ..... = " << nequ  << "\n";
  COUT << "neqx ..... = " << neqx  << "\n\n";

  return;
}










void InterfaceMatch::generateIdux(void)
{
  char fct[] = "InterfaceMatch::generateIdux";
  
  //cout << fct << "\n";

  int d, i, k, l, m, nn;

  idu.setDim(numnp,ndf,true);
  idu.zero();

  idx.setDim(numnp,ndm,true);
  idx.zero();

  for (d=0; d<domType.n; d++)
  {
    for (i=0; i<freeNode(d).n; i++)
    {
      l = freeNodeToIntfNode[d][i];

      m = freeNode(d)[i].dofU->n;

      for (k=0; k<m; k++) idu.x[ndf*l+k] = 1;

      m = freeNode(d)[i].dofX->n;

      for (k=0; k<m; k++) idx.x[ndm*l+k] = 1;
    }
  }

  nn = ndf + ndm + 3;

  m = bcTmp.n / nn;

  for (i=0; i<m; i++)
  {
    d = 0; while (d < domType.n)
      if (domType[d] == bcTmp[i*nn] && domId[d] == bcTmp[i*nn+1]) break; else d++;
    if (d == domType.n) prgError(1,fct,"invalid domain specified in boundary conditions!");

    l = bcTmp[i*nn+2];

    k = 0; while (k < freeNode(d).n)
      { if (freeNode(d)[k].nd == l) break; else k++; }
    if (k == freeNode(d).n) prgError(2,fct,"invalid node specified in boundary conditions!");

    for (l=0; l<ndf; l++) idu.x[freeNodeToIntfNode[d][k]*ndf+l] = bcTmp[i*nn+3+l];

    for (l=0; l<ndm; l++) idx.x[freeNodeToIntfNode[d][k]*ndm+l] = bcTmp[i*nn+3+ndf+l];
  }

  bcTmp.free();

  return;
}












void InterfaceMatch::setLagrangianNodes(void)
{
  char fct[] = "InterfaceMatch::setLagrangianFreeNodes";

  int dm, i, j, nDOF, *DOF, *BNDND, *IDX = idx.x;

  VectorArray<int> flag, isBnd, isLagrDom;

  // determine Lagrangian / non Lagrangian nodes

  flag.setDim(numnp);

  flag.zero();

  for (dm=0; dm<domType.n; dm++)
  {
    if (isFreeSurface(*domPtr[dm]))
    {
      isBnd.setDim(domPtr[dm]->numnp);

      isBnd.zero();

      BNDND = domPtr[dm]->bndNd;

      for (i=0; i<domPtr[dm]->nBndNd; i++) isBnd[BNDND[i]-1] = 1;

      F2I = freeNodeToIntfNode[dm].x;

      for (i=0; i<freeNode(dm).n; i++)
      {
        if (isBnd[i] == 1)
        {
          if (flag[F2I[i]] == 0) flag[F2I[i]] = -1;

          else if (flag[F2I[i]] == 1) prgError(1,fct,"freeSurface inconsistency!");
        }
        else
        {
          if (flag[F2I[i]] == 0) flag[F2I[i]] = 1;

          else if (flag[F2I[i]] == -1) prgError(2,fct,"freeSurface inconsistency!");
        }
      }
    }
  } // flag is now 1 for freeNodes that are subject to FreeSurface mesh motion

  isLagr.setDim(numnp);

  for (i=0; i<numnp; i++) isLagr[i] = (flag[i] != 1);

  flag.free();

  // adjust idx if required

  for (i=0; i<numnp; i++)
  {
    j = 0; while (j < ndm) if (IDX[i*ndm+j] == 1) break; else j++;

    if (j > 0 && !isLagr[i])

      prgError(3,fct,"freeSurface / boundary condition inconsistency!");

    else if (isLagr[i]) for (j=0; j<ndm; j++) IDX[i*ndm+j] = 0;
  }

  // if non Lagrangian nodes exist, uType has to be VELOCITY

  i = 0; while (i < numnp) if (!isLagr[i]) break; else i++;

  if (i < numnp && uType != VELOCITY) 

    prgError(1,fct,"set uType = VELOCITY for problems with non Lagrangian interface nodes!");

  // set Lagrangian nodes in each domain

  for (dm=0; dm<domType.n; dm++)
  {
    F2I = freeNodeToIntfNode[dm].x;

    isLagrDom.setDim(freeNode(dm).n);

    for (i=0; i<freeNode(dm).n; i++) if (isLagr[F2I[i]]) isLagrDom[i] = 1; else isLagrDom[i] = 0;

    domPtr[dm]->setLagrangianFreeNodes(isLagrDom);
  }

  return;
}













void InterfaceMatch::initKinData(void)
{
  char fct[] = "InterfaceMatch::initKinData";

  int i, j, k, d, n;

  x.setDim(numnp,ndm,true);
  xn.setDim(numnp,ndm,true);
  x0.setDim(numnp,ndm,true);

  u.setDim(numnp,ndf,true);    u.zero();   // this will one day have to be replaced
  un.setDim(numnp,ndf,true);   un.zero();  // by getting the u values from the 
  du.setDim(numnp,ndf,true);   du.zero();  // domains
  dun.setDim(numnp,ndf,true);  dun.zero(); //
  ddu.setDim(numnp,ndf,true);  ddu.zero(); //
  ddun.setDim(numnp,ndf,true); ddun.zero();//

  for (i=0; i<numnp; i++)
  {
    d = intfNodeToFreeNode[i][0];
    n = intfNodeToFreeNode[i][1];

    for (k=0; k<freeNode(d)[n].x.n; k++) 
    {
       x.x[i*ndm+k] = freeNode(d)[n].x.x[k];
      xn.x[i*ndm+k] = freeNode(d)[n].x.x[k];
      x0.x[i*ndm+k] = freeNode(d)[n].x.x[k];
    }

    for (j=2; j<intfNodeToFreeNode[i].n; j+=2)
    {
      d = intfNodeToFreeNode[i][j];
      n = intfNodeToFreeNode[i][j+1];

      for (k=0; k<freeNode(d)[n].x.n; k++)
      {
        if (abs(freeNode(d)[n].x.x[k] - x.x[i*ndm+k]) > 1.e-8)
          prgError(1,fct,"interface coordinate mismatch!");
      }
    }
  }

  return;
}









bool InterfaceMatch::converged(void)
{
  if (rNorm < tol) return true; else return false;
}









bool InterfaceMatch::diverging(double factor)
{
  if (rNormPrev > -0.1 && rNorm / rNormPrev > factor) return true;  

  if (prgNAN(rNorm)) return true;
  
  return false;
}










void InterfaceMatch::setTimeParam(void)
{
  char fct[] = "InterfaceMatch::setTimeParam";

  double dt = mpapTime.dt, alpf, alpm, beta, gamm, rho;
	
  for (int i=10; i<TD_DIM; i++) td[i] = 0.;

  td[5-1] = dt;

  if (uType == DISPLACEMENT)
  {
    switch (tis)
    {
      case  0: // quasi static

               td[5] = 1.0;
  	   
               prgError(1,fct,"quasi static interface !?!?");
   
  	       break;

      case  1: // generalised midpoint rule

               gamm = td[0];

               td[5]  = gamm;
               td[6]  = gamm;
               td[7]  = gamm;
  	     
               td[10] = - (1. - gamm) / gamm;      // dU_n    in    dU_n+1
               td[11] = 1. / (gamm * dt);          // U_n+1

               td[8]  = td[11] * td[11];  // coefficient of U_n+1 in ddU_n+1 
               td[9]  = td[11];           // coefficient of U_n+1 in dU_n+1 
 
               dxdu = 1. / td[9];

               break;

      case  2: // generalised alpha-method

               rho  = td[0];
  	     
               alpf = 1. / (1. + rho);
               alpm = (2. - rho) / (1. + rho);
               beta = .25 * (1. + alpm - alpf)*(1. + alpm - alpf);
               gamm = .5 + alpm - alpf;
  	     
               td[5]  = alpf;
               td[6]  = alpf;
               td[7]  = alpm;
  	     
               td[10] = gamm / (beta * dt);             // U_n+1   in  dU_n+1
               td[11] = - td[10];                       // U_n
               td[12] = 1. - gamm / beta;               // dU_n    
               td[13] = dt * (1. - gamm / (2. * beta)); // ddU_n
  	     
               td[14] = 1. / (beta * dt * dt);          // U_n+1   in  ddU_n+1
               td[15] = - td[14];                       // U_n
               td[16] = - 1. / (beta * dt);             // dU_n
               td[17] = 1. - 1. / (2. * beta);          // ddU_n

               td[8]  = td[14];  // coefficient of U_n+1 in ddU_n+1 
               td[9]  = td[10];  // coefficient of U_n+1 in dU_n+1 
  	     
               dxdu = 1. / td[9];

  	       break;
  	     
      default: prgError(1,fct,"invalid value of tis!");
    }
  }
  else if (uType == VELOCITY)
  {
    switch (tis)
    {
      case  0: // quasi static

               td[6-1] = 1.0;
  	   
               prgError(1,fct,"quasi static interface !?!?");
   
  	       break;

      case  1: // generalised midpoint rule

               gamm = td[0];

               td[5]  = gamm;
               td[6]  = gamm;
               td[7]  = gamm;
  	     
               td[10] = dt * gamm;              // dU_n+1   in   U_n+1
               td[11] = dt * (1. - gamm);       // dU_n

               td[12] = 1. / (dt * gamm);       // dU_n+1  in  ddU_n+1
               td[13] = - td[12];               // dU_n
               td[14] = (gamm - 1.) / gamm;     // ddU_n

               td[8]  = td[12];        // coefficient of dU_n+1 in ddU_n+1 
               td[9]  = td[10];        // coefficient of dU_n+1 in U_n+1 
  	     
               dxdu   = td[9];

               break;

      case  2: // generalised alpha-method

               rho  = td[0];
  	     
               alpf = 1. / (1. + rho);
               alpm = (2. - rho) / (1. + rho);
               beta = .25 * (1. + alpm - alpf)*(1. + alpm - alpf);
               gamm = .5 + alpm - alpf;
  	     
               td[5]  = alpf;
               td[6]  = alpf;
               td[7]  = alpm;
  	     
               td[10] = beta * dt / gamm;               // dU_n+1   in  U_n+1
               td[11] = 1.;                             // U_n
               td[12] = dt * (1. - beta / gamm);        // dU_n    
               td[13] = dt * dt * (.5 - beta / gamm);   // ddU_n
  	     
               td[14] = 1. / (gamm * dt);               // dU_n+1   in  ddU_n+1
               td[15] = 0.;                             // U_n
               td[16] = - td[14];                       // dU_n
               td[17] = 1. - 1. / gamm;                 // ddU_n

               td[8]  = td[14];   // coefficient of dU_n+1 in ddU_n+1 
               td[9]  = td[10];   // coefficient of dU_n+1 in U_n+1 
  	     
               dxdu   = td[9];

  	       break;
  	     
      default: prgError(2,fct,"invalid value of tis!");
    }
  }
  else prgError(4,fct,"invalid uType!");

  // set dxdu for induvidual domains

  for (int dm=0; dm<domType.n; dm++) domPtr[dm]->dxdu = dxdu;

  return;
}












void InterfaceMatch::timeUpdate(void)
{
  char fct[] = "InterfaceMatch::timeUpdate";

  int i, indf = 0, indm = 0, numnpXndf = numnp * ndf, numnpXndm = numnp * ndm;

  double *X = x.x, *Xn = xn.x, *X0 = x0.x,
	 *U = u.x, *Un = un.x, *dU = du.x, 
	 *dUn = dun.x, *ddU = ddu.x, *ddUn = ddun.x,
         dt, dt2, dt3, dt4, dt5, fact[6];

  // xn <- x
 
  for (i=0; i<numnpXndm; i++) Xn[i] = X[i];

  // new initial guess
  //
  // and
  //
  // un <- u
  // dun <- du
  // ddun <- ddu

  for (i=0; i<numnpXndf; i++) 
  {
    Un[i] = U[i];
    dUn[i] = dU[i];
    ddUn[i] = ddU[i];
  }

  if      (initGuess == KEEP_VELOCITY && uType == DISPLACEMENT)
  {
    switch (tis)
    {
      case 1: for (i=0; i<numnpXndf; i++) U[i] += (1.-td[10])*dUn[i] / td[11]; break;

      case 2: for (i=0; i<numnpXndf; i++) 

                U[i] += ((1.-td[12])*dUn[i] + td[13]*ddUn[i]) / td[10]; break;
    }
  }
  else if (initGuess == KEEP_DISPLACEMENT && uType == VELOCITY)
  {
    switch (tis)
    {
      case 1: for (i=0; i<numnpXndf; i++) dU[i] = -td[11]*dUn[i] / td[10]; break;

      case 2: for (i=0; i<numnpXndf; i++) dU[i] = -(td[12]*dUn[i]+td[13]*ddUn[i])/td[10]; break;
    }
  }

  // update interface traction forces

  double fact1, fact2, *IR, *IRn;

  int dm, n, j;

  for (dm=0; dm<domType.n; dm++)
  {
    fact1 = 1. / domReacTime[dm],
    fact2 = - (1.-domReacTime[dm]) * fact1;
         
    for (i=0; i<freeNode(dm).n; i++)  
    {
      IR  = freeNode(dm)[i].reac.x;
      IRn = freeNode(dm)[i].reacn.x;
      n   = freeNode(dm)[i].reacn.n;

      for (j=0; j<n; j++) IRn[j] = fact1 * IR[j] + fact2 * IRn[j];
    }
  }

  // set iteration flag
  
  firstIter = true;
 
  // update X, dU, ddU    !important!
  
  updateIterStep();
	 
  return;   
}











void InterfaceMatch::updateIterStep(void)
{
  char fct[] = "InterfaceMatch::updateIterStep";

  int i, j, im, indf = 0, indm = 0, numnpXndf = numnp * ndf, *IDX = idx.x;

  double *X = x.x, *X0 = x0.x, *U = u.x, 
	 *Un = un.x, *dU = du.x, *dUn = dun.x, *ddU = ddu.x, *ddUn = ddun.x;

  //cout << fct << "\n";

  if (uType == DISPLACEMENT)
  {
    im = 0;

    for (i=0; i<numnp; i++)
    {
      im += ndf;

      if (isLagr[i]) // node is Lagrangian

        switch (tis)
        {
          case  0: // quasi static
    
                   // nothing to be done here
    
                   break;
    
          case  1: // generalised midpoint rule
    
                   for (j=im-ndf; j<im; j++)
                   {
                      dU[j] = td[10]* dUn[j] + td[11]*( U[j]- Un[j]);
                     ddU[j] = td[10]*ddUn[j] + td[11]*(dU[j]-dUn[j]);
                   }
      	         break;
    
          case  2: // generalised alpha-method
   
                   for (j=im-ndf; j<im; j++)
      	           {
                      dU[j] = td[10]*(U[j]-Un[j]) + td[12]*dUn[j] + td[13]*ddUn[j];
                     ddU[j] = td[14]*(U[j]-Un[j]) + td[16]*dUn[j] + td[17]*ddUn[j];
      	           }
                   break;
      	     
          default: prgError(1,fct,"invalid value of tis!");
        }
    }
  }
  else if (uType == VELOCITY)
  {
    im = 0;

    for (i=0; i<numnp; i++)
    {
      im += ndf;

      switch (tis)
      {
        case  0: // quasi static
    
    	     // nothing to be done here
    
    	     break;
    
        case  1: // generalised midpoint rule

                 if (isLagr[i]) // node is Lagrangian

                   for (j=im-ndf; j<im; j++)
                   {
                       U[j] = Un[j] + td[10]*dU[j] + td[11]*dUn[j];
    
                     ddU[j] = td[14]*ddUn[j] + td[12]*(dU[j]-dUn[j]);
                   }

                 else

                   for (j=im-ndf; j<im; j++)
    
                     ddU[j] = td[14]*ddUn[j] + td[12]*(dU[j]-dUn[j]);

    	         break;
    
        case  2: // generalised alpha-method
    
                 if (isLagr[i]) // node is Lagrangian

                   for (j=im-ndf; j<im; j++)
                   {
                       U[j] = Un[j] + td[10]*dU[j] + td[12]*dUn[j] + td[13]*ddUn[j];
  
                     ddU[j] = td[14]*(dU[j]-dUn[j]) + td[17]*ddUn[j];
                   }

                 else

                   for (j=im-ndf; j<im; j++)
  
                     ddU[j] = td[14]*(dU[j]-dUn[j]) + td[17]*ddUn[j];
    	     
    	         break;
    	     
        default: prgError(1,fct,"invalid value of tis!");
      }
    }
  }
  else prgError(1,fct,"invalid uType!");

  // update coordinates
    
  for (i=0; i<numnp; i++)
  {
    for (j=0; j<ndm; j++) {  X[indm+j] = X0[indm+j] + U[indf+j]; }

    indf += ndf;
    indm += ndm;
  }

  return;
}



