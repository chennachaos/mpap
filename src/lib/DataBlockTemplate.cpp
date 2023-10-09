
#include <iostream>

#include "DataBlockTemplate.h"
#include "FunctionsProgram.h"

using namespace std;


extern bool readIfileInfo;


DataBlockTemplate::DataBlockTemplate(void)
{
  list = NULL;
	
  return;
}





DataBlockTemplate::~DataBlockTemplate()
{
  free();
	
  return;
}





bool DataBlockTemplate::initialise(char *tmplStr, char **l1, char **l2, char **l3, char **l4)
{
  MyString tmp; tmp = tmplStr;

  return initialise(tmp,l1,l2,l3,l4);
}




bool DataBlockTemplate::initialise(MyString &tmplStr, char **l1, char **l2, char **l3, char **l4)
{
  //cout << tmplStr << "\n";

  char **tmpl, fct[] = "DataBlockTemplate::DataBlockTemplate";
 
  int i, j, l, n = 0, cl, nl = 0;

  bool count = false;
 
  free();

  if (l1 != NULL) { nl++; 
    if (l2 != NULL) { nl++; 
      if (l3 != NULL) { nl++; 
        if (l4 != NULL) { nl++; }}}}

  if (nl > 0) { list = new char** [nl]; list[0] = l1; 
    if (nl > 1) { list[1] = l2; 
      if (nl > 2) { list[2] = l3; 
        if (nl > 3) { list[3] = l4; }}}}
  
  int nw = tmplStr.split(&tmpl);

  if (nw == 0) return true;
  
  // get type sequence

  cl = 0;
 
  for (i=0; i<nw; i++)
  {
    if (tmpl[i][0] == 'i') 
    { if (n>0 && t[n-1]==0) c[n-1]++; else { c.append(1); t.append(0); n++; } }
    
      else if (tmpl[i][0] == 'f') 
      { if (n>0 && t[n-1]==1) c[n-1]++; else { c.append(1); t.append(1); n++; } }

        else if (tmpl[i][0] == 's')
	{ if (n>0 && t[n-1]==2) c[n-1]++; else { c.append(1); t.append(2); n++; } }
      
          else if (tmpl[i][0] == 'l')
          { if (n>0 && t[n-1]==3) c[n-1]++; else { c.append(1); t.append(3); n++; } cl++; }
      
            else if (tmpl[i][strlen(tmpl[i])-1] == 'i') 
            { l = strlen(tmpl[i]);
              tmpl[i][l-1] = '\0'; 
              if (!scanInt(tmpl[i],&j)) prgError(1,fct,"template!"); 
              tmpl[i][l-1] = 'i'; 
              if (n>0 && t[n-1]==0) c[n-1]+=j; else { c.append(j); t.append(0); n++; }
            }
              else if (tmpl[i][strlen(tmpl[i])-1] == 'f') 
              { l = strlen(tmpl[i]);
                tmpl[i][l-1] = '\0'; 
                if (!scanInt(tmpl[i],&j)) prgError(2,fct,"template!"); 
                tmpl[i][l-1] = 'f'; 
                if (n>0 && t[n-1]==1) c[n-1]+=j; else { c.append(j); t.append(1); n++; }
              }
                else if (tmpl[i][strlen(tmpl[i])-1] == 's') 
                { l = strlen(tmpl[i]);
                  tmpl[i][l-1] = '\0'; 
                  if (!scanInt(tmpl[i],&j)) prgError(3,fct,"template!"); 
                  tmpl[i][l-1] = 's'; 
	          if (n>0 && t[n-1]==2) c[n-1]+=j; else { c.append(j); t.append(2); n++; }
                }
                  else if (tmpl[i][strlen(tmpl[i])-1] == 'l') 
                  { l = strlen(tmpl[i]);
                    tmpl[i][l-1] = '\0'; 
                    if (!scanInt(tmpl[i],&j)) prgError(4,fct,"template!"); 
                    tmpl[i][l-1] = 'l'; 
                    if (n>0 && t[n-1]==3) c[n-1]+=j; else { c.append(j); t.append(3); n++; }
                    cl++;
                  }
  }

  if (cl > nl) prgError(1,fct,"specify more lists!");

  // get details
  
  cl = 0;

  for (i=0; i<nw; i++)
  {
    if (tmpl[i][0] == 'i') 
    { type.append(0); 
      if (!scanInt(&(tmpl[i][1]),iDef.append())) prgError(4,fct,"template!"); 
      expr.append(1 - iDef.n);
    }
      else if (tmpl[i][0] == 'f') 
      { type.append(1); 
        if (!scanDbl(&(tmpl[i][1]),dDef.append())) prgError(5,fct,"template!"); 
        expr.append(1 - dDef.n);
      }
        else if (tmpl[i][0] == 's') 
        { type.append(2); 
          sDef.addNew(&(tmpl[i][1])); 
          expr.append(1 - sDef.n); 
        }
          else if (tmpl[i][0] == 'l') 
          { type.append(3); 
            lDef.append(which(&(tmpl[i][1]),list[cl])); cl++;
            expr.append(1 - lDef.n);
          }
            else if (tmpl[i][strlen(tmpl[i])-1] == 'i') 
            { type.append(0);
              l = strlen(tmpl[i]);
              tmpl[i][l-1] = '\0'; 
              if (!scanInt(tmpl[i],expr.append())) prgError(6,fct,"template!"); 
              tmpl[i][l-1] = 'i';
            }
              else if (tmpl[i][strlen(tmpl[i])-1] == 'f') 
              { type.append(1);
    	        l = strlen(tmpl[i]);
                tmpl[i][l-1] = '\0'; 
                if (!scanInt(tmpl[i],expr.append())) prgError(7,fct,"template!"); 
                tmpl[i][l-1] = 'f';
              }
                else if (tmpl[i][strlen(tmpl[i])-1] == 's') 
                { type.append(2);
                  l = strlen(tmpl[i]);
                  tmpl[i][l-1] = '\0'; 
                  if (!scanInt(tmpl[i],expr.append())) prgError(8,fct,"template!"); 
                  tmpl[i][l-1] = 's';
                }
                  else if (tmpl[i][strlen(tmpl[i])-1] == 'l')
                  { type.append(3);
                    l = strlen(tmpl[i]);
                    tmpl[i][l-1] = '\0'; 
                    if (!scanInt(tmpl[i],expr.append())) prgError(8,fct,"template!"); 
                    tmpl[i][l-1] = 'l';
                    cl++;
                  }
                    else if (strcmp(tmpl[i],"123")==0) // counter
                    { if (count) prgError(9,fct,"template!"); count = true;
            	      type.append(4);
                      expr.append(0);
                    }
                      else if (strcmp(tmpl[i],"012")==0) // counter
                      { if (count) prgError(10,fct,"template!"); count = true;
            	        type.append(5);
                        expr.append(0);
                      }
  }

  i = 0; if (tmpl != NULL) { while (tmpl[i] != NULL) delete [] tmpl[i++]; delete [] tmpl; }
 
  if (readIfileInfo)
  {
    cout << c << "\n";
    cout << t << "\n\n";

    cout << expr << "\n";
    cout << type << "\n\n";

    cout << iDef << "\n"; 
    cout << dDef << "\n"; 
    cout << sDef << "\n"; 
    cout << lDef << "\n\n"; 
  }
  
  return true;
}





bool DataBlockTemplate::expandToMatch(DataBlockTemplate &t2)
{
  int n, i, j;
	
  char fct[] = "DataBlockTemplate::expandToMatch";
 
  if (c.dim() > t2.c.dim())  prgError(1,fct,"templates inconsistent!");

  for (i=0; i<c.dim()-1; i++) 
    if (c[i] != t2.c[i] || t[i] != t2.t[i]) prgError(2,fct,"templates inconsistent!");

  if (c[i] != t2.c[i]) 
  {
    if (t[i] != t2.t[i]) prgError(3,fct,"templates inconsistent!");
    if (c[i] >  t2.c[i]) prgError(4,fct,"templates inconsistent!");
    
    n = t2.c[i]-c[i];
    if (n > 0) append(n, t[i]);
  }
    
  for (i=c.dim(); i<t2.c.dim(); i++)  append(t2.c[i],t2.t[i]);
  
  if (readIfileInfo)
  { 
    cout << "  after expanding : \n\n";

    cout << c << "\n";
    cout << t << "\n\n";

    cout << expr << "\n";
    cout << type << "\n\n";

    cout << iDef << "\n"; 
    cout << dDef << "\n"; 
    cout << sDef << "\n"; 
    cout << lDef << "\n\n"; 
  }
  
  return true;
}





void DataBlockTemplate::append(int n, int tp)
{
  for (int i=0; i<n; i++)
  {
    type.append(tp);
    
    switch (tp)
    {
      case 0: expr.append(-iDef.n);
	      iDef.append(0);
	      break;

      case 1: expr.append(-dDef.n);
	      dDef.append(0.);
	      break;

      case 2: expr.append(-sDef.n);
	      sDef.addNew(""); sDef[sDef.n-1].free();
	      break;

      case 3: expr.append(-lDef.n);
	      lDef.append(-1);
	      break;
    }
  }
	
  return;
}





bool DataBlockTemplate::readLine(ifstream &Ifile, MyString &line, 
		                 myVector<int> &iTmp, 
                                 myVector<double> &dTmp, 
                                 MyStringList &sTmp, 
                                 myVector<int> &lTmp)
{
  int i, j, l, cw, nw, cl = 0;

  char **word;
 
  bool fail = false;
  
  if (!Ifile) return false;
  
  line.getNextLine(Ifile);
    
  if (readIfileInfo) cout << line << "\n";
    
  nw = line.split(&word);   cw = 0;
   
  for (i=0; i<expr.n; i++)
  {    
    switch (type[i])
    {
      case 0: // integer

              if (expr[i] < 1) { iTmp.append(iDef[-expr[i]]); break; }
	
	      for (j=0; j<expr[i]; j++)
	      {	
                if (cw < nw)     
		{ if (!scanInt(word[cw++],iTmp.append(),false)) 
		  { fail = true; if (readIfileInfo) cout << "\nStop A\n\n"; break; }
		}
		else iTmp.append(0); 	
              }
              break;

      case 1: // double

	      if (expr[i] < 1) { dTmp.append(dDef[-expr[i]]); break; }
	
	      for (j=0; j<expr[i]; j++)
	      {	
                if (cw < nw)     
	        { if (!scanDbl(word[cw++],dTmp.append(),false)) 
		  { fail = true; if (readIfileInfo) cout << "\nStop A\n\n"; break; }
		}
		else dTmp.append(0.);
	      }
              break;

      case 2: // MyString

	      if (expr[i] < 1) { sTmp.addNew(sDef[-expr[i]].asCharArray()); break; }
	
              for (j=0; j<expr[i]; j++)
	      {
	        if (cw < nw)  sTmp.addNew(word[cw++]); 
		  
	          else { sTmp.addNew(""); sTmp[sTmp.n-1].free(); }
	      }
	      break;

      case 3: // which string from list

              if (expr[i] < 1) { lTmp.append(lDef[-expr[i]]); cl++; break; }

              j = 0; 
              fail = false;
              while (!fail && j < expr[i])
	      {
	        if (cw < nw)  
                {
                  l = which(word[cw++],list[cl++]); 
                  if (l < 0) fail = true; else lTmp.append(l);
                }
                else fail = true;
                j++;
	      }

              break;

      case 4: // counter 123

  	      if (cw >= nw) 
  	        { fail = true; if (readIfileInfo) cout << "\nStop C\n\n"; break; }
  	      if (!scanInt(word[cw++],&j,false))
  	        { fail = true; if (readIfileInfo) cout << "\nStop D\n\n"; break; }
  	      cnt++; if (j != cnt) 
  	        { fail = true; if (readIfileInfo) cout << "\nStop E\n\n"; } 
	      break; 
  	    
      case 5: // counter 012

  	      if (cw >= nw) 
  	        { fail = true; if (readIfileInfo) cout << "\nStop C\n\n"; break; }
  	      if (!scanInt(word[cw++],&j,false))
  	        { fail = true; if (readIfileInfo) cout << "\nStop D\n\n"; break; }
  	      cnt++; if (j != cnt-1) 
  	        { fail = true; if (readIfileInfo) cout << "\nStop E\n\n"; } 
	      break; 
  	    
      default: fail = true; if (readIfileInfo) cout << "\nStop F\n\n"; break;

    }
    if (fail) break;
  }
		 
  if (word != NULL) 
    { j = 0; while (word[j] != NULL) delete [] word[j++]; delete [] word; word = NULL; }

  if (i<expr.n) return false;
  
  return true;
}




bool DataBlockTemplate::readBlock(ifstream &Ifile, MyString &line, 
		                  myVector<int> &iTmp, 
                                  myVector<double> &fTmp, 
                                  MyStringList &sTmp, 
                                  myVector<int> &lTmp)
{
  int i, j,
      iDim, fDim, sDim, lDim,
      iDimOld = iTmp.n,
      fDimOld = fTmp.n,
      sDimOld = sTmp.n,
      lDimOld = lTmp.n,
      ni = 0, nf = 0, ns = 0, nl = 0;

  for (i=0; i<expr.n; i++) 
    switch (type[i])
    {
      case 0: if (expr[i]<1) ni++; else ni+=expr[i]; break;
      case 1: if (expr[i]<1) nf++; else nf+=expr[i]; break;
      case 2: if (expr[i]<1) ns++; else ns+=expr[i]; break;
      case 3: if (expr[i]<1) nl++; else nl+=expr[i]; break;
    }

  cnt = 0;

  while (Ifile && readLine(Ifile,line,iTmp,fTmp,sTmp,lTmp));
  
  iDim = iTmp.n - iDimOld;
  fDim = fTmp.n - fDimOld;
  sDim = sTmp.n - sDimOld;
  lDim = lTmp.n - lDimOld;

  i = intDiv(iDim,ni);
  j = intDiv(fDim,nf); if (i==0 || (j<i && j!=0)) i = j;
  j = intDiv(sDim,ns); if (i==0 || (j<i && j!=0)) i = j;
  j = intDiv(lDim,nl); if (i==0 || (j<i && j!=0)) i = j;
  
  iTmp.trunc(iDimOld+ni*i);
  fTmp.trunc(fDimOld+nf*i);
  sTmp.trunc(sDimOld+ns*i);
  lTmp.trunc(lDimOld+nl*i);

  if (iTmp.dim() <= iDimOld && ni>0) return false;
  if (fTmp.dim() <= fDimOld && nf>0) return false;
  if (sTmp.n     <= sDimOld && ns>0) return false;
  if (lTmp.n     <= lDimOld && nl>0) return false;
	  
  if (readIfileInfo) cout << iTmp << "\n\n";
  if (readIfileInfo) cout << fTmp << "\n\n";
  if (readIfileInfo) cout << sTmp << "\n\n";
  if (readIfileInfo) cout << lTmp << "\n\n";

  return true;
}




void DataBlockTemplate::free(void)
{
  iDef.free();
  dDef.free();
  sDef.free();
  lDef.free();

  c.free();
  t.free();

  expr.free();
  type.free();

  if (list != NULL) delete [] list;

  return;
}


