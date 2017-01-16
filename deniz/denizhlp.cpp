
#include <iostream>

#include "FunctionsProgram.h"
#include "FunctionsSupport.h"
#include "FunctionsDeniz.h"


int micro4to1(int i, int j, int k, int l, int isw)
{
  int index, m, n;


switch (isw)
{
case 1:

switch (i)
{
case 1:
	switch (j)
	{
	case 1: m=2; break;
	default: m=3; break;
	}
break;

default:
	switch (j)
	{
	case 1: m=3; break;
	default: m=1; break;
	}
break;
}

switch (k)
{
case 1:
	switch (l)
	{
	case 1: n=2; break;
	default: n=3; break;
	}
break;

default:
	switch (l)
	{
	case 1: n=3; break;
	default: n=1; break;
	}
break;
}

switch (m)
{

case 3:
	switch (n)
	{
	case 3: index=8; break;
	case 2: index=5; break;
	default: index=2; break;
	}
break;


case 2:
	switch (n)
	{
	case 3: index=7; break;
	case 2: index=4; break;
	default: index=1; break;
	}
break;

default:
	switch (n)
	{
	case 3: index=6; break;
	case 2: index=3; break;
	default: index=0; break;
	}
break;
}

break;

default:

switch (i)
{
case 1:
	switch (j)
	{
	case 1: m=2; break;
	default: m=4; break;
	}
break;

default:
	switch (j)
	{
	case 1: m=3; break;
	default: m=1; break;
	}
break;
}

switch (k)
{
case 1:
	switch (l)
	{
	case 1: n=2; break;
	default: n=4; break;
	}
break;

default:
	switch (l)
	{
	case 1: n=3; break;
	default: n=1; break;
	}
break;
}

switch (m)
{

case 4:
	switch (n)
	{
	case 4: index=15; break;
	case 3: index=11; break;
	case 2: index=7; break;
	default: index=3; break;
	}
break;


case 3:
	switch (n)
	{
	case 4: index=14; break;
	case 3: index=10; break;
	case 2: index=6; break;
	default: index=2; break;
	}
break;


case 2:
	switch (n)
	{
	case 4: index=13; break;
	case 3: index=9; break;
	case 2: index=5; break;
	default: index=1; break;
	}
break;

default:
	switch (n)
	{
	case 4: index=12; break;
	case 3: index=8; break;
	case 2: index=4; break;
	default: index=0; break;
	}
break;
}

break;
}

  return index;
}









 
 
int micro2to1(int i, int j, int isw)
{
  int index;

switch (isw)
{
case 1:

switch (i)
{
case 1:
	switch (j)
	{
	case 1: index=1; break;
	default: index=2; break;
	}
break;

default:
	switch (j)
	{
	case 1: index=2; break;
	default: index=0; break;
	}
break;
}

break;

default:

switch (i)
{
case 1:
	switch (j)
	{
	case 1: index=3; break;
	default: index=1; break;
	}
break;

default:
	switch (j)
	{
	case 1: index=2; break;
	default: index=0; break;
	}
break;
}

break;
}
 
  
  return index;
}









 
 
int integrateQuadraticBoundary 
(double *Integralx, double *Integraly, double *X, int edge, int nen){
double xi[6], shp[8], Jmtx[4], product,
	  nx, ny,
	  	  tx, ty,tsize,
		  	  nodex, nodey,nodesize,
	  d=1.0e-3, xA[3], yA[3], localtol=1.0e-6,length,
	  r1=1.e0, r2=2.e0, r0=0.e0, r06=6.e-1, wgp[3];

double sq06=sqrt(r06), sq2d2=sqrt(r2)/r2, r1d2=r1/r2;

  int l, i, nxi, Norder[3],jmtxc=0;

  for(i=0; i<3; i++) {Integralx[i]=0.0; Integraly[i]=0.0;}

wgp[0]=5.0/9.0;
wgp[1]=8.0/9.0;
wgp[2]=5.0/9.0;

  if(edge==2 || edge==4) jmtxc=2;


switch(nen)
{
		case  6:  
			prgError(1,"integrateQuadraticBoundary","not tested yet!");
		switch(edge)
		{
		case  2:
		   xi[0]= (r1-sq06)*r1d2;
		   xi[1]= r1-(r1-sq06)*r1d2;
		   xi[2]= r1d2;
		   xi[3]= r1d2;
		   xi[4]= (r1+sq06)*r1d2;
		   xi[5]= r1-(r1+sq06)*r1d2;
		   for(i=0; i<3; i++) wgp[i]*=sq2d2; 
		   Norder[0]=1; Norder[1]=4; Norder[2]=2;break;   

		case  3:
		   xi[0]= r0;
		   xi[1]= (r1-sq06)*r1d2;
		   xi[2]= r0;
		   xi[3]= r1d2;
		   xi[4]= r0;
		   xi[5]= (r1+sq06)*r1d2;
		   for(i=0; i<3; i++) wgp[i]*=r1d2;
		   Norder[0]=2; Norder[1]=5; Norder[2]=0;break; 

		default: 
		   xi[0]= (r1-sq06)*r1d2;
		   xi[1]= r0;
		   xi[2]= r1d2;
		   xi[3]= r0;
		   xi[4]= (r1+sq06)*r1d2;
		   xi[5]= r0;
		   for(i=0; i<3; i++) wgp[i]*=r1d2;
 			Norder[0]=0; Norder[1]=3; Norder[2]=1;break; 
		}
			
			
break;

		default:     

		for(i=0; i<3; i++) yA[i] = -r1;

		xA[0]= -sqrt(r06);
		xA[1]= r0;
		xA[2]= sqrt(r06);


		switch(edge)
		{
		case  2:     for(i=0; i<3; i++) {xi[2*i]=-yA[i]; xi[2*i+1]=+xA[i];} 
					 Norder[0]=1; Norder[1]=5; Norder[2]=2;break;
		case  3:     for(i=0; i<3; i++) {xi[2*i]=-xA[i]; xi[2*i+1]=-yA[i];}
					 Norder[0]=2; Norder[1]=6; Norder[2]=3;break;
		case  4:     for(i=0; i<3; i++) {xi[2*i]=+yA[i]; xi[2*i+1]=-xA[i];}
 					 Norder[0]=3; Norder[1]=7; Norder[2]=0;break;
		default:     for(i=0; i<3; i++) {xi[2*i]=+xA[i]; xi[2*i+1]=+yA[i];}
 					 Norder[0]=0; Norder[1]=4; Norder[2]=1;break;
		}
			
break;
}
    

//now calculate the integral
  for(l=0; l<3; l++)
  {

	  		switch(nen)
		{
			case  3:   compshpbt_(shp,Jmtx,&(xi[l+l]),X,&nen); break;
			default:   compshpbq_(shp,Jmtx,&(xi[l+l]),X,&nen); break;
		}

   // Choose  Jmtx[0]&[1](x) or Jmtx[2]&[3](y)
	tx=Jmtx[jmtxc];  ty=Jmtx[jmtxc+1];
	tsize=sqrt(tx*tx+ty*ty);
	tx/=tsize;  ty/=tsize;

	nodex=X[2*Norder[2]]-X[2*Norder[0]];
	nodey=X[2*Norder[2]+1]-X[2*Norder[0]+1];
	nodesize=sqrt(nodex*nodex+nodey*nodey);
	nodex/=nodesize;  nodey/=nodesize;



product=tx*nodex+ty*nodey;
if(product<localtol) {tx=-tx; ty=-ty;}

	nx=ty; ny=-tx;

	for(nxi=0; nxi<3; nxi++)
	{
	Integralx[nxi]+=nx*shp[Norder[nxi]]*wgp[l];
	Integraly[nxi]+=ny*shp[Norder[nxi]]*wgp[l];
	}
  }
  return 0;  
}
 
  
int integrateLinearBoundary 
(double *Integralx, double *Integraly, double *X, int edge, int nen){
double xi[4], shp[4], Jmtx[4],  product,
	  coorx[4], coory[4], m, b, 
	  nx, ny,
	  	  tx, ty,tsize,
		  	  nodex, nodey,nodesize,
	  xA[2], yA[2], localtol=1.0e-6, 
	  r1=1.e0, r0=0.e0, r3=3.e0,  wgp[]={1.0, 1.0}, r1d2=0.5e0;

  int l, i, nxi, Norder[2],jmtxc=0;

    for(i=0; i<2; i++) {Integralx[i]=0.0; Integraly[i]=0.0;}

  if(edge==2 || edge==4) jmtxc=2;

		xA[0]= -r1/sqrt(r3);
		xA[1]= +r1/sqrt(r3);


		switch(nen)
		{
			case  3:
			prgError(1,"integrateLinearBoundary","not tested yet!");
		xA[0]= (r1-r1/sqrt(r3))*r1d2; wgp[0]=r1d2;
		xA[1]= (r1+r1/sqrt(r3))*r1d2; wgp[1]=r1d2;



		switch(edge)
		{
		case  2:     for(i=0; i<2; i++) {xi[2*i]=xA[i]; xi[2*i+1]=r1-xA[i];} 
					 Norder[0]=1; Norder[1]=2;  break;
		case  3:     for(i=0; i<2; i++) {xi[2*i]=r0; xi[2*i+1]=xA[i];}
					 Norder[0]=2; Norder[1]=0;  break; 
		default:     for(i=0; i<2; i++) {xi[2*i]=+xA[i]; xi[2*i+1]=r0;}
 					 Norder[0]=0; Norder[1]=1; 	break;
		}
break;


			default:
				     for(i=0; i<2; i++) yA[i] = -r1;
		switch(edge)
		{
		case  2:     for(i=0; i<2; i++) {xi[2*i]=-yA[i]; xi[2*i+1]=+xA[i];} 
					 Norder[0]=1; Norder[1]=2;  break;
		case  3:     for(i=0; i<2; i++) {xi[2*i]=-xA[i]; xi[2*i+1]=-yA[i];}
					 Norder[0]=2; Norder[1]=3; break;
		case  4:     for(i=0; i<2; i++) {xi[2*i]=+yA[i]; xi[2*i+1]=-xA[i];}
 					 Norder[0]=3; Norder[1]=0;  break;
		default:     for(i=0; i<2; i++) {xi[2*i]=+xA[i]; xi[2*i+1]=+yA[i];}
 					 Norder[0]=0; Norder[1]=1; 
		break;
		}
break;
		}


//now calculate the integral
  for(l=0; l<2; l++)
  {

	  		switch(nen)
		{
			case  3:   compshpbt_(shp,Jmtx,&(xi[l+l]),X,&nen); break;
			default:   compshpbq_(shp,Jmtx,&(xi[l+l]),X,&nen); break;
		}

   // Choose  Jmtx[0]&[1](x) or Jmtx[2]&[3](y)
	tx=Jmtx[jmtxc];  ty=Jmtx[jmtxc+1];
	tsize=sqrt(tx*tx+ty*ty);
	tx/=tsize;  ty/=tsize;

	nodex=X[2*Norder[1]]-X[2*Norder[0]];
	nodey=X[2*Norder[1]+1]-X[2*Norder[0]+1];
	nodesize=sqrt(nodex*nodex+nodey*nodey);
	nodex/=nodesize;  nodey/=nodesize;



product=tx*nodex+ty*nodey;
if(product<localtol) {tx=-tx; ty=-ty;}

	nx=ty; ny=-tx;

	for(nxi=0; nxi<2; nxi++)
	{
	Integralx[nxi]+=nx*shp[Norder[nxi]]*wgp[l];
	Integraly[nxi]+=ny*shp[Norder[nxi]]*wgp[l];
	}
  }
  return 0;  
}




