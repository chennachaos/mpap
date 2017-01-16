
#include "NurbsShapeFunctions.h"

bool checkClosedness(int ndim, int dir, int ngbf1, int ngbf2, int ngbf3, int* ctrlpoints)
{
    int  ii, ind1, ind2;
    bool  flag = false;
    if(ndim == 1)
    {
      cout << " 'checkClosedness' ... Not yet implemented for Curves " << endl;
    }
    else if(ndim == 2)
    {
       if(dir == 1)
       {
          ind1 = ngbf2*(ngbf1-1);
          for(ii=0;ii<ngbf2;ii++)
          {
             //cout << ctrlpoints[ii] << '\t' << ctrlpoints[ind1+ii] << endl;
             if( ctrlpoints[ii] == ctrlpoints[ind1+ii] )
               flag = true;
             else
             {
                flag = false;
                break;
             }
          }
       }
       else
       {
          for(ii=0;ii<ngbf1;ii++)
          {
             //cout << ctrlpoints[ngbf1*ii] << '\t' << ctrlpoints[ngbf1*(ii+1)-1] << endl;
             if( ctrlpoints[ngbf1*ii] == ctrlpoints[ngbf1*(ii+1)-1] )
               flag = true;
             else
             {
                flag = false;
                break;
             }
          }
       }

    }
    else
    {  
       int jj, nn;
       
       nn = ngbf1*ngbf2;
       if(dir == 1)
       {          
          for(jj=0;jj<ngbf3;jj++)
          {
             ind1 = nn*jj;
             for(ii=0;ii<ngbf2;ii++)
             {
                 if( ctrlpoints[ind1+ngbf1*ii] == ctrlpoints[ind1+ngbf1*(ii+1)-1] )
                   flag = true;
                 else
                 {
                    flag = false;
                    break;
                 }
             }
          }
       }
       else if(dir == 2)
       {
       }
       else
       {
       }
    }


   return flag;
}

