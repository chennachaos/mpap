
#include <iostream>

#include "IsogeometricFEM.h"
#include "FunctionsProgram.h"
#include "DataBlockTemplate.h"
#include "PropertyTypeEnum.h"
#include "MathGeom.h"
#include "SolverMA41.h"
//#include "SolverPARDISO.h"
#include "NurbsShapeFns.h"
#include "NurbsMiscFunctions.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "Plot.h"
#include <assert.h>
#include "UnixGUI.h"
#include "PlotVTK.h"

#include <fstream>
#include <string.h>
#include <iomanip>
#include <iostream>
#include <fstream>

#include "util.h"

using namespace std;


#include <Xm/PushB.h>
#include <Xm/Form.h>

extern Plot plot;
extern PlotVTK  plotvtk;
extern ComputerTime computerTime;
extern MpapTime     mpapTime;
extern UnixGUI    unixGUI;




void IsogeometricFEM::checkInputData1()
{
  char fct[] = "IsogeometricFEM::checkInputData1";

  int i, j;
  int index=0;

  for (i=0; i<patch.n; i++)
  {

    if( patch[i][0] > patchgrpdata.n)
    {
       prgError(1,fct,"Patches and Patch group data Inconsistenct!");
    }

    index  = patch[i][0] - 1;

    if(patch[i][1] == 1)
    {
      if((patchgrpdata[index][1] > kv.n) )
        prgError(2,fct,"Knot Vector index in patch group data is out of range!");

      DEGREE deg = patchgrpdata[index][0];
      int m = kv[patchgrpdata[index][1]-1].n - 1;
      int nn = m - deg -1;
      if((patch[i].n - 3) != (nn+1))
      {
        prgError(3,fct,"Curve Degree, Knot Vector and # of Control Points in Patches not consistenct!");
      }
    }
    else if(patch[i][1] == 2)
    {
      if((patchgrpdata[index][2] > kv.n) || (patchgrpdata[index][3] > kv.n) )
        prgError(4,fct,"Knot Vector index in patch group data is out of range!");

      DEGREE deg1 = patchgrpdata[index][0];
      DEGREE deg2 = patchgrpdata[index][1];
      int m1 = kv[patchgrpdata[index][2]-1].n - 1;
      int m2 = kv[patchgrpdata[index][3]-1].n - 1;
      int nn1 = m1 - deg1 -1;
      int nn2 = m2 - deg2 -1;

      if((patch[i].n - 3) != ((nn1+1)*(nn2+1)))
      {
        prgError(5,fct,"Surface Degrees, Knot Vector and # of Control Points in Patches not consistenct!");
      }
    }
    else if(patch[i][1] == 3)
    {
      //cout << patchgrpdata[index] << endl;
      if((patchgrpdata[index][3] > kv.n) || (patchgrpdata[index][4] > kv.n) || (patchgrpdata[index][5] > kv.n) )
        prgError(6,fct,"Knot Vector index in patch group data is out of range!");

      DEGREE deg1 = patchgrpdata[index][0];
      DEGREE deg2 = patchgrpdata[index][1];
      DEGREE deg3 = patchgrpdata[index][2];
      int m1 = kv[patchgrpdata[index][3]-1].n - 1;
      int m2 = kv[patchgrpdata[index][4]-1].n - 1;
      int m3 = kv[patchgrpdata[index][5]-1].n - 1;
      int nn1 = m1 - deg1 -1;
      int nn2 = m2 - deg2 -1;
      int nn3 = m3 - deg3 -1;
      
      int  nn = ((nn1+1)*(nn2+1)*(nn3+1));
      
      //cout << patch[i].n << endl;
      //cout << patch[i] << endl;
      
      if((patch[i].n - 3) != nn)
      {
        prgError(7,fct,"Solid Degrees, Knot Vector and # of Control Points in Patches not consistenct!");
      }
    }

  }


  return;
}




void IsogeometricFEM::checkInputData2()
{
//
  char fct[] = " IsogeometricFEM::checkInputData2 ";

  int i, j;
  
  cout << " patchElemProp[i].id  = " << patchElemProp[i].id << endl;
  
  for (i=0; i<patchElemProp.n; i++)
  {
    // check d compatibility of element type and ndf specified
    bool flag=0;
    if( ndf == 1)
    {
       switch (patchElemProp[i].id)
       {
          case  0:               //  NurbsElem1DAdvectionDiffusion
          case  1: flag = true;  //  NurbsElem1DElasticBar
                   break;
          case  3: flag = true;  // NurbsElem1DElasticBarLSFEM
                   break;
          case  13:  flag = true;  //  NurbsElem2DAdvectionDiffusion
                   break;

          case  19:  flag = true;  //  NurbsElem2DHeatTransfer
                   break;

          default: flag = false;
       }
    }
    else if( ndf == 2)
    {
       switch (patchElemProp[i].id)
       {
          case  2:               // NurbsElem1DEulerBeam
          case  4: flag = true;  // NurbsElem2DStructSolid
                   break;
          case  5: flag = true;  // NurbsElem2DStructBbarSolid
                   break;
          case  6: flag = true;  // NurbsElem2DStructFbarSolid
                   break;

          case  7: flag = true;  // NurbsElem2DStructMixed2field

                   mixedSolverFlag = 7;
                   break;

          case  8: flag = true;  // NurbsElem2DStructMixed3field

                  mixedSolverFlag = 8;
                  break;

          case  17: flag = true;  // NurbsElem2DStructSolidLSFEM2dof

                  break;

          default: flag = false;
                   break;
       }
    }
    else if( ndf == 3)
    {
       switch (patchElemProp[i].id)
       {
          case  9:  flag = true;  // NurbsElemKirchhoffPlate

                   break;

          case  10:  flag = true;  // NurbsElemMindlinPlate

                   break;

          case  11: flag = true;  // NurbsElem3DStructSolid

                   break;

          case  12: flag = true;  // NurbsElem3DStructMixed2field;

                   mixedSolverFlag = 12;
                   break;

          case  14: flag = true;  // NurbsElem2DStokes;

                   break;

          case  15: flag = true;  // NurbsElem2DNavierStokes3dof;

                   break;

          case  18: flag = true;  // NurbsElem2DStructSolidLSFEM3dof;

                   break;

          case  21: flag = true;  // NurbsElem2DStructMixed2fieldStabilised

                  break;

          default: flag = false;
                   break;
       }
    }
    else if( ndf == 4)
    {
       switch (patchElemProp[i].id)
       {
          case  16: flag = true;  // NurbsElem2DNavierStokes4dof;

                   break;

          case  20: flag = true;  // NurbsElem2DTempCoupled4dof

                   break;

          case  22: flag = true;  // NurbsElem3DStructMixed2fieldStabilised

                   break;

          default: flag = false;
                   break;
       }      
    }

      if(!flag)
         prgError(1,fct, " Element Type is not compatible with the specified DOF ");
  }
//
  return;
}






void IsogeometricFEM::plotGeom(int val1, bool flag2, int col, bool PLOT_KNOT_LINES, int* resln)
{
  if(!plotvtk.ActiveFlag && patchGrp[0].ndom == 3)
  {
     plotvtk.set();
  }
  else
    plot.setColour(col);

  char fct[] = "IsogeometricFEM::plotgeometry";

  cout << "     ISOGEOMETRICFEM: plotgeometry ...\n\n";
  
  int i;
  
  if(val1 == 1) // initial geometry
  {
    if(flag2) // plot elements
    {
      for(i=0;i<patchGrp.n;i++)
         NurbsBaseOriginal[i]->PlotElements(col, PLOT_KNOT_LINES, resln);
    }
    else // plot control points
    {
      for(i=0;i<patchGrp.n;i++)
         NurbsBaseOriginal[i]->PlotControlPoints(col);
    }
  }
  else if(val1 == 2) // final geometry
  {
    if(flag2) // plot elements
    {
      for(i=0;i<patchGrp.n;i++)
         NurbsBaseFinal[i]->PlotElements(col, PLOT_KNOT_LINES, resln);
    }
    else // plot control points
    {
      for(i=0;i<patchGrp.n;i++)
         NurbsBaseFinal[i]->PlotControlPoints(col);
    }
  }
  else// if(val1 == 3) // result
  {
    if(flag2) // plot elements
    {
      for(i=0;i<patchGrp.n;i++)
         NurbsBaseResult[i]->PlotElements(col, PLOT_KNOT_LINES, resln);
    }
    else // plot control points
    {
      for(i=0;i<patchGrp.n;i++)
         NurbsBaseResult[i]->PlotControlPoints(col);
    }
  }

  if(patchGrp[0].ndom == 3)
  {
     plotvtk.rendr->ResetCamera();
     plotvtk.renWindow->Render();
  }

  return ;
}




void IsogeometricFEM::plotSolverMatrixPattern(char *fileName)
{
  if (solver != NULL) ((SolverSparse*)solver)->mtx.plotPattern(fileName);

  return;
}





void IsogeometricFEM::prepareInteractions(void)
{
  // go and inherit from ancestors

  Domain::prepareInteractions();

  cout << "     ISOGEOMETRICFEM: preparing interactions ...\n\n";

  char fct[] = "IsogeometricFEM::prepareInteractions";


  return;
}





void IsogeometricFEM::printInfo(void)
{
   printf("\t IsogeometricFEM Analysis Model Statistics\n ");
   printf("\t ---------------------------------------------\n ");
   printf("\n");
   printf("\t Problem dimensionality        : %10dD\n ", ndm);
   printf("\n");
   printf("\t No. of Patches                : %10d\n ", Npatch);
   printf("\n");
   printf("\t No. of Global Basis Functions : %10d\n ", ntotgbf);
   printf("\n");
   printf("\t Total number of elements      : %10d\n ", totnumel);
   printf("\n");
   printf("\t No. of DOF  per CP            : %10d\n ", ndf);
   printf("\n");
   printf("\t No. of Equations              : %10d\n ", ntoteqs);
   printf("\n");
   printf("\t --------------------------------------------- ");
   printf("\n");
   printf("\n");

  return;
}








bool IsogeometricFEM::converged()
{
  char fct[] = "IsogeometricFEM::converged";

  //prgWarning(1,fct,"not yet implemented!");

  if (rNorm < tol && localStiffnessError == 0)     return true;

  return false;
}






bool IsogeometricFEM::diverging(double factor)
{
  if (rNormPrev > -0.1 && (rNorm / rNormPrev) > factor) return true;

  if (localStiffnessError != 0) return true;

  if (prgNAN(rNorm)) return true;

  return false;
}





void IsogeometricFEM::setTimeParam()
{
  for(int iii=0;iii<Npatch;iii++)
    NurbsBaseResult[iii]->setTimeParam();

  return;
}





void IsogeometricFEM::timeUpdate()
{
  cout << " IsogeometricFEM::timeUpdate() ... STARTED " << endl;

  firstIter = true;
  localStiffnessError = 0;
  filecount++;

  for(int iii=0;iii<Npatch;iii++)
    NurbsBaseResult[iii]->geomToVector(&(solnPrev[0]));

  cout << " aaaaaaaaaaaaa  " << endl;

  int  ii, jj;

    for(ii=0; ii<NurbsBaseResult[0]->ndof; ii++)
    {
      for(jj=0; jj<NurbsBaseResult[0]->ngbf; jj++)
      {
        NurbsBaseResult[0]->ValuesPrev[ii][jj]        =  NurbsBaseResult[0]->Values[ii][jj]       ;
        NurbsBaseResult[0]->ValuesDotPrev[ii][jj]     =  NurbsBaseResult[0]->ValuesDot[ii][jj]    ;
        NurbsBaseResult[0]->ValuesDotDotPrev[ii][jj]  =  NurbsBaseResult[0]->ValuesDotDot[ii][jj] ;
      }
    }

  cout << " bbbbbbbbbbbbbb  " << endl;

   if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
   {
      for(int iii=0;iii<Npatch;iii++)
        secVarPrev[iii]  =  NurbsBaseSecondVar[iii]->Values[0] ;
   }
   if(mixedSolverFlag == 8)
   {
      for(int iii=0;iii<Npatch;iii++)
      {
        secVarPrev[iii]   =  NurbsBaseSecondVar[iii]->Values[0];
        secVarPrev2[iii]  =  NurbsBaseSecondVar[iii]->Values[1];
      }
   }

   if(patchElemProp[0].id >= 14)
   {
      for(int ii=0;ii<ndf;ii++)
        Values[ii] = NurbsBaseResult[0]->Values[ii] ;
   }
    // update internal variables

  cout << " aaaaaaaaaaaaa  " << endl;

    if(intVarFlag)
    {
       double *intVar1, *intVar2;
       int i, e, nivEl;

       nivEl = elem[0]->nivGP * elem[0]->nGP;

       for (e=0; e<totnumel; e++)
       {
          intVar1 = elem[e]->intVar1;
          intVar2 = elem[e]->intVar2;

          for (i=0; i<nivEl; i++)
            intVar1[i] = intVar2[i];
       }

       intVar1 = NULL;
       intVar2 = NULL;
    }

   updateIterStep();

    cout << " IsogeometricFEM::timeUpdate() ... FINISHED " << endl;

  return;
}






void IsogeometricFEM::updateIterStep()
{
  //cout << " IsogeometricFEM::updateIterStep() ... STARTED " << endl;
  
  char fct[] = "IsogeometricFEM::updateIterStep";

  // update the coordinates
  // SurfaceResult[0].updateCoordinates(0, uinter);

  for(int iii=0;iii<Npatch;iii++)
    NurbsBaseResult[iii]->updateIterStep();

  //cout << " IsogeometricFEM::updateIterStep() ... FINISHED " << endl;
  
  return;
}




void IsogeometricFEM::reset()
{
//    cout << " Uprev " << endl;        printVector(&(Uprev[0]), Uprev.n);

    for(int iii=0;iii<patch.n;iii++)
       NurbsBaseResult[iii]->resetGeometry(&(solnPrev[0]));

    if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
    {
       for(int iii=0;iii<Npatch;iii++)
          NurbsBaseSecondVar[iii]->Values[0] = secVarPrev[iii];
    }
    if(mixedSolverFlag == 8)
    {
       for(int iii=0;iii<Npatch;iii++)
       {
          NurbsBaseSecondVar[iii]->Values[0] = secVarPrev[iii];
          NurbsBaseSecondVar[iii]->Values[1] = secVarPrev2[iii];
       }
    }

   for(int kk=0;kk<ndf;kk++)
      NurbsBaseResult[0]->Values[kk] = Values[kk];

    if(intVarFlag)
      copyElemInternalVariables();

  return;
}





void IsogeometricFEM::printData(int index, int patchnum)
{
//   char fct[] = "IsogeometricFEM::printData";

     int ii;

     if(index == 1) // print disp vector
     {
         printf("\n\n");
         printf("     Full solution  Vector ...:  \n");
         for(ii=0;ii<solnFull.n;ii++)
            printf("\t%5d\t%14.8f\n", ii, solnFull[ii]);
         printf("\n\n");
     }

     if(index == 2) // print disp vector
     {
         printf("\n\n");
         printf("     Initial solution  Vector ...:  \n");
         for(ii=0;ii<solnInit.n;ii++)
            printf("\t%5d\t%14.8f\n", ii, solnInit[ii]);
         printf("\n\n");
     }

     if(index == 3) // print residue vector
     {
         printf("\n\n");
         printf("      rhsVec and ForceVec Vectors ...:  \n");
         for(ii=0;ii<ForceVec.n;ii++)
            printf("\t%5d\t%14.12f\t%14.12f\n", ii, solver->rhsVec[ii], ForceVec[ii]);
         printf("\n\n");
     }

     if(index == 4) // global stiffness matrix
     {
//              int n1, n2, ind1, ind2, ii, jj, iii, kk;

/*
             ofstream  fout1("globalstiff.dat");

	      if(fout1.fail())
	      {
		 cout << " Could not open the Output file" << endl;
		 exit(1);
	      }

	      fout1.setf(ios::fixed);
	      fout1.setf(ios::showpoint);
	      fout1.precision(10);

             //ostream  os(fout5)

             ((SolverSparse*)solver)->mtx.print(fout1);
             
             fout1.close();
*/

             cout << "      global stiffness matrix "  << endl;
             cout << endl;

             cout << endl;

             ((SolverSparse*)solver)->mtx.print();
             cout << endl;
             cout << endl;
     }

     if(index == 5) // connectivity arrays
     {
          for(ii=0;ii<Npatch;ii++)
          {
              printf("      Connectivity Arrays for patch ...# %5d \n\n", (ii+1));
              NurbsBaseFinal[ii]->printConnectivityArrays();
          }
         printf("\n\n");
     }

     if(index == 6) // print reac vector
     {
         printf("\n\n");
         printf("      reac Vector ...:  \n");

         for(ii=0;ii<reac.n;ii++)
            printf("\t%5d\t%14.8f\n", ii, reac[ii]);
         printf("\n\n");
     }

     if(index == 7) // print knot vectors
     {
         printf("\n\n");
         for(ii=0;ii<Npatch;ii++)
         {
            cout << "      Knot Vectors for patch #" << (ii+1) << endl;
            cout << endl;
            cout << '\t' << SurfaceListFinal[ii].U << endl;
            cout << endl;
            cout << endl;
            cout << '\t' << SurfaceListFinal[ii].V << endl;
            cout << endl;
            cout << endl;
         }
     }
//
     if(index == 8) // print oup value for contour plots
     {
         int n1, n2, ind1, ind2, jj, iii, kk;

         ofstream fout("contour-data.dat");

         if(fout.fail())
         {
           cout << " Could not open the Output file" << endl;
           exit(1);
         }

         fout.setf(ios::fixed);
         fout.setf(ios::showpoint);
         fout.precision(10);

         n1 = uu[0].n;
         n2 = vv[0].n;

         for(iii=1;iii<Npatch;iii++)
           n2 += (vv[iii].n);

         fout << n1 << '\t' << n2 << '\t' << 0.0 << endl;
         if(defUndefFlag)
         {
           for(iii=0;iii<Npatch;iii++)
           {
             n1 = uu[iii].n;
             n2 = vv[iii].n;
             
             if(iii == 0)
               kk = 0;
             else
               kk = 1;

             for(jj=0;jj<n2;jj++)
             {
                ind1 = jj * n1;
                for(ii=0;ii<n1;ii++)
                {
                  ind2 = ind1+ii;
                  fout << Sdef[iii][ind2].x << setw(25) << Sdef[iii][ind2].y << setw(25) << outp[iii][ind2] << endl;
                }
                fout << endl;
             }
           }
         }
         else
         {
           for(jj=0;jj<n2;jj++)
           {
             for(iii=0;iii<Npatch;iii++)
             {
               n1 = uu[iii].n;

               if(iii == 0)
                 kk = 0;
               else
                 kk = 1;

               ind1 = jj * n1;

               for(ii=kk;ii<n1;ii++)
               {
                 ind2 = ind1+ii;
                 fout << Sorig[iii][ind2].x << setw(25) << Sorig[iii][ind2].y << setw(25) << outp[iii][ind2] << endl;
               }
               fout << endl;
             }
          }
        }

	      
	      /*
	      n1 = uu[0].n;
	      n2 = vv[0].n;

	      if(Npatch > 1)
		for(iii=1;iii<Npatch;iii++)
		   n1 += (uu[iii].n - 1);

              fout << n2 << '\t' << n1 << '\t' << 0.0 << endl;
              if(defUndefFlag)
	      {
	          for(jj=0;jj<n2;jj++)
	          {
          	      for(iii=0;iii<Npatch;iii++)
          	      {
          	          n1 = uu[iii].n;

          	          if(iii == 0)
          	            kk = 0;
          	          else
          	            kk = 1;

          	          ind1 = jj * n1;

	                  for(ii=kk;ii<n1;ii++)
	                  {
	                     ind2 = ind1+ii;
	                     fout << Sdef[iii][ind2].x << setw(25) << Sdef[iii][ind2].y << setw(25) << outp[iii][ind2] << endl;
	                  }
	                  fout << endl;
	              }
                  }
	     }
             else
	     {
                 for(jj=0;jj<n2;jj++)
	         {
                     for(iii=0;iii<Npatch;iii++)
          	     {
          	          n1 = uu[iii].n;

          	          if(iii == 0)
          	            kk = 0;
          	          else
          	            kk = 1;

          	          ind1 = jj * n1;

	                  for(ii=kk;ii<n1;ii++)
	                  {
	                     ind2 = ind1+ii;
	                     fout << Sorig[iii][ind2].x << setw(25) << Sorig[iii][ind2].y << setw(25) << outp[iii][ind2] << endl;
	                  }
	                  fout << endl;
	              }
	         }
             }
*/
             fout.close();
     }


/*
     if(index == 8) // print oup value for contour plots
     {
         ofstream fout("contour-data.dat");

	      if(fout.fail())
	      {
		 cout << " Could not open the Output file" << endl;
		 exit(1);
	      }

	      fout.setf(ios::fixed);
	      fout.setf(ios::showpoint);
	      fout.precision(10);
	      
	      VectorArray<double> u1;
              create_vector2(CurveResult[0].U, 10, u1);
              //findunique(CurveResult[0].U, u1);
              
	      EPOINT  EP;

              for(ii=0;ii<u1.n;ii++)
              {
                 EP  =  CurveResult[0].CurvePoint(u1[ii]).CalcEuclid();
                 fout << EP.x << setw(16) << EP.y << endl;
              }

             fout.close();
     }
*/
     if(index == 9) // write matrix pattern to an output file 
     {
         int ind2, *rr, *cc;

         ofstream fout("matrix-pattern.dat");

         if(fout.fail())
	 {
	    cout << " Could not open the Output file" << endl;
	    exit(1);
         }

	 fout.setf(ios::fixed);
	 fout.precision(1);

         fout << ((SolverSparse*)solver)->mtx.nRow << setw(10) << ((SolverSparse*)solver)->mtx.nCol << endl;

         ind2 = ((SolverSparse*)solver)->mtx.x.n;

         rr   = &(((SolverSparse*)solver)->mtx.row[0]);
         cc   = &(((SolverSparse*)solver)->mtx.col[0]);

         for(ii=0;ii<ind2;ii++)
            fout << rr[ii]-1 << setw(10) << cc[ii]-1 << endl;


         fout.close();
     }
     
     if(index == 10)
     {
        NurbsBaseFinal[0]->createAndWriteBasisFunctions(1);
     }


  return;
}





void IsogeometricFEM::printComputerTime(bool reset, int detailFlg)
{

  printf("----------------------------------------------------\n");
  printf("IsogeometricFEM::calcStiffRes:%7.3f sec ->%5.1f %\n",
               ctimCalcStiffRes, ctimCalcStiffRes/ctimSinceLastCall*100.);

  printf("IsogeometricFEM::factSolvUpdt:%7.3f sec ->%5.1f %\n",
               ctimFactSolvUpdt, ctimFactSolvUpdt/ctimSinceLastCall*100.);

  if(reset)
  {
     ctimFactSolvUpdt = 0.;
     ctimCalcStiffRes = 0.;
  }

  return;
}






void IsogeometricFEM::preparePatchElemProp()
{
  char fct[] = "IsogeometricFEM::preparePatchElemProp";

  if (patchElemProp.n < 1) prgError(1,fct,"'patch element property data ' missing!");


  char *elmTypeNames[] = ISOGEOM_ELEMENT_TYPE_NAMES;

  for (int ii=0; ii<patchElemProp.n; ii++)
  {
    // assign correct id of the element type (as stored in the database) based on the element name
    patchElemProp[ii].id = patchElemProp[ii].name.which(elmTypeNames);

    if (patchElemProp[ii].id < 0)
       prgError(2,fct,"unknown element type name!");

  }

  return;
}





void IsogeometricFEM::preparePatchMatlProp()
{
  char fct[] = "IsogeometricFEM::preparePatchMatlProp";

  if (patchMatlProp.n < 1) prgError(1,fct,"'patch material property data ' missing!");


  char *matlTypeNames[] = MATERIAL_TYPE_NAMES;

  for (int ii=0; ii<patchMatlProp.n; ii++)
  {
    // assign correct id of the material type (as stored in the database) based on the material name
    patchMatlProp[ii].id = patchMatlProp[ii].name.which(matlTypeNames);

    if (patchMatlProp[ii].id < 0) prgError(2,fct,"unknown element type name!");
  }

  return;
}






void IsogeometricFEM::findMinMaxX(double *xmn, double *xmx, bool defFlg)
{
  double xmin, xmax, ymin, ymax, xx, yy;
  xmin = xmax = x.x[0];
  ymin = ymax = x.x[1];

  for(int ii=1;ii<numnp;ii++)
  {
     xx = x.x[(ndm+1)*ii + 0];
     yy = x.x[(ndm+1)*ii + 1];
     if(xx < xmin)       xmin = xx;
     if(xx > xmax)       xmax = xx;
     if(yy < ymin)       ymin = yy;
     if(yy > ymax)       ymax = yy;

  }
  xmn[0] = xmin;  xmn[1] = ymin;

  xmx[0] = xmax;  xmx[1] = ymax;

  return;
}






void IsogeometricFEM::findMinMaxResult(double *xmn, double *xmx, bool defFlg)
{
   double xmin, xmax, ymin, ymax, xx, yy;
   EPOINT EP;
   int ngbf1, ngbf2;
   ngbf1 = SurfaceResult[0].ngbf1;
   ngbf2 = SurfaceResult[0].ngbf2;

   EP = SurfaceResult[0].Pw[0][0].CalcEuclid();
   xmin = xmax = EP.x;
   ymin = ymax = EP.y;

   for(int jj=0;jj<ngbf2;jj++)
   {
      for(int ii=0;ii<ngbf1;ii++)
      {
         EP = SurfaceResult[0].Pw[ii][jj].CalcEuclid();
         xx = EP.x;
         yy = EP.y;

         if(xx < xmin)       xmin = xx;
         if(xx > xmax)       xmax = xx;

         if(yy < ymin)       ymin = yy;
         if(yy > ymax)       ymax = yy;
      }
   }

  xmn[0] = xmin;  xmn[1] = ymin;

  xmx[0] = xmax;  xmx[1] = ymax;

  return;
}






void IsogeometricFEM::PostprocessForCurves()
{
  char fct[] = "IsogeometricFEM::PostprocessForCurves";

  cout << "     Postprocessing for Curves .... " << endl;
  cout << endl;


  CurveResult.setDim(CurveListFinal.n);

  for(int ii=0;ii<CurveListFinal.n;ii++)
  {
     CurveResult[ii].Pw = CurveListFinal[ii].Pw;
     CurveResult[ii].U = CurveListFinal[ii].U;
     CurveResult[ii].p = CurveListFinal[ii].p;

     for(int jj=0;jj<CurveResult[ii].Pw.n;jj++)
     {
        CurveResult[ii].Pw[jj].y = solnFull[jj];
     }
  }

  cout << "     Postprocessing for Curves Over.... " << endl;
  cout << endl;

  return;
}




void IsogeometricFEM::PostprocessForSurfaces()
{
  char fct[] = "IsogeometricFEM::PostprocessForSurfaces";

  cout << "     Postprocessing for Surfaces .... " << endl;
  cout << endl;

 // CalcStrainStress();

  return;
}


void IsogeometricFEM::CalcStrainStress1(int rsys1)
{
  cerr << " Not implemented yet ..." << endl;

  return;
}



void IsogeometricFEM::setDifferentFlags(int index, int value)
{
   switch(index)
   {
       case 1:
              RSYS = value;

              break;
       case 2:
              defUndefFlag = !(value == 0);

              break;

   }

   return;
}



void IsogeometricFEM::printDifferentFlags()
{
   printf("\n");
   printf("\t  Current     'RSYS'       value : %2d \n", RSYS);
   printf("\t  Current  'defUndefFlag'  value : %2d \n", defUndefFlag);
   printf("\n");

   return;
}



void IsogeometricFEM::CalcStrainStress(int rsys1)
{
  char fct[] = "IsogeometricFEM::CalcStrainStress";

 // cout << "     ISOGEOMETRICFEM: calculating strains and stresses ...\n\n";
 // cout << endl;


//   cout << '\t' << epspeqv2[0] << endl;   cout << endl;
//   cout << '\t' << epspeqv2[1] << endl;   cout << endl;


  return;
}






void  IsogeometricFEM::printResultAtPoint(int patchnum, double u1, double v1, double w1)
{

  int  elemtype = patchGrp[0].eltype;

  if(!(elemtype == 14) && !(elemtype == 15) )
  {
     VectorArray<double>  vals;
     vals.setDim(3); vals.zero();

     EPOINT EP1, EP2;
     
     if(patchGrp[0].ndom == 1)
     {
        EP1 = CurveListFinal[patchnum].CurvePoint(u1).CalcEuclid();
        EP2 = CurveResult[patchnum].CurvePoint(u1).CalcEuclid();
     }
     if(patchGrp[0].ndom == 2)
     {
        EP1 = SurfaceListFinal[patchnum].SurfacePoint(u1, v1).CalcEuclid();
        EP2 = SurfaceResult[patchnum].SurfacePoint(u1, v1).CalcEuclid();
     }
     if(patchGrp[0].ndom == 3)
     {
        EP1 = SolidListFinal[patchnum].SolidPoint(u1, v1, w1).CalcEuclid();
        EP2 = SolidResult[patchnum].SolidPoint(u1, v1, w1).CalcEuclid();
     }

     vals[0] = EP2.x - EP1.x;
     vals[1] = EP2.y - EP1.y;
     vals[2] = EP2.z - EP1.z;

     printf("\n\n");
     printf("\t Patch Number = %5d \n", patchnum);
     printf("\t Output Values at parameters u = %8.6f \t v = %8.6f are \n", u1, v1);
     printf("\n\n");
     printf("\t X-Coord     = %12.6f\t%12.6f \n", EP1.x, EP2.x);
     printf("\t Y-Coord     = %12.6f\t%12.6f \n", EP1.y, EP2.y);
     //printf("\t Z-Coord     = %12.6f\t%12.6f \n", EP1.z, EP2.z);
     printf("\n\n");
     printf("\t X-disp      = %12.8f \n", vals[0]);
     printf("\t Y-disp      = %12.8f \n", vals[1]);
     //printf("\t Z-disp      = %12.8f \n", vals[2]);
     printf("\n\n");

    char        tmp[50];
    MyString    tmpStr, fname;

    sprintf(tmp,"\t %12.6f \t %12.6f ", vals[0], vals[1]);
    tmpStr.append(tmp);

    prgWriteToTFile(tmpStr);

    //printf("\t %12.6f\n", sqrt(EP1.x*EP1.x + EP1.y*EP1.y));
    /*
    int  dd, ii, jj, kk, ll, count, nlocal, ind1, ind2, dir, resln[2];
    
    resln[0] = 1;

    VectorArray<double>  U, uu;
    vector<vector<double> >  outp2;
    vector<double>  xx;

    outp2.resize(4);

    create_vector2(CurveListFinal[0].U, resln[0], uu);

    for(ii=0;ii<uu.n;ii++)
    {
        EP1 = CurveListFinal[0].CurvePoint(uu[ii]).CalcEuclid();

        outp2[0].push_back(CurveResult[0].computeValue(1,uu[ii]));

        xx.push_back(EP1.x);
    }

     // prepare and write a file to postprocess in matplotlib

      ofstream fout("post-process-curve.dat");

      if(fout.fail())
      {
         cout << " Could not open the Output file" << endl;
      exit(1);
      }

      fout.setf(ios::fixed);
      fout.setf(ios::showpoint);
      fout.precision(8);

      for(ii=0;ii<xx.size();ii++)
         fout << xx[ii] << '\t' << outp2[0][ii] << endl;

      fout.close();
      */
  }
  else
  {
    int  dd, ii, jj, kk, ll, count, nlocal, ind1, ind2, dir, resln[2];
    
    resln[0] = (int) u1;
    resln[1] = (int) v1;

    VectorArray<double>  U, V, uu, vv;
    vector<vector<double> >  outp2;
    vector<double>  xx, yy;
    
    outp2.resize(4);

    create_vector2(SurfaceListFinal[0].U, resln[0], uu);
    create_vector2(SurfaceListFinal[0].V, resln[1], vv);

    EPOINT EP1, EP2;

    for(ii=0;ii<uu.n;ii++)
    {
        EP1 = SurfaceListFinal[0].SurfacePoint(uu[ii], 0.5).CalcEuclid();
        //EP2 = SurfaceResult[0].SurfacePoint(uu[ii], 0.5).CalcEuclid();

        outp2[0].push_back(SurfaceResult[0].computeValue(1,uu[ii], 0.5));
        outp2[1].push_back(SurfaceResult[0].computeValue(2,uu[ii], 0.5));

        //outp2[0].push_back(EP2.x-EP1.x);
        //outp2[1].push_back(EP2.y-EP1.y);

        xx.push_back(EP1.x);
    }

    for(jj=0;jj<vv.n;jj++)
    {
        EP1 = SurfaceListFinal[0].SurfacePoint(0.5, vv[jj]).CalcEuclid();
        //EP2 = SurfaceResult[0].SurfacePoint(0.5, vv[jj]).CalcEuclid();

        //outp2[2].push_back(EP2.x-EP1.x);
        //outp2[3].push_back(EP2.y-EP1.y);

        outp2[2].push_back(SurfaceResult[0].computeValue(1, 0.5, vv[jj]));
        outp2[3].push_back(SurfaceResult[0].computeValue(2, 0.5, vv[jj]));


        yy.push_back(EP1.y);
   }
       
      // prepare and write a file to postprocess in matplotlib

      ofstream fout("post-process-curve.dat");

      if(fout.fail())
      {
         cout << " Could not open the Output file" << endl;
      exit(1);
      }

      fout.setf(ios::fixed);
      fout.setf(ios::showpoint);
      fout.precision(8);

      for(ii=0;ii<xx.size();ii++)
      {
         fout << xx[ii] << '\t' << outp2[0][ii] << '\t' << outp2[1][ii] << '\t' << yy[ii] << '\t' << outp2[2][ii] << '\t' << outp2[3][ii] << endl;
      }

      fout.close();
  }

/*
    VectorArray<double> paraincrs;

    double incr = 0.1, param;

    create_vector(SurfaceResult[0].V[0], SurfaceResult[0].V[SurfaceResult[0].V.n-1], incr, paraincrs);

    int nlbf = SurfaceListFinal[0].nlbf;

    VectorArray<double> N;
    N.setDim(nlbf);
 
    int startindex[2], ii, jj;

    double dN_dx[nlbf], dN_dy[nlbf], para[2];

    double F[4], detF, Jac, sxx, sxy, BULK, mu, lambda, volstr;

    double *matDat = &(SurfaceListFinal[patchnum].MatlProp.data[0]);

    BULK = matDat[0];
    mu = matDat[1];
    lambda = BULK - 2.0*mu/3.0;

    para[0] = u1;

    startindex[0] = FindSpan(&(SurfaceListFinal[0].U[0]), SurfaceListFinal[0].U.n, SurfaceListFinal[0].p, para[0]) - SurfaceListFinal[0].p;

    cout << para[0] << '\t' << startindex[0] << '\t' << startindex[1] << endl;

    for(ii=0;ii<paraincrs.n;ii++)
    {
      para[1] = paraincrs[ii];

      startindex[1] = FindSpan(&(SurfaceListFinal[0].V[0]), SurfaceListFinal[0].V.n, SurfaceListFinal[0].q, paraincrs[ii]) - SurfaceListFinal[0].q;

      SurfaceListFinal[0].ShapeFunDerivatives(&(startindex[0]), &(para[0]), dN_dx, dN_dy, Jac);
      SurfaceResult[0].deformationGradient(startindex[0], startindex[1], 1, dN_dx, dN_dy, F, detF);

      //for(jj=0;jj<nlbf;jj++)
        //printf("%12.8f \t %12.8f \n", dN_dx[jj], dN_dy[jj]);
      //printf("\n");

      volstr = F[0] + F[3] - 2.0;
      sxx = 2.0*mu*(F[0]-1.0) + lambda*volstr;
      sxy = mu*(F[1]+F[2]);

      printf("%12.8f \t %12.8f \t %12.8f \n", paraincrs[ii], sxx, sxy);
    }
    printf("\n\n");
*/

  return;
}






void IsogeometricFEM::printDispsAtParameter(int patchnum, int dir, double param, double incr)
{
    VectorArray<double> paraincrs;

 //   cout << " dir " << dir << endl;
 //   cout << " param " << param << endl;
 //   cout << " incr " << incr << endl;

    if(dir == 1)
       create_vector(SurfaceResult[patchnum].V[0], SurfaceResult[patchnum].V[SurfaceResult[patchnum].V.n-1], incr, paraincrs);
    else
       create_vector(SurfaceResult[patchnum].U[0], SurfaceResult[patchnum].U[SurfaceResult[patchnum].U.n-1], incr, paraincrs);

    int nlbf = SurfaceListFinal[patchnum].nlbf;

    ListArray<VectorArray<double> > vals;

    vals.setDim(2);
    for(int ii=0;ii<2;ii++)
      vals[ii].setDim(paraincrs.n);

    EPOINT EP1, EP2;

    for(int ii=0;ii<paraincrs.n;ii++)
    {
        if(dir == 1)
        {
           EP1 = SurfaceListFinal[patchnum].SurfacePoint(param, paraincrs[ii]).CalcEuclid();
           EP2 = SurfaceResult[patchnum].SurfacePoint(param, paraincrs[ii]).CalcEuclid();
        }
        else
        {
           EP1 = SurfaceListFinal[patchnum].SurfacePoint(paraincrs[ii], param).CalcEuclid();
           EP2 = SurfaceResult[patchnum].SurfacePoint(paraincrs[ii], param).CalcEuclid();
        }
        
        EP1.print2screen();
        cout << sqrt(EP1.x*EP1.x + EP1.y*EP1.y) << '\t' << atan(EP1.y/EP1.x)*180.0/PI << endl;

        vals[0][ii] = EP2.x - EP1.x;
        vals[1][ii] = EP2.y - EP1.y;

    }

       cout << '\t' << " Output Values at parameter = " << param << "  in direction = " << dir << endl;
       cout << endl;
       cout << endl;
       cout << '\t' << " S.No " << '\t' << " param value " << '\t' << " X-disp  " << '\t' << " Y-disp " << endl;
       cout << endl;
       for(int ii=0;ii<paraincrs.n;ii++)
       {
          cout << '\t' << ii << '\t' << paraincrs[ii] << '\t' << vals[0][ii] << '\t' << vals[1][ii] << endl;
       }
       cout << endl;
       cout << endl;


  return;
}



void IsogeometricFEM::printStrainsAtParameter(int patchnum, int dir, double param, double incr)
{

  return;
}


/*
void IsogeometricFEM::printStressesAtParameter(int patchnum, int dir, double param, double incr)
{
    VectorArray<double> paraincrs;

    if(dir == 1)
       create_vector(SurfaceResult[patchnum].V[0], SurfaceResult[patchnum].V[SurfaceResult[patchnum].V.n-1], incr, paraincrs);
    else
       create_vector(SurfaceResult[patchnum].U[0], SurfaceResult[patchnum].U[SurfaceResult[patchnum].U.n-1], incr, paraincrs);

    int nlbf = SurfaceListFinal[patchnum].nlbf;

    ListArray<VectorArray<double> > vals;

    VectorArray<double> N;
    VectorArray<int> ind;
    N.setDim(nlbf);
    ind.setDim(nlbf);
    vals.setDim(5);
    for(int ii=0;ii<5;ii++)
    {
       vals[ii].setDim(paraincrs.n);
       vals[ii].zero();
    }

   double *elmDat = &(SurfaceListFinal[patchnum].ElemProp.data[0]),
          *matDat = &(SurfaceListFinal[patchnum].MatlProp.data[0]);

   double F[4]={0.0, 0.0, 0.0, 0.0}, detF=0.0, F33=0.0, stre[4], cc[4][4], dvol=0.0, dt=0.0, detJ=0;

   double rad=0.0, rad0=0.0;

   int   finInt = (int) elmDat[2],
         sss    = (int) elmDat[3],
         matId  = SurfaceListFinal[patchnum].MatlProp.id + 1,
         err    = 0,
         isw    = 3,
         count  = 1,
         nivGP  = elem[0]->nivGP;


    bool   finite = (finInt >= 1),
           axsy   = (sss == 3);

    double *iv1, *iv2;
    iv1 = new double[nivGP];
    iv2 = new double[nivGP];


    for(int ii=0;ii<paraincrs.n;ii++)
    {
        if(intVarFlag)
        {
            double dmy[10];

            int  dmyI[3],
               isw   = 2,
               mDim  = matdim_(&matId);

            if (mDim == 1) matlib1d_(matDat,dmy,dmy,dmy,dmy,iv1,iv2,dmy,
       	                         &matId,&nivGP,&finInt,&sss,&isw,dmyI);

            else if (mDim == 2) matlib2d_(matDat,dmy,dmy,dmy,dmy,iv1,iv2,dmy,
                                    &matId,&nivGP,&finInt,&sss,&isw,dmyI);

            else if (mDim == 3) matlib3d_(matDat,dmy,dmy,dmy,iv1,iv2,dmy,
	                             &matId,&nivGP,&finInt,&isw,dmyI);

            else prgError(1,"CalcStrains--initialiseIntVar","invalid value of ndm!");
        }

//
        //  COMPUTE SHAPE FUNCTIONS, THEIR DERIVATIVES AND THE DEFORMATION GRADIENT "F"
        if(dir == 1)
        {
           NurbsShapeFunctions2DPost2(&SurfaceListFinal[patchnum], &SurfaceResult[patchnum], param, paraincrs[ii], N, ind, F, detF, detJ);
           //   for axisymmetric problems compute radius
           if(axsy)
           {
              rad0 = SurfaceListFinal[patchnum].SurfacePoint(param, paraincrs[ii]).CalcEuclid().x;
              rad  = SurfaceResult[patchnum].SurfacePoint(param, paraincrs[ii]).CalcEuclid().x;
           }
        }
        else
        {
           NurbsShapeFunctions2DPost2(&SurfaceListFinal[patchnum], &SurfaceResult[patchnum], paraincrs[ii], param, N, ind, F, detF, detJ);
           //   for axisymmetric problems compute radius
           if(axsy)
           {
              rad0 = SurfaceListFinal[patchnum].SurfacePoint(paraincrs[ii], param).CalcEuclid().x;
              rad  = SurfaceResult[patchnum].SurfacePoint(paraincrs[ii], param).CalcEuclid().x;
           }
        }
//
           //  ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

           if(sss == 1)  // plane stress
           {
             if(finite)
               F33 = 1.0/sqrt(detF);
             else
               F33 = 3.0 - F[0] - F[3];
           }
           else if(sss == 2)    // plane strain
             F33 = 1.0;
           else //(sss == 3)    // axisymmetric
             F33 = rad/rad0;

           // COMPUTE MATERIAL RESPONSE

           for(int pp=0;pp<4;pp++)
           {
              stre[pp] = 0.0;
              for(int qq=0;qq<4;qq++)
                cc[pp][qq] = 0.0;
           }

           matlib2d_(matDat, F, &F33, stre, cc[0], iv1, iv2, &dt, &matId, &nivGP, &finInt, &sss, &isw, &err, &count, (Element*) this);

           vals[0][ii] = stre[0]; // sxx
           vals[1][ii] = stre[1]; // syy
           vals[2][ii] = stre[3]; // sxy
           vals[3][ii] = stre[2]; // szz

           vals[4][ii] = sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6* stre[3]*stre[3])/2); // seqv
    }

    cout.precision(5);
    cout << '\t' << " Output Values at parameter = " << param << "  in direction = " << dir << endl;
    cout << endl;
    cout << '\t' << " S.No " << '\t' << " param value " << '\t' << " Stress-XX " << '\t' << " Stress-YY " << '\t' << " Stress-XY " << '\t' << " Stress-ZZ " << '\t' << " Stress-EQV "<< endl;
    cout << endl;
    for(int ii=0;ii<paraincrs.n;ii++)
       cout << '\t' << ii << '\t' << paraincrs[ii] << '\t' << vals[0][ii] << '\t' << vals[1][ii] << '\t' << vals[2][ii] << '\t' << vals[3][ii] << '\t' << vals[4][ii] << endl;
    cout << endl;
    cout << endl;


  return;
}
*/


void IsogeometricFEM::printStressesAtParameter(int patchnum, int dir, double param, double incr)
{
  cout << patchnum << '\t' << dir << endl;
  
  double *elmDat = &(SurfaceListFinal[patchnum].ElemProp.data[0]),
         *matDat = &(SurfaceListFinal[patchnum].MatlProp.data[0]);

  double F[4]={0.0, 0.0, 0.0, 0.0}, detF=0.0, F33=0.0, stre[4], cc[4][4], dvol=0.0, dt=0.0, detJ, Jac;

  int   finInt = (int) elmDat[2],
        sss    = (int) elmDat[3],
        matId  = SurfaceListFinal[patchnum].MatlProp.id + 1,
        err    = 0,
        isw    = 3,
        count  = 1,
        nivGP  = elem[0]->nivGP;

    double *iv1, *iv2;
    //iv1 = new double[nivGP];
    //iv2 = new double[nivGP];


    bool   finite = (finInt >= 1),
           axsy   = (sss == 3);

    int  dd, ii, jj, kk, ll, index, ind1, ind2, n1, n2, pp, qq, start[2];
    
    int  p = SurfaceListFinal[patchnum].p;
    int  q = SurfaceListFinal[patchnum].q;
    
    int  nel1 = SurfaceListFinal[patchnum].nelem1;
    int  nel2 = SurfaceListFinal[patchnum].nelem2;

    int nn=0;
    if(patchnum == 1)
      nn = SurfaceListFinal[0].nelem1 * SurfaceListFinal[0].nelem2;

    int  nlbf = (p+1)*(q+1);

    double   incr1, incr2, start1, start2, knots[2], dummy, pres;
    vector<double>   N(nlbf), dN_dx(nlbf), dN_dy(nlbf) ;


    VectorArray<double>  uu, vv;
    vector<vector<double> >  outp2, xx, yy;

    create_vector(0.0,1.0,0.01,uu);
    vv = uu;

    outp2.resize(vv.n);
    xx.resize(vv.n);
    yy.resize(vv.n);

    EPOINT EP1, EP2;

    for(jj=0;jj<vv.n;jj++)
    {
        outp2[jj].resize(uu.n);
        xx[jj].resize(uu.n);
        yy[jj].resize(uu.n);

        knots[1] = vv[jj];

        for(ii=0;ii<uu.n;ii++)
        {
          knots[0] = uu[ii];

          //EP1 = SurfaceListFinal[patchnum].SurfacePoint(uu[ii], vv[jj]).CalcEuclid();
          //if(finInt)
            EP1 = SurfaceResult[patchnum].SurfacePoint(uu[ii], vv[jj]).CalcEuclid();

          //outp2[0].push_back(EP2.x-EP1.x);
          //outp2[1].push_back(EP2.y-EP1.y);

          xx[jj][ii] = EP1.x;
          yy[jj][ii] = EP1.y;

          start[0] = FindSpan(&(SurfaceListFinal[patchnum].U[0]), SurfaceListFinal[patchnum].U.n, p, knots[0]) - p;
          start[1] = FindSpan(&(SurfaceListFinal[patchnum].V[0]), SurfaceListFinal[patchnum].V.n, q, knots[1]) - q;

          //cout << " indices ... " << start[0] << '\t' << start[1] << endl;
          
          iv1 = elem[nn + nel1*start[1] + start[0]]->intVar1;
          iv2 = elem[nn + nel1*start[1] + start[0]]->intVar2;

          SurfaceListFinal[patchnum].ShapeFunDerivatives(start, knots, &N[0], &dN_dx[0], &dN_dy[0], Jac);

          SurfaceResult[patchnum].deformationGradient(start[0], start[1], 1, &dN_dx[0], &dN_dy[0], F, detF);

          if(sss == 1)  // plane stress
          {
            if(finite)
              F33 = 1.0/sqrt(detF);
            else
              F33 = 3.0 - F[0] - F[3];
          }
          else if(sss == 2)    // plane strain
            F33 = 1.0;

          // COMPUTE MATERIAL RESPONSE

          for(pp=0;pp<4;pp++)
          {
            stre[pp] = 0.0;
            for(qq=0;qq<4;qq++)
              cc[pp][qq] = 0.0;
          }

          matlib2d_(matDat, F, &F33, stre, cc[0], iv1, iv2, &dt, &matId, &nivGP, &finInt, &sss, &isw, &err, &count, (Element*) this);

          pres = (stre[0]+stre[1]+stre[2])/3.0 ;

          if( SurfaceListFinal[patchnum].ElemProp.id == 7)
          {
            pres = surfSecondVar[patchnum].computeValue(1, knots[0], knots[1]);

            dummy = pres - (stre[0]+stre[1]+stre[2])/3.0 ;

            stre[0] += dummy;
            stre[1] += dummy;
            stre[2] += dummy;
          }

          if( SurfaceListFinal[patchnum].ElemProp.id == 21)
          {
            pres = SurfaceResult[patchnum].computeValue(3, knots[0], knots[1]);
            
            //pres = SurfaceListFinal[patchnum].computeValue(3, knots[0], knots[1]);

            dummy = pres - (stre[0]+stre[1]+stre[2])/3.0 ;

            stre[0] += dummy;
            stre[1] += dummy;
            stre[2] += dummy;
          }

          if(dir < 4)
            outp2[jj][ii] = stre[dir];
          else if (dir == 4)
            outp2[jj][ii] = pres;
          else
            outp2[jj][ii] = sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6.0* stre[3]*stre[3])/2.0); // seqv

        }
    }

      ofstream fout("post-process.dat");

      if(fout.fail())
      {
        cout << " Could not open the Output file" << endl;
        exit(1);
      }

      fout.setf(ios::fixed);
      fout.setf(ios::showpoint);
      fout.precision(8);
      
      fout << uu.n << '\t' << vv.n << '\t' << 0.0 << endl;

      for(jj=0;jj<xx.size();jj++)
      {
        for(ii=0;ii<xx[jj].size();ii++)
          fout << xx[ii][jj] << '\t' <<  yy[ii][jj] << '\t' << outp2[ii][jj] << endl;
        fout << endl; fout << endl;
      }

      fout.close();

  return;
}



void IsogeometricFEM::GenerateConnectivityArrays()
{
  if(ndm == 1)
    GenerateConnectivityArraysCurve();
  else if(ndm == 2)
    GenerateConnectivityArraysSurface();
  else
    GenerateConnectivityArraysSolid();

  return;
}



void IsogeometricFEM::GenerateConnectivityArraysCurve()
{
  cout << "     ISOGEOMETRICFEM: generating cinnectivity arrays for curve ...\n\n";
  char fct[] = "IsogeometricFEM::GenerateConnectivityArraysCurve";

   int iii, ii, jj, ngbf, count, ee;

   ntoteqs1=0;
   for(iii=0;iii<CurveListFinal.n;iii++)
     CurveListFinal[iii].GenerateConnectivityArrays1(ntoteqs1);

   solver->rhsVec.resize(ntoteqs1);
   solver->rhsVec.setZero();

   ForceVec.setDim(ntoteqs1);
   ForceVec.zero();
   assy4r = ForceVec;

   // copy Uinit values from  CurveListFinal -> CurveResult
   for(iii=0;iii<CurveListFinal.n;iii++)
      CurveResult[iii].Uinit = CurveListFinal[iii].Uinit;


   // calculate LM Arrays
   for(iii=0;iii<CurveListFinal.n;iii++)
      CurveListFinal[iii].GenerateConnectivityArrays2();

//   for(int iii=0;iii<CurveListFinal.n;iii++)      CurveListFinal[iii].printConnectivityArrays();

     // calculate gbf numbers for patch #1
     ngbf = CurveListFinal[0].ngbf;
     count = 0;
     for(jj=0;jj<ngbf;jj++)
        CurveListFinal[0].gbfnums[jj] = jj;


     if(Npatch > 1)
     {
        // compute gbfnums for other patches
        count = CurveListFinal[0].ngbf;
        int patchnum, side, ngbf31, ngbf32, index1, index2;
        for(iii=1;iii<Npatch;iii++)
        {        }
     }


        int index1, index2, *U1, *U2;
        for(iii=0;iii<Npatch;iii++)
        {
            U1 = &(CurveListFinal[iii].gbfnums[0]);
            U2 = &(CurveListFinal[iii].toUfull[0]);
            for(ii=0;ii<CurveListFinal[iii].ngbf;ii++)
            {
               index1 = U1[ii]*ndf;
               index2 = ii*ndf;
               for(jj=0;jj<ndf;jj++)
                 U2[index2+jj] = index1 + jj;
            }
        }


        solnInit.zero();
        int count1 = -1, *U4;
        double *U3;
        for(iii=0;iii<Npatch;iii++)
        {
            U4 = &(CurveListFinal[iii].toUfull[0]);
            U3 = &(CurveListFinal[iii].Uinit[0]);
            for(ii=0;ii<CurveListFinal[iii].Uinit.n;ii++)
            {
               if(U4[ii] > count1)
                 solnInit[U4[ii]] = U3[ii];
            }
           count1 += (CurveListFinal[iii-1].ngbf * ndf);
        }

        for(iii=0;iii<Npatch;iii++)
        {
           CurveResult[iii].toUfull = CurveListFinal[iii].toUfull;
           CurveResult[iii].gbfnums = CurveListFinal[iii].gbfnums;
        }


//  computerTime.stopAndPrint(fct);

  return;
}






void IsogeometricFEM::GenerateConnectivityArraysSurface()
{
   cout << "     ISOGEOMETRICFEM: generating cinnectivity arrays for surface ...\n\n";
   char fct[] = "IsogeometricFEM::GenerateConnectivityArraysSurface";

   int iii, ii, jj, kk;

   ntoteqs1=0;
   for(iii=0;iii<SurfaceListFinal.n;iii++)
     SurfaceListFinal[iii].GenerateConnectivityArrays1(ntoteqs1);

   cout << '\t' << " ntoteqs1  " << ntoteqs1 << endl;

//   for(iii=0;iii<SurfaceListFinal.n;iii++)     SurfaceListFinal[iii].printConnectivityArrays();

     // adjust ID arrays based on interface edges

     if(Npatch > 1)
     {
        int  patchnum, side, ngbf21, ngbf22, ngbf31, ngbf32, index1, index2, dir;

        for(iii=1;iii<Npatch;iii++)
        {
            ngbf21 = SurfaceListFinal[iii].ngbf1;
            ngbf22 = SurfaceListFinal[iii].ngbf2;

            for(kk=0;kk<4;kk++)
            {
               if(SurfaceListFinal[iii].edgedata[kk] != -1)
               {
                  patchnum = SurfaceListFinal[iii].intfdata[2*kk];
                  side     = SurfaceListFinal[iii].intfdata[2*kk+1];

                  ngbf31 = SurfaceListFinal[patchnum].ngbf1;
                  ngbf32 = SurfaceListFinal[patchnum].ngbf2;

                  switch(kk)
                  {
                      case 0: // left side

                               if(ngbf22 != ngbf32)
                                  prgError(1,fct," ngbf2 in the 2 patches do not match");

                               for(jj=0;jj<ngbf22;jj++)
                               {
                                  index1 = ngbf21*jj;
                                  index2 = ngbf31*(jj+1)-1;
                                  for(dir=0;dir<2;dir++)
                                  {
                                     SurfaceListFinal[iii].ID[dir][index1] = SurfaceListFinal[patchnum].ID[dir][index2];
                                     SurfaceListFinal[iii].dispBCs[index1][dir] = SurfaceListFinal[patchnum].dispBCs[index2][dir];

                                     SurfaceListFinal[iii].Uinit[ndf*index1 + dir] = SurfaceListFinal[patchnum].Uinit[ndf*index2 + dir] ;
                                  }
                               }

                              break;

                      case 1: // right side
                              break;

                               if(ngbf22 != ngbf32)
                                  prgError(1,fct," ngbf2 in the 2 patches do not match");

                               for(jj=0;jj<ngbf22;jj++)
                               {
                                  index2 = ngbf21*jj;
                                  index1 = ngbf31*(jj+1)-1;
                                  for(dir=0;dir<2;dir++)
                                  {
                                     SurfaceListFinal[iii].ID[dir][index1] = SurfaceListFinal[patchnum].ID[dir][index2];
                                     SurfaceListFinal[iii].dispBCs[index1][dir] = SurfaceListFinal[patchnum].dispBCs[index2][dir];

                                     SurfaceListFinal[iii].Uinit[ndf*index1 + dir] = SurfaceListFinal[patchnum].Uinit[ndf*index2 + dir] ;
                                  }
                               }

                      case 2: // bottom side

                               if(ngbf21 != ngbf31)
                                  prgError(2,fct," ngbf1 in the 2 patches do not match");

                               for(jj=0;jj<ngbf21;jj++)
                               {
                                  index1 = jj;
                                  index2 = ngbf31*(ngbf32-1)+jj;
                                  for(dir=0;dir<2;dir++)
                                  {
                                     SurfaceListFinal[iii].ID[dir][index1] = SurfaceListFinal[patchnum].ID[dir][index2];
                                     SurfaceListFinal[iii].dispBCs[index1][dir] = SurfaceListFinal[patchnum].dispBCs[index2][dir];

                                     SurfaceListFinal[iii].Uinit[ndf*index1 + dir] = SurfaceListFinal[patchnum].Uinit[ndf*index2 + dir] ;
                                  }
                               }

                              break;

                      case 3: // top side

                               if(ngbf21 != ngbf31)
                                  prgError(2,fct," ngbf1 in the 2 patches do not match");

                               for(jj=0;jj<ngbf21;jj++)
                               {
                                  index1 = ngbf21*(ngbf22-1)+jj;
                                  index2 = jj;
                                  for(dir=0;dir<2;dir++)
                                  {
                                     SurfaceListFinal[iii].ID[dir][index1] = SurfaceListFinal[patchnum].ID[dir][index2];
                                     SurfaceListFinal[iii].dispBCs[index1][dir] = SurfaceListFinal[patchnum].dispBCs[index2][dir];

                                     SurfaceListFinal[iii].Uinit[ndf*index1 + dir] = SurfaceListFinal[patchnum].Uinit[ndf*index2 + dir] ;
                                  }
                               }

                              break;

                  }
               }
            }
        }
     }


   // copy Uinit values from  SurfaceListFinal -> SurfaceResult
   for(iii=0;iii<SurfaceListFinal.n;iii++)
      SurfaceResult[iii].Uinit = SurfaceListFinal[iii].Uinit;


   // calculate LM Arrays
   for(iii=0;iii<SurfaceListFinal.n;iii++)
      SurfaceListFinal[iii].GenerateConnectivityArrays2();

   //for(iii=0;iii<SurfaceListFinal.n;iii++)
     //SurfaceListFinal[iii].printConnectivityArrays();

/*
     cout << "          CCCCCCCCCCCCC    " << endl;
     cout << "          CCCCCCCCCCCCC    " << endl;
     cout << "          CCCCCCCCCCCCC    " << endl;

     cout << endl;
     cout << "       .... displacement values for the final mesh ... " << endl;
     cout << endl;

     for(iii=0;iii<Npatch;iii++)
     {
         cout << "       patch... : " << (iii+1) << endl;
         cout << endl;
         for(ii=0;ii<SurfaceListFinal[iii].dispBCs.n;ii++)
         {
            cout << '\t' << ii << '\t' ;
            for(jj=0;jj<SurfaceListFinal[iii].dispBCs[0].n;jj++)
              cout << SurfaceListFinal[iii].dispBCs[ii][jj] << '\t' ;
            cout << endl;
         }
         cout << endl;
     }

     cout << endl;
     cout << "       .... Uinit values for the final mesh ... " << endl;
     cout << endl;
     for(iii=0;iii<Npatch;iii++)
     {
         cout << "       patch... : " << (iii+1) << endl;
         cout << endl;
         for(ii=0;ii<SurfaceListFinal[iii].ngbf;ii++)
         {
            cout << '\t' << ii << '\t' ;
            for(int dof=0;dof<ndf;dof++)
              cout << SurfaceListFinal[iii].Uinit[ndf*ii+dof] << '\t' ;
            cout << endl;
         }
         cout << endl;
     }

     cout << endl;
     cout << "       .... primary variable values for the elements ... " << endl;
     cout << endl;
     for(int ee=0;ee<totnumel;ee++)
     {
         cout << "       elem... : " << (ee+1) << endl;

         cout << endl;
         elem[ee]->printPrimVariable();
         cout << endl;
     }
*/

     int ngbf21, ngbf22, count;
     // calculate gbf numbers for patch #1
     ngbf21 = SurfaceListFinal[0].ngbf1;
     ngbf22 = SurfaceListFinal[0].ngbf2;
     count = 0;
     for(jj=0;jj<ngbf22;jj++)
     {
       for(kk=0;kk<ngbf21;kk++)
       {
          SurfaceListFinal[0].gbfnums[count] = count;
          count++;
       }
     }

     if(Npatch > 1)
     {
        count = SurfaceListFinal[0].ngbf;
        int patchnum, side, ngbf31, ngbf32, index1, index2;


        for(iii=1;iii<Npatch;iii++)
        {
            ngbf21 = SurfaceListFinal[iii].ngbf1;
            ngbf22 = SurfaceListFinal[iii].ngbf2;

            for(kk=0;kk<4;kk++)
            {
               if(SurfaceListFinal[iii].edgedata[kk] != -1)
               {
                  patchnum = SurfaceListFinal[iii].intfdata[2*kk];
                  side     = SurfaceListFinal[iii].intfdata[2*kk+1];

                  ngbf31 = SurfaceListFinal[patchnum].ngbf1;
                  ngbf32 = SurfaceListFinal[patchnum].ngbf2;

                  switch(kk)
                  {
                      case 0: // left

                              for(jj=0;jj<ngbf22;jj++)
                              {
                                 index1 = ngbf21*jj;
                                 index2 = ngbf31*(jj+1)-1;
                                 SurfaceListFinal[iii].gbfnums[index1] = SurfaceListFinal[patchnum].gbfnums[index2];
                              }
                         break;

                      case 1: // right

                              for(jj=0;jj<ngbf22;jj++)
                              {
                                 index2 = ngbf21*jj;
                                 index1 = ngbf31*(jj+1)-1;
                                 SurfaceListFinal[iii].gbfnums[index1] = SurfaceListFinal[patchnum].gbfnums[index2];
                              }
                         break;

                      case 2: // bottom

                              for(jj=0;jj<ngbf21;jj++)
                              {
                                 index1 = jj;
                                 index2 = ngbf31*(ngbf32-1)+jj;
                                 SurfaceListFinal[iii].gbfnums[index1] = SurfaceListFinal[patchnum].gbfnums[index2];
                              }
                         break;

                      case 3: // top

                              for(jj=0;jj<ngbf21;jj++)
                              {
                                 index1 = ngbf21*(ngbf22-1)+jj;
                                 index2 = jj;
                                 SurfaceListFinal[iii].gbfnums[index1] = SurfaceListFinal[patchnum].gbfnums[index2];
                              }
                         break;

                  }
               }
            }
            for(ii=0;ii<SurfaceListFinal[iii].ngbf;ii++)
            {
               if( SurfaceListFinal[iii].gbfnums[ii] == -5555)
                  SurfaceListFinal[iii].gbfnums[ii] = count++;
            }
        }
     }

//  cout << '\t' << "   SSSSSS    " << SurfaceListFinal[0].ngbf1 << '\t' << SurfaceListFinal[0].ngbf2 << endl;

        int index1, index2, *U1, *U2;
        for(iii=0;iii<Npatch;iii++)
        {
            U1 = &(SurfaceListFinal[iii].gbfnums[0]);
            U2 = &(SurfaceListFinal[iii].toUfull[0]);
            for(ii=0;ii<SurfaceListFinal[iii].ngbf;ii++)
            {
               index1 = U1[ii]*ndf;
               index2 = ii*ndf;
               for(jj=0;jj<ndf;jj++)
                 U2[index2+jj] = index1 + jj;
            }
        }

/*
     cout << endl;
     cout << "       .... gbfnumbers for the final mesh ... " << endl;
     cout << endl;
     for(iii=0;iii<Npatch;iii++)
     {
         ngbf21 = SurfaceListFinal[iii].ngbf1;

         cout << "       patch... : " << (iii+1) << endl;
         cout << endl;
         count = 0;
         for(ii=0;ii<SurfaceListFinal[iii].ngbf1;ii++)
         {
            cout << '\t' ;
            for(jj=0;jj<SurfaceListFinal[iii].ngbf2;jj++)
            {
              count = ngbf21*jj;
              cout << SurfaceListFinal[iii].gbfnums[count+ii] << '\t' ;
            }
            cout << endl;
         }
         cout << endl;
         cout << endl;
     }

     cout << endl;
     cout << "       .... toUfull values for the final mesh ... " << endl;
     cout << endl;

     for(iii=0;iii<Npatch;iii++)
     {
         cout << "       patch... : " << (iii+1) << endl;
         cout << endl;
         cout << '\t' << SurfaceListFinal[iii].toUfull << endl;
     }
         cout << endl;
         cout << endl;

    // printData(2);
*/

        solnInit.zero();
        int count1 = -1;
        double *U3;
        int *U4;
        for(iii=0;iii<Npatch;iii++)
        {
            U4 = &(SurfaceListFinal[iii].toUfull[0]);
            U3 = &(SurfaceListFinal[iii].Uinit[0]);
            for(ii=0;ii<SurfaceListFinal[iii].Uinit.n;ii++)
            {
               if(U4[ii] > count1)
                 solnInit[U4[ii]] = U3[ii];
            }
           count1 += (SurfaceListFinal[iii-1].ngbf * ndf);
        }

  //ufull = Uinit;

        for(iii=0;iii<Npatch;iii++)
        {
           SurfaceResult[iii].toUfull = SurfaceListFinal[iii].toUfull;
           SurfaceResult[iii].gbfnums = SurfaceListFinal[iii].gbfnums;
        }


  return;

}






void IsogeometricFEM::readSurfaceFromFile(MyString &fileName, int geom, int patchnum)
{
   char fct[] = "IsogeometricFEM::readSurfaceFromFile";

   if( patchnum > SurfaceListFinal.n )
      prgError(1,fct," invalid patch number ");


    SurfaceResult[patchnum-1].readSurfaceFromFile(fileName);

  return;
}







void IsogeometricFEM::writeGeomToFile(MyString &fileName, int geom, int index, int patchnum)
{
  char fct[] = "IsogeometricFEM::writeSurfaceToFile";

   if( patchnum > Npatch )
      prgError(1,fct," invalid patch number ");


  switch(geom)
  {
      case 1:

             SurfaceListOriginal[patchnum-1].writeToFile(fileName, index);
             break;

      case 2:

             SurfaceListFinal[patchnum-1].writeToFile(fileName, index);
             break;

      case 3:

             SurfaceResult[patchnum-1].writeToFile(fileName, index);
             break;

      default:

             prgError(2,fct," invalid geom number ");
	      break;

  }


  return;
}





void IsogeometricFEM::copyElemInternalVariables()
{
  // copy internal variables intVar1 to intVar2

  int e, i, nivEl;
  nivEl = elem[0]->nivGP * elem[0]->nGP;
  double *intVar1, *intVar2;

  for(e=0; e<totnumel; e++)
  {
    intVar1 = elem[e]->intVar1;
    intVar2 = elem[e]->intVar2;

    for (i=0; i<nivEl; i++) intVar2[i] = intVar1[i];
  }

  intVar1 = NULL;
  intVar2 = NULL;

  return;
}






int IsogeometricFEM::calcAndAssyLoadVector(double fact, double dt)
{
  int temp;
  ForceVec.zero();

  for(int e=0;e<totnumel;e++)
  {
    //cout << "       elem... : " << (e+1) << endl;  //          cout << endl;     
    temp = elem[e]->calcLoadVector();
    if(temp != 0)
    {
      cout << " local element failure --> load vector calculation. Elem # : " << e << endl;
      return temp;
    }
    elem[e]->AssembleElementVector(firstIter, 1, &(ForceVec[0]), &(reac[0]), 0, 0 );
  }

  return 0;
}







void IsogeometricFEM::writeNodalData()
{
  char fct[] = "IsogeometricFEM::writeNodalData";

  int         i, j, dof, patchnum, type;
  double      val, fact, uval, vval, temp1, temp2, wval;
  char        tmp[20];
  MyString    tmpStr, fname;

  for (i=0; i<outdparam.n; i++)
  {
    patchnum = (int) outdparam[i][0]-1;
    type     = (int) outdparam[i][1];
    dof      = (int) outdparam[i][2];
    fact     = outdparam[i][3];
    uval     = outdparam[i][4];
    vval     = outdparam[i][5];
    
    //cout << patchnum << '\t' << type << '\t' << dof << '\t' << fact << '\t' << uval << '\t' << vval << endl;

    val  = 0.0;
    switch(type)
    {
      case  1: // u
      
            if(patchGrp[0].ndom == 2)
            {
               if(dof == 1)
               {
                 temp1 = SurfaceResult[patchnum].SurfacePoint(uval, vval).CalcEuclid().x;
                 temp2 = SurfaceListFinal[patchnum].SurfacePoint(uval, vval).CalcEuclid().x;
               }
               else
               {
                 temp1 = SurfaceResult[patchnum].SurfacePoint(uval, vval).CalcEuclid().y;
                 temp2 = SurfaceListFinal[patchnum].SurfacePoint(uval, vval).CalcEuclid().y;
               }
               val = temp1 - temp2;
            }

            if(patchGrp[0].ndom == 3)
            {
               wval = outdparam[i][6];
               
               if(dof == 1)
               {
                 temp1 = SolidResult[patchnum].SolidPoint(uval, vval, wval).CalcEuclid().x;
                 temp2 = SolidListFinal[patchnum].SolidPoint(uval, vval, wval).CalcEuclid().x;
               }
               else if(dof == 2)
               {
                 temp1 = SolidResult[patchnum].SolidPoint(uval, vval, wval).CalcEuclid().y;
                 temp2 = SolidListFinal[patchnum].SolidPoint(uval, vval, wval).CalcEuclid().y;
               }
               else
               {
                 temp1 = SolidResult[patchnum].SolidPoint(uval, vval, wval).CalcEuclid().z;
                 temp2 = SolidListFinal[patchnum].SolidPoint(uval, vval, wval).CalcEuclid().z;
               }

               val = temp1 - temp2;
            }

               break;

      case  2: // reac at a particular point

               {
                  int nlbf = SurfaceResult[0].nlbf;

                  VectorArray<double> N;
                  VectorArray<int> ind;

                  //NurbsShapeFunctions2DPost3(&SurfaceResult[0], uval, vval, N, ind);

                  switch(dof)
                  {
                     case 1:

                            for(int kk=0;kk<nlbf;kk++)
                              val += N[kk] * reac[ndf*ind[kk]+0] ;
                            break;

                     case 2:

                            for(int kk=0;kk<nlbf;kk++)
                              val += N[kk] * reac[ndf*ind[kk]+0] ;
                            break;

                     default: prgError(1,fct,"invalid 'dof' number!");

                  }
               }

               break;

      case  3: // reaction on a complete face/edge

               val = computeReactions(1,i);

               break;

      case  4: // moment

               val = computeReactions(2,i);

               break;

      default: prgError(1,fct,"unknown wrndType!");

               break;
    }
    val *= fact;
    sprintf(tmp," %12.5f",val);
    tmpStr.append(tmp);
  }

cout << " val " << val << endl;

  prgWriteToTFile(tmpStr);

  return;
}


void IsogeometricFEM::plotGaussPoints(int num, bool defFlg)
{

   for(int e=0;e<totnumel;e++)
     elem[e]->plotGaussPoints(num, defFlg);

return;
}



void IsogeometricFEM::checkContSurfs(int geom, int val1, int val2, int dir)
{

  cout << '\t' << geom << '\t' << val1 << '\t' << val2 << '\t' << dir << endl;

  switch(geom)
  {
     case 1:
         checkContSurfaces(&(SurfaceListFinal[val1]), &(SurfaceListFinal[val2]), dir);
         break;

     case 2:
         checkContSurfaces(&(SurfaceResult[val1]), &(SurfaceResult[val2]), dir);
         break;

     default:
         cout << '\t' << " invalid value of 'geom' number  " << endl;
  }


  return;
}


void IsogeometricFEM::printElemInvVars(int patchnum, int elenum)
{
   if(elenum > SurfaceListFinal[patchnum].nelem)
   {
     cout << '\t' << " Element number for the given patch is out of range " << endl;
     cout << endl;
     return;
   }

   int ee=0, i, n, iii;
   double *intVar1, *intVar2;

   for(iii=0;iii<patchnum;iii++)
     ee += SurfaceListFinal[iii].nelem;

   ee += elenum;

   cout << '\t' << " global element number : " << ee << endl;

   intVar1 = elem[ee]->intVar1;
   intVar2 = elem[ee]->intVar2;

   n = (elem[ee]->nivGP) * ( elem[ee]->nGP);
   for(i=0; i<n; i++)
     cout << '\t' << i << '\t' << intVar1[i] << '\t' << intVar2[i] << endl;
   cout << endl;

  return;
}


void IsogeometricFEM::GenerateConnectivityArraysConstraintVariables()
{
  if(ndm == 1)
    GenerateConnectivityArraysConstraintVariablesCurve();
  else if(ndm == 2)
    GenerateConnectivityArraysConstraintVariablesSurface();
  else
    GenerateConnectivityArraysConstraintVariablesSolid();
  
  return;
}



void IsogeometricFEM::GenerateConnectivityArraysConstraintVariablesCurve()
{
  return;
}


void IsogeometricFEM::GenerateConnectivityArraysConstraintVariablesSurface()
{
   cout << "     ISOGEOMETRICFEM: generating cinnectivity arrays for Constraint Variables ...\n\n";
   char fct[] = "IsogeometricFEM::GenerateConnectivityArraysConstraintVariables";

     int ii, jj, iii, kk;

     ntotgbf2 = 0;
     // compute total number of independent global basis funtions for the constraint variable
     for(iii=0;iii<Npatch;iii++)
     {
        ntotgbf2 += surfSecondVar[iii].ngbf;
     }


   ntoteqs2=0;
   for(iii=0;iii<surfSecondVar.n;iii++)
     surfSecondVar[iii].GenerateConnectivityArrays1(ntoteqs2);

   cout << '\t' << " ntotgbf2 & ntoteqs2  : " << ntotgbf2 << '\t' << ntoteqs2 << endl;

   // calculate LM Arrays
   for(iii=0;iii<SurfaceListFinal.n;iii++)
      surfSecondVar[iii].GenerateConnectivityArrays2();

//   for(iii=0;iii<SurfaceListFinal.n;iii++)
//      surfSecondVar[iii].printConnectivityArrays();




     int ngbf21, ngbf22, count;
     // calculate gbf numbers for patch #1
     ngbf21 = surfSecondVar[0].ngbf1;
     ngbf22 = surfSecondVar[0].ngbf2;
     count = 0;
     for(jj=0;jj<ngbf22;jj++)
     {
       for(kk=0;kk<ngbf21;kk++)
       {
          surfSecondVar[0].gbfnums[count] = count;
          count++;
       }
     }

     if(Npatch > 1)
     {
        count = surfSecondVar[0].ngbf;
        int patchnum, side, ngbf31, ngbf32, index1, index2;


        for(iii=1;iii<Npatch;iii++)
        {
            for(ii=0;ii<surfSecondVar[iii].ngbf;ii++)
            {
               if( surfSecondVar[iii].gbfnums[ii] == -5555)
                  surfSecondVar[iii].gbfnums[ii] = count++;
            }
        }
     }


        int index1, index2, *U1, *U2, ndof;
        for(iii=0;iii<Npatch;iii++)
        {
            U1 = &(surfSecondVar[iii].gbfnums[0]);
            U2 = &(surfSecondVar[iii].toUfull[0]);
            ndof = surfSecondVar[iii].ndof;

            index2 = 0;
            for(ii=0;ii<surfSecondVar[iii].ngbf;ii++)
            {
               index1 = U1[ii]*ndof;

               for(jj=0;jj<ndof;jj++)
               {
                 U2[index2] = index1 + jj;
                 index2++;
               }
            }
        }

/*
        for(iii=0;iii<Npatch;iii++)
        {
           cout << surfSecondVar[iii].gbfnums << endl;
           cout << surfSecondVar[iii].toUfull << endl;
        }
*/

  return;

}




void IsogeometricFEM::GenerateConnectivityArraysConstraintVariablesSolid()
{
   cout << "     ISOGEOMETRICFEM: generating cinnectivity arrays for Constraint Variables ...\n\n";
   char fct[] = "IsogeometricFEM::GenerateConnectivityArraysConstraintVariables";

     int ii, jj, iii, kk;

     ntotgbf2 = 0;
     // compute total number of independent global basis funtions for the constraint variable
     for(iii=0;iii<patch.n;iii++)
     {
        ntotgbf2 += NurbsBaseSecondVar[iii]->ngbf;
     }


   ntoteqs2=0;
   for(iii=0;iii<patch.n;iii++)
     NurbsBaseSecondVar[iii]->GenerateConnectivityArrays1(ntoteqs2);

   cout << '\t' << " ntotgbf2 & ntoteqs2  : " << ntotgbf2 << '\t' << ntoteqs2 << endl;

   // calculate LM Arrays
   for(iii=0;iii<patch.n;iii++)
      NurbsBaseSecondVar[iii]->GenerateConnectivityArrays2();
      
//   for(iii=0;iii<patch.n;iii++)
//      NurbsBaseSecondVar[iii]->printConnectivityArrays();


     int ngbf21, ngbf22, ngbf32, count;

     ngbf21 = solidSecondVar[0].ngbf1;
     ngbf22 = solidSecondVar[0].ngbf2;
     ngbf32 = solidSecondVar[0].ngbf3;
     count = 0;
     
     for(kk=0;kk<ngbf32;kk++)
     {
        for(jj=0;jj<ngbf22;jj++)
        {
          for(ii=0;ii<ngbf21;ii++)
          {
             NurbsBaseSecondVar[0]->gbfnums[count] = count;
             count++;
          }
        }
     }

     if(patch.n > 1)
     {
        count = NurbsBaseSecondVar[0]->ngbf;
        int patchnum, side, ngbf31, ngbf32, index1, index2;


        for(iii=1;iii<patch.n;iii++)
        {
            for(ii=0;ii<NurbsBaseSecondVar[iii]->ngbf;ii++)
            {
               if( NurbsBaseSecondVar[iii]->gbfnums[ii] == -5555)
                  NurbsBaseSecondVar[iii]->gbfnums[ii] = count++;
            }
        }
     }


        int index1, index2, *U1, *U2, ndof;
        for(iii=0;iii<patch.n;iii++)
        {
            U1 = &(NurbsBaseSecondVar[iii]->gbfnums[0]);
            U2 = &(NurbsBaseSecondVar[iii]->toUfull[0]);
            ndof = NurbsBaseSecondVar[iii]->ndof;

            index2 = 0;
            for(ii=0;ii<NurbsBaseSecondVar[iii]->ngbf;ii++)
            {
               index1 = U1[ii]*ndof;

               for(jj=0;jj<ndof;jj++)
               {
                 U2[index2] = index1 + jj;
                 index2++;
               }
            }
        }

/*
        for(iii=0;iii<patch.n;iii++)
        {
           cout << NurbsBaseSecondVar[iii]->gbfnums << endl;
           cout << NurbsBaseSecondVar[iii]->toUfull << endl;
        }
*/

  return;
}
