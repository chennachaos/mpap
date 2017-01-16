
#include <iostream>

#include "IsogeometricFEM.h"
#include "FunctionsProgram.h"
#include "DataBlockTemplate.h"
#include "PropertyTypeEnum.h"
#include "MathGeom.h"
#include "SolverMA41.h"
#include "SolverPARDISO.h"
#include "SolverPardisoEigen.h"
#include "NurbsShapeFns.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "NurbsUtilitiesSOLID.h"
#include "PlotVTK.h"

#include "NurbsElem1DAdvectionDiffusion.h"
#include "NurbsElem2DAdvectionDiffusion.h"
#include "NurbsElem2DStokes.h"
#include "NurbsElem2DNavierStokes3dof.h"
#include "NurbsElem2DNavierStokes4dof.h"
#include "NurbsElem1DElasticBar.h"
#include "NurbsElem1DEulerBeam.h"
#include "NurbsElem1DElasticBarLSFEM.h"
#include "NurbsElem2DStructSolid.h"
#include "NurbsElem2DStructSolidLSFEM2dof.h"
#include "NurbsElem2DStructSolidLSFEM3dof.h"
#include "NurbsElemKirchhoffPlate.h"
#include "NurbsElemMindlinPlate.h"
#include "NurbsElem2DStructFbarSolid.h"
#include "NurbsElem2DStructBbarSolid.h"
#include "NurbsElem2DStructMixed2fieldStabilised.h"
#include "NurbsElem3DStructMixed2fieldStabilised.h"
#include "NurbsElem2DStructMixed2field.h"
#include "NurbsElem2DStructMixed3field.h"
#include "NurbsElem3DStructSolid.h"
#include "NurbsElem3DStructMixed2field.h"
#include "NurbsElem2DHeatTransfer.h"
#include "NurbsElem2DTempCoupled4dof.h"

#include "util.h"

#include "DistFunctions3D.h"

extern Plot               plot;
extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern MpapTime           mpapTime;
extern ComputerTime       computerTime;
extern PlotVTK  plotvtk;



IsogeometricFEM::IsogeometricFEM(void)
{
  // initialise stuff here

  totIntVar = 0;

  RSYS = 0;

  TSTEP = -1;

  numnp = -1;

  ndf = -1;

  tol = -2.0;

  solver = NULL;
  elem = NULL;

  NurbsBaseOriginal  = NULL;
  NurbsBaseFinal     = NULL;
  NurbsBaseResult    = NULL;
  NurbsBaseSecondVar = NULL;
  
  localStiffnessError = 0;

  ctimFactSolvUpdt = 0.;
  ctimCalcStiffRes = 0.;

  globalFirstIter = true;
  comprMtxFlg = false;

  ntoteqs1 = ntoteqs2 = ntoteqs = 0;
  ntotgbf1 = ntotgbf2 = ntotgbf = 0;

  Npatch = 0;

  kRefineFlag = true;
  defUndefFlag = true;
  mixedSolverFlag = -1;
  LSFEMflag = false;

  filecount = 0;

  // add new type

  DomainType *isogeometricFEM = domain.newType(ISOGEOMETRICFEM,ROOTDOMAIN);

  if (isogeometricFEM == NULL) return;  // domain type already exists

  isogeometricFEM->key.addNew("patches",
                              "control points",
                              "patch group data",
                              "knot vectors",
                              "displacement boundary conditions",
                              "force boundary conditions",
                              "traction boundary conditions",
                              "traction bc data",
                              "element type",
                              "material",
                              "analysis options",
                              "control parameters",
                              "data output",
                              "displacement bc extra data",
                              "interface data",
                              "bcs for constraint variable",
                              "time steps");
}





IsogeometricFEM::~IsogeometricFEM(void)
{
   cout << "     ISOGEOMETRICFEM: destructor ...\n\n";

   if(solver != NULL) delete solver;     solver = NULL;

   if (elem != NULL)
   {
      for(int ii=0;ii<totnumel;ii++)
        delete elem[ii];

      delete [] elem;
      elem = NULL;
   }

   if(NurbsBaseOriginal != NULL)
   {
      delete [] NurbsBaseOriginal;
      NurbsBaseOriginal = NULL;
   }

   if(NurbsBaseFinal != NULL)
   {
      delete [] NurbsBaseFinal;
      NurbsBaseFinal = NULL;
   }

   if(NurbsBaseResult != NULL)
   {
      delete [] NurbsBaseResult;
      NurbsBaseResult = NULL;
   }

   if(NurbsBaseSecondVar != NULL)
   {
      delete [] NurbsBaseSecondVar;
      NurbsBaseSecondVar = NULL;
   }

   if(plotvtk.ActiveFlag)
   {
     plotvtk.clearWindow();
     plotvtk.set();
   }
   
   cout << "     ISOGEOMETRICFEM: destructor ...\n\n";

}





void IsogeometricFEM::readInputData(std::ifstream &Ifile, MyString &line)
{
  char fct[] = "IsogeometricFEM::readInputData";

  MyString tmpl, *word;

  char tmp[30];

  int nw, i, j, k, n, nn;

  double xx[10];

  Vector<double> dTmp, dTmp2, outdfactTmp;
  Vector<int>    iTmp, iTmp2, lTmp, outdtypeTmp, outddofTmp;
  MyStringList   sTmp;
  List<Vector<int> > lviTmp, lviTmp2, lviTmp3, lviTmp4, lviTmp5;
  List<Vector<double> > lvdTmp, lvdTmp2, lvdTmp3, lvdTmp4;

  DataBlockTemplate t1, t2;

  switch (domain[ISOGEOMETRICFEM].key.whichBegins(line))
  {
    case  0: cout << "     ISOGEOMETRICFEM: reading patches ...\n\n";

             if (Npatch > 0) prgError(1,fct,"multiple definition of 'patches'");

             if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
               prgError(1,fct,"invalid node number or invalid keyword!");

             patch.setDim(lviTmp.n);

             Npatch = lviTmp.n;

             for (i=0; i<lviTmp.n; i++)
             {
               if (lviTmp[i][0] != i+1) prgError(2,fct,"counter error in 'patches'!");

               patch[i].setDim(lviTmp[i].n - 1);

               for(j=1; j<lviTmp[i].n; j++)
               {
                 if(lviTmp[i][j] < 0)
                   prgError(3,fct,"input data in patches can not be negative!");
                 patch[i][j-1] = lviTmp[i][j];
               }
             }
             break;

    case  1: cout << "     ISOGEOMETRICFEM: reading control points ...\n\n";

             if (numnp > 0) prgError(1,fct,"duplicate definition of 'control points'!");

             if (ndf < 1) prgError(1,fct,"'dimensions' must precede 'control points'!");

	      sprintf(tmp,"123 %df",ndm+2);

             if (!line.copyAfter('|',tmpl)) tmpl.free().append(tmp);

             t1.initialise(tmpl);
	      t2.initialise(tmp);
             t1.expandToMatch(t2);

             if (!t1.readBlock(Ifile,line,iTmp,dTmp,sTmp,lTmp))
                prgError(2,fct,"data error in 'control points'!");

             numnp = dTmp.dim() / (ndm+2);

             x.setDim(numnp,ndm+1,true);

             for (i=0; i<numnp; i++)
             {
               for (j=0; j<=ndm; j++)
                 x.x[i*(ndm+1)+j] = dTmp[i*(ndm+2)+j];
             }

	     break;

    case  2: cout << "     ISOGEOMETRICFEM: reading patch group data ...\n\n";

             if (Npatch < 1) prgError(1,fct,"'patches' must precede 'patches data'!");

             if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp2))
               prgError(1,fct,"invalid input in Patch Groups!");

             numPatchGrp = lviTmp2.n;

             patchgrpdata.setDim(lviTmp2.n);

             for(i=0; i<lviTmp2.n; i++)
             {
               if(lviTmp2[i].n < (4 * patch[i][1] + 3))
                 prgError(2,fct,"data error in 'patch group data'!");

               if (lviTmp2[i][0] != i+1) prgError(2,fct,"counter error in 'patch group data'!");

               patchgrpdata[i].setDim(lviTmp2[i].n - 1);

               for (j=1; j<lviTmp2[i].n; j++)
               {
                 if(lviTmp2[i][j] < 0)
                   prgError(3,fct,"input data in patch group data can not be negative!");
                 patchgrpdata[i][j-1] = lviTmp2[i][j];
               }
             }

             break;

    case  3: cout << "     ISOGEOMETRICFEM: reading knot vectors ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
               prgError(10,fct,"invalid knot vector or invalid keyword!");

             kv.setDim(lvdTmp.n);

             for(i=0;i<kv.n;i++)
             {
                 if (lvdTmp[i][0] != i+1) prgError(1,fct,"counter error in 'knot vectors'!");

                 kv[i].setDim(lvdTmp[i].n - 1);

                 for (j=1; j<lvdTmp[i].n; j++)
                 {
                    if(lvdTmp[i][j] < 0)
                      prgError(3,fct,"knot values in knot vectors can not be negative!");
                    kv[i][j-1] = lvdTmp[i][j];
                 }
             }

             break;


    case  4: cout << "     ISOGEOMETRICFEM: reading displacement boundary conditions ...\n\n";

             if (ndf < 1) prgError(1,fct,"'dimensions' must precede 'displacement boundary conditions'!");

             if (dispbc.n > 0) prgError(1,fct,"multiple definition of 'displacement boundary conditions'!");

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp2))
               prgError(2,fct,"invalid number or invalid keyword!");

             dispbc.setDim(lvdTmp2.n);

             for(i=0; i<lvdTmp2.n; i++)
             {
                 dispbc[i] = lvdTmp2[i];
             }

             /*
             for(i=0; i<lvdTmp2.n; i++)
             {
                 dispbc[i].setDim(4);

                 for(j=0; j<lvdTmp2[i].n; j++)
                 {
                    dispbc[i][j] = lvdTmp2[i][j];
                 }
             }
             */
             break;

    case  5: cout << "     ISOGEOMETRICFEM: reading force boundary conditions ...\n\n";

             if (ndf < 1) prgError(1,fct,"'dimensions' must precede ' force boundary conditions'!");

             if (forcebc.n > 0) prgError(1,fct,"multiple definition of ' force boundary conditions'!");


             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp3))
               prgError(2,fct,"invalid number or invalid keyword!");

             forcebc.setDim(lvdTmp3.n);

             for(i=0; i<lvdTmp3.n; i++)
             {
                 forcebc[i].setDim(4);

                 for(j=0; j<lvdTmp3[i].n; j++)
                 {
                    forcebc[i][j] = lvdTmp3[i][j];
                 }
             }

             break;

    case  6: cout << "     ISOGEOMETRICFEM: reading traction boundary conditions ...\n\n";

             if (ndf < 1) prgError(1,fct,"'dimensions' must precede ' traction boundary conditions'!");

             if (tracbc.n > 0) prgError(1,fct,"multiple definition of 'traction boundary conditions'!");

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp4))
               prgError(2,fct,"invalid number or invalid keyword!");

             tracbc.setDim(lvdTmp4.n);

             for(i=0; i<lvdTmp4.n; i++)
             {
               tracbc[i] = lvdTmp4[i];
               //cout << tracbc[i] << endl;
             }
             
             /*
             for(i=0; i<lvdTmp4.n; i++)
             {
               tracbc[i].setDim(4);

               for(j=0; j<lvdTmp4[i].n; j++)
                  tracbc[i][j] = lvdTmp4[i][j];
               //cout << tracbc[i] << endl;
             }
             */
   

             break;


    case  7: cout << "     ISOGEOMETRICFEM: reading traction bcs extra data ...\n\n";

             if (ndf < 1) prgError(1,fct,"'dimensions' must precede ' traction boundary conditions #2'!");

             if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp3))
               prgError(2,fct,"invalid number or invalid keyword!");

             tracbc2.setDim(lviTmp3.n);

             for (i=0; i<lviTmp3.n; i++)
             {
               tracbc2[i].setDim(4);

               for(j=0; j<lviTmp3[i].n; j++)
               {
                 if(lviTmp3[i][j] < 0)
                   prgError(3,fct,"input data in traction bcs extra data can not be negative!");
                 tracbc2[i][j] = lviTmp3[i][j];
               }
             }

             break;



    case  8: cout << "     ISOGEOMETRICFEM: reading element type ...\n\n";

             patchElemProp.add(new PropertyItem(ELEMENTTYPE));

             patchElemProp[patchElemProp.n-1].readInputData(Ifile,line,"input error in 'element type'!");

             break;

    case  9: cout << "     ISOGEOMETRICFEM: reading material ...\n\n";

             patchMatlProp.add(new PropertyItem(MATERIAL));

	      patchMatlProp[patchMatlProp.n-1].readInputData(Ifile,line,"input error in 'material'!");

	      break;

    case  10: cout << "     ISOGEOMETRICFEM: analysis options ...\n\n";

             if (analopts.n > 0)
                prgError(1,fct,"'analysis options' has already been read!");

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp2))
                prgError(2,fct,"invalid analysis options!");

             analopts.setDim(lvdTmp2[0].n);

             for (i=0; i<analopts.n; i++)
             {
                if( ((int) lvdTmp2[0][i]) < 0)
                   prgError(3,fct, " analysis options can not be negative!");

                analopts[i] = lvdTmp2[0][i];
             }

             break;

/*
    case  11: cout << "     ISOGEOMETRICFEM: reading control ...\n\n";

             if (tol > -1)         prgError(1,fct,"'control' has already been read!");

	     line.getNextLine(Ifile);

             nw = line.split(&word);

           //  if (nw < 3)                    prgError(1,fct,"input error in 'control'!");

             if (!word[0].toDbl(&tol))      prgError(2,fct,"input error in 'control'!");


             if (!word[1].toInt(&lumpType)) prgError(3,fct,"input error in 'control'!");
	     if (lumpType < 0) lumpType = 0;

	     if (!word[2].toInt(&tis))      prgError(4,fct,"input error in 'control'!");
	     if (tis < 0) tis = 0;

             for (i=0; i<10; i++) td[i] = 0.;

             nn = nw; if (nn > 13) nn = 13;

             for (i=3; i<nn; i++)
	       if (!word[i].toDbl(&(td[i-3]))) prgError(5,fct,"input error in 'control'!");


             for (i=0; i<nw; i++) word[i].free(); delete [] word;

       	     line.getNextLine(Ifile);

             break;
*/
    case  11: cout << "     IsogeometricFEM: reading 'control parameters' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'control parameters'!");

            if( lvdTmp[0].n < 3)
              cerr <<  " Error in (( IsogeometricFEM: reading 'control parameters' )) " << endl;

            tol      = lvdTmp[0][0];
            tis      = (int) lvdTmp[0][1];
            rhoInfty = lvdTmp[0][2];

            cout << tol << '\t' << tis << '\t' << rhoInfty << endl;

            break;


    case 12: cout << "     ISOGEOMETRICFEM: reading data output...\n\n";

             if (outdparam.n > 0) prgError(1,fct,"multiple definition of ' reading data output '!");

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp4))
               prgError(2,fct,"invalid number or invalid keyword!");

             outdparam.setDim(lvdTmp4.n);

             for(i=0; i<lvdTmp4.n; i++)
             {
               outdparam[i] = lvdTmp4[i];
               cout << outdparam[i] << endl;
             }

             break;

    case  13: cout << "     ISOGEOMETRICFEM: reading displacement bcs extra data ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp3))
               prgError(1,fct,"invalid number or invalid keyword!");

             dispbc2.setDim(lvdTmp3.n);

             for (i=0; i<lvdTmp3.n; i++)
             {
                dispbc2[i].setDim(6);

                for(j=0; j<lvdTmp3[i].n; j++)
                {
                   if(j != 5)
                   {
                     if(lvdTmp3[i][j] < 0)
                       prgError(2,fct,"input data in displacement bcs extra data can not be negative!");
                   }
                   dispbc2[i][j] = lvdTmp3[i][j];
                }
             }

             //cout << '\t' << " dispbc2 " << dispbc2[0] << endl;

             break;

    case  14: cout << "     ISOGEOMETRICFEM: reading interface data ...\n\n";

             if (intfdata.n > 0) prgError(1,fct," multiple definition of 'interface data' ");

             if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp5))
               prgError(2,fct,"invalid number or invalid keyword!");

             intfdata.setDim(lviTmp5.n);

             for (i=0; i<lviTmp5.n; i++)
             {
                 if (lviTmp5[i].n < 4) prgError(3,fct,"data error in 'interface data'!");

                 intfdata[i].setDim(4);

                 for(j=0; j<lviTmp5[i].n; j++)
                 {
                    if(lviTmp5[i][j] < 0)
                      prgError(4,fct,"input data in patches can not be negative!");
                    intfdata[i][j] = lviTmp5[i][j];
                 }
             }

        //     cout << '\t' << " intfdata " << intfdata[0] << endl;

             break;

    case  15: cout << "     ISOGEOMETRICFEM: reading BCs for constraint variable ...\n\n";

             if (ndf < 1) prgError(1,fct,"'dimensions' must precede 'bcs for constraint variable'!");

             if (ibc4.n > 0) prgError(1,fct,"multiple definition of 'bcs for constraint variable'!");

             if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp5))
               prgError(2,fct,"invalid number or invalid keyword!");

             ibc4.setDim(lviTmp5.n);

             for(i=0; i<lviTmp5.n; i++)
             {
                ibc4[i].setDim(2);

                for(j=0;j<2;j++)
                  ibc4[i][j] = lviTmp5[i][j];
             }

             break;

    case -1: // go and inherit from DOMAIN

             this->Domain::readInputData(Ifile,line);

             break;
  }

  return;
}








void IsogeometricFEM::prepareInputData(void)
{

  char fct[] = "IsogeometricFEM::prepareInputData";

  // go and inherit from ancestors

  Domain::prepareInputData();

  cout << "     ISOGEOMETRICFEM: prepare input data ...\n\n";

  int i, j, ndomT, index=0, tt1, tt2;

   // ==================================================
   //
   // Check the  consistency of input data
   //
   // ==================================================

   checkInputData1();

  // Prepare patchElemProp data. Assign suitable id(in the database) based on the element name

  if(patchElemProp.n > 0)
    preparePatchElemProp();

  // Prepare patchMatlProp data. Assign suitable id(in the database) based on the material name

  if(patchMatlProp.n > 0)
    preparePatchMatlProp();

   // generate patch groups

   if(patchGrp.n < 1)
   {
       if (Npatch < 1) prgError(1,fct,"define at least one patch group!");

       for (i=0; i<Npatch; i++)
         patchGrp.add(new PatchGroup);

       for (i=0; i<Npatch; i++)
       {
          index = patch[i][0] - 1;

          ndomT  = patch[i][1];

          //assign the data for member variables of PatchGroup class

          patchGrp[i].ndom = ndomT;
          patchGrp[i].reparmid = patch[i][2];

          //set dimensions for polydegree, kvector and subDiv member variable of PatchGroup Class

          patchGrp[i].polydegree.setDim(ndomT);
          patchGrp[i].kvindex.setDim(ndomT);
          patchGrp[i].subDiv.setDim(ndomT);
          patchGrp[i].elevdegree.setDim(ndomT);

          patchGrp[i].ctrlpoints.setDim(patch[i].n - 3);

          for(int ii=0;ii<ndomT;ii++)
          {
             patchGrp[i].polydegree[ii] = patchgrpdata[index][ii];
             patchGrp[i].kvindex[ii] = patchgrpdata[index][ndomT + ii];
             patchGrp[i].subDiv[ii] = patchgrpdata[index][2*ndomT + ii];
             patchGrp[i].elevdegree[ii] = patchgrpdata[index][3*ndomT + ii];
          }

          patchGrp[i].eltype = patchElemProp[patchgrpdata[index][4*ndomT]-1].id;

          for(int ii=0;ii<patchGrp[i].ctrlpoints.n;ii++)
            patchGrp[i].ctrlpoints[ii] = patch[i][3+ii];
       }
   }
   // generation of patch groups is over

   // Check the  consistency of element type and dof
   checkInputData2();

   // ==================================================
   //
   // compute the number of curves and surfaces in the input file
   //
   // ==================================================

   countcurve = 0;
   countsurf  = 0;
   countsolid = 0;

   for(i=0;i<patchGrp.n;i++)
   {
      if(patchGrp[i].ndom == 1)
        countcurve++;
      else if(patchGrp[i].ndom == 2)
        countsurf++;
      else
        countsolid++;
   }


   // assign dimensions to curve and surface lists
   if(countcurve > 0)
   {
      CurveListOriginal.setDim(countcurve);
      CurveListFinal.setDim(countcurve);
      CurveResult.setDim(countcurve);
   }
   if(countsurf > 0)
   {
      SurfaceListOriginal.setDim(countsurf);
      SurfaceListFinal.setDim(countsurf);
      SurfaceResult.setDim(countsurf);
   }
   if(countsolid > 0)
   {
      SolidListOriginal.setDim(countsolid);
      SolidListFinal.setDim(countsolid);
      SolidResult.setDim(countsolid);
   }

   countcurve = 0;
   countsurf  = 0;
   countsolid = 0;

   totnumel = 0;



   // ==================================================
   //
   // Generate NURBS Class Objects based on the input data in the patch groups
   //
   // ==================================================

  
   int   ii, jj, kk, ind1, ind2, ind3;
   double xx=0.0, yy=0.0, zz=0.0, wt=0.0;

  for(i=0;i<patchGrp.n;i++) // for loop A
  {
      if(patchGrp[i].ndom == 1)
      {
          KNOTVECTOR U; U = kv[patchGrp[i].kvindex[0] - 1];

          DEGREE p = patchGrp[i].polydegree[0];

          int mU = U.n -1;
          int nU = mU - p -1;

          CPOLYGON Pw;
          Pw.setDim(patchGrp[i].ctrlpoints.n);

          int ind=0;
          double xx=0.0, yy=0.0, zz=0.0, wt=0.0;

          // Generate the coordinates of the control points
          for(ii=0;ii<Pw.n;ii++)
          {
             ind1 = patchGrp[i].ctrlpoints[ii] - 1;
             ind2 = (ndm+1)*ind1;
             
             xx = x.x[ind2];
             yy = x.x[ind2 + 1];

             if(ndm == 3) // 3D Rational Curve
             {
               zz = x.x[ind2 + 2];
               wt = x.x[ind2 + 3];

               Pw[ii].x = xx * wt ;
               Pw[ii].y = yy * wt ;
               Pw[ii].z = zz * wt ;
               Pw[ii].w = wt ;
             }
             else // 1D/2D Rational Curve
             {
               wt = x.x[ind2 + 2];

               Pw[ii].x = xx * wt ;
               Pw[ii].y = yy * wt ;
               Pw[ii].z = 0.0 ;
               Pw[ii].w = wt ;
             }
         }

         CurveListOriginal[countcurve].Pw = Pw;
         CurveListOriginal[countcurve].U = U;
         CurveListOriginal[countcurve].p = p;

         CurveListOriginal[countcurve].ndof = ndf;

         CurveListOriginal[countcurve].initializeBCdata();

        NurbsCURVE curve_temp(Pw, U, p);

        kRefineFlag = true;

        if( patchgrpdata[i].n > (4 * patchGrp[i].ndom + 2) )
        {
           int temp = patchgrpdata[i][4 * patchGrp[i].ndom + 2];

           kRefineFlag = !(temp == 1);
        }

        cout << " kRefineFlag  " <<  kRefineFlag  << endl;

        if(kRefineFlag)
        {
           // Do the k-refinement
           // First Degree Elevate the Curve

           if(patchGrp[i].reparmid == 0)
           {
              if(patchGrp[i].elevdegree[0] > CurveListOriginal[countcurve].p)
              {
                 int tt = patchGrp[i].elevdegree[0] - CurveListOriginal[countcurve].p;
                 DegreeElevateCurve(&CurveListOriginal[countcurve], tt, &curve_temp);
              }

              // Then refine the knot vector of the degree elevated curve
              if(patchGrp[i].subDiv[0] > 0)
              {
                 KNOTVECTOR X;

                 GenKnotVecForRefining(curve_temp.U, patchGrp[i].subDiv[0], X);

                 RefineCurveKnotVector(&curve_temp, X, &CurveListFinal[countcurve]);
              }
              else
              {
                 CurveListFinal[countcurve].Pw = curve_temp.Pw;
                 CurveListFinal[countcurve].U = curve_temp.U;
                 CurveListFinal[countcurve].p = curve_temp.p;
              }
           }
           else if(patchGrp[i].reparmid == 1)
           {
              RepamCurveNelem(&CurveListOriginal[countcurve], patchGrp[i].elevdegree[0], patchGrp[i].subDiv[0], &CurveListFinal[countcurve]);
           }
           else if(patchGrp[i].reparmid == 2)
           {
              RepamCurvenCPs(&CurveListOriginal[countcurve], patchGrp[i].elevdegree[0], patchGrp[i].subDiv[0], &CurveListFinal[countcurve]);
           }
        }
        else
        {
              if(patchGrp[i].subDiv[0] > 0)
              {
                 KNOTVECTOR X;

                 GenKnotVecForRefining(curve_temp.U, patchGrp[i].subDiv[0], X);
                 
                 cout << X << endl;

                 RefineCurveKnotVector(&CurveListOriginal[countcurve], X, &curve_temp);
              }

              if(patchGrp[i].elevdegree[0] > CurveListOriginal[countcurve].p)
              {
                 int tt = patchGrp[i].elevdegree[0] - CurveListOriginal[countcurve].p;
                 DegreeElevateCurve(&curve_temp, tt, &CurveListFinal[countcurve]);
              }
              else
              {
                 CurveListFinal[countcurve].Pw = curve_temp.Pw;
                 CurveListFinal[countcurve].U = curve_temp.U;
                 CurveListFinal[countcurve].p = curve_temp.p;
              }
        }

         // initialize the data members of curve

         CurveListFinal[countcurve].ndof = ndf;
         CurveListFinal[countcurve].initializeBCdata();

         CurveListFinal[countcurve].ElemProp = patchElemProp[patchgrpdata[index][4*(patchGrp[i].ndom)]-1];

         if(patchMatlProp.n > 0)
           CurveListFinal[countcurve].MatlProp = patchMatlProp[patchgrpdata[index][4*(patchGrp[i].ndom)+1]-1];

         CurveResult[countcurve].Pw =  CurveListFinal[countcurve].Pw;
         CurveResult[countcurve].U  =  CurveListFinal[countcurve].U;
         CurveResult[countcurve].p  =  CurveListFinal[countcurve].p;

         CurveResult[countcurve].ndof = ndf;
         CurveResult[countcurve].initializeBCdata();

         countcurve++;

    }  // end of if condition for curves

    else if(patchGrp[i].ndom == 2)
    {
        KNOTVECTOR U; U = kv[patchGrp[i].kvindex[0] - 1];
        KNOTVECTOR V; V = kv[patchGrp[i].kvindex[1] - 1];

        DEGREE p = patchGrp[i].polydegree[0];
        DEGREE q = patchGrp[i].polydegree[1];

        int mU = U.n -1;
        int mV = V.n -1;
        int nU = mU - p;
        int nV = mV - q;
        bool  closed1=false, closed2=false;
        
        closed1 = checkClosedness(2, 1, nU, nV, 0, &(patchGrp[i].ctrlpoints[0]));
        //cout << endl;        cout << endl;
        closed2 = checkClosedness(2, 2, nU, nV, 0, &(patchGrp[i].ctrlpoints[0]));
        
        //cout << patchGrp[i].ctrlpoints << endl;
        cout <<  " closed1 & closed2 " << closed1 << '\t' << closed2 << endl;
        //cout << nU << '\t' << nV << endl;

        CNET Pw;
        Pw.setDim(nU);
        for(ii=0;ii<nU;ii++)
          Pw[ii].setDim(nV);
          
        ind3=0;
        for(ii=0;ii<nU;ii++)
        {
            for(jj=0;jj<nV;jj++)
            {
               ind1 = patchGrp[i].ctrlpoints[ind3++] - 1; 
               ind2 = (ndm+1)*ind1;
               
               xx = x.x[ind2];
               yy = x.x[ind2 + 1];

               if(ndm == 3) // 3D Rational object
               {
                  zz = x.x[ind2 + 2];
                  wt = x.x[ind2 + 3];

                  Pw[ii][jj].x = xx * wt ;
                  Pw[ii][jj].y = yy * wt ;
                  Pw[ii][jj].z = zz * wt ;
                  Pw[ii][jj].w = wt ;
               }
               else // 2D Rational object
               {
                  wt = x.x[ind2 + 2];

                  //if(wt < 0.9)
                    //wt = 0.5*sqrt(2.0);

                  //cout << " wt = " <<  wt << endl;

                  Pw[ii][jj].x = xx * wt ;
                  Pw[ii][jj].y = yy * wt ;
                  Pw[ii][jj].z = 0.0 ;
                  Pw[ii][jj].w = wt ;
               }
            }
        }

        SurfaceListOriginal[countsurf].Pw = Pw;
        SurfaceListOriginal[countsurf].U  = U;
        SurfaceListOriginal[countsurf].V  = V;
        SurfaceListOriginal[countsurf].p  = p;
        SurfaceListOriginal[countsurf].q  = q;

        SurfaceListOriginal[countsurf].ndof = ndf;

        SurfaceListOriginal[countsurf].initializeBCdata();
        
        SurfaceListOriginal[countsurf].closed[0] = closed1;
        SurfaceListOriginal[countsurf].closed[1] = closed2;

        kRefineFlag = true;

        if( patchgrpdata[i].n > (4 * patchGrp[i].ndom + 2) )
        {
           int temp = patchgrpdata[i][4 * patchGrp[i].ndom + 2];
           
           cout << " temp  " << temp << endl;

           kRefineFlag = !(temp == 1);
        }

        cout << " kRefineFlag  " <<  kRefineFlag  << '\t' << patchGrp[i].reparmid << endl;

        if(kRefineFlag)
        {
            cout << " performing k-refinement " << endl;
            //===============================================
            //
            // Do the k-refinement
            // 1.) Degree Elevate first
            // 2.) Refine the degree elevated surface
            //
            //===============================================

            //===============================================
            // 1.) Degree Elevate the Surface
            //===============================================


            NurbsSURFACE surface_temp1(Pw, U, V, p, q);

            tt1 = patchGrp[i].elevdegree[0] - p;
            tt2 = patchGrp[i].elevdegree[1] - q;


            NurbsSURFACE surface_tempD(Pw, U, V, p, q);

            //degree elevate the surface in U direction
            DegreeElevateSurf1D(&SurfaceListOriginal[countsurf], tt1, 1, &surface_tempD);

            //degree elevate the surface in V direction
            DegreeElevateSurf1D(&surface_tempD, tt2, 2, &surface_temp1);
            
            cout << surface_temp1.U << endl;
            cout << surface_temp1.V << endl;
            //cout << surface_temp1.Pw << endl;

            //=============================================================
            // 2.) Refine the knot vector of the degree elevated surface
            //=============================================================


            //   reparmid = 0 ---> no reparameterization in any direction
            //            = 1 ---> reparameterization only in the 1st direction
            //            = 2 ---> reparameterization only in the 2nd direction
            //            = 3 ---> reparameterization in both the directions
            //            = 5 ---> reparameterize elements in both the directions
            //            = 6 ---> reparameterize CPs in both the directions
            //

            //    cout << '\t' <<  " patch # : " << (i+1) << '\t' << patchGrp[i].reparmid << endl;

            if(patchGrp[i].reparmid < 5)
            {
               if(patchGrp[i].subDiv[0] > 0 && patchGrp[i].subDiv[1] == 0)
               {
                  KNOTVECTOR XU;

                  if(patchGrp[i].reparmid == 1 || patchGrp[i].reparmid == 3)
                     FindKnotVector(&SurfaceListOriginal[countsurf], 1, patchGrp[i].subDiv[0], XU);
                  else
                     GenKnotVecForRefining(surface_temp1.U, patchGrp[i].subDiv[0], XU);

                  RefineSurfKnotVector1D(&surface_temp1, XU, 1, &SurfaceListFinal[countsurf]);
               }
               else if(patchGrp[i].subDiv[0] == 0 && patchGrp[i].subDiv[1] > 0)
               {
                   KNOTVECTOR XV;

                   if(patchGrp[i].reparmid == 2 || patchGrp[i].reparmid == 3)
                      FindKnotVector(&SurfaceListOriginal[countsurf], 2, patchGrp[i].subDiv[1], XV);
                   else
                      GenKnotVecForRefining(surface_temp1.V, patchGrp[i].subDiv[1], XV);

                   RefineSurfKnotVector1D(&surface_temp1, XV, 2, &SurfaceListFinal[countsurf]);
               }
               else if(patchGrp[i].subDiv[0] > 0 && patchGrp[i].subDiv[1] > 0)
               {
                   KNOTVECTOR XU, XV;

                   switch(patchGrp[i].reparmid)
                   {
                      case  0:
                              GenKnotVecForRefining(surface_temp1.U, patchGrp[i].subDiv[0], XU);
                              GenKnotVecForRefining(surface_temp1.V, patchGrp[i].subDiv[1], XV);
                              break;

                      case  1:
                              FindKnotVector(&SurfaceListOriginal[countsurf], 1, patchGrp[i].subDiv[0], XU);
                              GenKnotVecForRefining(surface_temp1.V, patchGrp[i].subDiv[1], XV);
                              break;

                      case  2:
                              GenKnotVecForRefining(surface_temp1.U, patchGrp[i].subDiv[0], XU);
                              FindKnotVector(&SurfaceListOriginal[countsurf], 2, patchGrp[i].subDiv[1], XV);
                              break;

                      case  3:
                              FindKnotVector(&SurfaceListOriginal[countsurf], 1, patchGrp[i].subDiv[0], XU);
                              FindKnotVector(&SurfaceListOriginal[countsurf], 2, patchGrp[i].subDiv[1], XV);
                              break;
                  }
                  cout << XU << endl;
                  cout << XV << endl;

                  NurbsSURFACE surface_temp4(Pw, U, V, p, q);

                  //refine the knot vector in U direction
                  RefineSurfKnotVector1D(&surface_temp1, XU, 1, &surface_temp4);


                  //refine the knot vector in V direction
                  RefineSurfKnotVector1D(&surface_temp4, XV, 2, &SurfaceListFinal[countsurf]);
               }
               else
               {
                  SurfaceListFinal[countsurf].Pw = surface_temp1.Pw;
                  SurfaceListFinal[countsurf].U  = surface_temp1.U;
                  SurfaceListFinal[countsurf].V  = surface_temp1.V;
                  SurfaceListFinal[countsurf].p  = surface_temp1.p;
                  SurfaceListFinal[countsurf].q  = surface_temp1.q;
               }
               SurfaceListFinal[countsurf].ndof = ndf;
               //SurfaceListFinal[countsurf].initializeBCdata();
               //SurfaceListFinal[countsurf].computeNET();
            }
            else if(patchGrp[i].reparmid == 5)
            {
               //ProcessForceBoundaryConditions();
               SurfaceListFinal[countsurf].ndof = ndf;
               cout << " qwertyqwerty " << endl;
               RepamSurf2DNelem(&SurfaceListOriginal[countsurf], patchGrp[i].elevdegree[0], patchGrp[i].elevdegree[1], patchGrp[i].subDiv[0], patchGrp[i].subDiv[1], &SurfaceListFinal[countsurf]);
            }
            else if(patchGrp[i].reparmid == 6)
            {
               SurfaceListFinal[countsurf].ndof = ndf;
               RepamSurf2DnCPs(&SurfaceListOriginal[countsurf], patchGrp[i].elevdegree[0], patchGrp[i].elevdegree[1], patchGrp[i].subDiv[0], patchGrp[i].subDiv[1], &SurfaceListFinal[countsurf]);
            }
            else //discretization based on specified number of elements //equal element length in parametric space
            {
              cout << patchGrp[i].subDiv[0] << '\t' << patchGrp[i].subDiv[1] << endl;
            
               if(patchGrp[i].subDiv[0] > 0 && patchGrp[i].subDiv[1] == 0)
               {
                  KNOTVECTOR XU;
                  
                  double  incr = 1.0/patchGrp[i].subDiv[0];
                  
                  create_vector(incr, 1.0-incr, incr, XU);
                  
                  cout << XU << endl;

                  RefineSurfKnotVector1D(&surface_temp1, XU, 1, &SurfaceListFinal[countsurf]);
               }
               else if(patchGrp[i].subDiv[0] == 0 && patchGrp[i].subDiv[1] > 0)
               {
                   KNOTVECTOR XV;

                  double  incr = 1.0/patchGrp[i].subDiv[1];

                  create_vector(incr, 1.0-incr, incr, XV);
                  
                  cout << XV << endl;

                  RefineSurfKnotVector1D(&surface_temp1, XV, 2, &SurfaceListFinal[countsurf]);
               }
               else if(patchGrp[i].subDiv[0] > 0 && patchGrp[i].subDiv[1] > 0)
               {
                   KNOTVECTOR XU, XV, XU1, XV1;

                  double  incr;
                  
                  incr = 1.0/patchGrp[i].subDiv[0];
                  
                  create_vector(incr, 1.0001-incr, incr, XU1);

                  incr = 1.0/patchGrp[i].subDiv[1];
                  
                  create_vector(incr, 1.0001-incr, incr, XV1);

                  cout << XU1 << endl;
                  cout << XV1 << endl;

                  XU = XU1;
                  XV = XV1;

                  /*
                  for(ii=0;ii<XU1.n;ii++)
                  {
                    if(XU1[ii] < 0.5)
                      XU[ii] = 2.0*pow(XU1[ii],2.0);
                    else
                      XU[ii] = 1.0 - 2.0*pow((1.0-XU1[ii]),2.0);
                  }

                  XV = XU;

                  cout << XU << endl;
                  cout << XV << endl;
                  */

                  NurbsSURFACE surface_temp4(Pw, U, V, p, q);

                  //refine the knot vector in U direction
                  RefineSurfKnotVector1D(&surface_temp1, XU, 1, &surface_temp4);

                  //refine the knot vector in V direction
                  RefineSurfKnotVector1D(&surface_temp4, XV, 2, &SurfaceListFinal[countsurf]);
               }
               else
               {
                  SurfaceListFinal[countsurf].Pw = surface_temp1.Pw;
                  SurfaceListFinal[countsurf].U  = surface_temp1.U;
                  SurfaceListFinal[countsurf].V  = surface_temp1.V;
                  SurfaceListFinal[countsurf].p  = surface_temp1.p;
                  SurfaceListFinal[countsurf].q  = surface_temp1.q;
               }
               SurfaceListFinal[countsurf].ndof = ndf;
            }

        }
        else
        {

            //===============================================
            //
            // Do the refinement opposite to k-refinement
            // 1.) Refine the surface first
            // 1.) Then degree elevate refined surface
            //
            //===============================================


            //=============================================================
            // 1.) Refine the surface
            //=============================================================


            //   reparmid = 0 ---> no reparameterization in any direction
            //            = 1 ---> reparameterization only in the 1st direction
            //            = 2 ---> reparameterization only in the 2nd direction
            //            = 3 ---> reparameterization in both the directions
            //            = 5 ---> reparameterize elements in both the directions
            //            = 6 ---> reparameterize CPs in both the directions
            //


            //    cout << '\t' <<  " patch # : " << (i+1) << '\t' << patchGrp[i].reparmid << endl;

            NurbsSURFACE surface_temp1(Pw, U, V, p, q);

            if(patchGrp[i].reparmid < 5)
            {
               if(patchGrp[i].subDiv[0] > 0 && patchGrp[i].subDiv[1] == 0)
               {
                  KNOTVECTOR XU;

                  if(patchGrp[i].reparmid == 1 || patchGrp[i].reparmid == 3)
                     FindKnotVector(&SurfaceListOriginal[countsurf], 1, patchGrp[i].subDiv[0], XU);
                  else
                     GenKnotVecForRefining(SurfaceListOriginal[countsurf].U, patchGrp[i].subDiv[0], XU);

                  RefineSurfKnotVector1D(&SurfaceListOriginal[countsurf], XU, 1, &surface_temp1);
               }
               else if(patchGrp[i].subDiv[0] == 0 && patchGrp[i].subDiv[1] > 0)
               {
                   KNOTVECTOR XV;

                   if(patchGrp[i].reparmid == 2 || patchGrp[i].reparmid == 3)
                      FindKnotVector(&SurfaceListOriginal[countsurf], 2, patchGrp[i].subDiv[1], XV);
                   else
                      GenKnotVecForRefining(SurfaceListOriginal[countsurf].V, patchGrp[i].subDiv[1], XV);

                   RefineSurfKnotVector1D(&SurfaceListOriginal[countsurf], XV, 2, &surface_temp1);
               }
               else if(patchGrp[i].subDiv[0] > 0 && patchGrp[i].subDiv[1] > 0)
               {
                   KNOTVECTOR XU, XV;

                   switch(patchGrp[i].reparmid)
                   {
                      case  0:
                              GenKnotVecForRefining(SurfaceListOriginal[countsurf].U, patchGrp[i].subDiv[0], XU);
                              GenKnotVecForRefining(SurfaceListOriginal[countsurf].V, patchGrp[i].subDiv[1], XV);
                              break;

                      case  1:
                              FindKnotVector(&SurfaceListOriginal[countsurf], 1, patchGrp[i].subDiv[0], XU);
                              GenKnotVecForRefining(SurfaceListOriginal[countsurf].V, patchGrp[i].subDiv[1], XV);
                              break;

                      case  2:
                              GenKnotVecForRefining(SurfaceListOriginal[countsurf].U, patchGrp[i].subDiv[0], XU);
                              FindKnotVector(&SurfaceListOriginal[countsurf], 2, patchGrp[i].subDiv[1], XV);
                              break;

                      case  3:
                              FindKnotVector(&SurfaceListOriginal[countsurf], 1, patchGrp[i].subDiv[0], XU);
                              FindKnotVector(&SurfaceListOriginal[countsurf], 2, patchGrp[i].subDiv[1], XV);
                              break;
                  }

                  NurbsSURFACE surface_temp4(Pw, U, V, p, q);

                  //refine the knot vector in U direction
                  RefineSurfKnotVector1D(&SurfaceListOriginal[countsurf], XU, 1, &surface_temp4);

                  //refine the knot vector in V direction
                  RefineSurfKnotVector1D(&surface_temp4, XV, 2, &surface_temp1);
               }
            }

            //===============================================
            // 1.) Degree Elevate the Surface
            //===============================================

            tt1 = patchGrp[i].elevdegree[0] - p;
            tt2 = patchGrp[i].elevdegree[1] - q;


            NurbsSURFACE surface_tempD(Pw, U, V, p, q);

            //degree elevate the surface in U direction
            DegreeElevateSurf1D(&surface_temp1, tt1, 1, &surface_tempD);

            //degree elevate the surface in V direction
            DegreeElevateSurf1D(&surface_tempD, tt2, 2, &SurfaceListFinal[countsurf]);

        }

        SurfaceListFinal[countsurf].ndof = ndf;
        SurfaceListFinal[countsurf].initializeBCdata();
        SurfaceListFinal[countsurf].computeNET();

        //cout << closed1 << '\t' << closed2 << endl;

        SurfaceListFinal[countsurf].closed[0] = closed1;
        SurfaceListFinal[countsurf].closed[1] = closed2;
            
        //cout << SurfaceListFinal[countsurf].closed  << endl;

        index = patch[i][0] - 1;

        SurfaceListFinal[countsurf].ElemProp = patchElemProp[patchgrpdata[index][4*(patchGrp[i].ndom)]-1];
        SurfaceListFinal[countsurf].MatlProp = patchMatlProp[patchgrpdata[index][4*(patchGrp[i].ndom)+1]-1];


        SurfaceResult[countsurf].Pw =  SurfaceListFinal[countsurf].Pw;
        SurfaceResult[countsurf].U  =  SurfaceListFinal[countsurf].U;
        SurfaceResult[countsurf].V  =  SurfaceListFinal[countsurf].V;
        SurfaceResult[countsurf].p  =  SurfaceListFinal[countsurf].p;
        SurfaceResult[countsurf].q  =  SurfaceListFinal[countsurf].q;

        SurfaceResult[countsurf].ndof = ndf;
        SurfaceResult[countsurf].initializeBCdata();
        SurfaceResult[countsurf].computeNET();

        SurfaceResult[countsurf].closed[0] = closed1;
        SurfaceResult[countsurf].closed[1] = closed2;

        countsurf++;

    }  // end of if condition for surfaces
    else
    {
        KNOTVECTOR U; U = kv[patchGrp[i].kvindex[0] - 1];
        KNOTVECTOR V; V = kv[patchGrp[i].kvindex[1] - 1];
        KNOTVECTOR W; W = kv[patchGrp[i].kvindex[2] - 1];

        DEGREE p = patchGrp[i].polydegree[0];
        DEGREE q = patchGrp[i].polydegree[1];
        DEGREE r = patchGrp[i].polydegree[2];

        int mU = U.n -1;
        int mV = V.n -1;
        int mW = W.n -1;
        int nU = mU - p ;
        int nV = mV - q ;
        int nW = mW - r ;

        bool  closed1=false, closed2=false, closed3=false;
        
        //closed1 = checkClosedness(3, 1, nU, nV, nW, &(patchGrp[i].ctrlpoints[0]));
        //cout << endl;        cout << endl;
        //closed2 = checkClosedness(3, 2, nU, nV, nW, &(patchGrp[i].ctrlpoints[0]));
        //closed3 = checkClosedness(3, 3, nU, nV, nW, &(patchGrp[i].ctrlpoints[0]));
        
        //cout << patchGrp[i].ctrlpoints << endl;
        cout <<  " closed1 & closed2 & closed3 " << closed1 << '\t' << closed2 << '\t' << closed3 << endl;
        //cout << nU << '\t' << nV << endl;

        CNET3D Pw;
        Pw.setDim(nW);
        for(kk=0;kk<nW;kk++)
        {
           Pw[kk].setDim(nU);
           for(ii=0;ii<nU;ii++)
              Pw[kk][ii].setDim(nV);
        }

       
        ind3 = 0;
        for(kk=0;kk<nW;kk++)
        {
            for(ii=0;ii<nU;ii++)
            {
                for(jj=0;jj<nV;jj++)
                {
                   ind1 = patchGrp[i].ctrlpoints[ind3++] - 1;
                   ind2 = (ndm+1)*ind1;

                   xx = x.x[ind2];
                   yy = x.x[ind2 + 1];
                   zz = x.x[ind2 + 2];
                   wt = x.x[ind2 + 3];

                   Pw[kk][ii][jj].x = xx * wt ;
                   Pw[kk][ii][jj].y = yy * wt ;
                   Pw[kk][ii][jj].z = zz * wt ;
                   Pw[kk][ii][jj].w = wt ;
                }
            }
        }

        SolidListOriginal[countsolid].Pw = Pw;
        SolidListOriginal[countsolid].U  = U;
        SolidListOriginal[countsolid].V  = V;
        SolidListOriginal[countsolid].W  = W;
        SolidListOriginal[countsolid].p  = p;
        SolidListOriginal[countsolid].q  = q;
        SolidListOriginal[countsolid].r  = r;

        SolidListOriginal[countsolid].ndof = ndf;

        SolidListOriginal[countsolid].initializeBCdata();

            //===============================================
            // 1.) Degree Elevate the Solid
            //===============================================
            
            NurbsSOLID  solid_temp1(Pw, U, V, W, p, q, r);

            DegreeElevateSolid(&SolidListOriginal[countsolid], &(patchGrp[i].elevdegree[0]), &solid_temp1);

            //===============================================
            // 1.) Refine the Solid
            //===============================================

            if(patchGrp[i].reparmid > 5)
            {
               cerr << " ERROR in 'reparmid' for SOLIDs ... " << endl;
               return;
            }

            RefineSolidKnotVector(&solid_temp1, patchGrp[i].reparmid, &(patchGrp[i].subDiv[0]), &SolidListFinal[countsolid]);
            
/*
            cout << SolidListFinal[countsolid].U << endl;
            cout << SolidListFinal[countsolid].V << endl;
            cout << SolidListFinal[countsolid].W << endl;
            cout << SolidListFinal[countsolid].p << endl;
            cout << SolidListFinal[countsolid].q << endl;
            cout << SolidListFinal[countsolid].r << endl;
*/

            SolidListFinal[countsolid].ndof = ndf;
            SolidListFinal[countsolid].initializeBCdata();
            SolidListFinal[countsolid].computeNET();

            SolidResult[countsolid].Pw =  SolidListFinal[countsolid].Pw;
            SolidResult[countsolid].U  =  SolidListFinal[countsolid].U;
            SolidResult[countsolid].V  =  SolidListFinal[countsolid].V;
            SolidResult[countsolid].W  =  SolidListFinal[countsolid].W;
            SolidResult[countsolid].p  =  SolidListFinal[countsolid].p;
            SolidResult[countsolid].q  =  SolidListFinal[countsolid].q;
            SolidResult[countsolid].r  =  SolidListFinal[countsolid].r;

            SolidResult[countsolid].ndof = ndf;
            SolidResult[countsolid].initializeBCdata();
            SolidResult[countsolid].computeNET();

        SolidListFinal[countsolid].ElemProp = patchElemProp[patchgrpdata[index][4*(patchGrp[i].ndom)]-1];
        SolidListFinal[countsolid].MatlProp = patchMatlProp[patchgrpdata[index][4*(patchGrp[i].ndom)+1]-1];

        countsolid++;
    }

  } // end of for loop A

    subDivStab = false;
    //subDivStab = true;
    //eqOrder    = false;
    eqOrder    = true;

    if(mixedSolverFlag == 7 || mixedSolverFlag == 8)
    {
        cout << " mixed methods .... " << endl;
        int  e, iii, pnum, sizep;

        surfSecondVar.setDim(patch.n);

        for(iii=0;iii<patch.n;iii++)
        {
            KNOTVECTOR XU, XV;

            NurbsSURFACE surface_temp5(SurfaceListOriginal[iii].Pw, SurfaceListOriginal[iii].U, 
		    SurfaceListOriginal[iii].V, SurfaceListOriginal[iii].p, SurfaceListOriginal[iii].q);

            NurbsSURFACE surface_temp6(SurfaceListOriginal[iii].Pw, SurfaceListOriginal[iii].U, 
		    SurfaceListOriginal[iii].V, SurfaceListOriginal[iii].p, SurfaceListOriginal[iii].q);

            NurbsSURFACE surface_temp7(SurfaceListOriginal[iii].Pw, SurfaceListOriginal[iii].U, 
		    SurfaceListOriginal[iii].V, SurfaceListOriginal[iii].p, SurfaceListOriginal[iii].q);

            tt1 = patchGrp[iii].elevdegree[0] - SurfaceListOriginal[iii].p;
            tt2 = patchGrp[iii].elevdegree[1] - SurfaceListOriginal[iii].q;

            if(!eqOrder)
            {
              tt1 -= 1;
              tt2 -= 1;
            }
            cout << tt1 << '\t' << tt2 << endl;

            if(tt1 > 0)
	    {
             //degree elevate the surface in U direction
              DegreeElevateSurf1D(&SurfaceListOriginal[iii], tt1, 1, &surface_temp5);
	    }
	    else
	    {
	      surface_temp5.Pw = SurfaceListOriginal[iii].Pw;
	      surface_temp5.U  = SurfaceListOriginal[iii].U;
	      surface_temp5.V  = SurfaceListOriginal[iii].V;
	      surface_temp5.p  = SurfaceListOriginal[iii].p;
	      surface_temp5.q  = SurfaceListOriginal[iii].q;
	    }

            if(tt2 > 0)
	    {
              //degree elevate the surface in V direction
              DegreeElevateSurf1D(&surface_temp5, tt2, 2, &surface_temp6);
	    }
            else
            {
              surface_temp6.Pw = surface_temp5.Pw;
              surface_temp6.U = surface_temp5.U;
              surface_temp6.V = surface_temp5.V;
	      surface_temp6.p = surface_temp5.p;
	      surface_temp6.q = surface_temp5.q;
            }

            tt1 = patchGrp[iii].subDiv[0];
            tt2 = patchGrp[iii].subDiv[1];

            if(subDivStab)
            {
              tt1 -= 1;
              tt2 -= 1;
            }

            cout << tt1 << '\t' << tt2 << endl;

            //refine the knot vector in U direction
            if(tt1 > 0)
	    {
              GenKnotVecForRefining(surface_temp6.U, tt1, XU);
              cout << " LLLLLLLLLLL " << endl;
              RefineSurfKnotVector1D(&surface_temp6, XU, 1, &surface_temp7);
              cout << " cccccccccccccc " << endl;
	    }
	    else
            {
              surface_temp7.Pw = surface_temp6.Pw;
              surface_temp7.U = surface_temp6.U;
              surface_temp7.V = surface_temp6.V;
	      surface_temp7.p = surface_temp6.p;
	      surface_temp7.q = surface_temp6.q;
            }

            if(tt2 > 0)
	    {
              //refine the knot vector in V direction
              GenKnotVecForRefining(surface_temp7.V, tt2, XV);
              cout << " cccccccccccccc " << endl;
              RefineSurfKnotVector1D(&surface_temp7, XV, 2, &surfSecondVar[iii]);
              cout << " LLLLLLLLLLL " << endl;
	    }
	    else
            {
              surfSecondVar[iii].Pw = surface_temp7.Pw;
              surfSecondVar[iii].U = surface_temp7.U;
              surfSecondVar[iii].V = surface_temp7.V;
	      surfSecondVar[iii].p = surface_temp7.p;
	      surfSecondVar[iii].q = surface_temp7.q;
            }


            cout << SurfaceListFinal[iii].U << endl;
            cout << SurfaceListFinal[iii].V << endl;
            cout << " AAAAAAAAAAAAAA " << endl;
            cout << surfSecondVar[iii].U << endl;
            cout << surfSecondVar[iii].V << endl;
            cout << " LLLLLLLLLLL " << endl;
            cout << " wwwwwwwwwwwww " << endl;
            //

            surfSecondVar[iii].ndof = 1;

            surfSecondVar[iii].ElemProp.data = SurfaceListFinal[iii].ElemProp.data;
            cout << " AAAAAAAAAAAAAA " << endl;
            surfSecondVar[iii].initializeBCdata();
            cout << " wwwwwwwwwwwww " << endl;
        }
    }

    /*
    if(mixedSolverFlag == 12)
    {
        int e, iii, pnum, sizep;

        solidSecondVar.setDim(patch.n);

        for(iii=0;iii<patch.n;iii++)
        {
            getLowerOrderKnotVector(SolidListFinal[iii].U, SolidListFinal[iii].p, solidSecondVar[iii].U);
            getLowerOrderKnotVector(SolidListFinal[iii].V, SolidListFinal[iii].q, solidSecondVar[iii].V);
            getLowerOrderKnotVector(SolidListFinal[iii].W, SolidListFinal[iii].r, solidSecondVar[iii].W);
            solidSecondVar[iii].p = SolidListFinal[iii].p-1;
            solidSecondVar[iii].q = SolidListFinal[iii].q-1;
            solidSecondVar[iii].r = SolidListFinal[iii].r-1;
            solidSecondVar[iii].ndof = 1;
            
            solidSecondVar[iii].ElemProp.data = SolidListFinal[iii].ElemProp.data;

            solidSecondVar[iii].initializeBCdata();
         }
    }
    */

    if(mixedSolverFlag == 12)
    {
        int  e, iii, pnum, sizep, size2;

        solidSecondVar.setDim(patch.n);

        KNOTVECTOR XU, XV, XW;

        for(iii=0;iii<patch.n;iii++)
        {
          if(eqOrder)
          {
            solidSecondVar[iii].p = SolidListFinal[iii].p;
            solidSecondVar[iii].q = SolidListFinal[iii].q;
            solidSecondVar[iii].r = SolidListFinal[iii].r;
          }
          else
          {
            solidSecondVar[iii].p = SolidListFinal[iii].p-1;
            solidSecondVar[iii].q = SolidListFinal[iii].q-1;
            solidSecondVar[iii].r = SolidListFinal[iii].r-1;
          }

          if(subDivStab)
          {
            //if(patchGrp[iii].subDiv[0] == 1)
              //solidSecondVar[iii].U = SolidListOriginal[iii].U;
            //else
            {
              GenKnotVecForRefining(SolidListOriginal[iii].U, patchGrp[iii].subDiv[0]-1, XU);

              cout << " AAAAAAAAAAAAAA " << endl;
              cout << XU << endl;
              solidSecondVar[iii].U.setDim(XU.n + 2*(solidSecondVar[iii].p+1));
              for(ii=0;ii<=solidSecondVar[iii].p;ii++)
              {
                solidSecondVar[iii].U[ii] = 0.0;
                solidSecondVar[iii].U[solidSecondVar[iii].U.n-1-ii] = 1.0;
              }
 
              for(ii=0;ii<XU.n;ii++)
                solidSecondVar[iii].U[solidSecondVar[iii].p+1+ii] = XU[ii];
            }

            //if(patchGrp[iii].subDiv[1] == 1)
              //solidSecondVar[iii].V = SolidListOriginal[iii].V;
            //else
            {
              GenKnotVecForRefining(SolidListOriginal[iii].V, patchGrp[iii].subDiv[1]-1, XV);

              cout << " AAAAAAAAAAAAAA " << endl;
              cout << XV << endl;
              solidSecondVar[iii].V.setDim(XV.n + 2*(solidSecondVar[iii].q+1));
              for(ii=0;ii<=solidSecondVar[iii].q;ii++)
              {
                solidSecondVar[iii].V[ii] = 0.0;
                solidSecondVar[iii].V[solidSecondVar[iii].V.n-1-ii] = 1.0;
              }

              for(ii=0;ii<XV.n;ii++)
                solidSecondVar[iii].V[solidSecondVar[iii].q+1+ii] = XV[ii];
            }

            //if(patchGrp[iii].subDiv[2] == 1)
              //solidSecondVar[iii].W = SolidListOriginal[iii].W;
            //else
            {
              GenKnotVecForRefining(SolidListOriginal[iii].W, patchGrp[iii].subDiv[2]-1, XW);

              cout << " AAAAAAAAAAAAAA " << endl;
              cout << XW << endl;
              solidSecondVar[iii].W.setDim(XW.n + 2*(solidSecondVar[iii].r+1));
              for(ii=0;ii<=solidSecondVar[iii].r;ii++)
              {
                solidSecondVar[iii].W[ii] = 0.0;
                solidSecondVar[iii].W[solidSecondVar[iii].W.n-1-ii] = 1.0;
              }

              for(ii=0;ii<XW.n;ii++)
                solidSecondVar[iii].W[solidSecondVar[iii].r+1+ii] = XW[ii];
            }
            
            /*
            vector<double>  vectmp, vtmp2, vtmp3;
            
            for(ii=0;ii<=solidSecondVar[iii].r;ii++)
              vectmp.push_back(0.0);

            //vtmp2.push_back(0.125);
            //vtmp2.push_back(0.25);
            //vtmp2.push_back(0.375);

            //for(ii=0;ii<vtmp2.size();ii++)
              //vectmp.push_back(vtmp2[ii]);

            size2 = pow(2, patchGrp[iii].subDiv[2] - 1) - 1;
            for(ii=0;ii<size2;ii++)
              vectmp.push_back(XW[ii]);

            for(ii=0;ii<solidSecondVar[iii].r;ii++)
              vectmp.push_back(0.5);

            //vtmp3.push_back(0.625);
            //vtmp3.push_back(0.75);
            //vtmp3.push_back(0.875);
            
            //for(ii=0;ii<vtmp3.size();ii++)
              //vectmp.push_back(vtmp3[ii]);

            for(ii=0;ii<size2;ii++)
              vectmp.push_back(XW[size2+ii]);

            for(ii=0;ii<=solidSecondVar[iii].r;ii++)
              vectmp.push_back(1.0);

            solidSecondVar[iii].W.setDim(vectmp.size());

            for(ii=0;ii<vectmp.size();ii++)
              solidSecondVar[iii].W[ii] = vectmp[ii];
            */


            cout << SolidListOriginal[iii].U << endl;
            cout << XU << endl;
            cout << SolidListOriginal[iii].V << endl;
            cout << XV << endl;
            cout << SolidListOriginal[iii].W << endl;
            cout << XW << endl;
            cout << endl;cout << endl;cout << endl;
          }
          else
          {
            if(eqOrder)
            {
              solidSecondVar[iii].U  =  SolidListFinal[iii].U;
              solidSecondVar[iii].V  =  SolidListFinal[iii].V;
              solidSecondVar[iii].W  =  SolidListFinal[iii].W;
            }
            else
            {
              getLowerOrderKnotVector(SolidListFinal[iii].U, SolidListFinal[iii].p, solidSecondVar[iii].U);
              getLowerOrderKnotVector(SolidListFinal[iii].V, SolidListFinal[iii].q, solidSecondVar[iii].V);
              getLowerOrderKnotVector(SolidListFinal[iii].W, SolidListFinal[iii].r, solidSecondVar[iii].W);
            }
          }

            cout << " BBBBBBBBBBBBBB " << endl;
            cout << " BBBBBBBBBBBBBB " << endl;
            cout << SolidListFinal[iii].U << endl;
            cout << solidSecondVar[iii].U << endl;
            cout << " BBBBBBBBBBBBBB " << endl;
            cout << SolidListFinal[iii].V << endl;
            cout << solidSecondVar[iii].V << endl;
            cout << " BBBBBBBBBBBBBB " << endl;
            cout << SolidListFinal[iii].W << endl;
            cout << solidSecondVar[iii].W << endl;
            cout << " BBBBBBBBBBBBBB " << endl;

            solidSecondVar[iii].ndof = 1;
            //cout << " LLLLLLLLLLL " << endl;

            solidSecondVar[iii].ElemProp.data = SolidListFinal[iii].ElemProp.data;
            //cout << " LLLLLLLLLLL " << endl;

            solidSecondVar[iii].initializeBCdata();
            //cout << " LLLLLLLLLLL " << endl;
        }
    }

    /*
         cout << endl;
         cout << endl;
         cout << CurveListOriginal[0].U << endl;
         cout << CurveListFinal[0].U << endl;
         cout << CurveResult[0].U << endl;
         cout << endl;
         cout << endl;
    //
         cout << endl;
         cout << endl;
         cout << SurfaceListOriginal[0].U << endl;
         cout << SurfaceListFinal[0].U << endl;
         cout << SurfaceResult[0].U << endl;
         cout << endl;
         cout << endl;
    */
    
    cout << " ttttttttttttttt " << endl;
    cout << " LLLLLLLLLLL " << endl;

    NurbsBaseOriginal = new NurbsBASE* [patch.n];
    NurbsBaseFinal    = new NurbsBASE* [patch.n];
    NurbsBaseResult   = new NurbsBASE* [patch.n];

    if(mixedSolverFlag > 0)
    {
       solidSecondVar.setDim(patch.n);
       NurbsBaseSecondVar = new NurbsBASE* [patch.n];
    }
 
    for(i=0;i<patchGrp.n;i++)
    {
       if(patchGrp[i].ndom == 1)
       {
           NurbsBaseOriginal[i] =  new NurbsCURVE;
           NurbsBaseFinal[i]    =  new NurbsCURVE;
           NurbsBaseResult[i]   =  new NurbsCURVE;

           NurbsBaseOriginal[i] =  &CurveListOriginal[i];
           NurbsBaseFinal[i]    =  &CurveListFinal[i];
           NurbsBaseResult[i]   =  &CurveResult[i];
       }
       else if(patchGrp[i].ndom == 2)
       {
           NurbsBaseOriginal[i] =  new NurbsSURFACE;
           NurbsBaseFinal[i]    =  new NurbsSURFACE;
           NurbsBaseResult[i]   =  new NurbsSURFACE;

           NurbsBaseOriginal[i] =  &SurfaceListOriginal[i];
           NurbsBaseFinal[i]    =  &SurfaceListFinal[i];
           NurbsBaseResult[i]   =  &SurfaceResult[i];

           if(mixedSolverFlag > 0)
           {
              NurbsBaseSecondVar[i] = new NurbsSURFACE;
              NurbsBaseSecondVar[i] = &surfSecondVar[i];
           }
       }
       else
       {
           NurbsBaseOriginal[i] =  new NurbsSOLID;       
           NurbsBaseFinal[i]    =  new NurbsSOLID;
           NurbsBaseResult[i]   =  new NurbsSOLID;

           NurbsBaseOriginal[i] =  &SolidListOriginal[i];
           NurbsBaseFinal[i]    =  &SolidListFinal[i];
           NurbsBaseResult[i]   =  &SolidResult[i];

           if(mixedSolverFlag > 0)
           {
              NurbsBaseSecondVar[i] = new NurbsSOLID;
              NurbsBaseSecondVar[i] = &solidSecondVar[i];
           }
      }
    }

    totnumel = 0;
    for(i=0;i<patchGrp.n;i++)
      totnumel += NurbsBaseFinal[i]->nelem;

    //=============================================================
    // process patch interface data for multipatch domains
    //=============================================================
    
    if(Npatch > 1)
       processInterfaceData();

    //=============================================================
    // generate elements
    //=============================================================

    elem = new NurbsElement* [totnumel];

    int e=0, e1, e2, e3, e4, n0, n1, n2, nn;

    int  elemtype = patchGrp[0].eltype;

    //mixedSolverFlag = elemtype;

    if(elemtype == 9 || elemtype == 10)
    {
       for(i=0;i<Npatch;i++)
         SurfaceResult[i].PLATE_BENDING = true;
    }

    if(elemtype == 3 || ( (elemtype >= 14) && (elemtype != 21) && (elemtype != 22) ) )
      LSFEMflag = true;

    cout << " elemtype " << elemtype << endl;
    cout << " mixedSolverFlag  " << mixedSolverFlag  << endl;
    cout << " LSFEMflag  " << LSFEMflag  << endl;

    e=0;
    for(i=0;i<patch.n;i++)
    {
       if(patchGrp[i].ndom == 1)   // curves
       {
          VectorArray<double> uu1;
          findunique(CurveListFinal[i].U, uu1);

          n0 = uu1.n-1;
          for(e1=0;e1<n0;e1++)
          {
             elem[e] = newElement(elemtype);

             elem[e]->curve0 = &CurveListFinal[i];
             elem[e]->curve1 = &CurveResult[i];

             elem[e]->ndof     = ndf;
             elem[e]->elenum   = e1;
             elem[e]->patchnum = i;
             elem[e]->prepareElemData();

             e++;
          }
       }
       else if(patchGrp[i].ndom == 2) //surfaces
       {
          VectorArray<double> uu1, vv1;
          findunique(SurfaceListFinal[i].U, uu1);
          findunique(SurfaceListFinal[i].V, vv1);
          n0 = uu1.n-1;
          n1 = vv1.n-1;

          e3=0;
          for(e2=0;e2<n1;e2++)
          {
            for(e1=0;e1<n0;e1++)
            {
                elem[e] = newElement(elemtype);

                elem[e]->surf0 = &SurfaceListFinal[i];
                elem[e]->surf1 = &SurfaceResult[i];

                elem[e]->patchnum = i;
                elem[e]->ndof     = ndf;
                elem[e]->elenum   = e3;

                nn = n0/2*(e2/2)+e1/2;
                //cout << e3 << '\t' << nn << endl;

                if(subDivStab)
                  elem[e]->elenum2  = nn;
                else
                  elem[e]->elenum2  = e3;

                elem[e]->uvalues[0] = uu1[e1];
                elem[e]->uvalues[1] = uu1[e1+1];
                elem[e]->vvalues[0] = vv1[e2];
                elem[e]->vvalues[1] = vv1[e2+1];

                if(mixedSolverFlag > 0)
                   elem[e]->surf2 = &(surfSecondVar[i]);

                //cout << " uuuuuuuuuuu " << endl;
                elem[e]->prepareElemData();
                //cout << " uuuuuuuuuuu " << endl;

                e++; e3++;
            } // for(e1=0;e1<n0;e1++)
          } // for(e2=0;e2<n1;e2++)
          cout << " lllllllllll " << endl;
          /*
          if(subDivStab && mixedSolverFlag > 0)
	  {
            VectorArray<double> uu2, vv2;
            findunique(surfSecondVar[i].U, uu2);
            findunique(surfSecondVar[i].V, vv2);
            
            cout << uu2 << endl;
            cout << vv2 << endl;

            int m0 = uu2.n-1;
            int m1 = vv2.n-1;
            int abc[4], c1;

            for(e2=0;e2<m1;e2++)
            {
              for(e1=0;e1<m0;e1++)
              {
                nn = m0*e2+ e1;

                abc[0] = n0*m0*e2 + 2*e1;
                abc[1] = abc[0] + 1;
                abc[2] = abc[0] + n0;
                abc[3] = abc[2] + 1;
                cout << abc[0] << '\t' << abc[1] << '\t' << abc[2] << '\t' << abc[3] << endl;

                for(c1=0; c1<4; c1++)
		{
                  elem[abc[c1]]->elenum2  = nn;
                  //cout << " uuuuuuuuuuu " << endl;
                  elem[abc[c1]]->prepareElemData();
		}
                //cout << " uuuuuuuuuuu " << endl;

                //e++;
              } // for(e1=0;e1<n0;e1++)
            } // for(e2=0;e2<n1;e2++)
	  }
	  */
       }
       else        //solids
       {
          VectorArray<double> uu1, vv1, ww1;
          findunique(SolidListFinal[i].U, uu1);
          findunique(SolidListFinal[i].V, vv1);
          findunique(SolidListFinal[i].W, ww1);
          n0 = uu1.n-1;
          n1 = vv1.n-1;
          n2 = ww1.n-1;
          //cout << uu1.n << '\t' << vv1.n << '\t' << ww1.n << endl;

          e4=0;
          for(e3=0;e3<n2;e3++)
          {
              for(e2=0;e2<n1;e2++)
              {
                 for(e1=0;e1<n0;e1++)
                 {
                    //cout << e1 << '\t' << e2 << '\t' << e3 << '\t' << e << '\t' << e4 << endl;
                    elem[e] = newElement(elemtype);

                    elem[e]->solid0 = &SolidListFinal[i];
                    elem[e]->solid1 = &SolidResult[i];

                    elem[e]->ndof     = ndf;
                    elem[e]->patchnum = i;
                    elem[e]->elenum   = e4;

                    nn = (n1*n0/4)*(e3/2) + n0/2*(e2/2) + e1/2;
                    //cout << e4 << '\t' << nn << endl;
                    if(subDivStab)
                      elem[e]->elenum2  = nn;
                    else
                      elem[e]->elenum2  = e4;

                    elem[e]->uvalues[0] = uu1[e1];
                    elem[e]->uvalues[1] = uu1[e1+1];
                    elem[e]->vvalues[0] = vv1[e2];
                    elem[e]->vvalues[1] = vv1[e2+1];
                    elem[e]->wvalues[0] = ww1[e3];
                    elem[e]->wvalues[1] = ww1[e3+1];

                    if(mixedSolverFlag > 0)
                      elem[e]->solid2 = &(solidSecondVar[i]);

                    elem[e]->prepareElemData();
                    
                    e++; e4++;
                 }
              }
          }
       }
    }


  for(int iii=0;iii<Npatch;iii++)
    NurbsBaseResult[iii]->setTimeIntegrationParameters(tis, rhoInfty);


    cout << "   intVarFlag   " << intVarFlag << endl;

      // set the global intVarFlag

      if( (elem[0]->nivGP) > 0)
        intVarFlag = true;
      else
        intVarFlag = false;

      cout << "   intVarFlag   " << intVarFlag << endl;

//    EPOINT  EP;
//    EP = SurfaceListFinal[0].SurfacePoint(0.5,0.5).CalcEuclid();
//    EP.print2screen();

    //=============================================================
    //
    // assign point force BCs data
    //
    //=============================================================

    if(forcebc.n > 0) // check if force BCs are specified
    {
       cerr << " Point force data processing yet to be implemented ... " << endl;
       /*
       int jj=0, kk=0, ll=0;
       for(i=0;i<forcebc.n;i++)
       {
         jj = ibc2[i][0]-1;  // patch number
         kk = ibc2[i][1]-1;  // gbf number
         ll = ibc2[i][2]-1;  // dof number

           NurbsBaseOriginal[jj]->forceBCs[kk][ll] = forcebc[i];
       }
       */
    }
    
    //testFunction();

  return;
}






/*
int IsogeometricFEM::calcStiffnessAndResidual(int printRes, bool zeroMtx, bool zeroRes)
{
   //cout << "     ISOGEOMETRICFEM: generating coefficient Matrices ...\n\n";

   char fct[] = "IsogeometricFEM::calcStiffnessAndResidual";

   computerTime.go(fct);

   if(intVarFlag)
     copyElemInternalVariables();

   //cout << " PPPPPPPPPPPP " << endl;
   solver->zeroMtx();

   if(firstIter) rNorm = -1.0;

   rhsVec.zero();
   reac.zero();

   int start1,  start2, ee;
   start1 = start2  =  0;
   if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
      start1 = start2  =  ntoteqs1;
   if(mixedSolverFlag == 8)
   {
      start1 = ntoteqs1 ;
      start2 = ntoteqs1 + ntoteqs2;
   }

   //printData(4,0);
   for(ee=0;ee<totnumel;ee++)  // loop over all the elements
   {
      //cout << "       elem... : " << (ee+1) << endl;

      localStiffnessError = elem[ee]->calcStiffnessAndResidual();

      //cout << " AAAAAAAAAA " << endl;
      if (localStiffnessError != 0)
      {
         cout << '\t' << "local element failure!\n\n";
         return localStiffnessError;
      }

      //elem[ee]->printStiffnessMatrix();
      //elem[ee]->printForceVector();

      //cout << " AAAAAAAAAA " << endl;
      if(SOLVER_TYPE >= 4)
        elem[ee]->AssembleElementMatrix(1, ((SolverEigen*)solver)->mtx, start1, start2);
      else
        elem[ee]->AssembleElementMatrix(1, ((SolverSparse*)solver)->mtx);
      
      //cout << " AAAAAAAAAA " << endl;
      //((SolverEigen*)solver)->printMatrix(1,1,1,1,1);

      elem[ee]->AssembleElementVector(firstIter, 0, &(rhsVec[0]), &(reac[0]), start1, start2);
      //printData(3,0);
      //cout << " AAAAAAAAAA " << endl;
   }

   //cout << " rhsVec " << endl;        printVector(&(rhsVec[0]), rhsVec.n);

   //printData(4,0);
   //printData(3,0);

   printf("\n rhsVec norm = %12.6E \n", rhsVec.norm());

   int  elemtype = patchGrp[0].eltype;

  if(LSFEMflag)
  {
    int  temp;
    for(ee=0;ee<totnumel;ee++)  // loop over all the elements
    {
      if(elem[ee]->tracflag)
      {
         temp = elem[ee]->calcLoadVector();
         //elem[e]->AssembleElementMatrix(1, ((SolverSparse*)solver)->mtx);
         elem[ee]->AssembleElementVector(1, 0, &(rhsVec[0]), &(reac[0]), 0, 0);

         //temp = elem[e]->applyDirichletBCs();
         //elem[e]->AssembleElementMatrix(1, ((SolverSparse*)solver)->mtx);
         //elem[e]->AssembleElementVector(firstIter, 0, &(rhsVec[0]), &(reac[0]), 0, 0);
      }
    }
  }
  else
    addExternalForces();

  //printData(4,0);
  //printData(3,0);

  //cout << " rhsVec " << endl;        printVector(&(rhsVec[0]), rhsVec.n);

  printf("\n rhsVec norm = %12.6E \n", rhsVec.norm());

  if(firstIter)
  {
    int iii, ii, kk, ind1, ind2, ee;

    for(iii=0;iii<patch.n;iii++)
      NurbsBaseResult[iii]->addInitDOFvalues();

    for(ee=0;ee<totnumel;ee++)
      elem[ee]->initialiseDOFvalues();

    for(ii=0;ii<solnFull.n;ii++)
      solnFull[ii] = mpapTime.dt * solnInit[ii];

    ind1 = 0;
    for(iii=0;iii<Npatch;iii++)
    {
      for(kk=0;kk<NurbsBaseFinal[iii]->ngbf;kk++)
      {
        ind2 = NurbsBaseFinal[iii]->gbfnums[kk]*ndf;
        for(ii=0;ii<ndf;ii++)
        {
          //cout << iii << '\t' << kk << '\t' << ii << '\t' << ind1+ind2+ii << endl;
          NurbsBaseResult[iii]->Values[kk][ii]  += solnFull[ind2+ii];
        }
      }
    }
  }

  //cout << " kkkkkkkkkkkk " << endl;

  firstIter = false;
  rNormPrev = rNorm;
  rNorm     = rhsVec.norm();

  if (printRes > 1) { COUT << domain.name(this); printf("  %11.4e\n",rNorm);}

  solver->currentStatus = ASSEMBLY_OK;

  ctimCalcStiffRes += computerTime.stop(fct);

//   computerTime.stopAndPrint(fct);

  return 0;
}
*/




//
int IsogeometricFEM::calcStiffnessAndResidual(int printRes, bool zeroMtx, bool zeroRes)
{
  cout << "     ISOGEOMETRICFEM: generating coefficient Matrices ...\n\n";

  int iii, ii, kk, ind1, ind2, ee;
  
  double val[2], hmax=-1.0e-10;

  //for(ee=0;ee<totnumel;ee++)  // loop over all the elements
  //{
    //elem[ee]->computeBounds(val);
    //hmax = max(hmax, val[0]);
  //}
  
  //cout << " hmax = " << hmax << endl;

   char fct[] = "IsogeometricFEM::calcStiffnessAndResidual";
   
   //cout << " firstIter = " << firstIter << endl;

   computerTime.go(fct);

   if(intVarFlag)
     copyElemInternalVariables();

   //cout << " PPPPPPPPPPPP " << endl;
   solver->zeroMtx();

  if(firstIter)
  {
    for(iii=0;iii<patch.n;iii++)
      NurbsBaseResult[iii]->addInitDOFvalues();

    for(ee=0;ee<totnumel;ee++)
      elem[ee]->initialiseDOFvalues();

    for(ii=0;ii<solnFull.n;ii++)
      solnFull[ii] += mpapTime.dt * solnInit[ii];

  if(2 < 1)
  {
    cout << " NurbsBaseResult[0]->Values " << endl;         printVector(&(NurbsBaseResult[0]->Values[0][0]), NurbsBaseResult[0]->Values[0].n);
    cout << "\n\n" << endl;
    cout << " NurbsBaseResult[0]->Values2 " << endl;        printVector(&(NurbsBaseResult[0]->Values[1][0]), NurbsBaseResult[0]->Values[1].n);
    cout << "\n\n" << endl;
    cout << " NurbsBaseResult[0]->Values3 " << endl;        printVector(&(NurbsBaseResult[0]->Values[2][0]), NurbsBaseResult[0]->Values[2].n);
    cout << "\n\n" << endl;
  }
  }

   cout << " firstIter = " << firstIter << endl;

   if(firstIter) rNorm = -1.0;

   solver->rhsVec.setZero();
   reac.zero();

   int start1,  start2;

   start1 = start2  =  0;

   if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
      start1 = start2  =  ntoteqs1;
   if(mixedSolverFlag == 8)
   {
      start1 = ntoteqs1 ;
      start2 = ntoteqs1 + ntoteqs2;
   }

   //cout << " rhsVec " << endl;        printVector(&(solver->rhsVec[0]), solver->rhsVec.rows());

   //cout << " reacc " << endl;        printVector(&(reac[0]), reac.n);

  cout << " SOLVER_TYPE " << SOLVER_TYPE << endl;

   for(ee=0;ee<totnumel;ee++)  // loop over all the elements
   {
      //cout << "       elem... : " << (ee+1) << endl;

      localStiffnessError = elem[ee]->calcStiffnessAndResidual();

      //cout << " AAAAAAAAAA " << endl;
      if (localStiffnessError != 0)
      {
        cout << '\t' << "local element failure! " << localStiffnessError << "\n\n" << endl;
        return localStiffnessError;
      }

      //elem[ee]->printStiffnessMatrix();
      //elem[ee]->printForceVector();

      //cout << " AAAAAAAAAA " << endl;
      if(SOLVER_TYPE >= 4)
        elem[ee]->AssembleElementMatrix(1, solver->mtx, start1, start2);
      else
        elem[ee]->AssembleElementMatrix(1, ((SolverSparse*)solver)->mtx);
      
      //cout << " AAAAAAAAAA " << endl;
      //((SolverEigen*)solver)->printMatrix(1,1,1,1,1);
      
      //cout << " bbbbbbbbbbbbb " << endl;

      elem[ee]->AssembleElementVector(firstIter, 0, &(solver->rhsVec[0]), &(reac[0]), start1, start2);
      //printData(3,0);
      //cout << " AAAAAAAAAA " << endl;
   }

   //cout << " rhsVec " << endl;        printVector(&(rhsVec[0]), rhsVec.n);

   //printData(4,0);
   //printData(3,0);

   printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

   int  elemtype = patchGrp[0].eltype;

  if(LSFEMflag)
  {
    int  temp;
    for(ee=0;ee<totnumel;ee++)  // loop over all the elements
    {
      if(elem[ee]->tracflag)
      {
         temp = elem[ee]->calcLoadVector();
         //elem[e]->AssembleElementMatrix(1, ((SolverSparse*)solver)->mtx);
         elem[ee]->AssembleElementVector(1, 0, &(solver->rhsVec[0]), &(reac[0]), 0, 0);

         //temp = elem[e]->applyDirichletBCs();
         //elem[e]->AssembleElementMatrix(1, ((SolverSparse*)solver)->mtx);
         //elem[e]->AssembleElementVector(firstIter, 0, &(rhsVec[0]), &(reac[0]), 0, 0);
      }
    }
  }
  else
  {
  }
    addExternalForces();

  //printData(4,0);
  //printData(3,0);

  //cout << " rhsVec " << endl;        printVector(&(rhsVec[0]), rhsVec.n);

  //cout << solver->mtx << endl;
  
  printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

  //cout << " kkkkkkkkkkkk " << endl;

  firstIter = false;
  rNormPrev = rNorm;
  rNorm     = solver->rhsVec.norm();

  if (printRes > 1) { COUT << domain.name(this); printf("  %11.4e\n",rNorm);}

  solver->currentStatus = ASSEMBLY_OK;

  ctimCalcStiffRes += computerTime.stop(fct);

//   computerTime.stopAndPrint(fct);

  return 0;
}
//



int IsogeometricFEM::factoriseSolveAndUpdate()
{
  char fct[] = "IsogeometricFEM::factoriseSolveAndUpdate";
  computerTime.go(fct);

  //double *du;
  int    kk, iii, ii, ind1, ind2;

  int  elemtype = patchGrp[0].eltype;

  //cout << " residue_new " << endl;        printVector(&(rhsVec[0]), rhsVec.n);
  //cout << " rhsVec " << endl;        printVector(&(solver->rhsVec[0]), solver->rhsVec.rows());

  time_t tstart, tend;

  tstart = time(0);

  solver->factoriseAndSolve();

  tend = time(0); 
  printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );

  soln.zero();
  // update solution vector
  for(kk=0;kk<ntoteqs1;kk++)
    soln[assy4r[kk]]  = solver->soln[kk];

  //for(kk=0;kk<soln.n;kk++)
  //{
    //diff[kk] = soln[kk] - solnPrev[kk];
    //solnFull[kk] += soln[kk];
  //}

  //printf(" soln diff norm   =  %12.6E  \n", diff.norm()/solnFull.norm());

  //cout << " solnFull " << endl;        printVector(&(solnFull[0]), solnFull.n);
  //cout << " soln     " << endl;        printVector(&(soln[0]), soln.n);

    ind1 = 0;
    for(iii=0;iii<Npatch;iii++)
    {
      //cout << NurbsBaseFinal[iii]->gbfnums << endl;
      for(kk=0;kk<NurbsBaseFinal[iii]->ngbf;kk++)
      {
        ind2 = NurbsBaseFinal[iii]->gbfnums[kk]*ndf;
        for(ii=0;ii<ndf;ii++)
        {
          //cout << iii << '\t' << kk << '\t' << ii << '\t' << ind1+ind2+ii << endl;
          NurbsBaseResult[iii]->Values[ii][kk]  += soln[ind2+ii];
        }
        //printf("%5d \t %12.8f \t %12.8f \n", kk, NurbsBaseResult[iii]->Values[kk][0], NurbsBaseResult[iii]->Values[kk][1]);
      }
    }

  /*
  cout << " NurbsBaseResult[0]->Values " << endl;         printVector(&(NurbsBaseResult[0]->Values[0][0]), NurbsBaseResult[0]->Values[0].n);
  cout << "\n\n" << endl;
  cout << " NurbsBaseResult[0]->Values2 " << endl;        printVector(&(NurbsBaseResult[0]->Values[1][0]), NurbsBaseResult[0]->Values[1].n);
  cout << "\n\n" << endl;
  cout << " NurbsBaseResult[0]->Values3 " << endl;        printVector(&(NurbsBaseResult[0]->Values[2][0]), NurbsBaseResult[0]->Values[2].n);
  cout << "\n\n" << endl;
  */

  //cout << " solnFull " << endl;        printVector(&(solnFull[0]), solnFull.n);

  cout << " nnnnnnnnnnn " << endl;
  for(iii=0;iii<Npatch;iii++)
    NurbsBaseResult[iii]->updateCoordinates(&(soln[0]));

  if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
  {
    for(iii=0;iii<Npatch;iii++)
      NurbsBaseSecondVar[iii]->updateValues(1, &solver->soln[ntoteqs1]);
  }

  if(mixedSolverFlag == 8)
  {
    int ind = ntoteqs1 + ntoteqs2;

    for(iii=0;iii<Npatch;iii++)
    {
      NurbsBaseSecondVar[iii]->updateValues(1, &solver->soln[ntoteqs1]);
      NurbsBaseSecondVar[iii]->updateValues(2, &solver->soln[ind]);
    }
  }
  //du = NULL;

  ctimFactSolvUpdt += computerTime.stop(fct);

  return 0;
}




void IsogeometricFEM::computeElementErrors(int index)
{
  int  ii, count=0, e;

  totalError = 0.0;

  for(e=0;e<totnumel;e++)  // loop over all the elements
  {
    //cout << "       elem... : " << (e+1) << endl;  //          cout << endl;

    localStiffnessError = elem[e]->calcError(index);

    totalError += elem[e]->GetError();

    //totalError += ( elem[e]->GetError() * elem[e]->GetError() );
  }
    
  totalError = sqrt(totalError);

  printf(" \n\n \t totalError = %12.6E \n\n " , totalError);

  return;
}



void IsogeometricFEM::addExternalForces()
{
   if(globalFirstIter)
   {
      ForceVec.zero();
      calcAndAssyLoadVector(1.0, 0.0);
      globalFirstIter = false;
   }

  double fact, fact1;
  fact = timeFunction[0].prop;

  //for(int ii=0; ii<ForceVec.n; ii++)
    //printf(" %5d \t %12.8f \n", ii, ForceVec[ii] );

  //if(firstIter)     

  cout << "       fact .... : " << fact << endl;

   for(int kk=0;kk<solver->rhsVec.rows();kk++)
   {
      fact1 = fact * ForceVec[kk];
      solver->rhsVec[kk] += fact1;
      reac[assy4r[kk]] += fact1 ;
   }

/*
   assy4F2[0] = 0;
   assy4F2[1] = 0;

   for(int kk=0;kk<assy4F2.n;kk++)
   {
      fact1 = fact * ForceVec[assy4F2[kk]];
      rhsVec[assy4F2[kk]] += fact1;
      reac[assy4r[assy4F2[kk]]] += fact1 ;
   }
*/

  return;
}









NurbsElement* IsogeometricFEM::newElement(int type)
{
  switch (type)
  {
    case  0: return (NurbsElement*) new NurbsElem1DAdvectionDiffusion; break;

    case  1: return (NurbsElement*) new NurbsElem1DElasticBar; break;

    case  2: return (NurbsElement*) new NurbsElem1DEulerBeam; break;

    case  3: return (NurbsElement*) new NurbsElem1DElasticBarLSFEM; break;

    case  4: return (NurbsElement*) new NurbsElem2DStructSolid; break;

    case  5: return (NurbsElement*) new NurbsElem2DStructFbarSolid; break;

    case  6: return (NurbsElement*) new NurbsElem2DStructBbarSolid; break;

    case  7: return (NurbsElement*) new NurbsElem2DStructMixed2field; break;

    case  8: return (NurbsElement*) new NurbsElem2DStructMixed3field; break;

    case  9: return (NurbsElement*) new NurbsElemKirchhoffPlate; break;

    case 10: return (NurbsElement*) new NurbsElemMindlinPlate; break;

    case  11: return (NurbsElement*) new NurbsElem3DStructSolid; break;
    
    case  12: return (NurbsElement*) new NurbsElem3DStructMixed2field; break;

    case  13: return (NurbsElement*) new NurbsElem2DAdvectionDiffusion; break;
    
    case  14: return (NurbsElement*) new NurbsElem2DStokes; break;

    case  15: return (NurbsElement*) new NurbsElem2DNavierStokes3dof; break;

    case  16: return (NurbsElement*) new NurbsElem2DNavierStokes4dof; break;
    
    case  17: return (NurbsElement*) new NurbsElem2DStructSolidLSFEM2dof; break;

    case  18: return (NurbsElement*) new NurbsElem2DStructSolidLSFEM3dof; break;

    case  19: return (NurbsElement*) new NurbsElem2DHeatTransfer; break;

    case  20: return (NurbsElement*) new NurbsElem2DTempCoupled4dof; break;
    
    case  21: return (NurbsElement*) new NurbsElem2DStructMixed2fieldStabilised; break;
    
    case  22: return (NurbsElement*) new NurbsElem3DStructMixed2fieldStabilised; break;

    default: prgError(1,"IsogeometricFEM::newElement","unknown element type name!"); return NULL;
  }
}




void IsogeometricFEM::ProcessDispBoundaryConditions()
{
  if(ndm == 1)
    ProcessDispBCsCurve();
  else if(ndm == 2)
  {
    ProcessDispBCsSurface();
    ProcessDispBCsSurface2();
  }
  else
    ProcessDispBCsSolid();

  return;
}


void IsogeometricFEM::ProcessDispBCsCurve()
{
  cout << "     ISOGEOMETRICFEM: processing disp boundary conditions for Curves...\n\n";

  char fct[] = "IsogeometricFEM::ProcessDispBoundaryConditionsCurve";


  int iii;

     for(iii=0;iii<dispbc.n;iii++) // Loop A
     {
        int index, patchnum, val1, dof, dir, ngbf2;
        
        double  dispval;

        patchnum = (int) dispbc[iii][0]-1;

        ngbf2 = CurveListFinal[patchnum].ngbf;

        val1 = (int) dispbc[iii][1]-1;
        dof  = (int) dispbc[iii][2]-1;
        dir  = (int) dispbc[iii][3]-1; // direction in which traction is applied ( 0 -> 1st direction, 1-> 2nd direction)
        dispval = dispbc[iii][4];

        // search at the left boundary

        if( val1 == 0 )
        {
            CurveListFinal[patchnum].dispBCs[0][dir] = dispval;
            CurveListFinal[patchnum].Uinit[dir] = dispval;
        }

        // search at the right boundary
        else if( val1 == 1 )
        {
            CurveListFinal[patchnum].dispBCs[ngbf2-1][dir] = dispval;
            CurveListFinal[patchnum].Uinit[ndf*(ngbf2-1) + dir] = dispval;
        }
    }

     ntotgbf1 = 0;

     // compute total number of independent global basis funtions
     for(iii=0;iii<Npatch;iii++)
        ntotgbf1 += CurveListFinal[iii].ngbf;

  return;
}





void IsogeometricFEM::ProcessDispBCsSurface()
{
  cout << "     ISOGEOMETRICFEM: processing disp boundary conditions for Surfaces...\n\n";

  char fct[] = "IsogeometricFEM::ProcessDispBoundaryConditionsSurface";

  int ii, jj, iii;

     for(iii=0;iii<Npatch;iii++)
     {
         cout << SurfaceListFinal[iii].edgedata << endl;
         cout << SurfaceListFinal[iii].intfdata << endl;
         cout << endl;
         cout << endl;
     }

     // apply displacement BC values on interface edges if any specified

     if(Npatch > 1)
     {
        int patchnum, side, ngbf21, ngbf22, ngbf31, ngbf32, index1, index2, dir, kk, ndof;

        for(iii=1;iii<Npatch;iii++)
        {
            ngbf21 = SurfaceListFinal[iii].ngbf1;
            ngbf22 = SurfaceListFinal[iii].ngbf2;

            ndof = SurfaceListFinal[iii].ndof;

            for(kk=0;kk<4;kk++)
            {
               if(SurfaceListFinal[iii].edgedata[kk] != -1)
               {
                  patchnum = SurfaceListFinal[iii].intfdata[2*kk];
                  side     = SurfaceListFinal[iii].intfdata[2*kk+1];

                  ngbf31 = SurfaceListFinal[patchnum].ngbf1;
                  ngbf32 = SurfaceListFinal[patchnum].ngbf2;

                 // cout << '\t' << "  AAAA  " << '\t' << patchnum << '\t' << side << '\t' << ngbf31 << '\t' << ngbf32 << endl;

                  switch(kk)
                  {
                      case 0: // left side

                               if(ngbf22 != ngbf32)
                                  prgError(1,fct," ngbf2 in the 2 patches do not match");

                               for(jj=0;jj<ngbf22;jj++)
                               {
                                  index1 = ngbf21*jj;
                                  for(dir=0;dir<ndof;dir++)
                                    SurfaceListFinal[iii].dispBCs[index1][dir] = -9999;
                               }

                              break;

                      case 1: // right side
                              break;

                               if(ngbf22 != ngbf32)
                                  prgError(1,fct," ngbf2 in the 2 patches do not match");

                               for(jj=0;jj<ngbf22;jj++)
                               {
                                  index1 = ngbf21*(jj+1)-1;
                                  for(dir=0;dir<ndof;dir++)
                                    SurfaceListFinal[iii].dispBCs[index1][dir] = -9999;
                               }

                              break;

                      case 2: // bottom side

                               if(ngbf21 != ngbf31)
                                  prgError(1,fct," ngbf1 in the 2 patches do not match");

                               for(jj=0;jj<ngbf21;jj++)
                               {
                                  index1 = jj;
                                  for(dir=0;dir<ndof;dir++)
                                    SurfaceListFinal[iii].dispBCs[index1][dir] = -9999;
                               }

                              break;

                      case 3: // top side

                               if(ngbf21 != ngbf31)
                                  prgError(1,fct," ngbf1 in the 2 patches do not match");

                               for(jj=0;jj<ngbf21;jj++)
                               {
                                  index1 = ngbf31*(ngbf32-1)+jj;
                                  for(dir=0;dir<ndof;dir++)
                                    SurfaceListFinal[iii].dispBCs[index1][dir] = -9999;
                               }

                              break;

                  }
               }
            }
        }
     }

     for(iii=0;iii<dispbc.n;iii++) // Loop A
     {
        int index, patchnum, val1, val2, dir;
        int ngbf1, ngbf2, ngbf11, ngbf12, ngbf21, ngbf22;
        double dispval;

        patchnum = (int) dispbc[iii][0]-1;

        ngbf11 = SurfaceListOriginal[patchnum].ngbf1;
        ngbf12 = SurfaceListOriginal[patchnum].ngbf2;
        ngbf1  = SurfaceListOriginal[patchnum].ngbf;

        ngbf21 = SurfaceListFinal[patchnum].ngbf1;
        ngbf22 = SurfaceListFinal[patchnum].ngbf2;
        ngbf2  = SurfaceListFinal[patchnum].ngbf;

        val1 = (int) dispbc[iii][1]-1;
        val2 = (int) dispbc[iii][2]-1;
        dir  = (int) dispbc[iii][3]-1; // direction in which traction is applied ( 0 -> 1st direction, 1-> 2nd direction)
        dispval = dispbc[iii][4];

        // search at the top boundary
        if( (val1 == 0 && val2 == ngbf11*(ngbf12-1)) || (val2 == 0 && val1 == ngbf11*(ngbf12-1)) )
        {
             for(jj=0;jj<ngbf22;jj++)
             {
                index = ngbf21*jj;
                SurfaceListFinal[patchnum].dispBCs[index][dir] = dispval;
                SurfaceListFinal[patchnum].Uinit[ndf*index + dir] = dispval;
             }
        }

        // search at the bottom boundary
        else if( ((val1 == ngbf11-1) && (val2 == ngbf11*ngbf12-1)) || ((val2 == ngbf11-1) && (val1 == ngbf11*ngbf12-1)) )
        {
             for(jj=0;jj<ngbf22;jj++)
             {
                index = (ngbf21*(jj+1)-1);
                SurfaceListFinal[patchnum].dispBCs[index][dir] = dispval;
                SurfaceListFinal[patchnum].Uinit[ndf*index + dir] = dispval;
             }
        }

        // search at the left boundary
        else if( ((val1 == 0) && (val2 == ngbf11-1)) || ((val2 == 0) && (val1 == ngbf11-1)) )
        {
             for(jj=0;jj<ngbf21;jj++)
             {
                index = jj;
                SurfaceListFinal[patchnum].dispBCs[index][dir] = dispval;
                SurfaceListFinal[patchnum].Uinit[ndf*index + dir] = dispval;
             }
        }

        // search at the right boundary
        else if( ((val1 == ngbf11*(ngbf12-1)) && (val2 == ngbf11*ngbf12-1)) || ((val2 == ngbf11*(ngbf12-1)) && (val1 == ngbf11*ngbf12-1)) )
        {
             int temp = ngbf21*(ngbf22-1);
             for(jj=0;jj<ngbf21;jj++)
             {
                index = temp + jj;
                SurfaceListFinal[patchnum].dispBCs[index][dir] = dispval;
                SurfaceListFinal[patchnum].Uinit[ndf*index + dir] = dispval;
             }
        }
    }

     ntotgbf1 = 0;
     // compute total number of independent global basis funtions
     for(iii=0;iii<Npatch;iii++)
     {
        ntotgbf1 += SurfaceListFinal[iii].ngbf;

        for(int kk=0;kk<4;kk++)
        {
           if(SurfaceListFinal[iii].edgedata[kk] == 3)
           {
             if( kk <=1)
                ntotgbf1 -= SurfaceListFinal[iii].ngbf2;
             else
                ntotgbf1 -= SurfaceListFinal[iii].ngbf1;
           }
        }
     }

/*
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
           {
             cout << SurfaceListFinal[iii].dispBCs[ii][jj] << '\t' ;
           }
           cout << endl;
         }
         cout << endl;
     }
*/


  return;
}




void IsogeometricFEM::ProcessDispBCsSurface2()
{
  cout << "     ISOGEOMETRICFEM: processing disp bcs extra data ...\n\n";
  char fct[] = " IsogeometricFEM::ProcessDispBoundaryConditions2 ";

    int ngbf21=0, ngbf22=0;


    CPOINT CP1, CP2;

    int index, patchnum, side, ttt, index1, index2, dir, iii, ii, jj, kk;
    int val1, val2;
    double dispval;

    Vector<int> temp1;

    for(iii=0;iii<dispbc2.n;iii++)
    {
        patchnum = (int) (dispbc2[iii][0]-1);
            side = (int) (dispbc2[iii][1]-1);
            val1 = (int) (dispbc2[iii][2]-1);
            val2 = (int) (dispbc2[iii][3]-1);
             dir = (int) (dispbc2[iii][4]-1);
         dispval = dispbc2[iii][5];

        cout << '\t' << " disp value " << dispval << endl;

        ngbf21 = SurfaceListFinal[patchnum].ngbf1;
        ngbf22 = SurfaceListFinal[patchnum].ngbf2;

        CP1.x = x.x[(ndm+1)*val1+0];
        CP1.y = x.x[(ndm+1)*val1+1];
        CP2.x = x.x[(ndm+1)*val2+0];
        CP2.y = x.x[(ndm+1)*val2+1];

        if(ndm == 3)
        {
           CP1.z = x.x[(ndm+1)*val1 + 2];
           CP1.w = x.x[(ndm+1)*val1 + 3];
           CP2.z = x.x[(ndm+1)*val2 + 2];
           CP2.w = x.x[(ndm+1)*val2 + 3];
        }
        else
        {
           CP1.z = 0.0;
           CP1.w = x.x[(ndm+1)*val1 + 2];
           CP2.z = 0.0;
           CP2.w = x.x[(ndm+1)*val2 + 2];
        }

        CP1.print2screen();
        CP2.print2screen();

        cout << '\t' << " ngbf21 and ngbf22 " << '\t' << ngbf21 << '\t' << ngbf22 << endl;

        index1 = index2 = 0;

        switch(side)
        {
          case 0:  // search at the top boundary

                 for(jj=1;jj<(ngbf22-1);jj++)
                 {
                    if(SurfaceListFinal[patchnum].Pw[0][jj] == CP1)
                       index1 = jj;

                    if(SurfaceListFinal[patchnum].Pw[0][jj] == CP2)
                       index2 = jj;
                 }

                 if(!(index == 0 && index2 == 0) )
                 {
                    for(kk=index1;kk<=index2;kk++)
                      temp1.append(ngbf21*kk);
                 }

                 break;

          case 1:  // search at the bottom boundary

                 for(jj=1;jj<(ngbf22-1);jj++)
                 {
                    if(SurfaceListFinal[patchnum].Pw[ngbf21-1][jj] == CP1)
                       index1 = jj;

                    if(SurfaceListFinal[patchnum].Pw[ngbf21-1][jj] == CP2)
                       index2 = jj;
                 }

                 if(index != 0 && index2 != 0)
                 {
                    for(kk=index1;kk<=index2;kk++)
                      temp1.append(ngbf21*(kk+1) - 1);
                 }

                 break;

          case 2:  // search at the left boundary

                 for(jj=1;jj<(ngbf21-1);jj++)
                 {
                    if(SurfaceListFinal[patchnum].Pw[jj][0] == CP1)
                       index1 = jj;

                    if(SurfaceListFinal[patchnum].Pw[jj][0] == CP2)
                       index2 = jj;
                 }

                 //if(index != 0 && index2 != 0)
                 {
                    for(int kk=index1;kk<=index2;kk++)
                      temp1.append(kk);
                 }


                 break;

          case 3:  // search at the right boundary

                 ttt = ngbf21*(ngbf22-1);
                 for(jj=1;jj<(ngbf21-1);jj++)
                 {
                   // if(SurfaceListFinal[patchnum].Pw[jj][ngbf22-1] == CP1)
                   //    index1 = jj;

                    if(SurfaceListFinal[patchnum].Pw[jj][ngbf22-1].x <= CP2.x)
                       index2 = jj;
                 }

                 //if(index != 0 && index2 != 0)
                 {
                    for(kk=index1;kk<=index2;kk++)
                      temp1.append(ttt+kk);
                 }

                 break;
        }

        cout << '\t' << index1 << '\t' << index2 << '\t' << dir << endl;
        cout << endl;     cout << '\t' << " temp1 " << temp1 << endl;

        if(temp1.n > 0)
        {
           VectorArray<int> temp2, temp3;
           temp2 = temp1;

           finduniqueInt(temp2, temp3);
           cout << endl;     cout << '\t' << " temp2 " << temp2 << endl;

           SortArrayInt(temp3);
           cout << endl;     cout << '\t' << " temp3 " << temp3 << endl;

           if( dir == 0)
           {
              for(ii=0;ii<temp3.n;ii++)
              {
                SurfaceListFinal[patchnum].dispBCs[temp3[ii]][0] = dispval;
                SurfaceListFinal[patchnum].Uinit[ndf*temp3[ii]] = dispval;
              }
           }
           else if( dir == 1)
           {
              for(ii=0;ii<temp3.n;ii++)
              {
                SurfaceListFinal[patchnum].dispBCs[temp3[ii]][1] = dispval;
                SurfaceListFinal[patchnum].Uinit[ndf*temp3[ii]+1] = dispval;
              }
           }
           else
           {
              for(ii=0;ii<temp3.n;ii++)
              {
                SurfaceListFinal[patchnum].dispBCs[temp3[ii]][0]  = dispval;
                SurfaceListFinal[patchnum].dispBCs[temp3[ii]][1]  = dispval;
                SurfaceListFinal[patchnum].Uinit[ndf*temp3[ii]]   = dispval;
                SurfaceListFinal[patchnum].Uinit[ndf*temp3[ii]+1] = dispval;
              }
           }

        }
    }

  cout << "     ISOGEOMETRICFEM: processing disp bcs extra data ... DONE \n\n";

  return;
}




void IsogeometricFEM::ProcessDispBCsSolid()
{
  cout << "     ISOGEOMETRICFEM: processing disp boundary conditions for Solids...\n\n";

  char fct[] = "IsogeometricFEM::ProcessDispBoundaryConditionsSolid";

  int ii, jj, iii, kk;


     for(iii=0;iii<Npatch;iii++)
     {
         cout << SolidListFinal[iii].edgedata << endl;
         cout << SolidListFinal[iii].intfdata << endl;
         cout << endl;
         cout << endl;
     }

     if(Npatch > 1)
     {
        int patchnum, side, ngbf21, ngbf22, ngbf23, ngbf31, ngbf32, ngbf33, index1, index2, dir, kk, ndof, ngbf1m2, ind1, ind2, ll;

        for(iii=1;iii<Npatch;iii++)
        {
            ngbf21 = SolidListFinal[iii].ngbf1;
            ngbf22 = SolidListFinal[iii].ngbf2;
            ngbf23 = SolidListFinal[iii].ngbf3;
            ndof   = SolidListFinal[iii].ndof;

            ngbf1m2 = ngbf21*ngbf22;

            for(ll=0;ll<6;ll++)
            {
               if(SolidListFinal[iii].edgedata[ll] != -1)
               {
                  patchnum = SolidListFinal[iii].intfdata[2*ll];
                  side     = SolidListFinal[iii].intfdata[2*ll+1];

                  ngbf31 = SolidListFinal[patchnum].ngbf1;
                  ngbf32 = SolidListFinal[patchnum].ngbf2;
                  ngbf33 = SolidListFinal[patchnum].ngbf3;

                  //cout << '\t' << "  AAAA  " << '\t' << patchnum << '\t' << side << '\t' << ngbf31 << '\t' << ngbf32 << '\t' << ngbf33 << endl;

                  switch(ll)
                  {
                       case 0: // search at w = 0

                              if( (ngbf21 != ngbf31) && (ngbf22 != ngbf32) )
                                 prgError(1,fct," ngbf1 and ngbf2 in the 2 patches do not match");

                              ind1 = 0;
                              for(kk=0;kk<ngbf22;kk++)
                              {
                                  for(jj=0;jj<ngbf21;jj++)
                                  {
                                      ind2 = ind1++;
                                      for(dir=0;dir<ndof;dir++)
                                      {
                                           SolidListFinal[iii].dispBCs[ind2][dir] = -9999;
                                      }
                                  }
                              }
                         break;

                       case 1: // search at w = 1

                              if( (ngbf21 != ngbf31) && (ngbf22 != ngbf32) )
                                 prgError(1,fct," ngbf1 and ngbf2 in the 2 patches do not match");

                              ind1 = ngbf1m2 * (ngbf23-1);
                              for(kk=0;kk<ngbf22;kk++)
                              {
                                  for(jj=0;jj<ngbf21;jj++)
                                  {
                                      ind2 = ind1++;
                                      for(dir=0;dir<ndof;dir++)
                                      {
                                           SolidListFinal[iii].dispBCs[ind2][dir] = -9999;
                                      }
                                  }
                              }
                         break;

                       case 2: // search at v = 0

                              if( (ngbf21 != ngbf31) && (ngbf23 != ngbf33) )
                                  prgError(1,fct," ngbf1 and ngbf3 in the 2 patches do not match");

                              for(kk=0;kk<ngbf23;kk++)
                              {
                                  ind1 = ngbf1m2 * kk;
                                  for(jj=0;jj<ngbf21;jj++)
                                  {
                                      ind2 = ind1 + jj;
                                      for(dir=0;dir<ndof;dir++)
                                      {
                                           SolidListFinal[iii].dispBCs[ind2][dir] = -9999;
                                      }
                                  }
                              }
                         break;

                       case 3: // search at v = 1
            
                              if( (ngbf21 != ngbf31) && (ngbf23 != ngbf33) )
                                  prgError(1,fct," ngbf1 and ngbf3 in the 2 patches do not match");

                               for(kk=0;kk<ngbf23;kk++)
                               {
                                   ind1 = ngbf1m2 * (kk+1) - ngbf21 ;
                                   for(jj=0;jj<ngbf21;jj++)
                                   {
                                      ind2 = ind1 + jj;
                                      for(dir=0;dir<ndof;dir++)
                                      {
                                           SolidListFinal[iii].dispBCs[ind2][dir] = -9999;
                                      }
                                   }
                               }

                         break;

                       case 4: //

                               if( (ngbf22 != ngbf32) && (ngbf23 != ngbf33) )
                                  prgError(1,fct," ngbf2 and ngbf3 in the 2 patches do not match");

                               for(kk=0;kk<ngbf22;kk++)
                               {
                                   ind1 = ngbf21*kk;
                                   for(jj=0;jj<ngbf23;jj++)
                                   {
                                      ind2 = ind1 + ngbf1m2*jj;
                                      for(dir=0;dir<ndof;dir++)
                                      {
                                           SolidListFinal[iii].dispBCs[ind2][dir] = -9999;
                                      }
                                   }
                               }

                              break;

                       case 5: // search at u = 1

                               if( (ngbf22 != ngbf32) && (ngbf23 != ngbf33) )
                                  prgError(1,fct," ngbf2 and ngbf3 in the 2 patches do not match");

                              for(kk=0;kk<ngbf23;kk++)
                              {
                                  ind1 = ngbf1m2*kk;
                                  for(jj=0;jj<ngbf22;jj++)
                                  {
                                      ind2 = ind1 + ngbf21*(jj+1) - 1;
                                      for(dir=0;dir<ndof;dir++)
                                      {
                                           SolidListFinal[iii].dispBCs[ind2][dir] = -9999;
                                      }
                                  }
                              }
                         break;
                  }
               }
            }
        }
     }

/*
     cout << endl;
     cout << "       .... displacement values for the final mesh ... " << endl;
     cout << endl;
     for(int jjj=0;jjj<patch.n;jjj++)
     {
         cout << "       patch... : " << (jjj+1) << endl;
         cout << endl;
         for(ii=0;ii<SolidListFinal[jjj].dispBCs.n;ii++)
         {
           cout << '\t' << ii << '\t' ;
           for(jj=0;jj<SolidListFinal[jjj].dispBCs[0].n;jj++)
           {
             cout << SolidListFinal[jjj].dispBCs[ii][jj] << '\t' ;
           }
           cout << endl;
         }
         cout << endl;
     }
*/
     for(iii=0;iii<dispbc.n;iii++) // Loop A
     {
        int  ind1, ind2, ind3, patchnum, side, dir, ngbf1m2;
        int  ngbf1, ngbf2, ngbf3, ngbf;
        
        double  disp_val, angle, cs, sn, val1, val2;

        patchnum = (int) dispbc[iii][0]-1;

        ngbf1 = SolidListFinal[patchnum].ngbf1;
        ngbf2 = SolidListFinal[patchnum].ngbf2;
        ngbf3 = SolidListFinal[patchnum].ngbf3;
        ngbf  = SolidListFinal[patchnum].ngbf;

        ngbf1m2 = ngbf1*ngbf2;

        side = (int) dispbc[iii][1];
        dir  = (int) dispbc[iii][2]-1; // direction in which traction is applied ( 0 -> 1st direction, 1-> 2nd direction, 2-> 3rd direction)
        disp_val  = dispbc[iii][3];
        
        //cout << patchnum << '\t' << side << '\t' << dir << '\t' << disp_val << endl;

        if(side > 18) cerr << " side > 18 in IsogeometricFEM::ProcessDispBoundaryConditionsSolid() " << endl;

        switch(side)
        {
            case 1: // search at w = 0
            
                   ind1 = 0;
                   for(kk=0;kk<ngbf2;kk++)
                   {
                       for(jj=0;jj<ngbf1;jj++)
                       {
                          SolidListFinal[patchnum].dispBCs[ind1][dir]    = disp_val;
                          SolidListFinal[patchnum].Uinit[ndf*ind1 + dir] = disp_val;
                          ind1++;
                       }
                   }
              break;

            case 2: // search at w = 1
            
                   ind1 = ngbf1m2 * (ngbf3-1);
                   for(kk=0;kk<ngbf2;kk++)
                   {
                       for(jj=0;jj<ngbf1;jj++)
                       {
                          SolidListFinal[patchnum].dispBCs[ind1][dir]    = disp_val;
                          SolidListFinal[patchnum].Uinit[ndf*ind1 + dir] = disp_val;
                          ind1++;
                       }
                   }
              break;

            case 3: // search at v = 0
            
                   for(kk=0;kk<ngbf3;kk++)
                   {
                       ind1 = ngbf1m2 * kk;
                       for(jj=0;jj<ngbf1;jj++)
                       {
                          ind2 = ind1 + jj;
                          SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                          SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                       }
                   }
              break;

            case 4: // search at v = 1
            
                   //if(dir < 3)
                   //{
                       for(kk=0;kk<ngbf3;kk++)
                       {
                           ind1 = ngbf1m2 * (kk+1) - ngbf1 ;
                           for(jj=0;jj<ngbf1;jj++)
                           {
                              ind2 = ind1 + jj;
                              SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                              SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                           }
                       }
                   //}
                   //else
                   //{
                      /*
                       SolidListFinal[patchnum].face_rot_angle[3] = disp_val;
                       SolidResult[patchnum].face_rotation_flag = true;

                       angle = PI/180.0; // apply a rotation of 1 deg
        
                       cs = cos(angle);
                       sn = sin(angle);
                       
                       ind3 = ngbf2-1;
                       
                       EPOINT EP;

                       for(kk=0;kk<ngbf3;kk++)
                       {
                          ind1 = ngbf1m2 * (kk+1) - ngbf1 ;
                          for(jj=0;jj<ngbf1;jj++)
                          {
                              ind2 = ind1 + jj;
                              EP = SolidListFinal[patchnum].Pw[kk][jj][ind3].CalcEuclid();
               
                              //val1 = EP.x * cs - EP.z * sn - EP.x;
                              //val2 = EP.x * sn + EP.z * cs - EP.z;

                              val1 = EP.x * cs + EP.z * sn - EP.x;
                              val2 = -EP.x * sn + EP.z * cs - EP.z;

                              EP.print2screen();
                              cout << ind2 << '\t' << val1 << '\t' << val2 << endl;
                              dir=0;
                              SolidListFinal[patchnum].dispBCs[ind2][dir]    = val1;
                              SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = val1;

                              dir=2;
                              SolidListFinal[patchnum].dispBCs[ind2][dir]    = val2;
                              SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = val2;
                          }
                       }
                    */
                   //}
              break;

            case 5: // search at u = 0

                   for(kk=0;kk<ngbf3;kk++)
                   {
                       ind1 = ngbf1m2*kk;
                       for(jj=0;jj<ngbf2;jj++)
                       {
                          ind2 = ind1 + ngbf1*jj;
                          SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                          SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                       }
                   }
              break;

            case 6: // search at u = 1

                   for(kk=0;kk<ngbf3;kk++)
                   {
                       ind1 = ngbf1m2*kk;
                       for(jj=0;jj<ngbf2;jj++)
                       {
                          ind2 = ind1 + ngbf1*(jj+1) - 1;
                          SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                          SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                       }
                   }
              break;

            case 7:      // Edge 1

                 for(ii=0;ii<ngbf2;ii++)
                 {
                     ind2 = ngbf1*ii;
                     SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                     SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                 }

              break;


            case 8:        // Edge 2

                 for(ii=0;ii<ngbf1;ii++)
                 {
                      ind2 = ii;
                      SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                      SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                 }

              break;

            case 9:        // Edge 3

                 ind1 = 0;
                 for(ii=0;ii<ngbf2;ii++)
                 {
                     ind2 = ngbf1*(ii+1) - 1;
                     SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                     SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                 }

              break;

            case 10:        // Edge 4

                 ind1 = ngbf1m2 - ngbf1;
                 for(ii=0;ii<ngbf1;ii++)
                 {
                     ind2 = ind1 + ii;
                     SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                     SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                 }

              break;

            case 11:        // Edge 5

                 ind1 = ngbf1m2 * (ngbf3-1);
                 for(ii=0;ii<ngbf2;ii++)
                 {
                     ind2 = ind1 + ngbf1 * ii;
                     SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                     SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                 }

              break;

            case 12:        // Edge 6

                 ind1 = ngbf1m2 * (ngbf3-1);
                 for(ii=0;ii<ngbf1;ii++)
                 {
                     ind2 = ind1 + ii;
                     SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                     SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                 }

              break;

            case 13:        // Edge 7

                 ind1 = ngbf1m2 * (ngbf3-1);
                 for(ii=0;ii<ngbf2;ii++)
                 {
                     ind2 = ind1 + ngbf1*(ii+1) - 1;
                     SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                     SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                 }

              break;

            case 14:        // Edge 8

                 ind1 = ngbf - ngbf1;
                 for(ii=0;ii<ngbf1;ii++)
                 {
                     ind2 = ind1 + ii;
                     SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                     SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                 }

              break;

            case 15:        // Edge 9

                 ind1 = 0;
                 for(ii=0;ii<ngbf3;ii++)
                 {
                     ind2 = ind1 + ngbf1m2*ii;
                     SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                     SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                 }

              break;

            case 16:        // Edge 10

                 ind1 = ngbf1 -1;
                 for(ii=0;ii<ngbf3;ii++)
                 {
                     ind2 = ind1 + ngbf1m2*ii;
                     SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                     SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                 }

              break;

            case 17:        // Edge 11

                 ind1 = 0;
                 for(ii=0;ii<ngbf3;ii++)
                 {
                     ind2 = ind1 + ngbf1m2*(ii+1)-ngbf1;
                     SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                     SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                 }

              break;

            case 18:        // Edge 12

                 ind1 = 0;
                 for(ii=0;ii<ngbf3;ii++)
                 {
                     ind2 = ind1 + ngbf1m2*(ii+1)-1;
                     SolidListFinal[patchnum].dispBCs[ind2][dir]    = disp_val;
                     SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] = disp_val;
                 }

              break;
        }
    }



     int ngbf1, ngbf2, ngbf3, *tt;

     ntotgbf1 = 0;
     // compute total number of independent global basis funtions
     for(iii=0;iii<Npatch;iii++)
     {
        cout << " ngbf " << SolidListFinal[iii].ngbf  << endl;
     
        ntotgbf1 += SolidListFinal[iii].ngbf;

        ngbf1 = SolidListFinal[iii].ngbf1;
        ngbf2 = SolidListFinal[iii].ngbf2;
        ngbf3 = SolidListFinal[iii].ngbf3;
        
        //cout << ngbf1 << '\t' << ngbf2 << '\t' << ngbf3 << endl;
        
        tt = &(SolidListFinal[iii].edgedata[0]);

        for(kk=0;kk<6;kk++)
        {
           if(tt[kk] == 0 || tt[kk] == 1)
             ntotgbf1 -= ngbf1 * ngbf2;

           if(tt[kk] == 2 || tt[kk] == 3)
             ntotgbf1 -= ngbf1 * ngbf3;

           if(tt[kk] == 4 || tt[kk] == 5)
             ntotgbf1 -= ngbf2 * ngbf3;
        }
        
        if(tt[0] > -1 && tt[2] > -1)
          ntotgbf1 += ngbf1;
        if(tt[0] > -1 && tt[4] > -1)
          ntotgbf1 += ngbf2;
        if(tt[0] > -1 && tt[5] > -1)
          ntotgbf1 += ngbf2;
        if(tt[1] > -1 && tt[2] > -1)
          ntotgbf1 += ngbf1;
        if(tt[1] > -1 && tt[4] > -1)
          ntotgbf1 += ngbf2;
        if(tt[1] > -1 && tt[5] > -1)
          ntotgbf1 += ngbf2;

        if(tt[2] > -1 && tt[4] > -1)
          ntotgbf1 += ngbf3;
        if(tt[2] > -1 && tt[5] > -1)
          ntotgbf1 += ngbf3;

        if(tt[3] > -1 && tt[4] > -1)
          ntotgbf1 += ngbf3;
        if(tt[3] > -1 && tt[5] > -1)
          ntotgbf1 += ngbf3;
     }

/*
     cout << endl;
     cout << "       .... displacement values for the final mesh ... " << endl;
     cout << endl;
     for(int jjj=0;jjj<patch.n;jjj++)
     {
         cout << "       patch... : " << (jjj+1) << endl;
         cout << endl;
         for(ii=0;ii<SolidListFinal[jjj].dispBCs.n;ii++)
         {
           cout << '\t' << ii << '\t' ;
           for(jj=0;jj<SolidListFinal[jjj].dispBCs[0].n;jj++)
           {
             cout << SolidListFinal[jjj].dispBCs[ii][jj] << '\t' ;
           }
           cout << endl;
         }
         cout << endl;
     }
*/

  return;
}



void IsogeometricFEM::processInterfaceData()
{
  if(intfdata.n > 0)
  {
      int patchnum1, patchnum2, side1, side2;
      for(int iii=0;iii<intfdata.n;iii++)
      {
          patchnum1 = intfdata[iii][0] - 1;
          side1 = intfdata[iii][1] - 1;
          patchnum2 = intfdata[iii][2] - 1;
          side2 = intfdata[iii][3] - 1;

          NurbsBaseFinal[patchnum1]->edgedata[side1] = side1;

          NurbsBaseFinal[patchnum1]->intfdata[2*side1] = patchnum2;
          NurbsBaseFinal[patchnum1]->intfdata[2*side1+1] = side2;
      }

/*
   for(int iii=0;iii<Npatch;iii++)
   {
      cout << " Interface Data for patch #... : " << (iii+1) << endl;
      cout <<  NurbsBaseFinal[0]->edgedata << endl;
      cout <<  NurbsBaseFinal[0]->intfdata << endl;
      cout << endl;
      cout << endl;
   }
*/
  }

  return;
}



