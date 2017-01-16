#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_SUPERLU_SUPPORT
//#define EIGEN_UMFPACK_SUPPORT

#include <iostream>

#include "IsogeometricFEM.h"
#include "FunctionsProgram.h"
#include "DataBlockTemplate.h"
#include "PropertyTypeEnum.h"
#include "MathGeom.h"
#include "SolverMA41Eigen.h"
//#include "SolverPARDISO.h"
#include "NurbsShapeFns.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include <math.h>
#include "TimeFunction.h"

#include <fstream>
#include <string.h>
#include <iomanip>
#include <iostream>
#include <fstream>

//#include "GL/glut.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

//#include <Eigen/UmfPackSupport>
//#include <Eigen/SuperLUSupport>


#include "NurbsElem1DAdvectionDiffusion.h"
#include "NurbsElem1DElasticBar.h"
#include "NurbsElem1DEulerBeam.h"
#include "NurbsElem2DStructSolid.h"
#include "NurbsElem2DStructMixed2field.h"
#include "NurbsElem2DStructMixed3field.h"


using namespace std;
using namespace Eigen;

//typedef SparseMatrix<double> SparseMatrixXd;

//extern Plot plot;
extern ComputerTime computerTime;
extern MpapTime     mpapTime;
//extern Files        files;
extern List<TimeFunction> timeFunction;




inline void datafileToMatrix(MyString& fileName, MatrixXd& matx)
{
    ifstream Ffile;

    MyString line, tmpl, *word;

    int i, j, nw = 0, n = 0, n1, n2;

    double val=0.0;

    List<Vector<double> > list;

    Ffile.open(fileName.asCharArray());

    if(!Ffile)
      cerr << "failed to open the input file." << endl;

    while(!Ffile.eof())
    {
       list.add(new Vector<double>);

       line.getNextLine(Ffile);
       nw = line.split(&word);

       if(nw == 0)         break;

       i = 0;
       while(i < nw && word[i].toDbl(&val,false))
       {
          list[n].append(val);
          i++;
       }

       for(j=0; j<nw; j++)
          word[j].free();
       delete [] word;

       if(i < nw) break;
       n++;
    }

    list.del(n);

    n1 = list.n;
    n2 = list[0].n;

    matx.resize(n1,n2);

    for(i=0;i<n1;i++)
    {
       for(j=0;j<n2;j++)
          matx(i,j) = list[i][j];
    }

    Ffile.close();

  return;
}




/*
void IsogeometricFEM::ModalAnalysis()
{
    int ii, jj, e, ind2, *rr, *cc;

    double *val;

    globalK.resize(ntoteqs, ntoteqs);
    globalM.resize(ntoteqs, ntoteqs);

    globalK.setZero();
    globalM.setZero();

    /////////////////////////////////////////
    // STIFFNESS MATRIX
    /////////////////////////////////////////

    solverEigen->zeroMtx();

    for(e=0;e<totnumel;e++)
    {
       localStiffnessError = elem[e]->calcStiffnessAndResidual();

       if (localStiffnessError != 0)
       {
          cout << '\t' << "local element failure!\n\n";
          exit(0);
       }

       elem[e]->AssembleElementMatrix(1, solverEigen->mtx);
    }

    ind2 = solverEigen->mtx.x.n;

    rr   = &(solverEigen->mtx.row[0]);
    cc   = &(solverEigen->mtx.col[0]);
    val  = &(solverEigen->mtx.x[0]);

    for(ii=0;ii<ind2;ii++)
        globalK(rr[ii]-1,cc[ii]-1) = val[ii];


    /////////////////////////////////////////
    // MASS MATRIX
    /////////////////////////////////////////

    solverEigen->zeroMtx();

    for(e=0;e<totnumel;e++)
    {
       localStiffnessError = elem[e]->calcMassMatrix(1, 0.0);

       if (localStiffnessError != 0)
       {
          cout << '\t' << "local element failure!\n\n";
          exit(0);
       }

       elem[e]->AssembleElementMatrix(2, solverEigen->mtx);
    }

    ind2 = solverEigen->mtx.x.n;

    rr   = &(solverEigen->mtx.row[0]);
    cc   = &(solverEigen->mtx.col[0]);
    val  = &(solverEigen->mtx.x[0]);

    for(ii=0;ii<ind2;ii++)
        globalM(rr[ii]-1,cc[ii]-1) = val[ii];
//
    printf(" Global Stiffness matrix \n");
    printMatrix(globalK);
    printf("\n");
    printf("\n");
    printf(" Global Mass matrix \n");
    printMatrix(globalM);
    printf("\n");
//

    GeneralizedSelfAdjointEigenSolver<MatrixXd> es(globalK, globalM);

    eigen_values = es.eigenvalues();

    eigen_vectors = es.eigenvectors();

//
    cout << "The eigenvalues of the pencil (A,B) are:" << endl << eigen_values << endl;
    cout << "The matrix of eigenvectors, V, is:" << endl << eigen_vectors << endl << endl;

    double lambda = eigen_values[0];

    cout << "Consider the first eigenvalue, lambda = " << lambda << endl;

    VectorXd v = eigen_vectors.col(0);

    cout << " eigenvector v = " << endl << v << endl;
//

cout << "  Modal analysis successfully completed ... " << endl;

    cout << " The first " << min(10, (int) eigen_values.rows()) << " eigenvalues are ... " << endl;

   for(int ii=0;ii<min(10, (int) eigen_values.rows());ii++)
    cout << '\t' << (ii+1) << '\t' << sqrt(eigen_values(ii)) << endl;

//    plotFrequencyCurve(1, 0.0);


      ofstream fout("frequency-data.dat");

      if(fout.fail())
      {
	 cout << " Could not open the Output file" << endl;
	 exit(1);
      }

      fout.setf(ios::fixed);
      fout.setf(ios::showpoint);
      fout.precision(10);

      for(ii=0;ii<ntoteqs;ii++)
      {
         fout << (double) (ii+1)/ntoteqs << setw(16) << sqrt(eigen_values(ii))/PI/(ii+1) << endl;
      }

      fout.close();

    return;
}
*/


void IsogeometricFEM::ModalAnalysis()
{
    int ii, jj, e, ind2, *rr, *cc;

    double *val;
    
    if(LSFEMflag)
    {
      Vmat.resize(2*ntotgbf1, 2*ntotgbf1);
      Bmat.resize(2*ntotgbf1, ntotgbf1);
      Qmat.resize(ntotgbf1, ntotgbf1);
    
      globalK.resize(ntotgbf1, ntotgbf1);
      globalM.resize(ntotgbf1, ntotgbf1);
    }
    else
    {
      Vmat.resize(ntoteqs1, ntoteqs1);
      Bmat.resize(ntoteqs1, ntoteqs2);
      Qmat.resize(ntoteqs2, ntoteqs2);
    
      globalK.resize(ntoteqs2, ntoteqs2);
      globalM.resize(ntoteqs2, ntoteqs2);
    }

    globalK.setZero();
    globalM.setZero();
    
    Vmat.setZero();
    Bmat.setZero();
    Qmat.setZero();

    /////////////////////////////////////////
    // compute and assemble matrices
    /////////////////////////////////////////

    for(e=0;e<totnumel;e++)
    {
       //cout << "       elem... : " << (e+1) << endl;
       
       localStiffnessError = elem[e]-> toComputeInfSupCondition();

       if (localStiffnessError != 0)
       {
          cout << '\t' << "local element failure!\n\n";
          exit(0);
       }

       elem[e]->AssembleElementMatrix2(1, Vmat, globalK);
       elem[e]->AssembleElementMatrix2(2, Bmat, globalK);
       elem[e]->AssembleElementMatrix2(3, globalM, globalK);
    }

    globalK = Vmat.inverse();
    //printMatrix(Vmat);
    //printf("\n\n");
    //printMatrix(globalK);
    //printf("\n\n");
    //printMatrix(globalM);
    //printf("\n\n");
    globalK = (Bmat.transpose()*globalK)*Bmat;

/*
    printf(" Global Stiffness matrix \n");
    printMatrix(globalK);
    printf("\n");
    printf("\n");
    printf(" Global Mass matrix \n");
    printMatrix(globalM);
    printf("\n");
*/

    GeneralizedSelfAdjointEigenSolver<MatrixXd> es(globalK, globalM, EigenvaluesOnly);
 
    eigen_values = es.eigenvalues();

    //eigen_vectors = es.eigenvectors();

/*
    cout << "The eigenvalues of the pencil (A,B) are:" << endl << eigen_values << endl;
    cout << "The matrix of eigenvectors, V, is:" << endl << eigen_vectors << endl << endl;

    double lambda = eigen_values[0];

    cout << "Consider the first eigenvalue, lambda = " << lambda << endl;

    VectorXd v = eigen_vectors.col(0);

    cout << " eigenvector v = " << endl << v << endl;
*/

   cout << "  Modal analysis successfully completed ... " << endl;

   cout << " The first " << min(20, (int) eigen_values.rows()) << " eigenvalues are ... " << endl;

   for(int ii=0;ii<min(20, (int) eigen_values.rows());ii++)
     printf("\t %5d \t %12.10f \t %12.10f \n", (ii+1), eigen_values(ii), sqrt(abs(eigen_values(ii))));

    return;
}




void IsogeometricFEM::plotFrequencyCurve(int Nmode, double val)
{
    VectorXd   natural_freqs(ntoteqs), freq_ratio(ntoteqs), x_axis(ntoteqs);

    int  ii;

    double   fact ;

    for(ii=0;ii<ntoteqs;ii++)
    {
       fact = ii + 1.0;
       natural_freqs(ii) =  sqrt(eigen_values(ii));
       freq_ratio(ii)    =  natural_freqs(ii)/(fact*PI);
       x_axis(ii)        =  fact/ntoteqs;
    }


    ofstream   fout("freq-curve.dat");

    if(fout.fail())
    {
	cout << " Could not open the Output file" << endl;
	exit(1);
    }

    fout.setf(ios::fixed);
    fout.setf(ios::showpoint);
    fout.precision(10);

    for(ii=0;ii<ntoteqs;ii++)
       fout << x_axis(ii) << setw(16) << freq_ratio(ii) << endl;

    fout.close();

    return;
}


void IsogeometricFEM::plotModeShape(int Nmode, int flag)
{
    cout << "        Mode shape ... : " << Nmode+1 << endl;

    VectorXd  vec = eigen_vectors.col(Nmode);

    soln.zero();

    int ii, resln[3];
    double  fact = abs(vec.maxCoeff()) * 10.0;
    
    resln[0] = 1;    resln[1] = 1;    resln[2] = 1;

    for(ii=0;ii<ntoteqs;ii++)
      soln[assy4r[ii]] = vec(ii)/fact;

    if(patchGrp[0].ndom == 1)
    {
        //plot.setColour(1);

        CurveResult[0].Pw = CurveListFinal[0].Pw;

        CurveResult[0].updateCoordinates(&(soln[0]));

        CurveResult[0].PlotElements(1,0,resln);
    }
    else
    {
        resln[0] = resln[1] = resln[2] = 5;
    
        SurfaceResult[0].Pw = SurfaceListFinal[0].Pw;

        SurfaceResult[0].updateCoordinates(&(soln[0]));
        
        //NurbsBaseResult[0]->Values = soln;

        //NurbsBaseResult[0]->PlotElements(3, 1, resln);

        SurfaceResult[0].PlotValues(1);
    }

//    plotGeom(3,1);

	return;
}





int IsogeometricFEM::finalsolve()
{
//   char fct[] = "IsogeometricFEM::finalsolve";

   ModalAnalysis();

  //cout << " ERROR in finalsolve() for mixed problems using Eigen Solver " << endl;

  return 0;
}




int IsogeometricFEM::calcAndAssyTangentMatrix(bool flag, double dt)
{
       if(intVarFlag)  copyElemInternalVariables();

       //flag is true by default
       solverEigen->rhsVec.setZero();
       if(flag)  //calculate the new tangent matrix
       {
          //initialize the matrix
          solverEigen->mtx.setZero();

          for(int e=0;e<totnumel;e++)
          {
              if(e<0)
              {
                cout << "        elem # : " << e << endl;
                cout << endl;

                int nivEl = elem[e]->nivGP * elem[e]->nGP;
                for (int i=0; i<nivEl; i++)
                  cout <<  '\t' << i << '\t' << elem[e]->intVar1[i] << '\t' << elem[e]->intVar2[i] << endl; //fixed <<
              }

              localStiffnessError = elem[e]->calcStiffnessMatrix(dt);

              if (localStiffnessError != 0)
              {
                  cout << "local element failure --> tangent matrix calculation. Elem # : " << e << endl;
                  return localStiffnessError;
              }
              //if(e==0)
              // elem[e]->printStiffnessMatrix();   cout << endl;    cout << endl;

              elem[e]->AssembleElementMatrix(1, solverEigen->mtx);
              elem[e]->AssembleElementVector(firstIter, 0, &(solverEigen->rhsVec[0]), &(reac[0]), 0, 0);
          }
       }

  return 0;
}




int IsogeometricFEM::calcAndAssyInternalForceVector(double dt)
{
   cout << " firstIter = " << firstIter << endl;
   
   if(intVarFlag)
     copyElemInternalVariables();

   solverEigen->rhsVec.setZero();
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

  if(firstIter)
  {
    int iii, ii, kk, ind1, ind2, ee;

    for(iii=0;iii<patch.n;iii++)
      NurbsBaseResult[iii]->addInitDOFvalues();

    for(ee=0;ee<totnumel;ee++)
      elem[ee]->initialiseDOFvalues();

    for(ii=0;ii<solnFull.n;ii++)
      solnFull[ii] += mpapTime.dt * solnInit[ii];
  }


   for(int e=0;e<totnumel;e++)  // loop over all the elements
   {
      localStiffnessError = elem[e]->calcInternalForces();

      if (localStiffnessError != 0)
      {
         cout << '\t' << "local element failure!\n\n";
         return localStiffnessError;
      }
      
      //elem[e]->printForceVector();

      elem[e]->AssembleElementVector(firstIter, 0, &(solverEigen->rhsVec[0]), &(reac[0]), start1, start2);
   }
   
   firstIter = false;

//     printData(4,0);
//     printData(3,0);

//     addExternalForces();

  return 0;
}




void IsogeometricFEM::staticAnalysis(int nropt, int niter, double dt)
{


  return;
}




void IsogeometricFEM::quasiStaticAnalysis(int nropt, int niter, double dt)
{


  return;
}






void IsogeometricFEM::elementDiffStiffTest(double ddd, int elnum, int dig, int dig2, bool gfrmt)
{
  if(elnum > totnumel)
    cout << "  Error! Requested Element number is out of range " << endl;
  else
  {
     cout << endl;
     cout << endl;
     cout << "             ///////////////////////////////////// " << endl;
     cout << "             // " << endl;
     cout << "             //     diffStiffTest for Element # : " << elnum << endl;
     cout << "             // " << endl;
     cout << "             ///////////////////////////////////// " << endl;
     cout << endl;
     cout << endl;
     elem[elnum-1]->diffStiffTest(ddd, dig, dig2, 1);
     cout << endl;
     cout << endl;
  }

  return;
}







double IsogeometricFEM::computeReactions(int a1, int a2)
{
   char fct[] = " IsogeometricFEM::computeReactions ";

   int  jj, dof, patchnum, side, *UU;
   double   val = 0.0;


    patchnum = (int) outdparam[a2][0]-1;
    dof      = (int) outdparam[a2][2] - 1;
    side     = outdparam[a2][4] ;
    
    if(patchGrp[0].ndom == 2)
    {
         int  ngbf21, ngbf22, index1, index2;

        ngbf21 = SurfaceListFinal[patchnum].ngbf1;
        ngbf22 = SurfaceListFinal[patchnum].ngbf2;

        UU = &(SurfaceListFinal[patchnum].toUfull[0]);

        switch(side)
        {
            case 1: // top side

              for(jj=0;jj<ngbf22;jj++)
              {
                  index2 = ngbf21*jj*ndf + dof;
                  val += reac[UU[index2]];
              }

              break;

            case 2: // bottom side

              for(jj=0;jj<ngbf22;jj++)
              {
                  index2 = (ngbf21*(jj+1)-1)*ndf + dof;
                  val += reac[UU[index2]];
              }

              break;

            case 3: // left side

              for(jj=0;jj<ngbf21;jj++)
              {
                  index2 = jj*ndf + dof;
                  val += reac[UU[index2]];
              }

              break;

            case 4:  // right side

              index1 = ngbf21*(ngbf22-1);
              for(jj=0;jj<ngbf21;jj++)
              {
                  index2 = (index1+jj) * ndf + dof;
                  val += reac[UU[index2]];
              }

              break;

           default: prgError(1,fct," invalid patch side number");

              break;

        }
    }

    if(patchGrp[0].ndom == 3)
    {
         int  ngbf1, ngbf2, ngbf3, ind1, ind2, ind3, kk, ngbf1m2;

        ngbf1 = SolidListFinal[patchnum].ngbf1;
        ngbf2 = SolidListFinal[patchnum].ngbf2;
        ngbf3 = SolidListFinal[patchnum].ngbf3;
        
        ngbf1m2 = ngbf1*ngbf2;

        //cout << SolidListFinal[patchnum].toUfull << endl;

        UU = &(SolidListFinal[patchnum].toUfull[0]);
        
        if(a1 == 1) //reaction force
        {
            switch(side)
            {
                case 1: // Face #1
            
                    ind1=0;
                    for(kk=0;kk<ngbf2;kk++)
                    {
                         for(jj=0;jj<ngbf1;jj++)
                         {
                             ind2 = ind1++;
                             //cout << '\t' << ind2 << '\t' << UU[ind2] << endl;
                             ind3 = ind2*ndf + dof;
                             val += reac[UU[ind3]];
                         }
                    }

                  break;

                case 2: // Face #2
            
                    ind1=ngbf1*ngbf2*(ngbf3-1);
                    for(kk=0;kk<ngbf2;kk++)
                    {
                         for(jj=0;jj<ngbf1;jj++)
                         {
                             ind2 = ind1++;
                             //cout << '\t' << ind2 << '\t' << UU[ind2] << endl;
                             ind3 = ind2*ndf + dof;
                             val += reac[UU[ind3]];
                         }
                    }

                  break;

               default: prgError(1,fct," invalid patch side number");

                  break;
            } //end of switch
        } // end of if(a1 == 1)

        if(a1 == 2) //reaction moments
        {
            EPOINT  pos_vec, force_vec, moment_vec_local, moment_full(0.0,0.0,0.0);
        
            switch(side)
            {
                case 3: // Face #3
            
                   for(kk=0;kk<ngbf3;kk++)
                   {
                       ind1 = ngbf1m2 * kk;
                       for(jj=0;jj<ngbf1;jj++)
                       {
                           ind2 = ind1 + jj;
                           
                           pos_vec = SolidListFinal[patchnum].Pw[kk][jj][0].CalcEuclid();

                           ind3 = ind2*ndf;
                           
                           force_vec.x = reac[UU[ind3]];
                           force_vec.y = reac[UU[ind3+1]];
                           force_vec.z = reac[UU[ind3+2]];
                           
                           CrossProduct(pos_vec, force_vec, moment_vec_local);
                           
                           moment_full = moment_full + moment_vec_local;
                           
                       }
                   }
                   
                   if(dof == 1)
                     val = moment_full.x;
                   if(dof == 2)
                     val = moment_full.y;
                   if(dof == 3)
                     val = moment_full.z;


                  break;


               default: prgError(1,fct," invalid patch side number");

                  break;
            } //end of switch
        } // end of if(a1 == 1)

    }
//       cout << "          reaction  ... : " << val << endl; cout << endl;
//       cout << "          outdtype  ... : " << outdtype[0] << endl; cout << endl;

  return val;
}


void IsogeometricFEM::ProcessTractionBCs()
{
  if(ndm == 1)
    ProcessTractionBCsCurves();
  else if(ndm == 2)
  {
    ProcessTractionBCsSurface();
    ProcessTractionBCsSurface2();
  }
  else
    ProcessTractionBCsSolid();

  return;
}



void IsogeometricFEM::ProcessTractionBCsCurves()
{
    cout << "     ISOGEOMETRICFEM: processing Traction boundary conditions for curves ...\n\n";

    char fct[] = "IsogeometricFEM::ProcessTractionBCsForCurves";

    int ngbf1=0, ngbf2=0, ngbf11=0, ngbf12=0, ngbf21=0, ngbf22=0,  count1=0, start=0;
    double  spec_val;

    CPOINT CP1, CP2, CP3;

    int index, patchnum, val1, val2, dir, iii, ii, jj, kk, nelm;

    bool left, right, top, bottom;

    for(iii=0;iii<tracbc.n;iii++) // Loop A
    {
        //cout << tracbc[iii] << endl;
    
        patchnum  =  (int) tracbc[iii][0]-1; // patch number
        val1      =  (int) tracbc[iii][1]-1; // end
        dir       =  (int) tracbc[iii][2]-1; // dof
        spec_val  =  tracbc[iii][3];

        nelm = CurveListFinal[patchnum].nelem;

        left = right = false;

        // first search if the traction is applied on a complete side of either of the boundaries
        // if YES, apply the same amount of traction in the refined mesh on complete side of the specified boundary
        // if NO, check for the coordinate locations of the specified CPs and find their indices in the refined mesh
        // boundaries are as per the Surface Control Net matrix indices

        //cout << "       val1 && val2 .... : " << val1 << "   " << val2 << endl;
        //cout << "       ngbf21 && ngbf22 .... : " << ngbf21 << "   " << ngbf22 << endl;
        
        // search at the left boundary
        if( val1 == 0 )
        {
            left  = true;
            index = 0;
            elem[index]->tracflag = true;
            elem[index]->createTractionDataVariable();
            elem[index]->tracdata[0][dir] = spec_val;
        }
        // search at the right boundary
        else if(val1 == 1)
        {
            right = true;
            index = nelm-1;
            elem[index]->tracflag = true;
            elem[index]->createTractionDataVariable();
            elem[index]->tracdata[1][dir] = spec_val;
        }
        else
        {
           cerr << " wrong side in IsogeometricFEM::ProcessTractionBCsForCurves " << endl;
        }
    }

  return;
}




void IsogeometricFEM::ProcessTractionBCsSurface()
{
  cout << "     ISOGEOMETRICFEM: processing Traction boundary conditions ...\n\n";

  char fct[] = "IsogeometricFEM::ProcessTractionBoundaryConditions";
  //computerTime.go(fct);

    int ngbf1=0, ngbf2=0, ngbf11=0, ngbf12=0, ngbf21=0, ngbf22=0,  count1=0, start=0;
    double  spec_val;

    CPOINT CP1, CP2, CP3;

    int index, patchnum, val1, val2, dir, iii, ii, jj, kk, nelm1, nelm2;

    bool left, right, top, bottom;

    for(iii=0;iii<tracbc.n;iii++) // Loop A
    {
        cout << tracbc[iii] << endl;
    
        patchnum  =  (int) tracbc[iii][0]-1;
        val1      =  (int) tracbc[iii][1]-1;
        val2      =  (int) tracbc[iii][2]-1;
        dir       =  (int) tracbc[iii][3]-1; // direction in which traction is applied ( 0 -> 1st direction, 1-> 2nd direction)
        spec_val  =  tracbc[iii][4];

        ngbf11 = SurfaceListOriginal[patchnum].ngbf1;
        ngbf12 = SurfaceListOriginal[patchnum].ngbf2;

        ngbf21 = SurfaceListFinal[patchnum].ngbf1;
        ngbf22 = SurfaceListFinal[patchnum].ngbf2;

        nelm1 = SurfaceListFinal[patchnum].nelem1;
        nelm2 = SurfaceListFinal[patchnum].nelem2;

        left = right = top = bottom = false;


        // first search if the traction is applied on a complete side of either of the boundaries
        // if YES, apply the same amount of traction in the refined mesh on complete side of the specified boundary
        // if NO, check for the coordinate locations of the specified CPs and find their indices in the refined mesh
        // boundaries are as per the Surface Control Net matrix indices

        cout << "       val1 && val2 .... : " << val1 << "   " << val2 << endl;
        cout << "       ngbf11 && ngbf12 .... : " << ngbf11 << "   " << ngbf12 << endl;
        cout << "       ngbf21 && ngbf22 .... : " << ngbf21 << "   " << ngbf22 << endl;

        start = 0;
        for(ii=0;ii<patchnum;ii++)
          start += SurfaceListFinal[ii].nelem;

        // search at the top boundary
        if( (val1 == 0 && val2 == ngbf11*(ngbf12-1)) || (val2 == 0 && val1 == ngbf11*(ngbf12-1)) )
        {
             cout << " top boundary " << endl;
             top = true;
             for(jj=0;jj<nelm2;jj++)
             {
                index = start+nelm1*jj;

                elem[index]->tracflag = true;
                elem[index]->createTractionDataVariable();
                elem[index]->tracdata[3][dir] = spec_val;
             }
        }

        // search at the bottom boundary
        else if( ((val1 == ngbf11-1) && (val2 == ngbf11*ngbf12-1)) || ((val2 == ngbf11-1) && (val1 == ngbf11*ngbf12-1)) )
        {
            cout << " bottom boundary " << endl;
            bottom = true;
             for(jj=0;jj<nelm2;jj++)
             {
                index = start+nelm1*(jj+1)-1;
                elem[index]->tracflag = true;
                elem[index]->createTractionDataVariable();
                elem[index]->tracdata[1][dir] = spec_val;

             //cout << " AAAAAAAAAAAAAAA " << endl;
             }
        }

        // search at the left boundary
        else if( ((val1 == 0) && (val2 == ngbf11-1)) || ((val2 == 0) && (val1 == ngbf11-1)) )
        {
            cout << " left boundary " << endl;
            left = true;
             for(jj=0;jj<nelm1;jj++)
             {
                index = start+jj;
                elem[index]->tracflag = true;
                elem[index]->createTractionDataVariable();
                elem[index]->tracdata[0][dir] = spec_val;
             }
        }

        // search at the right boundary
        else if( ((val1 == ngbf11*(ngbf12-1)) && (val2 == ngbf11*ngbf12-1)) || ((val2 == ngbf11*(ngbf12-1)) && (val1 == ngbf11*ngbf12-1)) )
        {
            cout << " right boundary " << endl;
            right = true;
             for(jj=0;jj<nelm1;jj++)
             {
                index = start+nelm1*(nelm2-1)+jj;
                elem[index]->tracflag = true;
                elem[index]->createTractionDataVariable();
                elem[index]->tracdata[2][dir] = spec_val;
             }
        }

        SurfaceListFinal[patchnum].tracflag[0] = top;
        SurfaceListFinal[patchnum].tracflag[1] = bottom;
        SurfaceListFinal[patchnum].tracflag[2] = left;
        SurfaceListFinal[patchnum].tracflag[3] = right;

    } // end of for loop A


//    cout << '\t' << top << '\t' << bottom << '\t' << left << '\t' << right << endl;

//    cout << "  tracflag " << endl;
//    cout << '\t' << SurfaceListFinal[0].tracflag << endl;
/*
    Vector<int> temp1;
    for(iii=0;iii<Npatch;iii++)
    {
        int ngbf21, ngbf22;

        ngbf21 = SurfaceListFinal[iii].ngbf1;
        ngbf22 = SurfaceListFinal[iii].ngbf2;

       // compute the global basis numbers on which force is applied, to assemble to 'rhsVec' vector

       Vector<int> temp2;
       if(SurfaceListFinal[iii].tracflag[0])
       {
          for(ii=0;ii<ngbf22;ii++)
          {
              temp2.append(ngbf21 * ii);
          }
       }
       if(SurfaceListFinal[iii].tracflag[1])
       {
          for(ii=0;ii<ngbf22;ii++)
          {
              temp2.append(ngbf21*(ii+1)-1);
          }
       }
       if(SurfaceListFinal[iii].tracflag[2])
       {
          for(ii=0;ii<ngbf21;ii++)
          {
              temp2.append(ii);
          }
       }
       if(SurfaceListFinal[iii].tracflag[3])
       {
          int ttt = ngbf21*(ngbf22-1);
          for(ii=0;ii<ngbf21;ii++)
          {
              temp2.append(ttt+ii);
          }
       }

      cout << endl;  cout << '\t' << temp2 << endl;

       if(temp2.n > 0)
       {
          VectorArray<int> temp3, temp4;

          temp3.setDim(2*temp2.n);
          int ind1, ind2;
          for(ii=0;ii<temp2.n;ii++)
          {
             ind1 = ndf * ii;
             //ind2 = temp2[ii]*ndf;
             //for(jj=0;jj<SurfaceListFinal[iii].ID.n;jj++)
             for(jj=0;jj<2;jj++)
             {
                //temp3[ind1+jj] = SurfaceListFinal[iii].toUfull[ind2+jj];
                temp3[ind1+jj] = SurfaceListFinal[iii].ID[jj][temp2[ii]];
             }
          }

          SortArrayInt(temp3);
      //    cout << endl;     cout << '\t' << temp3 << endl;

          finduniqueInt(temp3, temp4);
          cout << endl;     cout << '\t' << temp4 << endl;

          if(temp4[0] == -1)
          {
             for(ii=1;ii<temp4.n;ii++)
               temp1.append(temp4[ii]);
          }
          else
          {
             for(ii=0;ii<temp4.n;ii++)
               temp1.append(temp4[ii]);
          }

       }
    }

    assy4F = temp1;

    cout << endl;     cout << '\t' << assy4F << endl;     cout << endl;
*/

/*
     cout << "       .... Traction BC data for the elements in final mesh ... " << endl;
     cout << endl;

  for(int e=0;e<totnumel;e++)
  {
    cout << '\t' << e << '\t' << elem[e]->tracflag << endl;
  }
    cout  << endl;
    cout  << endl;

  for(int e=0;e<totnumel;e++)
  {
    cout << '\t' << e << endl;
    for(int ii=0;ii<elem[e]->tracdata.n;ii++)
    {
      for(int jj=0;jj<elem[e]->tracdata[0].n;jj++)
      {
        cout << "   " << elem[e]->tracdata[ii][jj];
      }
    cout  << endl;
    }
    cout  << endl;
    cout  << endl;
  }
*/

//  computerTime.stopAndPrint(fct);

  return;
}




void IsogeometricFEM::ProcessTractionBCsSurface2()
{
  cout << "     ISOGEOMETRICFEM: processing Traction boundary conditions 2 ...\n\n";
  char fct[] = "IsogeometricFEM::ProcessTractionBoundaryConditions";
  //computerTime.go(fct);

    int ngbf21, ngbf22, iii, ii, jj, kk;

    int index, patchnum, side, ttt, index1, index2;
    int val1, val2;

    CPOINT CP1, CP2;

    Vector<int> temp1;

    for(iii=0;iii<tracbc2.n;iii++)
    {
        patchnum = tracbc2[iii][0]-1;
            side = tracbc2[iii][1]-1;
            val1 = tracbc2[iii][2]-1;
            val2 = tracbc2[iii][3]-1;

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

        index1 = index2 = 0;

        switch(side)
        {
          case 0:  // search at the top boundary

                 for(jj=1;jj<(ngbf21-1);jj++)
                 {
                    if(SurfaceListFinal[patchnum].Pw[0][jj] == CP1)
                       index1 = jj;

                    if(SurfaceListFinal[patchnum].Pw[0][jj] == CP2)
                       index2 = jj;
                 }

                 if(index != 0 && index2 != 0)
                 {
                    for(kk=0;kk<index1;kk++)
                      temp1.append(ngbf21*kk);
                    for(kk=(index2+1);kk<ngbf21;kk++)
                      temp1.append(ngbf21*kk);
                 }


                 break;

          case 1:  // search at the bottom boundary

                 for(jj=1;jj<(ngbf21-1);jj++)
                 {
                    if(SurfaceListFinal[patchnum].Pw[ngbf21-1][jj] == CP1)
                       index1 = jj;

                    if(SurfaceListFinal[patchnum].Pw[ngbf21-1][jj] == CP2)
                       index2 = jj;
                 }

                 if(index != 0 && index2 != 0)
                 {
                    for(kk=0;kk<index1;kk++)
                      temp1.append(ngbf21*(kk+1) - 1);
                    for(kk=(index2+1);kk<ngbf21;kk++)
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

                 if(index != 0 && index2 != 0)
                 {
                    for(kk=0;kk<index1;kk++)
                      temp1.append(kk);
                    for(kk=(index2+1);kk<ngbf21;kk++)
                      temp1.append(kk);
                 }

                 break;

          case 3:  // search at the right boundary

                 ttt = ngbf21*(ngbf22-1);
                 for(jj=1;jj<(ngbf21-1);jj++)
                 {
                    if(SurfaceListFinal[patchnum].Pw[jj][ngbf22-1] == CP1)
                       index1 = jj;

                    if(SurfaceListFinal[patchnum].Pw[jj][ngbf22-1] == CP2)
                       index2 = jj;
                 }

                 if(index != 0 && index2 != 0)
                 {
                    for(kk=0;kk<index1;kk++)
                      temp1.append(ttt+kk);
                    for(kk=(index2+1);kk<ngbf21;kk++)
                      temp1.append(ttt+kk);
                 }

                 break;
        }

    }

  //cout << '\t' << index1 << '\t' << index2 << endl;
  //cout << '\t' << temp1 << endl;
  //cout << endl;     cout << '\t' << " assy4F " << assy4F << endl;

  if(temp1.n > 0)
  {
     VectorArray<int> temp3, temp2;

     temp2.setDim(2*temp1.n);
     int ind;
     for(ii=0;ii<temp1.n;ii++)
     {
        ind = ndf * ii;
        for(jj=0;jj<SurfaceListFinal[patchnum].ID.n;jj++)
        {
           temp2[ind+jj] = SurfaceListFinal[patchnum].ID[jj][temp1[ii]];
        }
     }
  //  cout << endl;     cout << '\t' << " temp2 " << temp2 << endl;

    finduniqueInt(temp2, temp3);
  //  cout << endl;     cout << '\t' << " temp3 " << temp3 << endl;

    sub2VectorArraysInt(assy4F, temp3, assy4F2);

  }
  else
    assy4F2 = assy4F;

  //cout << endl;     cout << '\t' << " assy4F2 " << assy4F2 << endl;

  cout << "     ISOGEOMETRICFEM: processing Traction boundary conditions 2 ...DONE \n\n";

  return;
}



void  IsogeometricFEM::ProcessBCsConstraintVariable()
{
  cout << "     ISOGEOMETRICFEM: processing disp BCs for Surfaces for Constraint Variable...\n\n";

  char fct[] = "IsogeometricFEM::ProcessBCsConstraintVariable";

  int ii, jj, iii;

     for(iii=0;iii<ibc4.n;iii++) // Loop A
     {
        int index, patchnum, side, dir;
        int ngbf1, ngbf2, ngbf, temp, ndof;

        patchnum = ibc4[iii][0]-1;

        ngbf1 = surfSecondVar[patchnum].ngbf1;
        ngbf2 = surfSecondVar[patchnum].ngbf2;
        ngbf  = surfSecondVar[patchnum].ngbf;
        ndof  = surfSecondVar[patchnum].ndof;

        side = ibc4[iii][1]-1;

        dir = 0;

      //  dir  = ibc4[iii][2]-1; //( 0 -> 1st direction, 1-> 2nd direction)


        switch(side)
        {
             case 0:   // search at the top boundary

                    for(jj=0;jj<ngbf2;jj++)
                    {
                       index = ngbf1*jj;
                       surfSecondVar[patchnum].dispBCs[index][dir] = 0.0;
                       surfSecondVar[patchnum].Uinit[ndof*index + dir] = 0.0;
                    }

             break;

             case 1:        // search at the bottom boundary

                    for(jj=0;jj<ngbf2;jj++)
                    {
                       index = (ngbf1*(jj+1)-1);
                       surfSecondVar[patchnum].dispBCs[index][dir] = 0.0;
                       surfSecondVar[patchnum].Uinit[ndof*index + dir] = 0.0;
                    }

             break;

             case 2:        // search at the left boundary

                    for(jj=0;jj<ngbf1;jj++)
                    {
                       index = jj;
                       surfSecondVar[patchnum].dispBCs[index][dir] = 0.0;
                       surfSecondVar[patchnum].Uinit[ndof*index + dir] = 0.0;
                    }

             break;

             case 3:        // search at the right boundary

                    temp = ngbf1*(ngbf2-1);
                    for(jj=0;jj<ngbf1;jj++)
                    {
                       index = temp + jj;
                       surfSecondVar[patchnum].dispBCs[index][dir] = 0.0;
                       surfSecondVar[patchnum].Uinit[ndof*index + dir] = 0.0;
                    }

             break;
        }
    }

  return;
}






void IsogeometricFEM::GenerateConnectivityArraysSolid()
{
     cout << "     ISOGEOMETRICFEM: generating cinnectivity arrays for solid ...\n\n";
     char fct[] = "IsogeometricFEM::GenerateConnectivityArraysSolid";

     int  iii, ii, jj, kk, count;

     ntoteqs1=0;
     for(iii=0;iii<SolidListFinal.n;iii++)
       SolidListFinal[iii].GenerateConnectivityArrays1(ntoteqs1);

     cout << '\t' << " ntoteqs1  " << ntoteqs1 << endl;

     count = 0;
     SolidListFinal[0].computeGBFnumbers(1, count);


     // adjust ID arrays based on interface edges and also compute global basis function numbers

     if(Npatch > 1)
     {
        int  patchnum, side, ngbf21, ngbf22, ngbf23, ngbf31, ngbf32, ngbf33, ind4, ind2, dir, ll;

        for(iii=1;iii<Npatch;iii++)
        {
            ngbf21 = SolidListFinal[iii].ngbf1;
            ngbf22 = SolidListFinal[iii].ngbf2;
            ngbf23 = SolidListFinal[iii].ngbf3;
 
            for(ll=0;ll<6;ll++)
            {
               if(SolidListFinal[iii].edgedata[ll] != -1)
               {
                  patchnum = SolidListFinal[iii].intfdata[2*ll];
                  side     = SolidListFinal[iii].intfdata[2*ll+1];

                  ngbf31 = SolidListFinal[patchnum].ngbf1;
                  ngbf32 = SolidListFinal[patchnum].ngbf2;
                  ngbf33 = SolidListFinal[patchnum].ngbf3;

                  //cout << '\t' << ngbf21 << '\t' << ngbf22 << '\t' << ngbf23 << endl;
                  //cout << '\t' << ngbf31 << '\t' << ngbf32 << '\t' << ngbf33 << endl;

                  switch(ll)
                  {
                       case 0: // search at w = 0

                               if( (ngbf21 != ngbf31) && (ngbf22 != ngbf32) )
                                  prgError(1,fct," ngbf1 and ngbf2 in the 2 patches do not match");

                              for(kk=0;kk<ngbf22;kk++)
                              {
                                  for(jj=0;jj<ngbf21;jj++)
                                  {
                                      ind2 = SolidListFinal[iii].gbfNumFace1(jj, kk);
                                      ind4 = SolidListFinal[patchnum].gbfNumFace2(jj, kk);

                                      SolidListFinal[iii].gbfnums[ind2] = SolidListFinal[patchnum].gbfnums[ind4];

                                      for(dir=0;dir<ndf;dir++)
                                      {
                                         SolidListFinal[iii].ID[dir][ind2]         =  SolidListFinal[patchnum].ID[dir][ind4];
                                         SolidListFinal[iii].dispBCs[ind2][dir]    =  SolidListFinal[patchnum].dispBCs[ind4][dir];

                                         SolidListFinal[iii].Uinit[ndf*ind2 + dir] =  SolidListFinal[patchnum].Uinit[ndf*ind4 + dir] ;
                                      }
                                  }
                              }
                         break;

                       case 1: // search at w = 1

                               if( (ngbf21 != ngbf31) && (ngbf22 != ngbf32) )
                                  prgError(1,fct," ngbf1 and ngbf2 in the 2 patches do not match");

                              for(kk=0;kk<ngbf22;kk++)
                              {
                                  for(jj=0;jj<ngbf21;jj++)
                                  {
                                      ind2 = SolidListFinal[iii].gbfNumFace1(jj, kk);
                                      ind4 = SolidListFinal[patchnum].gbfNumFace2(jj, kk);

                                      SolidListFinal[iii].gbfnums[ind2] = SolidListFinal[patchnum].gbfnums[ind4];

                                      for(dir=0;dir<ndf;dir++)
                                      {
                                         SolidListFinal[iii].ID[dir][ind2]         =  SolidListFinal[patchnum].ID[dir][ind4];
                                         SolidListFinal[iii].dispBCs[ind2][dir]    =  SolidListFinal[patchnum].dispBCs[ind4][dir];

                                         SolidListFinal[iii].Uinit[ndf*ind2 + dir] =  SolidListFinal[patchnum].Uinit[ndf*ind4 + dir] ;
                                      }
                                  }
                              }
                         break;

                       case 2: // search at v = 0

                              if( (ngbf21 != ngbf31) && (ngbf23 != ngbf33) )
                                  prgError(1,fct," ngbf1 and ngbf3 in the 2 patches do not match");

                              for(kk=0;kk<ngbf23;kk++)
                              {
                                  for(jj=0;jj<ngbf21;jj++)
                                  {
                                      ind2 = SolidListFinal[iii].gbfNumFace3(jj, kk);
                                      ind4 = SolidListFinal[patchnum].gbfNumFace4(jj, kk);

                                      SolidListFinal[iii].gbfnums[ind2] = SolidListFinal[patchnum].gbfnums[ind4];

                                      for(dir=0;dir<ndf;dir++)
                                      {
                                         SolidListFinal[iii].ID[dir][ind2]         =  SolidListFinal[patchnum].ID[dir][ind4];
                                         SolidListFinal[iii].dispBCs[ind2][dir]    =  SolidListFinal[patchnum].dispBCs[ind4][dir];

                                         SolidListFinal[iii].Uinit[ndf*ind2 + dir] =  SolidListFinal[patchnum].Uinit[ndf*ind4 + dir] ;
                                      }
                                  }
                              }
                         break;

                       case 3: // search at v = 1
            
                              if( (ngbf21 != ngbf31) && (ngbf23 != ngbf33) )
                                  prgError(1,fct," ngbf1 and ngbf3 in the 2 patches do not match");

                              for(kk=0;kk<ngbf23;kk++)
                              {
                                  for(jj=0;jj<ngbf21;jj++)
                                  {
                                      ind2 = SolidListFinal[iii].gbfNumFace3(jj, kk);
                                      ind4 = SolidListFinal[patchnum].gbfNumFace4(jj, kk);

                                      SolidListFinal[iii].gbfnums[ind2] = SolidListFinal[patchnum].gbfnums[ind4];

                                      for(dir=0;dir<ndf;dir++)
                                      {
                                         SolidListFinal[iii].ID[dir][ind2]         =  SolidListFinal[patchnum].ID[dir][ind4];
                                         SolidListFinal[iii].dispBCs[ind2][dir]    =  SolidListFinal[patchnum].dispBCs[ind4][dir];

                                         SolidListFinal[iii].Uinit[ndf*ind2 + dir] =  SolidListFinal[patchnum].Uinit[ndf*ind4 + dir] ;
                                      }
                                  }
                              }

                         break;

                      case 4: // 

                               if( (ngbf22 != ngbf32) && (ngbf23 != ngbf33) )
                                  prgError(1,fct," ngbf2 and ngbf3 in the 2 patches do not match");

                               for(kk=0;kk<ngbf22;kk++)
                               {
                                   for(jj=0;jj<ngbf23;jj++)
                                   {
                                      ind2 = SolidListFinal[iii].gbfNumFace5(jj, kk);
                                      ind4 = SolidListFinal[patchnum].gbfNumFace6(jj, kk);

                                      SolidListFinal[iii].gbfnums[ind2] = SolidListFinal[patchnum].gbfnums[ind4];

                                      for(dir=0;dir<ndf;dir++)
                                      {
                                         SolidListFinal[iii].ID[dir][ind2]         =  SolidListFinal[patchnum].ID[dir][ind4];
                                         SolidListFinal[iii].dispBCs[ind2][dir]    =  SolidListFinal[patchnum].dispBCs[ind4][dir];

                                         SolidListFinal[iii].Uinit[ndf*ind2 + dir] =  SolidListFinal[patchnum].Uinit[ndf*ind4 + dir] ;
                                      }
                                   }
                               }

                          break;

                       case 5: // search at u = 1

                               if( (ngbf22 != ngbf32) && (ngbf23 != ngbf33) )
                                  prgError(1,fct," ngbf2 and ngbf3 in the 2 patches do not match");

                               for(kk=0;kk<ngbf22;kk++)
                               {
                                   for(jj=0;jj<ngbf23;jj++)
                                   {
                                      ind2 = SolidListFinal[iii].gbfNumFace5(jj, kk);
                                      ind4 = SolidListFinal[patchnum].gbfNumFace6(jj, kk);

                                      SolidListFinal[iii].gbfnums[ind2] = SolidListFinal[patchnum].gbfnums[ind4];

                                      for(dir=0;dir<ndf;dir++)
                                      {
                                         SolidListFinal[iii].ID[dir][ind4]         =  SolidListFinal[patchnum].ID[dir][ind2];
                                         SolidListFinal[iii].dispBCs[ind4][dir]    =  SolidListFinal[patchnum].dispBCs[ind2][dir];

                                         SolidListFinal[iii].Uinit[ndf*ind4 + dir] =  SolidListFinal[patchnum].Uinit[ndf*ind2 + dir] ;
                                      }
                                   }
                               }
                         break;

                  }
               }
            }
            SolidListFinal[iii].computeGBFnumbers(2, count);
        }
     }


     // calculate LM Arrays
     for(iii=0;iii<Npatch;iii++)
        SolidListFinal[iii].GenerateConnectivityArrays2();

     //for(iii=0;iii<Npatch;iii++)
       //SolidListFinal[iii].printConnectivityArrays();

     int count1, *U4;
     double *U3;

     for(iii=0;iii<Npatch;iii++)
        SolidListFinal[iii].computeToUfull();

/*
     for(iii=0;iii<Npatch;iii++)
     {
         cout << SolidListFinal[iii].gbfnums << endl;
         cout << SolidListFinal[iii].toUfull << endl;
         cout << endl;
         cout << endl;
     }
*/

     // set global 'Uinit' vector
     solnInit.zero();
     count1 = -1;
     for(iii=0;iii<Npatch;iii++)
     {
         U4 = &(SolidListFinal[iii].toUfull[0]);
         U3 = &(SolidListFinal[iii].Uinit[0]);
         for(ii=0;ii<SolidListFinal[iii].Uinit.n;ii++)
         {
            if(U4[ii] > count1)
               solnInit[U4[ii]] = U3[ii];
         }
         count1 += (SolidListFinal[iii].ngbf * ndf);
     }


     for(iii=0;iii<Npatch;iii++)
     {
        SolidResult[iii].toUfull = SolidListFinal[iii].toUfull;
        SolidResult[iii].gbfnums = SolidListFinal[iii].gbfnums;
        SolidResult[iii].Uinit = SolidListFinal[iii].Uinit;

        //SolidListFinal[iii].toUfull.free();
        //SolidListFinal[iii].gbfnums.free();
        //SolidListFinal[iii].Uinit.free();
     }

     cout << "     ISOGEOMETRICFEM: generating cinnectivity arrays for solid ...DONE ... \n\n";

  return;
}




void IsogeometricFEM::ProcessTractionBCsSolid()
{
    cout << "     ISOGEOMETRICFEM: processing Traction boundary conditions for SOLIDs...\n\n";

    int  ngbf1, ngbf2, ngbf3, ngbf, count1=0, start=0;
    
    int  nelm1, nelm2, nelm3, nelm1m2, nelm, ngbf1m2, ind1, ind2;

    CPOINT CP1, CP2, CP3;

    int  index, patchnum, type, dir, iii, ii, jj, kk;

    bool left, right, top, bottom, front, back, edges[12];
    
    double  value;

    for(iii=0;iii<tracbc.n;iii++) // Loop A
    {
        patchnum = (int) tracbc[iii][0]-1;
        type     = (int) tracbc[iii][1];
        dir      = (int) tracbc[iii][2]-1;
        value    = tracbc[iii][3];

        ngbf1 = SolidListFinal[patchnum].ngbf1;
        ngbf2 = SolidListFinal[patchnum].ngbf2;
        ngbf3 = SolidListFinal[patchnum].ngbf3;
        ngbf1m2 = SolidListFinal[patchnum].ngbf1m2;
        ngbf  = SolidListFinal[patchnum].ngbf;

        nelm1 = SolidListFinal[patchnum].nelem1;
        nelm2 = SolidListFinal[patchnum].nelem2;
        nelm3 = SolidListFinal[patchnum].nelem3;
        nelm  = SolidListFinal[patchnum].nelem;
        nelm1m2 = SolidListFinal[patchnum].nelem1m2;
        
        //cout << " Traction value " << value << endl;

        left = right = top = bottom = front = back = false;
        for(ii=0;ii<12;ii++)
          edges[ii] = false;

        cout << "       ngbf .... : " << ngbf1 << '\t' << ngbf2 << '\t' << ngbf3 << endl;

        //cout << "       val1 && val2 .... : " << patchnum << '\t' << type << '\t' << dir << '\t' << value << endl;
//        cout << "       nelem .... : " << nelm1 << '\t' << nelm2 << '\t' << nelm3 << endl;

        start = 0;
        for(ii=0;ii<patchnum;ii++)
          start += SolidListFinal[ii].nelem;

        // front boundary
        if(type == 1)
        {
             front = true;
             ind1 = start;
             for(jj=0;jj<nelm2;jj++)
             {
                for(ii=0;ii<nelm1;ii++)
                {
                   ind2 = ind1++;
                   elem[ind2]->tracflag = true;
                   elem[ind2]->createTractionDataVariable();
                   elem[ind2]->tracdata[0][dir] = value;
                }
             }
        }
        
        // back boundary
        if(type == 2)
        {
             back = true;
             ind1 = start + nelm1m2*(nelm3-1);
             for(jj=0;jj<nelm2;jj++)
             {
                for(ii=0;ii<nelm1;ii++)
                {
                   ind2 = ind1++;
                   elem[ind2]->tracflag = true;
                   elem[ind2]->createTractionDataVariable();
                   elem[ind2]->tracdata[1][dir] = value;
                }
             }
        }

        // bottom boundary
        if(type == 3)
        {
             bottom = true;
             for(jj=0;jj<nelm3;jj++)
             {
                ind1 = start + nelm1m2*jj;
                for(ii=0;ii<nelm1;ii++)
                {
                   ind2 = ind1 + ii;
                   elem[ind2]->tracflag = true;
                   elem[ind2]->createTractionDataVariable();
                   elem[ind2]->tracdata[2][dir] = value;
                }
             }
        }

        // top boundary
        if(type == 4)
        {
             top = true;
             for(jj=0;jj<nelm3;jj++)
             {
                ind1 = start + nelm1m2*(jj+1) - nelm1;
                for(ii=0;ii<nelm1;ii++)
                {
                   ind2 = ind1 + ii;
                   elem[ind2]->tracflag = true;
                   elem[ind2]->createTractionDataVariable();
                   elem[ind2]->tracdata[3][dir] = value;
                }
             }
        }

        // left boundary
        if(type == 5)
        {
             left = true;
             for(jj=0;jj<nelm2;jj++)
             {
                ind1 = start + nelm1*jj;
                for(ii=0;ii<nelm3;ii++)
                {
                   ind2 = ind1 + nelm1m2*ii;
                   elem[ind2]->tracflag = true;
                   elem[ind2]->createTractionDataVariable();
                   elem[ind2]->tracdata[4][dir] = value;
                }
             }
        }
        // right boundary
        if(type == 6)
        {
             right = true;
             for(jj=0;jj<nelm2;jj++)
             {
                ind1 = start + nelm1*(jj+1)-1;
                for(ii=0;ii<nelm3;ii++)
                {
                   ind2 = ind1 + nelm1m2*ii;
                   elem[ind2]->tracflag = true;
                   elem[ind2]->createTractionDataVariable();
                   elem[ind2]->tracdata[5][dir] = value;
                }
             }
        }

        // Edge 1
        if(type == 7)
        {
             edges[0] = true;
             ind1 = start;
             for(ii=0;ii<nelm2;ii++)
             {
                 ind2 = ind1 + nelm1*ii;
                 elem[ind2]->tracflag = true;
                 elem[ind2]->createTractionDataVariable();
                 elem[ind2]->tracdata[6][dir] = value;
             }
        }
        // Edge 2
        if(type == 8)
        {
             edges[1] = true;
             ind1 = start;
             for(ii=0;ii<nelm1;ii++)
             {
                 ind2 = ind1 + ii;
                 elem[ind2]->tracflag = true;
                 elem[ind2]->createTractionDataVariable();
                 elem[ind2]->tracdata[7][dir] = value;
             }
        }
        // Edge 3
        if(type == 9)
        {
             edges[2] = true;
             ind1 = start ;
             for(ii=0;ii<nelm2;ii++)
             {
                 ind2 = ind1 + nelm1*(ii+1) - 1;
                 elem[ind2]->tracflag = true;
                 elem[ind2]->createTractionDataVariable();
                 elem[ind2]->tracdata[8][dir] = value;
             }
        }

        // Edge 4
        if(type == 10)
        {
             edges[3] = true;
             ind1 = start + nelm1m2 - nelm1;
             for(ii=0;ii<nelm1;ii++)
             {
                 ind2 = ind1 + ii;
                 elem[ind2]->tracflag = true;
                 elem[ind2]->createTractionDataVariable();
                 elem[ind2]->tracdata[9][dir] = value;
             }
        }

        // Edge 5
        if(type == 11)
        {
             edges[4] = true;
             ind1 = start + nelm1m2 * (nelm3-1);
             for(ii=0;ii<nelm2;ii++)
             {
                 ind2 = ind1 + nelm1 * ii;
                 elem[ind2]->tracflag = true;
                 elem[ind2]->createTractionDataVariable();
                 elem[ind2]->tracdata[10][dir] = value;
             }
        }

        // Edge 6
        if(type == 12)
        {
             edges[5] = true;
             ind1 = start + nelm1m2 * (nelm3-1);
             for(ii=0;ii<nelm1;ii++)
             {
                 ind2 = ind1 + ii;
                 elem[ind2]->tracflag = true;
                 elem[ind2]->createTractionDataVariable();
                 elem[ind2]->tracdata[11][dir] = value;
             }
        }
        // Edge 7
        if(type == 13)
        {
             edges[6] = true;
             ind1 = start + nelm1m2 * (nelm3-1);
             for(ii=0;ii<nelm2;ii++)
             {
                 ind2 = ind1 + nelm1*(ii+1) - 1;
                 elem[ind2]->tracflag = true;
                 elem[ind2]->createTractionDataVariable();
                 elem[ind2]->tracdata[12][dir] = value;
             }
        }

        // Edge 8
        if(type == 14)
        {
             edges[7] = true;
             ind1 = start + nelm - nelm1;
             for(ii=0;ii<nelm1;ii++)
             {
                 ind2 = ind1 + ii;
                 elem[ind2]->tracflag = true;
                 elem[ind2]->createTractionDataVariable();
                 elem[ind2]->tracdata[13][dir] = value;
             }
        }

        // Edge 9
        if(type == 15)
        {
             edges[8] = true;
             ind1 = start;
             for(ii=0;ii<nelm3;ii++)
             {
                 ind2 = ind1 + nelm1m2*ii;
                 elem[ind2]->tracflag = true;
                 elem[ind2]->createTractionDataVariable();
                 elem[ind2]->tracdata[14][dir] = value;
             }
        }

        // Edge 10
        if(type == 16)
        {
             edges[9] = true;
             ind1 = start + nelm1 -1;
             for(ii=0;ii<nelm3;ii++)
             {
                 ind2 = ind1 + nelm1m2*ii;
                 elem[ind2]->tracflag = true;
                 elem[ind2]->createTractionDataVariable();
                 elem[ind2]->tracdata[15][dir] = value;
             }
        }

        // Edge 11
        if(type == 17)
        {
             edges[10] = true;
             ind1 = start;
             for(ii=0;ii<nelm3;ii++)
             {
                 ind2 = ind1 + nelm1m2*(ii+1)-nelm1;
                 elem[ind2]->tracflag = true;
                 elem[ind2]->createTractionDataVariable();
                 elem[ind2]->tracdata[16][dir] = value;
             }
        }

        // Edge 12
        if(type == 18)
        {
             edges[11] = true;
             ind1 = start;
             for(ii=0;ii<nelm3;ii++)
             {
                 ind2 = ind1 + nelm1m2*(ii+1)-1;
                 elem[ind2]->tracflag = true;
                 elem[ind2]->createTractionDataVariable();
                 elem[ind2]->tracdata[17][dir] = value;
             }
        }

        SolidListFinal[patchnum].tracflag[0] = front;
        SolidListFinal[patchnum].tracflag[1] = back;
        SolidListFinal[patchnum].tracflag[2] = bottom;
        SolidListFinal[patchnum].tracflag[3] = top;
        SolidListFinal[patchnum].tracflag[4] = left;
        SolidListFinal[patchnum].tracflag[5] = right;
        for(ii=0;ii<12;ii++)
           SolidListFinal[patchnum].tracflag[6+ii] = edges[ii];
        

    } // end of for loop A

/*
    for(ii=0;ii<totnumel;ii++)
       cout << ii << " ---> " << elem[ii]->tracflag << endl;

    for(ii=0;ii<totnumel;ii++)
    {
       cout << " elem ...: " << (ii+1) << endl;
       if(elem[ii]->tracdata.n > 0)
       {
          for(jj=0;jj<18;jj++)
             cout << '\t' << elem[ii]->tracdata[jj] ;
          cout << endl;
       }
       else
          cout << " ....... NO TRACTION FORCE ON THIS ELEMENT " << endl;
       cout << endl;
    }
*/

    cout << "  SolidListFinal  tracflag " << endl;
    for(iii=0;iii<patch.n;iii++)
      cout << '\t' << SolidListFinal[iii].tracflag << endl;

    Vector<int> temp1;
    for(iii=0;iii<patch.n;iii++)
    {
        ngbf1   = SolidListFinal[iii].ngbf1;
        ngbf2   = SolidListFinal[iii].ngbf2;
        ngbf3   = SolidListFinal[iii].ngbf3;
        ngbf1m2 = SolidListFinal[iii].ngbf1m2;
        ngbf    = SolidListFinal[iii].ngbf;

        Vector<int> temp2;
        if(SolidListFinal[iii].tracflag[0])
        {
           ind1 = 0;
           for(jj=0;jj<ngbf2;jj++)
           {
              for(ii=0;ii<ngbf1;ii++)
                 temp2.append(ind1++);
           }
        }
        if(SolidListFinal[iii].tracflag[1])
        {
           ind1 = ngbf1m2*(ngbf3-1);
           for(jj=0;jj<ngbf2;jj++)
           {
              for(ii=0;ii<ngbf1;ii++)
                 temp2.append(ind1++);
           }
        }
        if(SolidListFinal[iii].tracflag[2])
        {
           for(jj=0;jj<ngbf3;jj++)
           {
              ind1 = ngbf1m2*jj;
              for(ii=0;ii<ngbf1;ii++)
                 temp2.append(ind1+ii);
           }
        }
        if(SolidListFinal[iii].tracflag[3])
        {
           for(jj=0;jj<ngbf3;jj++)
           {
              ind1 = ngbf1m2*(jj+1) - ngbf1;
              for(ii=0;ii<ngbf1;ii++)
                 temp2.append(ind1+ii);
           }
        }
        if(SolidListFinal[iii].tracflag[4])
        {
           for(jj=0;jj<ngbf2;jj++)
           {
              ind1 = ngbf1*jj;
              for(ii=0;ii<ngbf3;ii++)
                 temp2.append(ind1 + ngbf1m2*ii);
           }
        }
        if(SolidListFinal[iii].tracflag[5])
        {
           for(jj=0;jj<ngbf2;jj++)
           {
              ind1 = ngbf1*(jj+1)-1;
              for(ii=0;ii<ngbf3;ii++)
                 temp2.append(ind1 + ngbf1m2*ii);
           }
        }

//       cout << endl;  cout << '\t' << temp2 << endl;

        if(temp2.n > 0)
        {
            VectorArray<int> temp3, temp4;
 
            temp3.setDim(3*temp2.n);

            for(ii=0;ii<temp2.n;ii++)
            {
               ind1 = ndf * ii;

               for(jj=0;jj<SolidListFinal[iii].ID.n;jj++)
               {
                  temp3[ind1+jj] = SolidListFinal[iii].ID[jj][temp2[ii]];
               }
            }

            SortArrayInt(temp3);
      //    cout << endl;     cout << '\t' << temp3 << endl;

            finduniqueInt(temp3, temp4);
      //    cout << endl;     cout << '\t' << temp4 << endl;

            if(temp4[0] == -1)
            {
               for(ii=1;ii<temp4.n;ii++)
                 temp1.append(temp4[ii]);
            }
            else
            {
               for(ii=0;ii<temp4.n;ii++)
                 temp1.append(temp4[ii]);
            }
        }
    }

    assy4F = temp1;
    
//    cout << " assy4F " << assy4F << endl;

  return;
}

/*
double  IsogeometricFEM::LineSearch(double gtol, VectorXd& prsd, VectorXd& du, double STOL)
{
     double  sa, ga, gb, sb, g, step, temp;
     int  j, linmax = 10;

           // Beginning of line search
           /////////////////////////////////////////////////////

           sa = 1.0;
           ga = gamma1(du, sa);
           sb = 0.0;
           gb = gtol;

           //     Find bracket on zero

           if(ga*gb > 0.0)
           {
               printf(" -> No line search - Energy positive at full step ... ");
           }
           else
           {
              temp = STOL*abs(gtol);

              //     Perform 'linmax' steps of line-search
              for(j=0; j <= linmax ; j++)
              {
                  step = sa - ga*(sa-sb)/(ga-gb);
                  g    = gamma1(du, step);

                  printf(" \t\t-> Iteration = %5d, Step Size = %12.5e, Energy = %12.5e,  sa= %8.4f, sb=%8.4f\n", j, step, g, sa, sb);

                  if( abs(g) < temp )
                    break;

                  //   Update postions for next iteration

                  gb = 0.5*gb;
                  if(g*ga < 0.0)
                  {
                    sb = sa;
                    gb = ga;
                  }
                  sa = step;
                  ga = g;
               }
          }

   return step;
}



double  IsogeometricFEM::LineSearch(double gtol, double* residue_old, double* du, double STOL)
{
     double  sa, ga, gb, sb, g, step, temp, step_old;
     int  j, linmax = 10, kk, iii;

     VectorArray<double>  geom_temp; geom_temp.setDim(ntotgbf1*ndf); geom_temp.zero();

     for(iii=0;iii<Npatch;iii++)
        NurbsBaseResult[iii]->geomToVector(&(geom_temp[0]));

           // Beginning of line search
           /////////////////////////////////////////////////////
           
          for(kk=0;kk<ntoteqs;kk++)
             uinter[assy4r[kk]] = du[kk];
             
          // update the geometry and compute new residuals

          for(iii=0;iii<Npatch;iii++)
             NurbsBaseResult[iii]->updateCoordinates(&(uinter[0]));

          sa = 1.0;
          ga = gamma1(&(du[0]), sa);
          sb = 0.0;
          gb = gtol;

          for(iii=0;iii<Npatch;iii++)
             NurbsBaseResult[iii]->resetGeometry(&(geom_temp[0]));


           printf("    Line search begins now -----------> gtol= %12.6e, ga=%12.6e \n\n", gtol, ga);

           step = 1.0;

           if(ga*gb > 0.0)
             printf(" -----------> No line search - Energy positive at full step ... \n");
           else
           {
              temp = STOL*abs(gtol);

              //     Perform 'linmax' steps of line-search
              
              for(j=0; j <= linmax ; j++)
              {
                  step_old = step;
                  step = sa - ga*(sa-sb)/(ga-gb);

                  for(kk=0;kk<ntoteqs;kk++)
                    uinter[assy4r[kk]] = step * du[kk];

                  for(iii=0;iii<Npatch;iii++)
                     NurbsBaseResult[iii]->updateCoordinates(&(uinter[0]));

                  g = gamma1(du, 1.0);

                  for(iii=0;iii<patch.n;iii++)
                     NurbsBaseResult[iii]->resetGeometry(&(geom_temp[0]));

                  printf(" \t\t-> Iteration = %5d, Step Size = %12.6f, Norm = %12.5e,  sa= %8.4f, sb=%8.4f\n", j, step, g, sa, sb);

                  //   Update postions for next iteration

                  gb = 0.5*gb;
                  if( (ga*g) < 0.0)
                  {
                    sb = sa;
                    gb = ga;
                  }
                  sa = step;
                  ga = g;

                  if( abs(g) < temp )
                    break;

               }
          }

          for(iii=0;iii<patch.n;iii++)
              NurbsBaseResult[iii]->resetGeometry(&(geom_temp[0]));

          printf(" \t\t-> Step size = %12.6f, \n", step);

   return step;
}
*/



double  IsogeometricFEM::LineSearch(double old, double* residue_old, double* du, double STOL)
{
     double  sa, ga, gb, sb, g, step, temp, step_old, val_new, dummy, val_old, step_out;
     int  j, linmax = 10, kk, iii;

     VectorArray<double>  disp_incr, geom_temp;
     
     ListArray<VectorArray<double> > secVarPrev_temp; secVarPrev_temp.setDim(Npatch);
    
     disp_incr.setDim(ntoteqs); disp_incr.zero();
     geom_temp.setDim(ntotgbf1*ndf); geom_temp.zero();

     for(iii=0;iii<Npatch;iii++)
        NurbsBaseResult[iii]->geomToVector(&(geom_temp[0]));
     if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
     {
        for(iii=0;iii<Npatch;iii++)
          secVarPrev_temp[iii]  =  NurbsBaseSecondVar[iii]->Values[0] ;
     }

     // Beginning of line search
     /////////////////////////////////////////////////////
           
     for(kk=0;kk<ntoteqs1;kk++)
        soln[assy4r[kk]] = du[kk];
             
     // update the geometry and compute new residuals

     for(iii=0;iii<Npatch;iii++)
        NurbsBaseResult[iii]->updateCoordinates(&(soln[0]));

     sa = 1.0;
     ga = gamma1(&(du[0]), sa);
     gb = dotProductVecs(residue_old, du, ntoteqs);
     sb = 0.0;
     //gb = gtol;
     
     val_old = old;
     val_new = solverEigen->rhsVec.norm();

     temp = STOL*abs(gb);

     for(iii=0;iii<Npatch;iii++)
        NurbsBaseResult[iii]->resetGeometry(&(geom_temp[0]));
     if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
     {
         for(iii=0;iii<Npatch;iii++)
            NurbsBaseSecondVar[iii]->Values[0] = secVarPrev_temp[iii];
     }


     printf("    Line search begins now -----------> Norm old= %12.6e,  Norm new= %12.6e,  ga= %12.6e, gb=%12.6e , temp=%12.6e \n\n", val_old, val_new, ga, gb, temp);

     step_out = step = 1.0;

           //if(ga*gb > 0.0)
             //printf(" -----------> No line search - Energy positive at full step ... \n");
           //else
           {
              //     Perform 'linmax' steps of line-search
              
              for(j=0; j <= linmax ; j++)
              {
                  step_old = step;

                  //printf(" \n\t\t sa = %12.6f, sb = %12.6f,  ga = %12.6e,  gb= %12.6e \n", sa, sb, ga, gb);

                  //step = sa - ga*(sa-sb)/(ga-gb);
                  
                  step = sa - val_new*(sa-sb)/(val_new-val_old);
                  /*
                  if(val_new < val_old)
                  {
                     step += 0.5;
                  }
                  else
                  {
                    step -= 0.5;
                  }
                  */
                  for(kk=0;kk<ntoteqs;kk++)
                    disp_incr[kk] = step * du[kk];

                  for(kk=0;kk<ntoteqs1;kk++)
                    soln[assy4r[kk]] = disp_incr[kk];

                  for(iii=0;iii<Npatch;iii++)
                     NurbsBaseResult[iii]->updateCoordinates(&(soln[0]));

                  if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
                  {
                     for(iii=0;iii<Npatch;iii++)
                        NurbsBaseSecondVar[iii]->updateValues(1, &disp_incr[ntoteqs1]);
                  }

                  g = gamma1(du, 1.0);
                  
                  dummy = solverEigen->rhsVec.norm();
                  
                  val_old = val_new;
                  val_new = solverEigen->rhsVec.norm();
                  sb = sa;
                  sa = step;


                  for(iii=0;iii<patch.n;iii++)
                     NurbsBaseResult[iii]->resetGeometry(&(geom_temp[0]));

                  if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
                  {
                     for(iii=0;iii<Npatch;iii++)
                        NurbsBaseSecondVar[iii]->Values[0] = secVarPrev_temp[iii];
                  }

                  printf(" \t\t-> Iteration = %5d, Step Size = %12.6f, Norm = %12.5e,  g = %12.5e,  sa= %8.4f, sb=%8.4f\n", j, step, solverEigen->rhsVec.norm(), g, sa, sb);

                  //   Update postions for next iteration
                  
                  if( val_new < old )
                    step_out = step;
                  
                  /*

                  //gb = 0.5*gb;
                  if( val_new > val_old )
                  {
                    sb = sa;
                    gb = ga;
                  }
                  sa = step;
                  ga = g;
                  */

                  if( abs(g) < temp )
                    break;

               }
          }

          for(iii=0;iii<patch.n;iii++)
              NurbsBaseResult[iii]->resetGeometry(&(geom_temp[0]));

          printf(" \t\t-> Step size = %12.6f, \n", step_out);

   return step_out;
}



double  IsogeometricFEM::gamma1(double* du, double step)
{
    double  val = 0.0;

    calcAndAssyInternalForceVector(1.0);

    for(int kk=0;kk<ntoteqs;kk++)
       val += (step * du[kk] * solverEigen->rhsVec[kk]);

    return  val;
}




int IsogeometricFEM::BFGSsolver(int nBFGS, bool LINE_SEARCH_FLAG, double STOL, double TOL)
{
    int ii, jj, iii, *rr, *cc, kk, ind2, linmax = 10;
    bool  UP=true;

    int nn=ntotgbf1*ndf;

    double  temp, gtol, step, g0, g, initNorm, curNorm, step_old, coeff;
    
    double  fact, delgam, dlkdl, stcond, condmx = 100000;

    ListArray<VectorArray<double> >   w_vec, v_vec;
    
    w_vec.setDim(nBFGS);
    v_vec.setDim(nBFGS);
    
    for(ii=0;ii<nBFGS;ii++)
    {
      w_vec[ii].setDim(ntoteqs);
      v_vec[ii].setDim(ntoteqs);
      w_vec[ii].zero();
      v_vec[ii].zero();
    }
    
    //VectorArray<double>  gamma, residue_new, disp_incr, residue_old, geom_temp;
    
    //gamma.setDim(ntoteqs);
    //gamma.zero();

    //geom_temp.setDim(ntotgbf1*ndf);
    //geom_temp.zero();
    
    VectorXd  gamma, residue_new, disp_incr, residue_old, geom_temp;

    gamma.resize(ntoteqs);
    gamma.setZero();

    residue_new = gamma;
    disp_incr   = gamma;
    residue_old = gamma;

    geom_temp.resize(ntotgbf1*ndf);
    geom_temp.setZero();

    solverEigen->currentStatus = ASSEMBLY_OK;
    solverEigen->factorise();

    step = 1.0;
    g0   = 0.0;
    //g = gamma1(&(disp_incr[0]), step);
    
    cout << " AAAAAAAAAA " << nBFGS << '\t' << LINE_SEARCH_FLAG << '\t' << STOL << '\t' << TOL << endl;

    //firstIter = true;
    calcAndAssyInternalForceVector(1.0);
    //calcStiffnessAndResidual();
    
    cout << " AAAAAAAAAA " << endl;

    //cout << " residue_new " << endl;        printVector(&(solverEigen->rhsVec[0]), rhsVec.n);
    
    initNorm = solverEigen->rhsVec.norm();
    
    rNormPrev = rNorm = initNorm;

    residue_new = solverEigen->rhsVec;
    residue_old = residue_new;
    
    for(ii=0;ii<nBFGS;ii++)
    {
        //cout << " residue_new " << endl;      printVector(residue_new);

        solverEigen->solve();

        for(kk=0;kk<ntoteqs;kk++)
           disp_incr[kk] = solverEigen->soln[kk];

        if(ii > 0)
        {
           for(jj=0;jj<ii;jj++)
           {
              coeff = dotProductVecs(&(v_vec[jj][0]), &(disp_incr[0]), ntoteqs);
              for(kk=0;kk<ntoteqs;kk++)
                 disp_incr[kk] += (coeff * w_vec[jj][kk] );
           }
        }
        
        //cout << " disp_incr " << endl;        printVector(&(disp_incr[0]), disp_incr.n);
        
        // perform line search

        if(LINE_SEARCH_FLAG)
        {
             //if(abs(g) > (STOL * abs(g0)) )

             //for(iii=0;iii<patch.n;iii++)
               // NurbsBaseResult[iii]->resetGeometry(&(geom_temp[0]));
               
             //g0 = dotProductVecs(&(residue_new[0]), &(disp_incr[0]), ntoteqs);

             step = LineSearch(residue_new.norm(), &(residue_new[0]), &(disp_incr[0]), STOL);

             for(kk=0;kk<ntoteqs;kk++)
                disp_incr[kk] = step * disp_incr[kk];
        }

        for(kk=0;kk<ntoteqs1;kk++)
           soln[assy4r[kk]] = disp_incr[kk];

        for(iii=0;iii<Npatch;iii++)
           NurbsBaseResult[iii]->updateCoordinates(&(soln[0]));

        if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
        {
           for(iii=0;iii<Npatch;iii++)
              NurbsBaseSecondVar[iii]->updateValues(1, &disp_incr[ntoteqs1]);
        }

        calcAndAssyInternalForceVector(1.0);
        //calcStiffnessAndResidual();
        
        //cout << " residue_new " << endl;        printVector(&(solverEigen->rhsVec[0]), rhsVec.n);

        residue_old = residue_new;
        residue_new = solverEigen->rhsVec;

        rNormPrev = rNorm;
        rNorm     = solverEigen->rhsVec.norm();

        COUT << domain.name(this); printf("%5d\t\t  %11.4e\n",ii, rNorm);

        //if( (rNorm/initNorm) < TOL ) break;
        if( rNorm < TOL )      break;
        
        //if( rNorm > rNormPrev ) break;

        for(kk=0;kk<ntoteqs;kk++)
           gamma[kk] = residue_old[kk] - residue_new[kk] ;

        temp = dotProductVecs(&(disp_incr[0]), &(gamma[0]), ntoteqs);
        
        cout << " temp " << temp << endl;

        coeff = 1.0/temp;

        for(kk=0;kk<ntoteqs;kk++)
           w_vec[ii][kk] = coeff * disp_incr[kk];

        temp /= (dotProductVecs(&(disp_incr[0]), &(residue_old[0]), ntoteqs));

        temp = sqrt(abs(temp));

        for(kk=0;kk<ntoteqs;kk++)
          v_vec[ii][kk] = - (temp * residue_old[kk]) - gamma[kk];

        for(jj=ii;jj>=0;jj--)
        {
          coeff = dotProductVecs(&(w_vec[jj][0]), &(solverEigen->rhsVec[0]), ntoteqs);
          for(kk=0;kk<ntoteqs;kk++)
            solverEigen->rhsVec[kk] += (coeff * v_vec[jj][kk] );
        }
   }


  return 0;
}



