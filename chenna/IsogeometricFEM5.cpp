
#include "IsogeometricFEM.h"
#include "FunctionsProgram.h"
#include "DataBlockTemplate.h"
#include "PropertyTypeEnum.h"
#include "MathGeom.h"
#include "SolverMA41Eigen.h"
#include "SolverPardisoEigen.h"
#include "NurbsShapeFns.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "NurbsUtilitiesSOLID.h"
#include "PlotVTK.h"

#include <set>

#include "util.h"

extern DomainTree         domain;
extern List<TimeFunction> timeFunction;
extern MpapTime           mpapTime;
extern ComputerTime       computerTime;
extern PlotVTK  plotvtk;


using namespace std;


void IsogeometricFEM::setSolver(int slv, int *parm, bool cIO)
{
   if(solverEigen != NULL)
     delete solverEigen;
   solverEigen = NULL;

    char fct[] = "HBSplineFEM::setSolver";
    
    int numProc;
    
    SOLVER_TYPE = slv;
    
    switch(slv)
    {
      case  4: // SolverEigen ..........................

             solverEigen = new SolverEigen;

             comprMtxFlg = false;

             printInfo();

             prepareMatrixPattern();

             cout << " kkkkkkkkkk " << endl;
             if(solverEigen->initialise(0, 0, ntoteqs) != 0)
                return;

             //solverEigen->SetSolverAndParameters();
             solverEigen->setAlgorithmType(parm[0]);

             solverEigen->printInfo();

             break;

        case  5: // PARDISO(sym) with Eigen
        case  6: // PARDISO(unsym) with Eigen

            //SOLVER_TYPE = SOLVER_TYPE_EIGEN;

            solverEigen = (SolverEigen*) new SolverPardisoEigen;

            if (parm != NULL) numProc = parm[0]; else numProc = 1;

            numProc = min(MAX_PROCESSORS,numProc);

            //cout << " numProc " <<  numProc << endl;

            solverEigen->STABILISED = true;

            //printInfo();
            //cout << " numProc " <<  numProc << endl;
            prepareMatrixPattern();

            if(slv == 5)
            {
              if(solverEigen->initialise(numProc, PARDISO_STRUCT_SYM, ntoteqs) != 0)
                return;
            }
            if(slv == 6)
            {
              if(solverEigen->initialise(numProc, PARDISO_UNSYM, ntoteqs) != 0)
                return;
            }

            solverEigen->printInfo();

        break;


    default: // invalid slv ...................

             cout << " this solver has not been implemented yet!\n\n";

             break;
    }

    solverOK = true;
    
    //cout << " cIO " << cIO << endl;

    if(solverEigen != NULL)
      solverEigen->checkIO = cIO;

    //if(tis > 0)
      //setInitialConditions();

    return;
}


void IsogeometricFEM::prepareMatrixPattern()
{
  prepareMatrixPattern1();
  cout << " jjjjjjjjjjjjjjj " << endl;
  prepareMatrixPattern3();
  cout << " jjjjjjjjjjjjjjj " << endl;
  
  if(SOLVER_TYPE >= 4)
  {
    if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
      secVarPrev.setDim(Npatch);

    if(mixedSolverFlag == 8)
    {
      secVarPrev.setDim(Npatch);
      secVarPrev2.setDim(Npatch);
    }
    
    prepareMatrixPatternPetsc();
  }
  else
    prepareMatrixPattern2();
  
  return;
}


void IsogeometricFEM::prepareMatrixPattern1()
{
    // generate connectivity arrays
    
    cout << ndm << endl;
    
    ProcessDispBoundaryConditions();

    //uinter.setDim(ntotgbf1*ndf+ntotgbf2);    uinter.zero();

    soln.setDim(ntotgbf1*ndf);
    soln.zero();
    
    solnFull = soln;
    solnInit = soln;
    solnPrev = soln;
    diff     = soln;
    reac     = soln;

    GenerateConnectivityArrays();
    cout << " jjjjjjjjjjjjjjj " << endl;
    //printData(5,0);
    ProcessTractionBCs();
    cout << " jjjjjjjjjjjjjjj " << endl;

    if(mixedSolverFlag > 0)
      GenerateConnectivityArraysConstraintVariables();

    ntoteqs = ntoteqs1;
    ntotgbf = ntotgbf1;

    if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
    {
        ntoteqs = ntoteqs1 + ntoteqs2;
        ntotgbf = ntotgbf1 + ntotgbf2;
    }
    if(mixedSolverFlag == 8)
    {
        ntoteqs = ntoteqs1 + 2*ntoteqs2;
        ntotgbf = ntotgbf1 + 2*ntotgbf2;
    }

    cout << " ntoteqs1, ntoteqs2 and ntoteqs....: " << ntoteqs1  << '\t' << ntoteqs2 << '\t' << ntoteqs << endl;
    cout << " ntotgbf1, ntotgbf2 and ntotgbf....: " << ntotgbf1  << '\t' << ntotgbf2 << '\t' << ntotgbf << endl;
    cout << endl;

    solverEigen->rhsVec.resize(ntoteqs);
    solverEigen->rhsVec.setZero();

    ForceVec.setDim(ntoteqs);
    ForceVec.zero();

    assy4r.setDim(ntoteqs);
    assy4r.zero();

    for(int ee=0;ee<totnumel;ee++)
      elem[ee]->initialiseDOFvalues();

/*
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

    return;
}






void IsogeometricFEM::prepareMatrixPattern3()
{
    int  ee, r, c, r1, c1, count=0, count1=0, count2=0, ii, jj, iii, e, val1, val2, n1, n2, kk, e1, a, b, ll, pp;
    int *tt, *tt1, *tt2;

   cout << " PPPPPPPPPPPPP " << endl;

     if(patchGrp[0].ndom == 1)
     {
        count2 = 0;
        for(ii=0;ii<CurveListFinal[0].ngbf;ii++)
        {
           for(jj=0;jj<ndf;jj++)
           {
              if( CurveListFinal[0].ID[jj][ii] != -1)
                 assy4r[count2++] = ii*ndf + jj;
           }
        }

     }
     else if(patchGrp[0].ndom == 2)
     {
        count2 = 0;
        int  ngbf1 = SurfaceListFinal[0].ngbf1;
        for(ii=0;ii<SurfaceListFinal[0].ngbf;ii++)
        {
           if( (ii+1)%ngbf1 == 0 )
           {
              for(jj=0;jj<ndf;jj++)
              {
                 if( SurfaceListFinal[0].ID[jj][ii] != -1 && (SurfaceListFinal[0].ID[jj][ii] != SurfaceListFinal[0].ID[jj][ii-ngbf1+1]) )
                    assy4r[count2++] = ii*ndf + jj;
              }
           }
           else
           {
              for(jj=0;jj<ndf;jj++)
              {
                 if( SurfaceListFinal[0].ID[jj][ii] != -1)
                    assy4r[count2++] = ii*ndf + jj;
              }
           }
        }
     }
     else if(patchGrp[0].ndom == 3)
     {
        count2 = 0;
        for(ii=0;ii<SolidListFinal[0].ngbf;ii++)
        {
           for(jj=0;jj<ndf;jj++)
           {
              if( SolidListFinal[0].ID[jj][ii] != -1)
                 assy4r[count2++] = ii*ndf + jj;
           }
        }
     }

     if(Npatch > 1)
     {
        int  *tt, ngbf, index1;
        for(iii=1;iii<Npatch;iii++)
        {
           ngbf = NurbsBaseFinal[iii-1]->ngbf;
           count1 = NurbsBaseFinal[iii-1]->gbfnums[ngbf-1];
           tt = &(NurbsBaseFinal[iii]->gbfnums[0]);

           for(ii=0;ii<NurbsBaseFinal[iii]->ngbf;ii++)
           {
              if(tt[ii] > count1)
              {
                 index1 = tt[ii]*ndf;
                 for(jj=0;jj<ndf;jj++)
                 {
                    if( NurbsBaseFinal[iii]->ID[jj][ii] != -1)
                    assy4r[count2++] = index1 + jj;
                 }
              }
           }
        }
     }


/*
         cout << endl;
         cout << "      assy4r Vector ...:  " << endl;
         cout << endl;
         for(ii=0;ii<ntoteqs1;ii++)
           cout << '\t' << ii << '\t' << assy4r[ii] << endl;
         cout << endl;
         cout << endl;
*/

   //Assign data to element properties

   for(e=0;e<totnumel;e++)
      elem[e]->initialiseKnotsAtGPs();


     // compute dimensions of variables for contour plots

     if(patchGrp[0].ndom == 1)
     {
         uu.setDim(Npatch);
         vv.setDim(Npatch);
         outp.setDim(Npatch);
         Sorig.setDim(Npatch);
         Sdef.setDim(Npatch);

         for(iii=0;iii<Npatch;iii++)
         {
             VectorArray<double> uu1;
             findunique(CurveListFinal[iii].U, uu1);

             uu[iii] = uu1;

             outp[iii].setDim(uu1.n);
             Sorig[iii].setDim(uu1.n);
             Sdef[iii].setDim(uu1.n);
         }

     }
     else if(patchGrp[0].ndom == 2)
     {
         uu.setDim(patch.n);
         vv.setDim(patch.n);
         outp.setDim(patch.n);
         Sorig.setDim(patch.n);
         Sdef.setDim(patch.n);

         int cnt,  ind1, nn1=0, nn2=0;

         for(iii=0;iii<patch.n;iii++)
         {
             VectorArray<double> uu1, vv1;
             findunique(SurfaceResult[iii].U, uu1);
             findunique(SurfaceResult[iii].V, vv1);

             uu[iii] = uu1;
             vv[iii] = vv1;

             cnt = uu1.n*vv1.n;

             outp[iii].setDim(cnt);
             Sorig[iii].setDim(cnt);
             Sdef[iii].setDim(cnt);
         }

         for(iii=0;iii<patch.n;iii++)
         {
            n1 = uu[iii].n;
            n2 = vv[iii].n;
            for(jj=0;jj<n2;jj++)
            {
                ind1 = n1*jj;
                for(ii=0;ii<n1;ii++)
                {
                   Sorig[iii][ind1+ii] = SurfaceListFinal[iii].SurfacePoint(uu[iii][ii], vv[iii][jj]).CalcEuclid();
                   Sdef[iii][ind1+ii]  = SurfaceListFinal[iii].SurfacePoint(uu[iii][ii], vv[iii][jj]).CalcEuclid();
                }
            }
            nn1 += n1*n2;
         }
     }
     else
     {
         uu.setDim(patch.n);
         vv.setDim(patch.n);
         ww.setDim(patch.n);
         outp.setDim(patch.n);
         Sorig.setDim(patch.n);
         Sdef.setDim(patch.n);

         int cnt,  ind1, nn1=0, nn2=0, nn3, n3;

         for(iii=0;iii<patch.n;iii++)
         {
             VectorArray<double> uu1, vv1, ww1;
             findunique(SolidResult[iii].U, uu1);
             findunique(SolidResult[iii].V, vv1);
             findunique(SolidResult[iii].W, ww1);

             uu[iii] = uu1;
             vv[iii] = vv1;
             ww[iii] = ww1;

             cnt = uu1.n*vv1.n*ww1.n;

             outp[iii].setDim(cnt);
             Sorig[iii].setDim(cnt);
             Sdef[iii].setDim(cnt);
         }
         
         EPOINT  EP;

         for(iii=0;iii<patch.n;iii++)
         {
            n1 = uu[iii].n;
            n2 = vv[iii].n;
            n3 = ww[iii].n;
            
            nn3 = n1*n2;
            
            for(kk=0;kk<n3;kk++)
            {
                for(jj=0;jj<n2;jj++)
                {
                   ind1 = nn3*kk + n1*jj;
                   for(ii=0;ii<n1;ii++)
                   {
                      EP = SolidListFinal[iii].SolidPoint(uu[iii][ii], vv[iii][jj], ww[iii][kk]).CalcEuclid();
                      Sorig[iii][ind1+ii] = EP;
                      Sdef[iii][ind1+ii]  = EP;
                   }
               }
            }
            nn1 += n1*n2;
         }
    }

    Values.setDim(ndf);
    
    for(kk=0;kk<ndf;kk++)
    {
      Values[kk].setDim(ntotgbf);
      Values[kk].zero();
    }

    return;
}






void IsogeometricFEM::prepareMatrixPatternPetsc()
{
    char fct[] = "HBSplineFEM::prepareMatrixPatternGrid";

    int nRow, nCol, size;
    VectorXd  nnzVec(ntoteqs);

    time_t tstart, tend; 

    computerTime.go(fct);

    int  r, c, r1, c1, ii, jj, e, ee, e1, e2, pnum;
    int *tt1, *tt2, val1, val2, kk, nnz;

    printf("\n element DOF values initialised \n\n");
    printf("\n Finding Global positions \n\n");

    vector<vector<int> >  DDconn;

    vector<int>::const_iterator location;
    set<int>::iterator it;

    DDconn.resize(ntoteqs);

    tstart = time(0);

    for(ee=0;ee<totnumel;ee++)
    {
      e1   = elem[ee]->elenum;
      pnum = elem[ee]->patchnum;

      val1 = elem[ee]->nsize;
      tt1 = &(NurbsBaseFinal[pnum]->LM[e1][0]);

      val1 = elem[ee]->nsize;
      for(ii=0;ii<val1;ii++)
      {
        r = tt1[ii];
        if(r != -1)
        {
          for(jj=0;jj<val1;jj++)
          {
            c = tt1[jj];
            if(c != -1)
              DDconn[r].push_back(c);
          }
        }
      }
    }

    if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
    {
      cout << " forassembly for Kup  " << endl;

      val1 = elem[0]->nsize;
      val2 = NurbsBaseSecondVar[0]->nsize;

      cout << val1 << '\t' << val2 << endl;

      //-----------------------------
      // forassy for Kup and Kpu
      //-----------------------------

      for(ee=0;ee<totnumel;ee++)
      {
        e1   =  elem[ee]->elenum;
        e2   =  elem[ee]->elenum2;
        pnum =  elem[ee]->patchnum;
        
        //cout << ee << '\t' << e1 << '\t' << e2 << '\t' << pnum << endl;

        tt1 = &(NurbsBaseFinal[pnum]->LM[e1][0]);
        tt2 = &(NurbsBaseSecondVar[pnum]->LM[e2][0]);

        for(ii=0;ii<val1;ii++) //ndof*nlbf
        {
          r = tt1[ii];
          if(r != -1)
          {
            for(jj=0;jj<val2;jj++) //ndof*nlbf
            {
              //cout << tt1[ii] << '\t' << tt2[jj] << endl;
              c = tt2[jj];
              if(c != -1)
              {
                c += ntoteqs1;
                DDconn[r].push_back(c);
                DDconn[c].push_back(r);
              }
            }//for(jj=0;...
          }
        }//for(ii=0;...
        //cout << " ppppppppppp " << endl;

        for(ii=0;ii<val2;ii++)
        {
          r = tt2[ii];
          if(r != -1)
          {
            r += ntoteqs1;
            for(jj=0;jj<val2;jj++)
            {
              c = tt2[jj];
              if(c != -1)
                DDconn[r].push_back(c+ntoteqs1);
            }//for(jj=0;...
          }
        }//for(ii=0;...
      }//for(ee=0;...
    }//if(mixedSolverFlag == 7 || mixedSolverFlag == 12)

    printf("\n Finding Global positions DONE \n\n");

    nnz = 0;
    for(ii=0;ii<ntoteqs;ii++)
    {
      findUnique(DDconn[ii]);
      nnzVec[ii] = DDconn[ii].size();
      nnz += nnzVec[ii];
    }
    cout << " nnz " << nnz << endl;

    nRow = nCol = ntoteqs;

    /*
    MatCreateSeqAIJ(PETSC_COMM_WORLD, nRow, nRow, 1000, nnzVec, &(((SolverEigen*)solver)->mtx));

    for(ii=0;ii<ntoteqs;ii++)
    {
      //for(jj=0;jj<DDconn[ii].size();jj++)
      for(it=DDconn[ii].begin(); it!=DDconn[ii].end(); ++it)
      {
        MatSetValues(((SolverEigen*)solver)->mtx, 1, &ii, 1 , &*it, &val, ADD_VALUES);
      }
    }
    */
    cout << " AAAAAAAAAA " << endl;
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, nRow, nRow, 500, nnzVec, &(((SolverEigen*)solver)->mtx));
    solverEigen->mtx.resize(ntoteqs, ntoteqs);
    solverEigen->mtx.reserve(nnz);
    solverEigen->mtx.reserve(nnzVec);
    cout << " AAAAAAAAAA " << endl;

    for(ii=0;ii<ntoteqs;ii++)
    {
      for(jj=0;jj<DDconn[ii].size();jj++)
      {
        //cout << ii << '\t' << DDconn[ii][jj] << endl;
        //ierr = MatSetValues(((SolverEigen*)solver)->mtx, 1, &ii, 1 , &DDconn[ii][jj], &val, ADD_VALUES);
        solverEigen->mtx.coeffRef(ii, DDconn[ii][jj]) = 0.0;
      }
    }

    solverEigen->mtx.makeCompressed();

  //((SolverEigen*)solver)->printMatrix();

  cout << " BBBBBBBBBBBBB " << endl;

  solverEigen->currentStatus = PATTERN_OK;

  //printf("\n     HBSplineFEM::prepareMatrixPattern()  .... FINISHED ...\n\n");

  return;
}




void IsogeometricFEM::prepareMatrixPattern2()
{
  cout << "     ISOGEOMETRICFEM: preparing Global Coefficient Matrix Pattern ...\n\n";
  char fct[] = "IsogeometricFEM::prepareMatrixPattern";

//  computerTime.go(fct);

    int  ee, r, c, r1, c1, count=0, count1=0, count2=0, ii, jj, iii, e, val1, val2, n1, n2, kk, a, b, ll, pp;
    int *tt, *tt1, *tt2, e1, e2;

    MatrixSparse<double> matxtemp;

    ListArray<Vector<int> > tmp1;
    ListArray<VectorArray<int> > tmp;

    tmp.setDim(ntoteqs1);
    tmp1.setDim(ntoteqs1);

    for(e=0;e<totnumel;e++)
    {
       e1 = elem[e]->elenum;
       pp = elem[e]->patchnum;

       tt = &(NurbsBaseFinal[pp]->LM[e1][0]);

       val1 = elem[e]->nsize;
       for(ii=0;ii<val1;ii++)
       {
          count1 = tt[ii];
          if(count1 != -1)
            tmp1[count1].append(e);
       }
    }

    for(ii=0;ii<ntoteqs1;ii++)
        tmp[ii] = tmp1[ii];

    tmp1.free();

/*
      cout << endl;
      cout << "      dof to element connectivity ...:  " << endl;
      cout << endl;
      for(ii=0;ii<ntoteqs1;ii++)
      {
         cout << " eqn # " << ii << " : ";
	 for(jj=0;jj<tmp[ii].n;jj++)
           cout << '\t' << tmp[ii][jj];
         cout << endl;
      }
      cout << endl;
*/

cout << " AAAAAAAAAAA " << endl;

      ee = a = b = 0;

      val1 = elem[0]->nsize;

      int  pnum, elnum;

      for(ee=0;ee<totnumel;ee++) 
      {
         elnum  =  elem[ee]->elenum;
         pnum   =  elem[ee]->patchnum;

         tt = &(NurbsBaseFinal[pnum]->LM[elnum][0]);

         for(ii=0;ii<val1;ii++) //ndof*nlbf
         {
            r = tt[ii];
            r1 = r+1;
            if(r1 != 0)
            {
                for(jj=0;jj<val1;jj++) //ndof*nlbf
                {
                    c = tt[jj];
                    c1 = c+1;
                    if(c1 != 0)
                    {
                        for(ll=0;ll<tmp[c].n;ll++)
                        {
                            e1 = tmp[c][ll];

                            a = b = -1;
                            if(e1 <= ee)
                            {
                               tt1 = &(NurbsBaseFinal[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                               for(kk=0;kk<val1;kk++)
                               {
                                  if(tt1[kk] == r)    a = kk;
                                  if(tt1[kk] == c)    b = kk;
                               }
                               if(a != -1 && b != -1)    break;
                            }
                        }
                        if(a == -1 && b == -1)
                        {
                            count1 = matxtemp.append(r1,c1)+1;
                            elem[ee]->forassembly[ii][jj] = count1;
                        }
                        else
                        {
                            count1 = elem[e1]->forassembly[a][b];
                            if( count1 == -1 )
                               count1 = matxtemp.append(r1,c1)+1;

                            elem[ee]->forassembly[ii][jj] = count1;
                        }
                    }
                 }
             }
         }
      }

   cout << " forassembly for Kup  " << endl;
   if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
   {
      secVarPrev.setDim(Npatch);

      ee = a = b = 0;

      val1 = elem[0]->nsize;
      val2 = NurbsBaseSecondVar[0]->nsize;

      int *row, *col;

      //-----------------------------
      // forassy for Kup
      //-----------------------------

      for(ee=0;ee<totnumel;ee++)
      {
         elnum  =  elem[ee]->elenum;
         pnum   =  elem[ee]->patchnum;

         row = &(NurbsBaseFinal[pnum]->LM[elnum][0]);
         col = &(NurbsBaseSecondVar[pnum]->LM[elnum][0]);
         
         //cout << NurbsBaseFinal[pnum]->LM[elnum] << endl;
         //cout << NurbsBaseSecondVar[pnum]->LM[elnum] << endl;
         //cout << row << '\t' << col << endl;

         for(ii=0;ii<val1;ii++) //ndof*nlbf
         {
             r = row[ii];
            r1 = r+1;

            if(r != -1)
            {
                for(jj=0;jj<val2;jj++) //ndof*nlbf
                {
                     c = col[jj];
                    c1 = c + 1 + ntoteqs1;

                    if(c != -1)
                    {
                        for(ll=0;ll<tmp[r].n;ll++)
                        {
                            e1 = tmp[r][ll];

                            a = b = -1;
                            if(e1 <= ee)
                            {
                               tt1 = &(NurbsBaseFinal[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                               tt2 = &(NurbsBaseSecondVar[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                               for(kk=0;kk<val1;kk++)
                               {
                                  if(tt1[kk] == r)    a = kk;
                               }
                               for(kk=0;kk<val2;kk++)
                               {
                                  if(tt2[kk] == c)    b = kk;
                               }

                               if(a != -1 && b != -1)    break;
                            }
                        }
                        if(a == -1 && b == -1)
                        {
                            count1 = matxtemp.append(r1,c1)+1;
                            elem[ee]->forassyKup[ii][jj] = count1;
                        }
                        else
                        {
                            count1 = elem[e1]->forassyKup[a][b];
                            if( count1 == -1 )
                               count1 = matxtemp.append(r1,c1)+1;

                            elem[ee]->forassyKup[ii][jj] = count1;
                        }
                    }
                 }
             }
         }
      }

      //-----------------------------
      // forassy for Kpu
      //-----------------------------

      for(ee=0;ee<totnumel;ee++)
      {
         elnum  =  elem[ee]->elenum;
         pnum   =  elem[ee]->patchnum;

         row = &(NurbsBaseSecondVar[pnum]->LM[elnum][0]);
         col = &(NurbsBaseFinal[pnum]->LM[elnum][0]);

         for(ii=0;ii<val2;ii++) //ndof*nlbf
         {
            r = row[ii];
            r1 = r + 1 + ntoteqs1;
            if(r != -1)
            {
                for(jj=0;jj<val1;jj++) //ndof*nlbf
                {
                    c = col[jj];
                    c1 = c + 1;
                    if(c != -1)
                    {
                        for(ll=0;ll<tmp[c].n;ll++)
                        {
                            e1 = tmp[c][ll];

                            a = b = -1;
                            if(e1 <= ee)
                            {
                               tt2 = &(NurbsBaseFinal[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                               tt1 = &(NurbsBaseSecondVar[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                               for(kk=0;kk<val2;kk++)
                               {
                                  if(tt1[kk] == r)    a = kk;
                               }
                               for(kk=0;kk<val1;kk++)
                               {
                                  if(tt2[kk] == c)    b = kk;
                               }

                               if(a != -1 && b != -1)    break;
                            }
                        }
                        if(a == -1 && b == -1)
                        {
                            count1 = matxtemp.append(r1,c1)+1;
                            elem[ee]->forassyKpu[ii][jj] = count1;
                        }
                        else
                        {
                            count1 = elem[e1]->forassyKpu[a][b];
                            if( count1 == -1 )
                               count1 = matxtemp.append(r1,c1)+1;

                            elem[ee]->forassyKpu[ii][jj] = count1;
                        }
                    }
                 }
             }

         }
      }

         ListArray<Vector<int> > tmp3;
         ListArray<VectorArray<int> > tmp2;

         tmp2.setDim(ntoteqs2);
         tmp3.setDim(ntoteqs2);

         for(e=0;e<totnumel;e++)
         {
            e1 = elem[e]->elenum;
            pp = elem[e]->patchnum;

            tt = &(NurbsBaseSecondVar[pp]->LM[e1][0]);

            val1 = NurbsBaseSecondVar[pp]->nsize;

            for(ii=0;ii<val1;ii++)
            {
               count1 = tt[ii];
               if(count1 != -1)
                 tmp3[count1].append(e);
            }
         }

         for(ii=0;ii<ntoteqs2;ii++)
            tmp2[ii] = tmp3[ii];

         tmp3.free();

         //-----------------------------
         // forassy for Ktt is used for Kpp Matrix
         //-----------------------------

         for(ee=0;ee<totnumel;ee++)
         {
            elnum  =  elem[ee]->elenum;
            pnum   =  elem[ee]->patchnum;

            tt = &(NurbsBaseSecondVar[pnum]->LM[elnum][0]);

            for(ii=0;ii<val2;ii++) //ndof*nlbf
            {
               r = tt[ii];
               r1 = r + 1 + ntoteqs1;
               if(r != -1)
               {
                  for(jj=0;jj<val2;jj++) //ndof*nlbf
                  {
                     c = tt[jj];
                     c1 = c + 1 + ntoteqs1;
                     if(c != -1)
                     {
                         for(ll=0;ll<tmp2[r].n;ll++)
                         {
                            e1 = tmp2[r][ll];

                            a = b = -1;
                            if(e1 <= ee)
                            {
                               tt1 = &(NurbsBaseSecondVar[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                               for(kk=0;kk<val2;kk++)
                               {
                                  if(tt1[kk] == r)    a = kk;
                               }
                               for(kk=0;kk<val2;kk++)
                               {
                                  if(tt1[kk] == c)    b = kk;
                               }

                               if(a != -1 && b != -1)    break;
                            }
                        }

                        if(a == -1 && b == -1)
                        {
                            count1 = matxtemp.append(r1,c1)+1;
                            elem[ee]->forassyKtt[ii][jj] = count1;
                        }
                        else
                        {
                            count1 = elem[e1]->forassyKtt[a][b];
                            if( count1 == -1 )
                               count1 = matxtemp.append(r1,c1)+1;

                            elem[ee]->forassyKtt[ii][jj] = count1;
                        }
                    }
                 }
              }
           }
        }
   }

   if(mixedSolverFlag == 8)
   {
       secVarPrev.setDim(Npatch);
       secVarPrev2.setDim(Npatch);

       // generate gbf and elem connectivity for constraint variables
       //---------------------------------------------------------------

       bool   flag = (elem[0]->elmDat[2] == 1);

       ListArray<Vector<int> > tmp3;
       ListArray<VectorArray<int> > tmp2;

       tmp2.setDim(ntoteqs2);
       tmp3.setDim(ntoteqs2);

       for(e=0;e<totnumel;e++)
       {
          e1 = elem[e]->elenum;
          pp = elem[e]->patchnum;

          tt = &(NurbsBaseSecondVar[pp]->LM[e1][0]);

          val1 = NurbsBaseSecondVar[pp]->nsize;

          for(ii=0;ii<val1;ii++)
          {
             count1 = tt[ii];
             if(count1 != -1)
               tmp3[count1].append(e);
          }
       }

       for(ii=0;ii<ntoteqs2;ii++)
           tmp2[ii] = tmp3[ii];

       tmp3.free();
/*
      cout << endl;
      cout << "      dof to element connectivity ...:  " << endl;
      cout << endl;
      for(ii=0;ii<ntoteqs2;ii++)
      {
         cout << " eqn # " << ii << " : ";
	 for(jj=0;jj<tmp2[ii].n;jj++)
           cout << '\t' << tmp2[ii][jj];
         cout << endl;
      }
      cout << endl;
*/
      ee = a = b = 0;

      val1 = elem[0]->nsize;
      val2 = surfSecondVar[0].nsize;

      int *row, *col, dummy;

      dummy = ntoteqs1 + ntoteqs2;

      //-----------------------------
      // forassy for Kut and Ktu
      //-----------------------------

      if(flag)
      {
           for(ee=0;ee<totnumel;ee++)
           {
               elnum  =  elem[ee]->elenum;
               pnum   =  elem[ee]->patchnum;

               row = &(NurbsBaseFinal[pnum]->LM[elnum][0]);

               for(ii=0;ii<val1;ii++) //ndof*nlbf
               {
                  r = row[ii];
                  r1 = r+1;
                  if(r != -1)
                  {
                      col = &(NurbsBaseSecondVar[pnum]->LM[elnum][0]);

                      for(jj=0;jj<val2;jj++) //ndof*nlbf
                      {
                          c = col[jj];
                          c1 = c + 1 + ntoteqs1;
                          if(c != -1)
                          {
                              for(ll=0;ll<tmp[r].n;ll++)
                              {
                                  e1 = tmp[r][ll];

                                  a = b = -1;
                                  if(e1 <= ee)
                                  {
                                     tt1 = &(NurbsBaseFinal[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                                     tt2 = &(NurbsBaseSecondVar[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                                     for(kk=0;kk<val1;kk++)
                                     {
                                        if(tt1[kk] == r)    a = kk;
                                     }
                                     for(kk=0;kk<val2;kk++)
                                     {
                                        if(tt2[kk] == c)    b = kk;
                                     }

                                     if(a != -1 && b != -1)    break;
                                  }
                              }
                              if(a == -1 && b == -1)
                              {
                                  count1 = matxtemp.append(r1,c1)+1;
                                  elem[ee]->forassyKut[ii][jj] = count1;
                              }
                              else
                              {
                                  count1 = elem[e1]->forassyKut[a][b];
                                  if( count1 == -1 )
                                     count1 = matxtemp.append(r1,c1)+1;

                                  elem[ee]->forassyKut[ii][jj] = count1;
                              }
                          }
                       }
                   }
               }
            }

            for(ee=0;ee<totnumel;ee++)
            {
               elnum  =  elem[ee]->elenum;
               pnum   =  elem[ee]->patchnum;

               row = &(NurbsBaseSecondVar[pnum]->LM[elnum][0]);

               for(ii=0;ii<val2;ii++) //ndof*nlbf
               {
                  r = row[ii];
                  r1 = r + 1 + ntoteqs1;
                  if(r != -1)
                  {
                      col = &(NurbsBaseFinal[pnum]->LM[elnum][0]);

                      for(jj=0;jj<val1;jj++) //ndof*nlbf
                      {
                          c = col[jj];
                          c1 = c + 1;
                          if(c != -1)
                          {
                              for(ll=0;ll<tmp[c].n;ll++)
                              {
                                  e1 = tmp[c][ll];

                                  a = b = -1;
                                  if(e1 <= ee)
                                  {
                                     tt2 = &(NurbsBaseFinal[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                                     tt1 = &(NurbsBaseSecondVar[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                                     for(kk=0;kk<val2;kk++)
                                     {
                                        if(tt1[kk] == r)    a = kk;
                                     }
                                     for(kk=0;kk<val1;kk++)
                                     {
                                        if(tt2[kk] == c)    b = kk;
                                     }

                                     if(a != -1 && b != -1)    break;
                                  }
                              }
                              if(a == -1 && b == -1)
                              {
                                  count1 = matxtemp.append(r1,c1)+1;
                                  elem[ee]->forassyKtu[ii][jj] = count1;
                              }
                              else
                              {
                                  count1 = elem[e1]->forassyKtu[a][b];
                                  if( count1 == -1 )
                                     count1 = matxtemp.append(r1,c1)+1;

                                  elem[ee]->forassyKtu[ii][jj] = count1;
                              }
                          }
                       }
                   }
      
               }
            }
            cout << " forassy for Kut and Ktu " << endl;
      }
      //-----------------------------
      // forassy for Kup and Kpu
      //-----------------------------

      for(ee=0;ee<totnumel;ee++)
      {
         elnum  =  elem[ee]->elenum;
         pnum   =  elem[ee]->patchnum;

         row = &(NurbsBaseFinal[pnum]->LM[elnum][0]);

         for(ii=0;ii<val1;ii++) //ndof*nlbf
         {
            r = row[ii];
            r1 = r+1;
            if(r != -1)
            {
                col = &(NurbsBaseSecondVar[pnum]->LM[elnum][0]);

                for(jj=0;jj<val2;jj++) //ndof*nlbf
                {
                    c = col[jj];
                    c1 = c + 1 + dummy;
                    if(c != -1)
                    {
                        for(ll=0;ll<tmp[r].n;ll++)
                        {
                            e1 = tmp[r][ll];

                            a = b = -1;
                            if(e1 <= ee)
                            {
                               tt1 = &(NurbsBaseFinal[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                               tt2 = &(NurbsBaseSecondVar[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                               for(kk=0;kk<val1;kk++)
                               {
                                  if(tt1[kk] == r)    a = kk;
                               }
                               for(kk=0;kk<val2;kk++)
                               {
                                  if(tt2[kk] == c)    b = kk;
                               }

                               if(a != -1 && b != -1)    break;
                            }
                        }
                        if(a == -1 && b == -1)
                        {
                            count1 = matxtemp.append(r1,c1)+1;
                            elem[ee]->forassyKup[ii][jj] = count1;
                        }
                        else
                        {
                            count1 = elem[e1]->forassyKup[a][b];
                            if( count1 == -1 )
                               count1 = matxtemp.append(r1,c1)+1;

                            elem[ee]->forassyKup[ii][jj] = count1;
                        }
                    }
                 }
             }
         }
      }

      for(ee=0;ee<totnumel;ee++)
      {
         elnum  =  elem[ee]->elenum;
         pnum   =  elem[ee]->patchnum;

         row = &(NurbsBaseSecondVar[pnum]->LM[elnum][0]);

         for(ii=0;ii<val2;ii++) //ndof*nlbf
         {
            r = row[ii];
            r1 = r + 1 + dummy;
            if(r != -1)
            {
                col = &(NurbsBaseFinal[pnum]->LM[elnum][0]);

                for(jj=0;jj<val1;jj++) //ndof*nlbf
                {
                    c = col[jj];
                    c1 = c + 1;
                    if(c != -1)
                    {
                        for(ll=0;ll<tmp[c].n;ll++)
                        {
                            e1 = tmp[c][ll];

                            a = b = -1;
                            if(e1 <= ee)
                            {
                               tt2 = &(NurbsBaseFinal[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                               tt1 = &(NurbsBaseSecondVar[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                               for(kk=0;kk<val2;kk++)
                               {
                                  if(tt1[kk] == r)    a = kk;
                               }
                               for(kk=0;kk<val1;kk++)
                               {
                                  if(tt2[kk] == c)    b = kk;
                               }

                               if(a != -1 && b != -1)    break;

                            }
                        }
                        if(a == -1 && b == -1)
                        {
                            count1 = matxtemp.append(r1,c1)+1;
                            elem[ee]->forassyKpu[ii][jj] = count1;
                        }
                        else
                        {
                            count1 = elem[e1]->forassyKpu[a][b];
                            if( count1 == -1 )
                               count1 = matxtemp.append(r1,c1)+1;

                            elem[ee]->forassyKpu[ii][jj] = count1;
                        }
                    }
                 }
             }
         }
      }
      cout << " forassy for Kup and Kpu " << endl;
      //-----------------------------
      // forassy for Ktt
      //-----------------------------

      for(ee=0;ee<totnumel;ee++)
      {
         elnum  =  elem[ee]->elenum;
         pnum   =  elem[ee]->patchnum;

         tt = &(NurbsBaseSecondVar[pnum]->LM[elnum][0]);

         for(ii=0;ii<val2;ii++) //ndof*nlbf
         {
            r = tt[ii];
            r1 = r + 1 + ntoteqs1;
            if(r != -1)
            {
                for(jj=0;jj<val2;jj++) //ndof*nlbf
                {
                    c = tt[jj];
                    c1 = c + 1 + ntoteqs1;
                    if(c != -1)
                    {
                        for(ll=0;ll<tmp2[r].n;ll++)
                        {
                            e1 = tmp2[r][ll];

                            //cout << " ll and e1 " << ll << '\t' << e1 << endl;

                            a = b = -1;
                            if(e1 <= ee)
                            {
                               tt1 = &(NurbsBaseSecondVar[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                               for(kk=0;kk<val2;kk++)
                               {
                                  if(tt1[kk] == r)    a = kk;
                               }
                               for(kk=0;kk<val2;kk++)
                               {
                                  if(tt1[kk] == c)    b = kk;
                               }

                               if(a != -1 && b != -1)    break;
                            }
                        }

                        if(a == -1 && b == -1)
                        {
                            count1 = matxtemp.append(r1,c1)+1;
                            elem[ee]->forassyKtt[ii][jj] = count1;
                        }
                        else
                        {
                            count1 = elem[e1]->forassyKtt[a][b];
                            if( count1 == -1 )
                               count1 = matxtemp.append(r1,c1)+1;

                            elem[ee]->forassyKtt[ii][jj] = count1;
                        }
                    }
                 }
             }
         }
      }
      cout << " forassy for Ktt " << endl;
      //-----------------------------
      // forassy for Ktp and Kpt
      //-----------------------------

      for(ee=0;ee<totnumel;ee++)
      {
         elnum  =  elem[ee]->elenum;
         pnum   =  elem[ee]->patchnum;

         tt = &(NurbsBaseSecondVar[pnum]->LM[elnum][0]);

         for(ii=0;ii<val2;ii++) //ndof*nlbf
         {
            r = tt[ii];
            r1 = r + 1 + ntoteqs1;
            if(r != -1)
            {
                for(jj=0;jj<val2;jj++) //ndof*nlbf
                {
                    c = tt[jj];
                    c1 = c + 1 + dummy;
                    if(c != -1)
                    {
                        for(ll=0;ll<tmp2[r].n;ll++)
                        {
                            e1 = tmp2[r][ll];

                            a = b = -1;
                            if(e1 <= ee)
                            {
                               tt1 = &(NurbsBaseSecondVar[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                               for(kk=0;kk<val2;kk++)
                               {
                                  if(tt1[kk] == r)    a = kk;
                               }
                               for(kk=0;kk<val2;kk++)
                               {
                                  if(tt1[kk] == c)    b = kk;
                               }

                               if(a != -1 && b != -1)    break;
                            }
                        }
                        if(a == -1 && b == -1)
                        {
                            count1 = matxtemp.append(r1,c1)+1;
                            elem[ee]->forassyKtp[ii][jj] = count1;
                        }
                        else
                        {
                            count1 = elem[e1]->forassyKtp[a][b];
                            if( count1 == -1 )
                               count1 = matxtemp.append(r1,c1)+1;

                            elem[ee]->forassyKtp[ii][jj] = count1;
                        }
                    }
                 }
             }
         }
      }

      for(ee=0;ee<totnumel;ee++)
      {
         elnum  =  elem[ee]->elenum;
         pnum   =  elem[ee]->patchnum;

         tt = &(NurbsBaseSecondVar[pnum]->LM[elnum][0]);

         for(ii=0;ii<val2;ii++) //ndof*nlbf
         {
            r = tt[ii];
            r1 = r + 1 + dummy;
            if(r != -1)
            {
                for(jj=0;jj<val2;jj++) //ndof*nlbf
                {
                    c = tt[jj];
                    c1 = c + 1 + ntoteqs1;
                    if(c != -1)
                    {
                        for(ll=0;ll<tmp2[r].n;ll++)
                        {
                            e1 = tmp2[r][ll];

                            a = b = -1;
                            if(e1 <= ee)
                            {
                               tt1 = &(NurbsBaseSecondVar[elem[e1]->patchnum]->LM[elem[e1]->elenum][0]);

                               for(kk=0;kk<val2;kk++)
                               {
                                  if(tt1[kk] == r)    a = kk;
                               }
                               for(kk=0;kk<val2;kk++)
                               {
                                  if(tt1[kk] == c)    b = kk;
                               }

                               if(a != -1 && b != -1)    break;
                            }
                        }
                        if(a == -1 && b == -1)
                        {
                            count1 = matxtemp.append(r1,c1)+1;
                            elem[ee]->forassyKpt[ii][jj] = count1;
                        }
                        else
                        {
                            count1 = elem[e1]->forassyKpt[a][b];
                            if( count1 == -1 )
                               count1 = matxtemp.append(r1,c1)+1;

                            elem[ee]->forassyKpt[ii][jj] = count1;
                        }
                    }
                 }
             }
         }
      }
      cout << " forassy for Ktp and Kpt " << endl;
/*
      double *utmp;
      for(iii=0;iii<Npatch;iii++)
      {
         count2 = surfSecondVar[iii].Values.n ;
         utmp   = &(surfSecondVar[iii].Values[0]);

         for(ii=0;ii<count2;ii++)
           utmp[ii] = 1.0;
      }
*/
   }
   

  if(comprMtxFlg)
  {
        VectorArray<int> perm, tmp4, compr1;

        // sort matrix entries

        COUT << "sorting matrix entries for compressed row format ...\n\n";

        computerTime.go("sorting matrix entries");

        matxtemp.quickSort(true,true,perm);

        computerTime.stopAndPrint("sorting matrix entries");

        // check for double entries

       for(ii=1;ii<matxtemp.row.n;ii++)
       {
            if (matxtemp.row[ii] == matxtemp.row[ii-1])
                if (matxtemp.col[ii] == matxtemp.col[ii-1])
                   prgError(100,fct,"double entries! Could the solver cope?!");
       }

       tmp4.setDim(perm.n);

       for(ii=0;ii<perm.n;ii++)
          tmp4[perm[ii]] = ii;

             //cout << " tmp4 \n " << tmp4 << endl;
             //cout << " perm \n " << perm << endl;

       int  size1, size2;
      for(ee=0;ee<totnumel;ee++)
      {
          for(ee=0;ee<totnumel;ee++)
          {
              size1  = elem[ee]->nsize;
              for(ii=0;ii<size1;ii++)
              {
                  for(jj=0;jj<size1;jj++)
                  {
                      if(elem[ee]->forassembly[ii][jj] != -1)
                         elem[ee]->forassembly[ii][jj] = tmp4.x[elem[ee]->forassembly[ii][jj]-1] +1;
                  }
              }
          }
      }

      if(mixedSolverFlag == 7 || mixedSolverFlag == 12)
      {
          bool  flg = (NurbsBaseFinal[0]->ElemProp.data[9] == 1);
          for(ee=0;ee<totnumel;ee++)
          {
              size1  = elem[ee]->nsize;
              if(mixedSolverFlag == 7)
                size2  = elem[ee]->surf2->nsize;

              if(mixedSolverFlag == 12)
                size2  = elem[ee]->solid2->nsize;


              for(ii=0;ii<size1;ii++)
              {
                  for(jj=0;jj<size2;jj++)
                  {
                      if(elem[ee]->forassyKup[ii][jj]  != -1)
                         elem[ee]->forassyKup[ii][jj] = tmp4.x[elem[ee]->forassyKup[ii][jj]-1] +1;

                       if(elem[ee]->forassyKpu[jj][ii]  != -1)
                         elem[ee]->forassyKpu[jj][ii] = tmp4.x[elem[ee]->forassyKpu[jj][ii]-1] +1;
                  }
              }
              if(flg)
              {
                 for(ii=0;ii<size2;ii++)
                 {
                     for(jj=0;jj<size2;jj++)
                     {
                         if(elem[ee]->forassyKtt[ii][jj]  != -1)
                            elem[ee]->forassyKtt[ii][jj] = tmp4.x[elem[ee]->forassyKtt[ii][jj]-1] +1;
                     }
                 }
              }
          }
      }
      if(mixedSolverFlag == 8)
      {
            bool   flag = (elem[0]->elmDat[2] == 1);
            for(ee=0;ee<totnumel;ee++)
            {
                  size1  = elem[ee]->nsize;
                  size2  = elem[ee]->surf2->nsize;

                  for(ii=0;ii<size1;ii++)
                 {
                       for(jj=0;jj<size2;jj++)
                      {
                            if(elem[ee]->forassyKup[ii][jj]  != -1)
                                 elem[ee]->forassyKup[ii][jj] = tmp4.x[elem[ee]->forassyKup[ii][jj]-1] +1;

                            if(elem[ee]->forassyKpu[jj][ii]  != -1)
                                 elem[ee]->forassyKpu[jj][ii] = tmp4.x[elem[ee]->forassyKpu[jj][ii]-1] +1;
                      }
                 }

                  for(ii=0;ii<size2;ii++)
                 {
                       for(jj=0;jj<size2;jj++)
                      {
                            if(elem[ee]->forassyKtt[ii][jj]  != -1)
                                 elem[ee]->forassyKtt[ii][jj] = tmp4.x[elem[ee]->forassyKtt[ii][jj]-1] +1;

                            if(elem[ee]->forassyKtp[ii][jj]  != -1)
                                 elem[ee]->forassyKtp[ii][jj] = tmp4.x[elem[ee]->forassyKtp[ii][jj]-1] +1;

                            if(elem[ee]->forassyKpt[ii][jj]  != -1)
                                 elem[ee]->forassyKpt[ii][jj] = tmp4.x[elem[ee]->forassyKpt[ii][jj]-1] +1;

                      }
                 }
                 if(flag)
                 {
                    for(ii=0;ii<size1;ii++)
                    {
                        for(jj=0;jj<size2;jj++)
                        {
                            if(elem[ee]->forassyKut[ii][jj]  != -1)
                                 elem[ee]->forassyKut[ii][jj] = tmp4.x[elem[ee]->forassyKut[ii][jj]-1] +1;

                            if(elem[ee]->forassyKtu[jj][ii]  != -1)
                                 elem[ee]->forassyKtu[jj][ii] = tmp4.x[elem[ee]->forassyKtu[jj][ii]-1] +1;
                        }
                    }
     		 }
            }
      }

     perm.free();
     tmp4.free();

     // generate compr
     compr1.setDim(matxtemp.nRow+1);
     jj = 0;
     compr1[jj++] = 1;
     for(ii=1;ii<matxtemp.row.n;ii++)
     {
          if(matxtemp.row[ii] != matxtemp.row[ii-1])
             compr1[jj++] = ii+1;
     }
     
     compr1[jj] = matxtemp.col.n+1;

//     cout << compr1 << endl;

     if(jj != matxtemp.nRow)
        prgError(200,fct,"fatal error!");

      //solverEigen->compr = compr1;

      compr1.free();

      //cout << matxtemp.col << "\n\n" << matxtemp.row << "\n\n" << compr << "\n\n";

      COUT << "matrix compression done.\n\n";
  }

/*
       cout << endl;
       cout << "      forassembly ...:  " << endl;
       cout << endl;
	for(e=0;e<totnumel;e++)
	{
	   val1 = elem[e]->nsize;

	   for(ii=0;ii<val1;ii++)
	   {
             for(jj=0;jj<val1;jj++)
                cout << '\t' << elem[e]->forassembly[ii][jj];
             cout << endl;
          }
          cout << endl;
          cout << endl;
	}


       cout << endl;
       cout << "      forassembly C ...:  " << endl;
       cout << endl;
	for(e=0;e<totnumel;e++)

	{
	   val1 = elem[e]->nsize;
	   val2 = surfSecondVar[0].nsize;

	   for(ii=0;ii<val1;ii++)
	   {
             for(jj=0;jj<val2;jj++)
                cout << '\t' << elem[e]->forassyKup[ii][jj];
             cout << endl;
          }

          cout << endl;
          cout << endl;
	}


       cout << endl;
       cout << "      forassembly Ct ...:  " << endl;
       cout << endl;
	for(e=0;e<totnumel;e++)
	{
	   val1 = elem[e]->nsize;
	   val2 = surfSecondVar[0].nsize;

	   for(ii=0;ii<val2;ii++)
	   {
             for(jj=0;jj<val1;jj++)
                cout << '\t' << elem[e]->forassyKpu[ii][jj];
             cout << endl;
          }
          cout << endl;
          cout << endl;
	}

*/



   if (ntoteqs1 < 30) matxtemp.printPattern();

   //MatrixSparseArray<double> &mtx2 = solverEigen->mtx;
   //mtx2 = matxtemp;

   // if(massMatrixflag)     M_global = matxtemp;

   matxtemp.free();

  // finalise
  solverEigen->currentStatus = PATTERN_OK;

//  computerTime.stopAndPrint(fct);

  return;
}


