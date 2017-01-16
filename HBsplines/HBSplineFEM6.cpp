

#include "HBSplineFEM.h"

#include "SolverPardisoEigen.h"
#include "SolverMA41Eigen.h"
#include "SolverEigen.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "ImmersedIntegrationElement.h"


extern ComputerTime       computerTime;
extern MpapTime           mpapTime;




int  HBSplineFEM::prepareMatrixPattern()
{
    char fct[] = "HBSplineFEM::prepareMatrixPattern";
    
    //cout << " HBSplineFEM::prepareMatrixPattern() " << endl;
    //cout << " GRID_CHANGED " << GRID_CHANGED << endl;
    //cout << " IB_MOVED " << IB_MOVED << endl;
    //cout << " STAGGERED " << STAGGERED << endl;
    
    //cout << " velDOF   = " << velDOF << endl;
    //cout << " fluidDOF = " << fluidDOF << endl;
    //cout << " IBDOF    = " << IBDOF << endl;
    //cout << " solidDOF = " << solidDOF << endl;
    //cout << " totalDOF = " << totalDOF << endl;

    int ii, jj, count, nnz, r, c;
    VectorXi  nnzVec(totalDOF);

    time_t tstart, tend;

    tstart = time(0);

    //printf("\n element DOF values initialised \n\n");
    //printf("\n Finding Global positions \n\n");

    ///////////////////////////////////////////////////////
    //
    // prepare the dof-to-dof connectivity
    //
    ///////////////////////////////////////////////////////

    if(GRID_CHANGED)
    {
      prepareMatrixPatternFluid();
      prepareMatrixPatternPostProcess();
      GRID_CHANGED = false;

      for(ii=0;ii<fluidDOF;ii++)
        findUnique(DDconn[ii]);
    }
    //cout << " AAAAAAAAAAAAA " << endl;

    if( (IBDOF > 0) && IB_MOVED )
    {
      prepareMatrixPatternLagrangeMultipliers();

      for(ii=0;ii<velDOF;ii++)
      {
        if( DDconnEt[ii].size() > 0)
          findUnique(DDconnEt[ii]);
      }
      for(ii=0;ii<IBDOF;ii++)
        findUnique(DDconnE[ii]);
    }
    //cout << " CCCCCCCCCCCCCCCCC " << endl;
    //cout << " DDconnE.size()   = " ;
    //cout <<  DDconnE.size() << endl;
    //cout << " DDconnEt.size()  = " ;
    //cout <<  DDconnEt.size() << endl;
    //cout << " CCCCCCCCCCCCCCCCC " << endl;

    if(!STAGGERED)
    {
      prepareMatrixPatternSolid();
      for(ii=0;ii<IBDOF;ii++)
        findUnique(DDconnHt[ii]);

      for(ii=0;ii<solidDOF;ii++)
      {
        findUnique(DDconnG[ii]);
        findUnique(DDconnH[ii]);
      }
    }

    //printf("\n Finding Global positions DONE \n\n");

    //for(ii=0;ii<fluidDOF;ii++)
      //if( DDconnEt[ii].size() > 0)
        //cout << ii << '\t' << DDconnEt[ii].size() << '\t' << DDconn[ii].size() << endl;

    ///////////////////////////////////////////////////////
    //
    // compute the number of nonzeros in each row
    //
    ///////////////////////////////////////////////////////

    nnzVec.setZero();
    for(ii=0;ii<fluidDOF;ii++)
      nnzVec[ii] = DDconn[ii].size();

    if(IBDOF > 0)
    {
      for(ii=0;ii<velDOF;ii++)
        nnzVec[ii] += DDconnEt[ii].size();

      for(ii=0;ii<IBDOF;ii++)
        nnzVec[fluidDOF+ii] += DDconnE[ii].size();
    }

    if(!STAGGERED)
    {
      for(ii=0;ii<IBDOF;ii++)
        nnzVec[fluidDOF+ii] += DDconnHt[ii].size();

      count = fluidDOF + IBDOF;
      for(ii=0;ii<solidDOF;ii++)
        nnzVec[count+ii] += DDconnH[ii].size() + DDconnG[ii].size();
    }

    nnz = 0;
    for(ii=0;ii<totalDOF;ii++)
      nnz += nnzVec[ii];

    //cout << " nnz " << nnz << endl;

    tend = time(0);

    cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;

    ///////////////////////////////////////////////////////
    //
    // allocate the memory for the sparse matrix
    //
    ///////////////////////////////////////////////////////

    //cout << " AAAAAAAAAA " << endl;
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, nRow, nRow, 500, nnzVec, &(((SolverEigen*)solver)->mtx));
    //cout << ii << '\t' << DDconn[ii][jj] << endl;
    //ierr = MatSetValues(((SolverEigen*)solver)->mtx, 1, &ii, 1 , &DDconn[ii][jj], &val, ADD_VALUES);

    solverEigen->mtx.setZero();
    solverEigen->mtx.uncompress();

    //solverEigen->mtx.resize(totalDOF, totalDOF);
    solverEigen->mtx.conservativeResize(totalDOF, totalDOF);
    //solverEigen->mtx.reserve(nnz);
    solverEigen->mtx.reserve(nnzVec);
    //cout << " AAAAAAAAAA " << endl;

    ///////////////////////////////////////////////////////
    //
    // set the entries in the sparse matrix
    //
    ///////////////////////////////////////////////////////

    typedef Eigen::Triplet<double> T;

    vector<T> tripletList;

    tripletList.reserve(nnz);

    //cout << " DDconn.size()    = " << DDconn.size() << endl;

    for(ii=0;ii<DDconn.size();ii++)
    {
      for(jj=0;jj<DDconn[ii].size();jj++)
        //solver->mtx.coeffRef(ii, DDconn[ii][jj]) = 0.0;
        tripletList.push_back(T(ii, DDconn[ii][jj], 0.0));
    }
    //cout << " BBBBBBBBBBBBBBB " << endl;

    //cout << " DDconnE.size()   = " ;
    //cout <<  DDconnE.size() << endl;
    //cout << " DDconnEt.size()  = " ;
    //cout <<  DDconnEt.size() << endl;
    //cout << " BBBBBBBBBBBBBBB " << endl;

    for(ii=0;ii<DDconnE.size();ii++)
    {
      r = fluidDOF + ii;
      for(jj=0;jj<DDconnE[ii].size();jj++)
      {
        c = DDconnE[ii][jj];
        //solver->mtx.coeffRef(r, c) = 0.0;
        //solver->mtx.coeffRef(c, r) = 0.0;

        tripletList.push_back(T(r, c, 0.0));
        tripletList.push_back(T(c, r, 0.0));
      }
    }

    //cout << " DDconnH.size()   = " << DDconnH.size() << endl;
    //cout << " DDconnHt.size()  = " << DDconnHt.size() << endl;

    count = fluidDOF + IBDOF;
    for(ii=0;ii<DDconnH.size();ii++)
    {
      r = count + ii;
      for(jj=0;jj<DDconnH[ii].size();jj++)
      {
        c = fluidDOF + DDconnH[ii][jj];
        //solver->mtx.coeffRef(r, c) = 0.0;
        //solver->mtx.coeffRef(c, r) = 0.0;

        tripletList.push_back(T(r, c, 0.0));
        tripletList.push_back(T(c, r, 0.0));
      }
    }

    //cout << " DDconnG.size()   = " << DDconnG.size() << endl;

    for(ii=0;ii<DDconnG.size();ii++)
    {
      r = count + ii;
      for(jj=0;jj<DDconnG[ii].size();jj++)
      {
        c = count + DDconnG[ii][jj];
        //solver->mtx.coeffRef(r, c) = 0.0;

        tripletList.push_back(T(r, c, 0.0));
      }
    }

    solverEigen->mtx.setFromTriplets(tripletList.begin(), tripletList.end());

    solverEigen->mtx.makeCompressed();

    solverEigen->currentStatus = PATTERN_OK;

    GRID_CHANGED = IB_MOVED = false;

    printf("\n     HBSplineFEM::prepareMatrixPattern()  .... FINISHED ...\n\n");

    return 1;
}
//



void HBSplineFEM::prepareMatrixPatternFluid()
{
    int  r, c, ii, jj, e, ee, val1, val2;
    int  *tt1, *tt2;

    DDconn.clear();
    DDconn.resize(fluidDOF);

    //cout << " activeElements.size() " << activeElements.size() << '\t' << totalDOF << endl;
    for(e=0;e<activeElements.size();e++)
    {
        ee = activeElements[e];

        val1 =  elems[ee]->forAssyVec.size();
        tt1  =  &(elems[ee]->forAssyVec[0]);
        //cout << e << '\t' << ee << '\t' << val1 << endl;
        //printVector(elems[ee]->forAssyVec);

        for(ii=0;ii<val1;ii++)
        {
          r = tt1[ii];

          for(jj=0;jj<val1;jj++)
            DDconn[r].push_back(tt1[jj]);
        }
       
        if( e % 1000 == 0)
        {
          for(ii=0;ii<fluidDOF;ii++)
            findUnique(DDconn[ii]);
        }
    }

    //cout << " solver->STABILISED  " << solver->STABILISED << endl;

    /*
    if( !LSFEM_FLAG && (ndof > 1) )
    {
      for(e=0;e<activeElements.size();e++)
      {
        ee = activeElements[e];

        val1 =  elems[ee]->forAssyVec.size();
        val2 =  elems[ee]->forAssyVec2.size();
        tt1  =  &(elems[ee]->forAssyVec[0]);
        tt2  =  &(elems[ee]->forAssyVec2[0]);
        //printVector(elems[ee]->forAssyVec);

        for(ii=0;ii<val2;ii++)
        {
          r = velDOF + tt2[ii];

          for(jj=0;jj<val1;jj++)
          {
            c = tt1[jj];
            DDconn[r].push_back(c);
            DDconn[c].push_back(r);
          }
        }

        if(solver->STABILISED)
        {
          for(ii=0;ii<val2;ii++)
          {
            r = velDOF + tt2[ii];

            for(jj=0;jj<val2;jj++)
             DDconn[r].push_back(velDOF + tt2[jj]);
          }
        }
      }
    }
    */
    
    return;
}

void  HBSplineFEM::prepareMatrixPatternPostProcess()
{
    ///////////////////////////// 
    // matrix pattern for the matrix to compute vorticity
    ///////////////////////////// 

    int  *tt1, r, c, ii, jj, e, ee, val1, val2, dd, nnz;

    dd = velDOF/ndof;

    VectorXi vectemp(dd);
    vector<vector<int> > conn2;

    SolnData.vorticity.resize(dd);
    conn2.resize(dd);

    for(e=0;e<activeElements.size();e++)
    {
        ee = activeElements[e];

        val1 =  elems[ee]->GlobalBasisFuncs.size();
        tt1  =  &(elems[ee]->GlobalBasisFuncs[0]);

        for(ii=0;ii<val1;ii++)
        {
          r = tt1[ii];

          for(jj=0;jj<val1;jj++)
            conn2[r].push_back(tt1[jj]);
        }
    }

    nnz = 0;
    for(ii=0;ii<dd;ii++)
    {
      findUnique(conn2[ii]);
      vectemp[ii] = conn2[ii].size();
      nnz += vectemp[ii];
    }

    //cout << " nnz " << nnz << endl;

    globalK2.setZero();

    globalK2.resize(dd, dd);
    globalK2.reserve(nnz);
    globalK2.reserve(vectemp);

    ///////////////////////////////////////////////////////
    //
    // set the entries in the sparse matrix
    //
    ///////////////////////////////////////////////////////

    for(ii=0;ii<conn2.size();ii++)
    {
      for(jj=0;jj<conn2[ii].size();jj++)
        globalK2.coeffRef(ii, conn2[ii][jj]) = 0.0;
    }

    return;
}


void HBSplineFEM::prepareMatrixPatternLagrangeMultipliers()
{
    //cout << " HBSplineFEM::prepareMatrixPatternLagrangeMultipliers() " << endl;

    int  r, c, ii, jj, ee, val1, val2, aa, bb, count, gp;
    int  *tt1, *tt2;

    //DDconnE.clear();
    //DDconnEt.clear();
    DDconnE.resize(IBDOF);
    DDconnEt.resize(velDOF);

    for(ii=0; ii<IBDOF; ii++)
      DDconnE[ii].clear();

    for(ii=0; ii<velDOF; ii++)
      DDconnEt[ii].clear();


    ImmersedIntegrationElement *lme;

    for(bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
      for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
      {
        lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

        //cout << bb << '\t' << aa << '\t' << ImmersedBodyObjects[bb]->ImmIntgElems.size() << '\t' << lme->gausspoints.size() << endl;

        for(gp=0;gp<lme->gausspoints.size();gp++)
        {
          lme->computePointAtGP(gp, geom);
          r = findCellNumber(geom);

          lme->elemNums[gp] = r;
          //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %5d, \n", geom[0], geom[1], geom[2], r);
        }
      }
    }
    //cout << " HBSplineFEM::prepareMatrixPatternLagrangeMultipliers() " << endl;

    for(bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
      if(ImmersedBodyObjects[bb]->IsBoundaryConditionTypeLagrange())
      {
        for(aa=0;aa<ImmersedBodyObjects[bb]->ImmIntgElems.size();aa++)
        {
          lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];

          val1 =  lme->posIndices.size();
          tt1  =  &(lme->posIndices[0]);
          
          for(ee=0;ee<lme->elemNums.size();ee++)
          {
            val2 =  elems[lme->elemNums[ee]]->forAssyVec.size();
            tt2  =  &(elems[lme->elemNums[ee]]->forAssyVec[0]);

            //cout << bb << '\t' << ee << '\t' << val1 << endl;
            //for(kk=0;kk<val1;kk++)
              //cout << tt[kk] << '\t' ;
            //cout << endl;

            for(ii=0;ii<val1;ii++)
            {
              r = tt1[ii];

              for(jj=0;jj<val2;jj++)
              {
                c = tt2[jj];
                DDconnE[r].push_back(c);
                DDconnEt[c].push_back(r);
              }
              // LagrangeMultiplier stabilisation terms
              //for(jj=0;jj<val1;jj++)
                //DDconnE[r].push_back(fluidDOF+tt1[jj]);
            }
          }
        } // for(aa=0...
      } // if(
    } // for(bb=0...

  return;
}



void HBSplineFEM::prepareMatrixPatternSolid()
{
    int  r, c, ii, jj, bb;

    //DDconnG.clear();
    //DDconnH.clear();
    //DDconnHt.clear();

    DDconnG.resize(solidDOF);
    DDconnH.resize(solidDOF);
    DDconnHt.resize(IBDOF);

    for(ii=0; ii<solidDOF; ii++)
    {
      DDconnG[ii].clear();
      DDconnH[ii].clear();
    }

    for(ii=0; ii<IBDOF; ii++)
      DDconnHt[ii].clear();


    for(bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
        for(ii=0;ii<ImmersedBodyObjects[bb]->forAssyMat.size();ii++)
        {
          for(jj=0;jj<ImmersedBodyObjects[bb]->forAssyMat[ii].size();jj++)
          {
            DDconnG[ii].push_back(ImmersedBodyObjects[bb]->forAssyMat[ii][jj]);
          }
        }

        // connectivity for coupled matrices

        for(ii=0;ii<ImmersedBodyObjects[bb]->forAssyCoupledHorz.size();ii++)
        {
          for(jj=0;jj<ImmersedBodyObjects[bb]->forAssyCoupledHorz[ii].size();jj++)
          {
            c = ImmersedBodyObjects[bb]->forAssyCoupledHorz[ii][jj];
            DDconnH[ii].push_back(c);
            DDconnHt[c].push_back(ii);
          }
        }
    }//for(bb=0;...

  for(bb=0;bb<contactElementObjects.size();bb++)
  {
    for(ii=0;ii<2;ii++)
    {
      for(jj=0;jj<2;jj++)
      {
        DDconnG[ii].push_back(jj);
      }
    }
  }

  return;
}


