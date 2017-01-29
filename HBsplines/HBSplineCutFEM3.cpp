

#include "HBSplineCutFEM.h"
#include "SolverPardisoEigen.h"
#include "SolverMA41Eigen.h"
#include "SolverEigen.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "conditionalOStream.h"
#include "my_types.h"
#include "mpi.h"
#include "metis.h"



extern ComputerTime       computerTime;
extern MpapTime           mpapTime;


void HBSplineCutFEM::setTimeParam()
{
  //cout << " HBSplineCutFEM::setTimeParam() " << endl;

  SolnData.setTimeParam();

  for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
    ImmersedBodyObjects[bb]->setTimeParam();

  //cout << " HBSplineCutFEM::setTimeParam() " << endl;

  return;
}



void HBSplineCutFEM::timeUpdate()
{
  //cout << " HBSplineCutFEM::timeUpdate() ... STARTED " << endl;

  IB_MOVED = false;

  firstIter = true;
  localStiffnessError = 0;
  iterCount = 1;
  filecount++;

  for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
  {
    ImmersedBodyObjects[bb]->timeUpdate();
  }

  //cout << " AAAAAAAAAAAAAAA " << endl;

  SolnData.timeUpdate();

  //cout << " AAAAAAAAAAAAAAA " << endl;

  if(STAGGERED)
  {
    solveSolidProblem();
    //updateImmersedPointPositions();
  }

  //cout << " AAAAAAAAAAAAAAA " << endl;

  updateIterStep();

  //cout << " HBSplineCutFEM::timeUpdate() ... FINISHED " << endl;

  return;
}



void HBSplineCutFEM::updateIterStep()
{
  //cout << " HBSplineCutFEM::updateIterStep() ... STARTED " << endl;

  int kk, bb, ee;

  SolnData.updateIterStep();

  //cout << " kkkkkkkkkkk " << endl;

  for(bb=0;bb<ImmersedBodyObjects.size();bb++)
  {
    computeTotalForce(bb);
    ImmersedBodyObjects[bb]->updateForce(&(totalForce(0)));
  }

  if(!STAGGERED)
  {
    for(bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
      //ImmersedBodyObjects[bb]->updateDisplacement(&(SolnData.var1(fluidDOF)));
      //cout << " AAAAAAAAAAA " << bb << endl;
      ImmersedBodyObjects[bb]->updateIterStep();
      //cout << " AAAAAAAAAAA " << bb << endl;
    }
    //cout << " AAAAAAAAAAA " << bb << endl;
    updateImmersedPointPositions();
    IB_MOVED = true;
  }

  //cout << " kkkkkkkkkkk " << endl;
  //cout << " PPPPPPPPPPPPPPPPP " << firstIter << '\t' << SOLVER_TYPE << '\t' << GRID_CHANGED << '\t' << IB_MOVED << endl;
  if(GRID_CHANGED || IB_MOVED)
  {
      solverEigen->free();

      prepareMatrixPattern();

      switch(slv_type)
      {
          case  1: // MA41 ..........................
          case  4: // SolverEigen ..........................

            //printInfo();

            if(solverEigen->initialise(0,0,totalDOF) != 0)
              return;

            //solverEigen->printInfo();

          break;

          case  5: // PARDISO(sym) with Eigen
          case  6: // PARDISO(unsym) with Eigen

            if(slv_type == 5)
            {
              if(solverEigen->initialise(numProc, PARDISO_STRUCT_SYM, totalDOF) != 0)
                return;
            }
            if(slv_type == 6)
            {
              if(solverEigen->initialise(numProc, PARDISO_UNSYM, totalDOF) != 0)
                return;
            }

          break;

        case  8: // SolverPetsc ..........................

            if(solverPetsc->initialise(0, 0, totalDOF) != 0)
              return;

            solverPetsc->SetSolverAndParameters();
            //cout << " kkkkkkkkkk " << endl;
            solverPetsc->printInfo();
            //cout << " aaaaaaaaa " << endl;

        break;

          default: // invalid SOLVER_TYPE ...................

            cout << " this solver has not been implemented yet!\n\n";

          break;
      }

      solverOK = true;
 
      if(solverEigen != NULL)
        solverEigen->checkIO = true;

      if(solverPetsc != NULL)
        solverPetsc->checkIO = true;
  }

  //cout << " HBSplineCutFEM::updateIterStep() ... FINISHED " << endl;

  return;
}




void HBSplineCutFEM::reset()
{
  SolnData.reset();

  for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
    ImmersedBodyObjects[bb]->reset();

  return;
}




int  HBSplineCutFEM::setCoveringUncovering()
{
  for(int dd=0; dd<forAssyCutFEM.size(); dd++)
  {
    if( (forAssyCutFEM[dd] != -1) && (forAssyCutFEMprev[dd] == -1) )
    {
      //if( dd % 3 == 1 )
        SolnData.var1Prev[dd] = 0.0;
    }
  }

  return  1;
}






int  HBSplineCutFEM::prepareMatrixPattern()
{
    assert(SOLVER_TYPE == SOLVER_TYPE_PETSC);

    time_t tstart, tend; 

    int  tempDOF, domTemp, npElem, ind;
    int  r, c, r1, c1, count=0, count1=0, count2=0, ii, jj, e, ee, dd, ind1, ind2;
    int  *tt1, *tt2, val1, val2, nnz, n1, n2, kk, e1, e2, a, b, ll, pp, aa, bb;
    int  side, NUM_NEIGHBOURS=2*DIM, start1, start2, nr1, nr2;

    node  *nd1, *nd2;

    std::vector<int>::iterator itInt;

    MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);

    //cout << " this_mpi_proc " << this_mpi_proc << endl;
    //cout << " n_mpi_procs " << n_mpi_procs << endl;

    int n_subdomains = n_mpi_procs, subdomain=0;


    tstart = time(0);

    prepareCutElements();

    tend = time(0);
    //printf("HBSplineCutFEM::prepareCutElements() took %8.4f second(s) \n ", difftime(tend, tstart) );

    nElem = activeElements.size();
    nNode = gridBF1;

    PetscSynchronizedPrintf(MPI_COMM_WORLD, " nNode = %8d \n", nNode);

    if(n_mpi_procs == 1)
    {
      elem_start = 0;
      elem_end   = nElem-1;

      nElem_local = nElem;

      totalDOF = nNode*ndof;

      row_start = 0;
      row_end   = totalDOF-1;

      ndofs_local = totalDOF;

      node_map_new_to_old.resize(nNode, 0);
      node_map_old_to_new.resize(nNode, 0);

      dof_map_new_to_old.resize(totalDOF, 0);
      dof_map_old_to_new.resize(totalDOF, 0);

      kk=0;
      for(ii=0; ii<nNode; ii++)
      {
        node_map_new_to_old[ii] = ii;
        node_map_old_to_new[ii] = ii;
        
        for(jj=0; jj<ndof; jj++)
        {
          dof_map_new_to_old[kk] = kk;
          dof_map_old_to_new[kk] = kk;
          kk++;
        }
      }
  }
  else
  {
      /////////////////////////////////////////////////////////////////////////////
      //
      // Partition the mesh. This can be done using software libraries
      // Chaco, Jostle, METIS and Scotch.
      // Here METIS is used.
      // 
      /////////////////////////////////////////////////////////////////////////////

      n2 = ceil(nElem/n_mpi_procs);

      elem_start = n2*this_mpi_proc;
      elem_end   = n2*(this_mpi_proc+1)-1;

      if(this_mpi_proc == (n_mpi_procs-1))
        elem_end = nElem-1;

      nElem_local = elem_end - elem_start + 1;

      cout << " elem_start = " << elem_start << '\t' << elem_end << '\t' << nElem_local << endl;

      idx_t eptr[nElem+1];
      idx_t epart[nElem];
      idx_t npart[nNode];


      eptr[0] = 0;

      idx_t npElem_total = 0;
      for(e1=0; e1<nElem; e1++)
      {
        npElem_total += elems[activeElements[e1]]->GlobalBasisFuncs.size();

        eptr[e1+1] = npElem_total;
      }

      cout << " npElem_total = " << npElem_total << endl;

      idx_t eind[npElem_total];

      //PetscInt  *eptr, *eind;

      //PetscInt  eptr[nElem_local+1];
      //PetscInt  eind[nElem_local*npElem];

      //ierr  = PetscMalloc1(nElem_local+1,  &eptr);CHKERRQ(ierr);
      //ierr  = PetscMalloc1(nElem_local*npElem,  &eind);CHKERRQ(ierr);

      //ierr  = PetscMalloc1(nElem+1,  &eptr);CHKERRQ(ierr);
      //ierr  = PetscMalloc1(nElem*npElem,  &eind);CHKERRQ(ierr);

      //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d \n", elem_start);

      vector<int>  vecTemp2;

      kk=0;
      //for(e1=0; e1<nElem_local; e1++)
      for(e1=0; e1<nElem; e1++)
      {
        //e2 = elem_start + e1;

        vecTemp2 = elems[activeElements[e1]]->GlobalBasisFuncs ;

        //if(this_mpi_proc == 0)
          //printVector(vecTemp2);
        //findUnique(vecTemp2);

        npElem = vecTemp2.size();

        for(ii=0; ii<npElem; ii++)
          eind[kk+ii] = vecTemp2[ii] ;

        kk += npElem;
      }

      //for(e1=0; e1<nElem; e1++)
      //{
        //kk = e1*npElem;
        //for(ii=0; ii<npElem; ii++)
          //cout << eind[kk+ii] << '\t' ;
        //cout << endl;
      //}

      //PetscInt  nodes_per_side;
      idx_t nodes_per_side;
      
      if(DIM == 2)
        nodes_per_side = 2;
      else
        nodes_per_side = 3;


      //
      idx_t  nWeights  = 1;
      idx_t  nParts = n_mpi_procs;
      idx_t  objval;
      idx_t  *xadj, *adjncy;
      idx_t  numflag=0;

      idx_t options[METIS_NOPTIONS];

      METIS_SetDefaultOptions(options);

      // Specifies the partitioning method.
      options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;    // Multilevel recursive bisectioning.
      options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;  // Multilevel k-way partitioning.

      //options[METIS_OPTION_NSEPS] = 10;

      options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;   // Edge-cut minimization
      //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL; // Total communication volume minimization

      options[METIS_OPTION_NUMBERING] = 0;  // C-style numbering is assumed that starts from 0.
      //options[METIS_OPTION_NUMBERING] = 1;  // Fortran-style numbering is assumed that starts from 1.


      // METIS partition routine
      int ret = METIS_PartMeshNodal(&nElem, &nNode, eptr, eind, NULL, NULL, &nParts, NULL, options, &objval, epart, npart);
      //int ret = METIS_PartMeshDual(&nElem, &nNode, eptr, eind, NULL, NULL, &nodes_per_side, &n_mpi_procs, NULL, options, &objval, epart, npart);

      //idx_t  wgtflag=0, numflag=0, ncon=0, ncommonnodes=2, nparts=n_mpi_procs;
      //idx_t  *elmdist;

      //int ret = ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, NULL, &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts, NULL, NULL, options, NULL, npart, MPI_COMM_WORLD);

      //pcout <<  " ccccccccccc " << endl;

      if(ret == METIS_OK)
        std::cout << " METIS partition routine successful "  << std::endl;
      else
        std::cout << " METIS partition routine FAILED "  << std::endl;

      //
      //if(this_mpi_proc == 0)
      //{
        //for(ee=0; ee<nElem; ee++)
          //cout << ee << '\t' << epart[ee] << endl;

        //cout << " \n\n\n\n " << endl;

        //for(ee=0; ee<nNode; ee++)
          //cout << ee << '\t' << npart[ee] << endl;
      //}
      //

      //cout << " nNode_local = " << nNode_local << endl;

      for(e1=0; e1<nElem; e1++)
        elems[activeElements[e1]]->set_subdomain_id(epart[e1]);


      std::vector<std::vector<int> >  locally_owned_nodes_total;
      std::vector<int>  locally_owned_nodes;

      locally_owned_nodes_total.resize(n_subdomains);

      for(ee=0; ee<nNode; ee++)
      {
        locally_owned_nodes_total[npart[ee]].push_back(ee);
      }

      locally_owned_nodes = locally_owned_nodes_total[this_mpi_proc];

      nNode_local = locally_owned_nodes.size();

      cout << " nNode_local = " << nNode_local << '\t' << this_mpi_proc << endl;

      node_map_new_to_old.resize(nNode, 0);
      node_map_old_to_new.resize(nNode, 0);

      ii=0;
      for(subdomain=0; subdomain<n_subdomains; subdomain++)
      {
        for(kk=0; kk<locally_owned_nodes_total[subdomain].size(); kk++)
        {
          node_map_new_to_old[ii++] = locally_owned_nodes_total[subdomain][kk];
        }
      }

      for(ii=0; ii<nNode; ii++)
        node_map_old_to_new[node_map_new_to_old[ii]] = ii;

      totalDOF = nNode*ndof;

      dof_map_new_to_old.resize(totalDOF, 0);
      dof_map_old_to_new.resize(totalDOF, 0);

      kk = 0;
      for(ii=0; ii<nNode; ii++)
      {
        ind1 = node_map_new_to_old[ii]*ndof;
        ind2 = node_map_old_to_new[ii]*ndof;

        for(jj=0; jj<ndof; jj++)
        {
          dof_map_new_to_old[kk] = ind1 + jj;
          dof_map_old_to_new[kk] = ind2 + jj;
          kk++;
        }
      }

      row_start = 0;
      row_end   = locally_owned_nodes_total[0].size()*ndof - 1;

      for(subdomain=1; subdomain<=this_mpi_proc; subdomain++)
      {
        row_start += locally_owned_nodes_total[subdomain-1].size()*ndof;
        row_end   += locally_owned_nodes_total[subdomain].size()*ndof;
      }

      ndofs_local = row_end - row_start + 1;


      /////////////////////////////////////////////////////////////
      // 
      // reorder element node numbers, and
      // the the associated DOFs
      /////////////////////////////////////////////////////////////

      MPI_Barrier(MPI_COMM_WORLD);

      for(e1=0; e1<nElem; e1++)
      {
        ee = activeElements[e1];
        //if(elems[ee]->get_subdomain_id() == this_mpi_proc)
        //{
          for(ii=0; ii<elems[ee]->GlobalBasisFuncs.size(); ii++)
          {
            elems[ee]->GlobalBasisFuncs[ii] = node_map_old_to_new[elems[ee]->GlobalBasisFuncs[ii]];
          }
        //}
      }
    } // if(n_mpi_procs > 1)


    for(ee=0; ee<nElem; ee++)
    {
      elems[activeElements[ee]]->initialiseDOFvalues();
    }

    PetscSynchronizedPrintf(MPI_COMM_WORLD, "\n    preparing matrix pattern    \n");
    //printf("\n element DOF values initialised \n\n");
    //printf("\n Finding Global positions \n\n");

    vector<vector<int> >  DDconnLoc;

    tempDOF  = gridBF1 * ndof;

    DDconnLoc.resize(tempDOF);

    //cout << " tempDOF =  " << tempDOF << endl;

    for(e=0; e<activeElements.size(); e++)
    {
        nd1 = elems[activeElements[e]];

        val1 =  nd1->GetNsize();
        tt1  =  &(nd1->forAssyVec[0]);

        if( nd1->domNums[0] == 0 )
        {
          for(ii=0; ii<val1; ii++)
          {
            r = tt1[ii];

            for(jj=0; jj<val1; jj++)
              DDconnLoc[r].push_back(tt1[jj]);
          }
        }

      // connectivity for ghost-penalty terms
      if( nd1->IsCutElement() )
      {
        for(side=0; side<NUM_NEIGHBOURS; side++)
        {
          nd2 = nd1->GetNeighbour(side);

          if( (nd2 != NULL) && !(nd2->IsGhost()) && nd2->IsLeaf() && (nd2->IsCutElement() || nd2->domNums[0] == 0) )
          {
              nr1 = nd1->forAssyVec.size();
              nr2 = nd2->forAssyVec.size();

              for(ii=0; ii<nr1; ii++)
              {
                r = nd1->forAssyVec[ii];

                for(jj=0; jj<nr2; jj++)
                {
                  c = nd2->forAssyVec[jj];

                  DDconnLoc[r].push_back(c);
                  DDconnLoc[c].push_back(r);
                } // for(jj=0; jj<nr2; jj++)
              } // for(ii=0; ii<nr1; ii++)
          }
        } //for(side=0; side<NUM_NEIGHBOURS; side++)
      } //if( nd1->IsCutElement() )
    } // for(e=0;e<activeElements.size();e++)

    PetscSynchronizedPrintf(MPI_COMM_WORLD, "\n    Global positions DONE for individual domains \n\n");

    forAssyCutFEMprev = forAssyCutFEM;

    forAssyCutFEM2.clear();
    forAssyCutFEM.assign(tempDOF, -1);


    fluidDOF = 0;
    for(ee=0; ee<tempDOF; ee++)
    {
      findUnique(DDconnLoc[ee]);

      if(DDconnLoc[ee].size() > 0)
      {
        forAssyCutFEM[ee]  = fluidDOF++;
        forAssyCutFEM2.push_back(ee);
      }
      else
        forAssyCutFEM[ee] = -1;
    }

    // for monolithic scheme add the DOF of the solids to the global DOF
    solidDOF=0;

    totalDOF = fluidDOF + solidDOF;

    //printf("\n \t   Fluid DOF in the model   =  %8d\n\n", fluidDOF);
    //printf("\n \t   Solid DOF in the model   =  %8d\n\n", solidDOF);
    PetscSynchronizedPrintf(MPI_COMM_WORLD, "\n \t   Total DOF in the model   =  %8d\n\n", totalDOF);


    vector<vector<int> > DDconn;

    DDconn.resize(totalDOF);

    for(ee=0; ee<DDconnLoc.size(); ee++)
    {
      if(forAssyCutFEM[ee] != -1)
      {
        val1 = DDconnLoc[ee].size();
        //cout << " val1 " << ee << '\t' << val1 << endl;
        r = forAssyCutFEM[ee];
        for(ii=0; ii<val1; ii++)
        {
          DDconn[r].push_back( forAssyCutFEM[DDconnLoc[ee][ii]] );
        }
      }
    } //for(ee=0; ee<DDconnLoc[dd].size(); ee++)

    //printf("\n Finding Global positions DONE \n\n");

    VectorXi  nnzVec(totalDOF);

    nnz = 0;
    for(ii=0;ii<totalDOF;ii++)
    {
      findUnique(DDconn[ii]);
      nnzVec[ii] = DDconn[ii].size();
      nnz += nnzVec[ii];
    }
    //cout << " nnz " << nnz << endl;

    tend = time(0); 
    //cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;


    SolnData.node_map_new_to_old = node_map_new_to_old;
    SolnData.node_map_old_to_new = node_map_old_to_new;

    GeomData.node_map_new_to_old = node_map_new_to_old;
    GeomData.node_map_old_to_new = node_map_old_to_new;


      ////////////////////////////////////////////////////////////////////////////////////////////
      // find the number of nonzeroes in the diagonal and off-diagonal portions of each processor
      ////////////////////////////////////////////////////////////////////////////////////////////

      cout << " local sizes " << row_start << '\t' << row_end << '\t' << ndofs_local << endl;

      PetscInt  count_diag=0, count_offdiag=0, tempInt;
      PetscInt  d_nnz[ndofs_local], o_nnz[ndofs_local];

      kk=0;
      for(ii=row_start; ii<=row_end; ii++)
      {
        count_diag=0, count_offdiag=0;
        for(jj=0; jj<DDconn[ii].size(); jj++)
        {
          tempInt = DDconn[ii][jj];

          if(tempInt >= row_start && tempInt <= row_end)
            count_diag++;
          else
            count_offdiag++;
        }
        d_nnz[kk] = count_diag;
        o_nnz[kk] = count_offdiag;
        kk++;
      }

      //Create parallel matrix, specifying only its global dimensions.
      //When using MatCreate(), the matrix format can be specified at
      //runtime. Also, the parallel partitioning of the matrix is
      //determined by Petsc at runtime.
      //Performance tuning note: For problems of substantial size,
      //preallocation of matrix memory is crucial for attaining good
      //performance. See the matrix chapter of the users manual for details.


      ierr = MatCreate(PETSC_COMM_WORLD, &solverPetsc->mtx);CHKERRQ(ierr);

      ierr = MatSetSizes(solverPetsc->mtx, ndofs_local, ndofs_local, totalDOF, totalDOF);CHKERRQ(ierr);

      ierr = MatSetFromOptions(solverPetsc->mtx);CHKERRQ(ierr);

      ierr = MatMPIAIJSetPreallocation(solverPetsc->mtx, 20, d_nnz, 50, o_nnz);CHKERRQ(ierr);

      ierr = MatSeqAIJSetPreallocation(solverPetsc->mtx, 50, d_nnz);CHKERRQ(ierr);

      for(ii=row_start; ii<=row_end; ii++)
      {
        for(jj=0; jj<DDconn[ii].size(); jj++)
        {
          ierr = MatSetValue(solverPetsc->mtx, ii, DDconn[ii][jj], 0.0, INSERT_VALUES);
        }
      }


      VecCreate(PETSC_COMM_WORLD, &solverPetsc->soln);
      VecCreate(PETSC_COMM_WORLD, &solverPetsc->solnPrev);
      VecCreate(PETSC_COMM_WORLD, &solverPetsc->rhsVec);
      VecCreate(PETSC_COMM_WORLD, &solverPetsc->reac);

      ierr = VecSetSizes(solverPetsc->soln, ndofs_local, totalDOF); CHKERRQ(ierr);
      //ierr = VecSetSizes(solnPrev, ndofs_local, totalDOF); CHKERRQ(ierr);
      //ierr = VecSetSizes(rhsVec, ndofs_local, totalDOF); CHKERRQ(ierr);
      ierr = VecSetSizes(solverPetsc->reac, ndofs_local, totalDOF); CHKERRQ(ierr);

      ierr = VecSetFromOptions(solverPetsc->soln);CHKERRQ(ierr);
      ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->rhsVec);CHKERRQ(ierr);
      ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->solnPrev);CHKERRQ(ierr);
      ierr = VecSetFromOptions(solverPetsc->reac);CHKERRQ(ierr);

      solverPetsc->currentStatus = PATTERN_OK;

    soln.resize(totalDOF);
    soln.setZero();
    solnInit = soln;

    GRID_CHANGED = IB_MOVED = false;

    //printf("\n     HBSplineCutFEM::prepareMatrixPattern()  .... FINISHED ...\n\n");

  return 1;
}
//





/*
int HBSplineCutFEM::prepareMatrixPattern()
{
    // subroutine for Eigen based solver

    cout << "  HBSplineCutFEM::prepareMatrixPattern() ... STARTED ... " << endl;

    time_t tstart, tend; 

    int  tempDOF, domTemp;
    int  r, c, r1, c1, count=0, count1=0, count2=0, ii, jj, e, ee, dd;
    int  *tt1, *tt2, val1, val2, nnz, n1, n2, kk, e1, a, b, ll, pp, aa, bb;
    int  side, NUM_NEIGHBOURS=2*DIM, start1, start2, nr1, nr2;

    nElem = activeElements.size();
    nNode = gridBF1;

    elem_start = 0;
    elem_end   = nElem-1;

    nElem_local = nElem;

    totalDOF = nNode*ndof;

    row_start = 0;
    row_end   = totalDOF-1;

    ndofs_local = totalDOF;

    node_map_new_to_old.resize(nNode, 0);
    node_map_old_to_new.resize(nNode, 0);

    dof_map_new_to_old.resize(totalDOF, 0);
    dof_map_old_to_new.resize(totalDOF, 0);

    kk=0;
    for(ii=0; ii<nNode; ii++)
    {
      node_map_new_to_old[ii] = ii;
      node_map_old_to_new[ii] = ii;

      for(jj=0; jj<ndof; jj++)
      {
        dof_map_new_to_old[kk] = kk;
        dof_map_old_to_new[kk] = kk;
        kk++;
      }
    }

    SolnData.node_map_new_to_old = node_map_new_to_old;
    SolnData.node_map_old_to_new = node_map_old_to_new;

    GeomData.node_map_new_to_old = node_map_new_to_old;
    GeomData.node_map_old_to_new = node_map_old_to_new;


    printf("\n Finding Global positions \n\n");

    vector<vector<int> >  DDconnLoc;

    tempDOF  = gridBF1 * ndof;

    DDconnLoc.resize(tempDOF);

    //cout << " activeElements.size() " << activeElements.size() << '\t' << tempDOF << endl;

    node  *nd1, *nd2;

    tstart = time(0);

    prepareCutElements();

    tend = time(0);
    printf("HBSplineCutFEM::prepareCutElements() took %8.4f second(s) \n ", difftime(tend, tstart) );

    for(e=0; e<activeElements.size(); e++)
    {
      elems[activeElements[e]]->initialiseDOFvalues();
    }

    for(e=0; e<activeElements.size(); e++)
    {
        nd1 = elems[activeElements[e]];

        val1 =  nd1->GetNsize();
        tt1  =  &(nd1->forAssyVec[0]);
        //cout << e << '\t' << val1 << endl;

        //printVector(nd1->forAssyVec);

        if( nd1->domNums[0] == 0 )
        {
          for(ii=0; ii<val1; ii++)
          {
            r = tt1[ii];

            for(jj=0; jj<val1; jj++)
              DDconnLoc[r].push_back(tt1[jj]);
          }
        }

      // connectivity for ghost-penalty terms
      //
      if( nd1->IsCutElement() )
      {
        for(side=0; side<NUM_NEIGHBOURS; side++)
        {
          nd2 = nd1->GetNeighbour(side);

          if( (nd2 != NULL) && !(nd2->IsGhost()) && nd2->IsLeaf() && (nd2->IsCutElement() || nd2->domNums[0] == 0) )
          {
              //cout << " side = " << side << endl;

              nr1 = nd1->forAssyVec.size();
              nr2 = nd2->forAssyVec.size();

              for(ii=0; ii<nr1; ii++)
              {
                r = nd1->forAssyVec[ii];

                for(jj=0; jj<nr2; jj++)
                {
                  c = nd2->forAssyVec[jj];

                  DDconnLoc[r].push_back(c);
                  DDconnLoc[c].push_back(r);
                } // for(jj=0; jj<nr2; jj++)
              } // for(ii=0; ii<nr1; ii++)
          }
        } //for(side=0; side<NUM_NEIGHBOURS; side++)
      } //if( nd1->IsCutElement() )
      //
    } // for(e=0;e<activeElements.size();e++)

    printf("\n Global positions DONE for individual domains \n\n");

    forAssyCutFEMprev = forAssyCutFEM;

    forAssyCutFEM2.clear();
    forAssyCutFEM.assign(tempDOF, -1);


    fluidDOF = 0;
    for(ee=0; ee<tempDOF; ee++)
    {
      findUnique(DDconnLoc[ee]);

      if(DDconnLoc[ee].size() > 0)
      {
        forAssyCutFEM[ee]  = fluidDOF++;
        forAssyCutFEM2.push_back(ee);
      }
      else
        forAssyCutFEM[ee] = -1;
    }

    bool pp1=false;
    //pp1=true;
    if(pp1)
    {
        printf("   dof to dof connectivity ...:  \n\n");
        for(ii=0;ii<totalDOF;ii++)
        {
          cout << " dof # " << ii << " : ";
          for(jj=0;jj<DDconnLoc[ii].size();jj++)
            cout << '\t' << DDconnLoc[ii][jj];
          cout << endl;
        }
        printf("\n\n\n\n\n\n");

        printf("   forAssyCutFEM ..... : \n\n");
        for(ii=0; ii<forAssyCutFEM.size(); ii++)
          cout << ii << '\t' << forAssyCutFEM[ii] << endl;

        printf("\n\n\n\n\n\n");

        printf("   forAssyCutFEM2 ..... : \n\n");
        for(ii=0;ii<forAssyCutFEM2.size();ii++)
          cout << ii << '\t' << forAssyCutFEM2[ii] << endl;

        printf("\n\n\n\n\n\n");
    }


    // for monolithic scheme add the DOF of the solids to the global DOF
    // 
    // hardcoded for a single DOF rigid body
    // 
    solidDOF=0;

    totalDOF = fluidDOF + solidDOF;

    //printf("\n \t   Fluid DOF in the model   =  %8d\n\n", fluidDOF);
    //printf("\n \t   Solid DOF in the model   =  %8d\n\n", solidDOF);
    printf("\n \t   Total DOF in the model   =  %8d\n\n", totalDOF);


    vector<vector<int> > DDconn;

    DDconn.resize(totalDOF);

    for(ee=0; ee<DDconnLoc.size(); ee++)
    {
      if(forAssyCutFEM[ee] != -1)
      {
        val1 = DDconnLoc[ee].size();
        //cout << " val1 " << ee << '\t' << val1 << endl;
        r = forAssyCutFEM[ee];
        for(ii=0; ii<val1; ii++)
        {
          DDconn[r].push_back( forAssyCutFEM[DDconnLoc[ee][ii]] );
        }
      }
    } //for(ee=0; ee<DDconnLoc[dd].size(); ee++)


    if(!STAGGERED)
    {
      start1 = fluidDOF;
      start2 = fluidDOF;
      domTemp = 0;
      nr1 = 1;

      for(e=0; e<activeElements.size(); e++)
      {
          nd1 = elems[activeElements[e]];

          if( nd1->IsCutElement() )
          {
              nr2 = nd1->forAssyVec.size();
              tt2 = &(nd1->forAssyVec[0]);

              for(ii=0; ii<nr1; ii++)
              {
                r = start1 + ii;

                for(jj=0; jj<nr2; jj++)
                {
                  c = tt2[jj];

                  DDconn[r].push_back(c);
                  DDconn[c].push_back(r);
                } // for(jj=0; jj<nr2; jj++)
              
                DDconn[r].push_back(r);
              } // for(ii=0; ii<nr1; ii++)
          } //if( nd1->IsCutElement() )
      } // for(e=0;e<activeElements.size();e++)
    } // if(!STAGGERED)

    //printf("\n Finding Global positions DONE \n\n");

    VectorXi  nnzVec(totalDOF);

    nnz = 0;
    for(ii=0;ii<totalDOF;ii++)
    {
      findUnique(DDconn[ii]);
      nnzVec[ii] = DDconn[ii].size();
      nnz += nnzVec[ii];
    }
    //cout << " nnz " << nnz << endl;

    tend = time(0); 
    //cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;

    pp1=false;
    //pp1=true;
    if(pp1)
    {
      printf("   Number of non-zeros = %5d \n\n", nnz);
      printf("   dof to dof connectivity ...:  \n\n");
      for(ii=0;ii<totalDOF;ii++)
      {
        cout << " dof # " << ii << " : ";
        for(jj=0;jj<DDconn[ii].size();jj++)
          cout << '\t' << DDconn[ii][jj];
        cout << endl;
      }
      printf("\n\n\n");
    }

    //cout << " AAAAAAAAAA " << endl;
    solverEigen->mtx.setZero();
    solverEigen->mtx.uncompress();

    //solverEigen->mtx.resize(totalDOF, totalDOF);
    solverEigen->mtx.conservativeResize(totalDOF, totalDOF);
    //solverEigen->mtx.reserve(nnz);
    solverEigen->mtx.reserve(nnzVec);
    //cout << " AAAAAAAAAA " << endl;

    FILE * pFile;
    char name [100];

    //pFile = fopen ("pattern.dat","w");

    typedef Eigen::Triplet<double> T;

    vector<T> tripletList;

    tripletList.reserve(nnz);

    for(ii=0;ii<DDconn.size();ii++)
    {
      for(jj=0;jj<DDconn[ii].size();jj++)
      {
        tripletList.push_back(T(ii, DDconn[ii][jj], 0.0));
      }
    }

    solverEigen->mtx.setFromTriplets(tripletList.begin(), tripletList.end());

    // mat is ready to go!

    //fclose(pFile);
    //solverEigen->mtx.data().squeeze();

    solverEigen->mtx.makeCompressed();

    soln.resize(totalDOF);
    soln.setZero();
    solnInit = soln;

    //cout << solverEigen->mtx.rows() << endl;
    //cout << solverEigen->mtx.cols() << endl;
    //cout << solverEigen->mtx.nonZeros() << endl;

    solverEigen->currentStatus = PATTERN_OK;

    GRID_CHANGED = IB_MOVED = false;

    printf("\n     HBSplineCutFEM::prepareMatrixPattern()  .... FINISHED ...\n\n");

  return 1;
}
*/




/*
void HBSplineCutFEM::prepareMatrixPatternCutFEM()
{
  // subroutine for multiple active domains

    time_t tstart, tend; 

    computerTime.go(fct);

    int  IBPsize, ndf2, ff, ind1, ind2, tempDOF;
    int  r, c, r1, c1, count=0, count1=0, count2=0, ii, jj, e, ee;
    int  *tt1, *tt2, val1, val2, nnz, n1, n2, kk, e1, a, b, ll, pp, aa, bb;

    printf("\n element DOF values initialised \n\n");
    printf("\n Finding Global positions \n\n");

    vector<vector<int> >  DDconn1, DDconn2;

    tempDOF = totalDOF;

    DDconn1.resize(tempDOF);
    DDconn2.resize(tempDOF);

    // var1 and var2 are set of size tempDOF so that it would be easier to compute solution variables
    // at Gausspoints while performing numerical intergration and also for postprocessing
    //

    SolnData.initialise(tempDOF, tempDOF, 0, 0);

    SolnData.SetSolidOrFluid(0);
    SolnData.SetTimeIncrementType(tis);
    SolnData.SetRho(td[0]);
    SolnData.STAGGERED = STAGGERED;

    tstart = time(0);

    cout << " activeElements.size() " << activeElements.size() << '\t' << totalDOF << endl;
    for(e=0;e<activeElements.size();e++)
    {
      ee = activeElements[e];

      val1 =  elems[ee]->GetNsize();
      tt1  =  &(elems[ee]->forAssyVec[0]);
      //cout << e << '\t' << ee << '\t' << val1 << endl;
      //printVector(elems[ee]->forAssyVec);

      elems[ee]->findCutElemType();

      if( elems[ee]->IsCutElement() )
      {
          for(ii=0;ii<val1;ii++)
          {
            r = tt1[ii];

            for(jj=0;jj<val1;jj++)
            {
              DDconn1[r].push_back(tt1[jj]);
              DDconn2[r].push_back(tt1[jj]);
            }
          }
      } // if( !elems[ee]->IsCutElement() )
      else
      {
        if( elems[ee]->GetDomainNumber() == 0 )
        {
          for(ii=0;ii<val1;ii++)
          {
            r = tt1[ii];

            for(jj=0;jj<val1;jj++)
              DDconn1[r].push_back(tt1[jj]);
          }
        }
        else
        {
          for(ii=0;ii<val1;ii++)
          {
            r = tt1[ii];

            for(jj=0;jj<val1;jj++)
              DDconn2[r].push_back(tt1[jj]);
          }
        }
      } // else
    } // for(e=0;e<activeElements.size();e++)

    printf("\n Global positions DONE for individual domains \n\n");

    domainTotalDOF.resize(2);
    forAssyCutFEM.resize(2);
    forAssyCutFEM2.resize(2);
    domainTotalDOF[0] = domainTotalDOF[1] = 0;

    for(ee=0; ee<tempDOF; ee++)
    {
      findUnique(DDconn1[ee]);
      if(DDconn1[ee].size() > 0)
      {
        forAssyCutFEM[0].push_back(domainTotalDOF[0]++);
        forAssyCutFEM2[0].push_back(ee);
      }
      else
        forAssyCutFEM[0].push_back(-1);

      findUnique(DDconn2[ee]);
      if(DDconn2[ee].size() > 0)
      {
        forAssyCutFEM[1].push_back(domainTotalDOF[1]++);
        forAssyCutFEM2[1].push_back(ee);
      }
      else
        forAssyCutFEM[1].push_back(-1);
    }

    cout << " domainDOF " << domainTotalDOF[0] << '\t' << domainTotalDOF[1] << endl;

    bool pp1=false;
    //pp1=true;
    if(pp1)
    {
       printf("   dof to dof connectivity -> Domain #1 ...:  \n\n");
       for(ii=0;ii<totalDOF;ii++)
       {
          cout << " dof # " << ii << " : ";
          for(jj=0;jj<DDconn1[ii].size();jj++)
            cout << '\t' << DDconn1[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("   dof to dof connectivity -> Domain #2 ...:  \n\n");
       for(ii=0;ii<totalDOF;ii++)
       {
          cout << " dof # " << ii << " : ";
          for(jj=0;jj<DDconn2[ii].size();jj++)
            cout << '\t' << DDconn2[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("   forAssyCutFEM[0] ..... : \n\n");
       for(ii=0;ii<forAssyCutFEM[0].size();ii++)
          cout << ii << '\t' << forAssyCutFEM[0][ii] << endl;

       printf("\n\n\n");

       printf("   forAssyCutFEM[1] ..... : \n\n");
       for(ii=0;ii<forAssyCutFEM[1].size();ii++)
          cout << ii << '\t' << forAssyCutFEM[1][ii] << endl;

       printf("\n\n\n");

       printf("   forAssyCutFEM2[0] ..... : \n\n");
       for(ii=0;ii<forAssyCutFEM2[0].size();ii++)
          cout << ii << '\t' << forAssyCutFEM2[0][ii] << endl;

       printf("\n\n\n");

       printf("   forAssyCutFEM2[1] ..... : \n\n");
       for(ii=0;ii<forAssyCutFEM2[1].size();ii++)
          cout << ii << '\t' << forAssyCutFEM2[1][ii] << endl;

       printf("\n\n\n");

    }

    velDOF  = domainTotalDOF[0];
    presDOF = domainTotalDOF[1];



    totalDOF = velDOF + presDOF;

    printf("\n \t   Total DOF in domain #1    =  %5d\n\n", velDOF);
    printf("\n \t   Total DOF in domain #2    =  %5d\n\n", presDOF);
    printf("\n \t   Total number of DOF       =  %5d\n\n", totalDOF);

    //totalDOF = domainTotalDOF[0] + domainTotalDOF[1];

    DDconn.resize(totalDOF);


    for(ee=0; ee<tempDOF; ee++)
    {
      val1 = DDconn1[ee].size();
      if(val1 > 0)
      {
        //cout << " val1 " << ee << '\t' << val1 << endl;
        r = forAssyCutFEM[0][ee];
        for(ii=0; ii<val1; ii++)
        {
          DDconn[r].push_back(forAssyCutFEM[0][DDconn1[ee][ii]]);
        }
      }

      val1 = DDconn2[ee].size();
      kk = domainTotalDOF[0];
      if(val1 > 0)
      {
        r = kk + forAssyCutFEM[1][ee] ;
        //cout << " val1 " << ee << '\t' << r << endl;
        for(ii=0; ii<val1; ii++)
        {
          DDconn[r].push_back(kk+forAssyCutFEM[1][DDconn2[ee][ii]]);
        }
      }
    } //for(ee=0; ee<tempDOF; ee++)


    // generate matrix pattern for the coupling terms
    //
    for(e=0;e<activeElements.size();e++)
    {
      ee = activeElements[e];

      if( elems[ee]->IsCutElement() )
      {
        val1 =  elems[ee]->GetNsize();
        tt1  =  &(elems[ee]->forAssyVec[0]);
        //cout << e << '\t' << ee << '\t' << val1 << endl;
        //printVector(elems[ee]->forAssyVec);

        kk = domainTotalDOF[0];

        for(ii=0;ii<val1;ii++)
        {
          r = forAssyCutFEM[0][tt1[ii]];

          for(jj=0;jj<val1;jj++)
          {
            c = kk + forAssyCutFEM[1][tt1[jj]] ;
            DDconn[r].push_back(c);
            DDconn[c].push_back(r);
          }
        }
      } // if( !elems[ee]->IsCutElement() )
    } //for(e=0;e<activeElements.size();e++)


    printf("\n Finding Global positions DONE \n\n");

    VectorXd  nnzVec(totalDOF);

    nnz = 0;
    for(ii=0;ii<totalDOF;ii++)
    {
      findUnique(DDconn[ii]);
      nnzVec[ii] = DDconn[ii].size();
      nnz += nnzVec[ii];
    }
    cout << " nnz " << nnz << endl;

    tend = time(0); 
    cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;

    pp1=false;
    //pp1=true;
    if(pp1)
    {
       printf("   Number of non-zeros = %5d \n\n", nnz);
       printf("   dof to dof connectivity ...:  \n\n");
       for(ii=0;ii<totalDOF;ii++)
       {
          cout << " dof # " << ii << " : ";
          for(jj=0;jj<DDconn[ii].size();jj++)
            cout << '\t' << DDconn[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");
    }

    cout << " AAAAAAAAAA " << endl;


    cout << " AAAAAAAAAA " << endl;
    //MatCreateSeqAIJ(PETSC_COMM_WORLD, totalDOF, totalDOF, 500, nnzVec, &(((SolverEigen*)solver)->mtx));
    //solverEigen->mtx.setZero();

    solverEigen->mtx.resize(totalDOF, totalDOF);
    solverEigen->mtx.reserve(nnz);
    solverEigen->mtx.reserve(nnzVec);
    cout << " AAAAAAAAAA " << endl;

    for(ii=0;ii<totalDOF;ii++)
    {
      for(jj=0;jj<DDconn[ii].size();jj++)
      {
        solverEigen->mtx.coeffRef(ii, DDconn[ii][jj]) = 0.0;
      }
    }

    solverEigen->mtx.makeCompressed();

    soln.resize(totalDOF);
    soln.setZero();
    solnInit = soln;

  //printVector(solver->rhsVec);
  //cout << " totalDOF = " << totalDOF << endl;
  //cout << solver->mtx << endl;

  //solverEigen->printMatrix();

  solverEigen->currentStatus = PATTERN_OK;

  GRID_CHANGED = IB_MOVED = false;

  //printf("\n     HBSplineCutFEM::prepareMatrixPattern()  .... FINISHED ...\n\n");

  return;
}
*/




void  HBSplineCutFEM::prepareMatrixPatternPostProcess()
{
    ///////////////////////////// 
    // matrix pattern for the matrix to compute vorticity
    ///////////////////////////// 

    int  *tt1, r, c, ii, jj, e, ee, val1, val2, dd, nnz;

    dd = gridBF1;

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




