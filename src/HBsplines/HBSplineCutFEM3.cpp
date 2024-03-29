
#include "HBSplineCutFEM.h"
#include "ImmersedIntegrationElement.h"
#include "QuadratureUtil.h"

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
  //PetscPrintf(MPI_COMM_WORLD, "   HBSplineCutFEM::timeUpdate() ... STARTED \n\n");

  //double tstart = MPI_Wtime();

  IB_MOVED = false;

  SOLID_SOLVED = 0;

  firstIter = true;
  localStiffnessError = 0;
  iterCount = 1;
  filecount++;

  for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
    ImmersedBodyObjects[bb]->timeUpdate();

  SolnData.timeUpdate();

  //if(STAGGERED)
    //solveSolidProblem();

  MPI_Barrier(MPI_COMM_WORLD);

  //updateIterStep();

  //double tend = MPI_Wtime();
  //PetscPrintf(MPI_COMM_WORLD, " HBSplineCutFEM::timeUpdate() took %f  milliseconds \n", (tend-tstart)*1000);

  //PetscPrintf(MPI_COMM_WORLD, "   HBSplineCutFEM::timeUpdate() ... FINISHED \n\n");

  return;
}



void HBSplineCutFEM::updateIterStep()
{
  //PetscPrintf(MPI_COMM_WORLD, "   HBSplineCutFEM::updateIterStep() ... STARTED \n\n");

  PetscLogDouble mem1, mem2, mem3, mem4;

  //the current resident set size (memory used) for the program.
  //ierr = PetscMemoryGetCurrentUsage(&mem1);   //CHKERRQ(ierr);
  //the maximum resident set size (memory used) for the program.
  //ierr = PetscMemoryGetMaximumUsage(&mem2);   //CHKERRQ(ierr);
  //the current amount of memory used that was PetscMalloc()ed
  //ierr = PetscMallocGetCurrentUsage(&mem3);   //CHKERRQ(ierr);
  //the maximum amount of memory used that was PetscMalloc()ed at any time during this run.
  //ierr = PetscMallocGetMaximumUsage(&mem4);   //CHKERRQ(ierr);

  //PetscPrintf(MPI_COMM_WORLD, " Petsc memory allocation details ... %12.8f \t %12.8f \t%12.8f \t%12.8f \n\n", mem1, mem2, mem3, mem4);

  int kk, bb, ee;

  SolnData.updateIterStep();


  if(!STAGGERED)
  {
    IB_MOVED = false;
    for(bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
      //cout << " AAAAAAAAAAA " << bb << endl;
      ImmersedBodyObjects[bb]->updateIterStep();
      //cout << " AAAAAAAAAAA " << bb << endl;
      IB_MOVED = (IB_MOVED || ImmersedBodyObjects[bb]->updatePointPositions() );
    }
    //cout << " AAAAAAAAAAA " << bb << endl;
  }

  //cout << " IB_MOVED = " << IB_MOVED << '\t' << slv_type << endl;
  //cout << " kkkkkkkkkkk " << endl;
  if(GRID_CHANGED || IB_MOVED)
  {
      if(SOLVER_TYPE == SOLVER_TYPE_PETSC)
        solverPetsc->free();
      else
        solverEigen->free();

      double tstart = MPI_Wtime();

      prepareMatrixPattern();

      double tend = MPI_Wtime();
      PetscPrintf(MPI_COMM_WORLD, " prepareMatrixPattern() took %f  milliseconds \n", (tend-tstart)*1000);

      switch(slv_type)
      {
        case  8: // SolverPetsc ..........................

            if(solverPetsc->initialise(0, 0, totalDOF) != 0)
              return;

            solverPetsc->setSolverAndParameters();
            //solverPetsc->printInfo();

        break;

        case  9: // PARDISO(sym) with Petsc
        case  10: // PARDISO(unsym) with Petsc

            if(slv_type == 9)
            {
              if(solverPetsc->initialise(numProc, PARDISO_STRUCT_SYM, totalDOF) != 0)
                return;
            }
            if(slv_type == 10)
            {
              if(solverPetsc->initialise(numProc, PARDISO_UNSYM, totalDOF) != 0)
                return;
            }

            solverPetsc->setSolverAndParameters();

        break;

          default: // invalid SOLVER_TYPE ...................

            cout << " this solver has not been implemented yet!\n\n";

          break;
      }

      solverOK = true;

      if(solverPetsc != NULL)
        solverPetsc->checkIO = true;
  }

  //PetscPrintf(MPI_COMM_WORLD, "   HBSplineCutFEM::updateIterStep() ... FINISHED \n\n");

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
    if(ImmersedBodyObjects.size() == 0)
      return 1;

    int  ttt=0;
    if(DIM == 1)
      ttt = setCoveringUncovering1D();
    else if(DIM == 2)
      ttt = setCoveringUncovering2D();
    else
      ttt = setCoveringUncovering3D();

  return ttt;
}



int  HBSplineCutFEM::setCoveringUncovering1D()
{
  return 1;
}


int  HBSplineCutFEM::setCoveringUncovering2D()
{
    int  aa, bb=0, dd, ee, ii;
    double  vel[3];

    int nlocal = (degree[0]+1) * (degree[1] + 1);

    vector<int>  bfTemp;
    VectorXd  NN(nlocal), N(nlocal);
    myPoint  knotIncr, knotBegin;

    node* ndTemp;


        for(aa=0; aa<ImmersedBodyObjects[bb]->GeomData.NodePosNew.size(); aa++)
        {
          geom    = ImmersedBodyObjects[bb]->GeomData.NodePosNew[aa];

          vel[0]  = ImmersedBodyObjects[bb]->GeomData.specValNew[aa][0];
          vel[1]  = ImmersedBodyObjects[bb]->GeomData.specValNew[aa][1];

          //printf(" %12.6f,  %12.6f \n\n", geom[0], geom[1]);
          //printf(" %12.6f,  %12.6f \n", vel[0], vel[1]);

          //cout << " uuuuuuuuuuu " << endl;

          ndTemp = elems[findCellNumber(geom)];

          if( ndTemp->getDomainNumber() > 0 )
          {
            geometryToParametric(geom, param);

            knotBegin = ndTemp->getKnotBegin();
            knotIncr  = ndTemp->getKnotIncrement();
            bfTemp    = ndTemp->GlobalBasisFuncs;

            GeomData.computeBasisFunctions2D(knotBegin, knotIncr, param, NN);

            if(ndTemp->getParent() == NULL)
            {
              N = NN;
            }
            else
            {
              N = ndTemp->SubDivMat*NN;
            }

            for(ii=0; ii<bfTemp.size(); ii++)
            {
              if( grid_to_cutfem_BF[bfTemp[ii]] == -1 )
              {
                for(dd=0; dd<DIM; dd++)
                {
                  SolnData.var1[bfTemp[ii]*ndof+dd] = vel[dd];
                }
              }
            }
          }
        }

  return  1;
}



int  HBSplineCutFEM::setCoveringUncovering3D()
{
    int  aa, bb=0, dd, ee, ii;
    double  vel[3];

    int nlocal = (degree[0]+1) * (degree[1] + 1) * (degree[2] + 1);

    vector<int>  bfTemp;
    VectorXd  NN(nlocal), N(nlocal);
    myPoint  knotIncr, knotBegin;

    node* ndTemp;


    for(aa=0; aa<ImmersedBodyObjects[bb]->GeomData.NodePosNew.size(); aa++)
    {
        geom    = ImmersedBodyObjects[bb]->GeomData.NodePosNew[aa];

        vel[0]  = ImmersedBodyObjects[bb]->GeomData.specValNew[aa][0];
        vel[1]  = ImmersedBodyObjects[bb]->GeomData.specValNew[aa][1];
        vel[2]  = ImmersedBodyObjects[bb]->GeomData.specValNew[aa][2];

        //printf(" %12.6f,  %12.6f \n\n", geom[0], geom[1]);
        //printf(" %12.6f,  %12.6f \n", vel[0], vel[1]);

        //cout << " uuuuuuuuuuu " << endl;

        ndTemp = elems[findCellNumber(geom)];

        if( ndTemp->getDomainNumber() > 0 )
        {
            geometryToParametric(geom, param);

            knotBegin = ndTemp->getKnotBegin();
            knotIncr  = ndTemp->getKnotIncrement();
            bfTemp    = ndTemp->GlobalBasisFuncs;

            GeomData.computeBasisFunctions3D(knotBegin, knotIncr, param, NN);

            if(ndTemp->getParent() == NULL)
            {
              N = NN;
            }
            else
            {
              N = ndTemp->SubDivMat*NN;
            }

            for(ii=0; ii<bfTemp.size(); ii++)
            {
              if( grid_to_cutfem_BF[bfTemp[ii]] == -1 )
              {
                for(dd=0; dd<DIM; dd++)
                {
                  SolnData.var1[bfTemp[ii]*ndof+dd] = vel[dd];
                }
              }
            }
        }
    }

    return  1;
}




int  HBSplineCutFEM::prepareMatrixPattern()
{
    //for PETSc - working

    assert(SOLVER_TYPE == SOLVER_TYPE_PETSC);

    int  tempDOF, domTemp, npElem, ind, size1;
    int  r, c, r1, c1, count=0, count1=0, count2=0, ii, jj, ee, dd, ind1, ind2;
    int  *tt1, *tt2, val1, val2, nnz, nnz_max_row, n1, n2, kk, e1, e2, ll, pp, aa, bb;
    int  side, NUM_NEIGHBOURS=2*DIM, start1, start2, nr1, nr2;

    node  *nd1, *nd2;

    int subdomain=0;

    PetscPrintf(MPI_COMM_WORLD, " preparing cut elements \n");

    double  tstart = MPI_Wtime();

    //
    // Find cut cells, without computng their quadrature points.
    // Because quadrature points information is not required to decompose the mesh
    // and prepare the matrix pattern.
    // Quadrature points can be computed can be computed "parallely" after
    // preparing the matrix pattern.
    // It has to be done this way because we don't know which elements belong
    // to which processor.
    ///////////////////////////////////////////////////////////////////////////

    //prepareCutElements();
    prepareCutElements2();

    MPI_Barrier(MPI_COMM_WORLD);

    double  tend = MPI_Wtime();

    PetscPrintf(MPI_COMM_WORLD, " preparing cut elements took %f  milliseconds \n", (tend-tstart)*1000);

    fluidDOF = nNode*ndof;

    PetscSynchronizedPrintf(MPI_COMM_WORLD, " nElem    = %8d \n", nElem);
    PetscSynchronizedPrintf(MPI_COMM_WORLD, " gridBF1  = %8d \n", gridBF1);
    PetscSynchronizedPrintf(MPI_COMM_WORLD, " nNode    = %8d \n", nNode);
    PetscSynchronizedPrintf(MPI_COMM_WORLD, " ndof     = %8d \n", ndof);
    PetscSynchronizedPrintf(MPI_COMM_WORLD, " fluidDOF = %8d \n", fluidDOF);
    PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT);

    //setCoveringUncovering();

    tstart = MPI_Wtime();


    // for monolithic scheme add the DOF of the solids to the global DOF
    solidDOF = 0;
    if(!STAGGERED)
    {
      for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
      {
        solidDOF += ImmersedBodyObjects[bb]->getTotalDOF();
      }
    }

    totalDOF = fluidDOF + solidDOF;

    PetscPrintf(MPI_COMM_WORLD, "\n \t   Number of processors     =  %8d\n\n", n_mpi_procs);
    PetscPrintf(MPI_COMM_WORLD, "\n \t   Fluid DOF in the model   =  %8d\n\n", fluidDOF);
    PetscPrintf(MPI_COMM_WORLD, "\n \t   Solid DOF in the model   =  %8d\n\n", solidDOF);
    PetscPrintf(MPI_COMM_WORLD, "\n \t   Total DOF in the model   =  %8d\n\n", totalDOF);

    proc_to_grid_DOF.resize(totalDOF);

    if(n_mpi_procs == 1)
    {
      elem_start = 0;
      elem_end   = nElem-1;

      nElem_local = nElem;

      row_start   = 0;
      row_end     = totalDOF-1;
      ndofs_local = totalDOF;

      bfs_start = 0;
      bfs_end   = nNode-1;
      bfs_local = nNode;

      node_map_get_old.resize(nNode, 0);
      node_map_get_new.resize(nNode, 0);

      dof_map_get_old.resize(totalDOF, 0);
      dof_map_get_new.resize(totalDOF, 0);

      kk=0;
      for(ii=0; ii<nNode; ii++)
      {
        node_map_get_old[ii] = ii;
        node_map_get_new[ii] = ii;

        for(jj=0; jj<ndof; jj++)
        {
          dof_map_get_old[kk] = kk;
          dof_map_get_new[kk] = kk;
          kk++;
        }
      }
    }
    else
    {
      /////////////////////////////////////////////////////////////////////////////
      //
      // Partition the mesh. This can be done using software libraries
      // Chaco, Jostle, METIS and Scotch, or PETSc's subroutines.
      // Here METIS is used.
      // 
      /////////////////////////////////////////////////////////////////////////////

      PetscSynchronizedPrintf(MPI_COMM_WORLD, " this_mpi_proc = %d \n", this_mpi_proc);
      PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT);

      n2 = ceil(nElem/n_mpi_procs);

      elem_start = n2*this_mpi_proc;
      elem_end   = n2*(this_mpi_proc+1)-1;

      if(this_mpi_proc == (n_mpi_procs-1))
        elem_end = nElem-1;

      nElem_local = elem_end - elem_start + 1;

      PetscInt  *elem_proc_id, *node_proc_id;

      ierr  = PetscMalloc1(nElem,    &elem_proc_id); CHKERRQ(ierr);
      ierr  = PetscMalloc1(nNode,    &node_proc_id); CHKERRQ(ierr);

      // paritioning is performed only by the first processor
      //
      if(this_mpi_proc == 0)
      {
        PetscInt  *eptr, *eind, *vwgtElem;

        ierr  = PetscMalloc1(nElem+1,  &eptr);  CHKERRQ(ierr);
        //ierr  = PetscMalloc1(nElem,    &vwgtElem); CHKERRQ(ierr);

        eptr[0] = 0;

        PetscInt npElem_total = 0;
        for(ee=0; ee<nElem; ee++)
        {
          npElem_total += elems[fluidElementIds[ee]]->GlobalBasisFuncs.size();

          eptr[ee+1] = npElem_total;

          //vwgtElem[ee] = elems[fluidElementIds[ee]]->getComputationalEffort();
          //vwgtElem[ee] = elems[fluidElementIds[ee]]->getNumberOfQuadraturePoints();
        }

        ierr  = PetscMalloc1(npElem_total,  &eind); CHKERRQ(ierr);

        vector<int>  vecTemp2;

        kk=0;
        for(e1=0; e1<nElem; e1++)
        {
          vecTemp2 = elems[fluidElementIds[e1]]->GlobalBasisFuncs ;

          npElem = vecTemp2.size();

          for(ii=0; ii<npElem; ii++)
            eind[kk+ii] = grid_to_cutfem_BF[vecTemp2[ii]] ;

          kk += npElem;
        }

        int  nodes_per_side=2;

        if(DIM == 3)
          nodes_per_side = 4;

        int  nWeights  = 1, numflag=0, objval;
        int  *xadj, *adjncy;
        int  options[METIS_NOPTIONS];

        METIS_SetDefaultOptions(options);

        // Specifies the partitioning method.
        //options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;    // Multilevel recursive bisectioning.
        options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;  // Multilevel k-way partitioning.

        //options[METIS_OPTION_NSEPS] = 10;

        //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;   // Edge-cut minimization
        options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL; // Total communication volume minimization

        options[METIS_OPTION_NUMBERING] = 0;  // C-style numbering is assumed that starts from 0.


        // METIS partition routine - based on the nodal graph of the mesh
        //int ret = METIS_PartMeshNodal(&nElem, &nNode, eptr, eind, NULL, NULL, &n_mpi_procs, NULL, options, &objval, elem_proc_id, node_proc_id);

        // METIS partition routine - based on the dual graph of the mesh
        int ret = METIS_PartMeshDual(&nElem, &nNode, eptr, eind, NULL, NULL, &nodes_per_side, &n_mpi_procs, NULL, options, &objval, elem_proc_id, node_proc_id);
        //int ret = METIS_PartMeshDual(&nElem, &nNode, eptr, eind, vwgtElem, NULL, &nodes_per_side, &n_mpi_procs, NULL, options, &objval, elem_proc_id, node_proc_id);


        if(ret == METIS_OK)
          PetscPrintf(MPI_COMM_WORLD, "   METIS partition routine successful \n");
        else
          PetscPrintf(MPI_COMM_WORLD, "   METIS partition routine FAILED \n");

        //for(ee=0; ee<nElem; ee++)
          //cout << ee << '\t' << elem_proc_id[ee] << endl;
        //cout << " \n\n\n\n " << endl;
        //for(ee=0; ee<nNode; ee++)
          //cout << ee << '\t' << node_proc_id[ee] << endl;

        ierr  = PetscFree(eptr);  CHKERRQ(ierr);
        ierr  = PetscFree(eind);  CHKERRQ(ierr);
        ierr  = PetscFree(vwgtElem); CHKERRQ(ierr);
      }  //if(this_mpi_proc == 0)

      ierr = PetscPrintf(MPI_COMM_WORLD, " After Metis \n\n\n");
      //PetscSynchronizedPrintf(MPI_COMM_WORLD, " this_mpi_proc = %d \n", this_mpi_proc);
      //PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT);
      //ierr = MPI_Barrier(MPI_COMM_WORLD);


      MPI_Request request;
      ierr = MPI_Bcast(node_proc_id, nNode, MPI_INT, rootproc, MPI_COMM_WORLD);
      ierr = PetscSynchronizedPrintf(MPI_COMM_WORLD, " After broadcasting node array on processor %d \n", this_mpi_proc);
      cout << "ierr = " << ierr << endl;
      PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT);
      ierr = MPI_Barrier(MPI_COMM_WORLD);

      ierr = MPI_Bcast(elem_proc_id, nElem, MPI_INT, rootproc, MPI_COMM_WORLD);
      ierr = PetscSynchronizedPrintf(MPI_COMM_WORLD, " After broadcasting element array on processor %d \n", this_mpi_proc);
      PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT);
      ierr = MPI_Barrier(MPI_COMM_WORLD);


      for(e1=0; e1<nElem; e1++)
        elems[fluidElementIds[e1]]->setSubdomainId(elem_proc_id[e1]);
      MPI_Barrier(MPI_COMM_WORLD);

      ierr = PetscPrintf(MPI_COMM_WORLD, " Finding local elements and nodes \n");

      // find nodes local to each processor

      nElem_local = 0;
      for(ii=0; ii<nElem; ii++)
      {
        if(elem_proc_id[ii] == this_mpi_proc)
          nElem_local++;
      }
      PetscSynchronizedPrintf(MPI_COMM_WORLD, "\t   nElem_local   =  %8d ... Processor = %d\n", nElem_local, this_mpi_proc);
      PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT);
      //MPI_Barrier(MPI_COMM_WORLD);


      std::vector<int>  nodelist_owned;

      for(ii=0; ii<nNode; ii++)
      {
        if( node_proc_id[ii] == this_mpi_proc )
          nodelist_owned.push_back(ii);
      }
      nNode_local = nodelist_owned.size();
      PetscSynchronizedPrintf(MPI_COMM_WORLD, "\t   nNode_local   =  %8d ... Processor = %d\n", nNode_local, this_mpi_proc);
      PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT);
      //MPI_Barrier(MPI_COMM_WORLD);

      // create the vector (of size n_mpi_procs)
      // consisting of nNode_owned from all the processors in the communication
      vector<int>  nNode_owned_vector(n_mpi_procs), nNode_owned_sum(n_mpi_procs);

      MPI_Allgather(&nNode_local, 1, MPI_INT, &nNode_owned_vector[0], 1, MPI_INT, MPI_COMM_WORLD);
      //printVector(nNode_owned_vector);

      // compute the numbers of first and last nodes in the local processor
      nNode_owned_sum = nNode_owned_vector;
      for(ii=1; ii<n_mpi_procs; ii++)
      {
        nNode_owned_sum[ii] += nNode_owned_sum[ii-1];
      }
      //printVector(nNode_owned_sum);

      bfs_start = 0;
      if(this_mpi_proc > 0)
        bfs_start = nNode_owned_sum[this_mpi_proc-1];
      bfs_end = nNode_owned_sum[this_mpi_proc]-1;

      PetscSynchronizedPrintf(MPI_COMM_WORLD, "\n    bfs_start = %d \t bfs_end = %d \n", bfs_start, bfs_end);
      PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT);

      row_start = bfs_start*ndof;
      row_end   = (bfs_end+1)*ndof-1;

      bfs_local   = bfs_end - bfs_start + 1;
      ndofs_local = row_end - row_start + 1;

      //cout << "row_start =" << row_start << '\t' << row_end << endl;
      PetscSynchronizedPrintf(MPI_COMM_WORLD, "\n    row_start = %d \t row_end = %d \n", row_start, row_end);
      PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT);

      MPI_Barrier(MPI_COMM_WORLD);

      std::vector<int>  displs(n_mpi_procs);

      displs[0] = 0;
      for(ii=0; ii<n_mpi_procs-1; ii++)
        displs[ii+1] = displs[ii] + nNode_owned_vector[ii];

      node_map_get_old.resize(nNode, -1);
      node_map_get_new.resize(nNode, -1);

      // create a global list of nodelist_owned
      // which will serve as a mapping from NEW node numbers to OLD node numbers
      ierr = MPI_Allgatherv(&nodelist_owned[0], nNode_local, MPI_INT, &node_map_get_old[0], &nNode_owned_vector[0], &displs[0], MPI_INT, MPI_COMM_WORLD);


      // create an array for mapping from OLD node numbers to NEW node numbers
      // Also, generate NodeTypeNew array for computing the local and global DOF size
      // as well as creating the element-wise array for element matrix/vector assembly
      for(ii=0; ii<nNode; ii++)
        node_map_get_new[node_map_get_old[ii]] = ii;


      ierr = PetscPrintf(MPI_COMM_WORLD, " Preparing DOF maps \n");

      tempDOF = nNode*ndof;

      dof_map_get_old.resize(tempDOF, 0);

      kk = 0;
      for(ii=0; ii<nNode; ii++)
      {
        ind1 = node_map_get_old[ii]*ndof;
        ind2 = node_map_get_new[ii]*ndof;

        for(jj=0; jj<ndof; jj++)
        {
          dof_map_get_old[kk] = ind1 + jj;
          kk++;
        }
      }

      ierr  = PetscFree(elem_proc_id); CHKERRQ(ierr);
      ierr  = PetscFree(node_proc_id); CHKERRQ(ierr);
      
      /////////////////////////////////////////////////////////////
      // 
      // reorder element node numbers, and
      // the the associated DOFs
      /////////////////////////////////////////////////////////////

      //for(e1=0; e1<nElem; e1++)
      //{
        //ee = fluidElementIds[e1];
        ////if(elems[ee]->getSubdomainId() == this_mpi_proc)
        ////{
          //for(ii=0; ii<elems[ee]->GlobalBasisFuncs.size(); ii++)
          //{
            //elems[ee]->GlobalBasisFuncs[ii] = node_map_get_new[grid_to_cutfem_BF[elems[ee]->GlobalBasisFuncs[ii]]];
          //}
        ////}
      //}
    } // if(n_mpi_procs > 1)

    /////////////////////////////////////////////////////
    // mesh partitioning and node reordering is done
    /////////////////////////////////////////////////////

    MPI_Barrier(MPI_COMM_WORLD);

    for(ee=0; ee<activeElements.size(); ee++)
    {
      elems[activeElements[ee]]->initialiseDOFvalues();
    }

    PetscPrintf(MPI_COMM_WORLD, "\n element DOF values initialised \n\n");
    //printf("\n Finding Global positions \n\n");

    grid_to_proc_BF.assign(gridBF1, -1);

    for(ii=0; ii<gridBF1; ii++)
    {
      //cout << ii << '\t' << grid_to_cutfem_BF[ii] << endl;
      if(grid_to_cutfem_BF[ii] != -1)
      {
        grid_to_proc_BF[ii]   =  node_map_get_new[grid_to_cutfem_BF[ii]];

        ind1 = ii*ndof;
        ind2 = grid_to_proc_BF[ii]*ndof;

        for(jj=0; jj<ndof; jj++)
          grid_to_proc_DOF[ind1+jj]  =  ind2+jj;
      }
    }

    for(ii=0; ii<totalDOF; ii++)
      proc_to_grid_DOF[ii]  =  cutfem_to_grid_DOF[dof_map_get_old[ii]];

    MPI_Barrier(MPI_COMM_WORLD);

    PetscSynchronizedPrintf(MPI_COMM_WORLD, "\n    preparing matrix pattern    \n");
    PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT);

    //vector<vector<int> > DDconn;
    vector<set<int> > DDconnLoc;
    set<int>::iterator it;

    DDconnLoc.resize(nNode);

    jj=200;
    if(totalDOF < 200)      jj= totalDOF;
    //for(ii=row_start; ii<=row_end; ii++)
      //DDconn[ii].reserve(jj);

    if(this_mpi_proc >= 90)
    {
        cout << "this_mpi_proc = " << this_mpi_proc << endl;
        cout << "grid_to_proc_BF : " << endl;
        printVector(grid_to_proc_BF);
        cout << endl;cout << endl;cout << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(this_mpi_proc == 90)
    {
        cout << "this_mpi_proc = " << this_mpi_proc << endl;
        cout << "grid_to_proc_DOF : " << endl;
        printVector(grid_to_proc_DOF);
        cout << endl;cout << endl;cout << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    vector<int>  vecIntTemp;

    for(ee=0; ee<nElem; ee++)
    {
        nd1 = elems[fluidElementIds[ee]];

        val1 =  nd1->getNumberOfBasisFunctions();
        vecIntTemp  =  nd1->GlobalBasisFuncs;

        //cout << " aa = " << ee << '\t' << this_mpi_proc << '\t' << nd1->domNums[0] << endl;
        //PetscSynchronizedPrintf(MPI_COMM_WORLD, " aa = %d \t %d \t %d \t %d \n", ee, fluidElementIds[ee], this_mpi_proc, nd1->domNums[0]);
        //if(this_mpi_proc==1) printVector(vecIntTemp);

        if( nd1->domNums[0] == 0 )
        {
          for(ii=0; ii<val1; ii++)
          {
            r = grid_to_proc_BF[vecIntTemp[ii]];

            // this saves a lot of memory
            // and prevents crashing due to SIGNAL 9
            if( r >= bfs_start && r <= bfs_end )
            {
              for(jj=0; jj<val1; jj++)
                DDconnLoc[r].insert(grid_to_proc_BF[vecIntTemp[jj]]);
            }
          }
        }

        //PetscSynchronizedPrintf(MPI_COMM_WORLD, " bb = %d \t %d \n", ee, fluidElementIds[ee]);
        //cout << " bb = " << ee << '\t' << this_mpi_proc << endl;

        // connectivity for ghost-penalty terms
        if( nd1->isCutElement() )
        {
          for(side=0; side<NUM_NEIGHBOURS; side++)
          {
            nd2 = nd1->getNeighbour(side);

            if( (nd2 != NULL) && !(nd2->isGhost()) && nd2->isLeaf() && (nd2->isCutElement() || nd2->domNums[0] == 0) )
            {
              nr1 = nd1->GlobalBasisFuncs.size();
              nr2 = nd2->GlobalBasisFuncs.size();

              for(ii=0; ii<nr1; ii++)
              {
                r = grid_to_proc_BF[nd1->GlobalBasisFuncs[ii]];

                if(r != -1)
                {
                  for(jj=0; jj<nr2; jj++)
                  {
                    c = grid_to_proc_BF[nd2->GlobalBasisFuncs[jj]];

                    if(c != -1)
                    {
                      DDconnLoc[r].insert(c);
                      DDconnLoc[c].insert(r);
                    }
                  }
                } // for(jj=0; jj<nr2; jj++)
              } // for(ii=0; ii<nr1; ii++)
            }
          } //for(side=0; side<NUM_NEIGHBOURS; side++)
        } //if( nd1->isCutElement() )

        //cout << " cc = " << ee << '\t' << this_mpi_proc << endl;
        //PetscSynchronizedPrintf(MPI_COMM_WORLD, " cc = %d \t %d \n", ee, fluidElementIds[ee]);
        //PetscSynchronizedFlush(MPI_COMM_WORLD, PETSC_STDOUT);
      //} //if( nd1->getSubdomainId() == this_mpi_proc )
    } // for(e=0;e<fluidElementIds.size();e++)

    MPI_Barrier(MPI_COMM_WORLD);

    tend = MPI_Wtime();

    PetscPrintf(MPI_COMM_WORLD, " Finding global positions ... FINISHED. Took %f  milliseconds \n", (tend-tstart)*1000);
    MPI_Barrier(MPI_COMM_WORLD);

    // Setup Petsc matrix and vectors
    //
    //////////////////////////////////////////////
    tstart = MPI_Wtime();

    SolnData.node_map_get_old = node_map_get_old;
    SolnData.node_map_get_new = node_map_get_new;

    GeomData.node_map_get_old = node_map_get_old;
    GeomData.node_map_get_new = node_map_get_new;

    ////////////////////////////////////////////////////////////////////////////////////////////
    // find the number of nonzeroes in the diagonal and off-diagonal portions of each processor
    ////////////////////////////////////////////////////////////////////////////////////////////

    MPI_Barrier(MPI_COMM_WORLD);

    PetscInt  count_diag=0, count_offdiag=0, tempInt;

    //PetscInt  d_nnz[ndofs_local], o_nnz[ndofs_local];
    // this leads to segmentation fault errors. So, allocate using Malloc functions. And, don't forget to free this memory.

    PetscInt  *d_nnz, *o_nnz;

    ierr  = PetscMalloc1(ndofs_local,  &d_nnz);CHKERRQ(ierr);
    ierr  = PetscMalloc1(ndofs_local,  &o_nnz);CHKERRQ(ierr);

    //cout << " bfs sizes " << this_mpi_proc << '\t' << bfs_start << '\t' << bfs_end << '\t' << bfs_local << endl;
    //cout << " dof sizes " << this_mpi_proc << '\t' << row_start << '\t' << row_end << '\t' << ndofs_local << endl;
    PetscPrintf(MPI_COMM_WORLD, " rank = %d  \t row_start = %d  \t row_end = %d  \t ndofs_local = %d  \n", this_mpi_proc,  row_start,  row_end,  ndofs_local);


    kk=0;
    nnz_max_row = 0;
    for(ii=bfs_start; ii<=bfs_end; ii++)
    {
      size1 = DDconnLoc[ii].size()*ndof;

      nnz_max_row = max(nnz_max_row, size1);

      count_diag=0, count_offdiag=0;
      for(it=DDconnLoc[ii].begin(); it!=DDconnLoc[ii].end(); ++it)
      {
        tempInt = *it;

        if(tempInt >= bfs_start && tempInt <= bfs_end)
          count_diag++;
        else
          count_offdiag++;
      }

      //cout << " count_diag ..." << ii << '\t' << count_diag << '\t' << count_offdiag << endl;

      ind1          = ndof*kk++;
      count_diag    = ndof*count_diag;
      count_offdiag = ndof*count_offdiag;
      for(jj=0; jj<ndof; jj++)
      {
        r=ind1+jj;
        d_nnz[r] = count_diag;
        o_nnz[r] = count_offdiag;
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    //PetscPrintf(MPI_COMM_WORLD, "count_diag... %d \t %d \t %d \n", count_diag, count_offdiag, nnz_max_row );
    //cout << "count_diag..." << this_mpi_proc << '\t' << count_diag << '\t' << count_offdiag << '\t' << nnz_max_row << endl;
    //Create parallel matrix, specifying only its global dimensions.
    //When using MatCreate(), the matrix format can be specified at
    //runtime. Also, the parallel partitioning of the matrix is
    //determined by Petsc at runtime.
    //Performance tuning note: For problems of substantial size,
    //preallocation of matrix memory is crucial for attaining good
    //performance. See the matrix chapter of the users manual for details.

    /*
    ierr = MatCreate(MPI_COMM_WORLD, &solverPetsc->mtx);CHKERRQ(ierr);
    ierr = MatSetSizes(solverPetsc->mtx, ndofs_local, ndofs_local, totalDOF, totalDOF);CHKERRQ(ierr);

    ierr = MatSetFromOptions(solverPetsc->mtx);CHKERRQ(ierr);

    ierr = MatMPIAIJSetPreallocation(solverPetsc->mtx, 20, d_nnz, 20, o_nnz);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(solverPetsc->mtx, 20, d_nnz);CHKERRQ(ierr);
    */

    ierr = MatCreateAIJ(MPI_COMM_WORLD, ndofs_local, ndofs_local, totalDOF, totalDOF, 50, d_nnz, 50, o_nnz, &(solverPetsc->mtx));
    MPI_Barrier(MPI_COMM_WORLD);

    //colTemp and arrayTemp are defined in HBSplineBase.h
    // and initially allocated of size 5000
    if(nnz_max_row > 5000)
    {
      ierr  = PetscMalloc1(nnz_max_row,  &colTemp);CHKERRQ(ierr);
      ierr  = PetscMalloc1(nnz_max_row,  &arrayTemp);CHKERRQ(ierr);
    }

    for(jj=0; jj<nnz_max_row; jj++)
      arrayTemp[jj] = 0.0;

    for(ii=bfs_start; ii<=bfs_end; ii++)
    {
      kk=0;
      for(it=DDconnLoc[ii].begin(); it!=DDconnLoc[ii].end(); ++it)
      {
        ind = (*it)*ndof;
        for(jj=0; jj<ndof; jj++)
        {
          colTemp[kk++] = ind+jj;
        }
      }

      ind1 = ndof*ii;
      ind2 = ndof*DDconnLoc[ii].size();

      for(jj=0; jj<ndof; jj++)
      {
        r = ind1+jj;
        ierr = MatSetValues(solverPetsc->mtx, 1, &r, ind2, colTemp, arrayTemp, INSERT_VALUES);
      }
    }

    for(ii=0; ii<DDconnLoc.size(); ii++)
      DDconnLoc[ii].clear();
    DDconnLoc.clear();

    MPI_Barrier(MPI_COMM_WORLD);

    ierr = MatSetOption(solverPetsc->mtx, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE);

    //VecCreate(MPI_COMM_WORLD, &solverPetsc->soln);
    ierr = VecCreateMPI(MPI_COMM_WORLD, ndofs_local, totalDOF, &solverPetsc->soln); CHKERRQ(ierr);
    ierr = VecSetOption(solverPetsc->soln, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
    ierr = VecSetFromOptions(solverPetsc->soln);CHKERRQ(ierr);
    ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->rhsVec);CHKERRQ(ierr);
    ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->solnPrev);CHKERRQ(ierr);
    ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->reac);CHKERRQ(ierr);

    solverPetsc->currentStatus = PATTERN_OK;

    soln.resize(totalDOF);
    soln.setZero();
    solnInit = soln;

    GRID_CHANGED = IB_MOVED = false;

    tend = MPI_Wtime();

    PetscPrintf(MPI_COMM_WORLD, " HBSplineCutFEM::prepareMatrixPattern()  .... FINISHED. Took %f  milliseconds \n", (tend-tstart)*1000);

    MPI_Barrier(MPI_COMM_WORLD);

    // compute quadrature points for cut cells
    //
    //////////////////////////////////////////////
    tstart = MPI_Wtime();


    for(ee=0; ee<fluidElementIds.size(); ee++)
    {
      nd1 = elems[fluidElementIds[ee]];

      if( nd1->getSubdomainId() == this_mpi_proc )
      {
        //nd1->clearSubtriangulation();

        if( nd1->isCutElement() )
        {
          if(cutFEMparams[0] == 1)
            nd1->computeGaussPointsSubTrias(cutFEMparams[2], false);
          else
            nd1->computeGaussPointsAdapIntegration(cutFEMparams[3], cutFEMparams[4], false, true);
        }
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    ierr  = PetscFree(d_nnz);   CHKERRQ(ierr);
    ierr  = PetscFree(o_nnz);   CHKERRQ(ierr);

    tend = MPI_Wtime();

    PetscPrintf(MPI_COMM_WORLD, " computing cutcell gps took %f  milliseconds \n", (tend-tstart)*1000);

    return 1;
}






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



