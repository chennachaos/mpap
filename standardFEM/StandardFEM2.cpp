
#include "SolverPardisoEigen.h"
#include "SolverMA41Eigen.h"
#include "StandardFEM.h"
#include "MpapTime.h"
#include "TimeFunction.h"
#include "SolverPetsc.h"
#include "LagrangeElement.h"
#include "SolutionData.h"
#include "conditionalOStream.h"
#include "my_types.h"
#include "mpi.h"
#include "metis.h"

extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;

typedef int idx_t;


void StandardFEM::setSolver(int slv, int *parm, bool cIO)
{
    //if(solver != NULL)
      //delete solver;
    //solver = NULL;

    Eigen::initParallel();
    
    int numProc=1;

    switch(slv)
    {
        case  1: // MA41 ..........................

            SOLVER_TYPE = SOLVER_TYPE_EIGEN;

            solverEigen = (SolverEigen*) new SolverMA41Eigen;

            prepareMatrixPattern();

            solverEigen->printInfo();

            if(solverEigen->initialise(0,0,totalDOF) != 0)
              return;

        break;

        case  4: // SolverEigen ..........................

            SOLVER_TYPE = SOLVER_TYPE_EIGEN;

            solverEigen = (SolverEigen*) new SolverEigen;

            //printInfo();
            solverEigen->setAlgorithmType(parm[0]);

            prepareMatrixPattern();

            if(solverEigen->initialise(0,0,totalDOF) != 0)
              return;
            //solver->setSolverAndParameters();

            solverEigen->printInfo();

        break;

        case  5: // PARDISO(sym) with Eigen
        case  6: // PARDISO(unsym) with Eigen

            SOLVER_TYPE = SOLVER_TYPE_EIGEN;

            solverEigen = (SolverEigen*) new SolverPardisoEigen;

            if (parm != NULL) numProc = parm[0]; else numProc = 1;

            numProc = min(MAX_PROCESSORS,numProc);

            //printInfo();
            prepareMatrixPattern();

            if(slv == 5)
            {
              if (solverEigen->initialise(numProc, PARDISO_STRUCT_SYM, totalDOF) != 0)
                return;
            }
            if(slv == 6)
            {
              if (solverEigen->initialise(numProc, PARDISO_UNSYM, totalDOF) != 0)
                return;
            }

            //solver->printInfo();

        break;

        case  8: // SolverPetsc ..........................

            SOLVER_TYPE = SOLVER_TYPE_PETSC;

            solverPetsc = (SolverPetsc*) new SolverPetsc;

            //printInfo();
            //solver2->setAlgorithmType(parm[0]);

            prepareMatrixPattern();

            if(solverPetsc->initialise(nNode*ndof, 0, totalDOF) != 0)
              return;

            //solverPetsc->setSolverAndParameters();
            solverPetsc->printInfo();

        break;

        default: // invalid slv ...................

             cout << " this solver has not been implemented yet!\n\n";

        break;
    }

    solverOK = true;

    if(solverEigen != NULL)
      solverEigen->checkIO = cIO;

    if(solverPetsc != NULL)
      solverPetsc->checkIO = cIO;

    //if( (tis > 0) )
      //setInitialConditions();

    //cout << " cIO " << cIO << endl;

    return;
}




int StandardFEM::prepareMatrixPattern()
{
    printf("\n     StandardFEM::prepareMatrixPattern()  .... STARTED ...\n");

    char fct[] = "StandardFEM::prepareMatrixPattern";

    int  r, c, r1, c1, count=0, count1=0, count2=0, iii, e, ind, nsize;
    int *tt, *tt1, *tt2,  val1, val2, n1, n2, a, b, ll, pp, nnz;
    int  nRow, nCol, ind1, ind2;
    int  ee, ii, jj, kk, e1, e2;

    ConditionalOStream  pcout(std::cout,  (this_mpi_proc == 0) );

    int n_subdomains = n_mpi_procs, subdomain=0;


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

      //cout << " 1111111111111111 "  << nElem_local << '\t' << npElem << endl;

      n2 = ceil(nElem/n_mpi_procs);

      elem_start = n2*this_mpi_proc;
      elem_end   = n2*(this_mpi_proc+1)-1;

      if(this_mpi_proc == (n_mpi_procs-1))
        elem_end = nElem-1;

      nElem_local = elem_end - elem_start + 1;

      //cout << " elem_start = " << elem_start << '\t' << elem_end << '\t' << nElem_local << endl;

      PetscInt  *eptr, *eind;

      //PetscInt  eptr[nElem_local+1];
      //PetscInt  eind[nElem_local*npElem];

      //ierr  = PetscMalloc1(nElem_local+1,  &eptr);CHKERRQ(ierr);
      //ierr  = PetscMalloc1(nElem_local*npElem,  &eind);CHKERRQ(ierr);

      ierr  = PetscMalloc1(nElem+1,  &eptr);CHKERRQ(ierr);
      ierr  = PetscMalloc1(nElem*npElem,  &eind);CHKERRQ(ierr);

      //PetscInt  eptr[nElem+1];
      //PetscInt  eind[nElem*npElem];

      //idx_t  eptr[nElem+1];
      //idx_t  eind[nElem*npElem];

      //PetscSynchronizedPrintf(PETSC_COMM_WORLD,"%d \n", elem_start);

      vector<int>  vecTemp2;

      ind = 0;
      //for(e1=0; e1<nElem_local; e1++)
      for(e1=0; e1<nElem; e1++)
      {
        eptr[e1] = e1*npElem;

        //e2 = elem_start + e1;
        e2 = e1;
        kk = e1*npElem;

        vecTemp2 = elems[e2]->nodeNums ;

        //findUnique(vecTemp2);

        for(ii=0; ii<npElem; ii++)
          eind[kk+ii] = vecTemp2[ii] ;

        ind++;
      }

      //eptr[nElem_local] = nElem_local*npElem;
      eptr[nElem] = nElem*npElem;

      //cout << " \n\n\n\n " << endl;

      //PetscInt  nodes_per_side;
      idx_t nodes_per_side;

      if(DIM == 2)
        nodes_per_side = 2;
      else
        nodes_per_side = 3;



    //
    idx_t nWeights  = 1;

    idx_t objval, numflag=0;
    idx_t *xadj, *adjncy;

    idx_t epart[nElem];
    idx_t npart[nNode];
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
      //int ret = METIS_PartMeshNodal(&nElem, &nNode, eptr, eind, NULL, NULL, &n_mpi_procs, NULL, options, &objval, epart, npart);
      int ret = METIS_PartMeshDual(&nElem, &nNode, eptr, eind, NULL, NULL, &nodes_per_side, &n_mpi_procs, NULL, options, &objval, epart, npart);

      //idx_t  wgtflag=0, numflag=0, ncon=0, ncommonnodes=2, nparts=n_mpi_procs;
      //idx_t  *elmdist;

      //int ret = ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, NULL, &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts, NULL, NULL, options, NULL, npart, MPI_COMM_WORLD);

      //pcout <<  " ccccccccccc " << endl;

      if(ret == METIS_OK)
        std::cout << " METIS partition routine successful "  << std::endl;
      else
        std::cout << " METIS partition routine FAILED "  << std::endl;

      //if(this_mpi_proc == 0)
      //{
        //for(ee=0; ee<nNode; ee++)
          //cout << ee << '\t' << npart[ee] << endl;
      //}

      for(ee=0; ee<nElem; ee++)
      {
        elems[ee]->setSubdomainId(epart[ee]);
      }

      ierr = PetscFree(eptr); CHKERRQ(ierr);
      ierr = PetscFree(eind); CHKERRQ(ierr);


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

      for(ee=0; ee<nElem; ee++)
      {
        //if(elems[ee]->getSubdomainId() == this_mpi_proc)
        //{
          for(ii=0; ii<elems[ee]->nodeNums.size(); ii++)
          {
            elems[ee]->nodeNums[ii] = node_map_old_to_new[elems[ee]->nodeNums[ii]];
          }
        //}
      }

      for(ii=0;ii<DirichletBCs.size();ii++)
      {
        DirichletBCs[ii][0] = node_map_old_to_new[DirichletBCs[ii][0]] ;
      }
  } // if(n_mpi_procs > 1)


    pcout << " preparing matrix pattern " << endl;

    /////////////////////////////////////////////////////////////
    //
    // prepare the matrix pattern
    /////////////////////////////////////////////////////////////

    IEN.resize(nElem);

    for(ee=0;ee<nElem;ee++)
    {
      for(ii=0;ii<npElem;ii++)
        IEN[ee].push_back(elems[ee]->nodeNums[ii]);
    }

    NodeType.resize(nNode);
    ID.resize(nNode);
    for(ii=0;ii<nNode;ii++)
    {
      NodeType[ii].resize(ndof);
      ID[ii].resize(ndof);
      for(jj=0;jj<ndof;jj++)
      {
        NodeType[ii][jj] = false;
        ID[ii][jj] = -1;
      }
    }

    for(ii=0;ii<DirichletBCs.size();ii++)
    {
      NodeType[DirichletBCs[ii][0]][DirichletBCs[ii][1]] = true;
    }

    int nn, dof;
    for(ii=0;ii<DirichletBCs.size();ii++)
    {
      nn  = DirichletBCs[ii][0];
      dof = DirichletBCs[ii][1];
      SolnData.var1applied[nn*ndof+dof] =  DirichletBCs[ii][2];
    }

    //cout << " BBBBBBBBBBBBBBBBB " << endl;

    soln.resize(nNode*ndof);
    soln.setZero();

    solnInit = soln;
    reac  =  soln;


    totalDOF = 0;
    for(ii=0;ii<nNode;ii++)
    {
      for(jj=0;jj<ndof;jj++)
      {
        //if(!NodeType[ii][jj])
          ID[ii][jj] = totalDOF++;
      }
    }

    cout << " totalDOF " << totalDOF << endl;
    cout << " nElem    " << nElem << endl;
    cout << " nNode    " << nNode << endl;
    cout << " npElem   " << npElem << endl;
    cout << " ndof     " << ndof << endl;

    LM.resize(nElem);

    for(ee=0;ee<nElem;ee++)
    {
      npElem = IEN[ee].size();

      //printVector(IEN[ee]);

      ind = ndof*npElem;
      LM[ee].resize(ind);

      for(ii=0;ii<npElem;ii++)
      {
        ind = ndof*ii;

        //itint = find(GlobalPointNumbers.begin(), GlobalPointNumbers.end(), IEN[ee][ii]);
        //kk = distance(GlobalPointNumbers.begin(), itint);
        //IEN2[ee][ii] = kk;
        kk = IEN[ee][ii];

        for(jj=0;jj<ndof;jj++)
        {
          //cout << ee << '\t' << ii << '\t' << jj << '\t' << ind << '\t' << ID[kk][jj] << endl;
          //cout << ee << '\t' << ii << '\t' << jj << '\t' << ind << '\t' << LM[ee][ind+jj] << '\t' << ID[IEN[ee][ii]][jj] << endl;
          LM[ee][ind+jj] = ID[kk][jj];
          //cout << " IIIIIIIII " << endl;
        }
      }
    }

    assy4r.resize(totalDOF);

    count = 0;
    for(ii=0;ii<nNode;ii++)
    {
      for(jj=0;jj<ndof;jj++)
      {
        //cout << ii << '\t' << jj << '\t' << ID[ii][jj] << endl;
        if( ID[ii][jj] != -1)
          assy4r[count++] = ii*ndof + jj;
      }
    }
    //cout << " totalDOF " << totalDOF << endl;
    for(ee=0;ee<nElem;ee++)
    {
      //elems[ee]->nodeNums   = IEN2[ee];
      //elems[ee]->nodeNums   = IEN[ee];
      elems[ee]->forAssyVec = LM[ee];
      elems[ee]->prepareElemData();
    }

    pp=false;
    //pp=true;
    if(pp)
    {
       printf("   IEN array \n\n");
       for(ii=0;ii<nElem;ii++)
       {
          for(jj=0;jj<npElem;jj++)
            cout << '\t' << IEN[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("   ID array \n\n");
       for(ii=0;ii<nNode;ii++)
       {
          for(jj=0;jj<ndof;jj++)
            cout << '\t' << ID[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("   LM array \n\n");
       for(ii=0;ii<nElem;ii++)
       {
          for(jj=0;jj<LM[ii].size();jj++)
            cout << '\t' << LM[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("  assy4r array \n\n");
       for(ii=0;ii<totalDOF;ii++)
       {
          cout << assy4r[ii] << endl;
       }
       printf("\n\n\n");
    }

    printf("\n element DOF values initialised \n\n");
    printf("\n Preparing matrix pattern \n\n");

    vector<int>::const_iterator location;
    //set<int>::iterator it;

    forAssyMat.resize(totalDOF);
    
    for(ee=0;ee<nElem;ee++)
    {
      tt = &(LM[ee][0]);
      nsize = LM[ee].size();

      for(ii=0;ii<nsize;ii++)
      {
        count1 = tt[ii];

        if(tt[ii] != -1)
        {
          for(jj=0;jj<nsize;jj++)
          {
            if(tt[jj] != -1)
              forAssyMat[count1].push_back(tt[jj]);
          }
        }
      }
    }

    printf("\n Preparing matrix pattern DONE \n\n");

    VectorXi  nnzVec(totalDOF);

    nnz = 0;
    for(ii=0;ii<totalDOF;ii++)
    {
      findUnique(forAssyMat[ii]);
      nnzVec[ii] = forAssyMat[ii].size();
      nnz += nnzVec[ii];
    }
    cout << " nnz " << nnz << endl;

    nRow = nCol = totalDOF;

    bool pp1=false;
    //pp1=true;
    if(pp1)
    {
       printf("   Number of non-zeros = %5d \n\n", nnz);
       printf("   dof to dof connectivity ...:  \n\n");
       for(ii=0;ii<totalDOF;ii++)
       {
          cout << " dof # " << ii << " : ";
          for(jj=0;jj<forAssyMat[ii].size();jj++)
            cout << '\t' << forAssyMat[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");
    }

    //cout << " AAAAAAAAAA " << endl;

    SolnData.node_map_new_to_old = node_map_new_to_old;
    SolnData.node_map_old_to_new = node_map_old_to_new;

    GeomData.node_map_new_to_old = node_map_new_to_old;
    GeomData.node_map_old_to_new = node_map_old_to_new;


    if(SOLVER_TYPE == SOLVER_TYPE_PETSC)
    {
      /////////////////////////////////////////
      //
      // Petsc based solver
      //
      /////////////////////////////////////////

      cout << " PetscSolver " << endl;

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
        for(jj=0; jj<forAssyMat[ii].size(); jj++)
        {
          tempInt = forAssyMat[ii][jj];

          if(tempInt >= row_start && tempInt <= row_end)
            count_diag++;
          else
            count_offdiag++;
        }
        d_nnz[kk] = count_diag;
        o_nnz[kk] = count_offdiag;
        kk++;
      }

      //MatCreateSeqAIJ(PETSC_COMM_WORLD, nRow, nCol, 500, &nnzVec(0), &(solver2->mtx));
      //PetscErrorCode ierr;

      //
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
        for(jj=0; jj<forAssyMat[ii].size(); jj++)
        {
          ierr = MatSetValue(solverPetsc->mtx, ii, forAssyMat[ii][jj], 0.0, INSERT_VALUES);
        }
      }

      VecCreate(PETSC_COMM_WORLD, &solverPetsc->soln);
      VecCreate(PETSC_COMM_WORLD, &solverPetsc->solnPrev);
      VecCreate(PETSC_COMM_WORLD, &solverPetsc->rhsVec);
      VecCreate(PETSC_COMM_WORLD, &solverPetsc->reac);

      ierr = VecSetSizes(solverPetsc->soln, ndofs_local, totalDOF); CHKERRQ(ierr);
      ierr = VecSetSizes(solverPetsc->reac, ndofs_local, totalDOF); CHKERRQ(ierr);

      ierr = VecSetFromOptions(solverPetsc->soln);CHKERRQ(ierr);
      ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->rhsVec);CHKERRQ(ierr);
      ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->solnPrev);CHKERRQ(ierr);
      ierr = VecSetFromOptions(solverPetsc->reac);CHKERRQ(ierr);

      solverPetsc->currentStatus = PATTERN_OK;
    }
    else
    {
      /////////////////////////////////////////
      //
      // Eigen based solver
      //
      /////////////////////////////////////////

      cout << " Eigen based solver " << totalDOF << endl;

      //solver->rhsVec.resize(nRow);

      solverEigen->mtx.setZero();

      solverEigen->mtx.resize(nRow, nCol);
      solverEigen->mtx.reserve(nnz);
      solverEigen->mtx.reserve(nnzVec);

      for(ii=0;ii<totalDOF;ii++)
      {
        for(jj=0;jj<forAssyMat[ii].size();jj++)
        {
          solverEigen->mtx.coeffRef(ii, forAssyMat[ii][jj]) = 0.0;
        }
      }

      solverEigen->mtx.makeCompressed();

      solverEigen->currentStatus = PATTERN_OK;
    }

    //SolnData.var1 += SolnData.var1applied;

    // remove data objects

    nodePosData.clear();
    elemConn.clear();
    IEN.clear();
    LM.clear();
    ID.clear();
    forAssyMat.clear();

    printf("\n     StandardFEM::prepareMatrixPattern()  .... FINISHED ...\n\n");

    return 1;
}




int StandardFEM::solveStep(int niter)
{
    for(int iter=0; iter<niter; iter++)
    {
      calcStiffnessAndResidual();

      if( converged() )
        break;

      factoriseSolveAndUpdate();

      updateIterStep();
    }
    
    if( !converged() )
      return 0;
  
  return 1;
}





int StandardFEM::calcStiffnessAndResidual(int printRes, bool zeroMtx, bool zeroRes)
{
    //cout << "     StandardFEM: generating coefficient Matrices ...\n\n";

    char fct[] = "StandardFEM::calcStiffnessAndResidual";

    //computerTime.go(fct);

    //if(solver == NULL || solver2 == NULL)
    //{
      //COUT << "You need to select a solver first!\n\n";
      //return 1;
    //}
  
    //cout << " firstIter = " << firstIter << endl;

    if(firstIter)
    {
      rNorm = -1.0;
    }

    //printVector(SolnData.var1applied);

    MatrixXd  Klocal;
    VectorXd  Flocal;

    SolnData.reac.setZero();

    if(SOLVER_TYPE == SOLVER_TYPE_PETSC)
    {
      //cout << " uuuuuuuuuuuuuu " << endl;
      solverPetsc->zeroMtx();
      //cout << " uuuuuuuuuuuuuu " << endl;
      for(int ee=0;ee<nElem;ee++)  // loop over all the elements
      {
        //cout << "       elem... : " << (ee+1) << endl;
        if(elems[ee]->getSubdomainId() == this_mpi_proc)
        {
          //elems[ee]->resetMatrixAndVector();
          //cout << " MMMMMMMMMMM " << endl;
          elems[ee]->calcStiffnessAndResidual(Klocal, Flocal);

          //elems[ee]->assembleElementMatrix(0, solver2->mtx);
          //cout << " ooooooooooooooo " << endl;
          //elems[ee]->assembleElementVector(0, 0, solverPetsc->rhsVec, solver2->reac, 0, 0);
          //cout << " MMMMMMMMMMM " << endl;
          solverPetsc->assembleMatrixAndVector(elems[ee]->forAssyVec, elems[ee]->forAssyVec, Klocal, Flocal);
        }
      }
    }
    else
    {
      solverEigen->zeroMtx();

      //printVector(solver->rhsVec);

      for(int ee=0;ee<nElem;ee++)  // loop over all the elements
      {
        //cout << "       elem... : " << (ee+1) << endl;

        //elems[ee]->resetMatrixAndVector();
        //cout << " MMMMMMMMMMM " << endl;
        elems[ee]->calcStiffnessAndResidual(Klocal, Flocal);

        //cout << " MMMMMMMMMMM " << endl;
        //elems[ee]->assembleElementMatrixAndVector(0, solver->mtx, &(solver->rhsVec(0)));

        //elems[ee]->assembleElementMatrix(0, solver->mtx);
        //cout << " aaaaaaaaaaaaa " << endl;
        //elems[ee]->assembleElementVector(0, 0, &(solverEigen->rhsVec(0)), &(SolnData.reac(0)), 0, 0);
        //cout << " MMMMMMMMMMM " << endl;
        solverEigen->assembleMatrixAndVector(0, 0, elems[ee]->forAssyVec, Klocal, Flocal);
      }
    }

    //cout << " MMMMMMMMMMM " << endl;

    //cout << solver->mtx << endl;

    //cout << " solver->rhsVec " << endl;        printVector(solver->rhsVec);

    //printf("\n solver->rhsVec norm = %12.6E \n", solver->rhsVec.norm());

    applyBoundaryConditions();

    applyExternalForces();
   
    //cout << solver->mtx << endl;

    // rhs due to external forces
    //rhsVec += rhsVec2;

    //cout << " rhsVec " << endl;        printVector(&(rhsVec[0]), totalDOF);

    //printf("\n rhsVec norm = %12.6E \n", solver->rhsVec.norm());

/*
  if(firstIter)
  {
    //SolnData.var1 = fact1*SolnData.var1applied;

    int ii, jj, nn, dof, aa, ind;

    double fact1=1.0, specVal;

    fact1 = timeFunction[0].prop;

    vector<int>::iterator itint;

    for(aa=0;aa<DirichletBCs.size();aa++)
    {
      nn  = (int) (DirichletBCs[aa][0]);
      dof = (int) (DirichletBCs[aa][1]);
      specVal = DirichletBCs[aa][2];

      //itint = find(GlobalPointNumbers.begin(), GlobalPointNumbers.end(), nn);
      //nn   = distance(GlobalPointNumbers.begin(), itint);

      ind = nn*ndof+dof;

      //SolnData.var1[ind] = fact1*SolnData.var1applied[ind];

      SolnData.var1[ind] = fact1*specVal ;
    }
  }
*/

    firstIter = false;
    rNormPrev = rNorm;

    if(SOLVER_TYPE == SOLVER_TYPE_PETSC)
    {
      VecAssemblyBegin(solverPetsc->rhsVec);
      VecAssemblyEnd(solverPetsc->rhsVec);

      VecNorm(solverPetsc->rhsVec, NORM_2, &rNorm);
      solverPetsc->currentStatus = ASSEMBLY_OK;
    }
    else
    {
      rNorm = solverEigen->rhsVec.norm();
      solverEigen->currentStatus = ASSEMBLY_OK;
    }

    //if(printRes > 1)
       COUT << "StandardFEM"; printf("  %11.4e\n",rNorm);

    //ctimCalcStiffRes += computerTime.stop(fct);

    //computerTime.stopAndPrint(fct);

    return 0;
}



int StandardFEM::factoriseSolveAndUpdate()
{
  //cout << " StandardFEM::factoriseSolveAndUpdate " << endl;
  //computerTime.go(fct);

  time_t tstart, tend;

  soln.setZero();

  if(SOLVER_TYPE == SOLVER_TYPE_PETSC)
  {
    //VecView(solver2->rhsVec, PETSC_VIEWER_STDOUT_WORLD);

    tstart = time(0);

    solverPetsc->factoriseAndSolve();

    tend = time(0);
    printf("HBSplineCutFEM::factoriseSolveAndUpdate() took %8.4f second(s) \n ", difftime(tend, tstart) );

    //VecView(solver2->soln, PETSC_VIEWER_STDOUT_WORLD);

    /////////////////////////////////////////////////////////////////////////////
    // get the solution vector onto all the processors
    /////////////////////////////////////////////////////////////////////////////
  
    Vec            vec_SEQ;
    VecScatter     ctx;
    PetscScalar *arrayTemp;

    VecScatterCreateToAll(solverPetsc->soln, &ctx, &vec_SEQ);
    VecScatterBegin(ctx, solverPetsc->soln, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx, solverPetsc->soln, vec_SEQ, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterDestroy(&ctx);

    VecGetArray(vec_SEQ, &arrayTemp);

    // update solution vector

    int ii, jj, dd, kk=0;
    for(ii=0; ii<nNode; ii++)
    {
      jj = ii*ndof;
      for(dd=0; dd<ndof; dd++)
        soln[kk++] = arrayTemp[jj+dd];
    }

    VecRestoreArray(vec_SEQ, &arrayTemp);
  }
  else
  {
    //cout << " residue_new " << endl;
    //printVector(&(solver->rhsVec[0]), totalDOF);

    //cout << solver->mtx << endl;

    cout << " Eigen Solver " << endl;
    solverEigen->factoriseAndSolve();

    //printVector(solver->soln); printf("\n\n\n");

    // update solution vector
    for(int kk=0;kk<totalDOF;kk++)
      soln[assy4r[kk]] = solverEigen->soln[kk];
  }

  //if( ElemProp[elemConn[0][1]].id == 7 )
    //SolnData.var1Dot += soln;
  //else
    SolnData.var1 += soln;

    //VectorXd  vecTemp(totalDOF);
    //vecTemp = soln - SolnData.var1;
    //double  dnorm = vecTemp.norm();
    //cout << " dnorm = " << dnorm << endl;
    //SolnData.var1Prev = SolnData.var1;
    //SolnData.var1 = soln;

  //printVector(SolnData.var1);

/*
  //if(firstIter)
  //{
    //SolnData.var1 = fact1*SolnData.var1applied;

    int ii, jj, nn, dof, aa, ind;

    double fact1=1.0, specVal;

    fact1 = timeFunction[0].prop;

    vector<int>::iterator itint;

    for(aa=0;aa<DirichletBCs.size();aa++)
    {
      nn  = (int) (DirichletBCs[aa][0]);
      dof = (int) (DirichletBCs[aa][1]);
      specVal = DirichletBCs[aa][2];

      //itint = find(GlobalPointNumbers.begin(), GlobalPointNumbers.end(), nn);
      //nn   = distance(GlobalPointNumbers.begin(), itint);

      ind = nn*ndof+dof;

      //SolnData.var1[ind] = fact1*SolnData.var1applied[ind];

      SolnData.var1[ind] = fact1*specVal ;
    }
  //}
*/

    //cout << " Max displacement = " << SolnData.var1.maxCoeff() << endl;
    //cout << " Tip displacement = " << SolnData.var1[104*2+1] << endl;

  //cout << " result " << endl;        printVector(SolnData.var1);

  //ctimFactSolvUpdt += computerTime.stop(fct);

  return 0;
}





void StandardFEM::applyBoundaryConditions()
{   
    cout <<  " applying boundary conditions .... " << endl;
    cout << " tis = " << tis << endl;
    if(PHYSICS_TYPE == PHYSICS_TYPE_FLUID)
      cout << " PHYSICS_TYPE = " << PHYSICS_TYPE_FLUID << endl;

    //printVector(SolnData.td);

    int ii, jj, nn, dof, aa, ind;
    double  specVal, PENALTY=1.0e6, val1, val2;

    double  af = SolnData.td(2);

    if(SOLVER_TYPE == SOLVER_TYPE_PETSC)
    {
      //int  row_start, row_end;
      //VecGetOwnershipRange(solver2->rhsVec, &row_start, &row_end);
      //row_end -= 1;

      cout << " row_start = " <<  row_start << '\t' << row_end << '\t' << af << endl;

      for(aa=0;aa<DirichletBCs.size();aa++)
      {
        nn  = (int) (DirichletBCs[aa][0]);
        dof = (int) (DirichletBCs[aa][1]);
        specVal = DirichletBCs[aa][2];

        ind = nn*ndof+dof;
        //specVal = SolnData.var1applied[ind];

        //if( ElemProp[elemConn[0][1]].id == 7 )
          //specVal  -=  SolnData.var1Cur[ind];
        //else
          specVal  -=  SolnData.var1Cur[ind];

        //cout << " values " << '\t' << nn << '\t' << ind << '\t' << specVal << '\t' << SolnData.var1Cur[ind] << endl;

        val1 = PENALTY*af;
        val2 = PENALTY*specVal;

        //ind += start;

        if( (ind >= row_start) && (ind <= row_end) )
        {
          MatSetValue(solverPetsc->mtx, ind, ind, val1, ADD_VALUES);
          VecSetValue(solverPetsc->rhsVec, ind, val2, ADD_VALUES);
        }
      }
    }
    else
    {
      for(aa=0;aa<DirichletBCs.size();aa++)
      {
        nn  = (int) (DirichletBCs[aa][0]);
        dof = (int) (DirichletBCs[aa][1]);
        specVal = DirichletBCs[aa][2];

        //itint = find(GlobalPointNumbers.begin(), GlobalPointNumbers.end(), nn);
        //nn   = distance(GlobalPointNumbers.begin(), itint);

        ind = nn*ndof+dof;
        //specVal = SolnData.var1applied[ind];

        //if( ElemProp[elemConn[0][1]].id == 7 )
          //specVal  -=  SolnData.var1Cur[ind];
        //else
          specVal  -=  SolnData.var1Cur[ind];

        //cout << start << '\t' << nn << '\t' << ind << '\t' << specVal << endl;

        val1 = PENALTY*af;
        val2 = PENALTY*specVal;

        //ind += start;

        solverEigen->mtx.coeffRef(ind, ind) += val1;
        solverEigen->rhsVec[ind]   += val2;
      }
    }

    return;
}



void StandardFEM::applyExternalForces()
{
    int  nn, dof, ii, ind;
    double specVal=0.0, fact=0.0, fact1=0.0;

    VectorXd  vecTemp(SolnData.forceCur.rows());
    vecTemp.setZero();
    for(ii=0;ii<nodeForcesData.size();ii++)
    {
      nn  = (int) (nodeForcesData[ii][0] - 1);
      dof = (int) (nodeForcesData[ii][1] - 1);
      specVal = nodeForcesData[ii][2];

      ind = nn*ndof+dof;

      //cout << nn << '\t' << dof << '\t' << ind << '\t' << specVal << endl;

      vecTemp[ind] += specVal;
    }
    //printVector(vecTemp);

    if(PHYSICS_TYPE == PHYSICS_TYPE_FLUID)
    {
      if(mpapTime.cur <= 5.0e-3)
      {
        //fact = ( 1.0-cos(628.3185*mpapTime.cur));
        //fact = sin(628.3185*mpapTime.cur);
        fact = 1.0;
        //fact = sin(20.0*mpapTime.cur);
      }
      else
        fact = 0.0;

      //cout << " fact = " << fact << endl;
    }

    //fact = 1.0;
    fact = timeFunction[0].prop;
    //fact = mpapTime.cur;
    //fact = sin(20.0*mpapTime.cur);
    //fact = 0.5*( 1.0-cos(628.3185*mpapTime.cur));

    if(SOLVER_TYPE == SOLVER_TYPE_PETSC)
    {
      for(ii=0;ii<totalDOF;ii++)
      {
        VecSetValue(solverPetsc->rhsVec, ii, SolnData.forceCur[assy4r[ii]], ADD_VALUES);

        fact1 = fact*vecTemp[assy4r[ii]];
        VecSetValue(solverPetsc->rhsVec, ii, fact1, ADD_VALUES);
      }
    }
    else
    {
      for(ii=0;ii<totalDOF;ii++)
        solverEigen->rhsVec[ii] += SolnData.forceCur[assy4r[ii]];

      for(ii=0;ii<totalDOF;ii++)
        solverEigen->rhsVec[ii] += (fact*vecTemp[assy4r[ii]]);
    }

    return;
}




void  StandardFEM::computeElementErrors(int index)
{
    int  ii, ee, count=0, dd, domTemp;

    //cout << " index = " << index << endl;

    totalError = 0.0;

    if(index < 4) // L2 or H1 norm based errors
    {
      for(int ee=0;ee<nElem;ee++)  // loop over all the elements
      {
        elems[ee]->calcError(index);

        totalError += elems[ee]->getError();
      }

      totalError = sqrt(totalError);

      if(index < 3)
        printf(" \n\n \t L2 Error = %12.6E \n\n " , totalError);
      else
        printf(" \n\n \t H1 Error = %12.6E \n\n " , totalError);
    }
    else if(index == 10) // total energy
    {
      VectorXd  energyElem(3), energyGlobal(3);
      energyGlobal.setZero();

      for(int ee=0;ee<nElem;ee++)  // loop over all the elements
      {
        elems[ee]->computeEnergy(0, 0, energyElem);

        energyGlobal += energyElem;
      }
      energyGlobal[2] = energyGlobal[0] + energyGlobal[1];

      char        tmp[200];
      MyString    tmpStr;

      sprintf(tmp," \t %12.6E \t %12.6E \t %12.6E", energyGlobal[0], energyGlobal[1], energyGlobal[2]);

      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);
    }
    else if(index == 11) // total momentum
    {
      VectorXd  momElem(6), momGlobal(6);
      momGlobal.setZero();

      for(int ee=0;ee<nElem;ee++)  // loop over all the elements
      {
        elems[ee]->computeMomentum(0, 0, momElem);

        momGlobal += momElem;
      }

      char        tmp[200];
      MyString    tmpStr;

      sprintf(tmp," \t %12.6E \t %12.6E \t %12.6E", momGlobal[0], momGlobal[1], momGlobal[5]);

      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);

    }
    else
    {
      cout << " Compute element errors not defined for thie index value " << endl;
    }
    //printf(" \n\n \t totalError = %12.6E \n\n " , totalError);
    
    return;
}




