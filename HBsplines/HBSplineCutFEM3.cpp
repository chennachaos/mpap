
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

#include <chrono>
typedef std::chrono::high_resolution_clock Clock;



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

  IB_MOVED = false;

  SOLID_SOLVED = 0;

  firstIter = true;
  localStiffnessError = 0;
  iterCount = 1;
  filecount++;

  for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
  {
    ImmersedBodyObjects[bb]->timeUpdate();
  }

  SolnData.timeUpdate();

  if(STAGGERED)
    solveSolidProblem();

  updateIterStep();

  //PetscPrintf(MPI_COMM_WORLD, "   HBSplineCutFEM::timeUpdate() ... FINISHED \n\n");

  return;
}



void HBSplineCutFEM::updateIterStep()
{
  //PetscPrintf(MPI_COMM_WORLD, "   HBSplineCutFEM::updateIterStep() ... STARTED \n\n");

  int kk, bb, ee;

  SolnData.updateIterStep();

  //cout << " kkkkkkkkkkk " << endl;

  for(bb=0;bb<ImmersedBodyObjects.size();bb++)
  {
    computeTotalForce(bb);
    //cout << " bbbbbbbbbbbbbbbb " << endl;
    ImmersedBodyObjects[bb]->updateForce(&(totalForce(0)));
  }

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

  //cout << " IB_MOVED = " << IB_MOVED << endl;
  //cout << " kkkkkkkkkkk " << endl;
  //cout << " PPPPPPPPPPPPPPPPP " << firstIter << '\t' << SOLVER_TYPE << '\t' << GRID_CHANGED << '\t' << IB_MOVED << endl;
  if(GRID_CHANGED || IB_MOVED)
  {
      if(SOLVER_TYPE == SOLVER_TYPE_PETSC)
        solverPetsc->free();
      else
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
            //solverPetsc->printInfo();
            //cout << " aaaaaaaaa " << endl;

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

            solverPetsc->SetSolverAndParameters();

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




int  HBSplineCutFEM::prepareMatrixPattern()
{
    assert(SOLVER_TYPE == SOLVER_TYPE_PETSC);

    //time_t tstart, tend;

    int  tempDOF, domTemp, npElem, ind, size1;
    int  r, c, r1, c1, count=0, count1=0, count2=0, ii, jj, ee, dd, ind1, ind2;
    int  *tt1, *tt2, val1, val2, nnz, nnz_max_row, n1, n2, kk, e1, e2, ll, pp, aa, bb;
    int  side, NUM_NEIGHBOURS=2*DIM, start1, start2, nr1, nr2;

    node  *nd1, *nd2;

    //MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);
    //MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);

    //cout << " this_mpi_proc " << this_mpi_proc << endl;
    //cout << " n_mpi_procs " << n_mpi_procs << endl;

    int n_subdomains = n_mpi_procs, subdomain=0;


    //tstart = time(0);
    auto tstart = Clock::now();

    PetscPrintf(MPI_COMM_WORLD, " preparing cut elements \n");

    prepareCutElements();

    auto tend = Clock::now();

    //cout << " preparing()   took " << std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count() << "   milliseconds " << endl;
    PetscPrintf(MPI_COMM_WORLD, " preparing cut elements took %12.6f  milliseconds \n", std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart).count());

    setCoveringUncovering();

    //cout << " activeElements.size () = "  << activeElements.size() << endl;
    //cout << " fluidElementIds.size () = "  << fluidElementIds.size() << endl;

    //PetscSynchronizedPrintf(MPI_COMM_WORLD, " nElem = %8d \n", nElem);
    //PetscSynchronizedPrintf(MPI_COMM_WORLD, " nNode = %8d \n", nNode);


    fluidDOF = nNode*ndof;

    // for monolithic scheme add the DOF of the solids to the global DOF
    solidDOF = 0;
    if(!STAGGERED)
    {
      for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
      {
        solidDOF += ImmersedBodyObjects[bb]->GetTotalDOF();
      }
    }

    totalDOF = fluidDOF + solidDOF;

    //PetscPrintf(MPI_COMM_WORLD, "\n \t   Fluid DOF in the model   =  %8d\n\n", fluidDOF);
    //PetscPrintf(MPI_COMM_WORLD, "\n \t   Solid DOF in the model   =  %8d\n\n", solidDOF);
    PetscPrintf(MPI_COMM_WORLD, "\n \t   Total DOF in the model   =  %8d\n\n", totalDOF);

    proc_to_grid_BF.resize(nNode);
    proc_to_grid_DOF.resize(totalDOF);


    if(n_mpi_procs == 1)
    {
      elem_start = 0;
      elem_end   = nElem-1;

      nElem_local = nElem;

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
      for(ee=0; ee<nElem; ee++)
      {
        npElem_total += elems[fluidElementIds[ee]]->GlobalBasisFuncs.size();

        eptr[ee+1] = npElem_total;
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

        vecTemp2 = elems[fluidElementIds[e1]]->GlobalBasisFuncs ;

        //if(this_mpi_proc == 0)
          //printVector(vecTemp2);
        //findUnique(vecTemp2);

        npElem = vecTemp2.size();

        for(ii=0; ii<npElem; ii++)
          eind[kk+ii] = grid_to_cutfem_BF[vecTemp2[ii]] ;

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
      //int ret = METIS_PartMeshNodal(&nElem, &nNode, eptr, eind, NULL, NULL, &nParts, NULL, options, &objval, epart, npart);
      int ret = METIS_PartMeshDual(&nElem, &nNode, eptr, eind, NULL, NULL, &nodes_per_side, &n_mpi_procs, NULL, options, &objval, epart, npart);

      //idx_t  wgtflag=0, numflag=0, ncon=0, ncommonnodes=2, nparts=n_mpi_procs;
      //idx_t  *elmdist;

      //int ret = ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, NULL, &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts, NULL, NULL, options, NULL, npart, MPI_COMM_WORLD);

      //pcout <<  " ccccccccccc " << endl;

      if(ret == METIS_OK)
        PetscPrintf(MPI_COMM_WORLD, "   METIS partition routine successful \n");
      else
        PetscPrintf(MPI_COMM_WORLD, "   METIS partition routine FAILED \n");

      //if(this_mpi_proc == 0)
      //{
        //for(ee=0; ee<nElem; ee++)
          //cout << ee << '\t' << epart[ee] << endl;

        //cout << " \n\n\n\n " << endl;

        //for(ee=0; ee<nNode; ee++)
          //cout << ee << '\t' << npart[ee] << endl;
      //}

      //cout << " nNode_local = " << nNode_local << endl;

      for(e1=0; e1<nElem; e1++)
        elems[fluidElementIds[e1]]->setSubdomainId(epart[e1]);


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

      tempDOF = nNode*ndof;

      dof_map_new_to_old.resize(tempDOF, 0);
      dof_map_old_to_new.resize(tempDOF, 0);

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

      //if(this_mpi_proc==0)
      //{
        //for(ii=0; ii<node_map_new_to_old.size(); ii++)
          //cout << ii << '\t' << node_map_new_to_old[ii] << '\t' << node_map_old_to_new[ii] << endl;
        //cout << " \n\n\n\n " << endl;
      //}

      //MPI_Barrier(MPI_COMM_WORLD);

      //for(e1=0; e1<nElem; e1++)
      //{
        //ee = fluidElementIds[e1];
        ////if(elems[ee]->getSubdomainId() == this_mpi_proc)
        ////{
          //for(ii=0; ii<elems[ee]->GlobalBasisFuncs.size(); ii++)
          //{
            //elems[ee]->GlobalBasisFuncs[ii] = node_map_old_to_new[grid_to_cutfem_BF[elems[ee]->GlobalBasisFuncs[ii]]];
          //}
        ////}
      //}
    } // if(n_mpi_procs > 1)

    /////////////////////////////////////////////////////
    // mesh partitioning and node reordering is done
    /////////////////////////////////////////////////////

    for(ee=0; ee<activeElements.size(); ee++)
    {
      elems[activeElements[ee]]->initialiseDOFvalues();
    }

    PetscSynchronizedPrintf(MPI_COMM_WORLD, "\n    preparing matrix pattern    \n");
    //printf("\n element DOF values initialised \n\n");
    //printf("\n Finding Global positions \n\n");


    for(ii=0; ii<gridBF1; ii++)
      grid_to_proc_BF[ii]   =  node_map_old_to_new[grid_to_cutfem_BF[ii]];

    for(ii=0; ii<nNode; ii++)
      proc_to_grid_BF[ii]   =  cutfem_to_grid_BF[node_map_new_to_old[ii]];

    kk = gridBF1*ndof;
    for(ii=0; ii<kk; ii++)
      grid_to_proc_DOF[ii]  =  dof_map_old_to_new[grid_to_cutfem_DOF[ii]];

    for(ii=0; ii<totalDOF; ii++)
      proc_to_grid_DOF[ii]  =  cutfem_to_grid_DOF[dof_map_new_to_old[ii]];

    //return 1;

    vector<vector<int> > DDconn;

    DDconn.resize(totalDOF);


    for(ee=0; ee<nElem; ee++)
    {
        nd1 = elems[fluidElementIds[ee]];

        //PetscPrintf(MPI_COMM_WORLD, " %d \t %d \t %d \n", ee, fluidElementIds[ee], nd1->getID());

        val1 =  nd1->getNsize2();
        tt1  =  &(nd1->forAssyVec[0]);

        //if(this_mpi_proc == 0)
          //printVector(nd1->forAssyVec);

        if( nd1->domNums[0] == 0 )
        {
          for(ii=0; ii<val1; ii++)
          {
            r = grid_to_proc_DOF[tt1[ii]];

            for(jj=0; jj<val1; jj++)
              DDconn[r].push_back(grid_to_proc_DOF[tt1[jj]]);
          }
        }

        //PetscPrintf(MPI_COMM_WORLD, " aaaaaaaaaaaaa \n");

      // connectivity for ghost-penalty terms
      if( nd1->isCutElement() )
      {
        for(side=0; side<NUM_NEIGHBOURS; side++)
        {
          nd2 = nd1->getNeighbour(side);

          if( (nd2 != NULL) && !(nd2->isGhost()) && nd2->isLeaf() && (nd2->isCutElement() || nd2->domNums[0] == 0) )
          {
              nr1 = nd1->forAssyVec.size();
              nr2 = nd2->forAssyVec.size();

              for(ii=0; ii<nr1; ii++)
              {
                r = grid_to_proc_DOF[nd1->forAssyVec[ii]];

                for(jj=0; jj<nr2; jj++)
                {
                  c = grid_to_proc_DOF[nd2->forAssyVec[jj]];

                  DDconn[r].push_back(c);
                  DDconn[c].push_back(r);
                } // for(jj=0; jj<nr2; jj++)
              } // for(ii=0; ii<nr1; ii++)
          }
        } //for(side=0; side<NUM_NEIGHBOURS; side++)
      } //if( nd1->isCutElement() )
      //PetscPrintf(MPI_COMM_WORLD, " aaaaaaaaaaaaa \n");
    } // for(e=0;e<fluidElementIds.size();e++)



    if(!STAGGERED)
    {
      start1 = fluidDOF;
      start2 = fluidDOF;
      domTemp = 0;
      nr1 = 1;

      ImmersedIntegrationElement  *lme;
      myPoly *poly;

      vector<double>  gausspoints, gaussweights;
      vector<int>  vecIntTemp;
      int  nlb, nGauss, elnum, gp;
      double  detJ;
      VectorXd  Nb, dNb, xx, yy;

      for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
      {
        // matrix pattern for coupling terms
        //////////////////////////////////////////////////////////////////
        
        nlb = ImmersedBodyObjects[bb]->ImmIntgElems[0]->pointNums.size();

        Nb.resize(nlb);
        dNb.resize(nlb);
        xx.resize(nlb);
        yy.resize(nlb);

        nGauss = (int) cutFEMparams[1];

        getGaussPoints1D(nGauss, gausspoints, gaussweights);

        for(aa=0; aa<ImmersedBodyObjects[bb]->ImmIntgElems.size(); aa++)
        {
          lme = ImmersedBodyObjects[bb]->ImmIntgElems[aa];
          poly = ImmersedBodyObjects[bb]->ImmersedFaces[aa];

          //cout << bb << '\t' << aa << '\t' << lme->isActive() << endl;

          nr1 = lme->pointNums.size();
          tt1 =  &(lme->pointNums[0]);

          if( lme->isActive() )
          {
            for(gp=0;gp<nGauss;gp++)
            {
              //computeLagrangeBFs1D2(nlb-1, gausspoints[gp], &xx(0), &yy(0), &Nb(0), &dNb(0), detJ);

              param[0] = gausspoints[gp];
              poly->computeBasisFunctions(param, geom, Nb, detJ);

              elnum = findCellNumber(geom);

              nd2 = elems[elnum];

              nr2  =  nd2->getNsize2();
              tt2  =  &(nd2->forAssyVec[0]);

              for(ii=0; ii<nr1; ii++)
              {
                r = start1 + tt1[ii]*2;

                for(jj=0; jj<nr2; jj++)
                {
                  c = grid_to_proc_DOF[tt2[jj]];

                  DDconn[r].push_back(c);
                  DDconn[r+1].push_back(c);

                  DDconn[c].push_back(r);
                  DDconn[c].push_back(r+1);
                } // for(jj=0; jj<nr2; jj++)
              } // for(ii=0; ii<nr1; ii++)
            }
          }
        } // for(aa=0; aa<ImmersedBodyObjects[bb]->ImmIntgElems.size(); aa++)

        // matrix pattern for solid elements
        //////////////////////////////////////////////////////////////////

        for(ii=0; ii<ImmersedBodyObjects[bb]->forAssyMat.size(); ii++)
        {
          tt2 = &(ImmersedBodyObjects[bb]->forAssyMat[ii][0]);
          for(jj=0; jj<ImmersedBodyObjects[bb]->forAssyMat[ii].size(); jj++)
          {
            DDconn[start1+ii].push_back(start2+tt2[jj]);
          }
        }

        start1 += ImmersedBodyObjects[bb]->GetTotalDOF();
        start2 = start1;

      } // for(bb=0; bb<ImmersedBodyObjects.size(); bb++)
    } // if(!STAGGERED)

    PetscSynchronizedPrintf(MPI_COMM_WORLD, "\n   Finding global positions DONE \n\n");

    //tend = time(0); 
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
    nnz_max_row = 0;
    for(ii=row_start; ii<=row_end; ii++)
    {
      findUnique(DDconn[ii]);

      size1 = DDconn[ii].size();

      nnz_max_row = max(nnz_max_row, jj);

      count_diag=0, count_offdiag=0;
      for(jj=0; jj<size1; jj++)
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

    //cout << " iiiiiiiiiiiii " << count_diag << '\t' << count_offdiag << endl;

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

    ierr = MatMPIAIJSetPreallocation(solverPetsc->mtx, 20, d_nnz, 20, o_nnz);CHKERRQ(ierr);
    ierr = MatSeqAIJSetPreallocation(solverPetsc->mtx, 20, d_nnz);CHKERRQ(ierr);

    //cout << " jjjjjjjjjjjjjjjj " << endl;

    PetscInt  colTemp[nnz_max_row];
    PetscScalar  arrayD[nnz_max_row];

    for(jj=0; jj<nnz_max_row; jj++)
      arrayD[jj] = 0.0;

    for(ii=row_start; ii<=row_end; ii++)
    {
      size1 = DDconn[ii].size();

      for(jj=0; jj<size1; jj++)
        colTemp[jj] = DDconn[ii][jj];

      ierr = MatSetValues(solverPetsc->mtx, 1, &ii, size1, colTemp, arrayD, INSERT_VALUES);
    }

    //cout << " iiiiiiiiiiiii " << endl;

    VecCreate(PETSC_COMM_WORLD, &solverPetsc->soln);
    VecCreate(PETSC_COMM_WORLD, &solverPetsc->solnPrev);
    VecCreate(PETSC_COMM_WORLD, &solverPetsc->rhsVec);
    VecCreate(PETSC_COMM_WORLD, &solverPetsc->reac);

    ierr = VecSetSizes(solverPetsc->soln, ndofs_local, totalDOF); CHKERRQ(ierr);
    //ierr = VecSetSizes(solnPrev, ndofs_local, totalDOF); CHKERRQ(ierr);
    //ierr = VecSetSizes(rhsVec, ndofs_local, totalDOF); CHKERRQ(ierr);
    //ierr = VecSetSizes(solverPetsc->reac, ndofs_local, totalDOF); CHKERRQ(ierr);

    ierr = VecSetFromOptions(solverPetsc->soln);CHKERRQ(ierr);
    ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->rhsVec);CHKERRQ(ierr);
    ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->solnPrev);CHKERRQ(ierr);
    ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->reac);CHKERRQ(ierr);

    solverPetsc->currentStatus = PATTERN_OK;

    soln.resize(totalDOF);
    soln.setZero();
    solnInit = soln;

    GRID_CHANGED = IB_MOVED = false;

    //printf("\n     HBSplineCutFEM::prepareMatrixPattern()  .... FINISHED ...\n\n");

  return 1;
}




/*
int  HBSplineCutFEM::prepareMatrixPattern()
{
    for petsc - working
    
    assert(SOLVER_TYPE == SOLVER_TYPE_PETSC);

    time_t tstart, tend; 

    int  tempDOF, domTemp, npElem, ind, size1;
    int  r, c, r1, c1, count=0, count1=0, count2=0, ii, jj, ee, dd, ind1, ind2;
    int  *tt1, *tt2, val1, val2, nnz, nnz_max_row, n1, n2, kk, e1, e2, ll, pp, aa, bb;
    int  side, NUM_NEIGHBOURS=2*DIM, start1, start2, nr1, nr2;

    node  *nd1, *nd2;

    //MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);
    //MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);

    //cout << " this_mpi_proc " << this_mpi_proc << endl;
    //cout << " n_mpi_procs " << n_mpi_procs << endl;

    int n_subdomains = n_mpi_procs, subdomain=0;


    tstart = time(0);

    prepareCutElements();

    tend = time(0);
    //printf("HBSplineCutFEM::prepareCutElements() took %8.4f second(s) \n ", difftime(tend, tstart) );

    //cout << " activeElements.size () = "  << activeElements.size() << endl;
    //cout << " fluidElementIds.size () = "  << fluidElementIds.size() << endl;

    PetscSynchronizedPrintf(MPI_COMM_WORLD, " nElem = %8d \n", nElem);
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
      //int ret = METIS_PartMeshNodal(&nElem, &nNode, eptr, eind, NULL, NULL, &nParts, NULL, options, &objval, epart, npart);
      int ret = METIS_PartMeshDual(&nElem, &nNode, eptr, eind, NULL, NULL, &nodes_per_side, &n_mpi_procs, NULL, options, &objval, epart, npart);

      //idx_t  wgtflag=0, numflag=0, ncon=0, ncommonnodes=2, nparts=n_mpi_procs;
      //idx_t  *elmdist;

      //int ret = ParMETIS_V3_PartMeshKway(elmdist, eptr, eind, NULL, &wgtflag, &numflag, &ncon, &ncommonnodes, &nparts, NULL, NULL, options, NULL, npart, MPI_COMM_WORLD);

      //pcout <<  " ccccccccccc " << endl;

      if(ret == METIS_OK)
        PetscPrintf(MPI_COMM_WORLD, "   METIS partition routine successful \n");
      else
        PetscPrintf(MPI_COMM_WORLD, "   METIS partition routine FAILED \n");

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
        elems[activeElements[e1]]->setSubdomainId(epart[e1]);


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
        //if(elems[ee]->getSubdomainId() == this_mpi_proc)
        //{
          for(ii=0; ii<elems[ee]->GlobalBasisFuncs.size(); ii++)
          {
            elems[ee]->GlobalBasisFuncs[ii] = node_map_old_to_new[elems[ee]->GlobalBasisFuncs[ii]];
          }
        //}
      }
    } // if(n_mpi_procs > 1)

    /////////////////////////////////////////////////////
    // mesh partitioning and node reordering is done
    /////////////////////////////////////////////////////

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

    for(ee=0; ee<activeElements.size(); ee++)
    {
        nd1 = elems[activeElements[ee]];

        val1 =  nd1->getNsize2();
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
      if( nd1->isCutElement() )
      {
        for(side=0; side<NUM_NEIGHBOURS; side++)
        {
          nd2 = nd1->getNeighbour(side);

          if( (nd2 != NULL) && !(nd2->isGhost()) && nd2->isLeaf() && (nd2->isCutElement() || nd2->domNums[0] == 0) )
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
      } //if( nd1->isCutElement() )
    } // for(e=0;e<activeElements.size();e++)

    PetscSynchronizedPrintf(MPI_COMM_WORLD, "\n    Global positions DONE for individual domains \n\n");

    grid_to_cutfem_DOFprev = grid_to_cutfem_DOF;

    cutfem_to_grid_DOF.clear();
    grid_to_cutfem_DOF.assign(tempDOF, -1);


    fluidDOF = 0;
    for(ee=0; ee<tempDOF; ee++)
    {
      findUnique(DDconnLoc[ee]);

      if(DDconnLoc[ee].size() > 0)
      {
        grid_to_cutfem_DOF[ee]  = fluidDOF++;
        cutfem_to_grid_DOF.push_back(ee);
      }
      else
        grid_to_cutfem_DOF[ee] = -1;
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
      if(grid_to_cutfem_DOF[ee] != -1)
      {
        val1 = DDconnLoc[ee].size();
        //cout << " val1 " << ee << '\t' << val1 << endl;
        r = grid_to_cutfem_DOF[ee];
        for(ii=0; ii<val1; ii++)
        {
          DDconn[r].push_back( grid_to_cutfem_DOF[DDconnLoc[ee][ii]] );
        }
      }
    } //for(ee=0; ee<DDconnLoc[dd].size(); ee++)

    //printf("\n Finding Global positions DONE \n\n");

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
    nnz_max_row = 0;
    for(ii=row_start; ii<=row_end; ii++)
    {
      size1 = DDconn[ii].size();

      nnz_max_row = max(nnz_max_row, jj);

      count_diag=0, count_offdiag=0;
      for(jj=0; jj<size1; jj++)
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

    PetscInt  colTemp[nnz_max_row];
    PetscScalar  arrayD[nnz_max_row];

    for(jj=0; jj<nnz_max_row; jj++)
      arrayD[jj] = 0.0;

    for(ii=row_start; ii<=row_end; ii++)
    {
      size1 = DDconn[ii].size();

      for(jj=0; jj<size1; jj++)
        colTemp[jj] = DDconn[ii][jj];

      ierr = MatSetValues(solverPetsc->mtx, 1, &ii, size1, colTemp, arrayD, INSERT_VALUES);
    }

    VecCreate(PETSC_COMM_WORLD, &solverPetsc->soln);
    VecCreate(PETSC_COMM_WORLD, &solverPetsc->solnPrev);
    VecCreate(PETSC_COMM_WORLD, &solverPetsc->rhsVec);
    VecCreate(PETSC_COMM_WORLD, &solverPetsc->reac);

    ierr = VecSetSizes(solverPetsc->soln, ndofs_local, totalDOF); CHKERRQ(ierr);
    //ierr = VecSetSizes(solnPrev, ndofs_local, totalDOF); CHKERRQ(ierr);
    //ierr = VecSetSizes(rhsVec, ndofs_local, totalDOF); CHKERRQ(ierr);
    //ierr = VecSetSizes(solverPetsc->reac, ndofs_local, totalDOF); CHKERRQ(ierr);

    ierr = VecSetFromOptions(solverPetsc->soln);CHKERRQ(ierr);
    ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->rhsVec);CHKERRQ(ierr);
    ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->solnPrev);CHKERRQ(ierr);
    ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->reac);CHKERRQ(ierr);

    solverPetsc->currentStatus = PATTERN_OK;

    soln.resize(totalDOF);
    soln.setZero();
    solnInit = soln;

    GRID_CHANGED = IB_MOVED = false;

    //printf("\n     HBSplineCutFEM::prepareMatrixPattern()  .... FINISHED ...\n\n");

  return 1;
}
*/





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

        val1 =  nd1->getNsize2();
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

      if( nd1->isCutElement() )
      {
        for(side=0; side<NUM_NEIGHBOURS; side++)
        {
          nd2 = nd1->getNeighbour(side);

          if( (nd2 != NULL) && !(nd2->isGhost()) && nd2->isLeaf() && (nd2->isCutElement() || nd2->domNums[0] == 0) )
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
      } //if( nd1->isCutElement() )
      //
    } // for(e=0;e<activeElements.size();e++)

    printf("\n Global positions DONE for individual domains \n\n");

    grid_to_cutfem_DOFprev = grid_to_cutfem_DOF;

    cutfem_to_grid_DOF.clear();
    grid_to_cutfem_DOF.assign(tempDOF, -1);


    fluidDOF = 0;
    for(ee=0; ee<tempDOF; ee++)
    {
      findUnique(DDconnLoc[ee]);

      if(DDconnLoc[ee].size() > 0)
      {
        grid_to_cutfem_DOF[ee]  = fluidDOF++;
        cutfem_to_grid_DOF.push_back(ee);
      }
      else
        grid_to_cutfem_DOF[ee] = -1;
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
      if(grid_to_cutfem_DOF[ee] != -1)
      {
        val1 = DDconnLoc[ee].size();
        //cout << " val1 " << ee << '\t' << val1 << endl;
        r = grid_to_cutfem_DOF[ee];
        for(ii=0; ii<val1; ii++)
        {
          DDconn[r].push_back( grid_to_cutfem_DOF[DDconnLoc[ee][ii]] );
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

          if( nd1->isCutElement() )
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
          } //if( nd1->isCutElement() )
      } // for(e=0;e<activeElements.size();e++)
    } // if(!STAGGERED)

    PetscPrintf(MPI_COMM_WORLD, "\n Finding Global positions DONE \n\n");

    VectorXi  nnzVec(totalDOF);

    nnz = 0;
    for(ii=0;ii<totalDOF;ii++)
    {
      findUnique(DDconn[ii]);
      nnzVec[ii] = DDconn[ii].size();
      nnz += nnzVec[ii];
    }
    cout << " nnz " << nnz << endl;

    tend = time(0); 
    //cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< endl;

    bool pp1=false;
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

      val1 =  elems[ee]->getNsize2();
      tt1  =  &(elems[ee]->forAssyVec[0]);
      //cout << e << '\t' << ee << '\t' << val1 << endl;
      //printVector(elems[ee]->forAssyVec);

      elems[ee]->findCutElemType();

      if( elems[ee]->isCutElement() )
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
      } // if( !elems[ee]->isCutElement() )
      else
      {
        if( elems[ee]->getDomainNumber() == 0 )
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
    grid_to_cutfem_DOF.resize(2);
    cutfem_to_grid_DOF.resize(2);
    domainTotalDOF[0] = domainTotalDOF[1] = 0;

    for(ee=0; ee<tempDOF; ee++)
    {
      findUnique(DDconn1[ee]);
      if(DDconn1[ee].size() > 0)
      {
        grid_to_cutfem_DOF[0].push_back(domainTotalDOF[0]++);
        cutfem_to_grid_DOF[0].push_back(ee);
      }
      else
        grid_to_cutfem_DOF[0].push_back(-1);

      findUnique(DDconn2[ee]);
      if(DDconn2[ee].size() > 0)
      {
        grid_to_cutfem_DOF[1].push_back(domainTotalDOF[1]++);
        cutfem_to_grid_DOF[1].push_back(ee);
      }
      else
        grid_to_cutfem_DOF[1].push_back(-1);
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

       printf("   grid_to_cutfem_DOF[0] ..... : \n\n");
       for(ii=0;ii<grid_to_cutfem_DOF[0].size();ii++)
          cout << ii << '\t' << grid_to_cutfem_DOF[0][ii] << endl;

       printf("\n\n\n");

       printf("   grid_to_cutfem_DOF[1] ..... : \n\n");
       for(ii=0;ii<grid_to_cutfem_DOF[1].size();ii++)
          cout << ii << '\t' << grid_to_cutfem_DOF[1][ii] << endl;

       printf("\n\n\n");

       printf("   cutfem_to_grid_DOF[0] ..... : \n\n");
       for(ii=0;ii<cutfem_to_grid_DOF[0].size();ii++)
          cout << ii << '\t' << cutfem_to_grid_DOF[0][ii] << endl;

       printf("\n\n\n");

       printf("   cutfem_to_grid_DOF[1] ..... : \n\n");
       for(ii=0;ii<cutfem_to_grid_DOF[1].size();ii++)
          cout << ii << '\t' << cutfem_to_grid_DOF[1][ii] << endl;

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
        r = grid_to_cutfem_DOF[0][ee];
        for(ii=0; ii<val1; ii++)
        {
          DDconn[r].push_back(grid_to_cutfem_DOF[0][DDconn1[ee][ii]]);
        }
      }

      val1 = DDconn2[ee].size();
      kk = domainTotalDOF[0];
      if(val1 > 0)
      {
        r = kk + grid_to_cutfem_DOF[1][ee] ;
        //cout << " val1 " << ee << '\t' << r << endl;
        for(ii=0; ii<val1; ii++)
        {
          DDconn[r].push_back(kk+grid_to_cutfem_DOF[1][DDconn2[ee][ii]]);
        }
      }
    } //for(ee=0; ee<tempDOF; ee++)


    // generate matrix pattern for the coupling terms
    //
    for(e=0;e<activeElements.size();e++)
    {
      ee = activeElements[e];

      if( elems[ee]->isCutElement() )
      {
        val1 =  elems[ee]->getNsize2();
        tt1  =  &(elems[ee]->forAssyVec[0]);
        //cout << e << '\t' << ee << '\t' << val1 << endl;
        //printVector(elems[ee]->forAssyVec);

        kk = domainTotalDOF[0];

        for(ii=0;ii<val1;ii++)
        {
          r = grid_to_cutfem_DOF[0][tt1[ii]];

          for(jj=0;jj<val1;jj++)
          {
            c = kk + grid_to_cutfem_DOF[1][tt1[jj]] ;
            DDconn[r].push_back(c);
            DDconn[c].push_back(r);
          }
        }
      } // if( !elems[ee]->isCutElement() )
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




