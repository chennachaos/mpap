

#include "HBSplineFEM.h"
#include "SolverPardisoEigen.h"
#include "SolverMA41Eigen.h"
#include "MpapTime.h"
#include "TimeFunction.h"


extern MpapTime           mpapTime;



void HBSplineFEM::setTimeParam()
{
  //cout << " HBSplineBase::setTimeParam() ... STARTED " << endl;

  SolnData.setTimeParam();

  for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
    ImmersedBodyObjects[bb]->setTimeParam();

  //cout << " HBSplineBase::setTimeParam() ... FINISHED " << endl;

  return;
}



void HBSplineFEM::timeUpdate()
{
  //cout << " HBSplineBase::timeUpdate() ... STARTED " << endl;

  IB_MOVED = false;

  firstIter = true;
  localStiffnessError = 0;
  iterCount = 1;
  filecount++;

  //ImmersedBoundaryBodyForceLSFEM();

  //if(!STAGGERED)
  //{
    for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
      ImmersedBodyObjects[bb]->timeUpdate();
  //}

  SolnData.timeUpdate();

  //cout << " zzzzzzzzzzzzzz " << endl;
  //if(LSFEM_FLAG)
    //ImmersedBoundaryBodyForceLSFEM();
  //else
  //{
    if(STAGGERED)
    {
      solveSolidProblem();
      updateImmersedPointPositions();
      //IB_MOVED = true;
    }
  //}
  //

  //ImmersedBoundaryBodyForceLSFEM();

  updateIterStep();

  //cout << " HBSplineBase::timeUpdate() ... FINISHED " << endl;

  return;
}



void HBSplineFEM::updateIterStep()
{
  int kk, bb, ee;

  SolnData.updateIterStep();

  //cout << " kkkkkkkkkkk " << endl;

  for(bb=0;bb<ImmersedBodyObjects.size();bb++)
    ImmersedBodyObjects[bb]->updateForce();

  //cout << " qqqqqqqqqqq " << endl;

  if(!STAGGERED)
  {
    for(bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
      //cout << " AAAAAAAAAAA " << bb << endl;
      ImmersedBodyObjects[bb]->updateDisplacement(&(SolnData.var4(0)));
      //cout << " AAAAAAAAAAA " << bb << endl;
      ImmersedBodyObjects[bb]->updateIterStep();
      //cout << " AAAAAAAAAAA " << bb << endl;
    }
    //cout << " AAAAAAAAAAA " << bb << endl;
    updateImmersedPointPositions();
  }

  //cout << " kkkkkkkkkkk " << GRID_CHANGED << '\t' << IB_MOVED << endl;

  //cout << " PPPPPPPPPPPPPPPPP " << firstIter << '\t' << SOLVER_TYPE << '\t' << GRID_CHANGED << '\t' << IB_MOVED << endl;
  if(GRID_CHANGED || IB_MOVED)
  {
      solverEigen->free();

      prepareMatrixPattern();

      switch(slv_type)
      {
          case  1: // MA41 ..........................
          case  4: // SolverEigen ..........................

            if(solverEigen->initialise(0,0,totalDOF) != 0)
              return;

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

          default: // invalid slv_type ...................

            cout << " this solver has not been implemented yet!\n\n";

          break;
      }

      solverOK = true;
 
      if(solverEigen != NULL)
        solverEigen->checkIO = true;
  }

  //cout << " qqqqqqqqqqq " << endl;

  return;
}




void HBSplineFEM::reset()
{
  SolnData.reset();

  for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
    //if(ImmersedBodyObjects[bb]->getNdof() > 0)
      ImmersedBodyObjects[bb]->reset();

  return;
}


void HBSplineFEM::printData(int index, int patchnum)
{
   if(index == 1)  // print global stiffness matrix
   {
      cout << "      global stiffness matrix "  << endl;
      printf("\n\n");
      solverEigen->printMatrix();
      printf("\n\n");
      return;
   }
   if(index == 2)  // print global stiffness matrix
   {
      printf("\n\n");
//      printMatrix(globalM);
      printf("\n\n");
      return;
   }
   if(index == 3)  // print global rhsVec vector
   {
      printf("\n\n");
      printVector(solverEigen->rhsVec);
      printf("\n\n");
      return;
   }
   if(index == 4)  // print solution vector
   {
      printf("\n\n");
      printVector(soln);
      printf("\n\n");
      return;
   }

   if(index == 5)
   {
      ofstream fout("ctrl-point-data.dat");

      if(fout.fail())
      {
         cout << " Could not open the Output file" << endl;
      exit(1);
      }

      fout.setf(ios::fixed);
      fout.setf(ios::showpoint);
      fout.precision(8);

      //for(int ii=0;ii<soln.rows();ii++)
        //fout << (ii+1) << '\t' <<  soln(ii) << '\t' <<  rhsVec(ii) << endl;;

      fout.close();
   }
   if(index == 6)
   {
      ofstream fout("ctrl-point-data.dat");

      if(fout.fail())
      {
         cout << " Could not open the Output file" << endl;
      exit(1);
      }

      fout.setf(ios::fixed);
      fout.setf(ios::showpoint);
      fout.precision(8);

      //for(int ii=0;ii<rhsVec.rows();ii++)
        //fout << (ii+1) << '\t' <<  rhsVec(ii) << endl;;

      fout.close();
   }

   //if(index == 8 || index == 9) // write matrix pattern to an output file 
     //solver->printMatrixPatternToFile();

   return;
}



void  HBSplineFEM::writeReadResult(int index, MyString &fileName)
{
    char fct[] = "HBSplineFEM::writeReadResult";

    int  ii, jj, ind, ee, bb, count;
    
    double  val1, val2;

   if(index == 1) // write result to a file
   {
      ofstream  fout;

      char tmp[500], ff[100];
      sprintf(ff,"%s%s%06d%s", fileName.asCharArray(),"-",filecount, ".dat");

      MyString  fName;
      
      fName = ff;

      //fout.open(fName.asCharArray());
      fout.open(ff);

      if(!fout)
      {
         prgWarning(1,"HBSplineFEM::writeReadResult","could not open file for writing!");
         return;
      }

      /*
      sprintf(tmp,"Model");
      fout << tmp << "\n";
      sprintf(tmp,"%10.6f X %10.6f", gridLEN[0], gridLEN[1]);
      //fout << tmp << "\t";
      sprintf(&(tmp[strlen(tmp)]),"\t %5d X %5d", nelem[0], nelem[1]);
      //fout << tmp << "\t";
      sprintf(&(tmp[strlen(tmp)]),"\t Q_%1d \t Level=%1d", degree[0], CURRENT_LEVEL);
      fout << tmp << "\n\n";
      */

      sprintf(tmp,"Time");
      fout << tmp << "\n";
      sprintf(tmp," %14.8f", mpapTime.cur);
      fout << tmp << "\n\n";
      sprintf(tmp,"Result-Fluid");
      fout << tmp << "\n";

      //for(ii=0;ii<(velDOF/ndof);ii++)
      for(ii=0;ii<(fluidDOF/ndof);ii++)
      {
         ind = ndof*ii;

         //sprintf(tmp," %14.8f", FluidSolnData.bodyForce(ind));
         //sprintf(&(tmp[strlen(tmp)])," %14.8f", FluidSolnData.bodyForce(ind+1));
         sprintf(tmp," %14.8f", SolnData.var1(ind));
         sprintf(&(tmp[strlen(tmp)])," %14.8f", SolnData.var1(ind+1));
         sprintf(&(tmp[strlen(tmp)])," %14.8f", SolnData.var1(ind+2));
         //sprintf(&(tmp[strlen(tmp)])," %14.8f", FluidSolnData.var2(ii));
         fout << tmp << "\n";
         //cout << ii << '\t' << ind << '\t' << FluidSolnData.soln(ind) << '\t' << FluidSolnData.soln(ind+1) << '\t' << FluidSolnData.soln(ind+2) << endl ;
      }
      fout << "\n\n" ;
      //cout << " AAAAAAAAAAAAA " << endl;

      sprintf(tmp,"Result-Solid");
      fout << tmp << "\n\n";

      ImmersedIntegrationElement  *lme;

      for(bb=0;bb<ImmersedBodyObjects.size();bb++)
      {
          if( ImmersedBodyObjects[bb]->IsBoundaryConditionTypeLagrange() )
          {
          for(ii=0;ii<ImmersedBodyObjects[bb]->GetNumNodes();ii++)
          {
            //ind = ImmersedBodyObjects[bb]->GlobalPointNumbers[ii];
            ind = ii;
            sprintf(tmp," %5d", ind);
            ind = ind*DIM;
            sprintf(&(tmp[strlen(tmp)])," %14.8f", SolnData.var3(ind));
            sprintf(&(tmp[strlen(tmp)])," %14.8f", SolnData.var3(ind+1));

            sprintf(&(tmp[strlen(tmp)])," %14.8f", ImmersedBodyObjects[bb]->GeomData.NodePosCur[ii][0]);
            sprintf(&(tmp[strlen(tmp)])," %14.8f", ImmersedBodyObjects[bb]->GeomData.NodePosCur[ii][1]);
            sprintf(&(tmp[strlen(tmp)])," %14.8f", ImmersedBodyObjects[bb]->GeomData.specValCur[ii][0]);
            sprintf(&(tmp[strlen(tmp)])," %14.8f", ImmersedBodyObjects[bb]->GeomData.specValCur[ii][1]);

            fout << tmp << "\n";
          }
          }
      }
      fout << "\n\n" ;
      sprintf(tmp,"End");
      fout << tmp ;

      fout.close();
   }
   else // read result from a file
   {
      ifstream  Ifile;

      MyString line, tmpl;

      MyStringList  names;
      names.addNew("Model", "Time", "Result-Fluid", "Result-Solid", "End");

      List<Vector<double> > lvdTmp;

      Ifile.open(fileName.asCharArray());

      if(!Ifile) prgError(1,fct,"failed to open the input file.");

      line.read(Ifile);

      do
      {
         lvdTmp.free();

         cout << " >" << line.stripToMin() << "<\n\n";

         switch(names.whichBegins(line))
         {
            case 0:  // read Model

               //if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
                 //prgError(1,fct,"invalid input in 'HBSplineFEM::writeReadResult'!");
               
               //line.read(Ifile);
               cout << " AAAAAAAAAAA " << endl;

               break;

            case 1:  // read time

               if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
                 prgError(1,fct,"invalid input in 'HBSplineFEM::writeReadResult'!");

               break;

            case 2:  //read result Fluid

               cout << line << endl;

               if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
                 prgError(1,fct,"invalid input in 'HBSplineFEM::writeReadResult'!");
               
               cout << lvdTmp.n << endl;

               if(lvdTmp.n >= fluidDOF/ndof)
               {
                 //for(ii=0;ii<lvdTmp.n;ii++)
                 for(ii=0;ii<fluidDOF/ndof;ii++)
                 {
                    ind = ii*ndof;
                    for(jj=0;jj<2;jj++)
                      SolnData.bodyForce(ind+jj) = lvdTmp[ii][jj]*1.0;
                    for(jj=0;jj<3;jj++)
                      SolnData.var1(ind+jj) = lvdTmp[ii][jj+2];
                 }
               }

               //printVector(FluidSolnData.soln);

               break;

            case 3:  //read result Solid

               cout << line << endl;

               if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
                 prgError(1,fct,"invalid input in 'HBSplineFEM::writeReadResult'!");

               /*
               if(lvdTmp.n == IBpoints.size())
               {
                 for(ii=0;ii<lvdTmp.n;ii++)
                 {
                    ind = ii*ndof;
                    for(jj=0;jj<2;jj++)
                      IBpoints[ii]->SetDisplacement(jj, lvdTmp[ii][jj]*1.0);
                 }
               }
               */
               break;


            case  4: //End

               break;
          }
      }
      while(prgNextDataKey(Ifile, line));

      Ifile.close();
   }
}




