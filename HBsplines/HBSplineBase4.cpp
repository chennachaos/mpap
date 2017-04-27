
#include "HBSplineBase.h"



void  HBSplineBase::computeTotalBodyForce(int index)
{
  double val, forceTotal[3];
  int dir, ee;
   
  for(dir=0;dir<DIM;dir++)
  {
    val = 0.0;
    for(ee=0;ee<activeElements.size();ee++)
      val += elems[activeElements[ee]]->computeTotalBodyForce(index, dir);

    forceTotal[dir] = val;
  }

  char        tmp[200];
  MyString    tmpStr;
    
  sprintf(tmp," \t %12.6E \t %12.6E ", forceTotal[0], forceTotal[1]);
  tmpStr.append(tmp);

  prgWriteToTFile(tmpStr);

  return;
}







void  HBSplineBase::updateImmersedPointPositions()
{
  IB_MOVED = false;
  for(int bb=0; bb<nImmSolids; bb++)
  {
    IB_MOVED = (IB_MOVED || ImmersedBodyObjects[bb]->updatePointPositions() );
  }
  //cout << " HBSplineFEM::updateImmersedPointPositions() " << endl;

  return;
}




void  HBSplineBase::writeNodalData()
{
  writeFluidOutput();
  writeImmersedSolidOutput();

  return;
}


void  HBSplineBase::writeFluidOutput()
{
  return;
}



void  HBSplineBase::writeImmersedSolidOutput()
{
  for(int bb=0; bb<nImmSolids; bb++)
    ImmersedBodyObjects[bb]->writeOutput();

  return;
}




void  HBSplineBase::printResultAtPoint(int index, double u1, double v1, double w1)
{
/*
    int  dd, ii, jj, kk, ll, count, nlocal, ind1, ind2, dir, resln[2];
    dir = index-1;
    
    resln[0] = resln[1] = 5;

    double   fact[2], incr1, incr2, start1, start2, val, *tmp1, *tmp2;
    double   N1[degree[0]+1], N2[degree[1]+1];

    vector<double>  U, V, uu, vv;
    vector<vector<double> >  outp2;
    vector<double>  xx, yy;
    
    outp2.resize(4);

    tmp1 = elems[0]->GetKnots(0);
    tmp2 = elems[0]->GetKnots(1);

    incr1 = tmp1[1] - tmp1[0];
    incr2 = tmp2[1] - tmp2[0];

    int  m1, m2, n1, n2, ngbf1;
    
    ngbf1 = nelem[0] + degree[0];

    m1 = nelem[0] + 2*degree[0];
    m2 = nelem[1] + 2*degree[1];

    U.resize(m1+1);
    V.resize(m2+1);
    
    start1 = -degree[0] * incr1;
    start2 = -degree[1] * incr2;

    for(ii=0;ii<U.size();ii++)
      U[ii] = start1 + ((double) ii) * incr1;
           
    for(ii=0;ii<V.size();ii++)
      V[ii] = start2 + ((double) ii) * incr2;

    incr1 = incr1/resln[0];
    incr2 = incr2/resln[1];

    create_vector(0.0, 1.0, incr1, uu);
    create_vector(0.0, 1.0, incr2, vv);


                 BasisFuns(&V[0], V.size(), degree[1], 0.5, N2);

                 n2 = FindSpan(&V[0], V.size(), degree[1], 0.5) - degree[1];
                 
                 for(ii=0;ii<uu.size();ii++)
                 {
                    BasisFuns(&U[0], U.size(), degree[0], uu[ii], N1);
                    
                    n1 = FindSpan(&U[0], U.size(), degree[0], uu[ii]) - degree[0];
                    fact[0] = fact[1] = 0.0;
                    for(ll=0;ll<=degree[1];ll++)
                    {
                       ind2 = (n2+ll)*ngbf1 + n1;
                       for(kk=0;kk<=degree[0];kk++)
                       {
                          val = N2[ll] * N1[kk];
                          ind1 = ind2 + kk;
                          
                          fact[0] += val * soln(ndf*ind1);
                          fact[1] += val * soln(ndf*ind1+1);
                       }
                    }

                    outp2[0].push_back(fact[0]);
                    outp2[1].push_back(fact[1]);

                    param[0] = uu[ii];
                    param[1] = 0.5;
                    ComputeGeometry(param, geom);

                    xx.push_back(geom[0]);
                 }


                 BasisFuns(&U[0], U.size(), degree[0], 0.5, N1);
                    
                 n1 = FindSpan(&U[0], U.size(), degree[0], 0.5) - degree[0];
                 
                 for(jj=0;jj<vv.size();jj++)
                 {
                    BasisFuns(&V[0], V.size(), degree[1], vv[jj], N2);

                    n2 = FindSpan(&V[0], V.size(), degree[1], vv[jj]) - degree[1];

                    fact[0] = fact[1] = 0.0;
                    for(ll=0;ll<=degree[1];ll++)
                    {
                       ind2 = (n2+ll)*ngbf1 + n1;
                       for(kk=0;kk<=degree[0];kk++)
                       {
                          val = N2[ll] * N1[kk];
                          ind1 = ind2 + kk;
                          
                          fact[0] += val * soln(ndf*ind1);
                          fact[1] += val * soln(ndf*ind1+1);
                       }
                    }

                    outp2[2].push_back(fact[0]);
                    outp2[3].push_back(fact[1]);

                    param[0] = 0.5;
                    param[1] = vv[jj];
                    ComputeGeometry(param, geom);

                    yy.push_back(geom[1]);
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
*/
   return;
}



bool HBSplineBase::converged()
{
  char fct[] = "HBSplineBase::converged";

  if (rNorm < tol && localStiffnessError == 0)
  {
    //if(STAGGERED && (SOLID_SOLVED == 0))
    //{
      //solveSolidProblem();
      //SOLID_SOLVED = 1;
    //}
    return true;
  }

  return false;
}




bool HBSplineBase::diverging(double factor)
{
  if (rNormPrev > -0.1 && (rNorm / rNormPrev) > factor) return true;

  if (localStiffnessError != 0) return true;

  if (prgNAN(rNorm)) return true;

  return false;
}


void HBSplineBase::printComputerTime(bool reset, int detailFlg)
{
  COUT << "----------------------------------------------------\n";

  COUT; printf("HBSplineBase::calcStiffnessRes:  %7.3f sec ->%5.1f %\n",
               ctimCalcStiffRes, ctimCalcStiffRes/ctimSinceLastCall*100.);

  COUT; printf("HBSplineBase::factoriseSolvUpdt: %7.3f sec ->%5.1f %\n",
               ctimFactSolvUpdt, ctimFactSolvUpdt/ctimSinceLastCall*100.);

  if (reset)
  {
    ctimFactSolvUpdt = 0.;
    ctimCalcStiffRes = 0.;
  }

  return;
}


void HBSplineBase::printData(int index, int patchnum)
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

  if(index == 8 || index == 9) // write matrix pattern to an output file 
  {
    if(SOLVER_TYPE == SOLVER_TYPE_EIGEN)
      solverEigen->printMatrixPatternToFile();
    //else
      //solverPetsc->printMatrixPatternToFile();
  }

  return;
}




