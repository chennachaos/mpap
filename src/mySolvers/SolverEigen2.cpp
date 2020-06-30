
#include "Debug.h"
#include "SolverEigen.h"
#include "SolverTime.h"
#include "ComputerTime.h"
#include "util.h"


#include <Eigen/SuperLUSupport>
#include <Eigen/SparseExtra>
#include <Eigen/IterativeSolvers>


extern SolverTime      solverTime;
extern ComputerTime    computerTime;

using namespace std;
using namespace Eigen;



void SolverEigen::computeConditionNumber()
{
/*
    MatrixXd  globalK;
    
    globalK.resize(nRow, nCol);
    globalK.setZero();

    int k, ii, jj;

    for(k=0; k<mtx.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
      {
        ii = it.row();
        jj = it.col();
        
        //cout << ii << '\t' << jj << '\t' << it.value() << endl;

        globalK.coeffRef(ii, jj) = it.value();
      }
    }


    VectorXd sing_vals = globalK.jacobiSvd().singularValues();

    //printf("\n Matrix condition number = %12.6E \n", sing_vals(0)/sing_vals(sing_vals.size()-1) );
    
    //printf("\n Minimum eigenvalue = %12.6f \n", sing_vals.minCoeff() );
    //printf("\n Minimum eigenvalue = %12.6f \n", sing_vals.maxCoeff() );
    printf("\n Matrix condition number = %12.6E \n", sing_vals.maxCoeff() / sing_vals.minCoeff() );
    printf("\n\n\n\n");
*/

    //myMatrix::OneNormEst estimator(mtx.rows(), 4);
    //double norm_a;
    
    //estimator.ANorm(mtx, mtx, norm_a);

    //cout << " norm_a = " << norm_a << endl;
//

  //VectorXd x0 = VectorXd::LinSpaced(nRow, 0.0, 1.0);

  //SuperLU<SparseMatrixXd>  solver;

  //double  cond = myCondNum(mtx, x0, 50, solver);

  //printf("\n Matrix condition number = %12.6E \n", cond);

  //myCondNumMatlab(mtx);

  return;
}




void SolverEigen::setupMatricesAndVectors()
{
    int nU, nP, nL, ii, jj, niter, k;
    double alpha, TOL, rho, beta, fact, f3norm;
    bool preCond;

    nP = nRow/3;
    nU = nP*2;

    //nU = 10918;
    //nP = 5459;
    //nL = 404;
    //nL = 0;
  
    cout << " nU and nP " << nU << '\t' << nP << '\t' << nL << endl;

    SparseMatrixXd  tempMat(nU, nU);

    matA.resize(nU,nU);
    matA.reserve(mtx.nonZeros()*2/3);
    matB.resize(nP,nU);
    matB.reserve(mtx.nonZeros()/3);
    matE.resize(nL,nU);
    matE.reserve(mtx.nonZeros()/3);

    PreCondSchur.resize(nP, nP);
    PreCondSchur.reserve(nP*nP);
    tempMat.reserve(nU);

    for(k=0; k<mtx.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
      {
        ii = it.row();
        jj = it.col();
        if( ii >= (nU+nP) )
        {
          if( jj < nU )
            matE.coeffRef(ii-nU-nP, jj) = it.value();
        }
        else if(ii >= nU && ii < (nU+nP) )
        {
          if( jj < nU )
            matB.coeffRef(ii-nU, jj) = it.value();
        }
        else
        {
          if( jj < nU )
            matA.coeffRef(ii, jj) = it.value();
        }
      }
    }

    //for(ii=0;ii<nU;ii++)
      //tempMat.coeffRef(ii,ii) = 1.0/matA.coeffRef(ii,ii);

    //PreCondSchur = matB*tempMat*matB.transpose();
    PreCondSchur = matB*matB.transpose();

    /*
    for(ii=0;ii<nP;ii++)
      PreCondSchur.coeffRef(ii,ii) = 0.0;

    for(k=0; k<matB.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(matB, k); it; ++it)
      {
        ii   = it.row();
        jj   = it.col();
        fact = it.value();

        PreCondSchur.coeffRef(ii,ii) += (fact*fact/matA.coeffRef(jj,jj));
      }
    }

    for(k=0; k<PreCondSchur.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(PreCondSchur, k); it; ++it)
        cout << it.row() << '\t' << it.col() << '\t' << it.value() << endl;
    }
    */

    if(STABILISED)
    {
      matC.resize(nP, nP);
      for(ii=0;ii<nP;ii++)
      {
        for(jj=0;jj<nP;jj++)
          matC.coeffRef(ii,jj) = mtx.coeffRef(nU+ii, nU+jj);
      }
    }

    f1.resize(nU);
    f2.resize(nP);
    f3.resize(nL);
    f1.setZero();
    f2.setZero();
    f3.setZero();

    for(ii=0;ii<nU;ii++)
      f1(ii) = rhsVec(ii);

    for(ii=0;ii<nP;ii++)
      f2(ii) = rhsVec(nU+ii);

    for(ii=0;ii<nL;ii++)
      f3(ii) = rhsVec(nU+nP+ii);

    var1.resize(nU);
    var2.resize(nP);
    var3.resize(nL);
    var1.setZero();
    var2.setZero();
    var3.setZero();

    var1Prev = var1;
    var2Prev = var2;
    var3Prev = var3;

    return;
}

void SolverEigen::solverSchurCG()
{
    setupMatricesAndVectors();
    
    time_t tstart, tend;

    int nU, nP, nL, ii, jj, niter, k;
    double alpha, TOL, rho, rhoPrev, beta, fact;
    bool preCond;

    nP = nRow/3;
    nU = nP*2;
    nL = 0;
  
    VectorXd  tempVec, f4, r, Ap, z, p;

    cout << " nU and nP " << nU << '\t' << nP << endl;

    //BiCGSTAB<SparseMatrixXd> solver;
  
    BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver;

    //GMRES<SparseMatrixXd, IncompleteLUT<double> > solver;

    //ConjugateGradient<SparseMatrixXd>  solver;

    solver.compute(matA);

    solver.setMaxIterations(1000);
    solver.setTolerance(1.0e-16);

    GMRES<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //ConjugateGradient<SparseMatrixXd> solverPCSchur;
  
    //solver.preconditioner().setDroptol(1.0e-15);
    //solver.preconditioner().setFillfactor(3);

    //solver.set_restart(50);
    //solver.setEigenv(10);

    solverPCSchur.compute(PreCondSchur);

    solverPCSchur.setMaxIterations(1000);
    solverPCSchur.setTolerance(1.0e-16);

    char        tmp[200];
    MyString    tmpStr;

    preCond = false;
    preCond = true;
    TOL = 1.0e-8;
    niter = 300;

    //for(ii=0;ii<nU;ii++)
      //var1(ii) = solnPrev(ii) ;

    //for(ii=0;ii<nP;ii++)
      //var2(ii) = solnPrev(nU+ii);

    tstart = time(0);

    //f4 = matB*solver.solve(f1) - f2;

    normRef = f4.norm();

    cout << f1.norm() << '\t' << f2.norm() << '\t' << normRef << endl;

    //r = f4 - matB*solver.solve(matB.transpose()*var2);

    if(preCond)
      z = solverPCSchur.solve(r);
    else
      z = r;
    
    p = z;

    for(ii=0;ii<niter;ii++)
    {
      rho = r.dot(z);

      //Ap = matB*solver.solve(matB.transpose()*p);
      //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;

      alpha = rho / (p.dot(Ap) );
      
      var2  += alpha * p; // update the solution variable (var2 here)
      r  -= alpha * Ap;   // update the residual

      fact = r.norm()/normRef;
      
      tmpStr.free();
      sprintf(tmp," \t %6d \t %12.6E \t %12.6E \t %12.6E \t %12.6E \t %12.6E \n", ii, fact, p.norm(), rho, alpha, beta);
      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);

      //printf("%5d \t %12.6f \n", ii, var2.norm());

      if( fact < TOL )
        break;
      
      if(preCond)
        z = solverPCSchur.solve(r);
      else
        z = r;

      beta = z.dot(r)/rho;
      p = z + beta*p;
    }

    cout << " Number of iterations = " << (ii+1) << endl;

    var1 = solver.solve(f1 - matB.transpose()*var2);
    cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;

    tend = time(0); 
    printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );
    
    var1Prev = var1;
    var2Prev = var2;
    var3Prev = var3;

    for(ii=0;ii<nU;ii++)
      soln(ii) = var1(ii);
  
    for(ii=0;ii<nP;ii++)
      soln(nU+ii) = var2(ii);

    for(ii=0;ii<nL;ii++)
      soln(nU+nP+ii) = var3(ii);

  return;
}


void SolverEigen::solverSchurBiCGSTAB()
{
    setupMatricesAndVectors();

    time_t tstart, tend;

    int nU, nP, nL, ii, jj, niter, k;
    double alpha, TOL, rho, beta, fact, omega, rhoPrev, error;
    bool preCond;

    nP = nRow/3;
    nU = nP*2;
    nL = 0;

    cout << " nU and nP " << nU << '\t' << nP << endl;

    VectorXd  tempVec, f4, r, rtilde, v, t, s, shat, z, p, phat;

    //BiCGSTAB<SparseMatrixXd> solver;
  
    BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver;

    //GMRES<SparseMatrixXd, IncompleteLUT<double> > solver;

    //ConjugateGradient<SparseMatrixXd>  solver;

    solver.compute(matA);

    solver.setMaxIterations(1000);
    solver.setTolerance(1.0e-16);

    GMRES<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //ConjugateGradient<SparseMatrixXd> solverPCSchur;
  
    //solver.preconditioner().setDroptol(1.0e-15);
    //solver.preconditioner().setFillfactor(3);

    //solver.set_restart(50);
    //solver.setEigenv(10);

    solverPCSchur.compute(PreCondSchur);

    solverPCSchur.setMaxIterations(1000);
    solverPCSchur.setTolerance(1.0e-16);

    char        tmp[200];
    MyString    tmpStr;

    preCond = false;
    //preCond = true;
    TOL = 1.0e-6;
    niter = 35;

    //for(ii=0;ii<nU;ii++)
      //var1(ii) = solnPrev(ii) ;

    //for(ii=0;ii<nP;ii++)
      //var2(ii) = solnPrev(nU+ii);

    p.resize(nP);
    p.setZero();
    v = p;
    t = p;

    tstart = time(0);

    ////f4 = matB*solver.solve(f1) - f2;

   // r = f4 - matB*solver.solve(matB.transpose()*var2);

    rho = rhoPrev = alpha = beta = omega = 1.0;

    normRef = f4.norm();

    rtilde = r;

    cout << f1.norm() << '\t' << f2.norm() << '\t' << f4.norm() << endl;

    for(ii=0;ii<niter;ii++)
    {
      rho = rtilde.dot(r);
      
      if(CompareDoubles(rho, 0.0))
        break;

      beta = (rho/rhoPrev)/(alpha/omega);
      p = r + beta*(p-omega*v);

      if(preCond)
        phat = solverPCSchur.solve(p);
      else
        phat = p;

      //v = matB*solver.solve(matB.transpose()*phat);
      //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;

      alpha = rho / (rtilde.dot(v) );
      
      s = r - alpha*v;

      if(s.norm() < TOL)
      {
        var2 = var2 + alpha*phat;
        break;
      }

      if(preCond)
        shat = solverPCSchur.solve(s);
      else
        shat = s;

      //t = matB*solver.solve(matB.transpose()*shat);
      //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;

      //printf("%5d \t %12.6f \t %12.6f \t %12.6f \t %12.6f \n", ii, s.norm(), t.norm(), p.norm(), var2.norm());
      printf("%5d \t %12.6E \t %12.6f \t %12.6f \t %12.6f \n", ii, rho, alpha, beta, omega);

      omega = t.dot(s)/t.dot(t);

      if( omega == 0.0 )
        break;

      //cout << rho << '\t' << alpha << '\t' << omega << endl;

      var2 = var2 + alpha*phat + omega*shat;
      r = s - omega*t;
      error = r.norm()/normRef;

      tmpStr.free();
      sprintf(tmp," \t %6d \t %12.6E \n", ii, error);
      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);
      
      if( error < TOL )
        break;

      rhoPrev = rho;
    }

    int flag;

    if( error <= TOL || s.norm() <= TOL )  // converged
    {
      flag =  0;
      if ( s.norm() <= TOL )
        error = s.norm() / normRef;
    }
    else if( omega == 0.0 ) // breakdown
      flag = -2;
    else if( rho == 0.0 )
      flag = -1;
    else        // no convergence
      flag = 1;

    cout << " Number of iterations = " << (ii+1) << endl;
    cout << " flag " << flag << endl;

    var1 = solver.solve(f1 - matB.transpose()*var2);
    //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;

    tend = time(0); 
    printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );
    
    var1Prev = var1;
    var2Prev = var2;
    var3Prev = var3;

    for(ii=0;ii<nU;ii++)
      soln(ii) = var1(ii);
  
    for(ii=0;ii<nP;ii++)
      soln(nU+ii) = var2(ii);

    for(ii=0;ii<nL;ii++)
      soln(nU+nP+ii) = var3(ii);

  return;
}


void SolverEigen::solverSchurGMRES()
{
    setupMatricesAndVectors();

    time_t tstart, tend;

    int nU, nP, nL, ii, jj, niter, k, m;
    double alpha, TOL, rho, beta, fact;
    bool preCond;

    nP = nRow/3;
    nU = nP*2;
    nL = 0;

    VectorXd  tempVec, f4, r, r1, r2, Ap, Ap1, Ap2, z, z1, z2, p, p1, p2;

    cout << " nU and nP " << nU << '\t' << nP << endl;

    //BiCGSTAB<SparseMatrixXd> solver;
  
    BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver;

    //GMRES<SparseMatrixXd, IncompleteLUT<double> > solver;

    //ConjugateGradient<SparseMatrixXd>  solver;

    solver.compute(matA);

    solver.setMaxIterations(1000);
    solver.setTolerance(1.0e-16);

    GMRES<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //ConjugateGradient<SparseMatrixXd> solverPCSchur;
  
    //solver.preconditioner().setDroptol(1.0e-15);
    //solver.preconditioner().setFillfactor(3);

    //solver.set_restart(50);
    //solver.setEigenv(10);

    solverPCSchur.compute(PreCondSchur);

    solverPCSchur.setMaxIterations(1000);
    solverPCSchur.setTolerance(1.0e-16);

    char        tmp[200];
    MyString    tmpStr;

    preCond = false;
    //preCond = true;
    TOL = 1.0e-7;
    niter = 400;

    //for(ii=0;ii<nU;ii++)
      //var1(ii) = solnPrev(ii) ;

    //for(ii=0;ii<nP;ii++)
      //var2(ii) = solnPrev(nU+ii);

    tstart = time(0);

    ////f4 = matB*solver.solve(f1) - f2;

    //r = f4 - matB*solver.solve(matB.transpose()*var2);
    
    rho = alpha = beta = 1.0;

    p.resize(nP);
    p.setZero();

    normRef = f4.norm();

    cout << f1.norm() << '\t' << f2.norm() << '\t' << f4.norm() << endl;
    
    m = 30;
    VectorXd  s(m+1), cs(m+1), sn(m+1);
    double  w;
    MatrixXd  v(nP, m+1);

    for(ii=0;ii<niter;ii++)
    {
      Ap1 = matA*p1 + matB.transpose()*p2;
      Ap2 = matB*p1;

      rho = r1.dot(z1) + r2.dot(z2);

      alpha = rho / (p1.dot(Ap1)+p2.dot(Ap2));

      var1  += alpha * p1;
      var2  += alpha * p2; // update the solution variable (var2 here)

      r1 -= alpha * Ap1;   // update the residual
      r2 -= alpha * Ap2;

      fact = var1.norm();

      tmpStr.free();
      sprintf(tmp," \t %6d \t %12.6E \n", ii, fact);
      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);

      //printf("%5d \t %12.6f \n", k, fact);

      if( fact < TOL )
        break;

      z1 = r1;
      //z2 = -matB*solver.solve(r1) + r2;

      beta = (z1.dot(r1)+z2.dot(r2))/rho;

      p1 = z1 + beta*p1;
      p2 = z2 + beta*p2;
    }

    cout << " Number of iterations = " << (ii+1) << endl;

    tend = time(0); 
    printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );
    
    var1Prev = var1;
    var2Prev = var2;
    var3Prev = var3;

    for(ii=0;ii<nU;ii++)
      soln(ii) = var1(ii);
  
    for(ii=0;ii<nP;ii++)
      soln(nU+ii) = var2(ii);

    for(ii=0;ii<nL;ii++)
      soln(nU+nP+ii) = var3(ii);

  return;
}

  
void SolverEigen::solverUzawaType1()
{
    setupMatricesAndVectors();

    time_t tstart, tend;

    int nU, nP, nL, ii, jj, niter, k;
    double alpha, TOL, rho, beta, normPrev, error;
    bool preCond;

    nP = nRow/3;
    nU = nP*2;
    
    //nU = 10918;
    //nP = 5459;
    //nL = 404;
    //nL = 0;

    VectorXd  tempVec, f4, r2, r2prev, z, q;

    cout << " nU and nP " << nU << '\t' << nP << endl;

    //BiCGSTAB<SparseMatrixXd> solver;
  
    BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver;

    //GMRES<SparseMatrixXd, IncompleteLUT<double> > solver;

    //ConjugateGradient<SparseMatrixXd>  solver;

    solver.compute(matA);

    solver.setMaxIterations(1000);
    solver.setTolerance(1.0e-16);

    GMRES<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //ConjugateGradient<SparseMatrixXd> solverPCSchur;
  
    //solver.preconditioner().setDroptol(1.0e-15);
    //solver.preconditioner().setFillfactor(3);

    //solver.set_restart(50);
    //solver.setEigenv(10);

    solverPCSchur.compute(PreCondSchur);

    solverPCSchur.setMaxIterations(1000);
    solverPCSchur.setTolerance(1.0e-10);

    char        tmp[200];
    MyString    tmpStr;

    //for(ii=0;ii<nU;ii++)
      //var1(ii) = solnPrev(ii) ;

    //for(ii=0;ii<nP;ii++)
      //var2(ii) = solnPrev(nU+ii);

    cout << f1.norm() << '\t' << f2.norm() << '\t' << f3.norm() << endl;

    preCond = false;
    //preCond = true;
    TOL = 1.0e-5;
    niter = 200;
    alpha = 30.0;
    beta  = 0.001;

    tstart = time(0);
  
/*
    r2 = f1 - matB.transpose()*var2 - matE.transpose()*var3;
    normPrev =  f2.norm() + 100.0;

    for(ii=0;ii<niter;ii++)
    {
      printf("%5d \t %12.6f \t %12.6f \t %12.6f \n", ii, var1.norm(), var2.norm(), var3.norm());

      var1 = solver.solveWithGuess(r2, var1Prev);
      //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;
      //var2 += alpha * solverPCSchur.solve(matB*var1 - f2);
      var2 += alpha * (matB*var1 - f2);
      var3 += beta  * (matE*var1 - f3);

      //r2 = f1 - matB.transpose()*var2 - matE.transpose()*var3;
      
      error = r2.norm();

      tempVec = var1 - var1Prev;
      //printf("%5d \t %12.6f \t %12.6f \n", k, var1.norm(), var2.norm());

      tmpStr.free();
      //sprintf(tmp," \t %6d \t %12.6E \n", ii, tempVec.norm());
      sprintf(tmp," \t %6d \t %12.6E \t %12.6E \n", ii, error, tempVec.norm());
      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);

      if( error < TOL )
        break;

      var1Prev = var1;
    }
*/
//
    for(ii=0;ii<niter;ii++)
    {
     // //var1 = solver.solve(f1 - matB.transpose()*var2);
      //var1 = solver.solveWithGuess(f1 - matB.transpose()*var2, var1Prev);
      cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;
      ////var2 += alpha * (matB*var1 - matC*var2);
      //var2 += alpha * (matB*var1 - f2);
      ////var2 += alpha * solverPCSchur.solve(matB*var1 - f2);

      tempVec = var1 - var1Prev;
      //printf("%5d \t %12.6f \t %12.6f \n", k, var1.norm(), var2.norm());

      tmpStr.free();
      sprintf(tmp," \t %6d \t %12.6E \n", ii, tempVec.norm());
      //sprintf(tmp," \t %6d \t %12.6E \n", ii, var1.norm());
      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);

      //if( tempVec.norm() < TOL )
        //break;

      var1Prev = var1;
      //printf("%5d \t %12.6f \t %12.6f \n", ii, u.norm(), p.norm());
    }
//

    VectorXd  var1b(nU), var1c(nU), var2b(nP), var2c(nP);
    var1b.setZero();
    var1c.setZero();

    var2b.setZero();
    var2c.setZero();

/*
    // with Aitken acceleration

    var1 = solver.solve(f1 - matB.transpose()*var2);
    var2 += alpha * (matB*var1 - f2);

    for(ii=0;ii<niter;ii++)
    {
      var1b = solver.solve(f1 - matB.transpose()*var2);
      //var1b = solver.solveWithGuess(f1 - matB.transpose()*var2, var1Prev);
      //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;
      //var2 += alpha * (matB*var1 - matC*var2);
      var2b += alpha * (matB*var1b - f2);
      //var2 += alpha * solverPCSchur.solve(matB*var1 - f2);

      //var1c = solver.solveWithGuess(f1 - matB.transpose()*var2b, var1Prev);
      var1c = solver.solve(f1 - matB.transpose()*var2b);
      var2c += alpha * (matB*var1c - f2);

      for(k=0;k<nU;k++)
        var1(k) = var1(k) - (var1b(k)-var1(k))*(var1b(k)-var1(k))/(var1c(k)-2.0*var1b(k)+var1(k));
      //var1 = var1c;

      for(k=0;k<nP;k++)
        var2(k) = var2(k) - (var2b(k)-var2(k))*(var2b(k)-var2(k))/(var2c(k)-2.0*var2b(k)+var2(k));

      //tempVec = var1 - var1Prev;
      //printf("%5d \t %12.6f \t %12.6f \n", k, var1.norm(), var2.norm());

      tmpStr.free();
      //sprintf(tmp," \t %6d \t %12.6E \n", ii, tempVec.norm());
      sprintf(tmp," \t %6d \t %12.6E \n", ii, var1c.norm());
      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);

      //if( tempVec.norm() < TOL )
        //break;

      var1Prev = var1;
      //printf("%5d \t %12.6f \t %12.6f \n", ii, u.norm(), p.norm());
    }
*/
    //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;
/*
    TOL = 1.0e-8;
    niter = 1000;
    alpha = 20.0;

    f4 = matB*solver.solve(f1) - f2;
    
    r2 = f4 - matB*solver.solve(matB.transpose()*var2);
    r2prev = r2;
    
    normPrev = r2.norm() + 10.0;

    for(ii=0;ii<niter;ii++)
    {
      //z = solverPCSchur.solve(r2);
      z = r2;

      var2 += alpha * z;

      r2 = f4 - matB*solver.solve(matB.transpose()*var2);
//
      if(r2.norm() < normPrev)
      {
        //r2prev = r2;
        alpha *= 1.1;
        normPrev = r2.norm();
      }
      else
      {
        r2 = r2prev;
        alpha *= 0.5;
      }
      //printf(" \t %6d \t %12.6E \n", ii, alpha);
//
//
      if(r2.norm() > 1.0e-4)
        var2 += alpha * r2;
      else
      {
        var2b = var2  + alpha * r2;
        var2c = var2b + alpha * (f4 - matB*solver.solve(matB.transpose()*var2b));

        for(k=0;k<nP;k++)
          var2(k) = var2(k) - (var2b(k)-var2(k))*(var2b(k)-var2(k))/(var2c(k)-2.0*var2b(k)+var2(k));
      }
//

      //printf("%5d \t %12.6f \t %12.6f \n", ii, var2b.norm(), var2c.norm());

      tmpStr.free();
      sprintf(tmp," \t %6d \t %12.6E \t %12.6E \n", ii, r2.norm(), alpha);
      //sprintf(tmp," \t %6d \t %12.6E \n", ii, r2.norm());
      //sprintf(tmp," \t %6d \t %12.6E \n", ii, var2.norm());
      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);

      //if( tempVec.norm() < TOL )
        //break;
    }

    var1 = solver.solve(f1 - matB.transpose()*var2);
    cout << " Number of iterations = " << (ii+1) << endl;

*/

    tend = time(0); 
    printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );

    for(ii=0;ii<nU;ii++)
      soln(ii) = var1(ii);

    for(ii=0;ii<nP;ii++)
      soln(nU+ii) = var2(ii);

    for(ii=0;ii<nL;ii++)
      soln(nU+nP+ii) = var3(ii);

    /*
    //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver2;
    GMRES<SparseMatrixXd, IncompleteLUT<double> > solver2;

    solver2.compute(mtx);
    solver2.setMaxIterations(100);
    solver2.setTolerance(1.0e-16);

    //solver2.preconditioner().setDroptol(1.0e-15);
    solver2.preconditioner().setFillfactor(4);

    soln = solver2.solveWithGuess(rhsVec, soln);
    */

  return;
}


void SolverEigen::solverUzawaType2()
{

  return;
}

int  SolverEigen::myBiCGSTAB()
{
    MatrixXd  globalK;
/*
    globalK.resize(nRow, nCol);
    globalK.setZero();

    int k, ii, jj;

    for(k=0; k<mtx.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
      {
        ii = it.row();
        jj = it.col();
        
        //cout << ii << '\t' << jj << '\t' << it.value() << endl;

        globalK.coeffRef(ii, jj) = it.value();
      }
    }
    SparseLU<SparseMatrixXd>  slu;
    slu.compute(mtx);
    slu.factorize(mtx);
    cout <<  slu.matrixL().nonZeros() << endl;
    //SparseMatrixXd  sL = slu.matrixL();
    //SparseMatrixXd  sU = slu.matrixU();
*/


  using std::sqrt;
  using std::abs;

  typedef double  RealScalar;
  typedef double  Scalar;
  typedef Matrix<Scalar,Dynamic,1>  VectorType;
  
  typedef  int Index;
  
  RealScalar  tol_error = 1.0e-8;
  Index iters = 500;
  Index iterTemp;
  
  RealScalar tol = tol_error;
  RealScalar tolTemp;
  Index maxIters = iters;

  
  Scalar rho    = 1.0;
  Scalar alpha  = 1.0;
  Scalar w      = 1.0;
  Scalar omega  = 0.1;

  //soln.setZero();
  soln = solnPrev;
  
  //for(Index ii=0; ii<nRow/3; ii++)
    //soln(3*ii) = 0.1;

  iterTemp = 100;
  tolTemp = 1.0e-3;
  Index  nTerms=20;
  //mySOR_Sparse(mtx, rhsVec, soln, iterTemp, tolTemp, omega);
  //myPolyprecond_Dense(globalK, rhsVec, soln, iterTemp, tolTemp, nTerms);
  //cout << iterTemp << '\t' << tolTemp << endl;


  Index n = mtx.cols();
  VectorType r  = rhsVec - mtx * soln;
  //VectorType r  = rhsVec ;
  VectorType r0 = r;
  
  RealScalar r0_sqnorm = r0.squaredNorm();
  RealScalar rhs_sqnorm = rhsVec.squaredNorm();

  if(rhs_sqnorm == 0.0)
  {
    soln.setZero();
    return true;
  }


  VectorType v = VectorType::Zero(n);
  VectorType p = VectorType::Zero(n);
  VectorType y = VectorType::Zero(n);
  VectorType z = VectorType::Zero(n);

  VectorType kt = VectorType::Zero(n);
  VectorType ks = VectorType::Zero(n);
  VectorType s = VectorType::Zero(n);
  VectorType t = VectorType::Zero(n);


  //cout << NumTraits<Scalar>::epsilon() << endl;

  RealScalar tol2 = tol*tol*rhs_sqnorm;
  RealScalar eps2 = NumTraits<Scalar>::epsilon()*NumTraits<Scalar>::epsilon();
  Index i = 0;
  Index restarts = 0;
  
  //cout << " Maximum coefficient = " << mtx.max() << endl;
  //cout << " Minimum coefficient = " << mtx.minCoeff() << endl;
  
  //if(update_precond <= 3)
    updatePreconditioner();
  
  //cout << " ooooooooooooooo  " << endl;

  //myIncompleteLUT2<double>  precond;

  //precond.setDroptol(1.0e-3);
  //precond.setFillfactor(10);
  //precond.compute(mtx);

  while( (r.squaredNorm() > tol2) && (i<maxIters) )
  {
    Scalar rho_old = rho;

    rho = r0.dot(r);
    if (abs(rho) < eps2*r0_sqnorm)
    {
      // The new residual vector became too orthogonal to the arbitrarily chosen direction r0
      // Let's restart with a new r0:
      r  = rhsVec - mtx * soln;
      r0 = r;
      rho = r0_sqnorm = r.squaredNorm();
      if(restarts++ == 0)
        i = 0;
    }
    Scalar beta = (rho/rho_old) * (alpha / w);
    p = r + beta * (p - w * v);


    //y = p;
    //iterTemp = 1;
    //tolTemp = 1.0e-4;
    //myJacobi_Sparse(mtx, p, y, iterTemp, tolTemp, omega);
    //mySOR_Sparse(mtx, p, y, iterTemp, tolTemp, omega);
    //myPolyprecond_Sparse(mtx, p, y, iterTemp, tolTemp, nTerms);
    //cout << iterTemp << '\t' << tolTemp << endl;

    v.noalias() = mtx * y;

    alpha = rho / r0.dot(v);
    s = r - alpha * v;

    //z = s;
    //iterTemp = 1;
    //tolTemp = 1.0e-4;
    //myJacobi_Sparse(mtx, s, z, iterTemp, tolTemp, omega);
    //mySOR_Sparse(mtx, s, z, iterTemp, tolTemp, omega);
    //myPolyprecond_Sparse(mtx, s, z, iterTemp, tolTemp, nTerms);
    //cout << iterTemp << '\t' << tolTemp << endl;
    t.noalias() = mtx * z;


    RealScalar tmp = t.squaredNorm();
    if( tmp > RealScalar(0.0) )
      w = t.dot(s) / tmp;
    else
      w = Scalar(0);
    soln += alpha * y + w * z;
    r = s - w * t;
    ++i;
  }

  tol_error = sqrt(r.squaredNorm()/rhs_sqnorm);

  iters = i;

  cout << " Iterative solver stats ... Tol = " << tol_error << "\t Niter = " << iters << endl;
  
  //soln = soln - solnPrev;
  soln;

  return 1;
}




void  SolverEigen::updatePreconditioner()
{
  ////////////////////// 
  //
  //  preconditioner
  //
  ////////////////////// 
 
  //precond1.setDroptol(1.0e-3);
  //precond1.setFillfactor(2);
  //precond1.compute(mtx);

    time_t tstart, tend;
    tstart = time(0);

  tend = time(0); 
  printf("ILUT preconditioner took %8.4f second(s) \n ", difftime(tend, tstart) );

  //mtx.makeCompressed();

  //DiagonalPreconditioner<double>  precond(mtx);

  update_precond++;

  return;

}






