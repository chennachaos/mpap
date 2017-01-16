
#ifndef MY_STABILISATION_ROUTINES_H
#define MY_STABILISATION_ROUTINES_H








/****
c============================================================================

c----------------------------------------------------------------------------
c  stabilisation parameters tau_u, tau_p and tau_c
c
c                         beta1i
c  zi(Re) = ---------------------------------- ,  i=1,2,3
c            sqrt(1 + (beta1i/(beta2i Re))^2)
c
c             h
c  tau_u = ------- z1(Re)
c           2 |u|
c
c               h
c  tau_p = ----------- z2(Re)
c           2 |u| rho
c
c  tau_c = h |u| z3(Re)
c
c----------------------------------------------------------------------------
**/

inline  void  evaluateStabParams_algo1(double* u, double h, double rho, double mu, double dt,  double* beta, double* tau)
{
  double  fact, fact2;

  double  uNorm = u[0]*u[0];
  for(int i=1; i<3; i++)
    uNorm += u[i]*u[i] ;

  uNorm  = sqrt(uNorm);    //   |u|

  double  Re  = rho*uNorm*h/(2.0*mu);   //  Re(u)  Reynolds number
  
  //double  C   = uNorm*dt/h ;                //  C(u)   Courant no.
  //double  nRC[3];
  //nRC[0]  = uNorm ;  nRC[1]  = Re ;  nRC[2]  = C ;

  if( Re < 1.0e-10 )
  {
    fact   = h*h/(4.0*mu/rho) ;
    tau[0] = fact * beta[0] ;      // tau_u
    tau[1] = fact * beta[2] ;      // tau_p
    tau[2] = h * uNorm * Re * beta[4] ;  // tau_c
  }
  else
  {
    fact2 = h/(2.0*uNorm*rho);

    //  tau_u
    fact = beta[0]/(beta[1]*Re);
    tau[0] = fact2*beta[0]/sqrt(1.0 + fact*fact) ;

    //  tau_p
    fact = beta[2]/(beta[3]*Re);
    tau[1] = fact2*beta[2]/sqrt(1.0 + fact*fact) ;

    //  tau_c
    fact2 = h*uNorm;
    fact = beta[4]/(beta[5]*Re);
    tau[2] = fact2*beta[4]/sqrt(1.0 + fact*fact) ;
  }

  return;
}



inline void  evaluateStabParams_algo2(double* u, double h, double rho, double mu, double dt,  double* beta, double* tau)
{
  double  uNorm = u[0]*u[0];
  for(int i=1; i<3; i++)
    uNorm += u[i]*u[i] ;

  uNorm  = sqrt(uNorm);    //   |u|

  double  t1 = 2.0/dt;
  double  t2 = 2.0*uNorm/h;
  double  t3 = 4.0*mu/rho/h/h;

  //  tau_u
  //tau[0] = t1*t1 + t2*t2 + t3*t3;

  tau[0] = t2*t2 + t3*t3;

  tau[0] = 1.0/sqrt( tau[0] );

  //  tau_c
  tau[2] = 1.0/(tau[0]*8.0/h/h);

  //tau[0] = tau[0]/rho ;

  tau[1] = h*h/12.0/mu;

  return;
}



inline void  evaluateStabParams_algo3(VectorXd& u, MatrixXd& G, double dt, double rho, double mu, double CI, double* tau)
{
  //tau[0] = 4.0/dt/dt + u.dot(G*u) + CI*(mu*mu/rho/rho)*G.cwiseAbs2().sum();

  tau[0] = u.dot(G*u) + CI*(mu*mu/rho/rho)*G.cwiseAbs2().sum();

  tau[0] = 1.0/sqrt(tau[0]); // tau_M

  tau[2] = 1.0/(tau[0]*G.trace());  // tau_C

  tau[0] = tau[0]/rho;

  tau[1] = tau[0];

  return;
}



#endif



