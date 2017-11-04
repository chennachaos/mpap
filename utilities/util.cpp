
#include "headersBasic.h"
#include "util.h"
#include "myConstants.h"

using namespace std;
using namespace Eigen;




bool doubleGreater(double left, double right, bool orequal)
{
  if (fabs(left - right) < DBL_EPSILON) 
    return (orequal);

  return (left > right);
}


bool doubleLess(double left, double right, bool orequal)
{
  if (fabs(left - right) < DBL_EPSILON)
    return (orequal);

  return (left < right);
}




void GenerateCoeffMatrices(int deg, MatrixXd& SL, MatrixXd& SR)
{
    int ii, jj, ind, size = deg+1;
    
    VectorXd  xx(3*deg+2);
    xx.setZero();
        
    double  denom = pow(2.0, deg);

    for(ii=0;ii<=size;ii++)
      xx(deg+ii) = Bin(size,ii)/denom;
        
    ind = xx.rows() - deg - 2;
    
    for(ii=0;ii<size;ii++)
    {
       for(jj=0;jj<size;jj++)
       {
          SL(ii,jj)  =  xx(ind+jj);
          SR(ii,jj)  =  xx(ind+jj+1);
       }
       ind -= 2;
    }
    
    return;
}


void TensorProduct(MatrixXd& A, MatrixXd& B, MatrixXd& C)
{
    int ii, jj, r1, r2, c1, c2;
    
    r1 = A.rows();
    c1 = A.cols();

    r2 = B.rows();
    c2 = B.cols();
    
    C.resize(r1*r2, c1*c2);
    for(ii=0;ii<r1;ii++)
    {
      for(jj=0;jj<c1;jj++)
      {
        C.block(r2*ii,c2*jj,r2,c2) = A(ii,jj)*B;
      }
    }
    
    return;
}



void printMatrix(MatrixXd& AA)
{
    int ii, jj;
    printf("\n\n");
    for(ii=0;ii<AA.rows();ii++)
    {
        for(jj=0;jj<AA.cols();jj++)
           printf("\t%14.8f", AA(ii,jj));
        printf("\n");
    }
    printf("\n\n");

    return;
}



void printVector(VectorXd& AA)
{
    printf("\n\n");
    for(int ii=0;ii<AA.rows();ii++)
        printf("\t%6d\t%12.8f\n", ii, AA(ii));
    printf("\n\n");

   return;
}



void printVector(vector<int>&  vec)
{
    printf("\n\n");
    for(int ii=0;ii<vec.size();ii++)
        printf("\t%6d ", vec[ii]);
    printf("\n\n");

   return;
}


void printVector(double* data, int nn)
{
    printf("\n\n");
    for(int ii=0;ii<nn;ii++)
      printf("\t%6d\t%12.8f\n", ii, data[ii]);
    printf("\n\n");

   return;
}



void printVector(vector<double>&  vec)
{
    printf("\n\n");
    for(int ii=0;ii<vec.size();ii++)
        printf("\t%12.6f ", vec[ii]);
    printf("\n\n");

   return;
}



double  HeavisideFunction(double uu, double ee)
{
   double  val=0.0;
   if(uu < -ee)
     return 0.0;
   else if(abs(uu) <= ee)
     return  (0.5 + 0.5*tanh(uu/ee/ee));
     //return  (0.5 + 0.5*uu/ee - sin(PI*uu/ee)/PI);
   else
     return  1.0;
}



double  DiracDelta1(double r, double alpha)
{
   if(abs(r) <= alpha)
     return  (0.5*(1.0+cos(PI*r/alpha))/alpha);
   else
     return  0.0;
}


double  DiracDelta2(double r, double alpha)
{
    //if( abs(r) >= 2.0 )
      //return 0.0;

    if( (r >= -2.0) && (r <= -1.0) )
      return  ((5.0+2.0*r-sqrt(-7.0-12.0*r-4.0*r*r))/8.0);
    else if ( (r >= -1.0) && (r <= 0.0) )
      return  ((3.0+2.0*r+sqrt(1.0-4.0*r-4.0*r*r))/8.0);
    else if ( (r >= 0.0) && (r <= 1.0) )
      return  ((3.0-2.0*r+sqrt(1.0+4.0*r-4.0*r*r))/8.0);
    else if ( (r >= 1.0) && (r <= 2.0) )
      return  ((5.0-2.0*r-sqrt(-7.0+12.0*r-4.0*r*r))/8.0);
    else
      return 0.0;
}


double  IntegralDoubleDiracDelta1(double beta, double gamma)
{
   if(gamma < 0.0)
     return  0.0;
   else
   {
     double  fact = PI*gamma/beta;
     return  ((gamma+0.5*gamma*cos(fact) - (1.5*beta/PI)* sin(fact))/(4*beta*beta));
   }
}




void SetTimeParametersFluid(int tis, double rho, double dt, VectorXd& td)
{
  td.setZero();
  
  double  alpf, alpm, beta, gamm;

  td[0] = dt;

  switch(tis)
  {
      case  0: // quasi static

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

       break;

      case  1: // generalised midpoint rule

            alpf = 1.0/(1.0 + rho);
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = alpf*dt; //af*dt
            td[6]  = alpm/gamm;
            td[7]  = 1.0-alpm/gamm;
            td[8]  = 0.0;

            td[9]  = 1.0/dt;  // v_{n+1}
            td[10] = -td[9];  // v_n
            td[11] = 0.0;     // v_{n-1}
            td[12] = 0.0;     // v_{n-2}
            td[15] = 0.0;     // a_n

        break;

      case  2: // generalised alpha-method

            alpf = 1.0/(1.0 + rho);
            alpm = 0.5*(3.0 - rho)/(1.0 + rho);
            gamm = 0.5 + alpm - alpf;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = alpf*dt;
            td[6]  = alpm/gamm;
            td[7]  = 1.0-alpm/gamm;
            td[8]  = alpm/gamm/dt;

            td[9]  = 1.0/gamm/dt;     // v_{n+1}
            td[10] = -td[9];          // v_n
            td[11] = 0.0;     // v_{n-1}
            td[12] = 0.0;     // v_{n-2}
            td[15] = 1.0 - 1.0/gamm;  // a_n

         break;

      case  3: // Backward Euler or BDF1

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = 1.0;
            td[2]  = 1.0;
            td[3]  = 1.0;
            td[4]  = 1.0;

            td[5]  = dt;
            td[6]  = 1.0;
            td[7]  = 0.0;
            td[8]  = 1.0/dt;

            td[9]  =  1.0/dt;    // v_{n+1}
            td[10] = -1.0/dt;    // v_n
            td[11] =  0.0;       // v_{n-1}
            td[12] =  0.0;       // v_{n-2}
            td[13] =  0.0;       // v_{n-3}
            td[15] =  0.0;       // a_n

         break;

      case  4: // generalised alpha-method

            alpf = 1.0/(1.0 + rho);
            alpm = 0.5*(3.0 - rho)/(1.0 + rho);
            gamm = 0.5 + alpm - alpf;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = alpf*dt;
            td[6]  = alpm/gamm;
            td[7]  = 1.0-alpm/gamm;
            td[8]  = alpm/gamm/dt;

            td[9]  = 1.0/gamm/dt;    // v_{n+1}
            td[10] = -td[9];         // v_n
            td[11] = 0.0;            // v_{n-1}
            td[12] = 0.0;            // v_{n-2}
            td[13] =  0.0;           // v_{n-3}
            td[15] = 1.0 - 1.0/gamm; // a_n

         break;

      case  5: // BDF2

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = 1.0;
            td[2]  = 1.0;
            td[3]  = 1.0;
            td[4]  = 1.0;

            td[5]  = dt;
            td[6]  = 1.0;
            td[7]  = 0.0;
            td[8]  = 1.0/dt;

            td[9]  =  1.5/dt;    // v_{n+1}
            td[10] = -2.0/dt;    // v_n
            td[11] =  0.5/dt;    // v_{n-1}
            td[12] =  0.0;       // v_{n-2}
            td[13] =  0.0;       // v_{n-3}
            td[15] =  0.0;       // a_n

         break;

      case  6: // BDF3

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = 1.0;
            td[2]  = 1.0;
            td[3]  = 1.0;
            td[4]  = 1.0;

            td[5]  = dt;
            td[6]  = 1.0;
            td[7]  = 0.0;
            td[8]  = 1.0/dt;

            td[9]  =  (11.0/6.0)/dt;    // v_{n+1}
            td[10] = -3.0/dt;           // v_n
            td[11] =  1.5/dt;           // v_{n-1}
            td[12] = -(1.0/3.0)/dt;     // v_{n-2}
            td[13] =  0.0;              // v_{n-3}
            td[15] =  0.0;              // a_n

         break;

      case  7: // BDF4

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = 1.0;
            td[2]  = 1.0;
            td[3]  = 1.0;
            td[4]  = 1.0;

            td[5]  = dt;
            td[6]  = 1.0;
            td[7]  = 0.0;
            td[8]  = 1.0/dt;

            td[9]  =  (25.0/12.0)/dt;    // v_{n+1}
            td[10] = -4.0/dt;           // v_n
            td[11] =  3.0/dt;           // v_{n-1}
            td[12] = -(4.0/3.0)/dt;     // v_{n-2}
            td[13] =  (1.0/4.0)/dt;     // v_{n-3}
            td[15] =  0.0;              // a_n

         break;

      default:
            cerr << " SetTimeParametersFluid ... invalid value of tis!" << endl;

         break;
  }

  return;
}



void SetTimeParametersSolid(int tis, double rho, double dt, VectorXd& td)
{
  td.setZero();
  
  double  alpf, alpm, beta, gamm;

  td[0] = dt;

  switch(tis)
  {
      case  0: // quasi static

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;
            td[7]  = alpf;

            // velocity is used as the primary variable for 
            // solid dynamics problem 
            // td[10] is the multiplication when converting from
            // displacement based formulation to velocity based formulation
            // It is set to ONE for static problem

            td[10] = 1.0; 

      break;

      case  1: // generalised midpoint rule

            alpf = 1.0/(1.0 + rho);
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = 1.0/alpf/dt/dt;
            td[6]  = 1.0/dt;
            td[7]  = alpf;

            //displacement as the primary variable
            //v_{n+1}  = td[10]*d_{n+1} + td[11]*d_n + td[12]*v_n + td[13]*a_n + td[14]*ddot_n;
            //a_{n+1}  = td[15]*d_{n+1} + td[16]*d_n + td[17]*v_n + td[18]*a_n + td[19]*ddot_n;

            td[10] = 1.0/alpf/dt;      // d_{n+1}
            td[11] = -td[10];          // d_n
            td[12] = (alpf-1.0)/alpf;  // v_n
            td[13] = 0.0;              // a_n
            td[14] = 0.0;              // ddot_n

            td[15] = 1.0/alpf/dt/dt;   // d_{n+1}
            td[16] = -td[15];          // d_n
            td[17] = -1.0/alpf/dt;     // vn
            td[18] = 0.0;              // an
            td[19] = 0.0;              // ddot_n

            //velocity as the primary variable
            //d_{n+1}  = td[20]*v_{n+1} + td[21]*d_n + td[22]*v_n + td[23]*a_n + td[24]*ddot_n ;
            //a_{n+1}  = td[25]*v_{n+1} + td[26]*d_n + td[27]*v_n + td[28]*a_n + td[29]*ddot_n ;

            td[40] = alpf*dt;          // v_{n+1}
            td[41] = 1.0;              // d_n
            td[42] = 1.0-alpf;         // v_n
            td[43] = 0.0;              // a_n
            td[44] = 0.0;              // ddot_n

            td[45] = 1.0/dt;           // v_{n+1}
            td[46] = 0.0;              // d_n
            td[47] = -td[45];          // v_n
            td[48] = 0.0;              // a_n
            td[49] = 0.0;              // ddot_n

      break;

      case  2: // generalised alpha-method --- similar to fluid one

            alpm = (2.0-rho)/(rho+1.0);
            alpf = 1.0/(rho+1.0);

            gamm = 0.5 + alpm - alpf;
            beta = 0.25*(1.0+alpm-alpf)*(1.0+alpm-alpf);

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = alpm/beta/dt/dt;
            td[6]  = alpf*gamm/beta/dt;
            td[7]  = alpf;

            //displacement as the primary variable
            //v_{n+1}  = td[10]*d_{n+1} + td[11]*d_n + td[12]*v_n + td[13]*a_n + td[14]*ddot_n;
            //a_{n+1}  = td[15]*d_{n+1} + td[16]*d_n + td[17]*v_n + td[18]*a_n + td[19]*ddot_n;

            td[10] = gamm/beta/dt;           // d_{n+1}
            td[11] = -td[10];                // d_n
            td[12] = 1.0-gamm/beta;          // v_n
            td[13] = dt*(1.0-gamm/2.0/beta); // a_n
            td[14] = 0.0;                    // ddot_n

            td[15] = 1.0/beta/dt/dt;         // d_{n+1}
            td[16] = -td[15];                // d_n
            td[17] = -1.0/beta/dt;           // v_n
            td[18] = 1.0-1.0/2.0/beta;       // a_n
            td[19] = 0.0;                    // ddot_n

            //velocity as the primary variable
            //d_{n+1}  = td[20]*v_{n+1} + td[21]*d_n + td[22]*v_n + td[23]*a_n + td[24]*ddot_n ;
            //a_{n+1}  = td[25]*v_{n+1} + td[26]*d_n + td[27]*v_n + td[28]*a_n + td[29]*ddot_n ;
            
            td[40] = dt*beta/gamm;                     // v_{n+1}
            td[41] = 1.0;                              // d_n
            td[42] = dt*(gamm-beta)/gamm;              // v_n
            td[43] = dt*dt*(gamm-2.0*beta)/(2.0*gamm); // a_n
            td[44] = 0.0;                              // ddot_n

            td[45] = 1.0/(gamm*dt);                    // v_{n+1}
            td[46] = 0.0;                              // d_n
            td[47] = -td[45];                          // v_n
            td[48] = (gamm-1.0)/gamm;                  // a_n
            td[49] = 0.0;                              // ddot_n

      break;

      case  3: // Backward Euler

            alpf = 1.0;
            alpm = 1.0;
            gamm = 1.0;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = 1.0/dt/dt;
            td[6]  = 1.0/dt;
            td[7]  = alpf;

            //displacement as the primary variable
            //v_{n+1}  = td[10]*d_{n+1} + td[11]*d_n + td[12]*v_n + td[13]*a_n + td[14]*ddot_n;
            //a_{n+1}  = td[15]*d_{n+1} + td[16]*d_n + td[17]*v_n + td[18]*a_n + td[19]*ddot_n;

            td[10] = 1.0/dt;    // d_{n+1}
            td[11] = -td[10];   // d_n
            td[12] = 0.0;       // v_n
            td[13] = 0.0;       // a_n
            td[14] = 0.0;       // ddot_n

            td[15] = 1.0/dt/dt; // d_{n+1}
            td[16] = -td[15];   // d_n
            td[17] = -1.0/dt;   // vn
            td[18] = 0.0;       // an
            td[19] = 0.0;       // ddot_n

            //velocity as the primary variable
            //d_{n+1}  = td[20]*v_{n+1} + td[21]*d_n + td[22]*v_n + td[23]*a_n + td[24]*ddot_n ;
            //a_{n+1}  = td[25]*v_{n+1} + td[26]*d_n + td[27]*v_n + td[28]*a_n + td[29]*ddot_n ;

            td[40] = dt;        // v_{n+1}
            td[41] = 1.0;       // d_n
            td[42] = 0.0;       // v_n
            td[43] = 0.0;       // a_n
            td[44] = 0.0;       // ddot_n

            td[45] = 1.0/dt;    // v_{n+1}
            td[46] = 0.0;       // d_n
            td[47] = -td[45];   // v_n
            td[48] = 0.0;       // a_n
            td[49] = 0.0;       // ddot_n

      break;

      case  4: // generalised alpha-method --- state-space formulation

            alpf = 1.0/(1.0 + rho);
            alpm = 0.5*(3.0 - rho)/(1.0 + rho);
            gamm = 0.5 + alpm - alpf;

            td[1]  = alpm;
            td[2]  = alpf;
            td[3]  = alpm;
            td[4]  = gamm;

            td[5]  = (alpm*alpm)/(alpf*gamm*gamm*dt*dt);
            td[6]  = alpm/gamm/dt;
            td[7]  = alpf;

            //displacement as the primary variable
            //v_{n+1}    = td[10]*d_{n+1} + td[11]*d_n + td[12]*v_n + td[13]*a_n + td[14]*ddot_n;
            //a_{n+1}    = td[15]*d_{n+1} + td[16]*d_n + td[17]*v_n + td[18]*a_n + td[19]*ddot_n;
            //ddot_{n+1} = td[20]*d_{n+1} + td[21]*d_n + td[22]*v_n + td[23]*a_n + td[24]*ddot_n;

            td[10] = alpm/(alpf*gamm*dt);             // d_{n+1}
            td[11] = -td[10];                         // d_n
            td[12] = (alpf-1.0)/alpf;                 // v_n
            td[13] = 0.0;                             // a_n
            td[14] = (gamm-alpm)/alpf/gamm;           // ddot_n

            td[15] = alpm/(alpf*gamm*gamm*dt*dt);     // d_{n+1}
            td[16] = -td[15];                         // d_n
            td[17] = -1.0/(alpf*gamm*dt);             // v_n
            td[18] = (gamm-1.0)/gamm;                 // a_n
            td[19] = (gamm-alpm)/(alpf*gamm*gamm*dt); // ddot_n

            td[20] = 1.0/(gamm*dt);                   // d_{n+1}
            td[21] = -td[20];                         // d_n
            td[22] = 0.0;                             // v_n
            td[23] = 0.0;                             // a_n
            td[24] = (gamm-1.0)/gamm;                 // ddot_n

            //velocity as the primary variable
            //d_{n+1}    = td[20]*v_{n+1} + td[21]*d_n + td[22]*v_n + td[23]*a_n + td[24]*ddot_n ;
            //a_{n+1}    = td[25]*v_{n+1} + td[26]*d_n + td[27]*v_n + td[28]*a_n + td[29]*ddot_n ;
            //ddot_{n+1} = td[30]*v_{n+1} + td[31]*d_n + td[32]*v_n + td[33]*a_n + td[34]*ddot_n;

            td[40] = alpf*gamm*dt/alpm;           // v_{n+1}
            td[41] = 1.0;                         // d_n
            td[42] = (1.0-alpf)*gamm*dt/alpm;     // v_n
            td[43] = 0.0;                         // a_n
            td[44] = (alpm-gamm)*dt/alpm;        // ddot_n

            td[45] = 1.0/gamm/dt;                 // v_{n+1}
            td[46] = 0.0;                         // d_n
            td[47] = -td[45];                     // v_n
            td[48] = 1.0-1.0/gamm;                // v_n
            td[49] = 0.0;                         // ddot_n

            td[50] = alpf/alpm;                   // v_{n+1}
            td[51] = 0.0;                         // d_n
            td[52] = (1.0-alpf)/alpm;             // v_n
            td[53] = 0.0;                         // a_n
            td[54] = (alpm-1.0)/alpm;             // ddot_n

      break;

      default:

            cerr << "SetTimeParametersSolid ... invalid value of tis!"  << endl;

      break;
  }

  return;
}



void create_vector(double start, double end, double incr, vector<double>& uuu)
{
  int steps = int (((end-start)/incr) + 1);
  uuu.resize(steps);

  uuu[0] = start;
  for(int ii=1;ii<steps;ii++)
    uuu[ii] = uuu[ii-1] + incr;

  return;
}



void map2DPointTo3DPoint(int side, myPoint& ptTemp, double val3)
{
    switch(side)
    {
        case 0:
        case 1:
                ptTemp[2] = ptTemp[1] ;
                ptTemp[1] = ptTemp[0] ;
                ptTemp[0] = val3 ;

        break;

        case 2:
        case 3:

                ptTemp[0] = ptTemp[0] ;
                ptTemp[2] = ptTemp[1] ;
                ptTemp[1] = val3 ;

        break;

        case 4:
        case 5:

                ptTemp[0] = ptTemp[0] ;
                ptTemp[1] = ptTemp[1] ;
                ptTemp[2] = val3 ;

        break;

        default :

            cout << " Invalid 'side' in map2DPointTo3DPoint in myMappings.h " << endl;
        break;
    } //switch(side)
  return;
}





double factorial(unsigned int nn)
{
  if(nn == 0 || nn == 1)
    return 1.0;
  else
  {
    double result=1.0;
    for(unsigned int ii=1;ii<=nn;ii++)
      result *= ii;

    return result;
  }
}


double Bin(unsigned int m, unsigned int n)
{
  if((m == n) || (n == 0))
    return 1;
  else if(n > m)
    return 0;
  else
  {
    double num=1.0;
    int decr = m;
    do
    {
      num = num * decr;
      decr--;
    }while(decr > (m-n));

    return num/factorial(n);
  }
}



//RED, BLUE, GREEN, YELLOW, CYAN, MAGENTA, LIGHTBLUE, WHITE, BLACK

void  getColorValue(int col, double* color)
{
    switch(col)
    {
       case 0 : //RED
           color[0] = 1.0;   color[1] = 0.0;   color[2] = 0.0;
         break;

       case 1 : //BLUE
           color[0] = 0.0;   color[1] = 0.0;   color[2] = 1.0;
         break;

       case 2 : //GREEN
           color[0] = 0.0;   color[1] = 1.0;   color[2] = 0.0;
         break;

       case 3 : //YELLOW
           color[0] = 1.0;   color[1] = 1.0;   color[2] = 0.0;
         break;

       case 4 : //CYAN
           color[0] = 0.0;   color[1] = 1.0;   color[2] = 1.0;
         break;

       case 5 : //MAGENTA
           color[0] = 1.0;   color[1] = 0.0;   color[2] = 1.0;
         break;

       case 6 : //LIGHT BLUE
           color[0] = 1.0;   color[1] = 0.0;   color[2] = 0.0;
         break;

       case 7 : //WHILE
           color[0] = 1.0;   color[1] = 1.0;   color[2] = 1.0;
         break;

       case 8 : //BLACK
           color[0] = 0.0;   color[1] = 0.0;   color[2] = 0.0;
         break;

    }

  return;
}




double dotProductVecs(double* vec1, double* vec2, int N)
{
   double  val=0.0;
   
   for(int ii=0;ii<N;ii++)
     val += (vec1[ii] * vec2[ii]);

   return val;
}



