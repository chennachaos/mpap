#include "QuadratureUtil.h"

#include <iostream>
#include <math.h>

using namespace std;


void getGaussPoints1D(int ngp, vector<double>& gausspoints, vector<double>& gaussweights)
{
  gausspoints.resize(ngp);
  gaussweights.resize(ngp);

  switch(ngp)
  {
     case 1:  // 1 Point quadrature rule

            gausspoints[0] = 0.0;        gaussweights[0] = 2.0;

            break;

     case 2:  //2 Point quadrature rule

            gausspoints[0] = -0.577350269189626;    gaussweights[0] = 1.0;
            gausspoints[1] =  0.577350269189626;    gaussweights[1] = 1.0;

            break;

     case 3:  //3 Point quadrature rule

            gausspoints[0] = -0.774596669241483;  gaussweights[0] = 0.555555555555556;
            gausspoints[1] =  0.0;                gaussweights[1] = 0.888888888888889;
            gausspoints[2] =  0.774596669241483;  gaussweights[2] = 0.555555555555556;

            break;

     case 4:  //4 Point quadrature rule

            gausspoints[0] = -0.861136311594953;    gaussweights[0] = 0.347854845137454;
            gausspoints[1] = -0.339981043584856;    gaussweights[1] = 0.652145154862546;
            gausspoints[2] =  0.339981043584856;    gaussweights[2] = 0.652145154862546;
            gausspoints[3] =  0.861136311594953;    gaussweights[3] = 0.347854845137454;

            break;

     case 5:  //5 Point quadrature rule

            gausspoints[0] = -0.906179845938664;    gaussweights[0] = 0.236926885056189;
            gausspoints[1] = -0.538469310105683;    gaussweights[1] = 0.478628670499366;
            gausspoints[2] =  0.0;                  gaussweights[2] = 0.568888888888889;
            gausspoints[3] =  0.538469310105683;    gaussweights[3] = 0.478628670499366;
            gausspoints[4] =  0.906179845938664;    gaussweights[4] = 0.236926885056189;

            break;

     case 6:  //6 Point quadrature rule

            gausspoints[0] = -0.932469514203152;     gaussweights[0] = 0.171324492379170;
            gausspoints[1] = -0.661209386466265;     gaussweights[1] = 0.360761573048139;
            gausspoints[2] = -0.238619186083197;     gaussweights[2] = 0.467913934572691;
            gausspoints[3] =  0.238619186083197;     gaussweights[3] = 0.467913934572691;
            gausspoints[4] =  0.661209386466265;     gaussweights[4] = 0.360761573048139;
            gausspoints[5] =  0.932469514203152;     gaussweights[5] = 0.171324492379170;

            break;

     case 7:  //7 Point quadrature rule

            gausspoints[0] = -0.9491079123427585245261897;	gaussweights[0] = 0.1294849661688696932706114 ;
            gausspoints[1] = -0.7415311855993944398638648; 	gaussweights[1] = 0.2797053914892766679014678 ;
            gausspoints[2] = -0.4058451513773971669066064; 	gaussweights[2] = 0.3818300505051189449503698 ;
            gausspoints[3] =  0.0 ;			        gaussweights[3] = 0.4179591836734693877551020 ;
            gausspoints[4] =  0.4058451513773971669066064; 	gaussweights[4] = 0.3818300505051189449503698 ;
            gausspoints[5] =  0.7415311855993944398638648; 	gaussweights[5] = 0.2797053914892766679014678 ;
            gausspoints[6] =  0.9491079123427585245261897; 	gaussweights[6] = 0.1294849661688696932706114 ;

            break;


     case 8:  //8 Point quadrature rule

            gausspoints[0] = -0.96028986; 	gaussweights[0] = 0.10122854 ;
            gausspoints[1] = -0.79666648; 	gaussweights[1] = 0.22238103 ;
            gausspoints[2] = -0.52553241; 	gaussweights[2] = 0.31370665 ;
            gausspoints[3] = -0.18343464; 	gaussweights[3] = 0.36268378 ;
            gausspoints[4] =  0.18343464; 	gaussweights[4] = 0.36268378 ;
            gausspoints[5] =  0.52553241; 	gaussweights[5] = 0.31370665 ;
            gausspoints[6] =  0.79666648; 	gaussweights[6] = 0.22238103 ;
            gausspoints[7] =  0.96028986; 	gaussweights[7] = 0.10122854 ;

            break;

     case 9:  //9 Point quadrature rule

            gaussweights[0] = 0.0812743883615744; 	gausspoints[0] = -0.9681602395076261;
            gaussweights[1] = 0.1806481606948574; 	gausspoints[1] = -0.8360311073266358;
            gaussweights[2] = 0.2606106964029354; 	gausspoints[2] = -0.6133714327005904;
            gaussweights[3] = 0.3123470770400029; 	gausspoints[3] = -0.3242534234038089;
            gaussweights[4] = 0.3302393550012598; 	gausspoints[4] = 0.0000000000000000;
            gaussweights[5] = 0.3123470770400029; 	gausspoints[5] = 0.3242534234038089;
            gaussweights[6] = 0.2606106964029354; 	gausspoints[6] = 0.6133714327005904;
            gaussweights[7] = 0.1806481606948574; 	gausspoints[7] = 0.8360311073266358;
            gaussweights[8] = 0.0812743883615744; 	gausspoints[8] = 0.9681602395076261;


     case 10:  //10 Point quadrature rule

            gaussweights[0] = 0.0666713443086881; 	gausspoints[0] = -0.9739065285171717;
            gaussweights[1] = 0.1494513491505806; 	gausspoints[1] = -0.8650633666889845;
            gaussweights[2] = 0.2190863625159820; 	gausspoints[2] = -0.6794095682990244;
            gaussweights[3] = 0.2692667193099963; 	gausspoints[3] = -0.4333953941292472;
            gaussweights[4] = 0.2955242247147529; 	gausspoints[4] = -0.1488743389816312;
            gaussweights[5] = 0.2955242247147529; 	gausspoints[5] = 0.1488743389816312;
            gaussweights[6] = 0.2692667193099963; 	gausspoints[6] = 0.4333953941292472;
            gaussweights[7] = 0.2190863625159820; 	gausspoints[7] = 0.6794095682990244;
            gaussweights[8] = 0.1494513491505806; 	gausspoints[8] = 0.8650633666889845;
            gaussweights[9] = 0.0666713443086881; 	gausspoints[9] = 0.9739065285171717;


            break;

     case 11:  //11 Point quadrature rule

            gaussweights[0] = 0.0556685671161737; 	gausspoints[0] = -0.9782286581460570;
            gaussweights[1] = 0.1255803694649046; 	gausspoints[1] = -0.8870625997680953;
            gaussweights[2] = 0.1862902109277343; 	gausspoints[2] = -0.7301520055740494;
            gaussweights[3] = 0.2331937645919905; 	gausspoints[3] = -0.5190961292068118;
            gaussweights[4] = 0.2628045445102467; 	gausspoints[4] = -0.2695431559523450;
            gaussweights[5] = 0.2729250867779006; 	gausspoints[5] = 0.0000000000000000;
            gaussweights[6] = 0.2628045445102467; 	gausspoints[6] = 0.2695431559523450;
            gaussweights[7] = 0.2331937645919905; 	gausspoints[7] = 0.5190961292068118;
            gaussweights[8] = 0.1862902109277343; 	gausspoints[8] = 0.7301520055740494;
            gaussweights[9] = 0.1255803694649046; 	gausspoints[9] = 0.8870625997680953;
            gaussweights[10] = 0.0556685671161737; 	gausspoints[10] = 0.9782286581460570;


            break;

     case 12:  //12 Point quadrature rule

            gaussweights[0] = 0.0471753363865118; 	gausspoints[0] = -0.9815606342467192;
            gaussweights[1] = 0.1069393259953184; 	gausspoints[1] = -0.9041172563704749;
            gaussweights[2] = 0.1600783285433462; 	gausspoints[2] = -0.7699026741943047;
            gaussweights[3] = 0.2031674267230659; 	gausspoints[3] = -0.5873179542866175;
            gaussweights[4] = 0.2334925365383548; 	gausspoints[4] = -0.3678314989981802;
            gaussweights[5] = 0.2491470458134028; 	gausspoints[5] = -0.1252334085114689;
            gaussweights[6] = 0.2491470458134028; 	gausspoints[6] = 0.1252334085114689;
            gaussweights[7] = 0.2334925365383548; 	gausspoints[7] = 0.3678314989981802;
            gaussweights[8] = 0.2031674267230659; 	gausspoints[8] = 0.5873179542866175;
            gaussweights[9] = 0.1600783285433462; 	gausspoints[9] = 0.7699026741943047;
            gaussweights[10] = 0.1069393259953184; 	gausspoints[10] = 0.9041172563704749;
            gaussweights[11] = 0.0471753363865118; 	gausspoints[11] = 0.9815606342467192;

            break;

     case 13:  //13 Point quadrature rule

            gaussweights[0] = 0.0404840047653159; 	gausspoints[7] = -0.9841830547185881;
            gaussweights[1] = 0.0921214998377285; 	gausspoints[7] = -0.9175983992229779;
            gaussweights[2] = 0.1388735102197872; 	gausspoints[7] = -0.8015780907333099;
            gaussweights[3] = 0.1781459807619457; 	gausspoints[7] = -0.6423493394403402;
            gaussweights[4] = 0.2078160475368885; 	gausspoints[7] = -0.4484927510364469;
            gaussweights[5] = 0.2262831802628972; 	gausspoints[7] = -0.2304583159551348;
            gaussweights[6] = 0.2325515532308739; 	gausspoints[7] = 0.0000000000000000;
            gaussweights[7] = 0.2262831802628972; 	gausspoints[7] = 0.2304583159551348;
            gaussweights[8] = 0.2078160475368885; 	gausspoints[7] = 0.4484927510364469;
            gaussweights[9] = 0.1781459807619457; 	gausspoints[7] = 0.6423493394403402;
            gaussweights[10] = 0.1388735102197872; 	gausspoints[10] = 0.8015780907333099;
            gaussweights[11] = 0.0921214998377285; 	gausspoints[11] = 0.9175983992229779;
            gaussweights[12] = 0.0404840047653159; 	gausspoints[12] = 0.9841830547185881;

            break;

     case 14:  //14 Point quadrature rule

            gaussweights[0] = 0.0351194603317519; 	gausspoints[0] = -0.9862838086968123;
            gaussweights[1] = 0.0801580871597602; 	gausspoints[1] = -0.9284348836635735;
            gaussweights[2] = 0.1215185706879032; 	gausspoints[2] = -0.8272013150697650;
            gaussweights[3] = 0.1572031671581935; 	gausspoints[3] = -0.6872929048116855;
            gaussweights[4] = 0.1855383974779378; 	gausspoints[4] = -0.5152486363581541;
            gaussweights[5] = 0.2051984637212956; 	gausspoints[5] = -0.3191123689278897;
            gaussweights[6] = 0.2152638534631578; 	gausspoints[6] = -0.1080549487073437;
            gaussweights[7] = 0.2152638534631578; 	gausspoints[7] = 0.1080549487073437;
            gaussweights[8] = 0.2051984637212956; 	gausspoints[8] = 0.3191123689278897;
            gaussweights[9] = 0.1855383974779378; 	gausspoints[9] = 0.5152486363581541;
            gaussweights[10] = 0.1572031671581935; 	gausspoints[10] = 0.6872929048116855;
            gaussweights[11] = 0.1215185706879032; 	gausspoints[11] = 0.8272013150697650;
            gaussweights[12] = 0.0801580871597602; 	gausspoints[12] = 0.9284348836635735;
            gaussweights[13] = 0.0351194603317519; 	gausspoints[13] = 0.9862838086968123;

            break;

     case 15:  //15 Point quadrature rule

            gaussweights[0] = 0.0307532419961173;     gausspoints[0] = -0.9879925180204854;
            gaussweights[1] = 0.0703660474881081;     gausspoints[1] = -0.9372733924007060;
            gaussweights[2] = 0.1071592204671719;     gausspoints[2] = -0.8482065834104272;
            gaussweights[3] = 0.1395706779261543;     gausspoints[3] = -0.7244177313601701;
            gaussweights[4] = 0.1662692058169939;     gausspoints[4] = -0.5709721726085388;
            gaussweights[5] = 0.1861610000155622;     gausspoints[5] = -0.3941513470775634;
            gaussweights[6] = 0.1984314853271116;     gausspoints[6] = -0.2011940939974345;
            gaussweights[7] = 0.2025782419255613;     gausspoints[7] =  0.0000000000000000;
            gaussweights[8] = 0.1984314853271116;     gausspoints[8] =  0.2011940939974345;
            gaussweights[9] = 0.1861610000155622;     gausspoints[9] =  0.3941513470775634;
            gaussweights[10] = 0.1662692058169939;    gausspoints[10] = 0.5709721726085388;
            gaussweights[11] = 0.1395706779261543;    gausspoints[11] = 0.7244177313601701;
            gaussweights[12] = 0.1071592204671719;    gausspoints[12] = 0.8482065834104272;
            gaussweights[13] = 0.0703660474881081;    gausspoints[13] = 0.9372733924007060;
            gaussweights[14] = 0.0307532419961173;    gausspoints[14] = 0.9879925180204854;

            break;

     default:
            cerr << " invalid value of 'ngp' ! " << endl;

          break;

  }
  return;
}



void getGaussPointsHex(int ngp, vector<double>& gp1, vector<double>& gp2, vector<double>& gp3, vector<double>& gws)
{
  ngp = pow(ngp,1.0/3.0);

  gp1.resize(ngp);
  gp2.resize(ngp);
  gp3.resize(ngp);
  gws.resize(ngp);

  vector<double> gpoints1, gweights1;

  getGaussPoints1D(ngp, gpoints1, gweights1);

  int ii, jj, kk, ind;
  
  ind=0;
  for(kk=0; kk<ngp; kk++)
  {
    for(jj=0; jj<ngp; jj++)
    {
      for(ii=0; ii<ngp; ii++)
      {
        gp1[ind] = gpoints1[ii];
        gp2[ind] = gpoints1[jj];
        gp3[ind] = gpoints1[kk];
        gws[ind] = gweights1[kk]*gweights1[jj]*gweights1[ii];
        ind++;
      }
    }
  }

  return;
}




void getGaussPointsQuad(int ngp, vector<double>& gp1, vector<double>& gp2, vector<double>& gws)
{
  //assert(ngp>0);
  
  gp1.resize(ngp);
  gp2.resize(ngp);
  gws.resize(ngp);
  
  int  nn=sqrt(ngp), ii, jj, kk;
  
  vector<double> gpoints1, gweights1;
  vector<double> gpoints2, gweights2;
  
  getGaussPoints1D(nn, gpoints1, gweights1);
  getGaussPoints1D(nn, gpoints2, gweights2);
  
  kk=0;
  for(jj=0; jj<nn; jj++)
  {
    for(ii=0; ii<nn; ii++)
    {
      gp1[kk] = gpoints1[ii];  gp2[kk] = gpoints2[jj];  gws[kk] = gweights2[jj]*gweights1[ii];
      kk++;
    }
  }

  return;
}





void getGaussPointsTet(int ngp, vector<double>& gps1, vector<double>& gps2, vector<double>& gps3, vector<double>& gws)
{
  char fct[] = "getGaussPointsTet";

  double  fact=1.0/6.0;

  gps1.resize(ngp);
  gps2.resize(ngp);
  gps3.resize(ngp);
  gws.resize(ngp);

  switch(ngp)
  {
     case 1:  // 1 Point quadrature rule

            gps1[0] = 0.25;
            gps2[0] = 0.25;
            gps3[0] = 0.25;
            gws[0] = fact;

            break;

     case 4:  // N=4

        gps1[0] = 0.5854101966249685;
        gps1[1] = 0.1381966011250105;
        gps1[2] = 0.1381966011250105;
        gps1[3] = 0.1381966011250105;

        gps2[0] = 0.1381966011250105;
        gps2[1] = 0.1381966011250105;
        gps2[2] = 0.1381966011250105;
        gps2[3] = 0.5854101966249685;

        gps3[0] = 0.1381966011250105;
        gps3[1] = 0.1381966011250105;
        gps3[2] = 0.5854101966249685;
        gps3[3] = 0.1381966011250105;

        gws[0]  = 0.2500000000000000*fact;
        gws[1]  = 0.2500000000000000*fact;
        gws[2]  = 0.2500000000000000*fact;
        gws[3]  = 0.2500000000000000*fact;

        break;

     case 5:  // N=5

        gps1[0] = 0.2500000000000000;
        gps1[1] = 0.5000000000000000;
        gps1[2] = 0.1666666666666667;
        gps1[3] = 0.1666666666666667;
        gps1[4] = 0.1666666666666667;
 
        gps2[0] = 0.2500000000000000;
        gps2[1] = 0.1666666666666667;
        gps2[2] = 0.1666666666666667;
        gps2[3] = 0.1666666666666667;
        gps2[4] = 0.5000000000000000;

        gps3[0] = 0.2500000000000000;
        gps3[1] = 0.1666666666666667;
        gps3[2] = 0.1666666666666667;
        gps3[3] = 0.5000000000000000;
        gps3[4] = 0.1666666666666667;
 
        gws[0]  = -0.8000000000000000*fact;
        gws[1]  =  0.4500000000000000*fact;
        gws[2]  =  0.4500000000000000*fact;
        gws[3]  =  0.4500000000000000*fact;
        gws[4]  =  0.4500000000000000*fact;

        break;

     case 11:  // N=11

        gps1[0] = 0.2500000000000000;
        gps1[1] = 0.7857142857142857;
        gps1[2] = 0.0714285714285714;
        gps1[3] = 0.0714285714285714;
        gps1[4] = 0.0714285714285714;
        gps1[5] = 0.1005964238332008;
        gps1[6] = 0.3994035761667992;
        gps1[7] = 0.3994035761667992;
        gps1[8] = 0.3994035761667992;
        gps1[9] = 0.1005964238332008;
        gps1[10] = 0.1005964238332008;

        gps2[0] = 0.2500000000000000;
        gps2[1] = 0.0714285714285714;
        gps2[2] = 0.0714285714285714;
        gps2[3] = 0.0714285714285714;
        gps2[4] = 0.7857142857142857;
        gps2[5] = 0.3994035761667992;
        gps2[6] = 0.1005964238332008;
        gps2[7] = 0.3994035761667992;
        gps2[8] = 0.1005964238332008;
        gps2[9] = 0.3994035761667992;
        gps2[10] = 0.1005964238332008;

        gps3[0] = 0.2500000000000000;
        gps3[1] = 0.0714285714285714;
        gps3[2] = 0.0714285714285714;
        gps3[3] = 0.7857142857142857;
        gps3[4] = 0.0714285714285714;
        gps3[5] = 0.3994035761667992;
        gps3[6] = 0.3994035761667992;
        gps3[7] = 0.1005964238332008;
        gps3[8] = 0.1005964238332008;
        gps3[9] = 0.1005964238332008;
        gps3[10] = 0.3994035761667992;
 
        gws[0]  = -0.0789333333333333*fact;
        gws[1]  =  0.0457333333333333*fact;
        gws[2]  =  0.0457333333333333*fact;
        gws[3]  =  0.0457333333333333*fact;
        gws[4]  =  0.0457333333333333*fact;
        gws[5]  =  0.1493333333333333*fact;
        gws[6]  =  0.1493333333333333*fact;
        gws[7]  =  0.1493333333333333*fact;
        gws[8]  =  0.1493333333333333*fact;
        gws[9]  =  0.1493333333333333*fact;
        gws[10]  =  0.1493333333333333*fact;

        break;

     case 15:  //  N=15

        gps1[0] = 0.2500000000000000;
        gps1[1] = 0.0000000000000000;
        gps1[2] = 0.3333333333333333;
        gps1[3] = 0.3333333333333333;
        gps1[4] = 0.3333333333333333;
        gps1[5] = 0.7272727272727273;
        gps1[6] = 0.0909090909090909;
        gps1[7] = 0.0909090909090909;
        gps1[8] = 0.0909090909090909;
        gps1[9] = 0.4334498464263357;
        gps1[10] = 0.0665501535736643;
        gps1[11] = 0.0665501535736643;
        gps1[12] = 0.0665501535736643;
        gps1[13] = 0.4334498464263357;
        gps1[14] = 0.4334498464263357;
        
        gps2[0] = 0.2500000000000000;
        gps2[1] = 0.3333333333333333;
        gps2[2] = 0.3333333333333333;
        gps2[3] = 0.3333333333333333;
        gps2[4] = 0.0000000000000000;
        gps2[5] = 0.0909090909090909;
        gps2[6] = 0.0909090909090909;
        gps2[7] = 0.0909090909090909;
        gps2[8] = 0.7272727272727273;
        gps2[9] = 0.0665501535736643;
        gps2[10] = 0.4334498464263357;
        gps2[11] = 0.0665501535736643;
        gps2[12] = 0.4334498464263357;
        gps2[13] = 0.0665501535736643;
        gps2[14] = 0.4334498464263357;

        gps3[0] = 0.2500000000000000;
        gps3[1] = 0.3333333333333333;
        gps3[2] = 0.3333333333333333;
        gps3[3] = 0.0000000000000000;
        gps3[4] = 0.3333333333333333;
        gps3[5] = 0.0909090909090909;
        gps3[6] = 0.0909090909090909;
        gps3[7] = 0.7272727272727273;
        gps3[8] = 0.0909090909090909;
        gps3[9] = 0.0665501535736643;
        gps3[10] = 0.0665501535736643;
        gps3[11] = 0.4334498464263357;
        gps3[12] = 0.4334498464263357;
        gps3[13] = 0.4334498464263357;
        gps3[14] = 0.0665501535736643;

        gws[0]  = 0.1817020685825351*fact;
        gws[1]  = 0.0361607142857143*fact;
        gws[2]  = 0.0361607142857143*fact;
        gws[3]  = 0.0361607142857143*fact;
        gws[4]  = 0.0361607142857143*fact;
        gws[5]  = 0.0698714945161738*fact;
        gws[6]  = 0.0698714945161738*fact;
        gws[7]  = 0.0698714945161738*fact;
        gws[8]  = 0.0698714945161738*fact;
        gws[9]  = 0.0656948493683187*fact;
        gws[10]  = 0.0656948493683187*fact;
        gws[11]  = 0.0656948493683187*fact;
        gws[12]  = 0.0656948493683187*fact;
        gws[13]  = 0.0656948493683187*fact;
        gws[14]  = 0.0656948493683187*fact;

        break;

     default:
            cerr << " getGaussPointsTet() .... invalid value of 'ngp' ! " << endl;
            break;

  }

  return;
}





void getGaussPointsTriangle(int ngp, vector<double>& gps1, vector<double>& gps2, vector<double>& gws)
{
  // weights are normalized to calculate the exact area of the triangle
  // i.e. each weight is divided by 2.0

  char fct[] = "getGaussPointsTriangle";

  double r1d3 = 1.0/3.0, a1, fact=0.5;

  gps1.resize(ngp);
  gps2.resize(ngp);
  gws.resize(ngp);

  switch(ngp)
  {
     case 1:  // 1 Point quadrature rule

            gps1[0] = r1d3;        gps2[0] = r1d3;        gws[0] = fact*1.0;

            break;

     case 3:  //3 Point quadrature rule

            a1 = 1.0/3.0;

            //gps1[0] = 0.5;        gps2[0] = 0.0;        gws[0] = fact*a1;
            //gps1[1] = 0.5;        gps2[1] = 0.5;        gws[1] = fact*a1;
            //gps1[2] = 0.0;        gps2[2] = 0.5;        gws[2] = fact*a1;

            gps1[0] = 1.0/6.0;        gps2[0] = 1.0/6.0;        gws[0] = fact*a1;
            gps1[1] = 1.0/6.0;        gps2[1] = 4.0/6.0;        gws[1] = fact*a1;
            gps1[2] = 4.0/6.0;        gps2[2] = 1.0/6.0;        gws[2] = fact*a1;

            break;

     case 4:  //4 Point quadrature rule

            a1 = 25.0/48.0;

            gps1[0] = r1d3;       gps2[0] = r1d3;       gws[0] = fact*(-27.0/48.0);
            gps1[1] = 0.6;        gps2[1] = 0.2;        gws[1] = fact*a1;
            gps1[2] = 0.2;        gps2[2] = 0.6;        gws[2] = fact*a1;
            gps1[3] = 0.2;        gps2[3] = 0.2;        gws[3] = fact*a1;

            break;

     case 6:  //6 Point quadrature rule

          gps1[0] = 0.10810301816807022736;    gps2[0] = 0.44594849091596488632;    gws[0] = fact*0.22338158967801146570;
          gps1[1] = 0.44594849091596488632;    gps2[1] = 0.10810301816807022736;    gws[1] = fact*0.22338158967801146570;
          gps1[2] = 0.44594849091596488632;    gps2[2] = 0.44594849091596488632;    gws[2] = fact*0.22338158967801146570;
          gps1[3] = 0.81684757298045851308;    gps2[3] = 0.09157621350977074346;    gws[3] = fact*0.10995174365532186764;
          gps1[4] = 0.09157621350977074346;    gps2[4] = 0.81684757298045851308;    gws[4] = fact*0.10995174365532186764;
          gps1[5] = 0.09157621350977074346;    gps2[5] = 0.09157621350977074346;    gws[5] = fact*0.10995174365532186764;

            break;

     case 7:  //7 Point quadrature rule

          gps1[0] = r1d3;                      gps2[0] = r1d3;                      gws[0] = fact*0.225;
          gps1[1] = 0.79742698535308732240;    gps2[1] = 0.10128650732345633880;    gws[1] = fact*0.12593918054482715260;
          gps1[2] = 0.10128650732345633880;    gps2[2] = 0.79742698535308732240;    gws[2] = fact*0.12593918054482715260;
          gps1[3] = 0.10128650732345633880;    gps2[3] = 0.10128650732345633880;    gws[3] = fact*0.12593918054482715260;
          gps1[4] = 0.05971587178976982045;    gps2[4] = 0.47014206410511508977;    gws[4] = fact*0.13239415278850618074;
          gps1[5] = 0.47014206410511508977;    gps2[5] = 0.05971587178976982045;    gws[5] = fact*0.13239415278850618074;
          gps1[6] = 0.47014206410511508977;    gps2[6] = 0.47014206410511508977;    gws[6] = fact*0.13239415278850618074;

            break;

     case 12:  //12 Point quadrature rule

          gps1[0]  = 0.87382197101699554332;   gps2[0]  = 0.06308901449150222834;   gws[0]  = fact*0.050844906370206816921;
          gps1[1]  = 0.06308901449150222834;   gps2[1]  = 0.87382197101699554332;   gws[1]  = fact*0.050844906370206816921;
          gps1[2]  = 0.06308901449150222834;   gps2[2]  = 0.06308901449150222834;   gws[2]  = fact*0.050844906370206816921;
          gps1[3]  = 0.50142650965817915742;   gps2[3]  = 0.24928674517091042129;   gws[3]  = fact*0.116786275726379366030;
          gps1[4]  = 0.24928674517091042129;   gps2[4]  = 0.50142650965817915742;   gws[4]  = fact*0.116786275726379366030;
          gps1[5]  = 0.24928674517091042129;   gps2[5]  = 0.24928674517091042129;   gws[5]  = fact*0.116786275726379366030;
          gps1[6]  = 0.05314504984481694735;   gps2[6]  = 0.31035245103378440542;   gws[6]  = fact*0.082851075618373575194;
          gps1[7]  = 0.31035245103378440542;   gps2[7]  = 0.05314504984481694735;   gws[7]  = fact*0.082851075618373575194;
          gps1[8]  = 0.05314504984481694735;   gps2[8]  = 0.63650249912139864723;   gws[8]  = fact*0.082851075618373575194;
          gps1[9]  = 0.31035245103378440542;   gps2[9]  = 0.63650249912139864723;   gws[9]  = fact*0.082851075618373575194;
          gps1[10] = 0.63650249912139864723;   gps2[10] = 0.05314504984481694735;   gws[10] = fact*0.082851075618373575194;
          gps1[11] = 0.63650249912139864723;   gps2[11] = 0.31035245103378440542;   gws[11] = fact*0.082851075618373575194;

            break;

     case 13:  // 13 point quadrature rule

          gps1[0]  = 0.33333333333333;   gps2[0]  =  0.33333333333333;   gws[0]  = fact*-0.14957004446768;
          gps1[1]  = 0.26034596607904;   gps2[1]  =  0.26034596607904;   gws[1]  = fact* 0.17561525743321;
          gps1[2]  = 0.26034596607904;   gps2[2]  =  0.47930806784192;   gws[2]  = fact* 0.17561525743321;
          gps1[3]  = 0.47930806784192;   gps2[3]  =  0.26034596607904;   gws[3]  = fact* 0.17561525743321;
          gps1[4]  = 0.06513010290222;   gps2[4]  =  0.06513010290222;   gws[4]  = fact* 0.05334723560884;
          gps1[5]  = 0.06513010290222;   gps2[5]  =  0.86973979419557;   gws[5]  = fact* 0.05334723560884;
          gps1[6]  = 0.86973979419557;   gps2[6]  =  0.06513010290222;   gws[6]  = fact* 0.05334723560884;
          gps1[7]  = 0.31286549600487;   gps2[7]  =  0.63844418856981;   gws[7]  = fact* 0.07711376089026;
          gps1[8]  = 0.63844418856981;   gps2[8]  =  0.04869031542532;   gws[8]  = fact* 0.07711376089026;
          gps1[9]  = 0.04869031542532;   gps2[9]  =  0.31286549600487;   gws[9]  = fact* 0.07711376089026;
          gps1[10] = 0.63844418856981;   gps2[10] =  0.31286549600487;   gws[10] = fact* 0.07711376089026;
          gps1[11] = 0.31286549600487;   gps2[11] =  0.04869031542532;   gws[11] = fact* 0.07711376089026;
          gps1[12] = 0.04869031542532;   gps2[12] =  0.63844418856981;   gws[12] = fact* 0.07711376089026;

            break;

     case 16:  // 16 point quadrature rule

          gps1[0]  = 0.33333333333333;   gps2[0]  =  0.33333333333333;   gws[0]  = fact* 0.14431560767779;
          gps1[1]  = 0.45929258829272;   gps2[1]  =  0.45929258829272;   gws[1]  = fact* 0.09509163426728;
          gps1[2]  = 0.45929258829272;   gps2[2]  =  0.08141482341455;   gws[2]  = fact* 0.09509163426728;
          gps1[3]  = 0.08141482341455;   gps2[3]  =  0.45929258829272;   gws[3]  = fact* 0.09509163426728;
          gps1[4]  = 0.17056930775176;   gps2[4]  =  0.17056930775176;   gws[4]  = fact* 0.10321737053472;
          gps1[5]  = 0.17056930775176;   gps2[5]  =  0.65886138449648;   gws[5]  = fact* 0.10321737053472;
          gps1[6]  = 0.65886138449648;   gps2[6]  =  0.17056930775176;   gws[6]  = fact* 0.10321737053472;
          gps1[7]  = 0.05054722831703;   gps2[7]  =  0.05054722831703;   gws[7]  = fact* 0.03245849762320;
          gps1[8]  = 0.05054722831703;   gps2[8]  =  0.89890554336594;   gws[8]  = fact* 0.03245849762320;
          gps1[9]  = 0.89890554336594;   gps2[9]  =  0.05054722831703;   gws[9]  = fact* 0.03245849762320;
          gps1[10] = 0.26311282963464;   gps2[10] =  0.72849239295540;   gws[10] = fact* 0.02723031417443;
          gps1[11] = 0.72849239295540;   gps2[11] =  0.00839477740996;   gws[11] = fact* 0.02723031417443;
          gps1[12] = 0.00839477740996;   gps2[12] =  0.26311282963464;   gws[12] = fact* 0.02723031417443;
          gps1[13] = 0.72849239295540;   gps2[13] =  0.26311282963464;   gws[13] = fact* 0.02723031417443;
          gps1[14] = 0.26311282963464;   gps2[14] =  0.00839477740996;   gws[14] = fact* 0.02723031417443;
          gps1[15] = 0.00839477740996;   gps2[15] =  0.72849239295540;   gws[15] = fact* 0.02723031417443;

            break;

     default:
            cerr << " invalid value of 'ngp' ! " << endl;
            break;

  }
  return;
}




void getLobattoPoints(int ngp, vector<double>& gausspoints, vector<double>& gaussweights)
{
  char fct[] = "getGaussPoints";

  switch(ngp)
  {
     case 2:  //2 Point quadrature rule

            gausspoints.resize(2);
            gaussweights.resize(2);
            gausspoints[0] = -1.0;    gaussweights[0] = 1.0;
            gausspoints[1] =  1.0;    gaussweights[1] = 1.0;

            break;

     case 3:  //3 Point quadrature rule

            gausspoints.resize(3);
            gaussweights.resize(3);
            gausspoints[0] = -1.0;  gaussweights[0] = 0.333333333333;
            gausspoints[1] =  0.0;  gaussweights[1] = 1.333333333333;
            gausspoints[2] =  1.0;  gaussweights[2] = 0.333333333333;

            break;

     case 4:  //4 Point quadrature rule

            gausspoints.resize(4);
            gaussweights.resize(4);
            gausspoints[0] = -1.0;               gaussweights[0] = 1.0/6.0;
            gausspoints[1] = -0.447213595500;    gaussweights[1] = 5.0/6.0;
            gausspoints[2] =  0.447213595500;    gaussweights[2] = 5.0/6.0;
            gausspoints[3] =  1.0;               gaussweights[3] = 1.0/6.0;

            break;

     case 5:  //5 Point quadrature rule

            gausspoints.resize(5);
            gaussweights.resize(5);
            gausspoints[0] = -1.0;               gaussweights[0] = 0.1;
            gausspoints[1] = -0.654653670708;    gaussweights[1] = 0.54444444444;
            gausspoints[2] =  0.0;               gaussweights[2] = 0.71111111111;
            gausspoints[3] =  0.654653670708;    gaussweights[3] = 0.54444444444;
            gausspoints[4] =  1.0;               gaussweights[4] = 0.1;

            break;

     default:
            cerr << " getGaussPoints1D()... invalid value of 'ngp' ! " << endl;
            break;

  }
  return;
}

