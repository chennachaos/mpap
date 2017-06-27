
#include "HBSplineCutFEM.h"
#include "ComputerTime.h"
#include "MpapTime.h"
#include "ImmersedIntegrationElement.h"
#include "QuadratureUtil.h"


extern MpapTime           mpapTime;

using namespace std;


void  HBSplineCutFEM::computeTotalBodyForce(int index)
{
  cout << " HBSplineCutFEM::computeTotalBodyForce() ... not implemented yet ... " << endl;

  return;
}


//
void  HBSplineCutFEM::solveSolidProblem()
{
  IB_MOVED = false;
  for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
  {
    //computeTotalForce(bb);
    //ImmersedBodyObjects[bb]->updateForce(&(totalForce(0)));

    cout << " kkkkkkkkkkk " << endl;
    ImmersedBodyObjects[bb]->SolveTimeStep();

    //IB_MOVED = (IB_MOVED || ImmersedBodyObjects[bb]->updatePointPositions() );
  }
  
  updateImmersedPointPositions();

  return;
}
//



/*
void  HBSplineCutFEM::solveSolidProblem()
{
  IB_MOVED = false;
  
  cout << " Solving solid " << endl;


  int  nPart=ImmersedBodyObjects.size(), nCont=1, bb, b1, b2, doftemp, ii, jj, kk, iter, r1;
  double  xc1, xc2, yc1, yc2, Rad=0.1;
  double  l0=3.0*gridLEN[0]/nelem[0], gap, gap2, gap3, lam, c, s, theta;
  MatrixXd  Klocal(3,3), Kcont(2,6), Kglobal, matM(3,3);
  VectorXd  Flocal(3), Fglobal, tdTemp, xc(nPart), yc(nPart);

  cout << " aaaaaaaaaaaaaaaa " << endl;

  double  cn=1.0e5, tol=1.0e-6, dx, dy, af, rhof=1.0, rhos=1.2*rhof, As, normTemp, Fx, Fy, lij;
  double  dlijdix, dlijdiy, dlijdjx, dlijdjy;

  //l0=0.2001;
  //xc[0] = 1.0; yc[0] = 6.8;

  xc[0] = 1.0; yc[0] = 6.8;
  xc[1] = 1.0; yc[1] = 7.2;

  //xc[0] = 0.8; yc[0] = 6.8;
  //xc[1] = 1.2; yc[1] = 6.8;
  //xc[2] = 0.8; yc[2] = 7.2;
  //xc[3] = 1.2; yc[3] = 7.2;

  doftemp = 3*nPart;

  cout << " doftemp  = " << doftemp << endl;
  cout << " l0 = " << l0 << endl;

  Kglobal.resize(doftemp, doftemp);
  Fglobal.resize(doftemp);


  tdTemp = ImmersedBodyObjects[0]->SolnData.td;

  af = tdTemp[2];
  
  cout << " af = " << af << endl;
  
  As = PI*Rad*Rad;

  matM.setZero();
  matM(0,0) = As*rhos;
  matM(1,1) = As*rhos;
  matM(2,2) = 0.5*PI*Rad*Rad*Rad*Rad*rhos;

  double  F0 = (rhos-rhof)*As*981.0;

  Klocal = tdTemp[5]*matM;
  
  vector<vector<int> >  forAssy1, contIBs;
  
  forAssy1.resize(nPart);
  kk=0;
  for(bb=0; bb<nPart; bb++)
  {
    forAssy1[bb].resize(3);
    for(ii=0; ii<3; ii++)
      forAssy1[bb][ii] = kk++;

    printVector(forAssy1[bb]);
  }

  contIBs.resize(nCont);
  for(bb=0; bb<nCont; bb++)
    contIBs[bb].resize(2);
  
  if( nCont > 0 )
  {
    contIBs[0][0] = 0; contIBs[0][1] = 1;
    //contIBs[1][0] = 0; contIBs[1][1] = 2;
    //contIBs[2][0] = 0; contIBs[2][1] = 3;
    //contIBs[3][0] = 1; contIBs[3][1] = 2;
    //contIBs[4][0] = 1; contIBs[4][1] = 3;
    //contIBs[5][0] = 2; contIBs[5][1] = 3;
  }

  //slnTemp = slnTempPrev;

  for(iter=0; iter<10; iter++)
  {
    cout << " iter = " << iter << endl;
    
    Kglobal.setZero();
    Fglobal.setZero();

    //slnTempCur = af*slnTemp + (1.0-af)*slnTempPrev;

    for(b1=0; b1<nPart; b1++)    // loop over the particles
    {
      Flocal = ImmersedBodyObjects[b1]->SolnData.forceCur - matM*ImmersedBodyObjects[b1]->SolnData.var1DotDotCur;
      Flocal(1) -= F0;

      for(ii=0;ii<3;ii++)
      {
        Fglobal[forAssy1[b1][ii]] += Flocal[ii];
      
        for(jj=0;jj<3;jj++)
        {
          Kglobal(forAssy1[b1][ii],forAssy1[b1][jj])  +=  Klocal(ii,jj);
        }
      }


      xc1 = xc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[0];
      yc1 = yc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[1];

      // for each particle check whether it is in contact with the boundary walls
      for(b2=0;b2<4;b2++)
      {
        switch(b2)
        {
          case 0:
            xc2 = origin[0] - Rad;
            yc2 = yc1;
        
          break;

          case 1:

            xc2 = origin[0] + gridLEN[0] + Rad;
            yc2 = yc1;

          break;

          case 2:
        
            xc2 = xc1;
            yc2 = origin[1] - Rad;
        
          break;

          case 3:

            xc2 = xc1;
            yc2 = origin[1] + gridLEN[1] + Rad;

          break;

          default:
          break;
        }

        dx = xc1-xc2;
        dy = yc1-yc2;
        
        lij = sqrt(dx*dx+dy*dy);

        theta = atan2(dy, dx);

        c = cos(theta);
        s = sin(theta);

        gap = 2.0*Rad + l0 - lij;
        
        gap2 = gap*gap;
        gap3 = gap2*gap;

        Fx = cn*dx*gap2;
        Fy = cn*dy*gap2;

        if( (lamn + cn*gn) > tol)
        {
            cout << " Contact " << b1 << " to " << b2 << " is active " << endl;

            Kglobal(forAssy1[b1][0], forAssy1[b1][0]) -= (  cn*gap2 + cn*dx*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b1][0], forAssy1[b1][1]) -= (    0.0   + cn*dx*2.0*gap*(-dlijdix) );

            Kglobal(forAssy1[b1][1], forAssy1[b1][0]) -= (    0.0   + cn*dy*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b1][1], forAssy1[b1][1]) -= (  cn*gap2 + cn*dy*2.0*gap*(-dlijdix) );

            Fglobal(forAssy1[b1][0])    += Fx;
            Fglobal(forAssy1[b1][1])    += Fy;
        }
        else
        {
            //cout << " Contact " << b1 << " to " << b2 << " is inactive " << endl;
        }

      } // for(b2=0;b2<4;b2++)
    } //     for(b1=0; b1<nPart; b1++)    // loop over the particles

      cout << " contact with other particles " << endl;

      // check contact with other particles
      for(bb=0; bb<nCont; bb++)
      {
        b1 = contIBs[bb][0];
        b2 = contIBs[bb][1];

        xc1 = xc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[0];
        yc1 = yc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[1];

        xc2 = xc[b2] + ImmersedBodyObjects[b2]->SolnData.var1Cur[0];
        yc2 = yc[b2] + ImmersedBodyObjects[b2]->SolnData.var1Cur[1];

        dx = xc1-xc2;
        dy = yc1-yc2;
        
        lij = sqrt(dx*dx+dy*dy);

        gap = 2.0*Rad + l0 - lij;

        cout << " xc1 = " << xc1 << '\t' << " yc1 = " << yc1 << endl;
        cout << " xc2 = " << xc2 << '\t' << " yc2 = " << yc2 << endl;
        cout << " dx = " << dx << '\t' << " dy = " << dy << endl;
        cout << " gap = " << gap << endl;

        gap2 = gap*gap;

        Fx = cn*dx*gap2;
        Fy = cn*dy*gap2;

        //cout << " Fx = " << Fx << endl;
        //cout << " Fy = " << Fy << endl;

        theta = atan2(dy, dx);

        c = cos(theta);
        s = sin(theta);
        
        dlijdix =  dx/lij;
        dlijdiy =  dy/lij;
        dlijdjx = -dx/lij;
        dlijdjy = -dy/lij;

        if( gap > tol)
        {
            cout << " Contact " << b1 << " to " << b2 << " is active " << endl;

            Kglobal(forAssy1[b1][0], forAssy1[b1][0]) -= (  cn*gap2 + cn*dx*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b1][0], forAssy1[b1][1]) -= (    0.0   + cn*dx*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b1][0], forAssy1[b2][0]) -= ( -cn*gap2 + cn*dx*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b1][0], forAssy1[b2][1]) -= (    0.0   + cn*dx*2.0*gap*(-dlijdix) );

            Kglobal(forAssy1[b1][1], forAssy1[b1][0]) -= (    0.0   + cn*dy*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b1][1], forAssy1[b1][1]) -= (  cn*gap2 + cn*dy*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b1][1], forAssy1[b2][0]) -= (    0.0   + cn*dy*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b1][1], forAssy1[b2][1]) -= ( -cn*gap2 + cn*dy*2.0*gap*(-dlijdix) );

            Kglobal(forAssy1[b2][0], forAssy1[b1][0]) -= ( -cn*gap2 + cn*-dx*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b2][0], forAssy1[b1][1]) -= (    0.0   + cn*-dx*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b2][0], forAssy1[b2][0]) -= (  cn*gap2 + cn*-dx*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b2][0], forAssy1[b2][1]) -= (    0.0   + cn*-dx*2.0*gap*(-dlijdix) );

            Kglobal(forAssy1[b2][1], forAssy1[b1][0]) -= (    0.0   + cn*-dy*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b2][1], forAssy1[b1][1]) -= ( -cn*gap2 + cn*-dy*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b2][1], forAssy1[b2][0]) -= (    0.0   + cn*-dy*2.0*gap*(-dlijdix) );
            Kglobal(forAssy1[b2][1], forAssy1[b2][1]) -= (  cn*gap2 + cn*-dy*2.0*gap*(-dlijdix) );

            Fglobal(forAssy1[b1][0])    += Fx;
            Fglobal(forAssy1[b1][1])    += Fy;
            Fglobal(forAssy1[b2][0])    -= Fx;
            Fglobal(forAssy1[b2][1])    -= Fy;
        }
        else
        {
            //cout << " Contact " << b1 << " to " << b2 << " is inactive " << endl;
        }
      } // for(bb=0; bb<nCont; bb++)

    //cout << " kkkkkkkkkkk " << endl;

    //printVector(Fglobal);
    
    cout << " norm = " << Fglobal.norm() << endl;
    //ImmersedBodyObjects[bb]->SolveTimeStep();

    if(Fglobal.norm() < 1.0e-6)
      break;

    slnTemp = Kglobal.fullPivLu().solve(Fglobal);
    
    printVector(slnTemp);

    for(bb=0; bb<nPart; bb++)
      ImmersedBodyObjects[bb]->updateDisplacement(&(slnTemp(3*bb)));

    cout << " kkkkkkkkkkk " << endl;

    //ImmersedBodyObjects[bb]->updateIterStep();
    //cout << " kkkkkkkkkkk " << endl;
  } // for(iter=0; iter<10; iter++)



  for(bb=0; bb<nPart; bb++)
  {
    ImmersedBodyObjects[bb]->updatePointPositions() ;
  }

  IB_MOVED = true;

  //slnTempPrev2 = slnTempPrev;
  //slnTempPrev  = slnTemp;

  return;
}
*/


/*
void  HBSplineCutFEM::solveSolidProblem()
{
  IB_MOVED = false;
  
  cout << " Solving solid " << endl;


  int  nPart=ImmersedBodyObjects.size(), nCont=1, bb, b1, b2, doftemp, ii, jj, kk, iter, r1;
  double  xc1, xc2, yc1, yc2, Rad=0.1;
  double  l0=2.5*gridLEN[0]/nelem[0], gap, gap2, gap3, lam, lam, c, s, theta;
  MatrixXd  Klocal(3,3), Kcont(2,6), Kglobal, matM(3,3);
  VectorXd  Flocal(3), Fglobal, tdTemp, xc(nPart), yc(nPart);

  cout << " aaaaaaaaaaaaaaaa " << endl;

  double  cn=1.0e3, tol=1.0e-10, dx, dy, af, rhof=1.0, rhos=1.2*rhof, As, normTemp, Fx, Fy, lij;

  //l0 = 0.21;
  //xc[0] = 1.0; yc[0] = 6.8;

  xc[0] = 1.0; yc[0] = 6.8;
  xc[1] = 1.0; yc[1] = 7.2;

  doftemp = 3*nPart;

  cout << " doftemp  = " << doftemp << endl;
  cout << " l0 = " << l0 << endl;

  Kglobal.resize(doftemp, doftemp);
  Fglobal.resize(doftemp);


  tdTemp = ImmersedBodyObjects[0]->SolnData.td;

  af = tdTemp[2];
  
  cout << " af = " << af << endl;
  
  As = PI*Rad*Rad;

  matM.setZero();
  matM(0,0) = As*rhos;
  matM(1,1) = As*rhos;
  matM(2,2) = 0.5*PI*Rad*Rad*Rad*Rad*rhos;

  double  F0 = (rhos-rhof)*As*981.0;

  Klocal = tdTemp[5]*matM;
  
  vector<vector<int> >  forAssy1, contIBs;
  
  forAssy1.resize(nPart);
  kk=0;
  for(bb=0; bb<nPart; bb++)
  {
    forAssy1[bb].resize(3);
    for(ii=0; ii<3; ii++)
      forAssy1[bb][ii] = kk++;

    printVector(forAssy1[bb]);
  }

  contIBs.resize(nCont);
  for(bb=0; bb<nCont; bb++)
    contIBs[bb].resize(2);
  
  if( nCont > 0 )
  {
    contIBs[0][0] = 0; contIBs[0][1] = 1;
  }

  //slnTemp = slnTempPrev;

  for(iter=0; iter<4; iter++)
  {
    cout << " iter = " << iter << endl;
    
    Kglobal.setZero();
    Fglobal.setZero();

    slnTempCur = af*slnTemp + (1.0-af)*slnTempPrev;

    for(b1=0; b1<nPart; b1++)    // loop over the particles
    {
      Flocal = ImmersedBodyObjects[b1]->SolnData.forceCur - matM*ImmersedBodyObjects[b1]->SolnData.var1DotDotCur;
      Flocal(1) -= F0;

      for(ii=0;ii<3;ii++)
      {
        Fglobal[forAssy1[b1][ii]] += Flocal[ii];
      
        for(jj=0;jj<3;jj++)
        {
          Kglobal(forAssy1[b1][ii],forAssy1[b1][jj])  +=  Klocal(ii,jj);
        }
      }


      xc1 = xc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[0];
      yc1 = yc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[1];

      // for each particle check whether it is in contact with the boundary walls
      for(b2=0;b2<4;b2++)
      {
        switch(b2)
        {
          case 0:
            xc2 = origin[0] - Rad;
            yc2 = yc1;
        
          break;

          case 1:

            xc2 = origin[0] + gridLEN[0] + Rad;
            yc2 = yc1;

          break;

          case 2:
        
            xc2 = xc1;
            yc2 = origin[1] - Rad;
        
          break;

          case 3:

            xc2 = xc1;
            yc2 = origin[1] + gridLEN[1] + Rad;

          break;

          default:
          break;
        }

        dx = xc1-xc2;
        dy = yc1-yc2;
        
        lij = sqrt(dx*dx+dy*dy);

        theta = atan2(dy, dx);

        c = cos(theta);
        s = sin(theta);

        gap = 2.0*Rad + l0 - lij;

        gaph = xc1*c + yc1*s - xc2*c - yc2*s;
        
        Fx = cn*gap*c;
        Fy = cn*gap*s;

        if( gap > tol )
        {
            cout << " Contact " << b1 << " to " << b2 << " is active " << endl;

            Kglobal(forAssy1[b1][0], forAssy1[b1][0]) += (cn*c*c*af);
            Kglobal(forAssy1[b1][0], forAssy1[b1][1]) += (cn*c*s*af);
            Kglobal(forAssy1[b1][1], forAssy1[b1][0]) += (cn*c*s*af);
            Kglobal(forAssy1[b1][1], forAssy1[b1][1]) += (cn*s*s*af);

            Fglobal(forAssy1[b1][0])    += Fx;
            Fglobal(forAssy1[b1][1])    += Fy;
        }
        else
        {
            //cout << " Contact " << b1 << " to " << b2 << " is inactive " << endl;
        }
      } // for(b2=0;b2<4;b2++)
    } //     for(b1=0; b1<nPart; b1++)    // loop over the particles

    cout << " norm = " << Fglobal.norm() << endl;

      cout << " contact with other particles " << endl;

      // check contact with other particles
      for(bb=0; bb<nCont; bb++)
      {
        b1 = contIBs[bb][0];
        b2 = contIBs[bb][1];

        xc1 = xc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[0];
        yc1 = yc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[1];

        xc2 = xc[b2] + ImmersedBodyObjects[b2]->SolnData.var1Cur[0];
        yc2 = yc[b2] + ImmersedBodyObjects[b2]->SolnData.var1Cur[1];

        dx = xc1-xc2;
        dy = yc1-yc2;
        
        lij = sqrt(dx*dx+dy*dy);

        theta = atan2(dy, dx);

        c = cos(theta);
        s = sin(theta);

        gap = 2.0*Rad + l0 - lij;

        gap2 = gap*gap;
        gap3 = gap*gap2;

        Fx = cn*gap*c;
        Fy = cn*gap*s;
        cout << " dx = " << dx << endl;
        cout << " dy = " << dy << endl;
        cout << " gap = " << gap << endl;
        cout << " Fx = " << Fx << endl;
        cout << " Fy = " << Fy << endl;

        if( gap > tol)
        {
            cout << " Contact " << b1 << " to " << b2 << " is active " << endl;

            Kglobal(forAssy1[b1][0], forAssy1[b1][0]) += (cn*c*c*af);
            Kglobal(forAssy1[b1][0], forAssy1[b1][1]) += (cn*c*s*af);
            Kglobal(forAssy1[b1][0], forAssy1[b2][0]) -= (cn*c*c*af);
            Kglobal(forAssy1[b1][0], forAssy1[b2][1]) -= (cn*c*s*af);

            Kglobal(forAssy1[b1][1], forAssy1[b1][0]) += (cn*c*s*af);
            Kglobal(forAssy1[b1][1], forAssy1[b1][1]) += (cn*s*s*af);
            Kglobal(forAssy1[b1][1], forAssy1[b2][0]) -= (cn*c*s*af);
            Kglobal(forAssy1[b1][1], forAssy1[b2][1]) -= (cn*s*s*af);

            Kglobal(forAssy1[b2][0], forAssy1[b1][0]) -= (cn*c*c*af);
            Kglobal(forAssy1[b2][0], forAssy1[b1][1]) -= (cn*c*s*af);
            Kglobal(forAssy1[b2][0], forAssy1[b2][0]) += (cn*c*c*af);
            Kglobal(forAssy1[b2][0], forAssy1[b2][1]) += (cn*c*s*af);

            Kglobal(forAssy1[b2][1], forAssy1[b1][0]) -= (cn*c*s*af);
            Kglobal(forAssy1[b2][1], forAssy1[b1][1]) -= (cn*s*s*af);
            Kglobal(forAssy1[b2][1], forAssy1[b2][0]) += (cn*c*s*af);
            Kglobal(forAssy1[b2][1], forAssy1[b2][1]) += (cn*s*s*af);

            Fglobal(forAssy1[b1][0])    += Fx;
            Fglobal(forAssy1[b1][1])    += Fy;
            Fglobal(forAssy1[b2][0])    -= Fx;
            Fglobal(forAssy1[b2][1])    -= Fy;

        }
        else
        {
            //cout << " Contact " << b1 << " to " << b2 << " is inactive " << endl;
        }

      } // for(bb=0; bb<nCont; bb++)

    cout << " kkkkkkkkkkk " << endl;
    cout << " norm = " << Fglobal.norm() << endl;
    //ImmersedBodyObjects[bb]->SolveTimeStep();

    if(Fglobal.norm() < 1.0e-6)
      break;

    slnTemp = Kglobal.fullPivLu().solve(Fglobal);

    for(bb=0; bb<nPart; bb++)
      ImmersedBodyObjects[bb]->updateDisplacement(&(slnTemp(3*bb)));

    cout << " kkkkkkkkkkk " << endl;

    //ImmersedBodyObjects[bb]->updateIterStep();
    //cout << " kkkkkkkkkkk " << endl;
  } // for(iter=0; iter<10; iter++)



  for(bb=0; bb<nPart; bb++)
  {
    ImmersedBodyObjects[bb]->updatePointPositions() ;
  }

  IB_MOVED = true;

  //slnTempPrev2 = slnTempPrev;
  //slnTempPrev  = slnTemp;

  return;
}
*/



/*
void  HBSplineCutFEM::solveSolidProblem()
{
  IB_MOVED = false;
  // Lagrange multipliers
  cout << " Solving solid " << endl;


  int  nPart=ImmersedBodyObjects.size(), nCont=300, bb, b1, b2, doftemp, ii, jj, kk, iter, r1;
  double  xc1, xc2, yc1, yc2, Rad=0.1;
  double  gap, gap2, gap3, lam, c, s, theta;
  MatrixXd  Klocal(3,3), Kcont(2,6), Kglobal, matM(3,3);
  VectorXd  Flocal(3), Fglobal, tdTemp, xc(nPart), yc(nPart);

  double  cn=1.0, tol=1.0e-10, dx, dy, af, As, normTemp, lij, fact, dx2, dy2, dxdy;

  //double  rhof=1.0, rhos=1.01*rhof;
  double  rhof=0.998, rhos=1.002*rhof;

  double  l0=2.0*gridLEN[0]/nelem[0];

  //xc[0] = 1.0; yc[0] = 7.2;

  xc[0] = 1.0; yc[0] = 6.8;
  xc[1] = 1.0; yc[1] = 7.2;

  xc[0] = 0.8; yc[0] = 6.8;
  xc[1] = 1.2; yc[1] = 6.8;
  xc[2] = 0.8; yc[2] = 7.2;
  xc[3] = 1.2; yc[3] = 7.2;

  xc[0]  = 0.7; yc[0]  = 3.0;
  xc[1]  = 1.1; yc[1]  = 3.0;
  xc[2]  = 1.5; yc[2]  = 3.0;
  xc[3]  = 1.9; yc[3]  = 3.0;
  xc[4]  = 2.3; yc[4]  = 3.0;
  xc[5]  = 0.7; yc[5]  = 3.4;
  xc[6]  = 1.1; yc[6]  = 3.4;
  xc[7]  = 1.5; yc[7]  = 3.4;
  xc[8]  = 1.9; yc[8]  = 3.4;
  xc[9]  = 2.3; yc[9]  = 3.4;
  xc[10] = 0.7; yc[10] = 3.8;
  xc[11] = 1.1; yc[11] = 3.8;
  xc[12] = 1.5; yc[12] = 3.8;
  xc[13] = 1.9; yc[13] = 3.8;
  xc[14] = 2.3; yc[14] = 3.8;
  xc[15] = 0.7; yc[15] = 4.2;
  xc[16] = 1.1; yc[16] = 4.2;
  xc[17] = 1.5; yc[17] = 4.2;
  xc[18] = 1.9; yc[18] = 4.2;
  xc[19] = 2.3; yc[19] = 4.2;
  xc[20] = 0.7; yc[20] = 4.6;
  xc[21] = 1.1; yc[21] = 4.6;
  xc[22] = 1.5; yc[22] = 4.6;
  xc[23] = 1.9; yc[23] = 4.6;
  xc[24] = 2.3; yc[24] = 4.6;



  doftemp = 3*nPart+4*nPart+nCont;

  cout << " doftemp  = " << doftemp << endl;
  //cout << " l0 = " << l0 << endl;

  Kglobal.resize(doftemp, doftemp);
  Fglobal.resize(doftemp);


  tdTemp = ImmersedBodyObjects[0]->SolnData.td;

  af = tdTemp[2];
  
  //cout << " af = " << af << endl;
  
  As = PI*Rad*Rad;

  matM.setZero();
  matM(0,0) = As*rhos;
  matM(1,1) = As*rhos;
  matM(2,2) = 0.5*(As*rhos)*Rad*Rad;

  double  F0 = (rhos-rhof)*As*981.0;

  Klocal = tdTemp[5]*matM;
  
  vector<vector<int> >  forAssy1, contIBs;
  
  forAssy1.resize(nPart);
  kk=0;
  for(bb=0; bb<nPart; bb++)
  {
    forAssy1[bb].resize(3);
    for(ii=0; ii<3; ii++)
      forAssy1[bb][ii] = kk++;

    //printVector(forAssy1[bb]);
  }

  contIBs.resize(nCont);
  for(bb=0; bb<nCont; bb++)
    contIBs[bb].resize(2);
  
  if( nCont > 0 )
  {
    contIBs[0][0] = 0; contIBs[0][1] = 1;
    contIBs[1][0] = 0; contIBs[1][1] = 2;
    contIBs[2][0] = 0; contIBs[2][1] = 3;
    contIBs[3][0] = 1; contIBs[3][1] = 2;
    contIBs[4][0] = 1; contIBs[4][1] = 3;
    contIBs[5][0] = 2; contIBs[5][1] = 3;
    
    bb=0;
    for(ii=0;ii<nPart-1;ii++)
    {
      for(jj=(ii+1);jj<nPart;jj++)
      {
	contIBs[bb][0] = ii;
	contIBs[bb][1] = jj;
	bb++;
      }
    }
  }

  //slnTemp = slnTempPrev;

  for(iter=0; iter<20; iter++)
  {
    cout << " iter = " << iter << endl;
    
    Kglobal.setZero();
    Fglobal.setZero();

    slnTempCur = af*slnTemp + (1.0-af)*slnTempPrev;

    for(b1=0; b1<nPart; b1++)    // loop over the particles
    {
      Flocal = ImmersedBodyObjects[b1]->SolnData.forceCur - matM*ImmersedBodyObjects[b1]->SolnData.var1DotDotCur;
      Flocal(1) -= F0;

      for(ii=0;ii<3;ii++)
      {
        Fglobal[forAssy1[b1][ii]] += Flocal[ii];
      
        for(jj=0;jj<3;jj++)
        {
          Kglobal(forAssy1[b1][ii],forAssy1[b1][jj])  +=  Klocal(ii,jj);
        }
      }


      xc1 = xc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[0];
      yc1 = yc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[1];

      // for each particle check whether it is in contact with the boundary walls
      for(b2=0;b2<4;b2++)
      {
        switch(b2)
        {
          case 0:
            xc2 = origin[0] - Rad;
            yc2 = yc1;

          break;

          case 1:

            xc2 = origin[0] + gridLEN[0] + Rad;
            yc2 = yc1;

          break;

          case 2:
        
            xc2 = xc1;
            yc2 = origin[1] - Rad;
        
          break;

          case 3:

            xc2 = xc1;
            yc2 = origin[1] + gridLEN[1] + Rad;

          break;

          default:
          break;
        }

        dx = xc2-xc1;
        dy = yc2-yc1;
        
        dx2  = dx*dx;
        dy2  = dy*dy;
        dxdy = dx*dy;

        lij = sqrt(dx2+dy2);

        c = dx/lij;
        s = dy/lij;

        gap = 2.0*Rad + l0 - lij;

        //cout << " gap = " << dx << '\t' << dy << '\t' << gap << endl;

        r1 = 3*nPart + 4*b1 + b2;

        lam = slnTempCur[r1];

        fact = lam/lij/lij/lij;
        
        if( (lam + cn*gap) > tol)
        //if( gap > tol)
        {
            //cout << " Contact " << b1 << " to " << b2 << " is active " << endl;

            Kglobal(forAssy1[b1][0], forAssy1[b1][0]) -= (af*fact*dy2);
            Kglobal(forAssy1[b1][0], forAssy1[b1][1]) += (af*fact*dxdy);
            Kglobal(forAssy1[b1][0], r1)              += (af*dx/lij);

            Kglobal(forAssy1[b1][1], forAssy1[b1][0]) += (af*fact*dxdy);
            Kglobal(forAssy1[b1][1], forAssy1[b1][1]) -= (af*fact*dx2);
            Kglobal(forAssy1[b1][1], r1)              += (af*dy/lij);

            Kglobal(r1, forAssy1[b1][0])              += (af*dx/lij);
            Kglobal(r1, forAssy1[b1][1])              += (af*dy/lij);

            Fglobal(forAssy1[b1][0])    -= lam*dx/lij;
            Fglobal(forAssy1[b1][1])    -= lam*dy/lij;
            Fglobal(r1)                 -= gap;
        }
        else
        {
            //cout << " Contact " << b1 << " to " << b2 << " is inactive " << endl;
    
            Kglobal(r1,   r1)   += af;
            Fglobal(r1)         -= lam;
        }
      } // for(b2=0;b2<4;b2++)
    } //     for(b1=0; b1<nPart; b1++)    // loop over the particles

      cout << " contact with other particles " << endl;

      // check contact with other particles
      for(bb=0; bb<nCont; bb++)
      {
        b1 = contIBs[bb][0];
        b2 = contIBs[bb][1];

        xc1 = xc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[0];
        yc1 = yc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[1];

        xc2 = xc[b2] + ImmersedBodyObjects[b2]->SolnData.var1Cur[0];
        yc2 = yc[b2] + ImmersedBodyObjects[b2]->SolnData.var1Cur[1];

        dx = xc2-xc1;
        dy = yc2-yc1;
        
        dx2  = dx*dx;
        dy2  = dy*dy;
        dxdy = dx*dy;

        lij = sqrt(dx2+dy2);

        c = dx/lij;
        s = dy/lij;

        gap = 2.0*Rad + l0 - lij;

        r1 = 3*nPart + 4*nPart + bb;

        lam = slnTempCur[r1];

        //cout << " gap = " << dx << '\t' << dy << '\t' << gap << '\t' << lam << endl;

        fact = lam/lij/lij/lij;

        if( (lam + cn*gap) > tol)
        //if( gap > tol)
        {
            cout << " Contact " << b1 << " to " << b2 << " is active " << endl;
            cout << " r1 = " << r1 << endl;

            Kglobal(forAssy1[b1][0], forAssy1[b1][0]) -= (af*fact*dy2);
            Kglobal(forAssy1[b1][0], forAssy1[b1][1]) += (af*fact*dxdy);
            Kglobal(forAssy1[b1][0], forAssy1[b2][0]) += (af*fact*dy2);
            Kglobal(forAssy1[b1][0], forAssy1[b2][1]) -= (af*fact*dxdy);
            Kglobal(forAssy1[b1][0], r1)              += (af*dx/lij);

            Kglobal(forAssy1[b1][1], forAssy1[b1][0]) += (af*fact*dxdy);
            Kglobal(forAssy1[b1][1], forAssy1[b1][1]) -= (af*fact*dx2);
            Kglobal(forAssy1[b1][1], forAssy1[b2][0]) -= (af*fact*dxdy);
            Kglobal(forAssy1[b1][1], forAssy1[b2][1]) += (af*fact*dx2);
            Kglobal(forAssy1[b1][1], r1)              += (af*dy/lij);

            Kglobal(forAssy1[b2][0], forAssy1[b1][0]) += (af*fact*dy2);
            Kglobal(forAssy1[b2][0], forAssy1[b1][1]) -= (af*fact*dxdy);
            Kglobal(forAssy1[b2][0], forAssy1[b2][0]) -= (af*fact*dy2);
            Kglobal(forAssy1[b2][0], forAssy1[b2][1]) += (af*fact*dxdy);
            Kglobal(forAssy1[b2][0], r1)              -= (af*dx/lij);

            Kglobal(forAssy1[b2][1], forAssy1[b1][0]) -= (af*fact*dxdy);
            Kglobal(forAssy1[b2][1], forAssy1[b1][1]) += (af*fact*dx2);
            Kglobal(forAssy1[b2][1], forAssy1[b2][0]) += (af*fact*dxdy);
            Kglobal(forAssy1[b2][1], forAssy1[b2][1]) -= (af*fact*dx2);
            Kglobal(forAssy1[b2][1], r1)              -= (af*dy/lij);

            Kglobal(r1, forAssy1[b1][0])              += (af*dx/lij);
            Kglobal(r1, forAssy1[b1][1])              += (af*dy/lij);
            Kglobal(r1, forAssy1[b2][0])              -= (af*dx/lij);
            Kglobal(r1, forAssy1[b2][1])              -= (af*dy/lij);

            Fglobal(forAssy1[b1][0])    -= lam*dx/lij;
            Fglobal(forAssy1[b1][1])    -= lam*dy/lij;
            Fglobal(forAssy1[b2][0])    += lam*dx/lij;
            Fglobal(forAssy1[b2][1])    += lam*dy/lij;
            Fglobal(r1)                 -= gap;
        }
        else
        {
            //cout << " Contact " << b1 << " to " << b2 << " is inactive " << endl;
            //cout << " r1 = " << r1 << endl;

            Kglobal(r1,   r1)   += af;
            Fglobal(r1)         -= lam;
        }
      } // for(bb=0; bb<nCont; bb++)

    //cout << " kkkkkkkkkkk " << endl;
    //printVector(Fglobal);
    //cout << " \n\n\n\n " << endl;
    cout << " norm = " << Fglobal.norm() << endl;

    if(Fglobal.norm() < 1.0e-6)
      break;

    slnTemp += Kglobal.fullPivLu().solve(Fglobal);

    for(bb=0; bb<nPart; bb++)
      ImmersedBodyObjects[bb]->updateDisplacement(&(slnTemp(3*bb)));

    //cout << " kkkkkkkkkkk " << endl;

    //ImmersedBodyObjects[bb]->updateIterStep();
    //cout << " kkkkkkkkkkk " << endl;
  } // for(iter=0; iter<10; iter++)



  for(bb=0; bb<nPart; bb++)
  {
    ImmersedBodyObjects[bb]->updatePointPositions() ;
  }

  IB_MOVED = true;

  slnTempPrev2 = slnTempPrev;
  slnTempPrev  = slnTemp;


  return;
}
*/



/*
void  HBSplineCutFEM::solveSolidProblem()
{
  IB_MOVED = false;
  // penalty method
  cout << " Solving solid " << endl;


  int  nPart=ImmersedBodyObjects.size(), nCont=300, bb, b1, b2, doftemp, ii, jj, kk, iter, r1;
  double  xc1, xc2, yc1, yc2, Rad=0.1;
  double  gap, gap2, gap3, lam, c, s, theta;
  MatrixXd  Klocal(3,3), Kcont(2,6), Kglobal, matM(3,3);
  VectorXd  Flocal(3), Fglobal, tdTemp, xc(nPart), yc(nPart);

  double  cn=1.0e4, tol=1.0e-10, dx, dy, af, As, normTemp, lij, fact, dx2, dy2, dxdy;

  double  rhof=1.0, rhos=1.1*rhof;

  double  l0=2.0*gridLEN[0]/nelem[0];

  //xc[0] = 1.0; yc[0] = 7.2;

  xc[0] = 1.0; yc[0] = 6.8;
  xc[1] = 1.0; yc[1] = 7.2;

  xc[0] = 0.8; yc[0] = 6.8;
  xc[1] = 1.2; yc[1] = 6.8;
  xc[2] = 0.8; yc[2] = 7.2;
  xc[3] = 1.2; yc[3] = 7.2;

  xc[0]  = 0.7; yc[0]  = 3.0;
  xc[1]  = 1.1; yc[1]  = 3.0;
  xc[2]  = 1.5; yc[2]  = 3.0;
  xc[3]  = 1.9; yc[3]  = 3.0;
  xc[4]  = 2.3; yc[4]  = 3.0;
  xc[5]  = 0.7; yc[5]  = 3.4;
  xc[6]  = 1.1; yc[6]  = 3.4;
  xc[7]  = 1.5; yc[7]  = 3.4;
  xc[8]  = 1.9; yc[8]  = 3.4;
  xc[9]  = 2.3; yc[9]  = 3.4;
  xc[10] = 0.7; yc[10] = 3.8;
  xc[11] = 1.1; yc[11] = 3.8;
  xc[12] = 1.5; yc[12] = 3.8;
  xc[13] = 1.9; yc[13] = 3.8;
  xc[14] = 2.3; yc[14] = 3.8;
  xc[15] = 0.7; yc[15] = 4.2;
  xc[16] = 1.1; yc[16] = 4.2;
  xc[17] = 1.5; yc[17] = 4.2;
  xc[18] = 1.9; yc[18] = 4.2;
  xc[19] = 2.3; yc[19] = 4.2;
  xc[20] = 0.7; yc[20] = 4.6;
  xc[21] = 1.1; yc[21] = 4.6;
  xc[22] = 1.5; yc[22] = 4.6;
  xc[23] = 1.9; yc[23] = 4.6;
  xc[24] = 2.3; yc[24] = 4.6;

  doftemp = 3*nPart;

  cout << " doftemp  = " << doftemp << endl;
  //cout << " l0 = " << l0 << endl;

  Kglobal.resize(doftemp, doftemp);
  Fglobal.resize(doftemp);


  tdTemp = ImmersedBodyObjects[0]->SolnData.td;

  af = tdTemp[2];
  
  //cout << " af = " << af << endl;
  
  As = PI*Rad*Rad;

  matM.setZero();
  matM(0,0) = As*rhos;
  matM(1,1) = As*rhos;
  matM(2,2) = 0.5*(As*rhos)*Rad*Rad;

  double  F0 = (rhos-rhof)*As*981.0;

  Klocal = tdTemp[5]*matM;
  
  vector<vector<int> >  forAssy1, contIBs;
  
  forAssy1.resize(nPart);
  kk=0;
  for(bb=0; bb<nPart; bb++)
  {
    forAssy1[bb].resize(3);
    for(ii=0; ii<3; ii++)
      forAssy1[bb][ii] = kk++;

    //printVector(forAssy1[bb]);
  }

  contIBs.resize(nCont);
  for(bb=0; bb<nCont; bb++)
    contIBs[bb].resize(2);
  
  if( nCont > 0 )
  {
    contIBs[0][0] = 0; contIBs[0][1] = 1;
    contIBs[1][0] = 0; contIBs[1][1] = 2;
    contIBs[2][0] = 0; contIBs[2][1] = 3;
    contIBs[3][0] = 1; contIBs[3][1] = 2;
    contIBs[4][0] = 1; contIBs[4][1] = 3;
    contIBs[5][0] = 2; contIBs[5][1] = 3;

    bb=0;
    for(ii=0;ii<nPart-1;ii++)
    {
      for(jj=(ii+1);jj<nPart;jj++)
      {
	contIBs[bb][0] = ii;
	contIBs[bb][1] = jj;
	//cout << bb << '\t' << contIBs[bb][0] << '\t' << contIBs[bb][1] << endl;
	bb++;
      }
    }
  }



  for(iter=0; iter<15; iter++)
  {
    cout << " iter = " << iter << endl;
    
    Kglobal.setZero();
    Fglobal.setZero();

    slnTempCur = af*slnTemp + (1.0-af)*slnTempPrev;

    for(b1=0; b1<nPart; b1++)    // loop over the particles
    {
      Flocal = ImmersedBodyObjects[b1]->SolnData.forceCur - matM*ImmersedBodyObjects[b1]->SolnData.var1DotDotCur;

      //if( (b1==9) || (b1==10) || (b1 == 11) )
        //printVector(Flocal);

      Flocal(1) -= F0;

      for(ii=0;ii<3;ii++)
      {
        Fglobal[forAssy1[b1][ii]] += Flocal[ii];
      
        for(jj=0;jj<3;jj++)
        {
          Kglobal(forAssy1[b1][ii],forAssy1[b1][jj])  +=  Klocal(ii,jj);
        }
      }


      xc1 = xc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[0];
      yc1 = yc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[1];

      // for each particle check whether it is in contact with the boundary walls
      for(b2=0;b2<4;b2++)
      {
        switch(b2)
        {
          case 0:
            xc2 = origin[0] - Rad;
            yc2 = yc1;

          break;

          case 1:

            xc2 = origin[0] + gridLEN[0] + Rad;
            yc2 = yc1;

          break;

          case 2:
        
            xc2 = xc1;
            yc2 = origin[1] - Rad;
        
          break;

          case 3:

            xc2 = xc1;
            yc2 = origin[1] + gridLEN[1] + Rad;

          break;

          default:
          break;
        }

        dx = xc2-xc1;
        dy = yc2-yc1;
        
        dx2  = dx*dx;
        dy2  = dy*dy;
        dxdy = dx*dy;

        lij = sqrt(dx2+dy2);

        c = dx/lij;
        s = dy/lij;

        gap = 2.0*Rad + l0 - lij;

        //cout << " gap = " << dx << '\t' << dy << '\t' << gap << endl;

        fact = gap/lij/lij/lij;

        //if( (lam + cn*gap) > tol)
        if( gap > tol)
        {
            //cout << " Contact " << b1 << " to " << b2 << " is active " << endl;

            Kglobal(forAssy1[b1][0], forAssy1[b1][0]) += af*cn*(-gap/lij + dx*dx/lij/lij + fact*dx*dx);
            Kglobal(forAssy1[b1][0], forAssy1[b1][1]) += af*cn*(    0.0  + dx*dy/lij/lij + fact*dx*dy);

            Kglobal(forAssy1[b1][1], forAssy1[b1][0]) += af*cn*(    0.0  + dx*dy/lij/lij + fact*dx*dy);
            Kglobal(forAssy1[b1][1], forAssy1[b1][1]) += af*cn*(-gap/lij + dy*dy/lij/lij + fact*dy*dy);

            Fglobal(forAssy1[b1][0])    -= cn*gap*dx/lij;
            Fglobal(forAssy1[b1][1])    -= cn*gap*dy/lij;
        }
        else
        {
            //cout << " Contact " << b1 << " to " << b2 << " is inactive " << endl;
        }
      } // for(b2=0;b2<4;b2++)
    } //     for(b1=0; b1<nPart; b1++)    // loop over the particles

      cout << " contact with other particles " << endl;

      // check contact with other particles
      for(bb=0; bb<nCont; bb++)
      {
        b1 = contIBs[bb][0];
        b2 = contIBs[bb][1];

        xc1 = xc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[0];
        yc1 = yc[b1] + ImmersedBodyObjects[b1]->SolnData.var1Cur[1];

        xc2 = xc[b2] + ImmersedBodyObjects[b2]->SolnData.var1Cur[0];
        yc2 = yc[b2] + ImmersedBodyObjects[b2]->SolnData.var1Cur[1];

        dx = xc2-xc1;
        dy = yc2-yc1;
        
        dx2  = dx*dx;
        dy2  = dy*dy;
        dxdy = dx*dy;

        lij = sqrt(dx2+dy2);

        c = dx/lij;
        s = dy/lij;

        gap = 2.0*Rad + l0 - lij;

        //cout << " gap = " << dx << '\t' << dy << '\t' << gap << endl;

        fact = gap/lij/lij/lij;

        //if( (lam + cn*gap) > tol)
        if( gap > tol)
        {
            cout << " Contact " << b1 << " to " << b2 << " is active " << endl;
            //cout << " r1 = " << r1 << endl;

            Kglobal(forAssy1[b1][0], forAssy1[b1][0]) += af*cn*(-gap/lij + dx*dx/lij/lij + fact*dx*dx);
            Kglobal(forAssy1[b1][0], forAssy1[b1][1]) += af*cn*(    0.0  + dx*dy/lij/lij + fact*dx*dy);
            Kglobal(forAssy1[b1][0], forAssy1[b2][0]) += af*cn*( gap/lij - dx*dx/lij/lij - fact*dx*dx);
            Kglobal(forAssy1[b1][0], forAssy1[b2][1]) += af*cn*(    0.0  - dx*dy/lij/lij - fact*dx*dy);

            Kglobal(forAssy1[b1][1], forAssy1[b1][0]) += af*cn*(    0.0  + dx*dy/lij/lij + fact*dx*dy);
            Kglobal(forAssy1[b1][1], forAssy1[b1][1]) += af*cn*(-gap/lij + dy*dy/lij/lij + fact*dy*dy);
            Kglobal(forAssy1[b1][1], forAssy1[b2][0]) += af*cn*(    0.0  - dx*dy/lij/lij - fact*dx*dy);
            Kglobal(forAssy1[b1][1], forAssy1[b2][1]) += af*cn*( gap/lij - dy*dy/lij/lij - fact*dy*dy);

            Kglobal(forAssy1[b2][0], forAssy1[b1][0]) -= af*cn*(-gap/lij + dx*dx/lij/lij + fact*dx*dx);
            Kglobal(forAssy1[b2][0], forAssy1[b1][1]) -= af*cn*(    0.0  + dx*dy/lij/lij + fact*dx*dy);
            Kglobal(forAssy1[b2][0], forAssy1[b2][0]) -= af*cn*( gap/lij - dx*dx/lij/lij - fact*dx*dx);
            Kglobal(forAssy1[b2][0], forAssy1[b2][1]) -= af*cn*(    0.0  - dx*dy/lij/lij - fact*dx*dy);

            Kglobal(forAssy1[b2][1], forAssy1[b1][0]) -= af*cn*(    0.0  + dx*dy/lij/lij + fact*dx*dy);
            Kglobal(forAssy1[b2][1], forAssy1[b1][1]) -= af*cn*(-gap/lij + dy*dy/lij/lij + fact*dy*dy);
            Kglobal(forAssy1[b2][1], forAssy1[b2][0]) -= af*cn*(    0.0  - dx*dy/lij/lij - fact*dx*dy);
            Kglobal(forAssy1[b2][1], forAssy1[b2][1]) -= af*cn*( gap/lij - dy*dy/lij/lij - fact*dy*dy);

            Fglobal(forAssy1[b1][0])    -= cn*gap*dx/lij;
            Fglobal(forAssy1[b1][1])    -= cn*gap*dy/lij;
            Fglobal(forAssy1[b2][0])    += cn*gap*dx/lij;
            Fglobal(forAssy1[b2][1])    += cn*gap*dy/lij;
        }
        else
        {
            //cout << " Contact " << b1 << " to " << b2 << " is inactive " << endl;
            //cout << " r1 = " << r1 << endl;
        }
      } // for(bb=0; bb<nCont; bb++)

    //cout << " kkkkkkkkkkk " << endl;
    cout << " norm = " << Fglobal.norm() << endl;
    //printVector(Fglobal);
    //cout << "\n\n\n\n" << endl;

    if(Fglobal.norm() < 1.0e-6)
      break;

    slnTemp += Kglobal.fullPivLu().solve(Fglobal);

    //printVector(slnTemp);
    //cout << "\n\n\n\n" << endl;

    for(bb=0; bb<nPart; bb++)
      ImmersedBodyObjects[bb]->updateDisplacement(&(slnTemp(3*bb)));

    //cout << " kkkkkkkkkkk " << endl;


    //cout << " kkkkkkkkkkk " << endl;
  } // for(iter=0; iter<10; iter++)



  for(bb=0; bb<nPart; bb++)
  {
    ImmersedBodyObjects[bb]->updatePointPositions() ;
  }

  IB_MOVED = true;

  slnTempPrev2 = slnTempPrev;
  slnTempPrev  = slnTemp;


  return;
}
*/




void  HBSplineCutFEM::writeReadResult(int index, MyString &fileName)
{
  char fct[] = "HBSplineCutFEM::writeReadResult";
  cout << "  HBSplineCutFEM::writeReadResult  " << endl;
}



void  HBSplineCutFEM::computeTotalForce(int bb)
{
  if(DIM == 2)
    computeTotalForce2D(bb);
  else if(DIM == 3)
    computeTotalForce3D(bb);  
  
  return;
}




void  HBSplineCutFEM::computeTotalForce2D(int bb)
{
    if( ImmersedBodyObjects.size() == 0 )
    {
      cout << " HBSplineCutFEM::computeTotalForce(int bb) ... There are no immersed bodies to compute the force " << endl;
      return;
    }

    int  aa, ii, jj, kk, nlb, nlf, nlocal, gp, nGauss;

    nlf = 1;
    for(ii=0; ii<DIM; ii++)
      nlf *= (degree[ii]+1);

    double  dvol, PENALTY, detJ, fact;

    VectorXd  Nb, fluidAcceOnSolid;
    myPoint   vel, velSpec, acceSpec, trac, acceFluid, centroid;
    MatrixXd  stress(DIM, DIM);

    bool   axsy = (GeomData.FluidProps[2] == 1);
    double  rho = GeomData.FluidProps[3];
    double  mu  = GeomData.FluidProps[4];
    double  rhoSolid = 1.0;

    vector<double>  gausspoints, gaussweights;

    normal.setZero();
    param.setZero();
    geom.setZero();

    node* ndTemp;

    ImmersedIntegrationElement  *lme;
    myPoly*  poly;

          PENALTY = ImmersedBodyObjects[bb]->GetPenaltyParameter();

          nlb = ImmersedBodyObjects[bb]->ImmIntgElems[0]->pointNums.size();

          nGauss = (int) cutFEMparams[1];

          Nb.resize(nlb);

          getGaussPoints1D(nGauss, gausspoints, gaussweights);

          bool rigidflag = ImmersedBodyObjects[bb]->IsRigidBody() ;

          if(rigidflag)
          {
            totalForce.resize(6);
            fluidAcceOnSolid.resize(6);
          }
          else
          {
            aa = ImmersedBodyObjects[bb]->SolnData.var1.rows();
            totalForce.resize(aa);
            fluidAcceOnSolid.resize(aa);
          }

          totalForce.setZero();
          fluidAcceOnSolid.setZero();

          centroid = ImmersedBodyObjects[bb]->GetCentroid(2);

          //cout << " centroid = " << centroid[0] << '\t' << centroid[1] << endl;
          //centroid[0] = 73.0;
          //centroid[1] = 60.0;

        for(aa=0; aa<ImmersedBodyObjects[bb]->ImmIntgElems.size(); aa++)
        {
          lme  = ImmersedBodyObjects[bb]->ImmIntgElems[aa];
          poly = ImmersedBodyObjects[bb]->ImmersedFaces[aa];

          //cout << " aa = " << aa << '\t' << lme->isActive() << endl;
          if( lme->isActive() )
          {
            for(gp=0; gp<nGauss; gp++)
            {
              param[0] = gausspoints[gp];

              poly->computeBasisFunctions(param, geom, Nb, detJ);

              dvol  = gaussweights[gp] * detJ;

              //printf(" %12.6f,  %12.6f, %12.6f \n", Nb[0], Nb[1], dvol);

              if(axsy)
              {
                dvol *= 2.0*PI*geom[0];
              }

              lme->computeVelocity(Nb, velSpec);
              lme->computeAcceleration(Nb, acceSpec);

              //cout << " uuuuuuuuuuu " << endl;

              ndTemp = elems[findCellNumber(geom)];

              geometryToParametric(geom, param);

              poly->computeNormal(geom, normal);

              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);
              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", param[0], param[1], param[2], dvol);

              ndTemp->computeVelocityAndStress(param, vel, stress, acceFluid);

              //printf("FluidAcce[0] = %12.6f, FluidAcce[1] = %12.6f \n", FluidAcce[0], FluidAcce[1]);

              trac[0] =  stress(0,0)*normal[0] + stress(0,1)*normal[1] ;
              trac[1] =  stress(1,0)*normal[0] + stress(1,1)*normal[1] ;

              trac[0] += PENALTY*(velSpec[0]-vel[0]);
              trac[1] += PENALTY*(velSpec[1]-vel[1]);

              //trac[0] += (stagParams[3]*rho)*acceFluid[0];
              //trac[1] += (stagParams[3]*rho)*acceFluid[1];

              //trac[0] += stagParams[3]*rho*(acceFluid[0] - acceSpec[0]);
              //trac[1] += stagParams[3]*rho*(acceFluid[1] - acceSpec[1]);

              //trac[0] += stagParams[2]*(rhoSolid*acceSpec[0] - rho*acceFluid[0]);
              //trac[1] += stagParams[2]*(rhoSolid*acceSpec[1] - rho*acceFluid[1]);

              trac *= dvol;

              //cout << trac[0] << '\t' << trac[1] << endl;

              fact = rho*dvol;

              if(rigidflag)
              {
                totalForce[0] += trac[0];
                totalForce[1] += trac[1];

                geom[0] -= centroid[0] ;
                geom[1] -= centroid[1] ;

                totalForce[5] +=  (geom[0]*trac[1]-geom[1]*trac[0]) ;
              }
              else
              {
                for(ii=0;ii<nlb;ii++)
                {
                  jj = lme->pointNums[ii];
                  kk = 2*jj;

                  totalForce[kk]   += (Nb[ii]*trac[0]);
                  totalForce[kk+1] += (Nb[ii]*trac[1]);

                  fluidAcceOnSolid[kk]   += (Nb[ii]*fact)*acceFluid[0];
                  fluidAcceOnSolid[kk+1] += (Nb[ii]*fact)*acceFluid[1];

                  //fluidAcceOnSolid[kk]   += (Nb[ii]*fact)*(acceFluid[0] - acceSpec[0]);
                  //fluidAcceOnSolid[kk+1] += (Nb[ii]*fact)*(acceFluid[1] - acceSpec[1]);
                }
              }

              //cout << " uuuuuuuuuuu " << endl;
            }//for(gp=0...
          } // if( lme->isActive() )
        }//for(aa=0...

  ImmersedBodyObjects[bb]->fluidAcce = stagParams[4]*fluidAcceOnSolid;

  return;
}




void  HBSplineCutFEM::computeTotalForce3D(int bb)
{
    if( ImmersedBodyObjects.size() == 0 )
    {
      cout << " HBSplineCutFEM::computeTotalForce(int bb) ... There are no immersed bodies to compute the force " << endl;
      return;
    }

    int  aa, ii, jj, nlb, nlocal, gp, nGauss;

    int nlf = (degree[0]+1)*(degree[1]+1)*(degree[2]+1);

    double  dvol, PENALTY, detJ;

    VectorXd  Nb;
    myPoint   vel, velSpec, trac, ptTemp, centroid, FluidAcce;
    MatrixXd  stress(DIM, DIM);

    double  rho = GeomData.FluidProps[3];
    double  mu  = GeomData.FluidProps[4];

    vector<double>  gausspoints1, gausspoints2, gaussweights;

    normal.setZero();
    param.setZero();
    geom.setZero();

    node* ndTemp;

    ImmersedIntegrationElement  *lme;
    myPoly*  poly;

          PENALTY = ImmersedBodyObjects[bb]->GetPenaltyParameter();

          nlb = ImmersedBodyObjects[bb]->ImmIntgElems[0]->pointNums.size();

          Nb.resize(nlb);

          nGauss = (int) cutFEMparams[1];

          if(nlb == 3)
            getGaussPointsTriangle(nGauss, gausspoints1, gausspoints2, gaussweights);
          else
            getGaussPointsQuad(nGauss, gausspoints1, gausspoints2, gaussweights);

          //printVector(gausspoints);          printf("\n\n\n\n");
          //printVector(gaussweights);          printf("\n\n\n\n");

          bool rigidflag = ImmersedBodyObjects[bb]->IsRigidBody() ;

          if(rigidflag)
            totalForce.resize(6);
          else
          {
            aa = ImmersedBodyObjects[bb]->GetNumNodes()*DIM;
            totalForce.resize(aa);
          }

          totalForce.setZero();

          centroid.setZero();
          //centroid = ImmersedBodyObjects[bb]->GetCentroid(2);

          //cout << " ImmersedBodyObjects[bb]->ImmIntgElems.size() = " << ImmersedBodyObjects[bb]->ImmIntgElems.size() << endl;

        for(aa=0; aa<ImmersedBodyObjects[bb]->ImmIntgElems.size(); aa++)
        {
          lme  = ImmersedBodyObjects[bb]->ImmIntgElems[aa];
          poly = ImmersedBodyObjects[bb]->ImmersedFaces[aa];

          //cout << " aa = " << aa << '\t' << lme->isActive() << '\t' << bb << endl;
          if( lme->isActive() )
          {
            for(gp=0; gp<nGauss; gp++)
            {
              param[0] = gausspoints1[gp];
              param[1] = gausspoints2[gp];

              //cout << " uuuuuuuuuuu " << endl;

              poly->computeBasisFunctions(param, geom, Nb, detJ);

              dvol  = gaussweights[gp] * detJ;

              lme->computeVelocity(Nb, velSpec);

              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", geom[0], geom[1], geom[2], dvol);

              //if(aa == 15734)
              //{
                //ptTemp = poly->GetPoint(0);
                //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f \n", ptTemp[0], ptTemp[1], ptTemp[2]);
                //ptTemp = poly->GetPoint(1);
                //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f \n", ptTemp[0], ptTemp[1], ptTemp[2]);
                //ptTemp = poly->GetPoint(2);
                //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f \n", ptTemp[0], ptTemp[1], ptTemp[2]);
              //}

              ndTemp = elems[findCellNumber(geom)];

              geometryToParametric(geom, param);

              poly->computeNormal(geom, normal);

              //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %12.6f, \n", param[0], param[1], param[2], dvol);

              ndTemp->computeVelocityAndStress(param, vel, stress, FluidAcce);

              //cout << " uuuuuuuuuuu " << endl;

              trac[0] =  stress(0,0)*normal[0] + stress(0,1)*normal[1] + stress(0,2)*normal[2] ;
              trac[1] =  stress(1,0)*normal[0] + stress(1,1)*normal[1] + stress(1,2)*normal[2] ;
              trac[2] =  stress(2,0)*normal[0] + stress(2,1)*normal[1] + stress(2,2)*normal[2] ;

              trac[0] += PENALTY*(velSpec[0]-vel[0]);
              trac[1] += PENALTY*(velSpec[1]-vel[1]);
              trac[2] += PENALTY*(velSpec[2]-vel[2]);

              trac *= dvol;

              if(rigidflag)
              {
                totalForce[0] += trac[0];
                totalForce[1] += trac[1];
                totalForce[2] += trac[2];

                geom[0] -= centroid[0] ;
                geom[1] -= centroid[1] ;
                geom[2] -= centroid[2] ;

                totalForce[3] +=  (geom[1]*trac[2]-geom[2]*trac[1]) ;
                totalForce[4] +=  (geom[2]*trac[0]-geom[0]*trac[2]) ;
                totalForce[5] +=  (geom[0]*trac[1]-geom[1]*trac[0]) ;
              }
              else
              {
                for(ii=0;ii<nlb;ii++)
                {
                  jj = lme->pointNums[ii];
                  jj = jj*3;

                  totalForce[jj]   += (Nb[ii]*trac[0]);
                  totalForce[jj+1] += (Nb[ii]*trac[1]);
                  totalForce[jj+2] += (Nb[ii]*trac[2]);
                }
              }

              //cout << " uuuuuuuuuuu " << endl;
            }//for(gp=0...
          } // if( lme->isActive() )
        }//for(aa=0...

    //printf("Force in X-direction = %12.6f \n", totalForce[0]);
    //printf("Force in Y-direction = %12.6f \n", totalForce[1]);
    //printf("Force in Z-direction = %12.6f \n", totalForce[2]);

  return;
}











