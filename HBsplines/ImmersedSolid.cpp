
#include "ImmersedSolid.h"
#include "SolutionData.h"
#include "GeomDataLagrange.h"
#include "ImmersedIntegrationElement.h"
#include "FunctionsProgram.h"
#include "MpapTime.h"
#include "TimeFunction.h"

extern MpapTime mpapTime;
extern List<TimeFunction> timeFunction;



ImmersedSolid::ImmersedSolid()
{
  id = count++;

  nImmInt = 0;
  PENALTY = 1.0e2;
  BC_ENFORCE_TYPE = BC_ENFORCE_TYPE_PENALTY;
  
  isNitsche = false;
  NitscheFact = 1.0;

  solver = NULL;

  uGrid    =  vtkSmartPointer<vtkUnstructuredGrid>::New();

  polyDataVTK =  vtkSmartPointer<vtkPolyData>::New();

  selectEnclosedPoints   =   vtkSmartPointer<vtkSelectEnclosedPoints>::New();
  selectEnclosedPoints2  =   vtkSmartPointer<vtkSelectEnclosedPoints>::New();

  selectEnclosedPoints->CheckSurfaceOn();
  //selectEnclosedPoints->InsideOutOff();
  selectEnclosedPoints->SetTolerance(0.000001);

  totalForce.resize(6);
  totalForce.setZero();
}




ImmersedSolid::~ImmersedSolid()
{
  for(vector<ImmersedIntegrationElement*>::iterator pObj = ImmIntgElems.begin();
       pObj != ImmIntgElems.end(); ++pObj)
  {
    delete *pObj; // Note that this is deleting what pObj points to, which is a pointer
  }
  ImmIntgElems.clear(); // Purge the contents so no one tries to delete them again

  if(solver != NULL)
    delete solver;
  
  --count;
}



void ImmersedSolid::setTimeParam()
{
  //SolnData.setTimeIncrementType(4);

  SolnData.setTimeParam();
  
  return;
}



void ImmersedSolid::timeUpdate()
{
  firstIter = true;
  SolnData.firstIter = firstIter;
  localStiffnessError = 0;
  iterCount = 1;

  //filecount++;

  fluidAccePrev = fluidAcce;
  
  SolnData.timeUpdate();

  //SolnData.var1 = SolnData.var1Extrap;
  //SolnData.var1.setZero();

  //SolnData.var1  = SolnData.var1Prev;
  //SolnData.var1 += mpapTime.dt * SolnData.var1DotPrev;
  //SolnData.var1 += (0.5*mpapTime.dt * mpapTime.dt)*SolnData.var1DotDotPrev;

  if(STAGGERED)
    updateIterStep();

  return;
}



void ImmersedSolid::updateIterStep()
{
  SolnData.updateIterStep();

  updatePointPositions();

  return;
}



void ImmersedSolid::reset()
{
  SolnData.reset();

  return;
}




bool ImmersedSolid::converged()
{
  tol=1.0e-6;

  if (rNorm < tol && localStiffnessError == 0)
  {
    return true;
  }

  return false;
}




bool ImmersedSolid::diverging(double factor)
{
    if(rNormPrev > -0.1 && (rNorm / rNormPrev) > factor)
      return true;

    if(localStiffnessError != 0)
      return true;

    if(prgNAN(rNorm))
      return true;

    return false;
}




void  ImmersedSolid::setImmersedIntegrationElements(vector<vector<int> >& datatemp)
{
    int ii=0, jj=0, ind=0;
    int kk = datatemp[0].size()-1 ;

    //cout << " datatemp.size() = " << datatemp.size() << '\t' << kk << endl;

    nImmInt = datatemp.size();

    ImmersedIntegrationElement* lme;

    for(ii=0; ii<datatemp.size(); ii++)
    {
      lme = new ImmersedIntegrationElement();

      lme->SolnData = &(SolnData);
      lme->GeomDataLag = &(GeomData);
      lme->setDimension(DIM);
      lme->prepareElemData();

      //cout << " ii " << ii << endl;
      for(jj=0; jj<kk; jj++)
      {
        ind = datatemp[ii][1+jj] ;
        lme->pointNums.push_back(ind);
      }

      lme->initialiseDOFvalues();

      if( datatemp[ii][0] == 0 )
      {
        lme->turnIsActiveOFF();
      }

      ImmIntgElems.push_back(lme);
    }

    //cout << " nImmInt " << nImmInt << endl;

    return;
}



void  ImmersedSolid::setDataForOutput(vector<vector<int> >& vectemp)
{
    OutputData.resize(vectemp.size());
    for(int ii=0; ii<vectemp.size(); ii++)
    {
      OutputData[ii] = vectemp[ii];
    }

    return;
}





void ImmersedSolid::setImmersedElemActiveFlag(vector<int>& datatemp)
{
    return;
}





void  ImmersedSolid::computeCentroid(int index)
{
    centroid.setZero();
  
    switch(index)
    {
      case 0: // original configuration

        for(int aa=0; aa<nNode; aa++)
          centroid += GeomData.NodePosOrig[aa];

      break;

      case 1: // configuration at t_{n+1}

        for(int aa=0; aa<nNode; aa++)
          centroid += GeomData.NodePosNew[aa];

      break;

      case 2: // configuration at t_{n+af}

        for(int aa=0; aa<nNode; aa++)
          centroid += GeomData.NodePosCur[aa];

      break;
    }

    centroid /= nNode;

    bool printOut = 0;
    if(printOut)
    {
      printf("\n  Centroid of the immersed body ... %5d \n", id);
      printf(" X-dir \t Y-dir \t Z-dir \n ");
      if(DIM == 2)
        printf(" %12.6f \t  %12.6f \t  %12.6f \n", centroid[0], centroid[1], 0.0);
      else
        printf(" %12.6f \t  %12.6f \t  %12.6f \n", centroid[0], centroid[1], centroid[2]);
    }

    return;
}




void  ImmersedSolid::computeAABB(int index)
{
    bbox.initialize();

    int  ii=0, aa=0;
    double  val=0.0;

    switch(index)
    {
      case 0: // original configuration

        for(aa=0; aa<nNode; aa++)
        {
          for(ii=0; ii<DIM; ii++)
          {
            val = GeomData.NodePosOrig[aa][ii];

            if( val < bbox.minBB[ii]  )
              bbox.minBB[ii] = val ;

            if( val > bbox.maxBB[ii]  )
              bbox.maxBB[ii] = val ;
          }
        }

      break;

      case 1: // configuration at t_{n+1}

        for(aa=0; aa<nNode; aa++)
        {
          for(ii=0; ii<DIM; ii++)
          {
            val = GeomData.NodePosNew[aa][ii];

            if( val < bbox.minBB[ii]  )
              bbox.minBB[ii] = val ;

            if( val > bbox.maxBB[ii]  )
              bbox.maxBB[ii] = val ;
          }
        }

      break;

      case 2: // configuration at t_{n+af}

        for(aa=0; aa<nNode; aa++)
        {
          for(ii=0; ii<DIM; ii++)
          {
            val = GeomData.NodePosCur[aa][ii];

            //cout << aa << '\t' << ii << '\t' << val << endl;

            if( val < bbox.minBB[ii]  )
              bbox.minBB[ii] = val ;

            if( val > bbox.maxBB[ii]  )
              bbox.maxBB[ii] = val ;
          }
        }

      break;
    }

    //bbox.printSelf();
  
    for(ii=0; ii<ImmersedFaces.size(); ii++)
      ImmersedFaces[ii]->computeAABB();

    return;
}



void  ImmersedSolid::computeTotalForce()
{
    totalForce.resize(6);
    totalForce.setZero();
    if(isBoundaryConditionTypeLagrange())
    {
      computeCentroid(1);

      for(int aa=0;aa<ImmIntgElems.size();aa++)
      {
        //cout << " aa " << aa << endl;
        ImmIntgElems[aa]->integrateForceAndMoment(totalForce, centroid);
        //printVector(vectemp);
      }
    }

    bool printOut = 0;
    if(printOut)
    {
      printf("\n  Sum of the forces on immersed body ... %5d \n", id);
      printf(" X-dir \t Y-dir \t Z-dir \n ");
      printf(" %12.6f \t  %12.6f \t  %12.6f \n\n\n", totalForce[0], totalForce[1], totalForce[2]);
      printf("\n  Sum of the moments on immersed body ... %5d \n", id);
      printf(" X-dir \t Y-dir \t Z-dir \n ");
      printf(" %12.6f \t  %12.6f \t  %12.6f \n\n\n", totalForce[3], totalForce[4], totalForce[5]);
    }

  return;
}



void ImmersedSolid::initialise_solid_state()
{
    SolnData.var1 = SolnData.var1Prev + mpapTime.dt*SolnData.var1DotPrev ;
    SolnData.var1Dot = SolnData.var1DotPrev;
  
    updateIterStep();
  
    SolnData.var1PrevIter = SolnData.var1;

    return;
}






