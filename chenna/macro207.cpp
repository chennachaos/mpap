
#include "Macro.h"
#include "DomainTree.h"
#include "HBSplineFEM.h"
#include "HBSplineCutFEM.h"
#include "StandardFEM.h"



extern DomainTree domain;


int macro207(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "pout";
    macro.type = "chen";
    macro.what = "print data onto screen";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    macro.db.frameButtonBox();

    macro.db.addRadioBox("*ufull","Uinit","resi","stif","CnAr","reac","kvts","cntr","mpat","bfns","elIV");
    
    macro.db.frameRadioBox();

    macro.db.frameButtonBox();

    macro.db.addTextField(" patch = ",1);

    macro.db.addTextField(" elem # = ",1);

    macro.db.frameButtonBox();


    // and other stuff

    return 0;	  
  }
//--------------------------------------------------------------------------------------------------

 // std::cout << "          " << macro << "\n\n";

  int  type, id, val, patch, elenum;

  type    = roundToInt(macro.p[0]);
  id      = roundToInt(macro.p[1]) - 1;
  val     = roundToInt(macro.p[2]);
  patch   = roundToInt(macro.p[3])-1;
  elenum  = roundToInt(macro.p[4])-1;

    //if(type == 26)
    //{
      //if(val<11)
        //isogeometricFEM(domain(type,id)).printData(val, patch);
      //else
        //isogeometricFEM(domain(type,id)).printElemInvVars(patch, elenum);
    //}

    if(type == 27)
    {
      if(val<11)
        hbsplineFEM(domain(type,id)).printData(val, patch);
      else
        printf("\n hbsplineFEM(domain(type,id)).printData .... Index out of range \n");
    }

    if(type == 28)
    {
      if(val<11)
        hbsplineCutFEM(domain(type,id)).printData(val, patch);
      else
        printf("\n hbsplineFEM(domain(type,id)).printData .... Index out of range \n");
    }

    //if(type == 30)
    //{
      //if(val<11)
        //hbscutFEMElasticity(domain(type,id)).printData(val, patch);
      //else
        //printf("\n hbsplineFEM(domain(type,id)).printData .... Index out of range \n");
    //}

//--------------------------------------------------------------------------------------------------
  return 0;  
}

