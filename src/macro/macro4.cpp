
#include "Macro.h"
#include "MacroQueue.h"
#include "If.h"
#include "DomainTree.h"
#include "FunctionsEssGrp.h"
#include "RunControl.h"
#include "MpapTime.h"
//#include "MicroCell.h"



extern MpapTime   mpapTime;
extern MacroQueue macroQueue;
extern RunControl runCtrl;
extern bool       noGUI;
extern DomainTree domains;
extern double     globalMaxIncrement;


int macro4(Macro &macro)
{
  if (!macro) 
  { 
    macro.name = "if";
    macro.type = "ctrl";
    macro.what = "begin if statement";

    macro.sensitivity[INTER] = true;
    macro.sensitivity[BATCH] = true;

    macro.db.selectDomain();

    macro.db.stringList("*user", "conv", "divg", "dtmn", "mous", "cinc", "poor", "merr");

    macro.db.addTextField("tol = ",0,10);

    return 0;
  }
//--------------------------------------------------------------------------------------------------

  MyString ch; 
  char   *yes[] = YESS, *no[] = NOO;
  bool   cond = false;
  char   *condStrg[] = { "user", "conv", "divg", "dtmn", "mous", "cinc", "poor", NULL };
  int    i    = roundToInt(macro.p[0]),
         type = roundToInt(macro.p[1]),
         id   = roundToInt(macro.p[2]) - 1,
         substepped = 0;

  double tol  = std::abs(macro.p[3]);

  If     *iff = &(macroQueue.iff[i]);

  switch (macro.strg.which(condStrg)) 
  {
    case  0: // user

             if (noGUI) 
             {
               while (ch.which(yes)<0 && ch.which(no)<0) 
               {
                 COUT << "\n";
                 COUT << "'if,user': condition = true ? (y/n) [y] ";

                 ch.free().append("y").inputKeepIfReturn();
               }
               COUT << "\n";

               if (ch.which(yes) > -1) cond = true;
             } 
             else 
             {
               COUT << "\n";
               COUT << "'if,user':   left mouse button -> true\n";
               COUT << "        any other mouse button -> false\n\n";

               //essGrpSetSensAllButDrawingArea(false);
               runCtrl.fixStatus(PRESSMOUSE);

               //if (essGrpWaitForMouseButtonPressed() == 1) cond = true; 
               //essGrpSetSensAllButDrawingArea(true);
               runCtrl.freeFixedStatus();
             }
        break;

    case  1: // conv

            cond = domain(type,id).converged();

        break;

    case  2: // divg

            cond = domain(type,id).diverging(macro.p[3]);

        break;

    case  3: // dtmn

            cond = (mpapTime.dt < macro.p[3]);

        break;

    case  4: // mous

            //cond = (essGrpWasMouseButtonPressed() != 0);

        break;

    case  5: // cinc

            COUT << "globalMaxIncrement = "; printf("%11.4e\n",globalMaxIncrement);

            cond = (globalMaxIncrement < tol);

            globalMaxIncrement = 0.;

        break;

    case  6: // poor

            cond = (!domain(type,id).elemSizeRatioOK(tol));

	     break;

    case  7: // merr

            //for (i=0; i<domain.nDomainOfType(MICROCELL); i++)

            //{ if(domain(MICROCELL,i).microError()==1) substepped++; }

            //cond =(intDiv(substepped*100,i)>roundToInt(macro.p[3]));

        break;

    default: // error

            COUT << "unknown condition! (->TRUE)\n\n";

            cond = true;

        break;
  }

  if (!cond) return iff->els + 2;

//--------------------------------------------------------------------------------------------------
  return 0;  
}

