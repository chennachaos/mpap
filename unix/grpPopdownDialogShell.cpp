
#include <Xm/Xm.h>


void grpPopdownDialogShell(Widget manager)
{
  XtUnmanageChild(manager);

  XtDestroyWidget(XtParent(manager));

  return;
}


