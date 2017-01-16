
#include "FunctionsEssGrp.h"
#include "MyString.h"



char prgUpdateDisplay(bool wait)
{
  essGrpCopyPixmap();
  essGrpUpdateDisplay();

  if (!wait) return '\n';

  MyString ch;

  ch.inputKeepIfReturn();

  if (ch.length() > 0) return ch.asCharArray()[0];

  return '\n';
}


