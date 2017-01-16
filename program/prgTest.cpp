

#include "MathVector.h"
#include "List.h"

#include "MyString.h"

//#include "RunControl.h"


void prgTest(void)
{
  cout << "\n  Put the stuff to be tested into prgTest.cpp.\n\n";



//  List<MyString> stringList;

//  Vector<double> x, y;

  Vector<int> vi, b;

  for (int i=0; i<10; i++) b.append(i+1);

  cout << b << "\n";

  b.reverseOrder();

  cout << b << "\n";

  //for (int i=9; i>-1; i--) cout << b[i] << "\n";

  b.move(9,0);

  cout << b << "\n";

  b.del(7);

  cout << b << "\n";

  b.del(7);

  cout << b << "\n";

  //cout << b[8] << "\n";

  cout << b.n << "\n";

  cout << b.lastCoeff() << "\n";

  return;

  int data[1000];

  for (int i=0; i<1000; i++) { data[i] = i+1; vi.append(i+1); }

  for (int i=0;   i<970; i+=3) { simpleSwap(data,i,i+10);   vi.swap(i,i+10);   }
  for (int i=0;   i<970; i+=4) { simpleSwap(data,i,i+20);   vi.swap(i,i+20);   }
  for (int i=999; i>30;  i-=7) { simpleSwap(data,i,i-8);    vi.swap(i,i-8);    }
  for (int i=1;   i<200; i+=9) { simpleSwap(data,i,1000-i); vi.swap(i,1000-i); }

  cout << vi << "\n";

  vi.quickSort(true,data);

//  for (int i=0; i<1000; i++) if (i+1 != vi[i]) cout << vi[i] << "!!!\n";

  cout << vi << "\n";

  for (int i=0; i<1000; i++) if (data[i] != vi[i]) cout << "HELLO!!!!!\n\n";
 
  return;

  vi.swap(1,9);

  cout << vi << "\n";

  vi.swap(0,3);

  cout << vi << "\n";

  vi.swap(9,0);

  cout << vi << "\n";

  vi.swap(vi.firstItem().next->prev,vi.firstItem().next->next->next->next);

  cout << vi << "A\n";

  //vi.swap(1,0);

  //cout << vi << "\n";


  vi.quickSort(true);

  cout << vi << "B\n";

  b = vi;

  //b = *((VectorBase<int>*)&vi);

  cout << b << "C\n";

  cout << "\n\n";

  return;
}

