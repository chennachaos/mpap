
#include "Plot.h"
#include "DAWindow.h"


extern Plot plot;





void prgPlotMatrixPattern(int nRow, int nCol, 
                          VectorBase<int> &row, VectorBase<int> &col,
                          char *fileName)
{
  DAWindow mtxPlot;

  double rdc = double(nRow)/double(nCol), r0, c0, dr, dc;

  if (!plot) plot.adjustToNewSize();

  if (plot.hdw < rdc)
  {
    dr = plot.h * .9;
    dc = plot.h / rdc * .9;
    r0 = plot.h * .05;
    c0 = .5 * (plot.w - dc);
  }
  else
  {
    dr = plot.w * rdc * .9;
    dc = plot.w * .9;
    r0 = .5 * (plot.h - dr);
    c0 = plot.w * .05;
  }
  //cout << int(c0) << "," << int(r0) << "; " << int(dc) << "," << int(dr) << "\n";

  mtxPlot.setup(int(c0),int(r0),int(dc),int(dr),double(nCol),double(nRow));

  mtxPlot.wipe();

  mtxPlot.frame();

  int i;

  dc /= double(nCol);
  dr /= double(nRow);

  for (i=0; i<row.n; i++) mtxPlot.fillRectangle(2,double(col[i]),double(row[i]),-1.,-1.,true);

  COUT << nRow << " x " << nCol << " matrix pattern plotted.\n";
  COUT << "matrix fill in: " << row.n << " / " << nRow * nCol 
       << " = " << double(row.n)/double(nRow*nCol)*100. << "%\n\n";

  return;
}


