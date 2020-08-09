

#include "ShowProgress.h"







template<typename Type> void MatrixSparse<Type>::quickSort(bool rowCompFlg, bool overwriteFlg, 
                                                           VectorArray<int> &perm)
{
  Vector<int> tmp1, tmp2;

  VectorArray<int> tmp3;

  int i, j, m = row.n, *P;

  perm.setDim(m); P = perm.x;

  for (i=0; i<m; i++) perm.x[i] = i;

  if (rowCompFlg) tmp1 = row; else tmp1 = col;

  //cout << " sorting rows ...   ";

  tmp1.quickSort(false,perm.x);

  //cout << " done.\n\n";

  if (rowCompFlg) { tmp3 = col; for (i=0; i<m; i++) tmp2.append(tmp3[*P++]); }
                                            
  else            { tmp3 = row; for (i=0; i<m; i++) tmp2.append(tmp3[*P++]); }

  tmp3.free();

  //cout << " sorting columns ...   ";

  i = 0; 
  while (i < m)
  {
    j = i++;
    while (i < m) if (tmp1[i] == tmp1[i-1]) i++; else break;

    tmp2.quickSort(false,perm.x,j,i-1);

    //cout << " column " << i << " sorted.\n";
  }

  //cout << " done.\n\n";

  if (overwriteFlg)
  {
    if (rowCompFlg) { row.takeOver(tmp1); col.takeOver(tmp2); }
    else            { col.takeOver(tmp1); row.takeOver(tmp2); }
  }

  return;
}













template<typename Type> void MatrixSparse<Type>::copySort(bool rowCompFlg, bool overwriteFlg, 
                                                          bool symFlg, VectorArray<int> &perm)
{
  if (!rowCompFlg)
  {
    std::cout << "   MatrixSparse<Type>::copySort: cannot yet handle rowCompFlg == false!\n\n";
    exit(1);
  }
  if (!overwriteFlg)
  {
    std::cout << "   MatrixSparse<Type>::copySort: cannot yet handle overwriteFlg == false!\n\n";
    exit(1);
  }

  ListArray< Vector<int> > tmpRow, tmpCol, tmpRowP, tmpColP;

  int c, i, j, k, m = row.n, z = 0, *P;

  ShowProgress<int> progress;

  // allocate memory

  tmpRow.setDim(this->nRow);
  tmpCol.setDim(this->nCol);

  tmpRowP.setDim(this->nRow);
  tmpColP.setDim(this->nCol);

  // initialise permutation vector

  perm.setDim(m); P = perm.x;

  for (i=0; i<m; i++) perm.x[i] = i;

  // sorting for columns

  VectorCoeff<int> *nvc, *v = &col.firstItem(), *w = &row.firstItem();

  progress.init("            sorting for columns:   ","%8d",m,true);
  for (k=0; k<m; k++)
  {
    progress.show(k+1);
    tmpCol [v->val-1].append(w->val-1);
    tmpColP[v->val-1].append(P[k]);
    v = (VectorCoeff<int>*)(v->next);
    w = (VectorCoeff<int>*)(w->next);
  }

  // sorting for rows

  progress.init("            sorting for rows:      ","%8d",tmpCol.n,true);
  for (i=0; i<tmpCol.n; i++)
  {
    progress.show(i+1);

    k = tmpCol[i].n;

    v = &tmpCol [i].firstItem();
    w = &tmpColP[i].firstItem();

    for (j=0; j<k; j++)
    {
      tmpRow [v->val].append(i);
      tmpRowP[v->val].append(w->val);
      v = (VectorCoeff<int>*)(v->next);
      w = (VectorCoeff<int>*)(w->next);
    }
  }
  tmpColP.free();

  // if required make symmetric

  if (symFlg)
  {
    // regenerate tmpCol (now with ordered rows in order)

    for (i=0; i<tmpCol.n; i++) tmpCol[i].free();

    progress.init("            regenerating columns:  ","%8d",tmpRow.n,true);
    for (i=0; i<tmpRow.n; i++)
    {
      progress.show(i+1);

      k = tmpRow[i].n;

      v = &tmpRow[i].firstItem();

      for (j=0; j<k; j++)
      {
        tmpCol[v->val].append(i);
        v = (VectorCoeff<int>*)(v->next);
      }
    }

    // identify and insert zero entries

    progress.init("            inserting zero entries:","%8d",tmpRow.n,true);
    for (i=0; i<tmpRow.n; i++)
    {
      progress.show(i+1);

      c = 0;
      j = 0;

      v = &tmpCol[i].firstItem();
      w = &tmpRow[i].firstItem();

      //cout << tmpCol[i].n << " - " << tmpRow[i].n << " = " << tmpCol[i].n - tmpRow[i].n << "\n";

      while (c < tmpCol[i].n)
      {
        k = v->val; v = (VectorCoeff<int>*)(v->next); c++;

        while (j < tmpRow[i].n)

          if (w->val < k) { w = (VectorCoeff<int>*)(w->next); j++; } else break;

        //cout << i << ": " << w->val << "," << k << "\n";

        if (j == tmpRow[i].n || w->val != k)
        {
          nvc           = new VectorCoeff<int>;
          nvc->val      = k;
          nvc->next     = (ListItem*)w;
          w->prev->next = (ListItem*)nvc;
          nvc->prev     = w->prev;
          w->prev       = (ListItem*)nvc;
          tmpRow[i].n++;

          tmpRowP[i].insert(perm.n+z++,j++);
        }
      }
    }

    // adjust memory allocation

    perm.setDim(perm.n + z);
    for (i=0; i<z; i++) x.append(0);
  }

  tmpCol. free();

  // regenerate row and col

  row.free();
  col.free();

  k = 0;

  progress.init("            regenerating matrix:   ","%8d",tmpRow.n,true);
  for (i=0; i<tmpRow.n; i++)
  {
    progress.show(i+1);

    m = tmpRow[i].n;

    v = &tmpRow [i].firstItem();
    w = &tmpRowP[i].firstItem();

    for (j=0; j<m; j++)
    {
      row.append(i+1);
      col.append(v->val+1);
      perm[k++] = w->val;
      v = (VectorCoeff<int>*)(v->next);
      w = (VectorCoeff<int>*)(w->next);
    }
  }
  cout << "\n";
  if (symFlg)
    cout << "            matrix pattern is now symmetric: " << z << " zero entries added.\n\n";

  tmpRow.free();
  tmpRowP.free();

  return;
}














template<typename Type> int MatrixSparse<Type>::makeSymmetric(bool rowCompFlg,
                                                              VectorArray<int> &perm)
{
  if (!rowCompFlg)
  {
    std::cout << "   MatrixSparse<Type>::makeSymmetric: cannot yet handle rowCompFlg == false!\n\n";
    exit(1);
  }

  int c, i, j, k, m = row.n;

  ListArray< Vector<int> > tmpCol, tmpRow, tmpRowP;

  tmpCol.setDim(this->nCol);
  tmpRow.setDim(this->nRow);
  tmpRowP.setDim(this->nRow);

  VectorCoeff<int> *nvc, *v = &col.firstItem(), *w = &row.firstItem();

  i = 0;

  // generate tmpCol and tmpRow

  for (k=0; k<m; k++)
  {
    if (i > w->val)
    {
      std::cout << "   MatrixSparse<Type>::makeSymmetric: matrix is not row compressed!\n\n";
      exit(1);
    }
    else if (w->val - i > 1) 
    {
      std::cout << "   MatrixSparse<Type>::makeSymmetric: empty row?!\n\n";
      exit(1);
    }
    i = w->val;
    tmpCol[v->val-1].append(i);
    tmpRow[i-1].append(v->val);
    tmpRowP[i-1].append(perm.x[k]);
    v = (VectorCoeff<int>*)(v->next);
    w = (VectorCoeff<int>*)(w->next);
  }

  // identify and insert zero entries

  m = 0;

  for (i=0; i<tmpRow.n; i++)
  {
    c = 0;
    j = 0;

    v = &tmpCol[i].firstItem();
    w = &tmpRow[i].firstItem();

    while (c < tmpCol[i].n)
    {
      k = v->val; v = (VectorCoeff<int>*)(v->next); c++;

      while (j < tmpRow[i].n)

        if (w->val < k) { w = (VectorCoeff<int>*)(w->next); j++; } else break;

      if (j == tmpRow[i].n || w->val != k)
      {
        nvc           = new VectorCoeff<int>;
        nvc->val      = k;
        nvc->next     = (ListItem*)w;
        w->prev->next = (ListItem*)nvc;
        nvc->prev     = w->prev;
        w->prev       = (ListItem*)nvc;
        tmpRow[i].n++;

        tmpRowP[i].insert(perm.n+m++,j++);
      }
    }
  }

  tmpCol.free();

  if (m > 0)
  {
    // adjust memory allocation

    perm.setDim(perm.n + m);

    for (i=0; i<m; i++) x.append(0);

    // regenerate row and col

    row.free();
    col.free();

    k = 0;

    for (i=0; i<tmpRow.n; i++)
    {
      v = &tmpRow [i].firstItem();
      w = &tmpRowP[i].firstItem();

      for (j=0; j<tmpRow[i].n; j++)
      {
        row.append(i + 1);
        col.append(v->val);
        perm[k++] = w->val;
        v = (VectorCoeff<int>*)(v->next);
        w = (VectorCoeff<int>*)(w->next);
      } 
    }
  }

  if (row.n != x.n || col.n != x.n) 
  {
    std::cout << "   MatrixSparse<Type>::makeSymmetric: row.n != x.n || col.n != x.n\n\n";
    exit(1);
  }

  tmpRow.free();
  tmpRowP.free();

  return m;
}

