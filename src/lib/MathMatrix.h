
#ifndef incl_MathMatrix_h
#define incl_MathMatrix_h

#include <iostream>

#include "MathBasic.h"
#include "MathVector.h"
#include "MyString.h"
#include "FunctionsProgram.h"



using namespace std;


extern int nVector;


namespace MatricesWulf
{


	//  matrix classes
	//
	//  MatrixFullArray       - full matrix stored as one long array, containing rows or columns
	//  MatrixSparse          - sparse matrix stored with linked list vectors
	//  MatrixSparseArray     - sparse matrix stored with arrays

	//
	//  TAKE CARE WHEN YOU DEAL WITH MATRICES WITH "DOUBLE ENTRIES"
	//
	//

	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//
	// abstract matrix base class

	template<typename Type> class MatrixBase
	{
	public:

		MatrixBase(void);

		virtual ~MatrixBase();

		int nRow, nCol;

		virtual Type &operator()(int, int, bool doubleEntries = false) = 0;

		virtual bool getEntry(int, Type &, int &, int &) = 0;

		virtual int  nNonZero(void) = 0;

		virtual void free(void) = 0;

		virtual bool nonZero(int, int) { return true; }

		virtual void print(std::ostream &);

		virtual void print(int dd = 12, int df = 5);

		virtual void showInteractive(int dd = 12, int df = 5, int dr = 10, int dc = 10);

		virtual void fmtChar(char *, int dd = 12, int df = 5);

		virtual void copyToSimpleMatrix(Type*&);

		virtual void printPattern(void);
	};


	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//
	// full matrix stored as one long array, containing rows or columns

	template<typename Type> class MatrixFullArray : public MatrixBase<Type>, public ListItem
	{
	public:

		MatrixFullArray(void);

		virtual ~MatrixFullArray();

		Type *x;

		bool transpose;

		void setDim(int, int, bool trnsp = false);

		virtual Type &operator()(int, int, bool doubleEntries = false);

		virtual MatrixFullArray<Type> &operator=(MatrixBase<Type> &);

		virtual MatrixFullArray<Type> &operator=(MatrixFullArray<Type> &);

		virtual void free(void);

		virtual bool getEntry(int, Type &, int &, int &);

		virtual int  nNonZero(void) { return this->nRow * this->nCol; }

		virtual void zero(void);

		virtual void takeOver(MatrixFullArray<Type> &mtx2);
	};


	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//
	// sparse matrix stored with linked list vectors

	template<typename Type> class MatrixSparse : public MatrixBase<Type>, public ListItem
	{
	public:

		MatrixSparse(void);

		virtual ~MatrixSparse();

		myVector<Type> x;

		myVector<int> row, col;

		virtual Type &operator()(int, int, bool doubleEntries = false);

		virtual Type &operator[](int);

		virtual MatrixSparse<Type> &operator=(MatrixBase<Type> &);

		virtual bool nonZero(int, int);

		virtual bool getEntry(int, Type &, int &, int &);

		virtual int  nNonZero(void) { return row.n; }

		virtual void free(void);

		virtual void printAll(void);

		virtual void zero(void);

		int append(int, int, Type val = (Type)0);

		void quickSort(bool, bool, VectorArray<int> &);

		void copySort(bool, bool, bool, VectorArray<int> &);

		int  makeSymmetric(bool, VectorArray<int> &);

		void append(int, int, MatrixBase<Type> &, bool trans = false);

		void plotPattern(char*);
	};


	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//
	// sparse matrix stored with arrays

	template<typename Type> class MatrixSparseArray : public MatrixBase<Type>, public ListItem
	{
	public:

		MatrixSparseArray(void);

		virtual ~MatrixSparseArray();

		VectorArray<Type> x;

		VectorArray<int> row, col;

		virtual Type &operator()(int, int, bool doubleEntries = false);

		virtual Type &operator[](int);

		virtual MatrixSparseArray<Type> &operator=(MatrixBase<Type> &);

		virtual bool nonZero(int, int);

		virtual bool getEntry(int, Type &, int &, int &);

		virtual int  nNonZero(void) { return row.n; }

		virtual void free(void);

		virtual void printAll(void);

		virtual void zero(void);

		void plotPattern(char*);
	};


#include "MathMatrix2.h"


	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	template<typename Type> MatrixBase<Type>::MatrixBase(void)
	{
		nRow = 0;
		nCol = 0;

		return;
	}




	template<typename Type> MatrixBase<Type>::~MatrixBase()
	{
		return;
	}




	template<typename Type> void MatrixBase<Type>::print(std::ostream &os)
	{
		int i, j;

		Type dmy;

		for (i = 0; i < nRow; i++)
		{
			for (j = 0; j < nCol; j++)
			{
				if (nonZero(i + 1, j + 1)) dmy = (*this)(i + 1, j + 1); else dmy = (Type)0;

				os << " " << dmy;
			}
			os << "\n";
		}
		os << "\n";

		return;
	}




	template<typename Type> void MatrixBase<Type>::print(int dd, int df)
	{
		int i, j;

		Type dmy;

		char fmt[30];

		fmtChar(fmt, dd, df);

		// print submatrix

		for (i = 1; i < nRow + 1; i++)
		{
			for (j = 1; j < nCol + 1; j++)
			{
				if (nonZero(i, j)) dmy = (*this)(i, j); else dmy = (Type)0;

				printf(fmt, dmy);
			}
			cout << "\n";
		}

		return;
	}




	template<typename Type> void MatrixBase<Type>::showInteractive(int dd, int df, int dr, int dc)
	{
		int i, j, ar = 1, ac = 1, nw = 0, lr, lc;

		Type dmy;

		MyString ch, *word = NULL;

		char *quit[] = QUIT, fmt[30];

		fmtChar(fmt, dd, df);

		while (1)
		{
			// print submatrix

			lr = ar + dr;  if (lr > nRow) lr = nRow + 1;
			lc = ac + dc;  if (lc > nCol) lc = nCol + 1;

			for (i = ar; i < lr; i++)
			{
				for (j = ac; j < lc; j++)
				{
					if (nonZero(i, j)) dmy = (*this)(i, j); else dmy = (Type)0;

					printf(fmt, dmy);
				}
				cout << "\n";
			}

			// select rows and columns

			while (1)
			{
				if (word != NULL) { for (i = 0; i < nw; i++) word[i].free(); delete[] word; word = NULL; }

				cout << "\n    input row and columns or 'q' to quit! ";

				ch.input();

				if (ch.which(quit) != -1) break;

				nw = ch.stripToMin().split(&word);

				if (nw > 1)  if (word[0].toInt(&ar))  if (word[1].toInt(&ac))  break;
			}
			cout << "\n";

			if (ch.which(quit) != -1) break;

			if (ar < 1) ar = 1; if (ar > nRow) ar = nRow;
			if (ac < 1) ac = 1; if (ac > nCol) ac = nCol;
		}

		if (word != NULL) { for (i = 0; i < nw; i++) word[i].free(); delete[] word; }

		return;
	}




	template<typename Type> void MatrixBase<Type>::fmtChar(char *fmt, int dd, int df)
	{
		sprintf(fmt, "  %d.%df", dd, df); fmt[1] = '%';

		return;
	}




	template<> inline void MatrixBase<int>::fmtChar(char *fmt, int dd, int df)
	{
		sprintf(fmt, "  %dd", dd); fmt[1] = '%';

		return;
	}




	template<typename Type> void MatrixBase<Type>::copyToSimpleMatrix(Type* &mtx)
	{
		int i, r, c;

		Type val;

		mtx = new Type[nRow*nCol];

		for (i = 0; i < nRow*nCol; i++) mtx[i] = (Type)0;

		i = 0; while (getEntry(i++, val, r, c)) mtx[c*nRow + r] += val;

		return;
	}




	template<typename Type> void MatrixBase<Type>::printPattern(void)
	{
		int i, j;

		cout << "\n   ";
		for (j = 1; j < nCol + 1; j++)  printf(" %1i", j % 10);

		for (i = 1; i < nRow + 1; i++)
		{
			printf("\n%3i ", i);
			for (j = 1; j < nCol + 1; j++)
			{
				if (nonZero(i, j)) cout << "x "; else cout << ". ";
			}
		}
		cout << "\n\n";

		return;
	}




	template<typename Type> std::ostream &operator<<(std::ostream &os, MatrixBase<Type> &mtx)
	{
		//if (os == cout) mtx.print();

		//else            mtx.print(os);

		return os;
	}



	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	template<typename Type> MatrixFullArray<Type>::MatrixFullArray(void)
	{
		x = NULL;

		return;
	}




	template<typename Type> MatrixFullArray<Type>::~MatrixFullArray()
	{
		free();

		return;
	}




	template<typename Type> void MatrixFullArray<Type>::setDim(int n, int m, bool trnsp)
	{
		free();

		this->nRow = n;
		this->nCol = m;

		transpose = trnsp;

		x = new Type[n*m];

		return;
	}




	template<typename Type> Type &MatrixFullArray<Type>::operator()(int i, int j, bool doubleEntries)
	{
		if (transpose) return x[(i - 1)*this->nCol + j - 1];

		else           return x[(j - 1)*this->nRow + i - 1];
	}




	template<typename Type>
	MatrixFullArray<Type> &MatrixFullArray<Type>::operator=(MatrixBase<Type> &b)
	{
		free();

		setDim(b.nRow, b.nCol);

		zero();

		int i = 0, r, c;

		Type val;

		while (b.getEntry(i++, val, r, c)) x[c*this->nRow + r] = val;

		return *this;
	}




	template<typename Type>
	MatrixFullArray<Type> &MatrixFullArray<Type>::operator=(MatrixFullArray<Type> &b)
	{
		free();

		setDim(b.nRow, b.nCol, b.transpose);

		int i, n = b.nRow * b.nCol;

		for (int i = 0; i < n; i++) x[i] = b.x[i];

		return *this;
	}




	template<typename Type> void MatrixFullArray<Type>::free(void)
	{
		if (x != NULL) delete[] x; x = NULL;

		this->nRow = 0;
		this->nCol = 0;

		return;
	}




	template<typename Type> bool MatrixFullArray<Type>::getEntry(int i, Type &val, int &r, int &c)
	{
		//           c
		//         0 1 2
		//
		//     0   0 3 6
		//   r 1   1 4 7  <- i
		//     2   2 5 8

		if (i + 1 > this->nRow*this->nCol) return false;

		val = x[i];

		r = i + 1 % this->nRow;

		c = (i + 1 - r) / this->nRow;

		r--;

		return true;
	}




	template<typename Type> void MatrixFullArray<Type>::zero(void)
	{
		int i, m = this->nRow * this->nCol;

		for (i = 0; i < m; i++)  x[i] = (Type)0;

		return;
	}




	template<typename Type> void MatrixFullArray<Type>::takeOver(MatrixFullArray<Type> &mtx2)
	{
		if (x != NULL) delete[] x;

		x = mtx2.x;           mtx2.x = NULL;
		this->nRow = mtx2.nRow;        mtx2.nRow = 0;
		this->nCol = mtx2.nCol;        mtx2.nCol = 0;
		transpose = mtx2.transpose;   mtx2.transpose = false;

		return;
	}




	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	template<typename Type> MatrixSparse<Type>::MatrixSparse(void)
	{
		return;
	}




	template<typename Type> MatrixSparse<Type>::~MatrixSparse()
	{
		free();

		return;
	}




	template<typename Type> Type &MatrixSparse<Type>::operator()(int ii, int jj, bool doubleEntries)
	{
		int i = 0, j, n = row.n;

		while (i < n) if (!(row[i] == ii && col[i] == jj)) i++; else break;

		if (i < n)
		{
			if (doubleEntries)

				for (j = i + 1; j<n; j++) if (row[j] == ii && col[j] == jj) { x[i] += x[j]; x[j] = (Type)0; }

			return x[i];
		}

		prgError(1, "MatrixSparse<Type>::operator()", "invalid indices!");

		return x[0];
	}




	template<typename Type> Type &MatrixSparse<Type>::operator[](int ii)
	{
		if (ii > x.n - 1 || ii < 0) prgError(1, "MatrixSparse<Type>::operator[]", "index out of bounds!");

		return x[ii];
	}




	template<typename Type> MatrixSparse<Type> &MatrixSparse<Type>::operator=(MatrixBase<Type> &b)
	{
		free();

		this->nRow = b.nRow;
		this->nCol = b.nCol;

		int i = 0, r, c;

		Type val;

		while (b.getEntry(i, val, r, c)) { append(r + 1, c + 1); x[i++] = val; }

		return *this;
	}




	template<typename Type> bool MatrixSparse<Type>::nonZero(int ii, int jj)
	{
		int i = 0, n = row.n;

		while (i < n && !(row[i] == ii && col[i] == jj)) i++;

		if (i < n) return true;

		return false;
	}




	template<typename Type> bool MatrixSparse<Type>::getEntry(int i, Type &val, int &r, int &c)
	{
		if (i + 1 > x.n) return false;

		val = x[i];
		r = row[i] - 1;
		c = col[i] - 1;

		return true;
	}




	template<typename Type> void MatrixSparse<Type>::free(void)
	{
		x.free();
		row.free();
		col.free();

		this->nRow = 0;
		this->nCol = 0;

		return;
	}




	template<typename Type> void MatrixSparse<Type>::printAll(void)
	{
		std::cout << " nRow = " << this->nRow << "\n nCol = " << this->nCol
			<< "\n rows = " << row << "\n cols = " << col << "\n coef = " << x << "\n\n";

		return;
	}




	template<typename Type> void MatrixSparse<Type>::zero(void)
	{
		for (int i = 0; i < x.n; i++) x[i] = (Type)0;

		return;
	}




	template<typename Type> int MatrixSparse<Type>::append(int ii, int jj, Type val)
	{
		int n = row.n;

		row.append(ii);
		col.append(jj);

		x.append(val);

		this->nRow = max(ii, this->nRow);
		this->nCol = max(jj, this->nCol);

		return n;
	}



	/*
	template<typename Type> int MatrixSparse<Type>::makeSymmetric(bool rowCompFlg,
	VectorArray<int> &perm)
	{
	if (!rowCompFlg)
	{
	std::cout << "   MatrixSparse<Type>::makeSymmetric: cannot yet handle rowCompFlg == false!\n\n";
	exit(1);
	}

	int c, i, j, k, n, m = row.n;

	ListArray< myVector<int> > tmpCol, tmpRow, tmpRowP;

	tmpCol.setDim(this->nCol);
	tmpRow.setDim(this->nRow);
	tmpRowP.setDim(this->nRow);

	i = 0;

	// generate tmpCol and tmpRow

	for (k=0; k<m; k++)
	{
	if (i > row[k])
	{
	std::cout << "   MatrixSparse<Type>::makeSymmetric: matrix is not row compressed!\n\n";
	exit(1);
	}
	else if (row[k] - i > 1)
	{
	std::cout << "   MatrixSparse<Type>::makeSymmetric: empty row?!\n\n";
	exit(1);
	}
	i = row[k];
	tmpCol[col[k]-1].append(i);
	tmpRow[i-1].append(col[k]);
	tmpRowP[i-1].append(perm.x[k]);
	}

	//

	n = 0;

	for (i=0; i<tmpRow.n; i++)
	{
	c = 0;

	j = 0;

	while (c < tmpCol[i].n)
	{
	k = tmpCol[i][c++];

	while (j < tmpRow[i].n)

	if (tmpRow[i][j] < k) j++; else break;

	if (j == tmpRow[i].n || tmpRow[i][j] != k)
	{
	tmpRow [i].insert( k,j);
	tmpRowP[i].insert(perm.n+n++,j++);
	}
	}
	}

	//

	if (n == 0) return 0;

	row.free();
	col.free();
	perm.setDim(perm.n + n);
	for (i=0; i<n; i++) x.append(0);

	k = 0;

	for (i=0; i<tmpRow.n; i++)
	{
	for (j=0; j<tmpRow[i].n; j++)
	{
	row.append(i + 1);
	col.append(tmpRow[i][j]);
	perm[k++] = tmpRowP[i][j];
	}
	}

	if (row.n != x.n || col.n != x.n)
	{
	std::cout << "   MatrixSparse<Type>::makeSymmetric: row.n != x.n || col.n != x.n\n\n";
	exit(1);
	}

	return n;
	}
	*/



	template<typename Type> void MatrixSparse<Type>::append(int iOff, int jOff,
		MatrixBase<Type> &mtx, bool trans)
	{
		Type val;

		this->nRow = max(this->nRow, iOff + mtx.nRow - 1);
		this->nCol = max(this->nCol, jOff + mtx.nCol - 1);

		int i = 0, c, r;

		if (!trans) while (mtx.getEntry(i++, val, r, c)) append(iOff + r, jOff + c);

		else        while (mtx.getEntry(i++, val, r, c)) append(iOff + c, jOff + r);

		return;
	}




	template<typename Type> void MatrixSparse<Type>::plotPattern(char *fileName)
	{
		prgPlotMatrixPattern(this->nRow, this->nCol, row, col, fileName);

		return;
	}




	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


	template<typename Type> MatrixSparseArray<Type>::MatrixSparseArray(void)
	{
		return;
	}




	template<typename Type> MatrixSparseArray<Type>::~MatrixSparseArray()
	{
		free();

		return;
	}




	template<typename Type> Type &MatrixSparseArray<Type>::operator()(int ii, int jj,
		bool doubleEntries)
	{
		int i = 0, j, n = row.n;

		while (i < n && !(row[i] == ii && col[i] == jj)) i++;

		if (i < n)
		{
			if (doubleEntries)
			{
				for (j = i + 1; j < n; j++)
				{
					if (row[j] == ii && col[j] == jj)
					{
						x[i] += x[j]; x[j] = (Type)0;
					}
				}
			}

			return x[i];
		}
		else
			prgError(1, "MatrixSparseArray<Type>::operator()", "invalid indices!");

		return x[0];
	}




	template<typename Type> Type &MatrixSparseArray<Type>::operator[](int ii)
	{
		if (ii > x.n - 1 || ii < 0)

			prgError(1, "MatrixSparseArray<Type>::operator[]", "index out of bounds!");

		return x[ii];
	}




	template<typename Type>
	MatrixSparseArray<Type> &MatrixSparseArray<Type>::operator=(MatrixBase<Type> &b)
	{
		free();

		this->nRow = b.nRow;
		this->nCol = b.nCol;

		int i = 0, m = 0, r, c;

		Type val;

		while (b.getEntry(m, val, r, c)) m++;

		x.setDim(m);
		row.setDim(m);
		col.setDim(m);

		while (b.getEntry(i, val, r, c)) { row[i] = r + 1; col[i] = c + 1; x[i++] = val; }

		return *this;
	}




	template<typename Type> bool MatrixSparseArray<Type>::nonZero(int ii, int jj)
	{
		int i = 0, n = row.n;

		while (i < n && !(row[i] == ii && col[i] == jj)) i++;

		if (i < n) return true;

		return false;
	}



	template<typename Type> bool MatrixSparseArray<Type>::getEntry(int i, Type &val, int &r, int &c)
	{
		if (i + 1 > x.n) return false;

		val = x[i];
		r = row[i] - 1;
		c = col[i] - 1;

		return true;
	}




	template<typename Type> void MatrixSparseArray<Type>::free(void)
	{
		x.free();
		row.free();
		col.free();

		this->nRow = 0;
		this->nCol = 0;

		return;
	}




	template<typename Type> void MatrixSparseArray<Type>::printAll(void)
	{
		std::cout << " nRow = " << this->nRow << "\n nCol = " << this->nCol
			<< "\n rows = " << row << "\n cols = " << col << "\n coef = " << x << "\n\n";

		return;
	}




	template<typename Type> void MatrixSparseArray<Type>::zero(void)
	{
		for (int i = 0; i < x.n; i++) x[i] = (Type)0;

		return;
	}




	template<typename Type> void MatrixSparseArray<Type>::plotPattern(char *fileName)
	{
		prgPlotMatrixPattern(this->nRow, this->nCol, row, col, fileName);

		return;
	}




	//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

}

#endif

