#define _MATRIX_H_
#ifdef _MATRIX_H_

#include <math.h>
#include <vector>

using namespace std;

const double EPS = 10e-9;

void print(const vector < vector <double> >&);
void copy(vector <double> &, const vector <double> &);
void copy(vector < vector <double> >&, const vector <vector <double>> & );

vector < vector <double> > multiplication(const vector < vector <double> >& l, const vector < vector <double> >&  r)
{
	if (l.size() == r[0].size()){

		vector <vector <double> > res;
		res.resize(r.size());
		for (int i = 0; i < l.size(); ++i){
			res[i].resize(l[0].size());
		}

		for (int i = 0; i < l[0].size(); ++i)
		{
			for (int j = 0; j < l[0].size(); ++j)
			{
				for (int k = 0; k < l[0].size(); ++k)
				{
					res[i][j] += l[i][k] * r[k][j];
				}
				if (fabs(res[i][j]) < EPS) { res[i][j] = 0.0; }
			}
		}

		return res;
	}
}

vector < vector <double> > addition(const vector < vector <double> >& l, const vector < vector <double> >&  r)
{
	vector <vector <double> > res;
	res.resize(r.size());
	for (int i = 0; i < l.size(); ++i){
		res[i].resize(l[0].size());
	}

	for (int i = 0; i < l.size(); ++i){
		for (int j = 0; j < l.size(); ++j){
			res[i][j] = l[i][j] + r[i][j];
		}
	}
	return res;
}

vector < vector <double> > subtraction (vector < vector <double> >& l, vector < vector <double> >&  r)
{
	vector <vector <double> > res;
	res.resize(r.size());
	for (int i = 0; i < l.size(); ++i){
		res[i].resize(l[0].size());
	}

	for (int i = 0; i < l.size(); ++i){
		for (int j = 0; j < l.size(); ++j){
			res[i][j] = l[i][j] - r[i][j];
		}
	}
	return res;
}

int calculateRank ( vector < vector <double> >& x )
{
	int r = 0;

	for (int i = 0; i < x[0].size(); ++i)
	{
		for (int j = 0; j < x.size(); ++j)
		{
			if (x[i][j]) {
				++r;
				break;
			}
		}
	}
	return r;
}

vector < vector <double> > identityMatrix(const int n){
	vector < vector <double> > res;
	res.resize(n);
	for (int i = 0; i < n; ++i){
		res[i].resize(n);
	}
	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j){
			res[i][j] = 0;
		}
		res[i][i] = 1;
	}
	return res;
}

void swapRows(vector < vector <double> >& x, const int r1, const int r2){
	int sz = x.size();

	if (r1 >= sz || r2 >= sz){
		cout << "Error in swapRows: one of numbers bigger than numbers of rows" << endl;
		return;
	}
	double tmp;

	for ( int i = 0; i < x.size(); ++i ){
		tmp = x[i][r1];
		x[i][r1] = x[i][r2];
		x[i][r2] = tmp;
	}
	return;
}

void swapColumns(vector < vector <double> >& x, const int r1, const int r2){
	int sz = x[0].size();

	if (r1 >= sz || r2 >= sz)
	{
		cout << "Error in swapRows: one of numbers bigger than numbers of rows" << endl;
		return;
	}
	double tmp;

	for (int i = 0; i < x.size(); ++i)
	{
		tmp = x[r1][i];
		x[r1][i] = x[r2][i];
		x[r2][i] = tmp;
	}
	return;
}

void PLU(const vector < vector <double> > &A, vector < vector <double> > &C, vector < vector <double> >& P) {
	//n - размерность исходной матрицы
	const int n = A.size();

	C = A;

	//загружаем в матрицу P единичную матрицу
	P = identityMatrix(A.size());

	for (int i = 0; i < n; ++i) {
		//поиск опорного элемента
		double pivotValue = 0;
		int pivot = -1;
		for (int row = i; row < n; row++) {
			if (fabs(C[row][i]) > pivotValue) {
				pivotValue = fabs(C[row][i]);
				pivot = row;
			}
		}
		if (pivotValue == 0) {
			cout << "Error: Matrix is singular" << endl;
			continue;
		}

		//меняем местами i-ю строку и строку с опорным элементом
		swapRows( P, pivot, i);
		swapRows( C, pivot, i);

		for (int j = i + 1; j < n; j++) {
			C[j][i] /= C[i][i];
			for (int k = i + 1; k < n; k++)
				C[j][k] -= C[j][i] * C[i][k];
		}
	}
	return;
}

void fromOnetoTwo( const vector < vector <double> > &A, vector < vector <double> > &L, vector < vector <double> > &U){
	const int sz = A.size();

	vector < vector <double> > E;
	E.resize( sz );
	for (int i = 0; i < sz; ++i){
		E[i].resize( sz );
	}
	E = identityMatrix ( sz );
	
	vector < vector <double> > TMP;
	TMP.resize(sz);
	for (int i = 0; i < sz; ++i){
		TMP[i].resize(sz);
	} 

	TMP = addition ( A, E );
	

	L = identityMatrix(sz);
	
	for (int i = 0; i < sz; i++){
		for (int j = 0; j < i; j++){
			L[i][j] = TMP[i][j];
		}
	}

	U = subtraction( TMP, L);

	cout << "Matrix L:" << endl;
	print(L);
	std::cout << "Matrix U:" << std::endl;
	print(U);
}

int GaussJordanElimination ( vector <vector < double > >& A, vector < double > & b){
	vector <vector < double > > TMP;
	TMP.resize ( A.size() );

	for (int i = 0; i < A.size(); ++i)
		TMP[i].resize(A[0].size() + 1);

	int sz = A[0].size();

	for (int i = 0; i < sz; ++i){
		for (int j = 0; j < sz; ++j){
			TMP[i][j] = A[i][j];
		}
	}

	for (int i = 0; i < sz; i++){
		TMP[i][sz] = b[i];
	}

	std::vector<int> where ( sz, -1);

	const int n = A[0].size();
	const int m = A.size();

	for (int col = 0, row = 0; col < m && row < n; ++col){
		int sel = row;
		for (int i = row; i < n; i++){
			if ( fabs(TMP[i][col]) > abs(TMP[sel][col])){
				sel = i;
			}
		}
		if ( fabs (TMP[sel][col]) < EPS)
			continue;
		for (int i = col; i <= m; i++){
			swapColumns( TMP, sel, i );
		}

		where[col] = row;
		
		for (int i = 0; i < n; i++){
			if (i != row){
				double tmp = TMP[i][col] / TMP[row][col];
				for (int j = col; j <= m; ++j){
					TMP[i][j] -= TMP[row][j] * tmp;
				}
			}
		}
		++row;
	}
	for (int i = 0; i < b.size(); ++i){
		b[i] = 0;
	}

	for (int i = 0; i < m; i++){
		if (where[i] != -1){
			b[i] = TMP[where[i]][m] / TMP[where[i]][i];
		}
	}

	for (int i = 0; i < A[0].size(); i++){
		for (int j = 0; j < A.size(); j++){
			A[i][j] = TMP[i][j];
		}
	}

	for (int i = 0; i < n; ++i) {
		double sum = 0.0;
		for (int j = 0; j < m; ++j)
			sum += b[j] * TMP[i][j];
		if ( abs (sum - TMP[i][m]) > EPS ){
			return 0;
		}
	}

	for (int i = 0; i < m; i++){
		if (where[i] == -1){
			return 2;
		}
	}
	return 1;
}

vector <double> solveLinerSystem (const vector < vector <double> > & P, const vector <vector <double> > & L, const vector <vector <double> > & U){
	vector <double> y;
	y.resize ( P.size() );

	for (int i = 0; i < P.size(); ++i){
		y[i] = P[i][0];
		for (int j = 0; j <= i - 1; ++j){
			y[i] -= L[i][j] * y[j];
		}
	}

	vector <double> x;
	x.resize(P.size());

	for (int i = P.size() - 1; i >= 0; --i){
		x[i] = y[i];
		for (int j = P.size() - 1; j > i; --j){
			x[i] -= U[i][j] * x[j];
		}
		x[i] /= U[i][i];
	}
	return x;
}

vector <vector <double> >& getInverseMatrix(const vector <vector <double> > & P, const vector < vector <double> >& L, const vector <vector <double>> & U)
{
	vector <vector <double> > Inv;
	Inv.resize(L[0].size());
	for (int i = 0; i < L[0].size(); ++i){
		Inv[i].resize(L[0].size());
	}
	
	vector < vector <double> > PL;
	PL.resize(L[0].size());
	for (int i = 0; i < L[0].size(); ++i){
		PL[i].resize(L[0].size());
	}
	copy (PL, multiplication( P, L));

	for ( int i = 0; i < L[0].size(); ++i )
	{
		vector <vector <double>> Ei;
		Ei.resize ( L[0].size());

		for (int j = 0; j < L[0].size(); ++j){
			Ei[j].resize(1);
			Ei[j][0] = 0;
		}

		Ei[i][0] = 1.0;

		vector <double> Ai;
		
		Ai.resize(L[0].size());

		for (int j = 0; j < L[0].size(); ++j){
			Ai[j] = 0;
		}
		
		copy( Ai, solveLinerSystem( Ei, L, U));

		for (int j = 0; j < L.size(); ++j){
			Inv[j][i] = Ai[j];
			if (fabs (Inv[j][i]) < EPS) Inv[j][i] = 0.0;
		}
	}

	std::cout << std::endl << "d) Inverse matrix:" << std::endl;
	print(multiplication(P,Inv));
	return Inv;
}

void print(const vector < vector <double> >& x){
	for (int i = 0; i < x.size(); ++i){
		for (int j = 0; j < x[0].size(); ++j){
			cout << x[i][j] << " ";
		}
		cout << endl;
	}
	return;
}

void copy(vector < vector <double> >& l, const vector <vector <double>> & r){
	for (int i = 0; i < l[0].size(); ++i){
		for (int j = 0; j < l.size(); ++j){
			l[i][j] = r[i][j];
		}
	}
}
void copy( vector <double> & l, const vector <double> & r){
	for (int j = 0; j < l.size(); ++j){
		l[j] = r[j];
	}
}

#endif	