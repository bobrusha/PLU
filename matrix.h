#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <iostream>

#define _MATRIX_
#ifdef _MATRIX_

class Matrix
{
public:
	int _row, _column;
	double _det;
	double **M = new double*[_row];
	char _name;
	bool _istriangle;
	double _sign;
	Matrix(const int n, const int m, char x)
	{
		_row = n;											// !!
		_column = m;
		_name = x;
		_det = 1.0;
		_istriangle = false;
		_sign = 1.0;
		for (int i = 0; i < _row; i++){
			M[i] = new double[_column];
		}

		for (int i = 0; i < _row; i++)
		{
			for (int j = 0; j < _column; j++)
			{
				M[i][j] = 0.0;
			}
		}
	}
	Matrix(const int n, const int m, double* A, char x)
	{
		_name = x;
		_row = n; 
		_column = m;

		for (int i = 0; i < _row; i++){
			M[i] = new double[_column];
		}

		for (int i = 0; i < n; i++){
			for (int j = 0; j < m; j++){
				M[i][j] = *(A + i*m + j);
			}
		}
	}
	Matrix ( const Matrix& X)									
	{
		if (_row == X._row && _column == X._column)
		{
			for (int i = 0; i < _row; i++)
			{
				for (int j = 0; j < _column; j++)
				{
					M[i][j] = X.M[i][j];
				}
			}
		}
	}

	~Matrix(){
		for (int count = 0; count < _row; count++)
			delete[] M[count];
		
		//delete[]M;  ??
	}
	//=============================================================================== перегрузка операторов

	Matrix& operator = (const Matrix& X) {
		//проверка на самоприсваивание
		if (this == &X) {
			return *this;
		}

		_row = X._row;
		_column = X._column;
		for (int count = 0; count < _row; count++)
			delete[] M[count];
		
		for (int i = 0; i < _row; i++){
			M[i] = new double[_column];
		}
		for (int i = 0; i < _row; i++){
			for (int j = 0; j < _column; j++){
				M[i][j] = X.M[i][j];
			}
		}
		return *this;
	}
	Matrix& operator * (const double x){									//ok
		for (int i = 0; i < _row; i++)
		{
			for (int j = 0; j < _column; j++)
			{
				M[i][j] *= x;
			}
		}
		return *this;
	}
	Matrix& operator * (const Matrix& X)
	{
		if (_column == X._row){
			Matrix T(_row, X._column, 'T');

			for (int i = 0; i < _row; i++)
			{
				for (int j = 0; j < X._column; j++)
				{
					for (int r = 0; r < _column; r++)
					{
						T.M[i][j] += M[i][r] * X.M[r][j];
					}
				}
			}
			*this = T;
		}
		return *this;
	}
	//============================================================
	void print(){
		for (int i = 0; i < _row; i++){
			for (int j = 0; j < _column; j++){
				std::cout<<M[i][j]<<" ";
			}
			std::cout << std::endl;
		}
		std::cout<<std::endl;
	}

	void full(const int k)
	{
		for (int i = 0; i < _column; i++)
		{
			for (int j = 0; j < _row; j++)
			{
				M[i][j] = k;
			}
		}
	}

	void IdentityMatrix()
	{
		if (_row != _column) return;
		this->full(0);
		for (int i = 0; i < _row; i++)
		{
			M[i][i] = 1.0;
		}
		return;
	}

	friend Matrix& operator + (Matrix&, const Matrix&);
	friend Matrix& operator - (Matrix&, const Matrix&);
	
	
	bool gaussian_elimination()							//возвращает true - если чётное количество перестановок false - если нечётное
	{
		//this->print();
		if (!_istriangle){
			bool sign = true;

			for (int i = 0; i < _row - 1; i++){				// поиск опорного элемента по всей матрице
				double tmp = 0.0;
				int num1 = 0, num2 = 0;

				for (int j = i; j < _column; j++)
				{
					for (int k = i; k < _row; k++){
						if (fabs(M[j][k]) > tmp)
						{
							tmp = fabs(M[j][k]);
							num1 = j;
							num2 = k;
						}
					}
				}
				if (i != num1){
					sign = !sign;
				}
				if (i != num2){
					sign = !sign;
				}
				swapRows(i, num1);
				swapColumn(i, num2);

				tmp = M[i + 1][i];
				for (int j = i + 1; j < _row; j++)
				{
					for (int k = i - 1; k < _column; k++){
						M[j][k] -= (tmp * M[i][k]) / M[i][i];
					}
				}
			}
			_istriangle = true;
			return sign;
		}
		return true;
	}
	
	void det(){
		if (_row != _column) { 
			std::cout << "Matrix isn't square!"; 
			return; 
		}
		if (_row == 1) {
			_det = M[0][0];
			return;
		}
		else {
			_det = 1.0;
			for (int i = 0; i < _row; i++){
				_det *= M[i][i];
			}
			return;
		}
	}

	int rank()
	{
		int r = 0;
		
		for (int i = 0; i < _column; i++)
		{
			for (int j = 0; j < _row; j++)
			{
				if (M[i][j]) { 
					++r; 
					break; 
				}
			}
		}
		return r;
	}
	void getInverseMatrix(Matrix& X){ 
		
	}
	friend void PLU(Matrix&, Matrix&, Matrix&, Matrix&, Matrix&);

	void swapRows(int x, int y){
		for (int i = 0; i < _column; i++){
			double tmp = M[x][i];
			M[x][i] = M[y][i];
			M[y][i] = tmp;
			//this->print();
		}
		return;
	}
	void swapColumn(int x, int y)
	{
		double tmp;
		for (int i = 0; i < _row; i++){
			tmp = M[i][x];
			M[i][x] = M[i][y];
			M[i][y] = tmp;
			//this->print();
		}
		return;
	}

};

Matrix& operator + (Matrix& M1, const Matrix& M2){						// make all const
	for (int i = 0; i < M1._row; i++)
	{
		for (int j = 0; j < M1._column; j++)
		{
				M1.M[i][j] += M2.M[i][j];
		}
	}
	return M1;
}
Matrix& operator - (Matrix& M1, const Matrix& M2){						// make all const
	for (int i = 0; i < M1._row; i++)
	{
		for (int j = 0; j < M1._column; j++)
		{
			M1.M[i][j] -= M2.M[i][j];
		}
	}
	return M1;
}

void PLU (Matrix& A, Matrix&B, Matrix& P, Matrix& L, Matrix& U)			//работает правильно вроде
{
	P.IdentityMatrix();
	B = A;
	for (int i = 0; i < A._row; i++){
		double tmp = 0.0;
		int num1 = 0, num2 = 0;

		for (int j = i; j < A._column; j++)
		{
			for (int k = i; k < A._row; k++){
				if (fabs(A.M[j][k]) > tmp)
				{
					tmp = fabs(A.M[j][k]);
					num1 = j;
					num2 = k;
				}
			}
		}
		B.swapRows( i, num1);
		B.swapColumn( i, num2);

		P.swapRows(i, num1);
		P.swapColumn(i, num2);
		
		if (i != num1){
			P._sign = - P._sign;
		}
		if (i != num2){
			P._sign = -P._sign;
		}
		if (num1 == 0 && num2 == 0){
			std::cout << "Matrix is singular"<<std::endl;
			return;
		}
		for (int j = i + 1; j < A._row; j++) {
			B.M[j][i] /= B.M[i][i];
			for (int k = i + 1; k < A._row; k++)
				B.M[j][k] -= B.M[j][i] * B.M[i][k];
		}
	}
	// B = L + U - E;
	// L - нижнетреугольная матрица с единицами на диагонали
	// U - верхнетреугольная матрица
	Matrix E(A._row, A._column, 'E');
	E.IdentityMatrix();
	//B.print();
	B = B + E;
	L.full(0); U.full(0);

	for (int i = 0; i < A._row; i++){ L.M[i][i] = 1; }
	for (int i = 0; i < A._row; i++){
		for (int j = 0; j < i; j++){
			L.M[i][j] = B.M[i][j];
		}
	}
	U = B - L;
	return;
}

//void getLandUMatrix()
#endif