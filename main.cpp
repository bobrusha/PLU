#include "matrix.h"
using namespace std;

int main(){
	/*
	double Arr1[][3] = { { 2, 7, -6 },
	{ 8, 2, 1 },
	{ 7, 4, 2 } }; // det = -171
	*/

	double Arr1[][3] = {{ 2, 1, 1 },
						{ 1, -1, 0},
						{ 3, -1, 2}};
	Matrix A(3, 3, &Arr1[0][0], 'A');
	Matrix B(3, 3, 'B');
	A.print();
	Matrix P(3, 3, 'P');
	Matrix L(3, 3, 'L');
	Matrix U(3, 3, 'U');

	PLU(A, B, P, L, U);

	//нахождение обратной матрицы
	/*Matrix E(3, 3, 'E');
	E.IdentityMatrix();

	Matrix U1(3, 3, 'U');
	U1.IdentityMatrix();

	//U*U1 = E;
	
	Matrix M(3, 4, 'M');
	for (int i = 0; i < U._column; i++){
		for (int j = 0; j < U._column; j++){
			double tmp = 0.0;
			for (int k = 0; k < U._column; k++){
				tmp +=  ;
			}
		}
	}
	*/

	L.det();
	U.det();
	
	cout << "a) det = " << P._sign * L._det * U._det<< endl;		// т.к. rank(A*B) = rank(A);
	cout << "b) rank = " << L.rank()<< endl;
	
	L.print();
	U.print();
	P.print();
	double arr2[3][1] = { { 2.0 }, { -2.0 }, { 2.0 } };
	
	Matrix b(3, 1, (& arr2)[0][0],'b');
	Matrix x(3, 1, 'x');
	Matrix y(3, 1, 'y');

	b.print();

	P = P*b;
	
	cout << "Matrix P:" << endl;
	P.print();

	for (int i = 0; i<b._row; i++){
		double tmp = 0.0;
		for (int j = 0; j < i; j++){
			tmp += y.M[i][0] * L.M[j][i];
		}
		y.M[i][0] = (b.M[i][0] - tmp)/L.M[i][i];
	}
	
	cout << "Matrix y:" << endl;
	y.print();

	for (int i = b._row - 1; i > 0; i--){
		double tmp = 0.0;
		for (int j = b._row-1; j > i; j--){
			tmp += x.M[i][0] * U.M[j][i];
		}
		x.M[i][0] = (y.M[i][0] - tmp) / U.M[i][i];
	}
	x.print();
	/*
	y.M[0][0] = P.M[0][0]/L.M[0][0];
	for (int i = 1; i < P._column; i++){
		double tmp = 0.0;
		for (int j = 0; j < i; j++){
			tmp += L.M[i][j] * y.M[0][j];
		}
		tmp -= P.M[0][i];
		tmp /= L.M[i][i];
		y.M[i][0] = tmp;
		y.print();
	}
	y.print();

	x.M[x._column][0] = y.M[0][x._column] / U.M[0][x._column];
	for (int i = x._column; i > 0; i--)
	{
		double tmp = 0.0;
		for (int j = U._column; j > i; j--){
			tmp += U.M[i][j] * x.M[j][0];
		}
		tmp -= y.M[i][0];
		tmp /= U.M[i][i];
		x.M[i][0] = tmp;
	}
	
	x.print();
	*/
	system("Pause");
	return 0;

}