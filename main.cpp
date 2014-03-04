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

	L.det();
	U.det();
	
	cout << "a) det = " << P._sign * L._det * U._det<< endl;		// ò.ê. rank(A*B) = rank(A);
	cout << "b) rank = " << L.rank()<< endl;
	
	L.print();
	U.print();
	P.print();

	double arr2[3][1] = { { 2.0 }, { -2.0 }, { 2.0 } };
	
	Matrix b(3, 1, (& arr2)[0][0],'b');
	Matrix x(3, 1, 'x');
	Matrix y(3, 1, 'y');

	b.print();
	//L = L*U;
	//A = P*A;
	//A.print();
	//L.print();

	P = P*b;
	
	cout << "Matrix P:" << endl;
	P.print();

	for (int i = 0; i < b._row; i++){
		double tmp = 0.0;
		for (int j = 1; j <= i; j++){
			cout << y.M[j-1][0] << " " << L.M[i][i] << " ";
			tmp += y.M[j-1][0] * L.M[j-1][i];
		}
		y.M[i][0] = (P.M[i][0] - tmp)/L.M[i][i];
		cout << y.M[i][0] << " ";
	}
	
	cout << "Matrix y:" << endl;
	y.print();

	for (int i = b._row - 1; i >= 0; i--){
		double tmp = 0.0;
		for (int j = b._column - 1; j >= i; j--){
			tmp += x.M[i][0] * U.M[i][j];
		}
		x.M[i][0] = (y.M[i][0] - tmp) / U.M[i][i];
	}
	cout << "Matrix x:" << endl;
	x.print();
	

	system("Pause");
	return 0;

}