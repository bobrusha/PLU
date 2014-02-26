#include "matrix.h"
using namespace std;

int main(){
	double Arr1[][3] = { { 2, 7, -6 },
	{ 8, 2, 1 },
	{ 7, 4, 2 } }; // det = -171

	Matrix A(3, 3, Arr1, 'A');
	Matrix B(3, 3, 'B');

	Matrix P(3, 3, 'P');
	Matrix L(3, 3, 'L');
	Matrix U(3, 3, 'U');

	PLU(A, B, P, L, U);
	P.print();
	L.det();
	U.det();

	cout << P._sign * L._det * U._det<<endl;

	L = L*U;
	A = P*A;
	P.print();
	A.print();
	
	system("Pause");
	return 0;

}