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
						{ 3, -1, 2}}; // det = -4
	Matrix A(3, 3, &Arr1[0][0], 'A');
	Matrix B(3, 3, 'B');
	A.print();
	Matrix P(3, 3, 'P');
	Matrix L(3, 3, 'L');
	Matrix U(3, 3, 'U');

	PLU(A, B, P, L, U);

	L.det();
	U.det();
	
	std::cout << "a) det = " << (P._kolvop / abs(P._kolvop))* L._det * U._det<< std::endl;		// ò.ê. rank(A*B) = rank(A);
	std::cout << "b) rank = " << L.rank()<< std::endl;
	
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

	//P = P*b;
	//cout << "Matrix P:" << endl;
	//P.print();

	//SLAU(P, L, U, x);
	Matrix Inv(3, 3, 'I');
	getInverseMatrix(L, U, Inv);
	P*Inv;
	Inv.print();
	/*
	Matrix a(3,1, 'a');
	Matrix A1( 3, 3,'A');
	for(int i = 0; i< A._row; i++){
		Matrix Ai( 3, 1, 'A');
		Matrix b( 3, 1, 'b');
		b.M[i][0]= 1;
		P = P*b;
		SLAU(P, L, U, a, Ai);
		for (int j = 0; j < A1._row; j++){
			A1.M[j][i] = Ai.M[j][0];
		}
	}
	*/
	system("Pause");
	return 0;

}