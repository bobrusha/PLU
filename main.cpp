#include "examples.h"
using namespace std;

int main(){

	initializeExample(2);
	Matrix A(3, 3, 'A');
	A = M;
	Matrix B(3, 3, 'B');
	Matrix P(3, 3, 'P');
	Matrix L(3, 3, 'L');
	Matrix U(3, 3, 'U');
	A.print();

	PLU(A, B, P, L, U);
	
	if (!A._issingular){
		L.det();
		U.det();

		cout << "a) det = " << (P._kolvop / abs(P._kolvop))* L._det * U._det << endl;

		cout << "b) rank = " << L.rank() << endl;										// т.к. rank(A*B) = rank(A);

		//Решение СЛАУ
		cout << "c) solution of liner system." << endl << "   vector-column b:" << endl;
		b.print();

		Matrix x(3, 1, 'x');
		Matrix tmp(3, 3, 't');
		tmp = P;
		tmp*b;
		solveLinerSystem(tmp, L, U, x);
		cout << "Result:" << endl;
		x.print();
		//Вычисление обратной матрицы
		Matrix Inv(3, 3, 'I');
		getInverseMatrix(L, U, Inv, A._issingular);
		Inv.swapColumn(0, 1);
		Inv.swapColumn(1, 2);
		Inv.print();

		cout << "   Test:" << endl;
		Matrix Test(3, 3, 'T');
		A.print();
		Test = A;
		Test*Inv;
		Test.print();


		//Число обусловленности матрицы
		A.getConditionNumber_1(Inv);
	}
	else{ 
		cout << " a) det = " << 0 << endl;
		
		Matrix tmp_b(3, 1, 'b');
		tmp_b = b;
		A.gaussian_elimination(b);
	}
	

	system("Pause");
	return 0;
}