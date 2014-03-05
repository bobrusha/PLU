#include "examples.h"
using namespace std;

int main(){

	initializeExample(0);
	Matrix A(3, 3, 'A');
	A = M;
	Matrix B(3, 3, 'B');
	Matrix P(3, 3, 'P');
	Matrix L(3, 3, 'L');
	Matrix U(3, 3, 'U');

	PLU(A, B, P, L, U);

	cout << "Test PLU" << endl;
	Matrix Tmp_l(3, 3, 't');
	Matrix Tmp_u(3, 3, 't');
	Matrix Tmp_p(3, 3, 't');
	Matrix Tmp_a(3, 3, 't');
	Tmp_a = A;
	Tmp_p = P;
	Tmp_l = L;
	Tmp_u = U;
	A.print();
	P.print();

	Tmp_p*Tmp_a;
	Tmp_l*Tmp_u;
	Tmp_p.print();
	Tmp_l.print();

	if (!A._issingular){

		L.det();
		U.det();
		cout << "a) det = " << (P._kolvop / abs(P._kolvop))* L._det * U._det << endl;
		cout << "b) rank = " << L.rank() << endl;				// т.к. rank(A*B) = rank(A);

		//Решение СЛАУ
		cout << "c) solution of liner system."<<endl<<"   vector-column b:" << endl;
		b.print();

		Matrix x(3, 1, 'x');

		solveLinerSystem(P, L, U, x);
		cout << "Result:" << endl;
		x.print();
		Matrix Tmp(3, 3, 'T');
		Tmp = P;
		Tmp*b;
		Tmp.print();

		solveLinerSystem(P, L, U, x);
		
		//Вычисление обратной матрицы
		Matrix Inv(3, 3, 'I');
		getInverseMatrix(L, U, Inv);

		cout << "   Test:" << endl;
		Matrix Test(3, 3, 'T');
		Test = A;
		Test*Inv;
		Test.print();

		//Число обусловленности матрицы
		A.getConditionNumber_1(Inv);
	}
	system("Pause");
	return 0;
}