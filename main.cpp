#include "examples.h"
using namespace std;


int main(){

	initializeExample(4);
	Matrix A(N, N, 'A');
	A = M;
	Matrix B(N, N, 'B');
	Matrix P(N, N, 'P');
	Matrix L(N, N, 'L');
	Matrix U(N, N, 'U');
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

		Matrix x(N, 1, 'x');
		Matrix tmp(N, N, 't');
		tmp = P;
		tmp*b;
		solveLinerSystem(tmp, L, U, x);
		cout << "Result:" << endl;
		x.print();
		//Вычисление обратной матрицы
		Matrix Inv(N, N, 'I');
		getInverseMatrix(L, U, Inv, A._issingular);

		//Inv.swapColumn(0, 1);
		//Inv.swapColumn(1, 2);
		//Inv.print();

		cout << "   Test:" << endl;
		Matrix Test(N, N, 'T');
		A.print();
		Test = A;
		P*Test*Inv;
		P.print();


		//Число обусловленности матрицы
		A.getConditionNumber_1(Inv);
	}
	else{ 
		int num_sol = GaussJordanElimination(A, b);
		cout << "Matrix A after Gaus-Jordan:" << endl;
		A.print();
		cout << "a) det = " << 0 << endl;
		cout << "b) rank = " << A.rank() << endl;
		cout << "c) solution:" << endl;
		if (num_sol == 0){
			cout << "System has not a solution" << endl;
		}
		else{
			cout << "System has endlessly many solutions" << endl;
			b.print();
		}
		cout << "d) Inverse matrix is not exist" << endl;
	}
	

	system("Pause");
	return 0;
}