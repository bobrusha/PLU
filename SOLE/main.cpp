#include <iostream>
#include "matrix.h"

int main(){
	vector < vector < double > > x, y;
	x.resize(3);

	for (int i = 0; i < 3; ++i){
		x[i].resize(3);
	}

	x[0][0] = 0;	x[0][1] = 1;	x[0][2] = 2;
	x[1][0] = 0;	x[1][1] = 1;	x[1][2] = 3;
	x[2][0] = 1;	x[2][1] = 2;	x[2][2] = 4;

	print(x);

	vector < vector <double> > P, L, U, LU;
	P.resize(3);
	LU.resize(3);
	L.resize(3);
	U.resize(3);

	for (int i = 0; i < 3; ++i){
		P[i].resize(3);
		LU[i].resize(3);
		L[i].resize(3);
		U[i].resize(3);
	}
	PLU ( x, LU, P);
	cout << "Matrix P:" << endl;
	print(P);
	fromOnetoTwo ( LU, L, U);
	
	LU = getInverseMatrix(P, L, U);


	system("Pause");
	return 0;
}

