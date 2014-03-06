#include "matrix.h"

Matrix M(3, 3, 'M');
Matrix b(3, 1, 'b');

void initializeExample(const int n){

	switch (n){
	case 0:{
			   double arr1[3][3] = { { 2, 1, 1 },			// det = -4
			   { 1, -1, 0 },
			   { 3, -1, 2 } };
			   Matrix A(3, 3, &arr1[0][0], 'M');
			   M = A;

			   double arr2[3][1] = { { 2.0 }, { -2.0 }, { 2.0 } };
			   Matrix x(3, 1, &arr2[0][0], 'x');
			   b = x;
	}break;
	case 1:{
			   double arr[3][3] = { { 2, 7, -6 },
			   { 8, 2, 1 },
			   { 7, 4, 2 } };			// det = -171
	}break;
	case 2:{
			   double arr1[3][3] = { { 2, 1, 1 },			// det = 
			   { 2, 1, 1 },
			   { 2, 1, 1 } };
			   Matrix A(3, 3, &arr1[0][0], 'M');
			   M = A;

			   double arr2[3][1] = { { 2.0 }, { 2.0 }, { 2.0 } };
			   Matrix x(3, 1, &arr2[0][0], 'x');
			   b = x;
	}break;
	case 3:{
			   double arr1[3][3] = { { 2, 1, 1 },			// det = -4
			   { 2, 1, 1 },
			   { 3, -1, 2 } };
			   Matrix A(3, 3, &arr1[0][0], 'M');
			   M = A;

			   double arr2[3][1] = { { 2.0 }, { -2.0 }, { 2.0 } };
			   Matrix x(3, 1, &arr2[0][0], 'x');
			   b = x;
	}
	}
	return;
}