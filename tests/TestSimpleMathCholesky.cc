#include <UnitTest++.h>

#include "SimpleMath/SimpleMath.h"
#include "SimpleMath/SimpleMathCholesky.h"

#include <iostream>

using namespace std;
using namespace SimpleMath;

typedef SimpleMath::Fixed::Matrix<double,3,3> Matrix33d;
typedef SimpleMath::Fixed::Matrix<double,3,1> Vector3d;
typedef SimpleMath::Dynamic::Matrix<double> MatrixNNd;
typedef SimpleMath::Dynamic::Matrix<double> VectorNd;

TEST (SimpleTestFixedLLTSimple) {
	Matrix33d test_matrix (Matrix33d::Identity());

	test_matrix <<
		4, 12, -16,
		12, 37, -43,
		-16, -43, 98;

	LLT<Matrix33d> llt = test_matrix.llt();	

	Matrix33d L_LT = llt.matrixL() * llt.matrixL().transpose();

	CHECK_ARRAY_EQUAL (test_matrix.data(), L_LT.data(), 9);
}

TEST (SimpleTestFixedLLTSolveSimple) {
	Matrix33d test_matrix (Matrix33d::Identity());

	test_matrix <<
		4, 12, -16,
		12, 37, -43,
		-16, -43, 98;

	Vector3d x;
	x <<
		1.,
		2.,
		3;

	Vector3d b = test_matrix * x;

	Vector3d x_dash = test_matrix.llt().solve(b);

	CHECK_ARRAY_EQUAL (x.data(), x_dash.data(), 3);
}
