#include <UnitTest++.h>

#include "SimpleMath/SimpleMath.h"

#include <iostream>

using namespace std;
using namespace SimpleMath;

typedef SimpleMath::Fixed::Matrix<double,3,3> Matrix33d;
typedef SimpleMath::Fixed::Matrix<double,3,1> Vector3d;
typedef SimpleMath::Dynamic::Matrix<double> MatrixNNd;
typedef SimpleMath::Dynamic::Matrix<double> VectorNd;

TEST (SimpleTestFixedQRSquare) {
	Matrix33d test_matrix (Matrix33d::Identity());

	test_matrix(0,0) = 10;
	test_matrix(1,1) = 2;
	test_matrix(1,2) = 4;
	test_matrix(2,1) = 3;
	test_matrix(2,2) = 1;

	Vector3d x;
	x[0] = 1.;
	x[1] = 2.;
	x[2] = 3.;

	Vector3d rhs = test_matrix * x;
	Vector3d x_qr =test_matrix.householderQr().solve (rhs);

	CHECK_ARRAY_CLOSE (x.data(), x_qr.data(), 3, 1.0e-14);
}

TEST (SimpleTestFixedPivotingQRSquare) {
	Matrix33d test_matrix;
	test_matrix <<
		1., 2., 3.,
		4., 4., 6.,
		8., 9., 7.;

	Vector3d x;
	x[0] = 1.;
	x[1] = 2.;
	x[2] = 3.;

	Vector3d rhs  = test_matrix * x;
	Vector3d x_qr =test_matrix.colPivHouseholderQr().solve (rhs);

	CHECK_ARRAY_CLOSE (x.data(), x_qr.data(), 3, 1.0e-14);
}

TEST (SimpleTestFixedQRRectangular) {
	SimpleMath::Fixed::Matrix<double, 6, 3> test_matrix (SimpleMath::Fixed::Matrix<double, 6, 3>::Zero());

	test_matrix(0,0) = 10;
	test_matrix(1,1) = 2;
	test_matrix(1,2) = 4;
	test_matrix(2,1) = 3;
	test_matrix(2,2) = 1;
	test_matrix(3,2) = 1;
	test_matrix(4,1) = 1;

	test_matrix(3,0) = 4.;
	test_matrix(3,1) = 3.;
	test_matrix(4,1) = 2.;
	test_matrix(4,2) = 1.;
	test_matrix(5,2) = 5.;

	Vector3d x;
	x[0] = 1.;
	x[1] = 2.;
	x[2] = 3.;

	SimpleMath::Fixed::Matrix<double, 6, 1> rhs = test_matrix * x;

	Vector3d x_qr = test_matrix.householderQr().solve (rhs);

	CHECK_ARRAY_CLOSE (x.data(), x_qr.data(), 3, 1.0e-14);
}

TEST (SimpleTestFixedColPivHouseholderQRSimple) {
	Matrix33d test_matrix (Matrix33d::Identity());

	test_matrix(0,0) = 10;
	test_matrix(1,1) = 2;
	test_matrix(1,2) = 4;
	test_matrix(2,1) = 3;
	test_matrix(2,2) = 1;

	Vector3d x;
	x[0] = 1.;
	x[1] = 2.;
	x[2] = 3.;

	Vector3d rhs = test_matrix * x;
	Vector3d x_qr =test_matrix.colPivHouseholderQr().solve (rhs);

	CHECK_ARRAY_CLOSE (x.data(), x_qr.data(), 3, 1.0e-14);
}

TEST (SimpleTestFixedColPivHouseholderQRRank) {
	Matrix33d test_matrix (Matrix33d::Identity());

	test_matrix(0,0) = 10;
	test_matrix(1,1) = 2;
	test_matrix(1,2) = 4;
	test_matrix(2,1) = 1;
	test_matrix(2,2) = 2.0001;

	ColPivHouseholderQR<Matrix33d> qr = test_matrix.colPivHouseholderQr();
	qr.setThreshold(1.0e-5);

	unsigned int rank = qr.rank();

	CHECK_EQUAL (2u, rank);
}

TEST (SimpleTestDynamicColPivHouseholderQRRank) {
	MatrixNNd test_matrix (3, 3);

	test_matrix(0,0) = 10;
	test_matrix(1,1) = 2;
	test_matrix(1,2) = 4;
	test_matrix(2,1) = 1;
	test_matrix(2,2) = 2.0001;

	ColPivHouseholderQR<MatrixNNd> qr = test_matrix.colPivHouseholderQr();
	qr.setThreshold(1.0e-5);

	unsigned int rank = qr.rank();

	CHECK_EQUAL (2u, rank);
}

TEST (SimpleMathQRInverseSimple) {
	MatrixNNd test_matrix (3,3);
	test_matrix <<
		1., 2., 3.,
		4., 4., 6.,
		8., 9., 7.;

	MatrixNNd inverse = test_matrix.householderQr().inverse();
	
	MatrixNNd product = inverse * test_matrix;

	MatrixNNd identity (MatrixNNd::Identity(3,3));

	CHECK_ARRAY_CLOSE (identity.data(), product.data(), 3 * 3, 1.0e-12);
}

TEST (SimpleMathPivotingQRInverseSimple) {
	MatrixNNd test_matrix (3,3);
	test_matrix <<
		1., 2., 3.,
		4., 4., 6.,
		8., 9., 7.;

	MatrixNNd inverse = test_matrix.colPivHouseholderQr().inverse();
	
	MatrixNNd product = inverse * test_matrix;

	MatrixNNd identity (MatrixNNd::Identity(3,3));

	CHECK_ARRAY_CLOSE (identity.data(), product.data(), 3 * 3, 1.0e-12);
}
