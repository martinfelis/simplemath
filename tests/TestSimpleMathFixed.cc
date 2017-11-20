#include <UnitTest++.h>
#include "SimpleMath/SimpleMathFixed.h"
#include <iostream>

using namespace std;
using namespace SimpleMath;

typedef SimpleMath::Fixed::Matrix<double, 6, 6> Matrix66d;
typedef SimpleMath::Fixed::Matrix<double, 4, 4> Matrix44d;
typedef SimpleMath::Fixed::Matrix<double, 3, 3> Matrix33d;
typedef SimpleMath::Fixed::Matrix<double, 4, 1> Vector4d;
typedef SimpleMath::Fixed::Matrix<double, 3, 1> Vector3d;

TEST (SimpleTestFixed) {
	Matrix44d mymatrix;
	Vector4d myvector;

	mymatrix.identity();
	myvector.random();

	Vector4d result = mymatrix * myvector;

	CHECK_ARRAY_EQUAL (myvector.data(), result, myvector.size());
}

TEST (FixedTestBlockFullMatrix) {
	Matrix44d mymatrix;
	mymatrix.identity();

	Matrix44d othermatrix = Matrix44d::Zero (4,4);

	mymatrix.block<4,4>(0,0) = othermatrix.block<4,4>(0,0);

	CHECK_EQUAL (mymatrix, othermatrix);
}

TEST (FixedTestBlockToSmaller) {
	Matrix66d matrix66;
	matrix66.identity();

	Matrix33d othermatrix = Matrix33d::Zero (3,3);
	Matrix33d identity33 = Matrix33d::Identity (3,3);

	othermatrix = matrix66.block<3,3>(1,1);

	CHECK_EQUAL (identity33, othermatrix);
}

Matrix33d assign_from_const (const Matrix66d &mat6) {
	Matrix33d result;
	result = mat6.block<3,3>(0,0);

	return result;
}

TEST (FixedTestBlockToSmallerFromConst) {
	Matrix66d matrix66;
	matrix66.identity();

	Matrix33d othermatrix = assign_from_const (matrix66);
	Matrix33d identity33 = Matrix33d::Identity (3,3);

	CHECK_EQUAL (identity33, othermatrix);
}


TEST (FixedTestBlockPartial) {
	Matrix44d mymatrix;
	mymatrix.setZero();
	mymatrix(0,0) = 1.;

	Matrix33d othermatrix = Matrix33d::Identity (3,3);

	mymatrix.block<3,3>(1,1) = othermatrix;

	CHECK_EQUAL (Matrix44d::Identity (4,4),mymatrix);
}

TEST (FixedTestBlock) {
	Matrix44d matrix;

	matrix(0,0) = 1.;
	matrix(0,1) = 2.;
	matrix(0,2) = 3.;
	matrix(0,3) = 4.;
	
	matrix(1,0) = 5.;
	matrix(1,1) = 6.;
	matrix(1,2) = 7.;
	matrix(1,3) = 8.;

	matrix(2,0) = 9.;
	matrix(2,1) = 10.;
	matrix(2,2) = 11.;
	matrix(2,3) = 12.;

	matrix(3,0) = 13.;
	matrix(3,1) = 14.;
	matrix(3,2) = 15.;
	matrix(3,3) = 16.;

	Matrix44d result_matrix;
	result_matrix(0,0) = 1.;
	result_matrix(0,1) = 2.;
	result_matrix(0,2) = 3.;
	result_matrix(0,3) = 4.;
	
	result_matrix(1,0) = 5.;
	result_matrix(1,1) = 6.;
	result_matrix(1,2) = 10.;
	result_matrix(1,3) = 8.;

	result_matrix(2,0) = 9.;
	result_matrix(2,1) = 7.;
	result_matrix(2,2) = 11.;
	result_matrix(2,3) = 12.;

	result_matrix(3,0) = 13.;
	result_matrix(3,1) = 14.;
	result_matrix(3,2) = 15.;
	result_matrix(3,3) = 16.;

	matrix.block<2,2>(1,1) = result_matrix.block<2,2>(1,1);

	CHECK_EQUAL (result_matrix, matrix);
}

TEST (FixedTestBlockTransposedSelf) {
	Matrix44d matrix;

	matrix(0,0) = 1.;
	matrix(0,1) = 2.;
	matrix(0,2) = 3.;
	matrix(0,3) = 4.;
	
	matrix(1,0) = 5.;
	matrix(1,1) = 6.;
	matrix(1,2) = 7.;
	matrix(1,3) = 8.;

	matrix(2,0) = 9.;
	matrix(2,1) = 10.;
	matrix(2,2) = 11.;
	matrix(2,3) = 12.;

	matrix(3,0) = 13.;
	matrix(3,1) = 14.;
	matrix(3,2) = 15.;
	matrix(3,3) = 16.;

	Matrix44d result_matrix;
	result_matrix(0,0) = 1.;
	result_matrix(0,1) = 2.;
	result_matrix(0,2) = 3.;
	result_matrix(0,3) = 4.;
	
	result_matrix(1,0) = 5.;
	result_matrix(1,1) = 6.;
	result_matrix(1,2) = 10.;
	result_matrix(1,3) = 8.;

	result_matrix(2,0) = 9.;
	result_matrix(2,1) = 7.;
	result_matrix(2,2) = 11.;
	result_matrix(2,3) = 12.;

	result_matrix(3,0) = 13.;
	result_matrix(3,1) = 14.;
	result_matrix(3,2) = 15.;
	result_matrix(3,3) = 16.;

	matrix.block<2,2>(1,1) = matrix.block<2,2>(1,1).transpose();

	CHECK_EQUAL (result_matrix, matrix);
}

TEST (FixedTestBlockDoubleTranspose) {
	Matrix44d matrix;

	matrix(0,0) = 1.;
	matrix(0,1) = 2.;
	matrix(0,2) = 3.;
	matrix(0,3) = 4.;
	
	matrix(1,0) = 5.;
	matrix(1,1) = 6.;
	matrix(1,2) = 7.;
	matrix(1,3) = 8.;

	matrix(2,0) = 9.;
	matrix(2,1) = 10.;
	matrix(2,2) = 11.;
	matrix(2,3) = 12.;

	matrix(3,0) = 13.;
	matrix(3,1) = 14.;
	matrix(3,2) = 15.;
	matrix(3,3) = 16.;

	Matrix44d result_matrix;

	result_matrix(0,0) = 1.;
	result_matrix(0,1) = 2.;
	result_matrix(0,2) = 3.;
	result_matrix(0,3) = 4.;
	
	result_matrix(1,0) = 5.;
	result_matrix(1,1) = 6.;
	result_matrix(1,2) = 7.;
	result_matrix(1,3) = 8.;

	result_matrix(2,0) = 9.;
	result_matrix(2,1) = 10.;
	result_matrix(2,2) = 11.;
	result_matrix(2,3) = 12.;

	result_matrix(3,0) = 13.;
	result_matrix(3,1) = 14.;
	result_matrix(3,2) = 15.;
	result_matrix(3,3) = 16.;

	matrix.block<2,2>(1,1) = matrix.block<2,2>(1,1).transpose().transpose();

	CHECK_EQUAL (result_matrix, matrix);
}


