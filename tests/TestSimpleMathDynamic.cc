#include <UnitTest++.h>
#include "SimpleMath/SimpleMathDynamic.h"
#include <iostream>

using namespace std;
using namespace SimpleMath;

typedef SimpleMath::Dynamic::Matrix<double> MatrixXd;
typedef SimpleMath::Dynamic::Matrix<double> VectorXd;

TEST (SimpleTestDynamic) {
	MatrixXd mymatrix (4,4);
	MatrixXd myvector (4);

	mymatrix.identity();
	myvector.random();

	MatrixXd result = MatrixXd::Zero (3);
  result = MatrixXd::Zero (3,3);
	result.resize(4);
	result = mymatrix * myvector;

	CHECK_EQUAL (4u, result.size());
	CHECK_EQUAL (4u, result.rows());
	CHECK_EQUAL (1u, result.cols());
	CHECK_ARRAY_EQUAL (myvector.data(), result, myvector.size());
}

TEST (SimpleTestAddScalarMult) {
	MatrixXd myvector (4), res_vector(4);
	myvector.random();

	res_vector = myvector + myvector;
	myvector = myvector * 2.;

	CHECK_ARRAY_CLOSE (myvector.data(), res_vector.data(), 4, 1.0e-12);
}
TEST (DynamicTestBlockFullMatrix) {
	MatrixXd mymatrix = MatrixXd::Identity(4,4);
	MatrixXd othermatrix = MatrixXd::Zero (4,4);

	mymatrix.block<4,4>(0,0) = othermatrix;

	CHECK_EQUAL (mymatrix, othermatrix);
}

TEST (DynamicTestBlockPartialToBorderMatrix) {
	MatrixXd mymatrix = MatrixXd::Zero(4,4);
	mymatrix(3,3) = 1.;

	MatrixXd othermatrix = MatrixXd::Identity (3,3);
	mymatrix.block(0,0,3,3) = othermatrix;

	CHECK_EQUAL (mymatrix, MatrixXd::Identity (4,4));
}
