#include <UnitTest++.h>
#include "SimpleMath/SimpleMath.h"
#include <iostream>

using namespace std;
using namespace SimpleMath;

typedef SimpleMath::Fixed::Matrix<double, 3, 3> Matrix33d;
typedef SimpleMath::Dynamic::Matrix<double> MatrixXd;

TEST (SimpleTestCommaInitializerFixed) {
	Matrix33d input;
	input << 
		1., 2., 3.,
		4., 5., 6.,
		7., 8., 9.;

	Matrix33d reference (
			1., 2., 3.,
			4., 5., 6.,
			7., 8., 9.
			);

	CHECK_ARRAY_EQUAL (reference.data(), input.data(), input.size());
}

TEST (SimpleTestCommaInitializerDynamic) {
	MatrixXd input (3, 3);
	input << 
		1., 2., 3.,
		4., 5., 6.,
		7., 8., 9.;

	Matrix33d reference (
			1., 2., 3.,
			4., 5., 6.,
			7., 8., 9.
			);

	CHECK_ARRAY_EQUAL (reference.data(), input.data(), input.size());
}


