#include <UnitTest++.h>
#include "SimpleMath/SimpleMathFixed.h"
#include "SimpleMath/SimpleMathCommaInitializer.h"
#include <iostream>

using namespace std;
using namespace SimpleMath;

typedef SimpleMath::Fixed::Matrix<double, 6, 6> Matrix66d;
typedef SimpleMath::Fixed::Matrix<double, 3, 3> Matrix33d;

TEST (SimpleTestBlockSimple) {
	Matrix66d input (
			1, 2, 3, 4, 5, 6,
			7, 8, 9, 0, 1, 2,
			3, 4, 5, 6, 7, 8,
			9, 0, 1, 2, 3, 4,
			5, 6, 7, 8, 9, 0,
			1, 2, 3, 4, 5, 6
	);

	Matrix33d block = input.block(2, 2,3,3);

	Matrix33d reference (
			5, 6, 7,
			1, 2, 3,
			7, 8, 9
			);

	CHECK_ARRAY_EQUAL (reference.data(), block.data(), 3 * 3);
}

TEST (SimpleTestBlockTranspose) {
	Matrix66d input (
			1, 0, 0, 0, 0, 0,
			0, 1, 0, 0, 0, 0,
			0, 0, 1, 0, 0, 0,
			-0, -0, -0, 1, 0, 0,
			-0, -0, 1.1, 0, 1, 0,
			-0, -1.1, -0, 0, 0, 1
	);

	Matrix66d res(input);

	res.block<3,3>(0,0) = input.block<3,3>(0,0).transpose();
	res.block<3,3>(3,0) = input.block<3,3>(3,0).transpose();
	res.block<3,3>(0,3) = input.block<3,3>(0,3).transpose();
	res.block<3,3>(3,3) = input.block<3,3>(3,3).transpose();

	Matrix66d reference (
			1, 0, 0, 0, 0, 0,
			0, 1, 0, 0, 0, 0,
			0, 0, 1, 0, 0, 0,
			0, 0, 0, 1, 0, 0,
			0, 0, -1.1, 0, 1, 0,
			0, 1.1, 0, 0, 0, 1
			);

	CHECK_ARRAY_EQUAL (res.data(), reference.data(), res.size());
}

TEST (SimpleMathBlockVectorTranspose) {
	Matrix33d input;
	
	input <<	1., 2., 3.,
			4., 5., 6.,
			7., 8., 9.;
			
	CHECK_EQUAL (3., (input.block<3, 1>(0,2))(0,0));
	CHECK_EQUAL (6., (input.block<3, 1>(0,2))(1,0));
	CHECK_EQUAL (9., (input.block<3, 1>(0,2))(2,0));

	CHECK_EQUAL (3., (input.block<3, 1>(0,2).transpose())(0,0));
	CHECK_EQUAL (6., (input.block<3, 1>(0,2).transpose())(0,1));
	CHECK_EQUAL (9., (input.block<3, 1>(0,2).transpose())(0,2));
}
