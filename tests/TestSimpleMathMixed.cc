#include <UnitTest++.h>

#include "SimpleMath/SimpleMath.h"

#include <iostream>

using namespace std;
using namespace SimpleMath;

typedef SimpleMath::Dynamic::Matrix<double> MatrixXd;
typedef SimpleMath::Dynamic::Matrix<double> VectorXd;
typedef SimpleMath::Fixed::Matrix<double, 3, 3> Matrix33d;
typedef SimpleMath::Fixed::Matrix<double, 3, 1> Vector3d;

TEST (SimpleTestMixed) {
	Fixed::Matrix<double, 1, 3> row_vector;
	Dynamic::Matrix<double> col_vector (3,1);

	row_vector[0] = 1.;
	row_vector[1] = 2.;
	row_vector[2] = 3.;

	col_vector[0] = 4.;
	col_vector[1] = 5.;
	col_vector[2] = 6.;

	Dynamic::Matrix<double> scalar_result = row_vector * col_vector;

	CHECK_EQUAL (1u, scalar_result.size());
	CHECK_EQUAL (1u, scalar_result.rows());
	CHECK_EQUAL (1u, scalar_result.cols());
	CHECK_EQUAL (32., scalar_result[0]);

	Dynamic::Matrix<double> outer_result = col_vector * row_vector;

	CHECK_EQUAL (9u, outer_result.size());
	CHECK_EQUAL (3u, outer_result.rows());
	CHECK_EQUAL (3u, outer_result.cols());

	Fixed::Matrix<double, 3, 3> outer_test_result;
	outer_test_result(0,0) =  4;
	outer_test_result(0,1) =  8;
	outer_test_result(0,2) = 12;

	outer_test_result(1,0) =  5;
	outer_test_result(1,1) = 10;
	outer_test_result(1,2) = 15;

	outer_test_result(2,0) =  6;
	outer_test_result(2,1) = 12;
	outer_test_result(2,2) = 18;

	CHECK_EQUAL (outer_test_result, outer_result);
}

TEST (ScalarTypeConversionFloatToDoubleDynamic) {
	SimpleMath::Dynamic::Matrix<float> fmatrix(3,3);

	fmatrix(0,0) = 1.1;
	fmatrix(0,1) = 2.2;
	fmatrix(0,2) = 3.3;
	fmatrix(1,0) = 4.4;
	fmatrix(1,1) = 5.5;
	fmatrix(1,2) = 6.6;
	fmatrix(2,0) = 7.7;
	fmatrix(2,1) = 8.8;
	fmatrix(2,2) = 9.9;

	SimpleMath::Dynamic::Matrix<double> dmatrix(fmatrix);

	CHECK_EQUAL (fmatrix.cols(), dmatrix.cols());
	CHECK_EQUAL (fmatrix.rows(), dmatrix.rows());

	CHECK_CLOSE (fmatrix(0,0), dmatrix(0,0), 1.0e-5);
	CHECK_CLOSE (fmatrix(0,1), dmatrix(0,1), 1.0e-5);
	CHECK_CLOSE (fmatrix(0,2), dmatrix(0,2), 1.0e-5);
	CHECK_CLOSE (fmatrix(1,0), dmatrix(1,0), 1.0e-5);
	CHECK_CLOSE (fmatrix(1,1), dmatrix(1,1), 1.0e-5);
	CHECK_CLOSE (fmatrix(1,2), dmatrix(1,2), 1.0e-5);
	CHECK_CLOSE (fmatrix(2,0), dmatrix(2,0), 1.0e-5);
	CHECK_CLOSE (fmatrix(2,1), dmatrix(2,1), 1.0e-5);
	CHECK_CLOSE (fmatrix(2,2), dmatrix(2,2), 1.0e-5);

	SimpleMath::Dynamic::Matrix<double> other_dmatrix;
	other_dmatrix = fmatrix;

	CHECK_EQUAL (dmatrix, other_dmatrix);
}

TEST (ScalarTypeConversionBlockDynamic) {
	SimpleMath::Dynamic::Matrix<float> fmatrix(3,3);

	fmatrix(0,0) = 1.1;
	fmatrix(0,1) = 2.2;
	fmatrix(0,2) = 3.3;
	fmatrix(1,0) = 4.4;
	fmatrix(1,1) = 5.5;
	fmatrix(1,2) = 6.6;
	fmatrix(2,0) = 7.7;
	fmatrix(2,1) = 8.8;
	fmatrix(2,2) = 9.9;

	SimpleMath::Dynamic::Matrix<double> dmatrix(2,2);
	dmatrix(0,0) = 1.;
	dmatrix(0,1) = 2.;
	dmatrix(1,0) = 4.;
	dmatrix(1,1) = 5.;

	fmatrix.block<2,2>(0,0) = dmatrix;

	CHECK_CLOSE (fmatrix(0,0), dmatrix(0,0), 1.0e-5);
	CHECK_CLOSE (fmatrix(0,1), dmatrix(0,1), 1.0e-5);
	CHECK_CLOSE (fmatrix(1,0), dmatrix(1,0), 1.0e-5);
	CHECK_CLOSE (fmatrix(1,1), dmatrix(1,1), 1.0e-5);
}

TEST (ScalarTypeConversionBlockFixed) {
	SimpleMath::Fixed::Matrix<float,3,3> fmatrix;

	fmatrix(0,0) = 1.1;
	fmatrix(0,1) = 2.2;
	fmatrix(0,2) = 3.3;
	fmatrix(1,0) = 4.4;
	fmatrix(1,1) = 5.5;
	fmatrix(1,2) = 6.6;
	fmatrix(2,0) = 7.7;
	fmatrix(2,1) = 8.8;
	fmatrix(2,2) = 9.9;

	SimpleMath::Fixed::Matrix<double,2,2> dmatrix;
	dmatrix(0,0) = 1.;
	dmatrix(0,1) = 2.;
	dmatrix(1,0) = 4.;
	dmatrix(1,1) = 5.;

	fmatrix.block<2,2>(0,0) = dmatrix;

	CHECK_CLOSE (fmatrix(0,0), dmatrix(0,0), 1.0e-5);
	CHECK_CLOSE (fmatrix(0,1), dmatrix(0,1), 1.0e-5);
	CHECK_CLOSE (fmatrix(1,0), dmatrix(1,0), 1.0e-5);
	CHECK_CLOSE (fmatrix(1,1), dmatrix(1,1), 1.0e-5);
}

TEST (ConversionFixedToDynamic) {
	Vector3d rhs_fixed;
	rhs_fixed[0] = 1.;
	rhs_fixed[1] = 2.;
	rhs_fixed[2] = 3.;

	VectorXd rhs_dynamic;

	rhs_dynamic = rhs_fixed * 2.;
	
	CHECK_EQUAL (3u, rhs_dynamic.rows());
	CHECK_EQUAL (1u, rhs_dynamic.cols());

	CHECK_EQUAL (2, rhs_dynamic[0]);
	CHECK_EQUAL (4, rhs_dynamic[1]);
	CHECK_EQUAL (6, rhs_dynamic[2]);

	VectorXd rhs_copy_constr (rhs_fixed * 2.);
	CHECK_EQUAL (rhs_dynamic, rhs_copy_constr);
}

TEST (ConversionDynamicToFixed) {
	VectorXd rhs_dynamic(3);
	rhs_dynamic[0] = 1.;
	rhs_dynamic[1] = 2.;
	rhs_dynamic[2] = 3.;

	Vector3d rhs_fixed;

	rhs_fixed = rhs_dynamic * 2.;

	CHECK_EQUAL (3u, rhs_fixed.rows());
	CHECK_EQUAL (1u, rhs_fixed.cols());

	CHECK_EQUAL (2, rhs_fixed[0]);
	CHECK_EQUAL (4, rhs_fixed[1]);
	CHECK_EQUAL (6, rhs_fixed[2]);
}

TEST (MultiplyFixedWithDynamic) {
	MatrixXd dynamic (3, 2);
	dynamic <<
		1., 2.,
		3., 4.,
		5., 6.;

	Matrix33d fixed;
	fixed <<
		1., 2., 3., 
		4., 5., 6.,
		7., 8., 9.;

	MatrixXd result = fixed * dynamic;

	MatrixXd reference (3, 2);
	reference <<
		22., 28.,
		49., 64.,
		76., 100;

	CHECK_ARRAY_EQUAL (reference.data(), result.data(), 3 * 2);
}

TEST (MultiplyDynamicWithFixed) {
	MatrixXd dynamic (2, 3);
	dynamic <<
		1., 2., 3.,
		4., 5., 6.;

	Matrix33d fixed;
	fixed <<
		1., 2., 3., 
		4., 5., 6.,
		7., 8., 9.;

	MatrixXd result = dynamic * fixed;

	MatrixXd reference (2, 3);
	reference <<
		30., 36., 42.,
		66., 81., 96.;

	CHECK_ARRAY_EQUAL (reference.data(), result.data(), 2 * 3);
}

TEST (MultiplyDynamicWithFixedStoreInFixed) {
	MatrixXd dynamic (3, 3);
	dynamic <<
		1., 2., 3., 
		4., 5., 6.,
		7., 8., 9.;

	Vector3d fixed;
	fixed <<
		1., 2., 3.;

	Vector3d result = dynamic * fixed;

	Vector3d reference;
	reference <<
		14., 32, 50;

	CHECK_ARRAY_EQUAL (reference.data(), result.data(), 3);
}
