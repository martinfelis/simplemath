#include "catch.hpp"
#include "../include/SimpleMath/SimpleMath.h"

#include <iostream>

using namespace std;
using namespace SimpleMath;

#define CHECK_ARRAY_CLOSE2 ( arr_a, arr_b, length, tol ) \
	do { \
		for (size_t array_check_i = 0; array_check_i < length; ++array_check_i) \
			REQUIRE( fabs (arr_a[array_check_i] - arr_b[array_check_i]) == Approx (0).epsilon(tol) ) \
	} while ( false )

template <typename ScalarType>
bool CHECK_ARRAY_CLOSE (const ScalarType* expected, const ScalarType* actual, size_t length, ScalarType tol) {
	for (size_t i = 0; i < length; i++) {
		cout << "i = " << i << " expected: " << expected[i] << " actual: " << actual[i] << endl;
		REQUIRE ( fabs(expected[i] - actual[i]) == Approx(0.0).epsilon(tol));
	}

	return true;
}

/*
TEST_CASE ("Basic SimpleMath works", "[SimpleMath]") {
	Matrix33f bla (Matrix33f::Identity());

	Vector3f x (1.0f, 2.0f, 3.0f);

	bla.block<3,1>(0,0) = x * (1.0 / 10.0f);
	
	cout << bla << endl;
}

TEST_CASE ("SimpleMatrixAdd", "[SimpleMath]") {
	Fixed<double, 3, 3> mat1;
	Dynamic<double> mat2(3,3);

	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			mat1(i,j) = (i + 1) * (j + 1);
		}
		mat2(i,i) = 1;
	}

	Fixed<double, 3, 3> sum_fixed_result = mat1 + mat2;
	double array_result[] = { 2., 2., 3., 2., 5., 6., 3., 6., 10.};
	CHECK_ARRAY_CLOSE(array_result, sum_fixed_result.data(), 9, 1.0e-12);

	Dynamic<double> sum_dynamic_result = mat1 + mat2;
	CHECK_ARRAY_CLOSE(array_result, sum_dynamic_result.data(), 9, 1.0e-12);
}

TEST_CASE ("SimpleMatrixValuesConstructor", "[SimpleMath]" ) {
	Fixed<double, 4, 1> vector (1.0, 2.0, 3.0, 4.0);

	double array_result[] = { 1.0, 2.0, 3.0, 4.0 };
	CHECK_ARRAY_CLOSE (array_result, vector.data(), 4, 1.0e-12);
}

*/

TEST_CASE ("SimpleMatrixMul", "[SimpleMath]") {
	Fixed<double, 3, 3> mat1;
	Dynamic<double> mat2(3,3);

	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			mat1(i,j) = (i + 1) * (j + 1);
		}
		mat2(i,i) = i + 1;
	}

	Fixed<double, 3, 3> sum_fixed_result = mat1 * mat2;
	double array_result[] = { 1., 4., 9., 2., 8., 18., 3., 12., 27.};
	CHECK_ARRAY_CLOSE(array_result, sum_fixed_result.data(), 9, 1.0e-12);

//	cout << "fixed result: " << endl << sum_fixed_result << endl;

	Dynamic<double> sum_dynamic_result = mat1 * mat2;
	CHECK_ARRAY_CLOSE(array_result, sum_dynamic_result.data(), 9, 1.0e-12);
}

TEST_CASE ("SimpleMatrixBlock", "[SimpleMath]") {
	Fixed<double, 6, 6> mat_fixed;
	Dynamic<double> mat_dynamic(6,6);

	for (size_t i = 0; i < mat_fixed.rows(); i++) {
		for (size_t j = 0; j < mat_fixed.cols(); j++) {
			mat_fixed(i,j) = (i + 1) * (j + 1);
			mat_dynamic(i,j) = (i + 1) * (j + 1);
		}
	}

	cout << mat_fixed << endl;
	cout << mat_fixed.block<3,3>(1,1) << endl;

	cout << mat_dynamic << endl;
	cout << mat_dynamic.block<2,2>(1,1) << endl;
	cout << mat_fixed.block<3,3>(1,1).block<2,2>(0,0) << endl;

	cout << "Assignment: dynamic = dynamic.block" << endl;
	Dynamic<double> mat_dynamic_block(2,2);
	mat_dynamic_block = mat_dynamic.block<2,2>(1,1);
	cout << mat_dynamic_block << endl;

	cout << "Assignment: dynamic = fixed.block" << endl;
	mat_dynamic_block = mat_fixed.block<2,2>(1,1);
	cout << mat_dynamic_block << endl;

	cout << "Add blocks:" << endl;
	cout << mat_dynamic.block<2,2>(1,1) + mat_fixed.block<2,2>(1,1) << endl;
}

TEST_CASE ("SimpleMatrixTranspose", "[SimpleMath]") {
	Fixed<double, 6, 6> mat_fixed;
	Dynamic<double> mat_dynamic(6,6);

	for (size_t i = 0; i < mat_fixed.rows(); i++) {
		for (size_t j = 0; j < mat_fixed.cols(); j++) {
			mat_fixed(i,j) = (i + 1) * (j + 3);
		}
	}

    CHECK(mat_fixed.transpose().transpose() == mat_fixed);
}

TEST_CASE ("SimpleMatrixAssignment", "[SimpleMath]") {
	Fixed<double, 6, 6> mat_fixed;
	Dynamic<double> mat_dynamic(3,3);

	for (size_t i = 0; i < mat_fixed.rows(); i++) {
		for (size_t j = 0; j < mat_fixed.cols(); j++) {
			mat_fixed(i,j) = (i + 1) * (j + 3);
		}
	}

	mat_dynamic(0,0) = 1.0; mat_dynamic(0,1) = 2.0; mat_dynamic(0,2) = 3.0;
	mat_dynamic(1,0) = 4.0; mat_dynamic(1,1) = 5.0; mat_dynamic(1,2) = 6.0;
	mat_dynamic(2,0) = 7.0; mat_dynamic(2,1) = 8.0; mat_dynamic(2,2) = 9.0;

	mat_fixed.block<3,3>(3,3).transpose().transpose() = mat_dynamic;
	REQUIRE (	(mat_fixed.block<3,3>(3,3).transpose().transpose()) == mat_dynamic);
}

TEST_CASE ("SimpleMatrixBlockComparison", "[SimpleMath]") {
	Fixed<double, 6, 6> mat_fixed;
	Dynamic<double> mat_dynamic(3,3);

	for (size_t i = 0; i < mat_fixed.rows(); i++) {
		for (size_t j = 0; j < mat_fixed.cols(); j++) {
			mat_fixed(i,j) = (i + 1) * (j + 3);
		}
	}

	mat_dynamic(0,0) = 1.0; mat_dynamic(0,1) = 2.0; mat_dynamic(0,2) = 3.0;
	mat_dynamic(1,0) = 4.0; mat_dynamic(1,1) = 5.0; mat_dynamic(1,2) = 6.0;
	mat_dynamic(2,0) = 7.0; mat_dynamic(2,1) = 8.0; mat_dynamic(2,2) = 9.0;

	// Blocks using block<nrows,ncols>(row,col) construction
	mat_fixed.block<3,3>(3,3).transpose().transpose() = mat_dynamic;
	REQUIRE (	(mat_fixed.block<3,3>(3,3).transpose().transpose()) == mat_dynamic);

	// Blocks using block(row,col,nrows,ncols)
	mat_fixed.block(3,3,3,3).transpose().transpose() = mat_dynamic;
	REQUIRE (	(mat_fixed.block(3,3,3,3).transpose().transpose()) == mat_dynamic);
}

TEST_CASE ("SimpleMatrixMultiplyScalar", "[SimpleMath]") {
	Dynamic<double> mat_dynamic(3,3);
  Fixed<double,3,3> mat_fixed;

	mat_dynamic(0,0) = 1.0; mat_dynamic(0,1) = 2.0; mat_dynamic(0,2) = 3.0;
	mat_dynamic(1,0) = 4.0; mat_dynamic(1,1) = 5.0; mat_dynamic(1,2) = 6.0;
	mat_dynamic(2,0) = 7.0; mat_dynamic(2,1) = 8.0; mat_dynamic(2,2) = 9.0;

  mat_fixed = mat_dynamic;

	Dynamic<double> mult = mat_dynamic * 3.;

	for (int i = 0, nr = mat_dynamic.rows(); i < nr; ++i)
		for (int j = 0, nc = mat_dynamic.cols(); j < nc; ++j)
			REQUIRE ((mult(i,j)) == (mat_dynamic(i,j) * 3.0));

	mult = 2.0 * mat_dynamic;

	for (int i = 0, nr = mat_dynamic.rows(); i < nr; ++i) {
		for (int j = 0, nc = mat_dynamic.cols(); j < nc; ++j) {
			REQUIRE ((mult(i,j)) == (mat_dynamic(i,j) * 2.0));
        }
    }

    // fixed * scalar
    Fixed<double,3,3> mult_fixed = mat_fixed * 3.;
    for (int i = 0, nr = mat_dynamic.rows(); i < nr; ++i)
        for (int j = 0, nc = mat_dynamic.cols(); j < nc; ++j)
                    REQUIRE ((mult_fixed(i,j)) == (mat_dynamic(i,j) * 3.0));

}

TEST_CASE ("SimpleMatrixCommaInitializer", "[SimpleMath]") {
	Dynamic<double> mat_dynamic(3,3);

	mat_dynamic(0,0) = 1.0; mat_dynamic(0,1) = 2.0; mat_dynamic(0,2) = 3.0;
	mat_dynamic(1,0) = 4.0; mat_dynamic(1,1) = 5.0; mat_dynamic(1,2) = 6.0;
	mat_dynamic(2,0) = 7.0; mat_dynamic(2,1) = 8.0; mat_dynamic(2,2) = 9.0;

	Dynamic<double> mat_dynamic_comma_initializer(3,3);
	
	mat_dynamic_comma_initializer << 
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0;

	for (int i = 0, nr = mat_dynamic.rows(); i < nr; ++i)
		for (int j = 0, nc = mat_dynamic.cols(); j < nc; ++j)
			REQUIRE ((mat_dynamic(i,j)) == (mat_dynamic_comma_initializer(i,j)));

	Fixed<double, 3, 3> mat_fixed_comma_initializer;
	mat_fixed_comma_initializer <<
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0;

	for (int i = 0, nr = mat_dynamic.rows(); i < nr; ++i)
		for (int j = 0, nc = mat_dynamic.cols(); j < nc; ++j)
			REQUIRE ((mat_dynamic(i,j)) == (mat_fixed_comma_initializer(i,j)));
}

TEST_CASE ("SimpleMathConstTranspose", "[SimpleMath]") {
	Dynamic<double> mat_dynamic(3,3);

	mat_dynamic(0,0) = 1.0; mat_dynamic(0,1) = 2.0; mat_dynamic(0,2) = 3.0;
	mat_dynamic(1,0) = 4.0; mat_dynamic(1,1) = 5.0; mat_dynamic(1,2) = 6.0;
	mat_dynamic(2,0) = 7.0; mat_dynamic(2,1) = 8.0; mat_dynamic(2,2) = 9.0;

	const Block<Dynamic<double>, double, 3, 1> const_block(&mat_dynamic, 0, 0);

	Fixed<double, 3, 1> vec1;
	vec1 = const_block;
	Dynamic<double> res1;

	cout << "Result:" << endl;
	res1 = vec1.transpose() * vec1;
	cout << res1 << endl;
	cout << vec1.transpose() * vec1 << endl;
	cout << const_block.transpose() * const_block << endl;
	float result = const_block.transpose() * const_block;
	cout << result << endl;
}


TEST_CASE ("SimpleMathUnifiedFixedDynamic", "[SimpleMath]") {
	Matrix<double, 3, 3>  fixed_mat33;
	Matrix<float> dynamic_mat33 (3, 3);

	dynamic_mat33 <<
			1.0, 2.0, 3.0,
			4.0, 5.0, 6.0,
			7.0, 8.0, 9.0;

	fixed_mat33 = dynamic_mat33;

	cout << "matrix " << sizeof(Matrix<double, 3, 3>) << endl;

	cout << fixed_mat33 << endl;
	cout << dynamic_mat33 << endl;

	cout << dynamic_mat33(0, 0) << ", " << dynamic_mat33(0, 1) << ", " << dynamic_mat33(0, 2) << endl
			<< dynamic_mat33(1, 0) << ", " << dynamic_mat33(1, 1) << ", " << dynamic_mat33(1, 2) << endl
			<< dynamic_mat33(2, 0) << ", " << dynamic_mat33(2, 1) << ", " << dynamic_mat33(2, 2) << endl;

	cout << sizeof(Matrix<float, 1, 1>) << endl;
    cout << sizeof(Matrix<float>) << endl;
}