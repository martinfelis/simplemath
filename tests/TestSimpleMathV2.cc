#include <UnitTest++.h>

#include "SimpleMath/SimpleMathBase.h"

#include <iostream>
#include <vector>

using namespace std;
using namespace SimpleMath;

TEST (SimpleMatrixAdd) {
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

TEST (SimpleMatrixMul) {
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

	Dynamic<double> sum_dynamic_result = mat1 * mat2;
	CHECK_ARRAY_CLOSE(array_result, sum_dynamic_result.data(), 9, 1.0e-12);
}

TEST (SimpleMatrixBlock) {
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

TEST (SimpleMatrixTranspose) {
	Fixed<double, 6, 6> mat_fixed;
	Dynamic<double> mat_dynamic(6,6);

	for (size_t i = 0; i < mat_fixed.rows(); i++) {
		for (size_t j = 0; j < mat_fixed.cols(); j++) {
			mat_fixed(i,j) = (i + 1) * (j + 3);
		}
	}

    CHECK(mat_fixed.transpose().transpose() == mat_fixed);
}

TEST (SimpleMatrixAssignment) {
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
	CHECK (	(mat_fixed.block<3,3>(3,3).transpose().transpose()) == mat_dynamic);
	CHECK_EQUAL (	(mat_fixed.block<3,3>(3,3).transpose().transpose()),  mat_dynamic);
}

TEST (SimpleMatrixBlockComparison) {
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
	CHECK (	(mat_fixed.block<3,3>(3,3).transpose().transpose()) == mat_dynamic);
	CHECK_EQUAL (	(mat_fixed.block<3,3>(3,3).transpose().transpose()),  mat_dynamic);

	// Blocks using block(row,col,nrows,ncols)
	mat_fixed.block(3,3,3,3).transpose().transpose() = mat_dynamic;
	CHECK (	(mat_fixed.block(3,3,3,3).transpose().transpose()) == mat_dynamic);
	CHECK_EQUAL (	(mat_fixed.block(3,3,3,3).transpose().transpose()),  mat_dynamic);
}

TEST (SimpleMatrixMultiplyScalar) {
	Dynamic<double> mat_dynamic(3,3);
    Fixed<double,3,3> mat_fixed(3,3);

	mat_dynamic(0,0) = 1.0; mat_dynamic(0,1) = 2.0; mat_dynamic(0,2) = 3.0;
	mat_dynamic(1,0) = 4.0; mat_dynamic(1,1) = 5.0; mat_dynamic(1,2) = 6.0;
	mat_dynamic(2,0) = 7.0; mat_dynamic(2,1) = 8.0; mat_dynamic(2,2) = 9.0;

    mat_fixed = mat_dynamic;

	Dynamic<double> mult = mat_dynamic * 3.;

	for (int i = 0, nr = mat_dynamic.rows(); i < nr; ++i)
		for (int j = 0, nc = mat_dynamic.cols(); j < nc; ++j)
			CHECK_EQUAL ((mult(i,j)), (mat_dynamic(i,j) * 3.0));

	mult = 2.0 * mat_dynamic;

	for (int i = 0, nr = mat_dynamic.rows(); i < nr; ++i) {
		for (int j = 0, nc = mat_dynamic.cols(); j < nc; ++j) {
			CHECK_EQUAL ((mult(i,j)), (mat_dynamic(i,j) * 2.0));
        }
    }

    // fixed * scalar
    Fixed<double,3,3> mult_fixed = mat_fixed * 3.;
    for (int i = 0, nr = mat_dynamic.rows(); i < nr; ++i)
        for (int j = 0, nc = mat_dynamic.cols(); j < nc; ++j)
                    CHECK_EQUAL ((mult_fixed(i,j)), (mat_dynamic(i,j) * 3.0));

}

TEST (SimpleMatrixCommaInitializer) {
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
			CHECK_EQUAL ((mat_dynamic(i,j)), (mat_dynamic_comma_initializer(i,j)));

	Fixed<double, 3, 3> mat_fixed_comma_initializer;
	mat_fixed_comma_initializer <<
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0;

	for (int i = 0, nr = mat_dynamic.rows(); i < nr; ++i)
		for (int j = 0, nc = mat_dynamic.cols(); j < nc; ++j)
			CHECK_EQUAL ((mat_dynamic(i,j)), (mat_fixed_comma_initializer(i,j)));
}

