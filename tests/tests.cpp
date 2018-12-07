#include "catch.hpp"
#include "../include/SimpleMath/SimpleMath.h"

#include <iostream>

using namespace std;
using namespace SimpleMath;

using Catch::Matchers::WithinULP;

#define CHECK_ARRAY_CLOSE2 ( arr_a, arr_b, length, tol ) \
	do { \
		for (size_t array_check_i = 0; array_check_i < length; ++array_check_i) \
			REQUIRE( fabs (arr_a[array_check_i] - arr_b[array_check_i]) == Approx (0).epsilon(tol) ) \
	} while ( false )

template <typename ScalarType>
bool CHECK_ARRAY_CLOSE (const ScalarType* expected, const ScalarType* actual, size_t length, const ScalarType& tol) {
	for (size_t i = 0; i < length; i++) {
//		cout << "i = " << i << " expected: " << expected[i] << " actual: " << actual[i] << endl;
		CHECK ( actual[i] == Approx(expected[i]).margin(tol));
	}

	return true;
}

TEST_CASE ("Basic SimpleMath works", "[SimpleMath]") {
	Matrix33f bla (Matrix33f::Identity());

	Vector3f x (1.0f, 2.0f, 3.0f);

	bla.block<3,1>(0,0) = x * (1.0 / 10.0f);
	
	cout << bla << endl;
}


TEST_CASE ("SimpleMatrixAdd", "[SimpleMath]") {
	Matrix<double, 3, 3> mat1;
	Matrix<double> mat2 (Matrix<double>::Zero(3, 3));

	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			mat1(i,j) = (i + 1) * (j + 1);
		}
		mat2(i,i) = 1;
	}

	Matrix<double, 3, 3> sum_fixed_result = mat1 + mat2;
	double array_result[] = { 2., 2., 3., 2., 5., 6., 3., 6., 10.};
	CHECK_ARRAY_CLOSE(array_result, sum_fixed_result.data(), 9, 1.0e-12);

	Matrix<double> sum_dynamic_result = mat1 + mat2;
	CHECK_ARRAY_CLOSE(array_result, sum_dynamic_result.data(), 9, 1.0e-12);
}

TEST_CASE ("SimpleMatrixValuesConstructor", "[SimpleMath]" ) {
	Matrix<double, 4, 1> vector (1.0, 2.0, 3.0, 4.0);

	double array_result[] = { 1.0, 2.0, 3.0, 4.0 };
	CHECK_ARRAY_CLOSE (array_result, vector.data(), 4, 1.0e-12);
}

TEST_CASE ("InitializeWithZeroMatrix", "[SimpleMath]") {
    Matrix<double, 3, 3> mat_fixed (Matrix<double, 3, 3>::Zero());
	Matrix<double> mat_dynamic (Matrix<double>::Zero(3, 3));

    double array_result[] = { 0., 0., 0., 0., 0., 0., 0., 0., 0.};
    CHECK_ARRAY_CLOSE(array_result, mat_fixed.data(), 9, 1.0e-12);
	CHECK_ARRAY_CLOSE(array_result, mat_dynamic.data(), 9, 1.0e-12);
}

TEST_CASE ("SimpleMatrixMul", "[SimpleMath]") {
	Matrix<double, 3, 3> mat1;
	Matrix<double> mat2(3,3);

	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			mat1(i,j) = (i + 1) * (j + 1);
		}
		mat2(i,i) = i + 1;
	}

	Matrix<double, 3, 3> sum_fixed_result = mat1 * mat2;
	double array_result[] = { 1., 4., 9., 2., 8., 18., 3., 12., 27.};
	CHECK_ARRAY_CLOSE(array_result, sum_fixed_result.data(), 9, 1.0e-12);

//	cout << "fixed result: " << endl << sum_fixed_result << endl;

	Matrix<double> sum_dynamic_result = mat1 * mat2;
	CHECK_ARRAY_CLOSE(array_result, sum_dynamic_result.data(), 9, 1.0e-12);
}


TEST_CASE ("SimpleMatrixBlock", "[SimpleMath]") {
	Matrix<double, 6, 6> mat_fixed;
	Matrix<double> mat_dynamic(6,6);

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
	Matrix<double> mat_dynamic_block(2,2);
	mat_dynamic_block = mat_dynamic.block<2,2>(1,1);
	cout << mat_dynamic_block << endl;

	cout << "Assignment: dynamic = fixed.block" << endl;
	mat_dynamic_block = mat_fixed.block<2,2>(1,1);
	cout << mat_dynamic_block << endl;

	cout << "Add blocks:" << endl;
	cout << mat_dynamic.block<2,2>(1,1) + mat_fixed.block<2,2>(1,1) << endl;
}

TEST_CASE ("SwapColumns", "[SimpleMath]") {
	Matrix<double, 3, 3> A;
	A <<
	  1., 2., 3.,
	  4., 5., 6.,
	  7., 8., 4.;

	Matrix<double, 3, 3> B;
	B <<
	  2., 1., 3.,
      5., 4., 6.,
	  8., 7., 4.;

	Matrix<double, 3, 1> column_0 = A.block(0, 0, A.rows(), 1);
	A.block(0, 0, A.rows(), 1) = A.block(0, 1, A.rows(), 1);
	A.block(0, 1, A.rows(), 1) = column_0;

	CHECK_ARRAY_CLOSE(B.data(), A.data(), 9, 1.0e-12);
}

TEST_CASE ("SimpleMatrixTranspose", "[SimpleMath]") {
	Matrix<double, 6, 6> mat_fixed;
	Matrix<double> mat_dynamic(6,6);

	for (size_t i = 0; i < mat_fixed.rows(); i++) {
		for (size_t j = 0; j < mat_fixed.cols(); j++) {
			mat_fixed(i,j) = (i + 1) * (j + 3);
		}
	}

    CHECK(mat_fixed.transpose().transpose() == mat_fixed);
}

TEST_CASE ("SimpleMatrixAssignment", "[SimpleMath]") {
	Matrix<double, 6, 6> mat_fixed;
	Matrix<double> mat_dynamic(3,3);

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
	Matrix<double, 6, 6> mat_fixed;
	Matrix<double> mat_dynamic(3,3);

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
	Matrix<double> mat_dynamic(3,3);
  	Matrix<double,3,3> mat_fixed;

	mat_dynamic(0,0) = 1.0; mat_dynamic(0,1) = 2.0; mat_dynamic(0,2) = 3.0;
	mat_dynamic(1,0) = 4.0; mat_dynamic(1,1) = 5.0; mat_dynamic(1,2) = 6.0;
	mat_dynamic(2,0) = 7.0; mat_dynamic(2,1) = 8.0; mat_dynamic(2,2) = 9.0;

  	mat_fixed = mat_dynamic;

	Matrix<double> mult = mat_dynamic * 3.;

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
    Matrix<double,3,3> mult_fixed = mat_fixed * 3.;
    for (int i = 0, nr = mat_dynamic.rows(); i < nr; ++i)
        for (int j = 0, nc = mat_dynamic.cols(); j < nc; ++j)
                    REQUIRE ((mult_fixed(i,j)) == (mat_dynamic(i,j) * 3.0));

}

TEST_CASE ("SimpleMatrixCommaInitializer", "[SimpleMath]") {
	Matrix<double> mat_dynamic(3,3);

	mat_dynamic(0,0) = 1.0; mat_dynamic(0,1) = 2.0; mat_dynamic(0,2) = 3.0;
	mat_dynamic(1,0) = 4.0; mat_dynamic(1,1) = 5.0; mat_dynamic(1,2) = 6.0;
	mat_dynamic(2,0) = 7.0; mat_dynamic(2,1) = 8.0; mat_dynamic(2,2) = 9.0;

	Matrix<double> mat_dynamic_comma_initializer(3,3);
	
	mat_dynamic_comma_initializer << 
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0;

	for (int i = 0, nr = mat_dynamic.rows(); i < nr; ++i)
		for (int j = 0, nc = mat_dynamic.cols(); j < nc; ++j)
			REQUIRE ((mat_dynamic(i,j)) == (mat_dynamic_comma_initializer(i,j)));

	Matrix<double, 3, 3> mat_fixed_comma_initializer;
	mat_fixed_comma_initializer <<
		1.0, 2.0, 3.0,
		4.0, 5.0, 6.0,
		7.0, 8.0, 9.0;

	for (int i = 0, nr = mat_dynamic.rows(); i < nr; ++i)
		for (int j = 0, nc = mat_dynamic.cols(); j < nc; ++j)
			REQUIRE ((mat_dynamic(i,j)) == (mat_fixed_comma_initializer(i,j)));
}

TEST_CASE ("SimpleMathConstTranspose", "[SimpleMath]") {
	Matrix<double> mat_dynamic(3,3);

	mat_dynamic(0,0) = 1.0; mat_dynamic(0,1) = 2.0; mat_dynamic(0,2) = 3.0;
	mat_dynamic(1,0) = 4.0; mat_dynamic(1,1) = 5.0; mat_dynamic(1,2) = 6.0;
	mat_dynamic(2,0) = 7.0; mat_dynamic(2,1) = 8.0; mat_dynamic(2,2) = 9.0;

	const Block<Matrix<double>, double, 3, 1> const_block(&mat_dynamic, 0, 0);

	Matrix<double, 3, 1> vec1;
	vec1 = const_block;
	Matrix<double> res1;

	cout << "Result:" << endl;
	res1 = vec1.transpose() * vec1;
	cout << res1 << endl;
	cout << vec1.transpose() * vec1 << endl;
	cout << const_block.transpose() * const_block << endl;
	float result = const_block.transpose() * const_block;
	cout << result << endl;
}

TEST_CASE ("SimpleMathUnifiedFixedDynamic", "[SimpleMath]") {
	Matrix<double, 3, 3> fixed_mat33;
	Matrix<float> dynamic_mat33 (3, 3);

	dynamic_mat33 <<
			1.0, 2.0, 3.0,
			4.0, 5.0, 6.0,
			7.0, 8.0, 9.0;

	dynamic_mat33(2, 2) = 1000.0;

	fixed_mat33 = dynamic_mat33;

	cout << fixed_mat33(1,1) << endl;
	cout << "matrix " << sizeof(Matrix<double, 3, 3>) << endl;

	cout << "fixed: " << fixed_mat33 << endl;
	cout << "dyn:" << dynamic_mat33 << endl;

	cout << dynamic_mat33(0, 0) << ", " << dynamic_mat33(0, 1) << ", " << dynamic_mat33(0, 2) << endl
			<< dynamic_mat33(1, 0) << ", " << dynamic_mat33(1, 1) << ", " << dynamic_mat33(1, 2) << endl
			<< dynamic_mat33(2, 0) << ", " << dynamic_mat33(2, 1) << ", " << dynamic_mat33(2, 2) << endl;

	cout << sizeof(Matrix<float, 1, 1>) << endl;
    cout << sizeof(Matrix<float>) << endl;
}

TEST_CASE ("MultFixedAndDynamic", "[SimpleMath]") {
	Matrix<double, 3, 3> A;
	A <<
				1., 2., 3.,
			4., 4., 6.,
			7., 8., 9.;

	Matrix<double> B (Matrix<double>::Identity(3,3));
	Matrix<double, 3, 3> res = A * B;

	CHECK_ARRAY_CLOSE (A.data(), res.data(), 9, 1.0e-14);

    Matrix<double> res2 = B * A;
    CHECK_ARRAY_CLOSE (A.data(), res2.data(), 9, 1.0e-14);
}

TEST_CASE ("HouseholderQRSimple", "[SimpleMath]") {
	Matrix<double, 3, 3> A;
	A <<
			1., 2., 3.,
			4., 5., 6.,
			7., 8., 4.;

	Matrix<double, 3, 1> x;
	x[0] = 1.;
	x[1] = 2.;
	x[2] = 3.;

	Matrix<double, 3, 1> b  = A * x;

	HouseholderQR<Matrix<double, 3, 3> > qr = A.householderQr();
	Matrix<double, 3, 3> Q = qr.householderQ();
	Matrix<double, 3, 3> R = qr.matrixR();
    Matrix<double, 3, 1> x_qr = qr.solve(b);

	CHECK_ARRAY_CLOSE (x.data(), x_qr.data(), 3, 1.0e-14);
}


TEST_CASE ("ColPivHouseholderQRSimple", "[SimpleMath]") {
	Matrix<double, 3, 3> A;
	A <<
	  1., 2., 3.,
	  4., 5., 6.,
	  7., 8., 4.;

	Matrix<double, 3, 1> x;
	x[0] = 1.;
	x[1] = 2.;
	x[2] = 3.;

	Matrix<double, 3, 1> b  = A * x;

	ColPivHouseholderQR<Matrix<double, 3, 3> > qr = A.colPivHouseholderQr();
	Matrix<double, 3, 3> Q = qr.householderQ();
	Matrix<double, 3, 3> R = qr.matrixR();
	Matrix<double, 3, 3> P = qr.matrixP();

	Matrix<double, 3, 1> x_qr = qr.solve(b);

	CHECK_ARRAY_CLOSE (x.data(), x_qr.data(), 3, 1.0e-14);
}

TEST_CASE ("ColPivHouseholderQRSimpleDynamic", "[SimpleMath]") {
	Matrix<double> A (3, 3);

	A <<
	  1., 2., 3.,
			4., 5., 6.,
			7., 8., 4.;

	Matrix<double> x (3, 1);
	x[0] = 1.;
	x[1] = 2.;
	x[2] = 3.;

	Matrix<double, 3, 1> b  = A * x;

	ColPivHouseholderQR<Matrix<double> > qr = A.colPivHouseholderQr();
	Matrix<double, 3, 3> Q = qr.householderQ();
	Matrix<double, 3, 3> R = qr.matrixR();
	Matrix<double, 3, 3> P = qr.matrixP();

	Matrix<double> x_qr = qr.solve(b);

	CHECK_ARRAY_CLOSE (x.data(), x_qr.data(), 3, 1.0e-14);
}

TEST_CASE ("InverseDynamic", "[SimpleMath]") {
	Matrix<float> A (3, 3);

	A <<
	  1., 2., 3.,
	  4., 5., 6.,
	  7., 8., 4.;

	Matrix<float> Ainv = A.inverse();
	Matrix<float, 3, 3> identity = Matrix<float, 3, 3>::Identity();
	Matrix<float> A_times_Ainv = A * Ainv;

	CHECK_ARRAY_CLOSE (identity.data(), A_times_Ainv.data(), 9, 1.0e-5f);
}

TEST_CASE("ScalarMultTranspose", "[SimpleMath]") {
	Matrix<float, 3, 1> x;
	Matrix<float, 1, 3> y;
	x << 1.f, 2.f, 3.f;

	cout << 2.0f * x.transpose() << endl;
	cout << x.transpose() * x << endl;
}


TEST_CASE("SpatialMatrix_Multiplication", "[SimpleMath]") {
	typedef Matrix<double, 6, 6> Matrix66d;
 	Matrix66d X_1 (
      1.,  2.,  3.,  4.,  5.,  6.,
      11., 12., 13., 14., 15., 16.,
      21., 22., 23., 24., 25., 26.,
      31., 32., 33., 34., 35., 36.,
      41., 42., 43., 44., 45., 46.,
      51., 52., 53., 54., 55., 56.
      );

  Matrix66d X_2 (X_1);

  X_2 *= 2.0;

  Matrix66d correct_result (
      1442,    1484,    1526,    1568,    1610,    1652,
      4562,    4724,    4886,    5048,    5210,    5372,
      7682,    7964,    8246,    8528,    8810,    9092,
      10802,   11204,   11606,   12008,   12410,   12812,
      13922,   14444,   14966,   15488,   16010,   16532,
      17042,   17684,   18326,   18968,   19610,   20252
      );

  Matrix66d test_result = X_1 * X_2;

  CHECK_ARRAY_CLOSE (correct_result.data(), test_result.data(), 6 * 6, 1.0e-12);

  // check the *= operator:
  test_result = X_1;
  test_result *= X_2;

  CHECK_ARRAY_CLOSE (correct_result.data(), test_result.data(), 6 * 6, 1.0e-12);
}

TEST_CASE("NegativeTransposeCreatesCopy", "[SimpleMath]") {
	typedef Matrix<double, 3, 3> Matrix33d;

	Matrix33d A (
			1., 2., 3.,
			4., 5., 6.,
			7., 8., 9.
			);
	Matrix33d Aref (A);
	
	Matrix33d Aneg = -A;
	Matrix33d ATneg = -A.transpose();
	Matrix33d AnegT = Aneg.transpose();

	CHECK (Aref(0, 2) == -AnegT(2, 0));
	CHECK_ARRAY_CLOSE (Aref.data(), A.data(), 3 * 3, 1.0e-12);
	CHECK_ARRAY_CLOSE (ATneg.data(), AnegT.data(), 3 * 3, 1.0e-12);
}

TEST_CASE("SpecifyDynamicSizeInConstructor", "[SimpleMath]") {
	typedef Matrix<double, -1, -1> MatrixNd;

	MatrixNd matrix45 (4, static_cast<unsigned int>(5));

	CHECK (matrix45.rows() == 4);
	CHECK (matrix45.cols() == 5);
}

TEST_CASE("SpecifyDynamicVectorSizeInConstructor", "[SimpleMath]") {
	typedef Matrix<double, -1, -1> MatrixNd;

	MatrixNd matrix45 (static_cast<unsigned int>(5));

	CHECK (matrix45.rows() == 5);
	CHECK (matrix45.cols() == 1);
}

TEST_CASE("SpecifySizeInConstructor", "[SimpleMath]") {
	typedef Matrix<double, 4, 5> MatrixNd;

	MatrixNd matrix45 (static_cast<size_t>(4), 5);

	CHECK (matrix45.rows() == 4);
	CHECK (matrix45.cols() == 5);
}

TEST_CASE("Initialize2DVector", "[SimpleMath]") {
	typedef Matrix<double, 2, 1> Vector2d;

	Vector2d vec(4.1, 1.2);

	CHECK (vec[0] == 4.1);
	CHECK (vec[1] == 1.2);
}
