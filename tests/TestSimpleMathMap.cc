#include <UnitTest++.h>
#include "SimpleMath/SimpleMathFixed.h"
#include "SimpleMath/SimpleMathDynamic.h"
#include "SimpleMath/SimpleMathMap.h"
#include <iostream>

using namespace std;
using namespace SimpleMath;

typedef SimpleMath::Fixed::Matrix<double, 3, 3> Matrix3d;
typedef SimpleMath::Fixed::Matrix<double, 3, 1> Vector3d;

typedef SimpleMath::Dynamic::Matrix<double> MatrixXd;
typedef SimpleMath::Dynamic::Matrix<double> VectorXd;

TEST (SimpleMathMapTestFixedSimple) {
	double data_3 [3] = {1.1, 1.2, 1.3};

	VectorXd vec = SimpleMath::Map<VectorXd> (data_3, 3, 1);

	CHECK_EQUAL (1.1, vec[0]);
	CHECK_EQUAL (1.2, vec[1]);
	CHECK_EQUAL (1.3, vec[2]);
}

TEST (SimpleMathMapTestFixedWrite) {
	double data_3 [3] = {1.1, 1.2, 1.3};

	VectorXd vec = SimpleMath::Map<VectorXd> (data_3, 3, 1);

	vec = vec * 2.;
	vec[2] = 999.;

	CHECK_EQUAL (2.2, data_3[0]);
	CHECK_EQUAL (2.4, data_3[1]);
	CHECK_EQUAL (999., data_3[2]);
}
