#include <iostream>
#include <iomanip>

#include "SimpleMath/SimpleMathDynamic.h"
#include "SimpleMath/SimpleMathFixed.h"

using namespace SimpleMath;
using namespace std;

typedef SimpleMath::Fixed::Matrix<double, 4, 4> Matrix44d;
typedef SimpleMath::Fixed::Matrix<double, 3, 3> Matrix33d;
typedef SimpleMath::Fixed::Matrix<double, 4, 1> Vector4d;
typedef SimpleMath::Fixed::Matrix<double, 3, 1> Vector3d;

typedef SimpleMath::Dynamic::Matrix<double> MatrixXd;
typedef SimpleMath::Dynamic::Matrix<double> VectorXd;

int main (int argc, char *argv[])
{
	cout << "= Fixed =" << endl;

	Matrix44d fix_matrix;
	Vector3d fix_vector;

	fix_matrix.random();
	fix_vector.random();

	cout << "vec             = " << fix_vector.transpose() << endl;
	cout << "matrix          = " << setw(15) << endl << fix_matrix << endl;
	//cout << "block<2,2>(2,1) = " << endl << fix_matrix.block<2,2>(2,1) << endl;

	MatrixXd dyn_matrix = MatrixXd::Zero(4,4);
	MatrixXd dyn_vector (4);

	dyn_matrix.random();
	dyn_vector.random();

	cout << endl << "= Dynamic =" << endl;

	cout << "vec             = " << dyn_vector.transpose() << endl;
	cout << "matrix          = " << endl << dyn_matrix << endl;
	//cout << "block<2,2>(2,1) = " << endl << dyn_matrix.block<2,2>(2,1) << endl;

}
