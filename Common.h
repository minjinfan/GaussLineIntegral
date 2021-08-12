#include <iostream>
#include <armadillo>
#include<vector>

using namespace std;
using namespace arma;

//typedef double real;

// X = a - b
inline arma::vec3 SubEq(const arma::vec3& X, const arma::vec3& Y)
{
	arma::vec3 res;
	res.zeros();
	for (int i = 0; i < 3; ++i) {
		res[i] = X[i] - Y[i];
	}
	return res;
}

inline void SubEq(arma::vec3& X, const arma::vec3& a, const arma::vec3& b)
{
	for (int i = 0; i < 3; ++i) {
		// for(int j = 0; j < 3; ++j){
		X[i] = a[i] - b[i];
		// }
	}
}

// Y = a * X + B
inline void faxpBy(arma::vec3& Y, double& a, arma::vec3& X, const arma::vec3& B)
{
	for (int i = 0; i < 3; ++i) {
		Y[i] = a * X[i] + B[i];
	}
}

// 点积
inline double fdot(arma::vec3 x, arma::vec3 y)
{
	double z = 0;
	for (size_t j = 0; j < 3; j++) {
		z += x[j] * y[j];
	}
	return z;
}

// 叉乘 X = a x b
inline void fcross(arma::vec3& X, const arma::vec3& a, const arma::vec3& b)
{
	X[0] = a[1] * b[2] - a[2] * b[1];
	X[1] = a[2] * b[0] - a[0] * b[2];
	X[2] = a[0] * b[1] - a[1] * b[0];
}

inline arma::vec3 fcross(const arma::vec3& a, const arma::vec3& b)
{
	arma::vec3 X;
	X.zeros();
	X[0] = a[1] * b[2] - a[2] * b[1];
	X[1] = a[2] * b[0] - a[0] * b[2];
	X[2] = a[0] * b[1] - a[1] * b[0];

	return X;
}

inline arma::vec3 faxpby(double& a, arma::vec3& X)
{
	arma::vec3 Y;
	Y.zeros();
	for (int i = 0; i < 3; ++i) {
		Y[i] = a * X[i];
	}
}

inline void faxpby(arma::vec3& Y, double& a, arma::vec3& X)
{
	for (int i = 0; i < 3; ++i) {
		Y[i] = a * X[i];
	}
}

