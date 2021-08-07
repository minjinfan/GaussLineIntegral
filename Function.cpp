#include"Function.h"

arma::vec3 Function::Getr0(const arma::vec3& A, const std::vector<arma::vec3>& VertexVec) {
	arma::vec3 temp1, temp2, n_, r0;

	r0.zeros();
	SubEq(temp1, VertexVec[0], VertexVec[1]);
	SubEq(temp2, VertexVec[0], VertexVec[2]);
	fcross(n_, temp1, temp2);
	n_ /= arma::norm(n_);
	double t = arma::dot(VertexVec[0] - A, n_) / pow(arma::norm(n_), 2);
	for (int i = 0; i < 3; ++i) {
		r0[i] = A[i] + n_[i] * t;
	}
	return r0;
}

arma::vec3 Function::OutNormal_edge(const arma::vec3& A, const std::vector<arma::vec3>& edge) 
{
	/*arma::vec3 temp1, temp2, n_, u_, r0;
	arma::vec3 n, u;
	r0.zeros();
	SubEq(temp1, edge[0], edge[1]);
	SubEq(temp2, edge[0], A);
	fcross(n_, temp1, temp2);
	n_ /= arma::norm(n_);
	fcross(u_, n_, temp1);
	u_ /= arma::norm(u_);
	double t = arma::dot(edge[0] - A, u_) / pow(arma::norm(u_), 2);
	for (int i = 0; i < 3; ++i) {
		r0[i] = A[i] + u_[i] * t;
	}*/

	arma::vec3 r0, u;
	r0 = GetDropFeet_edge(A, edge);

	SubEq(u, r0, A);
	u = u / arma::norm(u);
	 
	return u;

}

arma::vec3 Function::GetDropFeet_edge(const arma::vec3& A, const std::vector<arma::vec3>& edge)
{
	arma::vec3 temp1, temp2, n_, u_, r0;
	
	r0.zeros();
	SubEq(temp1, edge[0], edge[1]);
	SubEq(temp2, edge[0], A);
	fcross(n_, temp1, temp2);
	n_ /= arma::norm(n_);
	fcross(u_, n_, temp1);
	u_ /= arma::norm(u_);
	double t = arma::dot(edge[0] - A, u_) / pow(arma::norm(u_), 2);
	for (int i = 0; i < 3; ++i) {
		r0[i] = A[i] + u_[i] * t;
	}
	
	return r0;
}

double Function::GetArea(const std::vector<arma::vec3>& vertexVec)
{
	auto AreaOfTri = [](const arma::vec3& A, const arma::vec3& B, const arma::vec3& C)->double {
		double a = arma::norm(A - B);
		double b = arma::norm(A - C);
		double c = arma::norm(B - C);
		double p = (a + b + c) / 2;
		return sqrt(p * (p - a) * (p - b) * (p - c));
	};
	if (arma::norm(fcross(SubEq(vertexVec[0], vertexVec[1]), SubEq(vertexVec[2], vertexVec[1]))) == 0)
		return 0;
	return AreaOfTri(vertexVec[0], vertexVec[1], vertexVec[2]);
}