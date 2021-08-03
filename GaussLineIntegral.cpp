#include "GaussLineIntegral.h"


void GaussLineIntegral::GenerateGaussPointEdge(vector<arma::vec3>& edge, vector<arma::vec3>& list_gauss_point_s_edge, vector<double>& list_gauss_weight_s_edge)
{
	std::vector<double> X, W;
	int node_num = 12;
	if (node_num == 5) {
		X.push_back(0.0);
		X.push_back(0.5384693101);
		X.push_back(-0.5384693101);
		X.push_back(0.9061798459);
		X.push_back(-0.9061798459);

		W.push_back(128.0 / 225.0);
		W.push_back(0.4786286705);
		W.push_back(0.4786286705);
		W.push_back(0.2369268851);
		W.push_back(0.2369268851);
	}
	else if (node_num == 12) {
		X.push_back(-0.9815606342467190); W.push_back(0.0471753363865118);
		X.push_back(-0.9041172563704740); W.push_back(0.1069393259953180);
		X.push_back(-0.7699026741943040); W.push_back(0.1600783285433460);
		X.push_back(-0.5873179542866170); W.push_back(0.2031674267230650);
		X.push_back(-0.3678314989981800); W.push_back(0.2334925365383540);
		X.push_back(-0.1252334085114680); W.push_back(0.2491470458134020);
		X.push_back(0.1252334085114680); W.push_back(0.2491470458134020);
		X.push_back(0.3678314989981800); W.push_back(0.2334925365383540);
		X.push_back(0.5873179542866170); W.push_back(0.2031674267230650);
		X.push_back(0.7699026741943040); W.push_back(0.1600783285433460);
		X.push_back(0.9041172563704740); W.push_back(0.1069393259953180);
		X.push_back(0.9815606342467190); W.push_back(0.0471753363865118);
	}


	int X_size = X.size();
	arma::vec3 L_vec = SubEq(edge[1], edge[0]);
	double L_length = arma::norm(L_vec);
	arma::vec3 L_vec_n = L_vec / L_length;

	double common_term = L_length / 2.0;
	std::vector<double> new_x, new_w;
	new_x.resize(X_size);
	new_w.resize(X_size);
	for (int i = 0; i < X.size(); i++) {
		new_x[i] = common_term * X[i] + common_term;
		new_w[i] = common_term * W[i];
	}
	std::vector<arma::vec3> r_vec;
	for (int i = 0; i < X.size(); i++) {
		arma::vec3 r;
		faxpBy(r, new_x[i], L_vec_n, edge[0]);
		r_vec.push_back(r);
	}

	list_gauss_point_s_edge = r_vec;
	list_gauss_weight_s_edge = new_w;
}

arma::vec3 GaussLineIntegral::GetNormal_edge(const arma::vec3& A, const std::vector<arma::vec3>& edge) {
	arma::vec3 temp1, temp2, n_, u_, r0;
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
	}
	SubEq(u, r0, A);
	u = u / arma::norm(u);

	return u;

}