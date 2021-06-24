#include"main.h"
#include"GaussLineIntegral.h"

using namespace std;

int main(int argc, char** argv)
{
	arma::vec3 A = { 1,1,0 };
	arma::vec3 B = { 2,2,0 };
	arma::vec3 C = { 1,2,0 };
	vector<arma::vec3> VertexVec_s;
	VertexVec_s.push_back(C);
	VertexVec_s.push_back(A);
	VertexVec_s.push_back(B);

	GaussLineIntegral gs;

	size_t length = VertexVec_s.size();
	double res_total = 0.0;
	double tmp = 0.0;
	for (int i = 0; i < length; ++i) {
		double res = 0.0;
		arma::vec3 u;
		arma::vec3 A = VertexVec_s[i];
		vector<arma::vec3> edge;
		edge.push_back(VertexVec_s[(i + 1) % 3]);
		edge.push_back(VertexVec_s[(i + 2) % 3]);
		u = gs.GetNormal_edge(A, edge);

		std::vector<arma::vec3> list_gauss_point_s_edge;
		std::vector<double> list_gauss_weight_s_edge;
		gs.GenerateGaussPointEdge(edge, list_gauss_point_s_edge, list_gauss_weight_s_edge);

		int num_gauss_nodes_s_edge = list_gauss_point_s_edge.size();
		for (int j = 0; j < num_gauss_nodes_s_edge; ++j) {
			tmp = list_gauss_weight_s_edge[j];
			res += tmp;
		}
		cout << "Id:  " << i << "  edge" << "\t\t" << res << endl;
	}

	return 0;
}