#include"main.h"


using namespace std;

int main(int argc, char** argv)
{
#ifdef Test
	double x_tmp = 1e-3;
	double y_tmp = 1e-3;
	arma::vec3 A = { 0,0,0 };
	arma::vec3 B = { x_tmp,0,0 };
	arma::vec3 C = { 0,y_tmp,0 };
	vector<arma::vec3> VertexVec;
	VertexVec.push_back(A);
	VertexVec.push_back(B);
	VertexVec.push_back(C);

	arma::vec3 p = { 1.2e-3, 1.2e-3, 0.5e-3 };
#else
	double x_tmp = 1e-3;
	double y_tmp = 1e-3;
	arma::vec3 AS = { 0.001,0.052,0.00101 };
	arma::vec3 BS = { 0.001,0.051,0.00101 };
	arma::vec3 CS = { 0.00242875,0.051,0.00101 };
	vector<arma::vec3> VertexVec_s;
	VertexVec_s.push_back(AS);
	VertexVec_s.push_back(BS);
	VertexVec_s.push_back(CS);

	arma::vec3 AF = { 0.001,0.051,0.006 };
	arma::vec3 BF = { 0.001,0.051,0.003505 };
	arma::vec3 CF = { 0.00242875,0.051,0.006 };
	vector<arma::vec3> VertexVec_f;
	VertexVec_f.push_back(AF);
	VertexVec_f.push_back(BF);
	VertexVec_f.push_back(CF);

	arma::vec3 p = { 1.2e-3, 1.2e-3, 0.5e-3 };
#endif 

	{ // 测试高阶奇异
		cout << endl;
		double result = 0.0;
		fc.HighSingularity_EM(result, VertexVec_s, VertexVec_f);
		cout << "result:  " << result << endl;
		cout << endl;
	}

	//{
	//	cout << endl;
	//	arma::vec3 fp = p;
	//	double d = fp[2];

	//	double x_max = x_tmp;
	//	double y_max = y_tmp;

	//	vector<double> x_list;
	//	vector<double> y_list;

	//	int step = 2000;
	//	for (int i = 1; i < step; i++) {
	//		x_list.push_back(x_max / step * i);
	//		y_list.push_back(y_max / step * i);
	//	}

	//	vector<arma::vec3> a_list;
	//	vector<arma::vec3> DiscritPoint_list;
	//	DiscritPoint_list.clear();
	//	a_list.clear();
	//	for (auto x : x_list) {
	//		for (auto y : y_list) {
	//			if (x < x_max && y < y_max - x) {
	//				arma::vec3 a = { x, y, 0 };
	//				DiscritPoint_list.push_back(a);
	//			}
	//			/*arma::vec3 a = { x, y, 0 };
	//			a_list.push_back(a);*/
	//		}
	//	}
	//	cout << DiscritPoint_list.size() << endl;
	//	cout << a_list.size() << endl;

	//	/*vector<arma::vec3> DiscritPoint_list;
	//	for (auto a : a_list) {
	//		if (a[0] < x_max && a[1] < y_max - a[0]) {
	//			DiscritPoint_list.push_back(a);
	//		}
	//	}
	//	cout << DiscritPoint_list.size() << endl;*/

	//	double ds = x_max * y_max / 2.0 / DiscritPoint_list.size();
	//	double res1 = 0.0;
	//	double res2 = 0.0;

	//	/*for (int i = DiscritPoint_list.size() - 1; i > DiscritPoint_list.size() - 11; i--) {
	//		cout << DiscritPoint_list[i] << endl;
	//	}*/
	//	for (arma::vec3 sp : DiscritPoint_list) {
	//		double R = arma::norm(fp - sp);
	//		double one_over_R3 = d / pow(R, 3);

	//		res1 += 1.0 * ds;
	//		res2 += one_over_R3 * ds;
	//	}
	//	cout << "res1:  " << res1 << endl;
	//	cout << "res2:  " << res2 << endl;

	//	//cout << "res:  " << res2 * d << endl;
	//}
	{ // 测试函数区
#ifdef Test
		double result = 0.0;
		fc.InnerIntegralNonsingularFigure_R3_EM(result, p, VertexVec);
		cout << "result:  " << result << endl;

		double testArea = fc.GetArea(VertexVec);
		cout << "testArea:  " << testArea << endl;
#endif
	}

	/*vector<arma::vec3> VertexVec_s;
	VertexVec_s.push_back(C);
	VertexVec_s.push_back(A);
	VertexVec_s.push_back(B);

	fc.LineGaussIntegral(VertexVec_s);*/

	/*vector<arma::vec3> edge;
	edge.push_back(VertexVec_s[1]);
	edge.push_back(VertexVec_s[0]);
	arma::vec3 u = gs.GetNormal_edge(B, edge);
	std::vector<arma::vec3> list_gauss_point_s_edge;
	std::vector<double> list_gauss_weight_s_edge;
	gs.GenerateGaussPointEdge(edge, list_gauss_point_s_edge, list_gauss_weight_s_edge);

	int num_gauss_nodes_s_edge = list_gauss_point_s_edge.size();
	double res = 0.0;
	for (int j = 0; j < num_gauss_nodes_s_edge; ++j) {
		tmp = list_gauss_weight_s_edge[j];
		res += tmp;
	}*/

	return 0;
}