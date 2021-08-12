#include"Function.h"

void Function::LineGaussIntegral(vector<arma::vec3>& VertexVec)
{
	size_t length = VertexVec.size();
	double res_total = 0.0;
	double tmp = 0.0;
	for (int i = 0; i < length; ++i) {
		double res = 0.0;
		arma::vec3 u;
		arma::vec3 A = VertexVec[i];
		vector<arma::vec3> edge;
		edge.push_back(VertexVec[(i + 1) % 3]);
		edge.push_back(VertexVec[(i + 2) % 3]);
		u = GetNormal_edge(A, edge);

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
}

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

double Function::Getd(const arma::vec3& a, const std::vector<arma::vec3>& vec_b)
{
	double distance = 0.0;
	arma::vec3 temp1, temp2, n_;
	SubEq(temp1, vec_b[1], vec_b[0]);
	SubEq(temp2, vec_b[2], vec_b[0]);
	fcross(n_, temp1, temp2);

	distance = abs(n_[0] * a[0] + n_[1] * a[1] + n_[2] * a[2] - (n_[0] * vec_b[0][0] + n_[1] * vec_b[0][1] + n_[2] * vec_b[0][2])) /
		sqrt(pow(n_[0], 2) + pow(n_[1], 2) + pow(n_[2], 2));

	return distance;
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

arma::vec3 Function::GetNormal_edge(const arma::vec3& A, const std::vector<arma::vec3>& edge)
{
	
	arma::vec3 r0, u;
	r0 = GetDropFeet_edge(A, edge);

	SubEq(u, r0, A);
	u = u / arma::norm(u);

	return u;

}

arma::vec3 Function::GetDropFeet_edge(const arma::vec3& A, const std::vector<arma::vec3>& edge)
{
	/*arma::vec3 temp1, temp2, n_, u_, r0;
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
	return r0;*/


	arma::vec3 temp1, temp2, O;
	SubEq(temp1, A, edge[0]);
	SubEq(temp2, edge[0], edge[1]);

	double k = arma::dot(temp1, temp2) / pow(arma::norm(temp2), 2);
	for (int i = 0; i < 3; ++i) {
		O[i] = edge[0][i] + temp2[i] * k;
	}
	return O;
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
	arma::vec3 tmp1 = SubEq(vertexVec[0], vertexVec[1]);
	arma::vec3 tmp2 = SubEq(vertexVec[2], vertexVec[1]);
	arma::vec3 tmp3 = fcross(tmp1, tmp2);
	if (arma::norm(tmp3) == 0)
		return 0;

	double res = AreaOfTri(vertexVec[0], vertexVec[1], vertexVec[2]);
	return res;
}

bool Function::IsOnEdge(const arma::vec3& P, const std::vector<arma::vec3>& vertexVec)
{
	/*size_t length = vertexVec.size();
	double area = 0.0;

	for (int i = 0; i < length; ++i) {
		vector<arma::vec3> tri;
		tri.push_back(vertexVec[(i + 1) % 3]);
		tri.push_back(vertexVec[(i + 2) % 3]);

		tri.push_back(A);

		double tmp = GetArea(tri);
		if (tmp == 0) {
			return true;
			break;
		}
		area += tmp;
	}
	return false;*/

	bool Res = true;

	size_t length = vertexVec.size();
	for (int i = 0; i < length; ++i) {
		bool flag;
		vector<arma::vec3> tri;
		arma::vec3 A = vertexVec[(i + 1) % 3];
		arma::vec3 B = vertexVec[(i + 2) % 3];
		arma::vec3 C = vertexVec[i];

		arma::vec3 AB = SubEq(B, A);
		arma::vec3 AP = SubEq(P, A);
		arma::vec3 AC = SubEq(C, A);

		arma::vec3 tmp1 = fcross(AB, AP);
		arma::vec3 tmp2 = fcross(AB, AC);
		double res = fdot(tmp1, tmp2);

		if (res == 0) {
			flag = true;
			break;
		}
		else {
			flag = false;
		}

		Res = (Res && flag);
	}

	return Res;
}

bool Function::IsOnEdge_edge(const arma::vec3& r0, const std::vector<arma::vec3>& edge)
{
	bool Res = true;

	bool flag;
	vector<arma::vec3> tri;
	arma::vec3 A = edge[0];
	arma::vec3 B = edge[1];
	
	arma::vec3 tmp1 = SubEq(B, A);
	arma::vec3 tmp2 = SubEq(r0, A);
	double res = arma::norm(fcross(tmp1, tmp2));


	if (res == 0) {
		Res = true;
	}
	else {
		Res = false;
	}

	return Res;
}

bool Function::IsInTriangle(const arma::vec3& P, const std::vector<arma::vec3>& vertexVec)
{
	//double area = GetArea(vertexVec);

	//size_t length = vertexVec.size();
	//double areaAdd = 0.0;

	//for (int i = 0; i < length; ++i) {
	//	vector<arma::vec3> tri;
	//	tri.push_back(vertexVec[(i + 1) % 3]);
	//	tri.push_back(vertexVec[(i + 2) % 3]);
	//	tri.push_back(A);

	//	//double tmp = GetArea(tri);
	//	areaAdd += GetArea(tri);
	//}
	//if (area == areaAdd)
	//	return true;
	//return false;

	bool Res = true;

	size_t length = vertexVec.size();
	for (int i = 0; i < length; ++i) {
		bool flag;
		vector<arma::vec3> tri;
		arma::vec3 A = vertexVec[(i + 1) % 3];
		arma::vec3 B = vertexVec[(i + 2) % 3];
		arma::vec3 C = vertexVec[i];

		arma::vec3 AB = SubEq(B, A);
		arma::vec3 AP = SubEq(P, A);
		arma::vec3 AC = SubEq(C, A);

		arma::vec3 tmp1 = fcross(AB, AP);
		arma::vec3 tmp2 = fcross(AB, AC);
		double res = fdot(tmp1, tmp2);

		if (res < 0) {
			flag = false;
		}
		else {
			flag = true;
		}

		Res = (Res && flag);
	}
	return Res;

	// 以AB为直线 ， 判断 P、C 是否在AB的同一侧
	/*arma::vec3 AB = SubEq(vertexVec[1], vertexVec[0]);
	arma::vec3 AP = SubEq(P, vertexVec[0]);
	arma::vec3 AC = SubEq(vertexVec[2], vertexVec[0]);

	arma::vec3 tmp1 = fcross(AB, AP);
	arma::vec3 tmp2 = fcross(AB, AC);
	double res = fdot(tmp1, tmp2);

	if (res < 0) {
		return false;
	}
	else {
		return true;
	}*/

}

void Function::GetFactor(const arma::vec3& r, const std::vector<arma::vec3>& VertexVec)
{
	size_t length = VertexVec.size();

	arma::vec3 r0;
	r0 = Getr0(r, VertexVec);
	//cout << "r0:  " << r0 << endl;

	double d = Getd(r, VertexVec);
	double d2 = arma::norm(SubEq(r0, r));
	/*cout << "d: " << d << endl;
	cout << "d2: " << d2 << endl;*/
	bool IsInTri = IsInTriangle(r0, VertexVec);
	vector<arma::vec3> U;
	vector<arma::vec3> Pi;
	vector<arma::vec2> Li;
	vector<arma::vec2> Ri;
	vector<double> Ro_pow;
	U.resize(length);
	Pi.resize(length);
	Li.resize(length);
	Ri.resize(length);
	Ro_pow.resize(length);

	for (int i = 0; i < length; ++i) {
		bool flag;
		vector<arma::vec3> edge;
		arma::vec3 A = VertexVec[i];
		arma::vec3 B = VertexVec[(i + 1) % 3];
		arma::vec3 C = VertexVec[(i + 2) % 3];

		edge.push_back(B);
		edge.push_back(C);
		
		// Caculate U
		arma::vec3 u;
		arma::vec3 A_O = GetDropFeet_edge(A, edge);
		//cout << "A_O:  " << A_O << endl << endl;
		SubEq(u, A_O, A);
		u = u / arma::norm(u);
		U[i] = u;

		bool IsOnEdge = IsOnEdge_edge(r0, edge);
		arma::vec3 r0_O;
		if (IsOnEdge) {
			r0_O = r0;
		}
		else {
			r0_O = GetDropFeet_edge(r0, edge);
		}
		//cout << "r0_O:  " << r0_O << endl;

		// Caculate Pi
		arma::vec3 pi = SubEq(r0_O, r0);
		Pi[i] = SubEq(r0_O, r0);

		// Caculate Li
		bool IsInTri_O = IsInTriangle(r0_O, VertexVec);
		arma::vec2 li;
		if (IsInTri_O) {
			li[0] = arma::norm(SubEq(r0_O, B)) * (-1);
			li[1] = arma::norm(SubEq(r0_O, C));

			/*cout << SubEq(r0_O, B) << endl;
			cout << arma::norm(SubEq(r0_O, B)) << endl;
			cout << SubEq(r0_O, C) << endl;
			cout << arma::norm(SubEq(r0_O, C)) << endl;*/
		}
		else {
			double first = arma::norm(SubEq(r0_O, B));
			double second = arma::norm(SubEq(r0_O, C));
			if (first < second) {
				li[0] = first;
				li[1] - second;
			}
			else {
				li[0] = -1 * first;
				li[1] = -1 * second;
			}
		}
		Li[i] = li;

		// Caculate Ro 、Ri
		double Ro_pow_tmp = pow(arma::norm(pi), 2) + pow(d, 2);
		Ro_pow[i] = Ro_pow_tmp;
		arma::vec2 Ri_tmp;
		Ri_tmp[0] = sqrt(Ro_pow_tmp + pow(li[0], 2));
		Ri_tmp[1] = sqrt(Ro_pow_tmp + pow(li[1], 2));
		Ri[i] = Ri_tmp;

	}
	
	/*for (auto u : U) {
		cout << "U: " << u << endl;
	}
	cout << endl;
	for (auto p : Pi) {
		cout << "Pi:  " << p << endl;
	}
	cout << endl;
	for (auto l : Li) {
		cout << "Li:  " << l << endl;
	}*/

	cout << endl;
	double result = 0.0;
	for (int i = 0; i < length; ++i) {
		double pi = arma::norm(Pi[i]);
		double li_right = Li[i][1];
		double Ri_right = Ri[i][1];
		double res1 = atan2(pi * li_right, Ro_pow[i] + d * Ri_right);
		//cout << "atan2 res1:  " << res1 << endl;

		double li_left = Li[i][0];
		double Ri_left = Ri[i][0];
		double res2 = atan2(pi * li_left, Ro_pow[i] + d * Ri_left);
		//cout << "atan2 res2:  " << res2 << endl;

		double res = res1 - res2;
		//cout << "atan2 res:  " << res << endl;

		double result_tmp = fdot(U[i], Pi[i]) * res;
		cout << "result_tmp:  " << result_tmp << endl;
		result += result_tmp;
	}
	
	cout << endl;
}

void Function::InnerIntegralNonsingularFigure_R3_EM(double &result,const arma::vec3& r, const std::vector<arma::vec3>& VertexVec)
{
	size_t length = VertexVec.size();

	arma::vec3 r0;
	r0 = Getr0(r, VertexVec);
	cout << "r0:  " << r0 << endl;

	double d = Getd(r, VertexVec);
	double d2 = arma::norm(SubEq(r0, r));
	cout << "d: " << d << endl;
	//cout << "d2: " << d2 << endl;
	bool IsInTri = IsInTriangle(r0, VertexVec);
	vector<arma::vec3> U;
	vector<arma::vec3> Pi;
	vector<arma::vec2> Li;
	vector<arma::vec2> Ri;
	vector<double> Ro_pow;
	U.resize(length);
	Pi.resize(length);
	Li.resize(length);
	Ri.resize(length);
	Ro_pow.resize(length);

	for (int i = 0; i < length; ++i) {
		bool flag;
		vector<arma::vec3> edge;
		arma::vec3 A = VertexVec[i];
		arma::vec3 B = VertexVec[(i + 1) % 3];
		arma::vec3 C = VertexVec[(i + 2) % 3];

		edge.push_back(B);
		edge.push_back(C);

		// Caculate U
		arma::vec3 u;
		arma::vec3 A_O = GetDropFeet_edge(A, edge);
		//cout << "A_O:  " << A_O << endl << endl;
		SubEq(u, A_O, A);
		u = u / arma::norm(u);
		U[i] = u;

		bool IsOnEdge = IsOnEdge_edge(r0, edge);
		arma::vec3 r0_O;
		if (IsOnEdge) {
			r0_O = r0;
		}
		else {
			r0_O = GetDropFeet_edge(r0, edge);
		}
		cout << "ID " << i+1 << " Feet:  " << r0_O << endl;

		// Caculate Pi
		arma::vec3 pi = SubEq(r0_O, r0);
		Pi[i] = SubEq(r0_O, r0);

		// Caculate Li
		bool IsInTri_O = IsInTriangle(r0_O, VertexVec);
		arma::vec2 li;
		li.zeros();
		if (IsInTri_O) {
			li[0] = arma::norm(SubEq(r0_O, B)) * (-1);
			li[1] = arma::norm(SubEq(r0_O, C));

			/*cout << SubEq(r0_O, B) << endl;
			cout << arma::norm(SubEq(r0_O, B)) << endl;
			cout << SubEq(r0_O, C) << endl;
			cout << arma::norm(SubEq(r0_O, C)) << endl;*/
		}
		else {
			double first = arma::norm(SubEq(r0_O, B));
			double second = arma::norm(SubEq(r0_O, C));
			if (first < second) {
				li[0] = first;
				li[1] = second;
			}
			else {
				li[0] = -1 * first;
				li[1] = -1 * second;
			}
		}
		Li[i] = li;

		// Caculate Ro 、Ri
		double Ro_pow_tmp = pow(arma::norm(pi), 2) + pow(d, 2);
		Ro_pow[i] = Ro_pow_tmp;
		arma::vec2 Ri_tmp;
		Ri_tmp[0] = sqrt(Ro_pow_tmp + pow(li[0], 2));
		Ri_tmp[1] = sqrt(Ro_pow_tmp + pow(li[1], 2));
		Ri[i] = Ri_tmp;

	}

	for (auto u : U) {
		cout << "U: " << u << endl;
	}
	cout << endl;
	for (auto p : Pi) {
		cout << "Pi:  " << p << endl;
	}
	cout << endl;
	for (auto l : Li) {
		cout << "Li:  " << l << endl;
	}
	for (auto Ro2 : Ro_pow) {
		cout << "R02:  " << Ro2 << endl;
	}
	for (auto ri : Ri) {
		cout << "Ri:  " << ri << endl;
	}

	cout << endl;

	double Res = 0.0;
	for (int i = 0; i < length; ++i) {
		double pi = arma::norm(Pi[i]);
		double li_right = Li[i][1];
		double Ri_right = Ri[i][1];
		double res1 = atan2(pi * li_right, Ro_pow[i] + d * Ri_right);
		cout << "atan2 res1:  " << res1 << endl;

		double li_left = Li[i][0];
		double Ri_left = Ri[i][0];
		double res2 = atan2(pi * li_left, Ro_pow[i] + d * Ri_left);
		cout << "atan2 res2:  " << res2 << endl;

		double res = res1 - res2;
		cout << "atan2 res:  " << res << endl;

		arma::vec3 ppp;
		if (arma::norm(Pi[i]) != 0) {
			ppp = Pi[i] / arma::norm(Pi[i]);
		}
		else {
			ppp = { 0,0,0 };
		}
		double result_tmp = fdot(U[i], ppp) * res;
		cout << "result_tmp:  " << result_tmp << endl;
		Res += result_tmp;
	}
	result = Res ;
	cout << "Res:  " << Res << endl;

	cout << endl;
}

void Function::HighSingularity_EM(double& result,
	const std::vector<arma::vec3>& VertexVec_s, const std::vector<arma::vec3>& VertexVec_f)
{
	double NearSingular1_Yes = 0.0;
	size_t length = VertexVec_s.size();

	l_div_2area_f = 1.0 / 2.0 / GetArea(VertexVec_f);
	l_div_2area_s = 1.0 / 2.0 / GetArea(VertexVec_s);

	double sf = 1.0 / GetJacobi_f;
	double ss = 1.0 / GetJacobi_s;

	for (int i = 0; i < length; i++) {
		arma::vec3 u;
		arma::vec3 A = VertexVec_s[i];
		vector<arma::vec3> edge;
		edge.push_back(VertexVec_s[(i + 1) % 3]);
		edge.push_back(VertexVec_s[(i + 2) % 3]);
		u = GetNormal_edge(A, edge);

		double res_edge = 0.0;
		double integral_result_tmp1 = 0.0;
		double integral_result_tmp2 = 0.0;
		double integral_result_tmp3 = 0.0;

		std::vector<arma::vec3> list_gauss_point_s_edge;
		std::vector<double> list_gauss_weight_s_edge;
		gs.GenerateGaussPointEdge(edge, list_gauss_point_s_edge, list_gauss_weight_s_edge);
		// GenerateGaussPointEdge_new(edge, list_gauss_point_s_edge, list_gauss_weight_s_edge);
		int num_gauss_nodes_s_edge = list_gauss_point_s_edge.size();
		for (int j = 0; j < num_gauss_nodes_s_edge; j++) {

			arma::vec3 one_Rth_Vec_1;
			arma::vec3 one_Rth_Vec_2;
			double one_Rth = 0.0;

			arma::vec3 tmp1;
			arma::vec3 tmp2;
			arma::vec3 tmp3;

			InnerIntegralSingularFigureFreespace_R3_EM(one_Rth_Vec_1, one_Rth_Vec_2, one_Rth, VertexVec_f, list_gauss_point_s_edge[j]);

			fcross(tmp1, one_Rth_Vec_1, u);
			fcross(tmp2, one_Rth_Vec_2, u);
			faxpby(tmp3, one_Rth, u);

			integral_result_tmp1 = fdot(SubEq(Vn, Vm), tmp1);
			integral_result_tmp2 = fdot(SubEq(Vn, Vm), tmp2);
			integral_result_tmp3 = fdot(fcross(Vm, Vn), tmp3);

			NearSingular1_Yes += (integral_result_tmp1 + integral_result_tmp2 + integral_result_tmp3) * list_gauss_weight_s_edge[j];

			double res_tmp = (integral_result_tmp1 + integral_result_tmp2 + integral_result_tmp3) * list_gauss_weight_s_edge[j];
			cout << "(  " << i + 1 << " ,  " << j+1 << "  )  "<< res_tmp << endl;
			res_edge += res_tmp;
		}
		cout << "res_edge:  " << res_edge << endl;
	}
	//NearSingular1_Yes *= l_div_2area_f * l_div_2area_s;

	result = NearSingular1_Yes;

	cout << "result:  " << result << endl;
}

void Function::InnerIntegralSingularFigureFreespace_R3_EM(arma::vec3& one_Rth_Vec_1, arma::vec3& one_Rth_Vec_2,
									double& one_Rth, const vector<arma::vec3> &VertexVec_f, const arma::vec3& gauss_point_s)
{
	//arma::cx_vec3 one_Rth_Vec_tmp2, one_Rth_Vec_tmp3;
	arma::vec3 one_Rth_Vec_tmp1, one_Rth_Vec_tmp2, r0, SubVmr;

	one_Rth_Vec_tmp1 = GeneralIntegralOverROfTriVec(VertexVec_f, gauss_point_s) * (1.0);
	one_Rth = GeneralIntegralOverROfTri(VertexVec_f, gauss_point_s);

	
	faxpby(one_Rth_Vec_tmp2, one_Rth, Vm);

	one_Rth_Vec_1 = one_Rth_Vec_tmp1;
	one_Rth_Vec_2 = one_Rth_Vec_tmp2;

	one_Rth_Vec_1 *= GetJacobi_f;
	one_Rth_Vec_2 *= GetJacobi_f;
	one_Rth *= GetJacobi_f;
}

arma::vec3 Function::GeneralIntegralOverROfTriVec(const std::vector<arma::vec3>& vertexes,
	const arma::vec3& gauss_point)
{
	arma::vec3 integral({ 0, 0, 0 });

	if (vertexes.size() == 3) {
		integral = integral + compute_One_Over_R_vec(vertexes[0], vertexes[1], vertexes[2], gauss_point, vertexes[0]);
	}
	else {
		integral = integral + compute_One_Over_R_vec(vertexes[0], vertexes[1], vertexes[2], vertexes[3], gauss_point,
			vertexes[0], vertexes[1]) * 8.0;
	}
	return integral;
}

double Function::GeneralIntegralOverROfTri(const std::vector<arma::vec3>& vertexes, const arma::vec3& gauss_point)
{
	if (vertexes.size() == 3) {
		return compute_One_Over_R_part_Correct(vertexes[0], vertexes[1], vertexes[2], gauss_point);
	}
	else {
		return compute_One_Over_R_part_Correct(vertexes[0], vertexes[1], vertexes[2], vertexes[3], gauss_point) * 4.0;
	}
}