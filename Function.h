#pragma once
#include <iostream>
#include <armadillo>
#include<vector>
//#include<complex.h>

#include"GaussLineIntegral.h"
#include"One_Over_R.h"
//#include"Common.h"

using namespace std;
class Function : public GaussLineIntegral
{
public:
	GaussLineIntegral gs;

	double GetJacobi_f = 3.56473125 * 1e-6;
	double GetJacobi_s = 1.42875 * 1e-6;

	arma::vec3 Vm = { 0.00242875,0.051,0.006 };
	arma::vec3 Vn = { 0.001,0.051,0.00101 };

	double l_div_2area_f;
	double l_div_2area_s;

public:
	Function(){}
	~Function(){}

	arma::vec3 Getr0(const arma::vec3& A, const std::vector<arma::vec3>& VertexVec);

	double Getd(const arma::vec3& a, const std::vector<arma::vec3>& vec_b);

	arma::vec3 OutNormal_edge(const arma::vec3& A, const std::vector<arma::vec3>& edge);

	arma::vec3 GetNormal_edge(const arma::vec3& A, const std::vector<arma::vec3>& edge);

	arma::vec3 GetDropFeet_edge(const arma::vec3& A, const std::vector<arma::vec3>& edge);

	void GetFactor(const arma::vec3& A, const std::vector<arma::vec3>& VertexVec);

	void InnerIntegralNonsingularFigure_R3_EM(double& result, const arma::vec3& r, const std::vector<arma::vec3>& VertexVec);
	
	double GetArea(const std::vector<arma::vec3>& vertexVec);

	bool IsInTriangle(const arma::vec3& A, const std::vector<arma::vec3>& vertexVec);
	bool IsOnEdge(const arma::vec3& A, const std::vector<arma::vec3>& vertexVec);
	bool IsOnEdge_edge(const arma::vec3& r0, const std::vector<arma::vec3>& vertexVec);

	void LineGaussIntegral(vector<arma::vec3> & VertexVec);

	void HighSingularity_EM(double& result,const std::vector<arma::vec3>& VertexVec_s, 
		const std::vector<arma::vec3>& VertexVec_f);

	void InnerIntegralSingularFigureFreespace_R3_EM(arma::vec3& one_Rth_Vec_1, arma::vec3& one_Rth_Vec_2,
		double& one_Rth, const vector<arma::vec3>& VertexVec_f, const arma::vec3& gauss_point_s);

	arma::vec3 GeneralIntegralOverROfTriVec(const std::vector<arma::vec3>& vertexes,const arma::vec3& gauss_point);

	double GeneralIntegralOverROfTri(const std::vector<arma::vec3>& vertexes, const arma::vec3& gauss_point);
};

