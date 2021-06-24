#pragma once
#include <iostream>
#include <armadillo>
#include<vector>

#include"Common.h"

using namespace std;
using namespace arma;

typedef double real;

class GaussLineIntegral
{
public:
	GaussLineIntegral(){}
	~GaussLineIntegral(){}

	void GenerateGaussPointEdge(std::vector<arma::vec3>& edge, std::vector<arma::vec3>& list_gauss_point_s_edge, std::vector<double>& list_gauss_weight_s_edge);
	arma::vec3 GetNormal_edge(const arma::vec3& A, const std::vector<arma::vec3>& edge);
};

