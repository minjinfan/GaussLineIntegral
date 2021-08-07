#pragma once
#include <iostream>
#include <armadillo>
#include<vector>

#include"GaussLineIntegral.h"
//#include"Common.h"

class Function : public GaussLineIntegral
{
public:
	Function(){}
	~Function(){}

	arma::vec3 Getr0(const arma::vec3& A, const std::vector<arma::vec3>& VertexVec);

	arma::vec3 OutNormal_edge(const arma::vec3& A, const std::vector<arma::vec3>& edge);

	arma::vec3 GetDropFeet_edge(const arma::vec3& A, const std::vector<arma::vec3>& edge);

	double GetArea(const std::vector<arma::vec3>& vertexVec);
};

