#pragma once

#ifndef ONEOVERR_H_
#define ONEOVERR_H_
#include <armadillo>
#include <cmath>
#include <iostream>

	int compute_two_vectors_and_print();
	double compute_One_Over_R_part_simple(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3);
	double compute_One_Over_R_part(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r);
	double compute_One_Over_R(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r);
	arma::vec3 compute_One_Over_R_vec(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r4, arma::vec3 r, arma::vec3 rio, arma::vec3 unit_b);
	double compute_One_Over_R_part(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r4, arma::vec3 r);
	double compute_One_Over_R_part_Correct(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r);
	double compute_One_Over_R_part_Correct(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r4, arma::vec3 r);

	double compute_A_B_C_D_E_F(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r, double& A, double& B, double& C,
		double& D, double& E, double& F);

	arma::vec3 compute_One_Over_R_vec_part(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 rio);
	arma::vec3 compute_One_Over_R_vec(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r, arma::vec3 rio);

	double compute_Ip(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r);
	double compute_Iq(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r);

	// void swap(FPoint &p1, FPoint &p2);

#endif