/**********************************************************************
 *
 * FaradayEDA
 * http://www.faradynamics.com
 *
 * Copyright (C) 2017-2019 Faraday Dynamics, Ltd.
 * All rights reserved.
 *
 * Author: zhaopeng@faradynamics.com
 *
 * Version v201901
 * Date: 2019-03-07 09:44:10
 *
 **********************************************************************/

#include "One_Over_R.h"
#include "GaussLegendreQuadrature.h"
using namespace std;

#include <functional>

	// void swap(FPoint &p1, FPoint &p2){
	//     FPoint p;
	//     p = p1;
	//     p1 = p2;
	//     p2 = p;
	// }
	double norm_area(arma::vec3& p0, arma::vec3& p1, arma::vec3& p2)
	{
		return arma::norm(arma::cross(p1 - p0, p2 - p1));
	}

	static int n0[17] = { 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 32, 64, 96, 100, 128, 256, 512 };
	static double* x0[17] = { x2, x4, x6, x8, x10, x12, x14, x16, x18, x20, x32, x64, x96, x100, x128, x256, x512 };
	static double* w0[17] = { w2, w4, w6, w8, w10, w12, w14, w16, w18, w20, w32, w64, w96, w100, w128, w256, w512 };
	void plot(double x1, double x2, function<double(double)> fun)
	{
		for (int i = 0; i < 200; i++) {
			double x = i * (x2 - x1) / 20.0 + x1;
			cout << x << "\t" << fun(x) << endl;
		}
	}
	double AdaptiveIntegral(double x1, double x2, function<double(double)> fun)
	{
		// plot(x1,x2,fun);
		const double tol = 2e-3;
		double I0 = 0;
		for (int i = 0; i < 17; i++) {
			double I = 0;
			int NGauss = n0[i] / 2;
			double* x = x0[i];
			double* w = w0[i];
			for (int k = 0; k < NGauss; k++) {
				double u = x[k] * 0.5 + 0.5;
				I += w[k] * fun(u * x2 + (1 - u) * x1);
				u = -x[k] * 0.5 + 0.5;
				I += w[k] * fun(u * x2 + (1 - u) * x1);
			}
			if (std::abs((I - I0) / I) < tol && i != 0) {
				// cout <<"Done:"<< i <<":"<<I << endl;
				return I * 0.5 * (x2 - x1);
			}
			else {
				// cout << i << endl;
				I0 = I;
			}
		}
		cout << "Warning::Out of Gaussian Quadrature Precision Range" << endl;
	}
	double newton(double x1, double x2, function<double(double)> fun)
	{
		return fun(x2) - fun(x1);
	}
	int test_adap()
	{
		function<double(double)> fun = [](double x) -> double { return log(sqrt(x * x + 1) - 1); };
		function<double(double)> fun_ = [](double x) -> double { return (x / sqrt(x * x + 1)) / (sqrt(x * x + 1) - 1); };
		double x1 = 0;

		double x2 = 3.2;
		cout << AdaptiveIntegral(x1, x2, fun) << endl;
		cout << -AdaptiveIntegral(x1, x2, fun_) + newton(x1, x2, [](double x) -> double {
			return x * log(sqrt(x * x + 1) - 1);
			}) << endl;
		cout << newton(x1, x2, [](double x) -> double { return x * log(x) - x; }) << endl;

		// [Debug]
		// if(std::abs(Di)<1e-4){
		//             Ri = sqrt(Ai * pow(eta[j], 2.0) + Bi * eta[j] + Ci);
		//    Xi = Di * eta[j] + Ei;
		//    Ui = Xi + (bi) / (2 * (ai_1));
		//     tmp1 = Ui * log(Ri + Xi);
		//    integral_result += pow(-1.0, i)*log(sqrt(ai*Ei*Ei+bi*Ei+ci)+Ei)/A;
		//    continue;
		// }
		return 0;
	}

	// [Debug]
	double compute_One_Over_R_part_General(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r)
	{
		double integral_result = 0, integral_tmp = 0;
		double Ai1, Bi1, Ci1, Di1, Ei1, Ai2, Bi2, Ci2, Di2, Ei2;
		double A, A2;
		double eta[] = { 0.0, 1.0 }, Newton_Leibniz_formula[] = { -1.0, 1.0 };
		arma::vec3 r_tmp = { 0, 0, 0 };
		arma::vec3 ri = { 0, 0, 0 };
		// cout<<arma::dot((r2 - r1), (r1 - r3))<<endl;
		if (std::abs(arma::dot((r2 - r1), (r1 - r3))) / (arma::norm(r2 - r1) * arma::norm(r3 - r1)) < 1e-1) {
			r_tmp = r1;
			r1 = r3;
			r3 = r_tmp;
		}
		else if (std::abs(arma::dot((r2 - r1), (r2 - r3))) / (arma::norm(r2 - r1) * arma::norm(r3 - r2)) < 1e-1) {
			r_tmp = r2;
			r2 = r3;
			r3 = r_tmp;
		}
		double x1 = 0;
		double x2 = 1.0;
		A = arma::norm(r2 - r1);
		A2 = A * A;

		auto& ri1 = r1;
		auto& ri2 = r2;
		Ai1 = pow(arma::norm(r3 - ri1), 2.0) / A2;
		Bi1 = -2 * arma::dot((r3 - ri1), (r - ri1)) / A2;
		Ci1 = pow(arma::norm(r - ri1), 2.0) / A2;
		Di1 = (arma::dot((r2 - r1), (r3 - ri1))) / A2;
		Ei1 = (-1) * arma::dot((r2 - r1), (r - ri1)) / A2;

		Ai2 = pow(arma::norm(r3 - ri2), 2.0) / A2;
		Bi2 = -2 * arma::dot((r3 - ri2), (r - ri2)) / A2;
		Ci2 = pow(arma::norm(r - ri2), 2.0) / A2;
		Di2 = (arma::dot((r2 - r1), (r3 - ri2))) / A2;
		Ei2 = (-1) * arma::dot((r2 - r1), (r - ri2)) / A2;

		integral_result += (-AdaptiveIntegral(x1, x2, [Ai1, Bi1, Ci1, Di1, Ei1, Ai2, Bi2, Ci2, Di2,
			Ei2](double x) -> double {
				double sq1 = sqrt(Ai1 * x * x + Bi1 * x + Ci1);
				double sq2 = sqrt(Ai2 * x * x + Bi2 * x + Ci2);
				double x01 = (2 * Di1 * Ei1 - Bi1) / (Ai1 - Di1 * Di1) / 2.0;
				double x02 = (2 * Di2 * Ei2 - Bi2) / (Ai2 - Di2 * Di2) / 2.0;
				return -(x - x01) * ((Ai1 * x + 0.5 * Bi1) / sq1 + Di1) / (sq1 + Di1 * x + Ei1) +
					(x - x02) * ((Ai2 * x + 0.5 * Bi2) / sq2 + Di2) / (sq2 + Di2 * x + Ei2);
			}) + newton(x1, x2, [Ai1, Bi1, Ci1, Di1, Ei1, Ai2, Bi2, Ci2, Di2, Ei2](double x) -> double {
				double x01 = (2 * Di1 * Ei1 - Bi1) / (Ai1 - Di1 * Di1) / 2.0;
				double x02 = (2 * Di2 * Ei2 - Bi2) / (Ai2 - Di2 * Di2) / 2.0;
				return -(x - x01) * log(sqrt(Ai1 * x * x + Bi1 * x + Ci1) + Di1 * x + Ei1) +
					(x - x02) * log(sqrt(Ai2 * x * x + Bi2 * x + Ci2) + Di2 * x + Ei2);
				}));

		return integral_result / A;
	}

	int compute_two_vectors_and_print()
	{
		arma::vec3 v1({ 1, 2, 3 });
		arma::vec3 v2({ 2, 3, 4 });
		arma::cross(v1, v2).print();
		return 0;
	}
	double compute_One_Over_R_part_simple(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3)
	{
		double integral_result = 0;
		auto singular_integral = [&integral_result](arma::vec3 p1, arma::vec3 p2, arma::vec3 p3) {
			arma::vec3 v_1_3 = { 0, 0, 0 }, v_3_2 = { 0, 0, 0 };
			double jacob = 0, a = 0, b = 0, c = 0;
			v_1_3 = p1 - p3;
			v_3_2 = p3 - p2;
			a = arma::dot(v_1_3, v_1_3);
			b = 2 * arma::dot(v_3_2, v_1_3);
			c = arma::dot(v_3_2, v_3_2);

			integral_result += (log(2 * sqrt(c * (a + b + c)) + 2 * c + b) - log(2 * sqrt(c * a) + b)) / sqrt(c);
		};
		singular_integral(r1, r2, r3);
		return integral_result;
	}

	double compute_One_Over_R_part(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r)
	{
		//	arma::vec3 r = r1;

		double integral_result = 0, integral_tmp = 0;
		double Ai, Bi, Ci, Di, Ei, Ri, Ui, Xi, ai, bi, ci, di;
		double A, A2;
		double eta[] = { 0.0, 1.0 }, Newton_Leibniz_formula[] = { -1.0, 1.0 };
		arma::vec3 r_tmp = { 0, 0, 0 };
		arma::vec3 ri = { 0, 0, 0 };
		// cout<<arma::dot((r2 - r1), (r1 - r3))<<endl;
		if (std::abs(arma::dot((r2 - r1), (r1 - r3))) / (arma::norm(r2 - r1) * arma::norm(r3 - r1)) < 1e-12) {
			r_tmp = r1;
			r1 = r3;
			r3 = r_tmp;
		}
		else if (std::abs(arma::dot((r2 - r1), (r2 - r3))) / (arma::norm(r2 - r1) * arma::norm(r3 - r2)) < 1e-12) {
			r_tmp = r2;
			r2 = r3;
			r3 = r_tmp;
		}
		A = arma::norm(r2 - r1);
		A2 = A * A;
		for (int i = 1; i <= 2; i++) {
			ri = (i == 1 ? r1 : r2);
			Ai = pow(arma::norm(r3 - ri), 2.0) / A2;
			Bi = -2 * arma::dot((r3 - ri), (r - ri)) / A2;
			Ci = pow(arma::norm(r - ri), 2.0) / A2;
			Di = (arma::dot((r2 - r1), (r3 - ri))) / A2;
			Ei = (-1) * arma::dot((r2 - r1), (r - ri)) / A2;
			ai = Ai / pow(Di, 2.0);
			bi = (Bi * Di - 2 * Ai * Ei) / pow(Di, 2.0);
			ci = (Ci * pow(Di, 2.0) - Bi * Di * Ei + Ai * pow(Ei, 2.0)) / pow(Di, 2.0);
			double di2 = (ci) / (ai - 1) - pow(bi, 2.0) / (4 * pow(ai - 1, 2.0));
			// if (std::abs(di2) > 1e-12)
			// 	di = sqrt(std::abs(di2));
			// else
			// 	di = 0;
			di = std::sqrt(std::abs(di2));
			// if (di2 > 1e-12) {
			// } else if (di2 < -1e-3) {
			// 	di = 0;
			// getchar();

			// cout << "Warning: In One Over R calculation, the di2 is negative!"
			// 	 << "di2:" << di2 << endl;
			// if (abs(acos(sqrt(abs(pow(Di, 2.0) / Ai)))) < 1e-3) {
			// 	// Msg::Error("The triangle has a very small angle:%f rad", acos(sqrt(abs(pow(Di, 2.0) / Ai))));
			// 	Msg::Warning("The triangle has a very small angle.");
			// 	cout.precision(10);
			// 	r1.raw_print(cout, "r1:");
			// 	r2.raw_print(cout, "r2:");
			// 	r3.raw_print(cout, "r3:");
			// 	r.raw_print(cout, "r:");
			// 	// getchar();
			// }
			// } else {
			// 	di = 0;
			// }

			integral_tmp = 0;
			for (int j = 0; j < 2; j++) {
				Ri = sqrt(Ai * pow(eta[j], 2.0) + Bi * eta[j] + Ci);
				Xi = Di * eta[j] + Ei;
				Ui = Xi + (bi) / (2 * (ai - 1));
				auto fun = [&Ai, &Bi, &Di, &ai, &bi, &ci, &di, &Ri, &Ui, &Xi](double eta) {
					double integral_result = 0, tmp1 = 0, tmp2 = 0, tmp3 = 0;
					if (Ri + Xi > 0)
						tmp1 = Ui * log(Ri + Xi);
					else
						tmp1 = 0;
					if (std::abs(2 * sqrt(ai) * Ri + 2 * ai * Xi + bi) > 0)
						tmp2 = -(bi / (2 * sqrt(ai) * (ai - 1)) * log(std::abs(2 * sqrt(ai) * Ri + 2 * ai * Xi + bi)));
					else
						tmp2 = 0;

					tmp3 = di * atan2((Ui), (di)) - di * atan2((2 * di * Ri * (ai - 1)), (bi * Xi + 2 * ci));
					if (std::isnan(tmp1) || std::abs(tmp1) > numeric_limits<double>::max()) {
						tmp1 = 0;
					}
					if (std::isnan(tmp2) || std::abs(tmp2) > numeric_limits<double>::max()) {
						tmp2 = 0;
					}
					if (std::isnan(tmp3) || std::abs(tmp3) > numeric_limits<double>::max()) {
						tmp3 = 0;
					}
					integral_result = tmp1 + tmp2 + tmp3;
					return integral_result;
				};
				integral_tmp += fun(eta[j]) * Newton_Leibniz_formula[j];
			}
			double temp0 = pow(-1.0, i) * integral_tmp / Di;
			if (!std::isnan(temp0) && std::abs(temp0) < numeric_limits<double>::max()) {
				integral_result += temp0;
			}
		}
		integral_result /= A;
		if (std::isnan(integral_result) || std::abs(integral_result) > numeric_limits<double>::max()) {
			integral_result = 0;
		}
		return integral_result;
	}
	double max(double l1, double l2)
	{
		return l1 > l2 ? l1 : l2;
	}
	double min(double l1, double l2)
	{
		return l1 < l2 ? l1 : l2;
	}
	arma::vec3& max(double l1, double l2, arma::vec3& r1, arma::vec3& r2)
	{
		return l1 > l2 ? r1 : r2;
	}
	arma::vec3& min(double l1, double l2, arma::vec3& r1, arma::vec3& r2)
	{
		return l1 < l2 ? r1 : r2;
	}
	bool unequal(arma::vec3& r1, arma::vec3& r2)
	{
		return !((r1(0) == r2(0)) && (r1(1) == r2(1)) && (r1(2) == r2(2)));
	}
	double compute_One_Over_R_part_Correct(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r) {

		double l3 = arma::norm(r2 - r1);
		double l2 = arma::norm(r3 - r1);
		double l1 = arma::norm(r2 - r3);
		double l_max = max(max(l1, l2), l3);
		double l_min = min(min(l1, l2), l3);
		arma::vec3& r12_max = max(l1, l2, r1, r2);
		arma::vec3& r123_max = max(max(l1, l2), l3, r12_max, r3);
		arma::vec3& r12_min = min(l1, l2, r1, r2);
		arma::vec3& r123_min = min(min(l1, l2), l3, r12_min, r3);

		// if(l_min/l_max <1e-3){
		// arma::vec3 &r_0=r123_min;
		// arma::vec3 &r_temp = (unequal(r2,r123_min)&&(unequal(r2,r123_max))?r2:r3);
		// arma::vec3 &r_1=(unequal(r1,r123_min)&&unequal(r1,r123_max))?r1:r_temp;
		// arma::vec3 &r_3=r123_max;
		// arma::vec3 r_2 = (r_1-r_0)*l_min/l_max+r_1;
	  //     double I_s = compute_One_Over_R_part(r_0, r_2, r_3,r) * norm_area(r_0, r_2, r_3)
		//   		- compute_One_Over_R_part(r_2, r_3, r_1,r) * norm_area(r_2, r_3, r_1) ;
	  //     return  I_s / ( norm_area(r_0, r_2, r_3)- norm_area(r_2, r_3, r_1));
		// }else{
		arma::vec3& r_0 = r123_min;
		arma::vec3& r_temp = (unequal(r2, r123_min) && (unequal(r2, r123_max)) ? r2 : r3);
		arma::vec3& r_1 = (unequal(r1, r123_min) && unequal(r1, r123_max)) ? r1 : r_temp;
		arma::vec3& r_3 = r123_max;
		return compute_One_Over_R_part(r_0, r_1, r_3, r);
		// }
	}

	double compute_One_Over_R_part_Correct(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r4, arma::vec3 r)
	{
		double I_s = compute_One_Over_R_part_Correct(r1, r2, r3, r) +
			compute_One_Over_R_part_Correct(r3, r4, r1, r);
		return I_s;
	}
	double compute_One_Over_R_part(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r4, arma::vec3 r)
	{
		double I_s = compute_One_Over_R_part(r1, r2, r3, r) +
			compute_One_Over_R_part(r3, r4, r1, r);
		return I_s;
	}
	double compute_A_B_C_D_E_F(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r, double& A, double& B, double& C,
		double& D, double& E, double& F)
	{
		A = std::pow(arma::norm(r2 - r1), 2);
		B = std::pow(arma::norm(r3 - r1), 2);
		C = -2. * arma::dot((r - r1), (r2 - r1));
		D = -2. * arma::dot((r - r1), (r3 - r1));
		E = 2. * arma::dot((r2 - r1), (r3 - r1));
		F = std::pow(arma::norm(r - r1), 2);
		return 0;
	}

	arma::vec3 compute_One_Over_R_vec(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r, arma::vec3 rio)
	{
		arma::vec3 r_tmp = { 0, 0, 0 };
		if (std::abs(arma::dot((r2 - r1), (r1 - r3))) / (arma::norm(r2 - r1) * arma::norm(r3 - r1)) < 1e-12) {
			r_tmp = r1;
			r1 = r3;
			r3 = r_tmp;
		}
		else if (std::abs(arma::dot((r2 - r1), (r2 - r3))) / (arma::norm(r2 - r1) * arma::norm(r3 - r2)) < 1e-12) {
			r_tmp = r2;
			r2 = r3;
			r3 = r_tmp;
		}

		double I0 = compute_One_Over_R_part_Correct(r1, r2, r3, r);
		double A = 0, B = 0, C = 0, D = 0, E = 0, F = 0;
		// double jacob = arma::norm(arma::cross((r2 - r1), (r3 - r1)));
		compute_A_B_C_D_E_F(r1, r2, r3, r, A, B, C, D, E, F);

		double J1 = ((2 * B - C + D - E) * (sqrt(B + D + F)) + (2 * A + C - D - E) * (sqrt(A + C + F))) / (4 * (A + B - E));
		if (std::abs(2 * sqrt((A + B - E) * (A + C + F)) - (2 * A + C - D - E)) > 1e-12) {
			double tmp1 = (4 * (A + C) * (B + D + F) + 4 * F * (B - C - E) - pow(C + D + E, 2)) /
				(8 * pow(A + B - E, 1.5)) *
				log(std::abs((2 * sqrt((A + B - E) * (B + D + F)) + (2 * B - C + D - E)) /
					(2 * sqrt((A + B - E) * (A + C + F)) - (2 * A + C - D - E))));
			if (!isnan(tmp1) && std::abs(tmp1) < numeric_limits<double>::max())
				J1 += tmp1;
		}

		double J2 = ((2 * B + D) * (sqrt(B + D + F)) - D * sqrt(F)) / (4 * B);

		if (std::abs(2 * sqrt(B * F) + D) > 1e-12) {
			double tmp1 = (4 * B * F - D * D) / (8 * B * sqrt(B)) *
				log(std::abs((2 * sqrt(B * (B + D + F)) + 2 * B + D) / (2 * sqrt(B * F) + D)));
			if (!isnan(tmp1) && std::abs(tmp1) < numeric_limits<double>::max())
				J2 += tmp1;
		}
		double J3 = ((2 * A + C - D - E) * (sqrt(A + C + F)) + (2 * B - C + D - E) * (sqrt(B + D + F))) / (4 * (A + B - E));

		if (std::abs(2 * sqrt((A + B - E) * (B + D + F)) - (2 * B - C + D - E)) > 1e-12) {
			double tmp1 = (4 * (A + C) * (B + D + F) + 4 * F * (B - C - E) - pow(C + D + E, 2)) /
				(8 * pow(A + B - E, 1.5)) *
				log(std::abs((2 * sqrt((A + B - E) * (A + C + F)) + (2 * A + C - D - E)) /
					(2 * sqrt((A + B - E) * (B + D + F)) - (2 * B - C + D - E))));
			if (!isnan(tmp1) && std::abs(tmp1) < numeric_limits<double>::max())
				J3 += tmp1;
		}

		double J4 = ((2 * A + C) * (sqrt(A + C + F)) - C * sqrt(F)) / (4 * A);

		if (std::abs(2 * sqrt(A * F) + C) > 1e-12) {
			double tmp1 = (4 * A * F - C * C) / (8 * A * sqrt(A)) *
				log(std::abs((2 * sqrt(A * (A + C + F)) + 2 * A + C) / (2 * sqrt(A * F) + C)));
			if (!isnan(tmp1) && std::abs(tmp1) < numeric_limits<double>::max())
				J4 += tmp1;
		}

		double Ip = (4 * B * (J1 - J2) - 2 * E * (J3 - J4) - (2 * B * C - E * D) * I0) / (4 * A * B - E * E);
		double Iq = (4 * A * (J3 - J4) - 2 * E * (J1 - J2) - (2 * A * D - E * C) * I0) / (4 * A * B - E * E);
		if (isnan(I0) || std::abs(I0) > numeric_limits<double>::max())
			I0 = 0;
		if (isnan(Ip) || std::abs(Ip) > numeric_limits<double>::max())
			Ip = 0;
		if (isnan(Iq) || std::abs(Iq) > numeric_limits<double>::max())
			Iq = 0;
		return ((r1 - rio) * I0 + (r2 - r1) * Ip + (r3 - r1) * Iq);
	}

	arma::vec3 compute_One_Over_R_vec(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r, arma::vec3 rio,
		arma::vec3 unit_b)
	{
		unit_b = 1.0 / arma::norm(unit_b - rio) * (unit_b - rio);
		arma::vec3 r_tmp = { 0, 0, 0 };
		if (std::abs(arma::dot((r2 - r1), (r1 - r3))) / (arma::norm(r2 - r1) * arma::norm(r3 - r1)) < 1e-12) {
			r_tmp = r1;
			r1 = r3;
			r3 = r_tmp;
		}
		else if (std::abs(arma::dot((r2 - r1), (r2 - r3))) / (arma::norm(r2 - r1) * arma::norm(r3 - r2)) < 1e-12) {
			r_tmp = r2;
			r2 = r3;
			r3 = r_tmp;
		}

		double I0 = compute_One_Over_R_part(r1, r2, r3, r);
		double A = 0, B = 0, C = 0, D = 0, E = 0, F = 0;
		// double jacob = arma::norm(arma::cross((r2 - r1), (r3 - r1)));
		compute_A_B_C_D_E_F(r1, r2, r3, r, A, B, C, D, E, F);

		double J1 = ((2 * B - C + D - E) * (sqrt(B + D + F)) + (2 * A + C - D - E) * (sqrt(A + C + F))) / (4 * (A + B - E));
		if (std::abs(2 * sqrt((A + B - E) * (A + C + F)) - (2 * A + C - D - E)) > 1e-12) {
			double tmp1 = (4 * (A + C) * (B + D + F) + 4 * F * (B - C - E) - pow(C + D + E, 2)) /
				(8 * pow(A + B - E, 1.5)) *
				log(std::abs((2 * sqrt((A + B - E) * (B + D + F)) + (2 * B - C + D - E)) /
					(2 * sqrt((A + B - E) * (A + C + F)) - (2 * A + C - D - E))));
			if (!isnan(tmp1) && std::abs(tmp1) < numeric_limits<double>::max())
				J1 += tmp1;
		}

		double J2 = ((2 * B + D) * (sqrt(B + D + F)) - D * sqrt(F)) / (4 * B);

		if (std::abs(2 * sqrt(B * F) + D) > 1e-12) {
			double tmp1 = (4 * B * F - D * D) / (8 * B * sqrt(B)) *
				log(std::abs((2 * sqrt(B * (B + D + F)) + 2 * B + D) / (2 * sqrt(B * F) + D)));
			if (!isnan(tmp1) && std::abs(tmp1) < numeric_limits<double>::max())
				J2 += tmp1;
		}
		double J3 = ((2 * A + C - D - E) * (sqrt(A + C + F)) + (2 * B - C + D - E) * (sqrt(B + D + F))) / (4 * (A + B - E));

		if (std::abs(2 * sqrt((A + B - E) * (B + D + F)) - (2 * B - C + D - E)) > 1e-12) {
			double tmp1 = (4 * (A + C) * (B + D + F) + 4 * F * (B - C - E) - pow(C + D + E, 2)) /
				(8 * pow(A + B - E, 1.5)) *
				log(std::abs((2 * sqrt((A + B - E) * (A + C + F)) + (2 * A + C - D - E)) /
					(2 * sqrt((A + B - E) * (B + D + F)) - (2 * B - C + D - E))));
			if (!isnan(tmp1) && std::abs(tmp1) < numeric_limits<double>::max())
				J3 += tmp1;
		}

		double J4 = ((2 * A + C) * (sqrt(A + C + F)) - C * sqrt(F)) / (4 * A);

		if (std::abs(2 * sqrt(A * F) + C) > 1e-12) {
			double tmp1 = (4 * A * F - C * C) / (8 * A * sqrt(A)) *
				log(std::abs((2 * sqrt(A * (A + C + F)) + 2 * A + C) / (2 * sqrt(A * F) + C)));
			if (!isnan(tmp1) && std::abs(tmp1) < numeric_limits<double>::max())
				J4 += tmp1;
		}

		double Ip = (4 * B * (J1 - J2) - 2 * E * (J3 - J4) - (2 * B * C - E * D) * I0) / (4 * A * B - E * E);
		double Iq = (4 * A * (J3 - J4) - 2 * E * (J1 - J2) - (2 * A * D - E * C) * I0) / (4 * A * B - E * E);
		if (isnan(I0) || std::abs(I0) > numeric_limits<double>::max())
			I0 = 0;
		if (isnan(Ip) || std::abs(Ip) > numeric_limits<double>::max())
			Ip = 0;
		if (isnan(Iq) || std::abs(Iq) > numeric_limits<double>::max())
			Iq = 0;
		return arma::dot((r1 - rio) * I0 + (r2 - r1) * Ip + (r3 - r1) * Iq, unit_b) * unit_b;
	}
	arma::vec3 compute_One_Over_R_vec(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r4, arma::vec3 r,
		arma::vec3 rio, arma::vec3 unit_b)
	{
		arma::vec3 I_s = compute_One_Over_R_vec(r1, r2, r3, r, rio, unit_b) +
			compute_One_Over_R_vec(r3, r4, r1, r, rio, unit_b);
		return I_s;
	}

	double compute_Ip(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r)
	{
		double A, B, C, D, E, F;

		compute_A_B_C_D_E_F(r1, r2, r3, r, A, B, C, D, E, F);

		int Ni, Nj;
		double x1, x2, dx1, dx2;
		double Ip;

		Ni = 20000; // give the grid num for x1
		Nj = 20000; // give the grid num for x2
		dx1 = 1.0 / Ni;
		dx2 = 1.0 / Nj;
		Ip = 0; // initialze the integral

		for (int i = 1; i < Ni; i++) {
			x1 = i * dx1;
			for (int j = 1; j < Nj; j++) {
				x2 = j * dx2;

				if (x2 < 1 - x1) {
					Ip += x2 / sqrt(A * x2 * x2 + B * x1 * x1 + C * x2 + D * x1 + E * x2 * x1 + F);
				}
			}
		}
		Ip *= dx1 * dx2;
		return Ip;
	}
	double compute_Iq(arma::vec3 r1, arma::vec3 r2, arma::vec3 r3, arma::vec3 r)
	{
		double A, B, C, D, E, F;

		compute_A_B_C_D_E_F(r1, r2, r3, r, A, B, C, D, E, F);

		int Ni, Nj;
		double x1, x2, dx1, dx2;
		double Iq;

		Ni = 20000; // give the grid num for x1
		Nj = 20000; // give the grid num for x2
		dx1 = 1.0 / Ni;
		dx2 = 1.0 / Nj;
		Iq = 0; // initialze the integral

		for (int i = 1; i < Ni; i++) {
			x1 = i * dx1;

			for (int j = 1; j < Nj; j++) {
				x2 = j * dx2;

				if (x2 < 1 - x1) {
					Iq += x1 / sqrt(A * x2 * x2 + B * x1 * x1 + C * x2 + D * x1 + E * x2 * x1 + F);
				}
			}
		}
		Iq *= dx1 * dx2;
		return Iq;
	}

	// int main()
	// {
	// 	cout << "welcome" << endl;
	// 	compute_two_vectors_and_print();

	// 	arma::vec3 v1(0.2, 0.2, 0);
	// 	arma::vec3 v2(1, 0, 0);
	// 	arma::vec3 v3(0, 1, 0);
	// 	arma::vec3 v(0.2, 0.2, 0);

	// 	double result0 = compute_One_Over_R_part(v1, v2, v3);
	// 	double result1 =compute_One_Over_R_part_simple(v1, v2, v3);

	// 	cout << result0 << "," << result1 << endl;
	// 	//arma::vec3 result = compute_One_Over_R_vec(v1, v2, v3, v);
	// 	arma::vec3 result2 = compute_One_Over_R_vec_part(v1, v2, v3);

	// 	// cout << "One_Over_R_vec:"
	// 	// 	 << "(";
	// 	// cout << result.x << ",";
	// 	// cout << result.y << ",";
	// 	// cout << result.z;
	// 	// cout << ")" << endl;
	// 	cout << "Ip:";
	// 	cout << compute_Ip(v1, v2, v3, v) << endl;
	// 	cout << "Iq:";
	// 	cout << compute_Iq(v1, v2, v3, v) << endl;
	// 	return 1;

	// }

