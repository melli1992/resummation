#include <complex>
#include <vector>
#ifndef ttHsoftanom_H
#define ttHsoftanom_H

double beta34(double s34);
std::complex<double> log_beta34(double s34);
std::complex<double> Tij(std::complex<double> tij, std::complex<double> s);
std::complex<double> Omega3(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> Lambda3(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);

std::complex<double> lambda_qq_11(double s34);
std::complex<double> lambda_qq_12(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_qq_21(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_qq_22(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, double s34, std::complex<double> s);
std::complex<double> lambda_gg_11(double s34);
std::complex<double> lambda_gg_12();
std::complex<double> lambda_gg_13(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_gg_21();
std::complex<double> lambda_gg_22(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s, double s34);
std::complex<double> lambda_gg_23(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_gg_31(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_gg_32(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_gg_33(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s, double s34);

std::vector<std::complex<double>> eigenvalues_qq(std::complex<double> G11, std::complex<double> G12, std::complex<double> G21, std::complex<double> G22);
std::vector<std::complex<double>> eigenvalues_gg(std::complex<double> G11, std::complex<double> G13, std::complex<double> G22, std::complex<double> G23,  std::complex<double> G31, std::complex<double> G32, std::complex<double> G33);
std::vector<std::vector<std::complex<double>>> eigenvectors_qq(std::vector<std::complex<double>> lambda, std::complex<double> G11, std::complex<double> G12, std::complex<double> G21, std::complex<double> G22);
std::vector<std::vector<std::complex<double>>> eigenvectors_gg(std::vector<std::complex<double>> lambda, std::complex<double> G11, std::complex<double> G13, std::complex<double> G22, std::complex<double> G23,  std::complex<double> G31, std::complex<double> G32, std::complex<double> G33);
std::vector<std::vector<std::complex<double>>> R_qq(std::vector<std::vector<std::complex<double>>> ev);
std::vector<std::vector<std::complex<double>>> R_gg(std::vector<std::vector<std::complex<double>>> ev);
std::vector<std::vector<std::complex<double>>> Rinv_qq(std::vector<std::vector<std::complex<double>>> ev);
std::vector<std::vector<std::complex<double>>> Rinv_gg(std::vector<std::vector<std::complex<double>>> ev);
std::vector<std::vector<std::complex<double>>> Rdag_qq(std::vector<std::vector<std::complex<double>>> R);
std::vector<std::vector<std::complex<double>>> Rdag_gg(std::vector<std::vector<std::complex<double>>> R);
std::vector<std::vector<std::complex<double>>> mult_qq(std::vector<std::vector<std::complex<double>>> M1, std::vector<std::vector<std::complex<double>>> M2, std::vector<std::vector<std::complex<double>>> M3);
std::vector<std::vector<std::complex<double>>> mult_gg(std::vector<std::vector<std::complex<double>>> M1, std::vector<std::vector<std::complex<double>>> M2, std::vector<std::vector<std::complex<double>>> M3);
std::complex<double> omega_qq(std::vector<std::vector<std::complex<double>>> M1, std::vector<std::vector<std::complex<double>>> cM2, std::vector<std::vector<std::complex<double>>> M3, std::vector<std::vector<std::complex<double>>> M2);
std::complex<double> omega_gg(std::vector<std::vector<std::complex<double>>> M1, std::vector<std::vector<std::complex<double>>> cM2, std::vector<std::vector<std::complex<double>>> M3, std::vector<std::vector<std::complex<double>>> M2);
#endif
