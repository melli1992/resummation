#include <complex>
#include <vector>
#ifndef ttHres_H
#define ttHres_H

double SH_qq_LO(double s, double t13, double t14, double t23, double t24);
double SH_gg_LO(double s, double t13, double t14, double t23, double t24);
std::complex<double> SH_qq_LO_c(std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24);
std::complex<double> SH_gg_LO_c(std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24);

std::complex<double> delidelj_exp(std::complex<double> N, double A1);

std::complex<double> full_qq_res_abs(std::complex<double> N, double s, double t13, double t14, double t23, double t24);
std::complex<double> full_gg_res_abs(std::complex<double> N, double s, double t13, double t14, double t23, double t24);
std::complex<double> pT_qq_res_abs(std::complex<double> N, double pT2, std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24, double s34);
std::complex<double> pT_gg_res_abs(std::complex<double> N, double pT2, std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24, double s34);

std::complex<double> qq_res(std::complex<double> N, std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24, double s34);
std::complex<double> gg_res(std::complex<double> N, std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24, double s34);

#endif
