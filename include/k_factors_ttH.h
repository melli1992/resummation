#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

template <typename T>
struct identity_t { typedef T type; };

/// Make working with std::complex<> nubmers suck less... allow promotion.
#define COMPLEX_OPS(OP)                                                 \
  template <typename _Tp>                                               \
  std::complex<_Tp>                                                     \
  operator OP(std::complex<_Tp> lhs, const typename identity_t<_Tp>::type & rhs) \
  {                                                                     \
    return lhs OP rhs;                                                  \
  }                                                                     \
  template <typename _Tp>                                               \
  std::complex<_Tp>                                                     \
  operator OP(const typename identity_t<_Tp>::type & lhs, const std::complex<_Tp> & rhs) \
  {                                                                     \
    return lhs OP rhs;                                                  \
  }
COMPLEX_OPS(+)
COMPLEX_OPS(-)
COMPLEX_OPS(*)
COMPLEX_OPS(/)
#undef COMPLEX_OPS

#ifndef HARDFUNCttH_H
#define HARDFUNCttH_H

double H22_qq(double s, double t13, double t14, double t23, double t24);
double H11_gg(double s, double t13, double t14, double t23, double t24);
double H22_gg(double s, double t13, double t14, double t23, double t24);
double H23_gg(double s, double t13, double t14, double t23, double t24);
double H31_gg(double s, double t13, double t14, double t23, double t24);
double H33_gg(double s, double t13, double t14, double t23, double t24);
double S11_qq();
double S22_qq();
double S11_gg();
double S22_gg();
double S33_gg();
double SH_qq_LO(double s, double t13, double t14, double t23, double t24);
double SH_gg_LO(double s, double t13, double t14, double t23, double t24);
std::complex<double> SH_qq_res(std::complex<double> N, double s, double t13, double t14, double t23, double t24);
std::complex<double> SH_gg_res(std::complex<double> N, double s, double t13, double t14, double t23, double t24);

std::complex<double> H22_qq_c(std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24);
std::complex<double> H23_gg_c(std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24);
std::complex<double> H33_gg_c(std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24);
std::complex<double> SH_qq_LO_c(std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24);
std::complex<double> SH_gg_LO_c(std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24);
std::complex<double> SH_qq_res_c(std::complex<double> N, std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24);
std::complex<double> SH_gg_res_c(std::complex<double> N, std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24);

std::complex<double> g1_q(std::complex<double> lambda);
std::complex<double> g1_g(std::complex<double> lambda);
std::complex<double> g2_q(std::complex<double> lambda);
std::complex<double> g2_g(std::complex<double> lambda);
std::complex<double> log_delta_q(std::complex<double> N);
std::complex<double> log_delta_g(std::complex<double> N);
std::vector<double> D1_qq();
std::vector<double> D1_gg();
std::vector<std::complex<double> > h2_qqtt(std::complex<double> lambda);
std::vector<std::complex<double> > h2_ggtt(std::complex<double> lambda);
std::vector<std::complex<double> > log_delta_qq(std::complex<double> N);
std::vector<std::complex<double> > log_delta_gg(std::complex<double> N);
std::complex<double> omega_gg_res(std::complex<double> N, std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24, double s34);
std::complex<double> omega_qq_res(std::complex<double> N, std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24, double s34);

#endif
