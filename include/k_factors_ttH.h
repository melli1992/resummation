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

std::complex<double> H22_qq_c(std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24);
std::complex<double> H22_gg_c(std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24);
std::complex<double> H23_gg_c(std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24);
std::complex<double> H33_gg_c(std::complex<double> s, std::complex<double> t13, std::complex<double> t14, std::complex<double> t23, std::complex<double> t24);

#endif
