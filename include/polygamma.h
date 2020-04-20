#include <vector>

#ifndef POLYGAM_H
#define POLYGAM_H

double li2(double x);
std::complex<double> dilog(std::complex<double> x); // when it is complex!
double Li3(double x);
double S12(double x);
double clenshaw(std::vector<double> coeff, double x);
std::complex<double> Gamma(std::complex<double> x);

#endif
