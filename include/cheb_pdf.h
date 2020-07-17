#include <iostream>
#include <fstream>
#include <complex>
#include <string>

#ifndef CHEBPDF_H //need this otherwise compiler error that things were predefined (called a guard)
#define CHEBPDF_H

double easy_func (double x);
int cheb_approx_try(unsigned int Nch, const double umin, const double umax, double *c);
int ApproxLuminosity(double *coeff, double tau_min, unsigned int N, std::string LumChannel);
double Fueval(double u, double umin, double umax, unsigned int N, double *ccoeff);
std::complex<double> LumN(std::complex<double> MelMoment, unsigned int N, std::string LumChannel);


#endif
