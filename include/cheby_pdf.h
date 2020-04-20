#include <stdlib.h>
#include <gsl/gsl_math.h>
#include "LHAPDF/LHAPDF.h"

#ifndef CHEBYPDF_H
#define CHEBYPDF_H


double xfx_cheby(double x, double* ci, int size_ci);
std::vector<std::complex<double>> cheby_Nspace(double a, double b, std::complex<double> N);
std::vector<double> cheby_xspace(double a, double b, double x);
#endif
