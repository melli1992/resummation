#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>

#ifndef MELLIN_H
#define MELLIN_H

// test functions
double vegas_fofx2_full(double *k, size_t dim, void *params);
double vegas_fofx2_deriv(double *k, size_t dim, void *params);
double vegas_fofx2_defor(double *k, size_t dim, void *params);
double vegas_fofx2_Nspace(double *k, size_t dim, void *params);
#endif
