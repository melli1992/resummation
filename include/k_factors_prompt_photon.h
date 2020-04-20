#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTORPF_H
#define KFACTORPF_H

//needed for qqbar channel (identical quarks)
double vegas_FP_qqbartogg_LP(double *k, size_t dim, void *params);
double vegas_FP_qqbartogg_LP_corr(double *k, size_t dim, void *params);
double vegas_FP_qqbartogg_full(double *k, size_t dim, void *params);
double vegas_FP_qqbartogg_delta(double *k, size_t dim, void *params);
double vegas_FP_qqbartogg_power(double *k, size_t dim, void *params);
double FP_qqbartogg_LP(double v, double w);
double FP_qqbartogg_LP_corr(double v, double wlow);
double FP_qqbartogg_delta(double v);
double FP_qqbartogg_full(double v, double w);
double FP_qqbartogg_expansion(double v, double w, int power);

#endif
