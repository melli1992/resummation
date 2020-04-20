#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include "parameters.h"
#include <vector>

#ifndef MONTE_H
#define MONTE_H

struct results{double res; double err;};
struct functionint{gsl_monte_function G; std::vector<double> xl; std::vector<double> xu;};
void display_results(std::string title, double &result, double &error, double chisq = 0);
functionint init_vegas_pdfs(std::string lumi = "gg");
functionint init_vegas_ttH(std::string process = "LO");
functionint init_vegas_mellin(std::string process = "test", std::string order = "full");
functionint init_vegas_pf(std::string order, std::string power="full", std::string process = "qqbar", int power_number = 0, bool integrated = false);
results call_vegas(functionint integrand, lumni_params params, bool verbal = false,  bool high_prec = false);
double integrand(double *k, size_t dim, void *params);


#endif
