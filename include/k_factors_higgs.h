#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTORHIGGS_H
#define KFACTORHIGGS_H

//LO part
double higgs_LO_factor();
std::complex<double> AQ(double x);

//NLO part
//needed for gg
double higgs_NLO_gg_delta();
double higgs_NLO_gg_reg(double x);
double higgs_NLO_gg_plus(double x);
double higgs_NLO_gg_expansion(double x, int power);

//needed for qg
double higgs_NLO_qg_full(double x);
double higgs_NLO_qg_expansion(double x, int power);

//needed for qqbar
double higgs_NLO_qqbar_full(double x);
double higgs_NLO_qqbar_expansion(double x, int power);

#endif
