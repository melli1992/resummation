#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTOR_DY
#define KFACTOR_DY

//LO part
double DY_LO_factor();

//NLO part

//needed for gg
double DY_NLO_qqbar_delta();
double DY_NLO_qqbar_reg(double x);
double DY_NLO_qqbar_plus(double x);
double DY_NLO_qqbar_expansion(double x, int power);

//needed for qg
double DY_NLO_qg_full(double x);
double DY_NLO_qg_expansion(double x, int power);



//BSM part
double BSM_WR_LO_factor();

#endif
