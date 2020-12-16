#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTORHIGGSNNLO_H
#define KFACTORHIGGSNNLO_H

//needed for gg
double higgs_NNLO_gg_delta();
double higgs_NNLO_gg_reg(double x);
double higgs_NNLO_gg_plus(double x);
double higgs_NNLO_gg_expansion(double x, int power);
double logdep_gg_reg(double x);
double logdep_gg_constant();
double logdep_gg_plus(double x);
double logdep_gg_reg2(double x);
double logdep_gg_reg3(double x);

double higgs_scale_factor(double LO, double NLO);

double nnlo_reg(double z);
double FO_1_log_PL(double z);
double FO_1_nonlog(double z);
double FO_1_delta_const();
double FO_2_log(double z);
double FO_2_delta_const();
double FO_2_nonlog(double z);



// qg (+qbar g)
double higgs_NNLO_qg_reg(double x);
double higgs_NNLO_qg_expansion(double x, int power);
double logdep_qg(double x);

// qq (+qbar qbar)
double higgs_NNLO_qq_reg(double x);
double higgs_NNLO_qq_expansion(double x, int power);
double logdep_qq(double x);

// qq' (+qbar qbar')
double higgs_NNLO_qqp_reg(double x);
double higgs_NNLO_qqp_expansion(double x, int power);
double logdep_qqp(double x);

// qqbar
double higgs_NNLO_qqbar_reg(double x);
double higgs_NNLO_qqbar_expansion(double x, int power);
double logdep_qqbar(double x);


#endif
