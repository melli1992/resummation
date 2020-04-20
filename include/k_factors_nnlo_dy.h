#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTORNNLO_H
#define KFACTORNNLO_H
//needed for qqbar channel
double DY_NNLO_qqbar_delta();
double DY_NNLO_qqbar_plus(double x);
double DY_NNLO_qqbar_NS(double x);
double DY_NNLO_qqbar_NS_expansion(double x, int power);
double DY_NNLO_BB_full(double x);
double DY_NNLO_BB_expansion(double x, int power);
double DY_NNLO_BC_full(double x);
double DY_NNLO_BC_expansion(double x, int power);

//needed for qq, qbarqbar
double DY_NNLO_CC_full(double x);
double DY_NNLO_CC_expansion(double x, int power);
double DY_NNLO_CD_full(double x);
double DY_NNLO_CD_expansion(double x, int power);
double DY_NNLO_CE_full(double x);
double DY_NNLO_CE_expansion(double x, int power);
double DY_NNLO_CF_full(double x);
double DY_NNLO_CF_expansion(double x, int power);


//needed for gg
double DY_NNLO_gg_full(double x);
double DY_NNLO_gg_expansion(double x, int power);

//needed for qg
double DY_NNLO_qg_full(double x);
double DY_NNLO_qg_expansion(double x, int power);

#endif
