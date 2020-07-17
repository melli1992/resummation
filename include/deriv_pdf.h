#include <stdlib.h>
#include <gsl/gsl_math.h>
#include "LHAPDF/LHAPDF.h"

#ifndef DERIVPDF_H
#define DERIVPDF_H

// charge weighted (with photon)
double pdf_sum_qqbar_charge_weighted(double x, double tau_over_z);
double pdf_sum_qqbar_charge_unweighted(double x, double tau_over_z);
double pdf_sum_qq_charge_weighted_double(double x, double tau_over_z);
double pdf_sum_qq_charge_weighted_single(double x, double tau_over_z);
double pdf_sum_qq_charge_weighted_double_vivj(double x, double tau_over_z);
double pdf_sum_qq_charge_weighted_single_vivi(double x, double tau_over_z);
double pdf_sum_qg_charge_weighted(double x, double tau_over_z);
double pdf_sum_gg_charge_weighted(double x, double tau_over_z);

// no charge (without photon)
double pdf_sum_qqbar(double x, double tau_over_z);
double pdf_sum_qg(double x, double tau_over_z);
double pdf_sum_gg(double x, double tau_over_z);
double pdf_sum_qq(double x, double tau_over_z);
double pdf_sum_qqNI(double x, double tau_over_z);

// up and down seperately
double pdf_sum_qqbarUP(double x, double tau_over_z);
double pdf_sum_qqbarDOWN(double x, double tau_over_z);

//qqbar' for W
double pdf_sum_qqbarNI(double x, double tau_over_z);


//derivatives (only with photon)
double derivative_qg_pdf(double x, double z, double tau, double eps= 0.0001);
double derivative_gg_pdf(double x, double tau, double eps= 0.0001);
double derivative_qg_pdf_jac(double x, double z, double tau, double eps= 0.0001);
double vegas_pdf_up_minus_upbar(double *k, size_t dim, void *params);
double vegas_pdf_mom_consv(double *k, size_t dim, void *params);

// lumi functions
double vegas_lumi_gg(double *k, size_t dim, void *params);
double vegas_lumi_qg(double *k, size_t dim, void *params);
double vegas_lumi_qqbar(double *k, size_t dim, void *params);
double vegas_lumi_qq(double *k, size_t dim, void *params);
double vegas_lumi_qqNI(double *k, size_t dim, void *params);

// for ttH
double single_div_sophis(int i, double x);
double double_div_sophis(int i, double x);
double single_div(int i, double x);
double double_div(int i, double x);
double gluon_d0PDF(double x1, double x2);
double gluon_d1PDF(double x1, double x2);
double gluon_d2PDF(double x1, double x2);
double quark_d0PDF(double x1, double x2);
double quark_d1PDF(double x1, double x2);
double quark_d2PDF(double x1, double x2);

//for SCET
double lumni(std::string channel, double t, double y);
std::vector<double> deriv_to_y_luminosities(std::string channel, double x, double y);


#endif
