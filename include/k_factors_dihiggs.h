#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTORDIHIGGS_H
#define KFACTORDIHIGGS_H

//LO part
std::complex<double> Cab(double s, double mQ2);
std::complex<double> Cac(double t, double mHc2, double mQ2);
std::complex<double> Cad(double u, double mHd2, double mQ2);
std::complex<double> Cbc(double u, double mHc2, double mQ2);
std::complex<double> Cbd(double t, double mHd2, double mQ2);
std::complex<double> Ccd(double s, double mHc2, double mHd2, double mQ2);
std::complex<double> Dabc(double s, double t, double mHc2, double mHd2, double mQ2);
std::complex<double> Dacb(double t, double u, double mHc2, double mHd2, double mQ2);
std::complex<double> Dbac(double s, double t, double mHc2, double mHd2, double mQ2);

std::complex<double> Ftriangle_scalar_scalar(double s, double mQ2, std::complex<double> cab);
std::complex<double> Fbox_scalar_scalar(double s, double t, double u, double mQ2, double mHc2, double mHd2, std::complex<double> cab, std::complex<double> cac , std::complex<double> cbc, std::complex<double> cad, std::complex<double> cbd, std::complex<double> dabc, std::complex<double> dbac, std::complex<double> dacb);
std::complex<double> Gbox_scalar_scalar(double s, double t, double u, double mQ2, double mHc2, double mHd2, std::complex<double> cab, std::complex<double> cac , std::complex<double> cbc, std::complex<double> cad, std::complex<double> cbd, std::complex<double> ccd, std::complex<double> dabc, std::complex<double> dbac, std::complex<double> dacb);
std::complex<double> FtriangleA_pseudo_scalar(double mQ2, std::complex<double> cab);
std::complex<double> FtriangleZ_pseudo_scalar(double s, double mQ2, std::complex<double> cab, double mAc2, double mHd2);
std::complex<double> Fbox_pseudo_scalar(double s, double t, double u, double mQ2, double mAc2, double mHd2, std::complex<double> cab, std::complex<double> cac, std::complex<double> cbc, std::complex<double> cad, std::complex<double> cbd, std::complex<double> dabc, std::complex<double> dbac, std::complex<double> dacb);
std::complex<double> Gbox_pseudo_scalar(double s, double t, double u, double mQ2, double mAc2, double mHd2, std::complex<double> cab, std::complex<double> cac, std::complex<double> cbc, std::complex<double> cad, std::complex<double> cbd, std::complex<double> ccd, std::complex<double> dabc, std::complex<double> dbac, std::complex<double> dacb);
std::complex<double> Ftriangle_pseudo_pseudo(double s, double mQ2, std::complex<double> cab);
std::complex<double> Fbox_pseudo_pseudo(double s, double t, double u, double mQ2, double mAc2, double mAd2, std::complex<double> cab, std::complex<double> cac , std::complex<double> cbc, std::complex<double> cad, std::complex<double> cbd, std::complex<double> dabc, std::complex<double> dbac, std::complex<double> dacb);
std::complex<double> Gbox_pseudo_pseudo(double s, double t, double u, double mQ2, double mAc2, double mAd2, std::complex<double> cab, std::complex<double> cac , std::complex<double> cbc, std::complex<double> cad, std::complex<double> cbd, std::complex<double> ccd, std::complex<double> dabc, std::complex<double> dbac, std::complex<double> dacb);

double dihiggs_LO_factor_SM(double scale2, double ctheta1);
double dihiggs_LO_factor_approx(double scale2, double ctheta1);
double dihiggs_hh(double scale2, double ctheta);
double dihiggs_hH(double scale2, double ctheta);
double dihiggs_HH(double scale2, double ctheta);
double dihiggs_Ah(double scale2, double ctheta);
double dihiggs_AH(double scale2, double ctheta);
double dihiggs_AA(double scale2, double ctheta);
void test_LO(double scale2, double mC2, double mD2, double ctheta);
#endif
