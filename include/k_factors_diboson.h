#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>


#ifndef KFACTORDIBOSON_H
#define KFACTORDIBOSON_H
double cts(double s, double qi, double giL);
double css(double s, double qi, double giL, double giR);
double Fu(double s, double beta);
double Ju(double s, double beta);
double Ku(double s, double beta);
double Acoupling(double s, double qi, double gv, double ga);
double Icoupling(double s, double qi, double gv, double ga);
double Ecoupling();
double Astu(double s, double beta);
double Istu(double s, double beta);
double Estu(double s, double beta);
double partonic_up_wpwm(double s);
double partonic_down_wpwm(double s);
std::complex<double> cpartonic_up_wpwm(std::complex<double> s);
std::complex<double> cpartonic_down_wpwm(std::complex<double> s);
double partonic_up_zz(double s);
double partonic_down_zz(double s);
std::complex<double> cpartonic_up_zz(std::complex<double> s);
std::complex<double> cpartonic_down_zz(std::complex<double> s);
#endif
