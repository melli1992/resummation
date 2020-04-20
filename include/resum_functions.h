#ifndef RESUM_H
#define RESUM_H

// the resummation coefficients
std::complex<double> g1(double A1, std::complex<double>lambda);
std::complex<double> h1NLP(double A1, std::complex<double> N, std::complex<double>lambda);
std::complex<double> g2(double A1,double A2,std::complex<double>lambda);
std::complex<double> g3(double A1,double A2,double A3,std::complex<double>lambda);
std::complex<double> wideangle(double D2,std::complex<double>lambda);


std::complex<double> NLOmatch(std::complex<double> N, double plusone, double A1, double g0=0.);
std::complex<double> NNLOmatch(std::complex<double> N, double plusone, double A1, double A2, double D2, double g0=0.);

double higgs_g01();
double DY_g01();
double DY_g02();
double FDY();

// expansions of the resummed functions for large N at LL only (LP and NLP)
std::complex<double> LP_LL_function_expanded_LP(std::complex<double> N, double Col_Fac);
std::complex<double> LP_LL_function_expanded_NLP(std::complex<double> N, double Col_Fac);
std::complex<double> NLP_LL_function_expanded_NLP(std::complex<double> N, double Col_Fac);
std::complex<double> LP_LL_function_expanded(std::complex<double> N, double Col_Fac);
std::complex<double> NLP_LL_function_expanded(std::complex<double> N, double Col_Fac);
std::complex<double> LP_LL_function_full(std::complex<double> N, double Col_Fac);
std::complex<double> NLP_LL_function_full(std::complex<double> N, double Col_Fac);

#endif
