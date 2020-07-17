#include <complex>
#include <vector>
#ifndef ttHsoftanom_H
#define ttHsoftanom_H
double beta34(double s34);
std::complex<double> log_beta34(double s34);
std::complex<double> Tij(std::complex<double> tij, std::complex<double> s);
std::complex<double> Omega3(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> Lambda3(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);

std::complex<double> lambda_qq_11(double s34);
std::complex<double> lambda_qq_12(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_qq_21(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_qq_22(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, double s34, std::complex<double> s);
std::complex<double> lambda_gg_11(double s34);
std::complex<double> lambda_gg_12();
std::complex<double> lambda_gg_13(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_gg_21();
std::complex<double> lambda_gg_22(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s, double s34);
std::complex<double> lambda_gg_23(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_gg_31(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_gg_32(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s);
std::complex<double> lambda_gg_33(std::complex<double> t13, std::complex<double> t24, std::complex<double> t14, std::complex<double> t23, std::complex<double> s, double s34);

std::complex<double> eigen_qq_1(std::complex<double>  lambdaqq11, std::complex<double>  lambdaqq22, std::complex<double>  lambdaqq12, std::complex<double>  lambdaqq21);
std::complex<double> eigen_qq_2(std::complex<double>  lambdaqq11, std::complex<double>  lambdaqq22, std::complex<double>  lambdaqq12, std::complex<double>  lambdaqq21);
std::complex<double> R_qq_1i(std::complex<double>  lambdaqq11, std::complex<double>  lambdaqq22, std::complex<double>  lambdaqq12, std::complex<double>  lambdaqq21, std::complex<double> eigeni);
std::complex<double> R_qq_21();
std::complex<double> R_qq_22();
std::complex<double> Rinv_qq_const(std::complex<double>  lambdaqq11, std::complex<double>  lambdaqq22, std::complex<double>  lambdaqq12, std::complex<double>  lambdaqq21, std::complex<double> eigen1, std::complex<double> eigen2);

std::complex<double> a(std::complex<double> gamma11, std::complex<double> gamma22);
std::complex<double> b(std::complex<double> gamma11, std::complex<double> gamma22, std::complex<double> gamma31);
std::complex<double> c(std::complex<double> gamma11, std::complex<double> gamma22, std::complex<double> gamma31);
std::complex<double> collect(std::complex<double> a, std::complex<double> b, std::complex<double> c);
std::complex<double> eigen1_new(std::complex<double> a, std::complex<double> b, std::complex<double> Z);
std::complex<double> eigen2_new(std::complex<double> a, std::complex<double> b, std::complex<double> Z);
std::complex<double> eigen3_new(std::complex<double> a, std::complex<double> b, std::complex<double> Z);
std::complex<double> R_gg_1i(std::complex<double> lambda31, std::complex<double> eigeni, std::complex<double> lambda11);
std::complex<double> R_gg_2i(std::complex<double> lambda31, std::complex<double> eigeni, std::complex<double> lambda22);
std::complex<double> Rinv_gg_const(std::complex<double> v1_1, std::complex<double> v2_1, std::complex<double> v3_1, std::complex<double> v1_2, std::complex<double> v2_2, std::complex<double> v3_2);

std::complex<double> Check_Roots(std::complex<double> eigeni,std::complex<double> gamma11, std::complex<double> gamma31, std::complex<double> gamma22);

void Check_Eigen(std::complex<double> eigeni, std::complex<double> vi_1, std::complex<double> vi_2, std::complex<double> vi_3, std::complex<double> gamma31, std::complex<double> gamma11, std::complex<double> gamma22, std::complex<double> gamma13, std::complex<double> gamma23, std::complex<double> gamma32, std::complex<double> gamma33);


std::vector<std::complex<double>> eigenvalues_qq(std::complex<double> G11, std::complex<double> G12, std::complex<double> G21, std::complex<double> G22);
std::vector<std::complex<double>> eigenvalues_gg(std::complex<double> G11, std::complex<double> G13, std::complex<double> G22, std::complex<double> G23,  std::complex<double> G31, std::complex<double> G32, std::complex<double> G33);
std::vector<std::vector<std::complex<double>>> eigenvectors_qq(std::vector<std::complex<double>> lambda, std::complex<double> G11, std::complex<double> G12, std::complex<double> G21, std::complex<double> G22);
std::vector<std::vector<std::complex<double>>> eigenvectors_gg(std::vector<std::complex<double>> lambda, std::complex<double> G11, std::complex<double> G13, std::complex<double> G22, std::complex<double> G23,  std::complex<double> G31, std::complex<double> G32, std::complex<double> G33);
std::vector<std::vector<std::complex<double>>> R_qq(std::vector<std::vector<std::complex<double>>> ev);
std::vector<std::vector<std::complex<double>>> R_gg(std::vector<std::vector<std::complex<double>>> ev);
std::vector<std::vector<std::complex<double>>> Rinv_qq(std::vector<std::vector<std::complex<double>>> ev);
std::vector<std::vector<std::complex<double>>> Rinv_gg(std::vector<std::vector<std::complex<double>>> ev);
std::vector<std::vector<std::complex<double>>> Rdag_qq(std::vector<std::vector<std::complex<double>>> R);
std::vector<std::vector<std::complex<double>>> Rdag_gg(std::vector<std::vector<std::complex<double>>> R);
std::vector<std::vector<std::complex<double>>> mult_qq(std::vector<std::vector<std::complex<double>>> M1, std::vector<std::vector<std::complex<double>>> M2, std::vector<std::vector<std::complex<double>>> M3);
std::vector<std::vector<std::complex<double>>> mult_gg(std::vector<std::vector<std::complex<double>>> M1, std::vector<std::vector<std::complex<double>>> M2, std::vector<std::vector<std::complex<double>>> M3);
#endif
