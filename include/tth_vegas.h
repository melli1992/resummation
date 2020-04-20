#ifndef VEGASttH_H
#define VEGASttH_H

double vegas_ttH_Nspace_deform(double *k, size_t dim, void *params);
double vegas_ttH_Nspace(double *k, size_t dim, void *params);
double vegas_ttH_LO(double *k, size_t dim, void *params);
double vegas_ttH_LO_pT(double *k, size_t dim, void *params);
double vegas_ttH_LO_pT_N(double *k, size_t dim, void *params);
double vegas_ttH_pT_res(double *k, size_t dim, void *params);
double vegas_ttH_LO1st(double *k, size_t dim, void *params);
double vegas_ttH_LO2nd(double *k, size_t dim, void *params);
double vegas_ttH_resum_2nd(double *k, size_t dim, void *params);
double vegas_ttH_resum_2nd_z5(double *k, size_t dim, void *params);
double vegas_ttH_resum_2nd_def(double *k, size_t dim, void *params);
double vegas_ttH_resum_2nd_def_z5(double *k, size_t dim, void *params);
double vegas_ttH_resum_def(double *k, size_t dim, void *params);
double vegas_ttH_resum_def_z5(double *k, size_t dim, void *params);
double vegas_ttH_inv_mass_LO(double *k, size_t dim, void *params);
double vegas_ttH_inv_mass_res(double *k, size_t dim, void *params);
double vegas_ttH_inv_mass_LO_s34(double *k, size_t dim, void *params);
double vegas_ttH_inv_mass_res_s34(double *k, size_t dim, void *params);


double kallen(double x, double y, double z);
std::complex<double> kallen_c(std::complex<double> x, double y, double z);
std::vector<double> xsec_div(double rho,  double s34, double thetaHCM, double thetat, double phit);
std::vector<double> xsec_LO(double rho,  double s34, double thetaHCM, double thetat, double phit);
std::vector<double> xsec_LO_pT(double xT2, double pT2, double s34, double thetat, double phit);
std::vector<std::complex<double>> xsec_pT_res(std::complex<double> N, double xT2, double pT2, double s34, double thetat, double phit);
std::vector<std::complex<double>> xsec_LO_c(std::complex<double> rho, double s34, double thetaHCM, double thetat, double phit);
std::vector<std::complex<double>> xsec_res(std::complex<double> N, double rho, double s34, double thetaHCM, double thetat, double phit);
std::vector<std::complex<double>> xsec_res_c(std::complex<double> N, std::complex<double> rho, double s34, double thetaHCM, double thetat, double phit);
std::vector<std::complex<double>> xsec_res_c_z5(std::complex<double> N, std::complex<double> rho, double s34, double thetaHCM, double thetat, double phit);


#endif
