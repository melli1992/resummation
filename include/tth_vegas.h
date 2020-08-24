#ifndef VEGASttH_H
#define VEGASttH_H

double vegas_ttH_Nspace_deform(double *k, size_t dim, void *params);
double vegas_ttH_Nspace(double *k, size_t dim, void *params);
double vegas_ttH_LO(double *k, size_t dim, void *params);
double vegas_ttH_LO_pT(double *k, size_t dim, void *params);
double vegas_ttH_LO_pT_N(double *k, size_t dim, void *params);
double vegas_ttH_pT_res_Nfix(double *k, size_t dim, void *params);
double vegas_ttH_LO2nd(double *k, size_t dim, void *params);
double vegas_ttH_N1(double *k, size_t dim, void *params);
double vegas_ttH_N2(double *k, size_t dim, void *params);
double vegas_ttH_N3(double *k, size_t dim, void *params);
double vegas_ttH_N4(double *k, size_t dim, void *params);
double vegas_ttH_N5(double *k, size_t dim, void *params);
double vegas_ttH_inv_mass_N1(double *k, size_t dim, void *params);
double vegas_ttH_inv_mass_N2(double *k, size_t dim, void *params);
double vegas_ttH_inv_mass_N3(double *k, size_t dim, void *params);
double vegas_ttH_pT_N4(double *k, size_t dim, void *params);
double vegas_ttH_pT_N5(double *k, size_t dim, void *params);
double vegas_ttH_sttpT_N4(double *k, size_t dim, void *params);
double vegas_ttH_sttpT_N5(double *k, size_t dim, void *params);
double vegas_ttH_LO_pT_stt_dist(double *k, size_t dim, void *params);

double kallen(double x, double y, double z);
std::complex<double> kallen_c(std::complex<double> x, double y, double z);
std::vector<double> xsec_div(double rho,  double s34, double thetaHCM, double thetat, double phit);
std::vector<double> xsec_LO(double rho,  double s34, double thetaHCM, double thetat, double phit);
std::vector<double> xsec_LO_pT(double xT2, double pT2, double s34, double thetat, double phit);
std::vector<double> xsec_LO_pT_stt(double xT2, double pT2, double s34, double thetat, double phit);
std::vector<std::complex<double>> xsec_LO_pT_stt_c(std::complex<double> xT2, double s34, double thetaHCM, double thetat, double phit);
std::vector<std::complex<double>> xsec_pT_res_stt_c(std::complex<double> N,  std::complex<double> xT2, double s34, double thetaHCM, double thetat, double phit);
std::vector<std::complex<double>> xsec_pT_res(std::complex<double> N, double xT2, double pT2, double s34, double thetat, double phit);
std::vector<std::complex<double>> xsec_pT_res_abs(std::complex<double> N, double xT2, double pT2, double s34, double thetat, double phit);
std::vector<std::complex<double>> xsec_pT_res_Nfix(std::complex<double> N, double xT2, double pT2, double s34, double thetat, double phit);
std::vector<std::complex<double>> xsec_LO_c(std::complex<double> rho, double s34, double thetaHCM, double thetat, double phit);
std::vector<std::complex<double>> xsec_res(std::complex<double> N, double rho, double s34, double thetaHCM, double thetat, double phit);
std::vector<std::complex<double>> xsec_res_z5(std::complex<double> N, double rho, double s34, double thetaHCM, double thetat, double phit);
std::vector<std::complex<double>> xsec_res_c(std::complex<double> N, std::complex<double> rho, double s34, double thetaHCM, double thetat, double phit);
std::vector<std::complex<double>> xsec_res_c_z5(std::complex<double> N, std::complex<double> rho, double s34, double thetaHCM, double thetat, double phit);

int try_stt(double xT2, double pT2, double s34, double thetat, double phit);

#endif
