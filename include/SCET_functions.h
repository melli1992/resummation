#ifndef SCETRESUM_H
#define SCETRESUM_H


// solving lambdaQCD
std::complex<double> invasLambdaQCD(std::complex<double> LambdaQCD);
std::complex<double> DinvasLambdaQCD(std::complex<double> LambdaQCD);
void solveLambdaQCD();

// running alphas
std::complex<double> falphasQ2(std::complex<double> mu2);

// the hard functions
std::vector<std::complex<double>> hard_higgs(std::complex<double> mu2, std::complex<double> muh2);
double cT_higgs();
std::vector<std::complex<double>> hard_DY(std::complex<double> mu2, std::complex<double> muh2);

// acusp, S and tilde(s)
std::vector<std::complex<double>> acusp(double mu, double nu, std::vector<double> gamma);
std::vector<std::complex<double>> sudakov(double mu, double nu, std::vector<double> gamma);
std::vector<std::complex<double>> feta(double mu, double nu, std::vector<double> gamma);
std::vector<double> fsoft_higgs(double L);
std::vector<double> fsoft_DY(double L);

// C coefficients
std::vector<std::complex<double>> C_higgs();
std::vector<std::complex<double>> C_DY();

// eta derivatives
double Deta_z(double z, double eta);
double DDeta_z(double z, double eta);
double Deta_z_eta(double z, double eta, double BNsub);
double DDeta_z_eta(double z, double eta, double BNsub);
double Deta_tau(double eta);
double DDeta_tau(double eta);
double Deta_tau_2nd(double eta, double BNsub);
double DDeta_tau_2nd(double eta, double BNsub);

// cross sections
double diff_xsec_DY(double z, double x);
double xsec_higgs(double z, double x);

#endif
