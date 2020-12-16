#include "tth_softanom.h"
#include "resum_functions.h"
#include "k_factors_ttH.h"
#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "parameters.h"
#include "resum_tth.h"

using namespace std;

double SH_qq_LO(double s, double t13, double t14, double t23, double t24){
	return H22_qq(s, t13, t14, t23, t24)*S22_qq();
}
complex<double> SH_qq_LO_c(complex<double>  s, complex<double> t13, complex<double> t14, complex<double> t23, complex<double> t24){
	return H22_qq_c(s, t13, t14, t23, t24)*S22_qq();
}
double SH_gg_LO(double s, double t13, double t14, double t23, double t24){
	return H22_gg(s, t13, t14, t23, t24)*(S22_gg()+S11_gg()/(pow(CA,2))) + H33_gg(s, t13, t14, t23, t24)*S33_gg();
}
complex<double> SH_gg_LO_c(complex<double>  s, complex<double> t13, complex<double> t14, complex<double> t23, complex<double> t24){
	return H22_gg_c(s, t13, t14, t23, t24)*(S22_gg()+S11_gg()/(pow(CA,2))) + H33_gg_c(s, t13, t14, t23, t24)*S33_gg();
}

// LP LL function h0 (or g1) hep-ph/0306211 eqn 39 (note that gammaE is not part there of lambda and the factor 2 (as they have 2*h0 = g1!) or 1905.11771 eqn 6
complex<double> g1_M2(double A1,complex<double>lambda){
	return A1/(2.*M_PI*pow(b0,2))*(2.*lambda+(1.-2.*lambda)*log(1.-2.*lambda));
}

/////////////////////////////////////////////////////////////
// LP NLL function h1 (or g2) hep-ph/0306211 (note factor 2 and gammaE) eqn 40 or 1905.11771 eqn 61
// (checked with mathematica)
complex<double> g2_M2(double A1,double A2,complex<double>lambda){
	double INCeuler = 0.;
	if(INCEULER == 0) {INCeuler = 1.;}
	return 1./(2.*M_PI*b0)*(-A2/(M_PI*b0)+A1*log(M2/muR2))*(2.*lambda+log(1.-2.*lambda))
	+ A1*b1/(2.*M_PI*pow(b0,3))*(2.*lambda+log(1.-2.*lambda)+1./2.*pow(log(1.-2.*lambda),2))
	- A1/(M_PI*b0)*lambda*log(M2/muF2)
	- INCeuler*2.*M_gammaE*log(1.-2.*lambda)*A1/(2.*M_PI*b0);
}

//expanded version of the wide-angle
complex<double> delidelj_exp(complex<double> N, double A1){
	//if(!INCEULER)
	return alphas_muR*2.*A1/M_PI*log(N)*(log(N)-ISNLL*log(M2/muF2));//+ISNLL*INCEULER*M_gammaE);
}

complex<double> full_qq_res_abs(complex<double> N, double s, double t13, double t14, double t23, double t24){
	complex<double> lambda = b0*alphas_muR*log(N);
	complex<double> wide_soft = ISNLL*CA*(-1.)*log(1.-2.*lambda)/(2.*M_PI*b0);
	complex<double> H22qq = H22_qq_c(s, t13, t14, t23, t24);
	double S22qq = S22_qq();
	complex<double> result = H22qq*S22qq*exp(wide_soft)
														*exp(2.*(1./alphas_muR*ISLL*g1_M2(A1q,lambda)+ISNLL*g2_M2(A1q,A2q,lambda)));
	if(expansion){
		 	complex<double> wide_soft2 = 1.-ISNLL*CA*(-1.)*alphas_muR*log(N)/(M_PI);
			complex<double> exponent = delidelj_exp(N,A1q);
			//cout << "difference in result qq " << result << " " << (H22qq*S22qq*wide_soft2+H22qq*S22qq*delidelj_exp(N,A1q)) << endl;
			return result - (H22qq*S22qq*(wide_soft2+exponent));
	}
	return result;
}
complex<double> full_gg_res_abs(complex<double> N, double s, double t13, double t14, double t23, double t24)
{
	complex<double> lambda = b0*alphas_muR*log(N);
	complex<double> wide_soft = ISNLL*CA*(-1.)*log(1.-2.*lambda)/(2.*M_PI*b0);
	complex<double> H22gg = H22_gg_c(s, t13, t14, t23, t24);
	complex<double> H33gg = H33_gg_c(s, t13, t14, t23, t24);
	double S11gg = S11_gg();
	double S22gg = S22_gg();
	double S33gg = S33_gg();

	complex<double> result = (H22gg*(S22gg*exp(wide_soft)+S11gg/(pow(CA,2)))+ exp(wide_soft)*H33gg*S33gg)
														*exp(2.*(1./alphas_muR*ISLL*g1_M2(A1g,lambda)+ISNLL*g2_M2(A1g,A2g,lambda)));
	if(expansion){
		complex<double> wide_soft2 = 1.-ISNLL*CA*(-1.)*alphas_muR*log(N)/(M_PI);
		//cout << "difference in result gg " << result << " " << ((H22gg*(S22gg*wide_soft2+S11gg/(pow(CA,2)))+ wide_soft2*H33gg*S33gg)
		//									+ (H22gg*(S22gg+S11gg/(pow(CA,2)))+ H33gg*S33gg)*delidelj_exp(N,A1g)) << endl;
		complex<double> exponent = delidelj_exp(N,A1g);
		return result - ((H22gg*(S22gg*(exponent+wide_soft2)+S11gg/(pow(CA,2))*(1.+exponent)))+ (exponent+wide_soft2)*H33gg*S33gg);
					}
	return result;
}

complex<double> pT_qq_res_abs(complex<double> N, double pT2, complex<double> s, complex<double> t13, complex<double> t14, complex<double> t23, complex<double> t24, double s34){
	complex<double> lambda = b0*alphas_muR*log(N);
	complex<double> wide_soft = ISNLL*CA*(log(1.+pT2/(4.*mt2))-1.)*log(1.-2.*lambda)/(2.*M_PI*b0);
	complex<double> H22qq = H22_qq_c(s, t13, t14, t23, t24);
	double S22qq = S22_qq();
	complex<double> result = H22qq*S22qq*exp(wide_soft)
														*exp(2.*(1./alphas_muR*ISLL*g1_M2(A1q,lambda)+ISNLL*g2_M2(A1q,A2q,lambda)));
	if(expansion){
		 	complex<double> wide_soft2 = 1.-ISNLL*CA*(log(1.+pT2/(4.*mt2))-1.)*alphas_muR*log(N)/(M_PI);
			complex<double> exponent = delidelj_exp(N,A1q);
			//cout << "difference in result qq " << result << " " << (H22qq*S22qq*wide_soft2+H22qq*S22qq*delidelj_exp(N,A1q)) << endl;
			return result - (H22qq*S22qq*(wide_soft2+exponent));
	}
	return result;
}

complex<double> pT_gg_res_abs(complex<double> N, double pT2, complex<double> s, complex<double> t13, complex<double> t14, complex<double> t23, complex<double> t24, double s34){
	complex<double> lambda = b0*alphas_muR*log(N);
	complex<double> wide_soft = ISNLL*CA*(log(1.+pT2/(4.*mt2))-1.)*log(1.-2.*lambda)/(2.*M_PI*b0);

	complex<double> H22gg = H22_gg_c(s, t13, t14, t23, t24);
	complex<double> H33gg = H33_gg_c(s, t13, t14, t23, t24);
	double S11gg = S11_gg();
	double S22gg = S22_gg();
	double S33gg = S33_gg();

	complex<double> result = (H22gg*(S22gg*exp(wide_soft)+S11gg/(pow(CA,2)))+ exp(wide_soft)*H33gg*S33gg)
														*exp(2.*(1./alphas_muR*ISLL*g1_M2(A1g,lambda)+ISNLL*g2_M2(A1g,A2g,lambda)));
	if(expansion){
		complex<double> wide_soft2 = 1.-ISNLL*CA*(log(1.+pT2/(4.*mt2))-1.)*alphas_muR*log(N)/(M_PI);
		//cout << "difference in result gg " << result << " " << ((H22gg*(S22gg*wide_soft2+S11gg/(pow(CA,2)))+ wide_soft2*H33gg*S33gg)
		//									+ (H22gg*(S22gg+S11gg/(pow(CA,2)))+ H33gg*S33gg)*delidelj_exp(N,A1g)) << endl;
		complex<double> exponent = delidelj_exp(N,A1g);
		return result - ((H22gg*(S22gg*(exponent+wide_soft2)+S11gg/(pow(CA,2))*(1.+exponent)))+ (exponent+wide_soft2)*H33gg*S33gg);
					}
	return result;
}


complex<double> qq_res(complex<double> N, complex<double> s, complex<double> t13, complex<double> t14, complex<double> t23, complex<double> t24, double s34){

	vector<vector<complex<double> > > HR_IJ = {{0,0},{0,0}}; //hard matrix in R representatino
	vector<vector<complex<double> > > SR_IJ = {{0,0},{0,0}}; //soft matrix in R representatino
	vector<vector<complex<double> > > GammaR_IJ = {{0,0},{0,0}}; //anomalous dim in R representation (needs to be diagonal!)
	vector<vector<complex<double> > > H_IJ = {{0,0},{0,0}}; //hard matrix
	vector<vector<complex<double> > > S_IJ = {{0,0},{0,0}}; //soft matrix
	vector<vector<complex<double> > > Gamma_IJ = {{0,0},{0,0}}; //anomalous dim
// filling hard function
	H_IJ[0][0] = 0;
	H_IJ[0][1] = 0;
	H_IJ[1][0] = 0;
	H_IJ[1][1] = H22_qq_c(s, t13, t14, t23, t24);
// filling soft function
	S_IJ[0][0] = S11_qq();
	S_IJ[0][1] = 0;
	S_IJ[1][0] = 0;
	S_IJ[1][1] = S22_qq();
// filling anomalous dim
	complex<double> G11 = lambda_qq_11(s34);
	complex<double> G22 = lambda_qq_22(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s34, s);
	complex<double> G12 = lambda_qq_12(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s);
	complex<double> G21 = lambda_qq_21(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s);
	Gamma_IJ[0][0] = G11;
	Gamma_IJ[0][1] = G12;
	Gamma_IJ[1][0] = G21;
	Gamma_IJ[1][1] = G22;
	vector<complex<double>> eval = eigenvalues_qq(G11, G12, G21, G22);
	if(abs(real(G11)) == 0){eval[0] = G11;}
	vector<vector<complex<double>>> evec = eigenvectors_qq(eval, G11, G12, G21, G22);
	vector<vector<complex<double>>> resum = {{0,0},{0,0}};
	vector<vector<complex<double>>> resumC = {{0,0},{0,0}};
	vector<vector<complex<double>>> resumE = {{0,0},{0,0}};
	vector<vector<complex<double>>> resumEC = {{0,0},{0,0}};
	vector<vector<complex<double>>> ID = {{1.,0},{0,1.}};

	complex<double> lambda = b0*alphas_muR*log(N);
	complex<double> evol_factor = ISNLL*log(1.-2.*lambda)/(2.*M_PI*b0);
	resum[0][0] = exp(ISNLL*evol_factor*M_PI/alphas_muR*eval[0]);
	resum[1][1] = exp(ISNLL*evol_factor*M_PI/alphas_muR*eval[1]);
	resumC[0][0] = exp(ISNLL*evol_factor*M_PI/alphas_muR*conj(eval[0]));
	resumC[1][1] = exp(ISNLL*evol_factor*M_PI/alphas_muR*conj(eval[1]));
	resumE[0][0] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(eval[0]));
	resumE[1][1] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(eval[1]));
	resumEC[0][0] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(conj(eval[0])));
	resumEC[1][1] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(conj(eval[1])));

	if((real(eval[0])==0.)||(real(eval[1])==0.)){//works
			complex<double> omega = omega_qq(H_IJ,resumC,S_IJ,resum); //trace(S*Ubar*H*U)
			complex<double> result = omega*exp(2.*(1./alphas_muR*ISLL*g1_M2(A1q,lambda)+ISNLL*g2_M2(A1q,A2q,lambda)));
			if(expansion){
				complex<double> omegaE = omega_qq(H_IJ,resumEC,S_IJ,resumE);
				complex<double> omega0 = omega_qq(H_IJ,ID,S_IJ,ID);
				return result - (omegaE + omega0*delidelj_exp(N,A1q));
			}
			return result;
	}
	else{
		GammaR_IJ[0][0] = eval[0];
		GammaR_IJ[0][1] = 0;
		GammaR_IJ[1][0] = 0;
		GammaR_IJ[1][1] = eval[1];

		vector<vector<complex<double>>> R = R_qq(evec);
		vector<vector<complex<double>>> Rdag = Rdag_qq(R);
		vector<vector<complex<double>>> Rinv = Rinv_qq(evec);
		vector<vector<complex<double>>> Rinvdag = Rdag_qq(Rinv);

		HR_IJ = mult_qq(Rinv,H_IJ,Rinvdag);
		SR_IJ = mult_qq(Rdag,S_IJ,R);
		//GammaR_IJ = mult_qq(Rinv,Gamma_IJ,R);
		complex<double> omega = omega_qq(HR_IJ,resumC,SR_IJ,resum); //trace(S*H)
		complex<double> result = omega*exp(2.*(1./alphas_muR*ISLL*g1_M2(A1q,lambda)+ISNLL*g2_M2(A1q,A2q,lambda)));
		if(expansion){
			//cout << "exponent differences qq " << exp(2.*(1./alphas_muR*ISLL*g1(A1q,lambda)+ISNLL*g2(A1q,A2q,lambda))) << " " << 1.+delidelj_exp(N,A1q) << endl;
			complex<double> omegaE = omega_qq(HR_IJ,resumEC,SR_IJ,resumE);
			complex<double> omega0 = omega_qq(HR_IJ,ID,SR_IJ,ID);
			//cout << "omega differences qq " << omegaE << " " << omega0 << endl;
			//cout << "difference in result " << result << " " << (omegaE + omega0*delidelj_exp(N,A1q)) << endl;
			return result - (omegaE + omega0*delidelj_exp(N,A1q));
		}
		return result;
	}
	return 0.;

}

complex<double> gg_res(complex<double> N, complex<double> s, complex<double> t13, complex<double> t14, complex<double> t23, complex<double> t24, double s34){
	vector<vector<complex<double> > > HR_IJ = {{0,0,0},{0,0,0},{0,0,0}}; //hard matrix in R representatino
	vector<vector<complex<double> > > SR_IJ = {{0,0,0},{0,0,0},{0,0,0}}; //soft matrix in R representatino
	vector<vector<complex<double> > > H_IJ = {{0,0,0},{0,0,0},{0,0,0}}; //hard matrix
	vector<vector<complex<double> > > S_IJ = {{0,0,0},{0,0,0},{0,0,0}}; //soft matrix
	vector<vector<complex<double> > > Gamma_IJ = {{0,0,0},{0,0,0},{0,0,0}}; //anomalous dim
// filling hard function
	complex<double> H22 = H22_gg_c(s, t13, t14, t23, t24);
	complex<double> H32 = H23_gg_c(s, t13, t14, t23, t24);
	complex<double> H33 = H33_gg_c(s, t13, t14, t23, t24);
	H_IJ[0][0] = 1./pow(CA,2)*H22;
	H_IJ[0][1] = H22/CA;
	H_IJ[0][2] = 1./CA*H32;
	H_IJ[1][0] = H22/CA;
	H_IJ[1][1] = H22;
	H_IJ[1][2] = H32;
	H_IJ[2][0] = 1./CA*H32;
	H_IJ[2][1] = H32;
	H_IJ[2][2] = H33;
//filling soft function
  S_IJ[0][0] = S11_gg();
	S_IJ[0][1] = 0.;
	S_IJ[0][2] = 0.;
	S_IJ[1][0] = 0.;
	S_IJ[1][1] = S22_gg();
	S_IJ[1][2] = 0.;
	S_IJ[2][0] = 0.;
	S_IJ[2][1] = 0.;
	S_IJ[2][2] = S33_gg();

// filling anomalous dim
	complex<double> G11 = lambda_gg_11(s34);
	complex<double> G13 = lambda_gg_13(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s);
	complex<double> G22 = lambda_gg_22(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s, s34);
	complex<double> G23 = lambda_gg_23(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s);
	complex<double> G31 = lambda_gg_31(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s);
	complex<double> G32 = lambda_gg_32(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s);
	complex<double> G33 = lambda_gg_33(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s, s34);
	Gamma_IJ[0][0] = G11;
	Gamma_IJ[0][1] = 0.;
	Gamma_IJ[0][2] = G13;
	Gamma_IJ[1][0] = 0.;
	Gamma_IJ[1][1] = G22;
	Gamma_IJ[1][2] = G23;
	Gamma_IJ[2][0] = G31;
	Gamma_IJ[2][1] = G32;
	Gamma_IJ[2][2] = G33;

	vector<complex<double>> eval = eigenvalues_gg(G11, G13, G22, G23, G31, G32, G33);
	if(real(G11) == 0){eval[0] = G11;}
	vector<vector<complex<double>>> evec = eigenvectors_gg(eval, G11, G13, G22, G23, G31, G32, G33);

    vector<vector<complex<double>>> resum = {{0,0,0},{0,0,0},{0,0,0}};
	vector<vector<complex<double>>> resumE = {{0,0,0},{0,0,0},{0,0,0}};
	vector<vector<complex<double>>> resumC = {{0,0,0},{0,0,0},{0,0,0}};
	vector<vector<complex<double>>> resumEC = {{0,0,0},{0,0,0},{0,0,0}};
	vector<vector<complex<double>>> ID = {{1,0,0},{0,1,0},{0,0,1}};

	complex<double> lambda = b0*alphas_muR*log(N);
	complex<double> evol_factor = ISNLL*log(1.-2.*lambda)/(2.*M_PI*b0);
	resum[0][0] = exp(ISNLL*evol_factor*M_PI/alphas_muR*(eval[0]));
	resum[1][1] = exp(ISNLL*evol_factor*M_PI/alphas_muR*(eval[1]));
	resum[2][2] = exp(ISNLL*evol_factor*M_PI/alphas_muR*(eval[2]));
	resumC[0][0] = exp(ISNLL*evol_factor*M_PI/alphas_muR*(conj(eval[0])));
	resumC[1][1] = exp(ISNLL*evol_factor*M_PI/alphas_muR*(conj(eval[1])));
	resumC[2][2] = exp(ISNLL*evol_factor*M_PI/alphas_muR*(conj(eval[2])));
  //cout << eval[0] << " " << eval[1] << " " << eval[2] << endl;
	resumE[0][0] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(eval[0]));
	resumE[1][1] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(eval[1]));
	resumE[2][2] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(eval[2]));
	resumEC[0][0] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(conj(eval[0])));
	resumEC[1][1] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(conj(eval[1])));
	resumEC[2][2] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(conj(eval[2])));
	//cout << "SOFTgg " << resum[0][0] << " " << resum[1][1] << " " << resum[2][2] << endl;
	//cout << "SOFTggexp " << resumE[0][0] << " " << resumE[1][1] << " " << resumE[2][2] << endl;

	if(abs(real(eval[0])) == 0){//works
			complex<double> omega = omega_gg(H_IJ,resumC,S_IJ,resum); //trace(S*H)
			complex<double> result = omega*exp(2.*(1./alphas_muR*ISLL*g1_M2(A1g,lambda)+ISNLL*g2_M2(A1g,A2g,lambda)));
			if(expansion){
				vector<vector<complex<double>>> HtimesSE = mult_gg(H_IJ,S_IJ,resumE);
				vector<vector<complex<double>>> HtimesS0 = mult_gg(H_IJ,S_IJ,ID);
				complex<double> omegaE = omega_gg(H_IJ,resumEC,S_IJ,resumE);
				complex<double> omega0 = omega_gg(H_IJ,ID,S_IJ,ID);
				return result- (omegaE + omega0*delidelj_exp(N,A1g));
			}
			return result;
	}
	else{
		vector<vector<complex<double>>> R = R_gg(evec);
		vector<vector<complex<double>>> Rdag = Rdag_gg(R);
		vector<vector<complex<double>>> Rinv = Rinv_gg(evec);
		vector<vector<complex<double>>> Rinvdag = Rdag_gg(Rinv);
		HR_IJ = mult_gg(Rinv,H_IJ,Rinvdag);
		SR_IJ = mult_gg(Rdag,S_IJ,R);

		complex<double> omega = omega_gg(HR_IJ,resumC,SR_IJ,resum);
		complex<double> result =  omega*exp(2.*(1./alphas_muR*ISLL*g1_M2(A1g,lambda)+ISNLL*g2_M2(A1g,A2g,lambda)));
		if(expansion){
			//cout << "exponent differences gg " << exp(2.*(1./alphas_muR*ISLL*g1(A1g,lambda)+ISNLL*g2(A1g,A2g,lambda))) << " " << 1.+delidelj_exp(N,A1g) << endl;
			complex<double> omegaE = omega_gg(HR_IJ,resumEC,SR_IJ,resumE);
			complex<double> omega0 = omega_gg(HR_IJ,ID,SR_IJ,ID);
			//cout << "omega differences gg " << omegaE << " " << omega0 << endl;
			//cout << "difference in result " << result << " " << (omegaE + omega0*delidelj_exp(N,A1g)) << endl;
			return result- (omegaE + omega0*delidelj_exp(N,A1g));
		}
		return result;
	}
}
