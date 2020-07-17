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


vector<double> D1_qq() //normale array werkt niet, iets met geheugen.
{
	vector<double> D1(2);
	D1[0] = 0;
	D1[1] = -CA;
	return D1;
}

vector<std::complex<double> > log_delta_qq(std::complex<double> lambda)
{
	vector<std::complex<double> > h2_qq(D1_qq().size());
	complex<double> A = log(1.-2.*lambda)/(2.*M_PI*b0);
	for(int i = 0; i < D1_qq().size(); i++)
		h2_qq[i] = A*D1_qq()[i];
	return h2_qq;
}

vector<double> D1_gg()
{
	vector<double> D1(3);
	D1[0] = 0;
	D1[1] = -CA;
	D1[2] = -CA;
	return D1;
}

vector<std::complex<double> > log_delta_gg(std::complex<double> lambda)
{
	vector<complex<double> > h2_gg(D1_gg().size());
	complex<double> A = log(1.-2.*lambda)/(2.*M_PI*b0);
	for(int i = 0; i < D1_gg().size(); i++)
		h2_gg[i] = A*D1_gg()[i];
	return h2_gg;
}

complex<double> delidelj_exp(complex<double> N, double A1){
	return alphas_muR*2.*A1/M_PI*log(N)*(log(N)-ISNLL*log(Q2/muF2));
}

complex<double> SH_qq_res(complex<double> N, double s, double t13, double t14, double t23, double t24){
	complex<double> lambda = b0*alphas_muR*log(N);
	vector<complex<double>> wide_soft = log_delta_gg(lambda);
	return H22_qq(s, t13, t14, t23, t24)*S22_qq()*exp(ISNLL*wide_soft[1])
					*exp(2.*(1./alphas_muR*ISLL*g1(A1q,lambda)+ISNLL*g2(A1q,A2q,lambda)));

}
complex<double> SH_gg_res(complex<double> N, double s, double t13, double t14, double t23, double t24)
{
	complex<double> lambda = b0*alphas_muR*log(N);
	vector<complex<double>> wide_soft = log_delta_gg(lambda);
	return (H22_gg(s, t13, t14, t23, t24)*(S22_gg()*exp(ISNLL*wide_soft[1])+S11_gg()/(pow(CA,2))*exp(ISNLL*wide_soft[0]))
	        + exp(ISNLL*wide_soft[2])*H33_gg(s, t13, t14, t23, t24)*S33_gg())
					*exp(2.*(1./alphas_muR*ISLL*g1(A1g,lambda)+ISNLL*g2(A1g,A2g,lambda)));
}

complex<double> SH_qq_res_c(complex<double> N, complex<double>  s, complex<double> t13, complex<double> t14, complex<double> t23, complex<double> t24)
{
	complex<double> lambda = b0*alphas_muR*log(N);
	vector<complex<double>> wide_soft = log_delta_qq(lambda);
	return H22_qq_c(s, t13, t14, t23, t24)*S22_qq()*exp(ISNLL*wide_soft[1])
					*exp(2.*(1./alphas_muR*ISLL*g1(A1q,lambda)+ISNLL*g2(A1q,A2q,lambda)));

}

complex<double> SH_gg_res_c(complex<double> N, complex<double>  s, complex<double> t13, complex<double> t14, complex<double> t23, complex<double> t24)
{
	complex<double> lambda = b0*alphas_muR*log(N);
	vector<complex<double>> wide_soft = log_delta_gg(lambda);
	return (H22_gg_c(s, t13, t14, t23, t24)*(S22_gg()*exp(ISNLL*wide_soft[1])+S11_gg()/(pow(CA,2))*exp(ISNLL*wide_soft[0]))
	        + exp(ISNLL*wide_soft[2])*H33_gg_c(s, t13, t14, t23, t24)*S33_gg())
					*exp(2.*(1./alphas_muR*ISLL*g1(A1g,lambda)+ISNLL*g2(A1g,A2g,lambda)));
}

complex<double> pT_qq_res_abs(complex<double> N, double pT2, complex<double> s, complex<double> t13, complex<double> t14, complex<double> t23, complex<double> t24, double s34){
	complex<double> lambda = b0*alphas_muR*log(N);
	complex<double> wide_soft = ISNLL*CA*(log(1.+pT2/(4.*mt2))-1.)*log(1.-2.*lambda)/(2.*M_PI*b0);
	complex<double> H22qq = H22_qq_c(s, t13, t14, t23, t24);
	double S22qq = S22_qq();
	complex<double> result = H22qq*S22qq*exp(wide_soft)
														*exp(2.*(1./alphas_muR*ISLL*g1(A1q,lambda)+ISNLL*g2(A1q,A2q,lambda)));
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
														*exp(2.*(1./alphas_muR*ISLL*g1(A1g,lambda)+ISNLL*g2(A1g,A2g,lambda)));
	if(expansion){
		complex<double> wide_soft2 = 1.-ISNLL*CA*(log(1.+pT2/(4.*mt2))-1.)*alphas_muR*log(N)/(M_PI);
		//cout << "difference in result gg " << result << " " << ((H22gg*(S22gg*wide_soft2+S11gg/(pow(CA,2)))+ wide_soft2*H33gg*S33gg)
		//									+ (H22gg*(S22gg+S11gg/(pow(CA,2)))+ H33gg*S33gg)*delidelj_exp(N,A1g)) << endl;
		complex<double> exponent = delidelj_exp(N,A1g);
		return result - ((H22gg*(S22gg*(exponent+wide_soft2)+S11gg/(pow(CA,2))*(1.+exponent)))+ (exponent+wide_soft2)*H33gg*S33gg);
					}
	return result;
}

//checked the diagonalization, that works
complex<double> pT_qq_res(complex<double> N, complex<double> s, complex<double> t13, complex<double> t14, complex<double> t23, complex<double> t24, double s34){

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
	vector<vector<complex<double>>> resumE = {{0,0},{0,0}};
	vector<vector<complex<double>>> ID = {{1.,0},{0,1.}};

	complex<double> lambda = b0*alphas_muR*log(N);
	complex<double> evol_factor = ISNLL*log(1.-2.*lambda)/(2.*M_PI*b0);
	resum[0][0] = exp(ISNLL*evol_factor*M_PI/alphas_muR*(eval[0]+conj(eval[0])));
	resum[1][1] = exp(ISNLL*evol_factor*M_PI/alphas_muR*(eval[1]+conj(eval[1])));
	resumE[0][0] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(eval[0]+conj(eval[0])));
	resumE[1][1] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(eval[1]+conj(eval[1])));
	//cout << "SOFTqq " << resum[0][0] << " " << resum[1][1] << endl;
	//cout << "SOFTqqexp " << resumE[0][0] << " " << resumE[1][1] << endl;
	if((real(eval[0])==0.)||(real(eval[1])==0.)){//works
			vector<vector<complex<double>>> HtimesS = mult_qq(H_IJ,S_IJ,resum); //S*H
			complex<double> omega = HtimesS[0][0]+HtimesS[1][1]; //trace(S*H)
			complex<double> result = omega*exp(2.*(1./alphas_muR*ISLL*g1(A1q,lambda)+ISNLL*g2(A1q,A2q,lambda)));
			if(expansion){
				vector<vector<complex<double>>> HtimesSE = mult_qq(H_IJ,S_IJ,resumE); //S*H with S approximated
				vector<vector<complex<double>>> HtimesS0 = mult_qq(H_IJ,S_IJ,ID); //S*H with no exponent
				complex<double> omegaE = HtimesSE[0][0]+HtimesSE[1][1];
				complex<double> omega0 = HtimesS0[0][0]+HtimesS0[1][1];
				return /*result - */(omegaE + omega0*delidelj_exp(N,A1q));
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
		GammaR_IJ = mult_qq(Rinv,Gamma_IJ,R);
		vector<vector<complex<double>>> HtimesS = mult_qq(HR_IJ,SR_IJ,resum); //S*H in R representation

		complex<double> omega = HtimesS[0][0]+HtimesS[1][1]; //trace(S*H)
		complex<double> result = omega*exp(2.*(1./alphas_muR*ISLL*g1(A1q,lambda)+ISNLL*g2(A1q,A2q,lambda)));
		if(expansion){
			//cout << "exponent differences qq " << exp(2.*(1./alphas_muR*ISLL*g1(A1q,lambda)+ISNLL*g2(A1q,A2q,lambda))) << " " << 1.+delidelj_exp(N,A1q) << endl;
			vector<vector<complex<double>>> HtimesSE = mult_qq(H_IJ,S_IJ,resumE);
			vector<vector<complex<double>>> HtimesS0 = mult_qq(H_IJ,S_IJ,ID); //S*H with no exponent
			complex<double> omegaE = HtimesSE[0][0]+HtimesSE[1][1];
			complex<double> omega0 = HtimesS0[0][0]+HtimesS0[1][1];
			//cout << "omega differences qq " << omegaE << " " << omega0 << endl;
			//cout << "difference in result " << result << " " << (omegaE + omega0*delidelj_exp(N,A1q)) << endl;
			return result - (omegaE + omega0*delidelj_exp(N,A1q));
		}
		return result;
	}
	return 0.;

}

complex<double> pT_gg_res(complex<double> N, complex<double> s, complex<double> t13, complex<double> t14, complex<double> t23, complex<double> t24, double s34){
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
	vector<vector<complex<double>>> ID = {{1,0,0},{0,1,0},{0,0,1}};

	complex<double> lambda = b0*alphas_muR*log(N);
	complex<double> evol_factor = ISNLL*log(1.-2.*lambda)/(2.*M_PI*b0);
	resum[0][0] = exp(ISNLL*evol_factor*M_PI/alphas_muR*(eval[0]+conj(eval[0])));
	resum[1][1] = exp(ISNLL*evol_factor*M_PI/alphas_muR*(eval[1]+conj(eval[1])));
	resum[2][2] = exp(ISNLL*evol_factor*M_PI/alphas_muR*(eval[2]+conj(eval[2])));

	resumE[0][0] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(eval[0]+conj(eval[0])));
	resumE[1][1] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(eval[1]+conj(eval[1])));
	resumE[2][2] = 1. - ISNLL*alphas_muR/M_PI*log(N)*(M_PI/alphas_muR*(eval[2]+conj(eval[2])));
	//cout << "SOFTgg " << resum[0][0] << " " << resum[1][1] << " " << resum[2][2] << endl;
	//cout << "SOFTggexp " << resumE[0][0] << " " << resumE[1][1] << " " << resumE[2][2] << endl;

	if(abs(real(eval[0])) == 0){//works
			vector<vector<complex<double>>> HtimesS = mult_gg(H_IJ,S_IJ,resum); //S*H
			complex<double> omega = HtimesS[0][0]+HtimesS[1][1]+HtimesS[2][2]; //trace(S*H)
			complex<double> result = omega*exp(2.*(1./alphas_muR*ISLL*g1(A1g,lambda)+ISNLL*g2(A1g,A2g,lambda)));
			if(expansion){
				vector<vector<complex<double>>> HtimesSE = mult_gg(H_IJ,S_IJ,resumE);
				vector<vector<complex<double>>> HtimesS0 = mult_gg(H_IJ,S_IJ,ID);
				complex<double> omegaE = HtimesSE[0][0]+HtimesSE[1][1]+HtimesSE[2][2];
				complex<double> omega0 = HtimesS0[0][0]+HtimesS0[1][1]+HtimesS0[2][2];
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

		vector<vector<complex<double>>> HtimesS = mult_gg(HR_IJ,SR_IJ,resum); //S*H in R representation
    complex<double> omega = HtimesS[0][0]+HtimesS[1][1]+HtimesS[2][2]; //trace(S*H)
		complex<double> result =  omega*exp(2.*(1./alphas_muR*ISLL*g1(A1g,lambda)+ISNLL*g2(A1g,A2g,lambda)));
		if(expansion){
			//cout << "exponent differences gg " << exp(2.*(1./alphas_muR*ISLL*g1(A1g,lambda)+ISNLL*g2(A1g,A2g,lambda))) << " " << 1.+delidelj_exp(N,A1g) << endl;
			vector<vector<complex<double>>> HtimesSE = mult_gg(H_IJ,S_IJ,resumE);
			vector<vector<complex<double>>> HtimesS0 = mult_gg(H_IJ,S_IJ,ID);
			complex<double> omegaE = HtimesSE[0][0]+HtimesSE[1][1]+HtimesSE[2][2];
			complex<double> omega0 = HtimesS0[0][0]+HtimesS0[1][1]+HtimesS0[2][2];
			//cout << "omega differences gg " << omegaE << " " << omega0 << endl;
			//cout << "difference in result " << result << " " << (omegaE + omega0*delidelj_exp(N,A1g)) << endl;
			return result- (omegaE + omega0*delidelj_exp(N,A1g));
		}
		return result;
	}
}




complex<double> omega_qq_res(complex<double> N, complex<double> s, complex<double> t13, complex<double> t14, complex<double> t23, complex<double> t24, double s34)
{
	vector<vector<complex<double> > > HR_IJ(2);
	vector<vector<complex<double> > > SR_IJ(2);
	vector<vector<complex<double> > > Lambda(2);
	for ( int i = 0 ; i < 2 ; i++ ){HR_IJ[i].resize(2);SR_IJ[i].resize(2);Lambda[i].resize(2);}
	complex<double> Gamma[2][2];
	complex<double> lambdaqq11 = lambda_qq_11(s34), lambdaqq22 = lambda_qq_22(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s34, s), lambdaqq12 = lambda_qq_12(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s);
	complex<double> lambdaqq21 = lambda_qq_21(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s);
	Gamma[0][0] = lambdaqq11;
	Gamma[0][1] = lambdaqq12;
	Gamma[1][0] = lambdaqq21;
	Gamma[1][1] = lambdaqq22;


	complex<double> H[2][2]; //hard matrix
	H[0][0] = 0;
	H[0][1] = 0;
	H[1][0] = 0;
	H[1][1] = H22_qq_c(s, t13, t14, t23, t24);
	complex<double> S[2][2]; //soft matrix
	S[0][0] = S11_qq();
	S[0][1] = 0;
	S[1][0] = 0;
	S[1][1] = S22_qq();

	complex<double> eigen1 = eigen_qq_1(lambdaqq11, lambdaqq22, lambdaqq12, lambdaqq21);
	complex<double> eigen2 = eigen_qq_2(lambdaqq11, lambdaqq22, lambdaqq12, lambdaqq21);
	Lambda[0][0] =eigen1;
	Lambda[0][1] = 0;
	Lambda[1][0] = 0;
	Lambda[1][1] = eigen2;

	complex<double>  R[2][2], Rinv[2][2];
	complex<double> Rinv_const = Rinv_qq_const(lambdaqq11, lambdaqq22, lambdaqq12, lambdaqq21, eigen1, eigen2);
	complex<double> R11=Rinv_const*R_qq_1i(lambdaqq11, lambdaqq22, lambdaqq12, lambdaqq21, eigen1), R22=Rinv_const*1., R12=Rinv_const*R_qq_1i(lambdaqq11, lambdaqq22, lambdaqq12, lambdaqq21, eigen2), R21 = 0;

	Rinv[0][0] = R22;
	Rinv[0][1] = -R12;
	Rinv[1][0] = -R22;
	Rinv[1][1] = R11;

	R[0][0] = conj(R22);
	R[0][1] = conj(-R22);
	R[1][0] = conj(-R12);
	R[1][1] = conj(R11);

	complex<double> a = Rinv[0][0]*H[0][0]+Rinv[0][1]*H[1][0] , b =  Rinv[0][0]*H[0][1]+Rinv[0][1]*H[1][1]
			, c =  Rinv[1][0]*H[0][0]+Rinv[1][1]*H[1][0], d =  Rinv[1][0]*H[0][1]+Rinv[1][1]*H[1][1];

	HR_IJ[0][0] = a*R[0][0]+b*R[1][0];
	HR_IJ[0][1] = a*R[0][1]+b*R[1][1];
	HR_IJ[1][0] = c*R[0][0]+d*R[1][0];
	HR_IJ[1][1] = c*R[0][1]+d*R[1][1];


	complex<double> lambda1_2 = lambdaqq11-lambdaqq22;
	R11=R_qq_1i(lambdaqq11, lambdaqq22, lambdaqq12, lambdaqq21, eigen1), R22=1., R12=R_qq_1i(lambdaqq11, lambdaqq22, lambdaqq12, lambdaqq21, eigen2), R21=1.;
	R[0][0] = R11;
	R[0][1] = R12;
	R[1][0] = R21;
	R[1][1] = R22;
	Rinv[0][0] = conj(R11);
	Rinv[0][1] = conj(R21);
	Rinv[1][0] = conj(R12);
	Rinv[1][1] = conj(R22);
	a = Rinv[0][0]*S[0][0]+Rinv[0][1]*S[1][0] , b =  Rinv[0][0]*S[0][1]+Rinv[0][1]*S[1][1], c =  Rinv[1][0]*S[0][0]+Rinv[1][1]*S[1][0], d =  Rinv[1][0]*S[0][1]+Rinv[1][1]*S[1][1];

	SR_IJ[0][0] = a*R[0][0]+b*R[1][0];
	SR_IJ[0][1] = a*R[0][1]+b*R[1][1];
	SR_IJ[1][0] = c*R[0][0]+d*R[1][0];
	SR_IJ[1][1] = c*R[0][1]+d*R[1][1];

	complex<double> omega = 0;
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j ++){
				omega = HR_IJ[i][j]*SR_IJ[j][i]*exp(log(1.-2.*b0*alphas_muR*log(N))/(2.*M_PI*b0)*M_PI/alphas_muR*conj(Lambda[j][j]))*exp(log(1.-2.*b0*alphas_muR*log(N))/(2.*M_PI*b0)*M_PI/alphas_muR*(Lambda[i][i])) +omega;
			}
		}
	return omega*exp(2.*(1./alphas_muR*ISLL*g1(A1q,b0*alphas_muR*log(N))+ISNLL*g2(A1q,A2q,b0*alphas_muR*log(N))));

}



complex<double> omega_gg_res(complex<double> N, complex<double> s, complex<double> t13, complex<double> t14, complex<double> t23, complex<double> t24, double s34)
{
	vector<vector<complex<double> > > HR_IJ(3);
	vector<vector<complex<double> > > SR_IJ(3);
	vector<vector<complex<double> > > Lambda(3);
	for ( int i = 0 ; i < 3 ; i++ ){HR_IJ[i].resize(3);SR_IJ[i].resize(3);Lambda[i].resize(3);}

    complex<double> R[3][3], Rinv[3][3], Gamma[3][3], H[3][3], S[3][3];
	complex<double> H22 = H22_gg_c(s, t13, t14, t23, t24);
	complex<double> H32 = H23_gg_c(s, t13, t14, t23, t24);
	complex<double> H33 = H33_gg_c(s, t13, t14, t23, t24);

	H[0][0] = 1./pow(CA,2)*H22;
	H[0][1] = H22/CA;
	H[0][2] = 1./CA*H32;
	H[1][0] = H22/CA;
	H[1][1] = H22;
	H[1][2] = H32;
	H[2][0] = 1./CA*H32;
	H[2][1] = H32;
	H[2][2] = H33;

    S[0][0] = S11_gg();
	S[0][1] = 0.;
	S[0][2] = 0.;
	S[1][0] = 0.;
	S[1][1] = S22_gg();
	S[1][2] = 0.;
	S[2][0] = 0.;
	S[2][1] = 0.;
	S[2][2] = S33_gg();



	complex<double> lambda11 = lambda_gg_11(s34), lambda12 = 0;
	complex<double> lambda21 = 0, lambda22 = lambda_gg_22(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s, s34);
	complex<double> lambda31 = lambda_gg_31(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s);

	complex<double> a1 = a(lambda11, lambda22), b1 = b(lambda11,lambda22, lambda31), c1 = c(lambda11,lambda22,lambda31);
	complex<double> Z1 = collect(a1, b1, c1);
	complex<double> eigen1 = eigen1_new(a1, b1, Z1), eigen2 = eigen2_new(a1, b1, Z1), eigen3 = eigen3_new(a1, b1, Z1);

	Lambda[0][0] = eigen1; //to save time! a*R[0][0]+b*R[1][0]+c*R[2][0];
	Lambda[0][1] = 0;//a*R[0][1]+b*R[1][1]+c*R[2][1];
	Lambda[0][2] = 0;//a*R[0][2]+b*R[1][2]+c*R[2][2];
	Lambda[1][0] = 0;//d*R[0][0]+e*R[1][0]+f*R[2][0];
	Lambda[1][1] = eigen2; //d*R[0][1]+e*R[1][1]+f*R[2][1];
	Lambda[1][2] = 0;//d*R[0][2]+e*R[1][2]+f*R[2][2];
	Lambda[2][0] = 0;//g*R[0][0]+h*R[1][0]+i*R[2][0];
	Lambda[2][1] = 0;//g*R[0][1]+h*R[1][1]+i*R[2][1];
	Lambda[2][2] = eigen3; //g*R[0][2]+h*R[1][2]+i*R[2][2];

	complex<double> v1_1 = R_gg_1i(lambda31, eigen1, lambda11),  v2_1 = R_gg_1i(lambda31, eigen2, lambda11),  v3_1 = R_gg_1i(lambda31, eigen3, lambda11);
	complex<double> v1_2 = R_gg_2i(lambda31, eigen1, lambda22), v2_2 = R_gg_2i(lambda31, eigen2, lambda22), v3_2 = R_gg_2i(lambda31, eigen3, lambda22);
	complex<double> v1_3 = 1., v2_3 = 1., v3_3 = 1.;

	complex<double> Rinv_const = 1./Rinv_gg_const(v1_1, v2_1, v3_1, v1_2, v2_2, v3_2);
	Rinv[0][0] = Rinv_const*(v2_2-v3_2);
	Rinv[0][1] = Rinv_const*(v3_1-v2_1);
	Rinv[0][2] = Rinv_const*(v2_1*v3_2-v3_1*v2_2);
	Rinv[1][0] = Rinv_const*(v3_2-v1_2);
	Rinv[1][1] = Rinv_const*(v1_1-v3_1);
	Rinv[1][2] = Rinv_const*(v1_2*v3_1-v1_1*v3_2);
	Rinv[2][0] = Rinv_const*(v1_2-v2_2);
	Rinv[2][1] = Rinv_const*(v2_1-v1_1);
	Rinv[2][2] = Rinv_const*(v1_1*v2_2-v2_1*v1_2);


	R[0][0] = conj(Rinv[0][0]);
	R[0][1] = conj(Rinv[1][0]);
	R[0][2] = conj(Rinv[2][0]);
	R[1][0] = conj(Rinv[0][1]);
	R[1][1] = conj(Rinv[1][1]);
	R[1][2] = conj(Rinv[2][1]);
	R[2][0] = conj(Rinv[0][2]);
	R[2][1] = conj(Rinv[1][2]);
	R[2][2] = conj(Rinv[2][2]);

	complex<double> a = Rinv[0][0]*H[0][0]+Rinv[0][1]*H[1][0]+Rinv[0][2]*H[2][0],
				    b = Rinv[0][0]*H[0][1]+Rinv[0][1]*H[1][1]+Rinv[0][2]*H[2][1],
				    c = Rinv[0][0]*H[0][2]+Rinv[0][1]*H[1][2]+Rinv[0][2]*H[2][2],
				    d = Rinv[1][0]*H[0][0]+Rinv[1][1]*H[1][0]+Rinv[1][2]*H[2][0],
					e = Rinv[1][0]*H[0][1]+Rinv[1][1]*H[1][1]+Rinv[1][2]*H[2][1],
					f = Rinv[1][0]*H[0][2]+Rinv[1][1]*H[1][2]+Rinv[1][2]*H[2][2],
					g = Rinv[2][0]*H[0][0]+Rinv[2][1]*H[1][0]+Rinv[2][2]*H[2][0],
					h = Rinv[2][0]*H[0][1]+Rinv[2][1]*H[1][1]+Rinv[2][2]*H[2][1],
					i = Rinv[2][0]*H[0][2]+Rinv[2][1]*H[1][2]+Rinv[2][2]*H[2][2];

	HR_IJ[0][0] = a*R[0][0]+b*R[1][0]+c*R[2][0];
	HR_IJ[0][1] = a*R[0][1]+b*R[1][1]+c*R[2][1];
	HR_IJ[0][2] = a*R[0][2]+b*R[1][2]+c*R[2][2];
	HR_IJ[1][0] = d*R[0][0]+e*R[1][0]+f*R[2][0];
	HR_IJ[1][1] = d*R[0][1]+e*R[1][1]+f*R[2][1];
	HR_IJ[1][2] = d*R[0][2]+e*R[1][2]+f*R[2][2];
	HR_IJ[2][0] = g*R[0][0]+h*R[1][0]+i*R[2][0];
	HR_IJ[2][1] = g*R[0][1]+h*R[1][1]+i*R[2][1];
	HR_IJ[2][2] = g*R[0][2]+h*R[1][2]+i*R[2][2];

	R[0][0] = v1_1;
	R[0][1] = v2_1;
	R[0][2] = v3_1;
	R[1][0] = v1_2;
	R[1][1] = v2_2;
	R[1][2] = v3_2;
	R[2][0] = v1_3;
	R[2][1] = v2_3;
	R[2][2] = v3_3;

	Rinv[0][0] = conj(R[0][0]);
	Rinv[0][1] = conj(R[1][0]);
	Rinv[0][2] = conj(R[2][0]);
	Rinv[1][0] = conj(R[0][1]);
	Rinv[1][1] = conj(R[1][1]);
	Rinv[1][2] = conj(R[2][1]);
	Rinv[2][0] = conj(R[0][2]);
	Rinv[2][1] = conj(R[1][2]);
	Rinv[2][2] = conj(R[2][2]);

    a = Rinv[0][0]*S[0][0]+Rinv[0][1]*S[1][0]+Rinv[0][2]*S[2][0],
			b = Rinv[0][0]*S[0][1]+Rinv[0][1]*S[1][1]+Rinv[0][2]*S[2][1],
			c = Rinv[0][0]*S[0][2]+Rinv[0][1]*S[1][2]+Rinv[0][2]*S[2][2],
			d = Rinv[1][0]*S[0][0]+Rinv[1][1]*S[1][0]+Rinv[1][2]*S[2][0],
			e = Rinv[1][0]*S[0][1]+Rinv[1][1]*S[1][1]+Rinv[1][2]*S[2][1],
			f = Rinv[1][0]*S[0][2]+Rinv[1][1]*S[1][2]+Rinv[1][2]*S[2][2],
			g = Rinv[2][0]*S[0][0]+Rinv[2][1]*S[1][0]+Rinv[2][2]*S[2][0],
			h = Rinv[2][0]*S[0][1]+Rinv[2][1]*S[1][1]+Rinv[2][2]*S[2][1],
			i = Rinv[2][0]*S[0][2]+Rinv[2][1]*S[1][2]+Rinv[2][2]*S[2][2];

	SR_IJ[0][0] = a*R[0][0]+b*R[1][0]+c*R[2][0];
	SR_IJ[0][1] = a*R[0][1]+b*R[1][1]+c*R[2][1];
	SR_IJ[0][2] = a*R[0][2]+b*R[1][2]+c*R[2][2];
	SR_IJ[1][0] = d*R[0][0]+e*R[1][0]+f*R[2][0];
	SR_IJ[1][1] = d*R[0][1]+e*R[1][1]+f*R[2][1];
	SR_IJ[1][2] = d*R[0][2]+e*R[1][2]+f*R[2][2];
	SR_IJ[2][0] = g*R[0][0]+h*R[1][0]+i*R[2][0];
	SR_IJ[2][1] = g*R[0][1]+h*R[1][1]+i*R[2][1];
	SR_IJ[2][2] = g*R[0][2]+h*R[1][2]+i*R[2][2];

	complex<double> omega = 0;
	for (int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j ++){
				omega = HR_IJ[i][j]*SR_IJ[j][i]*exp(log(1.-2.*b0*alphas_muR*log(N))/(2.*M_PI*b0)*M_PI/alphas_muR*conj(Lambda[j][j]))*exp(log(1.-2.*b0*alphas_muR*log(N))/(2.*M_PI*b0)*M_PI/alphas_muR*(Lambda[i][i])) +omega;
			}
		}
	return omega*exp(2.*(1./alphas_muR*ISLL*g1(A1g,b0*alphas_muR*log(N))+ISNLL*g2(A1g,A2g,b0*alphas_muR*log(N))));

}
