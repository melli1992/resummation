#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include "parameters.h"
#include "resum_functions.h"
#include "polygamma.hpp"

//////////////////////////////////////////////////////////
///
/// contains all the resummation functions up to NNLL
/// note that no final state resummation effects are here
/// so code is only for processes with IS partons
/// matching for DY and Higgs is also included
///
//////////////////////////////////////////////////////////
using namespace std;

// LP LL function h0 (or g1) hep-ph/0306211 eqn 39 (note that gammaE is not part there of lambda and the factor 2 (as they have 2*h0 = g1!) or 1905.11771 eqn 6
// (checked with mathematica), seems that the I*0.0 is another side of the cut taken by mathematica and by c++ (maybe this is bad, have to find out)
complex<double> g1(double A1,complex<double>lambda){
	return A1/(2.*M_PI*pow(b0,2))*(2.*lambda+(1.-2.*lambda)*log(1.-2.*lambda));
}
/////////////////////////////////////////////////////////////
// NLP LL function
complex<double> h1NLP(double A1, complex<double> N, complex<double>lambda){
	return -A1/(2.*M_PI*b0)*(log(1.-2.*lambda))/N;
}

/////////////////////////////////////////////////////////////
// LP NLL function h1 (or g2) hep-ph/0306211 (note factor 2 and gammaE) eqn 40 or 1905.11771 eqn 61
// (checked with mathematica)
complex<double> g2(double A1,double A2,complex<double>lambda){
	double INCeuler = 0.;
	if(INCEULER == 0) INCeuler = 1.;
	return 1./(2.*M_PI*b0)*(-A2/(M_PI*b0)+A1*log(Q2/muR2))*(2.*lambda+log(1.-2.*lambda))
	+ A1*b1/(2.*M_PI*pow(b0,3))*(2.*lambda+log(1.-2.*lambda)+1./2.*pow(log(1.-2.*lambda),2))
	- A1/(M_PI*b0)*lambda*log(Q2/muF2)
	- INCeuler*2*M_gammaE*log(1.-2.*lambda)*A1/(2.*M_PI*b0);
}

////////////////////////////////////////////////////////////////////////
// LP NNLL function g3 hep-ph/0306211 eqn 41 (checked with mathematica)
// checked with Richards code as well (for Nbar => N)
complex<double> g3(double A1,double A2,double A3,complex<double>lambda){
	double INCeuler = 0.;
	if(INCEULER == 0) INCeuler = 1.;
	return 2.*A1/M_PI*zeta2*lambda/(1.-2.*lambda)
	+A1*pow(b1,2)/(2.*M_PI*pow(b0,4)*(1.-2.*lambda))*(2.*pow(lambda,2)+2.*lambda*log(1.-2.*lambda)+1./2.*pow(log(1.-2.*lambda),2))
	+A1*b2/(2.*M_PI*pow(b0,3))*(2.*lambda+log(1.-2.*lambda)+2.*pow(lambda,2)/(1.-2.*lambda))
	+A3/(pow(M_PI,3)*pow(b0,2))*pow(lambda,2)/(1.-2.*lambda)
	-A2*b1/(2.*pow(M_PI,2)*pow(b0,3))*(1./(1.-2.*lambda))*(2.*lambda+log(1.-2.*lambda)+2.*pow(lambda,2))
	-A2/(pow(M_PI,2)*b0)*lambda*log(Q2/muF2)
	-A1/(2.*M_PI)*lambda*pow(log(Q2/muF2),2)
	+A1/M_PI*lambda*log(Q2/muR2)*log(Q2/muF2)
	+1./(1.-2.*lambda)*(A1*b1/(2.*M_PI*pow(b0,2))*(2.*lambda+log(1.-2.*lambda))
	-2*A2/(pow(M_PI,2)*b0)*pow(lambda,2))*log(Q2/muR2)
	+A1/M_PI*pow(lambda,2)/(1.-2.*lambda)*pow(log(Q2/muR2),2)
	+ b0*((A1*b1/(M_PI*pow(b0,3))*INCeuler*M_gammaE)*(-2.*lambda/(1.-2.*lambda)-log(1.-2.*lambda)/(1.-2.*lambda))
- (A1/(M_PI*b0)*(INCeuler*M_gammaE*M_gammaE)+A2/(M_PI*M_PI*pow(b0,2))*INCeuler*M_gammaE)*-2.*lambda/(1.-2.*lambda)
+ ((A1/(M_PI*b0)*INCeuler*M_gammaE)*(-2.*lambda/(1.-2.*lambda)))*log(Q2/muR2));
}

// wide angle contribution
complex<double> wideangle(double D2,complex<double>lambda){
	return -D2/(pow(M_PI,2)*b0)*lambda/(1.-2.*lambda);
}

///////////////////////////////////////////////////////////////////
// matching functions to NLO and NNLO (checked with mathematica)

//https://arxiv.org/pdf/1009.5691.pdf (appendix D13 en 14)
complex<double> NLOmatch(complex<double> N, double plusone, double A1, double g0){
	return 1.
					+ g0
					+ 2.*alphas_muR/M_PI
						*(A1*log(N*exp(M_gammaE*INCEULER)+plusone)
								*(ISNLP/N+ISLL*log(N*exp(M_gammaE*INCEULER)+plusone)-ISNLL*log(Q2/muF2)));
}
complex<double> NNLOmatch(complex<double> N, double plusone, double A1, double A2, double D2, double g0){
	return 1.
				+ g0
				+ 2.*alphas_muR/(M_PI)
						*(A1*log(N*exp(M_gammaE*INCEULER)+plusone)
								*(ISNLP/N+ISLL*log(N*exp(M_gammaE*INCEULER)+plusone)-ISNLL*log(Q2/muF2)))
				+ pow(alphas_muR,2)/(3.*pow(M_PI,2))*log(N*exp(M_gammaE*INCEULER)+plusone)
						*(6.*pow(A1,2)*log(N*exp(M_gammaE*INCEULER)+plusone)
								*pow(ISNLP/N + ISLL*log(N*exp(M_gammaE*INCEULER)+plusone) - ISNLL*log(Q2/muF2),2)
							+ -3.*D2*ISNNLL+6.*A2*ISNLL*log(N*exp(M_gammaE*INCEULER)+plusone)
							+ 6.*ISNLP*A1*b0*log(N*exp(M_gammaE*INCEULER)+plusone)*M_PI/N
							+ 4.*A1*b0*ISLL*pow(log(N*exp(M_gammaE*INCEULER)+plusone),2)*M_PI
							+ 2.*A1*b0*ISNNLL*pow(M_PI,3)
							- 3.*ISNNLL*A1*b0*M_PI*pow(log(Q2/muF2),2)
							- 6.*ISNLL*A1*b0*log(N*exp(M_gammaE*INCEULER)+plusone)*M_PI*log(Q2/muR2)
							- 6.*ISNNLL*log(Q2/muF2)*(A2 - A1*b0*M_PI*log(Q2/muR2)))
				;
}

/////////////////////
/// hard functions
/////////////////////

//DY https://arxiv.org/pdf/1009.5691.pdf (appendix D)
//note that their beta0 is our b0, but beta1 = b1/b0
double DY_g01(){
	double INCeuler = 0.;
	if(INCEULER == 0) INCeuler = 1.;
	return alphas_muR*CF/M_PI*(4.*zeta2-4.+2.*INCeuler*pow(M_gammaE,2)+(3./2.-2.*INCeuler*M_gammaE)*log(Q2/muF2));
}
// not needed for NNLL
double DY_g02(){
	double INCeuler = 0.;
	if(INCEULER == 0) INCeuler = 1.;
	return pow(alphas_muR,2)*(CF/(16.*pow(M_PI,2))*(CF*(511./4.-198.*zeta2-60.*zeta3+552./5.*pow(zeta2,2)+INCeuler*(-128.*pow(M_gammaE,2)*(1.-zeta2)+32.*pow(M_gammaE,4)))
								+CA*(-1535./12.+376./3.*zeta2+604./9.*zeta3-92./5.*pow(zeta2,2)+INCeuler*(1616./27.*M_gammaE-56.*M_gammaE*zeta3+536./9.*pow(M_gammaE,2)-16.*pow(M_gammaE,2)*zeta2+176./9.*pow(M_gammaE,3)))
								+nF*(127./6.-64./3.*zeta2+8./9.*zeta3+INCeuler*(-224./27.*M_gammaE-80./9.*pow(M_gammaE,2)-32./9.*pow(M_gammaE,3)))
								+pow(log(Q2/muF2),2)*(CF*(INCeuler*(32.*pow(M_gammaE,2)-48.*M_gammaE)+18.)+CA*(44./3.*INCeuler*M_gammaE-11.)+nF*(2.-8./3.*INCeuler*M_gammaE))
								+log(Q2/muF2)*(CF*(48.*zeta3+72.*zeta2-93.+INCeuler*(-128.*zeta2*M_gammaE+128.*M_gammaE+48.*pow(M_gammaE,2)-64.*pow(M_gammaE,3)))
												+CA*(193./3.-24.*zeta3-88./3.*zeta2+INCeuler*(16.*M_gammaE*zeta2-536./9.*M_gammaE-88./3.*pow(M_gammaE,2)))
												+nF*(16./3.*zeta2-34./3.+INCeuler*(80./9.*M_gammaE+16./3.*pow(M_gammaE,2)))))
								-b0*CF/(M_PI)*(4.*zeta2-4.+2.*INCeuler*pow(M_gammaE,2)+(3./2.-2.*INCeuler*M_gammaE)*log(Q2/muF2))*log(muF2/muR2));
}
// higgs
// see https://arxiv.org/pdf/hep-ph/0306211.pdf, eqn 44 and
// or https://arxiv.org/pdf/1206.4133.pdf eqn 2
double higgs_g01(){
	double INCeuler = 0.;
	if(INCEULER == 0) INCeuler = 1.;
	return alphas_muR/M_PI*CA/3.*(
					11./2.+6.*zeta2+(33.-2.*nF)/6.*log(muR2/muF2)
					+6.*pow(INCeuler*M_gammaE,2)+pow(M_PI,2)-6.*INCeuler*M_gammaE*log(Q2/muF2));
}
//dihiggs: https://arxiv.org/pdf/1807.03704.pdf eqn 11
double Cgg1_dihiggs(double Q2)
{
	complex<double> CLO = 3.*mH2/(Q2-mH2+I*mH*GammaH)-1.;
	double sigma1fin_sigma0 = 1./norm(CLO)*(11.*norm(CLO)+4./3.*real(CLO));
	return CA*pow(M_PI,2)*4./3.+sigma1fin_sigma0;
}
//https://arxiv.org/pdf/1807.03704.pdf eqn 11
// see also  https://arxiv.org/pdf/1505.07122.pdf eqn 14-16
// and https://arxiv.org/pdf/1305.5206.pdf eqn 12, 13 for I2, V2
// and https://arxiv.org/pdf/1408.2422.pdf for R2, F2 (full scale dependence, checked that F2 is the same, R2 cannot be checked)
// check whether this is now fully correct
/*double Cgg2_dihiggs(double Q2, double ctheta)
{
	double t = -1./2.*(Q2-2.*mH2-sqrt(Q2*(Q2-4.*mH2))*ctheta);
	double u = -1./2.*(Q2-2.*mH2+sqrt(Q2*(Q2-4.*mH2))*ctheta);
	complex<double> CLO = 3.*mH2/(Q2-mH2+I*mH*GammaH)-1.;
	double sigma1fin_sigma0 = 1./norm(CLO)*(11.*norm(CLO)+4./3.*real(CLO));
	double tplus = -1./2.*(Q2-2.*mH2-Q*sqrt(Q2-4.*mH2));
	double tmin = -1./2.*(Q2-2.*mH2+Q*sqrt(Q2-4.*mH2));
	double Lm = log(muR2/mt2);
	double Ls = log(muR2/Q2);
	double Lu = log(muR2/(-u));
	double Lt = log(muR2/(-t));
	double V2 = 1./pow(3.*Q2*t*u,2)*(pow(mH2,4)*pow(t+u,2)-2.*mH2*mH2*t*u*pow(t+u,2)+pow(t*u,2)*(4.*Q2+pow(t+u,2)));
	double I2 = 4.*M_PI*(1.+2.*pow(mH2,2)/pow(Q2,2)*log((mH2-t)*(mH2-u)/(t*u));
	double F2 = pow(CA,2)*(23827./648.-83./6.*zeta2-253./36.*zeta3+5./8.*zeta4+7./2.*Lm+89./3.*Ls+121./12.*pow(Ls,2))
				+9.*pow(CF,2)+CA*CF*(-145./6.-11./2.*Lm-11.*Ls)+pow(nF,2)*pow(TF,2)*(4./3.*pow(Ls,2)-22./9.*zeta2)
				-5./24.*CA-1./3.*CF-1./3.*nF*TF*CA*(2255./54.+40.*Ls+22.*pow(Ls,2)-217./6.*zeta2+49./3.*zeta3)
				-1./3.*nF*TF*CF*(41.-12.*Lm-24.*zeta3);
	double R2 = -7.*pow(CA,2)+11.*CA*CF-8.*nF*CF*TF+1./3.*CA*(476./9.+11./3.*(4.*Ls+Lt+Lu+4.*mH2/Q2)
			-8.*CF-4./9.*TF*nF*(10./3.+4.*Ls+Lt+Lu))
			-CA/3.*(1.+2.*pow(mH2,2)/pow(Q2,2))*(2.*li2(1.-pow(mH2,2)/(t*u))+4.*li2(mH2/t)+4.*li2(mH2/u)+4.*log(1.-mH2/t)*log(-mH2/t)
					+4.*log(1.-mH2/u)*log(-mH2/u)-8.*zeta2-pow(log(t/u),2));
	double sigma2fin_sigma0 = 1./norm(CLO)*1./(tplus-tmin)*(norm(CLO)*F2+real(CLO)*R2+im(CLO)*I2+V2);
	return pow(CA,2)*(-55./36.*zeta3+607./81.+67./16.*pow(M_PI,2)+91./144.*pow(M_PI,4))+CA*nF*(5.*zeta3/18.-82./81-5./8.*pow(M_PI,2))+pow(b0,2)*11./3.*pow(M_PI,4)+CA*sigmafin1/sigma0*(4./3.*pow(M_PI,2))+sigmafin2/sigma0;
}*/



// expansions of the functions
/// direct N-space result
complex<double> LP_LL_function_expanded_LP(complex<double> N, double Col_Fac){
 return	alphas_muR/M_PI*Col_Fac*(2.*pow(log(1./N),2.) - 4.*M_gammaE*log(1./N)
					+pow(M_PI,2.)/3. +2.*pow(M_gammaE,2)
				);
}
complex<double> LP_LL_function_expanded_NLP(complex<double> N, double Col_Fac){
 return	alphas_muR/M_PI*Col_Fac*(- (2.*(-log(1./N) + M_gammaE + 1.))/N
          + 2.*pow(log(1./N),2.) - 4.*M_gammaE*log(1./N)
					+pow(M_PI,2.)/3. +2.*pow(M_gammaE,2)
				);
}
complex<double> NLP_LL_function_expanded_NLP(complex<double> N, double Col_Fac){
return alphas_muR/M_PI*Col_Fac*2./N*(
				2.*log(N)+ 2.*M_gammaE
				);
}
complex<double> LP_LL_function_expanded(complex<double> N, double Col_Fac){
 return	alphas_muR/M_PI*Col_Fac*(
		 			((1./33.)*log(1./N) - M_gammaE/33. - 23./15120.)/pow(N,10.)
	 				+ 7./(120.*pow(N,9.))
					+(-((1./60.)*log(1./N)) + M_gammaE/60. + 221./151200)/pow(N,8.)
					- 5./(126.*pow(N,7.))
					+ (40.*log(1./N) - 40.*M_gammaE - 7.)/(2520.*pow(N,6.))
					+ 1./(20.*pow(N,5.))
					+(-(12.*log(1./N)) + 12.*M_gammaE + 5.)/(360.*pow(N,4.))
					- 1./(6.*pow(N,3.))
					+(2.*log(1./N) - 2.*M_gammaE - 3.)/(6.*pow(N,2.))
					- (2.*(-log(1./N) + M_gammaE + 1.))/N
					+ 2.*pow(log(1./N),2.) - 4.*M_gammaE*log(1./N)
					+pow(M_PI,2.)/3. +2.*pow(M_gammaE,2)
				);
}
complex<double> NLP_LL_function_expanded(complex<double> N, double Col_Fac){
return alphas_muR/M_PI*Col_Fac*2./N*(
				- 1./(66.*pow(N,10))
				+ 1./(120.*pow(N,8))
				- 1./(126.*pow(N,6))
				+ 1./(60.*pow(N,4))
				- 1./(6.*pow(N,2))
				+ 1./N
				+ 2.*log(N)+ 2.*M_gammaE
				);
}
complex<double> LP_LL_function_full(complex<double> N, double Col_Fac){
return alphas_muR/M_PI*Col_Fac*1./3.*(
				6.*pGamma(0,N)*(pGamma(0,N)+2.*M_gammaE)-6.*pGamma(1,N)+pow(M_PI,2)+6.*pow(M_gammaE,2)
				);
}
complex<double> NLP_LL_function_full(complex<double> N, double Col_Fac){
return alphas_muR/M_PI*Col_Fac*2./N*(
				2.*pGamma(0,N+1.)+ 2.*M_gammaE
				);
}
