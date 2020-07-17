#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_psi.h>
#include "parameters.h"
#include "k_factors_dy.h"
#include "k_factors_higgs.h"
#include "SCET_functions.h"
#include "resum_functions.h"
#include "deriv_pdf.h"

//////////////////////////////////////////////////////////
///
/// contains all the SCET functions up to NNLL
///
//////////////////////////////////////////////////////////
using namespace std;

///////////////////////
/// solve Landau pole
///////////////////////
complex<double> invasLambdaQCD(complex<double> lambdaQCD){
	double b1 = beta1/beta0;
	double b2 = beta2/beta0;
	complex<double> LLambda = log((mZ2-I*1.E-16)/(lambdaQCD*lambdaQCD));
	return 4.*M_PI*(1./(beta0*LLambda)-1./(pow(beta0*LLambda,2))*b1*log(LLambda)
												+ISNNLL*1./(pow(beta0*LLambda,3))*(pow(b1,2)*(pow(log(LLambda),2)-log(LLambda)-1.)+b2))
								-pdfs[use_member]->alphasQ(sqrt(mZ2)) ;//+1./(pow(beta0*LLambda,4))*(pow(b1,3)*(-pow(log(LLambda),3)+5./2.*pow(log(LLambda),2)+2.*log(LLambda)-1./2.)-3.*b1*b2*log(LLambda)+b3/2.);
}

complex<double> DinvasLambdaQCD(complex<double> lambdaQCD){
	complex<double> delta = lambdaQCD*1.E-3;
	complex<double> fdp = invasLambdaQCD(lambdaQCD+delta);
	complex<double> fdm = invasLambdaQCD(lambdaQCD-delta);
	complex<double> fd2p = invasLambdaQCD(lambdaQCD+2.0*delta);
	complex<double> fd2m = invasLambdaQCD(lambdaQCD-2.0*delta);

	// derivative
	return (-fd2p+8.*fdp-8.*fdm+fd2m)/(12.*delta);
}

void solveLambdaQCD(){
	int i = 1;
	complex<double> lambda(0.18,0.);
	while((abs(invasLambdaQCD(lambda)) > 1.E-16) && (i <=20)){
		lambda = lambda - invasLambdaQCD(lambda)/DinvasLambdaQCD(lambda);
		i++;
	}
	//cout << lambda << endl;
	LambdaQCD = real(lambda);
}
///////////////////////
/// evolution as
///////////////////////

// https://arxiv.org/pdf/hep-ph/0408244.pdf eqn (2.5)
// Checked with Leonardo
complex<double> fourthKG(complex<double> DLR){
	complex<double> as = 0.;
	for(int i = 1; i < 20.; i++)
	 {
		complex<double> xk0 = DLR*betaF(alphas_muR/(4.*M_PI));
		complex<double> xk1 = DLR*betaF(alphas_muR/(4.*M_PI)+0.5*xk0);
		complex<double> xk2 = DLR*betaF(alphas_muR/(4.*M_PI)+0.5*xk1);
		complex<double> xk3 = DLR*betaF(alphas_muR/(4.*M_PI)+xk2);
		as = alphas_muR/(4.*M_PI) + 1./6.*(xk0+2.*xk1+2.*xk2+xk3);
	 }
	 return as;
}


complex<double> falphasQ2(complex<double> mu2){
	complex<double> as = 0.;
	complex<double> LLambda = log((mu2-I*1.E-16)/(LambdaQCD*LambdaQCD));
	double b1 = beta1/beta0;
	double b2 = beta2/beta0;
/*  complex<double> as = 1./(beta0*LLambda)-1./(pow(beta0*LLambda,2))*b1*log(LLambda)
												+1./(pow(beta0*LLambda,3))*(pow(b1,2)*(pow(log(LLambda),2)-log(LLambda)-1.)+b2)
												;//+1./(pow(beta0*LLambda,4))*(pow(b1,3)*(-pow(log(LLambda),3)+5./2.*pow(log(LLambda),2)+2.*log(LLambda)-1./2.)-3.*b1*b2*log(LLambda)+b3/2.);
*/
	if(ISNLL){
		as = 1./(beta0*LLambda)-1./(pow(beta0*LLambda,2))*b1*log(LLambda)
				+ISNNLL*1./(pow(beta0*LLambda,3))*(pow(b1,2)*(pow(log(LLambda),2)-log(LLambda)-1.)+b2);
		return 4.*M_PI*as;}
	else if(ISLL){
		LLambda = log((mu2-I*1.E-16)/(muR*muR));
	 	as = alphas_muR/(1.+alphas_muR*b0*LLambda);
		// https://arxiv.org/pdf/hep-ph/0408244.pdf eqn (2.5)
		return as;}
	return as;
/*
 complex<double> as = 0.;
 complex<double> LLambda = log((mu2-I*1.E-16)/(muR*muR));
 complex<double> as0 = alphas_muR/(4.*M_PI);
 complex<double> DLR = LLambda/20.;
 if(ISNLL){
	 as = fourthKG(DLR);
 }
 else{
	 as = as0/(1.+beta0*as0*LLambda);
 }
   return 4.*M_PI*as;*/
}

complex<double> falphasQ2_LO(complex<double> mu2){
	complex<double> LLambda = log((mu2-I*1.E-16)/(muR*muR));
	complex<double>	as = alphas_muR/(1.+alphas_muR*b0*LLambda);
	return as;
}

complex<double> betaF(complex<double> as){
	//return -alphas*alphas*(ISLL*b0+ISNLL*b1*alphas+ISNNLL*b2*alphas*alphas);
	return -as*as*(ISLL*beta0+ISNLL*beta1*as+ISNNLL*beta2*as*as);
}

///////////////////////
/// hard functions
///////////////////////

//https://arxiv.org/pdf/0809.4283.pdf, eqn (17). Kijk nog even goed naar de schalen die hier gebruikt worden
// Checked with Leonardo
vector<complex<double>> hard_higgs(complex<double> mu2, complex<double> muh2){
	complex<double> alphas_muh2 = falphasQ2(muh2);
	complex<double> logH(log(mu2/muh2));
	complex<double> c1top = CA*4.*11./6.;
	complex<double> c1L = CA*(-logH*logH + pow(M_PI,2)/6.);
	complex<double> c2L = pow(CA,2)*(pow(logH,4)/2. + 11./9.*pow(logH,3)+(-67./9.+pow(M_PI,2)/6.)*pow(logH,2)
																		+ (80./27.-11.*pow(M_PI,2)/9.-2.*zeta3)*logH
																		+ 5105./162.+67.*pow(M_PI,2)/36.+pow(M_PI,4)/72.-143./9.*zeta3)
												+CF*TF*nF*(4.*logH-67./3.+16.*zeta3)
												+CA*TF*nF*(-4./9.*pow(logH,3)+20./9.*pow(logH,2)+(104./27.+4.*pow(M_PI,2)/9.)*logH-1832./81.-5.*pow(M_PI,2)/9.-92./9.*zeta3);
	vector<complex<double>> hardH = {1.,ISNNLL*(alphas_muh2)/(4.*M_PI)*(c1top+c1L+conj(c1L)),pow((alphas_muh2)/(4.*M_PI),2)*(c1L*conj(c1L)+c2L+conj(c2L))};
	return hardH;
}

//https://arxiv.org/pdf/0809.4283.pdf, eqn (12), includes scale evolution
double cT_higgs(){
	complex<double> alphas_mut2 = falphasQ2((complex<double>)mut*mut);
	complex<double> alphas_mus2 = falphasQ2((complex<double>)mus*mus);
	complex<double> alphas_muf2 = falphasQ2((complex<double>)muF*muF);
	complex<double> Lt = log(mt*mt/(mut*mut));
	complex<double> c1L = 5.*CA-3.*CF;
	complex<double> c2L = (27./2.*CF*CF+(11.*Lt-100./3.)*CF*CA-(7.*Lt-1063./36.)*CA*CA
													-4./3.*CF*TF-5./6.*CA*TF-(8.*Lt+5.)*CF*TF*nF-47./9.*CA*TF*nF);
	//if(!ISNNLL){alphas_mus2 = alphas_muf2;}
	// cT at scale mT (O(1),O(alpha),O(alpha2))
	//complex<double> cT = 1.+(alphas_mut2)/(4.*M_PI)*(c1L+conj(c1L))+pow((alphas_mut2)/(4.*M_PI),2)*(c1L*conj(c1L)+c2L+conj(c2L));
	//Leonardo:
	complex<double> cT = (1.+(alphas_mut2)/(4.*M_PI)*(c1L)+pow((alphas_mut2)/(4.*M_PI),2)*(c2L));
	double cT2 = abs(cT*cT);
	// Is this consistent with counting??

	// beta(alphas)/alphas^2, the evolution factors
	complex<double> beta_ms_over_alphasms2 = -(beta0/(4.*M_PI) + ISNNLL*alphas_mus2/pow(4.*M_PI,2)*beta1);// + alphas_mus2*alphas_mus2/(64.*pow(M_PI,3))*beta2);
	complex<double> beta_mt_over_alphasmt2 = -(beta0/(4.*M_PI) + ISNNLL*alphas_mut2/pow(4.*M_PI,2)*beta1);// + alphas_mut2*alphas_mut2/(64.*pow(M_PI,3))*beta2);

	// evolution of cT
	// returns LL, NLL or NNLL (according to leonardo)
	double result = cT2*abs(pow(beta_ms_over_alphasms2/beta_mt_over_alphasmt2,2.)*pow(alphas_mus2/alphas_muf2,2*ISNLL));

	return result;
}

complex<double> cT_higgs_wo_heavy_top(){
	complex<double> alphas_mus2 = falphasQ2((complex<double>)mus*mus);
	complex<double> alphas_muf2 = falphasQ2((complex<double>)muF*muF);
	complex<double> alphas_muh2 = falphasQ2((complex<double>)muh*muh);
	complex<double> beta_ms_over_alphasms2 = -(beta0/(4.*M_PI) + ISNNLL*1./(4.*M_PI)*alphas_mus2/(4.*M_PI)*beta1);//+ ISNNLL*1./(4.*M_PI)*pow(alphas_mus2/(4.*M_PI),2)*beta2);
	complex<double> beta_mh_over_alphasmh2 = -(beta0/(4.*M_PI) + ISNNLL*1./(4.*M_PI)*alphas_muh2/(4.*M_PI)*beta1);//+ ISNNLL*1./(4.*M_PI)*pow(alphas_mus2/(4.*M_PI),2)*beta2);

	complex<double> result = pow(beta_ms_over_alphasms2/beta_mh_over_alphasmh2,2.*ISNLL)*pow(alphas_mus2/alphas_muf2,2*ISNLL);

	return result;
}

//https://arxiv.org/pdf/0710.0680.pdf, eqn (75) and (76), see also eqn (50) on how to combine
// Checked with Leonardo
vector<complex<double>> hard_DY(complex<double> mu2, complex<double> muh2){
	complex<double> alphas_muh2 = falphasQ2(muh2);
	complex<double> logH(log(mu2/muh2));
	complex<double> c1L = CF*(-logH*logH + 3.*logH-8.+pow(M_PI,2)/6.);
	complex<double> c2L = pow(CF,2)*(pow(logH,4)/2. - 3.*pow(logH,3)+(25./2.-pow(M_PI,2)/6.)*pow(logH,2)+ (-45./2.-3.*pow(M_PI,2)/2.+24.*zeta3)*logH
																		+ 255./8.+7.*pow(M_PI,2)/2.-83.*pow(M_PI,4)/360.-30.*zeta3)
												+CA*CF*(11./9.*pow(logH,3)+(-233./18.+pow(M_PI,2)/3.)*pow(logH,2)+(2545./54.+11.*pow(M_PI,2)/9.-26*zeta3)*logH
																										-51157./648.-337.*pow(M_PI,2)/108.+11.*pow(M_PI,4)/45.+313./9.*zeta3)
												+CF*TF*nF*(-4./9.*pow(logH,3)+38./9.*pow(logH,2)+(-418./27.-4.*pow(M_PI,2)/9.)*logH+4085./162.+23.*pow(M_PI,2)/27.+4./9.*zeta3);
	vector<complex<double>> hardDY = {1.,ISNNLL*(alphas_muh2)/(4.*M_PI)*(c1L+conj(c1L)),pow((alphas_muh2)/(4.*M_PI),2)*(c1L*conj(c1L)+c2L+conj(c2L))};
	return hardDY;
}

////////////////////
/// //acusp and S
/////////////////////
// Checked with Leonardo
vector<complex<double>> acusp(double nu, double mu, vector<double> gamma){
	complex<double> alphamu2(falphasQ2((complex<double>)mu*mu));
	complex<double> alphanu2(falphasQ2((complex<double>)nu*nu));
	complex<double> acusp_NLL(0.);
	complex<double> acusp_NNLL(0.);
	//complex<double> acusp_NNNLL(0.);
	complex<double> Las(log(alphamu2/alphanu2));
	acusp_NLL = ISNLL*gamma[0]/(2.*beta0)*Las;
	if((gamma[0] != 0)) acusp_NNLL = ISNNLL*gamma[0]/(2.*beta0)*((gamma[1]/gamma[0]-beta1/beta0)*(alphamu2-alphanu2)/(4.*M_PI));
	if((gamma[0] == 0)) acusp_NNLL = ISNNLL*1./(2.*beta0)*((gamma[1])*(alphamu2-alphanu2)/(4.*M_PI));
	//if(ISNNNLL) acusp_NNNLL = gamma[0]/(2.*beta0)*(gamma[2]/gamma[0]-beta2/beta0-beta1/beta0*(gamma[1]/gamma[0]-beta1/beta0))*(pow(alphamu2,2)-pow(alphanu2,2))/(32.*pow(M_PI,2)));

	return {0.,acusp_NLL,acusp_NNLL};
}

// Checked with Leonardo
vector<complex<double>> sudakov(double nu, double mu, vector<double> GammaC){
	complex<double> alphamu2(falphasQ2((complex<double>)mu*mu));
	complex<double> alphanu2(falphasQ2((complex<double>)nu*nu));
	complex<double> sud_LL(0.);
	complex<double> sud_NLL(0.);
	complex<double> sud_NNLL(0.);
	complex<double> r = alphamu2/alphanu2;
	sud_LL = ISLL*GammaC[0]/(4.*pow(beta0,2))*(4.*M_PI/alphanu2*(1.-1./r-log(r)));
	sud_NLL = ISNLL*GammaC[0]/(4.*pow(beta0,2))*((GammaC[1]/GammaC[0]-beta1/beta0)*(1.-r+log(r))+beta1/(2.*beta0)*pow(log(r),2));
	sud_NNLL = ISNNLL*GammaC[0]/(4.*pow(beta0,2))*alphanu2/(4.*M_PI)*((beta1*GammaC[1]/GammaC[0]/beta0-beta2/beta0)*(1.-r+r*log(r))
																															+(pow(beta1,2)/pow(beta0,2)-beta2/beta0)*(1.-r)*log(r)
																															-(pow(beta1,2)/pow(beta0,2)-beta2/beta0-beta1/beta0*GammaC[1]/GammaC[0]+GammaC[2]/GammaC[0])*pow(1.-r,2)/(2.));

	//cout << sud_LL << " " << sud_NLL << " " << sud_NNLL << endl;
	return {sud_LL,sud_NLL,sud_NNLL};
}

// Checked with Leonardo
vector<complex<double>> feta(double nu, double mu, vector<double> gamma)
{//nu=mus,mu=muF
		complex<double> alphamu2(falphasQ2((complex<double>)mu*mu));
		complex<double> alphanu2(falphasQ2((complex<double>)nu*nu));
		complex<double> eta_NLL(0.);
		complex<double> eta_NNLL(0.);
		//complex<double> acusp_NNNLL(0.);
		complex<double> Las(log(alphamu2/alphanu2));
		eta_NLL = 2.*gamma[0]/(2.*beta0)*Las;
		if((gamma[0] != 0)) eta_NNLL = 2.*ISNNLL*gamma[0]/(2.*beta0)*((gamma[1]/gamma[0]-beta1/beta0)*(alphamu2-alphanu2)/(4.*M_PI));
		if((gamma[0] == 0)) eta_NNLL = 2.*ISNNLL*1./(2.*beta0)*((gamma[1])*(alphamu2-alphanu2)/(4.*M_PI));

		return {eta_NLL,eta_NNLL};

}

//https://arxiv.org/pdf/0809.4283.pdf eqn. (28)
// L argument should be log(Q2/mus2)
// Checked with Leonardo
vector<double> fsoft_higgs(double L){ //contains vector with (non deriv, non deriv O(alpha2), derivative part, 2 derivatives)
	double alphas_mus2 = real(falphasQ2((complex<double>)mus*mus));
	//cout << "Higgs: as(mus2)=" << alphas_mus2<< " ,L0 " << CA*alphas_mus2/(4.*M_PI)*(pow(M_PI,2)/3.) << ", L2 " << CA*alphas_mus2/(4.*M_PI)*2.*pow(L,2) << endl;
	vector<double> results = {1.+ISNNLL*CA*alphas_mus2/(4.*M_PI)*(2.*pow(L,2)+pow(M_PI,2)/3.),
														2.*ISNNLL*CA*alphas_mus2/(4.*M_PI)*(2.*L),
														ISNNLL*CA*alphas_mus2/(4.*M_PI)*2.};
	return results;
}
// Checked with Leonardo
vector<double> fsoft_DY(double L){ //contains vector with (non deriv O(1), non deriv O(alpha2), derivative part, 2 derivatives)
	double alphas_mus2 = real(falphasQ2((complex<double>)mus*mus));
	//cout << "DY: as(mus2)=" << alphas_mus2<< ", L0=" << CF*alphas_mus2/(4.*M_PI)*(pow(M_PI,2)/3.) << ", L2=" << CF*alphas_mus2/(4.*M_PI)*2.*pow(L,2) << endl;
	vector<double> results = {1.+ISNNLL*CF*alphas_mus2/(4.*M_PI)*(2.*pow(L,2)+pow(M_PI,2)/3.),
														2.*ISNNLL*CF*alphas_mus2/(4.*M_PI)*(2.*L),
														ISNNLL*CF*alphas_mus2/(4.*M_PI)*2.};
	return results;
}


///////////////////
// eta derivatives
///////////////////

double Deta_z(double z, double eta){
	if(eta > -0.3)// derivative of pow(z,-BN*eta)*2.*eta*pow(1-z,-1+2*eta)*exp[-2*gam_E*eta]/Gamma[1+2eta]
		{return (-2.*pow(1. - z,-1. + 2.*eta)
			*(-1. + 2.*eta*M_gammaE - 2.*eta*log(1. - z)
				+ BN*eta*log(z) + 2.*eta*gsl_sf_psi_n(0,1. + 2.*eta)))/
   (exp(2.*eta*M_gammaE)*pow(z,BN*eta)*tgamma(1. + 2.*eta));}
  else if(eta > -0.9) // derivative of z^(-BN*eta)*(2*eta)*(1 + 2 eta)/Gamma[2 + 2 eta]*(1 - z)^(-1 + 2 eta)*Exp[-2*EulerGamma*eta]
 		{return ((-(1./tgamma(2.*eta + 2.)))*2.*pow(1. - z,2.*eta - 1.)*
      (2.*BN*eta*eta*log(z) + BN*eta*log(z) + 4.*M_gammaE*eta*eta -
         2.*(2.*eta + 1.)*eta*log(1. - z) + 2.*M_gammaE*eta - 4.*eta +
         2.*(2.*eta + 1.)*eta*gsl_sf_psi_n(0, 2.*eta + 2.) - 1.))/
   (exp(2.*M_gammaE*eta)*pow(z,BN*eta));}
	else{
		cout << "Eta is too large " << eta << " exiting program" << endl;
		exit(0);
		return 0;
	}
}
double DDeta_z(double z, double eta){ //w.o. z^-eta term
	// double derivative of pow(z,-BN*eta)*2.*eta*pow(1-z,-1+2*eta)*exp[-2*gam_E*eta]/Gamma[1+2eta]
	if(eta > -0.3)
	{return (2.*pow(1. - z,-1. + 2.*eta)
						*((2.*M_gammaE - 2.*log(1. - z) + BN*log(z))*(-2. + 2.*eta*M_gammaE - 2.*eta*log(1 - z)
									+ BN*eta*log(z))
									+ 4.*(-1. + 2.*eta*M_gammaE - 2.*eta*log(1. - z) + BN*eta*log(z))
											*gsl_sf_psi_n(0,1. + 2.*eta)
									+ 4.*eta*pow(gsl_sf_psi_n(0,1. + 2.*eta),2)
									- 4.*eta*gsl_sf_psi_n(1,1. + 2.*eta)))
					/(exp(2.*eta*M_gammaE)*pow(z,BN*eta)*tgamma(1. + 2.*eta));
				}
	else if(eta > -0.9) // derivative of z^(-BN*eta)*(2*eta)*(1 + 2 eta)/Gamma[2 + 2 eta]*(1 - z)^(-1 + 2 eta)*Exp[-2*EulerGamma*eta]
		{return ((1./tgamma(2.*eta + 2.))*2.*pow(1. - z,2.*eta - 1.)*
      (-(4.*log(1. - z)*(BN*eta*(2.*eta + 1.)*log(z) +
                 2.*eta*(2.*M_gammaE*eta + M_gammaE - 2.) - 1.)) +
         BN*log(z)*(BN*eta*(2.*eta + 1.)*log(z) +
              4.*eta*(2.*M_gammaE*eta + M_gammaE - 2.) - 2.) +
         4.*gsl_sf_psi_n(0, 2.*eta + 2.)*(eta*(2.*eta + 1.)*
                (BN*log(z) - 2.*log(1. - z)) + 2.*eta*(2.*M_gammaE*eta +
                   M_gammaE - 2.) - 1.) +
     4.*eta*(2.*eta + 1.)*pow(log(1. - z),2) +
         4.*M_gammaE*(eta*(2.*M_gammaE*eta + M_gammaE - 4.) - 1.) +
         4.*eta*(2.*eta + 1.)*pow(gsl_sf_psi_n(0, 2.*eta + 2.),2.) -
         4.*eta*(2.*eta + 1.)*gsl_sf_psi_n(1, 2.*eta + 2.) + 4.))/
   (exp(2.*M_gammaE*eta)*pow(z,BN*eta));}
	else{
		cout << "Eta is too large " << eta << " exiting program" << endl;
		exit(0);
		return 0;
	}
}

double Deta_z_eta(double z, double eta, double BNsub){
	// derivative of (1 + BN*eta)/(1 - z)^(1 - 2 eta)*2*eta*Exp[-2*M_gammaE*eta]/Gamma[1 + 2 eta]
 return -1.*(2.*pow(1. - z,-1. + 2.*eta)
 					*(-1. - 2.*BNsub*eta + 2.*eta*M_gammaE + 2.*BNsub*pow(eta,2)*M_gammaE
						- 2.*eta*(1. + BNsub*eta)*log(1. - z)
						+ 2.*eta*(1. + BNsub*eta)*gsl_sf_psi_n(0,1. + 2.*eta)))
				/(exp(2.*eta*M_gammaE)*tgamma(1. + 2.*eta));
}
double DDeta_z_eta(double z, double eta, double BNsub){
	// double derivative of (1 + BN*eta)/(1 - z)^(1 - 2 eta)*2*eta*Exp[-2*M_gammaE*eta]/Gamma[1 + 2 eta]
 return -1.*(-4.*pow(1. - z,-1. + 2.*eta)
 					*(BNsub - 2.*M_gammaE - 4.*BNsub*eta*M_gammaE + 2.*eta*pow(M_gammaE,2)
						+ 2.*BNsub*pow(eta,2)*pow(M_gammaE,2) + 2.*log(1. - z) + 4.*BNsub*eta*log(1. - z)
						- 4.*eta*M_gammaE*log(1. - z) - 4.*BNsub*pow(eta,2)*M_gammaE*log(1. - z)
						+ 2.*eta*pow(log(1. - z),2) + 2.*BNsub*pow(eta,2)*pow(log(1. - z),2)
						+ (-2.+4.*eta*M_gammaE+4.*BNsub*eta*(-1.+eta*M_gammaE)-4.*eta*(1.+BNsub*eta)*log(1. - z))
								*gsl_sf_psi_n(0,1 + 2*eta)
						+ 2.*eta*(1. + BNsub*eta)*pow(gsl_sf_psi_n(0,1. + 2.*eta),2)
						- 2.*eta*(1. + BNsub*eta)*gsl_sf_psi_n(1,1. + 2.*eta)))
			 		/(exp(2.*eta*M_gammaE)*tgamma(1. + 2.*eta));
}

double Deta_z_eta2(double z, double eta, double BNsub, double INCsqrtZ){
	// derivative of (-1 - BN*eta+1/2*INCSQRTZ)*(1 - z)^(2 eta)*2*eta*(1+2eta)Exp[-2*M_gammaE*eta]/Gamma[2 + 2 eta]
 return ((1./tgamma(2.*eta + 2.))*pow(1. - z,2.*eta)*(8.*M_gammaE*BNsub*pow(eta,3.) -
         12.*BNsub*pow(eta,2) + 4.*M_gammaE*BNsub*pow(eta,2) -
         2.*(2.*eta + 1.)*eta*log(1. - z)*(2.*BNsub*eta - INCsqrtZ + 2.) +
         2.*(2.*eta + 1.)*eta*gsl_sf_psi_n(0, 2.*eta + 2.)*
           (2.*BNsub*eta - INCsqrtZ + 2.) - 4.*BNsub*eta -
         4.*M_gammaE*pow(eta,2)*INCsqrtZ + 8.*M_gammaE*pow(eta,2) -
         2.*M_gammaE*eta*INCsqrtZ + 4.*eta*INCsqrtZ +
     4.*M_gammaE*eta -  8.*eta + INCsqrtZ - 2.))/exp(2.*M_gammaE*eta);
}
double DDeta_z_eta2(double z, double eta, double BNsub, double INCsqrtZ){
// double derivative of (-1 - BN*eta+1/2*INCSQRTZ)*(1 - z)^(2 eta)*2*eta*(1+2eta)Exp[-2*M_gammaE*eta]/Gamma[2 + 2 eta]
 return ((1./tgamma(2.*eta + 2.))*4.*pow(1. - z,2.*eta)*
      ((2.*eta*(2.*M_gammaE*eta + M_gammaE - 2.) - 1.)*
           (log(1. - z)*(2.*BNsub*eta - INCsqrtZ + 2.) +
        gsl_sf_psi_n(0, 2.*eta + 2.)*(-(2.*BNsub*eta) + INCsqrtZ - 2.) + BNsub)
				- eta*(2.*eta + 1.)*
           (log(1. - z)*(log(1. - z)*(2.*BNsub*eta - INCsqrtZ + 2.) +  2.*BNsub)
			- 2*gsl_sf_psi_n(0, 2.*eta + 2.)*
                (log(1. - z)*(2.*BNsub*eta - INCsqrtZ + 2.) + BNsub) +
              pow(gsl_sf_psi_n(0, 2.*eta + 2.),2)*(2.*BNsub*eta - INCsqrtZ + 2.) +
        gsl_sf_psi_n(1, 2.*eta + 2.)*(-(2.*BNsub*eta) + INCsqrtZ - 2.)) -
         (M_gammaE*(eta*(2.*M_gammaE*eta + M_gammaE - 4.) - 1.) + 1.)*
           (2.*BNsub*eta - INCsqrtZ + 2.)))/exp(2*M_gammaE*eta);
}

double Deta_tau(double eta){
	// derivative of pow(1-tau,2*eta)/(2eta)*exp[-2*gam_E*eta]/Gamma[2eta]
	if(eta > -0.3){
		return (-2.*pow(1. - tau,2.*eta)*(M_gammaE - log(1. - tau) + gsl_sf_psi_n(0,1. + 2.*eta)))/
   (exp(2.*eta*M_gammaE)*tgamma(1. + 2.*eta));}
	else if (eta > -0.9){
		return -((2.*pow(1. - tau,2.*eta)*(-((2.*eta + 1.)*log(1. - tau)) +
            2.*M_gammaE*eta + (2.*eta + 1.)*gsl_sf_psi_n(0, 2.*eta + 2.) +
            M_gammaE - 1))/(exp(2.*M_gammaE*eta)*tgamma(2.*eta + 2.)));
	}
 else{
	 cout << "Eta is too large " << eta << " exiting program" << endl;
	 exit(0);
	 return 0;
 }
}

double DDeta_tau(double eta){
	// double derivative of pow(1-tau,2*eta)/(2eta)*exp[-2*gam_E*eta]/Gamma[2eta]
 if(eta > -0.3){
	 return (4.*pow(1. - tau,2*eta)
						*(pow(M_gammaE - log(1. - tau),2)
							+ 2.*(M_gammaE - log(1. - tau))*gsl_sf_psi_n(0,1. + 2.*eta)
							+ pow( gsl_sf_psi_n(0,1. + 2.*eta),2)
							- gsl_sf_psi_n(1,1. + 2.*eta)))
				 /(exp(2*eta*M_gammaE)*tgamma(1. + 2.*eta));
 }
	else if (eta > -0.9){
		 		return ((1./tgamma(2.*eta + 2.))*4.*pow(1. - tau,2.*eta)*
      ((M_gammaE - log(1. - tau))*(-((2.*eta + 1.)*log(1. - tau)) +
              2.*M_gammaE*eta + M_gammaE - 2.) +
         2.*gsl_sf_psi_n(0, 2*eta + 2)*(-((2.*eta + 1.)*log(1. - tau)) +
              2.*M_gammaE*eta + M_gammaE - 1.) +
         (2.*eta + 1.)*pow(gsl_sf_psi_n(0, 2.*eta + 2.),2.) -
         (2.*eta + 1.)*gsl_sf_psi_n(1, 2.*eta + 2.)))/exp(2.*M_gammaE*eta);
		 	}
	 	else{
	 		cout << "Eta is too large " << eta << " exiting program" << endl;
	 		exit(0);
	 		return 0;
	 	}
}

double Deta_tau_2nd(double eta, double BNsub, double INCsqrtZ){
	// derivative of (-BNsub*eta + INCsqrtZ*1/2 - 1)*(1 - tau)^(1 + 2 eta)* Exp[-2*EulerGamma*eta]*(2 eta)/(Gamma[2 + 2 eta])
	return ((1./tgamma(2.*eta + 2.))*pow(1. - tau,2.*eta + 1.)*
      (4.*M_gammaE*BNsub*pow(eta,2) + 2.*eta*log(1. - tau)*
           (-(2.*BNsub*eta) + INCsqrtZ - 2.) +  2.*eta*
      gsl_sf_psi_n(0, 2.*eta + 2.)*(2.*BNsub*eta - INCsqrtZ + 2.) -
         4.*BNsub*eta - 2.*M_gammaE*eta*INCsqrtZ +
     4.*M_gammaE*eta +
         INCsqrtZ - 2.))/exp(2.*M_gammaE*eta);
 }
 double DDeta_tau_2nd(double eta, double BNsub, double INCsqrtZ){
	 // double derivative of (-BNsub*eta + INCsqrtZ*1/2 - 1)*(1 - tau)^(1 + 2 eta)* Exp[-2*EulerGamma*eta]*(2 eta)/(Gamma[2 + 2 eta])
 	return ((-(1./tgamma(2.*eta + 2.)))*4.*pow(1. - tau,2.*eta + 1.)*
      (log(1. - tau)*(eta*log(1. - tau)*(2.*BNsub*eta - INCsqrtZ + 2.) -
        4.*BNsub*eta*(M_gammaE*eta - 1.) + (2.*M_gammaE*eta - 1.)*
                (INCsqrtZ - 2.)) + gsl_sf_psi_n(0, 2.*eta + 2.)*
           (2*eta*log(1. - tau)*(-(2.*BNsub*eta) + INCsqrtZ - 2.) +
              4.*eta*(M_gammaE*BNsub*eta - BNsub + M_gammaE) -
              2.*M_gammaE*eta*INCsqrtZ + INCsqrtZ - 2.) +
         eta*pow(gsl_sf_psi_n(0, 2.*eta + 2),2)*(2.*BNsub*eta - INCsqrtZ + 2.) +
         eta*gsl_sf_psi_n(1, 2.*eta + 2.)*(-(2*BNsub*eta) + INCsqrtZ - 2.) +
         BNsub*(2.*M_gammaE*eta*(M_gammaE*eta - 2.) + 1.) -
         M_gammaE*(M_gammaE*eta - 1.)*(INCsqrtZ - 2.)))/
   exp(2.*M_gammaE*eta);
  }
////////////////////
/// C coefficients
////////////////////

vector<complex<double>> C_higgs(){
	complex<double> Ufunc(0.);
	vector<complex<double>> Ctot;
	vector<complex<double>> Cs = hard_higgs(-Q*Q-I*1.E-16, muh*muh);
	complex<double> cT_H = cT_higgs_wo_heavy_top();
	//cout << "warning cT_h is 1" << endl;
	vector<complex<double>> Ssud = sudakov(muh, mus, GammaCg);
	vector<complex<double>> aC = acusp(muh, mus, GammaCg);
	vector<complex<double>> aphiS = acusp(muh, mus, gammaS);
  vector<complex<double>> aphiB = acusp(mus, muF, gammaB);
	complex<double> SSsud = Ssud[0]+Ssud[1]+Ssud[2], SaC = aC[0]+aC[1]+aC[2], SaphiS = aphiS[0]+aphiS[1]+aphiS[2], SaphiB = aphiB[0]+aphiB[1]+aphiB[2];
	//cout << "Ssud: " << Ssud[0]<<" " <<Ssud[1]<<" " <<Ssud[2] << ", tot=" << SSsud << endl;
	//cout << "aC: " << aC[0]<<" " <<aC[1]<<" " <<aC[2] << ", tot=" << SaC << endl;
	//cout << "aphiS: " << aphiS[0]<<" " <<aphiS[1]<<" " <<aphiS[2] << ", tot=" << SaphiS << endl;
	//cout << "aphiB: " << aphiB[0]<<" " <<aphiB[1]<<" " <<aphiB[2] << ", tot=" << SaphiB << endl;
	//complex<double> correction_factor = pow(M_PI,2)*beta0*CA*(falphasQ2(Q2)-falphasQ2(mus*mus))/(4.*M_PI*beta0);
	Ufunc = (abs(pow((-Q*Q-I*1.E-16)/(muh*muh),-2.*SaC)))*abs(exp(4.*SSsud+-2.*SaphiS+4.*SaphiB));
	Ctot = {ISLL*Cs[0]*cT_H*Ufunc, //NLL is the same as LL, but then is Ufunc different
					ISNNLL*(Cs[1])*cT_H*Ufunc};
	return Ctot;
}

vector<complex<double>> C_DY(){
	complex<double> Ufunc(0.);
	vector<complex<double>> Ctot;
	vector<complex<double>> Cv = hard_DY(-Q*Q-I*1.E-16, muh*muh);
	vector<complex<double>> Ssud = sudakov(muh, mus, GammaCq);
	vector<complex<double>> aC = acusp(muh, mus, GammaCq);
	vector<complex<double>> aV = acusp(muh, mus, gammaV);
	vector<complex<double>> aphi = acusp(mus, muF, gammaphi);
	//cout << "Cv = " << Cv << ", Ssud =" << Ssud << ", aC =" << aC << ", aV =" << aV << ", aphi =" << aphi << endl;
	complex<double> SSsud = Ssud[0]+Ssud[1]+Ssud[2], SaC = aC[0]+aC[1]+aC[2], SaV = aV[0]+aV[1]+aV[2], Saphi = aphi[0]+aphi[1]+aphi[2];
	//cout << "hard " << Cv[0] << " " << Cv[1] << endl;
	//cout << "Ssud: " << Ssud[0]<<" " <<Ssud[1]<<" " <<Ssud[2] << ", tot=" << SSsud << endl;
	//cout << "aC: " << aC[0]<<" " <<aC[1]<<" " <<aC[2] << ", tot=" << SaC << endl;
	//cout << "aV: " << aV[0]<<" " <<aV[1]<<" " <<aV[2] << ", tot=" << SaV << endl;
	//cout << "aphi: " << aphi[0]<<" " <<aphi[1]<<" " <<aphi[2] << ", tot=" << Saphi << endl;

	Ufunc = abs(pow((-Q*Q-I*1.E-16)/(muh*muh),-2.*SaC))*abs(exp(4.*SSsud-2.*SaV+4.*Saphi));

	Ctot = {ISLL*Cv[0]*Ufunc, //NLL is the same as LL, but then is Ufunc different
					ISNNLL*Cv[1]*Ufunc};
	return Ctot;
}

// see notes leonardo for formula
complex<double> C_higgs_NLP(){
	complex<double> Ufunc(0.);
	complex<double> Ctot;
	vector<complex<double>> Cs = hard_higgs(-Q*Q-I*1.E-16, muh*muh);
	vector<complex<double>> SsudH = sudakov(muh, Q, GammaCg);
	vector<complex<double>> SsudS = sudakov(mus, Q, GammaCg);
	complex<double> alpha_Q2 = falphasQ2(Q*Q);
	complex<double> alpha_mus2 = falphasQ2(mus*mus);
	complex<double> SSsud = SsudH[0]-SsudS[0];
	Ufunc = abs(exp(4.*SSsud));
	Ctot = (Cs[0]+ISNNLL*INCHARD*Cs[1]+ISNNLL*INCHARD*CA*alpha_mus2/(2.*M_PI)*zeta2)*Ufunc*(-8.*CA/beta0)*real(log(alpha_Q2/alpha_mus2));
	return Ctot;
}
complex<double> C_DY_NLP(){
	complex<double> Ufunc(0.);
	complex<double> Ctot;

	vector<complex<double>> Cv = hard_DY(-Q*Q-I*1.E-16, muh*muh);
	vector<complex<double>> SsudH = sudakov(muh, Q, GammaCq);
	vector<complex<double>> SsudS = sudakov(mus, Q, GammaCq);
	complex<double> alpha_Q2 = falphasQ2(Q*Q);
	complex<double> alpha_mus2 = falphasQ2(mus*mus);
	complex<double> SSsud = SsudH[0]-SsudS[0];
	Ufunc = abs(exp(4.*SSsud));
	Ctot = (Cv[0]+ISNNLL*INCHARD*Cv[1]+ISNNLL*INCHARD*CF*alpha_mus2/(2.*M_PI)*zeta2)*Ufunc*(-8.*CF/beta0)*real(log(alpha_Q2/alpha_mus2));
	return Ctot;
}

double diff_xsec_DY(double z, double t){
	  // overal z dependence is pow(z,-eta-1.)*pow(1-z,-1+2*eta)
		// or pow(z,-1.)*pow(1-z,-1+2*eta) without eta improvement

// jacobian for integration from z = tau to z = 1, x = tau/z (or tau) to x = 1 (x = exp(t*log(tau/z)))
		double jacobian = (1.-tau);

// setting up	the resummation functions
		vector<complex<double>> res_eta = feta(mus,muF,GammaCq);
		double eta = real(res_eta[0]+res_eta[1]);
		double eulerGam = 2.*eta*exp(-2.*M_gammaE*INCEULER*eta)/tgamma(1.+2.*eta); //to avoid issues when eta = 0.;
	  if(eulerGam == 0) eulerGam = 1.;
		vector<complex<double>> CDY = C_DY();

		vector<double> DeltaDY = {DY_LO_factor()*real(CDY[0]),
				                      DY_LO_factor()*real(CDY[1])}	; //tau*simgaLO*resummedC (no z dependence)
		vector<double> fsoft = fsoft_DY(log(Q*Q/(mus*mus))); // log should be added, but to compare turn it off
		double BNfactor = pow(z,-eta*BN);
		if(INCSQRTZ){BNfactor = BNfactor*pow(z,1./2.);}
// luminosity
		//cout << DeltaDY[0] << " " << DeltaDY[1] << endl;
		//cout << DeltaDY[0] << " " << DeltaDY[1] << endl;
		vector<double> lumi_at_z = deriv_to_y_luminosities("qqbar",t,tau/z);
		vector<double> lumi_at_tau = deriv_to_y_luminosities("qqbar",t,tau);
		//cout << " eulergam " << eulerGam << endl;
// z integration terms at LP without eta derivative
		double wo_eta_deriv = jacobian*(DeltaDY[0]*(fsoft[0])+DeltaDY[1])
																	*eulerGam
																	*(BNfactor*lumi_at_z[0]/z-lumi_at_tau[0])
													*1./pow(1.-z,1.-2.*eta); //from z = tau to z = 1
// addition term, note that here only x integration is needed (from tau to 1)
		double add_correction = (DeltaDY[0]*(fsoft[0])+DeltaDY[1])
																*eulerGam
																*lumi_at_tau[0]
															*pow(1.-tau,2.*eta)/(2.*eta); //no integration for z
		if(eta < -0.3){
			wo_eta_deriv = jacobian*(DeltaDY[0]*(fsoft[0])+DeltaDY[1])
													*eulerGam
													*((BNfactor*lumi_at_z[0]/z-lumi_at_tau[0])/pow(1.-z,1.-2.*eta)
														+(-(1.+BN*eta)*lumi_at_tau[0]-tau*lumi_at_tau[1])/pow(1.-z,-2.*eta)); //from z = tau to z = 1
			add_correction = (DeltaDY[0]*(fsoft[0])+DeltaDY[1])
																*eulerGam
																*(pow(1.-tau,2.*eta)/(2.*eta)*lumi_at_tau[0]
																	- pow(1.-tau,1.+2.*eta)/(1.+2.*eta)*(-(1.+BN*eta)*lumi_at_tau[0]-tau*lumi_at_tau[1])); //no integration for z
		}
// total result without eta derivatives at LP
		double tot_no_eta_deriv = wo_eta_deriv+add_correction;

// z integration terms with one derivative of eta
		double tot_eta_deriv = 0;
		if(ISNNLL != 0){


				double one_eta_deriv = jacobian*DeltaDY[0]*(fsoft[1])*(lumi_at_z[0]/z*Deta_z(z,eta)-lumi_at_tau[0]*Deta_z_eta(z,eta,0.));
				double two_eta_deriv = jacobian*DeltaDY[0]*(fsoft[2])*(lumi_at_z[0]/z*DDeta_z(z,eta)-lumi_at_tau[0]*DDeta_z_eta(z,eta,0.));
				double one_eta_add_correction = DeltaDY[0]*(fsoft[1])*lumi_at_tau[0]*Deta_tau(eta); //no integration for z
				double two_eta_add_correction = DeltaDY[0]*(fsoft[2])*lumi_at_tau[0]*DDeta_tau(eta); //no integration for z

				if(eta < -0.3){
					one_eta_deriv = jacobian*DeltaDY[0]*fsoft[1]
															*(lumi_at_z[0]/z*Deta_z(z,eta)-lumi_at_tau[0]*Deta_z_eta(z,eta,0.)
																+lumi_at_tau[0]*Deta_z_eta2(z, eta, BN, INCSQRTZ)+tau*lumi_at_tau[1]*Deta_z_eta2(z, eta, 0,0)); //from z = tau to z = 1
					two_eta_deriv = jacobian*DeltaDY[0]*fsoft[2]
															*(lumi_at_z[0]/z*DDeta_z(z,eta)-lumi_at_tau[0]*DDeta_z_eta(z,eta,0.)
																+lumi_at_tau[0]*DDeta_z_eta2(z, eta, BN, INCSQRTZ)+tau*lumi_at_tau[1]*DDeta_z_eta2(z, eta, 0,0)); // -tau comes from dy/dz with y = tau/z
					one_eta_add_correction = DeltaDY[0]*fsoft[1]//finite
																	*(lumi_at_tau[0]*Deta_tau(eta)
																		-lumi_at_tau[0]*Deta_tau_2nd(eta,BN,INCSQRTZ)-tau*lumi_at_tau[1]*Deta_tau_2nd(eta,0,0)); //no integration for z
					two_eta_add_correction = DeltaDY[0]*fsoft[2]//finite
																	*(lumi_at_tau[0]*DDeta_tau(eta)
																		-lumi_at_tau[0]*DDeta_tau_2nd(eta,BN, INCSQRTZ)-tau*lumi_at_tau[1]*DDeta_tau_2nd(eta,0,0)); //no integration for z

				}


				tot_eta_deriv = one_eta_deriv+two_eta_deriv
													+ one_eta_add_correction+two_eta_add_correction;
			}
// NLP result

		double NLP_DeltaDY = ISNLP*DY_LO_factor()*real(C_DY_NLP()); //tau*simgaLO*resummedC (no z dependence)
		double NLP_wo_eta_deriv = jacobian*NLP_DeltaDY*lumi_at_z[0]/z;

		return tot_no_eta_deriv + tot_eta_deriv + NLP_wo_eta_deriv;
}


double xsec_higgs(double z, double t){
	  // overal z dependence is pow(z,-eta-1.)*pow(1-z,-1+2*eta)
		// or pow(z,-1.)*pow(1-z,-1+2*eta) without eta improvement

// jacobian for integration from z = tau to z = 1, x = tau/z (or tau) to x = 1 (x = exp(t*log(tau/z)))
		double jacobian = (1.-tau);

// setting up	the resummation functions
		vector<complex<double>> res_eta = feta(mus,muF,GammaCg);
		double eta = real(res_eta[0]+res_eta[1]);
		double eulerGam = 2.*eta*exp(-2.*M_gammaE*INCEULER*eta)/tgamma(1.+2.*eta); //to avoid issues when eta = 0.;
		if(eulerGam == 0) {eulerGam = 1.;}
		vector<complex<double>> CHiggs = C_higgs();
		vector<double> DeltaHiggs = {higgs_LO_factor()*real(CHiggs[0]), //O(1)
			 													 higgs_LO_factor()*real(CHiggs[1])}; //O(as) (at NNLL)
		vector<double> fsoft = fsoft_higgs(log(Q*Q/(mus*mus)));

		double BNfactor = pow(z,-eta*BN);
		if(INCSQRTZ){BNfactor = BNfactor*pow(z,1./2.);}
// luminosity
		vector<double> lumi_at_z = deriv_to_y_luminosities("gg", t,tau/z);
		vector<double> lumi_at_tau = deriv_to_y_luminosities("gg", t,tau);

// z integration terms at LP without eta derivative

		double wo_eta_deriv = jacobian*(DeltaHiggs[0]*(fsoft[0])+DeltaHiggs[1])
																	*eulerGam
																	*(BNfactor*lumi_at_z[0]/z-lumi_at_tau[0])
													*1./pow(1.-z,1.-2.*eta); //from z = tau to z = 1
		double add_correction = (DeltaHiggs[0]*(fsoft[0])+DeltaHiggs[1])
																	*eulerGam
																	*lumi_at_tau[0]
													*pow(1.-tau,2.*eta)/(2.*eta); //no integration for z
		if(eta < -0.3){
				wo_eta_deriv = jacobian*(DeltaHiggs[0]*(fsoft[0])+DeltaHiggs[1])
													*eulerGam
													*((BNfactor*lumi_at_z[0]/z-lumi_at_tau[0])/pow(1.-z,1.-2.*eta)
														+(-(1.+BN*eta)*lumi_at_tau[0]-tau*lumi_at_tau[1])/pow(1.-z,-2.*eta)); //from z = tau to z = 1
				add_correction = (DeltaHiggs[0]*(fsoft[0])+DeltaHiggs[1])
													*eulerGam
													*(pow(1.-tau,2.*eta)/(2.*eta)*lumi_at_tau[0]
														- pow(1.-tau,1.+2.*eta)/(1.+2.*eta)*(-(1.+BN*eta)*lumi_at_tau[0]-tau*lumi_at_tau[1])); //no integration for z
		}
		// total result without eta derivatives at LP
		double tot_no_eta_deriv = wo_eta_deriv+add_correction;

		// terms with eta derivatives
		double tot_eta_deriv = 0;
		if(ISNNLL != 0){
				double one_eta_deriv = jacobian*DeltaHiggs[0]*(fsoft[1])*(lumi_at_z[0]/z*Deta_z(z,eta)-lumi_at_tau[0]*Deta_z_eta(z,eta,0.));
				double two_eta_deriv = jacobian*DeltaHiggs[0]*(fsoft[2])*(lumi_at_z[0]/z*DDeta_z(z,eta)-lumi_at_tau[0]*DDeta_z_eta(z,eta,0.));
				double one_eta_add_correction = DeltaHiggs[0]*(fsoft[1])*lumi_at_tau[0]*Deta_tau(eta); //no integration for z
				double two_eta_add_correction = DeltaHiggs[0]*(fsoft[2])*lumi_at_tau[0]*DDeta_tau(eta); //no integration for z
				if(eta < -0.3){
					one_eta_deriv = jacobian*DeltaHiggs[0]*fsoft[1]
															*(lumi_at_z[0]/z*Deta_z(z,eta)-lumi_at_tau[0]*Deta_z_eta(z,eta,0.)
																+lumi_at_tau[0]*Deta_z_eta2(z, eta, BN, INCSQRTZ)+tau*lumi_at_tau[1]*Deta_z_eta2(z, eta, 0,0)); //from z = tau to z = 1
					two_eta_deriv = jacobian*DeltaHiggs[0]*fsoft[2]
															*(lumi_at_z[0]/z*DDeta_z(z,eta)-lumi_at_tau[0]*DDeta_z_eta(z,eta,0.)
																+lumi_at_tau[0]*DDeta_z_eta2(z, eta, BN, INCSQRTZ)+tau*lumi_at_tau[1]*DDeta_z_eta2(z, eta, 0,0)); // -tau comes from dy/dz with y = tau/z
					one_eta_add_correction = DeltaHiggs[0]*fsoft[1]//finite
																	*(lumi_at_tau[0]*Deta_tau(eta)
																		-lumi_at_tau[0]*Deta_tau_2nd(eta,BN,INCSQRTZ)-tau*lumi_at_tau[1]*Deta_tau_2nd(eta,0,0)); //no integration for z
					two_eta_add_correction = DeltaHiggs[0]*fsoft[2]//finite
																	*(lumi_at_tau[0]*DDeta_tau(eta)
																		-lumi_at_tau[0]*DDeta_tau_2nd(eta,BN, INCSQRTZ)-tau*lumi_at_tau[1]*DDeta_tau_2nd(eta,0,0)); //no integration for z

				}


				tot_eta_deriv = one_eta_deriv+two_eta_deriv
													+ one_eta_add_correction+two_eta_add_correction;
			}

// NLP result
		double NLP_DeltaHiggs = ISNLP*higgs_LO_factor()*real(C_higgs_NLP()); //tau*simgaLO*resummedC (no z dependence)
		double NLP_wo_eta_deriv = jacobian*NLP_DeltaHiggs*lumi_at_z[0]/z;

		return tot_no_eta_deriv + tot_eta_deriv + NLP_wo_eta_deriv;
}
