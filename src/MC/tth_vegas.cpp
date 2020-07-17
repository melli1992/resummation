#include "k_factors_ttH.h"
#include "tth_softanom.h"
#include "resum_tth.h"
#include "tth_vegas.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "parameters.h"
#include <limits>
#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using namespace std;

double vegas_ttH_Nspace_deform(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	rho	s34, 	angles,
	//k[0]	k[1]	k[2]	k[3]-k[5],
	double s34 = 4.*pow(mt,2)+(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*k[2];
	Q2 = pow(mH+sqrt(s34),2);
	tau = Q2/S2;
	double xvar = k[1]-1.;
	double wr = xvar/(1.+xvar);
	double wjac = 1./(pow(1.+xvar,2));
	complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
    complex<double> rhoin = exp(wr/N);
	if(k[0] == 1){return 0;}
	else if(tau > 1){return 0;}
   	else if(k[1] == 0){return 0;}
   	else{
		vector<complex<double>> partonic = xsec_LO_c(rhoin, s34, k[3], k[4], k[5]) ;
		double result = 2*(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_mellin_pdf_sum_gg(N)*partonic[1] + fit_mellin_pdf_sum_qqbar(N)*partonic[0])
					*wjac*exp(wr)/N*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}

double vegas_ttH_Nspace(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	rho	s34, 	angles,
	//k[0]	k[1]	k[2]	k[3]-k[5],

	double s34 = 4.*pow(mt,2)+(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*k[2];
	Q2 = pow(mH+sqrt(s34),2);
	tau = Q2/S2;
	if(k[0] == 1){return 0;}
	else if(tau > 1){return 0;}
   	else if(k[1] == 0){return 0;}
   	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<double> partonic = xsec_LO(k[1], s34, k[3], k[4], k[5]) ;
		double result = 2*((pow(sqrt(S2)-mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_mellin_pdf_sum_gg(N)*partonic[1] + fit_mellin_pdf_sum_qqbar(N)*partonic[0])
					*pow(k[1], N - 1.)*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}

double vegas_ttH_LO(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double s34 = 4.*pow(mt,2)+(pow(sqrt(k[0]*k[1]*S2)-mH,2)-4*pow(mt,2))*k[2];
	Q2 = pow(mH+sqrt(s34),2);
	if(k[0]*k[1]*S2  < Q2) {return 0;}
   	else{
		vector<double> partonic = xsec_LO(Q2/(k[0]*k[1]*S2), s34, k[3], k[4], k[5]) ;
		double result = 0;
		if(!fitPDF) result = ((pow(sqrt(k[0]*k[1]*S2)-mH,2)-4*pow(mt,2)))*(gluon_d0PDF(k[0],  k[1])*partonic[1] + quark_d0PDF(k[0],  k[1])*partonic[0]);
		if(fitPDF) result = ((pow(sqrt(k[0]*k[1]*S2)-mH,2)-4*pow(mt,2)))*(fit_gluon_d0PDF(k[0],  k[1])*partonic[1] + fit_quark_d0PDF(k[0],  k[1])*partonic[0]);
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}


double vegas_ttH_LO1st(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 x1, 	x2, 	rho	s34, 	angles, 	gluon/quark : x1*x2*S = tau/x*S
	//k[0]	 k[1]	k[2]	k[3]	k[4]	k[5]-k[7], 	0 or 1
	double s34 = 4.*pow(mt,2)+(pow(sqrt(Q2/k[3]) - mH,2)-4*pow(mt,2))*k[4];
	if(k[0] == 1){return 0;} //otherwise we get 1/0
	else if(k[3] == 0){return 0;}
	else if((Q2/k[3] > k[1]*k[2]*S2) && boundary == 0){return 0;}
	else if((Q2/k[3] <= k[1]*k[2]*S2) && boundary == 1){return 0;} // for the boundary term
	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<double> partonic = xsec_LO(k[3], s34, k[5], k[6], k[7]) ;
		double result = 0;
		if(!fitPDF) result = 2*((pow(sqrt(Q2/k[3]) - mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(gluon_d1PDF(k[1],  k[2])*partonic[1] + quark_d1PDF(k[1],  k[2])*partonic[0])
					*pow(k[3], N - 1.)*pow(N,-2.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if(fitPDF) result = 2*((pow(sqrt(Q2/k[3]) - mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_gluon_d1PDF(k[1],  k[2])*partonic[1] + fit_quark_d1PDF(k[1],  k[2])*partonic[0])
					*pow(k[3], N - 1.)*pow(N,-2.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}


		}
}



double vegas_ttH_LO2nd(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 x1, 	x2, 	rho	s34, 	angles, 	gluon/quark : x1*x2*S = tau/x*S
	//k[0]	 k[1]	k[2]	k[3]	k[4]	k[5]-k[7], 	0 or 1
	double s34 = 4.*pow(mt,2)+(pow(sqrt(Q2/k[3]) - mH,2)-4*pow(mt,2))*k[4];
	if(k[0] == 1){return 0;}
	else if(k[3] == 0){return 0;}
	else if((Q2/k[3] > k[1]*k[2]*S2) && boundary == 0){return 0;}
	else if((Q2/k[3] <= k[1]*k[2]*S2) && boundary == 1){return 0;} // for the boundary term
	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<double> partonic = xsec_LO(k[3], s34, k[5], k[6], k[7]) ;
		double result = 0;
		if(!fitPDF) result = 2*((pow(sqrt(Q2/k[3]) - mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(gluon_d2PDF(k[1],  k[2])*partonic[1] + quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(k[3], N - 1.)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if(fitPDF) result = 2*((pow(sqrt(Q2/k[3]) - mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_gluon_d2PDF(k[1],  k[2])*partonic[1] + fit_quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(k[3], N - 1.)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}


		}
}


double vegas_ttH_resum_2nd(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 x1, 	x2, 	rho	s34, 	angles, 	gluon/quark : x1*x2*S = tau/x*S
	//k[0]	 k[1]	k[2]	k[3]	k[4]	k[5]-k[7], 	0 or 1
	//double s34 = 4.*pow(mt,2)+(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*k[4];
	double s34 = 4.*pow(mt,2)+(pow(sqrt(Q2/k[3]) - mH,2)-4*pow(mt,2))*k[4];
	if(k[0] == 1){return 0;}
	else if(k[3] == 0){return 0;}
	else if((Q2/k[3] > k[1]*k[2]*S2) && boundary == 0){return 0;}
	else if((Q2/k[3] <= k[1]*k[2]*S2) && boundary == 1){return 0;} // for the boundary term
	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<complex<double>> partonic = xsec_res(N,k[3], s34, k[5], k[6], k[7]) ;
		double result = 0;
		if(!fitPDF) result = 2*((pow(sqrt(Q2/k[3]) - mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(gluon_d2PDF(k[1],  k[2])*partonic[1] + quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(k[3], N - 1.)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if(fitPDF) result = 2*((pow(sqrt(Q2/k[3]) - mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_gluon_d2PDF(k[1],  k[2])*partonic[1] + fit_quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(k[3], N - 1.)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}


		}
}


double vegas_ttH_resum_2nd_z5(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 x1, 	x2, 	rho	s34, 	angles, 	gluon/quark : x1*x2*S = tau/x*S
	//k[0]	 k[1]	k[2]	k[3]	k[4]	k[5]-k[7], 	0 or 1
	//double s34 = 4.*pow(mt,2)+(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*k[4];
	double s34 = 4.*pow(mt,2)+(pow(sqrt(Q2/k[3]) - mH,2)-4*pow(mt,2))*k[4];
	if(k[0] == 1){return 0;}
	else if(k[3] == 0){return 0;}
	else if((Q2/k[3] > k[1]*k[2]*S2) && boundary == 0){return 0;}
	else if((Q2/k[3] <= k[1]*k[2]*S2) && boundary == 1){return 0;} // for the boundary term
	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<complex<double>> partonic = xsec_res_c_z5(N,k[3], s34, k[5], k[6], k[7]) ;
		double result = 0;
		if(!fitPDF) result = 2*((pow(sqrt(Q2/k[3]) - mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(gluon_d2PDF(k[1],  k[2])*partonic[1] + quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(k[3], N - 1.)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if(fitPDF) result = 2*((pow(sqrt(Q2/k[3]) - mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_gluon_d2PDF(k[1],  k[2])*partonic[1] + fit_quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(k[3], N - 1.)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}


		}
}



double vegas_ttH_resum_2nd_def(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 x1, 	x2, 	rho	s34, 	angles, 	gluon/quark : x1*x2*S = tau/x*S
	//k[0]	 k[1]	k[2]	k[3]	k[4]	k[5]-k[7], 	0 or 1
	double s34 = 4.*pow(mt,2)+(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*k[4];
	if(k[0] == 1){return 0;} //otherwise we get 1/0
	Q2 = pow(mH+sqrt(s34),2);
	tau = Q2/S2;
	if(k[0] == 1){return 0;} //otherwise we get 1/0
	else if(k[1]*k[2]*S2  < Q2) {return 0;}
	else if(k[3] == 0){return 0;} //you get s = W2/rho so that is then infinity
	else if(Q2/k[3] > k[1]*k[2]*S2){return 0;} // s = W2/rho, cannot be bigger than S
	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<complex<double>> partonic = xsec_res(N,k[3], s34, k[5], k[6], k[7]) ;
		double result = 0;
		if(!fitPDF) result = 2*((pow(sqrt(S2)-mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(gluon_d2PDF(k[1],  k[2])*partonic[1] + quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(k[3], N - 1.)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if(fitPDF) result = 2*((pow(sqrt(S2)-mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_gluon_d2PDF(k[1],  k[2])*partonic[1] + fit_quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(k[3], N - 1.)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}


		}
}


double vegas_ttH_resum_2nd_def_z5(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 x1, 	x2, 	rho	s34, 	angles, 	gluon/quark : x1*x2*S = tau/x*S
	//k[0]	 k[1]	k[2]	k[3]	k[4]	k[5]-k[7], 	0 or 1
	double s34 = 4.*pow(mt,2)+(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*k[4];
	if(k[0] == 1){return 0;} //otherwise we get 1/0
	Q2 = pow(mH+sqrt(s34),2);
	tau = Q2/S2;
	if(k[0] == 1){return 0;} //otherwise we get 1/0
	else if(k[1]*k[2]*S2  < Q2) {return 0;}
	else if(k[3] == 0){return 0;} //you get s = W2/rho so that is then infinity
	else if(Q2/k[3] > k[1]*k[2]*S2){return 0;} // s = W2/rho, cannot be bigger than S
	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<complex<double>> partonic = xsec_res_c_z5(N,k[3], s34, k[5], k[6], k[7]) ;
		double result = 0;
		if(!fitPDF) result = 2*((pow(sqrt(S2)-mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(gluon_d2PDF(k[1],  k[2])*partonic[1] + quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(k[3], N - 1.)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if(fitPDF) result = 2*((pow(sqrt(S2)-mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_gluon_d2PDF(k[1],  k[2])*partonic[1] + fit_quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(k[3], N - 1.)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}


		}
}


double vegas_ttH_resum_def(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	rho	s34, 	angles,
	//k[0]	k[1]	k[2]	k[3]-k[5],
	double s34 = 4.*pow(mt,2)+(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*k[2];
	Q2 = pow(mH+sqrt(s34),2);
	tau = Q2/S2;
	double xvar = k[1]-1.;
	double wr = xvar/(1.+xvar);
	double wjac = 1./(pow(1.+xvar,2));
	complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
    complex<double> rhoin = exp(wr/N);
	if(k[0] == 1){return 0;}
	else if(tau > 1){return 0;}
   	else if(k[1] == 0){return 0;}
   	else{
		vector<complex<double>> partonic = xsec_res_c(N,rhoin, s34, k[3], k[4], k[5]) ;
		double result = 2*(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_mellin_pdf_sum_gg(N)*partonic[1] + fit_mellin_pdf_sum_qqbar(N)*partonic[0])
					*wjac*exp(wr)/N*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}


double vegas_ttH_resum_def_z5(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	rho	s34, 	angles,
	//k[0]	k[1]	k[2]	k[3]-k[5],
	double s34 = 4.*pow(mt,2)+(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*k[2];
	Q2 = pow(mH+sqrt(s34),2);
	tau = Q2/S2;
	double xvar = k[1]-1.;
	double wr = xvar/(1.+xvar);
	double wjac = 1./(pow(1.+xvar,2));
	complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
    complex<double> rhoin = exp(wr/N);
	if(k[0] == 1){return 0;}
	else if(tau > 1){return 0;}
   	else if(k[1] == 0){return 0;}
   	else{
		vector<complex<double>> partonic = xsec_res_c_z5(N,rhoin, s34, k[3], k[4], k[5]) ;
		double result = 2*(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_mellin_pdf_sum_gg(N)*partonic[1] + fit_mellin_pdf_sum_qqbar(N)*partonic[0])
					*wjac*exp(wr)/N*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}


double vegas_ttH_inv_mass_LO(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 x1, 	x2, 	rho	s34, 	angles, 	gluon/quark : x1*x2*S = tau/x*S
	//k[0]	 k[1]	k[2]	k[3]	k[4]	k[5]-k[7], 	0 or 1
	//double s34 = 4.*pow(mt,2)+(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*k[4];
	//Q2 = pow(2.*mt+mH,2)+(S2-pow(2.*mt+mH,2))*k[7];
	double s34 = 4.*pow(mt,2)+(pow(sqrt(Q2) - mH,2)-4*pow(mt,2))*k[3];
	if(k[0] == 1){return 0;}
	else if((Q2 > k[1]*k[2]*S2) && boundary == 0){return 0;}
	else if((pow(2.*mt+mH,2) > k[1]*k[2]*S2)){return 0;}
	else if((Q2 <= k[1]*k[2]*S2) && boundary == 1){return 0;} // for the boundary term
	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<double> partonic = xsec_LO(1., s34, k[4], k[5], k[6]) ;
		double result = 0;
		if(!fitPDF) result = 2./sqrt(Q2)*2.*(pow(sqrt(Q2) - mH,2)-4*pow(mt,2))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(gluon_d2PDF(k[1],  k[2])*partonic[1] + quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(N,-4.)*pow(k[1]*k[2], N)*pow(Q2/S2, -N));
					//(pow(sqrt(Q2) - mH,2)-4*pow(mt,2)) is jacobian factor from stt integration. If you want to integrate over Q too than you need this jacobian factor to be included as well!!
		if(fitPDF) result = 2./sqrt(Q2)*2.*(pow(sqrt(Q2) - mH,2)-4*pow(mt,2))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_gluon_d2PDF(k[1],  k[2])*partonic[1] + fit_quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(N,-4.)*pow(k[1]*k[2], N)*pow(Q2/S2, -N));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}

double vegas_ttH_inv_mass_res(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 x1, 	x2, 	rho	s34, 	angles, 	gluon/quark : x1*x2*S = tau/x*S
	//k[0]	 k[1]	k[2]	k[3]	k[4]	k[5]-k[7], 	0 or 1
	//double s34 = 4.*pow(mt,2)+(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*k[4];
	double s34 = 4.*pow(mt,2)+(pow(sqrt(Q2) - mH,2)-4*pow(mt,2))*k[3];
	if(k[0] == 1){return 0;}
	else if((Q2 > k[1]*k[2]*S2) && boundary == 0){return 0;}
	else if((pow(2.*mt+mH,2) > k[1]*k[2]*S2)){return 0;}
	else if((Q2 <= k[1]*k[2]*S2) && boundary == 1){return 0;} // for the boundary term
	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<complex<double>> partonic = xsec_res_c_z5(N,1., s34, k[4], k[5], k[6]) ;
		double result = 0;
		if(!fitPDF) result = 2./sqrt(Q2)*2.*(pow(sqrt(Q2) - mH,2)-4*pow(mt,2))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(gluon_d2PDF(k[1],  k[2])*partonic[1] + quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(N,-4.)*pow(k[1]*k[2], N)*pow(Q2/S2, -N));
		if(fitPDF) result = 2./sqrt(Q2)*(pow(sqrt(Q2) - mH,2)-4*pow(mt,2))*2*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_gluon_d2PDF(k[1],  k[2])*partonic[1] + fit_quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(N,-4.)*pow(k[1]*k[2], N)*pow(Q2/S2, -N));
		if (isnan(result)){return 0;}
		else{ return result;}


		}
}

double vegas_ttH_inv_mass_res_s34(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 x1, 	x2, 	rho	s34, 	angles, 	gluon/quark : x1*x2*S = tau/x*S
	//k[0]	 k[1]	k[2]	k[3]	k[4]	k[5]-k[7], 	0 or 1
	//double s34 = 4.*pow(mt,2)+(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*k[4];
	double s34 = 4.*pow(mt,2)+(pow(sqrt(s2)-mH,2)-4*pow(mt,2))*k[3];
	if(k[0] == 1){return 0;} //otherwise we get 1/0
	Q2 = pow(mH+sqrt(s34),2);
	tau = Q2/S2;
	double rho_set = Q2/s2;
	if(k[0] == 1){return 0;} //otherwise we get 1/0
	else if(k[1]*k[2]*S2  < Q2) {return 0;}
	else if(rho_set == 0){return 0;} //you get s = W2/rho so that is then infinity
	else if(Q2/rho_set > k[1]*k[2]*S2){return 0;} // s = W2/rho, cannot be bigger than S
	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<complex<double>> partonic = xsec_res_c_z5(N,rho_set, s34, k[4], k[5], k[6]) ;
		double result = 0;
		if(!fitPDF) result = 4*1./sqrt(s2)*((pow(sqrt(s2)-mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(gluon_d2PDF(k[1],  k[2])*partonic[1] + quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(rho_set, N)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if(fitPDF) result = 4*1./sqrt(s2)*((pow(sqrt(s2)-mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_gluon_d2PDF(k[1],  k[2])*partonic[1] + fit_quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(rho_set, N)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}


double vegas_ttH_inv_mass_LO_s34(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	//N,	 x1, 	x2, 	rho	s34, 	angles, 	gluon/quark : x1*x2*S = tau/x*S
	//k[0]	 k[1]	k[2]	k[3]	k[4]	k[5]-k[7], 	0 or 1
	//double s34 = 4.*pow(mt,2)+(pow(sqrt(S2)-mH,2)-4*pow(mt,2))*k[4];
	double s34 = 4.*pow(mt,2)+(pow(sqrt(s2)-mH,2)-4*pow(mt,2))*k[3];
	if(k[0] == 1){return 0;} //otherwise we get 1/0
	Q2 = pow(mH+sqrt(s34),2);
	tau = Q2/S2;
	double rho_set = Q2/s2;
	if(k[0] == 1){return 0;} //otherwise we get 1/0
	else if(k[1]*k[2]*S2  < Q2) {return 0;}
	else if(rho_set == 0){return 0;} //you get s = W2/rho so that is then infinity
	else if(Q2/rho_set > k[1]*k[2]*S2){return 0;} // s = W2/rho, cannot be bigger than S
	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<double> partonic = xsec_LO(rho_set, s34, k[4], k[5], k[6]) ;
		double result = 0;
		if(!fitPDF) result = 4*1./sqrt(s2)*((pow(sqrt(s2)-mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(gluon_d2PDF(k[1],  k[2])*partonic[1] + quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(rho_set, N)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if(fitPDF) result = 4*1./sqrt(s2)*((pow(sqrt(s2)-mH,2)-4*pow(mt,2)))*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_gluon_d2PDF(k[1],  k[2])*partonic[1] + fit_quark_d2PDF(k[1],  k[2])*partonic[0])
					*pow(rho_set, N)*pow(N,-4.)*pow(k[1]*k[2], N)*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}

double vegas_ttH_LO_pT(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double mT = sqrt(pT2+mH2);
	double mTtt = sqrt(pT2+4.*mt2);
	double s = (k[0]*k[1]*S2);
	double jacs34 = (s+mH2-2.*sqrt(s)*mT-4*mt2);
	double s34 = 4.*mt2+jacs34*k[2];
	vector<double> partonic = xsec_LO_pT((pow(mT+mTtt,2))/(s), pT2,s34, k[3], k[4]) ;
	double result = 0;
	if(!fitPDF) result = jacs34*(gluon_d0PDF(k[0],  k[1])*partonic[1] + quark_d0PDF(k[0],  k[1])*partonic[0]);
	if(fitPDF) result = jacs34*(fit_gluon_d0PDF(k[0],  k[1])*partonic[1] + fit_quark_d0PDF(k[0],  k[1])*partonic[0]);
	if (isnan(result)){return 0;}
	else{ return result;}
}

double vegas_ttH_LO_pT_stt_dist(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double mT = sqrt(pT2+mH2);
	double mTtt = sqrt(pT2+4.*mt2);
	double s = (k[0]*k[1]*S2);
	double jacs34 = 1.;//(s+mH2-2.*sqrt(s)*mT-4*mt2);
	double s34 = s34in;//4.*mt2+jacs34*k[2];
	if(s34 > s+mH2-2.*sqrt(s)*mT){return 0;}
	vector<double> partonic = xsec_LO_pT((pow(mT+mTtt,2))/(s), pT2,s34, k[2], k[3]) ;
	double result = 0;
	if(!fitPDF) result = jacs34*(gluon_d0PDF(k[0],  k[1])*partonic[1] + quark_d0PDF(k[0],  k[1])*partonic[0]);
	if(fitPDF) result = jacs34*(fit_gluon_d0PDF(k[0],  k[1])*partonic[1] + fit_quark_d0PDF(k[0],  k[1])*partonic[0]);
	if (isnan(result)){return 0;}
	else{ return result;}
}

double vegas_ttH_LO_pTfull(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double s = (k[0]*k[1]*S2);
	double pT2 = S2*k[2];
	double pT = sqrt(pT2);
	double mT = sqrt(pT2+mH2);
	double mTtt = sqrt(pT2+4.*mt2);
	double jacs34 = (s+mH2-2.*sqrt(s)*mT-4*mt2);
	double s34 = 4.*mt2+jacs34*k[3];
	//Q2 = pT2 here!

	vector<double> partonic = xsec_LO_pT((pow(mT+mTtt,2))/(s), pT2,s34, k[4], k[5]) ;
	double result = 0;
	if(!fitPDF) result = S2/(2.*pT)*jacs34*(gluon_d0PDF(k[0],  k[1])*partonic[1] + quark_d0PDF(k[0],  k[1])*partonic[0]);
	if(fitPDF) result = S2/(2.*pT)*jacs34*(fit_gluon_d0PDF(k[0],  k[1])*partonic[1] + fit_quark_d0PDF(k[0],  k[1])*partonic[0]);
	if (isnan(result)){return 0;}
	else{ return result;}
}

double vegas_ttH_LO_pT_stt(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double mT = sqrt(pT2+mH2);
	double jacs34 = (S2+mH2-2.*sqrt(S2)*mT-4*mt2);
	double s34 = 4.*mt2+jacs34*k[2];
	double mTtt = sqrt(pT2+s34);
	double s = (k[0]*k[1]*S2);
	//Q2 = pT2 here!
	if(k[0]*k[1]*S2  < Q2) {return 0;}
   	else{
		vector<double> partonic = xsec_LO_pT_stt((pow(mT+mTtt,2))/(s), pT2,s34, k[3], k[4]) ;
		double result = 0;
		if(!fitPDF) result = jacs34*(gluon_d0PDF(k[0],  k[1])*partonic[1] + quark_d0PDF(k[0],  k[1])*partonic[0]);
		if(fitPDF) result = jacs34*(fit_gluon_d0PDF(k[0],  k[1])*partonic[1] + fit_quark_d0PDF(k[0],  k[1])*partonic[0]);
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}

double vegas_ttH_LO_pT_N(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double mT = sqrt(pT2+mH2);
	double mTtt = sqrt(pT2+4.*mt2);
	double s = pow(mT+mTtt,2)/k[1];
	double jacs34 = (s+mH2-2.*sqrt(s)*mT-4*mt2);
	double s34 = 4.*mt2+jacs34*k[2];
	Q2 = pow(mT+mTtt,2); Q = sqrt(Q2);
  double tau = Q2/S2;
	if(k[0] == 1){return 0;}
	else if(tau > 1){return 0;}
   	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<double> partonic = xsec_LO_pT(k[1], pT2,s34, k[3], k[4]) ;
		double result = 2*jacs34*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_mellin_pdf_sum_gg(N)*partonic[1] + fit_mellin_pdf_sum_qqbar(N)*partonic[0])
					*pow(k[1], N - 1.)*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}


double vegas_ttH_LO_pT_stt_N(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double mT = sqrt(pT2+mH2);
	double jacs34 = (S2+mH2-2.*sqrt(S2)*mT-4*mt2);
	double s34 = 4.*mt2+jacs34*k[2];
	double mTtt = sqrt(pT2+s34);
	Q2 = pow(mT+mTtt,2); Q = sqrt(Q2);
  double tau = Q2/S2;
	if(k[0] == 1){return 0;}
	else if(tau > 1){return 0;}
   	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<double> partonic = xsec_LO_pT_stt(k[1], pT2,s34, k[3], k[4]) ;
		double result = 2*jacs34*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_mellin_pdf_sum_gg(N)*partonic[1] + fit_mellin_pdf_sum_qqbar(N)*partonic[0])
					*pow(k[1], N - 1.)*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}


double vegas_ttH_LO_pT_stt_defN(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double mT = sqrt(pT2+mH2);
	double jacs34 = (S2+mH2-2.*sqrt(S2)*mT-4*mt2);
	double s34 = 4.*mt2+jacs34*k[2];
	double mTtt = sqrt(pT2+s34);
	Q2 = pow(mT+mTtt,2); Q = sqrt(Q2);
  double tau = Q2/S2;
	double xvar = k[1]-1.;
	double wr = xvar/(1.+xvar);
	double wjac = 1./(pow(1.+xvar,2));
	complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
  complex<double> xT2 = exp(wr/N);

	if(k[0] == 1){return 0;}
	else if(tau > 1){return 0;}
	else if(k[1] == 0){return 0;}
   	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<complex<double>> partonic = xsec_LO_pT_stt_c(xT2, pT2,s34, k[3], k[4]) ;
		double result = 2*jacs34*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_mellin_pdf_sum_gg(N)*partonic[1] + fit_mellin_pdf_sum_qqbar(N)*partonic[0])
				*wjac*exp(wr)/N*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}

double vegas_ttH_pT_res(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double mT = sqrt(pT2+mH2);
	double mTtt = sqrt(pT2+4.*mt2);
	Q2 = pow(mT+mTtt,2); Q = sqrt(Q2);
	double s = Q2/k[1];
	double jacs34 = (s+mH2-2.*sqrt(s)*mT-4*mt2);
	double s34 = 4.*mt2+jacs34*k[2];
  double tau = Q2/S2;
	if(k[0] == 1){return 0;}
	else if(tau > 1){return 0;}
   	else{
		complex<double> Nint = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		complex<double> Njac = 1./pow(1.-k[0],2)*exp(phiMP*I);
		vector<complex<double>> partonic = xsec_pT_res(Nint*exp(M_gammaE*INCEULER), k[1], pT2,s34, k[3], k[4]);
		double result = 2*jacs34*imag(1./(2.*M_PI)*Njac*
				(fit_mellin_pdf_sum_gg(Nint)*partonic[1] + fit_mellin_pdf_sum_qqbar(Nint)*partonic[0])
					*pow(k[1], Nint - 1.)*pow(tau, -Nint));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}

double vegas_ttH_pT_res_abs(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double mT = sqrt(pT2+mH2);
	double mTtt = sqrt(pT2+4.*mt2);
	Q2 = pow(mT+mTtt,2); Q = sqrt(Q2);
	double tau = Q2/S2;
	double s = Q2/k[1];
	double jacs34 = (s+mH2-2.*sqrt(s)*mT-4*mt2);
	double s34 = 4.*mt2+jacs34*k[2];
  if(k[0] == 1){return 0;}
	else if(tau > 1){return 0;}
   	else{
		complex<double> Nint = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		complex<double> Njac = 1./pow(1.-k[0],2)*exp(phiMP*I);
		vector<complex<double>> partonic = xsec_pT_res_abs(Nint*exp(M_gammaE*INCEULER), k[1], pT2,s34, k[3], k[4]);
		double result = 2*jacs34*imag(1./(2.*M_PI)*Njac*
				(fit_mellin_pdf_sum_gg(Nint)*partonic[1] + fit_mellin_pdf_sum_qqbar(Nint)*partonic[0])
					*pow(k[1], Nint - 1.)*pow(tau, -Nint));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}

double vegas_ttH_pT_stt_res(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double mT = sqrt(pT2+mH2);
	double jacs34 = (S2+mH2-2.*sqrt(S2)*mT-4*mt2);
	double s34 = 4.*mt2+jacs34*k[2];
	double mTtt = sqrt(pT2+s34);
	Q2 = pow(mT+mTtt,2); Q = sqrt(Q2);
  double tau = Q2/S2;
	double xvar = k[1]-1.;
	double wr = xvar/(1.+xvar);
	double wjac = 1./(pow(1.+xvar,2));
	complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
  complex<double> xT2 = exp(wr/N);

	if(k[0] == 1){return 0;}
	else if(tau > 1){return 0;}
	else if(k[1] == 0){return 0;}

   	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<complex<double>> partonic = xsec_pT_res_stt_c(N*exp(M_gammaE*INCEULER), xT2, pT2,s34, k[3], k[4]);
		double result = 2*jacs34*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_mellin_pdf_sum_gg(N)*partonic[1] + fit_mellin_pdf_sum_qqbar(N)*partonic[0])
					*wjac*exp(wr)/N*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}


double vegas_ttH_pT_res_Nfix(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);

	double mT = sqrt(pT2+mH2);
	double mTtt = sqrt(pT2+4.*mt2);
	Q2 = pow(mT+mTtt,2); Q = sqrt(Q2);
	double s = Q2/k[1];
	double jacs34 = (s+mH2-2.*sqrt(s)*mT-4*mt2);
	double s34 = 4.*mt2+jacs34*k[2];
  double tau = Q2/S2;
	if(k[0] == 1){return 0;}
	else if(tau > 1){return 0;}
   	else{
		complex<double> Nint = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		complex<double> Njac = 1./pow(1.-k[0],2)*exp(phiMP*I);
		vector<complex<double>> partonic = xsec_pT_res_Nfix(Nint*exp(M_gammaE*INCEULER), k[1], pT2,s34, k[3], k[4]);
		double result = 2*jacs34*imag(1./(2.*M_PI)*Njac*
				(fit_mellin_pdf_sum_gg(Nint)*partonic[1] + fit_mellin_pdf_sum_qqbar(Nint)*partonic[0])
					*pow(k[1], Nint - 1.)*pow(tau, -Nint));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}

double kallen(double x, double y, double z){
return (pow(x,2)+pow(y,2)+pow(z,2)-2*x*y-2*x*z-2*y*z);
}

std::complex<double> kallen_c(std::complex<double> x, double y, double z){
return (pow(x,2)+pow(y,2)+pow(z,2)-2*x*y-2*x*z-2*y*z);
}

std::vector<double> xsec_div(double rho, double s34, double thetaHCM, double thetat, double phit){
	std::vector<double> xsec = xsec_LO(rho, s34, thetaHCM, thetat, phit);
	std::vector<double> res = {0,0};
	double eps = xsec[0]*1.E-3;
	std::vector<double> xsecpeps = xsec_LO(rho+eps, s34, thetaHCM, thetat, phit);
	std::vector<double> xsecmeps = xsec_LO(rho-eps, s34, thetaHCM, thetat, phit);
	res[0] = (xsecpeps[0]-xsecmeps[0])/(2.*eps);
	res[1] = (xsecpeps[1]-xsecmeps[1])/(2.*eps);
	return res;
}

template <class T>
std::vector<T> get_invariants(T s, T s34, T thetaHCM, T thetat, T phit){
  vector<T> invariant{0,0,0,0};
	T sqrts = sqrt(s);
	T pt_tt = sqrt(s34/4.-mt2);
	T EttCM = (sqrts*sqrts+s34-mH2)/2./sqrts;
	T pttCM = sqrt(EttCM*EttCM-s34);

	//boost matrix
	T gamma = EttCM/sqrt(s34);
	T betax = 0.;
	T betay = -pttCM*sin(thetaHCM)/EttCM;
	T betaz = -pttCM*cos(thetaHCM)/EttCM;
	T beta2 = pow(betay,2)+pow(betaz,2);
	T Et_tt = sqrt(s34)/2;
	T Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
	T ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
	T Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
	T ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);

	//invariants
	invariant[0] = mt2-sqrts*Et_CM+sqrts*ptz_CM; //t13
	invariant[1] = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM; //t14
	invariant[2] = mt2-sqrts*Et_CM-sqrts*ptz_CM; //t23
	invariant[3] = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM; //t24
	return invariant;
}

std::vector<double> xsec_LO(double rho, double s34, double thetaHCM, double thetat, double phit){

double s = Q2/rho;
double sqrts = sqrt(s);
double pt_tt = sqrt(s34/4.-mt2);
double EttCM = (sqrts*sqrts+s34-mH2)/2./sqrts;
double pttCM = sqrt(EttCM*EttCM-s34);

//boost matrix
double gamma = EttCM/sqrt(s34);
double betax = 0.;
double betay = -pttCM*sin(thetaHCM)/EttCM;
double betaz = -pttCM*cos(thetaHCM)/EttCM;
double beta2 = pow(betay,2)+pow(betaz,2);
double Et_tt = sqrt(s34)/2;
double Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
double ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
double Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
double ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);

//invariants
double t13 = mt2-sqrts*Et_CM+sqrts*ptz_CM;
double t14 = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM;
double t23 = mt2-sqrts*Et_CM-sqrts*ptz_CM;
double t24 = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM;

double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<double> result = {0,0};
double s34kal = pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8.*s34);
double skal = pow(kallen(s, s34, pow(mH,2)),0.5)/(8*s);
result[0] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(s34kal)
																*(skal)
																*(SH_qq_LO(s, t13, t14, t23, t24)/36.)
																*sin(thetaHCM)*sin(thetat);
result[1] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(s34kal)
																*(skal)
																*(SH_gg_LO(s, t13, t14, t23, t24)/256.)
																*sin(thetaHCM)*sin(thetat);
return result;

}


std::vector<double> xsec_LO_pT(double xT2, double pT2, double s34, double thetat, double phit){

double mT = sqrt(pT2+mH2);
double mTtt = sqrt(pT2+4.*mt2);
double s = pow(mT+mTtt,2)/xT2;
double sqrts = sqrt(s);
double cosheta = (s+mH2-s34)/2./sqrts/mT; //cosh eta in CM frame
double sinheta = sqrt(pow(cosheta,2)-1.);
	if(isnan(sinheta)){sinheta = 0;}
	if(pow(cosheta,2)-1. < 1.E-15){sinheta = 0;}
double pz5 = mT*sinheta; //pZ of Higgs boson in CM frame BUT! There are two branches: eta < 0 and eta > 0, so multiply result by 2
double E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;

//boost matrix
double gamma = E34CM/sqrt(s34);
double betax = 0.;
double betay = -sqrt(pT2)/E34CM;
double betaz = -pz5/E34CM;
double beta2 = pow(betay,2)+pow(betaz,2);
double Et_tt = sqrt(s34)/2;
double pt_tt = sqrt(s34/4.-mt2); //in ttbar CM frame we have Etop = sqrt(stt)/2 so ptop = sqrt(Etop^2-mt^2)
double Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
double ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
double Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
double ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);

//invariants
double t13 = mt2-sqrts*Et_CM+sqrts*ptz_CM;
double t14 = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM;
double t23 = mt2-sqrts*Et_CM-sqrts*ptz_CM;
double t24 = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM;

double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<double> result = {0,0};
result[0] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																*(SH_qq_LO(s, t13, t14, t23, t24)/36.) // = 1./(4*NC*NC)
																*sin(thetat);
result[1] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																*(SH_gg_LO(s, t13, t14, t23, t24)/256.) // = 1./(4*(NC^2-1)*(NC^2-1))
																*sin(thetat);
return result;

}


std::vector<double> xsec_LO_pT_stt(double xT2, double pT2, double s34, double thetat, double phit){

double mT = sqrt(pT2+mH2);
double mTtt = sqrt(pT2+s34);
double s = pow(mT+mTtt,2)/xT2;
double sqrts = sqrt(s);
double cosheta = (s+mH2-s34)/2./sqrts/mT; //cosh eta in CM frame
double sinheta = sqrt(pow(cosheta,2)-1.);
	if(isnan(sinheta)){sinheta = 0;}
	if(pow(cosheta,2)-1. < 1.E-15){sinheta = 0;}
double pz5 = mT*sinheta; //pZ of Higgs boson in CM frame BUT! There are two branches: eta < 0 and eta > 0, so multiply result by 2
double E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;

//boost matrix
double gamma = E34CM/sqrt(s34);
double betax = 0.;
double betay = -sqrt(pT2)/E34CM;
double betaz = -pz5/E34CM;
double beta2 = pow(betay,2)+pow(betaz,2);
double Et_tt = sqrt(s34)/2;
double pt_tt = sqrt(s34/4.-mt2); //in ttbar CM frame we have Etop = sqrt(stt)/2 so ptop = sqrt(Etop^2-mt^2)
double Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
double ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
double Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
double ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);

//invariants
double t13 = mt2-sqrts*Et_CM+sqrts*ptz_CM;
double t14 = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM;
double t23 = mt2-sqrts*Et_CM-sqrts*ptz_CM;
double t24 = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM;

double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<double> result = {0,0};
result[0] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(SH_qq_LO(s, t13, t14, t23, t24)/36.)
																													*sin(thetat);
result[1] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(SH_gg_LO(s, t13, t14, t23, t24)/256.)
																													*sin(thetat);
return result;

}


int try_stt(double xT2, double pT2, double s34, double thetat, double phit){

	double mT = sqrt(pT2+mH2);
	double mTtt = sqrt(pT2+s34);
	double s = pow(mT+mTtt,2)/xT2;
	double sqrts = sqrt(s);
	double cosheta = (s+mH2-s34)/2./sqrts/mT; //cosh eta in CM frame
	double sinheta = sqrt(pow(cosheta,2)-1.);
	if(isnan(sinheta)){sinheta = 0;}
	if(pow(cosheta,2)-1. < 1.E-15){sinheta = 0;}
cout << "===========================" << endl;
//cout << cosheta << " " << pow(cosheta, 2) << " " << (pow(cosheta,2)-1.)<< " " << sinheta << " " << sqrt(pT2) << endl;
double pz5 = mT*sinheta; //pZ of Higgs boson in CM frame
double E34CM = (s+s34-mH2)/2./sqrts;
double coshetap = (s+s34-mH2)/2./sqrts/mTtt;
double sinhetap = sqrt(pow(coshetap,2)-1.);
double pz5p = mTtt*sinhetap;

//cout << " pz5 " << pz5 << " pz5p " << pz5p << endl;

//boost matrix
double gamma = E34CM/sqrt(s34);
double betax = 0.;
double betay = -sqrt(pT2)/E34CM;
double betaz = -pz5/E34CM;
double beta2 = pow(betay,2)+pow(betaz,2);
double Et_tt = sqrt(s34)/2;
double pt_tt = sqrt(s34/4.-mt2); //in ttbar CM frame we have Etop = sqrt(stt)/2 so ptop = sqrt(Etop^2-mt^2)
double Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
double ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
double Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
double ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
//cout << "/.mh2->" << mH2 << "/.mt2->" <<mt2 << endl;
//cout << "/.gamma->" << gamma << "/.betaz->" << betaz << "/.betay->" <<betay << "/.beta2->" << beta2 << endl;
//cout << "/.Ets->" << Et_tt << "/.pts->" << pt_tt << endl;
double Et_CMT = gamma*Et_tt - gamma*betay*pt_tt*sin(phit)*sin(thetat);
double ptz_CMT = pt_tt*cos(thetat);
double Etbar_CMT = gamma*Et_tt + gamma*betay*pt_tt*sin(phit)*sin(thetat);
double ptbarz_CMT = -pt_tt*cos(thetat);
//invariants
//cout << "/.mt2->" << mt2 << "/.s->" << s <<"/.Et->"<<Et_CM<<"/.ptz->"<<ptz_CM << endl;
//cout << "/.mt2->" << mt2 << "/.s->" << s <<"/.Etb->"<<Etbar_CM<<"/.ptbz->"<<ptbarz_CM << endl;
double t13 = mt2-sqrts*Et_CM+sqrts*ptz_CM;
double t14 = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM;
double t23 = mt2-sqrts*Et_CM-sqrts*ptz_CM;
double t24 = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM;
double t13T = mt2-sqrts*Et_CMT+sqrts*ptz_CMT;
double t14T = mt2-sqrts*Etbar_CMT+sqrts*ptbarz_CMT;
double t23T = mt2-sqrts*Et_CMT-sqrts*ptz_CMT;
double t24T = mt2-sqrts*Etbar_CMT-sqrts*ptbarz_CMT;
cout << "t13 " << t13 << " t23 " << t23 << " t14 " << t14 << " t24 " << t24 << endl;
cout << "t13T " << t13T << " t23T " << t23T << " t14T " << t14T << " t24T " << t24T << endl;
//cout << sqrts << " " << sqrt(s34) << " " << cosheta << " " << sinheta << " " << pow(cosheta,2)-pow(sinheta,2)<< endl;
cout << " lambda qq 11" << lambda_qq_11(s34)/(alphas_muR/(2.*M_PI)) << endl
		<< " lambda qq 12" <<lambda_qq_12(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s)/(alphas_muR/(2.*M_PI))  << endl
		<< " lambda qq 21" <<  lambda_qq_21(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s) /(alphas_muR/(2.*M_PI)) << endl
		<< " lambda qq 22" << lambda_qq_22(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s34, s)/(alphas_muR/(2.*M_PI))  << endl;;

cout  << " lambda gg 11" << lambda_gg_11(s34)/(alphas_muR/(2.*M_PI))  << endl
   << " lambda gg 22" << lambda_gg_22(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s, s34)/(alphas_muR/(2.*M_PI))  << endl
	 << " lambda gg 23" << lambda_gg_32(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s) /(alphas_muR/(2.*M_PI)) << endl
	 << " lambda gg 31 " << lambda_gg_31(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s)/(alphas_muR/(2.*M_PI))  << endl
	 << " lambda gg 32 " << lambda_gg_32(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s)/(alphas_muR/(2.*M_PI))  << endl
	 << " lambda gg 33 " << lambda_gg_33(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s, s34)/(alphas_muR/(2.*M_PI))  << endl;

	 cout << "Lambda " << Lambda3(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s) << endl;
	 cout << "log(1.+pT2/4mt2)-1 " << log(1.+pT2/(4.*mt2))-1. << endl;
	 cout << "Omega " << Omega3(t13-mt2, t24-mt2, t14-mt2, t23-mt2, s) << endl;

	 complex<double> N = 0.5;
	 cout << "qq " << pT_qq_res(N, s, t13, t14, t23, t24, s34) << " " << pT_qq_res_abs(N, pT2, s, t13, t14, t23, t24, s34) << endl;
   cout << "gg " << pT_gg_res(N, s, t13, t14, t23, t24, s34) << " " << pT_gg_res_abs(N, pT2, s, t13, t14, t23, t24, s34) << endl;
return 0;

}

std::vector<complex<double>> xsec_LO_c(complex<double> rho, double s34, double thetaHCM, double thetat, double phit){

complex<double> s = Q2/rho;
complex<double> sqrts = sqrt(s);
complex<double> E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;
complex<double> p5CM = sqrt(E34CM*E34CM-s34);
//boost matrix
complex<double> gamma = E34CM/sqrt(s34);
complex<double> betax = 0.;
complex<double> betay = p5CM*sin(thetaHCM)/E34CM;
complex<double> betaz = p5CM*cos(thetaHCM)/E34CM;
complex<double> beta2 = pow(betay,2)+pow(betaz,2);
complex<double> Et_tt = sqrt(s34)/2.;
complex<double> pt_tt = sqrt(s34/4.-mt2);
complex<double> Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
complex<double> ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
complex<double> Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
complex<double> ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);

//invariants
complex<double> t13 = mt2-sqrts*Et_CM+sqrts*ptz_CM;
complex<double> t14 = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM;
complex<double> t23 = mt2-sqrts*Et_CM-sqrts*ptz_CM;
complex<double> t24 = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM;

double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<complex<double>> result = {0,0};
    result[0] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*
					(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))*
						(pow(kallen_c(s, s34, pow(mH,2)),0.5)/(8*s))*
						(SH_qq_LO_c(s, t13, t14, t23, t24)/36.)*sin(thetaHCM)*sin(thetat);
	result[1] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*
					(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))*
						(pow(kallen_c(s, s34, pow(mH,2)),0.5)/(8*s))*
						(SH_gg_LO_c(s, t13, t14, t23, t24)/256.)*sin(thetaHCM)*sin(thetat);
 	return result;
}



std::vector<complex<double>> xsec_LO_pT_stt_c(complex<double> xT2, double pT2, double s34, double thetat, double phit){
double mT = sqrt(pT2+mH2);
double mTtt = sqrt(pT2+s34);
complex<double> s = pow(mT+mTtt,2)/xT2;
complex<double> sqrts = sqrt(s);
complex<double> cosheta = (s+mH2-s34)/2./sqrts/mT; //cosh eta in CM frame
complex<double> sinheta = sqrt(pow(cosheta,2)-1.);
//	if(isnan(real(sinheta))){sinheta = 0;}
//	if(real(pow(cosheta,2)-1.) < 1.E-15){sinheta = 0;}
complex<double> pz5 = mT*sinheta; //pZ of Higgs boson in CM frame BUT! There are two branches: eta < 0 and eta > 0, so multiply result by 2
complex<double> E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;

//boost matrix
complex<double> gamma = E34CM/sqrt(s34);
complex<double> betax = 0.;
complex<double> betay = -sqrt(pT2)/E34CM;
complex<double> betaz = -pz5/E34CM;
complex<double> beta2 = pow(betay,2)+pow(betaz,2);
complex<double> Et_tt = sqrt(s34)/2.;
complex<double> pt_tt = sqrt(s34/4.-mt2);
complex<double> Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
complex<double> ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
complex<double> Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
complex<double> ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);

//invariants
complex<double> t13 = mt2-sqrts*Et_CM+sqrts*ptz_CM;
complex<double> t14 = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM;
complex<double> t23 = mt2-sqrts*Et_CM-sqrts*ptz_CM;
complex<double> t24 = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM;

complex<double> coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<complex<double>> result = {0,0};
result[0] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(SH_qq_LO_c(s, t13, t14, t23, t24)/36.)
																													*sin(thetat);
result[1] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(SH_gg_LO_c(s, t13, t14, t23, t24)/256.)
																													*sin(thetat);


return result;

}


std::vector<complex<double>> xsec_pT_res(complex<double> N, double xT2, double pT2, double s34, double thetat, double phit){

double mT = sqrt(pT2+mH2);
double mTtt = sqrt(pT2+4.*mt2);
double s = pow(mT+mTtt,2)/xT2;
double sqrts = sqrt(s);
double cosheta = (s+mH2-s34)/2./sqrts/mT; //cosh eta in CM frame
double sinheta = sqrt(pow(cosheta,2)-1.);
	if(isnan(sinheta)){sinheta = 0;}
	if(pow(cosheta,2)-1. < 1.E-15){sinheta = 0;}
double pz5 = mT*sinheta; //pZ of Higgs boson in CM frame BUT! There are two branches: eta < 0 and eta > 0, so multiply result by 2
double E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;

//boost matrix
double gamma = E34CM/sqrt(s34);
double betax = 0.;
double betay = -sqrt(pT2)/E34CM;
double betaz = -pz5/E34CM;
double beta2 = pow(betay,2)+pow(betaz,2);
double Et_tt = sqrt(s34)/2;
double pt_tt = sqrt(s34/4.-mt2); //in ttbar CM frame we have Etop = sqrt(stt)/2 so ptop = sqrt(Etop^2-mt^2)
double Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
double ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
double Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
double ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);

//invariants
double t13 = mt2-sqrts*Et_CM+sqrts*ptz_CM;
double t14 = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM;
double t23 = mt2-sqrts*Et_CM-sqrts*ptz_CM;
double t24 = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM;
//cout << s << " " << t13<< " " << t14<< " " <<  t23<< " " <<  t24<< " " <<  s34 << endl;
double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<complex<double>> result = {0,0};
result[0] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(pT_qq_res(N+1., s, t13, t14, t23, t24,s34)/36.)
																													*sin(thetat);
result[1] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(pT_gg_res(N+1., s, t13, t14, t23, t24,s34)/256.)
																													*sin(thetat);
return result;

}

std::vector<complex<double>> xsec_pT_res_Nfix(complex<double> N, double xT2, double pT2, double s34, double thetat, double phit){

double mT = sqrt(pT2+mH2);
double mTtt = sqrt(pT2+4.*mt2);
double s = pow(mT+mTtt,2)/xT2;
double sqrts = sqrt(s);
double cosheta = (s+mH2-s34)/2./sqrts/mT; //cosh eta in CM frame
double sinheta = sqrt(pow(cosheta,2)-1.);
	if(isnan(sinheta)){sinheta = 0;}
	if(pow(cosheta,2)-1. < 1.E-15){sinheta = 0;}
double pz5 = mT*sinheta; //pZ of Higgs boson in CM frame BUT! There are two branches: eta < 0 and eta > 0, so multiply result by 2
double E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;

//boost matrix
double gamma = E34CM/sqrt(s34);
double betax = 0.;
double betay = -sqrt(pT2)/E34CM;
double betaz = -pz5/E34CM;
double beta2 = pow(betay,2)+pow(betaz,2);
double Et_tt = sqrt(s34)/2;
double pt_tt = sqrt(s34/4.-mt2); //in ttbar CM frame we have Etop = sqrt(stt)/2 so ptop = sqrt(Etop^2-mt^2)
double Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
double ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
double Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
double ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);

//invariants
double t13 = mt2-sqrts*Et_CM+sqrts*ptz_CM;
double t14 = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM;
double t23 = mt2-sqrts*Et_CM-sqrts*ptz_CM;
double t24 = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM;

double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<complex<double>> result = {0,0};
double Ngluon = sgg*exp(INCEULER*M_gammaE);//newton_raphson_gg(1., 0., tau, xT2);
double Nquark = sqqbar*exp(INCEULER*M_gammaE);//newton_raphson_qqbar(1., 0., tau, xT2);
/*cout << Ngluon << endl;
cout << Nquark << endl;*/
result[0] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(pT_qq_res(Nquark, s, t13, t14, t23, t24,s34)/36.)
																													*sin(thetat);
result[1] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(pT_gg_res(Ngluon, s, t13, t14, t23, t24,s34)/256.)
																													*sin(thetat);
return result;

}

std::vector<complex<double>> xsec_pT_res_abs(complex<double> N, double xT2, double pT2, double s34, double thetat, double phit){

double mT = sqrt(pT2+mH2);
double mTtt = sqrt(pT2+4.*mt2);
double s = pow(mT+mTtt,2)/xT2;
double sqrts = sqrt(s);
double cosheta = (s+mH2-s34)/2./sqrts/mT; //cosh eta in CM frame
double sinheta = sqrt(pow(cosheta,2)-1.);
	if(isnan(sinheta)){sinheta = 0;}
	if(pow(cosheta,2)-1. < 1.E-15){sinheta = 0;}
double pz5 = mT*sinheta; //pZ of Higgs boson in CM frame BUT! There are two branches: eta < 0 and eta > 0, so multiply result by 2
double E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;

//boost matrix
double gamma = E34CM/sqrt(s34);
double betax = 0.;
double betay = -sqrt(pT2)/E34CM;
double betaz = -pz5/E34CM;
double beta2 = pow(betay,2)+pow(betaz,2);
double Et_tt = sqrt(s34)/2;
double pt_tt = sqrt(s34/4.-mt2); //in ttbar CM frame we have Etop = sqrt(stt)/2 so ptop = sqrt(Etop^2-mt^2)
double Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
double ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
double Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
double ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);

//invariants
double t13 = mt2-sqrts*Et_CM+sqrts*ptz_CM;
double t14 = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM;
double t23 = mt2-sqrts*Et_CM-sqrts*ptz_CM;
double t24 = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM;

double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<complex<double>> result = {0,0};

result[0] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(pT_qq_res_abs(N+1., pT2, s, t13, t14, t23, t24,s34)/36.)
																													*sin(thetat);
result[1] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(pT_gg_res_abs(N+1., pT2, s, t13, t14, t23, t24,s34)/256.)
																													*sin(thetat);
return result;

}

std::vector<complex<double>> xsec_pT_res_stt_c(complex<double> N, complex<double> xT2, double pT2, double s34, double thetat, double phit){

	double mT = sqrt(pT2+mH2);
	double mTtt = sqrt(pT2+s34);
	complex<double> s = pow(mT+mTtt,2)/xT2;
	complex<double> sqrts = sqrt(s);
	complex<double> cosheta = (s+mH2-s34)/2./sqrts/mT; //cosh eta in CM frame
	complex<double> sinheta = sqrt(pow(cosheta,2)-1.);
	//	if(isnan(real(sinheta))){sinheta = 0;}
	//	if(real(pow(cosheta,2)-1.) < 1.E-15){sinheta = 0;}
	complex<double> pz5 = mT*sinheta; //pZ of Higgs boson in CM frame BUT! There are two branches: eta < 0 and eta > 0, so multiply result by 2
	complex<double> E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;

	//boost matrix
	complex<double> gamma = E34CM/sqrt(s34);
	complex<double> betax = 0.;
	complex<double> betay = -sqrt(pT2)/E34CM;
	complex<double> betaz = -pz5/E34CM;
	complex<double> beta2 = pow(betay,2)+pow(betaz,2);
	complex<double> Et_tt = sqrt(s34)/2.;
	complex<double> pt_tt = sqrt(s34/4.-mt2);
	complex<double> Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
	complex<double> ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
	complex<double> Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
	complex<double> ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);

	//invariants
	complex<double> t13 = mt2-sqrts*Et_CM+sqrts*ptz_CM;
	complex<double> t14 = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM;
	complex<double> t23 = mt2-sqrts*Et_CM-sqrts*ptz_CM;
	complex<double> t24 = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM;

	double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
	std::vector<complex<double>> result = {0,0};

	result[0] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																														*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																														*(pT_qq_res(N+1., s, t13, t14, t23, t24,s34)/36.)
																														*sin(thetat);
	result[1] = 2.*coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																														*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																														*(pT_gg_res(N+1., s, t13, t14, t23, t24,s34)/256.)
																														*sin(thetat);
	return result;

}

std::vector<complex<double>> xsec_res(complex<double> N, double rho, double s34, double thetaHCM, double thetat, double phit){
double s = Q2/rho;
double sqrts = sqrt(s);
double E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;
double p5CM = sqrt(E34CM*E34CM-s34);

//boost matrix
double gamma = E34CM/sqrt(s34);
double betax = 0.;
double betay = p5CM*sin(thetaHCM)/E34CM;
double betaz = p5CM*cos(thetaHCM)/E34CM;
double beta2 = pow(betay,2)+pow(betaz,2);
double Et_tt = sqrt(s34)/2;
double pt_tt = sqrt(s34/4.-mt2);
double Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
double ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
double Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
double ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);

//invariants
double t13 = mt2-sqrts*Et_CM+sqrts*ptz_CM;
double t14 = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM;
double t23 = mt2-sqrts*Et_CM-sqrts*ptz_CM;
double t24 = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM;

double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<complex<double>> result = {0,0};
result[0] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*
					(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))*
						(pow(kallen_c(s, s34, pow(mH,2)),0.5)/(8*s))*
						(SH_qq_res(N+1., s, t13, t14, t23, t24)/36.)*sin(thetaHCM)*sin(thetat);
result[1] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*
					(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))*
						(pow(kallen_c(s, s34, pow(mH,2)),0.5)/(8*s))*
						(SH_gg_res(N+1., s, t13, t14, t23, t24)/256.)*sin(thetaHCM)*sin(thetat);
 	return result;
}


std::vector<complex<double>> xsec_res_c(complex<double> N, complex<double> rho, double s34, double thetaHCM, double thetat, double phit){

complex<double> s = Q2/rho;
complex<double> sqrts = sqrt(s);
complex<double> E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;
complex<double> p5CM = sqrt(E34CM*E34CM-s34);

//boost matrix
complex<double> gamma = E34CM/sqrt(s34);
complex<double> betax = 0.;
complex<double> betay = p5CM*sin(thetaHCM)/E34CM;
complex<double> betaz = p5CM*cos(thetaHCM)/E34CM;
complex<double> beta2 = pow(betay,2)+pow(betaz,2);
complex<double> Et_tt = sqrt(s34)/2.;
complex<double> pt_tt = sqrt(s34/4.-mt2);
complex<double> Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
complex<double> ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
complex<double> Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
complex<double> ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);

//invariants
complex<double> t13 = mt2-sqrts*Et_CM+sqrts*ptz_CM;
complex<double> t14 = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM;
complex<double> t23 = mt2-sqrts*Et_CM-sqrts*ptz_CM;
complex<double> t24 = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM;

double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<complex<double>> result = {0,0};
    result[0] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*
					(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))*
						(pow(kallen_c(s, s34, pow(mH,2)),0.5)/(8*s))*
						(SH_qq_res_c(N+1., s, t13, t14, t23, t24)/36.)*sin(thetaHCM)*sin(thetat);
	result[1] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*
					(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))*
						(pow(kallen_c(s, s34, pow(mH,2)),0.5)/(8*s))*
						(SH_gg_res_c(N+1., s, t13, t14, t23, t24)/256.)*sin(thetaHCM)*sin(thetat);
 	return result;
}




std::vector<complex<double>> xsec_res_c_z5(complex<double> N, complex<double> rho, double s34, double thetaHCM, double thetat, double phit){

complex<double> s = Q2/rho;
complex<double> sqrts = sqrt(s);
complex<double> E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;
complex<double> p5CM = sqrt(E34CM*E34CM-s34);
//boost matrix
complex<double> gamma = E34CM/sqrt(s34);
complex<double> betax = 0.;
complex<double> betay = p5CM*sin(thetaHCM)/E34CM;
complex<double> betaz = p5CM*cos(thetaHCM)/E34CM;
complex<double> beta2 = pow(betay,2)+pow(betaz,2);
complex<double> Et_tt = sqrt(s34)/2.;
complex<double> pt_tt = sqrt(s34/4.-mt2);
complex<double> Et_CM = gamma*Et_tt - gamma*betax*pt_tt*cos(phit)*sin(thetat) - gamma*betay*pt_tt*sin(phit)*sin(thetat) - gamma*betaz*pt_tt*cos(thetat);
complex<double> ptz_CM = -gamma*betaz*Et_tt + (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) + (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) + ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);
complex<double> Etbar_CM = gamma*Et_tt + gamma*betax*pt_tt*cos(phit)*sin(thetat) + gamma*betay*pt_tt*sin(phit)*sin(thetat) + gamma*betaz*pt_tt*cos(thetat);
complex<double> ptbarz_CM = -gamma*betaz*Et_tt - (gamma-1.)*betax*betaz/beta2*pt_tt*cos(phit)*sin(thetat) - (gamma-1.)*betay*betaz/beta2*pt_tt*sin(phit)*sin(thetat) - ((gamma-1.)*betaz*betaz/beta2+1.)*pt_tt*cos(thetat);

//invariants
complex<double> t13 = mt2-sqrts*Et_CM+sqrts*ptz_CM;
complex<double> t14 = mt2-sqrts*Etbar_CM+sqrts*ptbarz_CM;
complex<double> t23 = mt2-sqrts*Et_CM-sqrts*ptz_CM;
complex<double> t24 = mt2-sqrts*Etbar_CM-sqrts*ptbarz_CM;

double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<complex<double>> result = {0,0};
    result[0] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*
					(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))*
						(pow(kallen_c(s, s34, pow(mH,2)),0.5)/(8*s))*
						(omega_qq_res(N+1., s, t13, t14, t23, t24,s34)/36.)*sin(thetaHCM)*sin(thetat);
	result[1] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*
					(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))*
						(pow(kallen_c(s, s34, pow(mH,2)),0.5)/(8*s))*
						(omega_gg_res(N+1., s, t13, t14, t23, t24,s34)/256.)*sin(thetaHCM)*sin(thetat);
 	return result;
}
