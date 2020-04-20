#include "k_factors_ttH.h"
#include "tth_softanom.h"
#include "tth_vegas.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "parameters.h"
#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

using namespace std;

double vegas_ttH_Nspace_deform(double *k, size_t dim, void *params)
{
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

double vegas_ttH_Nspace(double *k, size_t dim, void *params)
{
	(void)(dim);
	(void)(params);
	//N,	rho	s34, 	angles,
	//k[0]	k[1]	k[2]	k[3]-k[5],
	//double s34 = 4.*pow(mt,2)+(S2-pow(mH,2)-4*pow(mt,2))*k[2];

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

double vegas_ttH_LO(double *k, size_t dim, void *params)
{
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


double vegas_ttH_LO1st(double *k, size_t dim, void *params)
{
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



double vegas_ttH_LO2nd(double *k, size_t dim, void *params)
{
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


double vegas_ttH_resum_2nd(double *k, size_t dim, void *params)
{
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


double vegas_ttH_resum_2nd_z5(double *k, size_t dim, void *params)
{
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



double vegas_ttH_resum_2nd_def(double *k, size_t dim, void *params)
{
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


double vegas_ttH_resum_2nd_def_z5(double *k, size_t dim, void *params)
{
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


double vegas_ttH_resum_def(double *k, size_t dim, void *params)
{
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


double vegas_ttH_resum_def_z5(double *k, size_t dim, void *params)
{
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


double vegas_ttH_inv_mass_LO(double *k, size_t dim, void *params)
{
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

double vegas_ttH_inv_mass_res(double *k, size_t dim, void *params)
{
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

double vegas_ttH_inv_mass_res_s34(double *k, size_t dim, void *params)
{
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


double vegas_ttH_inv_mass_LO_s34(double *k, size_t dim, void *params)
{
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

double vegas_ttH_LO_pT(double *k, size_t dim, void *params)
{
	(void)(dim);
	(void)(params);

	double mT = sqrt(Q2+mH2);
	double mTtt = sqrt(Q2+4.*mt2);
	double s = (k[0]*k[1]*S2);
	double jacs34 = (s+mH2-2.*sqrt(s)*mT-4*mt2);
	double s34 = 4.*mt2+jacs34*k[2];
	//Q2 = pT2 here!
	if(k[0]*k[1]*S2  < Q2) {return 0;}
   	else{
		vector<double> partonic = xsec_LO_pT((pow(mT+mTtt,2))/(k[0]*k[1]*S2), Q2,s34, k[3], k[4]) ;
		double result = 0;
		if(!fitPDF) result = jacs34*(gluon_d0PDF(k[0],  k[1])*partonic[1] + quark_d0PDF(k[0],  k[1])*partonic[0]);
		if(fitPDF) result = jacs34*(fit_gluon_d0PDF(k[0],  k[1])*partonic[1] + fit_quark_d0PDF(k[0],  k[1])*partonic[0]);
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}


double vegas_ttH_LO_pT_N(double *k, size_t dim, void *params)
{
	(void)(dim);
	(void)(params);

	double mT = sqrt(Q2+mH2);
	double mTtt = sqrt(Q2+4.*mt2);
	double s = pow(mT+mTtt,2)/k[1];
	double jacs34 = (s+mH2-2.*sqrt(s)*mT-4*mt2);
	double s34 = 4.*mt2+jacs34*k[2];
  double tau = pow(mT+mTtt,2)/S2;
	if(k[0] == 1){return 0;}
	else if(tau > 1){return 0;}
   	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<double> partonic = xsec_LO_pT(k[1], Q2,s34, k[3], k[4]) ;
		double result = 2*jacs34*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_mellin_pdf_sum_gg(N)*partonic[1] + fit_mellin_pdf_sum_qqbar(N)*partonic[0])
					*pow(k[1], N - 1.)*pow(tau, -N));
		if (isnan(result)){return 0;}
		else{ return result;}
		}
}


double vegas_ttH_pT_res(double *k, size_t dim, void *params)
{
	(void)(dim);
	(void)(params);

	//Q = 500.*k[5];
	//Q2 = pow(Q,2);

		double mT = sqrt(Q2+mH2);
		double mTtt = sqrt(Q2+4.*mt2);
		double s = pow(mT+mTtt,2)/k[1];
		double jacs34 = (s+mH2-2.*sqrt(s)*mT-4*mt2);
		double s34 = 4.*mt2+jacs34*k[2];
	  double tau = pow(mT+mTtt,2)/S2;
	if(k[0] == 1){return 0;}
	else if(tau > 1){return 0;}
   	else{
		complex<double> N = CMP + k[0]/(1-k[0])*exp(I*phiMP);
		vector<complex<double>> partonic = xsec_pT_res(N, k[1], Q2,s34, k[3], k[4]);
		double result = 2*jacs34*imag(1./(2.*M_PI)*exp(I*phiMP)/pow(1.-k[0], 2.)*
				(fit_mellin_pdf_sum_gg(N)*partonic[1] + fit_mellin_pdf_sum_qqbar(N)*partonic[0])
					*pow(k[1], N - 1.)*pow(tau, -N));
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

std::vector<double> xsec_LO(double rho, double s34, double thetaHCM, double thetat, double phit){

double s = Q2/rho;
double sqrts = sqrt(s);
double p3tt = sqrt(s34/4.-mt2);
double E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;
double p5CM = sqrt(E34CM*E34CM-s34);
double pz5 = p5CM*cos(thetaHCM);
double pz3tt = p3tt*(cos(thetaHCM)*cos(thetat)-sin(thetaHCM)*sin(thetat)*cos(phit));

double t13 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5-sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double t24 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5-sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double t14 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5+sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double t23 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5+sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<double> result = {0,0};
result[0] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))*(pow(kallen(s, s34, pow(mH,2)),0.5)/(8*s))*(SH_qq_LO(s, t13, t14, t23, t24)/36.)*sin(thetaHCM)*sin(thetat);
result[1] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))*(pow(kallen(s, s34, pow(mH,2)),0.5)/(8*s))*(SH_gg_LO(s, t13, t14, t23, t24)/256.)*sin(thetaHCM)*sin(thetat);
return result;

}


std::vector<double> xsec_LO_pT(double xT2, double pT2, double s34, double thetat, double phit){

double mT = sqrt(pT2+mH2);
double mTtt = sqrt(pT2+4.*mt2);
double s = pow(mT+mTtt,2)/xT2;
double sqrts = sqrt(s);
double cosheta = (s+mH2-s34)/2./sqrts/mT; //cosh eta in CM frame
double sinheta = sqrt(pow(cosheta,2)-1.);
double p3tt = sqrt(s34/4.-mt2); //in ttbar CM frame we have Etop = sqrt(stt)/2 so ptop = sqrt(Etop^2-mt^2)
double pz5 = mT*sinheta; //pZ of Higgs boson in CM frame
double p5CM = sqrt(pT2+pow(pz5,2)); // |p| of Higgs boson in CM frame
double E34CM = sqrt(s34+pow(p5CM,2)); // energy of ttbar in CM frame (E^2 - ptt^2 = stt with ptt = pH)
double thetaHCM = acos(pz5/p5CM);
double pz3tt = p3tt*(cos(thetaHCM)*cos(thetat)-sin(thetaHCM)*sin(thetat)*cos(phit));

double t13 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5-sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double t24 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5-sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double t14 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5+sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double t23 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5+sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<double> result = {0,0};
result[0] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(SH_qq_LO(s, t13, t14, t23, t24)/36.)
																													*sin(thetat);
result[1] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(SH_gg_LO(s, t13, t14, t23, t24)/256.)
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
double p3tt = sqrt(s34/4.-mt2); //in ttbar CM frame we have Etop = sqrt(stt)/2 so ptop = sqrt(Etop^2-mt^2)
double pz5 = mT*sinheta; //pZ of Higgs boson in CM frame
double p5CM = sqrt(pT2+pow(pz5,2)); // |p| of Higgs boson in CM frame
double E34CM = sqrt(s34+pow(p5CM,2)); // energy of ttbar in CM frame (E^2 - ptt^2 = stt with ptt = pH)
double thetaHCM = acos(pz5/p5CM);
double pz3tt = p3tt*(cos(thetaHCM)*cos(thetat)-sin(thetaHCM)*sin(thetat)*cos(phit));

double t13 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5-sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double t24 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5-sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double t14 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5+sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double t23 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5+sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double coupling = pow(mt/v,2)*pow(alphas_muR,2)*pow(4*M_PI,2);
std::vector<complex<double>> result = {0,0};
result[0] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(omega_qq_res(N+1., s, t13, t14, t23, t24,s34)/36.)
																													*sin(thetat);
result[1] = coupling*pbunits*(1/(2*s))*(1/(pow(2*M_PI,4)))*(pow(kallen(s34, pow(mt,2), pow(mt,2)), 0.5)/(8*s34))
																													*(sqrt(pT2)/(4.*sqrts)*1./pz5)
																													*(omega_gg_res(N+1., s, t13, t14, t23, t24,s34)/256.)
																													*sin(thetat);
return result;

}

std::vector<complex<double>> xsec_LO_c(complex<double> rho, double s34, double thetaHCM, double thetat, double phit){

complex<double> s = Q2/rho;
complex<double> sqrts = sqrt(s);
double p3tt = sqrt(s34/4.-mt2);
complex<double> E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;
complex<double> p5CM = sqrt(E34CM*E34CM-s34);
complex<double> pz5 = p5CM*cos(thetaHCM);
double pz3tt = p3tt*(cos(thetaHCM)*cos(thetat)-sin(thetaHCM)*sin(thetat)*cos(phit));

complex<double> t13 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5-sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
complex<double> t24 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5-sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
complex<double> t14 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5+sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
complex<double> t23 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5+sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
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


std::vector<complex<double>> xsec_res(complex<double> N, double rho, double s34, double thetaHCM, double thetat, double phit){
double s = Q2/rho;
double sqrts = sqrt(s);
double p3tt = sqrt(s34/4.-mt2);
double E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;
double p5CM = sqrt(E34CM*E34CM-s34);
double pz5 = p5CM*cos(thetaHCM);
double pz3tt = p3tt*(cos(thetaHCM)*cos(thetat)-sin(thetaHCM)*sin(thetat)*cos(phit));

double t13 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5-sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double t24 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5-sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double t14 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5+sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
double t23 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5+sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
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
double p3tt = sqrt(s34/4.-mt2);
complex<double> E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;
complex<double> p5CM = sqrt(E34CM*E34CM-s34);
complex<double> pz5 = p5CM*cos(thetaHCM);
double pz3tt = p3tt*(cos(thetaHCM)*cos(thetat)-sin(thetaHCM)*sin(thetat)*cos(phit));

complex<double> t13 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5-sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
complex<double> t24 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5-sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
complex<double> t14 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5+sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
complex<double> t23 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5+sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
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
double p3tt = sqrt(s34/4.-mt2);
complex<double> E34CM = (sqrts*sqrts+s34-mH2)/2./sqrts;
complex<double> p5CM = sqrt(E34CM*E34CM-s34);
complex<double> pz5 = p5CM*cos(thetaHCM);
double pz3tt = p3tt*(cos(thetaHCM)*cos(thetat)-sin(thetaHCM)*sin(thetat)*cos(phit));

complex<double> t13 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5-sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
complex<double> t24 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5-sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
complex<double> t14 = mt2-E34CM*sqrts/2.-sqrts/2.*pz5+sqrts*pz3tt-2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(-sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
complex<double> t23 = mt2-E34CM*sqrts/2.+sqrts/2.*pz5+sqrts*pz3tt+2.*p3tt/sqrt(s34)*p5CM*cos(thetat)*(sqrts/2.*pz5/(E34CM+sqrt(s34))-sqrts/2);
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
