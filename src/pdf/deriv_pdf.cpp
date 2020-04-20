#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "LHAPDF/LHAPDF.h"
#include "deriv_pdf.h"
#include "parameters.h"
using namespace std;

////////////////////////////////////////////////////////////////////////
///
/// this file contains all pdf related stuff, like sums of qqbar
/// and derivatives of pdfs
/// note that the xfxQ function of LHAPDF returns the x*f(x) value
///
////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////
/// this is the qqbar sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////

double pdf_sum_qqbar_charge_weighted(double x, double tau_over_z){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		sum_pdf+= eq[i-1]*eq[i-1]*1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
		//cout << "in the sum " << sum_pdf << endl;
	}
	return 1./x*sum_pdf;
}
double pdf_sum_qqbar_charge_unweighted(double x, double tau_over_z){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	double charge_factor = 0;
	for(int i = 1; i <=5; i++){
		charge_factor += eq[i-1]*eq[i-1];
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
	}
	return 1./x*charge_factor*sum_pdf;
}
//qq + qbarqbar +qqbar (identical+non-identical!), use it for DeltaNNLOC2
double pdf_sum_qq_charge_weighted_double(double x, double tau_over_z){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	for(int i = -5; i <=5; i++){
		for(int j = -5; j <=5; j++){
			if(j==0){continue;}
			if(i==0){continue;}
			sum_pdf+= eq[abs(i)-1]*eq[abs(i)-1]*1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(j,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(j,x,muF));
		}
	}
	return 1./x*sum_pdf;
}
//qq + qbarqbar +qqbar (identical+non-identical!), use it for DeltaNNLOCE
double pdf_sum_qq_charge_weighted_single(double x, double tau_over_z){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	for(int i = -5; i <=5; i++){
		if(i==0){continue;}
		sum_pdf+= eq[abs(i)-1]*eq[abs(i)-1]*1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(i,x,muF));
	}
	return 1./x*sum_pdf;
}


//qq + qbarqbar + qqbar (identical+nonidentical!), use it for DeltaNNLOCD
double pdf_sum_qq_charge_weighted_double_vivj(double x, double tau_over_z){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	for(int i = -5; i <=5; i++){
		for(int j = -5; j <=5; j++){
			if(j==0){continue;}
			if(i==0){continue;}
			sum_pdf+= eq[abs(i)-1]*eq[abs(j)-1]*1./x*1./(tau_over_z/x)*(-(float) i*j)/((float)abs(i*j))*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(j,tau_over_z/x,muF));
		}
	}
	return 1./x*sum_pdf;
}


//qq + qbarqbar + qqbar (identical only), use it for DeltaNNLOCF
double pdf_sum_qq_charge_weighted_single_vivi(double x, double tau_over_z){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	for(int i = -5; i <=5; i++){
		if(i==0){continue;}
		sum_pdf+= eq[abs(i)-1]*eq[abs(i)-1]*1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(i,tau_over_z/x,muF));
	}
	return 1./x*sum_pdf;
}

/////////////////////////////////////////////////////////////////////////
/// this is the qg sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////

double pdf_sum_qg_charge_weighted(double x, double tau_over_z){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		sum_pdf+= eq[i-1]*eq[i-1]*1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(0,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(0,x,muF));
		sum_pdf+= eq[i-1]*eq[i-1]*1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(-i,x,muF)*pdfs[use_member]->xfxQ(0,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(0,x,muF));
	}
	return 1./x*sum_pdf;
}


/////////////////////////////////////////////////////////////////////////
/// this is the gg sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////

double pdf_sum_gg_charge_weighted(double x, double tau_over_z){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	double charge_factor = 0;
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		charge_factor += eq[i-1]*eq[i-1]; //still need the charge factor
			}
	sum_pdf = 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(0,x,muF)*pdfs[use_member]->xfxQ(0,tau_over_z/x,muF));
	return 1./x*sum_pdf*charge_factor;
}


///////////////////////////////////////
/// HIGGS STUFF
//////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/// this is the qqbar sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////
double pdf_sum_qqbar(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
		//cout << "in the sum " << sum_pdf << endl;
	}
	return 1./x*sum_pdf;
}
/////////////////////////////////////////////////////////////////////////
/// this is the qg sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////
double pdf_sum_qg(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(0,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(0,x,muF));
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(-i,x,muF)*pdfs[use_member]->xfxQ(0,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(0,x,muF));
	}
	return 1./x*sum_pdf;
}
/////////////////////////////////////////////////////////////////////////
/// this is the gg sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////
double pdf_sum_gg(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	sum_pdf = 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(0,x,muF)*pdfs[use_member]->xfxQ(0,tau_over_z/x,muF));
	return 1./x*sum_pdf;
}
/////////////////////////////////////////////////////////////////////////
/// this is the qq (+qbarqbar) sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////
double pdf_sum_qq(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	for(int i = 1; i <=5; i++){
		sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(-i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF));
		//cout << "in the sum " << sum_pdf << endl;
	}
	return 1./x*sum_pdf;
}
/////////////////////////////////////////////////////////////////////////
/// this is the qq (+qbarqbar) sum with normal x and tau/z inputs for the integral
/// int[1/x*f_1(x,Q)*f_2(tau/(z*x),Q),{x,tau/z,1}]
/////////////////////////////////////////////////////////////////////////
double pdf_sum_qqNI(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	//int l = 0;
	for(int i = -5; i <=5; i++){
		for(int j = -5; j <=5; j++){
			if(abs(i)==abs(j)){continue;}
			if(i==0){continue;}
			if(j==0){continue;}
			//l += 1;
			//cout << "f_i(x1) i=" << i << ", f_j(x2), j=" << j << endl;//to avoid double counting!
			sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(j,tau_over_z/x,muF));
			//sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(-i,x,muF)*pdfs[use_member]->xfxQ(-j,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-j,x,muF));
			//sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(-i,x,muF)*pdfs[use_member]->xfxQ(j,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-j,x,muF));
		}
	}
	//cout << l << endl;
	return 1./x*sum_pdf;
}

double pdf_sum_qqbarNI(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(1,x,muF)*pdfs[use_member]->xfxQ(-2,tau_over_z/x,muF));
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(-2,x,muF)*pdfs[use_member]->xfxQ(1,tau_over_z/x,muF));
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(-1,x,muF)*pdfs[use_member]->xfxQ(2,tau_over_z/x,muF));
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(2,x,muF)*pdfs[use_member]->xfxQ(-1,tau_over_z/x,muF));
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(3,x,muF)*pdfs[use_member]->xfxQ(-4,tau_over_z/x,muF));
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(-4,x,muF)*pdfs[use_member]->xfxQ(3,tau_over_z/x,muF));
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(-3,x,muF)*pdfs[use_member]->xfxQ(4,tau_over_z/x,muF));
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(4,x,muF)*pdfs[use_member]->xfxQ(-3,tau_over_z/x,muF));
	return 1./x*sum_pdf;
}

//////////////////////////////////////////////////////
/// relevant for W+W-
//////////////////////////////////////////////////////
double pdf_sum_qqbarUP(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	int i = 2;
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
	i = 4;
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
	return 1./x*sum_pdf;
}
double pdf_sum_qqbarDOWN(double x, double tau_over_z){
	double sum_pdf(0);
	if(x < tau_over_z){return 0;}
	int i = 1;
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
	i = 3;
	sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
	//i = 5;
	//sum_pdf+= 1./x*1./(tau_over_z/x)*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,tau_over_z/x,muF)+pdfs[use_member]->xfxQ(i,tau_over_z/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
	return 1./x*sum_pdf;
}


///////////////////////////////////////////////////////
/// numerical derivative of the weighted pdf sum for qg
///////////////////////////////////////////////////////
double derivative_qg_pdf(double x, double z, double tau, double eps){
	double tau_over_zpeps(tau/(z+eps)); //this is the argument of z + eps
	double tau_over_zmeps(tau/(z-eps)); //argument z - eps
	return (1./(z+eps)*pdf_sum_qg_charge_weighted(x,tau_over_zpeps)-1./(z-eps)*pdf_sum_qg_charge_weighted(x,tau_over_zmeps))/(2.*eps); //numerical derivative
	//return (jac_zpeps/(z+eps)*pdf_sum_charge_weighted(x,tau_over_zpeps)-jac_zmeps/(z)*pdf_sum_charge_weighted(x,tau/z))/(eps); //numerical derivative
}


////////////////////////////////////////////////////////////////////////////////////////
/// numerical derivative of the weighted pdf sum including a z-dependent jacobean for qg
////////////////////////////////////////////////////////////////////////////////////////
double derivative_qg_pdf_jac(double x, double z, double tau, double eps){

	double zpeps = z + eps;
	double zmeps = z - eps;
	double jacp = 1.-tau/zpeps, jacm = 1.-tau/zmeps;
	double xpeps = tau/zpeps+x*jacp, xmeps = tau/zmeps+x*jacm;

	if(xmeps < tau/(zpeps)){return 0;}
	if(xpeps < tau/(zmeps)){return 0;}
	return (jacp/(zpeps)*pdf_sum_qg_charge_weighted(xpeps,tau/(zpeps))-jacm/(zmeps)*pdf_sum_qg_charge_weighted(xmeps,tau/(zmeps)))/(2.*eps); //numerical derivative
}



////////////////////////////////////////////////////////////////////////////////////////
/// numerical derivative of gg (for di-higgs)
////////////////////////////////////////////////////////////////////////////////////////
double derivative_gg_pdf(double x, double tau, double eps){
	eps = 1.E-6;
	double taupeps = tau + eps;
	double taumeps = tau - eps;
	return (pdf_sum_gg(x, taupeps) - pdf_sum_gg(x,taumeps))/(2.*eps);
}


///////////////////////////////////
/// checks for conservation charge
///////////////////////////////////
double vegas_pdf_up_minus_upbar(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	return 1./k[0]*(pdfs[use_member]->xfxQ(2,k[0],muF)-pdfs[use_member]->xfxQ(-2,k[0],muF));

}

///////////////////////////////////
/// checks for momentum cons.
///////////////////////////////////
double vegas_pdf_mom_consv(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	double sum_pdf(0);
	for(int i = 1; i <=5; i++){
		sum_pdf+= (pdfs[use_member]->xfxQ(i,k[0],muF)+pdfs[use_member]->xfxQ(-i,k[0],muF));
	}
	sum_pdf+=pdfs[use_member]->xfxQ(21,k[0],muF);
	return sum_pdf;
}

/////////////////////////////////////
/// lumi checks
/////////////////////////////////////
double vegas_lumi_gg(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	return pdf_sum_gg(k[0],tau);
}
double vegas_lumi_qg(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	return pdf_sum_qg(k[0],tau);
}
double vegas_lumi_qqbar(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	return pdf_sum_qqbar(k[0],tau);
}
double vegas_lumi_qq(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	return pdf_sum_qq(k[0],tau);
}
double vegas_lumi_qqNI(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	return pdf_sum_qqNI(k[0],tau);
}


//////////////////////////////////////////////////////////////
/// derivatives for ttH
//////////////////////////////////////////////////////////////

double single_div_sophis(int i, double x)
{
	double epsilon = 1E-5*x;
	if((x-2.*epsilon ) < 0){return 0;}
	else if((x+2.*epsilon) > 1){return 0;}
	else {return (-pdfs[use_member]->xfxQ(i,x+2.*epsilon,muF) + 8*pdfs[use_member]->xfxQ(i,x+epsilon,muF)-8.*pdfs[use_member]->xfxQ(i,x-epsilon,muF)+pdfs[use_member]->xfxQ(i,x-2.*epsilon,muF))/(12.*epsilon);}
}
double double_div_sophis(int i, double x)
{
	double epsilon = 1E-5*x;
	if((x-2.*epsilon ) < 0){return 0;}
	else if((x+2.*epsilon) > 1){return 0;}
	else {return (-single_div_sophis(i, x + 2.*epsilon) + 8*single_div_sophis(i, x + epsilon)-8.*single_div_sophis(i, x - epsilon)+single_div_sophis(i, x - 2.*epsilon))/(12.*epsilon);}
}
double single_div(int i, double x)
{
	double epsilon = 1E-3*x;
	if((x-epsilon/2. ) < 0){return 0;}
	else if((x+epsilon/2.) > 1){return 0;}
	else {return (pdfs[use_member]->xfxQ(i,x+epsilon/2.,muF) - pdfs[use_member]->xfxQ(i,x-epsilon/2.,muF))/(epsilon);}
}
double double_div(int i, double x)
{
	double epsilon = 1E-6*x;
	if((x-epsilon ) < 0){return 0;}
	else if((x+epsilon) > 1){return 0;}
	else {return (pdfs[use_member]->xfxQ(i,x+epsilon,muF) -2.*pdfs[use_member]->xfxQ(i,x,muF)+pdfs[use_member]->xfxQ(i,x-epsilon,muF))/(epsilon*epsilon);}

}
double gluon_d0PDF(double x1, double x2)
{
	//no div
	return (pdfs[use_member]->xfxQ(0,x1,muF))*(pdfs[use_member]->xfxQ(0,x2,muF))/(x1*x2);
}
double gluon_d1PDF(double x1, double x2)
{
	//one div
	return single_div(0, x1)*single_div(0, x2);
}
double gluon_d2PDF(double x1, double x2)
{
	//double div
	return (single_div(0, x1)+x1*double_div(0,x1))*(single_div(0, x2)+x2*double_div(0,x2));

}
double quark_d0PDF(double x1, double x2)
{
	double factor = 0;
	for(int i =1; i<6; i++)//it says so in the .info file!
			{
				// we can have qqbar or qbarq, so add those together!
				//no div
				factor = factor+((pdfs[use_member]->xfxQ(i,x1,muF))*(pdfs[use_member]->xfxQ(-i,x2,muF)) + (pdfs[use_member]->xfxQ(-i,x1,muF))*(pdfs[use_member]->xfxQ(i,x2,muF)))/(x1*x2);
			}
	return factor;

}
double quark_d1PDF(double x1, double x2)
{
	double factor = 0;
	for(int i =1; i<6; i++)
			{
				//one div
				factor = (single_div(i, x1)*single_div(-i, x2)+single_div(i, x2)*single_div(-i, x1)) + factor;
			}
	return factor;

}
double quark_d2PDF(double x1, double x2)
{
	double factor = 0;
	for(int i =1; i<6; i++)
			{
				//double div
				factor = (single_div(i, x1)+x1*double_div(i,x1))*(single_div(-i, x2)+x2*double_div(-i,x2)) +  (single_div(-i, x1)+x1*double_div(-i,x1))*(single_div(i, x2)+x2*double_div(i,x2)) + factor;
			}
	return factor;

}


double lumni(string channel, double t, double y){ //checked that this gives same after integration of pdf_sum_qqbar_charge_weighted
	double x = exp(t*log(y));
	double sum_pdf(0);
	if(channel == "qqbar"){
			double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
			if(x < y){return 0;}
			if(x >= 1){return 0;}
			for(int i = 1; i <=5; i++){
				sum_pdf+= eq[i-1]*eq[i-1]*1./y*(pdfs[use_member]->xfxQ(i,x,muF)*pdfs[use_member]->xfxQ(-i,y/x,muF)+pdfs[use_member]->xfxQ(i,y/x,muF)*pdfs[use_member]->xfxQ(-i,x,muF));
		//cout << "in the sum " << sum_pdf << endl;
		}
		return -log(y)*sum_pdf; //log(y) is jacobian left after dx/x -> dt, -1* to integrate from 0 to 1
	}
	else if(channel=="gg"){
		if(x < y){return 0;}
		if(x >= 1){return 0;}
		sum_pdf = 1./y*(pdfs[use_member]->xfxQ(0,x,muF)*pdfs[use_member]->xfxQ(0,y/x,muF));
		return -log(y)*sum_pdf; //log(y) is jacobian left after dx/x -> dt, -1* to integrate from 0 to 1
	}
}

vector<double> deriv_to_y_luminosities(string channel, double x, double y){ // returns vector with (lumni, d lumni / dy, d^2 lumni / dy^2 )
   double delta = y*1.E-3;
	 double lumi0 = lumni(channel,x,y);
   double lumidp = lumni(channel,x,y+delta);
   double lumidm = lumni(channel,x,y-delta);
   double lumid2p = lumni(channel,x,y+2.0*delta);
   double lumid2m = lumni(channel,x,y-2.0*delta);
   vector<double> results = {lumi0, (-lumid2p+8.*lumidp-8*lumidm+lumid2m)/(12.*delta), (-lumid2p+16.*lumidp-30.*lumi0+16.*lumidm-lumid2m)/(12.*delta*delta)};
   return results;
}
