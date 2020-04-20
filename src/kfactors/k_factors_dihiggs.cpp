#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "clooptools.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "k_factors_higgs.h"
#include "k_factors_dihiggs.h"
using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
/// Contains the LO dihiggs for SM and SUSY
/////////////////////////////////////////////////////////////////////////////////////////

/// all the triangle (C) and box (D) integrals
complex<double> Cab(double s, double mQ2){
	return C0(0,0,s,mQ2,mQ2,mQ2);
}
complex<double> Cac(double t, double mHc2, double mQ2){
	return C0(0,mHc2,t,mQ2,mQ2,mQ2);
}
complex<double> Cad(double u, double mHd2, double mQ2){
	return C0(0,mHd2,u,mQ2,mQ2,mQ2);
}
complex<double> Cbc(double u, double mHc2, double mQ2){
	return C0(0,mHc2,u,mQ2,mQ2,mQ2);
}
complex<double> Cbd(double t, double mHd2, double mQ2){
	return C0(0,mHd2,t,mQ2,mQ2,mQ2);
}
complex<double> Ccd(double s, double mHc2, double mHd2, double mQ2){
	return C0(mHc2,mHd2,s,mQ2,mQ2,mQ2);
}
complex<double> Dabc(double s, double u, double mHc2, double mHd2, double mQ2){
	return D0(0,0,mHc2,mHd2,s,u,mQ2,mQ2,mQ2,mQ2);
}
complex<double> Dacb(double t, double u, double mHc2, double mHd2, double mQ2){
	return D0(0,mHc2,0,mHd2,t,u,mQ2,mQ2,mQ2,mQ2);
}
complex<double> Dbac(double s, double t, double mHc2, double mHd2, double mQ2){
	return D0(0,0,mHc2,mHd2,s,t,mQ2,mQ2,mQ2,mQ2);
}


////////////////////////////////
/// form factors
////////////////////////////////

/// approximated
// https://arxiv.org/pdf/1603.00385.pdf formula 19 with TF = 1./2. includedx
double Ftriangle_approx(double s){
	return 1./2.*(4./3.+7./90.*s/mt2+1./126.*pow(s,2)/pow(mt2,2)+13./12600.*pow(s,3)/pow(mt2,3)+8./51975.*pow(s,4)/pow(mt2,4));
}
// https://arxiv.org/pdf/1603.00385.pdf formula 20 with TF = 1./2. included, but see fig3: poor approximation
// pT2 = (t*u-pow(mH2,2))/s
// also you see an extra factor of 1./2. there in the full cross section (11), which is definately due to identical particles
double Fbox_approx(double s,double pT2){
	return 1./2.*(-4./3.-7./14.*mH2/mt2-(45.*pow(mH2,2)-14.*mH2*s+6.*pow(s,2))/(315.*pow(mt2,2))+(13./630*pT2*s/pow(mt2,2))
		-(780.*pow(mH2,3)-620.*pow(mH2,2)*s+355.*mH2*pow(s,2)-16.*pow(s,3))/(18900.*pow(mt2,3))-pT2*(11.*pow(s,2)-36.*mH2*s)/(1890.*pow(mt2,3))
		-(2400.*pow(mH2,4)-3480.*pow(mH2,3)*s+2955.*pow(mH2,2)*pow(s,2)-704.*mH2*pow(s,3)+120.*pow(s,4))/(207000.*pow(mt2,4))
		+ pT2*s*(114.*pow(mH2,2)-85.*mH2*s+16.*pow(s,2)-8.*pT2*s)/(10395.*pow(mt2,4)));

}
// https://arxiv.org/pdf/1603.00385.pdf formula 21 with TF = 1./2. included
double Gbox_approx(double s,double pT2){
	return 1./2.*pT2/mt2*(-11./45.-(62.*mH2-5.*s)/(630.*mt2)-(400.*pow(mH2,2)-156.*mH2*s+49.*pow(s,2))/(12600.*pow(mt2,2))+103./18900.*pT2*s/pow(mt2,2)-(980.*pow(mH2,3)-867.*pow(mH2,2)*s+469.*mH2*pow(s,2)-34.*pow(s,3))/(103950.*pow(mt2,3))+pT2*s*(24.*mH2-7.*s)/(4950.*pow(mt2,3)));

}


/// scalar - scalar
complex<double> Ftriangle_scalar_scalar(double s, double mQ2, complex<double> cab)
{
	return 2./(s/mQ2)*(2.+(4.-s/mQ2)*mQ2*cab);
}
complex<double> Fbox_scalar_scalar(double s, double t, double u, double mQ2, double mHc2, double mHd2, complex<double> cab, complex<double> cac , complex<double> cbc, complex<double> cad, complex<double> cbd, complex<double> dabc, complex<double> dbac, complex<double> dacb)
{
	return 1./(pow(s/mQ2,2))*(4.*s/mQ2+8.*s*cab
														-2.*s*(s+mHc2+mHd2-8.*mQ2)*(dabc+dbac+dacb)
														+(mHc2+mHd2-8.*mQ2)*((t-mHc2)/mQ2*cac+(u-mHc2)/mQ2*cbc
														+(u-mHd2)/mQ2*cad+(t-mHd2)/mQ2*cbd
														-(t*u-mHc2*mHd2)/mQ2*dacb));
}
complex<double> Gbox_scalar_scalar(double s, double t, double u, double mQ2, double mHc2, double mHd2, complex<double> cab, complex<double> cac , complex<double> cbc, complex<double> cad, complex<double> cbd, complex<double> ccd, complex<double> dabc, complex<double> dbac, complex<double> dacb)
{
	return 1./(s/mQ2*(t*u-mHc2*mHd2)/mQ2/mQ2)*((t*t+mHc2*mHd2-8.*t*mQ2)/mQ2*(s/mQ2*cab+(t-mHc2)/mQ2*cac+(t-mHd2)/mQ2*cbd-s*t/mQ2*dbac)
																						+(u*u+mHc2*mHd2-8.*u*mQ2)/mQ2*(s/mQ2*cab+(u-mHc2)/mQ2*cbc+(u-mHd2)/mQ2*cad-s*u/mQ2*dabc)
																						-(t*t+u*u-2.*mHc2*mHd2)*(t+u-8.*mQ2)/mQ2/mQ2*ccd
																						-2.*(t+u-8.*mQ2)/mQ2*(t*u-mHd2*mHc2)*(dabc+dbac+dacb));
}

//pseudoscalar (C) - scalar (D)
complex<double> FtriangleA_pseudo_scalar(double mQ2, complex<double> cab)
{
	return -2.*mQ2*cab;
}
complex<double> FtriangleZ_pseudo_scalar(double s, double mQ2, complex<double> cab, double mAc2, double mHd2)
{
	return (1.-s/mZ2)*(mAc2-mHd2)/(s)*(1.+2.*mQ2*cab);
}
complex<double> Fbox_pseudo_scalar(double s, double t, double u, double mQ2, double mAc2, double mHd2, complex<double> cab, complex<double> cac, complex<double> cbc, complex<double> cad, complex<double> cbd,  complex<double> dabc, complex<double> dbac, complex<double> dacb)
{
	return 1./(pow(s/mQ2,2))*(-2.*s*(s+mAc2-mHd2)*(dabc+dbac+dacb)
														+(mAc2-mHd2)*((t-mAc2)/mQ2*cac+(u-mAc2)/mQ2*cbc+(u-mHd2)/mQ2*cad+(t-mHd2)/mQ2*cbd-(t*u-mHd2*mAc2)/mQ2*dacb));
}
complex<double> Gbox_pseudo_scalar(double s, double t, double u, double mQ2, double mAc2, double mHd2, complex<double> cab, complex<double> cac, complex<double> cbc, complex<double> cad, complex<double> cbd, complex<double> ccd,  complex<double> dabc, complex<double> dbac, complex<double> dacb)
{
	return 1./(s/mQ2*(t*u-mAc2*mHd2)/mQ2/mQ2)*((pow(u,2)-mAc2*mHd2)/mQ2*(s/mQ2*cab+(u-mAc2)/mQ2*cbc+(u-mHd2)/mQ2*cad-s*u/mQ2*dabc)
																							-(t*t-mAc2*mHd2)/mQ2*(s/mQ2*cab+(t-mAc2)/mQ2*cac+(t-mHd2)/mQ2*cbd-s*t/mQ2*dbac)
																							+(pow(t+u,2)-4.*mAc2*mHd2)/mQ2*(t-u)/mQ2*ccd
																							+2.*(t-u)*(t*u-mAc2*mHd2)/mQ2*(dabc+dbac+dacb));
}

// pseudo-pseudo
complex<double> Ftriangle_pseudo_pseudo(double s, double mQ2, complex<double> cab)
{
	return 2./(s/mQ2)*(2.+(4.-s/mQ2)*mQ2*cab);
}
complex<double> Fbox_pseudo_pseudo(double s, double t, double u, double mQ2, double mAc2, double mAd2, complex<double> cab, complex<double> cac , complex<double> cbc, complex<double> cad, complex<double> cbd, complex<double> dabc, complex<double> dbac, complex<double> dacb)
{
	return 1./(pow(s/mQ2,2))*(4.*s/mQ2+8.*s*cab-2.*s*(t+u)*(dabc+dbac+dacb)
														+(mAc2+mAd2)*((t-mAc2)/mQ2*cac+(u-mAc2)/mQ2*cbc+(u-mAd2)/mQ2*cad+(t-mAd2)/mQ2*cbd-(t*u-mAc2*mAd2)/mQ2*dacb));
}
complex<double> Gbox_pseudo_pseudo(double s, double t, double u, double mQ2, double mAc2, double mAd2, complex<double> cab, complex<double> cac , complex<double> cbc, complex<double> cad, complex<double> cbd, complex<double> ccd, complex<double> dabc, complex<double> dbac, complex<double> dacb)
{
	return 1./(s/mQ2*(t*u-mAc2*mAd2)/mQ2/mQ2)*((t*t+mAc2*mAd2)/mQ2*(s/mQ2*cab+(t-mAc2)/mQ2*cac+(t-mAd2)/mQ2*cbd-s*t/mQ2*dbac)
																						+(u*u+mAc2*mAd2)/mQ2*(s/mQ2*cab+(u-mAc2)/mQ2*cbc+(u-mAd2)/mQ2*cad-s*u/mQ2*dabc)
																						-(t*t+u*u-2.*mAc2*mAd2)*(t+u)/mQ2/mQ2*ccd
																						-2.*(t+u)/mQ2*(t*u-mAd2*mAc2)*(dabc+dbac+dacb)
																				);
}

///////////////////////////////////
// dsigma/dt
///////////////////////////////////

//SM
double dihiggs_LO_factor_SM(double scale2, double ctheta){
	double Cbox = 1;
	double Cdelta = 3.*mH2/(scale2-mH2);

	double s = scale2;
	double t = (1./2.)*ctheta*sqrt(scale2*scale2 - 4*mH2*scale2) + mH2 - scale2/2.;
	double u = 2.*mH2-t-s;
	//mQ2 = pow(100000.,2);
	//mH2-1./2.*(scale2+sqrt(Q2*Q2+mH2*mH2-2.*Q2*mH2)*ctheta);
	complex<double> cab = Cab(s,mt2);
	complex<double> cac = Cac(t,mH2,mt2);
	complex<double> cad = Cad(u,mH2,mt2);
	complex<double> cbc = cad; //Cbc(u,mH2,mQ2);
	complex<double> cbd = cac; //Cbd(t,mH2,mQ2);
	complex<double> ccd = Ccd(s,mH2,mH2,mt2);
	complex<double> dabc = Dabc(s,u,mH2,mH2,mt2);
	complex<double> dacb = Dacb(t,u,mH2,mH2,mt2);
	complex<double> dbac = Dbac(s,t,mH2,mH2,mt2);
	complex<double> fdel = Ftriangle_scalar_scalar(s, mt2, cab);
	complex<double> fbox = Fbox_scalar_scalar(s, t, u, mt2, mH2, mH2, cab, cac, cbc, cad, cbd, dabc, dbac, dacb);
	complex<double> gbox = Gbox_scalar_scalar(s, t, u, mt2, mH2, mH2, cab, cac, cbc, cad, cbd, ccd, dabc, dbac, dacb);

	// extra 0.5 for idential final states!
	double result = 0.5*fbunits*0.5*sqrt(((pow(scale2-mH2-mH2,2)-4.*mH2*mH2)))*pow(GF,2)*pow(alphas_muR,2)/(256.*pow(2.*M_PI,3))*(norm(Cdelta*fdel+Cbox*fbox)+norm(Cbox*gbox));
	clearcache();
	return result;
}

//approximate SM factor
double dihiggs_LO_factor_approx(double scale2, double ctheta){
  //ctheta from -1 to 1, scale2 from (2mH)^2 to S
	double Cbox = 1;
	double Cdelta = 3.*mH2/(scale2-mH2);
	complex<double> fdel = 2./3.;
	complex<double> fbox = -2./3.;
	complex<double> gbox = 0.;
	// approximation of Fbox and Gbox, implementation does not work yet
	double result = fbunits*0.5*0.5*scale2*sqrt(1.-4.*mH2/scale2)*pow(GF,2)*pow(alphas_muR,2)/(256.*pow(2.*M_PI,3))*(norm(Cdelta*fdel+Cbox*fbox)+norm(Cbox*gbox));
	clearcache();
	return result;
}

// SUSY cases
double  dihiggs_hh(double scale2, double ctheta){

	double lambdah = lambdahhh, lambdaH = lambdaHhh, Cbox_t = ght*ght, Cbox_b = ghb*ghb;
	complex<double> Cdeltah_t = lambdah*mZ2/(scale2-mH2+I*sqrt(mH2)*GammaH)*ght;
	complex<double> CdeltaH_t = lambdaH*mZ2/(scale2-mHeavy2+I*sqrt(mHeavy2)*GammaHeavy)*gHt;
	complex<double> Cdelta_t = Cdeltah_t+CdeltaH_t;
	complex<double> Cdeltah_b = lambdah*mZ2/(scale2-mH2+I*sqrt(mH2)*GammaH)*ghb;
	complex<double> CdeltaH_b = lambdaH*mZ2/(scale2-mHeavy2+I*sqrt(mHeavy2)*GammaHeavy)*gHb;
	complex<double> Cdelta_b = Cdeltah_b+CdeltaH_b;
	double mC2 = mH2;
	double mD2 = mH2;
	double s = scale2;
	double lamba2 = sqrt(scale2*scale2+mC2*mC2+mD2*mD2-2.*scale2*(mC2+mD2)-2.*mC2*mD2);
	double jac = 1./2.*lamba2;
	double t = (1./2.)*(-ctheta*lamba2 + mC2 + mD2 - scale2);
	double u = mD2+mC2-t-s;
	complex<double> cab_t = Cab(s,mt2);
	complex<double> cac_t = Cac(t,mC2,mt2);
	complex<double> cad_t = Cad(u,mD2,mt2);
	complex<double> cbc_t = Cbc(u,mC2,mt2);
	complex<double> cbd_t = Cbd(t,mD2,mt2);
	complex<double> ccd_t = Ccd(s,mC2,mD2,mt2);
	complex<double> dabc_t = Dabc(s,u,mC2,mD2,mt2);
	complex<double> dacb_t = Dacb(t,u,mC2,mD2,mt2);
	complex<double> dbac_t = Dbac(s,t,mC2,mD2,mt2);
	complex<double> ftriangle_ss_t = Ftriangle_scalar_scalar(s, mt2,cab_t);
	complex<double> fbox_ss_t = Fbox_scalar_scalar(s, t, u, mt2, mC2, mD2, cab_t, cac_t, cbc_t, cad_t, cbd_t, dabc_t, dbac_t, dacb_t);
	complex<double> gbox_ss_t = Gbox_scalar_scalar(s, t, u, mt2, mC2, mD2, cab_t, cac_t, cbc_t, cad_t, cbd_t, ccd_t, dabc_t, dbac_t, dacb_t);
	clearcache();
	complex<double> cab_b = Cab(s,mb2);
	complex<double> cac_b = Cac(t,mC2,mb2);
	complex<double> cad_b = Cad(u,mD2,mb2);
	complex<double> cbc_b = Cbc(u,mC2,mb2);
	complex<double> cbd_b = Cbd(t,mD2,mb2);
	complex<double> ccd_b = Ccd(s,mC2,mD2,mb2);
	complex<double> dabc_b = Dabc(s,u,mC2,mD2,mb2);
	complex<double> dacb_b = Dacb(t,u,mC2,mD2,mb2);
	complex<double> dbac_b = Dbac(s,t,mC2,mD2,mb2);
	clearcache();
	complex<double> ftriangle_ss_b = Ftriangle_scalar_scalar(s, mb2, cab_b);
	complex<double> fbox_ss_b = Fbox_scalar_scalar(s, t, u, mb2, mC2, mD2, cab_b, cac_b, cbc_b, cad_b, cbd_b, dabc_b, dbac_b, dacb_b);
	complex<double> gbox_ss_b = Gbox_scalar_scalar(s, t, u, mb2, mC2, mD2, cab_b, cac_b, cbc_b, cad_b, cbd_b, ccd_b, dabc_b, dbac_b, dacb_b);
	double scale = 0.5;
	double result = scale*fbunits*0.5*jac*pow(GF,2)*pow(alphas_muR,2)/(256.*pow(2.*M_PI,3))*(norm(Cdelta_t*ftriangle_ss_t+Cbox_t*fbox_ss_t+Cdelta_b*ftriangle_ss_b+Cbox_b*fbox_ss_b)+norm(Cbox_t*gbox_ss_t+Cbox_b*gbox_ss_b));
	return result;
}


double  dihiggs_hH(double scale2, double ctheta){
	double lambdah = lambdaHhh, lambdaH = lambdaHHh, Cbox_t = ght*gHt, Cbox_b = ghb*gHb;
	double Cdeltah_t = lambdah*mZ2/(scale2-mH2)*ght;
	double CdeltaH_t = lambdaH*mZ2/(scale2-mHeavy2)*gHt;
	double Cdelta_t = Cdeltah_t+CdeltaH_t;
	double Cdeltah_b = lambdah*mZ2/(scale2-mH2)*ghb;
	double CdeltaH_b = lambdaH*mZ2/(scale2-mHeavy2)*gHb;
	double Cdelta_b = Cdeltah_b+CdeltaH_b;
	double mC2 = mH2;
	double mD2 = mHeavy2;
	double s = scale2;
	double lamba2 = sqrt(scale2*scale2+mC2*mC2+mD2*mD2-2.*scale2*(mC2+mD2)-2.*mC2*mD2);
	double jac = 1./2.*lamba2;
	double t = (1./2.)*(-ctheta*lamba2 + mC2 + mD2 - scale2);
	double u = mD2+mC2-t-s;
	complex<double> cab_t = Cab(s,mt2);
	complex<double> cac_t = Cac(t,mC2,mt2);
	complex<double> cad_t = Cad(u,mD2,mt2);
	complex<double> cbc_t = Cbc(u,mC2,mt2);
	complex<double> cbd_t = Cbd(t,mD2,mt2);
	complex<double> ccd_t = Ccd(s,mC2,mD2,mt2);
	complex<double> dabc_t = Dabc(s,u,mC2,mD2,mt2);
	complex<double> dacb_t = Dacb(t,u,mC2,mD2,mt2);
	complex<double> dbac_t = Dbac(s,t,mC2,mD2,mt2);
	complex<double> ftriangle_ss_t = Ftriangle_scalar_scalar(s,mt2, cab_t);
	complex<double> fbox_ss_t = Fbox_scalar_scalar(s, t, u, mt2, mC2, mD2, cab_t, cac_t, cbc_t, cad_t, cbd_t, dabc_t, dbac_t, dacb_t);
	complex<double> gbox_ss_t = Gbox_scalar_scalar(s, t, u, mt2, mC2, mD2, cab_t, cac_t, cbc_t, cad_t, cbd_t, ccd_t, dabc_t, dbac_t, dacb_t);
	clearcache();
	complex<double> cab_b = Cab(s,mb2);
	complex<double> cac_b = Cac(t,mC2,mb2);
	complex<double> cad_b = Cad(u,mD2,mb2);
	complex<double> cbc_b = Cbc(u,mC2,mb2);
	complex<double> cbd_b = Cbd(t,mD2,mb2);
	complex<double> ccd_b = Ccd(s,mC2,mD2,mb2);
	complex<double> dabc_b = Dabc(s,u,mC2,mD2,mb2);
	complex<double> dacb_b = Dacb(t,u,mC2,mD2,mb2);
	complex<double> dbac_b = Dbac(s,t,mC2,mD2,mb2);
	clearcache();
	complex<double> ftriangle_ss_b = Ftriangle_scalar_scalar(s,mb2, cab_b);
	complex<double> fbox_ss_b = Fbox_scalar_scalar(s, t, u, mb2, mC2, mD2, cab_b, cac_b, cbc_b, cad_b, cbd_b, dabc_b, dbac_b, dacb_b);
	complex<double> gbox_ss_b = Gbox_scalar_scalar(s, t, u, mb2, mC2, mD2, cab_b, cac_b, cbc_b, cad_b, cbd_b, ccd_b, dabc_b, dbac_b, dacb_b);
	double scale = 1;
	double result = scale*fbunits*0.5*jac*pow(GF,2)*pow(alphas_muR,2)/(256.*pow(2.*M_PI,3))*(norm(Cdelta_t*ftriangle_ss_t+Cbox_t*fbox_ss_t+Cdelta_b*ftriangle_ss_b+Cbox_b*fbox_ss_b)+norm(Cbox_t*gbox_ss_t+Cbox_b*gbox_ss_b));
	return result;
}


double  dihiggs_HH(double scale2, double ctheta){
	double lambdah = lambdaHHh, lambdaH = lambdaHHH, Cbox_t = gHt*gHt, Cbox_b = gHb*gHb;
	double Cdeltah_t = lambdah*mZ2/(scale2-mH2)*ght;
	double CdeltaH_t = lambdaH*mZ2/(scale2-mHeavy2)*gHt;
	double Cdelta_t = Cdeltah_t+CdeltaH_t;
	double Cdeltah_b = lambdah*mZ2/(scale2-mH2)*ghb;
	double CdeltaH_b = lambdaH*mZ2/(scale2-mHeavy2)*gHb;
	double Cdelta_b = Cdeltah_b+CdeltaH_b;
	double mC2 = mHeavy2;
	double mD2 = mHeavy2;
	double s = scale2;
	double lamba2 = sqrt(scale2*scale2+mC2*mC2+mD2*mD2-2.*scale2*(mC2+mD2)-2.*mC2*mD2);
	double jac = 1./2.*lamba2;
	double t = (1./2.)*(-ctheta*lamba2 + mC2 + mD2 - scale2);
	double u = mD2+mC2-t-s;
	complex<double> cab_t = Cab(s,mt2);
	complex<double> cac_t = Cac(t,mC2,mt2);
	complex<double> cad_t = Cad(u,mD2,mt2);
	complex<double> cbc_t = Cbc(u,mC2,mt2);
	complex<double> cbd_t = Cbd(t,mD2,mt2);
	complex<double> ccd_t = Ccd(s,mC2,mD2,mt2);
	complex<double> dabc_t = Dabc(s,u,mC2,mD2,mt2);
	complex<double> dacb_t = Dacb(t,u,mC2,mD2,mt2);
	complex<double> dbac_t = Dbac(s,t,mC2,mD2,mt2);
	complex<double> ftriangle_ss_t = Ftriangle_scalar_scalar(s,mt2, cab_t);
	complex<double> fbox_ss_t = Fbox_scalar_scalar(s, t, u, mt2, mC2, mD2, cab_t, cac_t, cbc_t, cad_t, cbd_t, dabc_t, dbac_t, dacb_t);
	complex<double> gbox_ss_t = Gbox_scalar_scalar(s, t, u, mt2, mC2, mD2, cab_t, cac_t, cbc_t, cad_t, cbd_t, ccd_t, dabc_t, dbac_t, dacb_t);
	clearcache();
	complex<double> cab_b = Cab(s,mb2);
	complex<double> cac_b = Cac(t,mC2,mb2);
	complex<double> cad_b = Cad(u,mD2,mb2);
	complex<double> cbc_b = Cbc(u,mC2,mb2);
	complex<double> cbd_b = Cbd(t,mD2,mb2);
	complex<double> ccd_b = Ccd(s,mC2,mD2,mb2);
	complex<double> dabc_b = Dabc(s,u,mC2,mD2,mb2);
	complex<double> dacb_b = Dacb(t,u,mC2,mD2,mb2);
	complex<double> dbac_b = Dbac(s,t,mC2,mD2,mb2);
	clearcache();
	complex<double> ftriangle_ss_b = Ftriangle_scalar_scalar(s,mb2, cab_b);
	complex<double> fbox_ss_b = Fbox_scalar_scalar(s, t, u, mb2, mC2, mD2, cab_b, cac_b, cbc_b, cad_b, cbd_b, dabc_b, dbac_b, dacb_b);
	complex<double> gbox_ss_b = Gbox_scalar_scalar(s, t, u, mb2, mC2, mD2, cab_b, cac_b, cbc_b, cad_b, cbd_b, ccd_b, dabc_b, dbac_b, dacb_b);
	double result = 0.5*fbunits*0.5*jac*pow(GF,2)*pow(alphas_muR,2)/(256.*pow(2.*M_PI,3))*(norm(Cdelta_t*ftriangle_ss_t+Cbox_t*fbox_ss_t+Cdelta_b*ftriangle_ss_b+Cbox_b*fbox_ss_b)+norm(Cbox_t*gbox_ss_t+Cbox_b*gbox_ss_b));
	return result;
}


//mixed
double dihiggs_Ah(double scale2, double ctheta){

	double lambdaA = lambdahAA, lambdaZ = lambdaZAh, Cbox_t = gAt*ght, Cbox_b = gAb*ghb;
	double mC2 = mA2;
	double mD2 = mH2;
	double CdeltaA_t = lambdaA*mZ2/(scale2-mC2)*gAt;
	double CdeltaZ_t = lambdaZ*mZ2/(scale2-mZ2);
	double CdeltaA_b = lambdaA*mZ2/(scale2-mC2)*gAb;
	double CdeltaZ_b = lambdaZ*mZ2/(scale2-mZ2)*-1.;

	double s = scale2;
	double lamba2 = sqrt(scale2*scale2+mC2*mC2+mD2*mD2-2.*scale2*(mC2+mD2)-2.*mC2*mD2);
	double jac = 1./2.*lamba2;
	double t = (1./2.)*(-ctheta*lamba2 + mC2 + mD2 - scale2);
	double u = mD2+mC2-t-s;
	complex<double> cab_t = Cab(s,mt2);
	complex<double> cac_t = Cac(t,mC2,mt2);
	complex<double> cad_t = Cad(u,mD2,mt2);
	complex<double> cbc_t = Cbc(u,mC2,mt2);
	complex<double> cbd_t = Cbd(t,mD2,mt2);
	complex<double> ccd_t = Ccd(s,mC2,mD2,mt2);
	complex<double> dabc_t = Dabc(s,u,mC2,mD2,mt2);
	complex<double> dacb_t = Dacb(t,u,mC2,mD2,mt2);
	complex<double> dbac_t = Dbac(s,t,mC2,mD2,mt2);
	clearcache();
	complex<double> ftriangleA_ps_t = FtriangleA_pseudo_scalar(mt2, cab_t);
	complex<double> ftriangleZ_ps_t = FtriangleZ_pseudo_scalar(s, mt2, cab_t, mC2, mD2);
	complex<double> fbox_ps_t = Fbox_pseudo_scalar(s, t, u, mt2, mC2, mD2, cab_t, cac_t, cbc_t, cad_t, cbd_t, dabc_t, dbac_t, dacb_t);
	complex<double> gbox_ps_t = Gbox_pseudo_scalar(s, t, u, mt2, mC2, mD2, cab_t, cac_t, cbc_t, cad_t, cbd_t, ccd_t, dabc_t, dbac_t, dacb_t);

	complex<double> cab_b = Cab(s,mb2);
	complex<double> cac_b = Cac(t,mC2,mb2);
	complex<double> cad_b = Cad(u,mD2,mb2);
	complex<double> cbc_b = Cbc(u,mC2,mb2);
	complex<double> cbd_b = Cbd(t,mD2,mb2);
	complex<double> ccd_b = Ccd(s,mC2,mD2,mb2);
	complex<double> dabc_b = Dabc(s,u,mC2,mD2,mb2);
	complex<double> dacb_b = Dacb(t,u,mC2,mD2,mb2);
	complex<double> dbac_b = Dbac(s,t,mC2,mD2,mb2);
	clearcache();
	complex<double> ftriangleA_ps_b = FtriangleA_pseudo_scalar(mb2,cab_b);
	complex<double> ftriangleZ_ps_b = FtriangleZ_pseudo_scalar(s,mb2,cab_b, mC2, mD2);
	complex<double> fbox_ps_b = Fbox_pseudo_scalar(s, t, u,mb2, mC2, mD2, cab_b, cac_b, cbc_b, cad_b, cbd_b, dabc_b, dbac_b, dacb_b);
	complex<double> gbox_ps_b = Gbox_pseudo_scalar(s, t, u,mb2, mC2, mD2, cab_b, cac_b, cbc_b, cad_b, cbd_b, ccd_b, dabc_b, dbac_b, dacb_b);

	double result = fbunits*0.5*jac*pow(GF,2)*pow(alphas_muR,2)/(256.*pow(2.*M_PI,3))*(norm(CdeltaA_t*ftriangleA_ps_t+CdeltaZ_t*ftriangleZ_ps_t+CdeltaA_b*ftriangleA_ps_b+CdeltaZ_b*ftriangleZ_ps_b+Cbox_t*fbox_ps_t+Cbox_b*fbox_ps_b)+norm(Cbox_t*gbox_ps_t+Cbox_b*gbox_ps_b));

	return result;
}


double dihiggs_AH(double scale2, double ctheta){
	double lambdaA = lambdaHAA, lambdaZ = lambdaZAH, Cbox_t = gAt*gHt, Cbox_b = gAb*gHb;
	double mC2 = mA2;
	double mD2 = mHeavy2;
	double CdeltaA_t = lambdaA*mZ2/(scale2-mC2)*gAt;
	double CdeltaZ_t = lambdaZ*mZ2/(scale2-mZ2);
	double CdeltaA_b = lambdaA*mZ2/(scale2-mC2)*gAb;
	double CdeltaZ_b = lambdaZ*mZ2/(scale2-mZ2)*-1.;

	double s = scale2;
	double lamba2 = sqrt(scale2*scale2+mC2*mC2+mD2*mD2-2.*scale2*(mC2+mD2)-2.*mC2*mD2);
	double jac = 1./2.*lamba2;
	double t = (1./2.)*(-ctheta*lamba2 + mC2 + mD2 - scale2);
	double u = mD2+mC2-t-s;
	complex<double> cab_t = Cab(s,mt2);
	complex<double> cac_t = Cac(t,mC2,mt2);
	complex<double> cad_t = Cad(u,mD2,mt2);
	complex<double> cbc_t = Cbc(u,mC2,mt2);
	complex<double> cbd_t = Cbd(t,mD2,mt2);
	complex<double> ccd_t = Ccd(s,mC2,mD2,mt2);
	complex<double> dabc_t = Dabc(s,u,mC2,mD2,mt2);
	complex<double> dacb_t = Dacb(t,u,mC2,mD2,mt2);
	complex<double> dbac_t = Dbac(s,t,mC2,mD2,mt2);
	clearcache();
	complex<double> ftriangleA_ps_t = FtriangleA_pseudo_scalar(mt2,cab_t);
	complex<double> ftriangleZ_ps_t = FtriangleZ_pseudo_scalar(s,mt2,cab_t, mC2, mD2);
	complex<double> fbox_ps_t = Fbox_pseudo_scalar(s, t, u, mt2,mC2, mD2, cab_t, cac_t, cbc_t, cad_t, cbd_t, dabc_t, dbac_t, dacb_t);
	complex<double> gbox_ps_t = Gbox_pseudo_scalar(s, t, u, mt2,mC2, mD2, cab_t, cac_t, cbc_t, cad_t, cbd_t, ccd_t, dabc_t, dbac_t, dacb_t);

	complex<double> cab_b = Cab(s,mb2);
	complex<double> cac_b = Cac(t,mC2,mb2);
	complex<double> cad_b = Cad(u,mD2,mb2);
	complex<double> cbc_b = Cbc(u,mC2,mb2);
	complex<double> cbd_b = Cbd(t,mD2,mb2);
	complex<double> ccd_b = Ccd(s,mC2,mD2,mb2);
	complex<double> dabc_b = Dabc(s,u,mC2,mD2,mb2);
	complex<double> dacb_b = Dacb(t,u,mC2,mD2,mb2);
	complex<double> dbac_b = Dbac(s,t,mC2,mD2,mb2);
	clearcache();
	complex<double> ftriangleA_ps_b = FtriangleA_pseudo_scalar(mb2,cab_b);
	complex<double> ftriangleZ_ps_b = FtriangleZ_pseudo_scalar(s, mb2, cab_b, mC2, mD2);
	complex<double> fbox_ps_b = Fbox_pseudo_scalar(s, t, u, mb2, mC2, mD2, cab_b, cac_b, cbc_b, cad_b, cbd_b, dabc_b, dbac_b, dacb_b);
	complex<double> gbox_ps_b = Gbox_pseudo_scalar(s, t, u, mb2, mC2, mD2, cab_b, cac_b, cbc_b, cad_b, cbd_b, ccd_b, dabc_b, dbac_b, dacb_b);

	double result = fbunits*0.5*jac*pow(GF,2)*pow(alphas_muR,2)/(256.*pow(2.*M_PI,3))*(norm(CdeltaA_t*ftriangleA_ps_t+CdeltaZ_t*ftriangleZ_ps_t+CdeltaA_b*ftriangleA_ps_b+CdeltaZ_b*ftriangleZ_ps_b+Cbox_t*fbox_ps_t+Cbox_b*fbox_ps_b)+norm(Cbox_t*gbox_ps_t+Cbox_b*gbox_ps_b));

	return result;
}


double dihiggs_AA(double scale2, double ctheta){

	double mC2 = mA2;
	double mD2 = mA2;
	double lambdah = lambdahAA, lambdaH = lambdaHAA, Cbox_t = gAt*gAt, Cbox_b = gAb*gAb ;
	double Cdeltah_t = lambdah*mZ2/(scale2-mH2)*ght;
	double CdeltaH_t = lambdaH*mZ2/(scale2-mHeavy2)*gHt;
	double Cdelta_t = Cdeltah_t+CdeltaH_t;
	double Cdeltah_b = lambdah*mZ2/(scale2-mH2)*ghb;
	double CdeltaH_b = lambdaH*mZ2/(scale2-mHeavy2)*gHb;
	double Cdelta_b = Cdeltah_b+CdeltaH_b;

	double s = scale2;
	double lamba2 = sqrt(scale2*scale2+mC2*mC2+mD2*mD2-2.*scale2*(mC2+mD2)-2.*mC2*mD2);
	double jac = 1./2.*lamba2;
	double t = (1./2.)*(-ctheta*lamba2 + mC2 + mD2 - scale2);
	double u = mD2+mC2-t-s;
	complex<double> cab_t = Cab(s,mt2);
	complex<double> cac_t = Cac(t,mC2,mt2);
	complex<double> cad_t = Cad(u,mD2,mt2);
	complex<double> cbc_t = Cbc(u,mC2,mt2);
	complex<double> cbd_t = Cbd(t,mD2,mt2);
	complex<double> ccd_t = Ccd(s,mC2,mD2,mt2);
	complex<double> dabc_t = Dabc(s,u,mC2,mD2,mt2);
	complex<double> dacb_t = Dacb(t,u,mC2,mD2,mt2);
	complex<double> dbac_t = Dbac(s,t,mC2,mD2,mt2);
	clearcache();
	complex<double> ftriangle_pp_t = Ftriangle_pseudo_pseudo(s,mt2,cab_t);
	complex<double> fbox_pp_t = Fbox_pseudo_pseudo(s, t, u,mt2, mC2, mD2, cab_t, cac_t, cbc_t, cad_t, cbd_t, dabc_t, dbac_t, dacb_t);
	complex<double> gbox_pp_t = Gbox_pseudo_pseudo(s, t, u,mt2, mC2, mD2, cab_t, cac_t, cbc_t, cad_t, cbd_t, ccd_t, dabc_t, dbac_t, dacb_t);


	complex<double> cab_b = Cab(s,mb2);
	complex<double> cac_b = Cac(t,mC2,mb2);
	complex<double> cad_b = Cad(u,mD2,mb2);
	complex<double> cbc_b = Cbc(u,mC2,mb2);
	complex<double> cbd_b = Cbd(t,mD2,mb2);
	complex<double> ccd_b = Ccd(s,mC2,mD2,mb2);
	complex<double> dabc_b = Dabc(s,u,mC2,mD2,mb2);
	complex<double> dacb_b = Dacb(t,u,mC2,mD2,mb2);
	complex<double> dbac_b = Dbac(s,t,mC2,mD2,mb2);
	clearcache();
	complex<double> ftriangle_pp_b = Ftriangle_pseudo_pseudo(s,mb2,cab_b);
	complex<double> fbox_pp_b = Fbox_pseudo_pseudo(s, t, u, mb2, mC2, mD2, cab_b, cac_b, cbc_b, cad_b, cbd_b, dabc_b, dbac_b, dacb_b);
	complex<double> gbox_pp_b = Gbox_pseudo_pseudo(s, t, u, mb2, mC2, mD2, cab_b, cac_b, cbc_b, cad_b, cbd_b, ccd_b, dabc_b, dbac_b, dacb_b);

	double result = 0.5*fbunits*0.5*jac*pow(GF,2)*pow(alphas_muR,2)/(256.*pow(2.*M_PI,3))*(norm(Cdelta_t*ftriangle_pp_t+Cbox_t*fbox_pp_t+Cdelta_b*ftriangle_pp_b+Cbox_b*fbox_pp_b)+norm(Cbox_t*gbox_pp_t+Cbox_b*gbox_pp_b));
	return result;
}
