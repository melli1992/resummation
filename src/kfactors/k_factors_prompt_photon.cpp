#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "pf_pdf.h"
#include "k_factors_prompt_photon.h"
using namespace std;


// needed every time for the transformations
// integrals should be from 0 to 1 and xT, epeta, emeta should be defined
double vtrans(double xT,double epeta,double emeta,double k0){
  return xT/2.*epeta+(1.-xT/2.*(emeta+epeta))*k0;
}
double wtrans(double xT,double epeta,double k0, double v){
  return xT/2.*epeta/v+(1.-xT/2.*(epeta)/v)*k0;
}
double jacv(double xT,double epeta,double emeta){
  return 1.-xT/2.*(epeta+emeta);
}
double jacw(double xT,double epeta, double v){
  return 1.-xT/2.*(epeta)/v;
}

//////////////////////////////////////
/// qqbar to gamma gg channel
//////////////////////////////////////
double vegas_FP_qqbartogg_LP(double *k, size_t dim, void *params){
  (void)(dim);/* here we need xT and exp[eta]=epeta and exp[-eta]=emeta*/
  struct lumni_params * pf = (struct lumni_params *)params;
  double v = vtrans(pf->xT,pf->epeta,pf->emeta,k[0]);
  double w = wtrans(pf->xT,pf->epeta,k[1],v);
	double jac = jacv(pf->xT,pf->epeta,pf->emeta)*jacv(pf->xT,pf->epeta,v);
  double x1 = (pf->xT)*(pf->epeta)/(2.*v*w);
  double x1w1 = (pf->xT)*(pf->epeta)/(2.*v);
  double x2 = (pf->xT)*(pf->emeta)/(2.*(1-v));
  return pbunits*1./(M_PI)*1./(pow(pf->pT,4))*alphaEM*pow(alphas_Q,2)*(pf_sum_qqbar_charge_weighted(x1,x2)-pf_sum_qqbar_charge_weighted(x1w1,x2))*FP_qqbartogg_LP(v,w);
}


double vegas_FP_qqbartogg_LP_corr(double *k, size_t dim, void *params){
	(void)(dim);
  struct lumni_params * pf = (struct lumni_params *)params;
  double v = vtrans(pf->xT,pf->epeta,pf->emeta,k[0]);
  double wlow = (pf->xT)/2.*(pf->epeta)/v;
  double jac = jacv(pf->xT,pf->epeta,pf->emeta)*jacv(pf->xT,pf->epeta,v);
  double x1 = (pf->xT)*(pf->epeta)/(2.*v);
  double x2 = (pf->xT)*(pf->emeta)/(2.*(1-v));
  return pbunits*1./(M_PI)*1./(pow(pf->pT,4))*alphaEM*pow(alphas_Q,2)*(pf_sum_qqbar_charge_weighted(x1,x2))*FP_qqbartogg_LP_corr(v,wlow);
}


double vegas_FP_qqbartogg_full(double *k, size_t dim, void *params){
  (void)(dim);/* here we need xT and exp[eta]=epeta and exp[-eta]=emeta*/
  struct lumni_params * pf = (struct lumni_params *)params;
  double v = vtrans(pf->xT,pf->epeta,pf->emeta,k[0]);
  double w = wtrans(pf->xT,pf->epeta,k[1],v);
	double jac = jacv(pf->xT,pf->epeta,pf->emeta)*jacv(pf->xT,pf->epeta,v);
  double x1 = (pf->xT)*(pf->epeta)/(2.*v*w);
  double x2 = (pf->xT)*(pf->emeta)/(2.*(1-v));
  return pbunits*1./(M_PI)*1./(pow(pf->pT,4))*alphaEM*pow(alphas_Q,2)*(pf_sum_qqbar_charge_weighted(x1,x2))*FP_qqbartogg_full(v,w);
}

double vegas_FP_qqbartogg_delta(double *k, size_t dim, void *params){
  (void)(dim);/* here we need xT and exp[eta]=epeta and exp[-eta]=emeta*/
  struct lumni_params * pf = (struct lumni_params *)params;
  double v = vtrans(pf->xT,pf->epeta,pf->emeta,k[0]);
  double w = 1.;
	double jac = jacv(pf->xT,pf->epeta,pf->emeta)*jacv(pf->xT,pf->epeta,v);
  double x1 = (pf->xT)*(pf->epeta)/(2.*v*w);
  double x2 = (pf->xT)*(pf->emeta)/(2.*(1-v));
  return pbunits*1./(M_PI)*1./(pow(pf->pT,4))*alphaEM*pow(alphas_Q,2)*(pf_sum_qqbar_charge_weighted(x1,x2))*FP_qqbartogg_delta(v);
}

double vegas_FP_qqbartogg_power(double *k, size_t dim, void *params)
{
  (void)(dim);/* here we need xT and exp[eta]=epeta and exp[-eta]=emeta*/
  struct lumni_params * pf = (struct lumni_params *)params;
  double v = vtrans(pf->xT,pf->epeta,pf->emeta,k[0]);
  double w = 1.;
	double jac = jacv(pf->xT,pf->epeta,pf->emeta)*jacv(pf->xT,pf->epeta,v);
  double x1 = (pf->xT)*(pf->epeta)/(2.*v*w);
  double x2 = (pf->xT)*(pf->emeta)/(2.*(1-v));
  return pbunits*1./(M_PI)*1./(pow(pf->pT,4))*alphaEM*pow(alphas_Q,2)*(pf_sum_qqbar_charge_weighted(x1,x2))*FP_qqbartogg_expansion(v,w,pf->power);
}

double FP_qqbartogg_LP(double v, double w){
  return (CF*(1 - 2*v + 2*pow(v,2))*(11*CA - 2*nF + 24*CF*log(muR2/Q2) + 24*CF*log(-1 + 1/v) - 12*CA*log(1 - v) + 12*CA*log(1 - w) - 48*CF*log(1 - w)))/(6.*CA*(-1 + w));
}

double FP_qqbartogg_LP_corr(double v, double wlow){
  return -(CF*(1 - 2*v + 2*pow(v,2))*log(1 - wlow)*(11*CA - 2*nF + 24*CF*log(muR2/Q2) + 24*CF*log(-1 + 1/v) - 12*CA*log(1 - v) + 6*CA*log(1 - wlow) - 24*CF*log(1 - wlow)))/(6.*CA);
}

double FP_qqbartogg_delta(double v){
  return (CF*(-((-67 + 6*pow(M_PI,2))*(1 - 2*v + 2*pow(v,2))) + (18*beta0/2.*(1 - 2*v + 2*pow(v,2))*log(muF2/Q2))/CA - (18*CF*(1 - 2*v + 2*pow(v,2))*log(muR2/Q2)*(-3 + 2*log(-1 + 1/v)))/CA + 18*(-1 + v)*v*log(1 - v) - 9*(3 - 4*v + 3*pow(v,2))*pow(log(1 - v),2) - 3*(11 - 16*v + 16*pow(v,2))*log(v) + 9*(-2 + v)*v*pow(log(v),2) + (2*nF*(1 - 2*v + 2*pow(v,2))*(-5 + 3*log(v)))/CA + (18*CF*((-7 + pow(M_PI,2))*(1 - 2*v + 2*pow(v,2)) + (3 - 4*v + 3*pow(v,2))*pow(log(1 - v),2) + (3 - 4*v + pow(v,2))*log(v) + (2 - 2*v + 3*pow(v,2))*pow(log(v),2) + log(1 - v)*(v*(2 + v) + (-2 + 4*v - 4*pow(v,2))*log(v))))/CA))/18.;
}
double FP_qqbartogg_full(double v, double w){
  return (CF*(-2 + 11*pow(-1 + v,2) + (3 - 12*v + 11*pow(v,2))*w - (3*(-1 + v)*v)/pow(-1 + v*w,2) - (3*v*(-1 + 4*v))/(-1 + v*w) - (2*nF*v*(-2 + v + v*w))/CA + (6*CF*(4*(-1 + v)*v + (1 + v - 4*pow(v,2))*w - ((-1 + v)*v)/pow(-1 + v*w,2) + ((-1 + v)*v)/(-1 + v*w)))/CA + (6*CF*(pow(1 - 2*v,2) + (1 - 2*v + 4*pow(v,2))*w - ((-1 + v)*v)/pow(-1 + v*w,2) + (v - 2*pow(v,2))/(-1 + v*w))*log(muR2/Q2))/CA - (6*(2*CF*(1 - 2*v + pow(v,2)*(1 + pow(w,2))) + CA*(-1 + 2*v + pow(v,3)*w*(1 + w) - pow(v,2)*(2 + w + pow(w,2))))*log(1 - v))/(CA*(-1 + v*w)) - (6*CF*(pow(1 - 2*v,2) + (1 - 2*v + 4*pow(v,2))*w - ((-1 + v)*v)/pow(-1 + v*w,2) + (v - 2*pow(v,2))/(-1 + v*w))*log(v))/CA - (6*(CA*(-1 - 2*pow(v,4)*pow(w,2)*(1 + w) + v*(4 + w) + pow(v,3)*w*(5 + 6*w + pow(w,2)) - pow(v,2)*(3 + 8*w + pow(w,2))) + CF*(3 + w + 8*pow(v,4)*pow(w,2)*(1 + w) - 2*v*(6 + 3*w + pow(w,2)) - 4*pow(v,3)*w*(5 + 6*w + pow(w,2)) + pow(v,2)*(11 + 29*w + 7*pow(w,2) + pow(w,3))))*log(1 - w))/(CA*pow(-1 + v*w,2)) - (6*(1 - 2*v + 2*pow(v,2))*log(w))/(-1 + w) - (6*(2*CF*(1 - 2*v + pow(v,2)*(1 + pow(w,2))) + CA*(-1 + 2*v + pow(v,3)*w*(1 + w) - pow(v,2)*(2 + w + pow(w,2))))*log(w))/(CA*(-1 + v*w)) + (24*CF*v*(-1 + v + v*w + CF*(-2 + v + v*w))*log(1 - v*w))/CA + (6*(CA - 4*CF)*(1 - 2*v + 2*pow(v,2))*log((-1 + v)/(-1 + v*w)))/(CA*(-1 + w))))/6.;
}
double FP_qqbartogg_expansion(double v, double w, int power){
  if(power==1){
  	return -(CF*(3*CA + 3*CF - 11*CA*v - 18*CF*v + 2*nF*v + 16*CA*pow(v,2) + 30*CF*pow(v,2) - 4*nF*pow(v,2) + CA*pow(v,3) - 24*CF*pow(v,3) + 2*nF*pow(v,3) - 6*CF*(-1 + 4*v - 8*pow(v,2) + 4*pow(v,3))*log(muR2/Q2) + 3*(CA*(-1 + 2*v - 4*pow(v,2) + 2*pow(v,3)) - 2*CF*(-1 + 4*(1 + CF)*v - 8*(1 + CF)*pow(v,2) + 4*(1 + CF)*pow(v,3)))*log(1 - v) - 6*CF*log(v) + 24*CF*v*log(v) - 48*CF*pow(v,2)*log(v) + 24*CF*pow(v,3)*log(v) + 3*CA*log(1 - w) - 12*CF*log(1 - w) - 12*CA*v*log(1 - w) + 48*CF*v*log(1 - w) + 24*CA*pow(v,2)*log(1 - w) - 96*CF*pow(v,2)*log(1 - w) - 12*CA*pow(v,3)*log(1 - w) + 48*CF*pow(v,3)*log(1 - w)))/(3.*CA*(-1 + v));
  }
  if(power==2){
  	return (CF*(-1 + w)*(18*CF - 12*CA*v - 42*CF*v + 29*CA*pow(v,2) + 48*CF*pow(v,2) + 48*pow(CF,2)*pow(v,2) - 2*nF*pow(v,2) - 10*CA*pow(v,3) - 24*CF*pow(v,3) - 96*pow(CF,2)*pow(v,3) + 4*nF*pow(v,3) + 11*CA*pow(v,4) + 48*pow(CF,2)*pow(v,4) - 2*nF*pow(v,4) + 6*CF*(1 - 4*v + 10*pow(v,2) - 8*pow(v,3) + 4*pow(v,4))*log(muR2/Q2) - 6*v*(CA*(1 + v - 2*pow(v,2) + pow(v,3)) - 2*CF*(1 + 2*(1 + CF)*v - 4*(1 + CF)*pow(v,2) + 2*(1 + CF)*pow(v,3)))*log(1 - v) - 6*CF*log(v) + 24*CF*v*log(v) - 60*CF*pow(v,2)*log(v) + 48*CF*pow(v,3)*log(v) - 24*CF*pow(v,4)*log(v) - 6*CF*log(1 - w) + 6*CA*v*log(1 - w) + 12*CF*v*log(1 - w) + 12*CA*pow(v,2)*log(1 - w) - 84*CF*pow(v,2)*log(1 - w) - 24*CA*pow(v,3)*log(1 - w) + 96*CF*pow(v,3)*log(1 - w) + 12*CA*pow(v,4)*log(1 - w) - 48*CF*pow(v,4)*log(1 - w)))/(6.*CA*pow(-1 + v,2));
  }
  if(power==3){
  	return -(CF*pow(-1 + w,2)*(CA - 6*CF - 8*CA*v + 36*CF*v + 5*CA*pow(v,2) - 54*CF*pow(v,2) + 16*CA*pow(v,3) + 40*CF*pow(v,3) - 2*CA*pow(v,4) + 10*CF*pow(v,4) + 8*CA*pow(v,5) - 16*CF*pow(v,5) + 12*CF*pow(v,3)*(1 + v)*log(muR2/Q2) - 6*(CA - 2*CF)*pow(v,2)*(2 - 2*v + pow(v,2))*log(1 - v) - 12*CF*pow(v,3)*log(v) - 12*CF*pow(v,4)*log(v) + 12*CA*pow(v,2)*log(1 - w) - 24*CF*pow(v,2)*log(1 - w) - 12*CA*pow(v,3)*log(1 - w) + 12*CF*pow(v,3)*log(1 - w) + 6*CA*pow(v,4)*log(1 - w) - 24*CF*pow(v,4)*log(1 - w)))/(6.*CA*pow(-1 + v,3));
  }
  if(power==4){
  	return (CF*pow(-1 + w,3)*(-CA + 8*CF + 8*CA*v - 52*CF*v - 34*CA*pow(v,2) + 160*CF*pow(v,2) + 42*CA*pow(v,3) - 212*CF*pow(v,3) + 4*CA*pow(v,4) + 176*CF*pow(v,4) + 8*pow(CF,2)*pow(v,4) + 16*CA*pow(v,5) - 28*CF*pow(v,5) - 16*pow(CF,2)*pow(v,5) + 10*CA*pow(v,6) - 16*CF*pow(v,6) + 8*pow(CF,2)*pow(v,6) + 12*CF*pow(v,4)*(3 + 2*v)*log(muR2/Q2) - 12*(CA - 2*CF)*pow(v,3)*(2 - 2*v + pow(v,2))*log(1 - v) - 36*CF*pow(v,4)*log(v) - 24*CF*pow(v,5)*log(v) + 24*CA*pow(v,3)*log(1 - w) - 48*CF*pow(v,3)*log(1 - w) - 24*CA*pow(v,4)*log(1 - w) + 12*CF*pow(v,4)*log(1 - w) + 12*CA*pow(v,5)*log(1 - w) - 48*CF*pow(v,5)*log(1 - w)))/(12.*CA*pow(-1 + v,4));
  }
  if(power==5){
  	return -(CF*pow(-1 + w,4)*(3*CA - 30*CF - 26*CA*v + 220*CF*v + 106*CA*pow(v,2) - 720*CF*pow(v,2) - 290*CA*pow(v,3) + 1440*CF*pow(v,3) + 345*CA*pow(v,4) - 1570*CF*pow(v,4) - 46*CA*pow(v,5) + 1192*CF*pow(v,5) + 40*pow(CF,2)*pow(v,5) + 122*CA*pow(v,6) - 224*CF*pow(v,6) - 80*pow(CF,2)*pow(v,6) + 38*CA*pow(v,7) - 56*CF*pow(v,7) + 40*pow(CF,2)*pow(v,7) + 120*CF*pow(v,5)*(2 + v)*log(muR2/Q2) - 60*(CA - 2*CF)*pow(v,4)*(2 - 2*v + pow(v,2))*log(1 - v) - 240*CF*pow(v,5)*log(v) - 120*CF*pow(v,6)*log(v) + 120*CA*pow(v,4)*log(1 - w) - 240*CF*pow(v,4)*log(1 - w) - 120*CA*pow(v,5)*log(1 - w) + 60*CA*pow(v,6)*log(1 - w) - 240*CF*pow(v,6)*log(1 - w)))/(60.*CA*pow(-1 + v,5));
  }
  if(power==6){
  	return (CF*pow(-1 + w,5)*(-2*CA + 24*CF + 19*CA*v - 198*CF*v - 83*CA*pow(v,2) + 728*CF*pow(v,2) + 224*CA*pow(v,3) - 1580*CF*pow(v,3) - 445*CA*pow(v,4) + 2320*CF*pow(v,4) + 477*CA*pow(v,5) - 2114*CF*pow(v,5) - 89*CA*pow(v,6) + 1456*CF*pow(v,6) + 36*pow(CF,2)*pow(v,6) + 148*CA*pow(v,7) - 272*CF*pow(v,7) - 72*pow(CF,2)*pow(v,7) + 31*CA*pow(v,8) - 44*CF*pow(v,8) + 36*pow(CF,2)*pow(v,8) + 60*CF*pow(v,6)*(5 + 2*v)*log(muR2/Q2) - 60*(CA - 2*CF)*pow(v,5)*(2 - 2*v + pow(v,2))*log(1 - v) - 300*CF*pow(v,6)*log(v) - 120*CF*pow(v,7)*log(v) + 120*CA*pow(v,5)*log(1 - w) - 240*CF*pow(v,5)*log(1 - w) - 120*CA*pow(v,6)*log(1 - w) - 60*CF*pow(v,6)*log(1 - w) + 60*CA*pow(v,7)*log(1 - w) - 240*CF*pow(v,7)*log(1 - w)))/(60.*CA*pow(-1 + v,6));
  }
  if(power==7){
  	return -(CF*pow(-1 + w,6)*(10*CA - 140*CF - 104*CA*v + 1288*CF*v + 496*CA*pow(v,2) - 5320*CF*pow(v,2) - 1442*CA*pow(v,3) + 13020*CF*pow(v,3) + 2891*CA*pow(v,4) - 21070*CF*pow(v,4) - 4480*CA*pow(v,5) + 24360*CF*pow(v,5) + 4284*CA*pow(v,6) - 18928*CF*pow(v,6) - 846*CA*pow(v,7) + 11884*CF*pow(v,7) + 224*pow(CF,2)*pow(v,7) + 1167*CA*pow(v,8) - 2138*CF*pow(v,8) - 448*pow(CF,2)*pow(v,8) + 184*CA*pow(v,9) - 256*CF*pow(v,9) + 224*pow(CF,2)*pow(v,9) + 840*CF*pow(v,7)*(3 + v)*log(muR2/Q2) - 420*(CA - 2*CF)*pow(v,6)*(2 - 2*v + pow(v,2))*log(1 - v) - 2520*CF*pow(v,7)*log(v) - 840*CF*pow(v,8)*log(v) + 840*CA*pow(v,6)*log(1 - w) - 1680*CF*pow(v,6)*log(1 - w) - 840*CA*pow(v,7)*log(1 - w) - 840*CF*pow(v,7)*log(1 - w) + 420*CA*pow(v,8)*log(1 - w) - 1680*CF*pow(v,8)*log(1 - w)))/(420.*CA*pow(-1 + v,7));
  }
  if(power==8){
  	return (CF*pow(-1 + w,7)*(-15*CA + 240*CF + 170*CA*v - 2440*CF*v - 886*CA*pow(v,2) + 11232*CF*pow(v,2) + 2816*CA*pow(v,3) - 30912*CF*pow(v,3) - 6118*CA*pow(v,4) + 56616*CF*pow(v,4) + 9730*CA*pow(v,5) - 72940*CF*pow(v,5) - 12320*CA*pow(v,6) + 69440*CF*pow(v,6) + 10536*CA*pow(v,7) - 46832*CF*pow(v,7) - 2046*CA*pow(v,8) + 26972*CF*pow(v,8) + 400*pow(CF,2)*pow(v,8) + 2538*CA*pow(v,9) - 4636*CF*pow(v,9) - 800*pow(CF,2)*pow(v,9) + 320*CA*pow(v,10) - 440*CF*pow(v,10) + 400*pow(CF,2)*pow(v,10) + 840*CF*pow(v,8)*(7 + 2*v)*log(muR2/Q2) - 840*(CA - 2*CF)*pow(v,7)*(2 - 2*v + pow(v,2))*log(1 - v) - 5880*CF*pow(v,8)*log(v) - 1680*CF*pow(v,9)*log(v) + 1680*CA*pow(v,7)*log(1 - w) - 3360*CF*pow(v,7)*log(1 - w) - 1680*CA*pow(v,8)*log(1 - w) - 2520*CF*pow(v,8)*log(1 - w) + 840*CA*pow(v,9)*log(1 - w) - 3360*CF*pow(v,9)*log(1 - w)))/(840.*CA*pow(-1 + v,8));
  }
  if(power==9){
  	return -(CF*pow(-1 + w,8)*(35*CA - 630*CF - 430*CA*v + 7020*CF*v + 2440*CA*pow(v,2) - 35700*CF*pow(v,2) - 8478*CA*pow(v,3) + 109536*CF*pow(v,3) + 20178*CA*pow(v,4) - 225876*CF*pow(v,4) - 34944*CA*pow(v,5) + 330288*CF*pow(v,5) + 45990*CA*pow(v,6) - 353220*CF*pow(v,6) - 49140*CA*pow(v,7) + 285600*CF*pow(v,7) + 37803*CA*pow(v,8) - 169686*CF*pow(v,8) - 7018*CA*pow(v,9) + 90176*CF*pow(v,9) + 1080*pow(CF,2)*pow(v,9) + 8114*CA*pow(v,10) - 14788*CF*pow(v,10) - 2160*pow(CF,2)*pow(v,10) + 850*CA*pow(v,11) - 1160*CF*pow(v,11) + 1080*pow(CF,2)*pow(v,11) + 5040*CF*pow(v,9)*(4 + v)*log(muR2/Q2) - 2520*(CA - 2*CF)*pow(v,8)*(2 - 2*v + pow(v,2))*log(1 - v) - 20160*CF*pow(v,9)*log(v) - 5040*CF*pow(v,10)*log(v) + 5040*CA*pow(v,8)*log(1 - w) - 10080*CF*pow(v,8)*log(1 - w) - 5040*CA*pow(v,9)*log(1 - w) - 10080*CF*pow(v,9)*log(1 - w) + 2520*CA*pow(v,10)*log(1 - w) - 10080*CF*pow(v,10)*log(1 - w)))/(2520.*CA*pow(-1 + v,9));
  }
  if(power==10){
  	return (CF*pow(-1 + w,9)*(-28*CA + 560*CF + 371*CA*v - 6790*CF*v - 2281*CA*pow(v,2) + 37840*CF*pow(v,2) + 8630*CA*pow(v,3) - 128280*CF*pow(v,3) - 22473*CA*pow(v,4) + 295056*CF*pow(v,4) + 42714*CA*pow(v,5) - 486276*CF*pow(v,5) - 61446*CA*pow(v,6) + 591528*CF*pow(v,6) + 68922*CA*pow(v,7) - 540540*CF*pow(v,7) - 63630*CA*pow(v,8) + 379680*CF*pow(v,8) + 44323*CA*pow(v,9) - 201206*CF*pow(v,9) - 7769*CA*pow(v,10) + 99188*CF*pow(v,10) + 980*pow(CF,2)*pow(v,10) + 8536*CA*pow(v,11) - 15532*CF*pow(v,11) - 1960*pow(CF,2)*pow(v,11) + 763*CA*pow(v,12) - 1036*CF*pow(v,12) + 980*pow(CF,2)*pow(v,12) + 2520*CF*pow(v,10)*(9 + 2*v)*log(muR2/Q2) - 2520*(CA - 2*CF)*pow(v,9)*(2 - 2*v + pow(v,2))*log(1 - v) - 22680*CF*pow(v,10)*log(v) - 5040*CF*pow(v,11)*log(v) + 5040*CA*pow(v,9)*log(1 - w) - 10080*CF*pow(v,9)*log(1 - w) - 5040*CA*pow(v,10)*log(1 - w) - 12600*CF*pow(v,10)*log(1 - w) + 2520*CA*pow(v,11)*log(1 - w) - 10080*CF*pow(v,11)*log(1 - w)))/(2520.*CA*pow(-1 + v,10));
  }
}
