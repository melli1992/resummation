#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "k_factors_dy.h"
using namespace std;

//////////////////////////////////////////////////////////
///
/// contains all K factors for drell yan
/// split up in LO, NLO and the power expansion
///
//////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////
/// LO
/////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////
/// prefactor as given in Nucl.Phys. B359 (1991) 343-405, (A.1), alpha should be alphaEM
/// (tau*sigma0!)
////////////////////////////////////////////////////////////////////////////////////////////////////
double DY_LO_factor(){
	return pbunits*4.*M_PI*alphaEM*alphaEM/(3.*Q2*S2)*1./CA;
	//return 1;
}

/////////////////////////////////////////////////////////////////////////////////////////
/// NLO
/////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////
/// qqbar channel
///////////////////////////
/////////////////////////////////////////////////////////////////
/// Nucl.Phys. B359 (1991) 343-405 - B.4 - non constant piece
/////////////////////////////////////////////////////////////////
double DY_NLO_qqbar_reg(double x){
	return (alphas_muR*CF*(-4*(1 + x)*log(Q2/muF2) - 8*(1 + x)*log(1 - x) + (4*(1 + pow(x,2))*log(x))/(-1 + x)))/(4.*M_PI);
}
///////////////////////////////////////////////////////////
/// Nucl.Phys. B359 (1991) 343-405 - B.3 - plus dist
///////////////////////////////////////////////////////////
double DY_NLO_qqbar_plus(double x){
	return (2*alphas_muR*CF*(2*log(1.-x)/(1.-x) + log(Q2/muF2)/(1.-x)))/M_PI;
}
///////////////////////////////////////////////////////////
/// Nucl.Phys. B359 (1991) 343-405 - B.3 - delta(1-z) piece
///////////////////////////////////////////////////////////
double DY_NLO_qqbar_delta(){
	return (alphas_muR*CF*(-24 + 2*pow(M_PI,2) + 9*log(Q2/muF2)))/(6.*M_PI);
}
// power expansions
double DY_NLO_qqbar_expansion(double x, int power){
if(power==1){
	return (-2*alphas_muR*CF*(-1 + log(Q2/muF2) + 2*log(1 - x)))/M_PI;
}
if(power==2){
	return -((alphas_muR*CF*(-1 + x)*(-1 + log(Q2/muF2) + 2*log(1 - x)))/M_PI);
}
if(power==3){
	return (2*alphas_muR*CF*pow(-1 + x,2))/(3.*M_PI);
}
if(power==4){
	return -(alphas_muR*CF*pow(-1 + x,3))/(3.*M_PI);
}
if(power==5){
	return (7*alphas_muR*CF*pow(-1 + x,4))/(30.*M_PI);
}
if(power==6){
	return (-11*alphas_muR*CF*pow(-1 + x,5))/(60.*M_PI);
}
if(power==7){
	return (16*alphas_muR*CF*pow(-1 + x,6))/(105.*M_PI);
}
if(power==8){
	return (-11*alphas_muR*CF*pow(-1 + x,7))/(84.*M_PI);
}
if(power==9){
	return (29*alphas_muR*CF*pow(-1 + x,8))/(252.*M_PI);
}
if(power==10){
	return (-37*alphas_muR*CF*pow(-1 + x,9))/(360.*M_PI);
}
if(power==11){
	return (46*alphas_muR*CF*pow(-1 + x,10))/(495.*M_PI);
}
if(power==12){
	return (-14*alphas_muR*CF*pow(-1 + x,11))/(165.*M_PI);
}
if(power==13){
	return (67*alphas_muR*CF*pow(-1 + x,12))/(858.*M_PI);
}
if(power==14){
	return (-79*alphas_muR*CF*pow(-1 + x,13))/(1092.*M_PI);
}
if(power==15){
	return (92*alphas_muR*CF*pow(-1 + x,14))/(1365.*M_PI);
}
if(power==16){
	return (-53*alphas_muR*CF*pow(-1 + x,15))/(840.*M_PI);
}
if(power==17){
	return (121*alphas_muR*CF*pow(-1 + x,16))/(2040.*M_PI);
}
if(power==18){
	return (-137*alphas_muR*CF*pow(-1 + x,17))/(2448.*M_PI);
}
if(power==19){
	return (154*alphas_muR*CF*pow(-1 + x,18))/(2907.*M_PI);
}
if(power==20){
	return (-43*alphas_muR*CF*pow(-1 + x,19))/(855.*M_PI);
}
if(power==21){
	return (191*alphas_muR*CF*pow(-1 + x,20))/(3990.*M_PI);
}
}


///////////////////////////
/// qg channel
///////////////////////////
///// NLO functions
double DY_NLO_qg_full(double x){
	return (alphas_muR*TF*(1 + 6*x - 7*pow(x,2) + 2*(1 - 2*x + 2*pow(x,2))*(log(Q2/muF2) + 2*log(1 - x) - log(x))))/(4.*M_PI);
}
double DY_NLO_qg_expansion(double x, int power){
if(power==1){
	return (alphas_muR*TF*(log(Q2/muF2) + 2*log(1 - x)))/(2.*M_PI);
}
if(power==2){
	return (alphas_muR*TF*(-1 + x)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)))/(2.*M_PI);
}
if(power==3){
	return (alphas_muR*TF*pow(-1 + x,2)*(-5 + 2*log(Q2/muF2) + 4*log(1 - x)))/(2.*M_PI);
}
if(power==4){
	return (-2*alphas_muR*TF*pow(-1 + x,3))/(3.*M_PI);
}
if(power==5){
	return (7*alphas_muR*TF*pow(-1 + x,4))/(24.*M_PI);
}
if(power==6){
	return (-11*alphas_muR*TF*pow(-1 + x,5))/(60.*M_PI);
}
if(power==7){
	return (2*alphas_muR*TF*pow(-1 + x,6))/(15.*M_PI);
}
if(power==8){
	return (-11*alphas_muR*TF*pow(-1 + x,7))/(105.*M_PI);
}
if(power==9){
	return (29*alphas_muR*TF*pow(-1 + x,8))/(336.*M_PI);
}
if(power==10){
	return (-37*alphas_muR*TF*pow(-1 + x,9))/(504.*M_PI);
}
if(power==11){
	return (23*alphas_muR*TF*pow(-1 + x,10))/(360.*M_PI);
}
if(power==12){
	return (-28*alphas_muR*TF*pow(-1 + x,11))/(495.*M_PI);
}
if(power==13){
	return (67*alphas_muR*TF*pow(-1 + x,12))/(1320.*M_PI);
}
if(power==14){
	return (-79*alphas_muR*TF*pow(-1 + x,13))/(1716.*M_PI);
}
if(power==15){
	return (23*alphas_muR*TF*pow(-1 + x,14))/(546.*M_PI);
}
if(power==16){
	return (-53*alphas_muR*TF*pow(-1 + x,15))/(1365.*M_PI);
}
if(power==17){
	return (121*alphas_muR*TF*pow(-1 + x,16))/(3360.*M_PI);
}
if(power==18){
	return (-137*alphas_muR*TF*pow(-1 + x,17))/(4080.*M_PI);
}
if(power==19){
	return (77*alphas_muR*TF*pow(-1 + x,18))/(2448.*M_PI);
}
if(power==20){
	return (-86*alphas_muR*TF*pow(-1 + x,19))/(2907.*M_PI);
}
if(power==21){
	return (191*alphas_muR*TF*pow(-1 + x,20))/(6840.*M_PI);
}
}

////////////////////////////////
/// BSM part
////////////////////////////////
double BSM_WR_LO_factor(){
	return fbunits*M_PI*gR2/(MWR2)*1./(4.*CA);
}
