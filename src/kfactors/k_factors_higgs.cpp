#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "k_factors_higgs.h"
using namespace std;

////////////////////////////////////////////////////////////
///
/// contains all K factors for higgs
/// split up in LO, NLO and power corrections
/// also have two routines: either fitted pdfs or real ones
///
////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////
/// LO
/////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
/// prefactor as given in https://arxiv.org/pdf/hep-ph/0207004.pdf, beneath eqn. 43
/// note that we replace v^2 with GF ( v^2 = 1/(Sqrt[2]*GF)) like in https://arxiv.org/pdf/0809.4283.pdf
/// this is sigma0*tau!
////////////////////////////////////////////////////////////////////////////////////////////////////
double higgs_LO_factor(){
	complex<double> AQtot = 0;
	for(int i = 0; i<2; i++){AQtot = AQtot + AQ(4.*pow(quarkmasses[i],2)/pow(mH,2));}
	AQtot = norm(AQtot);
	AQtot = 1;
	return pbunits*alphas_muR*alphas_muR*Q2*sqrt(2.)*GF/576./M_PI/S2;//*real(AQtot);
}
// this factor is checked, get the same pb as the results in 0809.4283 fig 1 (also changed higgs mass to check it)
complex<double> AQ(double x){
	if(x>=1){ return 3./2.*x*(1.+(1.-x)*pow(asin(1./sqrt(x)),2));}
	if(x<1){ return 3./2.*x*(1.+(1.-x)*-1./4.*pow(log((1.+sqrt(1.-x))/(1.-sqrt(1.-x))-I*M_PI),2));}
}

/////////////////////////////////////////////////////////////////////////////////////////
/// NLO
/////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////
/// gg channel
///////////////////////////
/// and the NLO functions
//https://arxiv.org/pdf/0809.4283.pdf eqn. 9 with extra 1/x (see eqn. 1)
double higgs_NLO_gg_reg(double x){
	return alphas_muR/M_PI*(6.*(1./x-2.+x-pow(x,2))*log(Q2/muF2)+(11.*pow(-1.+x,4)-24.*(-1.+x)*(-1.+x*(2.+(-1.+x)*x))*log(1.-x) + 12.*pow(1.+(-1.+x)*x,2)*log(x))/(2.*(-1.+x)*x));
}
//https://arxiv.org/pdf/0809.4283.pdf eqn. 7 and 8
double higgs_NLO_gg_plus(double x){
	return alphas_muR/M_PI*((12.*log(1.-x)+6.*log(Q2/muF2))/(1.-x));
}
//https://arxiv.org/pdf/0809.4283.pdf eqn. 7 (see eqn. 1)
double higgs_NLO_gg_delta(){
	return alphas_muR/M_PI*(11./2.+pow(M_PI,2));
}
double higgs_NLO_gg_expansion(double x, int power){
	if(power==1){
	return (-6*alphas_muR*(-1 + log(Q2/muF2) + 2*log(1 - x)))/M_PI;
}
if(power==2){
	return (-3*alphas_muR*(-1 + x)*(-1 + 4*log(Q2/muF2) + 8*log(1 - x)))/M_PI;
}
if(power==3){
	return (11*alphas_muR*pow(-1 + x,2))/M_PI;
}
if(power==4){
	return (-6*alphas_muR*pow(-1 + x,3)*(log(Q2/muF2) + 2*log(1 - x)))/M_PI;
}
if(power==5){
	return (3*alphas_muR*pow(-1 + x,4)*(7 + 10*log(Q2/muF2) + 20*log(1 - x)))/(5.*M_PI);
}
if(power==6){
	return (-3*alphas_muR*pow(-1 + x,5)*(21 + 20*log(Q2/muF2) + 40*log(1 - x)))/(10.*M_PI);
}
if(power==7){
	return (3*alphas_muR*pow(-1 + x,6)*(181 + 140*log(Q2/muF2) + 280*log(1 - x)))/(70.*M_PI);
}
if(power==8){
	return (-3*alphas_muR*pow(-1 + x,7)*(83 + 56*log(Q2/muF2) + 112*log(1 - x)))/(28.*M_PI);
}
if(power==9){
	return (alphas_muR*pow(-1 + x,8)*(4129 + 2520*log(Q2/muF2) + 5040*log(1 - x)))/(420.*M_PI);
}
if(power==10){
	return -(alphas_muR*pow(-1 + x,9)*(319 + 180*log(Q2/muF2) + 360*log(1 - x)))/(30.*M_PI);
}
if(power==11){
	return (alphas_muR*pow(-1 + x,10)*(22.671861471861472 + 12*log(Q2/muF2) + 24*log(1 - x)))/(2.*M_PI);
}
if(power==12){
	return (alphas_muR*pow(1 - x,11)*(23.923376623376623 + 12*log(Q2/muF2) + 24*log(1 - x)))/(2.*M_PI);
}
if(power==13){
	return (alphas_muR*pow(-1 + x,12)*(25.052514152514153 + 12*log(Q2/muF2) + 24*log(1 - x)))/(2.*M_PI);
}
if(power==14){
	return (alphas_muR*pow(1 - x,13)*(26.081684981684983 + 12*log(Q2/muF2) + 24*log(1 - x)))/(2.*M_PI);
}
if(power==15){
	return (alphas_muR*pow(-1 + x,14)*(27.02753912753913 + 12*log(Q2/muF2) + 24*log(1 - x)))/(2.*M_PI);
}
if(power==16){
	return (alphas_muR*pow(1 - x,15)*(27.902813852813853 + 12*log(Q2/muF2) + 24*log(1 - x)))/(2.*M_PI);
}
if(power==17){
	return (alphas_muR*pow(-1 + x,16)*(28.717487414546238 + 12*log(Q2/muF2) + 24*log(1 - x)))/(2.*M_PI);
}
if(power==18){
	return (alphas_muR*pow(1 - x,17)*(29.47953223247341 + 12*log(Q2/muF2) + 24*log(1 - x)))/(2.*M_PI);
}
if(power==19){
	return (alphas_muR*pow(-1 + x,18)*(30.195424905332025 + 12*log(Q2/muF2) + 24*log(1 - x)))/(2.*M_PI);
}
if(power==20){
	return (alphas_muR*pow(1 - x,19)*(30.87050230471283 + 12*log(Q2/muF2) + 24*log(1 - x)))/(2.*M_PI);
}
if(power==21){
	return (alphas_muR*pow(-1 + x,20)*(31.50921673785451 + 12*log(Q2/muF2) + 24*log(1 - x)))/(2.*M_PI);
}
return 0.;
}

///////////////////////////
/// qg channel
///////////////////////////
///// NLO functions
double higgs_NLO_qg_full(double x){
	return (alphas_muR*(-3 - (-6 + x)*x + 2*(2 + (-2 + x)*x)*log(Q2/muF2)
	+ 4*(2 + (-2 + x)*x)*log(1 - x) - 2*(2 + (-2 + x)*x)*log(x)))/(3.*M_PI*x);
}
double higgs_NLO_qg_expansion(double x, int power){
	if(power==1){
	return (2*alphas_muR*(1 + log(Q2/muF2) + 2*log(1 - x)))/(3.*M_PI);
}
if(power==2){
	return (-2*alphas_muR*(-1 + x)*(log(Q2/muF2) + 2*log(1 - x)))/(3.*M_PI);
}
if(power==3){
	return (4*alphas_muR*pow(-1 + x,2)*(log(Q2/muF2) + 2*log(1 - x)))/(3.*M_PI);
}
if(power==4){
	return (-4*alphas_muR*pow(-1 + x,3)*(2 + 3*log(Q2/muF2) + 6*log(1 - x)))/(9.*M_PI);
}
if(power==5){
	return (alphas_muR*pow(-1 + x,4)*(25 + 24*log(Q2/muF2) + 48*log(1 - x)))/(18.*M_PI);
}
if(power==6){
	return (alphas_muR*pow(1 - x,5)*(5.233333333333333 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==7){
	return (alphas_muR*pow(-1 + x,6)*(91 + 60*log(Q2/muF2) + 120*log(1 - x)))/(45.*M_PI);
}
if(power==8){
	return (alphas_muR*pow(1 - x,7)*(6.752380952380952 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==9){
	return (alphas_muR*pow(-1 + x,8)*(7.335714285714285 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==10){
	return (alphas_muR*pow(1 - x,9)*(7.843650793650793 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==11){
	return (alphas_muR*pow(-1 + x,10)*(8.293650793650794 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==12){
	return (alphas_muR*pow(1 - x,11)*(8.697691197691197 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==13){
	return (alphas_muR*pow(-1 + x,12)*(9.064357864357865 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==14){
	return (alphas_muR*pow(1 - x,13)*(9.4000222000222 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==15){
	return (alphas_muR*pow(-1 + x,14)*(9.70954600954601 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==16){
	return (alphas_muR*pow(1 - x,15)*(9.996725496725496 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==17){
	return (alphas_muR*pow(-1 + x,16)*(10.26458263958264 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==18){
	return (alphas_muR*pow(1 - x,17)*(10.515563031739502 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==19){
	return (alphas_muR*pow(-1 + x,18)*(10.751674142850613 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==20){
	return (alphas_muR*pow(1 - x,19)*(10.97458435956888 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
if(power==21){
	return (alphas_muR*pow(-1 + x,20)*(11.18569547067999 + 4*log(Q2/muF2) + 8*log(1 - x)))/(3.*M_PI);
}
return 0.;
}


///////////////////////////
/// qqbar channel
///////////////////////////
///// NLO functions
double higgs_NLO_qqbar_full(double x){
	return (-32*alphas_muR*pow(-1 + x,3))/(27.*M_PI*x);
}
double higgs_NLO_qqbar_expansion(double x, int power){
	if(power==1){
		return 0;
	}
	if(power==2){
		return 0;
	}
	if(power==3){
		return 0;
	}
	if(power==4){
		return (-32*alphas_muR*pow(-1 + x,3))/(27.*M_PI);
	}
	if(power==5){
		return (32*alphas_muR*pow(-1 + x,4))/(27.*M_PI);
	}
	if(power==6){
		return (-32*alphas_muR*pow(-1 + x,5))/(27.*M_PI);
	}
	if(power==7){
		return (32*alphas_muR*pow(-1 + x,6))/(27.*M_PI);
	}
	if(power==8){
		return (-32*alphas_muR*pow(-1 + x,7))/(27.*M_PI);
	}
	if(power==9){
		return (32*alphas_muR*pow(-1 + x,8))/(27.*M_PI);
	}
	if(power==10){
		return (-32*alphas_muR*pow(-1 + x,9))/(27.*M_PI);
	}
	if(power==11){
		return (32*alphas_muR*pow(-1 + x,10))/(27.*M_PI);
	}
	if(power==12){
		return (-32*alphas_muR*pow(-1 + x,11))/(27.*M_PI);
	}
	if(power==13){
		return (32*alphas_muR*pow(-1 + x,12))/(27.*M_PI);
	}
	if(power==14){
		return (-32*alphas_muR*pow(-1 + x,13))/(27.*M_PI);
	}
	if(power==15){
		return (32*alphas_muR*pow(-1 + x,14))/(27.*M_PI);
	}
	if(power==16){
		return (-32*alphas_muR*pow(-1 + x,15))/(27.*M_PI);
	}
	if(power==17){
		return (32*alphas_muR*pow(-1 + x,16))/(27.*M_PI);
	}
	if(power==18){
		return (-32*alphas_muR*pow(-1 + x,17))/(27.*M_PI);
	}
	if(power==19){
		return (32*alphas_muR*pow(-1 + x,18))/(27.*M_PI);
	}
	if(power==20){
		return (-32*alphas_muR*pow(-1 + x,19))/(27.*M_PI);
	}
	if(power==21){
		return (32*alphas_muR*pow(-1 + x,20))/(27.*M_PI);
	}
	return 0.;
}
