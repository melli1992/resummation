#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "clooptools.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "k_factors_diboson.h"
using namespace std;

////////////////////////////////////////////////////////////
///
/// contains all K factors for diboson (WpWm, ZZ)
/// split up in LO
///
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
/// W+W-
////////////////////////////////////////////////////////////

// needed for W+W- production
// Nucl Phys B.410 (1993) 280-324 Frixione
// couplings eqn 3.4 (and a part in parameters.cpp, 3.5 and 3.6)
double cts(double s, double qi, double giL){
	return pow(e,4)/(4.*s*pow(sinw,2))*(qi + 2.*ez*giL/pow(e,2)*s/(s-mZ2));
}
double css(double s, double qi, double giL, double giR){
	return pow(e,4)/pow(s,2)*(pow(qi+ez*(giL+giR)/pow(e,2)*s/(s-mZ2),2)+pow(ez*(giL-giR)/pow(e,2)*s/(s-mZ2),2));
}
// eqn 3.10 but error in the log((1+beta)/(1-beta)) term (compare with https://journals.aps.org/prd/pdf/10.1103/PhysRevD.19.922 eqn 3.12)
double Fu(double s, double beta){ //Fd = Fu
	return 16.*pow(s,3)*beta/(pow(mW2,2))*(1./24.+5.*mW2/(6.*s)-2.*pow(mW2,2)/pow(s,2))+16.*2.*s*(1.-2.*mW2/s)*log((1.+beta)/(1.-beta));
	//return pow(s,3)*beta/(pow(mW2,2))*(1./24.+5./6.*mW2/s-2.*pow(mW2,2)/pow(s,2))+2.*s*(1.-2.*mW2/s)*log((1+beta)/(1-beta));
}
double Ju(double s, double beta){ //Jd = -Ju
	return 16.*pow(s,3)*beta/(pow(mW2,2))*(s/24.+3.*mW2/(4.)-7.*pow(mW2,2)/(6.*s)-pow(mW2,3)/pow(s,2))-16.*4.*s*mW2*(1.+mW2/(2.*s))*log((1.+beta)/(1.-beta));
	// return pow(s,3)*beta/(pow(mW2,2))*(1./24.+3.*mW2/(4.*s)-7./6.*pow(mW2,2)/pow(s,2)-pow(mW2,3)/pow(s,3))-4.*mW2*(1.+1./2.*mW2/s)*log((1.+beta)/(1.-beta));
}
double Ku(double s, double beta){ // Kd = Ku
	return 8.*pow(s,3)*pow(beta,3)/(pow(mW2,2))*(pow(s,2)/24.+5.*s*mW2/(6.)+pow(mW2,2)/(2.));
}
double Acoupling(double s, double qi, double gv, double ga){
	return pow(qi+ez*(gv)/pow(e,2)*s/(s-mZ2),2)+ pow(ez*(ga)/pow(e,2)*s/(s-mZ2),2);
}
double Icoupling(double s, double qi, double gv, double ga){
	return 4.*(qi+ez*(gv+ga)/pow(e,2)*s/(s-mZ2))*pow(1./(2.*sqrt(2)*sinw),2);
}
double Ecoupling(){
	return 8.*pow(pow(1./(2.*sqrt(2)*sinw),2),2);
}
double Astu(double s, double beta){ //Fd = Fu
	return pow(s,3)*pow(beta,3)/(pow(mW2,2))*(1./24.+5.*mW2/(6.*s)+pow(mW2,2)/(2.*pow(s,2)));
}
double Istu(double s, double beta){ //Jd = -Ju
	return pow(s,3)*beta/(pow(mW2,2))*(1./24.+3.*mW2/(4.*s)-7./6.*pow(mW2,2)/pow(s,2)-pow(mW2,3)/pow(s,3))-4.*mW2*(1.+1./2.*mW2/s)*log((1.+beta)/(1.-beta));
}
double Estu(double s, double beta){ // Kd = Ku
	return pow(s,3)*beta/(pow(mW2,2))*(1./24.+5./6.*mW2/s-2.*pow(mW2,2)/pow(s,2))+2.*s*(1.-2.*mW2/s)*log((1+beta)/(1-beta));
}
// eqn 3.9
double partonic_up_wpwm(double s){
		if(s < 4.*mW2){return 0;}
		double beta = sqrt(1.-4.*mW2/s);
		//double result = pbunits/(64.*M_PI*CA*pow(s,2))*(ctt*Fu(s,beta)-cts(s, 2./3., guL)*Ju(s,beta) + css(s, 2./3., guL, guR)*Ku(s, beta));
		double result = pbunits*2.*M_PI*alphaEM*alphaEM/pow(s,2)*1./3.*(Ecoupling()*Estu(s, beta)-Icoupling(s, 2./3., gVu, gAu)*Istu(s, beta)+Astu(s, beta)*Acoupling(s, 2./3., gVu, gAu));
		return result; 
}
double partonic_down_wpwm(double s){
		if(s < 4.*mW2){return 0;}
		double beta = sqrt(1.-4.*mW2/s);
		//double result = pbunits/(64.*M_PI*CA*pow(s,2))*(ctt*Fu(s,beta)+cts(s, -1./3., gdL)*Ju(s, beta)+css(s, -1./3., gdL, gdR)*Ku(s, beta));
		double result =  pbunits*2.*M_PI*alphaEM*alphaEM/pow(s,2)*1./3.*(Ecoupling()*Estu(s, beta) + Icoupling(s, -1./3., gVd, gAd)*Istu(s, beta) +Astu(s, beta)*Acoupling(s, -1./3., gVd, gAd));
		return result;
}



complex<double>  cpartonic_up_wpwm(complex<double>  s){
		complex<double>  beta = sqrt(1.-4.*mW2/s);
		complex<double>  Ecoup = Ecoupling();
		complex<double>  ESTU = pow(s,3)*beta/(pow(mW2,2))*(1./24.+5./6.*mW2/s-2.*pow(mW2,2)/pow(s,2))+2.*s*(1.-2.*mW2/s)*log((1.+beta)/(1.-beta));
		complex<double>  Icoup = 4.*(2./3.+ez*(gVu+gAu)/pow(e,2)*s/(s-mZ2))*pow(1./(2.*sqrt(2)*sinw),2);
		complex<double>  ISTU = pow(s,3)*beta/(pow(mW2,2))*(1./24.+3.*mW2/(4.*s)-7./6.*pow(mW2,2)/pow(s,2)-pow(mW2,3)/pow(s,3))-4.*mW2*(1.+1./2.*mW2/s)*log((1.+beta)/(1.-beta));
		complex<double>  Acoup = pow(2./3.+ez*(gVu)/pow(e,2)*s/(s-mZ2),2)+ pow(ez*(gAu)/pow(e,2)*s/(s-mZ2),2);
		complex<double>  ASTU = pow(s,3)*pow(beta,3)/(pow(mW2,2))*(1./24.+5.*mW2/(6.*s)+pow(mW2,2)/(2.*pow(s,2)));
		
		complex<double>  result = pbunits*2.*M_PI*alphaEM*alphaEM/pow(s,2)*1./3.*(Ecoup*ESTU-Icoup*ISTU+ASTU*Acoup);
		return result; 
}
complex<double>  cpartonic_down_wpwm(complex<double> s){
		complex<double>  beta = sqrt(1.-4.*mW2/s);
		complex<double>  Ecoup = Ecoupling();
		complex<double>  ESTU = pow(s,3)*beta/(pow(mW2,2))*(1./24.+5./6.*mW2/s-2.*pow(mW2,2)/pow(s,2))+2.*s*(1.-2.*mW2/s)*log((1.+beta)/(1.-beta));
		complex<double>  Icoup = 4.*(-1./3.+ez*(gVd+gAd)/pow(e,2)*s/(s-mZ2))*pow(1./(2.*sqrt(2)*sinw),2);
		complex<double>  ISTU = pow(s,3)*beta/(pow(mW2,2))*(1./24.+3.*mW2/(4.*s)-7./6.*pow(mW2,2)/pow(s,2)-pow(mW2,3)/pow(s,3))-4.*mW2*(1.+1./2.*mW2/s)*log((1.+beta)/(1.-beta));
		complex<double>  Acoup = pow(-1./3.+ez*(gVd)/pow(e,2)*s/(s-mZ2),2)+ pow(ez*(gAd)/pow(e,2)*s/(s-mZ2),2);
		complex<double>  ASTU = pow(s,3)*pow(beta,3)/(pow(mW2,2))*(1./24.+5.*mW2/(6.*s)+pow(mW2,2)/(2.*pow(s,2)));
		complex<double>  result = pbunits*2.*M_PI*alphaEM*alphaEM/pow(s,2)*1./3.*(Ecoup*ESTU+Icoup*ISTU+ASTU*Acoup);
		return result;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
/// For ZZ production (see eqn 3.7 of https://journals.aps.org/prd/pdf/10.1103/PhysRevD.19.922)
/// factor of 1/3 is because of their 1/CA normalization (see eqn 5.10)
/////////////////////////////////////////////////////////////////////////////////////////////////
double partonic_up_zz(double s){
		if(s < 4.*mZ2){return 0;}
		double beta = sqrt(1.-4.*mZ2/s);
		double result = 1./3.*pbunits*4.*M_PI*alphaEM*alphaEM/s*(pow(gVu,4)+pow(gAu,4)+6.*pow(gVu,2)*pow(gAu,2))/pow(e,4)*((1.+4.*pow(mZ2,2)/pow(s,2))/(1.-2.*mZ2/s)*log((1.+beta)/(1.-beta))-beta);
		return result; 
}
double partonic_down_zz(double s){
		if(s < 4.*mZ2){return 0;}
		double beta = sqrt(1.-4.*mZ2/s);
		double result = 1./3.*pbunits*4.*M_PI*alphaEM*alphaEM/s*(pow(gVd,4)+pow(gAd,4)+6.*pow(gVd,2)*pow(gAd,2))/pow(e,4)*((1.+4.*pow(mZ2,2)/pow(s,2))/(1.-2.*mZ2/s)*log((1+beta)/(1-beta))-beta);
		return result; 
}


complex<double> cpartonic_up_zz(complex<double> s){
		complex<double> beta = sqrt(1.-4.*mZ2/s);
		complex<double> result = 1./3.*pbunits*4.*M_PI*alphaEM*alphaEM/s*(pow(gVu,4)+pow(gAu,4)+6.*pow(gVu,2)*pow(gAu,2))/pow(e,4)*((1.+4.*pow(mZ2,2)/pow(s,2))/(1.-2.*mZ2/s)*log((1.+beta)/(1.-beta))-beta);
		return result; 
}
complex<double> cpartonic_down_zz(complex<double> s){
		complex<double> beta = sqrt(1.-4.*mZ2/s);
		complex<double> result = 1./3.*pbunits*4.*M_PI*alphaEM*alphaEM/s*(pow(gVd,4)+pow(gAd,4)+6.*pow(gVd,2)*pow(gAd,2))/pow(e,4)*((1.+4.*pow(mZ2,2)/pow(s,2))/(1.-2.*mZ2/s)*log((1.+beta)/(1.-beta))-beta);
		return result; 
}
