#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include "parameters.h"
#include "SCET_functions.h"
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

complex<double> g1_TROLL(double A1,complex<double>lambda){
	A1 = A1/M_PI;
	complex<double> z = -2.*lambda;
	return (log(z+1.)*(z+1.)-z)*A1/b0/b0;
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
	if(INCEULER == 0) {INCeuler = 1.;}
	return 1./(2.*M_PI*b0)*(-A2/(M_PI*b0)+A1*log(Q2/muR2))*(2.*lambda+log(1.-2.*lambda))
	+ A1*b1/(2.*M_PI*pow(b0,3))*(2.*lambda+log(1.-2.*lambda)+1./2.*pow(log(1.-2.*lambda),2))
	- A1/(M_PI*b0)*lambda*log(Q2/muF2)
	- INCeuler*2.*M_gammaE*log(1.-2.*lambda)*A1/(2.*M_PI*b0);
}

complex<double> g2_TROLL(double A1,double A2,complex<double>lambda){
A1 = A1/M_PI;
A2 = A2/M_PI/M_PI;
double D1 = 0.;
double EG = M_gammaE;
if(INCEULER==1) EG = 0.;
complex<double> z = -2.*lambda;
return ((pow(log(z+1.),2)*0.5 - z + log(z+1.) ) * ( A1*b1/b0/b0/b0 ) +
	    (z - log(z+1.) ) * A2/b0/b0
	    - (log(z+1.) * (EG - log(Q/muR)) + z*log(muF/muR) )*2.*A1/b0);
}
////////////////////////////////////////////////////////////////////////
// LP NNLL function g3 hep-ph/0306211 eqn 41 (checked with mathematica)
// checked with Richards code as well (for Nbar => N)
// checked also the output numerically with mathematica
complex<double> g3(double A1,double A2,double A3,complex<double>lambda){
	double INCeuler = 0.;
	double Lqr = log(Q2/muR2);
	double Lfr = log(muF2/muR2);
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
	-2.*A2/(pow(M_PI,2)*b0)*pow(lambda,2))*log(Q2/muR2)
	+A1/M_PI*pow(lambda,2)/(1.-2.*lambda)*pow(log(Q2/muR2),2)
	+ b0*((A1*b1/(M_PI*pow(b0,3))*INCeuler*M_gammaE)*(-2.*lambda/(1.-2.*lambda)-log(1.-2.*lambda)/(1.-2.*lambda))
- (A1/(M_PI*b0)*(INCeuler*M_gammaE*M_gammaE)+A2/(M_PI*M_PI*pow(b0,2))*INCeuler*M_gammaE)*-2.*lambda/(1.-2.*lambda)
+ ((A1/(M_PI*b0)*INCeuler*M_gammaE)*(-2.*lambda/(1.-2.*lambda)))*log(Q2/muR2));
}

complex<double> g3_TROLL(double A1,double A2,double A3,double D2, complex<double>lambda){
// note their b1 = our b1/b0 and their b2 = our b2/b0

A1 = A1/M_PI;
A2 = A2/M_PI/M_PI;
A3 = A3/M_PI/M_PI/M_PI;
D2 = D2/M_PI/M_PI;
double D1 = 0.;
double EG = M_gammaE;
if(INCEULER==1) EG = 0.;
complex<double> z = -2.*lambda;
     return ( (z*z/(z+1.)) * (A3 - A1*b2/b0 + A1*b1/b0*b1/b0 - A2*b1/b0)/2./b0/b0/b0
	     + ( log(z+1.)/(z+1.) + pow(log(z+1.),2.)/(z+1.)/2. ) * A1*b1/b0*b1/b0/b0/b0/b0
	     + log(z+1.) * A1*(b2/b0-b1/b0*b1/b0)/b0/b0/b0
	     + ( z/(z+1.) - log(z+1.)/(z+1.) ) * ( 2.*A1*b1/b0*EG + A2*b1/b0/b0 - D1/2.*b1/b0 )/b0/b0
	     - ( z/(z+1.) ) * ( A1*b2/b0/b0/b0 + 2.*A1*(EG*EG + zeta2) + 2.*A2*EG/b0 )/b0
	     + ( z/(z+1.) ) * ( D2/2./b0/b0 + D1*EG/b0 )
	     + ( ( z/(z+1.) ) * ( 2.*A1*EG + A2/b0 -A1*b1/b0/b0 -D1/2. )/b0
	      + log(z+1.)/(z+1.) * A1*b1/b0/b0/b0 ) * 2.*log(Q/muR)
	     - z * A2/b0/b0 * 2.*log(muF/muR)
	     + ( z * pow(log(muF/muR),2.) - z/(z+1.) * pow(log(Q/muR),2.) ) * 2.*A1/b0
	     )*b0;
}

//taken from TROLL
complex<double> g4(double A1, double A2, double A3, double A4, double D2, double D3, complex<double> lambda){
	A1 = A1/M_PI;
	A2 = A2/M_PI/M_PI;
	A3 = A3/M_PI/M_PI/M_PI;
	D2 = D2/M_PI/M_PI;
	D3 = D3/M_PI/M_PI/M_PI;
	double D1 = 0.;
	double EG = M_gammaE;
	if(INCEULER==1) EG = 0.;
	complex<double> z = -2.*lambda;
	return -(-6*A3*b1/b0*z + 6*A2*pow(b1/b0,2)*z - 6*A1*b1/b0*b2/b0*z + 6*A1*b3/b0*z + 6*b0*b1/b0*D2*z - 6*b0*D3*z +
      24*A3*b0*EG*z - 24*A2*b0*b1/b0*EG*z - 24*pow(b0,2)*D2*EG*z +
      48*A2*pow(b0,2)*pow(EG,2)*z - 24*pow(b0,3)*D1*pow(EG,2)*z +
      32*A1*pow(b0,3)*pow(EG,3)*z - 6*A4*pow(z,2) + 3*A3*b1/b0*pow(z,2) +
      3*A2*pow(b1/b0,2)*pow(z,2) - 9*A1*b1/b0*b2/b0*pow(z,2) + 9*A1*b3/b0*pow(z,2) -
      3*b0*pow(b1/b0,2)*D1*pow(z,2) + 3*b0*b2/b0*D1*pow(z,2) + 3*b0*b1/b0*D2*pow(z,2) -
      3*b0*D3*pow(z,2) + 12*A3*b0*EG*pow(z,2) - 12*A2*b0*b1/b0*EG*pow(z,2) +
      12*A1*b0*pow(b1/b0,2)*EG*pow(z,2) - 12*A1*b0*b2/b0*EG*pow(z,2) -
      12*pow(b0,2)*D2*EG*pow(z,2) + 24*A2*pow(b0,2)*pow(EG,2)*pow(z,2) -
      12*pow(b0,3)*D1*pow(EG,2)*pow(z,2) +
      16*A1*pow(b0,3)*pow(EG,3)*pow(z,2) - 4*A4*pow(z,3) + 4*A3*b1/b0*pow(z,3) -
      4*A2*pow(b1/b0,2)*pow(z,3) + 4*A1*pow(b1/b0,3)*pow(z,3) + 4*A2*b2/b0*pow(z,3) -
      8*A1*b1/b0*b2/b0*pow(z,3) + 4*A1*b3/b0*pow(z,3) + 48*A2*pow(b0,2)*z*zeta2 -
      24*pow(b0,3)*D1*z*zeta2 + 96*A1*pow(b0,3)*EG*z*zeta2 +
      24*A2*pow(b0,2)*pow(z,2)*zeta2 - 12*pow(b0,3)*D1*pow(z,2)*zeta2 +
      48*A1*pow(b0,3)*EG*pow(z,2)*zeta2 + 64*A1*pow(b0,3)*z*zeta3 +
      32*A1*pow(b0,3)*pow(z,2)*zeta3 - 16*A1*pow(b0,3)*z*(2.+z)*pow(log(Q/muR),3) +
      24*A3*b0*z*pow(1.+z,2)*log(muF/muR) -
      24*pow(b0,2)*(2*A2 + A1*b1/b0)*z*pow(1.+z,2)*pow(log(muF/muR),2) +
      32*A1*pow(b0,3)*z*pow(1.+z,2)*pow(log(muF/muR),3) + 6*A3*b1/b0*log(1.+z) -
      6*A2*pow(b1/b0,2)*log(1.+z) + 6*A1*b1/b0*b2/b0*log(1.+z) - 6*A1*b3/b0*log(1.+z) -
      6*b0*b1/b0*D2*log(1.+z) + 24*A2*b0*b1/b0*EG*log(1.+z) -
      12*pow(b0,2)*b1/b0*D1*EG*log(1.+z) +
      24*A1*pow(b0,2)*b1/b0*pow(EG,2)*log(1.+z) + 12*A1*b1/b0*b2/b0*z*log(1.+z) -
      12*A1*b3/b0*z*log(1.+z) - 6*A1*pow(b1/b0,3)*pow(z,2)*log(1.+z) +
      12*A1*b1/b0*b2/b0*pow(z,2)*log(1.+z) - 6*A1*b3/b0*pow(z,2)*log(1.+z) +
      24*A1*pow(b0,2)*b1/b0*zeta2*log(1.+z) - 6*A2*pow(b1/b0,2)*pow(log(1.+z),2) +
      3*b0*pow(b1/b0,2)*D1*pow(log(1.+z),2) - 12*A1*b0*pow(b1/b0,2)*EG*pow(log(1.+z),2) +
      2*A1*pow(b1/b0,3)*pow(log(1.+z),3) +
      12*pow(b0,2)*pow(log(Q/muR),2)*((2*A2 - b0*D1 + 4*A1*b0*EG)*z*(2.+z) +
         2*A1*b1/b0*log(1.+z)) - 12*b0*log(Q/muR)*
       (z*(-2*b0*D2 - 4*pow(b0,2)*D1*EG + 8*A1*pow(b0,2)*pow(EG,2) +
            A1*pow(b1/b0,2)*z - A1*b2/b0*z - b0*D2*z - 2*pow(b0,2)*D1*EG*z +
            4*A1*pow(b0,2)*pow(EG,2)*z + A3*(2.+z) - A2*(b1/b0 - 4*b0*EG)*(2.+z) +
            8*A1*pow(b0,2)*zeta2 + 4*A1*pow(b0,2)*z*zeta2) +
         b1/b0*(2*A2 - b0*D1 + 4*A1*b0*EG)*log(1.+z) - A1*pow(b1/b0,2)*pow(log(1.+z),2)))/
   (12.*pow(b0,2)*pow(1.+z,2));
}

// wide angle contribution
complex<double> wideangle(double D2,complex<double>lambda){
	return -D2/(pow(M_PI,2)*b0)*lambda/(1.-2.*lambda);
}

///////////////////////////////////////////////////////////////////
// matching functions to NLO and NNLO (checked with mathematica)
// obtained from expanding D12 like in ref below
// see https://arxiv.org/pdf/1009.5691.pdf (appendix D13 en 14)
// checked also the output numerically with mathematica

complex<double> NLOmatch(complex<double> N, double A1, double g01){
double as = alphas_muR;
double LQmuF2 = log(Q2/muF2);
complex<double> lnN = log(N*exp(M_gammaE*INCEULER));
double EG = 0.;
if(INCEULER == 0) EG = M_gammaE;
complex<double> LPpiece = 4.*A1/M_PI*as*EG*ISNLL*lnN
													+ as*((2.*A1*ISLL*pow(lnN,2))/M_PI
													- (2.*A1*ISNLL*lnN*LQmuF2)/M_PI)
													+ ISNNLL*g01
 												  + 1.;
complex<double> NLPpiece = 2.*((A1*as*ISNLP*lnN)/(M_PI*N));
//cout <<
return LPpiece + NLPpiece;
}

complex<double> NNLOmatch(complex<double> N, double A1, double A2, double D2, double g01, double g02){
double as = alphas_muR;
double LQmuF2 = log(Q2/muF2);
double LQmuR2 = log(Q2/muR2);
double LmuF2muR2 = log(muF2/muR2);
complex<double> lnN = log(N*exp(M_gammaE*INCEULER));
double EG = 0.;
if(INCEULER == 0) EG = M_gammaE;
complex<double> asS1 = (2.*A1*as*ISLL*lnN*lnN)/M_PI - (2.*A1*as*ISNLL*lnN*(-2.*EG + LQmuF2))/M_PI;
complex<double> as2S2 = 4.*pow(as,2)*pow(b0,2)*ISNLL*pow(lnN,2)*(A2/(2.*pow(b0,2)*pow(M_PI,2)) + (A1*(EG - LQmuR2/2.))/(b0*M_PI))
								+ (4.*A1*pow(as,2)*b0*ISLL*pow(lnN,3))/(3.*M_PI)
								- 2.*pow(as,2)*b0*ISNNLL*lnN*(D2/(2.*b0)/pow(M_PI,2) - (2.*A2*EG)/(b0*pow(M_PI,2)) + (A2*LQmuF2)/(b0*pow(M_PI,2)) + (2.*A1*EG*LQmuR2)/M_PI - (2.*A1*(-pow(LmuF2muR2,2)/4. + pow(LQmuR2,2)/4.))/M_PI - (2.*A1*(pow(EG,2) + zeta2))/M_PI);
//cout << "./lnN->" << lnN << "/.A1->"<< A1<<"/.A2->"<< A2<< "/.D2->"<< D2<<"/.b0->"<< b0<<"/.LQmuF2->" << LQmuF2 << "/.g01->" << g01/as*M_PI << "/.g02->" << g02/as/as*M_PI*M_PI << "/.as->" << as << endl;
complex<double> LPpiece = pow(asS1,2)/2. + as2S2 + g01*asS1 + g02
													+ asS1 + g01
													+ 1.;
//fix this with EG pieces!
complex<double> NLPpiece =
ISNLP*(as*as*((4.*A1*A1*ISLL*pow(lnN,3))/(pow(M_PI,2)*N)
						- (4.*A1*A1*ISNLL*pow(lnN,2)*LQmuF2)/(pow(M_PI,2)*N)
						+ (2.*A1*b0*pow(lnN,2))/(M_PI*N)
						+ (2.*A1*g01/as*M_PI*ISNNLL*lnN)/(pow(M_PI,2)*N)
					  + (2.*A1*A1*pow(lnN,2))/(pow(M_PI,2)*pow(N,2)))
						+ (2.*A1*as*lnN)/(M_PI*N));
return LPpiece+NLPpiece;
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

complex<double> DY_NLP_LL_qg(int power, complex<double> n){
	//double INCeuler = 0.;
	//if(INCEULER == 0) INCeuler = 1.;
	double as = alphas_muR;
	complex<double> lnN = log(n*exp(INCEULER*M_gammaE));
	if(power==0){
		return 0;
	}
	if(power==1){
		return -((as*lnN*TF)/(n*M_PI));
	}
	if(power==2){
		return (pow(as,2)*(-13*CA - 35*CF)*pow(lnN,3)*TF)/(24.*n*pow(M_PI,2));
	}
	if(power==3){
		return (pow(as,3)*(-27*pow(CA,2) - 50*CA*CF - 115*pow(CF,2))*pow(lnN,5)*TF)/(96.*n*pow(M_PI,3));
	}
	if(power==4){
		return (pow(as,4)*(-pow(CA - CF,3)/34560. + ((CA - CF)*pow(CF,2))/24. + pow(CF,4)/(3.*(CA - CF)) + (pow(CF,3)*(-CA + CF))/(3.*(CA - CF)) - pow(3*CA + CF,4)/(768.*(CA - CF)))*pow(lnN,7)*TF)/(n*pow(M_PI,4));
	}
	if(power==5){
		return (pow(as,5)*(-(pow(CA - CF,3)*CF)/17280. + ((CA - CF)*pow(CF,3))/36. + (2*pow(CF,5))/(15.*(CA - CF)) + (pow(CF,4)*(-CA + CF))/(6.*(CA - CF)) - pow(3*CA + CF,5)/(7680.*(CA - CF)))*pow(lnN,9)*TF)/(n*pow(M_PI,5));
	}
	if(power==6){
		return (pow(as,6)*(pow(CA - CF,5)/4.35456e7 - (pow(CA - CF,3)*pow(CF,2))/17280. + ((CA - CF)*pow(CF,4))/72. + (2*pow(CF,6))/(45.*(CA - CF)) + (pow(CF,5)*(-CA + CF))/(15.*(CA - CF)) - pow(3*CA + CF,6)/(92160.*(CA - CF)))*pow(lnN,11)*TF)/(n*pow(M_PI,6));
	}
	if(power==7){
		return (pow(as,7)*((pow(CA - CF,5)*CF)/2.17728e7 - (pow(CA - CF,3)*pow(CF,3))/25920. + ((CA - CF)*pow(CF,5))/180. + (4*pow(CF,7))/(315.*(CA - CF)) + (pow(CF,6)*(-CA + CF))/(45.*(CA - CF)) - pow(3*CA + CF,7)/(1.29024e6*(CA - CF)))*pow(lnN,13)*TF)/(n*pow(M_PI,7));
	}
	if(power==8){
		return (pow(as,8)*(-pow(CA - CF,7)/9.7542144e10 + (pow(CA - CF,5)*pow(CF,2))/2.17728e7 - (pow(CA - CF,3)*pow(CF,4))/51840. + ((CA - CF)*pow(CF,6))/540. + pow(CF,8)/(315.*(CA - CF)) + (2*pow(CF,7)*(-CA + CF))/(315.*(CA - CF)) - pow(3*CA + CF,8)/(2.064384e7*(CA - CF)))*pow(lnN,15)*TF)/(n*pow(M_PI,8));
	}
	if(power==9){
		return (pow(as,9)*(-(pow(CA - CF,7)*CF)/4.8771072e10 + (pow(CA - CF,5)*pow(CF,3))/3.26592e7 - (pow(CA - CF,3)*pow(CF,5))/129600. + ((CA - CF)*pow(CF,7))/1890. + (2*pow(CF,9))/(2835.*(CA - CF)) + (pow(CF,8)*(-CA + CF))/(630.*(CA - CF)) - pow(3*CA + CF,9)/(3.7158912e8*(CA - CF)))*pow(lnN,17)*TF)/(n*pow(M_PI,9));
	}
	if(power==10){
		return (pow(as,10)*(pow(CA - CF,9)/3.47640201216e14 - (pow(CA - CF,7)*pow(CF,2))/4.8771072e10 + (pow(CA - CF,5)*pow(CF,4))/6.53184e7 - (pow(CA - CF,3)*pow(CF,6))/388800. + ((CA - CF)*pow(CF,8))/7560. + (2*pow(CF,10))/(14175.*(CA - CF)) + (pow(CF,9)*(-CA + CF))/(2835.*(CA - CF)) - pow(3*CA + CF,10)/(7.4317824e9*(CA - CF)))*pow(lnN,19)*TF)/(n*pow(M_PI,10));
	}
	return 0.;
}


complex<double> higgs_NLP_LL_qg(int power, complex<double> n){
	complex<double> lnN = log(n*exp(INCEULER*M_gammaE));
	double as = alphas_muR;
	if(power==0){
		return 0.;
	}
	if(power==1){
		return -((as*CF*lnN)/(n*M_PI));
	}
	if(power==2){
		return -(pow(as,2)*CF*(35*CA + 13*CF)*pow(lnN,3))/(24.*n*pow(M_PI,2));
	}
	if(power==3){
		return -(pow(as,3)*CF*(115*pow(CA,2) + 50*CA*CF + 27*pow(CF,2))*pow(lnN,5))/(96.*n*pow(M_PI,3));
	}
	if(power==4){
		return (pow(as,4)*((pow(CA,4)*CF)/(3.*(-CA + CF)) + (pow(CA,3)*(CA - CF)*CF)/(3.*(-CA + CF)) + (pow(CA,2)*CF*(-CA + CF))/24. - (CF*pow(-CA + CF,3))/34560. - (CF*pow(CA + 3*CF,4))/(768.*(-CA + CF)))*pow(lnN,7))/(n*pow(M_PI,4));
	}
	if(power==5){
		return (pow(as,5)*((2*pow(CA,5)*CF)/(15.*(-CA + CF)) + (pow(CA,4)*(CA - CF)*CF)/(6.*(-CA + CF)) + (pow(CA,3)*CF*(-CA + CF))/36. - (CA*CF*pow(-CA + CF,3))/17280. - (CF*pow(CA + 3*CF,5))/(7680.*(-CA + CF)))*pow(lnN,9))/(n*pow(M_PI,5));
	}
	if(power==6){
		return (pow(as,6)*((2*pow(CA,6)*CF)/(45.*(-CA + CF)) + (pow(CA,5)*(CA - CF)*CF)/(15.*(-CA + CF)) + (pow(CA,4)*CF*(-CA + CF))/72. - (pow(CA,2)*CF*pow(-CA + CF,3))/17280. + (CF*pow(-CA + CF,5))/4.35456e7 - (CF*pow(CA + 3*CF,6))/(92160.*(-CA + CF)))*pow(lnN,11))/(n*pow(M_PI,6));
	}
	if(power==7){
		return (pow(as,7)*((4*pow(CA,7)*CF)/(315.*(-CA + CF)) + (pow(CA,6)*(CA - CF)*CF)/(45.*(-CA + CF)) + (pow(CA,5)*CF*(-CA + CF))/180. - (pow(CA,3)*CF*pow(-CA + CF,3))/25920. + (CA*CF*pow(-CA + CF,5))/2.17728e7 - (CF*pow(CA + 3*CF,7))/(1.29024e6*(-CA + CF)))*pow(lnN,13))/(n*pow(M_PI,7));
	}
	if(power==8){
		return (pow(as,8)*((pow(CA,8)*CF)/(315.*(-CA + CF)) + (2*pow(CA,7)*(CA - CF)*CF)/(315.*(-CA + CF)) + (pow(CA,6)*CF*(-CA + CF))/540. - (pow(CA,4)*CF*pow(-CA + CF,3))/51840. + (pow(CA,2)*CF*pow(-CA + CF,5))/2.17728e7 - (CF*pow(-CA + CF,7))/9.7542144e10 - (CF*pow(CA + 3*CF,8))/(2.064384e7*(-CA + CF)))*pow(lnN,15))/(n*pow(M_PI,8));
	}
	if(power==9){
		return (pow(as,9)*((2*pow(CA,9)*CF)/(2835.*(-CA + CF)) + (pow(CA,8)*(CA - CF)*CF)/(630.*(-CA + CF)) + (pow(CA,7)*CF*(-CA + CF))/1890. - (pow(CA,5)*CF*pow(-CA + CF,3))/129600. + (pow(CA,3)*CF*pow(-CA + CF,5))/3.26592e7 - (CA*CF*pow(-CA + CF,7))/4.8771072e10 - (CF*pow(CA + 3*CF,9))/(3.7158912e8*(-CA + CF)))*pow(lnN,17))/(n*pow(M_PI,9));
	}
	if(power==10){
		return (pow(as,10)*((2*pow(CA,10)*CF)/(14175.*(-CA + CF)) + (pow(CA,9)*(CA - CF)*CF)/(2835.*(-CA + CF)) + (pow(CA,8)*CF*(-CA + CF))/7560. - (pow(CA,6)*CF*pow(-CA + CF,3))/388800. + (pow(CA,4)*CF*pow(-CA + CF,5))/6.53184e7 - (pow(CA,2)*CF*pow(-CA + CF,7))/4.8771072e10 + (CF*pow(-CA + CF,9))/3.47640201216e14 - (CF*pow(CA + 3*CF,10))/(7.4317824e9*(-CA + CF)))*pow(lnN,19))/(n*pow(M_PI,10));
	}
	return 0.;
}

complex<double> NLP_expansions(int power, double A1, complex<double> n){
	complex<double> lnN = log(n*exp(INCEULER*M_gammaE));
	double as = alphas_muR;
	if(power==1){
	return 1;
}
if(power==2){
	return (2*A1*as*lnN)/(n*M_PI);
}
if(power==3){
	return (2*A1*pow(as,2)*pow(lnN,2)*(A1 + b0*n*M_PI))/(pow(n,2)*pow(M_PI,2));
}
if(power==4){
	return (4*A1*pow(as,3)*pow(lnN,3)*(A1 + b0*n*M_PI)*(A1 + 2*b0*n*M_PI))/(3.*pow(n,3)*pow(M_PI,3));
}
if(power==5){
	return (2*A1*pow(as,4)*pow(lnN,4)*(A1 + b0*n*M_PI)*(A1 + 2*b0*n*M_PI)*(A1 + 3*b0*n*M_PI))/(3.*pow(n,4)*pow(M_PI,4));
}
if(power==6){
	return (4*A1*pow(as,5)*pow(lnN,5)*(A1 + b0*n*M_PI)*(A1 + 2*b0*n*M_PI)*(A1 + 3*b0*n*M_PI)*(A1 + 4*b0*n*M_PI))/(15.*pow(n,5)*pow(M_PI,5));
}
if(power==7){
	return (4*A1*pow(as,6)*pow(lnN,6)*(A1 + b0*n*M_PI)*(A1 + 2*b0*n*M_PI)*(A1 + 3*b0*n*M_PI)*(A1 + 4*b0*n*M_PI)*(A1 + 5*b0*n*M_PI))/(45.*pow(n,6)*pow(M_PI,6));
}
if(power==8){
	return (8*A1*pow(as,7)*pow(lnN,7)*(A1 + b0*n*M_PI)*(A1 + 2*b0*n*M_PI)*(A1 + 3*b0*n*M_PI)*(A1 + 4*b0*n*M_PI)*(A1 + 5*b0*n*M_PI)*(A1 + 6*b0*n*M_PI))/(315.*pow(n,7)*pow(M_PI,7));
}
if(power==9){
	return (2*A1*pow(as,8)*pow(lnN,8)*(A1 + b0*n*M_PI)*(A1 + 2*b0*n*M_PI)*(A1 + 3*b0*n*M_PI)*(A1 + 4*b0*n*M_PI)*(A1 + 5*b0*n*M_PI)*(A1 + 6*b0*n*M_PI)*(A1 + 7*b0*n*M_PI))/(315.*pow(n,8)*pow(M_PI,8));
}
if(power==10){
	return (4*A1*pow(as,9)*pow(lnN,9)*(A1 + b0*n*M_PI)*(A1 + 2*b0*n*M_PI)*(A1 + 3*b0*n*M_PI)*(A1 + 4*b0*n*M_PI)*(A1 + 5*b0*n*M_PI)*(A1 + 6*b0*n*M_PI)*(A1 + 7*b0*n*M_PI)*(A1 + 8*b0*n*M_PI))/(2835.*pow(n,9)*pow(M_PI,9));
}
if(power==11){
	return (4*A1*pow(as,10)*pow(lnN,10)*(A1 + b0*n*M_PI)*(A1 + 2*b0*n*M_PI)*(A1 + 3*b0*n*M_PI)*(A1 + 4*b0*n*M_PI)*(A1 + 5*b0*n*M_PI)*(A1 + 6*b0*n*M_PI)*(A1 + 7*b0*n*M_PI)*(A1 + 8*b0*n*M_PI)*(A1 + 9*b0*n*M_PI))/(14175.*pow(n,10)*pow(M_PI,10));
}

	return 0;

}
// higgs
// see https://arxiv.org/pdf/hep-ph/0306211.pdf, eqn 44 and
// or https://arxiv.org/pdf/1206.4133.pdf eqn 2
// checked with TROLL implementation - it is the same
double higgs_g01(){
	double ZETA2 = zeta2;
	double A1 = A1g/M_PI;
  double lf  = log(Q2/muF2);
  double lr = log(muR2/Q2);
  double EulerGamma = M_gammaE;
  if(INCEULER==1) EulerGamma=0.;
  double g01 = (2*(b0*M_PI)*lf+2.*CA * ZETA2 + 11./2.+2*(b0*M_PI)*lr)/M_PI + 2*A1* (ZETA2+EulerGamma*EulerGamma - EulerGamma*lf);
  return alphas_muR*g01;
}

//TAKEN from TROLL
double higgs_g02(){
	double EulerGamma = M_gammaE;
	if(INCEULER==1) EulerGamma=0.;
	double ZETA2 = zeta2;
	double ZETA3 = zeta3;
	double ZETA4 = pow(M_PI,4)/90.;
  double A1 = A1g/M_PI;
  double A2 = A2g/M_PI/M_PI;
  double D2 = D2higgs/M_PI/M_PI;
	double lHt = log(Q2/mt2);
	double lf  = log(Q2/muF2);
  double lrf = log(muR2/muF2);
  double lr = log(muR2/Q2);
	double TROLL_g02 =
				pow(alphas_muR,2)*(((2*(b1/b0*b0*M_PI*M_PI)+3*(b0*M_PI)*(2*b0*M_PI*lf+2.*CA * ZETA2 + 11./2.))*lr + 3*(b0*M_PI)*(b0*M_PI)*lr*lr+(-2.*(CA)*(CA)*ZETA2 + b0*M_PI*(2*b0*M_PI+b0*M_PI))*lf*lf + (2*(CA*CA*(2./3.+3./4.*ZETA3) - (CA/6.+CF/8.)*nF)+2*b0*M_PI*(2.*CA * ZETA2 + 11./2.)+2*4.*CA*(CA)*ZETA3)*lf
				            +(837./16. + 67./2.*ZETA2 - 9./20.*ZETA2*ZETA2 - 165./4.*ZETA3 + nF * ( 5./6.*ZETA3 - 5./3.*ZETA2 - 247./36. ))
										+ 2.*CA * ZETA2*11./2. + 121./16. + 2777./144. + 19./8.*lHt + (-67./48. + 2./3.*lHt)*nF)/(M_PI*M_PI)
										+ (2*b0*M_PI*lf+2.*CA * ZETA2 + 11./2.+2*b0*M_PI*lr)/M_PI*(2*A1* (ZETA2+EulerGamma*EulerGamma - EulerGamma*lf))
										+ -4*A1*A1*((pow(EulerGamma,4) -6*EulerGamma*EulerGamma*ZETA2 +8*EulerGamma*ZETA3 +3*ZETA2*ZETA2 -6*ZETA4) -2*(2*EulerGamma-lf)*(pow(EulerGamma,3) -3*EulerGamma*ZETA2 +2*ZETA3) +pow(2*EulerGamma-lf,2)*(pow(EulerGamma,2) -ZETA2))/2. + 4*A1*b0/3.*(pow(EulerGamma,3) -3*EulerGamma*ZETA2 +2*ZETA3) -2*( A1*b0*(2*EulerGamma+lr) + A2 )*(pow(EulerGamma,2) -ZETA2)
										- EulerGamma* ( D2 - 4.*A2*EulerGamma - 4.*A1*b0*(EulerGamma*EulerGamma+ZETA2)
												- 4.*A1*b0*EulerGamma*lr +2*A2*lf
												- A1*b0*(lr*lr - lrf*lrf) ) + (2*A1* (ZETA2+EulerGamma*EulerGamma - EulerGamma*lf))*(2*A1* (ZETA2+EulerGamma*EulerGamma - EulerGamma*lf)));
	return TROLL_g02;
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

// without alphas expansion, see notes
complex<double> g1_alphas(double A1,complex<double> N){
	complex<double> N2 = pow(N,2);
	complex<double> alphasQ2 = falphasQ2(Q2);
	complex<double> alphasQ2_N2 = falphasQ2(Q2/N2);
	return A1/(2.*M_PI*pow(b0,2))*(1./alphasQ2-1./alphasQ2_N2+1./alphasQ2_N2*log(alphasQ2/alphasQ2_N2));
}


complex<double> h1_alphas(double A1,complex<double> N){
	//complex<double> N2 = pow(N,2);
	//complex<double> alphasQ2_N2 = falphasQ2(Q2/N2);
	//derivative from g1 //works (returns same as h1NLP), its 1/2*D[g1,N].
	complex<double> delta = 1.E-4;
	/*complex<double> g1p = 1./alphas_muR*g1(A1,alphas_muR*b0*log(N+delta));
	complex<double> g1m = 1./alphas_muR*g1(A1,alphas_muR*b0*log(N-delta));;
	return 1./2.*(g1p-g1m)/(2.*delta);*/
	//derivative from g1_alphas (returns same as h1NLP), its 1/N*beta(alphas)*D[g1_alphas,alphas].
	complex<double> N2 = pow(N,2);
	complex<double> alphasQ2 = falphasQ2(Q2);
	complex<double> alphasQ2_N2 = falphasQ2(Q2/N2);
	complex<double> alphasQ2_N2p = alphasQ2_N2*(1.+delta);
	complex<double> alphasQ2_N2m = alphasQ2_N2*(1.-delta);
	complex<double> g1p =  A1/(2.*M_PI*pow(b0,2))*(1./alphasQ2-1./alphasQ2_N2p+1./alphasQ2_N2p*log(alphasQ2/alphasQ2_N2p));
	complex<double> g1m =  A1/(2.*M_PI*pow(b0,2))*(1./alphasQ2-1./alphasQ2_N2m+1./alphasQ2_N2m*log(alphasQ2/alphasQ2_N2m));
	return 1./N*b0*alphasQ2_N2*alphasQ2_N2*(g1p-g1m)/(2.*alphasQ2_N2*delta);
}
