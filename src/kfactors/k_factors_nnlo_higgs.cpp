#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "polygamma.h"
#include "gsl/gsl_sf_zeta.h"
#include "k_factors_higgs.h"
#include "k_factors_nnlo_higgs.h"
using namespace std;

//////////////////////////////////////////////////////////
///
/// contains all K factors for higgs
/// split up in NNLO and LP, NLP, NNLP, NNNLP, full
///
/// extracted log dependent parts from ihixs
/// note that the power expansions do not contain log(muF/Q) deps
///
//////////////////////////////////////////////////////////

// log dep parts extracted from ihixs
// (log(z)) = log(z)
// (li2(-z)+log(z)*log(z+1.)) = li2(-z)+log(z)*log(1.+z)
// (-li2(z)-log(1.-z)*log(z)) = -li2(z)-log(z)*log(1.-z)
// (li2(z)) = li2(z)



///////////////////////////
/// gg channel
///////////////////////////
// the NNLO fuCAtions for the gg channel
// https://arxiv.org/pdf/hep-ph/0207004.pdf eqn. 47 and 48 (note that there is no scale dependeCAe now) also we multiply with 1/x and take into account the extra 1/x = (1-x)/x + 1 term from the plus distribution
// checked with mathematica
double higgs_NNLO_gg_reg(double x){
	return logdep_gg_reg(x)+pow(alphas_muR/M_PI,2)*(
	(1./x*((133.-90.*zeta2)*log(1.-x)+(-101./3.+33.*zeta2+351./2.*zeta3)-33.*pow(log(1.-x),2)+72.*pow(log(1.-x),3)) // this is the same
	+1./x*(9.*(38.*pow(x,2)-20.*pow(x,3)+18.*x-39.*pow(x,4)+14.+7.*pow(x,5))/(1.-pow(x,2))*Li3(x)-18.*pow((pow(x,2)+x+1),2)/(1.+x)*S12(pow(x,2))
			+9.*(4.*pow(x,4)+8.*pow(x,3)+21.*pow(x,2)+14.*x+7.)/(1.+x)*S12(-x)-9./2.*(5.*pow(x,5)-51.*pow(x,4)-57.*pow(x,3)+53.*pow(x,2)+59.*x-11.)/(1.-pow(x,2))*S12(x)
			-9./2.*(8.*pow(x,4)+8.*pow(x,3)-3.*pow(x,2)-2.*x-1.)/(1.+x)*Li3(-x)-9./2.*(16.+13.*pow(x,5)-40.*pow(x,3)-67.*pow(x,4)+64.*pow(x,2)+36.*x)/(1.-pow(x,2))*li2(x)*log(x)
			+9./2.*(2.*pow(x,4)-15.*pow(x,2)-10.*x-5.)/(1.+x)*li2(-x)*log(x)-9./4.*(59.+177.*pow(x,2)-116.*pow(x,3)+59.*pow(x,4)-118.*x)/(1.-x)*log(x)*pow(log(1.-x),2)
			+27.*(3.*pow(x,2)+2.*x+1)/(1.+x)*li2(-x)*log(1.+x)+9.*(6.-11.*pow(x,3)+18.*pow(x,2)-12.*x+6*pow(x,4))/(1.-x)*pow(log(x),2)*log(1.-x)
			+9./2.*(3.-8.*pow(x,3)+3.*pow(x,4)-6.*x+9.*pow(x,2))/(1.-x)*li2(x)*log(1.-x)-3./2.*(7.*x-7.*pow(x,3)+4.+18.*pow(x,2)-17.*pow(x,4)+9.*pow(x,5))/(1.-pow(x,2))*pow(log(x),3)
			+9./2.*(8.*pow(x,4)+16.*pow(x,3)+33.*pow(x,2)+22.*x+11.)/(1.+x)*zeta2*log(1.+x)-36.*pow(pow(x,2)+x+1.,2)/(1.+x)*li2(x)*log(1.+x)
			-9./4.*(4.*pow(x,4)+8.*pow(x,3)+27.*pow(x,2)+18.*x+9.)/(1.+x)*log(1.+x)*pow(log(x),2)+(-21.+63/2.*pow(x,2)-18.*x+33./2.*pow(x,3))*log(1.+x)*log(x)
			+27./2.*(3.*pow(x,2)+2.*x+1.)/(1.+x)*pow(log(1.+x),2)*log(x)-3./4.*(-280.*pow(x,3)+143.*pow(x,4)+394.*x-289+21.*pow(x,2))/(1.-x)*li2(x)
			+(-21.+63./2.*pow(x,2)-18.*x+33./2.*pow(x,3))*li2(-x)+(-2559./4.*pow(x,3)+1079./2.*pow(x,2)-2687./4.*x+2027./4.)*log(1.-x)
			-3./8.*(374.*pow(x,4)-389.*x+154.+699.*pow(x,2)-827.*pow(x,3))/(1.-x)*pow(log(x),2)+(330.*pow(x,3)-348.*pow(x,2)+381.*x-297.)*pow(log(1.-x),2)
			+3./4.*(-1180.*pow(x,3)+641.-1238.*x+1227.*pow(x,2)+605.*pow(x,4))/(1.-x)*log(x)*log(1.-x)-72.*(2.-x+pow(x,2))*x*pow(log(1.-x),3)
			-1./8.*(4318.*pow(x,4)-6955.*pow(x,3)+6447.*pow(x,2)-5611.*x+2333.)/(1.-x)*log(x)+3./4.*(495.*pow(x,4)-886.*pow(x,3)+564.*pow(x,2)-200.*x+16.)/(1.-x)*zeta2
			+9.*(6.*x+18.*pow(x,2)+2.+10.*pow(x,5)-6.*pow(x,3)-19.*pow(x,4))/(1.-pow(x,2))*zeta2*log(x)-9./2.*(-48.*pow(x,3)+23.*pow(x,4)-46.*x+3.+69.*pow(x,2))/(1.-x)*zeta2*log(1.-x)
			+9./2.*(-36.-15.*pow(x,4)-52.*x+19.*pow(x,2)+13.*pow(x,3)+33.*pow(x,5))/(1.-pow(x,2))*zeta3+7539./16.*pow(x,3)-24107./48.*pow(x,2)+22879./48.*x-18157./48.))

	+ nF*(1./x*(-10./3.*log(1.-x)+(14./9.-2.*zeta2)+2.*pow(log(1.-x),2)) // this piece same as mathematica
		+1./x*((31./6.*x+1./6.+65./12.*pow(x,2))*S12(x)+(-31./12.*pow(x,2)+1./6.-17./6.*x)*Li3(x)
				+(47./12.*pow(x,2)+25./6.*x-1./6.)*li2(x)*log(x)+(-1./12.*pow(x,2)+1./6.*x-1./6.)*zeta2*log(1.-x)-4.*x*(1.+x)*zeta2*log(x)
				+(-1./6.*x+1./6.+1./12.*pow(x,2))*li2(x)*log(1.-x)+(1./12.-1./12.*x+1./24.*pow(x,2))*log(1.-x)*log(x)*log((1.-x)/x)
				+5./9.*x*(1.+x)*pow(log(x),3)+(-17./6.*pow(x,2)-7./3.*x-1./3.)*zeta3+(-34./9.*pow(x,3)+2./3.*pow(x,2)-8./3.*x+16./9.)*(pow(log(1.-x),2)-zeta2)
				-2./9.*(21.*pow(x,2)+7.*x+25.*pow(x,4)+17.-61.*pow(x,3))/(1.-x)*log(x)*log(1.-x)+(785./54.*pow(x,3)-83./36.*pow(x,2)+49./18.*x-461./54.)*log(1.-x)
				+1./72.*(-351.*pow(x,3)+117.*pow(x,2)+68.+132.*pow(x,4)+52.*x)/(1.-x)*pow(log(x),2)+1./36.*(227.*pow(x,3)+68.+4.*pow(x,4)-302.*x+21.*pow(x,2))/(1.-x)*li2(1.-x)
				+1./216.*(333.*pow(x,2)+2384.*pow(x,4)-598.*x-3041.*pow(x,3)+1282.)/(1.-x)*log(x)-8887./648.*pow(x,3)+1267./432.*pow(x,2)-497./216.*x+12923./1296.)));
}

//Wilson Coefficient part
//HEFT::gg_n2lo_lzbar0_lz0_L_1(z)+HEFT::gg_n2lo_lzbar0_lz1_L_1(z)	//+HEFT::gg_n2lo_lzbar0_lz2_L_1(z)+HEFT::gg_n2lo_lzbar1_L_1(z)//+HEFT::gg_n2lo_lzbar2_L_1(z)
//HEFT::gg_n2lo_lzbar1_L_2(z)+HEFT::gg_n2lo_lzbar0_lz0_L_2(z)+HEFT::gg_n2lo_lzbar0_lz1_L_2(z)
double logdep_gg_reg(double z){
	double L = log(muF2/Q2);
	double L2 = log(muF2/Q2)*log(muF2/Q2); //checked
	double z2 = z*z;
	double z3 = z*z*z;
	double LZp = log(1.+z);
	double LZm = log(1.-z);
	double LZ = log(z); //checked
	double LZ2 = pow(log(z),2); //checked
	double LZm2 = pow(log(1.-z),2); //checked
	double result = L*(-((8./3.)*nF*z*li2(1. - z)) - (8./3.)*nF*li2(1. - z) - (36.*z2*li2(-z))/(z + 1.) - (18.*z3*li2(-z))/(z + 1.) + 144.*z*li2(1. - z) + 144.*li2(1. - z) - (54.*z*li2(-z))/(z + 1.) - (18.*li2(-z))/(z*(z + 1.)) - (36.*li2(-z))/(z + 1.)
					+ LZ*(LZm*(-((8.*nF*z)/3.) - (8.*nF)/3. + (396.*z2)/(z-1.) - (126.*z3)/(z-1.) - (378.*z)/(z-1.) + 108./(z-1.) - 126./((z-1.)*z))
								- (36.*LZp*z2)/(z + 1) - (18.*LZp*z3)/(z + 1.) - (54.*LZp*z)/(z + 1.) - (18.*LZp)/(z*(z + 1.)) - (36.*LZp)/(z + 1.)
								+ (3.*nF*z2)/(z-1.) - (nF*z3)/(z-1.) + 2.*nF*z - (3.*nF*z)/(z-1.) + nF/(z-1.) + (8.*nF)/(9.*z) - nF/((z-1.)*z) - (16.*nF*z2)/9. + (8.*nF)/3. - (1059.*z2)/(2.*(z-1.)) + (561.*z3)/(2.*(z-1.)) + (837.*z)/(2*(z-1.)) - 537./(2.*(z-1.)) + 231./(2.*(z-1.)*z))
					+ LZ2*((4.*nF*z)/3. + (4.*nF)/3. - 18./(z*(1.-z2)) - 72.*z + 27.*z2 - 27./(1.-z2)) + LZm*(-((4.*nF*z2)/(z-1.)) + (2.*nF*z3)/(z-1.) + (6.*nF*z)/(z-1.) + (4.*nF*z)/3. - (4.*nF)/(z-1.) - (34.*nF)/(9.*z) + (16.*nF*z2)/9. - (4.*nF)/3. + 348.*z + 330./z - 330.*z2 - 381.)
					+ LZm2*(-(108.*z) - 108./z + 108.*z2 + 216.) + (32.*nF*z2)/(3.*(z-1.)) - (19.*nF*z3)/(4.*(z-1.)) - (27.*nF*z)/(2.*(z-1.)) - (44.*nF*z)/9. + (32.*nF)/(3.*(z-1.)) - (37.*nF)/(12.*(z-1.)*z) + (77.*nF)/(27.*z) - (68.*nF*z2)/27. + (56.*nF)/9. + 36.*z*zeta2 + (45.*zeta2)/z - (9.*zeta2)/(z*(z + 1.)) - (1513.*z)/8. - 2161./(8.*z) - 54.*z2*zeta2 + (2161.*z2)/8. - 108.*zeta2 + 1781./8.) +
   			L2*(LZ*((2.*nF*z)/3. + (2.*nF)/3. - (72.*z2)/(z-1.) + (18.*z3)/(z-1.) + (54.*z)/(z-1.) + 18./((z-1.)*z))
				+ LZm*(36.*z + 36./z - 36.*z2 - 72.) + (nF*z2)/(z-1.) - (nF*z3)/(2.*(z-1.)) - (3.*nF*z)/(2.*(z-1.)) - (nF*z)/3. + nF/(z-1.) + (17.*nF)/(18.*z) - (4.*nF*z2)/9. + nF/3. - (249.*z)/4. - 297./(4.*z) + (297.*z2)/4. + 141./2.);
	return pow(alphas_muR/M_PI,2)*(11./2.*6.*(1./z-2.+z-pow(z,2))*log(Q2/muF2)+result);
}

// https://arxiv.org/pdf/hep-ph/0207004.pdf eqn. 47 and 48 (note that there is no scale dependeCAe now)
// same as in Mathematica code
double higgs_NNLO_gg_plus(double x){
	return logdep_gg_plus(x)+pow(alphas_muR/M_PI,2)*((133.-90.*zeta2)*log(1.-x)/(1.-x)+(-101./3.+33.*zeta2+351./2.*zeta3)*1./(1.-x)-33.*pow(log(1.-x),2)/(1.-x)+72.*pow(log(1.-x),3)/(1.-x)+nF*(-10./3.*log(1.-x)/(1.-x)+(14./9.-2.*zeta2)*1./(1.-x)+2.*pow(log(1.-x),2)/(1.-x)));
}
//https://arxiv.org/pdf/hep-ph/0302135.pdf also checked with this (appendix has scale dep terms)
//n_NNLO_D0_L() +n_NNLO_D1_L()+ n_NNLO_D2_L()
// n_NNLO_D0_L2()+ n_NNLO_D1_L2()
double logdep_gg_plus(double z){
	double L = log(muF2/Q2);
	double L2 = L*L;
	double D0 = 1./(1.-z);
	double D1 = log(1.-z)/(1.-z);
	double D2 = log(1.-z)*log(1.-z)/(1.-z);
	double result = L*(nF*((5*D0)/3.- 2*D1) + 45*D0*zeta2 - (67*D0)/2.+ 33*D1 - 108*D2 /*should this be there:*/-33*D0)
									+ L2*((D0*nF)/2.- (33*D0)/4.+ 36*D1);
	return pow(alphas_muR/M_PI,2)*(result);
}

// https://arxiv.org/pdf/hep-ph/0207004.pdf eqn. 47 and 48 (note that there is no scale dependeCAe now)
// Lt = log(muR2/mt2)
// same as in Mathematica code
double higgs_NNLO_gg_delta(){
	Lt = log(Q2/mt2);
	return logdep_gg_constant()+ pow(alphas_muR/M_PI,2)*(
		(11399./144.+133./2.*zeta2-165./4.*zeta3-9./20.*pow(zeta2,2)+19./8.*Lt)
		+nF*(-1189./144.+5./6.*zeta3-5./3.*zeta2+2./3.*Lt));
}

// from HiggsEFTGGfast::EvaluateDeltaAndPlusCoeffs
// checked - agreement with Leonardo
// also https://arxiv.org/pdf/hep-ph/0302135.pdf appendix B!!
// n_NNLO_delta_L()*_log_muf_mh_sq
//+ n_NNLO_delta_L2()*pow(_log_muf_mh_sq,2.)

double logdep_gg_constant(){
	double L = log(muF2/Q2);
	double L2 = L*L;
	double result = L*((nF*(-zeta2 - 11/6.) + (33*zeta2)/2.- (171*zeta3)/2.+ 27/2.))
									+ L2*(- 18*zeta2);
	return pow(alphas_muR/M_PI,2)*(result);
}

// power expansion, checked the coefficients with the mathematica code
double higgs_NNLO_gg_expansion(double x, int power){
	if(power==1){
		return pow(alphas_muR/M_PI,2)*(3147 - 684*pow(M_PI,2) + 2*nF*(-49 + 6*pow(M_PI,2)) + 3*log(1. - x)*(-1983 + 61*nF + 180*pow(M_PI,2) - 6*log(1. - x)*(-345 + 4*nF + 144*log(1. - x))) - 6318*zeta3)/36.;
	}
	if(power==2){
		return pow(alphas_muR/M_PI,2)*(1. - x)*((-55572 + 989*nF + 4734*pow(M_PI,2) - 64*nF*pow(M_PI,2) + 6*log(1 - x)*(12291 - 209*nF - 360*pow(M_PI,2) + log(1 - x)*(-5085 + 64*nF + 1728*log(1 - x))))/72. + 351*zeta3);
	}
	if(power==3){
		return pow(alphas_muR/M_PI,2)*(pow(-1. + x,2)*(44828 - 965*nF - 2376*pow(M_PI,2) + 24*log(1 - x)*(-1682 + 28*nF + 1023*log(1 - x))))/96.;
	}
	if(power==4){
		return pow(alphas_muR/M_PI,2)*pow(1. - x,3)*(-757.1067708333334 + 29*pow(M_PI,2) + nF*(11.765432098765432 - (5*pow(M_PI,2))/9.) + (798 - (28*nF)/3. - 15*pow(M_PI,2))*log(1 - x) + (5*(-855 + 16*nF)*pow(log(1 - x),2))/24. + 72*pow(log(1 - x),3) + (351*zeta3)/2.);
	}
	if(power==5){
		return pow(alphas_muR/M_PI,2)*pow(-1. + x,4)*(-327.39152083333335 + (391*pow(M_PI,2))/20. + nF*(5.709187885802469 - (5*pow(M_PI,2))/9.) + log(1 - x)*(485.4 - (272*nF)/45. - 15*pow(M_PI,2) + (-80.475 + (10*nF)/3.)*log(1 - x) + 72*pow(log(1 - x),2)) + (351*zeta3)/2.);
	}
	if(power==6){
		return pow(alphas_muR/M_PI,2)*pow(1. - x,5)*(-233.03654166666666 + (609*pow(M_PI,2))/40. + nF*(4.4048391203703705 - (76*pow(M_PI,2))/135.) + (413.56 - (3197*nF)/675. - 15*pow(M_PI,2))*log(1 - x) + ((-6129 + 608*nF)*pow(log(1 - x),2))/180. + 72*pow(log(1 - x),3) + (351*zeta3)/2.);
	}
	if(power==7){
		return pow(alphas_muR/M_PI,2)*pow(-1. + x,6)*(-187.7184452137998 + (3457*pow(M_PI,2))/280. + nF*(3.8538812200806247 - (77*pow(M_PI,2))/135.) + (log(1 - x)*(3671379 - 36731*nF - 141750*pow(M_PI,2) + 60*log(1 - x)*(-405 + 539*nF + 11340*log(1 - x))))/9450. + (351*zeta3)/2.);
	}
	if(power==8){
		return pow(alphas_muR/M_PI,2)*pow(1. - x,7)*(-157.21313551536687 + (1135*pow(M_PI,2))/112. + nF*(3.5441843393442034 - (109*pow(M_PI,2))/189.) + log(1 - x)*(379.50204081632654 - (85373*nF)/26460. - 15*pow(M_PI,2) + (21.776785714285715 + (218*nF)/63.)*log(1 - x) + 72*pow(log(1 - x),2)) + (351*zeta3)/2.);
	}
	if(power==9){
		return pow(alphas_muR/M_PI,2)*pow(-1. + x,8)*((-948743967014 + 23810115581*nF - 5174400*(-11421. + 800*nF)*pow(M_PI,2))/7.112448e9 - ((-9994419 + 70717*nF + 396900*pow(M_PI,2))*log(1 - x))/26460. + (41.87321428571428 + (220*nF)/63.)*pow(log(1 - x),2) + 72*pow(log(1 - x),3) + (351*zeta3)/2.);
	}
	if(power==10){
		return pow(alphas_muR/M_PI,2)*pow(1. - x,9)*(-113.44591817050895 + (809*pow(M_PI,2))/120. + nF*(3.2171579041421463 - (95*pow(M_PI,2))/162.) + log(1 - x)*(379.74150793650796 - (7451*nF)/3402. - 15*pow(M_PI,2) + log(1 - x)*(59.1 + (95*nF)/27. + 72*log(1 - x))) + (351*zeta3)/2.);
	}
	if(power==11){
		return pow(alphas_muR/M_PI,2)*pow(-1. + x,10)*(-96.11099293000008 + (24767*pow(M_PI,2))/4620. + nF*(3.130178372772995 - (239*pow(M_PI,2))/405.) + log(1 - x)*(383.98531746031745 - (823616*nF)/467775. - 15*pow(M_PI,2) + log(1 - x)*(74.2353896103896 + (478*nF)/135. + 72*log(1 - x))) + (351*zeta3)/2.);
	}
	if(power==12){
		return pow(alphas_muR/M_PI,2)*pow(1. - x,11)*(-80.69721017667689 + (5079*pow(M_PI,2))/1232. + nF*(3.0740456666125624 - (881*pow(M_PI,2))/1485.) + log(1 - x)*(389.61004820936637 - (9418411*nF)/6.8607e6 - 15*pow(M_PI,2) + log(1 - x)*(87.76607142857142 + (1762*nF)/495. + 72*log(1 - x))) + (351*zeta3)/2.);
	}
	if(power==13){
		return pow(alphas_muR/M_PI,2)*pow(-1. + x,12)*(-66.77320867484431 + (240059*pow(M_PI,2))/80080. + nF*(3.0407719398265622 - (59*pow(M_PI,2))/99.) + log(1 - x)*((28264271679 - 72689996*nF)/7.135128e7 - 15*pow(M_PI,2) + log(1 - x)*(100.01956793206793 + (118*nF)/33. + 72*log(1 - x))) + (351*zeta3)/2.);
	}
	if(power==14){
		return pow(alphas_muR/M_PI,2)*pow(1. - x,13)*(-54.044223123583215 + (14311*pow(M_PI,2))/7280. + nF*(3.024951566494933 - (70*pow(M_PI,2))/117.) + log(1 - x)*((595047198321 - 1022640010*nF)/1.4756742e9 - 15*pow(M_PI,2) + log(1 - x)*(111.22843406593407 + (140*nF)/39. + 72*log(1 - x))) + (351*zeta3)/2.);
	}
	if(power==15){
		return pow(alphas_muR/M_PI,2)*pow(-1. + x,14)*(-42.29690552559648 + (81009*pow(M_PI,2))/80080. + nF*(3.022742912785512 - (1475*pow(M_PI,2))/2457.) + log(1 - x)*(410.7389213991137 - (11548297*nF)/2.9513484e7 - 15*pow(M_PI,2) + log(1 - x)*(121.56464785214786 + (2950*nF)/819. + 72*log(1 - x))) + (351*zeta3)/2.);
	}
	if(power==16){
		return pow(alphas_muR/M_PI,2)*pow(1. - x,15)*(-31.370825678652583 + (1523*pow(M_PI,2))/12320. + nF*(3.031318854190447 - (569*pow(M_PI,2))/945.) - ((-21113379141 + 5566076*nF + 756756000*pow(M_PI,2))*log(1 - x))/5.04504e7 + (131.15949675324674 + (1138*nF)/315.)*pow(log(1 - x),2) + 72*pow(log(1 - x),3) + (351*zeta3)/2.);
	}
	if(power==17){
		return pow(alphas_muR/M_PI,2)*pow(-1. + x,16)*(-21.14204430455084 - (385087*pow(M_PI,2))/544544. + nF*(3.048543842443933 - (163*pow(M_PI,2))/270.) + log(1 - x)*((156737846709 + 56066504*nF)/3.675672e8 - 15*pow(M_PI,2) + log(1 - x)*(140.1157911941 + (163*nF)/45. + 72*log(1 - x))) + (351*zeta3)/2.);
	}
	if(power==18){
		return pow(alphas_muR/M_PI,2)*pow(1. - x,17)*(-11.512822272425945 - (506427*pow(M_PI,2))/340340. + nF*(3.0727715659285475 - (1111*pow(M_PI,2))/1836.) + log(1 - x)*((28503721680426 + 26208824915*nF)/6.56107452e10 - 15*pow(M_PI,2) + log(1 - x)*(148.51573867309162 + (1111*nF)/306. + 72*log(1 - x))) + (351*zeta3)/2.);
	}
	if(power==19){
		return pow(alphas_muR/M_PI,2)*pow(-1. + x,18)*(-2.4047783740564834 - (172635469*pow(M_PI,2))/7.759752e7 + nF*(3.102711853767016 - (835*pow(M_PI,2))/1377.) + log(1 - x)*((150442559584749 + 214951253020*nF)/3.399829524e11 - 15*pow(M_PI,2) + log(1 - x)*(156.42622628919068 + (1670*nF)/459. + 72*log(1 - x))) + (351*zeta3)/2.);
	}
	if(power==20){
		return pow(alphas_muR/M_PI,2)*pow(1. - x,19)*(6.245853262257745 - (5335601*pow(M_PI,2))/1.825824e6 + nF*(3.137339768041842 - (935*pow(M_PI,2))/1539.) + log(1 - x)*((1883304752426037 + 3562789557670*nF)/4.1797904148e12 - 15*pow(M_PI,2) + log(1 - x)*(163.90248402912877 + (1870*nF)/513. + 72*log(1 - x))) + (351*zeta3)/2.);
	}
	if(power==21){
		return pow(alphas_muR/M_PI,2)*pow(-1. + x,20)*(14.491609203985492 - (185443119*pow(M_PI,2))/5.173168e7 + nF*(3.175831629546658 - (1561*pow(M_PI,2))/2565.) + log(1 - x)*((638993812780395 + 1478484240002*nF)/1.3932634716e12 - 15*pow(M_PI,2) + log(1 - x)*(170.9906972091376 + (3122*nF)/855. + 72*log(1 - x))) + (351*zeta3)/2.);
	}
	return 0.;
}


///////////////////////////
/// qg channel
///////////////////////////
//NNLO coefficients  qg channel
// klopt met de code
double higgs_NNLO_qg_reg(double x){
	return  logdep_qg(x)+pow(alphas_muR/M_PI,2)*
		(1./x*((170./3.*x+338./9.+119./3.*pow(x,2))*Li3(x)+(4.*x+4.+2.*pow(x,2))*Li3(-x)+(16.+8.*pow(x,2)+16.*x)*S12(-x)
				+(-614./9.*x-269./9.*pow(x,2)-74./9.)*S12(x)+(-2.*pow(x,2)-4.-4.*x)*S12(pow(x,2))+(367./27.+367./54.*pow(x,2)-367./27.*x)*pow(log(1.-x),3)
				+((2.+pow(x,2)-2.*x)*log(1.-x)-(446./9.*x+214./9.+281./9.*pow(x,2))*log(x)-(8.+4.*pow(x,2)+8.*x)*log(1.+x))*li2(x)
				+(8.+8.*x+4.*pow(x,2))*log((1.+x)/x)*li2(-x)+(-115./9.*pow(x,2)+230./9.*x-230./9.)*log(x)*pow(log(1.-x),2)
				+(107./9+107./18.*pow(x,2)-107./9.*x)*pow(log(x),2)*log(1.-x)+(-145./54.*pow(x,2)-71./27.*x-2.)*pow(log(x),3)
				+(-3.*pow(x,2)-6.-6.*x)*log(1.+x)*pow(log(x),2)+(4.*x+4.+2.*pow(x,2))*pow(log(1.+x),2)*log(x)
				+(-4./27.*pow(x,3)-74./9.*x-11./9.*pow(x,2)-166./27.)*li2(-x)+(2605./54.-146./9.*x+74./27.*pow(x,3)-79./6.*pow(x,2))*li2(x)
				+(1139./18.*x+37./12.*pow(x,2)+8.*pow(x,3)-72.)*pow(log(1.-x),2)+(-121./18.*pow(x,2)-326./27.*pow(x,3)-826./9.*x+5935./54.)*log(x)*log(1.-x)
				+(113./27.*pow(x,3)+244./9.*x-13./3.*pow(x,2)-31./2.)*pow(log(x),2)+(-4./27.*pow(x,3)-74./9.*x-11./9.*pow(x,2)-166./27.)*log(1.+x)*log(x)
				+zeta2*(-59./9.*pow(x,2)+118./9.*x-118./9)*log(1.-x)+zeta2*(140./9.*x+128./9.*pow(x,2)+52./9.)*log(x)
				+zeta2*(12.+12.*x+6.*pow(x,2))*log(1.+x)+(-392./81.*pow(x,3)-49./3.*pow(x,2)+23671./162.-106.*x)*log(1.-x)
				+(1985./108.*pow(x,2)+800./9.*x-12209./162.+616./81.*pow(x,3))*log(x)+(-292./27.*pow(x,3)-82./3.*x+16./3.*pow(x,2)+221./27.)*zeta2
				+(-18.*x+10.*pow(x,2)+92./9.)*zeta3-210115./1944.+1537./486.*pow(x,3)+16465./162.*x+2393./648.*pow(x,2))

		+nF*1./x*((1./18.*pow(x,2)-1./9.*x+1./9.)*pow(log(1.-x),2)+(-38./27.*x+19./27.*pow(x,2)+29./27.)*log(x)-209./81.*x+265./162. //klopt
				+((-4./9.+4./9.*x-2./9.*pow(x,2))*log(x)-pow(x,2)+16./9.*x-13./9.)*log(1.-x)+179./162.*pow(x,2)+(1./9.*pow(x,2)-2./9.*x+2./9.)*pow(log(x),2))
		);

}


double logdep_qg(double z){
return pow(alphas_muR/M_PI,2)*(
	 log(Q2/muF2)*11.*(2. + (-2. + z)*z)/(3.*z) //Wilson coefficient part
	 +log(muF2/Q2)*(( -((-1 + CA)*(1 + CA)*(1728*pow(CA,2) + 1728*pow(CA,2)*z + 864*pow(CA,2)*pow(z,2))*(li2(-z)+log(z)*log(1.+z)))/(3456.*pow(CA,2)*z) - ((-1 + CA)*(1 + CA)*(-864 - 4320*pow(CA,2) - 15552*pow(CA,2)*z - 7776*pow(CA,2)*pow(z,2))*(-li2(z)-log(z)*log(1.-z)))/(3456.*pow(CA,2)*z) - ((-1 + CA)*(1 + CA)*(-636 + 29700*pow(CA,2) - 1392*CA*nF + 648*z - 21336*pow(CA,2)*z + 1824*CA*nF*z + 216*pow(z,2) - 4560*pow(CA,2)*pow(z,2) - 912*CA*nF*pow(z,2) + 744*pow(z,3) - 1128*pow(CA,2)*pow(z,3) + 864*zeta2 - 6048*pow(CA,2)*zeta2 - 1728*z*zeta2 - 12096*pow(CA,2)*z*zeta2 + 864*pow(z,2)*zeta2 - 8640*pow(CA,2)*pow(z,2)*zeta2 + 13392*pow(CA,2)*log(z) - 576*CA*nF*log(z) + 432*z*log(z) - 16560*pow(CA,2)*z*log(z) + 576*CA*nF*z*log(z) - 2052*pow(z,2)*log(z) + 2340*pow(CA,2)*pow(z,2)*log(z) - 288*CA*nF*pow(z,2)*log(z) - 288*pow(z,3)*log(z) - 3168*pow(CA,2)*pow(z,3)*log(z) + 2592*pow(CA,2)*pow(log(z),2) + 432*z*pow(log(z),2) + 3024*pow(CA,2)*z*pow(log(z),2) - 216*pow(z,2)*pow(log(z),2) + 3240*pow(CA,2)*pow(z,2)*pow(log(z),2)))/(3456.*pow(CA,2)*z))
	     + log(1.-z)*(((-1 + pow(CA,2))*(-162 + 2766*pow(CA,2) - 48*CA*nF + 288*z - 2496*pow(CA,2)*z + 48*CA*nF*z - 126*pow(z,2) + 6*pow(CA,2)*pow(z,2) - 24*CA*nF*pow(z,2) - 288*pow(CA,2)*pow(z,3)))/(288.*pow(CA,2)*z) + ((-1 + pow(CA,2))*(-72 + 1224*pow(CA,2) + 144*z + 432*pow(CA,2)*z - 72*pow(z,2) + 1080*pow(CA,2)*pow(z,2))*log(z))/(288.*pow(CA,2)*z))
			 + pow(log(1.-z),2)*((-3*(-1 + pow(CA,2))*(-1 + 7*pow(CA,2))*(2 - 2*z + pow(z,2)))/(16.*pow(CA,2)*z))
		 ) //HEFT::qg_n2lo_r_lz0_L1(z)+HEFT::qg_n2lo_r_lz1_L1(z)+HEFT::qg_n2lo_r_lz2_L1(z)
//double log agrees
		+ pow(log(muF2/Q2),2)*(-((-1. + CA)*(1. + CA)*(8280*pow(CA,2) - 288*CA*nF + 216*z - 6984*pow(CA,2)*z + 288*CA*nF*z - 54*pow(z,2) + 198*pow(CA,2)*pow(z,2) - 144*CA*nF*pow(z,2) - 864*pow(CA,2)*pow(z,3) + 2592*pow(CA,2)*log(z) + 216*z*log(z) + 2376*pow(CA,2)*z*log(z) - 108*pow(z,2)*log(z) + 2700*pow(CA,2)*pow(z,2)*log(z)))/(3456.*pow(CA,2)*z)
                          +log(1.-z)*((-1. + pow(CA,2))*(-36 + 252*pow(CA,2) + 36*z - 252*pow(CA,2)*z - 18*pow(z,2) + 126*pow(CA,2)*pow(z,2)))/(288.*pow(CA,2)*z))
 );
}

// checked
double higgs_NNLO_qg_expansion(double x, int power){
	if(power==1){
		return pow(alphas_muR/M_PI,2)*(132 + 52*nF + 261*pow(M_PI,2) + 3*log(1. - x)*(2046 - 72*nF - 100*pow(M_PI,2) + log(1. - x)*(255 + 6*nF + 734*log(1. - x))) + 5598*zeta3)/324.;
	}
	if(power==2){
		return pow(alphas_muR/M_PI,2)*((1 - x)*(-17130 + 52*nF + 693*pow(M_PI,2) + 3*log(1 - x)*(7242 - 24*nF - 100*pow(M_PI,2) + log(1 - x)*(-1527 + 6*nF + 734*log(1 - x))) + 5598*zeta3))/324.;
	}
	if(power==3){
		return pow(alphas_muR/M_PI,2)*(pow(-1. + x,2)*(-53805 + 772*nF + 2016*pow(M_PI,2) + 6*log(1 - x)*(-2*(-5451. + 60*nF + 100*pow(M_PI,2)) + log(1 - x)*(-2223 + 12*nF + 1468*log(1 - x))) + 22392*zeta3))/648.;
	}
	if(power==4){
		return pow(alphas_muR/M_PI,2)*(pow(1 - x,3)*(-31941. + 308*nF + 1696*pow(M_PI,2) + 6*log(1 - x)*(9018 - 88*nF - 200*pow(M_PI,2) + log(1 - x)*(-887 + 12*nF + 1468*log(1 - x))) + 22392*zeta3))/648.;
	}
	if(power==5){
		return pow(alphas_muR/M_PI,2)*(pow(-1. + x,4)*(-140757 + 1264*nF + 8984*pow(M_PI,2) + 24*log(1 - x)*(15237 - 140*nF - 400*pow(M_PI,2) + log(1 - x)*(515 + 24*nF + 2936*log(1 - x))) + 179136*zeta3))/5184.;
	}
	if(power==6){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,5)*(-16.547468621399176 + (5*nF)/36. + (3527*pow(M_PI,2))/3240. + (log(1 - x)*(1087583 - 8580*nF - 30000*pow(M_PI,2) + 75*log(1 - x)*(2163 + 24*nF + 2936*log(1 - x))))/16200. + (311*zeta3)/9.);
	}
	if(power==7){
		return pow(alphas_muR/M_PI,2)*pow(-1. + x,6)*(-8.753997942386832 + (71*nF)/810. + (101*pow(M_PI,2))/180. + (log(1 - x)*(53977 - 354*nF - 1500*pow(M_PI,2) + 6*log(1 - x)*(2176 + 15*nF + 1835*log(1 - x))))/810. + (311*zeta3)/9.);
	}
	if(power==8){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,7)*(-2.3669754826092695 + (73*nF)/1134. + (409*pow(M_PI,2))/3780. + (log(1 - x)*(13392887 - 71610*nF - 367500*pow(M_PI,2) + 210*log(1 - x)*(20086 + 105*nF + 12845*log(1 - x))))/198450. + (311*zeta3)/9.);
	}
	if(power==9){
		return pow(alphas_muR/M_PI,2)*pow(-1. + x,8)*(3.1236286366288137 + (2603*nF)/45360. - (4399*pow(M_PI,2))/15120. + (log(1 - x)*(219331957 - 939960*nF - 5880000*pow(M_PI,2) + 210*log(1 - x)*(388701. + 1680*nF + 205520*log(1 - x))))/3.1752e6 + (311*zeta3)/9.);
	}
	if(power==10){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,9)*(7.980193194991841 + (1649*nF)/27216. - (88271*pow(M_PI,2))/136080. + (log(1 - x)*(225779997 - 760760*nF - 5880000*pow(M_PI,2) + 70*log(1 - x)*(1344743 + 5040*nF + 616560*log(1 - x))))/3.1752e6 + (311*zeta3)/9.);
	}
	if(power==11){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,10)*((32969616809 + 187979680*nF - 2595706400*pow(M_PI,2))/2.667168e9 + (log(1 - x)*(116548839 - 301000*nF - 2940000*pow(M_PI,2) + 70*log(1 - x)*(752539 + 2520*nF + 308280*log(1 - x))))/1.5876e6 + (311*zeta3)/9.);
	}
	if(power==12){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,11)*(16.371439691932636 + (317707*nF)/3.7422e6 - (950819*pow(M_PI,2))/748440. + (75.88826037176085 - (1805*nF)/12474. - (50*pow(M_PI,2))/27.)*log(1 - x) + (36.38926166426167 + nF/9.)*pow(log(1 - x),2) + (367*pow(log(1 - x),3))/27. + (311*zeta3)/9.);
	}
	if(power==13){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,12)*(20.08352419129137 + (255601*nF)/2.4948e6 - (289013*pow(M_PI,2))/187110. + (log(1 - x)*(27133493609 - 35947296*nF - 640332000*pow(M_PI,2) + 8316*log(1 - x)*(1635261 + 4620*nF + 565180*log(1 - x))))/3.4577928e8 + (311*zeta3)/9.);
	}
	if(power==14){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,13)*(23.55038318418005 + (3963493*nF)/3.24324e7 - (4376459*pow(M_PI,2))/2.43243e6 + (81.1151583287637 - (27026*nF)/405405. - (50*pow(M_PI,2))/27.)*log(1 - x) + (42.03979908979909 + nF/9.)*pow(log(1 - x),2) + (367*pow(log(1 - x),3))/27. + (311*zeta3)/9.);
	}
	if(power==15){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,14)*(26.811561774613985 + (24435937*nF)/1.702701e8 - (19817981*pow(M_PI,2))/9.72972e6 + (83.79221021799508 - (26167*nF)/810810. - (50*pow(M_PI,2))/27.)*log(1 - x) + (44.55755448255448 + nF/9.)*pow(log(1 - x),2) + (367*pow(log(1 - x),3))/27. + (311*zeta3)/9.);
	}
	if(power==16){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,15)*(29.897415938273976 + (3138797*nF)/1.89189e7 - (21985949*pow(M_PI,2))/9.72972e6 + (86.4809662338329 - (59*nF)/162162. - (50*pow(M_PI,2))/27.)*log(1 - x) + (46.90757483257483 + nF/9.)*pow(log(1 - x),2) + (367*pow(log(1 - x),3))/27. + (311*zeta3)/9.);
	}
	if(power==17){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,16)*(32.83175519637401 + (2119591*nF)/1.12112e7 - (96107237*pow(M_PI,2))/3.891888e7 + (89.16682840462067 + (19069*nF)/648648. - (50*pow(M_PI,2))/27.)*log(1 - x) + (49.111047054797055 + nF/9.)*pow(log(1 - x),2) + (367*pow(log(1 - x),3))/27. + (311*zeta3)/9.);
	}
	if(power==18){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,17)*(35.633603799843875 + (9851981053*nF)/4.63134672e10 - (103818941*pow(M_PI,2))/3.891888e7 + (91.83945709710807 + (3158401*nF)/5.513508e7 - (50*pow(M_PI,2))/27.)*log(1 - x) + (51.18533899379488 + nF/9.)*pow(log(1 - x),2) + (367*pow(log(1 - x),3))/27. + (311*zeta3)/9.);
	}
	if(power==19){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,18)*(38.31841006646755 + (16444736207*nF)/6.94702008e10 - (104952839*pow(M_PI,2))/3.675672e7 + (94.49151733334287 + (2302423*nF)/2.756754e7 - (50*pow(M_PI,2))/27.)*log(1 - x) + (53.1448751230614 + nF/9.)*pow(log(1 - x),2) + (367*pow(log(1 - x),3))/27. + (311*zeta3)/9.);
	}
	if(power==20){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,19)*(40.89889806916082 + (344370580973*nF)/1.3199338152e12 - (2118708421*pow(M_PI,2))/6.9837768e8 + (97.11782729925943 + (56718997*nF)/5.2378326e8 - (50*pow(M_PI,2))/27.)*log(1 - x) + (55.00177226788424 + nF/9.)*pow(log(1 - x),2) + (367*pow(log(1 - x),3))/27. + (311*zeta3)/9.);
	}
	if(power==21){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,20)*(43.385680767448534 + (188202536987*nF)/6.599669076e11 - (70654979*pow(M_PI,2))/2.2054032e7 + (99.7147677337292 + (23001757*nF)/1.7459442e8 - (50*pow(M_PI,2))/27.)*log(1 - x) + (56.76631092935273 + nF/9.)*pow(log(1 - x),2) + (367*pow(log(1 - x),3))/27. + (311*zeta3)/9.);
	}
	return 0.;
}

///////////////////////////
/// qq channel (stricktly identical quarks!)
///////////////////////////
//NNLO coefficients  qq channel
//klopt met code
double higgs_NNLO_qq_reg(double x){
	return logdep_qq(x)+pow(alphas_muR/M_PI,2)*
		(1./x*((368./27.*x+104./27.*pow(x,2)+400./27.)*Li3(x)-32./9.*pow(x+2.,2)*S12(x)-4./27.*(2.+pow(x,2)-2.*x)*pow(log(x),2)*log(1.-x)
				-4./27.*pow(x+2.,2)*pow(log(x),3)-16./27.*(19.+5.*pow(x,2)+17.*x)*li2(x)*log(x)-32./9.*(x+3.)*(1.-x)*pow(log(1.-x),2)
				+16./3.*(x+3.)*(1.-x)*log(x)*log(1.-x)+4./27.*(26.*x-18.+9.*pow(x,2))*pow(log(x),2)-8./9.*(-6.+pow(x,2)+4.*x)*li2(x)
				+4./3.*(5.*x+17.)*(1.-x)*log(1.-x)+(8./9.*pow(x+2.,2)*zeta2-118./9.+248./27.*x+46./9.*pow(x,2))*log(x)
				+(-8./27.*pow(x,2)-16./27.+16./27.*x)*zeta3+(16./3.-32./9.*x-8./3.*pow(x,2))*zeta2-4./27.*(27.*x+160.)*(1.-x))
		+nF*1./x*(0.)
		);

}
// checked
double higgs_NNLO_qq_expansion(double x, int power){

if(power==1){
	return 0;
}
if(power==2){
	return pow(alphas_muR/M_PI,2)*(1 - x)*(-0.7021050077301801 - 1.7777777777777777*log(1 - x) + 1.7777777777777777*pow(log(1 - x),2));
}
if(power==3){
	return pow(alphas_muR/M_PI,2)*pow(-1 + x,2)*(0.20610174766398903 - 2.6666666666666665*log(1 - x) + 2.6666666666666665*pow(log(1 - x),2));
}
if(power==4){
	return pow(alphas_muR/M_PI,2)*pow(1 - x,3)*(0.9607868633905992 - 5.037037037037036*log(1 - x) + 4.444444444444445*pow(log(1 - x),2));
}
if(power==5){
	return pow(alphas_muR/M_PI,2)*pow(-1 + x,4)*(0.6842081075929632 - 5.555555555555555*log(1 - x) + 5.777777777777778*pow(log(1 - x),2));
}
if(power==6){
	return pow(alphas_muR/M_PI,2)*pow(1 - x,5)*(0.4842658825325742 - 5.571358024691358*log(1 - x) + 6.903703703703704*pow(log(1 - x),2));
}
if(power==7){
	return pow(alphas_muR/M_PI,2)*pow(-1 + x,6)*(0.35934818315065997 - 5.310617283950617*log(1 - x) + 7.881481481481481*pow(log(1 - x),2));
}
if(power==8){
	return pow(alphas_muR/M_PI,2)*pow(1 - x,7)*(0.30279414276688704 - 4.882620307382211*log(1 - x) + 8.744973544973545*pow(log(1 - x),2));
}
if(power==9){
	return pow(alphas_muR/M_PI,2)*pow(-1 + x,8)*(0.3064180285392837 - 4.350949861426052*log(1 - x) + 9.517460317460317*pow(log(1 - x),2));
}
if(power==10){
	return pow(alphas_muR/M_PI,2)*pow(1 - x,9)*(0.3615157940336111 - 3.7547291509196272*log(1 - x) + 10.215873015873017*pow(log(1 - x),2));
}
if(power==11){
	return pow(alphas_muR/M_PI,2)*pow(-1 + x,10)*(0.4599514555519086 - 3.118867892836147*log(1 - x) + 10.852910052910053*pow(log(1 - x),2));
}
if(power==12){
	return pow(alphas_muR/M_PI,2)*pow(1 - x,11)*(0.5946094480978931 - 2.459645882772386*log(1 - x) + 11.438319704986371*pow(log(1 - x),2));
}
if(power==13){
	return pow(alphas_muR/M_PI,2)*pow(-1 + x,12)*(0.7594646209861272 - 1.7879228807463536*log(1 - x) + 11.979733846400514*pow(log(1 - x),2));
}
if(power==14){
	return pow(alphas_muR/M_PI,2)*pow(1 - x,13)*(0.9494906332033051 - 1.1110560542291774*log(1 - x) + 12.483230349897015*pow(log(1 - x),2));
}
if(power==15){
	return pow(alphas_muR/M_PI,2)*pow(-1 + x,14)*(1.1605186847493472 - 0.4340831063160585*log(1 - x) + 12.953722820389485*pow(log(1 - x),2));
}
if(power==16){
	return pow(alphas_muR/M_PI,2)*pow(1 - x,15)*(1.3890944657024042 + 0.23952666472960993*log(1 - x) + 13.395236861903529*pow(log(1 - x),2));
}
if(power==17){
	return pow(alphas_muR/M_PI,2)*pow(-1 + x,16)*(1.632350839894967 + 0.907383499253111*log(1 - x) + 13.811109877776543*pow(log(1 - x),2));
}
if(power==18){
	return pow(alphas_muR/M_PI,2)*pow(1 - x,17)*(1.8879006372793563 + 1.5678518673723827*log(1 - x) + 14.204138200216631*pow(log(1 - x),2));
}
if(power==19){
	return pow(alphas_muR/M_PI,2)*pow(-1 + x,18)*(2.153748715447829 + 2.219829741264905*log(1 - x) + 14.576687219824475*pow(log(1 - x),2));
}
if(power==20){
	return pow(alphas_muR/M_PI,2)*pow(1 - x,19)*(2.42822070681569 + 2.8625962364016266*log(1 - x) + 14.930775053788471*pow(log(1 - x),2));
}
if(power==21){
	return pow(alphas_muR/M_PI,2)*pow(-1 + x,20)*(2.7099055665905247 + 3.4957048454251773*log(1 - x) + 15.268136977115306*pow(log(1 - x),2));
}
return 0.;
}

double logdep_qq(double z){
	 return pow(alphas_muR/M_PI,2)*(
		 log(muF2/Q2)*(
			 (pow(-1 + CA,2)*pow(1 + CA,2)*(96*CA + 96*CA*z + 24*CA*pow(z,2))*(li2(z)))/(96.*pow(CA,3)*z) + (pow(-1 + CA,2)*pow(1 + CA,2)*(192*CA + 192*CA*z + 48*CA*pow(z,2))*(-li2(z)-log(1.-z)*log(z)))/(96.*pow(CA,3)*z) + (pow(-1 + CA,2)*pow(1 + CA,2)*(-153*CA + 108*CA*z + 45*CA*pow(z,2) + 96*CA*zeta2 + 96*CA*z*zeta2 + 24*CA*pow(z,2)*zeta2 - 72*CA*log(z) + 48*CA*z*log(z) + 30*CA*pow(z,2)*log(z) - 24*CA*pow(log(z),2) - 24*CA*z*pow(log(z),2) - 6*CA*pow(z,2)*pow(log(z),2)))/(96.*pow(CA,3)*z)
       + log(1.-z)*((pow(-1 + pow(CA,2),2)*(72 - 48*z - 24*pow(z,2)))/(48.*pow(CA,2)*z) + (pow(-1 + pow(CA,2),2)*(96 + 96*z + 24*pow(z,2))*log(z))/(48.*pow(CA,2)*z)
        )
		 )//HEFT::qq_n2lo_lz0_L1(z) +HEFT::qq_n2lo_lz1_L1(z)
     + pow(log(muF2/Q2),2)*((pow(-1 + CA,2)*pow(1 + CA,2)*(-36*CA + 24*CA*z + 12*CA*pow(z,2) - 24*CA*log(z) - 24*CA*z*log(z) - 6*CA*pow(z,2)*log(z)))/(96.*pow(CA,3)*z)
		 )//HEFT::qq_n2lo_lz0_L2(z)
	 );
 }

///////////////////////////
/// qq' channel
///////////////////////////
//NNLO coefficients  qq' channel
// klopt met code
double higgs_NNLO_qqp_reg(double x){
	return logdep_qqp(x)+pow(alphas_muR/M_PI,2)*
		(1./x*(32./9.*pow(x+2.,2)*(Li3(x)-S12(x))-8./3.*pow(x+2.,2)*log(x)*li2(x)-4./27.*pow(x+2.,2)*pow(log(x),3)
				-8./9.*(4.*x-6.+pow(x,2))*li2(x)-32./9.*(x+3.)*(1.-x)*pow(log(1.-x),2)+16./3.*(x+3.)*(1.-x)*log(x)*log(1.-x)
				+8./9.*(pow(x,2)+4.*x-3.)*pow(log(x),2)+8./9.*zeta2*pow(x+2.,2)*log(x)+4./3.*(5.*x+17.)*(1.-x)*log(1.-x)
				+2./9.*(29.*pow(x,2)+44.*x-59.)*log(x)+(16./3.-32./9.*x-8./3.*pow(x,2))*zeta2-2./9.*(11.*x+105.)*(1.-x))
		+nF*1./x*(0.)
		);

}
// checked
double higgs_NNLO_qqp_expansion(double x, int power){
	if(power==1){
		return 0;
	}
	if(power==2){
		return pow(alphas_muR/M_PI,2)*(1 - x)*(-0.7021050077301801 - 1.7777777777777777*log(1 - x) + 1.7777777777777777*pow(log(1 - x),2));
	}
	if(power==3){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,2)*(0.2801758217380631 - 2.6666666666666665*log(1 - x) + 2.6666666666666665*pow(log(1 - x),2));
	}
	if(power==4){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,3)*(1.1830090856128215 - 5.037037037037037*log(1 - x) + 4.444444444444445*pow(log(1 - x),2));
	}
	if(power==5){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,4)*(1.0144550211732115 - 5.555555555555555*log(1 - x) + 5.777777777777778*pow(log(1 - x),2));
	}
	if(power==6){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,5)*(0.8984634133967717 - 5.571358024691358*log(1 - x) + 6.903703703703704*pow(log(1 - x),2));
	}
	if(power==7){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,6)*(0.8402809677871477 - 5.310617283950618*log(1 - x) + 7.881481481481482*pow(log(1 - x),2));
	}
	if(power==8){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,7)*(0.8378127592669851 - 4.882620307382211*log(1 - x) + 8.744973544973545*pow(log(1 - x),2));
	}
	if(power==9){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,8)*(0.886078277412496 - 4.350949861426052*log(1 - x) + 9.517460317460317*pow(log(1 - x),2));
	}
	if(power==10){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,9)*(0.9786311346828995 - 3.7547291509196272*log(1 - x) + 10.215873015873017*pow(log(1 - x),2));
	}
	if(power==11){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,10)*(1.1089461318860911 - 3.118867892836147*log(1 - x) + 10.852910052910053*pow(log(1 - x),2));
	}
	if(power==12){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,11)*(1.2710779217824189 - 2.459645882772386*log(1 - x) + 11.438319704986371*pow(log(1 - x),2));
	}
	if(power==13){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,12)*(1.459868031932748 - 1.7879228807463536*log(1 - x) + 11.979733846400514*pow(log(1 - x),2));
	}
	if(power==14){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,13)*(1.6709440503851178 - 1.1110560542291774*log(1 - x) + 12.483230349897015*pow(log(1 - x),2));
	}
	if(power==15){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,14)*(1.9006391812148966 - 0.43408310631605856*log(1 - x) + 12.953722820389487*pow(log(1 - x),2));
	}
	if(power==16){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,15)*(2.1458906144273784 + 0.2395266647296099*log(1 - x) + 13.395236861903527*pow(log(1 - x),2));
	}
	if(power==17){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,16)*(2.4041409283737765 + 0.907383499253111*log(1 - x) + 13.811109877776543*pow(log(1 - x),2));
	}
	if(power==18){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,17)*(2.6732511540869317 + 1.5678518673723827*log(1 - x) + 14.204138200216631*pow(log(1 - x),2));
	}
	if(power==19){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,18)*(2.951427393939112 + 2.219829741264905*log(1 - x) + 14.576687219824473*pow(log(1 - x),2));
	}
	if(power==20){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,19)*(3.2371601917823867 + 2.8625962364016266*log(1 - x) + 14.930775053788471*pow(log(1 - x),2));
	}
	if(power==21){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,20)*(3.5291749497289544 + 3.4957048454251773*log(1 - x) + 15.268136977115306*pow(log(1 - x),2));
	}
	return 0.;
}

double logdep_qqp(double z){
	return pow(alphas_muR/M_PI,2)*(
		    log(muF2/Q2)*(-(pow(-1 + CA,2)*pow(1 + CA,2)*(-96 - 96*z - 24*pow(z,2))*li2(z))/(96.*pow(CA,2)*z) - (pow(-1 + CA,2)*pow(1 + CA,2)*(-192 - 192*z - 48*pow(z,2))*(-li2(z)-log(1.-z)*log(z)))/(96.*pow(CA,2)*z) - (pow(-1 + CA,2)*pow(1 + CA,2)*(153 - 108*z - 45*pow(z,2) - 96*zeta2 - 96*z*zeta2 - 24*pow(z,2)*zeta2 + 72*log(z) - 48*z*log(z) - 30*pow(z,2)*log(z) + 24*pow(log(z),2) + 24*z*pow(log(z),2) + 6*pow(z,2)*pow(log(z),2)))/(96.*pow(CA,2)*z)
        + log(1.-z)*((pow(-1 + pow(CA,2),2)*(72 - 48*z - 24*pow(z,2)))/(48.*pow(CA,2)*z) + (pow(-1 + pow(CA,2),2)*(96 + 96*z + 24*pow(z,2))*log(z))/(48.*pow(CA,2)*z)
        )
				) //HEFT::q1q2_n2lo_lz0_L1(z)+HEFT::q1q2_n2lo_lz1_L1(z)
	      + pow(log(muF2/Q2),2)*(-(pow(-1 + CA,2)*pow(1 + CA,2)*(36 - 24*z - 12*pow(z,2) + 24*log(z) + 24*z*log(z) + 6*pow(z,2)*log(z)))/(96.*pow(CA,2)*z)
        ) //	HEFT::q1q2_n2lo_lz0_L2(z)
 );
}
///////////////////////////
/// qqbar channel
///////////////////////////
//NNLO coefficients  qqbar channel
//klopt met code
double higgs_NNLO_qqbar_reg(double x){
	return logdep_qqbar(x)+pow(alphas_muR/M_PI,2)*
		(1./x*((-16./9.-16./9.*x-8./9.*pow(x,2))*Li3(-x)+(-16./27*pow(x,2)-32./27.-32./27.*x)*S12(-x)
				+32./9.*pow(x+2.,2)*Li3(x)-32./9.*pow(x+2.,2)*S12(x)-4./27.*pow(x+2.,2)*pow(log(x),3)+4./9.*(2.+2.*x+pow(x,2))*log(1.+x)*pow(log(x),2)
				+(-8./27.*(2.+2.*x+pow(x,2))*pow(log(1.+x),2)-8./3.*pow(x+2.,2)*li2(x)+8./9.*(2.+2.*x+pow(x,2))*li2(-x))*log(x)
				-16./27.*(2.+2.*x+pow(x,2))*li2(-x)*log(1.+x)+32./81.*(1.-x)*(13.*pow(x,2)-35.*x-14.)*pow(log(1.-x),2)
				+-16./81.*(1.-x)*(37.*pow(x,2)-101.*x-44.)*log(x)*log(1.-x)-8./81.*(44.*pow(x,3)+39.*x-81.*pow(x,2)+27.)*pow(log(x),2)
				+16./27.*x*(x+6.*pow(x,2)+2.)*log(1.+x)*log(x)+8./81.*(42.*x-87.*pow(x,2)+12.+10.*pow(x,3))*li2(x)
				+16./27.*x*(x+6.*pow(x,2)+2.)*li2(-x)-4./81.*(1.-x)*(384.*pow(x,2)-967.*x-75.)*log(1.-x)+(-16./27.*pow(x,2)-32./27.-32./27.*x)*zeta3
				+(8./9.*pow(x+2.,2)*zeta2+4222./81*pow(x,2)-2896./81.*x-512./27.*pow(x,3)-10./3.)*log(x)-8./27.*(2.+2.*x+pow(x,2))*zeta2*log(1.+x)
				+(752./81.*pow(x,3)-544./27.*pow(x,2)+80./81.+400./27.*x)*zeta2+4./81.*(1.-x)*(783.*pow(x,2)-1925.*x+373.))
		+nF*1./x*(32./81.*pow(1.-x,3)*log(1.-x)+(-64./27.*pow(x,2)+64./81.*pow(x,3)-16./27.+80./27.*x)*log(x)-8./243.*(1.-x)*(41.*pow(x,2)-88.*x+23.))
		);

}

double logdep_qqbar(double z){
return pow(alphas_muR/M_PI,2)*(
	   log(muF2/Q2)*(-(pow(-1 + CA,2)*pow(1 + CA,2)*(-864*CA - 864*CA*z - 216*CA*pow(z,2))*(-li2(z)-log(1.-z)*log(z)))/(864.*pow(CA,3)*z) - (pow(-1 + CA,2)*pow(1 + CA,2)*(1377*CA - 396*pow(CA,2) + 72*CA*nF - 144*z - 972*CA*z + 1332*pow(CA,2)*z - 216*CA*nF*z + 144*pow(z,2) - 405*CA*pow(z,2) - 1332*pow(CA,2)*pow(z,2) + 216*CA*nF*pow(z,2) + 396*pow(CA,2)*pow(z,3) - 72*CA*nF*pow(z,3) - 864*CA*zeta2 - 864*CA*z*zeta2 - 216*CA*pow(z,2)*zeta2 + 648*CA*log(z) - 216*z*log(z) - 432*CA*z*log(z) + 216*pow(CA,2)*z*log(z) + 216*pow(z,2)*log(z) - 270*CA*pow(z,2)*log(z) - 216*pow(CA,2)*pow(z,2)*log(z) - 144*pow(z,3)*log(z) + 144*pow(CA,2)*pow(z,3)*log(z) + 216*CA*pow(log(z),2) + 216*CA*z*pow(log(z),2) + 54*CA*pow(z,2)*pow(log(z),2)))/(864.*pow(CA,3)*z)
        +log(1.-z)*(-(pow(-1 + pow(CA,2),2)*(-24 - 216*CA + 24*pow(CA,2) + 72*z + 144*CA*z - 72*pow(CA,2)*z - 72*pow(z,2) + 72*CA*pow(z,2) + 72*pow(CA,2)*pow(z,2) + 24*pow(z,3) - 24*pow(CA,2)*pow(z,3)))/(144.*pow(CA,3)*z) - (pow(-1 + pow(CA,2),2)*(-144*CA - 144*CA*z - 36*CA*pow(z,2))*(log(z)))/(144.*pow(CA,3)*z))
			) //HEFT::qqb_nnlo_r_lz0_L1(z) +HEFT::qqb_nnlo_r_lz1_L1(z)
	   +pow(log(muF2/Q2),2)*(-(pow(-1 + CA,2)*pow(1 + CA,2)*(324*CA - 216*CA*z - 108*CA*pow(z,2) + 216*CA*log(z) + 216*CA*z*log(z) + 54*CA*pow(z,2)*log(z)))/(864.*pow(CA,3)*z)
	 ) // HEFT::qqb_nnlo_r_lz0_L2(z)
);
}

// checked
double higgs_NNLO_qqbar_expansion(double x, int power){
	if(power==1){
		return 0;
	}
	if(power==2){
		return pow(alphas_muR/M_PI,2)*((1 - x)*(60 - 8*pow(M_PI,2) + 48*(-1. + log(1 - x))*log(1 - x)))/27.;
	}
	if(power==3){
		return pow(alphas_muR/M_PI,2)*(-2*pow(-1. + x,2)*(-83 + 6*pow(M_PI,2) + 36*log(1 - x) - 36*pow(log(1 - x),2)))/27.;
	}
	if(power==4){
		return pow(alphas_muR/M_PI,2)*(4*pow(1 - x,3)*(10972 - 240*nF - 393*pow(M_PI,2) + 6*log(1 - x)*(-833 + 12*nF + 291*log(1 - x))))/729.;
	}
	if(power==5){
		return pow(alphas_muR/M_PI,2)*(pow(-1. + x,4)*(33835 - 384*nF - 1734*pow(M_PI,2) + 6*log(1 - x)*(-2819 + 48*nF + 1326*log(1 - x))))/729.;
	}
	if(power==6){
		return pow(alphas_muR/M_PI,2)*(pow(1 - x,5)*(7689391 - 25800*nF - 467700*pow(M_PI,2) + 60*log(1 - x)*(-64043 + 1200*nF + 36570*log(1 - x))))/182250.;
	}
	if(power==7){
		return pow(alphas_muR/M_PI,2)*(pow(-1. + x,6)*(4890253 + 13200*nF - 331600*pow(M_PI,2) + 40*log(1 - x)*(-59171. + 1200*nF + 39540*log(1 - x))))/121500.;
	}
	if(power==8){
		return pow(alphas_muR/M_PI,2)*(pow(1 - x,7)*(1634627677 + 12171600*nF - 119736400*pow(M_PI,2) + 280*log(1 - x)*(-2691197 + 58800*nF + 2065980*log(1 - x))))/4.16745e7;
	}
	if(power==9){
		return pow(alphas_muR/M_PI,2)*(pow(-1. + x,8)*(6445572073 + 72676800*nF - 500407600*pow(M_PI,2) + 560*log(1 - x)*(-4999009 + 117600*nF + 4361910*log(1 - x))))/1.66698e8;
	}
	if(power==10){
		return pow(alphas_muR/M_PI,2)*(pow(1 - x,9)*(172752517031. + 2493237600*nF - 14034913200*pow(M_PI,2) + 5040*log(1 - x)*(-13912987 + 352800*nF + 13709430*log(1 - x))))/4.500846e9;
	}
	if(power==11){
		return pow(alphas_muR/M_PI,2)*(pow(-1 + x,10)*(172261391123 + 2941999200*nF - 14512780800*pow(M_PI,2) + 5040*log(1 - x)*(-12879787 + 352800*nF + 14278320*log(1 - x))))/4.500846e9;
	}
	if(power==12){
		return pow(alphas_muR/M_PI,2)*(pow(1 - x,11)*(5096385094711 + 98492548000*nF - 442244584320*pow(M_PI,2) + 3696*log(1 - x)*(-479527319 + 14229600*nF + 596977920*log(1 - x))))/1.331250228e11;
	}
	if(power==13){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,12)*(38.37648393741604 + (15262*nF)/18711. - (319234*pow(M_PI,2))/93555. + (-12.24571419779693 + (32*nF)/81.)*log(1 - x) + (177916*pow(log(1 - x),2))/10395.);
	}
	if(power==14){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,13)*(38.531952179726524 + (1074266*nF)/1.216215e6 - (4252102*pow(M_PI,2))/1.216215e6 + (-11.2157782182106 + (32*nF)/81.)*log(1 - x) + (2380948*pow(log(1 - x),2))/135135.);
	}
	if(power==15){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,14)*(38.733713151627256 + (1148426*nF)/1.216215e6 - (4347472*pow(M_PI,2))/1.216215e6 + (-10.220308951801163 + (32*nF)/81.)*log(1 - x) + (2444528*pow(log(1 - x),2))/135135.);
	}
	if(power==16){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,15)*(38.970858093761805 + (1215922*nF)/1.216215e6 - (4436968*pow(M_PI,2))/1.216215e6 + (-9.256697824087471 + (32*nF)/81.)*log(1 - x) + (2504192*pow(log(1 - x),2))/135135.);
	}
	if(power==17){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,16)*(39.23542402660024 + (255566*nF)/243243. - (9042533*pow(M_PI,2))/2.43243e6 + (-8.322716990106638 + (32*nF)/81.)*log(1 - x) + (2560391*pow(log(1 - x),2))/135135.);
	}
	if(power==18){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,17)*(39.52145679654574 + (4538959*nF)/4.135131e6 - (156431767*pow(M_PI,2))/4.135131e7 + (-7.41641399183901 + (32*nF)/81.)*log(1 - x) + (44429549*pow(log(1 - x),2))/2.297295e6);
	}
	if(power==19){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,18)*(39.82441463982607 + (70791089*nF)/6.2026965e7 - (79499666*pow(M_PI,2))/2.0675655e7 + (-6.536049125289343 + (32*nF)/81.)*log(1 - x) + (45285404*pow(log(1 - x),2))/2.297295e6);
	}
	if(power==20){
		return pow(alphas_muR/M_PI,2)*pow(1 - x,19)*(40.140774730542404 + (1393018631*nF)/1.178512335e9 - (1533676814*pow(M_PI,2))/3.92837445e8 + (-5.680055298934234 + (32*nF)/81.)*log(1 - x) + (291959372*pow(log(1 - x),2))/1.4549535e7);
	}
	if(power==21){
		return pow(alphas_muR/M_PI,2)*pow(-1 + x,20)*(40.46776562277762 + (1437991559*nF)/1.178512335e9 - (311152976*pow(M_PI,2))/7.8567489e7 + (-4.847010393104761 + (32*nF)/81.)*log(1 - x) + (534362096*pow(log(1 - x),2))/2.6189163e7);
	}
	return 0.;
}
