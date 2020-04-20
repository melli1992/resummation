#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include "parameters.h"
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "k_factors_nnlo_dy.h"
#include "k_factors_dy.h"
#include "polygamma.h"
using namespace std;

//////////////////////////////////////////////////////////
///
/// contains all K factors for drell yan  at NNLO
/// split up in LP, NLP, NNLP
/// base is mathematica code
///
//////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////
///
/// qqbar channel
///
////////////////////////////////////////////////////////////
double DY_NNLO_qqbar_plus(double x){
	return 1./4.*pow(alphas_muR/M_PI,2)*beta0*log(muR2/muF2)*(2*CF*(2*log(1.-x)/(1.-x) + log(Q2/muF2)/(1.-x)))+(pow(alphas_muR,2)*CF*(9*pow(log(Q2/muF2),2)*(11*CA - 2*(18*CF + nF) - 48*CF*log(1 - x)) + 6*log(Q2/muF2)*(-67*CA + 144*CF + 10*nF + 18*CA*zeta2 + 72*CF*zeta2 + 6*(11*CA - 2*(9*CF + nF))*log(1 - x) - 216*CF*pow(log(1 - x),2)) + 2*(202*CA - 28*nF - 198*CA*zeta2 + 36*nF*zeta2 - 189*CA*zeta3 - 864*CF*zeta3 + 6*(CA*(-67 + 18*zeta2) + 2*(5*nF + 36*CF*(2 + zeta2)))*log(1 - x) + 18*(11*CA - 2*nF)*pow(log(1 - x),2) - 432*CF*pow(log(1 - x),3))))/(108.*pow(M_PI,2)*(-1 + x));
}
double DY_NNLO_qqbar_delta(){
	return pow(alphas_muR/M_PI,2)*beta0*log(muR2/muF2)*(CF*(-24 + 2*pow(M_PI,2) + 9*log(Q2/muF2)))/(4.*6.)+(pow(alphas_muR,2)*CF*(22995*CF + 3810*nF - 12600*CF*zeta2 - 2240*nF*zeta2 + 288*CF*pow(zeta2,2) - 10800*CF*zeta3 + 1440*nF*zeta3 + CA*(-23025 + 11840*zeta2 - 432*pow(zeta2,2) + 5040*zeta3) + 60*(-34*nF + CA*(193 - 72*zeta3) + 3*CF*(-93 + 24*zeta2 + 176*zeta3))*log(Q2/muF2) - 180*(11*CA - 2*(nF + CF*(9 - 16*zeta2)))*pow(log(Q2/muF2),2)))/(2880.*pow(M_PI,2));
}
double DY_NNLO_qqbar_NS(double x){
	return (pow(alphas_muR,2)*CF*(beta0*log(muR2/muF2)*(-4*(1 + x)*log(Q2/muF2) - 8*(1 + x)*log(1 - x) + (4*(1 + pow(x,2))*log(x))/(-1 + x)) + (2*nF*(94 - 206*x - 18*(1 + x)*pow(log(Q2/muF2),2) + 24*(-1 + 11*x)*log(1 - x) + 36*(2 - 3*x)*log(x) - (12*log(Q2/muF2)*(-1 + 12*x - 11*pow(x,2) + 6*(-1 + pow(x,2))*log(1 - x) - 6*(1 + pow(x,2))*log(x)))/(-1 + x) - (18*(1 + pow(x,2))*(log(x)*(5 - 8*log(1 - x) + 3*log(x)) - li2(1 - x)))/(-1 + x) + 3*(1 + x)*(4*pow(M_PI,2) - 24*pow(log(1 - x),2) + 3*pow(log(x),2) + 6*li2(1 - x))))/27. + 2*(-CA/2. + CF)*(94 - 78*x + 8*(-8 + 7*x)*log(1 - x) + 2*(22 - 9*x)*log(x) + (4*log(Q2/muF2)*(8 - 15*x + 7*pow(x,2) + (5 - 2*pow(x,2))*log(x) + (1 + pow(x,2))*pow(log(x),2) + 2*(1 + pow(x,2))*li2(1 - x)))/(-1 + x) + ((1 + x)*(log(x)*(-168*log(1 - x) + log(x)*(69 + 4*log(x))) + 12*(-13 + 2*log(x))*li2(1 - x) - 48*Li3(1 - x)))/6. - ((1 + pow(x,2))*(72*log(x) - 72*log(1 - x)*log(x) + 45*pow(log(x),2) - 156*log(1 - x)*pow(log(x),2) + 16*pow(log(x),3) - 12*(3 + 8*log(1 - x) + 6*log(x))*li2(1 - x) - 216*log(x)*li2(x) + 96*Li3(1 - x) + 216*Li3(x) - 216*zeta3))/(6.*(-1 + x))) + CF*(-72 + 48*x + 4*(64 + 3*x)*log(1 - x) - 64*(-1 + x)*(zeta2 - pow(log(1 - x),2)) - 16*log(x) + 8*(-4 + 13*x)*log(x) + 16*(7 - 6*x)*log(1 - x)*log(x) - 16*(3 + x)*log(1 - x)*log(x) + 16*(-2 + x)*pow(log(x),2) + 8*(3 + x)*pow(log(x),2) + pow(log(Q2/muF2),2)*(-8*(5 + x) + (16*(1 + pow(x,2))*log(x))/(-1 + x) + 8*(1 + x)*(-4*log(1 - x) + log(x))) + 8*(3 - 2*x)*li2(1 - x) - 16*(3 + x)*li2(1 - x) + (4*log(Q2/muF2)*(-30 + 26*x + 4*pow(x,2) - 8*zeta2 + 8*pow(x,2)*zeta2 - 24*(-1 + pow(x,2))*pow(log(1 - x),2) - 2*log(x) + 20*x*log(x) - 6*pow(x,2)*log(x) - 3*pow(log(x),2) - 9*pow(x,2)*pow(log(x),2) + 4*log(1 - x)*(7 - 8*x + pow(x,2) + (5 + 9*pow(x,2))*log(x)) + 4*(-3 + pow(x,2))*li2(1 - x)))/(-1 + x) - (2*(1 + x)*(192*zeta3 - 96*zeta2*log(1 - x) + 96*pow(log(1 - x),3) + 48*zeta2*log(x) - 48*pow(log(1 - x),2)*log(x) + 24*log(1 - x)*pow(log(x),2) - 7*pow(log(x),3) - 72*log(1 - x)*li2(1 - x) - 24*log(x)*li2(x) + 60*Li3(1 - x) + 24*Li3(x) - 24*zeta3))/3. + (4*(1 + pow(x,2))*(-14*log(x) - 16*zeta2*log(x) + 31*pow(log(1 - x),2)*log(x) - 14*log(1 - x)*pow(log(x),2) + 3*pow(log(x),3) - 6*(log(1 - x) - log(x))*li2(1 - x) + 8*log(x)*li2(x) + 2*Li3(1 - x) - 8*Li3(x) + 8*zeta3))/(-1 + x)) + CA*(-16.51851851851852 + (2278*x)/27. - (4*(19 + 25*x)*zeta2)/3. + (22*(1 + x)*pow(log(Q2/muF2),2))/3. - (4*(38 + 239*x)*log(1 - x))/9. + (2*(-26 + 57*x)*log(x))/3. + 4*(-3 + x)*log(1 - x)*log(x) - 16*x*log(1 - x)*log(x) + ((23 - 25*x)*pow(log(x),2))/6. + 8*x*pow(log(x),2) + log(Q2/muF2)*((-4*(19 + 124*x))/9. + (1 + x)*(8*zeta2 + (88*log(1 - x))/3. - 6*log(x)) + ((1 + pow(x,2))*((70*log(x))/3. - 8*li2(1 - x)))/(1 - x)) - 16*x*li2(1 - x) - (4*(7 + x)*li2(1 - x))/3. + ((1 + pow(x,2))*(208*log(x) - 48*zeta2*log(x) - 280*log(1 - x)*log(x) + 87*pow(log(x),2) + 12*log(1 - x)*pow(log(x),2) + 8*(-1 + 6*log(1 - x) - 6*log(x))*li2(1 - x) + 24*log(x)*li2(x) + 72*Li3(1 - x) - 24*Li3(x) + 24*zeta3))/(6.*(-1 + x)) + (1 + x)*(-28*zeta3 + 16*zeta2*log(1 - x) + (88*pow(log(1 - x),2))/3. + 8*log(1 - x)*li2(1 - x) - 12*Li3(1 - x) + 8*(log(1 - x)*pow(log(x),2) + 2*log(x)*li2(x) - 2*Li3(x) + 2*zeta3)))))/(16.*pow(M_PI,2));
}
double DY_NNLO_qqbar_NS_expansion(double x, int power){
	if(power==1){
		return (0.009259259259259259*pow(alphas_muR,2)*CF*(-45.564747186927704*nF + CF*(-864. - 864.*zeta2 - 1728.*zeta3) + CA*(701. - 504.*zeta2 - 378.*zeta3) - 54.*beta0*log(muR2/muF2)*(-1. + log(Q2/muF2) + 2.*log(1 - x)) + 3.*(pow(log(Q2/muF2),2)*(33.*CA - 36.*CF - 6.*nF - 144.*CF*log(1 - x)) + log(Q2/muF2)*(-266.*CA + 450.*CF + 44.*nF + 36.*CA*zeta2 + 144.*CF*zeta2 + (132.*CA + 288.*CF - 24.*nF)*log(1 - x) - 432.*CF*pow(log(1 - x),2)) + log(1 - x)*(-487.*CA + 639.*CF + 88.*nF + 72.*CA*zeta2 + 288.*CF*zeta2 + (132.*CA + 558.*CF - 24.*nF)*log(1 - x) - 288.*CF*pow(log(1 - x),2)))))/pow(M_PI,2);
	}
	if(power==2){
		return (0.004629629629629629*pow(alphas_muR,2)*CF*(1 - x)*(384.5647471869277*nF + 2592.*CF*(-0.7395833333333334 + 1.*zeta2 + 0.6666666666666666*zeta3) + CA*(-2921. + 558.*zeta2 + 378.*zeta3) + 54.*beta0*log(muR2/muF2)*(-1. + log(Q2/muF2) + 2.*log(1 - x)) - 99.*pow(log(Q2/muF2),2)*(1.*CA + 3.272727272727273*CF - 0.18181818181818182*nF - 4.363636363636363*CF*log(1 - x)) - 6.*log(Q2/muF2)*(-135.*CF + 34.*nF + 72.*CF*zeta2 + CA*(-235. + 18.*zeta2) + (66.*CA + 432.*CF - 12.*nF)*log(1 - x) - 216.*CF*pow(log(1 - x),2)) + 6.*log(1 - x)*(522.*CF - 68.*nF + CA*(470. - 36.*zeta2) - 144.*CF*zeta2 + (-66.*CA - 567.*CF + 12.*nF)*log(1 - x) + 144.*CF*pow(log(1 - x),2))))/pow(M_PI,2);
	}
	if(power==3){
		return (-2.6666666666666665*CF*pow(alphas_muR - alphas_muR*x,2)*(-1.953125*CA - 2.7760416666666665*CF + 0.2673611111111111*nF + 0.125*CA*zeta2 + 1.*CF*zeta2 - 0.0625*beta0*log(muR2/muF2) - 0.25*CF*pow(log(Q2/muF2),2) + 0.875*CA*log(1 - x) + 3.5104166666666665*CF*log(1 - x) - 0.16666666666666666*nF*log(1 - x) - 1.9375*CF*pow(log(1 - x),2) + log(Q2/muF2)*(0.4583333333333333*CA + 1.3125*CF - 0.08333333333333333*nF - 1.75*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==4){
		return (-1.*pow(alphas_muR,2)*CF*pow(1 - x,3)*(-1.1391782407407407*CA - 0.7841435185185185*CF + 0.09837962962962964*nF + 0.16666666666666666*CA*zeta2 + 1.*CF*zeta2 - 0.08333333333333333*beta0*log(muR2/muF2) - 0.25*CF*pow(log(Q2/muF2),2) + 0.9097222222222222*CA*log(1 - x) + 0.1597222222222222*CF*log(1 - x) - 0.2222222222222222*nF*log(1 - x) - 2.25*CF*pow(log(1 - x),2) + log(Q2/muF2)*(0.4861111111111111*CA + 0.09722222222222222*CF - 0.1111111111111111*nF - 2.*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==5){
		return (-0.6*pow(alphas_muR,2)*CF*pow(-1 + x,4)*(-0.794729938271605*CA - 0.879429012345679*CF - 0.027469135802469135*nF + 0.19444444444444445*CA*zeta2 + 1.*CF*zeta2 - 0.09722222222222222*beta0*log(muR2/muF2) - 0.25*CF*pow(log(Q2/muF2),2) + 0.9425925925925925*CA*log(1 - x) - 1.7087962962962964*CF*log(1 - x) - 0.25925925925925924*nF*log(1 - x) - 2.4583333333333335*CF*pow(log(1 - x),2) + log(Q2/muF2)*(0.5046296296296297*CA - 0.5046296296296297*CF - 0.12962962962962962*nF - 2.1666666666666665*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==6){
		return (-0.43333333333333335*pow(alphas_muR,2)*CF*pow(1 - x,5)*(-0.5449412393162393*CA - 1.3180555555555555*CF - 0.12585470085470085*nF + 0.21153846153846154*CA*zeta2 + 1.*CF*zeta2 - 0.10576923076923077*beta0*log(muR2/muF2) - 0.25*CF*pow(log(Q2/muF2),2) + 0.9423076923076923*CA*log(1 - x) - 2.910576923076923*CF*log(1 - x) - 0.28205128205128205*nF*log(1 - x) - 2.5865384615384617*CF*pow(log(1 - x),2) + log(Q2/muF2)*(0.5032051282051282*CA - 0.8532051282051282*CF - 0.14102564102564102*nF - 2.269230769230769*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==7){
		return (-0.34285714285714286*pow(alphas_muR,2)*CF*pow(-1 + x,6)*(-0.35370976631393297*CA - 1.790043934240363*CF - 0.20404541446208113*nF + 0.2222222222222222*CA*zeta2 + 1.*CF*zeta2 - 0.1111111111111111*beta0*log(muR2/muF2) - 0.25*CF*pow(log(Q2/muF2),2) + 0.921957671957672*CA*log(1 - x) - 3.7474206349206347*CF*log(1 - x) - 0.2962962962962963*nF*log(1 - x) - 2.6666666666666665*CF*pow(log(1 - x),2) + log(Q2/muF2)*(0.49074074074074076*CA - 1.074537037037037*CF - 0.14814814814814814*nF - 2.3333333333333335*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==8){
		return (-0.2857142857142857*pow(alphas_muR,2)*CF*pow(1 - x,7)*(-0.20691174473261525*CA - 2.226536192602041*CF - 0.26770419973544973*nF + 0.22916666666666666*CA*zeta2 + 1.*CF*zeta2 - 0.11458333333333333*beta0*log(muR2/muF2) - 0.25*CF*pow(log(Q2/muF2),2) + 0.8918402777777777*CA*log(1 - x) - 4.368229166666667*CF*log(1 - x) - 0.3055555555555556*nF*log(1 - x) - 2.71875*CF*pow(log(1 - x),2) + log(Q2/muF2)*(0.4732638888888889*CA - 1.2259920634920636*CF - 0.1527777777777778*nF - 2.375*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==9){
		return (-0.24603174603174605*pow(alphas_muR,2)*CF*pow(-1 + x,8)*(-0.09404584242597876*CA - 2.6163036110501547*CF - 0.3208226659839563*nF + 0.23387096774193547*CA*zeta2 + 1.*CF*zeta2 - 0.11693548387096774*beta0*log(muR2/muF2) - 0.25*CF*pow(log(Q2/muF2),2) + 0.8576164874551971*CA*log(1 - x) - 4.852441756272402*CF*log(1 - x) - 0.3118279569892473*nF*log(1 - x) - 2.754032258064516*CF*pow(log(1 - x),2) + log(Q2/muF2)*(0.4538978494623656*CA - 1.3361175115207373*CF - 0.15591397849462366*nF - 2.403225806451613*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==10){
		return (-0.21666666666666667*pow(alphas_muR,2)*CF*pow(1 - x,9)*(-0.007061538824931682*CA - 2.9625655117222576*CF - 0.3661239994573328*nF + 0.23717948717948717*CA*zeta2 + 1.*CF*zeta2 - 0.11858974358974358*beta0*log(muR2/muF2) - 0.25*CF*pow(log(Q2/muF2),2) + 0.8221916971916972*CA*log(1 - x) - 5.2451745014245015*CF*log(1 - x) - 0.3162393162393162*nF*log(1 - x) - 2.7788461538461537*CF*pow(log(1 - x),2) + log(Q2/muF2)*(0.4341727716727717*CA - 1.4201821326821327*CF - 0.1581196581196581*nF - 2.423076923076923*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==11){
		return (-0.19393939393939394*pow(alphas_muR,2)*CF*pow(-1 + x,10)*(0.06006358035997342*CA - 3.2714887306777753*CF - 0.4054776434984768*nF + 0.23958333333333334*CA*zeta2 + 1.*CF*zeta2 - 0.11979166666666667*beta0*log(muR2/muF2) - 0.25*CF*pow(log(Q2/muF2),2) + 0.7870242604617604*CA*log(1 - x) - 5.573573532948533*CF*log(1 - x) - 0.3194444444444444*nF*log(1 - x) - 2.796875*CF*pow(log(1 - x),2) + log(Q2/muF2)*(0.4148189484126984*CA - 1.4868964947089947*CF - 0.1597222222222222*nF - 2.4375*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==12){
		return (-0.17575757575757575*pow(alphas_muR,2)*CF*pow(1 - x,11)*(0.1117950450405271*CA - 3.549069881542592*CF - 0.44019334353386075*nF + 0.2413793103448276*CA*zeta2 + 1.*CF*zeta2 - 0.1206896551724138*beta0*log(muR2/muF2) - 0.25*CF*pow(log(Q2/muF2),2) + 0.7528256704980842*CA*log(1 - x) - 5.854837786734339*CF*log(1 - x) - 0.3218390804597701*nF*log(1 - x) - 2.810344827586207*CF*pow(log(1 - x),2) + log(Q2/muF2)*(0.3961685823754789*CA - 1.5415189456137732*CF - 0.16091954022988506*nF - 2.4482758620689653*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==13){
		return (8.975105362789164e-15*pow(alphas_muR,2)*CF*pow(-1 + x,12)*(-2.714151079371e12*CA + 6.8105519568095e13*CF + 8.44438257663e12*nF - 4.3502875416e12*CA*zeta2 - 1.79205874848e13*CF*zeta2 + 2.1751437708e12*beta0*log(muR2/muF2) + 4.4801468712e12*CF*pow(log(Q2/muF2),2) - 1.290137772924e13*CA*log(1 - x) + 1.0932214545252e14*CF*log(1 - x) + 5.8003833888e12*nF*log(1 - x) + 5.05477440468e13*CF*pow(log(1 - x),2) + log(Q2/muF2)*(-6.78033177822e12*CA + 2.844683623782e13*CF + 2.9001916944e12*nF + 4.40223127344e13*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==14){
		return (7.051868499334342e-15*pow(alphas_muR,2)*CF*pow(1 - x,13)*(-3.818983687881e12*CA + 8.4772727209601e13*CF + 1.050225475299e13*nF - 5.1294435192e12*CA*zeta2 - 2.10372113952e13*CF*zeta2 + 2.5647217596e12*beta0*log(muR2/muF2) + 5.2593028488e12*CF*pow(log(Q2/muF2),2) - 1.448255236428e13*CA*log(1 - x) + 1.3291338456396e14*CF*log(1 - x) + 6.8392580256e12*nF*log(1 - x) + 5.95080377892e13*CF*pow(log(1 - x),2) + log(Q2/muF2)*(-7.60302717174e12*CA + 3.42210855987e13*CF + 3.4196290128e12*nF + 5.18138725104e13*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==15){
		return (2.820747399733737e-15*pow(alphas_muR,2)*CF*pow(-1 + x,14)*(-9.957279257753e12*CA + 2.07035229997002e14*CF + 2.562192140508e13*nF - 1.19470583232e13*CA*zeta2 - 4.88271079296e13*CF*zeta2 + 5.9735291616e12*beta0*log(muR2/muF2) + 1.22067769824e13*CF*pow(log(Q2/muF2),2) - 3.21457119984e13*CA*log(1 - x) + 3.1802936701536e14*CF*log(1 - x) + 1.59294110976e13*nF*log(1 - x) + 1.384300453536e14*CF*pow(log(1 - x),2) + log(Q2/muF2)*(-1.686066926544e13*CA + 8.109990953064e13*CF + 7.9647055488e12*nF + 1.205094578688e14*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==16){
		return (5.729643155709154e-16*pow(alphas_muR,2)*CF*pow(1 - x,15)*(-4.9388383599927e13*CA + 9.95122143460498e14*CF + 1.2301122281448e14*nF - 5.50603557504e13*CA*zeta2 - 2.243969215488e14*CF*zeta2 + 2.75301778752e13*beta0*log(muR2/muF2) + 5.60992303872e13*CF*pow(log(Q2/muF2),2) - 1.4129732753136e14*CA*log(1 - x) + 1.50134181685488e15*CF*log(1 - x) + 7.34138076672e13*nF*log(1 - x) + 6.373495896768e14*CF*pow(log(1 - x),2) + log(Q2/muF2)*(-7.405747116768e13*CA + 3.795126439104e14*CF + 3.67069038336e13*nF + 5.547590560512e14*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==17){
		return (1.6327092240978657e-18*pow(alphas_muR,2)*CF*pow(-1 + x,16)*(-1.7073038825791624e16*CA + 3.408752955471495e17*CF + 4.208776537117152e16*nF - 1.81642035871296e16*CA*zeta2 - 7.38577534286592e16*CF*zeta2 + 9.0821017935648e15*beta0*log(muR2/muF2) + 1.84644383571648e16*CF*pow(log(Q2/muF2),2) - 4.448542608173664e16*CA*log(1 - x) + 5.061221546204023e17*CF*log(1 - x) + 2.42189381161728e16*nF*log(1 - x) + 2.100892803321312e17*CF*pow(log(1 - x),2) + log(Q2/muF2)*(-2.330236517040432e16*CA + 1.269140033417256e17*CF + 1.21094690580864e16*nF + 1.828429749514368e17*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==18){
		return (1.3605910200815549e-18*pow(alphas_muR,2)*CF*pow(1 - x,17)*(-1.9866648813263e16*CA + 3.992889224346774e17*CF + 4.924329073888224e16*nF - 2.05660817474112e16*CA*zeta2 - 8.34652660697856e16*CF*zeta2 + 1.02830408737056e16*beta0*log(muR2/muF2) + 2.08663165174464e16*CF*pow(log(Q2/muF2),2) - 4.809172727110752e16*CA*log(1 - x) + 5.844275982983037e17*CF*log(1 - x) + 2.74214423298816e16*nF*log(1 - x) + 2.377108791753696e17*CF*pow(log(1 - x),2) + log(Q2/muF2)*(-2.518008387790896e16*CA + 1.4546696638730256e17*CF + 1.37107211649408e16*nF + 2.068617565542528e17*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==19){
		return (2.2216969565408095e-20*pow(alphas_muR,2)*CF*pow(-1 + x,18)*(-1.1661189872426816e18*CA + 2.3873773814285324e19*CF + 2.941007826463381e18*nF - 1.1922322718097792e18*CA*zeta2 - 4.830863231229235e18*CF*zeta2 + 5.961161359048896e17*beta0*log(muR2/muF2) + 1.2077158078073088e18*CF*pow(log(Q2/muF2),2) - 2.662885928733751e18*CA*log(1 - x) + 3.449515786851667e19*CF*log(1 - x) + 1.5896430290797056e18*nF*log(1 - x) + 1.377260526980258e19*CF*pow(log(1 - x),2) + log(Q2/muF2)*(-1.393784569830613e18*CA + 8.527240987112391e18*CF + 7.948215145398528e17*nF + 1.198425686208791e19*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==20){
		return (2.6977748757995544e-21*pow(alphas_muR,2)*CF*pow(1 - x,19)*(-9.121360539772599e18*CA + 1.9200775595277523e20*CF + 2.362810857960521e19*nF - 9.32108867051282e18*CA*zeta2 - 3.771789368998211e19*CF*zeta2 + 4.66054433525641e18*beta0*log(muR2/muF2) + 9.429473422495527e18*CF*pow(log(Q2/muF2),2) - 1.9890079573488284e19*CA*log(1 - x) + 2.7419941874478272e20*CF*log(1 - x) + 1.2428118227350426e19*nF*log(1 - x) + 1.0762605871882825e20*CF*pow(log(1 - x),2) + log(Q2/muF2)*(-1.0408384601470216e19*CA + 6.735069677077678e19*CF + 6.214059113675213e18*nF + 9.364442571305902e19*CF*log(1 - x))))/pow(M_PI,2);
	}
	if(power==21){
		return (4.624756929942094e-21*pow(alphas_muR,2)*CF*pow(-1 + x,20)*(-5.016437580478922e18*CA + 1.0942400880300922e20*CF + 1.345176155413112e19*nF - 5.175371907174269e18*CA*zeta2 - 2.091825713266249e19*CF*zeta2 + 2.5876859535871345e18*beta0*log(muR2/muF2) + 5.229564283165622e18*CF*pow(log(Q2/muF2),2) - 1.0552343319674237e19*CA*log(1 - x) + 1.546014403994688e20*CF*log(1 - x) + 6.900495876232358e18*nF*log(1 - x) + 5.973354643646951e19*CF*pow(log(1 - x),2) + log(Q2/muF2)*(-5.52132764646467e18*CA + 3.774839763246948e19*CF + 3.450247938116179e18*nF + 5.19704885757081e19*CF*log(1 - x))))/pow(M_PI,2);
	}
}

double DY_NNLO_BB_full(double x){
	return -(pow(alphas_muR,2)*CF*((1 + x)*(15*(-1 + x) + pow(M_PI,2)*(1 + x)) - 3*pow(1 + x,2)*pow(log(x),2) + 3*log(x)*(-3 - 4*x - 3*pow(x,2) + 4*pow(1 + x,2)*log(1 + x)) + 12*pow(1 + x,2)*li2(-x)))/(18.*pow(M_PI,2));
}

double DY_NNLO_BB_expansion(double x, int power){
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
		return 0;
	}
	if(power==5){
		return 0;
	}
	if(power==6){
		return (-0.016666666666666666*pow(alphas_muR,2)*CF*pow(-1. + x,5))/pow(M_PI,2);
	}
	if(power==7){
		return (0.025*pow(alphas_muR,2)*CF*pow(-1. + x,6))/pow(M_PI,2);
	}
	if(power==8){
		return (-0.029365079365079365*pow(alphas_muR,2)*CF*pow(-1. + x,7))/pow(M_PI,2);
	}
	if(power==9){
		return (0.031746031746031744*pow(alphas_muR,2)*CF*pow(-1. + x,8))/pow(M_PI,2);
	}
	if(power==10){
		return (-0.03302469135802469*pow(alphas_muR,2)*CF*pow(-1. + x,9))/pow(M_PI,2);
	}
	if(power==11){
		return (0.03364197530864198*pow(alphas_muR,2)*CF*pow(-1. + x,10))/pow(M_PI,2);
	}
	if(power==12){
		return (-0.03384239217572551*pow(alphas_muR,2)*CF*pow(-1. + x,11))/pow(M_PI,2);
	}
	if(power==13){
		return (0.03377224627224627*pow(alphas_muR,2)*CF*pow(-1. + x,12))/pow(M_PI,2);
	}
	if(power==14){
		return (-0.033523883523883524*pow(alphas_muR,2)*CF*pow(-1. + x,13))/pow(M_PI,2);
	}
	if(power==15){
		return (0.03315781440781441*pow(alphas_muR,2)*CF*pow(-1. + x,14))/pow(M_PI,2);
	}
	if(power==16){
		return (-0.032714754381421046*pow(alphas_muR,2)*CF*pow(-1. + x,15))/pow(M_PI,2);
	}
	if(power==17){
		return (0.03222262305595639*pow(alphas_muR,2)*CF*pow(-1. + x,16))/pow(M_PI,2);
	}
	if(power==18){
		return (-0.03170082277925415*pow(alphas_muR,2)*CF*pow(-1. + x,17))/pow(M_PI,2);
	}
	if(power==19){
		return (0.031162956959525586*pow(alphas_muR,2)*CF*pow(-1. + x,18))/pow(M_PI,2);
	}
	if(power==20){
		return (-0.030618609785276453*pow(alphas_muR,2)*CF*pow(-1. + x,19))/pow(M_PI,2);
	}
	if(power==21){
		return (0.03007453965787299*pow(alphas_muR,2)*CF*pow(-1. + x,20))/pow(M_PI,2);
	}
}

double DY_NNLO_BC_full(double x){
	return (pow(alphas_muR,2)*CF*(-CA/2. + CF)*(81 - 42*x - 39*pow(x,2) + 6*(9 + 11*x)*log(x) - 3*(-6 + 8*x + 15*pow(x,2))*pow(log(x),2) + 2*(1 + 4*x + pow(x,2))*pow(log(x),3) - 54*(-1 + pow(x,2))*li2(1 - x) + 24*(1 + 3*x + pow(x,2))*(log(x)*li2(1 - x) + 2*S12(1 - x)) + pow(1 + x,2)*(3*pow(M_PI,2) + 2*pow(M_PI,2)*log(x) - 6*pow(M_PI,2)*log(1 + x) + 36*log(x)*log(1 + x) + 30*pow(log(x),2)*log(1 + x) - 36*log(x)*pow(log(1 + x),2) + 36*(1 + log(x) - 2*log(1 + x))*li2(-x) - 12*Li3(-x) - 72*S12(-x))))/(24.*pow(M_PI,2));
}

double DY_NNLO_BC_expansion(double x, int power){
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
		return (-0.020833333333333332*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,3))/pow(M_PI,2);
	}
	if(power==5){
		return (0.043402777777777776*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,4))/pow(M_PI,2);
	}
	if(power==6){
		return (-0.04756944444444444*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,5))/pow(M_PI,2);
	}
	if(power==7){
		return (0.048136574074074075*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,6))/pow(M_PI,2);
	}
	if(power==8){
		return (-0.04766203703703704*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,7))/pow(M_PI,2);
	}
	if(power==9){
		return (0.04678689531368103*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,8))/pow(M_PI,2);
	}
	if(power==10){
		return (-0.04573436318972033*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,9))/pow(M_PI,2);
	}
	if(power==11){
		return (0.04460646179768204*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,10))/pow(M_PI,2);
	}
	if(power==12){
		return (-0.04345828967108729*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,11))/pow(M_PI,2);
	}
	if(power==13){
		return (0.042321788096296574*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,12))/pow(M_PI,2);
	}
	if(power==14){
		return (-0.0412156193515546*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,13))/pow(M_PI,2);
	}
	if(power==15){
		return (0.04015027100625778*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,14))/pow(M_PI,2);
	}
	if(power==16){
		return (-0.039131069305184375*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,15))/pow(M_PI,2);
	}
	if(power==17){
		return (0.03816007655097998*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,16))/pow(M_PI,2);
	}
	if(power==18){
		return (-0.03723731753309373*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,17))/pow(M_PI,2);
	}
	if(power==19){
		return (0.03636157868514285*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,18))/pow(M_PI,2);
	}
	if(power==20){
		return (-0.03553092848292547*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,19))/pow(M_PI,2);
	}
	if(power==21){
		return (0.03474305433512218*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,20))/pow(M_PI,2);
	}
}


////////////////////////////////////////////////////////////
///
/// qq + qbarqbar channel
///
////////////////////////////////////////////////////////////
//NNLO functions
double DY_NNLO_CC_full(double x){
	return (pow(alphas_muR,2)*CF*TF*(1779 + 116/x - 2598*x + 703*pow(x,2) - 24*(39 - 22/x - 39*x + 22*pow(x,2))*log(1 - x) + (72*(-4 - 3*x + 3*pow(x,2) + 4*pow(x,3))*(zeta2 - pow(log(1 - x),2)))/x + 6*(345 - 48*x + 20*pow(x,2))*log(x) + 216*(3 + 6*x + 4*pow(x,2))*log(1 - x)*log(x) - 45*(3 + 15*x + 8*pow(x,2))*pow(log(x),2) + 18*pow(log(Q2/muF2),2)*(3 + 4/x - 3*x - 4*pow(x,2) + 6*(1 + x)*log(x)) + 36*(39 + 16/x + 15*x + 8*pow(x,2))*li2(1 - x) + 12*log(Q2/muF2)*(-39 + 22/x + 39*x - 22*pow(x,2) + 6*(3 + 4/x - 3*x - 4*pow(x,2))*log(1 - x) + 9*(3 + 6*x + 4*pow(x,2))*log(x) + 18*(1 + x)*((2*log(1 - x) - log(x))*log(x) + 2*li2(1 - x))) - 18*(1 + x)*(4*pow(M_PI,2)*log(x) - 24*pow(log(1 - x),2)*log(x) - 12*log(1 - x)*pow(log(x),2) - 9*pow(log(x),3) - 12*(4*log(1 - x) + log(x))*li2(1 - x) - 72*log(x)*li2(x) + 48*Li3(1 - x) + 72*Li3(x) - 72*zeta3)))/(432.*pow(M_PI,2));
}

double DY_NNLO_CC_expansion(double x, int power){
	if(power==1){
		return 0;
	}
	if(power==2){
		return (0.041666666666666664*pow(alphas_muR,2)*CF*TF*(1 - x)*(105.95683520871486 - 60.*zeta2 + 3.*pow(log(Q2/muF2),2) + 12.*log(Q2/muF2)*(-1. + log(1 - x)) + 12.*(-2. + log(1 - x))*log(1 - x)))/pow(M_PI,2);
	}
	if(power==3){
		return (0.375*pow(alphas_muR,2)*CF*TF*pow(-1. + x,2)*(-2. + log(Q2/muF2) + 2.*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==4){
		return (-0.6666666666666666*pow(alphas_muR,2)*CF*TF*pow(1 - x,3)*(-1.5362335167120569 + 1.*zeta2 - 0.1875*pow(log(Q2/muF2),2) + log(Q2/muF2)*(0.375 - 0.75*log(1 - x)) + 0.75*log(1 - x) - 0.75*pow(log(1 - x),2)))/pow(M_PI,2);
	}
	if(power==5){
		return (-0.6666666666666666*pow(alphas_muR,2)*CF*TF*pow(-1 + x,4)*(-0.8018585167120565 + 1.*zeta2 - 0.1875*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-0.09375 - 0.75*log(1 - x)) - 0.1875*log(1 - x) - 0.75*pow(log(1 - x),2)))/pow(M_PI,2);
	}
	if(power==6){
		return (-0.6666666666666666*pow(alphas_muR,2)*CF*TF*pow(1 - x,5)*(-1.040582387263073 + 1.*zeta2 - 0.19375*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-0.21791666666666668 - 0.775*log(1 - x)) - 0.43583333333333335*log(1 - x) - 0.775*pow(log(1 - x),2)))/pow(M_PI,2);
	}
	if(power==7){
		return (-0.6666666666666667*pow(alphas_muR,2)*CF*TF*pow(-1 + x,6)*(-1.247153480036312 + 1.*zeta2 - 0.2*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-0.28583333333333333 - 0.8*log(1 - x)) - 0.5716666666666667*log(1 - x) - 0.8*pow(log(1 - x),2)))/pow(M_PI,2);
	}
	if(power==8){
		return (-0.6666666666666667*pow(alphas_muR,2)*CF*TF*pow(1 - x,7)*(-1.406460225736989 + 1.*zeta2 - 0.20535714285714285*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-0.33562925170068025 - 0.8214285714285714*log(1 - x)) - 0.6712585034013605*log(1 - x) - 0.8214285714285714*pow(log(1 - x),2)))/pow(M_PI,2);
	}
	if(power==9){
		return (-0.6666666666666666*pow(alphas_muR,2)*CF*TF*pow(-1 + x,8)*(-1.5343671368367593 + 1.*zeta2 - 0.20982142857142858*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-0.37623299319727893 - 0.8392857142857143*log(1 - x)) - 0.7524659863945579*log(1 - x) - 0.8392857142857143*pow(log(1 - x),2)))/pow(M_PI,2);
	}
	if(power==10){
		return (-0.6666666666666667*pow(alphas_muR,2)*CF*TF*pow(1 - x,9)*(-1.6412970872136363 + 1.*zeta2 - 0.21354166666666666*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-0.41077215608465606 - 0.8541666666666666*log(1 - x)) - 0.8215443121693121*log(1 - x) - 0.8541666666666666*pow(log(1 - x),2)))/pow(M_PI,2);
	}
	if(power==11){
		return (-0.6666666666666667*pow(alphas_muR,2)*CF*TF*pow(-1 + x,10)*(-1.7333452863465395 + 1.*zeta2 - 0.21666666666666667*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-0.44076719576719575 - 0.8666666666666667*log(1 - x)) - 0.8815343915343915*log(1 - x) - 0.8666666666666667*pow(log(1 - x),2)))/pow(M_PI,2);
	}
	if(power==12){
		return (-0.6666666666666667*pow(alphas_muR,2)*CF*TF*pow(1 - x,11)*(-1.8142249639399601 + 1.*zeta2 - 0.21931818181818183*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-0.4671472845336482 - 0.8772727272727273*log(1 - x)) - 0.9342945690672964*log(1 - x) - 0.8772727272727273*pow(log(1 - x),2)))/pow(M_PI,2);
	}
	if(power==13){
		return (-0.6666666666666667*pow(alphas_muR,2)*CF*TF*pow(-1 + x,12)*(-1.8863441905550307 + 1.*zeta2 - 0.2215909090909091*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-0.49056145874327695 - 0.8863636363636364*log(1 - x)) - 0.9811229174865539*log(1 - x) - 0.8863636363636364*pow(log(1 - x),2)))/pow(M_PI,2);
	}
	if(power==14){
		return (1.480892384860212e-13*pow(alphas_muR,2)*CF*TF*pow(1 - x,13)*(8.784612346567078e12 - 4.5017900928e12*zeta2 + 1.0064098044e12*pow(log(Q2/muF2),2) + 4.60530398328e12*log(1 - x) + 4.0256392176e12*pow(log(1 - x),2) + log(Q2/muF2)*(2.30265199164e12 + 4.0256392176e12*log(1 - x))))/pow(M_PI,2);
	}
	if(power==15){
		return (1.9745231798136158e-13*pow(alphas_muR,2)*CF*TF*pow(-1 + x,14)*(6.78803657832386e12 - 3.3763425696e12*zeta2 + 7.606046448e11*pow(log(Q2/muF2),2) + 3.58118159256e12*log(1 - x) + 3.0424185792e12*pow(log(1 - x),2) + log(Q2/muF2)*(1.79059079628e12 + 3.0424185792e12*log(1 - x))))/pow(M_PI,2);
	}
	if(power==16){
		return (-0.6666666666666666*pow(alphas_muR,2)*CF*TF*pow(1 - x,15)*(-2.064578340428136 + 1.*zeta2 - 0.22678571428571428*pow(log(Q2/muF2),2) + log(Q2/muF2)*(-0.5473784746999033 - 0.9071428571428571*log(1 - x)) - 1.0947569493998066*log(1 - x) - 0.9071428571428571*pow(log(1 - x),2)))/pow(M_PI,2);
	}
	if(power==17){
		return (3.2086001671971257e-13*pow(alphas_muR,2)*CF*TF*pow(-1 + x,16)*(4.393163676202957e12 - 2.0777492736e12*zeta2 + 4.7398655304e11*pow(log(Q2/muF2),2) + 2.339033264568e12*log(1 - x) + 1.89594621216e12*pow(log(1 - x),2) + log(Q2/muF2)*(1.169516632284e12 + 1.89594621216e12*log(1 - x))))/pow(M_PI,2);
	}
	if(power==18){
		return (3.2654184481957314e-17*pow(alphas_muR,2)*CF*TF*pow(1 - x,17)*(4.410773024026173e16 - 2.04159643623936e16*zeta2 + 4.6817859452364e15*pow(log(Q2/muF2),2) + 2.356136798563548e16*log(1 - x) + 1.87271437809456e16*pow(log(1 - x),2) + log(Q2/muF2)*(1.178068399281774e16 + 1.87271437809456e16*log(1 - x))))/pow(M_PI,2);
	}
	if(power==19){
		return (1.5674008551339512e-15*pow(alphas_muR,2)*CF*TF*pow(-1 + x,18)*(9.3710789598575e14 - 4.253325908832e14*zeta2 + 9.79932929976e13*pow(log(Q2/muF2),2) + 5.0190664111588e14*log(1 - x) + 3.919731719904e14*pow(log(1 - x),2) + log(Q2/muF2)*(2.5095332055794e14 + 3.919731719904e14*log(1 - x))))/pow(M_PI,2);
	}
	if(power==20){
		return (4.855994776439198e-19*pow(alphas_muR,2)*CF*TF*pow(1 - x,19)*(3.079507919192709e18 - 1.3728735251142912e18*zeta2 + 3.176275370604336e17*pow(log(Q2/muF2),2) + 1.6528571271425014e18*log(1 - x) + 1.2705101482417344e18*pow(log(1 - x),2) + log(Q2/muF2)*(8.264285635712507e17 + 1.2705101482417344e18*log(1 - x))))/pow(M_PI,2);
	}
	if(power==21){
		return (9.711989552878396e-20*pow(alphas_muR,2)*CF*TF*pow(-1 + x,20)*(1.5653478268779334e19 - 6.864367625571456e18*zeta2 + 1.5941590604123185e18*pow(log(Q2/muF2),2) + 8.415888167436029e18*log(1 - x) + 6.376636241649274e18*pow(log(1 - x),2) + log(Q2/muF2)*(4.2079440837180145e18 + 6.376636241649274e18*log(1 - x))))/pow(M_PI,2);
	}

}
double DY_NNLO_CD_full(double x){
	return (pow(alphas_muR,2)*CF*TF*(160*(-1 + x) - 24*(-6 + 4/x + x)*zeta3 - 16*(5 + 4*x)*log(x) + 8*(10 + x)*zeta2*log(x) - 52*x*pow(log(x),2) - (16*x*pow(log(x),3))/3. + 8*(5 - 4*x)*li2(1 - x) - 8*(-10 + x)*log(x)*li2(1 - x) + 32*(5/x + 2*x)*log(x)*li2(-x) + 40*(1 + x)*(zeta2 + 2*log(x)*log(1 + x) + 2*li2(-x)) + 8*(-6 + 4/x + 3*x)*Li3(1 - x) - 16*(-10 + 10/x + 3*x)*Li3(-x) - (8*(2 + 2*x + pow(x,2))*(log(1 + x)*(6*zeta2 - 5*pow(log(x),2) + 6*log(x)*log(1 + x) + 12*li2(-x)) - 4*S12(1 - x) + 12*S12(-x)))/x))/(16.*pow(M_PI,2));
}

double DY_NNLO_CD_expansion(double x, int power){
	if(power==1){
		return (-5.397207708399179*pow(alphas_muR,2)*CF*TF*(-1.3108566230742145 + 1.*zeta2 - 0.27792148848851744*zeta3))/pow(M_PI,2);
	}
	if(power==2){
		return (2.5794415416798353*pow(alphas_muR,2)*CF*TF*(-1. + x)*(-3.7419988682065273 + 1.*zeta2 + 1.7445636690294675*zeta3))/pow(M_PI,2);
	}
	if(power==3){
		return (-3.0338830833596715*pow(alphas_muR,2)*CF*TF*pow(-1. + x,2)*(-4.16640282776246 + 1.*zeta2 + 1.9776635536514149*zeta3))/pow(M_PI,2);
	}
	if(power==4){
		return (1.742216416693005*pow(alphas_muR,2)*CF*TF*pow(-1. + x,3)*(-6.31747766229338 + 1.*zeta2 + 3.443889026937838*zeta3))/pow(M_PI,2);
	}
	if(power==5){
		return (-1.2578414166930052*pow(alphas_muR,2)*CF*TF*pow(-1. + x,4)*(-8.296920993543955 + 1.*zeta2 + 4.770076672920039*zeta3))/pow(M_PI,2);
	}
	if(power==6){
		return (0.9932580833596714*pow(alphas_muR,2)*CF*TF*pow(-1. + x,5)*(-10.150366084644945 + 1.*zeta2 + 6.040726071621934*zeta3))/pow(M_PI,2);
	}
	if(power==7){
		return (-0.8239872500263381*pow(alphas_muR,2)*CF*TF*pow(-1. + x,6)*(-11.908682565588611 + 1.*zeta2 + 7.281666069236162*zeta3))/pow(M_PI,2);
	}
	if(power==8){
		return (0.7054604643120527*pow(alphas_muR,2)*CF*TF*pow(-1. + x,7)*(-13.598025995701239 + 1.*zeta2 + 8.505083280394812*zeta3))/pow(M_PI,2);
	}
	if(power==9){
		return (-0.617407078895386*pow(alphas_muR,2)*CF*TF*pow(-1. + x,8)*(-15.238838605245213 + 1.*zeta2 + 9.71806155953817*zeta3))/pow(M_PI,2);
	}
	if(power==10){
		return (0.5492027138160207*pow(alphas_muR,2)*CF*TF*pow(-1. + x,9)*(-16.84576727564191 + 1.*zeta2 + 10.924927807275843*zeta3))/pow(M_PI,2);
	}
	if(power==11){
		return (-0.4947105263160207*pow(alphas_muR,2)*CF*TF*pow(-1. + x,10)*(-18.428774157559562 + 1.*zeta2 + 12.128304697052688*zeta3))/pow(M_PI,2);
	}
	if(power==12){
		return (-1.1737087191695113e-14*pow(alphas_muR,2)*CF*TF*pow(-1. + x,11)*(7.667943428672704e14 - 3.835032234021594e13*zeta2 - 5.11200087552e14*zeta3))/pow(M_PI,2);
	}
	if(power==13){
		return (7.335679494809446e-16*pow(alphas_muR,2)*CF*TF*pow(-1. + x,12)*(1.2129192621058598e16 - 5.62912258136406e14*zeta2 - 8.179201400832e15*zeta3))/pow(M_PI,2);
	}
	if(power==14){
		return (0.38143598212106644*pow(alphas_muR,2)*CF*TF*pow(-1. + x,13)*(-23.089928805927016 + 1.*zeta2 + 15.73003146330233*zeta3))/pow(M_PI,2);
	}
	if(power==15){
		return (-0.35440758042835796*pow(alphas_muR,2)*CF*TF*pow(-1. + x,14)*(-24.624589800017898 + 1.*zeta2 + 16.929660456889906*zeta3))/pow(M_PI,2);
	}
	if(power==16){
		return (-3.338952887942397e-19*pow(alphas_muR,2)*CF*TF*pow(-1. + x,15)*(2.592252919644216e19 - 9.912029198682685e17*zeta2 - 1.7969705477627904e19*zeta3))/pow(M_PI,2);
	}
	if(power==17){
		return (-0.31041930088309444*pow(alphas_muR,2)*CF*TF*pow(-1. + x,16)*(-27.674962481359795 + 1.*zeta2 + 19.328695035814256*zeta3))/pow(M_PI,2);
	}
	if(power==18){
		return (0.29228065473860143*pow(alphas_muR,2)*CF*TF*pow(-1. + x,17)*(-29.19245690426609 + 1.*zeta2 + 20.528214586648048*zeta3))/pow(M_PI,2);
	}
	if(power==19){
		return (-0.27614437934119485*pow(alphas_muR,2)*CF*TF*pow(-1. + x,18)*(-30.705677167827513 + 1.*zeta2 + 21.727764346731817*zeta3))/pow(M_PI,2);
	}
	if(power==20){
		return (0.2616961858904876*pow(alphas_muR,2)*CF*TF*pow(-1. + x,19)*(-32.21510386170526 + 1.*zeta2 + 22.927349818200366*zeta3))/pow(M_PI,2);
	}
	if(power==21){
		return (-0.24868434358668168*pow(alphas_muR,2)*CF*TF*pow(-1. + x,20)*(-33.72113251511357 + 1.*zeta2 + 24.126971217665876*zeta3))/pow(M_PI,2);
	}
}
double DY_NNLO_CE_full(double x){
	return (pow(alphas_muR,2)*CF*(-CA/2. + CF)*(2*(-9 + 7*x)*log(x) - 4*(1 + 3*x)*pow(log(x),2) + 8*(3 + x)*li2(1 - x) + 4*(1 + x)*(zeta2 + 4*log(1 - x)*log(x) + 2*log(x)*log(1 + x) + 2*li2(-x)) + log(Q2/muF2)*(-16*(-1 + x) + 8*(1 + x)*log(x) - (4*(1 + pow(x,2))*(2*zeta2 - pow(log(x),2) + 4*log(x)*log(1 + x) + 4*li2(-x)))/(1 + x)) + (1 - x)*(-34 + 8*zeta3 + 32*log(1 - x) + 4*zeta2*log(x) - (2*pow(log(x),3))/3. - (4*pow(M_PI,2)*log(1 + x))/3. + 4*pow(log(x),2)*log(1 + x) - 8*log(x)*pow(log(1 + x),2) - 16*log(1 + x)*li2(-x) + 8*Li3(-x) - 16*S12(-x)) - (4*(1 + pow(x,2))*(3*zeta3 + 12*zeta2*log(1 - x) - 9*zeta2*log(x) - 6*log(1 - x)*pow(log(x),2) + 2*pow(log(x),3) + 6*zeta2*log(1 + x) + 24*log(1 - x)*log(x)*log(1 + x) - 21*pow(log(x),2)*log(1 + x) + 6*log(x)*pow(log(1 + x),2) - 18*log(x)*li2(1 - x) + 24*log(1 - x)*li2(-x) - 24*log(x)*li2(-x) + 12*log(1 + x)*li2(-x) + 24*Li3(1 - x) + 6*Li3(-x) - 24*Li3((1 - x)/(1 + x)) + 24*Li3((-1 + x)/(1 + x)) - 24*S12(1 - x) + 12*S12(-x)))/(3.*(1 + x))))/(16.*pow(M_PI,2));
}

double DY_NNLO_CE_expansion(double x, int power){
	if(power==1){
		return (0.25*pow(alphas_muR,2)*(CA - 2.*CF)*CF*(-0.09627579537442348 - 0.3068528194400547*zeta2 + 0.5*zeta3 + (-1.6449340668482264 + 1.*zeta2)*log(Q2/muF2) + (-3.289868133696453 + 2.*zeta2)*log(1 - x)))/pow(M_PI,2);
	}
	if(power==2){
		return (0.125*pow(alphas_muR,2)*(CA - 2.*CF)*CF*(-1. + x)*(1.2894785320028408 - 2.3068528194400546*zeta2 + 2.5*zeta3 + (-1.6449340668482264 + 1.*zeta2)*log(Q2/muF2) + (-3.289868133696453 + 2.*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==3){
		return (0.0625*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,2)*(-5.853545029343216 + 3.1931471805599454*zeta2 + 0.5*zeta3 + (-1.6449340668482264 + 1.*zeta2)*log(Q2/muF2) + (-3.289868133696453 + 2.*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==4){
		return (-0.03125*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,3)*(-9.906702540938484 + 5.859813847226612*zeta2 + 0.5*zeta3 + (-1.6449340668482264 + 1.*zeta2)*log(Q2/muF2) + (-3.289868133696453 + 2.*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==5){
		return (0.015625*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,4)*(-14.290183313902078 + 9.943147180559945*zeta2 + 0.5*zeta3 + (-1.6449340668482264 + 1.*zeta2)*log(Q2/muF2) + (-3.289868133696453 + 2.*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==6){
		return (-0.0078125*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,5)*(-14.985250352872097 + 16.50981384722661*zeta2 + 0.5*zeta3 + (-3.2449340668482263 + 1.*zeta2)*log(Q2/muF2) + (-6.4898681336964525 + 2.*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==7){
		return (0.00390625*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,6)*(-12.753566520783075 + 27.443147180559944*zeta2 + 0.5*zeta3 + (-8.044934066848226 + 1.*zeta2)*log(Q2/muF2) + (-16.089868133696452 + 2.*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==8){
		return (-0.001953125*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,7)*(-9.512042751598178 + 46.12886146627423*zeta2 + 0.5*zeta3 + (-19.016362638276796 + 1.*zeta2)*log(Q2/muF2) + (-38.03272527655359 + 2.*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==9){
		return (0.0009765625*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,8)*(-10.309714993437217 + 78.71814718055995*zeta2 + 0.5*zeta3 + (-41.87350549541966 + 1.*zeta2)*log(Q2/muF2) + (-83.74701099083931 + 2.*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==10){
		return (-0.00048828125*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,9)*(-27.656854067562765 + 136.47608368849646*zeta2 + 0.5*zeta3 + (-87.6216536435678 + 1.*zeta2)*log(Q2/muF2) + (-175.2433072871356 + 2.*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==11){
		return (0.000244140625*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,10)*(-90.66403973677579 + 240.1705281329409*zeta2 + 0.5*zeta3 + (-177.45869068060486 + 1.*zeta2)*log(Q2/muF2) + (-354.91738136120966 + 2.*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==12){
		return (-0.0001220703125*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,11)*(-264.1988770872763 + 428.30790187031465*zeta2 + 0.5*zeta3 + (-352.23810386001804 + 1.*zeta2)*log(Q2/muF2) + (-704.4762077200361 + 2.*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==13){
		return (0.00006103515625*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,12)*(-688.6630058320999 + 772.6427503551631*zeta2 + 0.5*zeta3 + (-690.7639335858478 + 1.*zeta2)*log(Q2/muF2) + (-1381.5278671716956 + 2.*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==14){
		return (-0.000030517578125*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,13)*(-1661.9088653914168 + 1407.4784146908273*zeta2 + 0.5*zeta3 + (-1345.2021176240319 + 1.*zeta2)*log(Q2/muF2) + (-2690.4042352480633 + 2.*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==15){
		return (-1.3428851146193327e-17*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,14)*(4.325588989993079e15 - 2.9374621844384395e15*zeta2 - 5.68134567e11*zeta3 + (2.9653196642940245e15 - 1.136269134e12*zeta2)*log(Q2/muF2) + (5.930639328588049e15 - 2.272538268e12*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==16){
		return (2.238141857698888e-18*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,15)*(2.867658060532632e16 - 1.6299690914652118e16*zeta2 - 1.704403701e12*zeta3 + (1.7226349705365914e16 - 3.408807402e12*zeta2)*log(Q2/muF2) + (3.4452699410731828e16 - 6.817614804e12*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==17){
		return (-5.59535464424722e-19*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,16)*(1.2357278997739048e17 - 6.0656816989408584e16*zeta2 - 3.408807402e12*zeta3 + (6.667169839314207e16 - 6.817614804e12*zeta2)*log(Q2/muF2) + (1.3334339678628414e17 - 1.3635229608e13*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==18){
		return (9.680544367209724e-22*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,17)*(7.55545044232028e19 - 3.278415727429747e19*zeta2 - 9.85145339178e14*zeta3 + (3.728953042231441e19 - 1.970290678356e15*zeta2)*log(Q2/muF2) + (7.457906084462882e19 - 3.940581356712e15*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==19){
		return (-2.0974512795621068e-21*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,18)*(3.639559484346177e19 - 1.4211684308666585e19*zeta2 - 2.27341232118e14*zeta3 + (1.6661857537795895e19 - 4.54682464236e14*zeta2)*log(Q2/muF2) + (3.332371507559179e19 - 9.09364928472e14*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==20){
		return (2.2346593645451807e-25*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,19)*(3.529754959854911e23 - 1.2576871600411686e23*zeta2 - 1.066912402329774e18*zeta3 + (1.515113771617754e23 - 2.133824804659548e18*zeta2)*log(Q2/muF2) + (3.030227543235508e23 - 4.267649609319096e18*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}
	if(power==21){
		return (-1.8994604598634037e-24*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,20)*(4.258117266154896e22 - 1.3998211020194622e22*zeta2 - 6.2759553078222e16*zeta3 + (1.7283545239026034e22 - 1.25519106156444e17*zeta2)*log(Q2/muF2) + (3.456709047805207e22 - 2.51038212312888e17*zeta2)*log(1. - 1.*x)))/pow(M_PI,2);
	}

}

double DY_NNLO_CF_full(double x){
	return (pow(alphas_muR,2)*CF*(-CA/2. + CF)*(-30 + 56*x - 26*pow(x,2) + 4*(-7 + 6*x)*log(x) - (4*pow(-1 + x,2)*(9*pow(log(x),2) + 6*log(1 - x)*pow(log(x),2) + 2*pow(log(x),3) + 6*(3 + 2*log(x))*li2(1 - x) + 12*log(x)*li2(x) - 12*Li3(1 - x) - 12*Li3(x) + 12*zeta3))/3.))/(32.*pow(M_PI,2));
}

double DY_NNLO_CF_expansion(double x, int power){

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
		return (0.08333333333333333*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,3))/pow(M_PI,2);
	}
	if(power==5){
		return (-0.078125*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,4))/pow(M_PI,2);
	}
	if(power==6){
		return (0.0738425925925926*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,5))/pow(M_PI,2);
	}
	if(power==7){
		return (-0.07022569444444444*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,6))/pow(M_PI,2);
	}
	if(power==8){
		return (0.06710912698412698*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,7))/pow(M_PI,2);
	}
	if(power==9){
		return (-0.06438161375661376*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,8))/pow(M_PI,2);
	}
	if(power==10){
		return (0.06196469063816003*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,9))/pow(M_PI,2);
	}
	if(power==11){
		return (-0.05980101155045352*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,10))/pow(M_PI,2);
	}
	if(power==12){
		return (0.057847434171012214*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,11))/pow(M_PI,2);
	}
	if(power==13){
		return (-0.05607071744943769*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,12))/pow(M_PI,2);
	}
	if(power==14){
		return (0.05444473023323623*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,13))/pow(M_PI,2);
	}
	if(power==15){
		return (-0.05294857619741503*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,14))/pow(M_PI,2);
	}
	if(power==16){
		return (0.05156529665742008*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,15))/pow(M_PI,2);
	}
	if(power==17){
		return (-0.05028095020386583*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,16))/pow(M_PI,2);
	}
	if(power==18){
		return (0.04908394521543117*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,17))/pow(M_PI,2);
	}
	if(power==19){
		return (-0.047964546372080315*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,18))/pow(M_PI,2);
	}
	if(power==20){
		return (0.04691450355533483*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,19))/pow(M_PI,2);
	}
	if(power==21){
		return (-0.04592676852536536*pow(alphas_muR,2)*(CA - 2.*CF)*CF*pow(-1. + x,20))/pow(M_PI,2);
	}

}


////////////////////////////////////////////////////////////
///
/// gg channel
///
////////////////////////////////////////////////////////////
///// NNLO functions
double DY_NNLO_gg_full(double x){
	return (pow(alphas_muR,2)*(-32 - 66*x + 98*pow(x,2) + 4*(5 + 9*x - 12*pow(x,2))*zeta2 + 4*(-1 - 2*x + 2*pow(x,2))*zeta3 + 2*(7 + 60*x - 67*pow(x,2))*log(1 - x) + 16*(-1 - 2*x + 3*pow(x,2))*pow(log(1 - x),2) + (-23 - 64*x + 105*pow(x,2))*log(x) + 4*(3 + 10*x + 10*pow(x,2))*zeta2*log(x) + 4*(1 + 8*x - 4*pow(x,2))*log(1 - x)*log(x) - 2*(3 + 7*x + 4*pow(x,2))*pow(log(x),2) - (2*(3 + 8*x + 8*pow(x,2))*pow(log(x),3))/3. - 2*pow(log(Q2/muF2),2)*(2 + 4*x - 6*pow(x,2) + pow(1 + 2*x,2)*log(x)) + log(Q2/muF2)*(7 + 60*x - 67*pow(x,2) + 16*(-1 - 2*x + 3*pow(x,2))*log(1 - x) + 2*(1 + 8*x - 4*pow(x,2))*log(x) + 2*pow(1 + 2*x,2)*(log(x)*(-4*log(1 - x) + log(x)) - 4*li2(1 - x))) + 4*(-5 - 4*x + 14*pow(x,2))*li2(1 - x) + 8*(2 + 4*x + pow(x,2))*log(x)*li2(-x) + 8*(1 + x)*(log(x)*log(1 + x) + li2(-x)) - 4*pow(1 + 2*x,2)*(log(1 - x)*(2*log(1 - x) - log(x))*log(x) + (4*log(1 - x) + log(x))*li2(1 - x) - 4*Li3(1 - x)) + 8*(-1 - 2*x + pow(x,2))*Li3(-x) - 8*(1 + 10*x + 7*pow(x,2))*S12(1 - x) - 4*pow(1 + x,2)*(log(1 + x)*(2*zeta2 - 3*pow(log(x),2) + 2*log(x)*log(1 + x) + 4*li2(-x)) + 4*S12(-x)) + (pow(CA,2)*(-47 - 144*x + 191*pow(x,2) - 2*(6 + 38*x + 75*pow(x,2))*log(x) + 2*(-2 + 2*x + 25*pow(x,2))*pow(log(x),2) - 24*pow(-1 + x,2)*S12(1 - x) + 4*pow(1 + x,2)*(2*zeta2 + 12*zeta3 + 6*zeta2*log(1 + x) + 4*log(x)*log(1 + x) - 9*pow(log(x),2)*log(1 + x) + 6*log(x)*pow(log(1 + x),2) + (4 - 18*log(x) + 12*log(1 + x))*li2(-x) + 18*Li3(-x) + 12*S12(-x))))/(3.*(-1 + pow(CA,2)))))/(16.*pow(M_PI,2));
}

double DY_NNLO_gg_expansion(double x, int power){
	if(power==1){
		return 0;
	}
	if(power==2){
		return (0.125*pow(alphas_muR,2)*(1 - x)*(1.4202637326070946 + 1.*pow(log(Q2/muF2),2) - 8.*log(1 - x) + 4.*pow(log(1 - x),2) + log(Q2/muF2)*(-4. + 4.*log(1 - x))))/pow(M_PI,2);
	}
	if(power==3){
		return (-0.1875*pow(alphas_muR - alphas_muR*x,2)*(11.420263732607095 + 1.*pow(log(Q2/muF2),2) - 16.*log(1 - x) + 4.*pow(log(1 - x),2) + log(Q2/muF2)*(-8. + 4.*log(1 - x))))/pow(M_PI,2);
	}
	if(power==4){
		return (0.006944444444444444*pow(alphas_muR,2)*pow(1 - x,3)*(-461. + 455.*pow(CA,2) - 118.4352528130723*(-1. + pow(CA,2)) + 18.*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 72.*(-1. + pow(CA,2))*log(1 - x)*(-4.333333333333333 + 1.*log(1 - x)) + 72.*(-1. + pow(CA,2))*log(Q2/muF2)*(-2.1666666666666665 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==5){
		return (0.003472222222222222*pow(alphas_muR,2)*pow(-1 + x,4)*(9.*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 36.*(-1. + pow(CA,2))*log(Q2/muF2)*(0.08333333333333333 + 1.*log(1 - x)) + 2.*(33.60881320326807 - 36.60881320326807*pow(CA,2) + 18.*(-1. + pow(CA,2))*log(1 - x)*(0.16666666666666666 + 1.*log(1 - x)))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==6){
		return (4.6296296296296296e-6*pow(alphas_muR,2)*pow(1 - x,5)*(8623. - 11323.*pow(CA,2) - 23687.05056261446*(-1. + pow(CA,2)) + 3600.*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 14400.*(-1. + pow(CA,2))*log(Q2/muF2)*(1.1541666666666666 + 1.*log(1 - x)) + 14400.*(-1. + pow(CA,2))*log(1 - x)*(2.308333333333333 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==7){
		return (6.944444444444445e-6*pow(alphas_muR,2)*pow(-1 + x,6)*(-2951. + 1201.*pow(CA,2) - 11843.52528130723*(-1. + pow(CA,2)) + 1800.*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 7200.*(-1. + pow(CA,2))*log(Q2/muF2)*(1.5083333333333333 + 1.*log(1 - x)) + 7200.*(-1. + pow(CA,2))*log(1 - x)*(3.0166666666666666 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==8){
		return (4.049238743116294e-8*pow(alphas_muR,2)*pow(1 - x,7)*(-1.484437e6 + 1.170837e6*pow(CA,2) - 1.7409982163521626e6*(-1. + pow(CA,2)) + 264600.*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 1.0584e6*(-1. + pow(CA,2))*log(Q2/muF2)*(1.5876984126984126 + 1.*log(1 - x)) + 1.0584e6*(-1. + pow(CA,2))*log(1 - x)*(3.175396825396825 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==9){
		return (8.435914048158946e-9*pow(alphas_muR,2)*pow(-1 + x,8)*(-9.969413e6 + 8.406313e6*pow(CA,2) - 7.544325604192705e6*(-1. + pow(CA,2)) + 1.1466e6*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 4.5864e6*(-1. + pow(CA,2))*log(Q2/muF2)*(1.5917582417582417 + 1.*log(1 - x)) + 4.5864e6*(-1. + pow(CA,2))*log(1 - x)*(3.1835164835164833 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==10){
		return (1.6871828096317892e-8*pow(alphas_muR,2)*pow(1 - x,9)*(-5.863801e6 + 5.062161e6*pow(CA,2) - 3.4819964327043253e6*(-1. + pow(CA,2)) + 529200.*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 2.1168e6*(-1. + pow(CA,2))*log(Q2/muF2)*(1.582142857142857 + 1.*log(1 - x)) + 2.1168e6*(-1. + pow(CA,2))*log(1 - x)*(3.164285714285714 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==11){
		return (1.3122532963802805e-9*pow(alphas_muR,2)*pow(-1 + x,10)*(-8.2569722e7 + 7.2100067e7*pow(CA,2) - 4.178395719245191e7*(-1. + pow(CA,2)) + 6.3504e6*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 2.54016e7*(-1. + pow(CA,2))*log(Q2/muF2)*(1.5749007936507937 + 1.*log(1 - x)) + 2.54016e7*(-1. + pow(CA,2))*log(1 - x)*(3.1498015873015874 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==12){
		return (6.572768827349263e-13*pow(alphas_muR,2)*pow(1 - x,11)*(-1.74110550661e11 + 1.53027772021e11*pow(CA,2) - 7.836581171444354e10*(-1. + pow(CA,2)) + 1.19101752e10*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 4.76407008e10*(-1. + pow(CA,2))*log(Q2/muF2)*(1.573037401666434 + 1.*log(1 - x)) + 4.76407008e10*(-1. + pow(CA,2))*log(1 - x)*(3.146074803332868 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==13){
		return (1.4788729861535843e-12*pow(alphas_muR,2)*pow(-1 + x,12)*(-8.0049209985e10 + 7.0647801353e10*pow(CA,2) - 3.2863082331863422e10*(-1. + pow(CA,2)) + 4.9945896e9*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 1.99783584e10*(-1. + pow(CA,2))*log(Q2/muF2)*(1.5761391386391386 + 1.*log(1 - x)) + 1.99783584e10*(-1. + pow(CA,2))*log(1 - x)*(3.152278277278277 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==14){
		return (1.3462658044183744e-14*pow(alphas_muR,2)*pow(1 - x,13)*(-5.562311481447204e12 + 4.530214379411204e12*pow(CA,2) + 5.194373184e11*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 2.0777492736e12*(-1. + pow(CA,2))*log(Q2/muF2)*(1.5831102924852924 + 1.*log(1 - x)) + 2.0777492736e12*(-1. + pow(CA,2))*log(1 - x)*(3.166220584970585 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==15){
		return (3.525934249667171e-15*pow(alphas_muR,2)*pow(-1 + x,14)*(-2.233125975199249e13 + 1.840495308370849e13*pow(CA,2) + 1.8829602792e12*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 7.5318411168e12*(-1. + pow(CA,2))*log(Q2/muF2)*(1.592908863167484 + 1.*log(1 - x)) + 7.5318411168e12*(-1. + pow(CA,2))*log(1 - x)*(3.185817726334968 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==16){
		return (4.231121099600605e-15*pow(alphas_muR,2)*pow(1 - x,15)*(-1.9303996440811836e13 + 1.6051739384763836e13*pow(CA,2) + 1.4933822904e12*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 5.9735291616e12*(-1. + pow(CA,2))*log(Q2/muF2)*(1.6047029661160095 + 1.*log(1 - x)) + 5.9735291616e12*(-1. + pow(CA,2))*log(1 - x)*(3.209405932232019 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==17){
		return (3.4377858934254918e-15*pow(alphas_muR,2)*pow(-1 + x,16)*(-2.4412690136300938e13 + 2.0441301216516938e13*pow(CA,2) + 1.7531009496e12*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 7.0124037984e12*(-1. + pow(CA,2))*log(Q2/muF2)*(1.6178667012000345 + 1.*log(1 - x)) + 7.0124037984e12*(-1. + pow(CA,2))*log(1 - x)*(3.235733402400069 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==18){
		return (6.530836896391463e-18*pow(alphas_muR,2)*pow(1 - x,17)*(-1.3113421214013556e16 + 1.1042298481938016e16*pow(CA,2) + 8.819396369784e14*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 3.5277585479136e15*(-1. + pow(CA,2))*log(Q2/muF2)*(1.6319399033729447 + 1.*log(1 - x)) + 3.5277585479136e15*(-1. + pow(CA,2))*log(1 - x)*(3.2638798067458894 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==19){
		return (8.163546120489329e-18*pow(alphas_muR,2)*pow(-1 + x,18)*(-1.064990731734814e16 + 9.010232146994206e15*pow(CA,2) + 6.755282325792e14*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 2.7021129303168e15*(-1. + pow(CA,2))*log(Q2/muF2)*(1.6465872634990282 + 1.*log(1 - x)) + 2.7021129303168e15*(-1. + pow(CA,2))*log(1 - x)*(3.2931745269980564 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==20){
		return (4.76077919258745e-21*pow(alphas_muR,2)*pow(1 - x,19)*(-1.8465396653877895e19 + 1.5685446681067764e19*pow(CA,2) + 1.1109437078227488e18*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 4.443774831290995e18*(-1. + pow(CA,2))*log(Q2/muF2)*(1.661565214813629 + 1.*log(1 - x)) + 4.443774831290995e18*(-1. + pow(CA,2))*log(1 - x)*(3.323130429627258 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}
	if(power==21){
		return (1.3488874378997772e-21*pow(alphas_muR,2)*pow(-1 + x,20)*(-6.5693280641250935e19 + 5.59989012129939e19*pow(CA,2) + 3.7663701313990753e18*(-1. + pow(CA,2))*pow(log(Q2/muF2),2) + 1.5065480525596301e19*(-1. + pow(CA,2))*log(Q2/muF2)*(1.6766970404865198 + 1.*log(1 - x)) + 1.5065480525596301e19*(-1. + pow(CA,2))*log(1 - x)*(3.3533940809730396 + 1.*log(1 - x))))/((-1. + pow(CA,2))*pow(M_PI,2));
	}

}


////////////////////////////////////////////////////////////
///
/// qg channel
///
////////////////////////////////////////////////////////////
///// NNLO functions
double DY_NNLO_qg_full(double x){
	return (pow(alphas_muR,2)*TF*(beta0*log(muR2/muF2)*(1 + 6*x - 7*pow(x,2) + 2*(1 + 2*(-1 + x)*x)*(log(Q2/muF2) + 2*log(1 - x) - log(x))) + CA*(59.888888888888886 + 116/(27.*x) - (1226*x)/9. + (1837*pow(x,2))/27. + (2*pow(M_PI,2)*(15 - 8/x - 84*x + 107*pow(x,2)))/9. + 8*x*li2(1 - x) + (4*(33 + 16/x + 90*x + 44*pow(x,2))*li2(1 - x))/3. - 4*(15 + 34*x + 12*pow(x,2))*Li3(1 - x) - (16*pow(M_PI,2)*(1 - x + 2*pow(x,2))*log(1 - x))/3. + (2*(-210 + 88/x + 75*x + 74*pow(x,2))*log(1 - x))/9. + 8*(7 + 10*x + 5*pow(x,2))*li2(1 - x)*log(1 - x) + (4*(6 + 8/x + 63*x - 77*pow(x,2))*pow(log(1 - x),2))/3. + (26*(1 - 2*x + 2*pow(x,2))*pow(log(1 - x),3))/3. - 8*x*log(x) + (8*pow(M_PI,2)*x*(-5 + 2*x)*log(x))/3. - (2*(-354 + 12*x + 457*pow(x,2))*log(x))/9. + 8*(7 - 2*x)*x*li2(1 - x)*log(x) + 8*x*log(1 - x)*log(x) + 20*(1 - 2*x + 13*pow(x,2))*log(1 - x)*log(x) + 4*(1 + 22*x - 6*pow(x,2))*pow(log(1 - x),2)*log(x) - 4*x*pow(log(x),2) - (5 + (346*pow(x,2))/3.)*pow(log(x),2) + 4*(-3 - 14*x + 2*pow(x,2))*log(1 - x)*pow(log(x),2) + (2*(9 + 20*x)*pow(log(x),3))/3. + (2*pow(log(Q2/muF2),2)*(3 + 4/x + 24*x - 31*pow(x,2) + 6*(1 - 2*x + 2*pow(x,2))*log(1 - x) + 6*(1 + 4*x)*log(x)))/3. + 8*(1 + 5*x + 4*pow(x,2))*(li2(-x) + log(x)*log(1 + x)) - 4*(1 + 2*x + 2*pow(x,2))*(2*Li3(-x) - 4*Li3((1 - x)/(1 + x)) + 4*Li3((-1 + x)/(1 + x)) + 4*li2(-x)*(log(1 - x) - log(x)) + 4*log(1 - x)*log(x)*log(1 + x) - 3*pow(log(x),2)*log(1 + x)) + (2*log(Q2/muF2)*(-87 + 44/x - 12*x + 73*pow(x,2) - 12*pow(M_PI,2)*(1 - x + 2*pow(x,2)) + 36*(3 + 6*x + 2*pow(x,2))*li2(1 - x) + 6*(9 + 8/x + 54*x - 71*pow(x,2))*log(1 - x) + 54*(1 - 2*x + 2*pow(x,2))*pow(log(1 - x),2) + 18*(3 - 2*x + 28*pow(x,2))*log(x) + 36*(1 + 10*x - 2*pow(x,2))*log(1 - x)*log(x) - 36*(1 + 3*x)*pow(log(x),2) - 36*(1 + 2*x + 2*pow(x,2))*(li2(-x) + log(x)*log(1 + x))))/9. - 4*(1 + 4*x + 2*pow(x,2))*zeta3 + 4*(9 + 16*x + 4*pow(x,2))*(-2*Li3(x) + 2*li2(x)*log(x) + log(1 - x)*pow(log(x),2) + 2*zeta3)) + CF*(-90.5 - 12*(-1 + x) + 233*x - (305*pow(x,2))/2. + (2*pow(M_PI,2)*(5 + 2*x - 12*pow(x,2)))/3. + 8*(-3 + x)*li2(1 - x) + 2*(3 - 28*x + 40*pow(x,2))*li2(1 - x) + 4*(-1 + 2*x + 18*pow(x,2))*Li3(1 - x) + 24*(-1 + x)*log(1 - x) + 2*(38 - 147*x + 88*pow(x,2))*log(1 - x) - 4*(3 - 6*x + 26*pow(x,2))*li2(1 - x)*log(1 - x) - 2*(23 - 80*x + 63*pow(x,2))*pow(log(1 - x),2) + (28 - 44*x)*log(x) + 4*pow(M_PI,2)*(1 - 2*x + 4*pow(x,2))*log(x) - (59 - 245*x + 174*pow(x,2))*log(x) + 4*(1 - 2*x)*li2(1 - x)*log(x) + 8*(-3 + x)*log(1 - x)*log(x) + 4*(13 - 50*x + 48*pow(x,2))*log(1 - x)*log(x) - 6*(7 - 14*x + 22*pow(x,2))*pow(log(1 - x),2)*log(x) - 4*(-3 + x)*pow(log(x),2) - ((35 - 68*x + 4*pow(x,2))*pow(log(x),2))/2. + 8*(3 - 6*x + 10*pow(x,2))*log(1 - x)*pow(log(x),2) - ((17 - 34*x + 52*pow(x,2))*pow(log(x),3))/3. + 3*pow(log(Q2/muF2),2)*(-1 + 4*x + (4 - 8*x + 8*pow(x,2))*log(1 - x) + (-2 + 4*x - 8*pow(x,2))*log(x)) + log(Q2/muF2)*(24 - 68*x + 22*pow(x,2) - 48*pow(x,2)*li2(1 - x) - 4*(8 - 34*x + 23*pow(x,2))*log(1 - x) - (4*(1 - 2*x + 2*pow(x,2))*(pow(M_PI,2) - 27*pow(log(1 - x),2)))/3. + 2*(5 - 40*x + 46*pow(x,2))*log(x) - 8*(5 - 10*x + 16*pow(x,2))*log(1 - x)*log(x) + 8*(1 - 2*x + 4*pow(x,2))*pow(log(x),2)) - 16*(-1 + 2*x + 3*pow(x,2))*(li2(-x) + log(x)*log(1 + x)) - 2*(11 - 22*x + 34*pow(x,2))*(-2*Li3(x) + 2*li2(x)*log(x) + log(1 - x)*pow(log(x),2) + 2*zeta3) + (2*(1 - 2*x + 2*pow(x,2))*(48*Li3(-x) - 4*pow(M_PI,2)*log(1 - x) + 35*pow(log(1 - x),3) - 24*li2(-x)*log(x) + 150*zeta3))/3.)))/(16.*pow(M_PI,2));
}

double DY_NNLO_qg_expansion(double x, int power){
	if(power==1){
		return (0.25*pow(alphas_muR,2)*TF*(-0.3989715484202028*CA + 25.273883360576974*CF + (-8.369604401089358*CA - 17.079736267392907*CF + 1.*beta0*log(muR2/muF2))*log(1 - x) - 3.*CF*pow(log(1 - x),2) + (2.1666666666666665*CA + 5.833333333333333*CF)*pow(log(1 - x),3) + pow(log(Q2/muF2),2)*(2.25*CF + (1.*CA + 3.*CF)*log(1 - x)) + log(Q2/muF2)*(-3.934802200544679*CA - 8.789868133696453*CF + 0.5*beta0*log(muR2/muF2) + 3.*CF*log(1 - x) + (3.*CA + 9.*CF)*pow(log(1 - x),2))))/pow(M_PI,2);
	}
	if(power==2){
		return (-0.5*pow(alphas_muR,2)*TF*(1 - x)*(-1.989499013634391*CA + 8.813751494273427*CF + 2.880395598910642*CA*log(1 - x) + 18.670263732607093*CF*log(1 - x) - 8.*CA*pow(log(1 - x),2) - 22.75*CF*pow(log(1 - x),2) + 2.1666666666666665*CA*pow(log(1 - x),3) + 5.833333333333333*CF*pow(log(1 - x),3) + beta0*log(muR2/muF2)*(-1.25 + 0.5*log(Q2/muF2) + 1.*log(1 - x)) + pow(log(Q2/muF2),2)*(-1.*CA - 0.75*CF + (1.*CA + 3.*CF)*log(1 - x)) + log(Q2/muF2)*(1.065197799455321*CA + 2.460131866303547*CF + (-7.*CA - 17.*CF)*log(1 - x) + (3.*CA + 9.*CF)*pow(log(1 - x),2))))/pow(M_PI,2);
	}
	if(power==3){
		return (0.5*TF*pow(alphas_muR - alphas_muR*x,2)*(-5.093628129257088*CA - 32.48151223833367*CF + 15.630395598910644*CA*log(1 - x) + 64.7952637326071*CF*log(1 - x) - 10.75*CA*pow(log(1 - x),2) - 32.625*CF*pow(log(1 - x),2) + 2.1666666666666665*CA*pow(log(1 - x),3) + 5.833333333333333*CF*pow(log(1 - x),3) + beta0*log(muR2/muF2)*(-1.25 + 0.5*log(Q2/muF2) + 1.*log(1 - x)) + pow(log(Q2/muF2),2)*(-1.5*CA - 3.375*CF + (1.*CA + 3.*CF)*log(1 - x)) + log(Q2/muF2)*(5.565197799455322*CA + 24.585131866303545*CF + (-9.*CA - 28.*CF)*log(1 - x) + (3.*CA + 9.*CF)*pow(log(1 - x),2))))/pow(M_PI,2);
	}
	if(power==4){
		return (-0.25*pow(alphas_muR,2)*TF*pow(-1. + x,3)*(8.970157094843819*CA + 50.906541302009146*CF + 0.6666666666666666*beta0*log(muR2/muF2) + (1.*CA + 3.*CF)*pow(log(Q2/muF2),2) - 17.22222222222222*CA*log(1. - 1.*x) - 57.388888888888886*CF*log(1. - 1.*x) + 8.*CA*pow(log(1. - 1.*x),2) + 18.*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(-7.444444444444445*CA - 25.166666666666668*CF + (6.666666666666667*CA + 17.333333333333332*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==5){
		return (0.1875*pow(alphas_muR,2)*TF*pow(-1. + x,4)*(-0.9233169912733352*CA + 0.9723609374646732*CF + 0.3888888888888889*beta0*log(muR2/muF2) + (1.*CA + 1.5*CF)*pow(log(Q2/muF2),2) - 4.62962962962963*CA*log(1. - 1.*x) - 10.324074074074074*CF*log(1. - 1.*x) + 6.333333333333333*CA*pow(log(1. - 1.*x),2) + 9.5*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(-2.7037037037037037*CA - 5.055555555555555*CF + (5.555555555555555*CA + 9.11111111111111*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==6){
		return (-0.16666666666666669*pow(alphas_muR,2)*TF*pow(-1. + x,5)*(-3.4798358110258123*CA - 1.1813570538417006*CF + 0.275*beta0*log(muR2/muF2) + (1.*CA + 0.975*CF)*pow(log(Q2/muF2),2) + 0.29*CA*log(1. - 1.*x) - 0.5008333333333334*CF*log(1. - 1.*x) + 5.65*CA*pow(log(1. - 1.*x),2) + 6.375*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(-0.5466666666666666*CA - 0.87*CF + (5.1*CA + 6.1*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==7){
		return (0.15833333333333333*pow(alphas_muR,2)*TF*pow(-1. + x,6)*(-3.335237778733355*CA - 0.324546496269552*CF + 0.21052631578947367*beta0*log(muR2/muF2) + (1.*CA + 0.7105263157894737*CF)*pow(log(Q2/muF2),2) + 2.495614035087719*CA*log(1. - 1.*x) + 2.6289473684210525*CF*log(1. - 1.*x) + 5.2631578947368425*CA*pow(log(1. - 1.*x),2) + 4.7368421052631575*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(0.5087719298245614*CA + 0.48947368421052634*CF + (4.842105263157895*CA + 4.526315789473684*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==8){
		return (-0.15476190476190477*pow(alphas_muR,2)*TF*pow(-1. + x,7)*(-2.6293620636491237*CA + 0.562668852328387*CF + 0.16923076923076924*beta0*log(muR2/muF2) + (1.*CA + 0.5538461538461539*CF)*pow(log(Q2/muF2),2) + 3.624908424908425*CA*log(1. - 1.*x) + 3.754065934065934*CF*log(1. - 1.*x) + 5.015384615384615*CA*pow(log(1. - 1.*x),2) + 3.7384615384615385*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(1.0965567765567765*CA + 1.0059340659340659*CF + (4.676923076923077*CA + 3.5692307692307694*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==9){
		return (0.1532738095238095*pow(alphas_muR,2)*TF*pow(-1. + x,8)*(-1.8830033560432784*CA + 1.2124903115063144*CF + 0.1407766990291262*beta0*log(muR2/muF2) + (1.*CA + 0.45145631067961167*CF)*pow(log(Q2/muF2),2) + 4.273370319001387*CA*log(1. - 1.*x) + 4.137274618585298*CF*log(1. - 1.*x) + 4.844660194174757*CA*pow(log(1. - 1.*x),2) + 3.0728155339805827*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(1.4630605640314378*CA + 1.2079519186315302*CF + (4.563106796116505*CA + 2.9320388349514563*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==10){
		return (-0.1527777777777778*pow(alphas_muR,2)*TF*pow(-1. + x,9)*(-1.213749019110908*CA + 1.6623667404752127*CF + 0.12012987012987013*beta0*log(muR2/muF2) + (1.*CA + 0.37987012987012986*CF)*pow(log(Q2/muF2),2) + 4.683663162234591*CA*log(1. - 1.*x) + 4.213345186559472*CF*log(1. - 1.*x) + 4.720779220779221*CA*pow(log(1. - 1.*x),2) + 2.6006493506493507*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(1.7133065347351062*CA + 1.2769635126777983*CF + (4.48051948051948*CA + 2.4805194805194803*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==11){
		return (0.1527777777777778*pow(alphas_muR,2)*TF*pow(-1. + x,10)*(-0.6364863240390094*CA + 1.9709467346150398*CF + 0.10454545454545454*beta0*log(muR2/muF2) + (1.*CA + 0.32727272727272727*CF)*pow(log(Q2/muF2),2) + 4.965003607503608*CA*log(1. - 1.*x) + 4.155169552669553*CF*log(1. - 1.*x) + 4.627272727272727*CA*pow(log(1. - 1.*x),2) + 2.25*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(1.8966089466089466*CA + 1.2854329004329004*CF + (4.418181818181818*CA + 2.1454545454545455*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==12){
		return (-0.15303030303030302*pow(alphas_muR,2)*TF*pow(-1. + x,11)*(-0.14112774897492836*CA + 2.18253418023331*CF + 0.0924092409240924*beta0*log(muR2/muF2) + (1.*CA + 0.2871287128712871*CF)*pow(log(Q2/muF2),2) + 5.170930902614071*CA*log(1. - 1.*x) + 4.038303354144938*CF*log(1. - 1.*x) + 4.554455445544554*CA*pow(log(1. - 1.*x),2) + 1.9801980198019802*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(2.038234061501388*CA + 1.2654079693683653*CF + (4.36963696369637*CA + 1.8877887788778878*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==13){
		return (0.1534090909090909*pow(alphas_muR,2)*TF*pow(-1. + x,12)*(0.28643512463431964*CA + 2.3272672507860195*CF + 0.08271604938271605*beta0*log(muR2/muF2) + (1.*CA + 0.25555555555555554*CF)*pow(log(Q2/muF2),2) + 5.3296999982185165*CA*log(1. - 1.*x) + 3.8978923271515864*CF*log(1. - 1.*x) + 4.496296296296296*CA*pow(log(1. - 1.*x),2) + 1.7666666666666666*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(2.1521564854898187*CA + 1.2321094543316766*CF + (4.330864197530865*CA + 1.6839506172839507*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==14){
		return (-2.243776340697291e-15*pow(alphas_muR,2)*TF*pow(-1. + x,13)*(4.517940892497431e13*CA + 1.662931937653972e14*CF + 5.1294435192e12*beta0*log(muR2/muF2) + (6.85657260288e13*CA + 1.57779085464e13*CF)*pow(log(Q2/muF2),2) + 3.741778132092e14*CA*log(1. - 1.*x) + 2.5718326786152e14*CF*log(1. - 1.*x) + 3.050395652304e14*CA*pow(log(1. - 1.*x),2) + 1.092766258584e14*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(1.5404238217368e14*CA + 8.180770630032e13*CF + (2.94780678192e14*CA + 1.041471823392e14*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==15){
		return (1.0577802749001514e-14*pow(alphas_muR,2)*TF*pow(-1. + x,14)*(1.4392801329777062e13*CA + 3.6326743229788555e13*CF + 9.955881936e11*beta0*log(muR2/muF2) + (1.45875313584e13*CA + 3.0516942456e12*CF)*pow(log(Q2/muF2),2) + 8.114962767912e13*CA*log(1. - 1.*x) + 5.25970080636e13*CF*log(1. - 1.*x) + 6.43236545952e13*CA*pow(log(1. - 1.*x),2) + 2.11670707248e13*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(3.39430699608e13*CA + 1.681015303968e13*CF + (6.2332478208e13*CA + 2.01714825312e13*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==16){
		return (-1.4103736998668685e-15*pow(alphas_muR,2)*TF*pow(-1. + x,15)*(1.4020689634011588e14*CA + 2.7777555813196706e14*CF + 6.8825444688e12*beta0*log(muR2/muF2) + (1.09731133512e14*CA + 2.10372113952e13*CF)*pow(log(Q2/muF2),2) + 6.2029105394256e14*CA*log(1. - 1.*x) + 3.8032670147472e14*CF*log(1. - 1.*x) + 4.802198008608e14*CA*pow(log(1. - 1.*x),2) + 1.460917458e14*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(2.629394199432e14*CA + 1.2200372232048e14*CF + (4.664547119232e14*CA + 1.392092013312e14*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==17){
		return (2.864821577854577e-16*pow(alphas_muR,2)*TF*pow(-1. + x,16)*(8.33510852975707e14*CA + 1.3843758557336062e15*CF + 3.14259577632e13*beta0*log(muR2/muF2) + (5.417731230912e14*CA + 9.58361852448e13*CF)*pow(log(Q2/muF2),2) + 3.10476428694432e15*CA*log(1. - 1.*x) + 1.80613985079528e15*CF*log(1. - 1.*x) + 2.355648238944e15*CA*pow(log(1. - 1.*x),2) + 6.66178360848e14*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(1.33117061261184e15*CA + 5.8110422321952e14*CF + (2.2927963234176e15*CA + 6.347524030848e14*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==18){
		return (-8.163546120489329e-19*pow(alphas_muR,2)*TF*pow(-1. + x,17)*(3.381870301762885e17*CA + 4.892802681469229e17*CF + 1.02830408737056e16*beta0*log(muR2/muF2) + (1.90649078972352e17*CA + 3.12994747761696e16*CF)*pow(log(Q2/muF2),2) + 1.1056413475265559e18*CA*log(1. - 1.*x) + 6.118889226495645e17*CF*log(1. - 1.*x) + 8.242945611316416e17*CA*pow(log(1. - 1.*x),2) + 2.177452669680288e17*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(4.787331348418819e17*CA + 1.973462315256432e17*CF + (8.037284793842304e17*CA + 2.074622260943232e17*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==19){
		return (9.524137140570883e-18*pow(alphas_muR,2)*TF*pow(-1. + x,18)*(3.256858913170483e16*CA + 4.207614876632963e16*CF + 8.256456175968e14*beta0*log(muR2/muF2) + (1.63842403076352e16*CA + 2.5091048638656e15*CF)*pow(log(Q2/muF2),2) + 9.602013484525176e16*CA*log(1. - 1.*x) + 5.067987713752392e16*CF*log(1. - 1.*x) + 7.04908349361216e16*CA*pow(log(1. - 1.*x),2) + 1.74672300138336e16*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(4.193492659795296e16*CA + 1.637836795953072e16*CF + (6.8839543700928e16*CA + 1.66415843962368e16*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==20){
		return (-1.5869263975291498e-21*pow(alphas_muR,2)*TF*pow(-1. + x,19)*(2.1523025217471578e20*CA + 2.526261153899053e20*CF + 4.66054433525641e18*beta0*log(muR2/muF2) + (9.857593192827219e19*CA + 1.414421013374329e19*CF)*pow(log(Q2/muF2),2) + 5.831360952415876e20*CA*log(1. - 1.*x) + 2.9420491361135434e20*CF*log(1. - 1.*x) + 4.222669937246272e20*CA*pow(log(1. - 1.*x),2) + 9.852173955228084e19*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(2.5661156924165882e20*CA + 9.524219315406188e19*CF + (4.1294590505411445e20*CA + 9.386119521702443e19*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
	if(power==21){
		return (2.6977748757995544e-21*pow(alphas_muR,2)*TF*pow(-1. + x,20)*(1.3735126513108405e20*CA + 1.4833590242307667e20*CF + 2.5876859535871345e18*beta0*log(muR2/muF2) + (5.812132325072674e19*CA + 7.844346424748433e18*CF)*pow(log(Q2/muF2),2) + 3.4673025266897755e20*CA*log(1. - 1.*x) + 1.6756776488602645e20*CF*log(1. - 1.*x) + 2.4801140872442977e20*CA*pow(log(1. - 1.*x),2) + 5.466655928127794e19*CF*pow(log(1. - 1.*x),2) + log(Q2/muF2)*(1.5361173195609712e20*CA + 5.4326130954681025e19*CF + (2.428360368172555e20*CA + 5.207887332769081e19*CF)*log(1. - 1.*x))))/pow(M_PI,2);
	}
}
