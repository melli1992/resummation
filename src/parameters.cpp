#include <bits/stdc++.h>
#include <iostream> //for stirng
#include <fstream>
#include <sstream>
#include <unistd.h> //for getting options
#include <gsl/gsl_math.h>
#include <complex>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <stdlib.h> // for exit()
#include <getopt.h> //for getting options
#include <vector>
#include "LHAPDF/LHAPDF.h"
#include "cuba.h"
#include "cuba_integration.h"
#include "SCET_functions.h"
#include "fit_coefficients.h"
#include "k_factors_dy.h"
#include "k_factors_higgs.h"
#include "cheb_pdf.h"

using namespace std;

//////////////////////////////////////////////////////////////////
///
/// parameters.cpp contains the functions to read in the parameters
/// that are used in the code and unchanged
/// contains stuff like invariant mass, com energy, tau, color
/// factors and the electromagnetic coupling constant
///
//////////////////////////////////////////////////////////////////


complex<double> I(0,1);
double CMP(2.5), phiMP(3./4.*M_PI);
int boundary = 0;

double pT(500.);
double s34in(4.*175.*175.);
double pT2(pow(pT,2));
double Q(500.);
double Q2(pow(Q,2.));
double S(13000.);
double S2(pow(S,2.));
double tau(Q2/S2), s2(1.);
double muF(Q);
double muR(Q);
double muF2(muF*muF);
double muR2(muR*muR);
double muh(Q);
double mus(Q);
double mut(Q);
//(6.62607015*10^(-34)*299792458/(2*pi))^2/(1.602176634*10^(−19))^2
// 3.8937937e
double pbunits(0.38937966*pow(10.,9.));
double fbunits(0.38937966*pow(10.,12.));

//double pbunits(0.389379323*pow(10.,9.)); //check this again // Jort!
//double fbunits(0.389379323*pow(10.,12.));
double quarkmasses[2] = {4.18,173.0};
double mt(172.5);//mt(173.0);
double mt2(pow(mt,2));
double mb(4.7);// should it be in MSbar?
double mb2(pow(mb,2));
double mZ(91.188000000000002);
double mZ2(pow(mZ,2));
double mW(80.419002445756163);
double mW2(pow(mW,2));
double mH(125.);//mH(125.09);
double mH2(pow(mH,2));
double GammaH(4.07E-3);
double Lt(log(pow(mH,2)/pow(mt,2)));
double M_gammaE(0.57721566490153286060651209008240243104215933593992);
double INCEULER(1.);
double zeta2(1.6449340668482264);
double zeta3(1.2020569031595942);
double zeta5(1.03692775514337);

// couplings needed for DY
double CF(4./3.);
double CA(3.);
double TF(1./2.);
//double alphaEM(1./127.940);
double alphaEM(1./132.5070);
double GF(1.1663790E-5); //GeV-2
double alphas_muF(0);
double alphas_muR(0);
double alphas_Q(0);
double LambdaQCD(0.208364837218848);
double nF(5); //4 for Richard!

double b0((11*CA/3.-2*nF/3.)/(4.*M_PI));
double b1((17.*pow(CA,2)-10.*CA*TF*nF-6.*CF*TF*nF)/(24.*pow(M_PI,2)));
double b2(1./(pow(4.*M_PI,3))*(2857./54.*pow(CA,3)-1415./27.*TF*pow(CA,2)*nF-205./9*TF*CA*CF*nF+2.*TF*pow(CF,2)*nF+79./54.*CA*pow(nF,2)+11./9.*CF*pow(nF,2)));
double b3((149753./6. + (1093.*nF*nF*nF)/729. + 3564*zeta3+ nF*nF*(50065./162. + (6472.*zeta3)/81.) - nF*(1078361./162. + (6508.*zeta3)/27.))/256./(M_PI*M_PI*M_PI*M_PI));
double beta0(4.*M_PI*b0);
double beta1(pow(4.*M_PI,2)*b1);
double beta2(pow(4.*M_PI,3)*b2);

// Nucl Phys B.410 (1993) 280-324 Frixione
// couplings needed for W+W-
//double sinw(sqrt(0.2223));23126
double sinw(sqrt(1.-mW2/mZ2));
double cosw(sqrt(1.-pow(sinw,2)));
double e(sqrt(4.*M_PI*alphaEM));
double ez(e*cosw/sinw);
double guL(e/(2.*sinw*cosw)*(1./2.-2./3.*pow(sinw,2)));
double gdL(e/(2.*sinw*cosw)*(-1./2.+1./3.*pow(sinw,2)));
double guR(-e/(2.*sinw*cosw)*(2./3.*pow(sinw,2)));
double gdR(e/(2.*sinw*cosw)*(1./3.*pow(sinw,2)));
double gVu(e/(2.*sinw*cosw)*(1./2.-4./3.*pow(sinw,2)));
double gAu(e/(2.*sinw*cosw)*(1./2.));
double gVd(e/(2.*sinw*cosw)*(-1./2.+2./3.*pow(sinw,2)));
double gAd(e/(2.*sinw*cosw)*(-1./2.));
double ctt = pow(e,4)/(16.*pow(sinw,4));
double v = 246.21845810181637;

// BSM stuff
double MWR(0.);
double MWR2(0.);
double gR2(4.*M_PI*alphaEM/pow(sinw,2));

//switches

double ISLL(1);
double ISNLL(1);
double ISNNLL(1);
double ISNNNLL(1);
double ISNLP(1);
bool INCSQRTZ = false,highscale=false,setdym = false, SCET=false, BN = false, DY=true, ttH=false,higgs=false, hh=false, WW=false, ZZ=false, full=true, diff=false, LO=true, NLO=true, NNLO=true, RES=true, realPDF=false, fitPDF=true, chebPDF = false, SUSY = false, INCHARD=false;
bool expansion = true;

// anomalous dimensions for dQCD
double A1q(CF); // 1405.4827 eq. 12
double A1g(CA);// 1405.4827 eq. 12
double A2q(CF/2.*(CA*(67./18.-pow(M_PI,2)/6.)-10./9.*TF*nF));// 1405.4827 eq. 12
double A2g(CA/2.*(CA*(67./18.-pow(M_PI,2)/6.)-10./9.*TF*nF));// 1405.4827 eq. 12
double A3q(CF*((245./96.-67./216.*pow(M_PI,2)+11./720.*pow(M_PI,4)+11./24.*zeta3)*pow(CA,2)+(-209./432.+5./108.*pow(M_PI,2)-7./12.*zeta3)*CA*nF+(-55./96.+1./2.*zeta3)*CF*nF-1./108.*pow(nF,2)));// 1405.4827 eq. 12
double A3g(CA*((245./96.-67./216.*pow(M_PI,2)+11./720.*pow(M_PI,4)+11./24.*zeta3)*pow(CA,2)+(-209./432.+5./108.*pow(M_PI,2)-7./12.*zeta3)*CA*nF+(-55./96.+1./2.*zeta3)*CF*nF-1./108.*pow(nF,2)));// 1405.4827 eq. 12
//pade for nF = 5
double A4q(0.0622776); // = A4/Pi/Pi/Pi/Pi
double A4g(0.0622776/CF*CA);


double C2qgamma(pow(CF,2));  // 0807.4412 eq 35
double Dbar2q(3./4.*pow(CF,2)-11./12.*CA*CF+1./6.*nF*CF); // 0807.4412 eq 35
double B1q(-3./4*CF);
double B1g(-M_PI*b0);
double D1DY(0.);
double D2DY((-101./27.+11./3.*zeta2+7./2.*zeta3)*CA*CF+(14./27.-2./3.*zeta2)*nF*CF); //hep-ph/0508284 eq. 3.17
double D3DY(CF*pow(CA,2)*(-297029./23328.+6139./324.*zeta2-187./60*pow(zeta2,2)+2509./108.*zeta3-11./6.*zeta2*zeta3-6.*zeta5)
+ nF*CA*CF*(31313./11664.-1837./324.*zeta2+23./30.*pow(zeta2,2)-155./36.*zeta3)
+ nF*pow(CF,2)*(1711./864.-1./2.*zeta2-1./5.*pow(zeta2,2)-19./18.*zeta3)
+ pow(nF,2)*CF*(-58./729.+10./27.*zeta2+5./27.*zeta3)); //hep-ph/0508284 eq 4.7
double D1higgs(0.);
double D2higgs((-101./27.+11./3.*zeta2+7./2.*zeta3)*CA*CA+(14./27.-2./3.*zeta2)*nF*CA); //hep-ph/0306211 eq. 37 //hep-ph/0508265 eq. 36 (pulled out 1/16=1/4^2)
double D3higgs(pow(CA,3)*(-297029./23328.+6139./324.*zeta2-187./60*pow(zeta2,2)+2509./108.*zeta3-11./6.*zeta2*zeta3-6.*zeta5)
+ nF*pow(CA,2)*(31313./11664.-1837./324.*zeta2+23./30.*pow(zeta2,2)-155./36.*zeta3)
+ nF*CF*CA*(1711./864.-1./2.*zeta2+1./5.*pow(zeta2,2)-19./18.*zeta3)
+ pow(nF,2)*CA*(-58./729.+10./27.*zeta2+5./27.*zeta3)); //hep-ph/0508265 eq 35 (note that they pulled out a factor of 1/64=1/4^3)

// SCET anomalous dimensions
double GammaC0q(4.*CF);
double GammaC1q(4.*CF*((67./9.-pow(M_PI,2)/3.)*CA-20./9.*TF*nF));
double GammaC2q(4.*CF*(pow(CA,2)*(245./6.-134.*pow(M_PI,2)/27.+11.*pow(M_PI,4)/45.+22./3.*zeta3)
                      +CA*TF*nF*(-418./27.+40.*pow(M_PI,2)/27.-56./3.*zeta3)
                      +CF*TF*nF*(-55./3.+16.*zeta3)-16./27.*pow(TF*nF,2)));
double GammaC0g(4.*CA);
double GammaC1g(4.*CA*((67./9.-pow(M_PI,2)/3.)*CA-20./9.*TF*nF));
double GammaC2g(4.*CA*(pow(CA,2)*(245./6.-134.*pow(M_PI,2)/27.+11.*pow(M_PI,4)/45.+22./3.*zeta3)
                      +CA*TF*nF*(-418./27.+40.*pow(M_PI,2)/27.-56./3.*zeta3)
                      +CF*TF*nF*(-55./3.+16.*zeta3)-16./27.*pow(TF*nF,2)));
vector<double> GammaCq = {GammaC0q,GammaC1q,GammaC2q};
vector<double> GammaCg = {GammaC0g,GammaC1g,GammaC2g};
vector<double> GammaC = {GammaC0g/CA,GammaC1g/CA,GammaC2g/CA};

double gammaV0q(-6.*CF);
double gammaV1q(CF*CF*(-3.+4.*pow(M_PI,2)-48*zeta3)+CF*CA*(-961./27.-11./3.*pow(M_PI,2)+52.*zeta3)+CF*TF*nF*(260./27.+4.*pow(M_PI,2)/3.));
double gammaV2q(0.);
double gammaS0(0.);
double gammaS1(CA*CA*(-160./27.+11.*pow(M_PI,2)/9.+4.*zeta3)+CA*TF*nF*(-208./27.-4.*pow(M_PI,2)/9.)-8.*CF*TF*nF);
vector<double> gammaS = {gammaS0,gammaS1};
vector<double> gammaV = {gammaV0q,gammaV1q};

double gammaphi0q(3.*CF);
double gammaphi1q(CF*CF*(3./2.-2.*pow(M_PI,2)+24.*zeta3)+CF*CA*(17./6.+22./9.*pow(M_PI,2)-12.*zeta3)-CF*TF*nF*(2./3.+8.*pow(M_PI,2)/9.));
double gammaB0(beta0);
double gammaB1(4.*CA*CA*(8./3.+3.*zeta3)-16./3.*CA*TF*nF-4.*CF*TF*nF);
vector<double> gammaB = {gammaB0,gammaB1};
vector<double> gammaphi = {gammaphi0q,gammaphi1q};

//double gammaphi2(pow(CF,3)*(29./2.+3.*pow(M_PI,4)+68));

//PDF declarations
string setname("PDF4LHC15_nnlo_100");
std::vector<LHAPDF::PDF*> pdfs;
double xmin_pdfs(0.), xmax_pdfs(0.); //min x, max x and alphas
int use_member(0);
double s1(0.), sgg(0.), sqqbar(0.);
double  *fitcheb_coeff_gg, *fitcheb_coeff_qqbar, *fitcheb_coeff_qqbarU, *fitcheb_coeff_qqbarD;
std::unordered_map<double, std::vector<std::vector<double>>> fitcoeff;

///////////////////////////////////////
/// update defaults of the programm
///////////////////////////////////////
void update_defaults(bool printout , bool pdfset){
    Q2 = pow(Q,2);
    muF2 = muF*muF;
    muR2 = muR*muR;
    tau = Q2/S2;
	  Lt = log(muR2/mt2);
    if(pdfset){
		LHAPDF::PDFSet setk(setname);
		int nmem(0.); //number of members
		vector<int> pids; //number of flavors, span from -5 to 5 with 0 = 21 gluon
		nmem = setk.size()-1;
		if(use_member>nmem){cout << "PDF member not there, using default value of 0" << endl; use_member = 0;}
		pdfs = setk.mkPDFs();
		pids = pdfs[use_member]->flavors();
		xmin_pdfs = pdfs[use_member]->xMin();
		xmax_pdfs = pdfs[use_member]->xMax();
	}
	if(fitPDF){
    cout << "Using fitted PDFs" << endl;
		if (setname=="PDF4LHC15_nnlo_100"){fitcoeff = fitcoeff_PDF4LHC15_nnlo_100;}
		else {cout << "PDFset " << setname << " not implemented, using default" << endl; fitcoeff = fitcoeff_PDF4LHC15_nnlo_100;}
		if (fitcoeff.find(muF) == fitcoeff.end())
		{  cout << "Scale " << muF << " not present in the fitted coefficients!" << endl;
				 cout << "Exiting program!" << endl;
				exit(0);}
	}
	if(chebPDF){
  	double *coeff   = new double[11];
    if(higgs || hh){ string lumchan = "gg"; ApproxLuminosity(coeff, Q2/S2, 10, lumchan); fitcheb_coeff_gg = coeff;}
    if(DY){ string lumchan = "qqbar"; ApproxLuminosity(coeff, Q2/S2, 10, lumchan); fitcheb_coeff_qqbar = coeff;}
    if(WW|| ZZ){ string lumchan = "qqbarU"; ApproxLuminosity(coeff, Q2/S2, 10, lumchan); fitcheb_coeff_qqbarU = coeff;
                  double *coeff2   = new double[11];
                  lumchan = "qqbarD"; ApproxLuminosity(coeff2, Q2/S2, 10, lumchan); fitcheb_coeff_qqbarD = coeff2;}
    if(ttH){string lumchan = "gg"; ApproxLuminosity(coeff, Q2/S2, 10, lumchan); fitcheb_coeff_gg = coeff;
                   lumchan = "qqbarH"; ApproxLuminosity(coeff, Q2/S2, 10, lumchan); fitcheb_coeff_qqbar = coeff;}
		}
  tau = Q2/S2;
	alphas_muF = pdfs[use_member]->alphasQ(muF);
	alphas_Q = pdfs[use_member]->alphasQ(Q);
	alphas_muR = pdfs[use_member]->alphasQ(muR);
	if(printout){
		cout << "=========================================" << endl;
		cout << "PDFset: 				" << setname << endl
			<< "PDFmember: 				" << use_member << endl
			<< "Center of Mass energy [TeV]: 		" << S/1000. << endl
			<< "Momentum [GeV]: 			" << Q << endl
			<< "Renormalization scale [GeV]: 		" << muR << endl
			<< "Factorization scale [GeV]:		" << muF << endl;
		if(SCET && setdym){
			solveLambdaQCD();
			muh = Q; mut = sqrt(mt2);
			if(higgs) {
        cout << "Calculating Higgs soft scale" << endl;
        vector<results_c> mus_scale = call_set_scale("gg");
        s1 = mus_scale[0].res;
      }
			else if(DY) {
        cout << "Calculating DY soft scale" << endl;
        vector<results_c> mus_scale = call_set_scale("qqbar");
        s1 = mus_scale[0].res;
      }
      mus = Q/(s1*exp(M_gammaE));
			cout << "LambdaQCD [GeV]:			" << LambdaQCD << endl
				 << "Soft scale [GeV]: 			" << mus << endl
   				 << "s1(tau) : 		          " << s1 << endl
			     << "Hard scale [GeV]:			" << muh << endl
				 << "alphas_mZ (grid):		 	" << pdfs[use_member]->alphasQ(sqrt(mZ2)) << endl
				 << "alphas_mZ (own):		 	" << falphasQ2(mZ2) << endl
			     << "alphas_muR (grid): 	 		" << alphas_muR << endl
				 << "alphas_muR (own):		 	" << falphasQ2(muR2) << endl;
		}
		else{
			cout << "alphas_mZ:			 	" << pdfs[use_member]->alphasQ(sqrt(mZ2)) << endl;
			cout << "alphas_muR: 		 		" << alphas_muR << endl;
		}
    if(ttH && setdym) {
      cout << "Calculating ttH soft scale" << endl;
      //vector<results_c> mus_scale = call_set_scale("ggH");
      //sgg = mus_scale[0].res; cout << "sgg = " << sgg << endl;
      //mus_scale = call_set_scale("qqbarH");
      //sqqbar = mus_scale[0].res; cout << "sqqbar = " << sqqbar << endl;
    }
		if(higgs) cout << "Higgs channel, tau = " << tau << ", prefactor = " << higgs_LO_factor() << endl;
		if(DY) cout << "DY channel, tau = " << tau << ", prefactor = " << DY_LO_factor() << endl;
		if(fitPDF){cout << "Using the fitted PDFs" << endl; realPDF = false;}
		else{ cout << "Using the real PDFs" << endl; realPDF = true;}
		cout << "INCEULER (resumming gammaE) = " << INCEULER << endl;
		cout << "=========================================" << endl;
	}
}
