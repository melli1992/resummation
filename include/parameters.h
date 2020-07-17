#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <stdlib.h> // for exit()
#include <getopt.h> //for getting options
#include <vector>
#include "LHAPDF/LHAPDF.h"


#ifndef PARAM_H //need this otherwise compiler error that things were predefined (called a guard)
#define PARAM_H

extern std::complex<double> I;
extern double CMP, phiMP;
extern int boundary;

extern double tau, s2, s34in;
extern double S;
extern double S2;
extern double pT;
extern double pT2;
extern double Q;
extern double Q2;
extern double muF;
extern double muF2;
extern double muR;
extern double muR2;
// SCET scales
extern double muh;
extern double mus;
extern double mut;

extern double pbunits;
extern double fbunits;

extern double quarkmasses[2];
extern double mt;
extern double mt2;
extern double mb;
extern double mb2;
extern double mZ;
extern double mZ2;
extern double mW;
extern double mW2;
extern double mH;
extern double mH2;
extern double GammaH;
extern double Lt;
extern double v;

extern double M_gammaE;
extern double zeta2;
extern double zeta3;
extern double zeta5;

extern double CF;
extern double CA;
extern double TF;
extern double alphas_muF;
extern double alphas_muR;
extern double alphas_Q;
extern double alphaEM;
extern double LambdaQCD;
extern double GF;
extern double nF;

extern double b0;
extern double b1;
extern double b2;
extern double b3;
extern double beta0;
extern double beta1;
extern double beta2;
//for W+W-, ZZ
extern double sinw;
extern double cosw;
extern double e;
extern double ez;
extern double guL;
extern double guR;
extern double gdL;
extern double gdR;
extern double gVu;
extern double gAu;
extern double gVd;
extern double gAd;
extern double ctt;

// BSM stuff
extern double MWR;
extern double MWR2;
extern double gR2;

//pQCD flags and parameters
extern double INCEULER;
extern double ISLL;
extern double ISNLL;
extern double ISNNLL;
extern double ISNNNLL;
extern double ISNLP;

extern double A1q;
extern double A1g;
extern double A2q;
extern double A2g;
extern double A3q;
extern double A3g;
extern double A4q;
extern double A4g;
extern double B1q;
extern double B1g;
extern double C2qgamma;
extern double Dbar2q;
extern double D1DY;
extern double D2DY;
extern double D3DY;
extern double D1higgs;
extern double D2higgs;
extern double D3higgs;

//double SCET anomalous dimensionsdouble GammaC0q(4.*CF);
extern double GammaC1q;
extern double GammaC2q;
extern double GammaC0g;
extern double GammaC1g;
extern double GammaC2g;
extern std::vector<double> GammaCq;
extern std::vector<double> GammaCg;

extern double gammaV0q;
extern double gammaV1q;
extern double gammaV2q;
extern double gammaS0;
extern double gammaS1;
extern std::vector<double> gammaS;
extern std::vector<double> gammaV;

extern double gammaphi0q;
extern double gammaphi1q;
extern double gammaB0;
extern double gammaB1;
extern std::vector<double> gammaB;
extern std::vector<double> gammaphi;



extern std::string setname;
extern std::vector<LHAPDF::PDF*> pdfs; //pdf vector
extern double xmin_pdfs, xmax_pdfs; //min x, max x
extern int use_member; //the member that one needs to use
extern double *fitcheb_coeff_gg, *fitcheb_coeff_qqbar, *fitcheb_coeff_qqbarU, *fitcheb_coeff_qqbarD;
extern double s1, sgg, sqqbar;
struct lumni_params {double z; double pT; double xT; double epeta; double emeta; int power; int flavor; int coefficient;};

extern bool INCSQRTZ,SCET, BN, DY, higgs, hh, WW, ZZ, ttH,diff, full, PF, LO, NLO, NNLO, RES, realPDF, fitPDF, chebPDF, SUSY, setdym,highscale, INCHARD;
extern bool expansion;

void update_defaults(bool printout = true , bool pdfset = true);

extern std::unordered_map<double, std::vector<std::vector<double>>> fitcoeff;
#endif
