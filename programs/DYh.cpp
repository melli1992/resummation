#include <cmath>
#include <iostream>
#include <fstream>
#include "cuba.h"
#include "cuba_integration.h"
#include <sstream>
#include "deriv_pdf.h"
#include "mellin_pdf.h"
#include "monte_carlo.h"
#include <chrono>
#include <gsl/gsl_sf_dilog.h>
#include "k_factors_dy.h"
#include "polygamma.h"
#include "resum_functions.h"
#include "k_factors_higgs.h"
#include "k_factors_dihiggs.h"
#include "k_factors_nnlo_dy.h"
#include "k_factors_nnlo_higgs.h"
#include "parameters.h"
#include "LHAPDF/LHAPDF.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <string.h>
#include <sstream>
#include "inout.h"
#include "cheb_pdf.h"

using namespace std;


////////////////////////////////////////////////////////////////
///
/// resums DY and higgs and calculates fixed-order expansions
///
/// inputs and flags via DYh.cfg
///
////////////////////////////////////////////////////////////////


string to_string(double Q){
	ostringstream q_to_str;
	q_to_str << Q;
	return q_to_str.str();
}


string to_string2(string Q){
	ostringstream q_to_str;
	q_to_str << Q;
	return q_to_str.str();
}



int main(int argc, char* argv[]){
	//////////////////////////////////////////////
	/// predefinition of everything, setting it up
	//////////////////////////////////////////////
	configure(argc, argv, to_string2("DYh.cfg"), true);
  bool PDFmemuse = false;
	double z = 0.5;
	cout.precision(16);

	/*double jacx = 1.-tau;
	double x = 1.;
	cout <<  1./(tau/jacx+x)*xfit_pdfs(4, tau+jacx*x)/(tau+jacx*x)*xfit_pdfs(6, tau/(tau+jacx*x))/(tau/(tau+jacx*x)) << endl;
	vector<results_c> testit;
	Q = 0.5;
	cout << "z = " << Q << endl;
	testit = call_cuhre_dy("test","qqbar","1",fitPDF);
	cout << "O(1/N) " << testit[0].res << endl;
	testit = call_cuhre_dy("test","qqbar","2",fitPDF);
	cout << "O(1/N^2) " << testit[0].res << endl;
	testit = call_cuhre_dy("test","qqbar","3",fitPDF);
	cout << "O(1/N^3) " << testit[0].res << endl;
	testit = call_cuhre_dy("test","qqbar","4",fitPDF);
	cout << "O(1/N^5) " << testit[0].res << endl;
	Q = 0.99;
	cout << "z = " << Q << endl;
	testit = call_cuhre_dy("test","qqbar","1",fitPDF);
	cout << "O(1/N) " << testit[0].res << endl;
	testit = call_cuhre_dy("test","qqbar","2",fitPDF);
	cout << "O(1/N^2) " << testit[0].res << endl;
	testit = call_cuhre_dy("test","qqbar","3",fitPDF);
	cout << "O(1/N^3) " << testit[0].res<< endl;
	testit = call_cuhre_dy("test","qqbar","4",fitPDF);
	cout << "O(1/N^5) " << testit[0].res << endl;
  exit(0);
*/
	bool doexpansions = true;

	//////////////////////
	/// LO declarations
	//////////////////////
	// DY
	vector<results_c> DY_LO_qqbar_full;
	// higgs
	vector<results_c> higgs_LO_gg_full;
	////////////////////////////////////////////////////////////////////

	//////////////////////
	/// NLO declarations
	//////////////////////
	// DY
	vector<results_c> res_DY_NLO_qqbar_full, res_DY_NLO_qqbar_hard, res_DY_NLO_qqbar_LP, res_DY_NLO_qqbar_LP_part1, res_DY_NLO_qqbar_LP_cor,res_DY_NLO_qqbar_delta, res_DY_NLO_qqbar_accum, DY_NLO_qqbar_full;
	vector<results_c> res_DY_NLO_qg_full, res_DY_NLO_qg_accum;
	vector<results_c> DY_NLO_qg_powers, DY_NLO_qqbar_powers;
	// higgs
	vector<results_c> res_higgs_NLO_gg_full, res_higgs_NLO_gg_hard, res_higgs_NLO_gg_LP, res_higgs_NLO_gg_LP_part1, res_higgs_NLO_gg_LP_cor,res_higgs_NLO_gg_delta, res_higgs_NLO_gg_accum, higgs_NLO_gg_full;
	vector<results_c> res_higgs_NLO_qg_full, res_higgs_NLO_qg_accum;
	vector<results_c> res_higgs_NLO_qqbar_full,res_higgs_NLO_qqbar_accum;
	vector<results_c> higgs_NLO_qg_powers, higgs_NLO_gg_powers, higgs_NLO_qqbar_powers;
	////////////////////////////////////////////////////////////////////

	//////////////////////
	/// NNLO declarations
	//////////////////////
	// DY
	vector<results_c> res_DY_NNLO_qqbar_full, DY_NNLO_qqbar_full, res_DY_NNLO_qqbar_hard, res_DY_NNLO_qqbar_LP, res_DY_NNLO_qqbar_LP_part1, res_DY_NNLO_qqbar_LP_cor,res_DY_NNLO_qqbar_delta, res_DY_NNLO_qqbar_accum;
	vector<results_c> res_DY_NNLO_qg_full, res_DY_NNLO_qg_accum;
	vector<results_c> res_DY_NNLO_qq_full, res_DY_NNLO_qq_accum;
	vector<results_c> res_DY_NNLO_gg_full, res_DY_NNLO_gg_accum;
	vector<results_c> DY_NNLO_gg_powers, DY_NNLO_qg_powers, DY_NNLO_qq_powers, DY_NNLO_qqbar_powers;
	// higgs
	vector<results_c> res_higgs_NNLO_gg_full, higgs_NNLO_gg_full, res_higgs_NNLO_gg_hard, res_higgs_NNLO_gg_LP, res_higgs_NNLO_gg_LP_part1, res_higgs_NNLO_gg_LP_cor,res_higgs_NNLO_gg_delta, res_higgs_NNLO_gg_accum;
	vector<results_c> res_higgs_NNLO_qg_full, res_higgs_NNLO_qg_accum;
	vector<results_c> res_higgs_NNLO_qq_full, res_higgs_NNLO_qq_accum;
	vector<results_c> res_higgs_NNLO_qqbar_full, res_higgs_NNLO_qqbar_accum;
	vector<results_c> res_higgs_NNLO_qqp_full, res_higgs_NNLO_qqp_accum;
	vector<results_c> higgs_NNLO_gg_powers, higgs_NNLO_qg_powers, higgs_NNLO_qq_powers, higgs_NNLO_qqp_powers, higgs_NNLO_qqbar_powers;
	////////////////////////////////////////////////////////////////////

	///////////////////////
	/// resummed declarations
	///////////////////////

	//DY
	vector<results_c> resummed_DY_NLP_exp_LL, resummed_DY_NLP_exp_NLL, resummed_DY_NLP_exp_NNLL, resummed_DY_NLP_qg;
	vector<results_c> resummed_DY_LP_NNNLL, resummed_DY_LP_NNNLL_NLP_LL, resummed_DY_LP_NNLL, resummed_DY_LP_NNLL_NLP_LL, resummed_DY_LP_NLL, resummed_DY_LP_NLL_NLP_LL, resummed_DY_LP_LL_NLP_LL, resummed_DY_LP_LL;
	vector<results_c> resummed_DY_LP_NNNLL_exp_NNLO, resummed_DY_LP_NNNLL_NLP_LL_exp_NNLO,resummed_DY_LP_NNLL_exp_NNLO, resummed_DY_LP_NNLL_NLP_LL_exp_NNLO, resummed_DY_LP_NLL_exp_NNLO, resummed_DY_LP_NLL_NLP_LL_exp_NNLO, resummed_DY_LP_LL_NLP_LL_exp_NNLO, resummed_DY_LP_LL_exp_NNLO;
	vector<results_c> resummed_DY_LP_NNNLL_exp_NLO, resummed_DY_LP_NNNLL_NLP_LL_exp_NLO,resummed_DY_LP_NNLL_exp_NLO, resummed_DY_LP_NNLL_NLP_LL_exp_NLO, resummed_DY_LP_NLL_exp_NLO, resummed_DY_LP_NLL_NLP_LL_exp_NLO, resummed_DY_LP_LL_NLP_LL_exp_NLO, resummed_DY_LP_LL_exp_NLO;
	vector<results_c> resummed_DY_LP_NLLp, resummed_DY_LP_NLLp_NLP_LL;
	vector<results_c> resummed_DY_LP_NLLp_exp_NLO, resummed_DY_LP_NLLp_NLP_LL_exp_NLO;
	vector<results_c> SCET_DY_LP_NNLL, SCET_DY_LP_NNLL_NLP_LL, SCET_DY_LP_NLL, SCET_DY_LP_NLL_NLP_LL, SCET_DY_LP_LL_NLP_LL, SCET_DY_LP_LL;

	// higgs
	vector<results_c> resummed_higgs_NLP_exp_LL, resummed_higgs_NLP_exp_NLL, resummed_higgs_NLP_exp_NNLL, resummed_higgs_NLP_qg;
	vector<results_c> resummed_higgs_LP_NNNLL, resummed_higgs_LP_NNNLL_NLP_LL,resummed_higgs_LP_NNLL, resummed_higgs_LP_NNLL_NLP_LL, resummed_higgs_LP_NLL, resummed_higgs_LP_NLL_NLP_LL, resummed_higgs_LP_LL_NLP_LL, resummed_higgs_LP_LL;
	vector<results_c> resummed_higgs_LP_NNNLL_exp_NNLO, resummed_higgs_LP_NNNLL_NLP_LL_exp_NNLO,resummed_higgs_LP_NNLL_exp_NNLO, resummed_higgs_LP_NNLL_NLP_LL_exp_NNLO, resummed_higgs_LP_NLL_exp_NNLO, resummed_higgs_LP_NLL_NLP_LL_exp_NNLO, resummed_higgs_LP_LL_NLP_LL_exp_NNLO, resummed_higgs_LP_LL_exp_NNLO;
	vector<results_c> resummed_higgs_LP_NNNLL_exp_NLO, resummed_higgs_LP_NNNLL_NLP_LL_exp_NLO,resummed_higgs_LP_NNLL_exp_NLO, resummed_higgs_LP_NNLL_NLP_LL_exp_NLO, resummed_higgs_LP_NLL_exp_NLO, resummed_higgs_LP_NLL_NLP_LL_exp_NLO, resummed_higgs_LP_LL_NLP_LL_exp_NLO, resummed_higgs_LP_LL_exp_NLO;
	vector<results_c> resummed_higgs_LP_NLLp, resummed_higgs_LP_NLLp_NLP_LL;
	vector<results_c> resummed_higgs_LP_NLLp_exp_NLO, resummed_higgs_LP_NLLp_NLP_LL_exp_NLO;
	vector<results_c> SCET_higgs_LP_NNLL, SCET_higgs_LP_NNLL_NLP_LL, SCET_higgs_LP_NLL, SCET_higgs_LP_NLL_NLP_LL, SCET_higgs_LP_LL_NLP_LL, SCET_higgs_LP_LL;

	cout << "computing the resummed results" << endl;
	cout << "LP NNLL'" << endl;
	ISNNNLL = 0;
	ISNLP = 0;
	ISNNLL = 1;
	ISNLL = 1;
	ISLL = 1;
	resummed_higgs_LP_NLL = call_cuhre_higgs("resum","gg","pQCD_prime",fitPDF);
	cout << resummed_higgs_LP_NLL[0].res << endl;
	cout << "LP NNLL" << endl;
	ISNNNLL = 0;
	ISNLP = 0;
	ISNNLL = 1;
	ISNLL = 1;
	ISLL = 1;
	resummed_higgs_LP_NLL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
	cout << resummed_higgs_LP_NLL[0].res << endl;
	cout << "LP NLL'" << endl;
	ISNNNLL = 0;
	ISNLP = 0;
	ISNNLL = 0;
	ISNLL = 1;
	ISLL = 1;
	resummed_higgs_LP_NLL = call_cuhre_higgs("resum","gg","pQCD_prime",fitPDF);
	cout << resummed_higgs_LP_NLL[0].res << endl;
	cout << "LP NLL" << endl;
	ISNNNLL = 0;
	ISNLP = 0;
	ISNNLL = 0;
	ISNLL = 1;
	ISLL = 1;
	resummed_higgs_LP_NLL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
	cout << resummed_higgs_LP_NLL[0].res << endl;

	cout << "LP NNLL'+NLP" << endl;
	ISNNNLL = 0;
	ISNLP = 1;
	ISNNLL = 1;
	ISNLL = 1;
	ISLL = 1;
	resummed_higgs_LP_NLL = call_cuhre_higgs("resum","gg","pQCD_prime",fitPDF);
	cout << resummed_higgs_LP_NLL[0].res << endl;
	cout << "LP NNLL+NLP" << endl;
	ISNNNLL = 0;
	ISNLP = 1;
	ISNNLL = 1;
	ISNLL = 1;
	ISLL = 1;
	resummed_higgs_LP_NLL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
	cout << resummed_higgs_LP_NLL[0].res << endl;
	cout << "LP NLL'+NLP" << endl;
	ISNNNLL = 0;
	ISNLP = 1;
	ISNNLL = 0;
	ISNLL = 1;
	ISLL = 1;
	resummed_higgs_LP_NLL = call_cuhre_higgs("resum","gg","pQCD_prime",fitPDF);
	cout << resummed_higgs_LP_NLL[0].res << endl;
	cout << "LP NLL+NLP" << endl;
	ISNNNLL = 0;
	ISNLP = 1;
	ISNNLL = 0;
	ISNLL = 1;
	ISLL = 1;
	resummed_higgs_LP_NLL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
	cout << resummed_higgs_LP_NLL[0].res << endl;
 exit(0);
	////////////////////////////////////////////////////////////////////

	///////////////////
	/// output
	///////////////////
	/*complex<double> Nint = 2.5;
	string lumchan = "qqbar"; complex<double> lumni = (complex<double>) LumN(Nint, 10, lumchan);
	cout << "qqbar " << lumni << endl;
  lumchan = "qg_charge"; lumni = (complex<double>) LumN(Nint, 10, lumchan);
	cout << "qg_charge " << lumni << endl;
	Nint = 1.5;
	lumchan = "qqbar"; lumni = (complex<double>) LumN(Nint, 10, lumchan);
	cout << "qqbar " << lumni << endl;
	lumchan = "qg_charge"; lumni = (complex<double>) LumN(Nint, 10, lumchan);
	cout << "qg_charge " << lumni << endl;

	cout << "WARNING LOOK AT DEFINITION OF NNLL or NNLL' " << endl;
	cout << "EXIT TO MAKE SURE THAT YOU DO " << endl;
	exit(0);*/
	ofstream output;
	ostringstream x_convert; // need this for the output
	x_convert << Q;
	string Qstring  = x_convert.str();
	ostringstream x_convert2; // need this for the output
	x_convert2 << alphas_Q;
	string asstring  = x_convert2.str();
	ostringstream x_convert3; // need this for the output
	x_convert3 << muF;
	string mufstring  = x_convert3.str();
	ostringstream x_convert4; // need this for the output
	x_convert4 << muR;
	string murstring  = x_convert4.str();
	string q_str = "results/DYh_16122020_without_prime/output_Q" + Qstring +"_as"+asstring+"_"+setname+"_muF"+mufstring+"_muR"+murstring;
	if(DY) q_str = q_str+"_DY";
	else if(higgs) q_str = q_str+"_Higgs";
	else {cout << "process not well defined... exiting" << endl; exit(0);}
	if((DY) || (higgs)){
		if(LO) q_str = q_str+"_LO";
		if(NLO) q_str = q_str+"_NLO";
		if(NNLO) q_str = q_str+"_NNLO";
	}
	if(doexpansions){q_str = q_str+"_with_expansions";}
	if(PDFmemuse){q_str = q_str+"_members"; RES=false; fitPDF=false;}
	if(chebPDF){q_str = q_str+"_cheb_pdfs"; fitPDF=false;} //if real pdfs, the resummed cannot be used!
	else if(fitPDF){q_str = q_str+"_fitted_pdfs";}
	else if(toy_pdfs){q_str = q_str+"_toy_pdfs";}
	else if(realPDF){q_str = q_str+"_real_pdfs"; RES=false; fitPDF=false;} //if real pdfs, the resummed cannot be used!
	if(RES) q_str = q_str+"_resummed";
	q_str = q_str + ".txt";
	cout << q_str.c_str() << endl;
	////////////////////////////////////////////////////////////////////

	/////////////////////////
	/// computation
	/////////////////////////
	int powerexpansion = 20;
	int NLPexpansion = 10;
	if(higgs){
		cout << "working on higgs" << endl;
		if((PDFmemuse == true) && (fitPDF == false)){
			for(unsigned int i = 0; i < 101;i++)
				{
					use_member = i;
					update_defaults();

					higgs_LO_gg_full.push_back(call_cuhre_higgs("LO","gg","LP",fitPDF)[0]);
					res_higgs_NLO_gg_hard.push_back(call_cuhre_higgs("NLO","gg","reg",fitPDF)[0]);
					res_higgs_NLO_gg_LP_part1.push_back(call_cuhre_higgs("NLO","gg","LP",fitPDF)[0]);
					res_higgs_NLO_gg_LP_cor.push_back(call_cuhre_higgs("NLO","gg","LP_corr",fitPDF)[0]);
					res_higgs_NLO_gg_delta.push_back(call_cuhre_higgs("NLO","gg","delta",fitPDF)[0]);
					res_higgs_NLO_qg_full.push_back(call_cuhre_higgs("NLO","qg","full",fitPDF)[0]);
					res_higgs_NLO_qqbar_full.push_back(call_cuhre_higgs("NLO","qqbar","full",fitPDF)[0]);
					res_higgs_NLO_gg_LP.push_back({res_higgs_NLO_gg_LP_part1[i].res + res_higgs_NLO_gg_LP_cor[i].res, res_higgs_NLO_gg_LP_part1[i].err + res_higgs_NLO_gg_LP_cor[i].err, res_higgs_NLO_gg_LP_part1[i].prob + res_higgs_NLO_gg_LP_cor[i].prob});
					higgs_NLO_gg_full.push_back({res_higgs_NLO_gg_hard[i].res+res_higgs_NLO_gg_LP_cor[i].res+res_higgs_NLO_gg_delta[i].res, res_higgs_NLO_gg_hard[i].err+res_higgs_NLO_gg_LP_cor[i].err+res_higgs_NLO_gg_delta[i].err, res_higgs_NLO_gg_hard[i].prob+res_higgs_NLO_gg_LP_cor[i].prob+res_higgs_NLO_gg_delta[i].prob});
					res_higgs_NNLO_gg_hard.push_back(call_cuhre_higgs("NNLO","gg","reg",fitPDF)[0]);
					res_higgs_NNLO_gg_LP_part1.push_back(call_cuhre_higgs("NNLO","gg","LP",fitPDF)[0]);
					res_higgs_NNLO_gg_LP_cor.push_back(call_cuhre_higgs("NNLO","gg","LP_corr",fitPDF)[0]);
					res_higgs_NNLO_gg_delta.push_back(call_cuhre_higgs("NNLO","gg","delta",fitPDF)[0]);
					res_higgs_NNLO_qg_full.push_back(call_cuhre_higgs("NNLO","qg","full",fitPDF)[0]);
					res_higgs_NNLO_qq_full.push_back(call_cuhre_higgs("NNLO","qq","full",fitPDF)[0]);
					res_higgs_NNLO_qqp_full.push_back(call_cuhre_higgs("NNLO","qqp","full",fitPDF)[0]);
					res_higgs_NNLO_qqbar_full.push_back(call_cuhre_higgs("NNLO","qqbar","full",fitPDF)[0]);
					res_higgs_NNLO_gg_LP.push_back({res_higgs_NNLO_gg_LP_part1[i].res + res_higgs_NNLO_gg_LP_cor[i].res, res_higgs_NNLO_gg_LP_part1[i].err + res_higgs_NNLO_gg_LP_cor[i].err, res_higgs_NNLO_gg_LP_part1[i].prob + res_higgs_NNLO_gg_LP_cor[i].prob});
					higgs_NNLO_gg_full.push_back({res_higgs_NNLO_gg_hard[i].res+res_higgs_NNLO_gg_LP_cor[i].res+res_higgs_NNLO_gg_delta[i].res, res_higgs_NNLO_gg_hard[i].err+res_higgs_NNLO_gg_LP_cor[i].err+res_higgs_NNLO_gg_delta[i].err, res_higgs_NNLO_gg_hard[i].prob+res_higgs_NNLO_gg_LP_cor[i].prob+res_higgs_NNLO_gg_delta[i].prob});
				}
		}
		else{
			if(LO&&higgs){
				cout << "working on LO" << endl;
				higgs_LO_gg_full = call_cuhre_higgs("LO","gg","LP",fitPDF);
				cout << higgs_LO_gg_full[0].res << " " << higgs_LO_gg_full[0].err << endl;
			}
			if(NLO&&higgs){
				cout << "working on NLO" << endl;
				res_higgs_NLO_gg_hard = call_cuhre_higgs("NLO","gg","reg",fitPDF);
				res_higgs_NLO_gg_LP_part1 = call_cuhre_higgs("NLO","gg","LP",fitPDF);
				res_higgs_NLO_gg_LP_cor = call_cuhre_higgs("NLO","gg","LP_corr",fitPDF);
				res_higgs_NLO_gg_delta = call_cuhre_higgs("NLO","gg","delta",fitPDF);
				res_higgs_NLO_qg_full = call_cuhre_higgs("NLO","qg","full",fitPDF);
				res_higgs_NLO_qqbar_full = call_cuhre_higgs("NLO","qqbar","full",fitPDF);
				res_higgs_NLO_gg_LP.push_back({res_higgs_NLO_gg_LP_part1[0].res + res_higgs_NLO_gg_LP_cor[0].res, res_higgs_NLO_gg_LP_part1[0].err + res_higgs_NLO_gg_LP_cor[0].err, res_higgs_NLO_gg_LP_part1[0].prob + res_higgs_NLO_gg_LP_cor[0].prob});
				if(doexpansions){
					higgs_NLO_gg_powers = call_cuhre_higgs("NLO","gg","exp",fitPDF,powerexpansion);
					higgs_NLO_qg_powers = call_cuhre_higgs("NLO","qg","exp",fitPDF,powerexpansion);
					higgs_NLO_qqbar_powers = call_cuhre_higgs("NLO","qqbar","exp",fitPDF,powerexpansion);
				}
				higgs_NLO_gg_full.push_back({res_higgs_NLO_gg_hard[0].res+res_higgs_NLO_gg_LP_cor[0].res+res_higgs_NLO_gg_delta[0].res, res_higgs_NLO_gg_hard[0].err+res_higgs_NLO_gg_LP_cor[0].err+res_higgs_NLO_gg_delta[0].err, res_higgs_NLO_gg_hard[0].prob+res_higgs_NLO_gg_LP_cor[0].prob+res_higgs_NLO_gg_delta[0].prob});

			}
			if(NNLO&&higgs){
				cout << "working on NNLO" << endl;
				res_higgs_NNLO_gg_hard = call_cuhre_higgs("NNLO","gg","reg",fitPDF);
				res_higgs_NNLO_gg_LP_part1 = call_cuhre_higgs("NNLO","gg","LP",fitPDF);
				res_higgs_NNLO_gg_LP_cor = call_cuhre_higgs("NNLO","gg","LP_corr",fitPDF);
				res_higgs_NNLO_gg_delta = call_cuhre_higgs("NNLO","gg","delta",fitPDF);
				res_higgs_NNLO_qg_full = call_cuhre_higgs("NNLO","qg","full",fitPDF);
				res_higgs_NNLO_qq_full = call_cuhre_higgs("NNLO","qq","full",fitPDF);
				res_higgs_NNLO_qqp_full = call_cuhre_higgs("NNLO","qqp","full",fitPDF);
				res_higgs_NNLO_qqbar_full = call_cuhre_higgs("NNLO","qqbar","full",fitPDF);
				res_higgs_NNLO_gg_LP.push_back({res_higgs_NNLO_gg_LP_part1[0].res + res_higgs_NNLO_gg_LP_cor[0].res, res_higgs_NNLO_gg_LP_part1[0].err + res_higgs_NNLO_gg_LP_cor[0].err, res_higgs_NNLO_gg_LP_part1[0].prob + res_higgs_NNLO_gg_LP_cor[0].prob});
				if(doexpansions){
					higgs_NNLO_gg_powers = call_cuhre_higgs("NNLO","gg","exp",fitPDF,powerexpansion);
					higgs_NNLO_qg_powers = call_cuhre_higgs("NNLO","qg","exp",fitPDF,powerexpansion);
					higgs_NNLO_qq_powers = call_cuhre_higgs("NNLO","qq","exp",fitPDF,powerexpansion);
					higgs_NNLO_qqp_powers = call_cuhre_higgs("NNLO","qqp","exp",fitPDF,powerexpansion);
					higgs_NNLO_qqbar_powers = call_cuhre_higgs("NLO","qqbar","exp",fitPDF,powerexpansion);
				}
				double error = higgs_scale_factor(higgs_LO_gg_full[0].res,higgs_NLO_gg_full[0].res);
				higgs_NNLO_gg_full.push_back({error+res_higgs_NNLO_gg_hard[0].res+res_higgs_NNLO_gg_LP_cor[0].res+res_higgs_NNLO_gg_delta[0].res, res_higgs_NNLO_gg_hard[0].err+res_higgs_NNLO_gg_LP_cor[0].err+res_higgs_NNLO_gg_delta[0].err, res_higgs_NNLO_gg_hard[0].prob+res_higgs_NNLO_gg_LP_cor[0].prob+res_higgs_NNLO_gg_delta[0].prob});
				error = higgs_scale_factor(0.,res_higgs_NLO_qg_full[0].res);
				res_higgs_NNLO_qg_full[0].res += error;
				error = higgs_scale_factor(0.,res_higgs_NLO_qqbar_full[0].res);
				res_higgs_NNLO_qqbar_full[0].res += error;
			}
			if(RES&&higgs){
				cout << "computing the resummed results" << endl;
				cout << "LP NNNLL + NLP LL" << endl;
				ISNNNLL = 1;
				ISNNLL = 1;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 1;
				resummed_higgs_LP_NNNLL_NLP_LL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
				resummed_higgs_LP_NNNLL_NLP_LL_exp_NLO =call_cuhre_higgs("resum","gg","pQCDNLO",fitPDF);
				resummed_higgs_LP_NNNLL_NLP_LL_exp_NNLO = call_cuhre_higgs("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_higgs_LP_NNNLL_NLP_LL[0].res << " " << resummed_higgs_LP_NNNLL_NLP_LL[0].err << endl;

				cout << "LP NNLL + NLP LL" << endl;
				ISNNNLL = 0;
				ISNNLL = 1;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 1;
				resummed_higgs_LP_NNLL_NLP_LL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
				resummed_higgs_LP_NNLL_NLP_LL_exp_NLO =call_cuhre_higgs("resum","gg","pQCDNLO",fitPDF);
				resummed_higgs_LP_NNLL_NLP_LL_exp_NNLO = call_cuhre_higgs("resum","gg","pQCDNNLO",fitPDF);
				resummed_higgs_NLP_exp_NNLL = call_cuhre_higgs("resum","gg","NLPexp",fitPDF, NLPexpansion);
				resummed_higgs_NLP_qg = call_cuhre_higgs("resum","gg","offdiag",fitPDF, NLPexpansion);
				cout << resummed_higgs_LP_NNLL_NLP_LL[0].res << " " << resummed_higgs_LP_NNLL_NLP_LL[0].err << endl;
				cout << "LP NLL + NLP LL" << endl;
				ISNNNLL = 0;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 1;
				resummed_higgs_LP_NLLp_NLP_LL  = call_cuhre_higgs("resum","gg","pQCD_prime",fitPDF);
				resummed_higgs_LP_NLLp_NLP_LL_exp_NLO = call_cuhre_higgs("resum","gg","pQCDNLO_prime",fitPDF);
				resummed_higgs_LP_NLL_NLP_LL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
				resummed_higgs_LP_NLL_NLP_LL_exp_NLO = call_cuhre_higgs("resum","gg","pQCDNLO",fitPDF);
				resummed_higgs_LP_NLL_NLP_LL_exp_NNLO = call_cuhre_higgs("resum","gg","pQCDNNLO",fitPDF);
				resummed_higgs_NLP_exp_NLL = call_cuhre_higgs("resum","gg","NLPexp",fitPDF, NLPexpansion);
				cout << resummed_higgs_LP_NLL_NLP_LL[0].res << " " << resummed_higgs_LP_NLL_NLP_LL[0].err << endl;
				for(unsigned int i = 0; i < resummed_higgs_NLP_exp_NLL.size(); i++){
							cout << "	order as "<<i<<" : " << resummed_higgs_NLP_exp_NLL[i].res << endl;
						}
				resummed_higgs_NLP_exp_NLL = call_cuhre_higgs("resum","gg","NLPexp",fitPDF, NLPexpansion);
				cout << "NLO matched" << endl;
				for(unsigned int i = 0; i < resummed_higgs_NLP_exp_NLL.size(); i++){
							cout << "	order as "<<i<<" : " << resummed_higgs_NLP_exp_NLL[i].res << endl;
								}
				cout << "LP LL + NLP LL" << endl;
				ISNNNLL = 0;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 0;
				resummed_higgs_LP_LL_NLP_LL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
				resummed_higgs_LP_LL_NLP_LL_exp_NLO = call_cuhre_higgs("resum","gg","pQCDNLO",fitPDF);
				resummed_higgs_LP_LL_NLP_LL_exp_NNLO = call_cuhre_higgs("resum","gg","pQCDNNLO",fitPDF);
				resummed_higgs_NLP_exp_LL = call_cuhre_higgs("resum","gg","NLPexp",fitPDF, NLPexpansion);

				cout << resummed_higgs_LP_LL_NLP_LL[0].res << " " << resummed_higgs_LP_LL_NLP_LL[0].err << endl;

				cout << "LP NNNLL" << endl;
				ISNNNLL = 1;
				ISNNLL = 1;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 1;
				resummed_higgs_LP_NNNLL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
				resummed_higgs_LP_NNNLL_exp_NLO = call_cuhre_higgs("resum","gg","pQCDNLO",fitPDF);
				resummed_higgs_LP_NNNLL_exp_NNLO = call_cuhre_higgs("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_higgs_LP_NNNLL[0].res << " " << resummed_higgs_LP_NNNLL[0].err << endl;

				cout << "LP NNLL" << endl;
				ISNNNLL = 0;
				ISNNLL = 1;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 1;
				resummed_higgs_LP_NNLL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
				resummed_higgs_LP_NNLL_exp_NLO = call_cuhre_higgs("resum","gg","pQCDNLO",fitPDF);
				resummed_higgs_LP_NNLL_exp_NNLO = call_cuhre_higgs("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_higgs_LP_NNLL[0].res << " " << resummed_higgs_LP_NNLL[0].err << endl;

				cout << "LP NLL" << endl;
				ISNNNLL = 0;
				ISNNLL = 0;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 1;
				resummed_higgs_LP_NLLp  = call_cuhre_higgs("resum","gg","pQCD_prime",fitPDF);
				resummed_higgs_LP_NLLp_exp_NLO = call_cuhre_higgs("resum","gg","pQCDNLO_prime",fitPDF);
				resummed_higgs_LP_NLL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
				resummed_higgs_LP_NLL_exp_NLO = call_cuhre_higgs("resum","gg","pQCDNLO",fitPDF);
				resummed_higgs_LP_NLL_exp_NNLO = call_cuhre_higgs("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_higgs_LP_NLL[0].res << " " << resummed_higgs_LP_NLL[0].err << endl;


				cout << "LP LL" << endl;
				ISNNNLL = 0;
				ISNNLL = 0;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 0;
				resummed_higgs_LP_LL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
				resummed_higgs_LP_LL_exp_NLO = call_cuhre_higgs("resum","gg","pQCDNLO",fitPDF);
				resummed_higgs_LP_LL_exp_NNLO = call_cuhre_higgs("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_higgs_LP_LL[0].res << " " << resummed_higgs_LP_LL[0].err << endl;



		}
	}
}
  if(DY){
		cout << "working on DY" << endl;
		if((PDFmemuse == true) && (fitPDF == false)){
		for(unsigned int i = 0; i < 101;i++)
				{
					use_member = i;
					update_defaults();
					DY_LO_qqbar_full.push_back(call_cuhre_dy("LO","qqbar","LP",fitPDF)[0]);
					res_DY_NLO_qqbar_hard.push_back(call_cuhre_dy("NLO","qqbar","reg",fitPDF)[0]);
					res_DY_NLO_qqbar_LP_part1.push_back(call_cuhre_dy("NLO","qqbar","LP",fitPDF)[0]);
					res_DY_NLO_qqbar_LP_cor.push_back(call_cuhre_dy("NLO","qqbar","LP_corr",fitPDF)[0]);
					res_DY_NLO_qqbar_delta.push_back(call_cuhre_dy("NLO","qqbar","delta",fitPDF)[0]);
					res_DY_NLO_qg_full.push_back(call_cuhre_dy("NLO","qg","full",fitPDF)[0]);
					res_DY_NLO_qqbar_LP.push_back({res_DY_NLO_qqbar_LP_part1[i].res + res_DY_NLO_qqbar_LP_cor[i].res, res_DY_NLO_qqbar_LP_part1[i].err + res_DY_NLO_qqbar_LP_cor[i].err, res_DY_NLO_qqbar_LP_part1[i].prob + res_DY_NLO_qqbar_LP_cor[i].prob});
					DY_NLO_qqbar_full.push_back({res_DY_NLO_qqbar_hard[i].res+res_DY_NLO_qqbar_LP_cor[i].res+res_DY_NLO_qqbar_delta[i].res, res_DY_NLO_qqbar_hard[i].err+res_DY_NLO_qqbar_LP_cor[i].err+res_DY_NLO_qqbar_delta[i].err, res_DY_NLO_qqbar_hard[i].prob+res_DY_NLO_qqbar_LP_cor[i].prob+res_DY_NLO_qqbar_delta[i].prob});
					res_DY_NNLO_qqbar_hard.push_back(call_cuhre_dy("NNLO","qqbar","reg",fitPDF)[0]);
					res_DY_NNLO_qqbar_LP_part1.push_back(call_cuhre_dy("NNLO","qqbar","LP",fitPDF)[0]);
					res_DY_NNLO_qqbar_LP_cor.push_back(call_cuhre_dy("NNLO","qqbar","LP_corr",fitPDF)[0]);
					res_DY_NNLO_qqbar_delta.push_back(call_cuhre_dy("NNLO","qqbar","delta",fitPDF)[0]);
					res_DY_NNLO_qg_full.push_back(call_cuhre_dy("NNLO","qg","full",fitPDF)[0]);
					res_DY_NNLO_qq_full.push_back(call_cuhre_dy("NNLO","qq","full",fitPDF)[0]);
					res_DY_NNLO_gg_full.push_back(call_cuhre_dy("NNLO","gg","full",fitPDF)[0]);
					res_DY_NNLO_qqbar_LP.push_back({res_DY_NNLO_qqbar_LP_part1[i].res + res_DY_NNLO_qqbar_LP_cor[i].res, res_DY_NNLO_qqbar_LP_part1[i].err + res_DY_NNLO_qqbar_LP_cor[i].err, res_DY_NNLO_qqbar_LP_part1[i].prob + res_DY_NNLO_qqbar_LP_cor[i].prob});
					DY_NNLO_qqbar_full.push_back({res_DY_NNLO_qqbar_hard[i].res+res_DY_NNLO_qqbar_LP_cor[i].res+res_DY_NNLO_qqbar_delta[i].res, res_DY_NNLO_qqbar_hard[i].err+res_DY_NNLO_qqbar_LP_cor[i].err+res_DY_NNLO_qqbar_delta[i].err, res_DY_NNLO_qqbar_hard[i].prob+res_DY_NNLO_qqbar_LP_cor[i].prob+res_DY_NNLO_qqbar_delta[i].prob});
					cout << DY_NLO_qqbar_full[i].res << endl;
				}
		}
		else{
			if(LO&&DY){
				cout << "working on LO" << endl;
				DY_LO_qqbar_full = call_cuhre_dy("LO","qqbar","full",fitPDF);
				cout << DY_LO_qqbar_full[0].res << " " << DY_LO_qqbar_full[0].err << endl;
			}
			if(NLO&&DY){
				cout << "working on NLO" << endl;
				res_DY_NLO_qqbar_hard = call_cuhre_dy("NLO","qqbar","reg",fitPDF);
				res_DY_NLO_qqbar_LP_part1 = call_cuhre_dy("NLO","qqbar","LP",fitPDF);
				res_DY_NLO_qqbar_LP_cor = call_cuhre_dy("NLO","qqbar","LP_corr",fitPDF);
				res_DY_NLO_qqbar_delta = call_cuhre_dy("NLO","qqbar","delta",fitPDF);
				res_DY_NLO_qg_full = call_cuhre_dy("NLO","qg","full",fitPDF);
				res_DY_NLO_qqbar_LP.push_back({res_DY_NLO_qqbar_LP_part1[0].res + res_DY_NLO_qqbar_LP_cor[0].res, res_DY_NLO_qqbar_LP_part1[0].err + res_DY_NLO_qqbar_LP_cor[0].err, res_DY_NLO_qqbar_LP_part1[0].prob + res_DY_NLO_qqbar_LP_cor[0].prob});
				if(doexpansions){
					DY_NLO_qqbar_powers = call_cuhre_dy("NLO","qqbar","exp",fitPDF,powerexpansion);
					DY_NLO_qg_powers = call_cuhre_dy("NLO","qg","exp",fitPDF,powerexpansion);
				}
				DY_NLO_qqbar_full.push_back({res_DY_NLO_qqbar_hard[0].res+res_DY_NLO_qqbar_LP_cor[0].res+res_DY_NLO_qqbar_delta[0].res, res_DY_NLO_qqbar_hard[0].err+res_DY_NLO_qqbar_LP_cor[0].err+res_DY_NLO_qqbar_delta[0].err, res_DY_NLO_qqbar_hard[0].prob + res_DY_NLO_qqbar_LP_cor[0].prob+res_DY_NLO_qqbar_delta[0].prob});

			}
			if(NNLO&&DY){
				cout << "working on NNLO" << endl;
				res_DY_NNLO_qqbar_hard = call_cuhre_dy("NNLO","qqbar","reg",fitPDF);
				res_DY_NNLO_qqbar_LP_part1 = call_cuhre_dy("NNLO","qqbar","LP",fitPDF);
				res_DY_NNLO_qqbar_LP_cor = call_cuhre_dy("NNLO","qqbar","LP_corr",fitPDF);
				res_DY_NNLO_qqbar_delta = call_cuhre_dy("NNLO","qqbar","delta",fitPDF);
				res_DY_NNLO_qg_full = call_cuhre_dy("NNLO","qg","full",fitPDF);
				res_DY_NNLO_qq_full = call_cuhre_dy("NNLO","qq","full",fitPDF);
				res_DY_NNLO_gg_full = call_cuhre_dy("NNLO","gg","full",fitPDF);
				res_DY_NNLO_qqbar_LP.push_back({res_DY_NNLO_qqbar_LP_part1[0].res + res_DY_NNLO_qqbar_LP_cor[0].res, res_DY_NNLO_qqbar_LP_part1[0].err + res_DY_NNLO_qqbar_LP_cor[0].err, res_DY_NNLO_qqbar_LP_part1[0].prob + res_DY_NNLO_qqbar_LP_cor[0].prob});
				if(doexpansions){
					DY_NNLO_gg_powers = call_cuhre_dy("NNLO","gg","exp",fitPDF,powerexpansion);
					DY_NNLO_qg_powers = call_cuhre_dy("NNLO","qg","exp",fitPDF,powerexpansion);
					DY_NNLO_qq_powers = call_cuhre_dy("NNLO","qq","exp",fitPDF,powerexpansion);
				}
				DY_NNLO_qqbar_powers = call_cuhre_dy("NNLO","qqbar","exp",fitPDF,powerexpansion);
				DY_NNLO_qqbar_full.push_back({res_DY_NNLO_qqbar_hard[0].res+res_DY_NNLO_qqbar_LP_cor[0].res+res_DY_NNLO_qqbar_delta[0].res, res_DY_NNLO_qqbar_hard[0].err+res_DY_NNLO_qqbar_LP_cor[0].err+res_DY_NNLO_qqbar_delta[0].err, res_DY_NNLO_qqbar_hard[0].prob+res_DY_NNLO_qqbar_LP_cor[0].prob+res_DY_NNLO_qqbar_delta[0].prob});
			}
			if(RES&&DY){
				cout << "computing the resummed results" << endl;
				cout << "LP NNNLL + NLP LL" << endl;
				ISNNNLL = 1;
				ISNNLL = 1;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 1;
				resummed_DY_LP_NNNLL_NLP_LL = call_cuhre_dy("resum","gg","pQCD",fitPDF);
				resummed_DY_LP_NNNLL_NLP_LL_exp_NLO =call_cuhre_dy("resum","gg","pQCDNLO",fitPDF);
				resummed_DY_LP_NNNLL_NLP_LL_exp_NNLO = call_cuhre_dy("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_DY_LP_NNNLL_NLP_LL[0].res << " " << resummed_DY_LP_NNNLL_NLP_LL[0].err << endl;

				cout << "LP NNLL + NLP LL" << endl;
				ISNNNLL = 0;
				ISNNLL = 1;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 1;
				cout <<  ISNNNLL << endl;
				resummed_DY_LP_NNLL_NLP_LL = call_cuhre_dy("resum","gg","pQCD",fitPDF);
				resummed_DY_LP_NNLL_NLP_LL_exp_NLO =call_cuhre_dy("resum","gg","pQCDNLO",fitPDF);
				resummed_DY_LP_NNLL_NLP_LL_exp_NNLO = call_cuhre_dy("resum","gg","pQCDNNLO",fitPDF);
				resummed_DY_NLP_exp_NNLL = call_cuhre_dy("resum","gg","NLPexp",fitPDF, NLPexpansion);
				resummed_DY_NLP_qg = call_cuhre_dy("resum","gg","offdiag",fitPDF, NLPexpansion);
				cout << resummed_DY_LP_NNLL_NLP_LL[0].res << " " << resummed_DY_LP_NNLL_NLP_LL[0].err << endl;
				cout <<  ISNNNLL << endl;

				cout << "LP NLL + NLP LL" << endl;
				ISNNNLL = 0;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 1;

				resummed_DY_LP_NLLp_NLP_LL  = call_cuhre_dy("resum","gg","pQCD_prime",fitPDF);
				resummed_DY_LP_NLLp_NLP_LL_exp_NLO = call_cuhre_dy("resum","gg","pQCDNLO_prime",fitPDF);
				resummed_DY_LP_NLL_NLP_LL = call_cuhre_dy("resum","gg","pQCD",fitPDF);
				resummed_DY_LP_NLL_NLP_LL_exp_NLO = call_cuhre_dy("resum","gg","pQCDNLO",fitPDF);
				resummed_DY_LP_NLL_NLP_LL_exp_NNLO = call_cuhre_dy("resum","gg","pQCDNNLO",fitPDF);
				resummed_DY_NLP_exp_NLL = call_cuhre_dy("resum","gg","NLPexp",fitPDF, NLPexpansion);
				cout << resummed_DY_LP_NLL_NLP_LL[0].res << " " << resummed_DY_LP_NLL_NLP_LL[0].err << endl;
				for(unsigned int i = 0; i < resummed_DY_NLP_exp_NLL.size(); i++){
							cout << "	order as "<<i<<" : " << resummed_DY_NLP_exp_NLL[i].res << endl;
						}

				cout << "LP LL + NLP LL" << endl;
				ISNNNLL = 0;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 0;
				resummed_DY_LP_LL_NLP_LL = call_cuhre_dy("resum","gg","pQCD",fitPDF);
				resummed_DY_LP_LL_NLP_LL_exp_NLO = call_cuhre_dy("resum","gg","pQCDNLO",fitPDF);
				resummed_DY_LP_LL_NLP_LL_exp_NNLO = call_cuhre_dy("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_DY_LP_LL_NLP_LL[0].res << " " << resummed_DY_LP_LL_NLP_LL[0].err << endl;

				cout << "LP NNNLL" << endl;
				ISNNNLL = 1;
				ISNNLL = 1;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 1;
				resummed_DY_LP_NNNLL = call_cuhre_dy("resum","gg","pQCD",fitPDF);
				resummed_DY_LP_NNNLL_exp_NLO = call_cuhre_dy("resum","gg","pQCDNLO",fitPDF);
				resummed_DY_LP_NNNLL_exp_NNLO = call_cuhre_dy("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_DY_LP_NNNLL[0].res << " " << resummed_DY_LP_NNNLL[0].err << endl;


				cout << "LP NNLL" << endl;
				ISNNNLL = 0;
				ISNNLL = 1;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 1;
				resummed_DY_LP_NNLL = call_cuhre_dy("resum","gg","pQCD",fitPDF);
				resummed_DY_LP_NNLL_exp_NLO = call_cuhre_dy("resum","gg","pQCDNLO",fitPDF);
				resummed_DY_LP_NNLL_exp_NNLO = call_cuhre_dy("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_DY_LP_NNLL[0].res << " " << resummed_DY_LP_NNLL[0].err << endl;


				cout << "LP NLL" << endl;
				ISNNNLL = 0;
				ISNNLL = 0;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 1;
				resummed_DY_LP_NLLp  = call_cuhre_dy("resum","gg","pQCD_prime",fitPDF);
				resummed_DY_LP_NLLp_exp_NLO = call_cuhre_dy("resum","gg","pQCDNLO_prime",fitPDF);
				resummed_DY_LP_NLL = call_cuhre_dy("resum","gg","pQCD",fitPDF);
				resummed_DY_LP_NLL_exp_NLO = call_cuhre_dy("resum","gg","pQCDNLO",fitPDF);
				resummed_DY_LP_NLL_exp_NNLO = call_cuhre_dy("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_DY_LP_NLL[0].res << " " << resummed_DY_LP_NLL[0].err << endl;


				cout << "LP LL" << endl;
				ISNNNLL = 0;
				ISNNLL = 0;
				ISNLP = 0;
				ISLL = 1;
				ISNLL = 0;
				resummed_DY_LP_LL = call_cuhre_dy("resum","gg","pQCD",fitPDF);
				resummed_DY_LP_LL_exp_NLO = call_cuhre_dy("resum","gg","pQCDNLO",fitPDF);
				resummed_DY_LP_LL_exp_NNLO = call_cuhre_dy("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_DY_LP_LL[0].res << " " << resummed_DY_LP_LL[0].err << endl;



			}
		}
	}

	/////////////////////////
	/// printouts
	/////////////////////////
	output.open(q_str.c_str()); //.c_str() needed to do a constant string conversion
	//output.open(result_map.c_str());
	cout << "Making the printouts" << endl;
	cout << "muF = " << muF << " GeV ,muR = " << muR << " GeV, alphas(muR2) = " << alphas_muR << " ,alphas(Q2) = " << alphas_Q << " ,alphas(muF2) = " << alphas_muF << endl;
	output << "muF = " << muF << " GeV ,muR = " << muR << " GeV, alphas(muR2) = " << alphas_muR << " ,alphas(Q2) = " << alphas_Q << " ,alphas(muF2) = " << alphas_muF << endl;
	if(higgs){
		output << "======================================================" << endl;
		output << "Higgs results" << endl;
		output << "======================================================" << endl;
		if(PDFmemuse)
		{
			output << "LO, NLOgg, NLOqg, NLOqqbar, NNLOgg, NNLOqg, NNLOqqbar, NNLOqq, NNLOqqp \n" << endl;
			for(unsigned int i = 0; i < higgs_LO_gg_full.size(); i++){
			output << "	member "<<i<<" : " <<  higgs_LO_gg_full[i].res << "," <<  higgs_NLO_gg_full[i].res << "," <<  res_higgs_NLO_qg_full[i].res << ","
			<<  res_higgs_NLO_qqbar_full[i].res << "," <<  higgs_NNLO_gg_full[i].res << ","  <<  res_higgs_NNLO_qg_full[i].res << ","
			 <<  res_higgs_NNLO_qqbar_full[i].res << ","  <<  res_higgs_NNLO_qq_full[i].res << ","  <<  res_higgs_NNLO_qqp_full[i].res << endl;

			}
		}
		else
		{
		if(LO&&higgs){
			cout << "LO" << endl;
			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "LO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////

			output << "---------------------------------------" << endl;

			output << "Total xsec:					" << higgs_LO_gg_full[0].res << " pb +/- " << higgs_LO_gg_full[0].err << " " << higgs_LO_gg_full[0].prob << endl;
		}
		if(NLO&&higgs){

			cout << "NLO" << endl;
			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "NLO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////
			output << "gg channel" << endl;
			/////////////////////////////////////
			output << "Z-dependent, including LP (NLO):		" <<  res_higgs_NLO_gg_hard[0].res << " pb +/- " << res_higgs_NLO_gg_hard[0].err <<  endl;
			output << "Constant piece (NLO):				" <<  res_higgs_NLO_gg_delta[0].res << " pb +/- " << res_higgs_NLO_gg_delta[0].err <<  endl;
			output << "Fractional gg xsec (at NLO):			" <<  higgs_NLO_gg_full[0].res << " pb " <<  endl;
			if(doexpansions){
				res_higgs_NLO_gg_accum.push_back(res_higgs_NLO_gg_LP[0]);
				output << "	power "<<0<<" : " << res_higgs_NLO_gg_accum[0].res << " pb, fractional: "<< res_higgs_NLO_gg_accum[0].res/higgs_NLO_gg_full[0].res << endl;
				for(unsigned int i = 0; i < higgs_NLO_gg_powers.size(); i++){
							res_higgs_NLO_gg_accum[0].res = res_higgs_NLO_gg_accum[0].res+higgs_NLO_gg_powers[i].res;
							output << "	power "<<i+1<<" : " << res_higgs_NLO_gg_accum[0].res << " pb, fractional: "<< res_higgs_NLO_gg_accum[0].res/higgs_NLO_gg_full[0].res << "; increase is "  << higgs_NLO_gg_powers[i].res << " (+/-"<< higgs_NLO_gg_powers[i].err << ") pb" << endl;
						}
			}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qg channel" << endl;
			/////////////////////////////////////
			output << "Fractional qg xsec (at NLO):			" <<  res_higgs_NLO_qg_full[0].res << " pb +/-"  <<  res_higgs_NLO_qg_full[0].err  <<  endl;
			if(doexpansions){
				res_higgs_NLO_qg_accum.push_back({0,0,0});
				for(unsigned int i = 0; i < higgs_NLO_qg_powers.size(); i++){
						res_higgs_NLO_qg_accum[0].res = res_higgs_NLO_qg_accum[0].res+higgs_NLO_qg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NLO_qg_accum[0].res << " pb, fractional: "<< res_higgs_NLO_qg_accum[0].res/res_higgs_NLO_qg_full[0].res << "; increase is "  << higgs_NLO_qg_powers[i].res << " (+/-"<< higgs_NLO_qg_powers[i].err << ") pb" << endl;
					}
				}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qqbar channel" << endl;
			/////////////////////////////////////
			output << "Fractional qqbar xsec (at NLO):			" <<  res_higgs_NLO_qqbar_full[0].res << " pb +/-"  <<  res_higgs_NLO_qqbar_full[0].err  <<  endl;
			if(doexpansions){
				res_higgs_NLO_qqbar_accum.push_back({0,0,0});
				for(unsigned int i = 0; i < higgs_NLO_qqbar_powers.size(); i++){
						res_higgs_NLO_qqbar_accum[0].res = res_higgs_NLO_qqbar_accum[0].res+higgs_NLO_qqbar_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NLO_qqbar_accum[0].res << " pb, fractional: "<< res_higgs_NLO_qqbar_accum[0].res/res_higgs_NLO_qqbar_full[0].res << "; increase is "  << higgs_NLO_qqbar_powers[i].res << " (+/-"<< higgs_NLO_qqbar_powers[i].err << ") pb" << endl;
					}
			}
			output << "---------" << endl;
			output << "......................................................." << endl;
			output << "Total (LO+NLO): " << higgs_LO_gg_full[0].res + higgs_NLO_gg_full[0].res + res_higgs_NLO_qg_full[0].res + res_higgs_NLO_qqbar_full[0].res << endl;
			output << "Only gg channel (LO+NLO): " << higgs_LO_gg_full[0].res + higgs_NLO_gg_full[0].res << endl;
			output << "......................................................." << endl;
			////////////////////////////////////////////////////////////////////
		}
		if(NNLO&&higgs){
			cout << "NNLO" << endl;
			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "NNLO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////
			output << "gg channel" << endl;
			/////////////////////////////////////
			output << "Z-dependent, including LP (NNLO):		" <<  res_higgs_NNLO_gg_hard[0].res << " pb +/- " << res_higgs_NNLO_gg_hard[0].err <<  endl;
			output << "Constant piece (NNLO):				" <<  res_higgs_NNLO_gg_delta[0].res << " pb +/- " << res_higgs_NNLO_gg_delta[0].err <<  endl;
			output << "muR-dependent piece (NNLO): " << higgs_scale_factor(higgs_LO_gg_full[0].res,higgs_NLO_gg_full[0].res)<< endl;
			output << "Fractional gg xsec (at NNLO):			" <<  higgs_NNLO_gg_full[0].res << " pb "  <<  endl;
			if(doexpansions){
				res_higgs_NNLO_gg_accum.push_back(res_higgs_NNLO_gg_LP[0]);
				output << "	power "<<0<<" : " << res_higgs_NNLO_gg_accum[0].res << " pb, fractional: "<< res_higgs_NNLO_gg_accum[0].res/higgs_NNLO_gg_full[0].res << endl;
				for(unsigned int i = 0; i < higgs_NNLO_gg_powers.size(); i++){
						res_higgs_NNLO_gg_accum[0].res = res_higgs_NNLO_gg_accum[0].res+higgs_NNLO_gg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_gg_accum[0].res << " pb, fractional: "<< res_higgs_NNLO_gg_accum[0].res/higgs_NNLO_gg_full[0].res << "; increase is "  << higgs_NNLO_gg_powers[i].res << " (+/-"<< higgs_NNLO_gg_powers[i].err << ") pb" << endl;
					}
				}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qg channel" << endl;
			/////////////////////////////////////
			output << "muR-dependent piece (NNLO): " << higgs_scale_factor(0,res_higgs_NLO_qg_full[0].res)<< endl;
			output << "Fractional qqbar xsec (at NNLO):		" <<  res_higgs_NNLO_qg_full[0].res << " pb +/-"  <<  res_higgs_NNLO_qg_full[0].err  <<  endl;
			if(doexpansions){
				res_higgs_NNLO_qg_accum.push_back({0,0,0});
				for(unsigned int i = 0; i < higgs_NNLO_qg_powers.size(); i++){
						res_higgs_NNLO_qg_accum[0].res = res_higgs_NNLO_qg_accum[0].res+higgs_NNLO_qg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_qg_accum[0].res << " pb, fractional: "<< res_higgs_NNLO_qg_accum[0].res/res_higgs_NNLO_qg_full[0].res << "; increase is "  << higgs_NNLO_qg_powers[i].res << " (+/-"<< higgs_NNLO_qg_powers[i].err << ") pb" << endl;
					}
				}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qqbar channel" << endl;
			/////////////////////////////////////
			output << "muR-dependent piece (NNLO): " << higgs_scale_factor(0,res_higgs_NLO_qqbar_full[0].res)<< endl;
			output << "Fractional qqbar xsec (at NNLO):		" <<  res_higgs_NNLO_qqbar_full[0].res << " pb +/-"  <<  res_higgs_NNLO_qqbar_full[0].err  <<  endl;
			if(doexpansions){
				res_higgs_NNLO_qqbar_accum.push_back({0,0,0});
				for(unsigned int i = 0; i < higgs_NNLO_qqbar_powers.size(); i++){
						res_higgs_NNLO_qqbar_accum[0].res = res_higgs_NNLO_qqbar_accum[0].res+higgs_NNLO_qqbar_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_qqbar_accum[0].res << " pb, fractional: "<< res_higgs_NNLO_qqbar_accum[0].res/res_higgs_NNLO_qqbar_full[0].res << "; increase is "  << higgs_NNLO_qqbar_powers[i].res << " (+/-"<< higgs_NNLO_qqbar_powers[i].err << ") pb" << endl;
					}
				}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qq channel" << endl;
			/////////////////////////////////////
			output << "Fractional qq xsec (at NNLO):			" <<  res_higgs_NNLO_qq_full[0].res << " pb +/-"  <<  res_higgs_NNLO_qq_full[0].err  <<  endl;
			if(doexpansions){
				res_higgs_NNLO_qq_accum.push_back({0,0,0});
				for(unsigned int i = 0; i < higgs_NNLO_qq_powers.size(); i++){
						res_higgs_NNLO_qq_accum[0].res = res_higgs_NNLO_qq_accum[0].res+higgs_NNLO_qq_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_qq_accum[0].res << " pb, fractional: "<< res_higgs_NNLO_qq_accum[0].res/res_higgs_NNLO_qq_full[0].res << "; increase is "  << higgs_NNLO_qq_powers[i].res << " (+/-"<< higgs_NNLO_qq_powers[i].err << ") pb" << endl;
					}
				}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qq' channel" << endl;
			/////////////////////////////////////
			output << "Fractional qq' xsec (at NNLO):			" <<  res_higgs_NNLO_qqp_full[0].res << " pb +/-"  <<  res_higgs_NNLO_qqp_full[0].err  <<  endl;
			if(doexpansions){
				res_higgs_NNLO_qqp_accum.push_back({0,0,0});
				for(unsigned int i = 0; i < higgs_NNLO_qqp_powers.size(); i++){
						res_higgs_NNLO_qqp_accum[0].res = res_higgs_NNLO_qqp_accum[0].res+higgs_NNLO_qqp_powers[i].res;
						output << "	power "<<i+1<<" : " << res_higgs_NNLO_qqp_accum[0].res << " pb, fractional: "<< res_higgs_NNLO_qqp_accum[0].res/res_higgs_NNLO_qqp_full[0].res << "; increase is "  << higgs_NNLO_qqp_powers[i].res << " (+/-"<< higgs_NNLO_qqp_powers[i].err << ") pb" << endl;
					}
				}
			output << "---------" << endl;
			output << "......................................................." << endl;
			output << "Total (LO+NLO+NNLO): " << higgs_LO_gg_full[0].res + higgs_NLO_gg_full[0].res + res_higgs_NLO_qg_full[0].res + res_higgs_NLO_qqbar_full[0].res + higgs_NNLO_gg_full[0].res + res_higgs_NNLO_qg_full[0].res + res_higgs_NNLO_qqbar_full[0].res+ res_higgs_NNLO_qq_full[0].res+ res_higgs_NNLO_qqp_full[0].res << " pb" << endl;
			output << "Only gg channel (LO+NLO+NNLO): " << higgs_LO_gg_full[0].res + higgs_NLO_gg_full[0].res+ higgs_NNLO_gg_full[0].res << " pb" << endl;
			output << "......................................................." << endl;
		}
		if(RES&&higgs){
			/*cout << "My code: " << endl;
			cout << "LO:	" << higgs_LO_gg_full[0].res << endl;
			cout << "NLO:	" << higgs_LO_gg_full[0].res+higgs_NLO_gg_full[0].res << endl;
			cout << "NNLO:	" << higgs_LO_gg_full[0].res+higgs_NLO_gg_full[0].res+higgs_NNLO_gg_full[0].res << endl;
			cout << "LL_NLO: " <<  resummed_higgs_LP_LL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_LL_exp_NLO[0].err <<  endl;
			cout << "LL_NNLO: " <<  resummed_higgs_LP_LL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_LL_exp_NNLO[0].err <<  endl;
			cout << "NLL_NLO: " <<  resummed_higgs_LP_NLL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NLL_exp_NLO[0].err <<  endl;
			cout << "NLL_NNLO: " <<  resummed_higgs_LP_NLL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NLL_exp_NNLO[0].err <<  endl;
			cout << "NNLL_NLO: " <<  resummed_higgs_LP_NNLL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NNLL_exp_NLO[0].err <<  endl;
			cout << "NNLL_NNLO: " <<  resummed_higgs_LP_NNLL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NNLL_exp_NNLO[0].err <<  endl;
			cout << "NNNLL_NLO: " <<  resummed_higgs_LP_NNNLL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NNNLL_exp_NLO[0].err <<  endl;
			cout << "NNNLL_NNLO: " <<  resummed_higgs_LP_NNNLL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NNNLL_exp_NNLO[0].err <<  endl;
			*/
			/////////////////////////////////////
			///
			output << "=======================" << endl;
			output << "Resummed results (pQCD)" << endl;
			output << "=======================" << endl;
			///
			/////////////////////////////////////
			output << "Resummed (LP NNNLL + NLP LL): 		" << resummed_higgs_LP_NNNLL_NLP_LL[0].res << " pb +/- " << resummed_higgs_LP_NNNLL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NNNLL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NNNLL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NNNLL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NNNLL_NLP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;


			output << "Resummed (LP NNLL + NLP LL): 		" << resummed_higgs_LP_NNLL_NLP_LL[0].res << " pb +/- " << resummed_higgs_LP_NNLL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NNLL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NNLL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NNLL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NNLL_NLP_LL_exp_NNLO[0].err <<  endl;
			for(unsigned int i = 0; i < resummed_higgs_NLP_exp_NNLL.size(); i++){
						output << " NLP(NNLL) expanded to order as "<<i<<" : " << resummed_higgs_NLP_exp_NNLL[i].res << endl;
					}
			output << "--------------------------------" << endl;

			output << "Resummed (LP NLL + NLP LL): 		" << resummed_higgs_LP_NLL_NLP_LL[0].res << " pb +/- " << resummed_higgs_LP_NLL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NLL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NLL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NLL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NLL_NLP_LL_exp_NNLO[0].err <<  endl;
			for(unsigned int i = 0; i < resummed_higgs_NLP_exp_NLL.size(); i++){
						output << "	NLP(NLL) expanded to order as "<<i<<" : " << resummed_higgs_NLP_exp_NLL[i].res << endl;
					}
			output << "--------------------------------" << endl;
			output << "Resummed (LP NLLp + NLP LL): 		" << resummed_higgs_LP_NLLp_NLP_LL[0].res << " pb +/- " << resummed_higgs_LP_NLLp_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO_prime:				" <<  resummed_higgs_LP_NLLp_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NLLp_NLP_LL_exp_NLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP LL + NLP LL): 		" << resummed_higgs_LP_LL_NLP_LL[0].res << " pb +/- " << resummed_higgs_LP_LL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_LL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_LL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_LL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_LL_NLP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NNNLL): 		" << resummed_higgs_LP_NNNLL[0].res << " pb +/- " << resummed_higgs_LP_NNNLL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NNNLL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NNNLL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NNNLL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NNNLL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NNLL): 		" << resummed_higgs_LP_NNLL[0].res << " pb +/- " << resummed_higgs_LP_NNLL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NNLL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NNLL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NNLL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NNLL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NLL): 		" << resummed_higgs_LP_NLL[0].res << " pb +/- " << resummed_higgs_LP_NLL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NLL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NLL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NLL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NLL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NLLp): 		" << resummed_higgs_LP_NLLp[0].res << " pb +/- " << resummed_higgs_LP_NLLp[0].err <<  endl;
			output << "Expanded to NLO_prime:				" <<  resummed_higgs_LP_NLLp_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NLLp_exp_NLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP LL): 		" << resummed_higgs_LP_LL[0].res << " pb +/- " << resummed_higgs_LP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_LL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_LL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "NNLL NLP qg: 		" << resummed_higgs_NLP_qg[3].res << " pb +/- " << resummed_higgs_NLP_qg[3].err << endl;
			output << "NLL NLP qg: 		" << resummed_higgs_NLP_qg[2].res << " pb +/- " << resummed_higgs_NLP_qg[2].err << endl;

			for(unsigned int i = 0; i < resummed_higgs_NLP_qg.size(); i++){
						output << "	NLP qg expanded to order as "<<i<<" : " << resummed_higgs_NLP_qg[i].res << endl;
					}
			output << "======================================================" << endl;
		}
	}
}
	if(DY){
		output << "======================================================" << endl;
		output << "DY results" << endl;
		output << "======================================================" << endl;
		if(PDFmemuse)
		{
			output << "LO, NLOqqbar, NLOqg, NNLOqqbar, NNLOqg, NNLOgg, NNLOqq \n" << endl;
			for(unsigned int i = 0; i < DY_LO_qqbar_full.size(); i++){
			output << "	member "<<i<<" : " <<  DY_LO_qqbar_full[i].res << "," <<  DY_NLO_qqbar_full[i].res << "," <<  res_DY_NLO_qg_full[i].res << ","
		     <<  DY_NNLO_qqbar_full[i].res << ","  <<  res_DY_NNLO_qg_full[i].res << ","
			 <<  res_DY_NNLO_gg_full[i].res << ","  <<  res_DY_NNLO_qq_full[i].res << endl;

			}
		}
		else{
		if(LO&&DY){

			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "LO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////

			output << "---------------------------------------" << endl;

			output << "Total xsec:					" << DY_LO_qqbar_full[0].res << " pb +/- " << DY_LO_qqbar_full[0].err <<  endl;
		}
		if(NLO&&DY){

			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "NLO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////
			output << "qqbar channel" << endl;
			/////////////////////////////////////
			output << "Z-dependent, including LP (NLO):		" <<  res_DY_NLO_qqbar_hard[0].res + res_DY_NLO_qqbar_LP_cor[0].res << " pb +/- " << res_DY_NLO_qqbar_hard[0].err + res_DY_NLO_qqbar_LP_cor[0].err <<  endl;
			output << "Constant piece (NLO):				" <<  res_DY_NLO_qqbar_delta[0].res << " pb +/- " << res_DY_NLO_qqbar_delta[0].err <<  endl;
			output << "Fractional qqbar xsec (at NLO):			" <<  DY_NLO_qqbar_full[0].res << " pb " <<  endl;
			if(doexpansions){
				res_DY_NLO_qqbar_accum.push_back(res_DY_NLO_qqbar_LP[0]);
				output << "	power "<<0<<" : " << res_DY_NLO_qqbar_accum[0].res << " pb, fractional: "<< res_DY_NLO_qqbar_accum[0].res/DY_NLO_qqbar_full[0].res << endl;
				for(unsigned int i = 0; i < DY_NLO_qqbar_powers.size(); i++){
						res_DY_NLO_qqbar_accum[0].res = res_DY_NLO_qqbar_accum[0].res+DY_NLO_qqbar_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NLO_qqbar_accum[0].res << " pb, fractional: "<< res_DY_NLO_qqbar_accum[0].res/DY_NLO_qqbar_full[0].res << "; increase is "  << DY_NLO_qqbar_powers[i].res << " (+/-"<< DY_NLO_qqbar_powers[i].err << ") pb" << endl;
					}
				}
			output << "---------" << endl;
		/////////////////////////////////////
			output << "qg channel" << endl;
			/////////////////////////////////////
			output << "Fractional qg xsec (at NLO):			" <<  res_DY_NLO_qg_full[0].res << " pb +/-"  <<  res_DY_NLO_qg_full[0].err  <<  endl;
			if(doexpansions){
				res_DY_NLO_qg_accum.push_back({0,0,0});
				for(unsigned int i = 0; i < DY_NLO_qg_powers.size(); i++){
						res_DY_NLO_qg_accum[0].res = res_DY_NLO_qg_accum[0].res+DY_NLO_qg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NLO_qg_accum[0].res << " pb, fractional: "<< res_DY_NLO_qg_accum[0].res/res_DY_NLO_qg_full[0].res << "; increase is "  << DY_NLO_qg_powers[i].res << " (+/-"<< DY_NLO_qg_powers[i].err << ") pb" << endl;
					}
				}
			output << "---------" << endl;
			output << "......................................................." << endl;
			output << "Total (LO+NLO): " << DY_LO_qqbar_full[0].res + DY_NLO_qqbar_full[0].res + res_DY_NLO_qg_full[0].res << " pb" << endl;
			output << "Only qqbar channel (LO+NLO): " << DY_LO_qqbar_full[0].res + DY_NLO_qqbar_full[0].res << " pb" << endl;
			output << "......................................................." << endl;
			////////////////////////////////////////////////////////////////////
		}
		if(NNLO&&DY){
			/////////////////////////////////////
			///
			output << "===================" << endl;
			output << "NNLO results" << endl;
			output << "===================" << endl;
			///
			/////////////////////////////////////
			output << "qqbar channel" << endl;
			/////////////////////////////////////
			output << "Z-dependent, including LP (NNLO):		" <<  res_DY_NNLO_qqbar_hard[0].res + res_DY_NNLO_qqbar_LP_cor[0].res<< " pb +/- " << res_DY_NNLO_qqbar_hard[0].err + res_DY_NNLO_qqbar_LP_cor[0].err <<  endl;
			output << "Constant piece (NNLO):				" <<  res_DY_NNLO_qqbar_delta[0].res << " pb +/- " << res_DY_NNLO_qqbar_delta[0].err <<  endl;
			output << "Fractional qqbar xsec (at NNLO):			" <<  DY_NNLO_qqbar_full[0].res << " pb "  <<  endl;
			if(doexpansions){
				res_DY_NNLO_qqbar_accum.push_back(res_DY_NNLO_qqbar_LP[0]);
				output << "	power "<<0<<" : " << res_DY_NNLO_qqbar_accum[0].res << " pb, fractional: "<< res_DY_NNLO_qqbar_accum[0].res/DY_NNLO_qqbar_full[0].res << endl;
				for(unsigned int i = 0; i < DY_NNLO_qqbar_powers.size(); i++){
						res_DY_NNLO_qqbar_accum[0].res = res_DY_NNLO_qqbar_accum[0].res+DY_NNLO_qqbar_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NNLO_qqbar_accum[0].res << " pb, fractional: "<< res_DY_NNLO_qqbar_accum[0].res/DY_NNLO_qqbar_full[0].res << "; increase is "  << DY_NNLO_qqbar_powers[i].res << " (+/-"<< DY_NNLO_qqbar_powers[i].err << ") pb" << endl;
					}
				}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qg channel" << endl;
			/////////////////////////////////////
			output << "Fractional qg xsec (at NNLO):			" <<  res_DY_NNLO_qg_full[0].res << " pb +/-"  <<  res_DY_NNLO_qg_full[0].err  <<  endl;
			if(doexpansions){
				res_DY_NNLO_qg_accum.push_back({0,0,0});
				for(unsigned int i = 0; i < DY_NNLO_qg_powers.size(); i++){
						res_DY_NNLO_qg_accum[0].res = res_DY_NNLO_qg_accum[0].res+DY_NNLO_qg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NNLO_qg_accum[0].res << " pb, fractional: "<< res_DY_NNLO_qg_accum[0].res/res_DY_NNLO_qg_full[0].res << "; increase is "  << DY_NNLO_qg_powers[i].res << " (+/-"<< DY_NNLO_qg_powers[i].err << ") pb" << endl;
					}
				}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "gg channel" << endl;
			/////////////////////////////////////
			output << "Fractional gg xsec (at NNLO):		" <<  res_DY_NNLO_gg_full[0].res << " pb +/-"  <<  res_DY_NNLO_gg_full[0].err  <<  endl;
			res_DY_NNLO_gg_accum.push_back({0,0,0});
			for(unsigned int i = 0; i < DY_NNLO_gg_powers.size(); i++){
						res_DY_NNLO_gg_accum[0].res = res_DY_NNLO_gg_accum[0].res+DY_NNLO_gg_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NNLO_gg_accum[0].res << " pb, fractional: "<< res_DY_NNLO_gg_accum[0].res/res_DY_NNLO_gg_full[0].res << "; increase is "  << DY_NNLO_gg_powers[i].res << " (+/-"<< DY_NNLO_gg_powers[i].err << ") pb" << endl;
					}
			output << "---------" << endl;
			/////////////////////////////////////
			output << "qq channel" << endl;
			/////////////////////////////////////
			output << "Fractional qq xsec (at NNLO):			" <<  res_DY_NNLO_qq_full[0].res << " pb +/-"  <<  res_DY_NNLO_qq_full[0].err  <<  endl;
			if(doexpansions){
				res_DY_NNLO_qq_accum.push_back({0,0,0});
				for(unsigned int i = 0; i < DY_NNLO_qq_powers.size(); i++){
						res_DY_NNLO_qq_accum[0].res = res_DY_NNLO_qq_accum[0].res+DY_NNLO_qq_powers[i].res;
						output << "	power "<<i+1<<" : " << res_DY_NNLO_qq_accum[0].res << " pb, fractional: "<< res_DY_NNLO_qq_accum[0].res/res_DY_NNLO_qq_full[0].res << "; increase is "  << DY_NNLO_qq_powers[i].res << " (+/-"<< DY_NNLO_qq_powers[i].err << ") pb" << endl;
					}
				}
			output << "---------" << endl;
			output << "......................................................." << endl;
			output << "Total (LO+NLO+NNLO): " << DY_LO_qqbar_full[0].res + DY_NLO_qqbar_full[0].res + res_DY_NLO_qg_full[0].res + DY_NNLO_qqbar_full[0].res + + res_DY_NNLO_gg_full[0].res + res_DY_NNLO_qg_full[0].res+ res_DY_NNLO_qq_full[0].res << " pb" << endl;
			output << "Only qqbar channel (LO+NLO+NNLO): " << DY_LO_qqbar_full[0].res + DY_NLO_qqbar_full[0].res + DY_NNLO_qqbar_full[0].res << " pb" << endl;
			output << "......................................................." << endl;
		}
		if(RES&&DY){

			/////////////////////////////////////
			///
			output << "=======================" << endl;
			output << "Resummed results (pQCD)" << endl;
			output << "=======================" << endl;
			///
			/////////////////////////////////////
			output << "Resummed (LP NNNLL + NLP LL): 		" << resummed_DY_LP_NNNLL_NLP_LL[0].res << " pb +/- " << resummed_DY_LP_NNNLL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NNNLL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_NNNLL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NNNLL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_NNNLL_NLP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NNLL + NLP LL): 		" << resummed_DY_LP_NNLL_NLP_LL[0].res << " pb +/- " << resummed_DY_LP_NNLL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NNLL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_NNLL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NNLL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_NNLL_NLP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			for(unsigned int i = 0; i < resummed_DY_NLP_exp_NNLL.size(); i++){
						output << "	NLP(NNLL) expanded to order as "<<i<<" : " << resummed_DY_NLP_exp_NNLL[i].res << endl;
					}
			output << "Resummed (LP NLL + NLP LL): 		" << resummed_DY_LP_NLL_NLP_LL[0].res << " pb +/- " << resummed_DY_LP_NLL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NLL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_NLL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NLL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_NLL_NLP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;
			for(unsigned int i = 0; i < resummed_DY_NLP_exp_NLL.size(); i++){
						output << "	NLP(NLL) expanded to order as "<<i<<" : " << resummed_DY_NLP_exp_NLL[i].res << endl;
					}
			output << "--------------------------------" << endl;
			output << "Resummed (LP NLLp + NLP LL): 		" << resummed_DY_LP_NLLp_NLP_LL[0].res << " pb +/- " << resummed_DY_LP_NLLp_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO_prime:				" <<  resummed_DY_LP_NLLp_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_NLLp_NLP_LL_exp_NLO[0].err <<  endl;
			output << "--------------------------------" << endl;
			output << "Resummed (LP LL + NLP LL): 		" << resummed_DY_LP_LL_NLP_LL[0].res << " pb +/- " << resummed_DY_LP_LL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_LL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_LL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_LL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_LL_NLP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NNNLL): 		" << resummed_DY_LP_NNNLL[0].res << " pb +/- " << resummed_DY_LP_NNNLL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NNNLL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_NNNLL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NNNLL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_NNNLL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NNLL): 		" << resummed_DY_LP_NNLL[0].res << " pb +/- " << resummed_DY_LP_NNLL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NNLL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_NNLL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NNLL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_NNLL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NLL): 		" << resummed_DY_LP_NLL[0].res << " pb +/- " << resummed_DY_LP_NLL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NLL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_NLL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NLL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_NLL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP NLLp): 		" << resummed_DY_LP_NLLp[0].res << " pb +/- " << resummed_DY_LP_NLLp[0].err <<  endl;
			output << "Expanded to NLO_prime:				" <<  resummed_DY_LP_NLLp_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_NLLp_exp_NLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "Resummed (LP LL): 		" << resummed_DY_LP_LL[0].res << " pb +/- " << resummed_DY_LP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_LL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_LL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "NNLL NLP qg: 		" << resummed_DY_NLP_qg[3].res << " pb +/- " << resummed_DY_NLP_qg[3].err << endl;
			output << "NLL NLP qg: 		" << resummed_DY_NLP_qg[2].res << " pb +/- " << resummed_DY_NLP_qg[2].err << endl;

			for(unsigned int i = 0; i < resummed_DY_NLP_qg.size(); i++){
						output << "	NLP qg expanded to order as "<<i<<" : " << resummed_DY_NLP_qg[i].res << endl;
					}
			output << "======================================================" << endl;
		}
		}
}
	cout << Q << " " << Q2 << " " << muR << " " << muR2 << " " << muF << " " << muF2 << endl;
	output.close();
	return 0;
}
