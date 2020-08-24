#include <cmath>
#include <iostream>
#include <fstream>
#include "cuba.h"
#include "cuba_integration.h"
#include <sstream>
#include "deriv_pdf.h"
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

	cout.precision(16);
	/*complex<double> N = 2.1;
	Q2 = 500.*500.;
	muF2 = 500.*500./4.;
	muR2 = 500.*500./4.;
	cout << "/.Lt->"<< log(Q2/mt2) << "/.as->" << alphas_muR << "/.lnN->" << log(N*exp(M_gammaE)) << "/.x->0.1/.CA->3/.CF->4/3/.LQmuF2->"<< log(Q2/muF2) << endl;
	cout << "DY NLO: " << DY_NLO_qqbar_delta() << " " << DY_NLO_qqbar_plus(0.1) << endl;
	cout << "DY NNLO: " << DY_NNLO_qqbar_delta() << " " << DY_NNLO_qqbar_plus(0.1) << endl;
	cout << "DY g01 g02 " << DY_g01() << " " << DY_g02() << endl;
	cout << "Higgs NLO: " << higgs_NLO_gg_delta() << " " << higgs_NLO_gg_plus(0.1) << endl;
	cout << "Higgs NNLO: " << higgs_NNLO_gg_delta() << " " << higgs_NNLO_gg_plus(0.1) << endl;
	cout << "Higgs g01 g02 " << higgs_g01() << " " << higgs_g02() << endl;

	cout << "NNLO match at LP NNLL " << endl;
	ISLL = 1;
	ISNLL = 1;
	ISNNLL = 1;
	ISNLP = 0;
	cout << "DY " << NNLOmatch(N, A1q, A2q, D2DY, DY_g01(), DY_g02()) << endl;
	cout << "Higgs " << NNLOmatch(N, A1g, A2g, D2higgs, higgs_g01(), higgs_g02()) << endl;
	cout << "NNLO match at NLP NNLL " << endl;
	ISLL = 1;
	ISNLL = 1;
	ISNNLL = 1;
	ISNLP = 1;
	cout << "DY " << NNLOmatch(N, A1q, A2q, D2DY, DY_g01(), DY_g02()) << endl;
	cout << "Higgs " << NNLOmatch(N, A1g, A2g, D2higgs, higgs_g01(), higgs_g02()) << endl;

	exit(0);*/
	//////////////////////
	/// TEST FUNCTIONS
	//////////////////////
	//cout << "REMEMBER THAT M2 was introduced for tth!!!! " << endl;
	bool TEST = false;
	bool MAKECOEF = false;
	bool doexpansions = false;
	if(MAKECOEF){
		Q = 125.;
		muF = 125.;
		muR = 125.;
		int Nsteps = 10000;
		ofstream output2;
		string q_str2 = "COEFFICIENTS_DY.txt";
		output2.open(q_str2.c_str()); //.c_str() needed to do a constant string conversion
		cout << q_str2 << endl;

		update_defaults();
		output2 << "x NLOqqbar NLOqg NNLOqqbarNS NNLOqg NNLOgg NNLOB2 NNLOBC NNLOC2 NNLOCD NNLOCE NNLOCF" << endl;
		for(int i=1;i<Nsteps;i++)
		{
			double x = double(float(i)/double(Nsteps));
			output2 << x << " ";
			output2 << /*M_PI/(alphas_muR)**/(DY_NLO_qqbar_reg(x)+DY_NLO_qqbar_plus(x)) << "," << /*M_PI/(alphas_muR)**/(DY_NLO_qqbar_plus(x)) << "," << /*M_PI/(alphas_muR)**/(DY_NLO_qqbar_plus(x)+DY_NLO_qqbar_expansion(x, 1)) << "," << /*M_PI/(alphas_muR)**/(DY_NLO_qqbar_plus(x)+DY_NLO_qqbar_expansion(x, 1)+DY_NLO_qqbar_expansion(x, 2)) << " ";
			output2 << /*M_PI/(alphas_muR)**/(DY_NLO_qg_full(x)) << "," << /*M_PI/(alphas_muR)**/(0) << "," << /*M_PI/(alphas_muR)**/(DY_NLO_qg_expansion(x, 1)) << "," << /*M_PI/(alphas_muR)**/(DY_NLO_qg_expansion(x, 1)+DY_NLO_qg_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qqbar_NS(x)+DY_NNLO_qqbar_plus(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qqbar_plus(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qqbar_plus(x)+DY_NNLO_qqbar_NS_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qqbar_plus(x)+DY_NNLO_qqbar_NS_expansion(x, 1)+DY_NNLO_qqbar_NS_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qg_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qg_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_qg_expansion(x, 1)+DY_NNLO_qg_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_gg_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_gg_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_gg_expansion(x, 1)+DY_NNLO_gg_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_BB_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_BB_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_BB_expansion(x, 1)+DY_NNLO_BB_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_BC_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_BC_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_BC_expansion(x, 1)+DY_NNLO_BC_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CC_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CC_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CC_expansion(x, 1)+DY_NNLO_CC_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CD_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CD_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CD_expansion(x, 1)+DY_NNLO_CD_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CE_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CE_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CE_expansion(x, 1)+DY_NNLO_CE_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CF_full(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CF_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(DY_NNLO_CF_expansion(x, 1)+DY_NNLO_CF_expansion(x, 2)) << " ";
			output2 << endl;
		}
		output2.close();

		q_str2 = "COEFFICIENTS_HIGGS.txt";
		output2.open(q_str2.c_str()); //.c_str() needed to do a constant string conversion
		cout << q_str2 << endl;

		output2 << "x NLOgg NLOqg NLOqqbar NNLOgg NNLOqg NNLOqq NNLOqqp NNLOqqbar" << endl;
		for(int i=1;i<Nsteps;i++)
		{
			double x = double(float(i)/double(Nsteps));
			output2 << x << " ";
			output2 << /*M_PI/(alphas_muR)**/(higgs_NLO_gg_reg(x)+higgs_NLO_gg_plus(x)) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_gg_plus(x)) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_gg_plus(x)+higgs_NLO_gg_expansion(x, 1)) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_gg_plus(x)+higgs_NLO_gg_expansion(x, 1)+higgs_NLO_gg_expansion(x, 2)) << " ";
			output2 << /*M_PI/(alphas_muR)**/(higgs_NLO_qg_full(x)) << "," << /*M_PI/(alphas_muR)**/(0) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_qg_expansion(x, 1)) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_qg_expansion(x, 1)+higgs_NLO_qg_expansion(x, 2)) << " ";
			output2 << /*M_PI/(alphas_muR)**/(higgs_NLO_qqbar_full(x)) << "," << /*M_PI/(alphas_muR)**/(0) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_qqbar_expansion(x, 1)) << "," << /*M_PI/(alphas_muR)**/(higgs_NLO_qqbar_expansion(x, 1)+higgs_NLO_qqbar_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_gg_reg(x)+higgs_NNLO_gg_plus(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_gg_plus(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_gg_plus(x)+higgs_NNLO_gg_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_gg_plus(x)+higgs_NNLO_gg_expansion(x, 1)+higgs_NNLO_gg_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qg_reg(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qg_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qg_expansion(x, 1)+higgs_NNLO_qg_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qq_reg(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qq_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qq_expansion(x, 1)+higgs_NNLO_qq_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qqp_reg(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qqp_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qqp_expansion(x, 1)+higgs_NNLO_qqp_expansion(x, 2)) << " ";
			output2 << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qqbar_reg(x)) << "," << /*pow(M_PI/(alphas_muR),2)**/(0) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qqbar_expansion(x, 1)) << "," << /*pow(M_PI/(alphas_muR),2)**/(higgs_NNLO_qqbar_expansion(x, 1)+higgs_NNLO_qqbar_expansion(x, 2)) << " ";
			output2 << endl;
		}
		output2.close();
		exit(0);
	}
	if(TEST){
		double z = 0.5;
		ofstream output2;
		string q_str2 = "COMPARISON_CODES.txt";
		output2.open(q_str2.c_str()); //.c_str() needed to do a constant string conversion
		cout << q_str2 << endl;
		output2 << "============================================================" << endl;
		output2 << "DY coefficients" << endl;
		output2 << "============================================================" << endl;
		output2 << "z=" << z <<", Q=" << Q << " GeV, muF=" << muF << " GeV, muR=" << muR << " GeV, tau=" << tau << ", alphas(muR2)=" << alphas_muR << endl;
		output2 << "LO factor (tau*sigma0): " << DY_LO_factor() << endl;
		output2 << "===============NLO-qqbar==========================" << endl;
		output2 << "NLO qqbar LP (*alphas/(4pi)): " << DY_LO_factor()*DY_NLO_qqbar_plus(z) << endl;
		output2 << "NLO qqbar z-regular: " << DY_LO_factor()*DY_NLO_qqbar_reg(z) << endl;
		output2 << "NLO qqbar delta: " << DY_LO_factor()*DY_NLO_qqbar_delta() << endl;

		output2 << "===============NLO-qg==========================" << endl;
		output2 << "NLO qg full: " << DY_LO_factor()*DY_NLO_qg_full(z) << endl;

		output2 << "===============NNLO-qqbar/qq/qbarqbar==========================" << endl;
		output2 << "NNLO qqbar LP (*alphas/(4pi))**2: " << DY_LO_factor()*DY_NNLO_qqbar_plus(z) << endl;
		output2 << "NNLO qqbar NS: " << DY_LO_factor()*DY_NNLO_qqbar_NS(z) << endl;
		output2 << "NNLO qqbar delta: " << DY_LO_factor()*DY_NNLO_qqbar_delta() << endl;
		output2 << "NNLO delta BB: " << DY_LO_factor()*DY_NNLO_BB_full(z) << endl;
		output2 << "NNLO delta BC: " << DY_LO_factor()*DY_NNLO_BC_full(z) << endl;
		output2 << "NNLO delta CC: " << DY_LO_factor()*DY_NNLO_CC_full(z) << endl;
		output2 << "NNLO delta CD: " << DY_LO_factor()*DY_NNLO_CD_full(z) << endl;
		output2 << "NNLO delta CE: " << DY_LO_factor()*DY_NNLO_CE_full(z) << endl;
		output2 << "NNLO delta CF: " << DY_LO_factor()*DY_NNLO_CF_full(z) << endl;
		output2 << "===============NNLO-qg==========================" << endl;
		output2 << "NNLO qg full: " << DY_LO_factor()*DY_NNLO_qg_full(z) << endl;
		output2 << "===============NNLO-gg==========================" << endl;
		output2 << "NNLO gg full: " << DY_LO_factor()*DY_NNLO_gg_full(z) << endl;

		output2 << "============================================================" << endl;
		output2 << "PDF convultions - weighted with charges" << endl;
		output2 << "============================================================" << endl;
		output2 << "Sum_i(Q_q[i]Q_q[i]*1/x*q_i(z)*qbar_i(tau/z) + <-> ), for DeltaNNLONS, DeltaNNLOBC: " << pdf_sum_qqbar_charge_weighted(z, tau) << endl;
		output2 << "Sum_i(Q_q[i]Q_q[i])*1/x*Sum_i(q_i(z)*qbar_i(tau/z) + <-> ), for DeltaNNLOB2: " << pdf_sum_qqbar_charge_unweighted(z, tau) << endl;
		output2 << "Sum_(i,j)(Q_q[i]Q_q[i]*1/x*q_i(z)*q_j(tau/z) + <-> + ...), for qq + qbarqbar +qqbar (identical+non-identical), for DeltaNNLOC2: " << pdf_sum_qq_charge_weighted_double(z, tau) << endl;
		output2 << "Sum_(i)(Q_q[i]Q_q[i]*1/x*q_i(z)*q_i(tau/z) + <-> + ...), for qq + qbarqbar +qqbar (identical), for DeltaNNLOCE: " << pdf_sum_qq_charge_weighted_single(z, tau) << endl;
		output2 << "Sum_(i,j)(Q_q[i]Q_q[j]*1/x*q_i(z)*q_j(tau/z) + ...), for qq + qbarqbar +qqbar (identical+NI), for DeltaNNLOCD: " << pdf_sum_qq_charge_weighted_double_vivj(z, tau) << endl;
		output2 << "Sum_(i)(Q_q[i]Q_q[i]*1/x*q_i(z)*q_i(tau/z) + ...), for qq + qbarqbar +qqbar (identical), for DeltaNNLOCF: " << pdf_sum_qq_charge_weighted_single_vivi(z, tau) << endl;
		output2 << "Sum_(i)(Q_q[i]Q_q[i]*1/x*q_i(z)*g(tau/z) + <-> + ...), for qg + qbarg: " << pdf_sum_qg_charge_weighted(z, tau) << endl;
		output2 << "Sum_(i)(Q_q[i]Q_q[i])*1/x*g(z)*g(tau/z)), for gg: " << pdf_sum_gg_charge_weighted(z, tau) << endl;


		output2 << "============================================================" << endl;
		output2 << "Higgs coefficients" << endl;
		output2 << "============================================================" << endl;
		output2 << "z=" << z <<", Q=" << Q << " GeV, muF=" << muF << " GeV, muR=" << muR << " GeV, tau=" << tau << ", alphas(muR2)=" << alphas_muR << endl;
		output2 << "LO factor (tau*sigma0): " << higgs_LO_factor() << endl;
		output2 << "===============NLO-gg==========================" << endl;
		output2 << "NLO gg LP (*alphas/(pi)): " << higgs_LO_factor()*higgs_NLO_gg_plus(z) << endl;
		output2 << "NLO gg z-regular: " << higgs_LO_factor()*higgs_NLO_gg_reg(z) << endl;
		output2 << "NLO gg delta: " << higgs_LO_factor()*higgs_NLO_gg_delta() << endl;

		output2 << "===============NLO-qg==========================" << endl;
		output2 << "NLO qg full: " << higgs_LO_factor()*higgs_NLO_qg_full(z) << endl;
		output2 << "===============NLO-qqbar==========================" << endl;
		output2 << "NLO qqbar full: " << higgs_LO_factor()*higgs_NLO_qqbar_full(z) << endl;


		output2 << "===============NNLO-gg==========================" << endl;
		output2 << "NNLO gg LP (*alphas/(pi))**2: " << higgs_LO_factor()*(higgs_NNLO_gg_plus(z)+logdep_gg_plus(z)) << endl;
		output2 << "NNLO gg z-regular: " <<  higgs_LO_factor()*(higgs_NNLO_gg_reg(z) + logdep_gg_reg(z)) << endl;
		output2 << "NNLO gg delta: " << higgs_LO_factor()*(higgs_NNLO_gg_delta() + logdep_gg_constant()) << endl;
		output2 << "===============NNLO-qg==========================" << endl;
		output2 << "NNLO qg full: " << higgs_LO_factor()*(higgs_NNLO_qg_reg(z) +logdep_qg(z)) << endl;
		output2 << "===============NNLO-qqbar==========================" << endl;
		output2 << "NNLO qqbar full: " << higgs_LO_factor()*(higgs_NNLO_qqbar_reg(z)+logdep_qqbar(z)) << endl;
		output2 << "===============NNLO-qq==========================" << endl;
		output2 << "NNLO qq full: " << higgs_LO_factor()*(higgs_NNLO_qq_reg(z)+logdep_qq(z)) << endl;
		output2 << "===============NNLO-qq'==========================" << endl;
		output2 << "NNLO qq' full: " << higgs_LO_factor()*(higgs_NNLO_qqp_reg(z)+logdep_qq(z)) << endl;


		output2 << "============================================================" << endl;
		output2 << "PDF convultions" << endl;
		output2 << "============================================================" << endl;
		output2 << "Sum_i(1/x*q_i(z)*qbar_i(tau/z) + <-> ), for qqbar (identical): " << pdf_sum_qqbar(z, tau) << endl;
		output2 << "Sum_(i)(1/x*q_i(z)*g(tau/z) + <->), for qg (and qbarg): " << pdf_sum_qg(z, tau) << endl;
		output2 << "1/x*g(z)*g(tau/z)), for gg: " << pdf_sum_gg(z, tau) << endl;
		output2 << "Sum_i(1/x*q_i(z)*q_i(tau/z)+qbar_i(z)*qbar_i(tau/z)), for qq+qbarqbar (identical) : " << pdf_sum_qq(z, tau) << endl;
		output2 << "Sum_(i,j)(1/x*q_i(z)*q_j(tau/z)+<->+qbar_i(z)*qbar_j(tau/z)+<->), for qq+qbarqbar (non-identical) : " << pdf_sum_qqNI(z, tau) << endl;

		output2.close();
		exit(0);
	}

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
	vector<results_c> resummed_DY_LP_NNNLL, resummed_DY_LP_NNNLL_NLP_LL, resummed_DY_LP_NNLL, resummed_DY_LP_NNLL_NLP_LL, resummed_DY_LP_NLL, resummed_DY_LP_NLL_NLP_LL, resummed_DY_LP_LL_NLP_LL, resummed_DY_LP_LL;
	vector<results_c> resummed_DY_LP_NNNLL_exp_NNLO, resummed_DY_LP_NNNLL_NLP_LL_exp_NNLO,resummed_DY_LP_NNLL_exp_NNLO, resummed_DY_LP_NNLL_NLP_LL_exp_NNLO, resummed_DY_LP_NLL_exp_NNLO, resummed_DY_LP_NLL_NLP_LL_exp_NNLO, resummed_DY_LP_LL_NLP_LL_exp_NNLO, resummed_DY_LP_LL_exp_NNLO;
	vector<results_c> resummed_DY_LP_NNNLL_exp_NLO, resummed_DY_LP_NNNLL_NLP_LL_exp_NLO,resummed_DY_LP_NNLL_exp_NLO, resummed_DY_LP_NNLL_NLP_LL_exp_NLO, resummed_DY_LP_NLL_exp_NLO, resummed_DY_LP_NLL_NLP_LL_exp_NLO, resummed_DY_LP_LL_NLP_LL_exp_NLO, resummed_DY_LP_LL_exp_NLO;
	vector<results_c> SCET_DY_LP_NNLL, SCET_DY_LP_NNLL_NLP_LL, SCET_DY_LP_NLL, SCET_DY_LP_NLL_NLP_LL, SCET_DY_LP_LL_NLP_LL, SCET_DY_LP_LL;

	// higgs
	vector<results_c> resummed_higgs_LP_NNNLL, resummed_higgs_LP_NNNLL_NLP_LL,resummed_higgs_LP_NNLL, resummed_higgs_LP_NNLL_NLP_LL, resummed_higgs_LP_NLL, resummed_higgs_LP_NLL_NLP_LL, resummed_higgs_LP_LL_NLP_LL, resummed_higgs_LP_LL;
	vector<results_c> resummed_higgs_LP_NNNLL_exp_NNLO, resummed_higgs_LP_NNNLL_NLP_LL_exp_NNLO,resummed_higgs_LP_NNLL_exp_NNLO, resummed_higgs_LP_NNLL_NLP_LL_exp_NNLO, resummed_higgs_LP_NLL_exp_NNLO, resummed_higgs_LP_NLL_NLP_LL_exp_NNLO, resummed_higgs_LP_LL_NLP_LL_exp_NNLO, resummed_higgs_LP_LL_exp_NNLO;
	vector<results_c> resummed_higgs_LP_NNNLL_exp_NLO, resummed_higgs_LP_NNNLL_NLP_LL_exp_NLO,resummed_higgs_LP_NNLL_exp_NLO, resummed_higgs_LP_NNLL_NLP_LL_exp_NLO, resummed_higgs_LP_NLL_exp_NLO, resummed_higgs_LP_NLL_NLP_LL_exp_NLO, resummed_higgs_LP_LL_NLP_LL_exp_NLO, resummed_higgs_LP_LL_exp_NLO;
	vector<results_c> SCET_higgs_LP_NNLL, SCET_higgs_LP_NNLL_NLP_LL, SCET_higgs_LP_NLL, SCET_higgs_LP_NLL_NLP_LL, SCET_higgs_LP_LL_NLP_LL, SCET_higgs_LP_LL;

	////////////////////////////////////////////////////////////////////

	///////////////////
	/// output
	///////////////////
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
	string q_str = "results/DYh_24082020/output_Q" + Qstring +"_as"+asstring+"_"+setname+"_muF"+mufstring;
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
	else if(realPDF){q_str = q_str+"_real_pdfs"; RES=false; fitPDF=false;} //if real pdfs, the resummed cannot be used!
	if(RES) q_str = q_str+"_resummed";
	q_str = q_str + ".txt";
	cout << q_str.c_str() << endl;
	////////////////////////////////////////////////////////////////////

	/////////////////////////
	/// computation
	/////////////////////////
	int powerexpansion = 20;

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
				higgs_NNLO_gg_full.push_back({res_higgs_NNLO_gg_hard[0].res+res_higgs_NNLO_gg_LP_cor[0].res+res_higgs_NNLO_gg_delta[0].res, res_higgs_NNLO_gg_hard[0].err+res_higgs_NNLO_gg_LP_cor[0].err+res_higgs_NNLO_gg_delta[0].err, res_higgs_NNLO_gg_hard[0].prob+res_higgs_NNLO_gg_LP_cor[0].prob+res_higgs_NNLO_gg_delta[0].prob});
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
				cout << resummed_higgs_LP_NNLL_NLP_LL[0].res << " " << resummed_higgs_LP_NNLL_NLP_LL[0].err << endl;

				cout << "LP NLL + NLP LL" << endl;
				ISNNNLL = 0;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 1;
				resummed_higgs_LP_NLL_NLP_LL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
				resummed_higgs_LP_NLL_NLP_LL_exp_NLO = call_cuhre_higgs("resum","gg","pQCDNLO",fitPDF);
				resummed_higgs_LP_NLL_NLP_LL_exp_NNLO = call_cuhre_higgs("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_higgs_LP_NLL_NLP_LL[0].res << " " << resummed_higgs_LP_NLL_NLP_LL[0].err << endl;

				cout << "LP LL + NLP LL" << endl;
				ISNNNLL = 0;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 0;
				resummed_higgs_LP_LL_NLP_LL = call_cuhre_higgs("resum","gg","pQCD",fitPDF);
				resummed_higgs_LP_LL_NLP_LL_exp_NLO = call_cuhre_higgs("resum","gg","pQCDNLO",fitPDF);
				resummed_higgs_LP_LL_NLP_LL_exp_NNLO = call_cuhre_higgs("resum","gg","pQCDNNLO",fitPDF);
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
				resummed_DY_LP_NNLL_NLP_LL = call_cuhre_dy("resum","gg","pQCD",fitPDF);
				resummed_DY_LP_NNLL_NLP_LL_exp_NLO =call_cuhre_dy("resum","gg","pQCDNLO",fitPDF);
				resummed_DY_LP_NNLL_NLP_LL_exp_NNLO = call_cuhre_dy("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_DY_LP_NNLL_NLP_LL[0].res << " " << resummed_DY_LP_NNLL_NLP_LL[0].err << endl;

				cout << "LP NLL + NLP LL" << endl;
				ISNNNLL = 0;
				ISNNLL = 0;
				ISNLP = 1;
				ISLL = 1;
				ISNLL = 1;
				resummed_DY_LP_NLL_NLP_LL = call_cuhre_dy("resum","gg","pQCD",fitPDF);
				resummed_DY_LP_NLL_NLP_LL_exp_NLO = call_cuhre_dy("resum","gg","pQCDNLO",fitPDF);
				resummed_DY_LP_NLL_NLP_LL_exp_NNLO = call_cuhre_dy("resum","gg","pQCDNNLO",fitPDF);
				cout << resummed_DY_LP_NLL_NLP_LL[0].res << " " << resummed_DY_LP_NLL_NLP_LL[0].err << endl;

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
	cout << "Making the printouts" << endl;
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
			if(doexpansions){
				output << "Fractional qg xsec (at NNLO):			" <<  res_higgs_NNLO_qg_full[0].res << " pb +/-"  <<  res_higgs_NNLO_qg_full[0].err  <<  endl;
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
			output << "--------------------------------" << endl;

			output << "Resummed (LP NLL + NLP LL): 		" << resummed_higgs_LP_NLL_NLP_LL[0].res << " pb +/- " << resummed_higgs_LP_NLL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_NLL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_NLL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_NLL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_NLL_NLP_LL_exp_NNLO[0].err <<  endl;
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

			output << "Resummed (LP LL): 		" << resummed_higgs_LP_LL[0].res << " pb +/- " << resummed_higgs_LP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_higgs_LP_LL_exp_NLO[0].res << " pb +/- " << resummed_higgs_LP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_higgs_LP_LL_exp_NNLO[0].res << " pb +/- " << resummed_higgs_LP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

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
			output << "Z-dependent, including LP (NLO):		" <<  res_DY_NLO_qqbar_hard[0].res << " pb +/- " << res_DY_NLO_qqbar_hard[0].err <<  endl;
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
			output << "Z-dependent, including LP (NNLO):		" <<  res_DY_NNLO_qqbar_hard[0].res << " pb +/- " << res_DY_NNLO_qqbar_hard[0].err <<  endl;
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

			output << "Resummed (LP NLL + NLP LL): 		" << resummed_DY_LP_NLL_NLP_LL[0].res << " pb +/- " << resummed_DY_LP_NLL_NLP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_NLL_NLP_LL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_NLL_NLP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_NLL_NLP_LL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_NLL_NLP_LL_exp_NNLO[0].err <<  endl;
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

			output << "Resummed (LP LL): 		" << resummed_DY_LP_LL[0].res << " pb +/- " << resummed_DY_LP_LL[0].err <<  endl;
			output << "Expanded to NLO:				" <<  resummed_DY_LP_LL_exp_NLO[0].res << " pb +/- " << resummed_DY_LP_LL_exp_NLO[0].err <<  endl;
			output << "Expanded to NNLO:				" <<  resummed_DY_LP_LL_exp_NNLO[0].res << " pb +/- " << resummed_DY_LP_LL_exp_NNLO[0].err <<  endl;
			output << "--------------------------------" << endl;

			output << "======================================================" << endl;
		}
		}
}
	output.close();
	return 0;
}
