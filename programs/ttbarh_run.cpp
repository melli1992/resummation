#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include "monte_carlo.h"
#include "mellin_pdf.h"
#include "tth_vegas.h"
#include "resum_tth.h"
#include "tth_softanom.h"
#include "parameters.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <string.h>
#include <sstream>
#include "inout.h"
using namespace std;


string to_string_round(double Q){
	ostringstream q_to_str;
	q_to_str << Q;
	return q_to_str.str();
}


string to_string2(string Q){
	ostringstream q_to_str;
	q_to_str << Q;
	return q_to_str.str();
}

double closest(vector<double> const& vec, double value) {
    auto it = lower_bound(vec.begin(), vec.end(), value);
		auto itm = prev(it);

		double itlow = abs(*itm - value);
		double ithigh = abs(*it - value);
		if (itlow < ithigh){ return *itm;}
    return *it;
}

vector<double> setscales(double M, vector<double> muF_values_pdf){
	vector<double> muvar{0,0,0};
	if(setdym){
			muF = M;
			//muF = sqrt(mT*mTtt);
			muF = closest(muF_values_pdf,muF);
			muR = muF;
			muvar[0] = closest(muF_values_pdf,muF/2.);
			muvar[1] = muF;
			muvar[2] = closest(muF_values_pdf,muF*2);
			}
	else{
		if(highscale){
		muF = 470.;
		muR = 470.;
		muvar[0] = closest(muF_values_pdf,muF/2.);
		muvar[1] = muF;
		muvar[2] = closest(muF_values_pdf,muF*2);
		}
		else{
			muF = 235.;
			muR = 235.;
			muvar[0] = closest(muF_values_pdf,muF/2.);
			muvar[1] = muF;
			muvar[2] = closest(muF_values_pdf,muF*2);
		}
	}
	return muvar;
}



int main(int argc, char* argv[]){
	//////////////////////////////////////////////
	/// predefinition of everything, setting it up
	//////////////////////////////////////////////
	configure(argc,argv, to_string2("ttH.cfg"), false);
	//cout << (int)argv[3] << endl;
	Q = (2.*mt+mH);
	muF = 235.;
	muR = 235.;
	update_defaults();
	vector<double> muF_values_pdf;
	for (std::unordered_map<double, std::vector<std::vector<double>>>::iterator it=fitcoeff.begin(); it!=fitcoeff.end(); ++it)
	muF_values_pdf.push_back((double) it->first);
	sort(muF_values_pdf.begin(),muF_values_pdf.end());
	bool totalxsec=false, pTdist=false, Qdist=false,fixedN=false, sttpTdist=true;
	results higgs1;
	double z =0.1;
	lumni_params params = {z, Q, 2*Q/S, 0, 0, 0,0,0};

	//////////////////////////////////////////////
	/// creating the output file
	//////////////////////////////////////////////
	ofstream output;
	string homedir = "ttH_22072020";
	string q_str = "results/"+homedir+"/output_ttH_res_CMP_"+to_string_round(CMP)+"_phiMP_"+to_string_round(phiMP);
	if(LO){q_str = "results/"+homedir+"/output_ttH_LO_CMP_"+to_string_round(CMP)+"_phiMP_"+to_string_round(phiMP);}

	//distributions
	if(totalxsec){q_str = q_str + "_totalxsec";}
	else if(pTdist){q_str = q_str + "_pTdist";}
	else if(Qdist){q_str = q_str + "_Qdist";}
	else if(sttpTdist){ostringstream number;
											number << argv[3];
											q_str = q_str + "_sttpTdist"+"_"+number.str();}

	//scale choices
	if(setdym){
		if(totalxsec){
				 cout << "Scale choice not available, exiting" << endl;
				 exit(0);
			 }
		if(highscale){q_str = q_str + "_dymscale_Q";}
		else{
			 if(Qdist){
				 cout << "Scale choice not available, exiting" << endl;
				 exit(0);
			 }
			q_str = q_str + "_dymscale_HT";
		}
	}
	if(!setdym){
		if(sttpTdist){
				 cout << "Scale choice not available, exiting" << endl;
				 exit(0);
			 }
		if(highscale){q_str = q_str + "_fixscale_high";}
		else {q_str = q_str + "_fixscale";}}

	//match to NLO
	if(!LO){
		expansion=INCSQRTZ;
		if(expansion){q_str = q_str + "_tomatch";}

	//euler
		if(INCEULER == 1){q_str = q_str + "_withEuler.txt";}
		if(INCEULER == 0){q_str = q_str + "_woEuler.txt";}
	}
	else{
		q_str = q_str + ".txt";
	}

	output.open(q_str.c_str()); //.c_str() needed to do a constant string conversion
	cout << q_str << endl;


	//////////////////////////////////////////////
	/// code for total xsec
	//////////////////////////////////////////////
	if(totalxsec){
		cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
		output << "CMP = " << CMP << " phiMP = " << phiMP << endl;
		Q2 = pow(mH+2.*mt,2); Q = sqrt(Q2);
		vector<double> muvar = setscales(Q,muF_values_pdf);
		for(int j = 0; j < 3; j++){
			muF = muvar[j];
			muR = muvar[j];
			tau = Q2/S2;
			update_defaults();
			if(j == 0){output << "LOW VALUE" << endl;}
			if(j == 1){output << "CENTRAL VALUE" << endl;}
			if(j == 2){output << "HIGH VALUE" << endl;}
			cout << " muF[j] = " << muF << ", j=" << j << " muR[j] = " << muR << ", j=" << j << endl;
			output << " muF[j] = " << muF << ", j=" << j << " muR[j] = " << muR << ", j=" << j << endl;
			LO = true;
			if(LO){
				higgs1 = call_vegas(init_vegas_ttH("tot_N1_2nd"),params, true, true);
				cout << "LO_N1, " << higgs1.res << " " << higgs1.err << endl;
				output << "LO_N1, " << higgs1.res << " " << higgs1.err << endl;
				higgs1 = call_vegas(init_vegas_ttH("tot_N1"),params, true, true);
				cout << "LO_N1, " << higgs1.res << " " << higgs1.err << endl;
				output << "LO_N1, " << higgs1.res << " " << higgs1.err << endl;
				higgs1 = call_vegas(init_vegas_ttH("tot_N2"),params, true, true);
				cout << "LO_N2, " << higgs1.res << " " << higgs1.err << endl;
				output << "LO_N2, " << higgs1.res << " " << higgs1.err << endl;
				higgs1 = call_vegas(init_vegas_ttH("tot_N3"),params, true, true);
				cout << "LO_N3, " << higgs1.res << " " << higgs1.err << endl;
				output << "LO_N3, " << higgs1.res << " " << higgs1.err << endl;
				higgs1 = call_vegas(init_vegas_ttH("tot_N5"),params, true, true);
				cout << "LO_N5, " << higgs1.res << " " << higgs1.err << endl;
				output << "LO_N5, " << higgs1.res << " " << higgs1.err << endl;
				fitPDF = false;
				higgs1 = call_vegas(init_vegas_ttH("LO"),params, true, true);
				cout << "LO_real_PDF, " << higgs1.res << " " << higgs1.err << endl;
				output << "LO_real_PDF, " << higgs1.res << " " << higgs1.err << endl;
				fitPDF = true;
			}
			else{
				ISLL = 1; ISNLL = 0;
				higgs1 = call_vegas(init_vegas_ttH("tot_N1_2nd"),params, true, true);
				cout << "LO_N1, " << higgs1.res << " " << higgs1.err << endl;
				output << "LO_N1, " << higgs1.res << " " << higgs1.err << endl;
			  higgs1 = call_vegas(init_vegas_ttH("tot_N1"),params, true, true);
				cout << "N1(abs), " << higgs1.res << " " << higgs1.err << endl;
				output << "N1(abs), " << higgs1.res << " " << higgs1.err << endl;
				higgs1 = call_vegas(init_vegas_ttH("tot_N2"),params, true, true);
				cout << "N2(stt), " << higgs1.res << " " << higgs1.err << endl;
				output << "N2(stt), " << higgs1.res << " " << higgs1.err << endl;
				higgs1 = call_vegas(init_vegas_ttH("tot_N3"),params, true, true);
				cout << "N3(Q), " << higgs1.res << " " << higgs1.err << endl;
				output << "N3(Q), " << higgs1.res << " " << higgs1.err << endl;
				higgs1 = call_vegas(init_vegas_ttH("tot_N5"),params, true, true);
				cout << "N5(sttpT), " << higgs1.res << " " << higgs1.err << endl;
				output << "N5(sttpT), " << higgs1.res << " " << higgs1.err << endl;
				ISLL = 1; ISNLL = 1;
				diagsoft = true;
				higgs1 = call_vegas(init_vegas_ttH("tot_N1"),params, true, true);
				cout << "N1(abs)_diag, " << higgs1.res << " " << higgs1.err << endl;
				output << "N1(abs)_diag, " << higgs1.res << " " << higgs1.err << endl;
				diagsoft = false;
				higgs1 = call_vegas(init_vegas_ttH("tot_N1"),params, true, true);
				cout << "N1(abs), " << higgs1.res << " " << higgs1.err << endl;
				output << "N1(abs), " << higgs1.res << " " << higgs1.err << endl;
				higgs1 = call_vegas(init_vegas_ttH("tot_N2"),params, true, true);
				cout << "N2(stt), " << higgs1.res << " " << higgs1.err << endl;
				output << "N2(stt), " << higgs1.res << " " << higgs1.err << endl;
				higgs1 = call_vegas(init_vegas_ttH("tot_N3"),params, true, true);
				cout << "N3(Q), " << higgs1.res << " " << higgs1.err << endl;
				output << "N3(Q), " << higgs1.res << " " << higgs1.err << endl;
				higgs1 = call_vegas(init_vegas_ttH("tot_N5"),params, true, true);
				cout << "N5(sttpT), " << higgs1.res << " " << higgs1.err << endl;
				output << "N5(sttpT), " << higgs1.res << " " << higgs1.err << endl;

			}
		}
	}
	//////////////////////////////////////////////
	/// code for pT dist
	//////////////////////////////////////////////

	if(pTdist){
		double probescales_pt[31] = {10,20,25,30,40,50,60,70,75,80,90,100,110,120,130,140,150,160,170,180,200,210,220,240,260,280,300,350,400,450,500};
		cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
		output << "CMP = " << CMP << " phiMP = " << phiMP << endl;
		vector<double> muvar = setscales(Q,muF_values_pdf);
		for(int i=0; i<31; i++){
			pT = probescales_pt[i];
			pT2 = pow(probescales_pt[i],2);
			double mT = sqrt(pT2+mH2);
			double mTtt = sqrt(pT2+4.*mt2);
			Q2 = pow(mT+mTtt,2); Q = sqrt(Q2);
			tau = Q2/S2;
			if(highscale){muvar = setscales(Q,muF_values_pdf);}
			else{muvar = setscales(sqrt(mT*mTtt),muF_values_pdf);}
			cout << "Values for scales: " << muvar[0] << " " << muvar[1] << " " << muvar[2] << endl;
			for(int j = 0; j < 3; j++){
				muF = muvar[j];
				muR = muvar[j];
				Q2 = pow(mT+mTtt,2); Q = sqrt(Q2);
				tau = Q2/S2;
				update_defaults();
				if(j == 0){output << "LOW VALUE" << endl;}
				if(j == 1){output << "CENTRAL VALUE" << endl;}
				if(j == 2){output << "HIGH VALUE" << endl;}
				cout << "pT = " << pT << " Q = " << Q << " muF[j] = " << muF << ", j=" << j << " muR[j] = " << muR << ", j=" << j << endl;
				output << "pT = " << pT << " Q = " << Q << " muF[j] = " << muF << ", j=" << j << " muR[j] = " << muR << ", j=" << j << endl;
				if(LO){
					higgs1 = call_vegas(init_vegas_ttH("pT_N4"),params, true, true);
					cout << "LO_N4, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
					output << "LO_N4, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
					fitPDF = false;
					higgs1 = call_vegas(init_vegas_ttH("LOpT"),params, true, true);
					cout << "LO_real_PDF, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
					output << "LO_real_PDF, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
					fitPDF = true;
				}
				else{
					if(fixedN){
						/*ISLL = 1; ISNLL = 0;
					sgg = newton_raphson_gg(1.,0.,tau,1.);
					sqqbar = newton_raphson_qqbar(1.,0.,tau,1.);
					cout << "mT = " << mT << ", mTtt =" << mTtt << ", Q2 = " << Q2 << ", tau =" << tau << ", N_gluon = " << sgg*exp(INCEULER*M_gammaE) << ", N_quark = " << sqqbar*exp(INCEULER*M_gammaE) << endl;
					if(INCEULER) output << "mT = " << mT << ", mTtt =" << mTtt << ", tau =" << tau << ", Nbar_gluon = " << sgg*exp(INCEULER*M_gammaE) << ", Nbar_quark = " << sqqbar*exp(INCEULER*M_gammaE) << endl;
					else output << "mT = " << mT << ", mTtt =" << mTtt << ", tau =" << tau << ", N_gluon = " << sgg*exp(INCEULER*M_gammaE) << ", N_quark = " << sqqbar*exp(INCEULER*M_gammaE) << endl;
					higgs1 = call_vegas(init_vegas_ttH("pTresNfix"),params, true, true);
					cout << "LLResummedNfix, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
					output << "LLResummedNfix, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
					*/
					}
					else{
						ISLL = 1; ISNLL = 0;
						higgs1 = call_vegas(init_vegas_ttH("pT_N4"),params, true, true);
						cout << "N4(4mt2), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "N4(4mt2), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						higgs1 = call_vegas(init_vegas_ttH("pT_N5"),params, true, true);
						cout << "N5(stt), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "N5(stt), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						ISLL = 1; ISNLL = 1;
						diagsoft = true;
						higgs1 = call_vegas(init_vegas_ttH("pT_N4"),params, true, true);
						cout << "N4(4mt2)_diagsoft, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "N4(4mt2)_diagsoft, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						diagsoft = false;
						higgs1 = call_vegas(init_vegas_ttH("pT_N4"),params, true, true);
						cout << "N4(4mt2), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "N4(4mt2), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						higgs1 = call_vegas(init_vegas_ttH("pT_N5"),params, true, true);
						cout << "N5(stt), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "N5(stt), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
					}
				}
			}
		}
	}

	//////////////////////////////////////////////
	/// code for Q dist
	//////////////////////////////////////////////

	if(Qdist){
		double probescales_Q[31] = {490,500,510,520,530,540,550,560,570,580,590,600,610,620,630,640,650,660,670,680,690,700,710,725,750,775,800,850,900,950,1000};
		cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
		output << "CMP = " << CMP << " phiMP = " << phiMP << endl;
		vector<double> muvar = setscales(Q,muF_values_pdf);
		for(int i=0; i<31; i++){
			Q = probescales_Q[i];
			Q2 = pow(probescales_Q[i],2);
			tau = Q2/S2;
			muvar = setscales(Q,muF_values_pdf);
			cout << "Values for scales: " << muvar[0] << " " << muvar[1] << " " << muvar[2] << endl;
			for(int j = 0; j < 3; j++){
				muF = muvar[j];
				muR = muvar[j];
				update_defaults();
				if(j == 0){output << "LOW VALUE" << endl;}
				if(j == 1){output << "CENTRAL VALUE" << endl;}
				if(j == 2){output << "HIGH VALUE" << endl;}
				cout << " Q = " << Q << " muF[j] = " << muF << ", j=" << j << " muR[j] = " << muR << ", j=" << j << endl;
				output << " Q = " << Q << " muF[j] = " << muF << ", j=" << j << " muR[j] = " << muR << ", j=" << j << endl;
				if(LO){
					higgs1 = call_vegas(init_vegas_ttH("inv_mass_N1"),params, true, true);
					cout << "N1(abs), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
					output << "N1(abs), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
					higgs1 = call_vegas(init_vegas_ttH("inv_mass_N2"),params, true, true);
					cout << "N2(stt), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
					output << "N2(stt), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
					higgs1 = call_vegas(init_vegas_ttH("inv_mass_N3"),params, true, true);
					cout << "N3(Q), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
					output << "N3(Q), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
				}
				else{
					if(fixedN){
						/*ISLL = 1; ISNLL = 0;
					sgg = newton_raphson_gg(1.,0.,tau,1.);
					sqqbar = newton_raphson_qqbar(1.,0.,tau,1.);
					cout << "mT = " << mT << ", mTtt =" << mTtt << ", Q2 = " << Q2 << ", tau =" << tau << ", N_gluon = " << sgg*exp(INCEULER*M_gammaE) << ", N_quark = " << sqqbar*exp(INCEULER*M_gammaE) << endl;
					if(INCEULER) output << "mT = " << mT << ", mTtt =" << mTtt << ", tau =" << tau << ", Nbar_gluon = " << sgg*exp(INCEULER*M_gammaE) << ", Nbar_quark = " << sqqbar*exp(INCEULER*M_gammaE) << endl;
					else output << "mT = " << mT << ", mTtt =" << mTtt << ", tau =" << tau << ", N_gluon = " << sgg*exp(INCEULER*M_gammaE) << ", N_quark = " << sqqbar*exp(INCEULER*M_gammaE) << endl;
					higgs1 = call_vegas(init_vegas_ttH("pTresNfix"),params, true, true);
					cout << "LLResummedNfix, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
					output << "LLResummedNfix, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
					*/
					}
					else{
						ISLL = 1; ISNLL = 0;
						higgs1 = call_vegas(init_vegas_ttH("inv_mass_N1"),params, true, true);
						cout << "N1(abs), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "N1(abs), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						higgs1 = call_vegas(init_vegas_ttH("inv_mass_N2"),params, true, true);
						cout << "N2(stt), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "N2(stt), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						higgs1 = call_vegas(init_vegas_ttH("inv_mass_N3"),params, true, true);
						cout << "N3(Q), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "N3(Q), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						ISLL = 1; ISNLL = 1;
						diagsoft = true;
						higgs1 = call_vegas(init_vegas_ttH("inv_mass_N1"),params, true, true);
						cout << "N1(abs)_diagsoft, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "N1(abs)_diagsoft, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						diagsoft = false;
						higgs1 = call_vegas(init_vegas_ttH("inv_mass_N1"),params, true, true);
						cout << "N1(abs), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "N1(abs), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						higgs1 = call_vegas(init_vegas_ttH("inv_mass_N2"),params, true, true);
						cout << "N2(stt), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "N2(stt), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						higgs1 = call_vegas(init_vegas_ttH("inv_mass_N3"),params, true, true);
						cout << "N3(Q), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "N3(Q), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
					}
				}
			}
		}
	}

	//////////////////////////////////////////////
	/// code for pT stt dist
	//////////////////////////////////////////////

	if(sttpTdist){
		ostringstream number;
		number << argv[3];
		int k = stoi(number.str());
		double probescales_pt[31] = {10,20,25,30,40,50,60,70,75,80,90,100,110,120,130,140,150,160,170,180,200,210,220,240,260,280,300,350,400,450,500};
		//double probescales_stt[7] = {4.5*mt2,4.75*mt2,5.*mt2,5.5*mt2,6.*mt2,8.*mt2,12*mt2};
		double probescales_stt[7] = {4.2*mt2,5.*mt2,6.*mt2,8.*mt2,12.*mt2,16.*mt2,20.*mt2};
		cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
		output << "CMP = " << CMP << " phiMP = " << phiMP << endl;
		vector<double> muvar = setscales(Q,muF_values_pdf);
		s34in = probescales_stt[k];
		for(int i=0; i<31; i++){
					pT = probescales_pt[i];
					pT2 = pow(probescales_pt[i],2);
					double mT = sqrt(pT2+mH2);
					double mTtt = sqrt(pT2+4.*mt2);
					Q2 = pow(mT+mTtt,2); Q = sqrt(Q2);
					tau = Q2/S2;
					vector<double> muvar1 = setscales(Q,muF_values_pdf);
					if(!highscale){ muvar1 = setscales(sqrt(mT*mTtt),muF_values_pdf);}
					mTtt = sqrt(pT2+s34in);
					Q2 = pow(mT+mTtt,2); Q = sqrt(Q2);
					tau = Q2/S2;
					vector<double> muvar2 = setscales(Q,muF_values_pdf);
					if(!highscale){ muvar2 = setscales(sqrt(mT*mTtt),muF_values_pdf);}
					muvar[0] = muvar1[1];
					muvar[1] = muvar2[1];
					cout << "Values for scales: " << muvar[0] << " " << muvar[1] << endl;
					for(int j = 0; j < 2; j++){
						muF = muvar[j];
						muR = muvar[j];
						Q2 = pow(mT+mTtt,2); Q = sqrt(Q2);
						tau = Q2/S2;
						update_defaults();
						if(j == 0){output << "mTtt_4mt2" << endl;cout << "mTtt_4mt2" << endl;}
						if(j == 1){output << "mTtt_stt" << endl;cout << "mTtt_stt" << endl;}
						//fitPDF = false;
						cout << "stt = " << s34in << " pT = " << pT << " Q = " << Q << " muF[j] = " << muF << ", j=" << j << " muR[j] = " << muR << ", j=" << j << endl;
						output <<  "stt = " << s34in << " pT = " << pT << " Q = " << Q << " muF[j] = " << muF << ", j=" << j << " muR[j] = " << muR << ", j=" << j << endl;
						//higgs1 = call_vegas(init_vegas_ttH("pTstt"),params, true, true);
						//cout << "LO, stt= " << s34in << " pT= " << pT  << " " << higgs1.res << " " << higgs1.err << endl;
						//output << "LO, stt= " << s34in << " pT= " << pT  << " " << higgs1.res << " " << higgs1.err << endl;
						if(LO){
					  higgs1 = call_vegas(init_vegas_ttH("sttpT_N4"),params, true, true);
						cout << "N4, stt= " << s34in << " pT= " << pT  << " " << higgs1.res << " " << higgs1.err << endl;
						output << "N4, stt= " << s34in << " pT= " << pT  << " " << higgs1.res << " " << higgs1.err << endl;
						higgs1 = call_vegas(init_vegas_ttH("sttpT_N5"),params, true, true);
						cout << "N5, stt= " << s34in << " pT= " << pT  << " " << higgs1.res << " " << higgs1.err << endl;
						output << "N5, stt= " << s34in << " pT= " << pT  << " " << higgs1.res << " " << higgs1.err << endl;
						}
						diagsoft = false;
						higgs1 = call_vegas(init_vegas_ttH("sttpT_N4"),params, true, true);
						cout << "N4res, stt= " << s34in << " pT= " << pT  << " " << higgs1.res << " " << higgs1.err << endl;
						output << "N4res, stt= " << s34in << " pT= " << pT  << " " << higgs1.res << " " << higgs1.err << endl;
						higgs1 = call_vegas(init_vegas_ttH("sttpT_N5"),params, true, true);
						cout << "N5res, stt= " << s34in << " pT= " << pT  << " " << higgs1.res << " " << higgs1.err << endl;
						output << "N5res, stt= " << s34in << " pT= " << pT  << " " << higgs1.res << " " << higgs1.err << endl;
					}
				}
		}
}
