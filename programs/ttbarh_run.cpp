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
			if(fitPDF) muF = closest(muF_values_pdf,muF);
			muR = muF;
			}
	else{
		if(highscale){
		muF = 470.;
		muR = 470.;
		}
		else{
			muF = 235.;
			muR = 235.;
		}
	}
	if(fitPDF){
		muvar[0] = closest(muF_values_pdf,muF/2.);
		muvar[1] = muF;
		muvar[2] = closest(muF_values_pdf,muF*2);
	}
	else{
		muvar[0] = muF/2.;
		muvar[1] = muF;
		muvar[2] = muF*2;
	}
	return muvar;
}



int main(int argc, char* argv[]){
	//////////////////////////////////////////////
	/// predefinition of everything, setting it up
	//////////////////////////////////////////////
	configure(argc,argv, to_string2("ttH.cfg"), false);
	results higgs1;
	double z =0.1;
	lumni_params params = {z, Q, 2*Q/S, 0, 0, 0,0,0};

//cout << (int)argv[3] << endl;
  //double Qvals[26] = {510.0, 530.0, 550.0, 570.0, 590.0, 610.0, 630.0, 650.0, 670.0, 690.0, 710.0, 730.0, 750.0, 770.0, 790.0, 810.0, 830.0, 850.0, 870.0, 890.0, 910.0, 930.0, 950.0, 970.0, 990.0, 1010.0};
	//double pTvals[50] = {5.0,15.0,25.0,35.0,45.0,55.0,65.0,75.0,85.0,95.0,105.0,115.0,125.0,135.0,145.0,155.0,165.0,175.0,185.0,195.0,205.0,215.0,225.0,235.0,245.0,255.0,265.0,275.0,285.0,295.0,305.0,315.0,325.0,335.0,345.0,355.0,365.0,375.0,385.0,395.0,405.0,415.0,425.0,435.0,445.0,455.0,465.0,475.0,485.0,495.0};
//{10.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0, 190.0, 210.0, 230.0, 250.0, 270.0, 290.0, 310.0, 330.0, 350.0, 370.0, 390.0, 410.0, 430.0, 450.0, 470.0, 490.0, 510.0};

  //double Qvals = {2281.65,1137.22,1140.83,568.61,1010,505,2020,1020,1104.96,552.481,275.248,2209.93,1100.99,1980,980,1069.37,532.494,534.685,266.247,2138.74,1064.99,970,485,1940,225,999.094,496.86,499.547,248.43,1998.19,993.719,930,465,215,1860,860,1928.94,958.495,482.236,964.472,910,455,205,1820,820,1860.48,923.574,930.239,461.787,465.12,445,195,1780,780,890,896.435,444.491,448.217,222.246,1792.87,888.983,870,435,1740,740,175,830.296,410.452,415.148,205.226,1660.59,820.904,830,415,1660,660,798.069,393.741,399.034,196.87,1596.14,787.481,620,1620,405,810,766.487,377.261,383.243,188.63,1532.97,754.522,1580,395,790,1471.25,722.07,735.624,361.035,367.812,1540,172.545,1411.13,676.411,329.456,338.206,164.728,1352.82,658.911,365,1460,730,1420,710,355,621.29,299.292,310.645,149.646,1242.58,598.585,1340,1380,335,670,1020,1060,1140,305,295,1180,610,1220,630,1260}
	//Q = (2.*mt+mH);
	muF = 235.;
	muR = 235.;
	update_defaults();
	vector<double> muF_values_pdf;
	for (std::unordered_map<double, std::vector<std::vector<double>>>::iterator it=fitcoeff.begin(); it!=fitcoeff.end(); ++it)
	muF_values_pdf.push_back((double) it->first);
	sort(muF_values_pdf.begin(),muF_values_pdf.end());
/*
	expansion=false;
	deform = true;
	cout << "START COUNT" << endl;
	higgs1 = call_vegas(init_vegas_ttH("tot_N2"),params, true, true);
	cout << "END COUNT" << endl;
	cout << "N2, deform, " << higgs1.res << " " << higgs1.err << endl;
	deform = false;
	cout << "START COUNT" << endl;
	higgs1 = call_vegas(init_vegas_ttH("tot_N2"),params, true, true);
	cout << "END COUNT" << endl;
	cout << "N2, fit PDF, " << higgs1.res << " " << higgs1.err << endl;
	fitPDF = false; boundary = false;
	cout << "START COUNT" << endl;
	higgs1 = call_vegas(init_vegas_ttH("tot_N2"),params, true, true);
	cout << "END COUNT" << endl;
	cout << "N2, true PDF, " << higgs1.res << " " << higgs1.err << endl;
	fitPDF = false; boundary = true;
	cout << "START COUNT" << endl;
	higgs1 = call_vegas(init_vegas_ttH("tot_N2"),params, true, true);
	cout << "END COUNT" << endl;
	cout << "N2, true PDF boundary term, " << higgs1.res << " " << higgs1.err << endl;

	exit(0);
*/
	//muF = 470.;
	//muR = 470.;
	//update_defaults();
/*
		for(int i = 0; i < 26; i++){
			//cout << "=====" << endl;
			//cout << Qvals[i] << "," << pTvals[i] << "," << Qvals[i]/2. << "," << pTvals[i]/2. << "," << Qvals[i]*2. << "," << pTvals[i]*2. << endl;
			//cout << closest(muF_values_pdf,Qvals[i]) << "," << closest(muF_values_pdf,pTvals[i]) << "," << closest(muF_values_pdf,Qvals[i]/2.) << "," << closest(muF_values_pdf,pTvals[i]/2.) << "," << closest(muF_values_pdf,Qvals[i]*2) << "," << closest(muF_values_pdf,pTvals[i]*2) << endl;
			//cout << "------" << endl;
			double mTt = sqrt(pTvals[i]*pTvals[i]+mH2);
			double mTttt = sqrt(pTvals[i]*pTvals[i]+4.*mt2);
			Q2 = pow(mTt+mTttt,2); Q = sqrt(Q2);
			//cout << Q << "," << sqrt(mTt*mTttt) << "," << Q/2. << "," << sqrt(mTt*mTttt)/2. << "," << Q*2. << "," << sqrt(mTt*mTttt)*2. << "," << endl;
			//cout << closest(muF_values_pdf,Q) << "," << closest(muF_values_pdf,sqrt(mTt*mTttt)) << "," << closest(muF_values_pdf,Q/2.) << "," << closest(muF_values_pdf,sqrt(mTt*mTttt)/2.) << "," << closest(muF_values_pdf,Q*2.) << "," << closest(muF_values_pdf,sqrt(mTt*mTttt)*2.) << "," << endl;
			//if(abs(pTvals[i] - closest(muF_values_pdf,pTvals[i])) > 1) cout << pTvals[i] << endl;
			//if(abs(pTvals[i]/2. - closest(muF_values_pdf,pTvals[i]/2.)) > 1) cout << pTvals[i]/2. << endl;
			//if(abs(pTvals[i]*2. - closest(muF_values_pdf,pTvals[i]*2.)) > 1) cout << pTvals[i]*2. << endl;
			if(abs(Q - closest(muF_values_pdf,Q)) > 1) cout << Q << ",";
			if(abs(Q/2. - closest(muF_values_pdf,Q/2.)) > 1) cout << Q/2. << ",";
			if(abs(Q*2. - closest(muF_values_pdf,Q*2)) > 1) cout << Q*2 << ",";

			if(abs(sqrt(mTt*mTttt) - closest(muF_values_pdf,sqrt(mTt*mTttt))) > 1) cout << sqrt(mTt*mTttt) << ",";//<< " " << closest(muF_values_pdf,sqrt(mTt*mTttt)) << endl;
			if(abs(sqrt(mTt*mTttt)/2. - closest(muF_values_pdf,sqrt(mTt*mTttt)/2.)) > 1) cout << sqrt(mTt*mTttt)/2. << ",";
			if(abs(sqrt(mTt*mTttt)*2. - closest(muF_values_pdf,sqrt(mTt*mTttt)*2.)) > 1) cout << sqrt(mTt*mTttt)*2 << ",";

		}

	exit(0);
	*/
	bool totalxsec=false, pTdist=false, Qdist=false, sttpTdist=false;
	if(observable == (string)"pT"){pTdist = true;}
	else if(observable == (string)"Qinv"){Qdist = true;}
	else if(observable == (string)"stt"){sttpTdist = true;}
	else if(observable == (string)"total"){totalxsec = true;}
	else{
		cout << "Error: observalbe not specified - exiting" <<endl;
		exit(0);
	}

	//////////////////////////////////////////////
	/// creating the output file
	//////////////////////////////////////////////
	ofstream output;
	//string homedir = "ttH_03112020_fixN";
	string homedir = "ttH_26112020_dddd";
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
				//fitPDF = false;
				//higgs1 = call_vegas(init_vegas_ttH("LO"),params, true, true);
				//cout << "LO_real_PDF, " << higgs1.res << " " << higgs1.err << endl;
				//output << "LO_real_PDF, " << higgs1.res << " " << higgs1.err << endl;
				//fitPDF = true;
			}
			else{
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
		//double pTvals[50] = {5.0,15.0,25.0,35.0,45.0,55.0,65.0,75.0,85.0,95.0,105.0,115.0,125.0,135.0,145.0,155.0,165.0,175.0,185.0,195.0,205.0,215.0,225.0,235.0,245.0,255.0,265.0,275.0,285.0,295.0,305.0,315.0,325.0,335.0,345.0,355.0,365.0,375.0,385.0,395.0,405.0,415.0,425.0,435.0,445.0,455.0,465.0,475.0,485.0,495.0};
		//double probescales_pt[50] = {5.0,15.0,25.0,35.0,45.0,55.0,65.0,75.0,85.0,95.0,105.0,115.0,125.0,135.0,145.0,155.0,165.0,175.0,185.0,195.0,205.0,215.0,225.0,235.0,245.0,255.0,265.0,275.0,285.0,295.0,305.0,315.0,325.0,335.0,345.0,355.0,365.0,375.0,385.0,395.0,405.0,415.0,425.0,435.0,445.0,455.0,465.0,475.0,485.0,495.0};
		double probescales_pt[26] = {10.0, 30.0, 50.0, 70.0, 90.0, 110.0, 130.0, 150.0, 170.0, 190.0, 210.0, 230.0, 250.0, 270.0, 290.0, 310.0, 330.0, 350.0, 370.0, 390.0, 410.0, 430.0, 450.0, 470.0, 490.0};

		cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
		output << "CMP = " << CMP << " phiMP = " << phiMP << endl;
		vector<double> muvar = setscales(Q,muF_values_pdf);
		for(int i=0; i<25; i++){
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
				if(Nfixed && (j==0)) continue;
				if(Nfixed && (j==2)) continue;
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
					cout << "LO, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
					output << "LO, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
					//fitPDF = false;
					//higgs1 = call_vegas(init_vegas_ttH("LOpT"),params, true, true);
					//cout << "LO_real_PDF, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
					//output << "LO_real_PDF, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
					//fitPDF = true;
				}
				else{
					if(Nfixed){
						tau = pow(mT+mTtt,2)/S2;
						sgg = newton_raphson_gg(1.,0.,tau,1.);
						sqqbar = newton_raphson_qqbar(1.,0.,tau,1.);
						diagsoft = false;
				  	cout << "N4: mT = " << mT << ", mTtt =" << mTtt << ", M2 = " << pow(mT+mTtt,2) << ", tau =" << tau << ", N_gluon = " << sgg*exp(INCEULER*M_gammaE) << ", N_quark = " << sqqbar*exp(INCEULER*M_gammaE) << endl;
						output << "N4: mT = " << mT << ", mTtt =" << mTtt << ", M2 = " << pow(mT+mTtt,2) << ", tau =" << tau << ", N_gluon = " << sgg*exp(INCEULER*M_gammaE) << ", N_quark = " << sqqbar*exp(INCEULER*M_gammaE) << endl;
						higgs1 = call_vegas(init_vegas_ttH("pT_N4"),params, true, true);
						cout << "Resummed_Nfix N4(4mt2), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "Resummed_Nfix N4(4mt2), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;

					}
					else{
						//ISLL = 1; ISNLL = 0;
						//higgs1 = call_vegas(init_vegas_ttH("pT_N4"),params, true, true);
						//cout << "LL N4(4mt2), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						//output << "LL N4(4mt2), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						//higgs1 = call_vegas(init_vegas_ttH("pT_N5"),params, true, true);
						//cout << "LL N5(stt), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						//output << "LL N5(stt), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						//ISLL = 1; ISNLL = 1;
						//diagsoft = true;
						//higgs1 = call_vegas(init_vegas_ttH("pT_N4"),params, true, true);
						//cout << "NLL N4(4mt2)_diagsoft, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						//output << "NLL N4(4mt2)_diagsoft, pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						diagsoft = false;
						higgs1 = call_vegas(init_vegas_ttH("pT_N4"),params, true, true);
						cout << "NLL N4(4mt2), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "NLL N4(4mt2), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						higgs1 = call_vegas(init_vegas_ttH("pT_N5"),params, true, true);
						cout << "NLL N5(stt), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "NLL N5(stt), pT=" << pT  << ", " << higgs1.res << " " << higgs1.err << endl;

					}
				}
			}
		}
	}

	//////////////////////////////////////////////
	/// code for Q dist
	//////////////////////////////////////////////

	if(Qdist){
		double probescales_Q[26] = {990.0, 530.0,550.0, 570.0, 590.0, 610.0, 630.0, 650.0, 670.0, 690.0, 710.0, 730.0, 750.0, 770.0, 790.0, 810.0, 830.0, 850.0, 870.0, 890.0, 910.0, 930.0, 950.0, 970.0, 990.0, 1010.0};
		cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
		output << "CMP = " << CMP << " phiMP = " << phiMP << endl;
		vector<double> muvar = setscales(Q,muF_values_pdf);
		cout << muvar[0] << endl;
		for(int i=0; i<26; i++){
			Q = probescales_Q[i];
			Q2 = pow(probescales_Q[i],2);
			tau = Q2/S2;
			muvar = setscales(Q,muF_values_pdf);
			cout << "Values for scales: " << muvar[0] << " " << muvar[1] << " " << muvar[2] << endl;
			for(int j = 1; j < 3; j++){
				if(Nfixed && (j==0)) continue;
				if(Nfixed && (j==2)) continue;
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
					cout << "LO, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
					output << "LO, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
				}
				else{
					if(Nfixed){
					ISLL = 1; ISNLL = 1;
					tau = pow((2.*mt+mH),2)/S2;
					sgg = newton_raphson_gg(1.,0.,tau,1.);
					sqqbar = newton_raphson_qqbar(1.,0.,tau,1.);
					cout << "N1: Q2 = " << Q2 << ", tau =" << tau << ", N_gluon = " << sgg*exp(INCEULER*M_gammaE) << ", N_quark = " << sqqbar*exp(INCEULER*M_gammaE) << endl;
					output << "N1: Q2 = " << Q2 << ", tau =" << tau << ", N_gluon = " << sgg*exp(INCEULER*M_gammaE) << ", N_quark = " << sqqbar*exp(INCEULER*M_gammaE) << endl;
					diagsoft = false;
					higgs1 = call_vegas(init_vegas_ttH("inv_mass_N1"),params, true, true);
					cout << "Resummed_Nfix N1(abs), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
					output << "Resummed_Nfix N1(abs), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
					tau = Q2/S2;
					sgg = newton_raphson_gg(1.,0.,tau,1.);
					sqqbar = newton_raphson_qqbar(1.,0.,tau,1.);
					cout << "N3: Q2 = " << Q2 << ", tau =" << tau << ", N_gluon = " << sgg*exp(INCEULER*M_gammaE) << ", N_quark = " << sqqbar*exp(INCEULER*M_gammaE) << endl;
					output << "N3: Q2 = " << Q2 << ", tau =" << tau << ", N_gluon = " << sgg*exp(INCEULER*M_gammaE) << ", N_quark = " << sqqbar*exp(INCEULER*M_gammaE) << endl;
					higgs1 = call_vegas(init_vegas_ttH("inv_mass_N3"),params, true, true);
					cout << "Resummed_Nfix N3(Q), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
					output << "Resummed_Nfix N3(Q), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
					}
					else{
						//ISLL = 1; ISNLL = 0;
						//higgs1 = call_vegas(init_vegas_ttH("inv_mass_N1"),params, true, true);
						//cout << "LL, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						//output << "LL, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						/*cout << "LL N1(abs), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "LL N1(abs), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						higgs1 = call_vegas(init_vegas_ttH("inv_mass_N2"),params, true, true);
						cout << "LL N2(stt), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "LL N2(stt), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						higgs1 = call_vegas(init_vegas_ttH("inv_mass_N3"),params, true, true);
						cout << "LL N3(Q), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "LL N3(Q), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						*/
						//ISLL = 1; ISNLL = 1;
						//diagsoft = true;
						//higgs1 = call_vegas(init_vegas_ttH("inv_mass_N1"),params, true, true);
						//cout << "NLL N1(abs)_diagsoft, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						//output << "NLL N1(abs)_diagsoft, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						//diagsoft = false;
						higgs1 = call_vegas(init_vegas_ttH("inv_mass_N1"),params, true, true);
						cout << "NLL N1(abs), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						output << "NLL N1(abs), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						exit(0);
						//higgs1 = call_vegas(init_vegas_ttH("inv_mass_N2"),params, true, true);
						//cout << "NLL N2(stt), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						//output << "NLL N2(stt), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						//higgs1 = call_vegas(init_vegas_ttH("inv_mass_N3"),params, true, true);
						//cout << "NLL N3(Q), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
						//output << "NLL N3(Q), Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
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
