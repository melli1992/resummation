#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include "monte_carlo.h"
#include "tth_vegas.h"
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

int main(int argc, char* argv[]){
	//////////////////////////////////////////////
	/// predefinition of everything, setting it up
	//////////////////////////////////////////////
	configure(argc,argv, to_string2("ttH.cfg"), false);
	results higgs1;
	ofstream output;
	string q_str = "results/ttH_fixed/output_ttH_Qscan_CMP_"+to_string_round(CMP)+"_phiMP_"+to_string_round(phiMP)+".txt";
	if(LO){q_str = "results/ttH_fixed/output_ttH_Qscan_LO.txt";}
	output.open(q_str.c_str()); //.c_str() needed to do a constant string conversion
	cout << q_str << endl;

	double z =0.1;
	lumni_params params = {z, Q, 2*Q/S, 0, 0, 0,0,0};
	double probescales_pt[29] = {10.,20.,25.,30.,40.,50.,60.,70.,75.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,200.,220.,260.,280.,300.,350.,400.,500.,600.};
	if(!LO) cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
	if(!LO) output << "CMP = " << CMP << " phiMP = " << phiMP << endl;
	/*muF = 240.;
	muR = 240.;
	update_defaults();
	higgs1 = call_vegas(init_vegas_ttH("LOpTfull"),params, true, true);
	cout << "LOpTfull = " << higgs1.res << " " << higgs1.err << endl;
	higgs1 = call_vegas(init_vegas_ttH("LO"),params, true, true);
	cout << "LO = " << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
	*/
	for(int i=0; i<29; i++){
	Q = probescales_pt[i];
	muF = 240.;
	muR = 240.;
	update_defaults();
	if(LO){
		higgs1 = call_vegas(init_vegas_ttH("LOpT"),params, true, true);
		cout << "LO, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
		output << "LO, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
	}
	else{
		higgs1 = call_vegas(init_vegas_ttH("LOpTN"),params, true, true);
		cout << "LONspace, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
		output << "LONspace, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
		higgs1 = call_vegas(init_vegas_ttH("LOpTsttN"),params, true, true);
		cout << "LOsttN, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
		output << "LOsttN, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
		higgs1 = call_vegas(init_vegas_ttH("LOpTsttdefN"),params, true, true);
		cout << "LOsttdef, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
		output << "LOsttdef, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
		higgs1 = call_vegas(init_vegas_ttH("pTres"),params, true, true);
		cout << "Resummed, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
		output << "Resummed, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
		higgs1 = call_vegas(init_vegas_ttH("pTresstt"),params, true, true);
		cout << "Resummedstt, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
		output << "Resummedstt, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
	}
	}

}
