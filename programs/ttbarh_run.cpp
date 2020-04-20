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
	configure(argc,argv, to_string2("ttH.cfg"), false);
	results higgs1;
	CMP = 1.5;
	phiMP = 3./4.*M_PI;
	SCET = false;
	double z =0.1;
	lumni_params params = {z, Q, 2*Q/S, 0, 0, 0,0,0};
	double probescales_pt[29] = {10.,20.,25.,30.,40.,50.,60.,70.,75.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,200.,220.,260.,280.,300.,350.,400.,500.,600.};
	cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
	for(int i=0; i<29; i++){
	Q = probescales_pt[i];
	muF = Q;
	muR = Q;
	update_defaults();
	higgs1 = call_vegas(init_vegas_ttH("LOpT"),params, true, true);
	cout << "LO, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
	higgs1 = call_vegas(init_vegas_ttH("LOpTres"),params, true, true);
	cout << "Resummed, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
	}

	phiMP = 1./2.*M_PI;
	cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
	for(int i=0; i<29; i++){
	Q = probescales_pt[i];
	muF = Q;
	muR = Q;
	update_defaults();
	higgs1 = call_vegas(init_vegas_ttH("LOpTres"),params, true, true);
	cout << "Resummed, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
	}
	
	CMP = 2.1;
	phiMP = 1./2.*M_PI;
	cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
	for(int i=0; i<29; i++){
	Q = probescales_pt[i];
	muF = Q;
	muR = Q;
	update_defaults();
	higgs1 = call_vegas(init_vegas_ttH("LOpTres"),params, true, true);
	cout << "Resummed, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
	}

	CMP = 2.1;
	phiMP = 3./4.*M_PI;
	cout << "CMP = " << CMP << " phiMP = " << phiMP << endl;
	for(int i=0; i<29; i++){
	Q = probescales_pt[i];
	muF = Q;
	muR = Q;
	update_defaults();
	higgs1 = call_vegas(init_vegas_ttH("LOpTres"),params, true, true);
	cout << "Resummed, Q=" << Q  << ", " << higgs1.res << " " << higgs1.err << endl;
	}

}
