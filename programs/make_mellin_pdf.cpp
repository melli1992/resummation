#include <cmath>
#include <iostream>
#include <fstream>
#include "parameters.h"
#include "inout.h"
#include "LHAPDF/LHAPDF.h"
#include <string>
#include <sstream>
#include <unordered_map>
#include "mellin_pdf.h"

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
	configure(argc, argv, to_string2("input.cfg"), false);
	cout << muF << endl;
	LHAPDF::PDFSet setk(setname);
	int nmem(0.); //number of members
	vector<int> pids; //number of flavors, span from -5 to 5 with 0 = 21 gluon
	nmem = setk.size()-1;
	pdfs = setk.mkPDFs();
	pids = pdfs[0]->flavors();
	xmin_pdfs = pdfs[0]->xMin();
	xmax_pdfs = pdfs[0]->xMax();
	///////////////
	ofstream output, output2;
	ostringstream x_convert;
	x_convert << muF;
	string Qstring  = x_convert.str();
	x_convert << alphas_muF;
	string asstring  = x_convert.str();
	string q2_str = "fit_pdfs/"+setname+"/muF" + Qstring +"_"+setname;
	//string q2_str = "fit_pdfs/"+setname+"/muF" + muF +"_"+setname;
	cout << q2_str << endl;
	q2_str = q2_str + "_pdfoutput.txt";
	
	output2.open(q2_str.c_str()); //.c_str() needed to do a constant string conversion
	output2 << "x xf(x) f(x)" << endl;
	//////////////////////////////////////////////
    cout << q2_str << endl;
    for(int flav = -5; flav < 6; flav++)
	{
		cout << flav << endl;
		cout << endl;
		cout << "ErrorType: " << setk.errorType() << endl;
		cout << "ErrorConfLevel: " << setk.errorConfLevel() << endl;


		double MAX = 1;
		output2 << "muF=" << muF <<" GeV, flavor=" << flav << endl;
		for (double i = 1E-10, increment = 1E-11, counter = 1; i <= MAX; i += increment) {
			vector<double> xfx;
			for (signed int imem = 0; imem <= nmem; imem++) {
			xfx.push_back(pdfs[imem]->xfxQ(flav,i,muF));
			//cout << pdfs[imem]->type() << endl;
			}
			
			const LHAPDF::PDFUncertainty xErr = setk.uncertainty(xfx, -1);
			output2 << i << " " << (xErr.central) << " " << (xErr.errsymm) << " " << xErr.errplus << " " << xErr.errminus << " " << (xfx[0]) << endl;

			if (counter == 100) {
				increment *= 10;
				counter = 1;
			}
			++counter;
		}
	}
	output2.close();
	return 0;
}
