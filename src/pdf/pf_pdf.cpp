#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "LHAPDF/LHAPDF.h"
#include "deriv_pdf.h"
#include "parameters.h"
using namespace std;

////////////////////////////////////////////////////////////////////////
///
/// this file contains all pdf related stuff, like sums of qqbar
/// and derivatives of pdfs
/// note that the xfxQ function of LHAPDF returns the x*f(x) value
/// use this only for the prompt photon stuff
///
////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////
/// this is the qqbar sum with x1 and x2 inputs
/////////////////////////////////////////////////////////////////////////

double pf_sum_qqbar_charge_weighted(double x1, double x2){
	double sum_pdf(0);
	double eq[5] = {-1./3.,2./3,-1./3.,2./3.,-1./3.}; //these are the charges
	for(int i = 1; i <=5; i++){
		sum_pdf+= eq[i-1]*eq[i-1]*(pdfs[0]->xfxQ(i,x1,muF)*pdfs[0]->xfxQ(-i,x2,muF)+pdfs[0]->xfxQ(i,x2,muF)*pdfs[0]->xfxQ(-i,x1,muF));
		//cout << "in the sum " << sum_pdf << endl;
	}
	return sum_pdf;
}
