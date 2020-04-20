#include <iostream>
#include <fstream>
#include <gsl/gsl_math.h>
#include "parameters.h"
#include "mellin_functions.h"
#include "resum_functions.h"
#include "polygamma.hpp"
#include "deriv_pdf.h"
#include "mellin_pdf.h"

using namespace std;
// contains only test functions
///////////////////////////
/// test functions
///////////////////////////
/// full transformation
double vegas_fofx2_full(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
  //cout << k[1] << " " << pow(k[1],Nint+1.) << endl;
  double result = 2.*imag(1./(2*M_PI)*Njac*pow(0.5,-Nint)*pow(k[1],Nint)*(1.-k[1]));
  if (isnan(result)){return 0;}
	else{return result;}
}
/// derivative approach
double vegas_fofx2_deriv(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
  double result = -2.*imag(1./(2*M_PI)*Njac*pow(0.5,-Nint)*pow(k[1],Nint)/Nint*(1.-2.*k[1]));
  if (isnan(result)){return 0;}
	else{return result;}
}
/// deformation of the contour
double vegas_fofx2_defor(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*I;
  complex<double> Njac = 1./pow(1.-k[0],2)*I;
  double wr = k[1]/(1.+k[1]);
  double wjac = 1./(pow(1.+k[1],2));
  complex<double> x = exp(wr/Nint);
  //cout << "x " << x << endl;
	double result = 2.*imag(1./(2*M_PI)*Njac*pow(0.5,-Nint)*wjac*exp(wr)/Nint*x*(1.-x));
  if (isnan(result)){return 0;}
	else{return result;}
}
/// direct N-space result
double vegas_fofx2_Nspace(double *k, size_t dim, void *params){
	(void)(dim);
	(void)(params);
	complex<double> Nint = CMP+k[0]/(1.-k[0])*exp(I*phiMP);
  complex<double> Njac = 1./pow(1.-k[0],2)*exp(I*phiMP);
	return 2.*imag(1./(2*M_PI)*Njac*pow(0.5,-Nint)*1./((Nint*Nint+3.*Nint+2.)));
}
