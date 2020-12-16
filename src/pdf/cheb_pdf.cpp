#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <string>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include "cuba.h"
#include "parameters.h"
#include "cuba_integration.h"
#include "deriv_pdf.h"


double Lum_func (double x, std::string LumChannel)
{
	tau = exp(x);
	std::vector<results_c> lumni;
	lumni = call_lumni(LumChannel);
    return exp(x)*lumni[0].res;
}

double easy_func (double x)
{
	//return pow(x,2);
	if(x==1){return 0;}
	return log10(pdfs[use_member]->xfxQ(1,pow(10,x),muF));//(exp(x)*(1.-exp(x)));
}

int **Tnk;
unsigned int current_Nmax=0;

// coefficients T_nk of Chebyshev expansion
int T(int n, int k) {
  int value=1;
  if(k>n || n<0 || k<0) {
    value = 0;
  } else if(((n-k)%2)==1) {
    value = 0;
  } else if(n>0) {
    if(k==0) {
      value= -Tnk[n-2][k];
    } else if(k==n && n!=1) {
      value = 2*Tnk[n-1][k-1];
    } else if(k<n) {
      value = 2*Tnk[n-1][k-1]-Tnk[n-2][k];
    }
  }
  return value;
}

int initT(unsigned int N) {
  if(N<=current_Nmax) return 1;
  else {
    current_Nmax = N;
    Tnk = new int*[N+1];
    for(unsigned int n=0; n<=N; n++) {
      Tnk[n]= new int[N+1];
      for(unsigned int k=0; k<=N; k++) {
	Tnk[n][k]=T(n,k);
      }
    }
    return 1;
  }
  return 0;
}

double ctilde(double *c, int k, int N) {
  double value=0;
  for(int n=k; n<=N; n++){
    value += c[n]*Tnk[n][k];
  }
  if(k==0) value -= 0.5*c[0];
  return value;
}
double cbar(double *c, int p, int N, double umin) {
  double value=0;
  for(int k=p; k<=N; k++){
    value += (c[k]*gsl_sf_fact(k))/gsl_sf_fact(k-p);
  }
  return value*pow(2./umin,p);
}

// routine adapted from GSL library
int cheb_approx_try(unsigned int Nch, const double umin, const double umax, double *c) {
  if(umin >= umax) {
    std::cout << "ERROR: in Chebyshev approx, wrong interval [umin,umax]" << std::endl;
    abort();
  }
  //
  double bma = 0.5*(umax-umin);
  double bpa = 0.5*(umax+umin);
  double fac = 2./(Nch+1.);
  double *f = new double[Nch+1];
  //
  for(unsigned int k=0; k<=Nch; k++) {
    double y = cos(M_PI*(k+0.5)/(Nch+1.));
    f[k] = easy_func(y*bma + bpa);
  }
  //
  for(unsigned int j=0; j<=Nch; j++) {
    double sum = 0.;
    for(unsigned int k=0; k<=Nch; k++)
      sum += f[k]*cos(M_PI*j*(k+0.5)/(Nch+1.));
    c[j] = fac * sum;
  }
  return 1;
}

// routine adapted from GSL library
int cheb_approx(unsigned int N, const double a, const double b, double *c, std::string LumChannel) {
  if(a >= b) {
    std::cout << "ERROR: in Chebyshev approx, wrong interval [a,b]" << std::endl;
    abort();
  }
  //
  double bma = 0.5*(b-a);
  double bpa = 0.5*(b+a);
  double fac = 2./(N+1.);
  double *f = new double[N+1];
  //
  for(unsigned int k=0; k<=N; k++) {
    double y = cos(M_PI*(k+0.5)/(N+1.));
    f[k] = Lum_func(y*bma + bpa,LumChannel);
  }
  //
  for(unsigned int j=0; j<=N; j++) {
    double sum = 0.;
    for(unsigned int k=0; k<=N; k++)
      sum += f[k]*cos(M_PI*j*(k+0.5)/(N+1.));
    c[j] = fac * sum;
  }
  return 1;
}

double Fueval(double u, double umin, double umax, unsigned int N, double *ccoeff)
{
  size_t i;
  double d1 = 0.0;
  double d2 = 0.0;

  double y = (2.0 * u - umin - umax) / (umax - umin);
  double y2 = 2.0 * y;

  for (i = N; i >= 1; i--)
    {
      double temp = d1;
      d1 = y2 * d1 - d2 + ccoeff[i];
      d2 = temp;
    }

  return y * d1 - d2 + 0.5 * ccoeff[0];
}

int ApproxLuminosity(double *coeff, double tau_min, unsigned int N, std::string LumChannel) {
  initT(N);
  //
	std::cout << "Calculating cheb PDF " << std::endl;
  double *coeff_tmp   = new double[N+1];
  double *coeff_tilde = new double[N+1];
  cheb_approx(N, log(tau_min), 0., coeff_tmp, LumChannel);

  for(unsigned int k=0; k<=N; k++){
	std::cout << coeff_tmp[k] << std::endl;
	}
  std::cout << "========== CHECK ==========" << std::endl;
  std::cout << Lum_func(log(tau_min), LumChannel) <<  " " << Fueval(log(tau_min), log(tau_min), 0., N, coeff_tmp) << std::endl;

	for(unsigned int k=0; k<=N; k++) {
    coeff_tilde[k] = ctilde(coeff_tmp,k,N);
  }
  for(unsigned int k=0; k<=N; k++) {
    coeff[k] = cbar(coeff_tilde,k,N,log(tau_min));
  }
  return 1;
}



std::complex<double> LumN(std::complex<double> MelMoment, unsigned int N, std::string LumChannel) {
 std::complex<double> lumniN = 0.;
 double *fitcheb_coeff;
 if (LumChannel == "gg"){fitcheb_coeff = fitcheb_coeff_gg;}
 else if (LumChannel == "qg"){fitcheb_coeff = fitcheb_coeff_qg;}
 else if (LumChannel == "qg_charge"){fitcheb_coeff = fitcheb_coeff_qg_charge;}
 else if (LumChannel == "qqbar"){fitcheb_coeff = fitcheb_coeff_qqbar;}
 else if (LumChannel == "qqbarU"){fitcheb_coeff = fitcheb_coeff_qqbarU;}
 else if (LumChannel == "qqbarD"){fitcheb_coeff = fitcheb_coeff_qqbarD;}
 else if (LumChannel == "qqbarH"){fitcheb_coeff = fitcheb_coeff_qqbar;}
 else{ std::cout << "Channel not implemented " << std::endl; exit(0);}
	for(unsigned int p = 0; p<=N; p++){
		lumniN += fitcheb_coeff[p]/pow(MelMoment-1.,p+1.);
	}
  return lumniN;
}
