#include "tth_softanom.h"
#include "k_factors_ttH.h"
#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "parameters.h"

using namespace std;
//FEED t wiggle to all the functions!
double beta34(double s34)
{
	return sqrt(1.-4.*mt2/s34);
}//check
complex<double> log_beta34(double s34)
{
	double betatt = beta34(s34);
	double result = (1.+pow(betatt,2))/(2.*betatt)*(log((1.-betatt)/(1.+betatt)));
	if(isnan(result) || isinf(result)){return -1.+I*M_PI;}
	return result +I*M_PI;
}//check

complex<double> Tij(complex<double> tij, complex<double> s)
{
	return log(-tij/(mt*sqrt(s)))+(I*M_PI-1.)/2.;
} //check

complex<double> Omega3(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return (Tij(t13, s)+Tij(t24,s)-Tij(t14,s)-Tij(t23,s))/2.;
}//LET OP DIT MOET t13 = (p1+p3)^2-mt^2 zijn!! Zorg dat dat klopt

complex<double> Lambda3(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return (Tij(t13, s)+Tij(t24,s)+Tij(t14,s)+Tij(t23,s))/2.;
}//LET OP DIT MOET t13 = (p1+p3)^2-mt^2 zijn!! Zorg dat dat klopt

complex<double> lambda_qq_11(double s34)
{
	//cout << log_beta34(s34) << endl;
	return -alphas_muR/(M_PI)*CF*(log_beta34(s34)+1.);//+alphas_muR/(M_PI)*CF*(1.-I*M_PI); //klopt dit?
}

complex<double> lambda_qq_12(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return alphas_muR/(M_PI)*(CF)/CA*Omega3(t13,t24,t14,t23,s);
}

complex<double> lambda_qq_21(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return alphas_muR/(M_PI)*2.*Omega3(t13,t24,t14,t23,s);
}

complex<double> lambda_qq_22(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, double s34, complex<double> s)
{
	return alphas_muR/(2.*M_PI)*((CA-2.*CF)*(log_beta34(s34)+1.)+CA*Lambda3(t13,t24,t14,t23,s)+(8.*CF-3.*CA)*Omega3(t13,t24,t14,t23,s));//+alphas_muR/(M_PI)*CF*(1.-I*M_PI);
}

complex<double> lambda_gg_11(double s34)
{
	return -alphas_muR/(M_PI)*CF*(log_beta34(s34)+1.);//+alphas_muR/(M_PI)*CA*(1.-I*M_PI);
}

complex<double> lambda_gg_12()
{
	return 0;
}

complex<double> lambda_gg_13(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return alphas_muR/(M_PI)*Omega3(t13,t24,t14,t23,s);
}

complex<double> lambda_gg_21()
{
	return 0;
}

complex<double> lambda_gg_22(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s, double s34)
{
	return alphas_muR/(2.*M_PI)*(1./CA*(log_beta34(s34)+1.)+CA*Lambda3(t13,t24,t14,t23,s));//+alphas_muR/(M_PI)*CA*(1.-I*M_PI);
}

complex<double> lambda_gg_23(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return alphas_muR/(2.*M_PI)*CA*Omega3(t13,t24,t14,t23,s);
}

complex<double> lambda_gg_31(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return alphas_muR/(M_PI)*2.*Omega3(t13,t24,t14,t23,s);
}

complex<double> lambda_gg_32(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s)
{
	return alphas_muR*(pow(CA,2)-4.)/(2.*CA*M_PI)*Omega3(t13,t24,t14,t23,s);
}

complex<double> lambda_gg_33(complex<double> t13, complex<double> t24, complex<double> t14, complex<double> t23, complex<double> s, double s34)
{
	return alphas_muR/(2.*M_PI)*(1./CA*(log_beta34(s34)+1.)+CA*Lambda3(t13,t24,t14,t23,s));//+alphas_muR/(M_PI)*CA*(1.-I*M_PI);
}

vector<complex<double>> eigenvalues_qq(complex<double> G11, complex<double> G12, complex<double> G21, complex<double> G22){
	vector<complex<double>> lambdapm = {0,0};
	complex<double> sqrtD = sqrt(pow(G11-G22,2)+4.*G21*G12);
	lambdapm[0] = 1./2.*(G11+G22+sqrtD);
	lambdapm[1] = 1./2.*(G11+G22-sqrtD);
	return lambdapm;
}

vector<vector<complex<double>>> eigenvectors_qq(vector<complex<double>> lambdapm, complex<double> G11, complex<double> G12, complex<double> G21, complex<double> G22){
	//vector<complex<double>> lambdapm = eigenvalues_qq(G11, G12, G21, G22);
	vector<complex<double>> v1 = {0,0}, v2 = {0,0};
	v1[0] = G12/(lambdapm[0]-G11); v1[1] = 1.;
	v2[0] = G12/(lambdapm[1]-G11); v2[1] = 1.;
	if(isnan(real(v1[0])) || isinf(real(v1[0]))){v1[0]=0.;}
	if(isnan(real(v2[0])) || isinf(real(v2[0]))){v2[0]=0.;}
	vector<vector<complex<double>>> eigenvectors = {v1,v2};
	return eigenvectors;
}

vector<vector<complex<double>>> R_qq(vector<vector<complex<double>>> eigenvectors){ //eigenvectors_qq(G11,G12,G21,G22)
	vector<vector<complex<double>>> R = {{0,0},{0,0}};
	R[0][0] = eigenvectors[0][0];
	R[0][1] = eigenvectors[1][0];
	R[1][0] = eigenvectors[0][1];
	R[1][1] = eigenvectors[1][1];
	return R;
}

vector<vector<complex<double>>> Rdag_qq(vector<vector<complex<double>>> R){
	vector<vector<complex<double>>> Rdag = {{0,0},{0,0}};
	Rdag[0][0] = conj(R[0][0]);
	Rdag[1][0] = conj(R[0][1]);
	Rdag[0][1] = conj(R[1][0]);
	Rdag[1][1] = conj(R[1][1]);
	return Rdag;
}

vector<vector<complex<double>>> Rinv_qq(vector<vector<complex<double>>> eigenvectors){ //eigenvectors_qq(G11,G12,G21,G22)
	complex<double> Cinv = 1./(eigenvectors[0][0] - eigenvectors[1][0]);
	vector<vector<complex<double>>> R = {{0,0},{0,0}};
	R[0][0] = Cinv*1.; R[0][1] = -Cinv*eigenvectors[1][0]; R[1][0] = -Cinv*1.; R[1][1] = Cinv*eigenvectors[0][0];
	return R;
}


vector<vector<complex<double>>> mult_qq(vector<vector<complex<double>>> M1, vector<vector<complex<double>>> M2, vector<vector<complex<double>>> M3){
	vector<vector<complex<double>>> mult = {{0,0},{0,0}};
  mult[0][0] = M1[0][0]*M2[0][0]*M3[0][0] + M1[0][1]*M2[1][0]*M3[0][0] + M1[0][0]*M2[0][1]*M3[1][0] + M1[0][1]*M2[1][1]*M3[1][0];
	mult[0][1] = M1[0][0]*M2[0][0]*M3[0][1] + M1[0][1]*M2[1][0]*M3[0][1] + M1[0][0]*M2[0][1]*M3[1][1] + M1[0][1]*M2[1][1]*M3[1][1];
	mult[1][0] = M1[1][0]*M2[0][0]*M3[0][0] + M1[1][1]*M2[1][0]*M3[0][0] + M1[1][0]*M2[0][1]*M3[1][0] + M1[1][1]*M2[1][1]*M3[1][0];
	mult[1][1] = M1[1][0]*M2[0][0]*M3[0][1] + M1[1][1]*M2[1][0]*M3[0][1] + M1[1][0]*M2[0][1]*M3[1][1] + M1[1][1]*M2[1][1]*M3[1][1];
	return mult;
}


complex<double> omega_qq(vector<vector<complex<double>>> HR_IJ, vector<vector<complex<double>>> resumC,  vector<vector<complex<double>>> SR_IJ, vector<vector<complex<double>>> resum){
	return resumC[0][0]*HR_IJ[0][0]*resum[0][0]*SR_IJ[0][0]+resumC[0][0]*HR_IJ[1][0]*resum[1][1]*SR_IJ[0][1]+resumC[1][1]*HR_IJ[0][1]*resum[0][0]*SR_IJ[1][0]+resumC[1][1]*HR_IJ[1][1]*resum[1][1]*SR_IJ[1][1];
}

complex<double> omega_gg(vector<vector<complex<double>>> H, vector<vector<complex<double>>> RC, vector<vector<complex<double>>> S, vector<vector<complex<double>>> R){
	return H[0][0]*R[0][0]*RC[0][0]*S[0][0] + H[0][1]*R[0][0]*RC[1][1]*S[1][0] + H[0][2]*R[0][0]*RC[2][2]*S[2][0] 
			+ H[1][0]*R[1][1]*RC[0][0]*S[0][1] + H[1][1]*R[1][1]*RC[1][1]*S[1][1] + H[1][2]*R[1][1]*RC[2][2]*S[2][1] 
			+ H[2][0]*R[2][2]*RC[0][0]*S[0][2] + H[2][1]*R[2][2]*RC[1][1]*S[1][2] + H[2][2]*R[2][2]*RC[2][2]*S[2][2];
}

//correct
vector<complex<double>> eigenvalues_gg(complex<double> G11, complex<double> G13, complex<double> G22, complex<double> G23, complex<double> G31, complex<double> G32, complex<double> G33){
	vector<complex<double>> lambda = {0,0,0};
	complex<double> b = -G11-G22-G33, c = G11*G22+G11*G33+G22*G33-G31*G13-G32*G23, d = (G11*G23*G32+G22*G13*G31-G11*G22*G33);
	complex<double> Delta0 = pow(b,2)-3.*c, Delta1 = 2.*pow(b,3)-9.*b*c+27.*d;
	complex<double> Cin = pow((Delta1+sqrt(pow(Delta1,2)-4.*pow(Delta0,3)))/2.,1./3.);
	complex<double> xi = -1./2.+I*sqrt(3.)/2.;
	lambda[2] = -1./3.*(b+Cin+Delta0/Cin);
	lambda[0] = -1./3.*(b+Cin*xi+Delta0/(xi*Cin));
	lambda[1] = -1./3.*(b+Cin*xi*xi+Delta0/(xi*xi*Cin));
	return lambda;
}
//correct
vector<vector<complex<double>>> eigenvectors_gg(vector<complex<double>> lambda, complex<double> G11, complex<double> G13, complex<double> G22, complex<double> G23,  complex<double> G31, complex<double> G32, complex<double> G33){
	//vector<complex<double>> lambda = eigenvalues_gg(G11, G13, G22, G23, G31, G32, G33);
	vector<complex<double>> v1 = {0,0,0}, v2 = {0,0,0}, v3 = {0,0,0};
	v1[0] = G13/(lambda[0]-G11); v1[1] = G23/(lambda[0]-G22); v1[2] = 1.;
	v2[0] = G13/(lambda[1]-G11); v2[1] = G23/(lambda[1]-G22); v2[2] = 1.;
	v3[0] = G13/(lambda[2]-G11); v3[1] = G23/(lambda[2]-G22); v3[2] = 1.;
	vector<vector<complex<double>>> eigenvectors = {v1,v2,v3};

	if(isnan(real(v1[0])) || isinf(real(v1[0]))){v1[0]=0.;}
	if(isnan(real(v2[0])) || isinf(real(v2[0]))){v2[0]=0.;}
	if(isnan(real(v3[0])) || isinf(real(v2[0]))){v3[0]=0.;}
	if(isnan(real(v1[1])) || isinf(real(v1[1]))){v1[1]=0.;}
	if(isnan(real(v2[1])) || isinf(real(v2[1]))){v2[1]=0.;}
	if(isnan(real(v3[1])) || isinf(real(v2[1]))){v3[1]=0.;}
	return eigenvectors;
}

vector<vector<complex<double>>> R_gg(vector<vector<complex<double>>> eigenvectors){ //eigenvectors_qq(G11,G12,G21,G22)
	vector<vector<complex<double>>> R = {{0,0,0},{0,0,0},{0,0,0}};
	R[0][0] = eigenvectors[0][0];
	R[0][1] = eigenvectors[1][0];
	R[0][2] = eigenvectors[2][0];
	R[1][0] = eigenvectors[0][1];
	R[1][1] = eigenvectors[1][1];
	R[1][2] = eigenvectors[2][1];
	R[2][0] = eigenvectors[0][2];
	R[2][1] = eigenvectors[1][2];
	R[2][2] = eigenvectors[2][2];
	return R;
}

vector<vector<complex<double>>> Rinv_gg(vector<vector<complex<double>>> ev){ //eigenvectors_qq(G11,G12,G21,G22)
	complex<double> Cinv = 1./(ev[0][0]*(ev[1][1]-ev[2][1])-ev[0][1]*(ev[1][0]-ev[2][0])+ev[1][0]*ev[2][1]-ev[1][1]*ev[2][0]);
	vector<vector<complex<double>>> R = {{0,0,0},{0,0,0},{0,0,0}};
	R[0][0] = Cinv*(ev[1][1]-ev[2][1]);
	R[0][1] = Cinv*(ev[2][0]-ev[1][0]);
	R[0][2] = Cinv*(ev[1][0]*ev[2][1]-ev[1][1]*ev[2][0]);
	R[1][0] = Cinv*(ev[2][1]-ev[0][1]);
	R[1][1] = Cinv*(ev[0][0]-ev[2][0]);
	R[1][2] = Cinv*(ev[0][1]*ev[2][0]-ev[0][0]*ev[2][1]);
	R[2][0] = Cinv*(ev[0][1]-ev[1][1]);
	R[2][1] = Cinv*(ev[1][0]-ev[0][0]);
	R[2][2] = Cinv*(ev[0][0]*ev[1][1]-ev[0][1]*ev[1][0]);
	return R;
}

vector<vector<complex<double>>> Rdag_gg(vector<vector<complex<double>>> R){
	vector<vector<complex<double>>> Rdag = {{0,0,0},{0,0,0},{0,0,0}};
	Rdag[0][0] = conj(R[0][0]);
	Rdag[0][1] = conj(R[1][0]);
	Rdag[0][2] = conj(R[2][0]);
	Rdag[1][0] = conj(R[0][1]);
	Rdag[1][1] = conj(R[1][1]);
	Rdag[1][2] = conj(R[2][1]);
	Rdag[2][0] = conj(R[0][2]);
	Rdag[2][1] = conj(R[1][2]);
	Rdag[2][2] = conj(R[2][2]);
	return Rdag;
}


vector<vector<complex<double>>> mult_gg(vector<vector<complex<double>>> M1, vector<vector<complex<double>>> M2, vector<vector<complex<double>>> M3){
	vector<vector<complex<double>>> mult = {{0,0,0},{0,0,0},{0,0,0}};
    mult[0][0] = M1[0][0]*(M2[0][0]*M3[0][0] + M2[0][1]*M3[1][0] + M2[0][2]*M3[2][0]) + M1[0][1]*(M2[1][0]*M3[0][0] + M2[1][1]*M3[1][0] + M2[1][2]*M3[2][0]) + M1[0][2]*(M2[2][0]*M3[0][0] + M2[2][1]*M3[1][0] + M2[2][2]*M3[2][0]);
	mult[0][1] = M1[0][0]*(M2[0][0]*M3[0][1] + M2[0][1]*M3[1][1] + M2[0][2]*M3[2][1]) + M1[0][1]*(M2[1][0]*M3[0][1] + M2[1][1]*M3[1][1] + M2[1][2]*M3[2][1]) + M1[0][2]*(M2[2][0]*M3[0][1] + M2[2][1]*M3[1][1] + M2[2][2]*M3[2][1]);
	mult[0][2] = M1[0][0]*(M2[0][0]*M3[0][2] + M2[0][1]*M3[1][2] + M2[0][2]*M3[2][2]) + M1[0][1]*(M2[1][0]*M3[0][2] + M2[1][1]*M3[1][2] + M2[1][2]*M3[2][2]) + M1[0][2]*(M2[2][0]*M3[0][2] + M2[2][1]*M3[1][2] + M2[2][2]*M3[2][2]);
	mult[1][0] = M1[1][0]*(M2[0][0]*M3[0][0] + M2[0][1]*M3[1][0] + M2[0][2]*M3[2][0]) + M1[1][1]*(M2[1][0]*M3[0][0] + M2[1][1]*M3[1][0] + M2[1][2]*M3[2][0]) + M1[1][2]*(M2[2][0]*M3[0][0] + M2[2][1]*M3[1][0] + M2[2][2]*M3[2][0]);
	mult[1][1] = M1[1][0]*(M2[0][0]*M3[0][1] + M2[0][1]*M3[1][1] + M2[0][2]*M3[2][1]) + M1[1][1]*(M2[1][0]*M3[0][1] + M2[1][1]*M3[1][1] + M2[1][2]*M3[2][1]) + M1[1][2]*(M2[2][0]*M3[0][1] + M2[2][1]*M3[1][1] + M2[2][2]*M3[2][1]);
	mult[1][2] = M1[1][0]*(M2[0][0]*M3[0][2] + M2[0][1]*M3[1][2] + M2[0][2]*M3[2][2]) + M1[1][1]*(M2[1][0]*M3[0][2] + M2[1][1]*M3[1][2] + M2[1][2]*M3[2][2]) + M1[1][2]*(M2[2][0]*M3[0][2] + M2[2][1]*M3[1][2] + M2[2][2]*M3[2][2]);
	mult[2][0] = M1[2][0]*(M2[0][0]*M3[0][0] + M2[0][1]*M3[1][0] + M2[0][2]*M3[2][0]) + M1[2][1]*(M2[1][0]*M3[0][0] + M2[1][1]*M3[1][0] + M2[1][2]*M3[2][0]) + M1[2][2]*(M2[2][0]*M3[0][0] + M2[2][1]*M3[1][0] + M2[2][2]*M3[2][0]);
	mult[2][1] = M1[2][0]*(M2[0][0]*M3[0][1] + M2[0][1]*M3[1][1] + M2[0][2]*M3[2][1]) + M1[2][1]*(M2[1][0]*M3[0][1] + M2[1][1]*M3[1][1] + M2[1][2]*M3[2][1]) + M1[2][2]*(M2[2][0]*M3[0][1] + M2[2][1]*M3[1][1] + M2[2][2]*M3[2][1]);
	mult[2][2] = M1[2][0]*(M2[0][0]*M3[0][2] + M2[0][1]*M3[1][2] + M2[0][2]*M3[2][2]) + M1[2][1]*(M2[1][0]*M3[0][2] + M2[1][1]*M3[1][2] + M2[1][2]*M3[2][2]) + M1[2][2]*(M2[2][0]*M3[0][2] + M2[2][1]*M3[1][2] + M2[2][2]*M3[2][2]);
	return mult;
}
