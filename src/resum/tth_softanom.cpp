#include "tth_softanom.h"
#include <bits/stdc++.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "parameters.h"

using namespace std;
//FEED t wiggle to all the functions!
complex<double> beta34(double s34)
{
	return sqrt(1.-(4.*pow(mt,2))/s34);
}//check
complex<double> log_beta34(double s34)
{
	return (1.+pow(beta34(s34),2))/(2.*beta34(s34))*(log((1.-beta34(s34))/(1.+beta34(s34)))+I*M_PI);
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
	return -alphas_muR/(M_PI)*CF*(log_beta34(s34)+1.)+alphas_muR/(M_PI)*CF*(1.-I*M_PI); //klopt dit?
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
	return alphas_muR/(2.*M_PI)*((CA-2.*CF)*(log_beta34(s34)+1.)+CA*Lambda3(t13,t24,t14,t23,s)+(8.*CF-3.*CA)*Omega3(t13,t24,t14,t23,s))+alphas_muR/(M_PI)*CF*(1.-I*M_PI);
}

complex<double> lambda_gg_11(double s34)
{
	return -alphas_muR/(M_PI)*CF*(log_beta34(s34)+1.)+alphas_muR/(M_PI)*CA*(1.-I*M_PI);
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
	return alphas_muR/(2.*M_PI)*(1./CA*(log_beta34(s34)+1.)+CA*Lambda3(t13,t24,t14,t23,s))+alphas_muR/(M_PI)*CA*(1.-I*M_PI);
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
	return alphas_muR/(2.*M_PI)*(1./CA*(log_beta34(s34)+1.)+CA*Lambda3(t13,t24,t14,t23,s))+alphas_muR/(M_PI)*CA*(1.-I*M_PI);
}


complex<double> eigen_qq_1(complex<double>  lambdaqq11, complex<double>  lambdaqq22, complex<double>  lambdaqq12, complex<double>  lambdaqq21)
{
	return 1./2.*(lambdaqq11+lambdaqq22+sqrt(pow(lambdaqq11-lambdaqq22, 2)+4.*lambdaqq12*lambdaqq21));
}

complex<double> eigen_qq_2(complex<double>  lambdaqq11, complex<double>  lambdaqq22, complex<double>  lambdaqq12, complex<double>  lambdaqq21)
{
	return 1./2.*(lambdaqq11+lambdaqq22-sqrt(pow(lambdaqq11-lambdaqq22, 2)+4.*lambdaqq12*lambdaqq21));
}
complex<double> R_qq_1i(complex<double>  lambdaqq11, complex<double>  lambdaqq22, complex<double>  lambdaqq12, complex<double>  lambdaqq21, complex<double> eigeni)
{
	return lambdaqq12/(eigeni-lambdaqq11);
} //=Rinv_qq_22  //=-Rinv_qq_12

complex<double> R_qq_21()
{
	return 1.;
}//=-Rinv_qq_21

complex<double> R_qq_22()
{
	return 1.;
} //=Rinv_qq_11

complex<double> Rinv_qq_const(complex<double>  lambdaqq11, complex<double>  lambdaqq22, complex<double>  lambdaqq12, complex<double>  lambdaqq21, complex<double> eigen1, complex<double> eigen2)
{
	return ((eigen1-lambdaqq11)*(eigen2-lambdaqq11))/(lambdaqq12*(eigen2-eigen1));

}

complex<double> a(complex<double> gamma11, complex<double> gamma22)
{
	return -(gamma11+2.*gamma22);
}

complex<double> b(complex<double> gamma11, complex<double> gamma22, complex<double> gamma31)
{
	return -(-2.*gamma11*gamma22-pow(gamma22,2)+1./16.*pow(gamma31,2)*pow(CA,2)+1./4.*pow(gamma31,2));
}

complex<double> c(complex<double> gamma11, complex<double> gamma22, complex<double> gamma31)
{
	return -(gamma11*pow(gamma22,2)-1./16.*gamma11*pow(gamma31,2)*pow(CA,2)+1./4.*gamma11*pow(gamma31,2)-1./2.*pow(gamma31,2)*gamma22);
}

complex<double> collect(complex<double> a, complex<double> b, complex<double> c)
{
	return 36.*a*b - 108.*c - 8.*pow(a,3)+12.*sqrt(12.*pow(a,3)*c-3.*pow(a,2)*pow(b,2)-54.*a*b*c+12.*pow(b,3)+81.*pow(c,2));
}


complex<double> eigen1_new(complex<double> a, complex<double> b, complex<double> Z)
{
	return 1./6.*pow(Z,1./3.)-6.*(1./3.*b-1./9.*pow(a,2))/pow(Z,1./3.)-1./3.*a;

}

complex<double> eigen2_new(complex<double> a, complex<double> b, complex<double> Z)
{
	return -1./12.*pow(Z,1./3.)+3.*(1./3.*b-1./9.*pow(a,2))/pow(Z,1./3.)-1./3.*a+1./2.*I*sqrt(3)*(1./6.*pow(Z,1./3.)+pow(Z,-1./3.)*6.*(1./3.*b-1./9.*pow(a,2)));

}

complex<double> eigen3_new(complex<double> a, complex<double> b, complex<double> Z)
{
	return -1./12.*pow(Z,1./3.)+pow(Z,-1./3.)*3.*(1./3.*b-1./9.*pow(a,2))-1./3.*a-1./2.*I*sqrt(3)*(1./6.*pow(Z,1./3.)+pow(Z,-1./3.)*6.*(1./3.*b-1./9.*pow(a,2)));

}

complex<double> R_gg_1i(complex<double> lambda31, complex<double> eigeni, complex<double> lambda11) //=vi_1
{
	return lambda31/(2.*(eigeni-lambda11));
} //check


complex<double> R_gg_2i(complex<double> lambda31, complex<double> eigeni, complex<double> lambda22) //=vi_2
{
	return CA*lambda31/(4.*(eigeni-lambda22));
} 
//check

complex<double> Rinv_gg_const(complex<double> v1_1, complex<double> v2_1, complex<double> v3_1, complex<double> v1_2, complex<double> v2_2, complex<double> v3_2)
{
	return (v1_1*v2_2-v1_1*v3_2-v1_2*v2_1+v1_2*v3_1+v2_1*v3_2-v2_2*v3_1);
} //check


complex<double> Check_Roots(complex<double> eigeni,complex<double> gamma11, complex<double> gamma31, complex<double> gamma22)
{
	pow(eigeni, 3)-(gamma11-2.*gamma22)*pow(eigeni,2)-((pow(CA,2)-4.)/16.*pow(gamma31,2)+pow(gamma31,2)/2.-pow(gamma22,2)-2.*gamma11*gamma22)*eigeni+gamma11*pow(gamma22,2)-(pow(CA,2)-4.)/16.*pow(gamma31,2)-pow(gamma31,2)*gamma22/2.;
} //works

void Check_Eigen(complex<double> eigeni, complex<double> vi_1, complex<double> vi_2, complex<double> vi_3, complex<double> gamma31, complex<double> gamma11, complex<double> gamma22, complex <double> gamma13, complex<double> gamma23, complex<double> gamma32, complex<double> gamma33)
{
	cout << "Eigenvalue i: " << eigeni << " gamma11 " << gamma11 << " gamma31 " << gamma31 << " gamma22 " << gamma22 << " CA " << CA << endl;
	cout << "Eigenvectors i: " << vi_1 << "," << vi_2 << "," << vi_3 << endl;
} //worked!

