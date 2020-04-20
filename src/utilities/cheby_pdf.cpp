#include <bits/stdc++.h>
#include <stdlib.h>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include "LHAPDF/LHAPDF.h"
#include "cheby_pdf.h"
#include "polygamma.h"
#include "parameters.h"
#include "chebychev.h"
using namespace std;

/*
 * 
double f(double x, void *p)
{ (void)(p);
   return pow(x,4)/(pow(1.-x,2))*pdfs[use_member]->xfxQ(2,x,muF);
}
 *     int n=100;
    gsl_cheb_series *cs = gsl_cheb_alloc (14);

    gsl_function F;

    F.function = f;
    F.params = 0;

    gsl_cheb_init (cs, &F, 0.0, 1.0);
	cout << "order " << gsl_cheb_order(cs) << endl;
	cout << "size " << gsl_cheb_size(cs) << endl;
    double *coeff = gsl_cheb_coeffs(cs);
    for (int i = 0; i < gsl_cheb_size(cs); i++)
    cout << coeff[i] << endl;
    for (int i = 0; i < n; i++)
    {
      double x = i / (10.*(double)n);
	  
	  double tot = xfx_cheby(x, coeff, gsl_cheb_size(cs));
	  // double r10 = gsl_cheb_eval_n (cs, 10, x);
      double r40 = gsl_cheb_eval (cs, x);
      double abserr;
      gsl_cheb_eval_err(cs, x, &r40, &abserr);
      printf ("%g %g %g %g %g %g\n",
              x, pow(x,-4)/(pow(1.-x,-2))*r40, pdfs[use_member]->xfxQ(2,x,muF), tot, real(xfit_pdfs(7, x)),abserr);
    }

    gsl_cheb_free (cs);
    exit(0);
*/
double xfx_cheby(double x, double* ci, int size_ci){
	double tot = 0;
	vector<double> chebs = cheby_xspace(4.,-2.,x);
	for(int i = 0; i < size_ci; i++)
	{tot = tot + ci[i]*chebs[i]*x;}
	return tot;	
}

std::vector<std::complex<double>> cheby_Nspace(double a, double b, complex<double> N){
	return {1./2.*(Gamma(1. - b)*Gamma(-1. - a + N))/Gamma(-a - b + N),
   ((-2. - a + b + N)*Gamma(1. - b)*Gamma(-1. - a + N))/Gamma(1. - a - b + N),
   ((8. + pow(a,2) + pow(b,2) + a*(7. - 6.*b - 2.*N) + (-7. + N)*N + b*(-9. + 6.*N))*
      Gamma(1. - b)*Gamma(-1. - a + N))/Gamma(2. - a - b + N),
   Gamma(1. - b)*(-(Gamma(-1. - a + N)/Gamma(-a - b + N)) + 
      (2.*(18. + pow(a,2) - 27.*b + a*(5. - 6.*b - 2.*N) - 5.*N + pow(3.*b + N,2))*
         Gamma(-a + N))/Gamma(3. - a - b + N)),
   Gamma(1. - b)*(Gamma(-1. - a + N)/Gamma(-a - b + N) + 
      (32.*(-1.+ b)*(6. + pow(a,2) - 5.*b + a*(3. - 2.*b - 2.*N) - 3.*N + 
           pow(b + N,2))*Gamma(-a + N))/Gamma(4. - a - b + N)),
   Gamma(1. - b)*(((-50.- 49.*a + b + 49.*N)*Gamma(-1. - a + N))/
       Gamma(1. - a - b + N) - (400.*Gamma(1. - a + N))/Gamma(2. - a - b + N) + 
      (1120.*Gamma(2. - a + N))/Gamma(3. - a - b + N) - 
      (1280.*Gamma(3. - a + N))/Gamma(4. - a - b + N) + 
      (512.*Gamma(4. - a + N))/Gamma(5. - a - b + N)),
   Gamma(1. - b)*(((72. + 769.*pow(a,2) + pow(b,2) + a*(839.- 70.*b - 1538.*N) - 
           839.*N + 769.*pow(N,2) + b*(-73. + 70.*N))*Gamma(-1. - a + N))/
       Gamma(2. - a - b + N) - (3584.*Gamma(2. - a + N))/Gamma(3. - a - b + N) + 
      (6912.*Gamma(3. - a + N))/Gamma(4. - a - b + N) - 
      (6144.*Gamma(4. - a + N))/Gamma(5. - a - b + N) + 
      (2048.*Gamma(5. - a + N))/Gamma(6. - a - b + N)),
   Gamma(1. - b)*(-(((98. + 97.*a - b - 97.*N)*Gamma(-1. - a + N))/
         Gamma(1. - a - b + N)) + 
      32.*((-49.*Gamma(1. - a + N))/Gamma(2. - a - b + N) + 
         (294.*Gamma(2. - a + N))/Gamma(3. - a - b + N) - 
         (840.*Gamma(3. - a + N))/Gamma(4. - a - b + N) + 
         (1232.*Gamma(4. - a + N))/Gamma(5. - a - b + N) - 
         (896.*Gamma(5. - a + N))/Gamma(6. - a - b + N) + 
         (256.*Gamma(6. - a + N))/Gamma(7. - a - b + N))),
   Gamma(1. - b)*(Gamma(-1. - a + N)/Gamma(-a - b + N) + 
      128.*(-(Gamma(-a + N)/Gamma(1. - a - b + N)) + 
         (21.*Gamma(1. - a + N))/Gamma(2. - a - b + N) - 
         (168.*Gamma(2. - a + N))/Gamma(3. - a - b + N) + 
         (660.*Gamma(3. - a + N))/Gamma(4. - a - b + N) + 
         128.*((-11.*Gamma(4. - a + N))/Gamma(5. - a - b + N) + 
            (13.*Gamma(5. - a + N))/Gamma(6. - a - b + N) - 
            (8.*Gamma(6. - a + N))/Gamma(7. - a - b + N) + 
            (2.*Gamma(7. - a + N))/Gamma(8.- a - b + N)))),
   Gamma(1. - b)*(-(((162. + 161.*a - b - 161.*N)*Gamma(-1. - a + N))/
         Gamma(1. - a - b + N)) + 
      32.*((-135.*Gamma(1. - a + N))/Gamma(2. - a - b + N) + 
         (1386.*Gamma(2. - a + N))/Gamma(3. - a - b + N) + 
         8.*((-891.*Gamma(3. - a + N))/Gamma(4. - a - b + N) + 
            (2574.*Gamma(4. - a + N))/Gamma(5. - a - b + N) - 
            (4368.*Gamma(5. - a + N))/Gamma(6. - a - b + N) + 
            (4320.*Gamma(6. - a + N))/Gamma(7. - a - b + N) - 
            (2304.*Gamma(7. - a + N))/Gamma(8.- a - b + N) + 
            (512.*Gamma(8.- a + N))/Gamma(9.- a - b + N)))),
   Gamma(1. - b)*(Gamma(-1. - a + N)/Gamma(-a - b + N) + 
      8.*((-25.*Gamma(-a + N))/Gamma(1. - a - b + N) + 
         (825.*Gamma(1. - a + N))/Gamma(2. - a - b + N) + 
         32.*((-330.*Gamma(2. - a + N))/Gamma(3. - a - b + N) + 
            (2145.*Gamma(3. - a + N))/Gamma(4. - a - b + N) - 
            (8008.*Gamma(4. - a + N))/Gamma(5. - a - b + N) + 
            (18200.*Gamma(5. - a + N))/Gamma(6. - a - b + N) + 
            256.*((-100.*Gamma(6. - a + N))/Gamma(7. - a - b + N) + 
               (85.*Gamma(7. - a + N))/Gamma(8.- a - b + N) - 
               (40.*Gamma(8.- a + N))/Gamma(9.- a - b + N) + 
               (8.*Gamma(9.- a + N))/Gamma(10.- a - b + N))))),
   Gamma(1. - b)*(-(((242. + 241.*a - b - 241.*N)*Gamma(-1. - a + N))/
         Gamma(1. - a - b + N)) + 
      16.*((-605.*Gamma(1. - a + N))/Gamma(2. - a - b + N) + 
         (9438.*Gamma(2. - a + N))/Gamma(3. - a - b + N) + 
         16.*((-4719.*Gamma(3. - a + N))/Gamma(4. - a - b + N) + 
            (22022.*Gamma(4. - a + N))/Gamma(5. - a - b + N) - 
            (64064.*Gamma(5. - a + N))/Gamma(6. - a - b + N) + 
            128.*((935.*Gamma(6. - a + N))/Gamma(7. - a - b + N) - 
               (1122.*Gamma(7. - a + N))/Gamma(8.- a - b + N) + 
               (836.*Gamma(8.- a + N))/Gamma(9.- a - b + N) - 
               (352.*Gamma(9.- a + N))/Gamma(10.- a - b + N) + 
               (64.*Gamma(10.- a + N))/Gamma(11. - a - b + N))))),
   Gamma(1. - b)*(Gamma(-1. - a + N)/Gamma(-a - b + N) + 
      32.*((-9.*Gamma(-a + N))/Gamma(1. - a - b + N) + 
         (429.*Gamma(1. - a + N))/Gamma(2. - a - b + N) - 
         (8008.*Gamma(2. - a + N))/Gamma(3. - a - b + N) + 
         (77220.*Gamma(3. - a + N))/Gamma(4. - a - b + N) + 
         1024.*((-429.*Gamma(4. - a + N))/Gamma(5. - a - b + N) + 
            (1547.*Gamma(5. - a + N))/Gamma(6. - a - b + N) - 
            (3672.*Gamma(6. - a + N))/Gamma(7. - a - b + N) + 
            (5814.*Gamma(7. - a + N))/Gamma(8.- a - b + N) - 
            (6080.*Gamma(8.- a + N))/Gamma(9.- a - b + N) + 
            (4032.*Gamma(9.- a + N))/Gamma(10.- a - b + N) - 
            (1536.*Gamma(10.- a + N))/Gamma(11. - a - b + N) + 
            (256.*Gamma(11. - a + N))/Gamma(12. - a - b + N)))),
   Gamma(1. - b)*(-(Gamma(-1. - a + N)/Gamma(-a - b + N)) + 
      (338.*Gamma(-a + N))/Gamma(1. - a - b + N) - 
      (18928.*Gamma(1. - a + N))/Gamma(2. - a - b + N) + 
      (416416.*Gamma(2. - a + N))/Gamma(3. - a - b + N) + 
      512.*((-9295.*Gamma(3. - a + N))/Gamma(4. - a - b + N) + 
         (63206.*Gamma(4. - a + N))/Gamma(5. - a - b + N) + 
         32.*((-8619.*Gamma(5. - a + N))/Gamma(6. - a - b + N) + 
            (25194.*Gamma(6. - a + N))/Gamma(7. - a - b + N) - 
            (50388.*Gamma(7. - a + N))/Gamma(8.- a - b + N) + 
            (69160.*Gamma(8.- a + N))/Gamma(9.- a - b + N) + 
            64.*((-1001.*Gamma(9.- a + N))/Gamma(10.- a - b + N) + 
               (598.*Gamma(10.- a + N))/Gamma(11. - a - b + N) - 
               (208.*Gamma(11. - a + N))/Gamma(12. - a - b + N) + 
               (32.*Gamma(12. - a + N))/Gamma(13. - a - b + N))))),
   Gamma(1. - b)*(Gamma(-1. - a + N)/Gamma(-a - b + N) + 
      8.*((-49.*Gamma(-a + N))/Gamma(1. - a - b + N) + 
         (3185.*Gamma(1. - a + N))/Gamma(2. - a - b + N) + 
         64.*((-1274.*Gamma(2. - a + N))/Gamma(3. - a - b + N) + 
            (17017.*Gamma(3. - a + N))/Gamma(4. - a - b + N) + 
            8.*((-17017.*Gamma(4. - a + N))/Gamma(5. - a - b + N) + 
               (88179.*Gamma(5. - a + N))/Gamma(6. - a - b + N) + 
               16.*((-19380.*Gamma(6. - a + N))/Gamma(7. - a - b + N) + 
                  (47481.*Gamma(7. - a + N))/Gamma(8.- a - b + N) - 
                  (81928.*Gamma(8.- a + N))/Gamma(9.- a - b + N) + 
                  (99176.*Gamma(9.- a + N))/Gamma(10.- a - b + N) + 
                  256.*((-322.*Gamma(10.- a + N))/Gamma(11. - a - b + N) + 
                     (175.*Gamma(11. - a + N))/Gamma(12. - a - b + N) - 
                     (56.*Gamma(12. - a + N))/Gamma(13. - a - b + N) + 
                     (8.*Gamma(13. - a + N))/Gamma(14. - a - b + N)))))))};
}


std::vector<double> cheby_xspace(double a, double b, double x){
	return {1./2.*pow(x,-1. - a)/pow(1. - x,b),(pow(x,-1. - a)*(-1.+ 2.*x))/pow(1. - x,b),
   (pow(x,-1. - a)*(-1.+ 2.*pow(-1.+ 2.*x,2)))/pow(1. - x,b),
   (pow(x,-1. - a)*(-3.*(-1.+ 2.*x) + 4.*pow(-1.+ 2.*x,3)))/pow(1. - x,b),
   (pow(x,-1. - a)*(1. - 8.*pow(-1.+ 2.*x,2) + 8.*pow(-1.+ 2.*x,4)))/
    pow(1. - x,b),(pow(x,-1. - a)*
      (5.*(-1.+ 2.*x) - 20.*pow(-1.+ 2.*x,3) + 16.*pow(-1.+ 2.*x,5)))/pow(1. - x,b)
    ,(pow(x,-1. - a)*(-1.+ 18.*pow(-1.+ 2.*x,2) - 48.*pow(-1.+ 2.*x,4) + 
        32.*pow(-1.+ 2.*x,6)))/pow(1. - x,b),
   (pow(x,-1. - a)*(-7.*(-1.+ 2.*x) + 56.*pow(-1.+ 2.*x,3) - 
        112.*pow(-1.+ 2.*x,5) + 64.*pow(-1.+ 2.*x,7)))/pow(1. - x,b),
   (pow(x,-1. - a)*(1. - 32.*pow(-1.+ 2.*x,2) + 160.*pow(-1.+ 2.*x,4) - 
        256.*pow(-1.+ 2.*x,6) + 128.*pow(-1.+ 2.*x,8)))/pow(1. - x,b),
   (pow(x,-1. - a)*(9.*(-1.+ 2.*x) - 120.*pow(-1.+ 2.*x,3) + 
        432.*pow(-1.+ 2.*x,5) - 576.*pow(-1.+ 2.*x,7) + 256.*pow(-1.+ 2.*x,9)))/
    pow(1. - x,b),(pow(x,-1. - a)*
      (-1.+ 50.*pow(-1.+ 2.*x,2) - 400.*pow(-1.+ 2.*x,4) + 
        1120.*pow(-1.+ 2.*x,6) - 1280.*pow(-1.+ 2.*x,8) + 512.*pow(-1.+ 2.*x,10)))
     /pow(1. - x,b),(pow(x,-1. - a)*
      (-11.*(-1.+ 2.*x) + 220.*pow(-1.+ 2.*x,3) - 1232.*pow(-1.+ 2.*x,5) + 
        2816.*pow(-1.+ 2.*x,7) - 2816.*pow(-1.+ 2.*x,9) + 1024.*pow(-1.+ 2.*x,11))
      )/pow(1. - x,b),(pow(x,-1. - a)*
      (1. - 72.*pow(-1.+ 2.*x,2) + 840.*pow(-1.+ 2.*x,4) - 
        3584.*pow(-1.+ 2.*x,6) + 6912.*pow(-1.+ 2.*x,8) - 
        6144.*pow(-1.+ 2.*x,10) + 2048.*pow(-1.+ 2.*x,12)))/pow(1. - x,b),
   (pow(x,-1. - a)*(13.*(-1.+ 2.*x) - 364.*pow(-1.+ 2.*x,3) + 
        2912.*pow(-1.+ 2.*x,5) - 9984.*pow(-1.+ 2.*x,7) + 
        16640.*pow(-1.+ 2.*x,9) - 13312.*pow(-1.+ 2.*x,11) + 
        4096.*pow(-1.+ 2.*x,13)))/pow(1. - x,b),
   (pow(x,-1. - a)*(-1.+ 98.*pow(-1.+ 2.*x,2) - 1568.*pow(-1.+ 2.*x,4) + 
        9408.*pow(-1.+ 2.*x,6) - 26880.*pow(-1.+ 2.*x,8) + 
        39424.*pow(-1.+ 2.*x,10) - 28672.*pow(-1.+ 2.*x,12) + 
        8192.*pow(-1.+ 2.*x,14)))/pow(1. - x,b)};
	}
