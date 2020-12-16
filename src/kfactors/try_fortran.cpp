
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
 
using namespace std;
 
//extern "C" double c2regAM_(double *z,int *nf,double *muf_mu_rat,double *mur_muf_rat);

extern "C" {
  double c2reghk_qq_(double *z, int *nf, double *muf_mu_rat, double *mur_muf_rat);
}

int main(int argc, char ** argv)
{
	double sum = 0;
	double z = 0.5;
	int nF = 5;
	double muf_mu_ratio = 2.;
	double mur_muf_ratio = 1.;
    //sum = c2regAM_(&z, &nF, &muf_mu_rat, &mur_muf_rat);
    sum = c2reghk_qq_(&z, &nF, &muf_mu_ratio, &mur_muf_ratio);
    cout << sum << endl;
 
}
