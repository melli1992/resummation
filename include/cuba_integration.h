#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>

#ifndef CUHRE_H
#define CUHRE_H


struct results_c{cubareal res; cubareal err; cubareal prob;};
std::vector<results_c> call_cuhre_dy(std::string orde, std::string chan, std::string pow, bool fitted=false, int maxpower = 1, int verbose=0);
std::vector<results_c> call_cuhre_higgs(std::string orde, std::string chan, std::string pow, bool fitted=false, int maxpower = 1, int verbose=0);
std::vector<results_c> call_cuhre_dihiggs(std::string orde, std::string chan, bool fitted=false, int verbose=0);
std::vector<results_c> call_cuhre_diboson(std::string orde, std::string chan, bool fitted=false, int verbose=0);
std::vector<results_c> call_cuhre_dy_bsm(std::string orde, std::string chan, bool fitted=false, int verbose=0);
std::vector<results_c> call_cuhre_test(std::string orde, std::string chan, int k, bool fitted=false, int verbose=0);
std::vector<results_c> call_set_scale(std::string channel, bool fitted=false, int verbose=0);
std::vector<results_c> call_lumni(std::string channel, bool fitted=false, int verbose=0);
#endif
