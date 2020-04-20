#include <vector>

#ifndef SUSY_H
#define SUSY_H

//for SUSY dihiggs;
extern double SUSYset;
extern double mHeavy2;
extern double mA2;
extern double GammaHeavy;
extern double tb;
extern double ght;
extern double gHt;
extern double gAt;
extern double ghb;
extern double gHb;
extern double gAb;
extern double lambdaZAh;
extern double lambdaZAH;
extern double lambdahhh;
extern double lambdaHhh;
extern double lambdaHHh;
extern double lambdaHHH;
extern double lambdahAA;
extern double lambdaHAA;

extern std::unordered_map<int, std::vector< std::vector<double> >> SUSYparameters;
void set_SUSY(float mA, float TANB, bool printout = true);

#endif
