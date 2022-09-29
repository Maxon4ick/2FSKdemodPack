#ifndef FILTER_H
#define FILTER_H
#include <cmath>
#include <complex>
#include <cstdlib>
#include <vector>
#define PI 3.1415926535897932384626433832795
//#define N 10

namespace SignalProccesing {
double *ComputeLP(int FilterOrder);
double *ComputeHP(int FilterOrder);
double *TrinomialMultiply(int FilterOrder, double *b, double *c);

double *ComputeNumCoeffs(int FilterOrder, double Lcutoff, double Ucutoff,
                         double *DenC);
double *ComputeDenCoeffs(int FilterOrder, double Lcutoff, double Ucutoff);
void filter(int ord, double *a, double *b, int np, double *x, double *y);
std::vector<double> filtfilt(std::vector<double> &inRawRe,
                             std::vector<double> &inRawIm, int oldsize,
                             int FiltOrd, double FrequencyBandsLo,
                             double FrequencyBandsHi);
std::vector<double> getOp();
} // namespace SignalProccesing

#endif // FILTER_H
