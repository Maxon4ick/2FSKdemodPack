#ifndef HILBERT_H
#define HILBERT_H

#include "fftw3.h"
#include <complex>
#include <cstring>
#include <math.h>
#include <vector>
#define REAL 0
#define IMAG 1
void hilbert(const double *in, fftw_complex *out, int n);
std::vector<double> filtOgib(std::vector<short> &inRaw, int oldsize,
                             int filtOrd, double FrequencyBandsLol,
                             double FrequencyBandsLoh, double FrequencyBandsHil,
                             double FrequencyBandsHih);

#endif // HILBERT_H
