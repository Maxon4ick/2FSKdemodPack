#include "hilbert.h"
#include "../Pll/filter.h"
#include "../anotherFun/anotherFun.h"
void hilbert(const double *in, fftw_complex *out, int n) {
  // copy the data to the complex array
  for (int i = 0; i < n; ++i) {
    out[i][REAL] = in[i];
    out[i][IMAG] = 0;
  }
  // creat a DFT plan and execute it
  fftw_plan plan = fftw_plan_dft_1d(n, out, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  // destroy a plan to prevent memory leak
  fftw_destroy_plan(plan);
  int hN = n >> 1; // half of the length (N/2)
  int numRem = hN; // the number of remaining elements
  // multiply the appropriate value by 2
  //(those should multiplied by 1 are left intact because they wouldn't change)
  for (int i = 1; i < hN; ++i) {
    out[i][REAL] *= 2;
    out[i][IMAG] *= 2;
  }
  // if the length is even, the number of the remaining elements decrease by 1
  if (n % 2 == 0)
    numRem--;
  else if (n > 1) {
    out[hN][REAL] *= 2;
    out[hN][IMAG] *= 2;
  }
  // set the remaining value to 0
  // (multiplying by 0 gives 0, so we don't care about the multiplicands)
  memset(&out[hN + 1][REAL], 0, numRem * sizeof(fftw_complex));
  // creat a IDFT plan and execute it
  plan = fftw_plan_dft_1d(n, out, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  // do some cleaning
  fftw_destroy_plan(plan);
  fftw_cleanup();
  // scale the IDFT output
  for (int i = 0; i < n; ++i) {
    out[i][REAL] /= n;
    out[i][IMAG] /= n;
  }
}

std::vector<double> filtOgib(std::vector<short> &inRaw, int oldsize,
                             int filtOrd, double FrequencyBandsLol,
                             double FrequencyBandsLoh, double FrequencyBandsHil,
                             double FrequencyBandsHih) {
  int n = HelpFun::retPow(inRaw.size());
  // inRaw.resize(n);
  double *in = new double[n];
  for (int i = 0; i < oldsize; i++) {
    in[i] = static_cast<double>(inRaw[i]);
  }
  fftw_complex *out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * n);

  hilbert(in, out, n);
  std::vector<double> outRe(n);
  std::vector<double> outIm(n);

  for (int i = 0; i < n; i++) {
    outRe[i] = out[i][REAL];
    outIm[i] = out[i][IMAG];
  }
  delete[] in;
  fftw_free(out);
  std::vector<double> outLo = SignalProccesing::filtfilt(
      outRe, outIm, oldsize, filtOrd, FrequencyBandsLol, FrequencyBandsLoh);
  std::vector<double> outHi = SignalProccesing::filtfilt(
      outRe, outIm, oldsize, filtOrd, FrequencyBandsHil, FrequencyBandsHih);
  outRe.clear();
  outIm.clear();
  std::vector<double> filout = HelpFun::minus(outHi, outLo);
  return filout;
}
