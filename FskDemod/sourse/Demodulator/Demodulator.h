#ifndef DEMOD_H
#define DEMOD_H
#include <cmath>
#include <vector>
class Demodulate {
public:
  virtual std::vector<double> Demod(std::vector<short> &input, double &fs,
                                    double &ubd, double &wLo, double &wHi) = 0;
};
class DemodulatorFSK : public Demodulate {
public:
  std::vector<double> Demod(std::vector<short> &input, double &fs, double &ubd,
                            double &wLo, double &wHi) override;
};
std::vector<double> rovn(std::vector<double> &filt, int cik);
std::vector<double> sinhr(std::vector<double> &sigRovn, int corSize);
std::vector<double> sinhrPoint(std::vector<double> &synx, int corSize);
std::vector<double> decod(std::vector<double> &sigRovn,
                          std::vector<double> &trPoint, int size);
#endif // DEMOD_H
