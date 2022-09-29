#include "Demodulator.h"
#include "../Hilbert/hilbert.h"
#include "../Pll/filter.h"
#include "../anotherFun/anotherFun.h"
std::vector<double> rovn(std::vector<double> &filt, int cik) {
  std::vector<double> strob(1200);
  int nach = 0;
  double delta = 0;
  std::vector<double> sigRovn(1200 * cik);
  for (int i = 0; i < cik; i++) {
    nach = i * 1200;
    strob = HelpFun::getStrob(filt, nach, 1200);
    delta = HelpFun::sum(strob) / 1200;
    strob = HelpFun::minus(strob, delta);
    for (int i = 0; i < 1200; i++) {
      sigRovn[i + nach] = strob[i];
    }
  }
  return sigRovn;
}
std::vector<double> sinhr(std::vector<double> &sigRovn, int corSize) {
  std::vector<double> strobOp(4800);
  std::vector<double> synx(corSize);
  std::vector<double> op = SignalProccesing::getOp();
  for (int i = 0; i < corSize; i++) {
    strobOp = HelpFun::getStrob(sigRovn, i, 4800);
    synx[i] = HelpFun::cor(strobOp, op);
  }
  return synx;
}
std::vector<double> sinhrPoint(std::vector<double> &synx, int corSize) {
  int nach = 0;
  int noc = corSize / 1200;
  std::vector<double> strob(1200);
  std::vector<double> maxAr(noc);
  for (int i = 0; i < noc; i++) {
    nach = i * 1200;
    strob = HelpFun::getStrob(synx, nach, 1200);
    maxAr[i] = HelpFun::max(strob);
  }
  std::vector<double> point(noc, 0);
  for (int i = 0; i < noc; i++) {
    nach = i * 1200;
    strob = HelpFun::getStrob(synx, nach, 1200);
    for (int j = 0; j < 1200; j++) {
      if ((strob[j] == maxAr[i]) & (strob[j] > 10e5))
        point[i] = i * 1200 + j;
    }
  }
  return point;
}
std::vector<double> decod(std::vector<double> &sigRovn,
                          std::vector<double> &trPoint, int size) {
  std::vector<double> strob(1200);
  std::vector<double> sigDecod(5 * size, 0);
  std::vector<double> strobMa(160);
  double su = 0;
  int b = 0;
  int sn = 0;
  for (int i = 0; i < size; i++) {
    strob = HelpFun::getStrob(sigRovn, trPoint[i], 1200);
    for (int j = 0; j < 5; j++) {
      b = (j + 1) * 160;
      sn = i * 5;
      strobMa = HelpFun::getStrob(strob, b, 160);
      su = HelpFun::sum(strobMa);
      if (su > 0)
        sigDecod[sn + j] = 1;
    }
  }
  return sigDecod;
}
std::vector<double> DemodulatorFSK::Demod(std::vector<short> &input, double &fs,
                                          double &ubd, double &wLo,
                                          double &wHi) {
  double FrequencyBandsLo[2] = {0.15, 0.25};
  double FrequencyBandsHi[2] = {0.27, 0.37};
  int oldsize = input.size();
  int filOrd = 3;
  std::vector<double> filt =
      filtOgib(input, oldsize, filOrd, FrequencyBandsLo[0], FrequencyBandsLo[1],
               FrequencyBandsHi[0], FrequencyBandsHi[1]);
  input.clear();
  int cik = std::round(oldsize / 1200);
  std::vector<double> sigRovn = rovn(filt, cik);
  for (int i = 0; i < sigRovn.size(); i++) {
    if (sigRovn[i] < 0)
      sigRovn[i] = -1000;
    else
      sigRovn[i] = 1000;
  }
  filt.clear();
  int corSize = sigRovn.size() - 4800;
  std::vector<double> synx = sinhr(sigRovn, corSize);
  std::vector<double> point = sinhrPoint(synx, corSize);
  std::vector<double> trPoint = HelpFun::nonZero(point);
  std::vector<double> sigDecod = decod(sigRovn, trPoint, trPoint.size());
  return sigDecod;
}
