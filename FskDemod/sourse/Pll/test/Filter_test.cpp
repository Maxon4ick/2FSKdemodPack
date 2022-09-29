#include "../../FileManager/FileManager.h"
#include "../filter.h"
#include <fstream>
#include <gtest/gtest.h>
#include <vector>
template <typename t> void EQUAL_VEC(t *test, t *ex, int size) {
  t num = 0;
  for (int i = 0; i < size; i++) {
    if (test[i] == ex[i]) {
      num += 1;
    }
  }
  assert(num == size);
}
template <typename t> void EQUAL_VEC(std::vector<t> &test, std::vector<t> &ex) {
  int num = 0;
  for (int i = 0; i < test.size(); i++) {
    if (test[i] == ex[i]) {
      num += 1;
    }
  }
  assert(num == test.size());
}
TEST(FilterTest, CompDenomTest) {
  double FrequencyBands[2] = {0.15, 0.25};
  int FiltOrd = 3;
  double *DenC = 0;
  DenC = SignalProccesing::ComputeDenCoeffs(FiltOrd, FrequencyBands[0],
                                            FrequencyBands[1]);
  const std::string fileNameTest = "TestFilDemC.txt";
  std::vector<double> test = FileManager::readBin<double>(fileNameTest);
  double *t = &test[0];
  EQUAL_VEC(DenC, t, 7);
}
TEST(FilterTest, CompNumTest) {
  double FrequencyBands[2] = {0.15, 0.25};
  int FiltOrd = 3;
  double *DenC = 0;
  DenC = SignalProccesing::ComputeDenCoeffs(FiltOrd, FrequencyBands[0],
                                            FrequencyBands[1]);
  double *NumC = 0;
  NumC = SignalProccesing::ComputeNumCoeffs(FiltOrd, FrequencyBands[0],
                                            FrequencyBands[1], DenC);
  const std::string fileNameTest = "TestFilNum.txt";
  std::vector<double> test = FileManager::readBin<double>(fileNameTest);
  double *t = &test[0];
  EQUAL_VEC(NumC, t, 7);
}
TEST(FilterTest, FiltTest) {
  double FrequencyBands[2] = {0.15, 0.25};
  int FiltOrd = 3;
  double *DenC = 0;
  DenC = SignalProccesing::ComputeDenCoeffs(FiltOrd, FrequencyBands[0],
                                            FrequencyBands[1]);
  double *NumC = 0;
  NumC = SignalProccesing::ComputeNumCoeffs(FiltOrd, FrequencyBands[0],
                                            FrequencyBands[1], DenC);
  double inRe[5] = {1, 2, 3, 4, 5};
  double oRe[5];
  SignalProccesing::filter(3, DenC, NumC, 4, inRe, oRe);
  const std::string fileNameTest = "TestFilFil.txt";
  std::vector<double> test = FileManager::readBin<double>(fileNameTest);
  double *t = &test[0];
  EQUAL_VEC(oRe, t, 4);
}
TEST(strob, GenStr) {

  std::vector<double> out = SignalProccesing::getOp();
  const std::string fileNameTest = "TestFilStrob.txt";
  std::vector<double> test = FileManager::readBin<double>(fileNameTest);
  EQUAL_VEC(test, out);
}

TEST(Filt, FiltCom) {
  double FrequencyBands[2] = {0.15, 0.25};
  int FiltOrd = 3;
  std::vector<double> rez(100, 1);
  std::vector<double> imz(100, 1);
  int oldSize = 7;
  std::vector<double> out = SignalProccesing::filtfilt(
      rez, imz, oldSize, FiltOrd, FrequencyBands[0], FrequencyBands[1]);
  const std::string fileNameTest = "TestFilBig.txt";
  std::vector<double> test = FileManager::readBin<double>(fileNameTest);
  EQUAL_VEC(test, out);
}
