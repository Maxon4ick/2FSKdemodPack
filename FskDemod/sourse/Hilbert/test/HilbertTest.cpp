#include "../hilbert.h"
#include "../../FileManager/FileManager.h"
#include <fstream>
#include <gtest/gtest.h>
#include <vector>
template <typename t> void EQUAL_VEC(std::vector<t> &test, std::vector<t> &ex) {
  int num = 0;
  for (int i = 0; i < test.size(); i++) {
    if (test[i] == ex[i]) {
      num += 1;
    }
  }
  assert(num == test.size());
}
TEST(HilbertTest, FilOgib) {
  std::vector<short> input(120, 1);
  double FrequencyBandsLo[2] = {0.15, 0.25};
  double FrequencyBandsHi[2] = {0.27, 0.37};
  int oldsize = input.size();
  int filOrd = 3;
  std::vector<double> filt =
      filtOgib(input, oldsize, filOrd, FrequencyBandsLo[0], FrequencyBandsLo[1],
               FrequencyBandsHi[0], FrequencyBandsHi[1]);
  const std::string fileNameTest = "TestFilHil.txt";
  std::vector<double> test = FileManager::readBin<double>(fileNameTest);
  EQUAL_VEC(filt, test);
}
