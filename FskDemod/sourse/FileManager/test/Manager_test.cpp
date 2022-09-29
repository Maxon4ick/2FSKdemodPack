#include "../FileManager.h"
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

TEST(FileTests, Size) {
  const std::string fileName = "FSK2_8000_50Bd_7_5st5.wav";
  std::vector<short> a = FileManager::readWav<short>(fileName);
  assert(a.size() == 968358);
}
TEST(WavTest, Data) {
  const std::string fileNameData = "FSK2_8000_50Bd_7_5st5.wav";
  const std::string fileNameTest = "Test1.txt";
  std::vector<short> test = FileManager::readBin<short>(fileNameTest);
  std::vector<short> ex = FileManager::readWav<short>(fileNameData);
  EQUAL_VEC(test, ex);
}

TEST(FileTests, WorkingFloat) {
  std::vector<float> sam1(1000, 0);
  FileManager::saveSignal(sam1, "test.bin");
  std::vector<float> sam2 = FileManager::readBin<float>("test.bin");
  EQUAL_VEC(sam1, sam2);
}
