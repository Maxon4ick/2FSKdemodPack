#include "../anotherFun.h"
#include <fstream>
#include <gtest/gtest.h>
#include <vector>

template <typename t> void EQUAL_VEC(std::vector<t> &test, std::vector<t> &ex) {
  t num = 0;
  int size = test.size();
  for (int i = 0; i < size; i++) {
    if (test[i] == ex[i]) {
      num++;
    }
  }
  assert(num == size);
}

TEST(AnoTests, getPow) { assert(HelpFun::retPow(960) == 1024); }
TEST(AnoTests, sum) {
  std::vector<short> a(100, 1);
  assert(HelpFun::sum(a) == 100);
}
TEST(AnoTests, cor) {
  std::vector<short> a(100, 1);
  std::vector<short> b(100, 1);
  assert(HelpFun::cor(a, b) == 100);
}
TEST(AnoTests, max) {
  std::vector<short> a = {1, 2, 3, 4, 5};
  assert(HelpFun::max(a) == 5);
}
TEST(AnoTests, minusVec) {
  std::vector<int> a = {1, 2, 3, 4, 5};
  std::vector<int> b = {0, 1, 2, 3, 4};
  std::vector<int> ex(5, 1);
  std::vector<int> test = HelpFun::minus(a, b);
  EQUAL_VEC(test, ex);
}
TEST(AnoTests, minusDig) {
  std::vector<int> a = {1, 2, 3, 4, 5};
  int b = 1;
  std::vector<int> ex = {0, 1, 2, 3, 4};
  std::vector<int> test = HelpFun::minus(a, b);
  EQUAL_VEC(test, ex);
}
TEST(AnoTests, strob) {
  std::vector<int> a(100, 1);
  std::vector<int> ex(10, 1);
  std::vector<int> test = HelpFun::getStrob(a, 0, 10);
  EQUAL_VEC(test, ex);
}
TEST(AnoTests, nonZero) {
  std::vector<int> a = {1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0};
  std::vector<int> ex(5, 1);
  std::vector<int> test = HelpFun::nonZero(a);
  EQUAL_VEC(test, ex);
}
