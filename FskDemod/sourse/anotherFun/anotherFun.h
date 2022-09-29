#ifndef ANOTHERFUN_H
#define ANOTHERFUN_H
#include <vector>
namespace HelpFun {
template <typename T> T sum(std::vector<T> &in) {
  double sum = 0;
  for (int i = 0; i < in.size(); i++) {
    sum += in[i];
  }
  return sum;
}
template <typename T> T cor(std::vector<T> &x, std::vector<T> &y) {
  T r = 0;
  for (int i = 0; i < x.size(); i++) {
    r += x[i] * y[i];
  }
  return r;
}
template <typename T> std::vector<T> nonZero(std::vector<T> &in) {
  int size = 0;
  for (int i = 0; i < in.size(); i++) {
    if (in[i] != 0)
      size++;
  }
  std::vector<T> out(size);
  int m = 0;
  for (int i = 0; i < in.size(); i++) {
    if (in[i] != 0) {
      out[m] = in[i];
      m++;
    }
  }
  return out;
}
template <typename T> T max(std::vector<T> &in) {
  int max = 0;
  for (int i = 0; i < in.size(); i++) {
    if (in[i] > max)
      max = in[i];
  }
  return max;
}
template <typename T>
std::vector<T> minus(std::vector<T> &Hi, std::vector<T> &Lo) {
  for (int i = 0; i < Hi.size(); i++) {
    Hi[i] -= Lo[i];
  }
  return Hi;
}
template <typename T> std::vector<T> minus(std::vector<T> &Hi, T &Lo) {
  for (int i = 0; i < Hi.size(); i++) {
    Hi[i] -= Lo;
  }
  return Hi;
}
template <typename T> T retPow(T size) {
  int p = 1;
  while (p < size) {
    p *= 2;
  }
  return p;
}
template <typename T>
std::vector<T> getStrob(std::vector<T> &in, int n, int m) {
  /*n = это сдвиг, m длина строба*/
  std::vector<T> out(m);
  for (int i = 0; i < out.size(); i++) {
    out[i] = in[i + n];
  }
  return out;
}

} // namespace HelpFun
#endif // ANOTHERFUN_H
