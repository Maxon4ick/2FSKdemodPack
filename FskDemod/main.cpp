#include "sourse/Demodulator/Demodulator.h"
#include "sourse/FileManager/FileManager.h"
#include <string>

int main() {
  const std::string fileName = "FSK2_8000_50Bd_7_5st5.wav";
  double fs = 8e3;
  double ubd = 50;
  double wLo = 789;
  double wHi = 1234;
  std::vector<short> inRaw = FileManager::readWav<short>(fileName);
  Demodulate *p = new DemodulatorFSK;
  std::vector<double> out = p->Demod(inRaw, fs, ubd, wLo, wHi);
  FileManager::saveSignal(out, "out.bin");
  return 0;
}
