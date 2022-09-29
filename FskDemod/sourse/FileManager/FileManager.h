#ifndef FILEMANAGER_H
#define FILEMANAGER_H
#include <fstream>
#include <vector>
struct WavFormat {
  WavFormat() { this->sampleRate = 44100; }
  WavFormat(uint16_t fs) { this->sampleRate = fs; }
  void wrireByte(std::ofstream &file, uint16_t &val, uint16_t byteSize) {
    file.write(reinterpret_cast<char *>(&val), byteSize);
  }
  /*Head*/
  std::string chunckId = "RIFF";
  std::string chuckSize = "----";
  std::string format = "WAVE";
  /*ftm sub*/
  std::string subChunck1Id = "ftm ";
  uint16_t subChunck1Suze = 16;
  uint16_t audioFormat = 1;
  uint16_t channels = 2;
  uint16_t sampleRate = 0;
  uint16_t byteRate = sampleRate * audioFormat * subChunck1Suze / 8;
  uint16_t blockAlign = channels * subChunck1Suze / 8;
  uint16_t bitsPerSample = subChunck1Suze;
  /*data sub*/
  std::string subChunck2Id = "data";
  std::string subChunck2Suze = "----";
};
namespace FileManager {
template <typename T>
void saveSignal(std::vector<T> &savedSignal, const std::string &fileName) {
  std::ofstream fout;
  const uint32_t len = savedSignal.size();
  fout.open(fileName, std::ios::trunc | std::ios::binary);
  fout.write(reinterpret_cast<char *>(&savedSignal[0]), (len * sizeof(T)));
  fout.close();
}
template <typename T> std::vector<T> readBin(const std::string &fileName) {
  std::ifstream fin;
  fin.open(fileName, std::ios::binary);
  fin.seekg(0, std::ios::end);
  uint32_t actualLength = fin.tellg() / sizeof(T);
  std::vector<T> out(actualLength);
  fin.seekg(0, std::ios::beg);
  fin.read(reinterpret_cast<char *>(&out[0]), actualLength * sizeof(T));
  fin.close();
  return out;
}

template <typename T> std::vector<T> readWav(const std::string &fileName) {
  std::ifstream fin;
  fin.open(fileName, std::ios::binary);
  fin.seekg(-44, std::ios::end);
  uint32_t actualLength = fin.tellg() / sizeof(T);
  std::vector<T> out(actualLength);
  fin.seekg(44, std::ios::beg);
  fin.read(reinterpret_cast<char *>(&out[0]), actualLength * sizeof(T));
  fin.close();
  return out;
}

void saveWav(std::uint32_t fs, const std::string &fileName,
             std::vector<float> &savedSignal) {
  std::ofstream wav;
  wav.open(fileName, std::ios::trunc | std::ios::binary);
  WavFormat file(fs);
  /*Golova*/
  wav << file.chunckId;
  wav << file.chuckSize;
  wav << file.format;
  /*sub */
  wav << file.subChunck1Id;
  file.wrireByte(wav, file.subChunck1Suze, 4);
  file.wrireByte(wav, file.audioFormat, 2);
  file.wrireByte(wav, file.channels, 2);
  file.wrireByte(wav, file.sampleRate, 4);
  file.wrireByte(wav, file.byteRate, 4);
  file.wrireByte(wav, file.blockAlign, 2);
  file.wrireByte(wav, file.bitsPerSample, 2);
  /*data*/
  wav << file.subChunck2Id;
  wav << file.subChunck2Suze;
  uint16_t start = wav.tellp();
  const uint32_t len = savedSignal.size();
  wav.write(reinterpret_cast<char *>(&savedSignal[0]), (len * sizeof(float)));
  uint16_t end = wav.tellp();
  uint16_t dur = end - start;
  wav.seekp(start - 4);
  wav.write(reinterpret_cast<char *>(&dur), 4);
  wav.seekp(4, std::ios::beg);
  uint16_t temp = end - 8;
  wav.write(reinterpret_cast<char *>(&temp), 4);
  wav.close();
}
};     // namespace FileManager
#endif // FILEMANAGER_H
