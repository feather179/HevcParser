#ifndef HEVC_DCT_H
#define HEVC_DCT_H

#include <vector>

void partialButterfly4(std::vector<std::vector<int>> &src, std::vector<std::vector<int>> &dst, int shift);
void partialButterflyInverse4(std::vector<std::vector<int>> &src,
                              std::vector<std::vector<int>> &dst,
                              int shift,
                              const int outputMinimum,
                              const int outputMaximum);
void partialButterfly8(std::vector<std::vector<int>> &src, std::vector<std::vector<int>> &dst, int shift);
void partialButterflyInverse8(std::vector<std::vector<int>> &src,
                              std::vector<std::vector<int>> &dst,
                              int shift,
                              const int outputMinimum,
                              const int outputMaximum);
void partialButterfly16(std::vector<std::vector<int>> &src, std::vector<std::vector<int>> &dst, int shift);
void partialButterflyInverse16(std::vector<std::vector<int>> &src,
                               std::vector<std::vector<int>> &dst,
                               int shift,
                               const int outputMinimum,
                               const int outputMaximum);
void partialButterfly32(std::vector<std::vector<int>> &src, std::vector<std::vector<int>> &dst, int shift);
void partialButterflyInverse32(std::vector<std::vector<int>> &src,
                               std::vector<std::vector<int>> &dst,
                               int shift,
                               const int outputMinimum,
                               const int outputMaximum);

#endif
