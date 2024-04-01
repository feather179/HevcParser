#ifndef HEVC_DST_H
#define HEVC_DST_H

#include <vector>

void fastForwardDst(std::vector<std::vector<int>> &src, std::vector<std::vector<int>> &dst, int shift);

void fastInverseDst(std::vector<std::vector<int>> &src,
                    std::vector<std::vector<int>> &dst,
                    int shift,
                    const int outputMinimum,
                    const int outputMaximum);

#endif
