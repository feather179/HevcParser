#ifndef HEVC_FILTER_H
#define HEVC_FILTER_H

#include "Typedef.h"

#include <memory>
#include <vector>

class HevcFrame;
class HevcCtu;
class HevcCu;

class HevcFilter {
private:
    void deblockingFilter(std::shared_ptr<HevcCtu> ctu, DeblockEdgeDir edgeDir);
    void getTbBoundary(std::shared_ptr<HevcCu> cu,
                       int xB0,
                       int yB0,
                       int log2TrafoSize,
                       int trafoDepth,
                       int filterEdgeFlag,
                       DeblockEdgeDir edgeDir,
                       std::vector<std::vector<uint8_t>> &edgeFlags);
    void getPbBoundary(std::shared_ptr<HevcCu> cu,
                       DeblockEdgeDir edgeDir,
                       std::vector<std::vector<uint8_t>> &edgeFlags);
    void getBoundaryStrength(std::shared_ptr<HevcFrame> frame,
                             std::shared_ptr<HevcCu> cu,
                             DeblockEdgeDir edgeDir,
                             std::vector<std::vector<uint8_t>> &edgeFlags,
                             std::vector<std::vector<uint8_t>> &bS);
    void filterLumaBlockEdge(std::shared_ptr<HevcFrame> frame,
                             std::shared_ptr<HevcCu> cu,
                             int xBl,
                             int yBl,
                             DeblockEdgeDir edgeDir,
                             uint8_t bS);
    void filterChromaBlockEdge(std::shared_ptr<HevcFrame> frame,
                               std::shared_ptr<HevcCu> cu,
                               int xBl,
                               int yBl,
                               DeblockEdgeDir edgeDir,
                               int cIdx,
                               int qpOffset);
    void sao(std::shared_ptr<HevcCtu> ctu, int cIdx);

public:
    void deblockingFilter(std::shared_ptr<HevcFrame> frame);
    void sao(std::shared_ptr<HevcFrame> frame);
};

#endif
