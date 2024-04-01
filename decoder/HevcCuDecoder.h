#ifndef HEVC_CU_DECODER_H
#define HEVC_CU_DECODER_H

#include "Typedef.h"

#include <cstdint>
#include <memory>

class CabacReader;
class HevcCtu;
class HevcCu;
class HevcTu;
class HevcFrame;

class HevcCuDecoder {
private:
    std::shared_ptr<CabacReader> mCabacReader;

    bool checkNbBlockAvail(std::shared_ptr<HevcFrame> frame, int xAddr, int yAddr, int xNb, int yNb);
    bool checkNbPbAvail(std::shared_ptr<HevcFrame> frame,
                        int xCb,
                        int yCb,
                        int nCbS,
                        int xPb,
                        int yPb,
                        int nPbW,
                        int nPbH,
                        int partIdx,
                        int xNb,
                        int yNb);

    void reconstructIntraCu(std::shared_ptr<HevcCtu> ctu, std::shared_ptr<HevcCu> cu, int cIdx);
    void reconstructInterCu(std::shared_ptr<HevcCtu> ctu, std::shared_ptr<HevcCu> cu, int cIdx);

public:
    void setCabacReader(std::shared_ptr<CabacReader> cabacReader) { mCabacReader = cabacReader; }
    void decodeCodingQuadtree(std::shared_ptr<HevcCtu> ctu, int x0, int y0, int log2CtuSize, int depth);
    void decodeCodingUnit(std::shared_ptr<HevcCtu> ctu, int x0, int y0, int log2CuSize);
    void decodePredictionUnit(
        std::shared_ptr<HevcCtu> ctu, std::shared_ptr<HevcCu> cu, int x0, int y0, int nPbW, int nPbH, int partIdx);
    void decodeTransformTree(std::shared_ptr<HevcCtu> ctu,
                             std::shared_ptr<HevcCu> cu,
                             int x0,
                             int y0,
                             int xBase,
                             int yBase,
                             int log2TrafoSize,
                             int trafoDepth,
                             int blkIdx);
    void decodeTransformUnit(std::shared_ptr<HevcCtu> ctu,
                             std::shared_ptr<HevcCu> cu,
                             std::shared_ptr<HevcTu> tu,
                             int x0,
                             int y0,
                             int xBase,
                             int yBase,
                             int log2TrafoSize,
                             int trafoDepth,
                             int blkIdx);
    void decodeResidualCoding(std::shared_ptr<HevcCtu> ctu,
                              std::shared_ptr<HevcCu> cu,
                              std::shared_ptr<HevcTu> tu,
                              int x0,
                              int y0,
                              int log2TrafoSize,
                              int cIdx);
};

#endif
