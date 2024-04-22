#include "HevcFrame.h"

#include "ParameterSet.h"
#include "HevcSlice.h"
#include "HevcCtu.h"
#include "HevcCu.h"
#include "HevcPu.h"

#include <assert.h>

#include <algorithm>

using std::vector;

void HevcFrame::init(std::shared_ptr<ParameterSet> paramSet) {
    auto sliceHeader = paramSet->getSliceHeader();
    auto pps = mPps = paramSet->getPps(sliceHeader->slice_pic_parameter_set_id);
    auto sps = mSps = paramSet->getSps(pps->pps_seq_parameter_set_id);
    mWidth = sps->width;
    mHeight = sps->height;

    auto slice = std::make_shared<HevcSlice>(sliceHeader);

    mCtus.clear();
    for (int i = 0; i < (int)sps->ctbHeight; i++) {
        for (int j = 0; j < (int)sps->ctbWidth; j++) {
            int rsAddr = i * sps->ctbWidth + j;
            auto ctu =
                std::make_shared<HevcCtu>(j << sps->log2CtbSize, i << sps->log2CtbSize, rsAddr, shared_from_this());
            ctu->setSlice(slice);
            mCtus.push_back(ctu);
        }
    }

    uint32_t minCtbWidth = sps->ctbWidth << (sps->log2CtbSize - sps->log2MinCtbSize);
    uint32_t minCtbHeight = sps->ctbHeight << (sps->log2CtbSize - sps->log2MinCtbSize);
    mCtDepth.clear();
    mCtDepth.resize(minCtbHeight, vector<uint8_t>(minCtbWidth, 0));

    int subWidthC = mSps->subWidthC, subHeightC = mSps->subHeightC;
    mReconY.clear();
    mReconCb.clear();
    mReconCr.clear();
    mReconY.resize(mHeight, vector<uint16_t>(mWidth, 0));
    mReconCb.resize(mHeight / subHeightC, vector<uint16_t>(mWidth / subWidthC, 0));
    mReconCr.resize(mHeight / subHeightC, vector<uint16_t>(mWidth / subWidthC, 0));

    mSaoY.clear();
    mSaoCb.clear();
    mSaoCr.clear();
    mSaoY.resize(mHeight, vector<uint16_t>(mWidth, 0));
    mSaoCb.resize(mHeight / subHeightC, vector<uint16_t>(mWidth / subWidthC, 0));
    mSaoCr.resize(mHeight / subHeightC, vector<uint16_t>(mWidth / subWidthC, 0));

    mSlices.clear();
    mSliceStartCtbRsAddr.clear();
    for (int i = 0; i < 2; i++) {
        mRefPicList[i].clear();
        mRefPocList[i].clear();
    }
}

std::shared_ptr<HevcCtu> HevcFrame::getCtu(int rsAddr) {
    assert(rsAddr < mCtus.size());
    return mCtus[rsAddr];
}

std::shared_ptr<HevcCtu> HevcFrame::getCtu(int xAddr, int yAddr) {
    int x = xAddr >> mSps->log2CtbSize;
    int y = yAddr >> mSps->log2CtbSize;
    return mCtus[y * mSps->ctbWidth + x];
}

std::shared_ptr<HevcCu> HevcFrame::getCu(int xAddr, int yAddr) {
    xAddr = std::clamp(xAddr, 0, mWidth - 1);
    yAddr = std::clamp(yAddr, 0, mHeight - 1);
    int x = xAddr >> mSps->log2CtbSize;
    int y = yAddr >> mSps->log2CtbSize;
    int ctuRsAddr = y * mSps->ctbWidth + x;
    auto ctu = mCtus[ctuRsAddr];
    return ctu->getCu(xAddr, yAddr);
}

std::shared_ptr<HevcPu> HevcFrame::getPu(int xAddr, int yAddr) {
    auto cu = getCu(xAddr, yAddr);
    return cu->getPu(xAddr, yAddr);
}

void HevcFrame::setCtDepth(int xAddr, int yAddr, int log2CbSize, uint8_t depth) {
    int x = xAddr >> mSps->log2MinCtbSize;
    int y = yAddr >> mSps->log2MinCtbSize;
    int count = 1 << (log2CbSize - mSps->log2MinCtbSize);

    for (int i = 0; i < count; i++)
        for (int j = 0; j < count; j++) mCtDepth[y + i][x + j] = depth;
}

uint8_t HevcFrame::getCtDepth(int xAddr, int yAddr) {
    int x = xAddr >> mSps->log2MinCtbSize;
    int y = yAddr >> mSps->log2MinCtbSize;

    return mCtDepth[y][x];
}

// top-left point address of the slice contains (xAddr, yAddr)
PointAddr HevcFrame::getSliceAddr(int xAddr, int yAddr) {
    int x = xAddr >> mSps->log2CtbSize;
    int y = yAddr >> mSps->log2CtbSize;
    int rsAddr = mCtus[mSps->ctbWidth * y + x]->getSliceRsAddr();
    auto ctu = mCtus[rsAddr];

    return {ctu->getX(), ctu->getY()};
}

// top-left point address of the tile contains (xAddr, yAddr)
PointAddr HevcFrame::getTileAddr(int xAddr, int yAddr) {
    int x = xAddr >> mSps->log2CtbSize;
    int y = yAddr >> mSps->log2CtbSize;

    int tileX = 0, tileY = 0;
    for (int i = 0; i < (mPps->tileColumnBoundary.size() - 1); i++) {
        if (x >= (int)mPps->tileColumnBoundary[i] && x < (int)mPps->tileColumnBoundary[i + 1]) {
            tileX = mPps->tileColumnBoundary[i];
            break;
        }
    }
    for (int i = 0; i < (mPps->tileRowBoundary.size() - 1); i++) {
        if (y >= (int)mPps->tileRowBoundary[i] && y < (int)mPps->tileRowBoundary[i + 1]) {
            tileY = mPps->tileRowBoundary[i];
            break;
        }
    }

    x = tileX << mSps->log2CtbSize;
    y = tileY << mSps->log2CtbSize;

    return {x, y};
}

// top-left point address of the quantization group contains (xAddr, yAddr)
PointAddr HevcFrame::getQGAddr(int xAddr, int yAddr) {
    int shift = mPps->log2MinCuQpDeltaSize;

    return {(xAddr >> shift) << shift, (yAddr >> shift) << shift};
}

// clang-format off
static const int gZscanToRaster2x2[] = {
    0, 1,
    2, 3,
};

static const int gZscanToRaster4x4[] = {
    0,  1,  4,  5,
    2,  3,  6,  7,
    8,  9,  12, 13,
    10, 11, 14, 15,
};

static const int gZscanToRaster8x8[] = {
    0,  1,  8,  9,  2,  3,  10, 11,
    16, 17, 24, 25, 18, 19, 26, 27,
    4,  5,  12, 13, 6,  7,  14, 15,
    20, 21, 28, 29, 22, 23, 30, 31,
    32, 33, 40, 41, 34, 35, 42, 43,
    48, 49, 56, 57, 50, 51, 58, 59,
    36, 37, 44, 45, 38, 39, 46, 47,
    52, 53, 60, 61, 54, 55, 62, 63,
};

static const int gRasterToZscan2x2[] = {
    0, 1,
    2, 3,
};

static const int gRasterToZscan4x4[] = {
    0,  1,  4,  5,
    2,  3,  6,  7,
    8,  9,  12, 13,
    10, 11, 14, 15,
};

static const int gRasterToZscan8x8[] = {
    0,  1,  4,  5,  16, 17, 20, 21,
    2,  3,  6,  7,  18, 19, 22, 23,
    8,  9,  12, 13, 24, 25, 28, 29,
    10, 11, 14, 15, 26, 27, 30, 31,
    32, 33, 36, 37, 48, 49, 52, 53,
    34, 35, 38, 39, 50, 51, 54, 55,
    40, 41, 44, 45, 56, 57, 60, 61,
    42, 43, 46, 47, 58, 59, 62, 63,
};

// clang-format on

// previous QG in decoding order
PointAddr HevcFrame::getPrevQGAddr(int xAddr, int yAddr) {
    auto getZscanAddr = [](int log2Size, int rsAddr) -> int {
        switch (log2Size) {
            case 1:
                return gRasterToZscan2x2[rsAddr];
            case 2:
                return gRasterToZscan4x4[rsAddr];
            case 3:
                return gRasterToZscan8x8[rsAddr];
            default:
                return -1;
        }
    };

    auto getRscanAddr = [](int log2Size, int zsAddr) -> int {
        switch (log2Size) {
            case 1:
                return gZscanToRaster2x2[zsAddr];
            case 2:
                return gZscanToRaster4x4[zsAddr];
            case 3:
                return gZscanToRaster8x8[zsAddr];
            default:
                return -1;
        }
    };

    auto ctu = getCtu(xAddr, yAddr);
    int x = (xAddr - ctu->getX()) >> mPps->log2MinCuQpDeltaSize;
    int y = (yAddr - ctu->getY()) >> mPps->log2MinCuQpDeltaSize;
    int rsAddr = (y << (mSps->log2CtbSize - mPps->log2MinCuQpDeltaSize)) + x;

    if (rsAddr == 0) {
        // get previous ctu
        if (ctu->getRsAddr() == 0) return {-1, -1};
        int tsAddr = mPps->ctbAddrRsToTs[ctu->getRsAddr()];
        auto prevCtu = mCtus[mPps->ctbAddrTsToRs[tsAddr - 1]];
        // get last QG
        x = prevCtu->getX() + (1 << mSps->log2CtbSize) - (1 << mPps->log2MinCuQpDeltaSize);
        y = prevCtu->getY() + (1 << mSps->log2CtbSize) - (1 << mPps->log2MinCuQpDeltaSize);
    } else {
        int zsAddr = getZscanAddr(mSps->log2CtbSize - mPps->log2MinCuQpDeltaSize, rsAddr);
        rsAddr = getRscanAddr(mSps->log2CtbSize - mPps->log2MinCuQpDeltaSize, zsAddr - 1);
        x = rsAddr & ((1 << (mSps->log2CtbSize - mPps->log2MinCuQpDeltaSize)) - 1);
        y = rsAddr >> (mSps->log2CtbSize - mPps->log2MinCuQpDeltaSize);
        x = ctu->getX() + (x << mPps->log2MinCuQpDeltaSize);
        y = ctu->getY() + (y << mPps->log2MinCuQpDeltaSize);
    }

    return {x, y};
}

void HevcFrame::setRefPicList(vector<std::shared_ptr<HevcFrame>> list[2]) {
    for (int i = 0; i < 2; i++)
        for (auto frame : list[i]) {
            mRefPicList[i].push_back(frame);
            mRefPocList[i].push_back(frame->getPoc());
        }
}
