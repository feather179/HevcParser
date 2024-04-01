#include "HevcCu.h"
#include "HevcTu.h"
#include "HevcPu.h"

#include <assert.h>

using std::vector;

HevcCu::HevcCu(int xAddr, int yAddr, int indexInCtu, int log2Size, int depth)
    : mXAddr(xAddr), mYAddr(yAddr), mIndexInCtu(indexInCtu), mLog2Size(log2Size), mDepth(depth) {
    mPredY.resize(1 << log2Size, vector<uint16_t>(1 << log2Size, 0));
    mPredCb.resize(1 << log2Size, vector<uint16_t>(1 << log2Size, 0));
    mPredCr.resize(1 << log2Size, vector<uint16_t>(1 << log2Size, 0));
    mResidualY.resize(1 << log2Size, vector<int>(1 << log2Size, 0));
    mResidualCb.resize(1 << log2Size, vector<int>(1 << log2Size, 0));
    mResidualCr.resize(1 << log2Size, vector<int>(1 << log2Size, 0));
}

void HevcCu::setIntraPredModeY(int blkIdx, IntraPredMode mode) {
    assert(mPredMode == PRED_MODE_INTRA);

    int puNum = 1;
    if (mPartMode == PART_MODE_NxN) puNum = 4;

    assert(blkIdx < puNum);

    if (mIntraPredModeY.size() != puNum) mIntraPredModeY.resize(puNum);
    mIntraPredModeY[blkIdx] = mode;
}

IntraPredMode HevcCu::getIntraPredModeY(int x, int y) {
    assert(mPredMode == PRED_MODE_INTRA);
    assert((x >= mXAddr) && (x < (mXAddr + (1 << mLog2Size))));
    assert((y >= mYAddr) && (y < (mYAddr + (1 << mLog2Size))));

    if (mPartMode == PART_MODE_NxN) {
        int blkIdx = ((y - mYAddr) >> (mLog2Size - 1) << 1) + ((x - mXAddr) >> (mLog2Size - 1));
        return mIntraPredModeY[blkIdx];
    } else
        return mIntraPredModeY[0];
}

IntraPredMode HevcCu::getIntraPredModeY(int blkIdx) {
    return mIntraPredModeY[blkIdx];
}

void HevcCu::setIntraPredModeC(int blkIdx, IntraPredMode mode) {
    assert(mPredMode == PRED_MODE_INTRA);

    int puNum = 1;
    if (mPartMode == PART_MODE_NxN && mChromaArrayType == 3) puNum = 4;

    assert(blkIdx < puNum);

    if (mIntraPredModeC.size() != puNum) mIntraPredModeC.resize(puNum);
    mIntraPredModeC[blkIdx] = mode;
}

IntraPredMode HevcCu::getIntraPredModeC(int x, int y) {
    assert(mPredMode == PRED_MODE_INTRA);
    assert((x >= mXAddr) && (x < (mXAddr + (1 << mLog2Size))));
    assert((y >= mYAddr) && (y < (mYAddr + (1 << mLog2Size))));

    if (mPartMode == PART_MODE_NxN && mChromaArrayType == 3) {
        int blkIdx = ((y - mYAddr) >> (mLog2Size - 1) << 1) + ((x - mXAddr) >> (mLog2Size - 1));
        return mIntraPredModeC[blkIdx];
    } else
        return mIntraPredModeC[0];
}

IntraPredMode HevcCu::getIntraPredModeC(int blkIdx) {
    return mIntraPredModeC[blkIdx];
}

void HevcCu::addTu(std::shared_ptr<HevcTu> tu) {
    if (tu->getDepth() == 0)
        mTu = tu;
    else
        mTu->addTu(tu);
}

std::shared_ptr<HevcTu> HevcCu::getTu(int x, int y, int depth) {
    if (mTu) return mTu->getTu(x, y, depth);
    return nullptr;
}

void HevcCu::copyResidual(std::shared_ptr<HevcTu> tu) {
    int rowBaseY = tu->getY() - mYAddr;
    int colBaseY = tu->getX() - mXAddr;

    int rowBaseC = rowBaseY >> 1;
    int colBaseC = colBaseY >> 1;

    if (tu->getLog2Size() == 2) {
        rowBaseC = (rowBaseY - 4) >> 1;
        colBaseC = (colBaseY - 4) >> 1;
    }

    for (int row = 0; row < tu->mResidualY.size(); row++)
        for (int col = 0; col < tu->mResidualY[row].size(); col++)
            mResidualY[row + rowBaseY][col + colBaseY] = tu->mResidualY[row][col];

    for (int row = 0; row < tu->mResidualCb.size(); row++)
        for (int col = 0; col < tu->mResidualCb[row].size(); col++)
            mResidualCb[row + rowBaseC][col + colBaseC] = tu->mResidualCb[row][col];

    for (int row = 0; row < tu->mResidualCr.size(); row++)
        for (int col = 0; col < tu->mResidualCr[row].size(); col++)
            mResidualCr[row + rowBaseC][col + colBaseC] = tu->mResidualCr[row][col];
}

std::shared_ptr<HevcPu> HevcCu::getPu(int xAddr, int yAddr) {
    for (auto pu : mPus) {
        if (xAddr >= pu->getX() && xAddr < (pu->getX() + pu->getWidth()) && yAddr >= pu->getY() &&
            yAddr < (pu->getY() + pu->getHeight()))
            return pu;
    }

    return nullptr;
}

uint8_t HevcCu::getMergeFlag(int xAddr, int yAddr) {
    auto pu = getPu(xAddr, yAddr);
    if (pu) return pu->getMergeFlag();

    return 0;
}
