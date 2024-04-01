#include "HevcTu.h"

#include <assert.h>

using std::shared_ptr;
using std::vector;

HevcTu::HevcTu(int xAddr, int yAddr, int log2Size, int depth)
    : mXAddr(xAddr), mYAddr(yAddr), mLog2Size(log2Size), mDepth(depth) {}

void HevcTu::setComponentLog2Size(int cIdx, int log2Size) {
    assert(cIdx >= 0 && cIdx <= 2);

    auto &transCoeff = (cIdx == 0) ? mTransCoeffY : ((cIdx == 1) ? mTransCoeffCb : mTransCoeffCr);
    auto &residual = (cIdx == 0) ? mResidualY : ((cIdx == 1) ? mResidualCb : mResidualCr);

    int count = 1 << log2Size;

    transCoeff.resize(count, vector<int>(count, 0));
    residual.resize(count, vector<int>(count, 0));
}

void HevcTu::addTu(shared_ptr<HevcTu> tu) {
    if (!mSplitTransformFlag) return;
    if (mSubTus.size() < 4) mSubTus.resize(4);

    int x = tu->getX() - mXAddr;
    int y = tu->getY() - mYAddr;
    int index = ((y >> (mLog2Size - 1)) << 1) + (x >> (mLog2Size - 1));
    if (tu->getDepth() == (mDepth + 1))
        mSubTus[index] = tu;
    else
        mSubTus[index]->addTu(tu);
}

// depth < 0 means get the Tu with max depth
shared_ptr<HevcTu> HevcTu::getTu(int xAddr, int yAddr, int depth) {
    if (depth < 0) {
        if (!mSplitTransformFlag) return shared_from_this();
    } else {
        if (xAddr == mXAddr && yAddr == mYAddr && depth == mDepth) return shared_from_this();
        if (!mSplitTransformFlag) return nullptr;
    }

    int x = xAddr - mXAddr;
    int y = yAddr - mYAddr;
    int index = ((y >> (mLog2Size - 1)) << 1) + (x >> (mLog2Size - 1));
    return mSubTus[index]->getTu(xAddr, yAddr, depth);
}
