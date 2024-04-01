#include "HevcCtu.h"
#include "HevcFrame.h"
#include "HevcCu.h"

using std::vector;

HevcCtu::HevcCtu(int xAddr, int yAddr, int rsAddr, std::shared_ptr<HevcFrame> frame)
    : mXAddr(xAddr), mYAddr(yAddr), mRsAddr(rsAddr), mFrame(frame) {
    auto sps = frame->getSps();
    auto pps = frame->getPps();
    mLog2Size = sps->log2CtbSize;
    mLog2MinCuQpDeltaSize = pps->log2MinCuQpDeltaSize;

    int cuIndexCount = 1 << (sps->log2CtbSize - sps->log2MinCtbSize);
    mCuIndex.resize(cuIndexCount, vector<int>(cuIndexCount, 0));
    mCus.reserve(cuIndexCount * cuIndexCount);
}

void HevcCtu::setSaoParam(SaoParam saoParam[3]) {
    for (int i = 0; i < 3; i++) {
        if (saoParam[i].mode == SAO_MODE_MERGE_LEFT)
            mSaoParam[i] = mFrame.lock()->getCtu(mXAddr - 1, mYAddr)->getSaoParam(i);
        else if (saoParam[i].mode == SAO_MODE_MERGE_ABOVE)
            mSaoParam[i] = mFrame.lock()->getCtu(mXAddr, mYAddr - 1)->getSaoParam(i);
        else
            mSaoParam[i] = saoParam[i];
    }
}

void HevcCtu::addCu(std::shared_ptr<HevcCu> cu) {
    if ((cu->getX() < mXAddr) || ((cu->getX() + (1 << cu->getLog2Size())) > (mXAddr + (1 << mLog2Size)))) return;
    if ((cu->getY() < mYAddr) || ((cu->getY() + (1 << cu->getLog2Size())) > (mYAddr + (1 << mLog2Size)))) return;

    auto sps = mFrame.lock()->getSps();
    int x = (cu->getX() - mXAddr) >> sps->log2MinCtbSize;
    int y = (cu->getY() - mYAddr) >> sps->log2MinCtbSize;
    int count = 1 << (cu->getLog2Size() - sps->log2MinCtbSize);

    for (int i = 0; i < count; i++)
        for (int j = 0; j < count; j++) mCuIndex[y + i][x + j] = cu->getIndexInCtu();

    if (mCus.size() <= cu->getIndexInCtu()) mCus.resize(cu->getIndexInCtu() + 1);
    mCus[cu->getIndexInCtu()] = cu;
}

std::shared_ptr<HevcCu> HevcCtu::getCu(int xAddr, int yAddr) const {
    auto sps = mFrame.lock()->getSps();
    int x = (xAddr - mXAddr) >> sps->log2MinCtbSize;
    int y = (yAddr - mYAddr) >> sps->log2MinCtbSize;
    int indexInCtu = mCuIndex[y][x];
    if (indexInCtu >= mCus.size()) return nullptr;
    return mCus[indexInCtu];
}

std::shared_ptr<HevcCu> HevcCtu::getCu(int index) const {
    return mCus[index];
}
