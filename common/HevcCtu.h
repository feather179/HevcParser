#ifndef HEVC_CTU_H
#define HEVC_CTU_H

#include "Typedef.h"

#include <cstdint>
#include <memory>
#include <vector>

class CabacReader;
class HevcFrame;
class HevcSlice;
class HevcCu;

class HevcCtu {
private:
    int mRsAddr;
    int mXAddr, mYAddr;
    int mLog2Size;
    int mLog2MinCuQpDeltaSize;
    int mSliceRsAddr = -1;
    int mTileId = -1;
    SaoParam mSaoParam[3];
    std::weak_ptr<HevcFrame> mFrame;
    std::shared_ptr<HevcSlice> mSlice;
    std::vector<std::shared_ptr<HevcCu>> mCus;
    std::vector<std::vector<int>> mCuIndex;

    uint8_t mIsCuQpDeltaCoded;
    int mCuQpDeltaVal;

public:
    HevcCtu(int xAddr, int yAddr, int rsAddr, std::shared_ptr<HevcFrame> frame);
    int getX() { return mXAddr; }
    int getY() { return mYAddr; }
    int getLog2Size() { return mLog2Size; }
    int getCuCount() { return mCus.size(); }
    int getRsAddr() { return mRsAddr; }
    void setSliceRsAddr(int addr) { mSliceRsAddr = addr; }
    int getSliceRsAddr() { return mSliceRsAddr; }
    void setTileId(int tileId) { mTileId = tileId; }
    int getTileId() { return mTileId; }
    void setSaoParam(SaoParam saoParam[3]);
    SaoParam getSaoParam(int cIdx) { return mSaoParam[cIdx]; }
    std::shared_ptr<HevcFrame> getFrame() { return mFrame.lock(); }
    void setSlice(std::shared_ptr<HevcSlice> slice) { mSlice = slice; }
    std::shared_ptr<HevcSlice> getSlice() { return mSlice; }
    void addCu(std::shared_ptr<HevcCu> cu);
    std::shared_ptr<HevcCu> getCu(int xAddr, int yAddr) const;
    std::shared_ptr<HevcCu> getCu(int index) const;
    void setIsCuQpDeltaCoded(uint8_t flag) { mIsCuQpDeltaCoded = flag; }
    uint8_t getIsCuQpDeltaCoded() { return mIsCuQpDeltaCoded; }
    void setCuQpDeltaVal(int value) { mCuQpDeltaVal = value; }
    int getCuQpDeltaVal() { return mCuQpDeltaVal; }
};

#endif
