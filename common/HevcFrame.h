#ifndef HEVC_FRAME_H
#define HEVC_FRAME_H

#include "Typedef.h"

#include <cstdint>
#include <vector>
#include <memory>

class ParameterSet;
class HevcSlice;
class HevcCtu;
class HevcCu;
class HevcPu;

class HevcFrame : public std::enable_shared_from_this<HevcFrame> {
    friend class HevcDecoder;
    friend class HevcCuDecoder;
    friend class HevcFilter;

private:
    int mWidth, mHeight;
    int mPoc, mTid;
    bool mUsedForReference = false;
    std::shared_ptr<HevcSps> mSps;
    std::shared_ptr<HevcPps> mPps;

    std::vector<std::shared_ptr<HevcSlice>> mSlices;
    std::vector<int> mSliceStartCtbRsAddr;
    std::vector<std::shared_ptr<HevcCtu>> mCtus;
    std::vector<std::vector<uint8_t>> mCtDepth;

    std::vector<std::shared_ptr<HevcFrame>> mRefPicList[2];
    std::vector<int> mRefPocList[2];

    std::vector<std::vector<uint16_t>> mReconY;
    std::vector<std::vector<uint16_t>> mReconCb;
    std::vector<std::vector<uint16_t>> mReconCr;

    std::vector<std::vector<uint16_t>> mSaoY;
    std::vector<std::vector<uint16_t>> mSaoCb;
    std::vector<std::vector<uint16_t>> mSaoCr;

public:
    void setPoc(int poc) { mPoc = poc; }
    int getPoc() { return mPoc; }
    void setTid(int tid) { mTid = tid; }
    int getTid() { return mTid; }
    int getCtuCount() { return (int)mCtus.size(); }
    void setUsedForReference(bool usedForReference) { mUsedForReference = usedForReference; }
    bool getUsedForReference() { return mUsedForReference; }
    void init(std::shared_ptr<ParameterSet> paramSet);
    void addSlice(int sliceStartCtbRsAddr) { mSliceStartCtbRsAddr.push_back(sliceStartCtbRsAddr); }
    std::shared_ptr<HevcSps> getSps() { return mSps; }
    std::shared_ptr<HevcPps> getPps() { return mPps; }
    std::shared_ptr<HevcCtu> getCtu(int rsAddr);
    std::shared_ptr<HevcCtu> getCtu(int xAddr, int yAddr);
    std::shared_ptr<HevcCu> getCu(int xAddr, int yAddr);
    std::shared_ptr<HevcPu> getPu(int xAddr, int yAddr);
    void setCtDepth(int xAddr, int yAddr, int log2CbSize, uint8_t depth);
    uint8_t getCtDepth(int xAddr, int yAddr);
    PointAddr getSliceAddr(int xAddr, int yAddr);
    PointAddr getTileAddr(int xAddr, int yAddr);
    PointAddr getQGAddr(int xAddr, int yAddr);      // QuantizationGroup
    PointAddr getPrevQGAddr(int xAddr, int yAddr);  // QuantizationGroup
    void setRefPicList(std::vector<std::shared_ptr<HevcFrame>> list[2]);
};

#endif
