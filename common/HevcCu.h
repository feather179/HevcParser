#ifndef HEVC_CU_H
#define HEVC_CU_H

#include "Typedef.h"

#include <vector>
#include <memory>

class HevcTu;
class HevcPu;

class HevcCu {
    friend class HevcCuDecoder;

private:
    int mXAddr, mYAddr;
    int mIndexInCtu;
    int mLog2Size;
    int mDepth;
    int mChromaArrayType;
    bool mCuSkipFlag;
    bool mPcmFlag = false;
    bool mCuTransquantBypassFlag;
    PredMode mPredMode;
    PartMode mPartMode;
    int mQpY, mQpYPred;
    int mQpCb, mQpCr;

    std::vector<IntraPredMode> mIntraPredModeY;
    std::vector<IntraPredMode> mIntraPredModeC;
    std::shared_ptr<HevcTu> mTu;
    std::vector<std::shared_ptr<HevcPu>> mPus;

    std::vector<std::vector<uint16_t>> mPredY;
    std::vector<std::vector<uint16_t>> mPredCb;
    std::vector<std::vector<uint16_t>> mPredCr;
    std::vector<std::vector<int>> mResidualY;
    std::vector<std::vector<int>> mResidualCb;
    std::vector<std::vector<int>> mResidualCr;

public:
    HevcCu(int xAddr, int yAddr, int indexInCtu, int log2Size, int depth);
    int getX() { return mXAddr; }
    int getY() { return mYAddr; }
    int getIndexInCtu() { return mIndexInCtu; }
    int getLog2Size() { return mLog2Size; }
    int getDepth() { return mDepth; }
    void setChromaArrayType(int type) { mChromaArrayType = type; }
    void setCuSkipFlag(bool flag) { mCuSkipFlag = flag; }
    bool getCuSkipFlag() { return mCuSkipFlag; }
    void setPcmFlag(bool flag) { mPcmFlag = flag; }
    bool getPcmFlag() { return mPcmFlag; }
    void setCuTransquantBypassFlag(bool flag) { mCuTransquantBypassFlag = flag; }
    bool getCuTransquantBypassFlag() { return mCuTransquantBypassFlag; }
    void setPredMode(PredMode mode) { mPredMode = mode; }
    PredMode getPredMode() { return mPredMode; }
    void setPartMode(PartMode mode) { mPartMode = mode; }
    PartMode getPartMode() { return mPartMode; }
    void setQpY(int qpY) { mQpY = qpY; }
    int getQpY() { return mQpY; }
    void setQpYPred(int qpYPred) { mQpYPred = qpYPred; }
    int getQpYPred() { return mQpYPred; }
    void setQpCb(int qpCb) { mQpCb = qpCb; }
    int getQpCb() { return mQpCb; }
    void setQpCr(int qpCr) { mQpCr = qpCr; }
    int getQpCr() { return mQpCr; }
    void setIntraPredModeY(int blkIdx, IntraPredMode mode);
    IntraPredMode getIntraPredModeY(int x, int y);
    IntraPredMode getIntraPredModeY(int blkIdx);
    void setIntraPredModeC(int blkIdx, IntraPredMode mode);
    IntraPredMode getIntraPredModeC(int x, int y);
    IntraPredMode getIntraPredModeC(int blkIdx);

    void addTu(std::shared_ptr<HevcTu> tu);
    std::shared_ptr<HevcTu> getTu(int x, int y, int depth = -1);
    void copyResidual(std::shared_ptr<HevcTu> tu);

    void addPu(std::shared_ptr<HevcPu> pu) { mPus.push_back(pu); }
    int getPuCount() { return (int)mPus.size(); }
    std::shared_ptr<HevcPu> getPu(int xAddr, int yAddr);
    std::shared_ptr<HevcPu> getPu(int index) { return mPus[index]; }
    uint8_t getMergeFlag(int xAddr, int yAddr);
};

#endif
