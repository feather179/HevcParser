#ifndef HEVC_TU_H
#define HEVC_TU_H

#include <cstdint>
#include <vector>
#include <memory>

class HevcTu : public std::enable_shared_from_this<HevcTu> {
    friend class HevcCuDecoder;
    friend class HevcCu;

private:
    int mXAddr, mYAddr;
    int mLog2Size;
    int mDepth;
    bool mSplitTransformFlag;
    std::vector<std::shared_ptr<HevcTu>> mSubTus;

    uint8_t mCbfLuma;
    uint8_t mCbfCb[2], mCbfCr[2];

    bool mTransformSkipFlag[3];
    bool mExplicitRdpcmFlag[3];
    bool mExplicitRdpcmDirFlag[3];

    uint8_t mCodedSubBlockFlag[8][8] = {{0}};

    std::vector<std::vector<int>> mTransCoeffY;
    std::vector<std::vector<int>> mTransCoeffCb;
    std::vector<std::vector<int>> mTransCoeffCr;
    std::vector<std::vector<int>> mResidualY;
    std::vector<std::vector<int>> mResidualCb;
    std::vector<std::vector<int>> mResidualCr;

public:
    HevcTu(int xAddr, int yAddr, int log2Size, int depth);
    int getX() { return mXAddr; }
    int getY() { return mYAddr; }
    int getLog2Size() { return mLog2Size; }
    int getDepth() { return mDepth; }
    void setSplitTransformFlag(bool flag) { mSplitTransformFlag = flag; }
    bool getSplitTransformFlag() { return mSplitTransformFlag; }
    void setCbfLuma(uint8_t cbfLuma) { mCbfLuma = cbfLuma; }
    uint8_t getCbfLuma() { return mCbfLuma; }
    void setCbfCb(const uint8_t cbfCb[2]) {
        mCbfCb[0] = cbfCb[0];
        mCbfCb[1] = cbfCb[1];
    }
    uint8_t getCbfCb(int idx) { return mCbfCb[idx]; }
    void setCbfCr(const uint8_t cbfCr[2]) {
        mCbfCr[0] = cbfCr[0];
        mCbfCr[1] = cbfCr[1];
    }
    uint8_t getCbfCr(int idx) { return mCbfCr[idx]; }
    void setTransformSkipFlag(int cIdx, bool flag) { mTransformSkipFlag[cIdx] = flag; }
    bool getTransformSkipFlag(int cIdx) { return mTransformSkipFlag[cIdx]; }
    void setExplicitRdpcmFlag(int cIdx, bool flag) { mExplicitRdpcmFlag[cIdx] = flag; }
    bool getExplicitRdpcmFlag(int cIdx) { return mExplicitRdpcmFlag[cIdx]; }
    void setExplicitRdpcmDirFlag(int cIdx, bool flag) { mExplicitRdpcmDirFlag[cIdx] = flag; }
    bool getExplicitRdpcmDirFlag(int cIdx) { return mExplicitRdpcmDirFlag[cIdx]; }
    void setCodedSubBlockFlag(int xAddr, int yAddr, uint8_t flag) {
        mCodedSubBlockFlag[yAddr][xAddr] = flag;
    }
    uint8_t getCodedSubBlockFlag(int xAddr, int yAddr) { return mCodedSubBlockFlag[yAddr][xAddr]; }
    void resetCodedSubBlockFlag() {
        for (int y = 0; y < 8; y++)
            for (int x = 0; x < 8; x++) mCodedSubBlockFlag[y][x] = 0;
    }
    void setComponentLog2Size(int cIdx, int log2Size);

    void addTu(std::shared_ptr<HevcTu> tu);
    std::shared_ptr<HevcTu> getTu(int xAddr, int yAddr, int depth);
};

#endif
