#ifndef HEVC_PU_H
#define HEVC_PU_H

#include "Typedef.h"

#include <cstdint>

class HevcPu {
private:
    int mXAddr, mYAddr;
    int mWidth, mHeight;

    uint8_t mMergeFlag;

    MvField mMvField;

public:
    HevcPu(int xAddr, int yAddr, int width, int height);
    int getX() { return mXAddr; }
    int getY() { return mYAddr; }
    int getWidth() { return mWidth; }
    int getHeight() { return mHeight; }
    void setMergeFlag(uint8_t flag) { mMergeFlag = flag; }
    uint8_t getMergeFlag() { return mMergeFlag; }
    void getMvField(MvField &mvField) { mvField = mMvField; }
    void setMvField(MvField &mvField) { mMvField = mvField; }
};

#endif
