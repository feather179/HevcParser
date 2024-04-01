#ifndef HEVC_SLICE_H
#define HEVC_SLICE_H

#include "Typedef.h"

#include <memory>

class HevcSlice {
private:
    std::shared_ptr<HevcSliceHeader> mSliceHeader;

public:
    HevcSlice(std::shared_ptr<HevcSliceHeader> sliceHeader);
    std::shared_ptr<HevcSliceHeader> getSliceHeader() { return mSliceHeader; }
};

#endif
