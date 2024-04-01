#ifndef HEVC_SLICE_DECODER_H
#define HEVC_SLICE_DECODER_H

#include "Typedef.h"

#include <cstdint>
#include <memory>
#include <vector>

class HevcFrame;

class HevcSliceDecoder {
private:
    std::shared_ptr<HevcSps> mSps;
    std::shared_ptr<HevcPps> mPps;
    // std::shared_ptr<HevcSliceHeader> mSliceHeader;

    uint32_t getSubStreamIdx(uint32_t ctuRsAddrInSlice);

public:
    void decodeSlice(std::shared_ptr<HevcFrame> frame,
                     std::shared_ptr<HevcSliceHeader> sliceHeader,
                     std::vector<std::vector<uint8_t>> &subStreams);
};

#endif
