#ifndef HEVC_DECODER_H
#define HEVC_DECODER_H

#include <cstdint>
#include <cstdio>
#include <vector>
#include <memory>

class ParameterSet;
class HevcFrame;
class HevcSliceDecoder;
class HevcFilter;

class HevcDecoder {
private:
    FILE *mOutputFile = nullptr;

    std::vector<uint32_t> mEmulationPreventionByteLocation;
    std::shared_ptr<ParameterSet> mParamSet;
    bool mFirstSliceInPicture = true;
    std::vector<std::shared_ptr<HevcFrame>> mFrames;
    std::shared_ptr<HevcFrame> mCurrFrame;
    std::shared_ptr<HevcSliceDecoder> mSliceDecoder;
    std::shared_ptr<HevcFilter> mFilter;

    int mPrevTid0Poc = 0;

    void extractRbsp(std::vector<uint8_t> &nalu);
    void decodePacket(std::vector<uint8_t> &buffer);
    void decodeNalu(std::vector<uint8_t> &nalu);

    std::shared_ptr<HevcFrame> getNewFrame(uint8_t tid);

public:
    HevcDecoder();
    // ~HevcDecoder();
    bool decode(const char *inputFilePath, const char *outputFilePath);
};

#endif
