#ifndef CAVLC_READER_H
#define CAVLC_READER_H

#include "BitReader.h"
#include "Typedef.h"

#include <cstdint>
#include <memory>
#include <vector>

// parse SEI/SPS/PPS/VPS/Slice Header

class HevcFrame;
class ParameterSet;

class CavlcReader : public BitReader {
private:
    NalUnitType mNaluType;
    uint8_t mTid;

    void parseVps(std::shared_ptr<ParameterSet> ps);
    void parseSps(std::shared_ptr<ParameterSet> ps);
    void parsePps(std::shared_ptr<ParameterSet> ps);

    void parseSei(std::shared_ptr<ParameterSet> ps);
    void parseSliceHeader(std::shared_ptr<ParameterSet> ps,
                          int prevTid0Poc,
                          std::vector<std::shared_ptr<HevcFrame>> &dpb);

    void parsePtl(HevcProfileTierLevel *ptl, bool profilePresentFlag, uint8_t maxNumSubLayersMinus1);
    void parseStRps(HevcStRefPicSet *stRps, int stRpsIdx, int numStRps);

public:
    CavlcReader(const std::vector<uint8_t> &rbsp);
    ~CavlcReader();

    void parse(std::shared_ptr<ParameterSet> ps, int prevTid0Poc, std::vector<std::shared_ptr<HevcFrame>> &dpb);
    int getNaluType() { return mNaluType; }
    uint8_t getTid() { return mTid; }
};

#endif
