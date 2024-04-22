#ifndef BIT_READER_H
#define BIT_READER_H

#include <cstdint>

class BitReader {
private:
    const uint8_t *mpData;
    uint32_t mSize;
    uint32_t mReservoir;
    uint32_t mNumBitsLeft;
    uint32_t mReadCount = 0;

    bool fillReservoir();
    bool getBitsGraceful(uint32_t n, uint32_t *out);

public:
    BitReader(const uint8_t *data, uint32_t size);
    BitReader(const BitReader &) = delete;
    BitReader &operator=(const BitReader &) = delete;
    virtual ~BitReader();

    bool skipBits(uint32_t n);
    uint32_t numBytesLeft() const;
    uint32_t numBitsUntilByteAligned() const;
    uint32_t numBytesRead() const;
    const uint8_t *data() const;

    uint32_t getUInt(uint32_t n);
    int32_t getSInt(uint32_t n);
    uint32_t getFlag();
    uint32_t getUVlc();
    int32_t getSVlc();
};

#endif
