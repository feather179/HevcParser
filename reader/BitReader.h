#ifndef BIT_READER_H
#define BIT_READER_H

#include <cstdint>

class BitReader {
private:
    const uint8_t *mpData;
    size_t mSize;
    uint32_t mReservoir;
    uint32_t mNumBitsLeft;
    uint32_t mReadCount = 0;

    bool fillReservoir();
    bool getBitsGraceful(size_t n, uint32_t *out);

public:
    BitReader(const uint8_t *data, uint32_t size);
    BitReader(const BitReader &) = delete;
    BitReader &operator=(const BitReader &) = delete;
    virtual ~BitReader();

    bool skipBits(size_t n);
    size_t numBytesLeft() const;
    size_t numBitsUntilByteAligned() const;
    size_t numBytesRead() const;
    const uint8_t *data() const;

    uint32_t getUInt(size_t n);
    int32_t getSInt(size_t n);
    uint32_t getFlag();
    uint32_t getUVlc();
    int32_t getSVlc();
};

#endif
