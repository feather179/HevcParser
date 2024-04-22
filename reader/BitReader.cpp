#include "BitReader.h"

#include <cmath>

BitReader::BitReader(const uint8_t *data, uint32_t size) : mpData(data), mSize(size), mReservoir(0), mNumBitsLeft(0) {}

BitReader::~BitReader() {}

bool BitReader::fillReservoir() {
    if (mSize == 0) return false;

    mReservoir = 0;
    int i = 0;
    for (; i < 4 && mSize > 0; i++) {
        mReservoir = (mReservoir << 8) | *mpData;
        mpData++;
        mSize--;
    }

    mNumBitsLeft = 8 * i;
    mReservoir <<= (32 - 8 * i);
    mReadCount += i;
    return true;
}

bool BitReader::getBitsGraceful(uint32_t n, uint32_t *out) {
    if (n > 32) return false;

    uint32_t ret = 0;
    while (n > 0) {
        if (mNumBitsLeft == 0)
            if (!fillReservoir()) return false;

        uint32_t m = n;
        if (m > mNumBitsLeft) m = mNumBitsLeft;
        ret = (ret << m) | (mReservoir >> (32 - m));
        mReservoir <<= m;
        n -= m;
        mNumBitsLeft -= m;
    }

    *out = ret;
    return true;
}

bool BitReader::skipBits(uint32_t n) {
    uint32_t ret;
    while (n > 32) {
        if (!getBitsGraceful(32, &ret)) return false;
        n -= 32;
    }

    if (n > 0) return getBitsGraceful(n, &ret);

    return true;
}

uint32_t BitReader::numBytesLeft() const {
    return mSize + ((mNumBitsLeft + 7) >> 3);
}

uint32_t BitReader::numBitsUntilByteAligned() const {
    return mNumBitsLeft & 0x07;
}

uint32_t BitReader::numBytesRead() const {
    return mReadCount - ((mNumBitsLeft + 7) >> 3);
}

const uint8_t *BitReader::data() const {
    return mpData - ((mNumBitsLeft + 7) >> 3);
}

uint32_t BitReader::getUInt(uint32_t n) {
    uint32_t ret = 0;
    getBitsGraceful(n, &ret);
    return ret;
}

int32_t BitReader::getSInt(uint32_t n) {
    uint32_t value = 0;
    int32_t ret = 0;

    getBitsGraceful(n, &value);
    if (value & (1 << (n - 1))) {
        value &= ((1 << (n - 1)) - 1);
        ret = ~value + 1;
    } else {
        ret = value;
    }

    return ret;
}

uint32_t BitReader::getFlag() {
    return getUInt(1);
}

uint32_t BitReader::getUVlc() {
    uint32_t ret = 0;
    uint32_t numZero = 0;

    while (!getUInt(1)) numZero++;

    ret = getUInt(numZero);
    ret += ((1 << numZero) - 1);

    return ret;
}

int32_t BitReader::getSVlc() {
    uint32_t value = 0;
    int32_t ret = 0;

    value = getUVlc();
    ret = (int32_t)std::ceil(value / 2.0f);
    if (!(value & 0x01)) ret = -ret;

    return ret;
}
