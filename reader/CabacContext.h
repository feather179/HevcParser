#ifndef CABAC_CONTEXT_H
#define CABAC_CONTEXT_H

#include <cstdint>

class CabacContext {
    friend class CabacReader;

private:
    uint8_t mStateIdx = 0;
    uint8_t mMps = 0;

public:
    CabacContext();
    ~CabacContext();

    void init(int qp, int initValue);
    void updateLps();
    void updateMps();

    uint8_t getStateIdx() { return mStateIdx; }
    uint8_t getMps() { return mMps; }
    uint8_t getLps() { return 1 - mMps; }
};

#endif
