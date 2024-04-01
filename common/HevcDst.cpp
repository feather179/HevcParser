#include "HevcDst.h"

#include <algorithm>

using std::vector;
using std::clamp;

template <typename ValueType>
inline ValueType leftShift(const ValueType value, const int shift) {
    return (shift >= 0) ? (value << shift) : (value >> -shift);
}
template <typename ValueType>
inline ValueType rightShift(const ValueType value, const int shift) {
    return (shift >= 0) ? (value >> shift) : (value << -shift);
}

// clang-format off
#define DEFINE_DST4x4_MATRIX(a,b,c,d) \
    {  a,  b,  c,  d }, \
    {  c,  c,  0, -c }, \
    {  d, -a, -c,  b }, \
    {  b, -d,  c, -a }, \

const static int gDstMatrix4[4][4] = {
    DEFINE_DST4x4_MATRIX(29, 55, 74, 84)
};

// clang-format on

void fastForwardDst(vector<vector<int>> &src, vector<vector<int>> &dst, int shift) {
    int i;
    int c[4];
    int rndFactor = (shift > 0) ? (1 << (shift - 1)) : 0;
    for (i = 0; i < 4; i++) {
        // Intermediate Variables
        c[0] = src[i][0];
        c[1] = src[i][1];
        c[2] = src[i][2];
        c[3] = src[i][3];

        for (int row = 0; row < 4; row++) {
            int result = 0;
            for (int column = 0; column < 4; column++) {
                result += c[column] * gDstMatrix4[row][column];
            }

            dst[row][i] = rightShift((result + rndFactor), shift);
        }
    }
}

void fastInverseDst(
    vector<vector<int>> &src, vector<vector<int>> &dst, int shift, const int outputMinimum, const int outputMaximum) {
    int i;
    int c[4];
    int rndFactor = (shift > 0) ? (1 << (shift - 1)) : 0;
    for (i = 0; i < 4; i++) {
        // Intermediate Variables
        c[0] = src[0][i];
        c[1] = src[1][i];
        c[2] = src[2][i];
        c[3] = src[3][i];

        for (int column = 0; column < 4; column++) {
            int result = 0;
            for (int row = 0; row < 4; row++) {
                result += c[row] * gDstMatrix4[row][column];
            }

            dst[i][column] = clamp(rightShift((result + rndFactor), shift), outputMinimum, outputMaximum);
        }
    }
}
