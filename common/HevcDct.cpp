#include "HevcDct.h"

#include <algorithm>

using std::vector;
using std::clamp;

// clang-format off
#define DEFINE_DCT4x4_MATRIX(a,b,c) \
    { a,  a,  a,  a}, \
    { b,  c, -c, -b}, \
    { a, -a, -a,  a}, \
    { c, -b,  b, -c}, \

#define DEFINE_DCT8x8_MATRIX(a,b,c,d,e,f,g) \
    { a,  a,  a,  a,  a,  a,  a,  a}, \
    { d,  e,  f,  g, -g, -f, -e, -d}, \
    { b,  c, -c, -b, -b, -c,  c,  b}, \
    { e, -g, -d, -f,  f,  d,  g, -e}, \
    { a, -a, -a,  a,  a, -a, -a,  a}, \
    { f, -d,  g,  e, -e, -g,  d, -f}, \
    { c, -b,  b, -c, -c,  b, -b,  c}, \
    { g, -f,  e, -d,  d, -e,  f, -g}, \

#define DEFINE_DCT16x16_MATRIX(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) \
    { a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a}, \
    { h,  i,  j,  k,  l,  m,  n,  o, -o, -n, -m, -l, -k, -j, -i, -h}, \
    { d,  e,  f,  g, -g, -f, -e, -d, -d, -e, -f, -g,  g,  f,  e,  d}, \
    { i,  l,  o, -m, -j, -h, -k, -n,  n,  k,  h,  j,  m, -o, -l, -i}, \
    { b,  c, -c, -b, -b, -c,  c,  b,  b,  c, -c, -b, -b, -c,  c,  b}, \
    { j,  o, -k, -i, -n,  l,  h,  m, -m, -h, -l,  n,  i,  k, -o, -j}, \
    { e, -g, -d, -f,  f,  d,  g, -e, -e,  g,  d,  f, -f, -d, -g,  e}, \
    { k, -m, -i,  o,  h,  n, -j, -l,  l,  j, -n, -h, -o,  i,  m, -k}, \
    { a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a}, \
    { l, -j, -n,  h, -o, -i,  m,  k, -k, -m,  i,  o, -h,  n,  j, -l}, \
    { f, -d,  g,  e, -e, -g,  d, -f, -f,  d, -g, -e,  e,  g, -d,  f}, \
    { m, -h,  l,  n, -i,  k,  o, -j,  j, -o, -k,  i, -n, -l,  h, -m}, \
    { c, -b,  b, -c, -c,  b, -b,  c,  c, -b,  b, -c, -c,  b, -b,  c}, \
    { n, -k,  h, -j,  m,  o, -l,  i, -i,  l, -o, -m,  j, -h,  k, -n}, \
    { g, -f,  e, -d,  d, -e,  f, -g, -g,  f, -e,  d, -d,  e, -f,  g}, \
    { o, -n,  m, -l,  k, -j,  i, -h,  h, -i,  j, -k,  l, -m,  n, -o}, \

#define DEFINE_DCT32x32_MATRIX(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,A,B,C,D,E) \
    { a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a,  a}, \
    { p,  q,  r,  s,  t,  u,  v,  w,  x,  y,  z,  A,  B,  C,  D,  E, -E, -D, -C, -B, -A, -z, -y, -x, -w, -v, -u, -t, -s, -r, -q, -p}, \
    { h,  i,  j,  k,  l,  m,  n,  o, -o, -n, -m, -l, -k, -j, -i, -h, -h, -i, -j, -k, -l, -m, -n, -o,  o,  n,  m,  l,  k,  j,  i,  h}, \
    { q,  t,  w,  z,  C, -E, -B, -y, -v, -s, -p, -r, -u, -x, -A, -D,  D,  A,  x,  u,  r,  p,  s,  v,  y,  B,  E, -C, -z, -w, -t, -q}, \
    { d,  e,  f,  g, -g, -f, -e, -d, -d, -e, -f, -g,  g,  f,  e,  d,  d,  e,  f,  g, -g, -f, -e, -d, -d, -e, -f, -g,  g,  f,  e,  d}, \
    { r,  w,  B, -D, -y, -t, -p, -u, -z, -E,  A,  v,  q,  s,  x,  C, -C, -x, -s, -q, -v, -A,  E,  z,  u,  p,  t,  y,  D, -B, -w, -r}, \
    { i,  l,  o, -m, -j, -h, -k, -n,  n,  k,  h,  j,  m, -o, -l, -i, -i, -l, -o,  m,  j,  h,  k,  n, -n, -k, -h, -j, -m,  o,  l,  i}, \
    { s,  z, -D, -w, -p, -v, -C,  A,  t,  r,  y, -E, -x, -q, -u, -B,  B,  u,  q,  x,  E, -y, -r, -t, -A,  C,  v,  p,  w,  D, -z, -s}, \
    { b,  c, -c, -b, -b, -c,  c,  b,  b,  c, -c, -b, -b, -c,  c,  b,  b,  c, -c, -b, -b, -c,  c,  b,  b,  c, -c, -b, -b, -c,  c,  b}, \
    { t,  C, -y, -p, -x,  D,  u,  s,  B, -z, -q, -w,  E,  v,  r,  A, -A, -r, -v, -E,  w,  q,  z, -B, -s, -u, -D,  x,  p,  y, -C, -t}, \
    { j,  o, -k, -i, -n,  l,  h,  m, -m, -h, -l,  n,  i,  k, -o, -j, -j, -o,  k,  i,  n, -l, -h, -m,  m,  h,  l, -n, -i, -k,  o,  j}, \
    { u, -E, -t, -v,  D,  s,  w, -C, -r, -x,  B,  q,  y, -A, -p, -z,  z,  p,  A, -y, -q, -B,  x,  r,  C, -w, -s, -D,  v,  t,  E, -u}, \
    { e, -g, -d, -f,  f,  d,  g, -e, -e,  g,  d,  f, -f, -d, -g,  e,  e, -g, -d, -f,  f,  d,  g, -e, -e,  g,  d,  f, -f, -d, -g,  e}, \
    { v, -B, -p, -C,  u,  w, -A, -q, -D,  t,  x, -z, -r, -E,  s,  y, -y, -s,  E,  r,  z, -x, -t,  D,  q,  A, -w, -u,  C,  p,  B, -v}, \
    { k, -m, -i,  o,  h,  n, -j, -l,  l,  j, -n, -h, -o,  i,  m, -k, -k,  m,  i, -o, -h, -n,  j,  l, -l, -j,  n,  h,  o, -i, -m,  k}, \
    { w, -y, -u,  A,  s, -C, -q,  E,  p,  D, -r, -B,  t,  z, -v, -x,  x,  v, -z, -t,  B,  r, -D, -p, -E,  q,  C, -s, -A,  u,  y, -w}, \
    { a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a,  a, -a, -a,  a}, \
    { x, -v, -z,  t,  B, -r, -D,  p, -E, -q,  C,  s, -A, -u,  y,  w, -w, -y,  u,  A, -s, -C,  q,  E, -p,  D,  r, -B, -t,  z,  v, -x}, \
    { l, -j, -n,  h, -o, -i,  m,  k, -k, -m,  i,  o, -h,  n,  j, -l, -l,  j,  n, -h,  o,  i, -m, -k,  k,  m, -i, -o,  h, -n, -j,  l}, \
    { y, -s, -E,  r, -z, -x,  t,  D, -q,  A,  w, -u, -C,  p, -B, -v,  v,  B, -p,  C,  u, -w, -A,  q, -D, -t,  x,  z, -r,  E,  s, -y}, \
    { f, -d,  g,  e, -e, -g,  d, -f, -f,  d, -g, -e,  e,  g, -d,  f,  f, -d,  g,  e, -e, -g,  d, -f, -f,  d, -g, -e,  e,  g, -d,  f}, \
    { z, -p,  A,  y, -q,  B,  x, -r,  C,  w, -s,  D,  v, -t,  E,  u, -u, -E,  t, -v, -D,  s, -w, -C,  r, -x, -B,  q, -y, -A,  p, -z}, \
    { m, -h,  l,  n, -i,  k,  o, -j,  j, -o, -k,  i, -n, -l,  h, -m, -m,  h, -l, -n,  i, -k, -o,  j, -j,  o,  k, -i,  n,  l, -h,  m}, \
    { A, -r,  v, -E, -w,  q, -z, -B,  s, -u,  D,  x, -p,  y,  C, -t,  t, -C, -y,  p, -x, -D,  u, -s,  B,  z, -q,  w,  E, -v,  r, -A}, \
    { c, -b,  b, -c, -c,  b, -b,  c,  c, -b,  b, -c, -c,  b, -b,  c,  c, -b,  b, -c, -c,  b, -b,  c,  c, -b,  b, -c, -c,  b, -b,  c}, \
    { B, -u,  q, -x,  E,  y, -r,  t, -A, -C,  v, -p,  w, -D, -z,  s, -s,  z,  D, -w,  p, -v,  C,  A, -t,  r, -y, -E,  x, -q,  u, -B}, \
    { n, -k,  h, -j,  m,  o, -l,  i, -i,  l, -o, -m,  j, -h,  k, -n, -n,  k, -h,  j, -m, -o,  l, -i,  i, -l,  o,  m, -j,  h, -k,  n}, \
    { C, -x,  s, -q,  v, -A, -E,  z, -u,  p, -t,  y, -D, -B,  w, -r,  r, -w,  B,  D, -y,  t, -p,  u, -z,  E,  A, -v,  q, -s,  x, -C}, \
    { g, -f,  e, -d,  d, -e,  f, -g, -g,  f, -e,  d, -d,  e, -f,  g,  g, -f,  e, -d,  d, -e,  f, -g, -g,  f, -e,  d, -d,  e, -f,  g}, \
    { D, -A,  x, -u,  r, -p,  s, -v,  y, -B,  E,  C, -z,  w, -t,  q, -q,  t, -w,  z, -C, -E,  B, -y,  v, -s,  p, -r,  u, -x,  A, -D}, \
    { o, -n,  m, -l,  k, -j,  i, -h,  h, -i,  j, -k,  l, -m,  n, -o, -o,  n, -m,  l, -k,  j, -i,  h, -h,  i, -j,  k, -l,  m, -n,  o}, \
    { E, -D,  C, -B,  A, -z,  y, -x,  w, -v,  u, -t,  s, -r,  q, -p,  p, -q,  r, -s,  t, -u,  v, -w,  x, -y,  z, -A,  B, -C,  D, -E}, \

const static int gDctMatrix4[4][4] = {
    DEFINE_DCT4x4_MATRIX(64, 83, 36)
};

const static int gDctMatrix8[8][8] = {
    DEFINE_DCT8x8_MATRIX(64, 83, 36, 89, 75, 50, 18)
};

const static int gDctMatrix16[16][16] = {
    DEFINE_DCT16x16_MATRIX(64, 83, 36, 89, 75, 50, 18, 90, 87, 80, 70, 57, 43, 25, 9)
};

const static int gDctMatrix32[32][32] = {
    DEFINE_DCT32x32_MATRIX(64, 83, 36, 89, 75, 50, 18, 90, 87, 80, 70, 57, 43, 25, 9, 90, 90, 88, 85, 82, 78, 73, 67, 61, 54, 46, 38, 31, 22, 13, 4)
};

// clang-format on

void partialButterfly4(vector<vector<int>> &src, vector<vector<int>> &dst, int shift) {
    int E[2], O[2];
    int add = (shift > 0) ? (1 << (shift - 1)) : 0;

    for (int i = 0; i < 4; i++) {
        /* E and O */
        E[0] = src[i][0] + src[i][3];
        O[0] = src[i][0] - src[i][3];
        E[1] = src[i][1] + src[i][2];
        O[1] = src[i][1] - src[i][2];

        dst[0][i] = (gDctMatrix4[0][0] * E[0] + gDctMatrix4[0][1] * E[1] + add) >> shift;
        dst[2][i] = (gDctMatrix4[2][0] * E[0] + gDctMatrix4[2][1] * E[1] + add) >> shift;
        dst[1][i] = (gDctMatrix4[1][0] * O[0] + gDctMatrix4[1][1] * O[1] + add) >> shift;
        dst[3][i] = (gDctMatrix4[3][0] * O[0] + gDctMatrix4[3][1] * O[1] + add) >> shift;
    }
}

void partialButterflyInverse4(
    vector<vector<int>> &src, vector<vector<int>> &dst, int shift, const int outputMinimum, const int outputMaximum) {
    int E[2], O[2];
    int add = (shift > 0) ? (1 << (shift - 1)) : 0;

    for (int i = 0; i < 4; i++) {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        O[0] = gDctMatrix4[1][0] * src[1][i] + gDctMatrix4[3][0] * src[3][i];
        O[1] = gDctMatrix4[1][1] * src[1][i] + gDctMatrix4[3][1] * src[3][i];
        E[0] = gDctMatrix4[0][0] * src[0][i] + gDctMatrix4[2][0] * src[2][i];
        E[1] = gDctMatrix4[0][1] * src[0][i] + gDctMatrix4[2][1] * src[2][i];

        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial
         * domain vector */
        dst[i][0] = clamp((E[0] + O[0] + add) >> shift, outputMinimum, outputMaximum);
        dst[i][1] = clamp((E[1] + O[1] + add) >> shift, outputMinimum, outputMaximum);
        dst[i][2] = clamp((E[1] - O[1] + add) >> shift, outputMinimum, outputMaximum);
        dst[i][3] = clamp((E[0] - O[0] + add) >> shift, outputMinimum, outputMaximum);
    }
}

void partialButterfly8(vector<vector<int>> &src, vector<vector<int>> &dst, int shift) {
    int k;
    int E[4], O[4];
    int EE[2], EO[2];
    int add = (shift > 0) ? (1 << (shift - 1)) : 0;

    for (int i = 0; i < 8; i++) {
        /* E and O*/
        for (k = 0; k < 4; k++) {
            E[k] = src[i][k] + src[i][7 - k];
            O[k] = src[i][k] - src[i][7 - k];
        }
        /* EE and EO */
        EE[0] = E[0] + E[3];
        EO[0] = E[0] - E[3];
        EE[1] = E[1] + E[2];
        EO[1] = E[1] - E[2];

        dst[0][i] = (gDctMatrix8[0][0] * EE[0] + gDctMatrix8[0][1] * EE[1] + add) >> shift;
        dst[4][i] = (gDctMatrix8[4][0] * EE[0] + gDctMatrix8[4][1] * EE[1] + add) >> shift;
        dst[2][i] = (gDctMatrix8[2][0] * EO[0] + gDctMatrix8[2][1] * EO[1] + add) >> shift;
        dst[6][i] = (gDctMatrix8[6][0] * EO[0] + gDctMatrix8[6][1] * EO[1] + add) >> shift;

        dst[1][i] = (gDctMatrix8[1][0] * O[0] + gDctMatrix8[1][1] * O[1] + gDctMatrix8[1][2] * O[2] +
                     gDctMatrix8[1][3] * O[3] + add) >>
                    shift;
        dst[3][i] = (gDctMatrix8[3][0] * O[0] + gDctMatrix8[3][1] * O[1] + gDctMatrix8[3][2] * O[2] +
                     gDctMatrix8[3][3] * O[3] + add) >>
                    shift;
        dst[5][i] = (gDctMatrix8[5][0] * O[0] + gDctMatrix8[5][1] * O[1] + gDctMatrix8[5][2] * O[2] +
                     gDctMatrix8[5][3] * O[3] + add) >>
                    shift;
        dst[7][i] = (gDctMatrix8[7][0] * O[0] + gDctMatrix8[7][1] * O[1] + gDctMatrix8[7][2] * O[2] +
                     gDctMatrix8[7][3] * O[3] + add) >>
                    shift;
    }
}

void partialButterflyInverse8(
    vector<vector<int>> &src, vector<vector<int>> &dst, int shift, const int outputMinimum, const int outputMaximum) {
    int k;
    int E[4], O[4];
    int EE[2], EO[2];
    int add = (shift > 0) ? (1 << (shift - 1)) : 0;

    for (int i = 0; i < 8; i++) {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        for (k = 0; k < 4; k++) {
            O[k] = gDctMatrix8[1][k] * src[1][i] + gDctMatrix8[3][k] * src[3][i] + gDctMatrix8[5][k] * src[5][i] +
                   gDctMatrix8[7][k] * src[7][i];
        }

        EO[0] = gDctMatrix8[2][0] * src[2][i] + gDctMatrix8[6][0] * src[6][i];
        EO[1] = gDctMatrix8[2][1] * src[2][i] + gDctMatrix8[6][1] * src[6][i];
        EE[0] = gDctMatrix8[0][0] * src[0][i] + gDctMatrix8[4][0] * src[4][i];
        EE[1] = gDctMatrix8[0][1] * src[0][i] + gDctMatrix8[4][1] * src[4][i];

        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial
         * domain vector */
        E[0] = EE[0] + EO[0];
        E[3] = EE[0] - EO[0];
        E[1] = EE[1] + EO[1];
        E[2] = EE[1] - EO[1];
        for (k = 0; k < 4; k++) {
            dst[i][k] = clamp((E[k] + O[k] + add) >> shift, outputMinimum, outputMaximum);
            dst[i][k + 4] = clamp((E[3 - k] - O[3 - k] + add) >> shift, outputMinimum, outputMaximum);
        }
    }
}

void partialButterfly16(vector<vector<int>> &src, vector<vector<int>> &dst, int shift) {
    int k;
    int E[8], O[8];
    int EE[4], EO[4];
    int EEE[2], EEO[2];
    int add = (shift > 0) ? (1 << (shift - 1)) : 0;

    for (int i = 0; i < 16; i++) {
        /* E and O*/
        for (k = 0; k < 8; k++) {
            E[k] = src[i][k] + src[i][15 - k];
            O[k] = src[i][k] - src[i][15 - k];
        }
        /* EE and EO */
        for (k = 0; k < 4; k++) {
            EE[k] = E[k] + E[7 - k];
            EO[k] = E[k] - E[7 - k];
        }
        /* EEE and EEO */
        EEE[0] = EE[0] + EE[3];
        EEO[0] = EE[0] - EE[3];
        EEE[1] = EE[1] + EE[2];
        EEO[1] = EE[1] - EE[2];

        dst[0][i] = (gDctMatrix16[0][0] * EEE[0] + gDctMatrix16[0][1] * EEE[1] + add) >> shift;
        dst[8][i] = (gDctMatrix16[8][0] * EEE[0] + gDctMatrix16[8][1] * EEE[1] + add) >> shift;
        dst[4][i] = (gDctMatrix16[4][0] * EEO[0] + gDctMatrix16[4][1] * EEO[1] + add) >> shift;
        dst[12][i] = (gDctMatrix16[12][0] * EEO[0] + gDctMatrix16[12][1] * EEO[1] + add) >> shift;

        for (k = 2; k < 16; k += 4) {
            dst[k][i] = (gDctMatrix16[k][0] * EO[0] + gDctMatrix16[k][1] * EO[1] + gDctMatrix16[k][2] * EO[2] +
                         gDctMatrix16[k][3] * EO[3] + add) >>
                        shift;
        }

        for (k = 1; k < 16; k += 2) {
            dst[k][i] = (gDctMatrix16[k][0] * O[0] + gDctMatrix16[k][1] * O[1] + gDctMatrix16[k][2] * O[2] +
                         gDctMatrix16[k][3] * O[3] + gDctMatrix16[k][4] * O[4] + gDctMatrix16[k][5] * O[5] +
                         gDctMatrix16[k][6] * O[6] + gDctMatrix16[k][7] * O[7] + add) >>
                        shift;
        }
    }
}

void partialButterflyInverse16(
    vector<vector<int>> &src, vector<vector<int>> &dst, int shift, const int outputMinimum, const int outputMaximum) {
    int k;
    int E[8], O[8];
    int EE[4], EO[4];
    int EEE[2], EEO[2];
    int add = (shift > 0) ? (1 << (shift - 1)) : 0;

    for (int i = 0; i < 16; i++) {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        for (k = 0; k < 8; k++) {
            O[k] = gDctMatrix16[1][k] * src[1][i] + gDctMatrix16[3][k] * src[3][i] + gDctMatrix16[5][k] * src[5][i] +
                   gDctMatrix16[7][k] * src[7][i] + gDctMatrix16[9][k] * src[9][i] + gDctMatrix16[11][k] * src[11][i] +
                   gDctMatrix16[13][k] * src[13][i] + gDctMatrix16[15][k] * src[15][i];
        }
        for (k = 0; k < 4; k++) {
            EO[k] = gDctMatrix16[2][k] * src[2][i] + gDctMatrix16[6][k] * src[6][i] + gDctMatrix16[10][k] * src[10][i] +
                    gDctMatrix16[14][k] * src[14][i];
        }
        EEO[0] = gDctMatrix16[4][0] * src[4][i] + gDctMatrix16[12][0] * src[12][i];
        EEE[0] = gDctMatrix16[0][0] * src[0][i] + gDctMatrix16[8][0] * src[8][i];
        EEO[1] = gDctMatrix16[4][1] * src[4][i] + gDctMatrix16[12][1] * src[12][i];
        EEE[1] = gDctMatrix16[0][1] * src[0][i] + gDctMatrix16[8][1] * src[8][i];

        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial
         * domain vector */
        for (k = 0; k < 2; k++) {
            EE[k] = EEE[k] + EEO[k];
            EE[k + 2] = EEE[1 - k] - EEO[1 - k];
        }
        for (k = 0; k < 4; k++) {
            E[k] = EE[k] + EO[k];
            E[k + 4] = EE[3 - k] - EO[3 - k];
        }
        for (k = 0; k < 8; k++) {
            dst[i][k] = clamp((E[k] + O[k] + add) >> shift, outputMinimum, outputMaximum);
            dst[i][k + 8] = clamp((E[7 - k] - O[7 - k] + add) >> shift, outputMinimum, outputMaximum);
        }
    }
}

void partialButterfly32(vector<vector<int>> &src, vector<vector<int>> &dst, int shift) {
    int k;
    int E[16], O[16];
    int EE[8], EO[8];
    int EEE[4], EEO[4];
    int EEEE[2], EEEO[2];
    int add = (shift > 0) ? (1 << (shift - 1)) : 0;

    for (int i = 0; i < 32; i++) {
        /* E and O*/
        for (k = 0; k < 16; k++) {
            E[k] = src[i][k] + src[i][31 - k];
            O[k] = src[i][k] - src[i][31 - k];
        }
        /* EE and EO */
        for (k = 0; k < 8; k++) {
            EE[k] = E[k] + E[15 - k];
            EO[k] = E[k] - E[15 - k];
        }
        /* EEE and EEO */
        for (k = 0; k < 4; k++) {
            EEE[k] = EE[k] + EE[7 - k];
            EEO[k] = EE[k] - EE[7 - k];
        }
        /* EEEE and EEEO */
        EEEE[0] = EEE[0] + EEE[3];
        EEEO[0] = EEE[0] - EEE[3];
        EEEE[1] = EEE[1] + EEE[2];
        EEEO[1] = EEE[1] - EEE[2];

        dst[0][i] = (gDctMatrix32[0][0] * EEEE[0] + gDctMatrix32[0][1] * EEEE[1] + add) >> shift;
        dst[16][i] = (gDctMatrix32[16][0] * EEEE[0] + gDctMatrix32[16][1] * EEEE[1] + add) >> shift;
        dst[8][i] = (gDctMatrix32[8][0] * EEEO[0] + gDctMatrix32[8][1] * EEEO[1] + add) >> shift;
        dst[24][i] = (gDctMatrix32[24][0] * EEEO[0] + gDctMatrix32[24][1] * EEEO[1] + add) >> shift;
        for (k = 4; k < 32; k += 8) {
            dst[k][i] = (gDctMatrix32[k][0] * EEO[0] + gDctMatrix32[k][1] * EEO[1] + gDctMatrix32[k][2] * EEO[2] +
                         gDctMatrix32[k][3] * EEO[3] + add) >>
                        shift;
        }
        for (k = 2; k < 32; k += 4) {
            dst[k][i] = (gDctMatrix32[k][0] * EO[0] + gDctMatrix32[k][1] * EO[1] + gDctMatrix32[k][2] * EO[2] +
                         gDctMatrix32[k][3] * EO[3] + gDctMatrix32[k][4] * EO[4] + gDctMatrix32[k][5] * EO[5] +
                         gDctMatrix32[k][6] * EO[6] + gDctMatrix32[k][7] * EO[7] + add) >>
                        shift;
        }
        for (k = 1; k < 32; k += 2) {
            dst[k][i] = (gDctMatrix32[k][0] * O[0] + gDctMatrix32[k][1] * O[1] + gDctMatrix32[k][2] * O[2] +
                         gDctMatrix32[k][3] * O[3] + gDctMatrix32[k][4] * O[4] + gDctMatrix32[k][5] * O[5] +
                         gDctMatrix32[k][6] * O[6] + gDctMatrix32[k][7] * O[7] + gDctMatrix32[k][8] * O[8] +
                         gDctMatrix32[k][9] * O[9] + gDctMatrix32[k][10] * O[10] + gDctMatrix32[k][11] * O[11] +
                         gDctMatrix32[k][12] * O[12] + gDctMatrix32[k][13] * O[13] + gDctMatrix32[k][14] * O[14] +
                         gDctMatrix32[k][15] * O[15] + add) >>
                        shift;
        }
    }
}

void partialButterflyInverse32(
    vector<vector<int>> &src, vector<vector<int>> &dst, int shift, const int outputMinimum, const int outputMaximum) {
    int k;
    int E[16], O[16];
    int EE[8], EO[8];
    int EEE[4], EEO[4];
    int EEEE[2], EEEO[2];
    int add = (shift > 0) ? (1 << (shift - 1)) : 0;

    for (int i = 0; i < 32; i++) {
        /* Utilizing symmetry properties to the maximum to minimize the number of multiplications */
        for (k = 0; k < 16; k++) {
            O[k] = gDctMatrix32[1][k] * src[1][i] + gDctMatrix32[3][k] * src[3][i] + gDctMatrix32[5][k] * src[5][i] +
                   gDctMatrix32[7][k] * src[7][i] + gDctMatrix32[9][k] * src[9][i] + gDctMatrix32[11][k] * src[11][i] +
                   gDctMatrix32[13][k] * src[13][i] + gDctMatrix32[15][k] * src[15][i] +
                   gDctMatrix32[17][k] * src[17][i] + gDctMatrix32[19][k] * src[19][i] +
                   gDctMatrix32[21][k] * src[21][i] + gDctMatrix32[23][k] * src[23][i] +
                   gDctMatrix32[25][k] * src[25][i] + gDctMatrix32[27][k] * src[27][i] +
                   gDctMatrix32[29][k] * src[29][i] + gDctMatrix32[31][k] * src[31][i];
        }
        for (k = 0; k < 8; k++) {
            EO[k] = gDctMatrix32[2][k] * src[2][i] + gDctMatrix32[6][k] * src[6][i] + gDctMatrix32[10][k] * src[10][i] +
                    gDctMatrix32[14][k] * src[14][i] + gDctMatrix32[18][k] * src[18][i] +
                    gDctMatrix32[22][k] * src[22][i] + gDctMatrix32[26][k] * src[26][i] +
                    gDctMatrix32[30][k] * src[30][i];
        }
        for (k = 0; k < 4; k++) {
            EEO[k] = gDctMatrix32[4][k] * src[4][i] + gDctMatrix32[12][k] * src[12][i] +
                     gDctMatrix32[20][k] * src[20][i] + gDctMatrix32[28][k] * src[28][i];
        }
        EEEO[0] = gDctMatrix32[8][0] * src[8][i] + gDctMatrix32[24][0] * src[24][i];
        EEEO[1] = gDctMatrix32[8][1] * src[8][i] + gDctMatrix32[24][1] * src[24][i];
        EEEE[0] = gDctMatrix32[0][0] * src[0][i] + gDctMatrix32[16][0] * src[16][i];
        EEEE[1] = gDctMatrix32[0][1] * src[0][i] + gDctMatrix32[16][1] * src[16][i];

        /* Combining even and odd terms at each hierarchy levels to calculate the final spatial
         * domain vector */
        EEE[0] = EEEE[0] + EEEO[0];
        EEE[3] = EEEE[0] - EEEO[0];
        EEE[1] = EEEE[1] + EEEO[1];
        EEE[2] = EEEE[1] - EEEO[1];
        for (k = 0; k < 4; k++) {
            EE[k] = EEE[k] + EEO[k];
            EE[k + 4] = EEE[3 - k] - EEO[3 - k];
        }
        for (k = 0; k < 8; k++) {
            E[k] = EE[k] + EO[k];
            E[k + 8] = EE[7 - k] - EO[7 - k];
        }
        for (k = 0; k < 16; k++) {
            dst[i][k] = clamp((E[k] + O[k] + add) >> shift, outputMinimum, outputMaximum);
            dst[i][k + 16] = clamp((E[15 - k] - O[15 - k] + add) >> shift, outputMinimum, outputMaximum);
        }
    }
}
