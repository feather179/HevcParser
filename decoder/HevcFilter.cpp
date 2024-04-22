#include "HevcFilter.h"

#include "HevcFrame.h"
#include "HevcSlice.h"
#include "HevcCtu.h"
#include "HevcCu.h"
#include "HevcTu.h"
#include "HevcPu.h"

#include <algorithm>

using std::shared_ptr;
using std::vector;
using std::clamp;
using std::abs;

const static uint8_t betaTable[] = {
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15,
    16, 17, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 64,
};

const static uint8_t tcTable[] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,  1,  1,  1,  1,  1,  1,  1,  1,
    2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6, 6, 7, 8, 9, 10, 11, 13, 14, 16, 18, 20, 22, 24,
};

const static int qpiTable[] = {
    29, 30, 31, 32, 33, 33, 34, 34, 35, 35, 36, 36, 37, 37,
};

void HevcFilter::deblockingFilter(shared_ptr<HevcCtu> ctu, DeblockEdgeDir edgeDir) {
    auto sliceHeader = ctu->getSlice()->getSliceHeader();
    if (sliceHeader->slice_deblocking_filter_disabled_flag) return;

    auto frame = ctu->getFrame();
    auto sps = frame->getSps();
    auto pps = frame->getPps();

    auto filterCu = [this, frame, sliceHeader, sps, pps](shared_ptr<HevcCu> cu, DeblockEdgeDir edgeDir) {
        int x = cu->getX(), y = cu->getY();
        int log2Size = cu->getLog2Size();
        int nCbS = 1 << log2Size;
        int filterEdgeFlag = 1;
        auto edgeFlags = vector<vector<uint8_t>>(nCbS, vector<uint8_t>(nCbS, 0));
        auto bS = vector<vector<uint8_t>>(nCbS, vector<uint8_t>(nCbS, 0));

        if (edgeDir == EDGE_VER) {
            if (x == 0) filterEdgeFlag = 0;
            if (!pps->loop_filter_across_tiles_enabled_flag && (x == frame->getTileAddr(x, y).x)) filterEdgeFlag = 0;
            if (!sliceHeader->slice_loop_filter_across_slices_enabled_flag && (x == frame->getSliceAddr(x, y).x))
                filterEdgeFlag = 0;
        } else {
            if (y == 0) filterEdgeFlag = 0;
            if (!pps->loop_filter_across_tiles_enabled_flag && (y == frame->getTileAddr(x, y).y)) filterEdgeFlag = 0;
            if (!sliceHeader->slice_loop_filter_across_slices_enabled_flag && (y == frame->getSliceAddr(x, y).y))
                filterEdgeFlag = 0;
        }

        getTbBoundary(cu, 0, 0, log2Size, 0, filterEdgeFlag, edgeDir, edgeFlags);
        getPbBoundary(cu, edgeDir, edgeFlags);
        getBoundaryStrength(frame, cu, edgeDir, edgeFlags, bS);

        int nD = 1 << (log2Size - 3);
        // Luma
        {
            if (edgeDir == EDGE_VER) {
                for (int m = 0; m < nD * 2; m++) {
                    int yD = m << 2;
                    for (int k = 0; k < nD; k++) {
                        int xD = k << 3;
                        if (bS[yD][xD] > 0) filterLumaBlockEdge(frame, cu, xD, yD, edgeDir, bS[yD][xD]);
                    }
                }
            } else {
                for (int k = 0; k < nD * 2; k++) {
                    int xD = k << 2;
                    for (int m = 0; m < nD; m++) {
                        int yD = m << 3;
                        if (bS[yD][xD] > 0) filterLumaBlockEdge(frame, cu, xD, yD, edgeDir, bS[yD][xD]);
                    }
                }
            }
        }

        // Chroma
        {
            uint8_t subWidthC = sps->subWidthC, subHeightC = sps->subHeightC;
            int edgeSpacing = 8 / subWidthC;
            int edgeSections = nD * (2 / subHeightC);

            if (edgeDir == EDGE_VER) {
                for (int m = 0; m < edgeSections; m++) {
                    int yD = m << 2;
                    for (int k = 0; k < nD; k++) {
                        int xD = k * edgeSpacing;
                        if (bS[yD * subHeightC][xD * subHeightC] == 2 && ((x / subWidthC + xD) & 7) == 0) {
                            filterChromaBlockEdge(frame, cu, xD, yD, edgeDir, 1, pps->pps_cb_qp_offset);
                            filterChromaBlockEdge(frame, cu, xD, yD, edgeDir, 2, pps->pps_cr_qp_offset);
                        }
                    }
                }
            } else {
                for (int k = 0; k < edgeSections; k++) {
                    int xD = k << 2;
                    for (int m = 0; m < nD; m++) {
                        int yD = m * edgeSpacing;
                        if (bS[yD * subHeightC][xD * subHeightC] == 2 && ((y / subHeightC + yD) & 7) == 0) {
                            filterChromaBlockEdge(frame, cu, xD, yD, edgeDir, 1, pps->pps_cb_qp_offset);
                            filterChromaBlockEdge(frame, cu, xD, yD, edgeDir, 2, pps->pps_cr_qp_offset);
                        }
                    }
                }
            }
        }
    };

    int cuCount = ctu->getCuCount();
    for (int i = 0; i < cuCount; i++) filterCu(ctu->getCu(i), edgeDir);
}

void HevcFilter::getTbBoundary(shared_ptr<HevcCu> cu,
                               int xB0,
                               int yB0,
                               int log2TrafoSize,
                               int trafoDepth,
                               int filterEdgeFlag,
                               DeblockEdgeDir edgeDir,
                               vector<vector<uint8_t>> &edgeFlags) {
    int xCb = cu->getX(), yCb = cu->getY();
    auto tu = cu->getTu(xCb + xB0, yCb + yB0, trafoDepth);
    bool splitTransformFlag;
    if (tu)
        splitTransformFlag = tu->getSplitTransformFlag();
    else
        splitTransformFlag = false;

    if (splitTransformFlag) {
        int xB1 = xB0 + (1 << (log2TrafoSize - 1));
        int yB1 = yB0 + (1 << (log2TrafoSize - 1));

        getTbBoundary(cu, xB0, yB0, log2TrafoSize - 1, trafoDepth + 1, filterEdgeFlag, edgeDir, edgeFlags);
        getTbBoundary(cu, xB1, yB0, log2TrafoSize - 1, trafoDepth + 1, filterEdgeFlag, edgeDir, edgeFlags);
        getTbBoundary(cu, xB0, yB1, log2TrafoSize - 1, trafoDepth + 1, filterEdgeFlag, edgeDir, edgeFlags);
        getTbBoundary(cu, xB1, yB1, log2TrafoSize - 1, trafoDepth + 1, filterEdgeFlag, edgeDir, edgeFlags);
    } else {
        if (edgeDir == EDGE_VER) {
            for (int k = 0; k < (1 << log2TrafoSize); k++) {
                if (xB0)
                    edgeFlags[yB0 + k][xB0] = 1;
                else
                    edgeFlags[yB0 + k][xB0] = filterEdgeFlag;
            }
        } else {
            for (int k = 0; k < (1 << log2TrafoSize); k++) {
                if (yB0)
                    edgeFlags[yB0][xB0 + k] = 1;
                else
                    edgeFlags[yB0][xB0 + k] = filterEdgeFlag;
            }
        }
    }
}

void HevcFilter::getPbBoundary(shared_ptr<HevcCu> cu, DeblockEdgeDir edgeDir, vector<vector<uint8_t>> &edgeFlags) {
    PartMode partMode = cu->getPartMode();
    int log2Size = cu->getLog2Size();

    if (edgeDir == EDGE_VER) {
        int x = -1;
        if (partMode == PART_MODE_Nx2N || partMode == PART_MODE_NxN)
            x = 1 << (log2Size - 1);
        else if (partMode == PART_MODE_nLx2N)
            x = 1 << (log2Size - 2);
        else if (partMode == PART_MODE_nRx2N)
            x = 3 * (1 << (log2Size - 2));

        if (x >= 0)
            for (int k = 0; k < (1 << log2Size); k++) edgeFlags[k][x] = 1;
    } else {
        int y = -1;
        if (partMode == PART_MODE_2NxN || partMode == PART_MODE_NxN)
            y = 1 << (log2Size - 1);
        else if (partMode == PART_MODE_2NxnU)
            y = 1 << (log2Size - 2);
        else if (partMode == PART_MODE_2NxnD)
            y = 3 * (1 << (log2Size - 2));

        if (y >= 0)
            for (int k = 0; k < (1 << log2Size); k++) edgeFlags[y][k] = 1;
    }
}

void HevcFilter::getBoundaryStrength(shared_ptr<HevcFrame> frame,
                                     shared_ptr<HevcCu> cu,
                                     DeblockEdgeDir edgeDir,
                                     vector<vector<uint8_t>> &edgeFlags,
                                     vector<vector<uint8_t>> &bS) {
    auto &recPicture = frame->mReconY;
    int xCb = cu->getX(), yCb = cu->getY();
    int log2Size = cu->getLog2Size();
    int xN = (edgeDir == EDGE_VER) ? (1 << (log2Size - 3)) : (1 << (log2Size - 2));
    int yN = (edgeDir == EDGE_VER) ? (1 << (log2Size - 2)) : (1 << (log2Size - 3));

    auto getInterStrength = [frame](shared_ptr<HevcPu> puP, shared_ptr<HevcPu> puQ) -> uint8_t {
        MvField mvFieldP, mvFieldQ;
        puP->getMvField(mvFieldP);
        puQ->getMvField(mvFieldQ);

        int mvCountP = mvFieldP.predFlagLX[0] + mvFieldP.predFlagLX[1];
        int mvCountQ = mvFieldQ.predFlagLX[0] + mvFieldQ.predFlagLX[1];

        if (mvCountP != mvCountQ) return 1;

        // same motion vector count
        if (mvCountP == 2) {
            int refPocP0 = -1, refPocP1 = -1, refPocQ0 = -1, refPocQ1 = -1;
            if (mvFieldP.refIdxLX[0] >= 0) refPocP0 = frame->mRefPicList[0][mvFieldP.refIdxLX[0]]->getPoc();
            if (mvFieldP.refIdxLX[1] >= 0) refPocP1 = frame->mRefPicList[1][mvFieldP.refIdxLX[1]]->getPoc();
            if (mvFieldQ.refIdxLX[0] >= 0) refPocQ0 = frame->mRefPicList[0][mvFieldQ.refIdxLX[0]]->getPoc();
            if (mvFieldQ.refIdxLX[1] >= 0) refPocQ1 = frame->mRefPicList[1][mvFieldQ.refIdxLX[1]]->getPoc();

            if ((refPocP0 == refPocQ0 && refPocP1 == refPocQ1) || (refPocP0 == refPocQ1 && refPocP1 == refPocQ0)) {
                if (refPocP0 != refPocP1) {
                    if (refPocP0 == refPocQ0) {
                        if (abs(mvFieldP.mvLX[0].x - mvFieldQ.mvLX[0].x) >= 4) return 1;
                        if (abs(mvFieldP.mvLX[0].y - mvFieldQ.mvLX[0].y) >= 4) return 1;
                        if (abs(mvFieldP.mvLX[1].x - mvFieldQ.mvLX[1].x) >= 4) return 1;
                        if (abs(mvFieldP.mvLX[1].y - mvFieldQ.mvLX[1].y) >= 4) return 1;
                    } else {
                        if (abs(mvFieldP.mvLX[0].x - mvFieldQ.mvLX[1].x) >= 4) return 1;
                        if (abs(mvFieldP.mvLX[0].y - mvFieldQ.mvLX[1].y) >= 4) return 1;
                        if (abs(mvFieldP.mvLX[1].x - mvFieldQ.mvLX[0].x) >= 4) return 1;
                        if (abs(mvFieldP.mvLX[1].y - mvFieldQ.mvLX[0].y) >= 4) return 1;
                    }
                } else {
                    if (((abs(mvFieldP.mvLX[0].x - mvFieldQ.mvLX[0].x) >= 4) ||
                         (abs(mvFieldP.mvLX[0].y - mvFieldQ.mvLX[0].y) >= 4) ||
                         (abs(mvFieldP.mvLX[1].x - mvFieldQ.mvLX[1].x) >= 4) ||
                         (abs(mvFieldP.mvLX[1].y - mvFieldQ.mvLX[1].y) >= 4)) &&
                        ((abs(mvFieldP.mvLX[0].x - mvFieldQ.mvLX[1].x) >= 4) ||
                         (abs(mvFieldP.mvLX[0].y - mvFieldQ.mvLX[1].y) >= 4) ||
                         (abs(mvFieldP.mvLX[1].x - mvFieldQ.mvLX[0].x) >= 4) ||
                         (abs(mvFieldP.mvLX[1].y - mvFieldQ.mvLX[0].y) >= 4)))
                        return 1;
                }
            } else  // different ref picture
                return 1;
        } else {
            int refPocP0, refPocQ0;
            Mv mvP0, mvQ0;
            int X = 0;
            if (mvFieldP.predFlagLX[1]) X = 1;
            refPocP0 = frame->mRefPicList[X][mvFieldP.refIdxLX[X]]->getPoc();
            mvP0 = mvFieldP.mvLX[X];

            X = 0;
            if (mvFieldQ.predFlagLX[1]) X = 1;
            refPocQ0 = frame->mRefPicList[X][mvFieldQ.refIdxLX[X]]->getPoc();
            mvQ0 = mvFieldQ.mvLX[X];

            if (refPocP0 != refPocQ0) return 1;
            if (abs(mvP0.x - mvQ0.x) >= 4) return 1;
            if (abs(mvP0.y - mvQ0.y) >= 4) return 1;
        }

        return 0;
    };

    for (int j = 0; j < yN; j++) {
        int yD = (edgeDir == EDGE_VER) ? (j << 2) : (j << 3);
        for (int i = 0; i < xN; i++) {
            int xD = (edgeDir == EDGE_VER) ? (i << 3) : (i << 2);

            if (edgeFlags[yD][xD]) {
                int xP0 = (edgeDir == EDGE_VER) ? xCb + xD - 1 : xCb + xD;
                int yP0 = (edgeDir == EDGE_VER) ? yCb + yD : yCb + yD - 1;
                int xQ0 = xCb + xD, yQ0 = yCb + yD;
                uint16_t p0 = recPicture[yP0][xP0];
                uint16_t q0 = recPicture[yQ0][xQ0];
                auto cuP0 = frame->getCu(xP0, yP0);
                auto cuQ0 = frame->getCu(xQ0, yQ0);

                if (cuP0->getPredMode() == PRED_MODE_INTRA || cuQ0->getPredMode() == PRED_MODE_INTRA)
                    bS[yD][xD] = 2;
                else {
                    auto tuP0 = cuP0->getTu(xP0, yP0);
                    auto tuQ0 = cuQ0->getTu(xQ0, yQ0);
                    if ((tuP0 && tuP0->getCbfLuma()) || (tuQ0 && tuQ0->getCbfLuma()))
                        bS[yD][xD] = 1;
                    else {
                        auto puP = cuP0->getPu(xP0, yP0);
                        auto puQ = cuQ0->getPu(xQ0, yQ0);
                        bS[yD][xD] = getInterStrength(puP, puQ);
                    }
                }
            }
        }
    }
}

void HevcFilter::filterLumaBlockEdge(
    shared_ptr<HevcFrame> frame, shared_ptr<HevcCu> cu, int xBl, int yBl, DeblockEdgeDir edgeDir, uint8_t bS) {
    int x = cu->getX() + xBl;
    int y = cu->getY() + yBl;
    auto sliceHeader = frame->getCtu(x, y)->getSlice()->getSliceHeader();
    auto sps = frame->getSps();
    int tcOffset = sliceHeader->slice_tc_offset_div2, betaOffset = sliceHeader->slice_beta_offset_div2;
    auto &recPicture = frame->mReconY;
    int qpL;
    {
        auto cuQ = frame->getCu(x, y);
        auto cuP = (edgeDir == EDGE_VER) ? frame->getCu(x - 1, y) : frame->getCu(x, y - 1);
        qpL = (cuQ->getQpY() + cuP->getQpY() + 1) >> 1;
    }
    uint16_t p[4][4], q[4][4];
    if (edgeDir == EDGE_VER) {
        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < 4; k++) {
                p[i][k] = recPicture[y + k][x - i - 1];
                q[i][k] = recPicture[y + k][x + i];
            }
        }
    } else {
        for (int i = 0; i < 4; i++) {
            for (int k = 0; k < 4; k++) {
                p[i][k] = recPicture[y - i - 1][x + k];
                q[i][k] = recPicture[y + i][x + k];
            }
        }
    }

    int beta = betaTable[clamp(qpL + (betaOffset << 1), 0, 51)] * (1 << (sps->bitDepthY - 8));
    int tc = tcTable[clamp(qpL + 2 * (bS - 1) + (tcOffset << 1), 0, 53)] * (1 << (sps->bitDepthY - 8));
    int dE = 0, dEp = 0, dEq = 0;

    auto decideLumaSample = [](uint16_t p0, uint16_t p3, uint16_t q0, uint16_t q3, int dpq, int beta, int tc) -> int {
        if (dpq >= (beta >> 2)) return 0;
        if (abs(p3 - p0) + abs(q0 - q3) >= (beta >> 3)) return 0;
        if (abs(p0 - q0) >= ((5 * tc + 1) >> 1)) return 0;

        return 1;
    };

    {
        int dp0 = abs(p[2][0] - 2 * p[1][0] + p[0][0]);
        int dp3 = abs(p[2][3] - 2 * p[1][3] + p[0][3]);
        int dq0 = abs(q[2][0] - 2 * q[1][0] + q[0][0]);
        int dq3 = abs(q[2][3] - 2 * q[1][3] + q[0][3]);
        int dpq0 = dp0 + dq0;
        int dpq3 = dp3 + dq3;
        int dp = dp0 + dp3;
        int dq = dq0 + dq3;
        int d = dpq0 + dpq3;

        if (d < beta) {
            int dSam0 = decideLumaSample(p[0][0], p[3][0], q[0][0], q[3][0], 2 * dpq0, beta, tc);
            int dSam3 = decideLumaSample(p[0][3], p[3][3], q[0][3], q[3][3], 2 * dpq3, beta, tc);

            dE = 1;
            if (dSam0 == 1 && dSam3 == 1) dE = 2;
            if (dp < ((beta + (beta >> 1)) >> 3)) dEp = 1;
            if (dq < ((beta + (beta >> 1)) >> 3)) dEq = 1;
        }
    }

    auto filterLumaSample = [&recPicture, sps, edgeDir, p, q, dE, dEp, dEq, tc](int k, vector<uint16_t> &filterP,
                                                                                vector<uint16_t> &filterQ) {
        uint16_t p0 = p[0][k], p1 = p[1][k], p2 = p[2][k], p3 = p[3][k];
        uint16_t q0 = q[0][k], q1 = q[1][k], q2 = q[2][k], q3 = q[3][k];

        if (dE == 2) {
            filterP.push_back(clamp((p2 + 2 * p1 + 2 * p0 + 2 * q0 + q1 + 4) >> 3, p0 - 2 * tc, p0 + 2 * tc));
            filterP.push_back(clamp((p2 + p1 + p0 + q0 + 2) >> 2, p1 - 2 * tc, p1 + 2 * tc));
            filterP.push_back(clamp((2 * p3 + 3 * p2 + p1 + p0 + q0 + 4) >> 3, p2 - 2 * tc, p2 + 2 * tc));
            filterQ.push_back(clamp((p1 + 2 * p0 + 2 * q0 + 2 * q1 + q2 + 4) >> 3, q0 - 2 * tc, q0 + 2 * tc));
            filterQ.push_back(clamp((p0 + q0 + q1 + q2 + 2) >> 2, q1 - 2 * tc, q1 + 2 * tc));
            filterQ.push_back(clamp((p0 + q0 + q1 + 3 * q2 + 2 * q3 + 4) >> 3, q2 - 2 * tc, q2 + 2 * tc));
        } else {
            int delta = (9 * (q0 - p0) - 3 * (q1 - p1) + 8) >> 4;
            if (abs(delta) < tc * 10) {
                delta = clamp(delta, -tc, tc);
                filterP.push_back(clamp(p0 + delta, 0, (1 << sps->bitDepthY) - 1));
                filterQ.push_back(clamp(q0 - delta, 0, (1 << sps->bitDepthY) - 1));

                if (dEp == 1) {
                    int deltaP = clamp((((p2 + p0 + 1) >> 1) - p1 + delta) >> 1, -(tc >> 1), tc >> 1);
                    filterP.push_back(clamp(p1 + deltaP, 0, (1 << sps->bitDepthY) - 1));
                }
                if (dEq == 1) {
                    int deltaQ = clamp((((q2 + q0 + 1) >> 1) - q1 - delta) >> 1, -(tc >> 1), tc >> 1);
                    filterQ.push_back(clamp(q1 + deltaQ, 0, (1 << sps->bitDepthY) - 1));
                }
            }
        }
    };

    if (dE != 0) {
        for (int k = 0; k < 4; k++) {
            vector<uint16_t> filterP, filterQ;
            filterLumaSample(k, filterP, filterQ);

            auto cuP = (edgeDir == EDGE_VER) ? frame->getCu(x - 1, y + k) : frame->getCu(x + k, y - 1);
            auto cuQ = (edgeDir == EDGE_VER) ? frame->getCu(x, y + k) : frame->getCu(x + k, y);

            if (sps->pcm_loop_filter_disabled_flag && cuP->getPcmFlag()) filterP.clear();
            if (cuP->getCuTransquantBypassFlag()) filterP.clear();
            // TODO: palette mode
            if (sps->pcm_loop_filter_disabled_flag && cuQ->getPcmFlag()) filterQ.clear();
            if (cuQ->getCuTransquantBypassFlag()) filterQ.clear();
            // TODO: palette mode

            if (edgeDir == EDGE_VER) {
                for (int i = 0; i < filterP.size(); i++) recPicture[y + k][x - i - 1] = filterP[i];
                for (int j = 0; j < filterQ.size(); j++) recPicture[y + k][x + j] = filterQ[j];
            } else {
                for (int i = 0; i < filterP.size(); i++) recPicture[y - i - 1][x + k] = filterP[i];
                for (int j = 0; j < filterQ.size(); j++) recPicture[y + j][x + k] = filterQ[j];
            }
        }
    }
}

void HevcFilter::filterChromaBlockEdge(shared_ptr<HevcFrame> frame,
                                       shared_ptr<HevcCu> cu,
                                       int xBl,
                                       int yBl,
                                       DeblockEdgeDir edgeDir,
                                       int cIdx,
                                       int qpOffset) {
    auto sps = frame->getSps();
    uint8_t subWidthC = sps->subWidthC, subHeightC = sps->subHeightC;
    int x = (cu->getX() / subWidthC) + xBl;
    int y = (cu->getY() / subHeightC) + yBl;
    auto sliceHeader = frame->getCtu(x, y)->getSlice()->getSliceHeader();
    int tcOffset = sliceHeader->slice_tc_offset_div2;
    auto &recPicture = (cIdx == 1) ? frame->mReconCb : frame->mReconCr;
    int qpi, qpC;
    {
        auto cuQ = frame->getCu(x * subWidthC, y * subHeightC);
        auto cuP = (edgeDir == EDGE_VER) ? frame->getCu(x * subWidthC - 1, y * subHeightC)
                                         : frame->getCu(x * subWidthC, y * subHeightC - 1);
        qpi = ((cuQ->getQpY() + cuP->getQpY() + 1) >> 1) + qpOffset;

        if (sps->chromaArrayType == 1) {
            if (qpi < 30)
                qpC = qpi;
            else if (qpi >= 30 && qpi <= 43)
                qpC = qpiTable[qpi - 30];
            else
                qpC = qpi - 6;
        } else
            qpC = std::min(qpi, 51);
    }
    int tc = tcTable[clamp(qpC + 2 + (tcOffset << 1), 0, 53)] * (1 << (sps->bitDepthC - 8));

    uint16_t p[2][4], q[2][4];
    if (edgeDir == EDGE_VER) {
        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < 4; k++) {
                p[i][k] = recPicture[y + k][x - i - 1];
                q[i][k] = recPicture[y + k][x + i];
            }
        }
    } else {
        for (int i = 0; i < 2; i++) {
            for (int k = 0; k < 4; k++) {
                p[i][k] = recPicture[y - i - 1][x + k];
                q[i][k] = recPicture[y + i][x + k];
            }
        }
    }

    auto filterChromaSample = [sps, p, q](int k, int tc, uint16_t &filterP, uint16_t &filterQ) {
        int delta = clamp((((q[0][k] - p[0][k]) << 2) + p[1][k] - q[1][k] + 4) >> 3, -tc, tc);
        filterP = clamp(p[0][k] + delta, 0, (1 << sps->bitDepthC) - 1);
        filterQ = clamp(q[0][k] - delta, 0, (1 << sps->bitDepthC) - 1);
    };

    for (int k = 0; k < 4; k++) {
        uint16_t filterP, filterQ;
        filterChromaSample(k, tc, filterP, filterQ);

        auto cuP = (edgeDir == EDGE_VER) ? frame->getCu((x - 1) * subWidthC, (y + k) * subHeightC)
                                         : frame->getCu((x + k) * subWidthC, (y - 1) * subHeightC);
        auto cuQ = (edgeDir == EDGE_VER) ? frame->getCu(x * subWidthC, (y + k) * subHeightC)
                                         : frame->getCu((x + k) * subWidthC, y * subHeightC);

        if (sps->pcm_loop_filter_disabled_flag && cuP->getPcmFlag()) filterP = p[0][k];
        if (cuP->getCuTransquantBypassFlag()) filterP = p[0][k];
        // TODO: palette mode
        if (sps->pcm_loop_filter_disabled_flag && cuQ->getPcmFlag()) filterQ = q[0][k];
        if (cuQ->getCuTransquantBypassFlag()) filterQ = q[0][k];
        // TODO: palette mode

        if (edgeDir == EDGE_VER) {
            recPicture[y + k][x - 1] = filterP;
            recPicture[y + k][x] = filterQ;
        } else {
            recPicture[y - 1][x + k] = filterP;
            recPicture[y][x + k] = filterQ;
        }
    }
}

void HevcFilter::deblockingFilter(shared_ptr<HevcFrame> frame) {
    int ctuCount = frame->getCtuCount();
    for (int i = 0; i < ctuCount; i++) deblockingFilter(frame->getCtu(i), EDGE_VER);
    for (int i = 0; i < ctuCount; i++) deblockingFilter(frame->getCtu(i), EDGE_HOR);
}

void HevcFilter::sao(shared_ptr<HevcCtu> ctu, int cIdx) {
    auto saoParam = ctu->getSaoParam(cIdx);

    auto frame = ctu->getFrame();
    auto sps = frame->getSps();
    auto pps = frame->getPps();
    auto &recPicture = (cIdx == 0) ? frame->mReconY : (cIdx == 1) ? frame->mReconCb : frame->mReconCr;
    auto &saoPicture = (cIdx == 0) ? frame->mSaoY : (cIdx == 1) ? frame->mSaoCb : frame->mSaoCr;
    int bitDepth = (cIdx == 0) ? sps->bitDepthY : sps->bitDepthC;
    int subWidthC = sps->subWidthC, subHeightC = sps->subHeightC;

    auto copyPixel = [&recPicture, &saoPicture, subWidthC, subHeightC, cIdx](int x, int y, int log2Size) {
        int width = 1 << log2Size, height = 1 << log2Size;
        if (cIdx > 0) {
            x = x / subWidthC;
            y = y / subHeightC;
            width = width / subWidthC;
            height = height / subHeightC;
        }

        if (x + width > saoPicture[0].size()) width = (int)saoPicture[0].size() - x;
        if (y + height > saoPicture.size()) height = (int)saoPicture.size() - y;

        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) saoPicture[y + j][x + i] = recPicture[y + j][x + i];
        }
    };

    if (saoParam.mode == SAO_MODE_OFF) {
        copyPixel(ctu->getX(), ctu->getY(), ctu->getLog2Size());
        return;
    }

    uint8_t bandTable[32] = {0};
    if (saoParam.mode == SAO_MODE_BO) {
        int saoLeftClass = saoParam.bandPosition;
        for (int k = 0; k < 4; k++) bandTable[(k + saoLeftClass) & 31] = k + 1;
    }

    auto saoCuBO = [&recPicture, &saoPicture, saoParam, bitDepth, subWidthC, subHeightC, cIdx,
                    bandTable](shared_ptr<HevcCu> cu) {
        int log2Size = cu->getLog2Size();
        int bandShift = bitDepth - 5;
        int maxValue = (1 << bitDepth) - 1;
        int x = cu->getX(), y = cu->getY();
        int width = 1 << log2Size, height = 1 << log2Size;
        if (cIdx > 0) {
            x = x / subWidthC;
            y = y / subHeightC;
            width = width / subWidthC;
            height = height / subHeightC;
        }

        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                uint16_t pixValue = recPicture[y + j][x + i];
                uint8_t bandIdx = bandTable[pixValue >> bandShift];
                int offset = bandIdx > 0 ? saoParam.offset[bandIdx - 1] : 0;
                saoPicture[y + j][x + i] = clamp(pixValue + offset, 0, maxValue);
            }
        }
    };

    int hPos[2], vPos[2];
    switch (saoParam.mode) {
        case SAO_MODE_EO_0:
            hPos[0] = -1;
            hPos[1] = 1;
            vPos[0] = vPos[1] = 0;
            break;
        case SAO_MODE_EO_45:
            hPos[0] = vPos[1] = 1;
            hPos[1] = vPos[0] = -1;
            break;
        case SAO_MODE_EO_90:
            hPos[0] = hPos[1] = 0;
            vPos[0] = -1;
            vPos[1] = 1;
            break;
        case SAO_MODE_EO_135:
            hPos[0] = vPos[0] = -1;
            hPos[1] = vPos[1] = 1;
            break;
        default:
            break;
    }

    auto saoCuEdge = [frame, ctu, sps, pps, &recPicture, &saoPicture, saoParam, bitDepth, subWidthC, subHeightC, cIdx,
                      hPos, vPos](shared_ptr<HevcCu> cu) {
        int log2Size = cu->getLog2Size();
        int log2MinTbSize = sps->log2MinTbSize;
        int maxValue = (1 << bitDepth) - 1;
        int frameWidth = (int)saoPicture[0].size(), frameHeight = (int)saoPicture.size();
        int x = cu->getX(), y = cu->getY();
        int width = 1 << log2Size, height = 1 << log2Size;
        if (cIdx > 0) {
            x = x / subWidthC;
            y = y / subHeightC;
            width = width / subWidthC;
            height = height / subHeightC;
        }

        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                int xS = x + i, yS = y + j;
                int xY = (cIdx == 0) ? xS : xS * subWidthC;
                int yY = (cIdx == 0) ? yS : yS * subHeightC;
                int edgeIdx = -1;

                for (int k = 0; k < 2; k++) {
                    int xSk = xS + hPos[k], ySk = yS + vPos[k];
                    int xYk = (cIdx == 0) ? xSk : xSk * subWidthC;
                    int yYk = (cIdx == 0) ? ySk : ySk * subHeightC;

                    if (!(xSk >= 0 && xSk < frameWidth && ySk >= 0 && ySk < frameHeight)) {
                        edgeIdx = 0;
                        break;
                    }

                    auto ctuK = frame->getCtu(xYk, yYk);

                    if (ctu->getSliceRsAddr() != ctuK->getSliceRsAddr()) {
                        if (pps->minTbAddrZs[yYk >> log2MinTbSize][xYk >> log2MinTbSize] <
                                pps->minTbAddrZs[yY >> log2MinTbSize][xY >> log2MinTbSize] &&
                            !frame->getCtu(xY, yY)
                                 ->getSlice()
                                 ->getSliceHeader()
                                 ->slice_loop_filter_across_slices_enabled_flag) {
                            edgeIdx = 0;
                            break;
                        }
                        if (pps->minTbAddrZs[yY >> log2MinTbSize][xY >> log2MinTbSize] <
                                pps->minTbAddrZs[yYk >> log2MinTbSize][xYk >> log2MinTbSize] &&
                            !frame->getCtu(xYk, yYk)
                                 ->getSlice()
                                 ->getSliceHeader()
                                 ->slice_loop_filter_across_slices_enabled_flag) {
                            edgeIdx = 0;
                            break;
                        }
                    }

                    if (!pps->loop_filter_across_tiles_enabled_flag && ctu->getTileId() != ctuK->getTileId()) {
                        edgeIdx = 0;
                        break;
                    }
                }

                if (edgeIdx < 0) {
                    edgeIdx = 2;
                    for (int k = 0; k < 2; k++) {
                        if (recPicture[yS][xS] < recPicture[yS + vPos[k]][xS + hPos[k]])
                            edgeIdx--;
                        else if (recPicture[yS][xS] > recPicture[yS + vPos[k]][xS + hPos[k]])
                            edgeIdx++;
                    }
                    if (edgeIdx >= 0 && edgeIdx <= 2) edgeIdx = (edgeIdx == 2) ? 0 : edgeIdx + 1;
                }

                int offset = edgeIdx > 0 ? saoParam.offset[edgeIdx - 1] : 0;
                saoPicture[yS][xS] = clamp(recPicture[yS][xS] + offset, 0, maxValue);
            }
        }
    };

    int cuCount = ctu->getCuCount();
    for (int i = 0; i < cuCount; i++) {
        auto cu = ctu->getCu(i);

        if ((sps->pcm_loop_filter_disabled_flag && cu->getPcmFlag()) || cu->getCuTransquantBypassFlag()) {
            copyPixel(cu->getX(), cu->getY(), cu->getLog2Size());
            continue;
        }

        if (saoParam.mode == SAO_MODE_BO)
            saoCuBO(cu);
        else
            saoCuEdge(cu);
    }
}

void HevcFilter::sao(shared_ptr<HevcFrame> frame) {
    if (!frame->getSps()->sample_adaptive_offset_enabled_flag) {
        for (int i = 0; i < frame->mSaoY.size(); i++)
            frame->mSaoY[i].assign(frame->mReconY[i].begin(), frame->mReconY[i].end());
        for (int i = 0; i < frame->mSaoCb.size(); i++)
            frame->mSaoCb[i].assign(frame->mReconCb[i].begin(), frame->mReconCb[i].end());
        for (int i = 0; i < frame->mSaoCr.size(); i++)
            frame->mSaoCr[i].assign(frame->mReconCr[i].begin(), frame->mReconCr[i].end());

        return;
    }

    int ctuCount = frame->getCtuCount();
    for (int i = 0; i < ctuCount; i++) {
        auto ctu = frame->getCtu(i);
        sao(ctu, 0);  // Y
        sao(ctu, 1);  // Cb
        sao(ctu, 2);  // Cr
    }
}
