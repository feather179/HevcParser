#include "HevcCuDecoder.h"

#include "CabacReader.h"
#include "HevcFrame.h"
#include "HevcSlice.h"
#include "HevcCtu.h"
#include "HevcCu.h"
#include "HevcTu.h"
#include "HevcPu.h"
#include "HevcDst.h"
#include "HevcDct.h"

#include <algorithm>

using std::shared_ptr;
using std::vector;
using std::clamp;
using std::abs;

bool HevcCuDecoder::checkNbBlockAvail(shared_ptr<HevcFrame> frame, int xAddr, int yAddr, int xNb, int yNb) {
    auto sps = frame->getSps();
    auto pps = frame->getPps();

    if (xNb < 0 || yNb < 0) return false;
    if (xNb >= sps->width || yNb >= sps->height) return false;

    uint32_t minBlockAddr = pps->minTbAddrZs[yAddr >> sps->log2MinTbSize][xAddr >> sps->log2MinTbSize];
    uint32_t minBlockAddrNb = pps->minTbAddrZs[yNb >> sps->log2MinTbSize][xNb >> sps->log2MinTbSize];
    if (minBlockAddrNb > minBlockAddr) return false;

    auto ctu = frame->getCtu(xAddr, yAddr);
    auto ctuNb = frame->getCtu(xNb, yNb);
    if (ctu->getSliceRsAddr() != ctuNb->getSliceRsAddr()) return false;
    if (ctu->getTileId() != ctuNb->getTileId()) return false;

    return true;
}

bool HevcCuDecoder::checkNbPbAvail(shared_ptr<HevcFrame> frame,
                                   int xCb,
                                   int yCb,
                                   int nCbS,
                                   int xPb,
                                   int yPb,
                                   int nPbW,
                                   int nPbH,
                                   int partIdx,
                                   int xNb,
                                   int yNb) {
    bool sameCb = false;
    bool availableN = true;

    if (xCb <= xNb && yCb <= yNb && (xCb + nCbS) > xNb && (yCb + nCbS) > yNb) sameCb = true;

    if (!sameCb)
        availableN = checkNbBlockAvail(frame, xPb, yPb, xNb, yNb);
    else if ((nPbW << 1) == nCbS && (nPbH << 1) == nCbS && partIdx == 1 && (yCb + nPbH) <= yNb && (xCb + nPbW) > xNb)
        availableN = false;

    if (availableN && frame->getCu(xNb, yNb)->getPredMode() == PRED_MODE_INTRA) availableN = false;

    return availableN;
}

void HevcCuDecoder::reconstructIntraCu(shared_ptr<HevcCtu> ctu, shared_ptr<HevcCu> cu, int cIdx) {
    auto frame = ctu->getFrame();
    auto sps = frame->getSps();
    auto pps = frame->getPps();
    int log2PuSize = (cIdx == 0) ? cu->getLog2Size() : cu->getLog2Size() - 1;
    int blkNum = 1;

    if (cu->getPartMode() == PART_MODE_NxN) {
        if (cIdx == 0 || (cIdx > 0 && log2PuSize > 2)) {
            log2PuSize--;
            blkNum = 4;
        }
    }

    auto &pred = (cIdx == 0) ? cu->mPredY : (cIdx == 1) ? cu->mPredCb : cu->mPredCr;
    auto &residual = (cIdx == 0) ? cu->mResidualY : (cIdx == 1) ? cu->mResidualCb : cu->mResidualCr;
    auto &recon = (cIdx == 0) ? frame->mReconY : (cIdx == 1) ? frame->mReconCb : frame->mReconCr;

    int subWidth = (cIdx == 0) ? 1 : 2, subHeight = (cIdx == 0) ? 1 : 2;
    int bitDepth = (cIdx == 0) ? sps->bitDepthY : sps->bitDepthC;
    int nPbS = 1 << log2PuSize;
    int xBase = cu->getX() / subWidth;
    int yBase = cu->getY() / subHeight;

    for (int blkIdx = 0; blkIdx < blkNum; blkIdx++) {
        int refSampleNum = nPbS * 4 + 1;
        vector<int> p(refSampleNum, -1);
        vector<int> filter(refSampleNum, -1);
        int firstAvailPos = -1;
        int rowBase = (blkIdx >> 1) * nPbS, colBase = (blkIdx & 1) * nPbS;
        int x = xBase + colBase;
        int y = yBase + rowBase;

        for (int k = 0; k < refSampleNum; k++) {
            int nbX = x + ((k < (nPbS << 1)) ? -1 : k - (nPbS << 1) - 1);
            int nbY = y + ((k < (nPbS << 1)) ? ((nPbS << 1) - k - 1) : -1);
            if (checkNbBlockAvail(frame, x * subWidth, y * subHeight, nbX * subWidth, nbY * subHeight)) {
                if (!(pps->constrained_intra_pred_flag &&
                      frame->getCu(nbX * subWidth, nbY * subHeight)->getPredMode() != PRED_MODE_INTRA)) {
                    p[k] = recon[nbY][nbX];
                    if (firstAvailPos < 0) firstAvailPos = k;
                }
            }
        }

        if (p[0] < 0) {
            if (firstAvailPos < 0)
                p[0] = 1 << (bitDepth - 1);
            else
                p[0] = p[firstAvailPos];
        }
        for (int j = 1; j < refSampleNum; j++)
            if (p[j] < 0) p[j] = p[j - 1];

        IntraPredMode predMode = (cIdx == 0) ? cu->getIntraPredModeY(blkIdx) : cu->getIntraPredModeC(blkIdx);
        bool filterFlag = false, biIntFlag = false;

        if (cIdx == 0 || sps->chromaArrayType == 3) {
            if (nPbS == 4 || predMode == INTRA_PRED_MODE_DC)
                filterFlag = false;
            else {
                int minDistVerHor = std::min(abs(predMode - 26), abs(predMode - 10));
                int intraHorVerDistThres = (nPbS == 8) ? 7 : ((nPbS == 16) ? 1 : 0);

                if (minDistVerHor > intraHorVerDistThres)
                    filterFlag = true;
                else
                    filterFlag = false;
            }
        }

        if (filterFlag && sps->strong_intra_smoothing_enabled_flag && cIdx == 0 && nPbS == 32 &&
            abs(p[64] + p[128] - 2 * p[96]) < (1 << (bitDepth - 5)) &&
            abs(p[64] + p[0] - 2 * p[32]) < (1 << (bitDepth - 5)))
            biIntFlag = true;

        if (filterFlag) {
            if (biIntFlag) {
                filter[0] = p[0];
                filter[64] = p[64];
                filter[128] = p[128];
                for (int i = 0; i < 63; i++) {
                    filter[63 - i] = ((63 - i) * p[64] + (i + 1) * p[0] + 32) >> 6;
                    filter[65 + i] = ((63 - i) * p[64] + (i + 1) * p[128] + 32) >> 6;
                }
            } else {
                filter[0] = p[0];
                filter[refSampleNum - 1] = p[refSampleNum - 1];
                for (int k = 1; k < (refSampleNum - 1); k++) filter[k] = (p[k - 1] + p[k] * 2 + p[k + 1] + 2) >> 2;
            }
        } else
            for (int k = 0; k < refSampleNum; k++) filter[k] = p[k];

        auto indexFromAddr = [nPbS](int x, int y) -> int { return x < 0 ? (nPbS * 2 - 1 - y) : (nPbS * 2 + 1 + x); };

        const static int intraPredAngle[] = {
            0,   0,   32,  26,  21,  17, 13, 9,  5, 2, 0, -2, -5, -9, -13, -17, -21, -26,
            -32, -26, -21, -17, -13, -9, -5, -2, 0, 2, 5, 9,  13, 17, 21,  26,  32,
        };
        const static int invAngle[] = {
            0,    0,    0,    0,    0,    0,    0,     0,     0, 0, 0, -4096, -1638, -910, -630, -482, -390, -315,
            -256, -315, -390, -482, -630, -910, -1638, -4096, 0, 0, 0, 0,     0,     0,    0,    0,    0,
        };

        if (predMode == INTRA_PRED_MODE_PLANAR) {
            for (int row = 0; row < nPbS; row++) {
                for (int col = 0; col < nPbS; col++)
                    pred[row + rowBase][col + colBase] = ((nPbS - 1 - col) * filter[indexFromAddr(-1, row)] +
                                                          (col + 1) * filter[indexFromAddr(nPbS, -1)] +
                                                          (nPbS - 1 - row) * filter[indexFromAddr(col, -1)] +
                                                          (row + 1) * filter[indexFromAddr(-1, nPbS)] + nPbS) >>
                                                         (log2PuSize + 1);
            }
        } else if (predMode == INTRA_PRED_MODE_DC) {
            int dcVal = 0;
            for (int j = nPbS; j < nPbS * 3 + 1; j++) {
                if (j == nPbS * 2) continue;
                dcVal += filter[j];
            }
            dcVal = (dcVal + nPbS) >> (log2PuSize + 1);

            for (int row = 0; row < nPbS; row++)
                for (int col = 0; col < nPbS; col++) pred[row + rowBase][col + colBase] = dcVal;

            if (cIdx == 0 && nPbS < 32) {
                pred[rowBase][colBase] =
                    (filter[indexFromAddr(-1, 0)] + 2 * dcVal + filter[indexFromAddr(0, -1)] + 2) >> 2;
                for (int row = 1; row < nPbS; row++)
                    pred[row + rowBase][colBase] = (filter[indexFromAddr(-1, row)] + 3 * dcVal + 2) >> 2;
                for (int col = 1; col < nPbS; col++)
                    pred[rowBase][col + colBase] = (filter[indexFromAddr(col, -1)] + 3 * dcVal + 2) >> 2;
            }
        } else if (predMode < 18) {
            vector<int> v(refSampleNum, 0);
            int *ref = &v[nPbS * 2];
            if (intraPredAngle[predMode] < 0) {
                for (int j = 0; j < (nPbS + 1); j++) ref[j] = filter[indexFromAddr(-1, j - 1)];
                if (((nPbS * intraPredAngle[predMode]) >> 5) < -1)
                    for (int j = -1; j >= ((nPbS * intraPredAngle[predMode]) >> 5); j--)
                        ref[j] = filter[indexFromAddr(-1 + ((j * invAngle[predMode] + 128) >> 8), -1)];
            } else {
                for (int j = 0; j < (2 * nPbS + 1); j++) ref[j] = filter[indexFromAddr(-1, j - 1)];
            }

            for (int row = 0; row < nPbS; row++) {
                for (int col = 0; col < nPbS; col++) {
                    int iIdx = ((col + 1) * intraPredAngle[predMode]) >> 5;
                    int iFact = ((col + 1) * intraPredAngle[predMode]) & 31;
                    if (iFact == 0)
                        pred[row + rowBase][col + colBase] = ref[row + iIdx + 1];
                    else
                        pred[row + rowBase][col + colBase] =
                            ((32 - iFact) * ref[row + iIdx + 1] + iFact * ref[row + iIdx + 2] + 16) >> 5;
                }
            }

            if (predMode == INTRA_PRED_MODE_ANGULAR_10 && cIdx == 0 && nPbS < 32) {
                for (int col = 0; col < nPbS; col++)
                    pred[rowBase][col + colBase] =
                        clamp(filter[indexFromAddr(-1, 0)] +
                                  ((filter[indexFromAddr(col, -1)] - filter[indexFromAddr(-1, -1)]) >> 1),
                              0, (1 << bitDepth) - 1);
            }
        } else {
            vector<int> v(refSampleNum, 0);
            int *ref = &v[nPbS * 2];
            if (intraPredAngle[predMode] < 0) {
                for (int j = 0; j < (nPbS + 1); j++) ref[j] = filter[indexFromAddr(j - 1, -1)];
                if (((nPbS * intraPredAngle[predMode]) >> 5) < -1)
                    for (int j = -1; j >= ((nPbS * intraPredAngle[predMode]) >> 5); j--)
                        ref[j] = filter[indexFromAddr(-1, -1 + ((j * invAngle[predMode] + 128) >> 8))];
            } else {
                for (int j = 0; j < (2 * nPbS + 1); j++) ref[j] = filter[indexFromAddr(j - 1, -1)];
            }

            for (int row = 0; row < nPbS; row++) {
                for (int col = 0; col < nPbS; col++) {
                    int iIdx = ((row + 1) * intraPredAngle[predMode]) >> 5;
                    int iFact = ((row + 1) * intraPredAngle[predMode]) & 31;
                    if (iFact == 0)
                        pred[row + rowBase][col + colBase] = ref[col + iIdx + 1];
                    else
                        pred[row + rowBase][col + colBase] =
                            ((32 - iFact) * ref[col + iIdx + 1] + iFact * ref[col + iIdx + 2] + 16) >> 5;
                }
            }

            if (predMode == INTRA_PRED_MODE_ANGULAR_26 && cIdx == 0 && nPbS < 32) {
                for (int row = 0; row < nPbS; row++)
                    pred[row + rowBase][colBase] =
                        clamp(filter[indexFromAddr(0, -1)] +
                                  ((filter[indexFromAddr(-1, row)] - filter[indexFromAddr(-1, -1)]) >> 1),
                              0, (1 << bitDepth) - 1);
            }
        }

        // reconstruction
        for (int row = 0; row < nPbS; row++)
            for (int col = 0; col < nPbS; col++)
                recon[row + y][col + x] =
                    clamp(pred[row + rowBase][col + colBase] + residual[row + rowBase][col + colBase], 0,
                          (1 << bitDepth) - 1);
    }
}

void HevcCuDecoder::reconstructInterCu(shared_ptr<HevcCtu> ctu, shared_ptr<HevcCu> cu, int cIdx) {
    const static int lumaFilter[][8] = {
        {0, 0, 0, 0, 0, 0, 0, 0},
        {-1, 4, -10, 58, 17, -5, 1, 0},
        {-1, 4, -11, 40, 40, -11, 4, -1},
        {0, 1, -5, 17, 58, -10, 4, -1},
    };
    const static int chromaFilter[][4] = {
        {0, 0, 0, 0},     {-2, 58, 10, -2}, {-4, 54, 16, -2}, {-6, 46, 28, -4},
        {-4, 36, 36, -4}, {-4, 28, 46, -6}, {-2, 16, 54, -4}, {-2, 10, 58, -2},
    };

    auto doLumaFilter = [](const int src[8], const int filter[8]) -> int {
        int value = 0;
        for (int i = 0; i < 8; i++) value += src[i] * filter[i];
        return value;
    };

    auto doChromaFilter = [](const int src[4], const int filter[4]) -> int {
        int value = 0;
        for (int i = 0; i < 4; i++) value += src[i] * filter[i];
        return value;
    };

    auto frame = ctu->getFrame();
    auto sps = frame->getSps();
    auto pps = frame->getPps();
    int subWidth = (cIdx == 0) ? 1 : sps->subWidthC, subHeight = (cIdx == 0) ? 1 : sps->subHeightC;
    int frameWidth = sps->width / subWidth, frameHeight = sps->height / subHeight;
    int bitDepth = (cIdx == 0) ? sps->bitDepthY : sps->bitDepthC;
    int maxValue = (1 << bitDepth) - 1;
    auto &pred = (cIdx == 0) ? cu->mPredY : (cIdx == 1) ? cu->mPredCb : cu->mPredCr;
    auto &residual = (cIdx == 0) ? cu->mResidualY : (cIdx == 1) ? cu->mResidualCb : cu->mResidualCr;
    auto &recon = (cIdx == 0) ? frame->mReconY : (cIdx == 1) ? frame->mReconCb : frame->mReconCr;

    // weighted sample prediction
    auto weightedPrediction = [=, &pred](shared_ptr<HevcPu> pu, vector<vector<int>> &predSamplesLX0,
                                         vector<vector<int>> &predSamplesLX1) {
        MvField mvField;
        pu->getMvField(mvField);
        int xPb = pu->getX() / subWidth, yPb = pu->getY() / subHeight;
        int nPbW = pu->getWidth() / subWidth, nPbH = pu->getHeight() / subHeight;
        int xCb = cu->getX() / subWidth, yCb = cu->getY() / subHeight;
        int xBl = xPb - xCb, yBl = yPb - yCb;

        int shift1 = std::max(2, 14 - bitDepth);
        int shift2 = std::max(3, 15 - bitDepth);
        int offset1 = 1 << (shift1 - 1);
        int offset2 = 1 << (shift2 - 1);
        auto sliceHeader = ctu->getSlice()->getSliceHeader();
        int weightedPredFlag = pps->weighted_pred_flag;
        if (sliceHeader->slice_type == B_SLICE) weightedPredFlag = pps->weighted_bipred_flag;

        if (!weightedPredFlag) {
            // Default weighted sample prediction
            if (mvField.predFlagLX[0] && !mvField.predFlagLX[1]) {
                for (int j = 0; j < nPbH; j++) {
                    for (int i = 0; i < nPbW; i++)
                        pred[yBl + j][xBl + i] = clamp((predSamplesLX0[j][i] + offset1) >> shift1, 0, maxValue);
                }
            } else if (!mvField.predFlagLX[0] && mvField.predFlagLX[1]) {
                for (int j = 0; j < nPbH; j++) {
                    for (int i = 0; i < nPbW; i++)
                        pred[yBl + j][xBl + i] = clamp((predSamplesLX1[j][i] + offset1) >> shift1, 0, maxValue);
                }
            } else if (mvField.predFlagLX[0] && mvField.predFlagLX[1]) {
                for (int j = 0; j < nPbH; j++) {
                    for (int i = 0; i < nPbW; i++)
                        pred[yBl + j][xBl + i] =
                            clamp((predSamplesLX0[j][i] + predSamplesLX1[j][i] + offset2) >> shift2, 0, maxValue);
                }
            }
        } else {
            // Explicit weighted sample prediction
            auto predWeightTable = sliceHeader->pred_weight_table;
            int wpOffsetBdShift = bitDepth - 8;
            int log2Wd, w0, w1, o0, o1;
            int refIdxL0 = mvField.refIdxLX[0], refIdxL1 = mvField.refIdxLX[1];

            if (cIdx == 0) {
                log2Wd = predWeightTable.luma_log2_weight_denom + shift1;
                w0 = 1 << predWeightTable.luma_log2_weight_denom;
                o0 = 0;
                if (predWeightTable.luma_weight_l0_flag[refIdxL0]) {
                    w0 += predWeightTable.delta_luma_weight_l0[refIdxL0];
                    o0 = predWeightTable.luma_offset_l0[refIdxL0] << wpOffsetBdShift;
                }

                if (sliceHeader->slice_type == B_SLICE) {
                    w1 = 1 << predWeightTable.luma_log2_weight_denom;
                    o1 = 0;
                    if (predWeightTable.luma_weight_l1_flag[refIdxL1]) {
                        w1 += predWeightTable.delta_luma_weight_l1[refIdxL1];
                        o1 = predWeightTable.luma_offset_l1[refIdxL1] << wpOffsetBdShift;
                    }
                }
            } else {
                log2Wd =
                    predWeightTable.luma_log2_weight_denom + predWeightTable.delta_chroma_log2_weight_denom + shift1;
                w0 = 1 << (predWeightTable.luma_log2_weight_denom + predWeightTable.delta_chroma_log2_weight_denom);
                o0 = 0;

                // TODO: if (predWeightTable.chroma_weight_l0_flag[refIdxL0])

                if (sliceHeader->slice_type == B_SLICE) {
                    w1 = 1 << (predWeightTable.luma_log2_weight_denom + predWeightTable.delta_chroma_log2_weight_denom);
                    o1 = 0;
                    // TODO: if (predWeightTable.chroma_weight_l1_flag[refIdxL1])
                }
            }

            if (mvField.predFlagLX[0] && !mvField.predFlagLX[1]) {
                for (int j = 0; j < nPbH; j++) {
                    for (int i = 0; i < nPbW; i++)
                        pred[yBl + j][xBl + i] =
                            clamp(((predSamplesLX0[j][i] * w0 + (1 << (log2Wd - 1))) >> log2Wd) + o0, 0, maxValue);
                }
            } else if (!mvField.predFlagLX[0] && mvField.predFlagLX[1]) {
                for (int j = 0; j < nPbH; j++) {
                    for (int i = 0; i < nPbW; i++)
                        pred[yBl + j][xBl + i] =
                            clamp(((predSamplesLX1[j][i] * w1 + (1 << (log2Wd - 1))) >> log2Wd) + o1, 0, maxValue);
                }
            } else if (mvField.predFlagLX[0] && mvField.predFlagLX[1]) {
                for (int j = 0; j < nPbH; j++) {
                    for (int i = 0; i < nPbW; i++)
                        pred[yBl + j][xBl + i] =
                            clamp((predSamplesLX0[j][i] * w0 + predSamplesLX1[j][i] * w1 + ((o0 + o1 + 1) << log2Wd)) >>
                                      (log2Wd + 1),
                                  0, maxValue);
                }
            }
        }
    };

    auto getPredSamples = [=](shared_ptr<HevcPu> pu) {
        MvField mvField;
        pu->getMvField(mvField);
        int xPb = pu->getX() / subWidth, yPb = pu->getY() / subHeight;
        int nPbW = pu->getWidth() / subWidth, nPbH = pu->getHeight() / subHeight;
        int shift1 = std::min(4, bitDepth - 8);
        int shift2 = 6;
        int shift3 = std::max(2, 14 - bitDepth);

        vector<vector<int>> predSamplesLX[2] = {vector<vector<int>>(nPbH, vector<int>(nPbW, 0)),
                                                vector<vector<int>>(nPbH, vector<int>(nPbW, 0))};

        for (int X = 0; X < 2; X++) {
            if (!mvField.predFlagLX[X]) continue;

            auto &predSamples = predSamplesLX[X];
            auto refFrame = frame->mRefPicList[X][mvField.refIdxLX[X]];
            auto &refSamples = (cIdx == 0) ? refFrame->mSaoY : (cIdx == 1) ? refFrame->mSaoCb : refFrame->mSaoCr;

            if (cIdx == 0) {
                // Luma
                int xInt = xPb + (mvField.mvLX[X].x >> 2);
                int yInt = yPb + (mvField.mvLX[X].y >> 2);
                int xFrac = mvField.mvLX[X].x & 3;
                int yFrac = mvField.mvLX[X].y & 3;

                if (xFrac == 0 && yFrac == 0) {
                    for (int j = 0; j < nPbH; j++) {
                        for (int i = 0; i < nPbW; i++) {
                            int xA = clamp(xInt + i, 0, frameWidth - 1);
                            int yA = clamp(yInt + j, 0, frameHeight - 1);
                            predSamples[j][i] = refSamples[yA][xA] << shift3;
                        }
                    }
                } else if (xFrac == 0) {
                    int A[8];
                    for (int j = 0; j < nPbH; j++) {
                        for (int i = 0; i < nPbW; i++) {
                            for (int k = 0; k < 8; k++)
                                A[k] = refSamples[clamp(yInt + j + k - 3, 0, frameHeight - 1)]
                                                 [clamp(xInt + i, 0, frameWidth - 1)];
                            predSamples[j][i] = doLumaFilter(A, lumaFilter[yFrac]) >> shift1;
                        }
                    }
                } else if (yFrac == 0) {
                    int A[8];
                    for (int j = 0; j < nPbH; j++) {
                        for (int i = 0; i < nPbW; i++) {
                            for (int k = 0; k < 8; k++)
                                A[k] = refSamples[clamp(yInt + j, 0, frameHeight - 1)]
                                                 [clamp(xInt + i + k - 3, 0, frameWidth - 1)];
                            predSamples[j][i] = doLumaFilter(A, lumaFilter[xFrac]) >> shift1;
                        }
                    }
                } else {
                    int A[8], abc[8];
                    for (int j = 0; j < nPbH; j++) {
                        for (int i = 0; i < nPbW; i++) {
                            for (int k = 0; k < 8; k++) {
                                for (int kk = 0; kk < 8; kk++)
                                    A[kk] = refSamples[clamp(yInt + j + k - 3, 0, frameHeight - 1)]
                                                      [clamp(xInt + i + kk - 3, 0, frameWidth - 1)];
                                abc[k] = doLumaFilter(A, lumaFilter[xFrac]) >> shift1;
                            }
                            predSamples[j][i] = doLumaFilter(abc, lumaFilter[yFrac]) >> shift2;
                        }
                    }
                }
            } else {
                // Chroma
                int xInt = xPb + (mvField.mvCLX[X].x >> 3);
                int yInt = yPb + (mvField.mvCLX[X].y >> 3);
                int xFrac = mvField.mvCLX[X].x & 7;
                int yFrac = mvField.mvCLX[X].y & 7;

                if (xFrac == 0 && yFrac == 0) {
                    for (int j = 0; j < nPbH; j++) {
                        for (int i = 0; i < nPbW; i++) {
                            int xA = clamp(xInt + i, 0, frameWidth - 1);
                            int yA = clamp(yInt + j, 0, frameHeight - 1);
                            predSamples[j][i] = refSamples[yA][xA] << shift3;
                        }
                    }
                } else if (xFrac == 0) {
                    int A[4];
                    for (int j = 0; j < nPbH; j++) {
                        for (int i = 0; i < nPbW; i++) {
                            for (int k = 0; k < 4; k++)
                                A[k] = refSamples[clamp(yInt + j + k - 1, 0, frameHeight - 1)]
                                                 [clamp(xInt + i, 0, frameWidth - 1)];
                            predSamples[j][i] = doChromaFilter(A, chromaFilter[yFrac]) >> shift1;
                        }
                    }
                } else if (yFrac == 0) {
                    int A[4];
                    for (int j = 0; j < nPbH; j++) {
                        for (int i = 0; i < nPbW; i++) {
                            for (int k = 0; k < 4; k++)
                                A[k] = refSamples[clamp(yInt + j, 0, frameHeight - 1)]
                                                 [clamp(xInt + i + k - 1, 0, frameWidth - 1)];
                            predSamples[j][i] = doChromaFilter(A, chromaFilter[xFrac]) >> shift1;
                        }
                    }
                } else {
                    int A[4], abc[4];
                    for (int j = 0; j < nPbH; j++) {
                        for (int i = 0; i < nPbW; i++) {
                            for (int k = 0; k < 4; k++) {
                                for (int kk = 0; kk < 4; kk++)
                                    A[kk] = refSamples[clamp(yInt + j + k - 1, 0, frameHeight - 1)]
                                                      [clamp(xInt + i + kk - 1, 0, frameWidth - 1)];
                                abc[k] = doChromaFilter(A, chromaFilter[xFrac]) >> shift1;
                            }
                            predSamples[j][i] = doChromaFilter(abc, chromaFilter[yFrac]) >> shift2;
                        }
                    }
                }
            }
        }

        weightedPrediction(pu, predSamplesLX[0], predSamplesLX[1]);
    };

    int puCount = cu->getPuCount();
    for (int k = 0; k < puCount; k++) {
        auto pu = cu->getPu(k);
        getPredSamples(pu);

        int xPb = pu->getX() / subWidth, yPb = pu->getY() / subHeight;
        int nPbW = pu->getWidth() / subWidth, nPbH = pu->getHeight() / subHeight;
        int xCb = cu->getX() / subWidth, yCb = cu->getY() / subHeight;
        int xBl = xPb - xCb, yBl = yPb - yCb;

        // reconstruction
        for (int j = 0; j < nPbH; j++) {
            for (int i = 0; i < nPbW; i++)
                recon[yPb + j][xPb + i] = clamp(pred[yBl + j][xBl + i] + residual[yBl + j][xBl + i], 0, maxValue);
        }
    }
}

void HevcCuDecoder::decodeCodingQuadtree(shared_ptr<HevcCtu> ctu, int x0, int y0, int log2CtuSize, int depth) {
    auto frame = ctu->getFrame();
    auto sps = frame->getSps();
    auto pps = frame->getPps();
    int nCbS = 1 << log2CtuSize;
    bool splitCuFlag = false;

    if (x0 >= sps->width || y0 >= sps->height) return;

    if ((x0 + nCbS <= sps->width) && (y0 + nCbS <= sps->height) && log2CtuSize > sps->log2MinCtbSize) {
        bool leftAvail = checkNbBlockAvail(frame, x0, y0, x0 - 1, y0);
        bool aboveAvail = checkNbBlockAvail(frame, x0, y0, x0, y0 - 1);
        int inc = 0;

        if (leftAvail && (frame->getCtDepth(x0 - 1, y0) > depth)) inc++;
        if (aboveAvail && (frame->getCtDepth(x0, y0 - 1) > depth)) inc++;
        splitCuFlag = mCabacReader->parseSplitCuFlag(inc);
    } else {
        if (log2CtuSize > sps->log2MinCtbSize)
            splitCuFlag = true;
        else
            splitCuFlag = false;
    }

    if (pps->cu_qp_delta_enabled_flag && log2CtuSize >= pps->log2MinCuQpDeltaSize) {
        // ctu->setIsCuQpDeltaCoded(x0, y0, 0);
        // ctu->setCuQpDeltaVal(x0, y0, 0);
        ctu->setIsCuQpDeltaCoded(0);
        ctu->setCuQpDeltaVal(0);
    }

    if (splitCuFlag) {
        int x1 = x0 + (1 << (log2CtuSize - 1));
        int y1 = y0 + (1 << (log2CtuSize - 1));
        decodeCodingQuadtree(ctu, x0, y0, log2CtuSize - 1, depth + 1);
        decodeCodingQuadtree(ctu, x1, y0, log2CtuSize - 1, depth + 1);
        decodeCodingQuadtree(ctu, x0, y1, log2CtuSize - 1, depth + 1);
        decodeCodingQuadtree(ctu, x1, y1, log2CtuSize - 1, depth + 1);
    } else
        decodeCodingUnit(ctu, x0, y0, log2CtuSize);
}

void HevcCuDecoder::decodeCodingUnit(shared_ptr<HevcCtu> ctu, int x0, int y0, int log2CuSize) {
    static int indexInCtu = 0;

    if (x0 == ctu->getX() && y0 == ctu->getY()) indexInCtu = 0;

    auto frame = ctu->getFrame();
    auto sps = frame->getSps();
    auto pps = frame->getPps();
    auto sliceHeader = ctu->getSlice()->getSliceHeader();
    int nCbS = 1 << log2CuSize;
    auto cu = std::make_shared<HevcCu>(x0, y0, indexInCtu, log2CuSize, sps->log2CtbSize - log2CuSize);
    ctu->addCu(cu);

    // get predicted qp
    {
        auto qgAddr = frame->getQGAddr(x0, y0);
        int xQg = qgAddr.x, yQg = qgAddr.y;
        int sliceQpY = sliceHeader->sliceQp;
        int qpYPrev, qpYA, qpYB;

        auto sliceAddr = frame->getSliceAddr(x0, y0);
        auto tileAddr = frame->getTileAddr(x0, y0);

        if ((xQg == sliceAddr.x && yQg == sliceAddr.y) || (xQg == tileAddr.x && yQg == tileAddr.y) ||
            (pps->entropy_coding_sync_enabled_flag && xQg == tileAddr.x && yQg == ctu->getY()))
            qpYPrev = sliceQpY;
        else {
            auto prevQgAddr = frame->getPrevQGAddr(x0, y0);
            auto prevCu = frame->getCu(prevQgAddr.x + (1 << pps->log2MinCuQpDeltaSize) - 1,
                                       prevQgAddr.y + (1 << pps->log2MinCuQpDeltaSize) - 1);
            qpYPrev = prevCu->getQpY();
        }

        auto getCtbAddr = [sps, pps](int xAddr, int yAddr) -> int {
            int x = xAddr >> sps->log2MinTbSize;
            int y = yAddr >> sps->log2MinTbSize;
            int minTbAddr = pps->minTbAddrZs[y][x];
            int ctbAddr = minTbAddr >> (2 * (sps->log2CtbSize - sps->log2MinTbSize));
            return ctbAddr;
        };

        if (!checkNbBlockAvail(frame, x0, y0, xQg - 1, yQg) ||
            (pps->ctbAddrRsToTs[ctu->getRsAddr()] != getCtbAddr(xQg - 1, yQg)))
            qpYA = qpYPrev;
        else
            qpYA = frame->getCu(xQg - 1, yQg)->getQpY();

        if (!checkNbBlockAvail(frame, x0, y0, xQg, yQg - 1) ||
            (pps->ctbAddrRsToTs[ctu->getRsAddr()] != getCtbAddr(xQg, yQg - 1)))
            qpYB = qpYPrev;
        else
            qpYB = frame->getCu(xQg, yQg - 1)->getQpY();

        cu->setQpYPred((qpYA + qpYB + 1) >> 1);
        // CuQpDeltaVal default 0, may be changed by other CU in the same QG
        // cu->setQpY(
        //     ((cu->getQpYPred() + ctu->getCuQpDeltaVal(x0, y0) + 52 + 2 * sps->qpBdOffsetY) % (52 + sps->qpBdOffsetY))
        //     - sps->qpBdOffsetY);
        cu->setQpY(((cu->getQpYPred() + ctu->getCuQpDeltaVal() + 52 + 2 * sps->qpBdOffsetY) % (52 + sps->qpBdOffsetY)) -
                   sps->qpBdOffsetY);
    }

    cu->setChromaArrayType(sps->chromaArrayType);

    if (pps->transquant_bypass_enabled_flag)
        cu->setCuTransquantBypassFlag(mCabacReader->parseCuTransquantBypassFlag());
    else
        cu->setCuTransquantBypassFlag(false);

    if (sliceHeader->slice_type != I_SLICE) {
        bool leftAvail = checkNbBlockAvail(frame, x0, y0, x0 - 1, y0);
        bool aboveAvail = checkNbBlockAvail(frame, x0, y0, x0, y0 - 1);
        int inc = 0;

        if (leftAvail && frame->getCu(x0 - 1, y0)->getCuSkipFlag()) inc++;
        if (aboveAvail && frame->getCu(x0, y0 - 1)->getCuSkipFlag()) inc++;
        cu->setCuSkipFlag(mCabacReader->parseCuSkipFlag(inc));
    } else
        cu->setCuSkipFlag(false);

    if (cu->getCuSkipFlag()) {
        cu->setPredMode(PRED_MODE_SKIP);
        cu->setPartMode(PART_MODE_2Nx2N);
        decodePredictionUnit(ctu, cu, x0, y0, nCbS, nCbS, 0);
    } else {
        if (sliceHeader->slice_type != I_SLICE) {
            if (mCabacReader->parsePredModeFlag())
                cu->setPredMode(PRED_MODE_INTRA);
            else
                cu->setPredMode(PRED_MODE_INTER);
        } else
            cu->setPredMode(PRED_MODE_INTRA);

        // TODO: palette mode, skip

        if (cu->getPredMode() != PRED_MODE_INTRA || log2CuSize == sps->log2MinCtbSize)
            cu->setPartMode(
                mCabacReader->parsePartMode(log2CuSize, sps->log2MinCtbSize, sps->amp_enabled_flag, cu->getPredMode()));
        else
            cu->setPartMode(PART_MODE_2Nx2N);

        if (cu->getPredMode() == PRED_MODE_INTRA) {
            // TODO"
            cu->setPcmFlag(false);
            // TODO: pcm releated, skip

            int puNum = cu->getPartMode() == PART_MODE_NxN ? 4 : 1;
            int pbOffset = cu->getPartMode() == PART_MODE_NxN ? (nCbS >> 1) : nCbS;
            vector<uint8_t> prevIntraLumaPredFlag(puNum, 0);
            vector<uint8_t> mpmIdx(puNum, 0);
            vector<uint8_t> remIntraLumaPredMode(puNum, 0);
            vector<uint8_t> intraChromaPredMode(puNum, 0);

            for (int i = 0; i < puNum; i++) prevIntraLumaPredFlag[i] = mCabacReader->parsePrevIntraLumaPredFlag();

            for (int i = 0; i < puNum; i++) {
                if (prevIntraLumaPredFlag[i])
                    mpmIdx[i] = mCabacReader->parseMpmIdx();
                else
                    remIntraLumaPredMode[i] = mCabacReader->parseRemIntraLumaPredMode();
            }

            int blkIdx = 0;
            for (int j = 0; j < nCbS; j += pbOffset) {
                for (int i = 0; i < nCbS; i += pbOffset) {
                    int x = x0 + i;
                    int y = y0 + j;
                    bool leftAvail = checkNbBlockAvail(frame, x, y, x - 1, y);
                    bool aboveAvail = checkNbBlockAvail(frame, x, y, x, y - 1);
                    IntraPredMode leftCandIntraPredMode, aboveCandIntraPredMode, intraPredModeY;

                    if (!leftAvail)
                        leftCandIntraPredMode = INTRA_PRED_MODE_DC;
                    else {
                        auto leftCu = frame->getCu(x - 1, y);
                        if (leftCu->getPredMode() != PRED_MODE_INTRA || leftCu->getPcmFlag() == 1)
                            leftCandIntraPredMode = INTRA_PRED_MODE_DC;
                        else
                            leftCandIntraPredMode = leftCu->getIntraPredModeY(x - 1, y);
                    }

                    if (!aboveAvail)
                        aboveCandIntraPredMode = INTRA_PRED_MODE_DC;
                    else {
                        auto aboveCu = frame->getCu(x, y - 1);
                        if (aboveCu->getPredMode() != PRED_MODE_INTRA || aboveCu->getPcmFlag() == 1)
                            aboveCandIntraPredMode = INTRA_PRED_MODE_DC;
                        else if ((y - 1) < ((y >> sps->log2CtbSize) << sps->log2CtbSize))
                            aboveCandIntraPredMode = INTRA_PRED_MODE_DC;
                        else
                            aboveCandIntraPredMode = aboveCu->getIntraPredModeY(x, y - 1);
                    }

                    IntraPredMode candModeList[3];
                    if (leftCandIntraPredMode == aboveCandIntraPredMode) {
                        if (leftCandIntraPredMode < 2) {
                            candModeList[0] = INTRA_PRED_MODE_PLANAR;
                            candModeList[1] = INTRA_PRED_MODE_DC;
                            candModeList[2] = INTRA_PRED_MODE_ANGULAR_26;
                        } else {
                            candModeList[0] = leftCandIntraPredMode;
                            candModeList[1] = (IntraPredMode)(((leftCandIntraPredMode + 29) % 32) + 2);
                            candModeList[2] = (IntraPredMode)(((leftCandIntraPredMode - 2 + 1) % 32) + 2);
                        }
                    } else {
                        candModeList[0] = leftCandIntraPredMode;
                        candModeList[1] = aboveCandIntraPredMode;
                        if (candModeList[0] != INTRA_PRED_MODE_PLANAR && candModeList[1] != INTRA_PRED_MODE_PLANAR)
                            candModeList[2] = INTRA_PRED_MODE_PLANAR;
                        else if (candModeList[0] != INTRA_PRED_MODE_DC && candModeList[1] != INTRA_PRED_MODE_DC)
                            candModeList[2] = INTRA_PRED_MODE_DC;
                        else
                            candModeList[2] = INTRA_PRED_MODE_ANGULAR_26;
                    }

                    if (prevIntraLumaPredFlag[blkIdx])
                        intraPredModeY = candModeList[mpmIdx[blkIdx]];
                    else {
                        intraPredModeY = (IntraPredMode)remIntraLumaPredMode[blkIdx];
                        if (candModeList[0] > candModeList[1]) std::swap(candModeList[0], candModeList[1]);
                        if (candModeList[0] > candModeList[2]) std::swap(candModeList[0], candModeList[2]);
                        if (candModeList[1] > candModeList[2]) std::swap(candModeList[1], candModeList[2]);
                        for (int i = 0; i < 3; i++)
                            if (intraPredModeY >= candModeList[i]) intraPredModeY = (IntraPredMode)(intraPredModeY + 1);
                    }

                    cu->setIntraPredModeY(blkIdx, intraPredModeY);
                    blkIdx++;
                }
            }

            if (sps->chromaArrayType == 3) {
                for (int i = 0; i < puNum; i++) intraChromaPredMode[i] = mCabacReader->parseIntraChromaPredMode();
            } else if (sps->chromaArrayType != 0) {
                intraChromaPredMode.resize(1);
                intraChromaPredMode[0] = mCabacReader->parseIntraChromaPredMode();
            }

            if (sps->chromaArrayType == 3) {
                blkIdx = 0;
                for (int j = 0; j < nCbS; j += pbOffset) {
                    for (int i = 0; i < nCbS; i += pbOffset) {
                        int x = x0 + i;
                        int y = y0 + j;
                        int modeIdx;
                        if (intraChromaPredMode[blkIdx] == 0)
                            modeIdx = cu->getIntraPredModeY(x, y) == 0 ? 34 : 0;
                        else if (intraChromaPredMode[blkIdx] == 1)
                            modeIdx = cu->getIntraPredModeY(x, y) == 26 ? 34 : 26;
                        else if (intraChromaPredMode[blkIdx] == 2)
                            modeIdx = cu->getIntraPredModeY(x, y) == 10 ? 34 : 10;
                        else if (intraChromaPredMode[blkIdx] == 3)
                            modeIdx = cu->getIntraPredModeY(x, y) == 1 ? 34 : 1;
                        else if (intraChromaPredMode[blkIdx] == 4)
                            modeIdx = cu->getIntraPredModeY(x, y);

                        cu->setIntraPredModeC(blkIdx, (IntraPredMode)modeIdx);
                        blkIdx++;
                    }
                }
            } else if (sps->chromaArrayType != 0) {
                int x = x0, y = y0;
                int modeIdx;
                if (intraChromaPredMode[0] == 0)
                    modeIdx = cu->getIntraPredModeY(x, y) == 0 ? 34 : 0;
                else if (intraChromaPredMode[0] == 1)
                    modeIdx = cu->getIntraPredModeY(x, y) == 26 ? 34 : 26;
                else if (intraChromaPredMode[0] == 2)
                    modeIdx = cu->getIntraPredModeY(x, y) == 10 ? 34 : 10;
                else if (intraChromaPredMode[0] == 3)
                    modeIdx = cu->getIntraPredModeY(x, y) == 1 ? 34 : 1;
                else if (intraChromaPredMode[0] == 4)
                    modeIdx = cu->getIntraPredModeY(x, y);

                if (sps->chromaArrayType == 2) {
                    const static int tabModeIdx[] = {
                        0,  1,  2,  2,  2,  2,  3,  5,  7,  8,  10, 12, 13, 15, 17, 18, 19, 20,
                        21, 22, 23, 23, 24, 24, 25, 25, 26, 27, 27, 28, 28, 29, 29, 30, 31,
                    };
                    cu->setIntraPredModeC(0, (IntraPredMode)tabModeIdx[modeIdx]);
                } else
                    cu->setIntraPredModeC(0, (IntraPredMode)modeIdx);
            }

        } else {
            switch (cu->getPartMode()) {
                case PART_MODE_2Nx2N:
                    decodePredictionUnit(ctu, cu, x0, y0, nCbS, nCbS, 0);
                    break;
                case PART_MODE_2NxN:
                    decodePredictionUnit(ctu, cu, x0, y0, nCbS, nCbS / 2, 0);
                    decodePredictionUnit(ctu, cu, x0, y0 + nCbS / 2, nCbS, nCbS / 2, 1);
                    break;
                case PART_MODE_Nx2N:
                    decodePredictionUnit(ctu, cu, x0, y0, nCbS / 2, nCbS, 0);
                    decodePredictionUnit(ctu, cu, x0 + nCbS / 2, y0, nCbS / 2, nCbS, 1);
                    break;
                case PART_MODE_2NxnU:
                    decodePredictionUnit(ctu, cu, x0, y0, nCbS, nCbS / 4, 0);
                    decodePredictionUnit(ctu, cu, x0, y0 + nCbS / 4, nCbS, nCbS * 3 / 4, 1);
                    break;
                case PART_MODE_2NxnD:
                    decodePredictionUnit(ctu, cu, x0, y0, nCbS, nCbS * 3 / 4, 0);
                    decodePredictionUnit(ctu, cu, x0, y0 + nCbS * 3 / 4, nCbS, nCbS / 4, 1);
                    break;
                case PART_MODE_nLx2N:
                    decodePredictionUnit(ctu, cu, x0, y0, nCbS / 4, nCbS, 0);
                    decodePredictionUnit(ctu, cu, x0 + nCbS / 4, y0, nCbS * 3 / 4, nCbS, 1);
                    break;
                case PART_MODE_nRx2N:
                    decodePredictionUnit(ctu, cu, x0, y0, nCbS * 3 / 4, nCbS, 0);
                    decodePredictionUnit(ctu, cu, x0 + nCbS * 3 / 4, y0, nCbS / 4, nCbS, 1);
                    break;
                case PART_MODE_NxN:
                    decodePredictionUnit(ctu, cu, x0, y0, nCbS / 2, nCbS / 2, 0);
                    decodePredictionUnit(ctu, cu, x0 + nCbS / 2, y0, nCbS / 2, nCbS / 2, 1);
                    decodePredictionUnit(ctu, cu, x0, y0 + nCbS / 2, nCbS / 2, nCbS / 2, 2);
                    decodePredictionUnit(ctu, cu, x0 + nCbS / 2, y0 + nCbS / 2, nCbS / 2, nCbS / 2, 3);
                    break;
                default:
                    break;
            }
        }

        if (!cu->getPcmFlag()) {
            bool rqtRootCbf = true;

            if (cu->getPredMode() != PRED_MODE_INTRA &&
                !(cu->getPartMode() == PART_MODE_2Nx2N && cu->getMergeFlag(x0, y0))) {
                rqtRootCbf = mCabacReader->parseRqtRootCbf();
            }

            if (rqtRootCbf) {
                decodeTransformTree(ctu, cu, x0, y0, x0, y0, log2CuSize, 0, 0);
            }
        }
    }

    // prediction
    if (cu->getPredMode() == PRED_MODE_INTRA) {
        reconstructIntraCu(ctu, cu, 0);
        if (sps->chromaArrayType != 0) {
            reconstructIntraCu(ctu, cu, 1);
            reconstructIntraCu(ctu, cu, 2);
        }
    } else {  // Inter
        reconstructInterCu(ctu, cu, 0);
        if (sps->chromaArrayType != 0) {
            reconstructInterCu(ctu, cu, 1);
            reconstructInterCu(ctu, cu, 2);
        }
    }

    frame->setCtDepth(x0, y0, log2CuSize, cu->getDepth());
    indexInCtu++;
}

void HevcCuDecoder::decodePredictionUnit(
    shared_ptr<HevcCtu> ctu, shared_ptr<HevcCu> cu, int x0, int y0, int nPbW, int nPbH, int partIdx) {
    auto frame = ctu->getFrame();
    auto sps = frame->getSps();
    auto pps = frame->getPps();
    auto sliceHeader = ctu->getSlice()->getSliceHeader();
    auto pu = std::make_shared<HevcPu>(x0, y0, nPbW, nPbH);
    MvField mvField;

    auto decodeMvdCoding = [this](Mv &mvd) {
        uint8_t absMvdGreater0Flag[2] = {0};
        uint8_t absMvdGreater1Flag[2] = {0};
        uint8_t mvdSignFlag[2] = {0};
        int absMvdMinus2[2] = {-1, -1};

        absMvdGreater0Flag[0] = mCabacReader->parseAbsMvdGreater0Flag();
        absMvdGreater0Flag[1] = mCabacReader->parseAbsMvdGreater0Flag();

        if (absMvdGreater0Flag[0]) absMvdGreater1Flag[0] = mCabacReader->parseAbsMvdGreater1Flag();
        if (absMvdGreater0Flag[1]) absMvdGreater1Flag[1] = mCabacReader->parseAbsMvdGreater1Flag();

        if (absMvdGreater0Flag[0]) {
            if (absMvdGreater1Flag[0]) absMvdMinus2[0] = mCabacReader->parseAbsMvdMinus2();
            mvdSignFlag[0] = mCabacReader->parseMvdSignFlag();
        }
        if (absMvdGreater0Flag[1]) {
            if (absMvdGreater1Flag[1]) absMvdMinus2[1] = mCabacReader->parseAbsMvdMinus2();
            mvdSignFlag[1] = mCabacReader->parseMvdSignFlag();
        }

        mvd.x = absMvdGreater0Flag[0] * (absMvdMinus2[0] + 2);
        if (mvdSignFlag[0]) mvd.x = -mvd.x;

        mvd.y = absMvdGreater0Flag[1] * (absMvdMinus2[1] + 2);
        if (mvdSignFlag[1]) mvd.y = -mvd.y;
    };

    uint8_t mergeIdx = 0;
    uint8_t mergeFlag = 0;
    InterPredIdc predIdc = INTER_PRED_L0;
    uint8_t refIdxLX[2] = {0};
    uint8_t mvpLXFlag[2] = {0};
    Mv mvdLX[2];

    if (cu->getCuSkipFlag()) {
        mergeFlag = 1;
        if (sliceHeader->maxNumMergeCand > 1) mergeIdx = mCabacReader->parseMergeIdx(sliceHeader->maxNumMergeCand - 1);
    } else {
        mergeFlag = mCabacReader->parseMergeFlag();
        if (mergeFlag) {
            if (sliceHeader->maxNumMergeCand > 1)
                mergeIdx = mCabacReader->parseMergeIdx(sliceHeader->maxNumMergeCand - 1);
        } else {
            if (sliceHeader->slice_type == B_SLICE)
                predIdc = (InterPredIdc)mCabacReader->parseInterPredIdc(nPbW, nPbH, cu->getDepth());

            if (predIdc != INTER_PRED_L1) {
                if (sliceHeader->num_ref_idx_l0_active_minus1 > 0)
                    refIdxLX[0] = mCabacReader->parseRefIdxLX(sliceHeader->num_ref_idx_l0_active_minus1);
                decodeMvdCoding(mvdLX[0]);
                mvpLXFlag[0] = mCabacReader->parseMvpLXFlag();
            }

            if (predIdc != INTER_PRED_L0) {
                if (sliceHeader->num_ref_idx_l1_active_minus1 > 0)
                    refIdxLX[1] = mCabacReader->parseRefIdxLX(sliceHeader->num_ref_idx_l1_active_minus1);
                if (sliceHeader->mvd_l1_zero_flag && predIdc == INTER_PRED_BI)
                    mvdLX[1].reset();
                else
                    decodeMvdCoding(mvdLX[1]);
                mvpLXFlag[1] = mCabacReader->parseMvpLXFlag();
            }
        }
    }

    pu->setMergeFlag(mergeFlag);

    shared_ptr<HevcFrame> colPic;
    if (sliceHeader->slice_type == B_SLICE && sliceHeader->collocated_from_l0_flag == 0)
        colPic = sliceHeader->refPicList[1][sliceHeader->collocated_ref_idx];
    else
        colPic = sliceHeader->refPicList[0][sliceHeader->collocated_ref_idx];

    auto scaleMv = [](int td, int tb, Mv mv) {
        Mv mvScaled;
        int tx = (16384 + (abs(td) >> 1)) / td;
        int distScaleFactor = clamp((tb * tx + 32) >> 6, -4096, 4095);
        {
            int value = distScaleFactor * mv.x;
            int sign = value < 0 ? -1 : 1;
            mvScaled.x = clamp(sign * ((abs(value) + 127) >> 8), -32768, 32767);
        }
        {
            int value = distScaleFactor * mv.y;
            int sign = value < 0 ? -1 : 1;
            mvScaled.y = clamp(sign * ((abs(value) + 127) >> 8), -32768, 32767);
        }

        return mvScaled;
    };

    auto getTemporalMvp = [scaleMv, sps, sliceHeader, colPic](int xPb, int yPb, int nPbW, int nPbH, int refIdxLX, int X,
                                                              Mv &mvLX) -> int {
        auto getCollocatedMv = [scaleMv, sliceHeader, colPic](int xColPb, int yColPb, int refIdxLX, int X,
                                                              Mv &mvLX) -> int {
            if (colPic->getCu(xColPb, yColPb)->getPredMode() == PRED_MODE_INTRA) {
                mvLX.reset();
                return 0;
            }

            MvField mvFieldCol;
            auto puCol = colPic->getPu(xColPb, yColPb);
            puCol->getMvField(mvFieldCol);

            int refIdxCol, listCol;
            Mv mvCol;
            if (mvFieldCol.predFlagLX[0] == 0) {
                mvCol = mvFieldCol.mvLX[1];
                refIdxCol = mvFieldCol.refIdxLX[1];
                listCol = 1;
            } else if (mvFieldCol.predFlagLX[0] == 1 && mvFieldCol.predFlagLX[1] == 0) {
                mvCol = mvFieldCol.mvLX[0];
                refIdxCol = mvFieldCol.refIdxLX[0];
                listCol = 0;
            } else {
                int noBackwardPredFlag = 1;
                for (int i = 0; (i < 2) && noBackwardPredFlag; i++) {
                    for (auto refPic : sliceHeader->refPicList[i]) {
                        if (refPic->getPoc() > sliceHeader->poc) {
                            noBackwardPredFlag = 0;
                            break;
                        }
                    }
                }

                if (noBackwardPredFlag) {
                    mvCol = mvFieldCol.mvLX[X];
                    refIdxCol = mvFieldCol.refIdxLX[X];
                    listCol = X;
                } else {
                    int i = sliceHeader->collocated_from_l0_flag;
                    mvCol = mvFieldCol.mvLX[i];
                    refIdxCol = mvFieldCol.refIdxLX[i];
                    listCol = i;
                }
            }

            int colPocDiff = colPic->getPoc() - colPic->mRefPocList[listCol][refIdxCol];
            int currPocDiff = sliceHeader->poc - sliceHeader->refPicList[X][refIdxLX]->getPoc();
            if (colPocDiff == currPocDiff)
                mvLX = mvCol;
            else {
                int td = clamp(colPocDiff, -128, 127);
                int tb = clamp(currPocDiff, -128, 127);
                mvLX = scaleMv(td, tb, mvCol);
            }

            return 1;
        };

        if (!sliceHeader->slice_temporal_mvp_enabled_flag) {
            mvLX.reset();
            return 0;
        }

        int available = 0;
        int xColBr = xPb + nPbW;
        int yColBr = yPb + nPbH;

        if ((yPb >> sps->log2CtbSize) == (yColBr >> sps->log2CtbSize) && yColBr < sps->height && xColBr < sps->width) {
            int xColPb = (xColBr >> 4) << 4;
            int yColPb = (yColBr >> 4) << 4;
            available = getCollocatedMv(xColPb, yColPb, refIdxLX, X, mvLX);
        } else {
            mvLX.reset();
            available = 0;
        }

        if (!available) {
            int xColCtr = xPb + (nPbW >> 1);
            int yColCtr = yPb + (nPbH >> 1);
            int xColPb = (xColCtr >> 4) << 4;
            int yColPb = (yColCtr >> 4) << 4;
            available = getCollocatedMv(xColPb, yColPb, refIdxLX, X, mvLX);
        }

        return available;
    };

    auto getMvp = [this, scaleMv, frame](int xCb, int yCb, int nCbS, int xPb, int yPb, int nPbW, int nPbH, int refIdxLX,
                                         int X, int partIdx, Mv &mvLXA, Mv &mvLXB, int &availableLXA,
                                         int &availableLXB) {
        int isScaled = 0;

        // A0 A1
        {
            availableLXA = 0;
            mvLXA.reset();
            int xNbA0 = xPb - 1, yNbA0 = yPb + nPbH;
            int xNbA1 = xNbA0, yNbA1 = yNbA0 - 1;
            bool availableA0 = checkNbPbAvail(frame, xCb, yCb, nCbS, xPb, yPb, nPbW, nPbH, partIdx, xNbA0, yNbA0);
            bool availableA1 = checkNbPbAvail(frame, xCb, yCb, nCbS, xPb, yPb, nPbW, nPbH, partIdx, xNbA1, yNbA1);
            int xNbTab[] = {xNbA0, xNbA1};
            int yNbTab[] = {yNbA0, yNbA1};
            bool availableTab[] = {availableA0, availableA1};

            if (availableA0 || availableA1) isScaled = 1;

            for (int i = 0; i < 2; i++) {
                if (availableTab[i] && !availableLXA) {
                    auto puNb = frame->getPu(xNbTab[i], yNbTab[i]);
                    MvField mvFieldNb;
                    puNb->getMvField(mvFieldNb);
                    if (mvFieldNb.predFlagLX[X] == 1 && (frame->mRefPicList[X][mvFieldNb.refIdxLX[X]]->getPoc() ==
                                                         frame->mRefPicList[X][refIdxLX]->getPoc())) {
                        availableLXA = 1;
                        mvLXA = mvFieldNb.mvLX[X];
                    } else if (mvFieldNb.predFlagLX[1 - X] == 1 &&
                               (frame->mRefPicList[1 - X][mvFieldNb.refIdxLX[1 - X]]->getPoc() ==
                                frame->mRefPicList[X][refIdxLX]->getPoc())) {
                        availableLXA = 1;
                        mvLXA = mvFieldNb.mvLX[1 - X];
                    }
                }
            }

            if (!availableLXA) {
                for (int i = 0; i < 2 && !availableLXA; i++) {
                    int refIdxA, refListA;
                    if (availableTab[i] && !availableLXA) {
                        auto puNb = frame->getPu(xNbTab[i], yNbTab[i]);
                        MvField mvFieldNb;
                        puNb->getMvField(mvFieldNb);
                        if (mvFieldNb.predFlagLX[X] == 1) {
                            availableLXA = 1;
                            mvLXA = mvFieldNb.mvLX[X];
                            refIdxA = mvFieldNb.refIdxLX[X];
                            refListA = X;
                        } else if (mvFieldNb.predFlagLX[1 - X] == 1) {
                            availableLXA = 1;
                            mvLXA = mvFieldNb.mvLX[1 - X];
                            refIdxA = mvFieldNb.refIdxLX[1 - X];
                            refListA = 1 - X;
                        }
                    }

                    if (availableLXA && (frame->mRefPicList[refListA][refIdxA]->getPoc() !=
                                         frame->mRefPicList[X][refIdxLX]->getPoc())) {
                        int td = clamp(frame->getPoc() - frame->mRefPicList[refListA][refIdxA]->getPoc(), -128, 127);
                        int tb = clamp(frame->getPoc() - frame->mRefPicList[X][refIdxLX]->getPoc(), -128, 127);
                        mvLXA = scaleMv(td, tb, mvLXA);
                    }
                }
            }
        }

        // B0 B1 B2
        {
            availableLXB = 0;
            mvLXB.reset();
            int xNbB0 = xPb + nPbW, yNbB0 = yPb - 1;
            int xNbB1 = xPb + nPbW - 1, yNbB1 = yPb - 1;
            int xNbB2 = xPb - 1, yNbB2 = yPb - 1;
            bool availableB0 = checkNbPbAvail(frame, xCb, yCb, nCbS, xPb, yPb, nPbW, nPbH, partIdx, xNbB0, yNbB0);
            bool availableB1 = checkNbPbAvail(frame, xCb, yCb, nCbS, xPb, yPb, nPbW, nPbH, partIdx, xNbB1, yNbB1);
            bool availableB2 = checkNbPbAvail(frame, xCb, yCb, nCbS, xPb, yPb, nPbW, nPbH, partIdx, xNbB2, yNbB2);
            int xNbTab[] = {xNbB0, xNbB1, xNbB2};
            int yNbTab[] = {yNbB0, yNbB1, yNbB2};
            bool availableTab[] = {availableB0, availableB1, availableB2};

            for (int i = 0; i < 3; i++) {
                if (availableTab[i] && !availableLXB) {
                    auto puNb = frame->getPu(xNbTab[i], yNbTab[i]);
                    MvField mvFieldNb;
                    puNb->getMvField(mvFieldNb);
                    if (mvFieldNb.predFlagLX[X] == 1 && (frame->mRefPicList[X][mvFieldNb.refIdxLX[X]]->getPoc() ==
                                                         frame->mRefPicList[X][refIdxLX]->getPoc())) {
                        availableLXB = 1;
                        mvLXB = mvFieldNb.mvLX[X];
                    } else if (mvFieldNb.predFlagLX[1 - X] == 1 &&
                               (frame->mRefPicList[1 - X][mvFieldNb.refIdxLX[1 - X]]->getPoc() ==
                                frame->mRefPicList[X][refIdxLX]->getPoc())) {
                        availableLXB = 1;
                        mvLXB = mvFieldNb.mvLX[1 - X];
                    }
                }
            }

            if (!isScaled && availableLXB) {
                availableLXA = 1;
                mvLXA = mvLXB;
            }

            if (!isScaled) {
                availableLXB = 0;
                for (int i = 0; i < 3 && !availableLXB; i++) {
                    int refIdxB, refListB;
                    if (availableTab[i] && !availableLXB) {
                        auto puNb = frame->getPu(xNbTab[i], yNbTab[i]);
                        MvField mvFieldNb;
                        puNb->getMvField(mvFieldNb);
                        if (mvFieldNb.predFlagLX[X] == 1) {
                            availableLXB = 1;
                            mvLXB = mvFieldNb.mvLX[X];
                            refIdxB = mvFieldNb.refIdxLX[X];
                            refListB = X;
                        } else if (mvFieldNb.predFlagLX[1 - X] == 1) {
                            availableLXB = 1;
                            mvLXB = mvFieldNb.mvLX[1 - X];
                            refIdxB = mvFieldNb.refIdxLX[1 - X];
                            refListB = 1 - X;
                        }
                    }

                    if (availableLXB && (frame->mRefPicList[refListB][refIdxB]->getPoc() !=
                                         frame->mRefPicList[X][refIdxLX]->getPoc())) {
                        int td = clamp(frame->getPoc() - frame->mRefPicList[refListB][refIdxB]->getPoc(), -128, 127);
                        int tb = clamp(frame->getPoc() - frame->mRefPicList[X][refIdxLX]->getPoc(), -128, 127);
                        mvLXB = scaleMv(td, tb, mvLXB);
                    }
                }
            }
        }
    };

    if (mergeFlag) {
        int xCb = cu->getX(), yCb = cu->getY();
        int xPb = x0, yPb = y0;
        int xOrigP = x0, yOrigP = y0;
        int nOrigPbW = nPbW, nOrigPbH = nPbH;
        int nCbS = 1 << cu->getLog2Size();
        int log2ParMrgLevel = pps->log2ParMrgLevel;
        auto partMode = cu->getPartMode();
        vector<MvField> mergeCandList;
        MvField mvFieldA0, mvFieldA1, mvFieldB0, mvFieldB1, mvFieldB2, mvFieldCol;
        bool availableA0, availableA1, availableB0, availableB1, availableB2;

        if (log2ParMrgLevel > 2 && nCbS == 8) {
            xPb = xCb;
            yPb = yCb;
            nPbW = nPbH = nCbS;
            partIdx = 0;
        }

        // A1
        {
            int xNb = xPb - 1, yNb = yPb + nPbH - 1;
            availableA1 = checkNbPbAvail(frame, xCb, yCb, nCbS, xPb, yPb, nPbW, nPbH, partIdx, xNb, yNb);
            if (((xPb >> log2ParMrgLevel) == (xNb >> log2ParMrgLevel) &&
                 (yPb >> log2ParMrgLevel) == (yNb >> log2ParMrgLevel)) ||
                ((partMode == PART_MODE_Nx2N || partMode == PART_MODE_nLx2N || partMode == PART_MODE_nRx2N) &&
                 partIdx == 1))
                availableA1 = false;
            if (availableA1) {
                auto puNb = frame->getPu(xNb, yNb);
                puNb->getMvField(mvFieldA1);
                mvFieldA1.availableFlag = 1;
            }
            if (mvFieldA1.availableFlag) mergeCandList.push_back(mvFieldA1);
        }
        // B1
        if (mergeCandList.size() < sliceHeader->maxNumMergeCand) {
            int xNb = xPb + nPbW - 1, yNb = yPb - 1;
            availableB1 = checkNbPbAvail(frame, xCb, yCb, nCbS, xPb, yPb, nPbW, nPbH, partIdx, xNb, yNb);
            if (((xPb >> log2ParMrgLevel) == (xNb >> log2ParMrgLevel) &&
                 (yPb >> log2ParMrgLevel) == (yNb >> log2ParMrgLevel)) ||
                ((partMode == PART_MODE_2NxN || partMode == PART_MODE_2NxnU || partMode == PART_MODE_2NxnD) &&
                 partIdx == 1))
                availableB1 = false;
            if (availableB1) {
                auto puNb = frame->getPu(xNb, yNb);
                puNb->getMvField(mvFieldB1);
                mvFieldB1.availableFlag = 1;
                if (availableA1 && mvFieldA1 == mvFieldB1) mvFieldB1.availableFlag = 0;
            }
            if (mvFieldB1.availableFlag) mergeCandList.push_back(mvFieldB1);
        }
        // B0
        if (mergeCandList.size() < sliceHeader->maxNumMergeCand) {
            int xNb = xPb + nPbW, yNb = yPb - 1;
            availableB0 = checkNbPbAvail(frame, xCb, yCb, nCbS, xPb, yPb, nPbW, nPbH, partIdx, xNb, yNb);
            if ((xPb >> log2ParMrgLevel) == (xNb >> log2ParMrgLevel) &&
                (yPb >> log2ParMrgLevel) == (yNb >> log2ParMrgLevel))
                availableB0 = false;
            if (availableB0) {
                auto puNb = frame->getPu(xNb, yNb);
                puNb->getMvField(mvFieldB0);
                mvFieldB0.availableFlag = 1;
                if (availableB1 && mvFieldB1 == mvFieldB0) mvFieldB0.availableFlag = 0;
            }
            if (mvFieldB0.availableFlag) mergeCandList.push_back(mvFieldB0);
        }
        // A0
        if (mergeCandList.size() < sliceHeader->maxNumMergeCand) {
            int xNb = xPb - 1, yNb = yPb + nPbH;
            availableA0 = checkNbPbAvail(frame, xCb, yCb, nCbS, xPb, yPb, nPbW, nPbH, partIdx, xNb, yNb);
            if ((xPb >> log2ParMrgLevel) == (xNb >> log2ParMrgLevel) &&
                (yPb >> log2ParMrgLevel) == (yNb >> log2ParMrgLevel))
                availableA0 = false;
            if (availableA0) {
                auto puNb = frame->getPu(xNb, yNb);
                puNb->getMvField(mvFieldA0);
                mvFieldA0.availableFlag = 1;
                if (availableA1 && mvFieldA1 == mvFieldA0) mvFieldA0.availableFlag = 0;
            }
            if (mvFieldA0.availableFlag) mergeCandList.push_back(mvFieldA0);
        }
        // B2
        if (mergeCandList.size() < sliceHeader->maxNumMergeCand) {
            int xNb = xPb - 1, yNb = yPb - 1;
            availableB2 = checkNbPbAvail(frame, xCb, yCb, nCbS, xPb, yPb, nPbW, nPbH, partIdx, xNb, yNb);
            if ((xPb >> log2ParMrgLevel) == (xNb >> log2ParMrgLevel) &&
                (yPb >> log2ParMrgLevel) == (yNb >> log2ParMrgLevel))
                availableB2 = false;
            if (availableB2) {
                auto puNb = frame->getPu(xNb, yNb);
                puNb->getMvField(mvFieldB2);
                mvFieldB2.availableFlag = 1;
                if ((availableA1 && mvFieldA1 == mvFieldB2) || (availableB1 && mvFieldB1 == mvFieldB2) ||
                    (mvFieldA0.availableFlag && mvFieldA1.availableFlag && mvFieldB0.availableFlag &&
                     mvFieldB1.availableFlag))
                    mvFieldB2.availableFlag = 0;
            }
            if (mvFieldB2.availableFlag) mergeCandList.push_back(mvFieldB2);
        }
        // Col
        if (mergeCandList.size() < sliceHeader->maxNumMergeCand) {
            mvFieldCol.refIdxLX[0] = 0;
            int availableFlagL0 = getTemporalMvp(xPb, yPb, nPbW, nPbH, mvFieldCol.refIdxLX[0], 0, mvFieldCol.mvLX[0]);
            mvFieldCol.availableFlag = availableFlagL0;
            mvFieldCol.predFlagLX[0] = availableFlagL0;
            mvFieldCol.predFlagLX[1] = 0;
            if (sliceHeader->slice_type == B_SLICE) {
                mvFieldCol.refIdxLX[1] = 0;
                int availableFlagL1 =
                    getTemporalMvp(xPb, yPb, nPbW, nPbH, mvFieldCol.refIdxLX[1], 1, mvFieldCol.mvLX[1]);
                mvFieldCol.availableFlag = (availableFlagL0 || availableFlagL1);
                mvFieldCol.predFlagLX[1] = availableFlagL1;
            }
            if (mvFieldCol.availableFlag) mergeCandList.push_back(mvFieldCol);
        }

        // combined bi-predictive merge candidates (applies for B slices)
        if (sliceHeader->slice_type == B_SLICE && mergeCandList.size() > 1 &&
            mergeCandList.size() < sliceHeader->maxNumMergeCand) {
            const static int l0CandIdxTab[] = {
                0, 1, 0, 2, 1, 2, 0, 3, 1, 3, 2, 3,
            };
            const static int l1CandIdxTab[] = {
                1, 0, 2, 0, 2, 1, 3, 0, 3, 1, 3, 2,
            };

            int numCurrMergeCand = mergeCandList.size();
            int numOrigMergeCand = mergeCandList.size();
            int numInputMergeCand = mergeCandList.size();
            int combIdx = 0;
            while (true) {
                int l0CandIdx = l0CandIdxTab[combIdx], l1CandIdx = l1CandIdxTab[combIdx];
                auto l0Cand = mergeCandList[l0CandIdx], l1Cand = mergeCandList[l1CandIdx];

                if (l0Cand.predFlagLX[0] == 1 && l1Cand.predFlagLX[1] == 1) {
                    int diffPoc = sliceHeader->refPicList[0][l0Cand.refIdxLX[0]]->getPoc() -
                                  sliceHeader->refPicList[1][l1Cand.refIdxLX[1]]->getPoc();
                    if (diffPoc || l0Cand.mvLX[0] != l1Cand.mvLX[1]) {
                        // int k = numCurrMergeCand - numInputMergeCand;
                        // auto &combCand = mergeCandList[k];
                        // mergeCandList.push_back(combCand);
                        MvField combCand;

                        combCand.refIdxLX[0] = l0Cand.refIdxLX[0];
                        combCand.refIdxLX[1] = l1Cand.refIdxLX[1];
                        combCand.predFlagLX[0] = 1;
                        combCand.predFlagLX[1] = 1;
                        combCand.mvLX[0] = l0Cand.mvLX[0];
                        combCand.mvLX[1] = l1Cand.mvLX[1];
                        mergeCandList.push_back(combCand);
                        numCurrMergeCand++;
                    }
                }
                combIdx++;

                if (combIdx == (numOrigMergeCand * (numOrigMergeCand - 1)) ||
                    numCurrMergeCand == sliceHeader->maxNumMergeCand)
                    break;
            }
        }

        // append zero motion vector candidates
        if (mergeCandList.size() < sliceHeader->maxNumMergeCand) {
            int numCurrMergeCand = mergeCandList.size();
            // int numInputMergeCand = mergeCandList.size();
            int zeroIdx = 0;
            int numRefIdx;
            if (sliceHeader->slice_type == P_SLICE)
                numRefIdx = sliceHeader->num_ref_idx_l0_active_minus1 + 1;
            else
                numRefIdx = std::min(sliceHeader->num_ref_idx_l0_active_minus1 + 1,
                                     sliceHeader->num_ref_idx_l1_active_minus1 + 1);

            while (true) {
                // int m = numCurrMergeCand - numInputMergeCand;
                // auto &zeroCand = mergeCandList[m];
                // mergeCandList.push_back(zeroCand);
                MvField zeroCand;

                zeroCand.refIdxLX[0] = (zeroIdx < numRefIdx) ? zeroIdx : 0;
                zeroCand.refIdxLX[1] = -1;
                zeroCand.predFlagLX[0] = 1;
                zeroCand.predFlagLX[1] = 0;
                zeroCand.mvLX[0].reset();
                zeroCand.mvLX[1].reset();
                if (sliceHeader->slice_type == B_SLICE) {
                    zeroCand.refIdxLX[1] = (zeroIdx < numRefIdx) ? zeroIdx : 0;
                    zeroCand.predFlagLX[1] = 1;
                }
                mergeCandList.push_back(zeroCand);
                numCurrMergeCand++;
                zeroIdx++;

                if (numCurrMergeCand >= sliceHeader->maxNumMergeCand) break;
            }
        }

        auto cand = mergeCandList[mergeIdx];
        mvField.refIdxLX[0] = cand.refIdxLX[0];
        mvField.refIdxLX[1] = cand.refIdxLX[1];
        mvField.predFlagLX[0] = cand.predFlagLX[0];
        mvField.predFlagLX[1] = cand.predFlagLX[1];
        if (sliceHeader->use_integer_mv_flag == 0) {
            mvField.mvLX[0] = cand.mvLX[0];
            mvField.mvLX[1] = cand.mvLX[1];
        } else {
            mvField.mvLX[0].x = cand.mvLX[0].x & 3;
            mvField.mvLX[0].y = cand.mvLX[0].y & 3;
            mvField.mvLX[1].x = cand.mvLX[1].x & 3;
            mvField.mvLX[1].y = cand.mvLX[1].y & 3;
        }

        if (mvField.predFlagLX[0] == 1 && mvField.predFlagLX[1] == 1 && (nOrigPbW + nOrigPbH) == 12) {
            mvField.refIdxLX[1] = -1;
            mvField.predFlagLX[1] = 0;
        }
    } else {
        int xCb = cu->getX(), yCb = cu->getY();
        int xPb = x0, yPb = y0;
        int nCbS = 1 << cu->getLog2Size();

        for (int i = 0; i < 2; i++) {
            if (predIdc == ((i == 0) ? INTER_PRED_L0 : INTER_PRED_L1) || predIdc == INTER_PRED_BI) {
                mvField.refIdxLX[i] = refIdxLX[i];
                mvField.predFlagLX[i] = 1;
            } else {
                mvField.refIdxLX[i] = -1;
                mvField.predFlagLX[i] = 0;
            }

            if (mvField.predFlagLX[i]) {
                int availableLXA, availableLXB, availableLXCol;
                Mv mvLXA, mvLXB, mvLXCol;

                getMvp(xCb, yCb, nCbS, xPb, yPb, nPbW, nPbH, mvField.refIdxLX[i], i, partIdx, mvLXA, mvLXB,
                       availableLXA, availableLXB);
                if (availableLXA && availableLXB && (mvLXA != mvLXB))
                    availableLXCol = 0;
                else
                    availableLXCol = getTemporalMvp(xPb, yPb, nPbW, nPbH, mvField.refIdxLX[i], i, mvLXCol);

                Mv mvpListLX[2];
                int k = 0;
                if (availableLXA) {
                    mvpListLX[k++] = mvLXA;
                    if (availableLXB && (mvLXA != mvLXB)) mvpListLX[k++] = mvLXB;
                } else if (availableLXB)
                    mvpListLX[k++] = mvLXB;
                if (k < 2 && availableLXCol) mvpListLX[k++] = mvLXCol;
                while (k < 2) mvpListLX[k++].reset();

                Mv mvpLX, uLX;
                mvpLX = mvpListLX[mvpLXFlag[i]];

                if (sliceHeader->use_integer_mv_flag) {
                    uLX.x = ((((mvpLX.x >> 2) + mvdLX[i].x) << 2) + 65536) & 65535;
                    uLX.y = ((((mvpLX.y >> 2) + mvdLX[i].y) << 2) + 65536) & 65535;
                } else {
                    uLX.x = (mvpLX.x + mvdLX[i].x + 65536) & 65535;
                    uLX.y = (mvpLX.y + mvdLX[i].y + 65536) & 65535;
                }
                mvField.mvLX[i].x = (uLX.x > 32768) ? (uLX.x - 65536) : uLX.x;
                mvField.mvLX[i].y = (uLX.y > 32768) ? (uLX.y - 65536) : uLX.y;
            }
        }
    }

    // TwoVersionsOfCurrDecPicFlag default0, skip

    // chroma
    if (sps->chromaArrayType != 0) {
        for (int i = 0; i < 2; i++) {
            if (mvField.predFlagLX[i]) {
                mvField.mvCLX[i].x = mvField.mvLX[i].x * 2 / sps->subWidthC;
                mvField.mvCLX[i].y = mvField.mvLX[i].y * 2 / sps->subHeightC;
            }
        }
    }

    pu->setMvField(mvField);
    cu->addPu(pu);
}

void HevcCuDecoder::decodeTransformTree(shared_ptr<HevcCtu> ctu,
                                        shared_ptr<HevcCu> cu,
                                        int x0,
                                        int y0,
                                        int xBase,
                                        int yBase,
                                        int log2TrafoSize,
                                        int trafoDepth,
                                        int blkIdx) {
    auto frame = ctu->getFrame();
    auto sps = frame->getSps();
    auto pps = frame->getPps();
    uint8_t intraSplitFlag = cu->getPredMode() == PRED_MODE_INTRA ? (cu->getPartMode() == PART_MODE_NxN ? 1 : 0) : 0;
    uint8_t maxTrafoDepth = cu->getPredMode() == PRED_MODE_INTRA
                                ? sps->max_transform_hierarchy_depth_intra + intraSplitFlag
                                : sps->max_transform_hierarchy_depth_inter;
    bool splitTransformFlag = false;

    uint8_t cbfCb[2] = {0, 0};
    uint8_t cbfCr[2] = {0, 0};

    auto tu = std::make_shared<HevcTu>(x0, y0, log2TrafoSize, trafoDepth);

    if (log2TrafoSize <= sps->log2MaxTbSize && log2TrafoSize > sps->log2MinTbSize && trafoDepth < maxTrafoDepth &&
        !(intraSplitFlag && (trafoDepth == 0)))
        splitTransformFlag = mCabacReader->parseSplitTransformFlag(5 - log2TrafoSize);
    else {
        bool interSplitFlag = sps->max_transform_hierarchy_depth_inter == 0 && cu->getPredMode() == PRED_MODE_INTER &&
                              cu->getPartMode() != PART_MODE_2Nx2N && trafoDepth == 0;
        splitTransformFlag =
            log2TrafoSize > sps->log2MaxTbSize || (intraSplitFlag && trafoDepth == 0) || interSplitFlag;
    }

    tu->setSplitTransformFlag(splitTransformFlag);

    if ((log2TrafoSize > 2 && sps->chromaArrayType != 0) || sps->chromaArrayType == 3) {
        if (trafoDepth == 0 || cu->getTu(xBase, yBase, trafoDepth - 1)->getCbfCb(0)) {
            cbfCb[0] = mCabacReader->parseCbfCbCr(trafoDepth);
            if (sps->chromaArrayType == 2 && (!splitTransformFlag || log2TrafoSize == 3))
                cbfCb[0] = mCabacReader->parseCbfCbCr(trafoDepth);
        }

        if (trafoDepth == 0 || cu->getTu(xBase, yBase, trafoDepth - 1)->getCbfCr(0)) {
            cbfCr[0] = mCabacReader->parseCbfCbCr(trafoDepth);
            if (sps->chromaArrayType == 2 && (!splitTransformFlag || log2TrafoSize == 3))
                cbfCr[0] = mCabacReader->parseCbfCbCr(trafoDepth);
        }
    }
    tu->setCbfCb(cbfCb);
    tu->setCbfCr(cbfCr);
    cu->addTu(tu);

    if (splitTransformFlag) {
        int x1 = x0 + (1 << (log2TrafoSize - 1));
        int y1 = y0 + (1 << (log2TrafoSize - 1));
        decodeTransformTree(ctu, cu, x0, y0, x0, y0, log2TrafoSize - 1, trafoDepth + 1, 0);
        decodeTransformTree(ctu, cu, x1, y0, x0, y0, log2TrafoSize - 1, trafoDepth + 1, 1);
        decodeTransformTree(ctu, cu, x0, y1, x0, y0, log2TrafoSize - 1, trafoDepth + 1, 2);
        decodeTransformTree(ctu, cu, x1, y1, x0, y0, log2TrafoSize - 1, trafoDepth + 1, 3);
    } else {
        uint8_t cbfLuma = 1;
        if (cu->getPredMode() == PRED_MODE_INTRA || trafoDepth != 0 || cbfCb[0] || cbfCr[0] ||
            (sps->chromaArrayType == 2 && (cbfCb[1] || cbfCr[1])))
            cbfLuma = mCabacReader->parseCbfLuma(trafoDepth == 0 ? 1 : 0);
        tu->setCbfLuma(cbfLuma);
        decodeTransformUnit(ctu, cu, tu, x0, y0, xBase, yBase, log2TrafoSize, trafoDepth, blkIdx);
    }
}

void HevcCuDecoder::decodeTransformUnit(shared_ptr<HevcCtu> ctu,
                                        shared_ptr<HevcCu> cu,
                                        shared_ptr<HevcTu> tu,
                                        int x0,
                                        int y0,
                                        int xBase,
                                        int yBase,
                                        int log2TrafoSize,
                                        int trafoDepth,
                                        int blkIdx) {
    auto frame = ctu->getFrame();
    auto sps = frame->getSps();
    auto pps = frame->getPps();
    auto sliceHeader = ctu->getSlice()->getSliceHeader();

    int log2TrafoSizeC = std::max<int>(2, log2TrafoSize - (sps->chromaArrayType == 3 ? 0 : 1));
    int cbfDepthC = trafoDepth - ((sps->chromaArrayType != 3 && log2TrafoSize == 2) ? 1 : 0);
    int xC = (sps->chromaArrayType != 3 && log2TrafoSize == 2) ? xBase : x0;
    int yC = (sps->chromaArrayType != 3 && log2TrafoSize == 2) ? yBase : y0;
    uint8_t cbfLuma = tu->getCbfLuma();
    auto tuC = cu->getTu(xC, yC, cbfDepthC);
    uint8_t cbfChroma =
        tuC->getCbfCb(0) || tuC->getCbfCr(0) || (sps->chromaArrayType == 2 && (tuC->getCbfCb(1) || tuC->getCbfCr(1)));

    if (cbfLuma || cbfChroma) {
        int xP = (x0 >> sps->log2MinCtbSize) << sps->log2MinCtbSize;
        int yP = (y0 >> sps->log2MinCtbSize) << sps->log2MinCtbSize;
        int nCbS = 1 << sps->log2MinCtbSize;
        // residual_adaptive_colour_transform_enabled_flag default 0
        // skip tu_residual_act_flag

        // if (pps->cu_qp_delta_enabled_flag && !ctu->getIsCuQpDeltaCoded(x0, y0)) {
        //     ctu->setIsCuQpDeltaCoded(x0, y0, 1);
        //     uint8_t cuQpDeltaAbs = mCabacReader->parseCuQpDeltaAbs();
        //     if (cuQpDeltaAbs) {
        //         uint8_t cuQpDeltaSignFlag = mCabacReader->parseCuQpDeltaSignFlag();
        //         int cuQpDeltaVal = cuQpDeltaAbs;
        //         if (cuQpDeltaSignFlag) cuQpDeltaVal = -cuQpDeltaVal;
        //         ctu->setCuQpDeltaVal(x0, y0, cuQpDeltaVal);
        //     }
        //     cu->setQpY(((cu->getQpYPred() + ctu->getCuQpDeltaVal(x0, y0) + 52 + 2 * sps->qpBdOffsetY) %
        //                 (52 + sps->qpBdOffsetY)) -
        //                sps->qpBdOffsetY);
        // }
        if (pps->cu_qp_delta_enabled_flag && !ctu->getIsCuQpDeltaCoded()) {
            ctu->setIsCuQpDeltaCoded(1);
            uint8_t cuQpDeltaAbs = mCabacReader->parseCuQpDeltaAbs();
            if (cuQpDeltaAbs) {
                uint8_t cuQpDeltaSignFlag = mCabacReader->parseCuQpDeltaSignFlag();
                int cuQpDeltaVal = cuQpDeltaAbs;
                if (cuQpDeltaSignFlag) cuQpDeltaVal = -cuQpDeltaVal;
                ctu->setCuQpDeltaVal(cuQpDeltaVal);
            }
            cu->setQpY(
                ((cu->getQpYPred() + ctu->getCuQpDeltaVal() + 52 + 2 * sps->qpBdOffsetY) % (52 + sps->qpBdOffsetY)) -
                sps->qpBdOffsetY);
        }

        // TODO: chroma_qp_offset -> skip

        // tu_residual_act_flag default 0
        if (sps->chromaArrayType != 0) {
            int qpiCb = clamp(
                cu->getQpY() + pps->pps_cb_qp_offset + sliceHeader->slice_cb_qp_offset /* + CuQpOffsetCb default0*/,
                -sps->qpBdOffsetC, 57);
            int qpiCr = clamp(
                cu->getQpY() + pps->pps_cr_qp_offset + sliceHeader->slice_cr_qp_offset /* + CuQpOffsetCr default0*/,
                -sps->qpBdOffsetC, 57);
            const static int qpiTable[] = {
                29, 30, 31, 32, 33, 33, 34, 34, 35, 35, 36, 36, 37, 37,
            };
            int qpCb, qpCr;

            if (sps->chromaArrayType == 1) {
                if (qpiCb < 30)
                    qpCb = qpiCb;
                else if (qpiCb >= 30 && qpiCb <= 43)
                    qpCb = qpiTable[qpiCb - 30];
                else
                    qpCb = qpiCb - 6;

                if (qpiCr < 30)
                    qpCr = qpiCr;
                else if (qpiCr >= 30 && qpiCr <= 43)
                    qpCr = qpiTable[qpiCr - 30];
                else
                    qpCr = qpiCr - 6;
            } else {
                qpCb = std::min(qpiCb, 51);
                qpCr = std::min(qpiCr, 51);
            }

            cu->setQpCb(qpCb);
            cu->setQpCr(qpCr);
        }

        if (cbfLuma) decodeResidualCoding(ctu, cu, tu, x0, y0, log2TrafoSize, 0);

        if (log2TrafoSize > 2 || sps->chromaArrayType == 3) {
            // TODO: cross_comp_pred -> skip
            if (tu->getCbfCb(0)) decodeResidualCoding(ctu, cu, tu, x0, y0, log2TrafoSizeC, 1);
            if (sps->chromaArrayType == 2 && tu->getCbfCb(1))
                decodeResidualCoding(ctu, cu, tu, x0, y0 + (1 << log2TrafoSizeC), log2TrafoSizeC, 1);

            if (tu->getCbfCr(0)) decodeResidualCoding(ctu, cu, tu, x0, y0, log2TrafoSizeC, 2);
            if (sps->chromaArrayType == 2 && tu->getCbfCr(1))
                decodeResidualCoding(ctu, cu, tu, x0, y0 + (1 << log2TrafoSizeC), log2TrafoSizeC, 2);
        } else if (blkIdx == 3) {
            auto tuBase = cu->getTu(xBase, yBase, trafoDepth - 1);
            if (tuBase->getCbfCb(0)) decodeResidualCoding(ctu, cu, tu, xBase, yBase, log2TrafoSize, 1);
            if (sps->chromaArrayType == 2 && tuBase->getCbfCb(1))
                decodeResidualCoding(ctu, cu, tu, xBase, yBase + (1 << log2TrafoSizeC), log2TrafoSize, 1);

            if (tuBase->getCbfCr(0)) decodeResidualCoding(ctu, cu, tu, xBase, yBase, log2TrafoSize, 2);
            if (sps->chromaArrayType == 2 && tuBase->getCbfCr(1))
                decodeResidualCoding(ctu, cu, tu, xBase, yBase + (1 << log2TrafoSizeC), log2TrafoSize, 2);
        }
    }

    cu->copyResidual(tu);
}

void HevcCuDecoder::decodeResidualCoding(shared_ptr<HevcCtu> ctu,
                                         shared_ptr<HevcCu> cu,
                                         shared_ptr<HevcTu> tu,
                                         int x0,
                                         int y0,
                                         int log2TrafoSize,
                                         int cIdx) {
    auto frame = ctu->getFrame();
    auto sps = frame->getSps();
    auto pps = frame->getPps();

    auto getScanPosition = [sps](int log2Size, int scanIdx, int scanPos) -> PointAddr {
        switch (scanIdx) {
            case 0: {
                switch (log2Size) {
                    case 0:
                        return sps->scanOrderDiagonal1x1[scanPos];
                    case 1:
                        return sps->scanOrderDiagonal2x2[scanPos];
                    case 2:
                        return sps->scanOrderDiagonal4x4[scanPos];
                    case 3:
                        return sps->scanOrderDiagonal8x8[scanPos];
                    default:
                        return {-1, -1};
                }
            } break;
            case 1: {
                switch (log2Size) {
                    case 0:
                        return sps->scanOrderHorizontal1x1[scanPos];
                    case 1:
                        return sps->scanOrderHorizontal2x2[scanPos];
                    case 2:
                        return sps->scanOrderHorizontal4x4[scanPos];
                    case 3:
                        return sps->scanOrderHorizontal8x8[scanPos];
                    default:
                        return {-1, -1};
                }
            } break;
            case 2: {
                switch (log2Size) {
                    case 0:
                        return sps->scanOrderVertical1x1[scanPos];
                    case 1:
                        return sps->scanOrderVertical2x2[scanPos];
                    case 2:
                        return sps->scanOrderVertical4x4[scanPos];
                    case 3:
                        return sps->scanOrderVertical8x8[scanPos];
                    default:
                        return {-1, -1};
                }
            } break;
            case 3: {
                switch (log2Size) {
                    case 2:
                        return sps->scanOrderTraverse4x4[scanPos];
                    case 3:
                        return sps->scanOrderTraverse8x8[scanPos];
                    case 4:
                        return sps->scanOrderTraverse16x16[scanPos];
                    case 5:
                        return sps->scanOrderTraverse32x32[scanPos];
                    default:
                        return {-1, -1};
                }
            } break;
            default:
                return {-1, -1};
        }
    };

    tu->setComponentLog2Size(cIdx, log2TrafoSize);

    int scanIdx = 0;
    if (cu->getPredMode() == PRED_MODE_INTRA) {
        if (log2TrafoSize == 2 || (log2TrafoSize == 3 && cIdx == 0) ||
            (log2TrafoSize == 3 && sps->chromaArrayType == 3)) {
            int predModeIntra = (cIdx == 0) ? cu->getIntraPredModeY(x0, y0) : cu->getIntraPredModeC(x0, y0);
            if (predModeIntra >= 6 && predModeIntra <= 14)
                scanIdx = 2;
            else if (predModeIntra >= 22 && predModeIntra <= 30)
                scanIdx = 1;
        }
    }

    if (pps->transform_skip_enabled_flag && !cu->getCuTransquantBypassFlag() &&
        (log2TrafoSize <= pps->log2MaxTransformSkipSize))
        tu->setTransformSkipFlag(cIdx, mCabacReader->parseTransformSkipFlag());
    else
        tu->setTransformSkipFlag(cIdx, false);

    // explicit_rdpcm_enabled_flag default 0
    // TODO: explicit_rdpcm_flag
    tu->setExplicitRdpcmFlag(cIdx, false);

    int lastSignificantCoeffX = 0, lastSignificantCoeffY = 0;
    uint8_t lastSigCoeffXPrefix = mCabacReader->parseLastSigCoeffXPrefix(cIdx, log2TrafoSize);
    uint8_t lastSigCoeffYPrefix = mCabacReader->parseLastSigCoeffYPrefix(cIdx, log2TrafoSize);
    if (lastSigCoeffXPrefix > 3) {
        uint8_t lastSigCoeffXSuffix = mCabacReader->parseLastSigCoeffXSuffix((lastSigCoeffXPrefix >> 1) - 1);
        lastSignificantCoeffX =
            (1 << ((lastSigCoeffXPrefix >> 1) - 1)) * (2 + (lastSigCoeffXPrefix & 1)) + lastSigCoeffXSuffix;
    } else
        lastSignificantCoeffX = lastSigCoeffXPrefix;

    if (lastSigCoeffYPrefix > 3) {
        uint8_t lastSigCoeffYSuffix = mCabacReader->parseLastSigCoeffYSuffix((lastSigCoeffYPrefix >> 1) - 1);
        lastSignificantCoeffY =
            (1 << ((lastSigCoeffYPrefix >> 1) - 1)) * (2 + (lastSigCoeffYPrefix & 1)) + lastSigCoeffYSuffix;
    } else
        lastSignificantCoeffY = lastSigCoeffYPrefix;

    if (scanIdx == 2) std::swap(lastSignificantCoeffX, lastSignificantCoeffY);

    int lastScanPos = 16;
    int lastSubBlock = (1 << (log2TrafoSize - 2)) * (1 << (log2TrafoSize - 2)) - 1;
    {
        int xC = 0, yC = 0;
        do {
            if (lastScanPos == 0) {
                lastScanPos = 16;
                lastSubBlock--;
            }
            lastScanPos--;
            auto scanPosition = getScanPosition(log2TrafoSize - 2, scanIdx, lastSubBlock);
            auto currPosition = getScanPosition(2, scanIdx, lastScanPos);
            xC = (scanPosition.x << 2) + currPosition.x;
            yC = (scanPosition.y << 2) + currPosition.y;
        } while ((xC != lastSignificantCoeffX) || (yC != lastSignificantCoeffY));
    }

    auto &transCoeff = (cIdx == 0) ? tu->mTransCoeffY : ((cIdx == 1) ? tu->mTransCoeffCb : tu->mTransCoeffCr);

    int greater1Ctx = 1;
    tu->resetCodedSubBlockFlag();
    for (int i = lastSubBlock; i >= 0; i--) {
        auto scanPosition = getScanPosition(log2TrafoSize - 2, scanIdx, i);
        int xS = scanPosition.x, yS = scanPosition.y;
        uint8_t inferSbDcSigCoeffFlag = 0;

        if (i < lastSubBlock && i > 0) {
            int csbfCtx = 0;
            if (xS < ((1 << (log2TrafoSize - 2)) - 1)) csbfCtx += tu->getCodedSubBlockFlag(xS + 1, yS);
            if (yS < ((1 << (log2TrafoSize - 2)) - 1)) csbfCtx += tu->getCodedSubBlockFlag(xS, yS + 1);

            int inc = 0;
            if (cIdx == 0)
                inc = std::min(csbfCtx, 1);
            else
                inc = 2 + std::min(csbfCtx, 1);
            tu->setCodedSubBlockFlag(xS, yS, mCabacReader->parseCodedSubBlockFlag(inc));

            inferSbDcSigCoeffFlag = 1;
        } else {
            if ((xS == 0 && yS == 0) || (xS == (lastSignificantCoeffX >> 2) && yS == (lastSignificantCoeffY >> 2)))
                tu->setCodedSubBlockFlag(xS, yS, 1);
            else
                tu->setCodedSubBlockFlag(xS, yS, 0);
        }

        uint8_t sigCoeffFlag[16] = {0};
        uint8_t coeffAbsLevelGreater1Flag[16] = {0};
        uint8_t coeffAbsLevelGreater2Flag[16] = {0};
        uint8_t coeffSignFlag[16] = {0};
        uint16_t coeffAbsLevelRemaining[16] = {0};
        uint8_t numNonZero = 0, signHidden = 0;

        for (int n = (i == lastSubBlock) ? lastScanPos : 15; n >= 0; n--) {
            auto currPosition = getScanPosition(2, scanIdx, n);
            int xC = (xS << 2) + currPosition.x, yC = (yS << 2) + currPosition.y;
            int xP = currPosition.x, yP = currPosition.y;
            if (xC == lastSignificantCoeffX && yC == lastSignificantCoeffY) {
                sigCoeffFlag[n] = 1;
                numNonZero++;
                continue;
            }

            if (tu->getCodedSubBlockFlag(xS, yS) && (n > 0 || !inferSbDcSigCoeffFlag)) {
                const static uint8_t ctxIdxMap[] = {
                    0, 1, 4, 5, 2, 3, 4, 5, 6, 6, 8, 8, 7, 7, 8,
                };
                int sigCtx = 0;
                // transform_skip_context_enabled_flag default 0
                if (log2TrafoSize == 2)
                    sigCtx = ctxIdxMap[(yC << 2) + xC];
                else if (xC + yC == 0)
                    sigCtx = 0;
                else {
                    int prevCsbf = 0;
                    if (xS < ((1 << (log2TrafoSize - 2)) - 1)) prevCsbf += tu->getCodedSubBlockFlag(xS + 1, yS);
                    if (yS < ((1 << (log2TrafoSize - 2)) - 1)) prevCsbf += (tu->getCodedSubBlockFlag(xS, yS + 1) << 1);

                    if (prevCsbf == 0)
                        sigCtx = (xP + yP == 0) ? 2 : (xP + yP < 3) ? 1 : 0;
                    else if (prevCsbf == 1)
                        sigCtx = (yP == 0) ? 2 : (yP == 1) ? 1 : 0;
                    else if (prevCsbf == 2)
                        sigCtx = (xP == 0) ? 2 : (xP == 1) ? 1 : 0;
                    else
                        sigCtx = 2;

                    if (cIdx == 0) {
                        if (xS + yS > 0) sigCtx += 3;
                        if (log2TrafoSize == 3)
                            sigCtx += (scanIdx == 0) ? 9 : 15;
                        else
                            sigCtx += 21;
                    } else {
                        if (log2TrafoSize == 3)
                            sigCtx += 9;
                        else
                            sigCtx += 12;
                    }
                }

                int inc = sigCtx;
                if (cIdx > 0) inc = sigCtx + 27;

                sigCoeffFlag[n] = mCabacReader->parseSigCoeffFlag(inc);

                if (sigCoeffFlag[n]) {
                    inferSbDcSigCoeffFlag = 0;
                    numNonZero++;
                }
            } else {
                if ((xC & 3) == 0 && (yC & 3) == 0 && inferSbDcSigCoeffFlag == 1 && tu->getCodedSubBlockFlag(xS, yS)) {
                    sigCoeffFlag[n] = 1;
                    numNonZero++;
                } else
                    sigCoeffFlag[n] = 0;
            }
        }

        if (numNonZero > 0) {
            int firstSigScanPos = 16;
            int lastSigScanPos = -1;
            int numGreater1Flag = 0;
            int lastGreater1ScanPos = -1;

            int ctxSet = (i > 0 && cIdx == 0) ? 2 : 0;
            if (!(i == lastSubBlock) && greater1Ctx == 0) ctxSet++;
            greater1Ctx = 1;

            for (int n = 15; n >= 0; n--) {
                if (sigCoeffFlag[n]) {
                    if (numGreater1Flag < 8) {
                        int inc = (ctxSet << 2) + greater1Ctx;
                        if (cIdx > 0) inc += 16;
                        coeffAbsLevelGreater1Flag[n] = mCabacReader->parseCoeffAbsLevelGreater1Flag(inc);
                        numGreater1Flag++;
                        if (coeffAbsLevelGreater1Flag[n]) {
                            greater1Ctx = 0;
                            if (lastGreater1ScanPos == -1) lastGreater1ScanPos = n;
                        } else if (greater1Ctx > 0 && greater1Ctx < 3)
                            greater1Ctx++;
                    }

                    if (lastSigScanPos == -1) lastSigScanPos = n;
                    firstSigScanPos = n;
                }
            }

            if (lastGreater1ScanPos != -1) {
                int inc = ctxSet;
                if (cIdx > 0) inc += 4;
                coeffAbsLevelGreater2Flag[lastGreater1ScanPos] = mCabacReader->parseCoeffAbsLevelGreater2Flag(inc);
            }

            // implicit_rdpcm_enabled_flag default 0
            if (cu->getCuTransquantBypassFlag() || tu->getExplicitRdpcmFlag(cIdx))
                signHidden = 0;
            else if ((lastSigScanPos - firstSigScanPos) > 3)
                signHidden = 1;

            for (int n = 15; n >= 0; n--) {
                if (sigCoeffFlag[n] && (!pps->sign_data_hiding_enabled_flag || !signHidden || (n != firstSigScanPos)))
                    coeffSignFlag[n] = mCabacReader->parseCoeffSignFlag();
            }

            int numSigCoeff = 0;
            int sumAbsLevel = 0;
            int cAbsLevel = 0, cRiceParam = 0;
            for (int n = 15; n >= 0; n--) {
                auto currPosition = getScanPosition(2, scanIdx, n);
                int xC = (xS << 2) + currPosition.x;
                int yC = (yS << 2) + currPosition.y;
                if (sigCoeffFlag[n]) {
                    int baseLevel = 1 + coeffAbsLevelGreater1Flag[n] + coeffAbsLevelGreater2Flag[n];
                    if (baseLevel == ((numSigCoeff < 8) ? ((n == lastGreater1ScanPos) ? 3 : 2) : 1)) {
                        cRiceParam = std::min(cRiceParam + (cAbsLevel > (3 * (1 << cRiceParam)) ? 1 : 0), 4);
                        coeffAbsLevelRemaining[n] = mCabacReader->parseCoeffAbsLevelRemaining(cRiceParam);
                        cAbsLevel = baseLevel + coeffAbsLevelRemaining[n];
                    }

                    int coeff = (coeffAbsLevelRemaining[n] + baseLevel) * (1 - 2 * coeffSignFlag[n]);
                    if (pps->sign_data_hiding_enabled_flag && signHidden) {
                        sumAbsLevel += (coeffAbsLevelRemaining[n] + baseLevel);
                        if ((n == firstSigScanPos) && (sumAbsLevel & 1)) coeff = -coeff;
                    }

                    transCoeff[yC][xC] = coeff;
                    numSigCoeff++;
                }
            }
        }
    }

    int qp;
    if (cIdx == 0)
        qp = clamp(cu->getQpY() + sps->qpBdOffsetY, 0, 51 + sps->qpBdOffsetY);
    else if (cIdx == 1)
        qp = cu->getQpCb() + sps->qpBdOffsetC;
    else
        qp = cu->getQpCr() + sps->qpBdOffsetC;

    int bitDepth = (cIdx == 0) ? sps->bitDepthY : sps->bitDepthC;
    int bdShift = std::max(20 - bitDepth, 0);
    int tsShift = 5 + log2TrafoSize;
    int rotateCoeffs = 0;
    int nTbS = 1 << log2TrafoSize;
    auto &residual = (cIdx == 0) ? tu->mResidualY : ((cIdx == 1) ? tu->mResidualCb : tu->mResidualCr);

    if (cu->getCuTransquantBypassFlag()) {
        for (int y = 0; y < nTbS; y++)
            for (int x = 0; x < nTbS; x++)
                residual[y][x] = rotateCoeffs ? transCoeff[nTbS - y - 1][nTbS - x - 1] : transCoeff[y][x];
    } else {
        vector<vector<int>> d(nTbS, vector<int>(nTbS));
        int coeffMin = -(1 << 15);
        int coeffMax = (1 << 15) - 1;

        {
            int log2TransformRange = 15;
            int bdShift = bitDepth + log2TrafoSize + 10 - log2TransformRange;
            const static int levelScale[] = {
                40, 45, 51, 57, 64, 72,
            };
            for (int y = 0; y < nTbS; y++)
                for (int x = 0; x < nTbS; x++)
                    d[y][x] = clamp(
                        ((transCoeff[y][x] * 16 * levelScale[qp % 6] << (qp / 6)) + (1 << (bdShift - 1))) >> bdShift,
                        coeffMin, coeffMax);
        }

        if (tu->getTransformSkipFlag(cIdx)) {
            for (int y = 0; y < nTbS; y++)
                for (int x = 0; x < nTbS; x++)
                    residual[y][x] = (((rotateCoeffs ? d[nTbS - y - 1][nTbS - x - 1] : d[y][x]) << tsShift) +
                                      (1 << (bdShift - 1))) >>
                                     bdShift;
        } else {
            // Intermediate variables
            vector<vector<int>> t(nTbS, vector<int>(nTbS));

            switch (nTbS) {
                case 4:
                    if (cIdx == 0 && cu->getPredMode() == PRED_MODE_INTRA) {
                        fastInverseDst(d, t, 7, coeffMin, coeffMax);
                        fastInverseDst(t, residual, bdShift, coeffMin, coeffMax);
                    } else {
                        partialButterflyInverse4(d, t, 7, coeffMin, coeffMax);
                        partialButterflyInverse4(t, residual, bdShift, coeffMin, coeffMax);
                    }
                    break;
                case 8:
                    partialButterflyInverse8(d, t, 7, coeffMin, coeffMax);
                    partialButterflyInverse8(t, residual, bdShift, coeffMin, coeffMax);
                    break;
                case 16:
                    partialButterflyInverse16(d, t, 7, coeffMin, coeffMax);
                    partialButterflyInverse16(t, residual, bdShift, coeffMin, coeffMax);
                    break;
                case 32:
                    partialButterflyInverse32(d, t, 7, coeffMin, coeffMax);
                    partialButterflyInverse32(t, residual, bdShift, coeffMin, coeffMax);
                    break;
                default:
                    break;
            }
        }
    }
}
