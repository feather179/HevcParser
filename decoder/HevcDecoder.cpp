#include "HevcDecoder.h"

#include "ParameterSet.h"
#include "HevcFrame.h"
#include "CavlcReader.h"
#include "HevcSliceDecoder.h"
#include "HevcFilter.h"

#include <vector>
#include <cstdio>
#include <cstdint>

#include <assert.h>

using std::vector;

// copy from ffmpeg libavformat/avc.c
static const uint8_t *findStartcodeInternal(const uint8_t *p, const uint8_t *end) {
    const uint8_t *a = p + 4 - ((intptr_t)p & 3);

    for (end -= 3; p < a && p < end; p++) {
        if (p[0] == 0 && p[1] == 0 && p[2] == 1) return p;
    }

    for (end -= 3; p < end; p += 4) {
        uint32_t x = *(const uint32_t *)p;
        //      if ((x - 0x01000100) & (~x) & 0x80008000) // little endian
        //      if ((x - 0x00010001) & (~x) & 0x00800080) // big endian
        if ((x - 0x01010101) & (~x) & 0x80808080) {  // generic
            if (p[1] == 0) {
                if (p[0] == 0 && p[2] == 1) return p;
                if (p[2] == 0 && p[3] == 1) return p + 1;
            }
            if (p[3] == 0) {
                if (p[2] == 0 && p[4] == 1) return p + 2;
                if (p[4] == 0 && p[5] == 1) return p + 3;
            }
        }
    }

    for (end += 3; p < end; p++) {
        if (p[0] == 0 && p[1] == 0 && p[2] == 1) return p;
    }

    return end + 3;
}

const uint8_t *findStartcode(const uint8_t *p, const uint8_t *end) {
    const uint8_t *out = findStartcodeInternal(p, end);
    if (p < out && out < end && !out[-1]) out--;
    return out;
}

HevcDecoder::HevcDecoder() {
    mSliceDecoder = std::make_shared<HevcSliceDecoder>();
    mFilter = std::make_shared<HevcFilter>();
    mParamSet = std::make_shared<ParameterSet>();
}

std::shared_ptr<HevcFrame> HevcDecoder::getNewFrame(uint8_t tid) {
    std::shared_ptr<HevcFrame> frame;

    auto sliceHeader = mParamSet->getSliceHeader();
    auto pps = mParamSet->getPps(sliceHeader->slice_pic_parameter_set_id);
    auto sps = mParamSet->getSps(pps->pps_seq_parameter_set_id);
    uint32_t dpbSize = sps->sps_max_dec_pic_buffering_minus1[tid] + 1;

    if (mFrames.size() < dpbSize) {
        frame = std::make_shared<HevcFrame>();
        frame->init(mParamSet);
        mFrames.push_back(frame);
    } else {
        int minPoc = INT_MAX;
        for (auto refFrame : mFrames) {
            if (!refFrame->getUsedForReference() && refFrame->getPoc() < minPoc) {
                frame = refFrame;
                minPoc = refFrame->getPoc();
            }
        }
        frame->init(mParamSet);
    }

    return frame;
}

void HevcDecoder::extractRbsp(vector<uint8_t> &nalu) {
    int zeroCount = 0;
    vector<uint8_t>::iterator iterRead, iterWrite;
    uint32_t position = 0;

    mEmulationPreventionByteLocation.clear();
    for (iterRead = iterWrite = nalu.begin(); iterRead != nalu.end(); iterRead++, iterWrite++, position++) {
        if (zeroCount == 2 && *iterRead == 3) {
            mEmulationPreventionByteLocation.push_back(position);
            position++;
            iterRead++;
            zeroCount = 0;
            if (iterRead == nalu.end()) break;
        }

        zeroCount = (*iterRead == 0) ? zeroCount + 1 : 0;
        *iterWrite = *iterRead;
    }
    nalu.resize(iterWrite - nalu.begin());
}

void HevcDecoder::decodePacket(vector<uint8_t> &buffer) {
    if (mParamSet == nullptr) mParamSet = std::make_shared<ParameterSet>();
    if (buffer.size() == 0) return;

    const uint8_t *pStart = buffer.data();
    const uint8_t *pEnd = buffer.data() + buffer.size();

    const uint8_t *pNaluStart = findStartcode(pStart, pEnd);
    const uint8_t *pNaluEnd = nullptr;
    for (;;) {
        while (pNaluStart && (pNaluStart < pEnd) && !*(pNaluStart++)) {}
        if (pNaluStart == pEnd) break;
        pNaluEnd = findStartcode(pNaluStart, pEnd);

        vector<uint8_t> nalu;
        nalu.insert(nalu.end(), pNaluStart, pNaluEnd);
        pNaluStart = pNaluEnd;

        decodeNalu(nalu);
    }
}

void HevcDecoder::decodeNalu(vector<uint8_t> &nalu) {
    // get RBSP
    extractRbsp(nalu);

    auto cavlcReader = std::make_unique<CavlcReader>(nalu);
    cavlcReader->parse(mParamSet, mPrevTid0Poc, mFrames);

    int naluType = cavlcReader->getNaluType();
    // slice segment
    if ((naluType >= NAL_UNIT_CODED_SLICE_TRAIL_N && naluType <= NAL_UNIT_CODED_SLICE_RASL_R) ||
        (naluType >= NAL_UNIT_CODED_SLICE_BLA_W_LP && naluType <= NAL_UNIT_CODED_SLICE_CRA)) {
        auto sliceHeader = mParamSet->getSliceHeader();

        // update prevTid0Poc
        {
            auto isSLNR = [](int naluType) -> bool {
                if (naluType == NAL_UNIT_CODED_SLICE_TRAIL_N || naluType == NAL_UNIT_CODED_SLICE_TSA_N ||
                    naluType == NAL_UNIT_CODED_SLICE_STSA_N || naluType == NAL_UNIT_CODED_SLICE_RADL_N ||
                    naluType == NAL_UNIT_RESERVED_VCL_N10 || naluType == NAL_UNIT_RESERVED_VCL_N12 ||
                    naluType == NAL_UNIT_RESERVED_VCL_N14)
                    return true;
                return false;
            };
            if (cavlcReader->getTid() == 0 && naluType != NAL_UNIT_CODED_SLICE_RASL_R &&
                naluType != NAL_UNIT_CODED_SLICE_RADL_R && !isSLNR(naluType))
                mPrevTid0Poc = sliceHeader->poc;
        }

        if (mFirstSliceInPicture || !mCurrFrame) {
            mCurrFrame = getNewFrame(cavlcReader->getTid());
            mCurrFrame->setPoc(sliceHeader->poc);
            mCurrFrame->setTid(cavlcReader->getTid());
            mCurrFrame->setRefPicList(sliceHeader->refPicList);
            // mFirstSliceInPicture = false;
        }

        uint32_t sliceHeaderSize = cavlcReader->numBytesRead();
        vector<uint32_t> subStreamSizes;
        vector<uint32_t> &entryPointOffsetMinus1 = sliceHeader->entry_point_offset_minus1;
        uint32_t currEntryPointOffset = 0;
        uint32_t prevEntryPointOffset = 0;
        for (int i = 0; i < entryPointOffsetMinus1.size(); i++) {
            uint32_t offset = entryPointOffsetMinus1[i] + 1;
            uint32_t count = 0;
            currEntryPointOffset += offset;
            for (int j = 0; j < mEmulationPreventionByteLocation.size(); j++) {
                if (mEmulationPreventionByteLocation[j] >= (prevEntryPointOffset + sliceHeaderSize) &&
                    mEmulationPreventionByteLocation[j] < (currEntryPointOffset + sliceHeaderSize))
                    count++;
            }
            prevEntryPointOffset = currEntryPointOffset;
            subStreamSizes.push_back(offset - count);
        }

        vector<vector<uint8_t>> subStreams(subStreamSizes.size() + 1);
        uint32_t offset = sliceHeaderSize;
        int i = 0;
        for (; i < subStreamSizes.size(); i++) {
            subStreams[i].insert(subStreams[i].end(), nalu.begin() + offset, nalu.begin() + offset + subStreamSizes[i]);
            offset += subStreamSizes[i];
        }
        subStreams[i].insert(subStreams[i].end(), nalu.begin() + offset, nalu.end());

        mSliceDecoder->decodeSlice(mCurrFrame, sliceHeader, subStreams);

        mFilter->deblockingFilter(mCurrFrame);
        mFilter->sao(mCurrFrame);

        if (mOutputFile) {
            int height = mCurrFrame->mSaoY.size();
            int width = mCurrFrame->mSaoY[0].size();
            uint8_t *buffer = (uint8_t *)malloc(width * height);
            int index = 0;
            for (int j = 0; j < height; j++)
                for (int i = 0; i < width; i++) buffer[index++] = mCurrFrame->mSaoY[j][i];
            fwrite(buffer, 1, width * height, mOutputFile);

            height = mCurrFrame->mSaoCb.size();
            width = mCurrFrame->mSaoCb[0].size();
            index = 0;
            for (int j = 0; j < height; j++)
                for (int i = 0; i < width; i++) buffer[index++] = mCurrFrame->mSaoCb[j][i];
            fwrite(buffer, 1, width * height, mOutputFile);

            index = 0;
            for (int j = 0; j < height; j++)
                for (int i = 0; i < width; i++) buffer[index++] = mCurrFrame->mSaoCr[j][i];
            fwrite(buffer, 1, width * height, mOutputFile);

            fflush(mOutputFile);
        }
    }
}

bool HevcDecoder::decode(const char *inputFilePath, const char *outputFilePath) {
    assert(inputFilePath != nullptr);
    assert(outputFilePath != nullptr);

    FILE *pInputFile = nullptr;
    if (!(pInputFile = fopen(inputFilePath, "rb"))) {
        printf("Can not open input file:%s\n", outputFilePath);
        return false;
    }

    if (!(mOutputFile = fopen(outputFilePath, "wb"))) {
        printf("Can not open output file:%s\n", outputFilePath);
        return false;
    }

#define BUFFER_SIZE (16 * 1024)

    uint8_t *pBuffer = (uint8_t *)malloc(BUFFER_SIZE);
    size_t readPosition = 0;
    vector<uint8_t> nalu;

    while (true) {
        nalu.clear();
        int readCount = 0;

    READ_PACKET:
        memset(pBuffer, 0, BUFFER_SIZE);
        fseek(pInputFile, readPosition, SEEK_SET);
        int readSize = fread(pBuffer, 1, BUFFER_SIZE, pInputFile);
        readCount++;
        if (readSize > 0) {
            const uint8_t *pStart = pBuffer;
            const uint8_t *pEnd = pBuffer + readSize;
            const uint8_t *pNaluStart = nullptr;
            const uint8_t *pNaluEnd = nullptr;
            int i = 0;

            if (readCount == 1) {
                pNaluStart = findStartcode(pStart, pEnd);
                while (pNaluStart && (pNaluStart < pEnd) && !*(pNaluStart++)) i++;
                i++;
                if (pNaluStart == pEnd) break;
                pNaluEnd = findStartcode(pNaluStart, pEnd);
            } else {
                pNaluStart = pStart;
                pNaluEnd = findStartcode(pNaluStart, pEnd);
            }

            if (pNaluEnd == pEnd && readSize == BUFFER_SIZE) {
                nalu.insert(nalu.end(), pNaluStart, pNaluEnd - 4);
                readPosition += (pNaluEnd - 4 - pNaluStart + i);
                goto READ_PACKET;
            } else {
                nalu.insert(nalu.end(), pNaluStart, pNaluEnd);
                readPosition += (pNaluEnd - pNaluStart + i);
            }
        } else
            break;

        printf("get nalu size:%zu\n", nalu.size());
        decodeNalu(nalu);
    }

    free(pBuffer);

    fclose(pInputFile);
    fclose(mOutputFile);

    return true;
}
