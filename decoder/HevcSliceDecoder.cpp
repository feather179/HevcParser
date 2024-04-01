#include "HevcSliceDecoder.h"

#include "HevcFrame.h"
#include "HevcSlice.h"
#include "HevcCtu.h"
#include "CabacReader.h"
#include "HevcCuDecoder.h"

using std::shared_ptr;
using std::vector;

uint32_t HevcSliceDecoder::getSubStreamIdx(uint32_t ctuRsAddrInSlice) {
    if (mPps->entropy_coding_sync_enabled_flag) return ctuRsAddrInSlice / mSps->ctbWidth;

    return 0;
}

void HevcSliceDecoder::decodeSlice(shared_ptr<HevcFrame> frame,
                                   shared_ptr<HevcSliceHeader> sliceHeader,
                                   vector<vector<uint8_t>> &subStreams) {
    auto sps = mSps = frame->getSps();
    auto pps = mPps = frame->getPps();

    vector<shared_ptr<CabacReader>> cabacReaders;
    for (int i = 0; i < subStreams.size(); i++) {
        auto reader = std::make_shared<CabacReader>(subStreams[i]);
        reader->reset(sliceHeader);
        cabacReaders.push_back(reader);
    }

    shared_ptr<CabacReader> cabacSyncContextState = std::make_shared<CabacReader>();

    uint32_t ctbRsAddr = sliceHeader->slice_segment_address;
    uint32_t ctbRsAddrStart = ctbRsAddr;
    uint32_t ctbTsAddr = pps->ctbAddrRsToTs[ctbRsAddr];
    bool endOfSliceSegmentFlag = false;

    auto cuDecoder = std::make_shared<HevcCuDecoder>();

    // loop over every ctu
    while (!endOfSliceSegmentFlag) {
        uint32_t idx = getSubStreamIdx(ctbRsAddr - ctbRsAddrStart);
        auto cabacReader = cabacReaders[idx];
        if (pps->entropy_coding_sync_enabled_flag && idx > 0)
            if (ctbRsAddr % sps->ctbWidth == 0) cabacReader->loadContexts(cabacSyncContextState);

        uint32_t ctuX = ctbRsAddr % sps->ctbWidth;
        uint32_t ctuY = ctbRsAddr / sps->ctbWidth;

        bool sliceSaoEnabled[3];
        uint8_t bitDepth[3];
        SaoParam saoParam[3];

        if (sliceHeader->slice_sao_luma_flag)
            sliceSaoEnabled[0] = true;
        else
            sliceSaoEnabled[0] = false;

        if (sliceHeader->slice_sao_chroma_flag)
            sliceSaoEnabled[1] = sliceSaoEnabled[2] = true;
        else
            sliceSaoEnabled[1] = sliceSaoEnabled[2] = false;

        bitDepth[0] = sps->bitDepthY;
        bitDepth[1] = bitDepth[2] = sps->bitDepthC;

        if (sliceHeader->slice_sao_luma_flag || sliceHeader->slice_sao_chroma_flag) {
            bool leftMergeAvail = false, aboveMergeAvail = false;

            if (ctuX > 0) {
                bool leftCtbInSliceSeg = (ctbRsAddr > ctbRsAddrStart);
                bool leftCtbInTile = (pps->ctbTileId[ctbTsAddr] == pps->ctbTileId[pps->ctbAddrRsToTs[ctbRsAddr - 1]]);
                if (leftCtbInSliceSeg && leftCtbInTile) leftMergeAvail = true;
            }
            if (ctuY > 0) {
                bool aboveCtbInSliceSeg = ((ctbRsAddr - sps->ctbWidth) >= ctbRsAddrStart);
                bool aboveCtbInTile =
                    (pps->ctbTileId[ctbTsAddr] == pps->ctbTileId[pps->ctbAddrRsToTs[ctbRsAddr - sps->ctbWidth]]);
                if (aboveCtbInSliceSeg && aboveCtbInTile) aboveMergeAvail = true;
            }
            cabacReader->parseSaoParam(saoParam, sliceSaoEnabled, leftMergeAvail, aboveMergeAvail, bitDepth);
        }

        auto ctu = frame->getCtu(ctbRsAddr);
        ctu->setSaoParam(saoParam);
        if (!sliceHeader->dependent_slice_segment_flag)
            ctu->setSliceRsAddr(sliceHeader->slice_segment_address);
        else {
            uint32_t addr = pps->ctbAddrTsToRs[pps->ctbAddrRsToTs[sliceHeader->slice_segment_address] - 1];
            ctu->setSliceRsAddr(frame->getCtu(addr)->getSliceRsAddr());
        }
        ctu->setTileId(pps->ctbTileId[ctbTsAddr]);
        cuDecoder->setCabacReader(cabacReader);
        cuDecoder->decodeCodingQuadtree(ctu, ctuX << sps->log2CtbSize, ctuY << sps->log2CtbSize, sps->log2CtbSize, 0);

        endOfSliceSegmentFlag = cabacReader->parseEndOfSliceSegmentFlag();
        ctbTsAddr++;

        if (ctbTsAddr >= sps->ctbCount) break;

        ctbRsAddr = pps->ctbAddrTsToRs[ctbTsAddr];

        if (pps->entropy_coding_sync_enabled_flag && (ctbTsAddr % sps->ctbWidth == 2))
            cabacSyncContextState->loadContexts(cabacReader);

        if (!endOfSliceSegmentFlag &&
            ((pps->tiles_enabled_flag && pps->ctbTileId[ctbTsAddr] != pps->ctbTileId[ctbTsAddr - 1]) ||
             (pps->entropy_coding_sync_enabled_flag &&
              (ctbRsAddr % sps->ctbWidth == 0 ||
               pps->ctbTileId[ctbTsAddr] != pps->ctbTileId[pps->ctbAddrRsToTs[ctbRsAddr - 1]])))) {
            // cabacReader->parseEndOfS
        }
    }
}
