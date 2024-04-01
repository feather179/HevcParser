#ifndef CABAC_READER_H
#define CABAC_READER_H

#include "BitReader.h"
#include "CabacContext.h"
#include "Typedef.h"

#include <cstdint>
#include <vector>
#include <memory>

#define HEVC_CONTEXTS 199

enum SyntaxElement {
    SAO_MERGE_FLAG = 0,
    SAO_TYPE_IDX,
    SAO_EO_CLASS,
    SAO_BAND_POSITION,
    SAO_OFFSET_ABS,
    SAO_OFFSET_SIGN,
    END_OF_SLICE_FLAG,
    SPLIT_CODING_UNIT_FLAG,
    CU_TRANSQUANT_BYPASS_FLAG,
    SKIP_FLAG,
    CU_QP_DELTA,
    PRED_MODE_FLAG,
    PART_MODE,
    PCM_FLAG,
    PREV_INTRA_LUMA_PRED_FLAG,
    MPM_IDX,
    REM_INTRA_LUMA_PRED_MODE,
    INTRA_CHROMA_PRED_MODE,
    MERGE_FLAG,
    MERGE_IDX,
    INTER_PRED_IDC,
    REF_IDX_L0,
    REF_IDX_L1,
    ABS_MVD_GREATER0_FLAG,
    ABS_MVD_GREATER1_FLAG,
    ABS_MVD_MINUS2,
    MVD_SIGN_FLAG,
    MVP_LX_FLAG,
    NO_RESIDUAL_DATA_FLAG,
    SPLIT_TRANSFORM_FLAG,
    CBF_LUMA,
    CBF_CB_CR,
    TRANSFORM_SKIP_FLAG,
    EXPLICIT_RDPCM_FLAG,
    EXPLICIT_RDPCM_DIR_FLAG,
    LAST_SIGNIFICANT_COEFF_X_PREFIX,
    LAST_SIGNIFICANT_COEFF_Y_PREFIX,
    LAST_SIGNIFICANT_COEFF_X_SUFFIX,
    LAST_SIGNIFICANT_COEFF_Y_SUFFIX,
    SIGNIFICANT_COEFF_GROUP_FLAG,
    SIGNIFICANT_COEFF_FLAG,
    COEFF_ABS_LEVEL_GREATER1_FLAG,
    COEFF_ABS_LEVEL_GREATER2_FLAG,
    COEFF_ABS_LEVEL_REMAINING,
    COEFF_SIGN_FLAG,
    LOG2_RES_SCALE_ABS,
    RES_SCALE_SIGN_FLAG,
    CU_CHROMA_QP_OFFSET_FLAG,
    CU_CHROMA_QP_OFFSET_IDX,
};

class CabacReader : public BitReader {
private:
    uint16_t mCurrRange;
    uint16_t mOffset;
    CabacContext mContexts[HEVC_CONTEXTS];

    uint8_t decodeBin(CabacContext &context);
    uint8_t decodeBypass();
    uint8_t decodeTerminate();
    void renormalization();

public:
    CabacReader();
    CabacReader(const std::vector<uint8_t> &data);
    ~CabacReader();
    void reset(std::shared_ptr<HevcSliceHeader> sliceHeader);
    void loadContexts(std::shared_ptr<CabacReader> cabacReader);

    void parseSaoParam(
        SaoParam *saoParam, bool sliceSaoEnabled[3], bool leftMergeAvail, bool aboveMergeAvail, uint8_t bitDepth[3]);
    bool parseSplitCuFlag(int inc);
    bool parseCuTransquantBypassFlag();
    bool parseCuSkipFlag(int inc);
    bool parsePredModeFlag();
    PartMode parsePartMode(int log2CbSize, int minLog2CbSize, bool ampEnabledFlag, PredMode predMode);
    uint8_t parsePrevIntraLumaPredFlag();
    uint8_t parseMpmIdx();
    uint8_t parseRemIntraLumaPredMode();
    uint8_t parseIntraChromaPredMode();
    uint8_t parseMergeIdx(int maxValue);
    uint8_t parseMergeFlag();
    uint8_t parseInterPredIdc(int nPbW, int nPbH, int ctDepth);
    uint8_t parseRefIdxLX(int maxValue);
    uint8_t parseMvpLXFlag();
    uint8_t parseAbsMvdGreater0Flag();
    uint8_t parseAbsMvdGreater1Flag();
    uint16_t parseAbsMvdMinus2();
    uint8_t parseMvdSignFlag();
    bool parseRqtRootCbf();
    bool parseSplitTransformFlag(int inc);
    uint8_t parseCbfCbCr(int inc);
    uint8_t parseCbfLuma(int inc);
    uint8_t parseCuQpDeltaAbs();
    uint8_t parseCuQpDeltaSignFlag();
    bool parseTransformSkipFlag();
    bool parseExplicitRdpcmFlag();
    bool parseExplicitRdpcmDirFlag();
    uint8_t parseLastSigCoeffXPrefix(int cIdx, int log2TrafoSize);
    uint8_t parseLastSigCoeffYPrefix(int cIdx, int log2TrafoSize);
    uint8_t parseLastSigCoeffXSuffix(int length);
    uint8_t parseLastSigCoeffYSuffix(int length);
    uint8_t parseCodedSubBlockFlag(int inc);
    uint8_t parseSigCoeffFlag(int inc);
    uint8_t parseCoeffAbsLevelGreater1Flag(int inc);
    uint8_t parseCoeffAbsLevelGreater2Flag(int inc);
    uint16_t parseCoeffAbsLevelRemaining(int riceParam);
    uint8_t parseCoeffSignFlag();
    bool parseEndOfSliceSegmentFlag();
};

#endif
