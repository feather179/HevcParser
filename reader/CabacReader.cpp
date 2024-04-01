#include "CabacReader.h"

#include <algorithm>

using std::shared_ptr;
using std::vector;

// clang-format off
static const uint8_t gRangeTabLps[][4] = {
    {128, 176, 208, 240},
    {128, 167, 197, 227},
    {128, 158, 187, 216},
    {123, 150, 178, 205},
    {116, 142, 169, 195},
    {111, 135, 160, 185},
    {105, 128, 152, 175},
    {100, 122, 144, 166},
    {95,  116, 137, 158},
    {90,  110, 130, 150},
    {85,  104, 123, 142},
    {81,  99,  117, 135},
    {77,  94,  111, 128},
    {73,  89,  105, 122},
    {69,  85,  100, 116},
    {66,  80,  95,  110},
    {62,  76,  90,  104},
    {59,  72,  86,  99},
    {56,  69,  81,  94},
    {53,  65,  77,  89},
    {51,  62,  73,  85},
    {48,  59,  69,  80},
    {46,  56,  66,  76},
    {43,  53,  63,  72},
    {41,  50,  59,  69},
    {39,  48,  56,  65},
    {37,  45,  54,  62},
    {35,  43,  51,  59},
    {33,  41,  48,  56},
    {32,  39,  46,  53},
    {30,  37,  43,  50},
    {29,  35,  41,  48},
    {27,  33,  39,  45},
    {26,  31,  37,  43},
    {24,  30,  35,  41},
    {23,  28,  33,  39},
    {22,  27,  32,  37},
    {21,  26,  30,  35},
    {20,  24,  29,  33},
    {19,  23,  27,  31},
    {18,  22,  26,  30},
    {17,  21,  25,  28},
    {16,  20,  23,  27},
    {15,  19,  22,  25},
    {14,  18,  21,  24},
    {14,  17,  20,  23},
    {13,  16,  19,  22},
    {12,  15,  18,  21},
    {12,  14,  17,  20},
    {11,  14,  16,  19},
    {11,  13,  15,  18},
    {10,  12,  15,  17},
    {10,  12,  14,  16},
    {9,   11,  13,  15},
    {9,   11,  12,  14},
    {8,   10,  12,  14},
    {8,   9,   11,  13},
    {7,   9,   11,  12},
    {7,   9,   10,  12},
    {7,   8,   10,  11},
    {6,   8,   9,   11},
    {6,   7,   9,   10},
    {6,   7,   8,   9},
    {2,   2,   2,   2},
};

/**
 * Offset to ctxIdx 0 in init_values and states, indexed by SyntaxElement.
 */
static const int gElementOffset[] = {
    0, // sao_merge_flag
    1, // sao_type_idx
    2, // sao_eo_class
    2, // sao_band_position
    2, // sao_offset_abs
    2, // sao_offset_sign
    2, // end_of_slice_flag
    2, // split_coding_unit_flag
    5, // cu_transquant_bypass_flag
    6, // skip_flag
    9, // cu_qp_delta
    12, // pred_mode
    13, // part_mode
    17, // pcm_flag
    17, // prev_intra_luma_pred_mode
    18, // mpm_idx
    18, // rem_intra_luma_pred_mode
    18, // intra_chroma_pred_mode
    20, // merge_flag
    21, // merge_idx
    22, // inter_pred_idc
    27, // ref_idx_l0
    29, // ref_idx_l1
    31, // abs_mvd_greater0_flag
    33, // abs_mvd_greater1_flag
    35, // abs_mvd_minus2
    35, // mvd_sign_flag
    35, // mvp_lx_flag
    36, // no_residual_data_flag
    37, // split_transform_flag
    40, // cbf_luma
    42, // cbf_cb, cbf_cr
    47, // transform_skip_flag[][]
    49, // explicit_rdpcm_flag[][]
    51, // explicit_rdpcm_dir_flag[][]
    53, // last_significant_coeff_x_prefix
    71, // last_significant_coeff_y_prefix
    89, // last_significant_coeff_x_suffix
    89, // last_significant_coeff_y_suffix
    89, // significant_coeff_group_flag
    93, // significant_coeff_flag
    137, // coeff_abs_level_greater1_flag
    161, // coeff_abs_level_greater2_flag
    167, // coeff_abs_level_remaining
    167, // coeff_sign_flag
    167, // log2_res_scale_abs
    175, // res_scale_sign_flag
    177, // cu_chroma_qp_offset_flag
    178, // cu_chroma_qp_offset_idx
};

#define CNU 154
/**
 * Indexed by init_type
 */
static const uint8_t gInitValues[3][HEVC_CONTEXTS] = {
    { // sao_merge_flag
      153,
      // sao_type_idx
      200,
      // split_coding_unit_flag
      139, 141, 157,
      // cu_transquant_bypass_flag
      154,
      // skip_flag
      CNU, CNU, CNU,
      // cu_qp_delta
      154, 154, 154,
      // pred_mode
      CNU,
      // part_mode
      184, CNU, CNU, CNU,
      // prev_intra_luma_pred_mode
      184,
      // intra_chroma_pred_mode
      63, 139,
      // merge_flag
      CNU,
      // merge_idx
      CNU,
      // inter_pred_idc
      CNU, CNU, CNU, CNU, CNU,
      // ref_idx_l0
      CNU, CNU,
      // ref_idx_l1
      CNU, CNU,
      // abs_mvd_greater1_flag
      CNU, CNU,
      // abs_mvd_greater1_flag
      CNU, CNU,
      // mvp_lx_flag
      CNU,
      // no_residual_data_flag
      CNU,
      // split_transform_flag
      153, 138, 138,
      // cbf_luma
      111, 141,
      // cbf_cb, cbf_cr
      94, 138, 182, 154, 154,
      // transform_skip_flag
      139, 139,
      // explicit_rdpcm_flag
      139, 139,
      // explicit_rdpcm_dir_flag
      139, 139,
      // last_significant_coeff_x_prefix
      110, 110, 124, 125, 140, 153, 125, 127, 140, 109, 111, 143, 127, 111,
       79, 108, 123,  63,
      // last_significant_coeff_y_prefix
      110, 110, 124, 125, 140, 153, 125, 127, 140, 109, 111, 143, 127, 111,
       79, 108, 123,  63,
      // significant_coeff_group_flag
      91, 171, 134, 141,
      // significant_coeff_flag
      111, 111, 125, 110, 110,  94, 124, 108, 124, 107, 125, 141, 179, 153,
      125, 107, 125, 141, 179, 153, 125, 107, 125, 141, 179, 153, 125, 140,
      139, 182, 182, 152, 136, 152, 136, 153, 136, 139, 111, 136, 139, 111,
      141, 111,
      // coeff_abs_level_greater1_flag
      140,  92, 137, 138, 140, 152, 138, 139, 153,  74, 149,  92, 139, 107,
      122, 152, 140, 179, 166, 182, 140, 227, 122, 197,
      // coeff_abs_level_greater2_flag
      138, 153, 136, 167, 152, 152,
      // log2_res_scale_abs
      154, 154, 154, 154, 154, 154, 154, 154,
      // res_scale_sign_flag
      154, 154,
      // cu_chroma_qp_offset_flag
      154,
      // cu_chroma_qp_offset_idx
      154,
    },
    { // sao_merge_flag
      153,
      // sao_type_idx
      185,
      // split_coding_unit_flag
      107, 139, 126,
      // cu_transquant_bypass_flag
      154,
      // skip_flag
      197, 185, 201,
      // cu_qp_delta
      154, 154, 154,
      // pred_mode
      149,
      // part_mode
      154, 139, 154, 154,
      // prev_intra_luma_pred_mode
      154,
      // intra_chroma_pred_mode
      152, 139,
      // merge_flag
      110,
      // merge_idx
      122,
      // inter_pred_idc
      95, 79, 63, 31, 31,
      // ref_idx_l0
      153, 153,
      // ref_idx_l1
      153, 153,
      // abs_mvd_greater1_flag
      140, 198,
      // abs_mvd_greater1_flag
      140, 198,
      // mvp_lx_flag
      168,
      // no_residual_data_flag
      79,
      // split_transform_flag
      124, 138, 94,
      // cbf_luma
      153, 111,
      // cbf_cb, cbf_cr
      149, 107, 167, 154, 154,
      // transform_skip_flag
      139, 139,
      // explicit_rdpcm_flag
      139, 139,
      // explicit_rdpcm_dir_flag
      139, 139,
      // last_significant_coeff_x_prefix
      125, 110,  94, 110,  95,  79, 125, 111, 110,  78, 110, 111, 111,  95,
       94, 108, 123, 108,
      // last_significant_coeff_y_prefix
      125, 110,  94, 110,  95,  79, 125, 111, 110,  78, 110, 111, 111,  95,
       94, 108, 123, 108,
      // significant_coeff_group_flag
      121, 140, 61, 154,
      // significant_coeff_flag
      155, 154, 139, 153, 139, 123, 123,  63, 153, 166, 183, 140, 136, 153,
      154, 166, 183, 140, 136, 153, 154, 166, 183, 140, 136, 153, 154, 170,
      153, 123, 123, 107, 121, 107, 121, 167, 151, 183, 140, 151, 183, 140,
      140, 140,
      // coeff_abs_level_greater1_flag
      154, 196, 196, 167, 154, 152, 167, 182, 182, 134, 149, 136, 153, 121,
      136, 137, 169, 194, 166, 167, 154, 167, 137, 182,
      // coeff_abs_level_greater2_flag
      107, 167, 91, 122, 107, 167,
      // log2_res_scale_abs
      154, 154, 154, 154, 154, 154, 154, 154,
      // res_scale_sign_flag
      154, 154,
      // cu_chroma_qp_offset_flag
      154,
      // cu_chroma_qp_offset_idx
      154,
    },
    { // sao_merge_flag
      153,
      // sao_type_idx
      160,
      // split_coding_unit_flag
      107, 139, 126,
      // cu_transquant_bypass_flag
      154,
      // skip_flag
      197, 185, 201,
      // cu_qp_delta
      154, 154, 154,
      // pred_mode
      134,
      // part_mode
      154, 139, 154, 154,
      // prev_intra_luma_pred_mode
      183,
      // intra_chroma_pred_mode
      152, 139,
      // merge_flag
      154,
      // merge_idx
      137,
      // inter_pred_idc
      95, 79, 63, 31, 31,
      // ref_idx_l0
      153, 153,
      // ref_idx_l1
      153, 153,
      // abs_mvd_greater1_flag
      169, 198,
      // abs_mvd_greater1_flag
      169, 198,
      // mvp_lx_flag
      168,
      // no_residual_data_flag
      79,
      // split_transform_flag
      224, 167, 122,
      // cbf_luma
      153, 111,
      // cbf_cb, cbf_cr
      149, 92, 167, 154, 154,
      // transform_skip_flag
      139, 139,
      // explicit_rdpcm_flag
      139, 139,
      // explicit_rdpcm_dir_flag
      139, 139,
      // last_significant_coeff_x_prefix
      125, 110, 124, 110,  95,  94, 125, 111, 111,  79, 125, 126, 111, 111,
       79, 108, 123,  93,
      // last_significant_coeff_y_prefix
      125, 110, 124, 110,  95,  94, 125, 111, 111,  79, 125, 126, 111, 111,
       79, 108, 123,  93,
      // significant_coeff_group_flag
      121, 140, 61, 154,
      // significant_coeff_flag
      170, 154, 139, 153, 139, 123, 123,  63, 124, 166, 183, 140, 136, 153,
      154, 166, 183, 140, 136, 153, 154, 166, 183, 140, 136, 153, 154, 170,
      153, 138, 138, 122, 121, 122, 121, 167, 151, 183, 140, 151, 183, 140,
      140, 140,
      // coeff_abs_level_greater1_flag
      154, 196, 167, 167, 154, 152, 167, 182, 182, 134, 149, 136, 153, 121,
      136, 122, 169, 208, 166, 167, 154, 152, 167, 182,
      // coeff_abs_level_greater2_flag
      107, 167, 91, 107, 107, 167,
      // log2_res_scale_abs
      154, 154, 154, 154, 154, 154, 154, 154,
      // res_scale_sign_flag
      154, 154,
      // cu_chroma_qp_offset_flag
      154,
      // cu_chroma_qp_offset_idx
      154,
    },
};
// clang-format on

CabacReader::CabacReader() : BitReader(nullptr, 0) {}

CabacReader::CabacReader(const vector<uint8_t> &data) : BitReader(data.data(), data.size()) {
    mCurrRange = 510;
    mOffset = getUInt(9);
}

CabacReader::~CabacReader() {}

void CabacReader::reset(shared_ptr<HevcSliceHeader> sliceHeader) {
    uint8_t initType = 0;
    if (sliceHeader->slice_type == I_SLICE)
        initType = 0;
    else if (sliceHeader->slice_type == P_SLICE)
        initType = sliceHeader->cabac_init_flag ? 2 : 1;
    else
        initType = sliceHeader->cabac_init_flag ? 1 : 2;

    int32_t qp = sliceHeader->sliceQp;
    for (int i = 0; i < HEVC_CONTEXTS; i++) mContexts[i].init(qp, gInitValues[initType][i]);
}

void CabacReader::loadContexts(shared_ptr<CabacReader> cabacReader) {
    for (int i = 0; i < HEVC_CONTEXTS; i++) {
        mContexts[i].mStateIdx = cabacReader->mContexts[i].mStateIdx;
        mContexts[i].mMps = cabacReader->mContexts[i].mMps;
    }
}

uint8_t CabacReader::decodeBin(CabacContext &context) {
    uint8_t binVal = 0;
    uint8_t rangeIdx = (mCurrRange >> 6) & 3;
    uint8_t lpsRange = gRangeTabLps[context.getStateIdx()][rangeIdx];
    mCurrRange -= lpsRange;

    if (mOffset >= mCurrRange) {
        // LPS path
        binVal = context.getLps();
        mOffset -= mCurrRange;
        mCurrRange = lpsRange;
        context.updateLps();
    } else {
        // MPS path
        binVal = context.getMps();
        context.updateMps();
    }

    renormalization();

    return binVal;
}

uint8_t CabacReader::decodeBypass() {
    uint8_t binVal = 0;

    mOffset = (mOffset << 1) | getFlag();
    if (mOffset >= mCurrRange) {
        binVal = 1;
        mOffset -= mCurrRange;
    }

    return binVal;
}

uint8_t CabacReader::decodeTerminate() {
    uint8_t binVal = 0;

    mCurrRange -= 2;
    if (mOffset >= mCurrRange)
        binVal = 1;
    else
        renormalization();

    return binVal;
}

void CabacReader::renormalization() {
    while (mCurrRange < 256) {
        mCurrRange <<= 1;
        mOffset = (mOffset << 1) | getFlag();
    }
}

void CabacReader::parseSaoParam(
    SaoParam *saoParam, bool sliceSaoEnabled[3], bool leftMergeAvail, bool aboveMergeAvail, uint8_t bitDepth[3]) {
    auto parseSaoMerge = [this]() -> bool {
        return decodeBin(mContexts[gElementOffset[SAO_MERGE_FLAG]]) ? true : false;
    };
    auto parseSaoTypeIdx = [this]() -> uint8_t {
        uint8_t code = decodeBin(mContexts[gElementOffset[SAO_TYPE_IDX]]);
        if (code == 0)
            return 0;
        else {
            code = decodeBypass();
            if (code == 0) return 1;
            return 2;
        }
    };
    auto parseSaoMaxUvlc = [this](uint32_t maxValue) -> uint8_t {
        if (maxValue == 0) return 0;

        uint8_t value = 0;
        while (true) {
            uint8_t code = decodeBypass();
            if (code == 0) break;
            value++;
            if (value == maxValue) break;
        }
        return value;
    };
    auto parseSaoSign = [this]() -> uint8_t { return decodeBypass(); };
    auto parseSaoUflc = [this](uint8_t bits) -> uint8_t {
        uint8_t value = 0;
        for (uint8_t i = 0; i < bits; i++) value = ((value << 1) | decodeBypass());
        return value;
    };

    bool isLeftMerge = false, isAboveMerge = false;
    uint8_t saoTypeIdx[3];

    if (leftMergeAvail) isLeftMerge = parseSaoMerge();
    if (!isLeftMerge && aboveMergeAvail) isAboveMerge = parseSaoMerge();

    if (isLeftMerge || isAboveMerge) {
        for (int cIdx = 0; cIdx < 3; cIdx++) {
            if (sliceSaoEnabled[cIdx]) {
                if (isLeftMerge)
                    saoParam[cIdx].mode = SAO_MODE_MERGE_LEFT;
                else
                    saoParam[cIdx].mode = SAO_MODE_MERGE_ABOVE;
            } else
                saoParam[cIdx].mode = SAO_MODE_OFF;
        }
    } else {
        for (int cIdx = 0; cIdx < 3; cIdx++) {
            if (sliceSaoEnabled[cIdx]) {
                if (cIdx == 0 || cIdx == 1)
                    saoTypeIdx[cIdx] = parseSaoTypeIdx();
                else
                    saoTypeIdx[cIdx] = saoTypeIdx[1];

                if (saoTypeIdx[cIdx] != 0) {
                    uint32_t maxValue = (1 << (std::min<uint8_t>(bitDepth[cIdx], 10) - 5)) - 1;
                    for (int i = 0; i < 4; i++) saoParam[cIdx].offset[i] = parseSaoMaxUvlc(maxValue);

                    if (saoTypeIdx[cIdx] == 1) {
                        saoParam[cIdx].mode = SAO_MODE_BO;
                        for (int i = 0; i < 4; i++) {
                            if (saoParam[cIdx].offset[i])
                                if (parseSaoSign()) saoParam[cIdx].offset[i] = -saoParam[cIdx].offset[i];
                        }
                        saoParam[cIdx].bandPosition = parseSaoUflc(5);
                    } else {
                        saoParam[cIdx].offset[2] = -saoParam[cIdx].offset[2];
                        saoParam[cIdx].offset[3] = -saoParam[cIdx].offset[3];
                        if (cIdx == 0 || cIdx == 1) {
                            uint8_t mode = parseSaoUflc(2);
                            if (mode == 0)
                                saoParam[cIdx].mode = SAO_MODE_EO_0;
                            else if (mode == 1)
                                saoParam[cIdx].mode = SAO_MODE_EO_90;
                            else if (mode == 2)
                                saoParam[cIdx].mode = SAO_MODE_EO_135;
                            else if (mode == 3)
                                saoParam[cIdx].mode = SAO_MODE_EO_45;
                        } else
                            saoParam[cIdx].mode = saoParam[1].mode;
                    }
                } else
                    saoParam[cIdx].mode = SAO_MODE_OFF;
            } else
                saoParam[cIdx].mode = SAO_MODE_OFF;
        }
    }
}

bool CabacReader::parseSplitCuFlag(int inc) {
    return decodeBin(mContexts[gElementOffset[SPLIT_CODING_UNIT_FLAG] + inc]) ? true : false;
}

bool CabacReader::parseCuTransquantBypassFlag() {
    return decodeBin(mContexts[gElementOffset[CU_TRANSQUANT_BYPASS_FLAG]]) ? true : false;
}

bool CabacReader::parseCuSkipFlag(int inc) {
    return decodeBin(mContexts[gElementOffset[SKIP_FLAG] + inc]) ? true : false;
}

bool CabacReader::parsePredModeFlag() {
    return decodeBin(mContexts[gElementOffset[PRED_MODE_FLAG]]) ? true : false;
}

PartMode CabacReader::parsePartMode(int log2CbSize, int minLog2CbSize, bool ampEnabledFlag, PredMode predMode) {
    if (decodeBin(mContexts[gElementOffset[PART_MODE]]))  // 1
        return PART_MODE_2Nx2N;
    if (log2CbSize == minLog2CbSize) {
        if (predMode == PRED_MODE_INTRA)  // 0
            return PART_MODE_NxN;
        if (decodeBin(mContexts[gElementOffset[PART_MODE] + 1]))  // 01
            return PART_MODE_2NxN;
        if (log2CbSize == 3)  // 00
            return PART_MODE_Nx2N;
        if (decodeBin(mContexts[gElementOffset[PART_MODE] + 2]))  // 001
            return PART_MODE_Nx2N;
        return PART_MODE_NxN;
    }

    if (!ampEnabledFlag) {
        if (decodeBin(mContexts[gElementOffset[PART_MODE] + 1]))  // 01
            return PART_MODE_2NxN;
        return PART_MODE_Nx2N;  // 00
    }

    if (decodeBin(mContexts[gElementOffset[PART_MODE] + 1])) {    // 01x
        if (decodeBin(mContexts[gElementOffset[PART_MODE] + 3]))  // 011
            return PART_MODE_2NxN;
        if (decodeBypass())  // 0101
            return PART_MODE_2NxnD;
        return PART_MODE_2NxnU;  // 0100
    }

    if (decodeBin(mContexts[gElementOffset[PART_MODE] + 3]))  // 001
        return PART_MODE_Nx2N;
    if (decodeBypass())  // 0001
        return PART_MODE_nRx2N;

    return PART_MODE_nLx2N;  // 0000
}

uint8_t CabacReader::parsePrevIntraLumaPredFlag() {
    return decodeBin(mContexts[gElementOffset[PREV_INTRA_LUMA_PRED_FLAG]]);
}

uint8_t CabacReader::parseMpmIdx() {
    if (decodeBypass()) {
        if (decodeBypass()) return 2;
        return 1;
    }
    return 0;
}

uint8_t CabacReader::parseRemIntraLumaPredMode() {
    uint8_t value = 0;
    for (int i = 0; i < 5; i++) value = ((value << 1) | decodeBypass());
    return value;
}

uint8_t CabacReader::parseIntraChromaPredMode() {
    if (!decodeBin(mContexts[gElementOffset[INTRA_CHROMA_PRED_MODE]])) return 4;

    return ((decodeBypass() << 1) | decodeBypass());
}

uint8_t CabacReader::parseMergeIdx(int maxValue) {
    uint8_t i = decodeBin(mContexts[gElementOffset[MERGE_IDX]]);
    if (i != 0)
        while (i < maxValue && decodeBypass()) i++;

    return i;
}

uint8_t CabacReader::parseMergeFlag() {
    return decodeBin(mContexts[gElementOffset[MERGE_FLAG]]);
}

uint8_t CabacReader::parseInterPredIdc(int nPbW, int nPbH, int ctDepth) {
    if (nPbW + nPbH == 12) return decodeBin(mContexts[gElementOffset[INTER_PRED_IDC] + 4]);
    if (decodeBin(mContexts[gElementOffset[INTER_PRED_IDC] + ctDepth])) return 2;
    return decodeBin(mContexts[gElementOffset[INTER_PRED_IDC] + 4]);
}

uint8_t CabacReader::parseRefIdxLX(int maxValue) {
    uint8_t i = 0;

    while (i < 2 && i < maxValue && decodeBin(mContexts[gElementOffset[REF_IDX_L0] + i])) i++;
    if (i == 2)
        while (i < maxValue && decodeBypass()) i++;

    return i;
}

uint8_t CabacReader::parseMvpLXFlag() {
    return decodeBin(mContexts[gElementOffset[MVP_LX_FLAG]]);
}

uint8_t CabacReader::parseAbsMvdGreater0Flag() {
    return decodeBin(mContexts[gElementOffset[ABS_MVD_GREATER0_FLAG]]);
}

uint8_t CabacReader::parseAbsMvdGreater1Flag() {
    return decodeBin(mContexts[gElementOffset[ABS_MVD_GREATER1_FLAG] + 1]);
}

uint16_t CabacReader::parseAbsMvdMinus2() {
    auto parseEGK = [this](int k) -> uint16_t {
        int numNonZero = 0;
        uint16_t value = 0;

        while (decodeBypass()) numNonZero++;

        for (int i = 0; i < (numNonZero + k); i++) value = ((value << 1) | decodeBypass());

        return value + (((1 << numNonZero) - 1) << k);
    };

    return parseEGK(1);
}

uint8_t CabacReader::parseMvdSignFlag() {
    return decodeBypass();
}

bool CabacReader::parseRqtRootCbf() {
    return decodeBin(mContexts[gElementOffset[NO_RESIDUAL_DATA_FLAG]]) ? true : false;
}

bool CabacReader::parseSplitTransformFlag(int inc) {
    return decodeBin(mContexts[gElementOffset[SPLIT_TRANSFORM_FLAG] + inc]) ? true : false;
}

uint8_t CabacReader::parseCbfCbCr(int inc) {
    return decodeBin(mContexts[gElementOffset[CBF_CB_CR] + inc]);
}

uint8_t CabacReader::parseCbfLuma(int inc) {
    return decodeBin(mContexts[gElementOffset[CBF_LUMA] + inc]);
}

uint8_t CabacReader::parseCuQpDeltaAbs() {
    uint8_t prefix = 0, suffix = 0;
    int inc = 0;

    auto parseEGK = [this](int k) -> uint8_t {
        int numNonZero = 0;
        uint8_t value = 0;

        while (decodeBypass()) numNonZero++;

        for (int i = 0; i < (numNonZero + k); i++) value = ((value << 1) | decodeBypass());

        return value + (((1 << numNonZero) - 1) << k);
    };

    while (prefix < 5 && decodeBin(mContexts[gElementOffset[CU_QP_DELTA] + inc])) {
        prefix++;
        inc = 1;
    }

    if (prefix > 4) suffix = parseEGK(0);

    return prefix + suffix;
}

uint8_t CabacReader::parseCuQpDeltaSignFlag() {
    return decodeBypass();
}

bool CabacReader::parseTransformSkipFlag() {
    return decodeBin(mContexts[gElementOffset[TRANSFORM_SKIP_FLAG]]) ? true : false;
}

bool CabacReader::parseExplicitRdpcmFlag() {
    return decodeBin(mContexts[gElementOffset[EXPLICIT_RDPCM_FLAG]]) ? true : false;
}

bool CabacReader::parseExplicitRdpcmDirFlag() {
    return decodeBin(mContexts[gElementOffset[EXPLICIT_RDPCM_DIR_FLAG]]) ? true : false;
}

uint8_t CabacReader::parseLastSigCoeffXPrefix(int cIdx, int log2TrafoSize) {
    uint8_t maxValue = (log2TrafoSize << 1) - 1;
    uint8_t ctxOffset = 15;
    uint8_t ctxShift = log2TrafoSize - 2;
    if (cIdx == 0) {
        ctxOffset = 3 * (log2TrafoSize - 2) + ((log2TrafoSize - 1) >> 2);
        ctxShift = (log2TrafoSize + 1) >> 2;
    }

    uint8_t i = 0;
    while (i < maxValue &&
           decodeBin(mContexts[gElementOffset[LAST_SIGNIFICANT_COEFF_X_PREFIX] + (i >> ctxShift) + ctxOffset]))
        i++;

    return i;
}

uint8_t CabacReader::parseLastSigCoeffYPrefix(int cIdx, int log2TrafoSize) {
    uint8_t maxValue = (log2TrafoSize << 1) - 1;
    uint8_t ctxOffset = 15;
    uint8_t ctxShift = log2TrafoSize - 2;
    if (cIdx == 0) {
        ctxOffset = 3 * (log2TrafoSize - 2) + ((log2TrafoSize - 1) >> 2);
        ctxShift = (log2TrafoSize + 1) >> 2;
    }

    uint8_t i = 0;
    while (i < maxValue &&
           decodeBin(mContexts[gElementOffset[LAST_SIGNIFICANT_COEFF_Y_PREFIX] + (i >> ctxShift) + ctxOffset]))
        i++;

    return i;
}

uint8_t CabacReader::parseLastSigCoeffXSuffix(int length) {
    uint8_t value = 0;
    for (int i = 0; i < length; i++) value = ((value << 1) | decodeBypass());

    return value;
}

uint8_t CabacReader::parseLastSigCoeffYSuffix(int length) {
    uint8_t value = 0;
    for (int i = 0; i < length; i++) value = ((value << 1) | decodeBypass());

    return value;
}

uint8_t CabacReader::parseCodedSubBlockFlag(int inc) {
    return decodeBin(mContexts[gElementOffset[SIGNIFICANT_COEFF_GROUP_FLAG] + inc]);
}

uint8_t CabacReader::parseSigCoeffFlag(int inc) {
    return decodeBin(mContexts[gElementOffset[SIGNIFICANT_COEFF_FLAG] + inc]);
}

uint8_t CabacReader::parseCoeffAbsLevelGreater1Flag(int inc) {
    return decodeBin(mContexts[gElementOffset[COEFF_ABS_LEVEL_GREATER1_FLAG] + inc]);
}

uint8_t CabacReader::parseCoeffAbsLevelGreater2Flag(int inc) {
    return decodeBin(mContexts[gElementOffset[COEFF_ABS_LEVEL_GREATER2_FLAG] + inc]);
}

uint16_t CabacReader::parseCoeffAbsLevelRemaining(int riceParam) {
    uint16_t prefix = 0, suffix = 0;

    auto parseEGK = [this](int k) -> uint16_t {
        int numNonZero = 0;
        uint16_t value = 0;

        while (decodeBypass()) numNonZero++;

        for (int i = 0; i < (numNonZero + k); i++) value = ((value << 1) | decodeBypass());

        return value + (((1 << numNonZero) - 1) << k);
    };

    while (prefix < 4 && decodeBypass()) prefix++;

    if (prefix < 4) {
        for (int i = 0; i < riceParam; i++) suffix = ((suffix << 1) | decodeBypass());
        return (prefix << riceParam) + suffix;
    } else {
        suffix = parseEGK(riceParam + 1);
        return (4 << riceParam) + suffix;
    }
}

uint8_t CabacReader::parseCoeffSignFlag() {
    return decodeBypass();
}

bool CabacReader::parseEndOfSliceSegmentFlag() {
    return decodeTerminate() ? true : false;
}
