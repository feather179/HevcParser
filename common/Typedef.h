#ifndef TYPEDEF_H
#define TYPEDEF_H

#include <cstdint>
#include <vector>
#include <memory>

enum NalUnitType {
    NAL_UNIT_CODED_SLICE_TRAIL_N = 0,  // 0
    NAL_UNIT_CODED_SLICE_TRAIL_R,      // 1

    NAL_UNIT_CODED_SLICE_TSA_N,  // 2
    NAL_UNIT_CODED_SLICE_TSA_R,  // 3

    NAL_UNIT_CODED_SLICE_STSA_N,  // 4
    NAL_UNIT_CODED_SLICE_STSA_R,  // 5

    NAL_UNIT_CODED_SLICE_RADL_N,  // 6
    NAL_UNIT_CODED_SLICE_RADL_R,  // 7

    NAL_UNIT_CODED_SLICE_RASL_N,  // 8
    NAL_UNIT_CODED_SLICE_RASL_R,  // 9

    NAL_UNIT_RESERVED_VCL_N10,
    NAL_UNIT_RESERVED_VCL_R11,
    NAL_UNIT_RESERVED_VCL_N12,
    NAL_UNIT_RESERVED_VCL_R13,
    NAL_UNIT_RESERVED_VCL_N14,
    NAL_UNIT_RESERVED_VCL_R15,

    NAL_UNIT_CODED_SLICE_BLA_W_LP,    // 16
    NAL_UNIT_CODED_SLICE_BLA_W_RADL,  // 17
    NAL_UNIT_CODED_SLICE_BLA_N_LP,    // 18
    NAL_UNIT_CODED_SLICE_IDR_W_RADL,  // 19
    NAL_UNIT_CODED_SLICE_IDR_N_LP,    // 20
    NAL_UNIT_CODED_SLICE_CRA,         // 21
    NAL_UNIT_RESERVED_IRAP_VCL22,
    NAL_UNIT_RESERVED_IRAP_VCL23,

    NAL_UNIT_RESERVED_VCL24,
    NAL_UNIT_RESERVED_VCL25,
    NAL_UNIT_RESERVED_VCL26,
    NAL_UNIT_RESERVED_VCL27,
    NAL_UNIT_RESERVED_VCL28,
    NAL_UNIT_RESERVED_VCL29,
    NAL_UNIT_RESERVED_VCL30,
    NAL_UNIT_RESERVED_VCL31,

    NAL_UNIT_VPS,                    // 32
    NAL_UNIT_SPS,                    // 33
    NAL_UNIT_PPS,                    // 34
    NAL_UNIT_ACCESS_UNIT_DELIMITER,  // 35
    NAL_UNIT_EOS,                    // 36
    NAL_UNIT_EOB,                    // 37
    NAL_UNIT_FILLER_DATA,            // 38
    NAL_UNIT_PREFIX_SEI,             // 39
    NAL_UNIT_SUFFIX_SEI,             // 40

    NAL_UNIT_RESERVED_NVCL41,
    NAL_UNIT_RESERVED_NVCL42,
    NAL_UNIT_RESERVED_NVCL43,
    NAL_UNIT_RESERVED_NVCL44,
    NAL_UNIT_RESERVED_NVCL45,
    NAL_UNIT_RESERVED_NVCL46,
    NAL_UNIT_RESERVED_NVCL47,
    NAL_UNIT_UNSPECIFIED_48,
    NAL_UNIT_UNSPECIFIED_49,
    NAL_UNIT_UNSPECIFIED_50,
    NAL_UNIT_UNSPECIFIED_51,
    NAL_UNIT_UNSPECIFIED_52,
    NAL_UNIT_UNSPECIFIED_53,
    NAL_UNIT_UNSPECIFIED_54,
    NAL_UNIT_UNSPECIFIED_55,
    NAL_UNIT_UNSPECIFIED_56,
    NAL_UNIT_UNSPECIFIED_57,
    NAL_UNIT_UNSPECIFIED_58,
    NAL_UNIT_UNSPECIFIED_59,
    NAL_UNIT_UNSPECIFIED_60,
    NAL_UNIT_UNSPECIFIED_61,
    NAL_UNIT_UNSPECIFIED_62,
    NAL_UNIT_UNSPECIFIED_63,
    NAL_UNIT_INVALID,
};

enum SliceType { B_SLICE = 0, P_SLICE = 1, I_SLICE = 2, NUMBER_OF_SLICE_TYPES = 3 };

namespace Profile {
enum Name { NONE = 0, MAIN = 1, MAIN10 = 2, MAINSTILLPICTURE = 3, MAINREXT = 4, HIGHTHROUGHPUTREXT = 5 };
}

namespace Level {
enum Tier { MAIN = 0, HIGH = 1, NUMBER_OF_TIERS = 2 };

enum Name {
    // code = (level * 30)
    NONE = 0,
    LEVEL1 = 30,
    LEVEL2 = 60,
    LEVEL2_1 = 63,
    LEVEL3 = 90,
    LEVEL3_1 = 93,
    LEVEL4 = 120,
    LEVEL4_1 = 123,
    LEVEL5 = 150,
    LEVEL5_1 = 153,
    LEVEL5_2 = 156,
    LEVEL6 = 180,
    LEVEL6_1 = 183,
    LEVEL6_2 = 186,
    LEVEL8_5 = 255,
};
}  // namespace Level

enum {
    // 7.4.2.1: vps_video_parameter_set_id is u(4).
    HEVC_MAX_VPS_COUNT = 16,
    // 7.4.3.2.1: sps_seq_parameter_set_id is in [0, 15].
    HEVC_MAX_SPS_COUNT = 16,
    // 7.4.3.3.1: pps_pic_parameter_set_id is in [0, 63].
    HEVC_MAX_PPS_COUNT = 64,
};

struct PointAddr {
    int x;
    int y;
};

enum DeblockEdgeDir {
    EDGE_VER = 0,
    EDGE_HOR = 1,
    EDGE_DIR_COUNT,
};

enum SaoMode {
    SAO_MODE_OFF,
    SAO_MODE_EO_0,
    SAO_MODE_EO_90,
    SAO_MODE_EO_135,
    SAO_MODE_EO_45,
    SAO_MODE_BO,
    SAO_MODE_MERGE_LEFT,
    SAO_MODE_MERGE_ABOVE,
};

struct SaoParam {
    SaoMode mode = SAO_MODE_OFF;
    int32_t offset[4];
    uint32_t bandPosition;

    SaoParam &operator=(const SaoParam &saoParam) {
        mode = saoParam.mode;
        bandPosition = saoParam.bandPosition;
        for (int i = 0; i < 4; i++) offset[i] = saoParam.offset[i];
        return *this;
    }
};

enum PredMode {
    PRED_MODE_INTER = 0,
    PRED_MODE_INTRA = 1,
    PRED_MODE_SKIP,
    PRED_MODE_COUNT,
};

enum PartMode {
    PART_MODE_2Nx2N = 0,  ///< symmetric motion partition,  2Nx2N
    PART_MODE_2NxN = 1,   ///< symmetric motion partition,  2Nx N
    PART_MODE_Nx2N = 2,   ///< symmetric motion partition,   Nx2N
    PART_MODE_NxN = 3,    ///< symmetric motion partition,   Nx N
    PART_MODE_2NxnU = 4,  ///< asymmetric motion partition, 2Nx( N/2) + 2Nx(3N/2)
    PART_MODE_2NxnD = 5,  ///< asymmetric motion partition, 2Nx(3N/2) + 2Nx( N/2)
    PART_MODE_nLx2N = 6,  ///< asymmetric motion partition, ( N/2)x2N + (3N/2)x2N
    PART_MODE_nRx2N = 7,  ///< asymmetric motion partition, (3N/2)x2N + ( N/2)x2N
    PART_MODE_COUNT = 8,
};

enum IntraPredMode {
    INTRA_PRED_MODE_PLANAR = 0,
    INTRA_PRED_MODE_DC,
    INTRA_PRED_MODE_ANGULAR_2,
    INTRA_PRED_MODE_ANGULAR_3,
    INTRA_PRED_MODE_ANGULAR_4,
    INTRA_PRED_MODE_ANGULAR_5,
    INTRA_PRED_MODE_ANGULAR_6,
    INTRA_PRED_MODE_ANGULAR_7,
    INTRA_PRED_MODE_ANGULAR_8,
    INTRA_PRED_MODE_ANGULAR_9,
    INTRA_PRED_MODE_ANGULAR_10,
    INTRA_PRED_MODE_ANGULAR_11,
    INTRA_PRED_MODE_ANGULAR_12,
    INTRA_PRED_MODE_ANGULAR_13,
    INTRA_PRED_MODE_ANGULAR_14,
    INTRA_PRED_MODE_ANGULAR_15,
    INTRA_PRED_MODE_ANGULAR_16,
    INTRA_PRED_MODE_ANGULAR_17,
    INTRA_PRED_MODE_ANGULAR_18,
    INTRA_PRED_MODE_ANGULAR_19,
    INTRA_PRED_MODE_ANGULAR_20,
    INTRA_PRED_MODE_ANGULAR_21,
    INTRA_PRED_MODE_ANGULAR_22,
    INTRA_PRED_MODE_ANGULAR_23,
    INTRA_PRED_MODE_ANGULAR_24,
    INTRA_PRED_MODE_ANGULAR_25,
    INTRA_PRED_MODE_ANGULAR_26,
    INTRA_PRED_MODE_ANGULAR_27,
    INTRA_PRED_MODE_ANGULAR_28,
    INTRA_PRED_MODE_ANGULAR_29,
    INTRA_PRED_MODE_ANGULAR_30,
    INTRA_PRED_MODE_ANGULAR_31,
    INTRA_PRED_MODE_ANGULAR_32,
    INTRA_PRED_MODE_ANGULAR_33,
    INTRA_PRED_MODE_ANGULAR_34,
};

enum InterPredIdc {
    INTER_PRED_L0 = 0,
    INTER_PRED_L1,
    INTER_PRED_BI,
};

struct Mv {
    int x = 0;
    int y = 0;

    bool operator==(const Mv &mv) const {
        if (x == mv.x && y == mv.y) return true;
        return false;
    }

    bool operator!=(const Mv &mv) const {
        if (x != mv.x || y != mv.y) return true;
        return false;
    }

    Mv &operator=(const Mv &mv) {
        x = mv.x;
        y = mv.y;
        return *this;
    }

    void reset() {
        x = 0;
        y = 0;
    }
};

struct MvField {
    int availableFlag = 0;
    Mv mvLX[2];
    Mv mvCLX[2];
    int refIdxLX[2] = {-1, -1};
    int predFlagLX[2] = {0, 0};

    bool operator==(const MvField &mvField) const {
        for (int i = 0; i < 2; i++) {
            if (mvLX[i] != mvField.mvLX[i]) return false;
            if (refIdxLX[i] != mvField.refIdxLX[i]) return false;
        }
        return true;
    }

    MvField &operator=(const MvField &mvField) {
        for (int i = 0; i < 2; i++) {
            mvLX[i] = mvField.mvLX[i];
            mvCLX[i] = mvField.mvCLX[i];
            refIdxLX[i] = mvField.refIdxLX[i];
            predFlagLX[i] = mvField.predFlagLX[i];
        }
        return *this;
    }

    void reset() {
        availableFlag = 0;
        for (int i = 0; i < 2; i++) {
            mvLX[i].reset();
            mvCLX[i].reset();
            refIdxLX[i] = -1;
            predFlagLX[i] = 0;
        }
    }
};

struct HevcCommonProfileTierLevel {
    uint8_t profile_space = 0;
    uint8_t tier_flag = 0;
    uint8_t profile_idc = 0;
    uint8_t profile_compatibility_flag[32] = {0};
    uint8_t progressive_source_flag = 0;
    uint8_t interlaced_source_flag = 0;
    uint8_t non_packed_constraint_flag = 0;
    uint8_t frame_only_constraint_flag = 0;
    uint8_t max_12bit_constraint_flag = 0;
    uint8_t max_10bit_constraint_flag = 0;
    uint8_t max_8bit_constraint_flag = 0;
    uint8_t max_422chroma_constraint_flag = 0;
    uint8_t max_420chroma_constraint_flag = 0;
    uint8_t max_monochrome_constraint_flag = 0;
    uint8_t intra_constraint_flag = 0;
    uint8_t one_picture_only_constraint_flag = 0;
    uint8_t lower_bit_rate_constraint_flag = 0;
    uint8_t inbld_flag = 0;
};

struct HevcProfileTierLevel {
    HevcCommonProfileTierLevel general_profile_tier_level;
    uint8_t general_level_idc = 0;
    std::vector<uint8_t> sub_layer_profile_present_flag;
    std::vector<uint8_t> sub_layer_level_present_flag;
    std::vector<HevcCommonProfileTierLevel> sub_layer_profile_tier_level;
    std::vector<uint8_t> sub_layer_level_idc;
};

struct HevcScalingList {};

struct HevcStRefPicSet {
    uint8_t inter_ref_pic_set_prediction_flag = 0;
    uint32_t delta_idx_minus1 = 0;
    uint8_t delta_rps_sign = 0;
    uint32_t abs_delta_rps_minus1 = 0;
    std::vector<uint8_t> used_by_curr_pic_flag;
    std::vector<uint8_t> use_delta_flag;

    uint32_t num_negative_pics = 0;
    uint32_t num_positive_pics = 0;
    std::vector<uint32_t> delta_poc_s0_minus1;
    std::vector<uint8_t> used_by_curr_pic_s0_flag;
    std::vector<uint32_t> delta_poc_s1_minus1;
    std::vector<uint8_t> used_by_curr_pic_s1_flag;

    // calculated
    uint32_t numNegativePics;
    uint32_t numPositivePics;
    uint32_t numDeltaPocs;
    std::vector<uint8_t> usedByCurrPicS0;
    std::vector<uint8_t> usedByCurrPicS1;
    std::vector<int32_t> deltaPocS0;
    std::vector<int32_t> deltaPocS1;
};

struct HevcVuiParameters {};

struct HevcVps {
    uint8_t vps_video_parameter_set_id = 0;
    uint8_t vps_base_layer_internal_flag = 0;
    uint8_t vps_base_layer_available_flag = 0;
    uint8_t vps_max_layers_minus1 = 0;
    uint8_t vps_max_sub_layers_minus1 = 0;
    uint8_t vps_temporal_id_nesting_flag = 0;
    uint16_t vps_reserved_0xffff_16bits = 0xffff;
    HevcProfileTierLevel profile_tier_level;
    uint8_t vps_sub_layer_ordering_info_present_flag = 0;
    std::vector<uint32_t> vps_max_dec_pic_buffering_minus1;
    std::vector<uint32_t> vps_max_num_reorder_pics;
    std::vector<uint32_t> vps_max_latency_increase_plus1;
    uint8_t vps_max_layer_id = 0;
    uint32_t vps_num_layer_sets_minus1 = 0;
    std::vector<std::vector<uint8_t>> layer_id_included_flag;
    uint8_t vps_timing_info_present_flag = 0;
    uint32_t vps_num_units_in_tick = 0;
    uint32_t vps_time_scale = 0;
    uint8_t vps_poc_proportional_to_timing_flag = 0;
    uint32_t vps_num_ticks_poc_diff_one_minus1 = 0;
    uint32_t vps_num_hrd_parameters = 0;
    uint8_t vps_extension_flag = 0;
};

struct HevcSps {
    uint8_t sps_video_parameter_set_id = 0;
    uint8_t sps_max_sub_layers_minus1 = 0;
    uint8_t sps_temporal_id_nesting_flag = 0;
    HevcProfileTierLevel profile_tier_level;
    uint8_t sps_seq_parameter_set_id = 0;
    uint8_t chroma_format_idc = 0;
    uint8_t separate_colour_plane_flag = 0;
    uint32_t pic_width_in_luma_samples = 0;
    uint32_t pic_height_in_luma_samples = 0;
    uint8_t conformance_window_flag = 0;
    uint32_t conf_win_left_offset = 0;
    uint32_t conf_win_right_offset = 0;
    uint32_t conf_win_top_offset = 0;
    uint32_t conf_win_bottom_offset = 0;
    uint8_t bit_depth_luma_minus8 = 0;
    uint8_t bit_depth_chroma_minus8 = 0;
    uint32_t log2_max_pic_order_cnt_lsb_minus4 = 0;
    uint8_t sps_sub_layer_ordering_info_present_flag = 0;
    std::vector<uint32_t> sps_max_dec_pic_buffering_minus1;
    std::vector<uint32_t> sps_max_num_reorder_pics;
    std::vector<uint32_t> sps_max_latency_increase_plus1;
    uint8_t log2_min_luma_coding_block_size_minus3 = 0;
    uint8_t log2_diff_max_min_luma_coding_block_size = 0;
    uint8_t log2_min_luma_transform_block_size_minus2 = 0;
    uint8_t log2_diff_max_min_luma_transform_block_size = 0;
    uint8_t max_transform_hierarchy_depth_inter = 0;
    uint8_t max_transform_hierarchy_depth_intra = 0;
    uint8_t scaling_list_enabled_flag = 0;
    uint8_t sps_scaling_list_data_present_flag = 0;
    HevcScalingList scaling_list_data;
    uint8_t amp_enabled_flag = 0;
    uint8_t sample_adaptive_offset_enabled_flag = 0;
    uint8_t pcm_enabled_flag = 0;
    uint8_t pcm_sample_bit_depth_luma_minus1 = 0;
    uint8_t pcm_sample_bit_depth_chroma_minus1 = 0;
    uint8_t log2_min_pcm_luma_coding_block_size_minus3 = 0;
    uint8_t log2_diff_max_min_pcm_luma_coding_block_size = 0;
    uint8_t pcm_loop_filter_disabled_flag = 0;
    uint8_t num_short_term_ref_pic_sets = 0;
    std::vector<HevcStRefPicSet> st_ref_pic_set;
    uint8_t long_term_ref_pics_present_flag = 0;
    uint8_t num_long_term_ref_pics_sps = 0;
    std::vector<uint32_t> lt_ref_pic_poc_lsb_sps;
    std::vector<uint8_t> used_by_curr_pic_lt_sps_flag;
    uint8_t sps_temporal_mvp_enabled_flag = 0;
    uint8_t strong_intra_smoothing_enabled_flag = 0;
    uint8_t vui_parameters_present_flag = 0;
    HevcVuiParameters vui_parameters;
    uint8_t sps_extension_present_flag = 0;
    uint8_t sps_range_extension_flag = 0;
    uint8_t sps_multilayer_extension_flag = 0;
    uint8_t sps_3d_extension_flag = 0;
    uint8_t sps_scc_extension_flag = 0;
    uint8_t sps_extension_4bits = 0;

    // calculated value
    uint8_t bitDepthY;
    uint8_t bitDepthC;
    uint8_t qpBdOffsetY;
    uint8_t qpBdOffsetC;
    uint8_t chromaArrayType;
    uint8_t subWidthC;
    uint8_t subHeightC;
    uint32_t width;
    uint32_t height;
    uint32_t log2MinCtbSize;
    uint32_t log2CtbSize;
    uint32_t ctbWidth;
    uint32_t ctbHeight;
    uint32_t ctbCount;
    uint32_t log2MinTbSize;
    uint32_t log2MaxTbSize;

    std::vector<PointAddr> scanOrderDiagonal1x1;
    std::vector<PointAddr> scanOrderDiagonal2x2;
    std::vector<PointAddr> scanOrderDiagonal4x4;
    std::vector<PointAddr> scanOrderDiagonal8x8;
    std::vector<PointAddr> scanOrderHorizontal1x1;
    std::vector<PointAddr> scanOrderHorizontal2x2;
    std::vector<PointAddr> scanOrderHorizontal4x4;
    std::vector<PointAddr> scanOrderHorizontal8x8;
    std::vector<PointAddr> scanOrderVertical1x1;
    std::vector<PointAddr> scanOrderVertical2x2;
    std::vector<PointAddr> scanOrderVertical4x4;
    std::vector<PointAddr> scanOrderVertical8x8;
    std::vector<PointAddr> scanOrderTraverse4x4;
    std::vector<PointAddr> scanOrderTraverse8x8;
    std::vector<PointAddr> scanOrderTraverse16x16;
    std::vector<PointAddr> scanOrderTraverse32x32;
};

struct HevcPps {
    uint8_t pps_pic_parameter_set_id = 0;
    uint8_t pps_seq_parameter_set_id = 0;
    uint8_t dependent_slice_segments_enabled_flag = 0;
    uint8_t output_flag_present_flag = 0;
    uint8_t num_extra_slice_header_bits = 0;
    uint8_t sign_data_hiding_enabled_flag = 0;
    uint8_t cabac_init_present_flag = 0;
    uint32_t num_ref_idx_l0_default_active_minus1 = 0;
    uint32_t num_ref_idx_l1_default_active_minus1 = 0;
    int32_t init_qp_minus26 = 0;
    uint8_t constrained_intra_pred_flag = 0;
    uint8_t transform_skip_enabled_flag = 0;
    uint8_t cu_qp_delta_enabled_flag = 0;
    uint32_t diff_cu_qp_delta_depth = 0;
    int32_t pps_cb_qp_offset = 0;
    int32_t pps_cr_qp_offset = 0;
    uint8_t pps_slice_chroma_qp_offsets_present_flag = 0;
    uint8_t weighted_pred_flag = 0;
    uint8_t weighted_bipred_flag = 0;
    uint8_t transquant_bypass_enabled_flag = 0;
    uint8_t tiles_enabled_flag = 0;
    uint8_t entropy_coding_sync_enabled_flag = 0;
    uint32_t num_tile_columns_minus1 = 0;
    uint32_t num_tile_rows_minus1 = 0;
    uint8_t uniform_spacing_flag = 0;
    std::vector<uint32_t> column_width_minus1;
    std::vector<uint32_t> row_height_minus1;
    uint8_t loop_filter_across_tiles_enabled_flag = 0;
    uint8_t pps_loop_filter_across_slices_enabled_flag = 0;
    uint8_t deblocking_filter_control_present_flag = 0;
    uint8_t deblocking_filter_override_enabled_flag = 0;
    uint8_t pps_deblocking_filter_disabled_flag = 0;
    int32_t pps_beta_offset_div2 = 0;
    int32_t pps_tc_offset_div2 = 0;
    uint8_t pps_scaling_list_data_present_flag = 0;
    HevcScalingList scaling_list_data;
    uint8_t lists_modification_present_flag = 0;
    uint32_t log2_parallel_merge_level_minus2 = 0;
    uint8_t slice_segment_header_extension_present_flag = 0;
    uint8_t pps_extension_present_flag = 0;
    uint8_t pps_range_extension_flag = 0;
    uint8_t pps_multilayer_extension_flag = 0;
    uint8_t pps_3d_extension_flag = 0;
    uint8_t pps_scc_extension_flag = 0;
    uint8_t pps_extension_4bits = 0;

    // calculated value
    uint32_t log2MinCuQpDeltaSize = 0;
    uint32_t log2MaxTransformSkipSize = 2;
    uint32_t log2ParMrgLevel;
    std::vector<uint32_t> tileColumnWidth;
    std::vector<uint32_t> tileRowHeight;
    std::vector<uint32_t> tileColumnBoundary;
    std::vector<uint32_t> tileRowBoundary;
    std::vector<uint32_t> ctbAddrRsToTs;
    std::vector<uint32_t> ctbAddrTsToRs;
    std::vector<uint32_t> ctbTileId;
    uint32_t minTbWidth;
    uint32_t minTbHeight;
    std::vector<std::vector<uint32_t>> minTbAddrZs;
};

struct HevcSei {};

// Weighted prediction parameters
struct HevcPredWeightTable {
    uint32_t luma_log2_weight_denom = 0;
    int32_t delta_chroma_log2_weight_denom = 0;
    std::vector<uint8_t> luma_weight_l0_flag;
    std::vector<uint8_t> chroma_weight_l0_flag;
    std::vector<int32_t> delta_luma_weight_l0;
    std::vector<int32_t> luma_offset_l0;
    std::vector<int32_t> delta_chroma_weight_l0;
    std::vector<int32_t> delta_chroma_offset_l0;
    std::vector<uint8_t> luma_weight_l1_flag;
    std::vector<uint8_t> chroma_weight_l1_flag;
    std::vector<int32_t> delta_luma_weight_l1;
    std::vector<int32_t> luma_offset_l1;
    std::vector<int32_t> delta_chroma_weight_l1;
    std::vector<int32_t> delta_chroma_offset_l1;
};

class HevcFrame;

struct HevcSliceHeader {
    uint8_t first_slice_segment_in_pic_flag = 0;
    uint8_t no_output_of_prior_pics_flag = 0;
    uint8_t slice_pic_parameter_set_id = 0;
    uint8_t dependent_slice_segment_flag = 0;
    uint32_t slice_segment_address = 0;
    std::vector<uint8_t> slice_reserved_flag;
    uint8_t slice_type = 0;
    uint8_t pic_output_flag = 0;
    uint8_t colour_plane_id = 0;
    uint32_t slice_pic_order_cnt_lsb = 0;
    uint8_t short_term_ref_pic_set_sps_flag = 0;
    HevcStRefPicSet st_ref_pic_set;
    uint32_t short_term_ref_pic_set_idx = 0;
    uint32_t num_long_term_sps = 0;
    uint32_t num_long_term_pics = 0;
    //

    uint8_t slice_temporal_mvp_enabled_flag = 0;
    uint8_t slice_sao_luma_flag = 0;
    uint8_t slice_sao_chroma_flag = 0;
    uint8_t num_ref_idx_active_override_flag = 0;
    uint32_t num_ref_idx_l0_active_minus1 = 0;
    uint32_t num_ref_idx_l1_active_minus1 = 0;
    // ref_pic_lists_modification
    uint8_t mvd_l1_zero_flag = 0;
    uint8_t cabac_init_flag = 0;
    uint8_t collocated_from_l0_flag = 1;
    uint32_t collocated_ref_idx = 0;
    HevcPredWeightTable pred_weight_table;
    uint32_t five_minus_max_num_merge_cand = 0;
    uint8_t use_integer_mv_flag = 0;

    int32_t slice_qp_delta = 0;
    int32_t slice_cb_qp_offset = 0;
    int32_t slice_cr_qp_offset = 0;
    int32_t slice_act_y_qp_offset = 0;
    int32_t slice_act_cb_qp_offset = 0;
    int32_t slice_act_cr_qp_offset = 0;
    uint8_t cu_chroma_qp_offset_enabled_flag = 0;
    uint8_t deblocking_filter_override_flag = 0;
    uint8_t slice_deblocking_filter_disabled_flag = 0;
    int32_t slice_beta_offset_div2 = 0;
    int32_t slice_tc_offset_div2 = 0;
    uint8_t slice_loop_filter_across_slices_enabled_flag = 0;
    uint32_t num_entry_point_offsets = 0;
    uint32_t offset_len_minus1 = 0;
    std::vector<uint32_t> entry_point_offset_minus1;
    uint32_t slice_segment_header_extension_length = 0;
    std::vector<uint8_t> slice_segment_header_extension_data_byte;

    // calculated value
    int32_t sliceQp;
    uint32_t maxNumMergeCand;
    uint32_t poc;
    std::vector<std::shared_ptr<HevcFrame>> refPicList[2];
};

#endif
