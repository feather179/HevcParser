#include "CavlcReader.h"

#include "Typedef.h"
#include "HevcFrame.h"
#include "ParameterSet.h"

#include <cstring>
#include <algorithm>

using std::shared_ptr;
using std::vector;

#define READ_SCODE(length, name) getSInt(length)
#define READ_CODE(length, name) getUInt(length)
#define READ_UVLC(name) getUVlc()
#define READ_SVLC(name) getSVlc()
#define READ_FLAG(name) getFlag()

static uint32_t ceilLog2(uint32_t value) {
    uint32_t ret = 0;
    while (value > (uint32_t)(1 << ret)) ret++;
    return ret;
}

CavlcReader::CavlcReader(const vector<uint8_t> &rbsp) : BitReader(rbsp.data(), (uint32_t)rbsp.size()) {}

CavlcReader::~CavlcReader() {}

void CavlcReader::parse(shared_ptr<ParameterSet> ps, int prevTid0Poc, vector<shared_ptr<HevcFrame>> &dpb) {
    READ_FLAG("forbidden_zero_bit");
    mNaluType = (NalUnitType)READ_CODE(6, "nal_unit_type");
    READ_CODE(6, "nuh_layer_id");
    mTid = READ_CODE(3, "nuh_temporal_id_plus1") - 1;

    switch (mNaluType) {
        case NAL_UNIT_VPS:
            parseVps(ps);
            break;
        case NAL_UNIT_SPS:
            parseSps(ps);
            break;
        case NAL_UNIT_PPS:
            parsePps(ps);
            break;
        case NAL_UNIT_PREFIX_SEI:
        case NAL_UNIT_SUFFIX_SEI:
            parseSei(ps);
            break;
        case NAL_UNIT_CODED_SLICE_TRAIL_R:
        case NAL_UNIT_CODED_SLICE_TRAIL_N:
        case NAL_UNIT_CODED_SLICE_TSA_R:
        case NAL_UNIT_CODED_SLICE_TSA_N:
        case NAL_UNIT_CODED_SLICE_STSA_R:
        case NAL_UNIT_CODED_SLICE_STSA_N:
        case NAL_UNIT_CODED_SLICE_BLA_W_LP:
        case NAL_UNIT_CODED_SLICE_BLA_W_RADL:
        case NAL_UNIT_CODED_SLICE_BLA_N_LP:
        case NAL_UNIT_CODED_SLICE_IDR_W_RADL:
        case NAL_UNIT_CODED_SLICE_IDR_N_LP:
        case NAL_UNIT_CODED_SLICE_CRA:
        case NAL_UNIT_CODED_SLICE_RADL_N:
        case NAL_UNIT_CODED_SLICE_RADL_R:
        case NAL_UNIT_CODED_SLICE_RASL_N:
        case NAL_UNIT_CODED_SLICE_RASL_R:
            parseSliceHeader(ps, prevTid0Poc, dpb);
            break;
        default:
            break;
    }
}

void CavlcReader::parseVps(shared_ptr<ParameterSet> ps) {
    auto vps = std::make_shared<HevcVps>();

    vps->vps_video_parameter_set_id = READ_CODE(4, "vps_video_parameter_set_id");
    vps->vps_base_layer_internal_flag = READ_FLAG("vps_base_layer_internal_flag");
    vps->vps_base_layer_available_flag = READ_FLAG("vps_base_layer_available_flag");
    vps->vps_max_layers_minus1 = READ_CODE(6, "vps_max_layers_minus1");
    vps->vps_max_sub_layers_minus1 = READ_CODE(3, "vps_max_sub_layers_minus1");
    vps->vps_temporal_id_nesting_flag = READ_FLAG("vps_temporal_id_nesting_flag");
    vps->vps_reserved_0xffff_16bits = READ_CODE(16, "vps_reserved_0xffff_16bits");
    parsePtl(&vps->profile_tier_level, true, vps->vps_max_sub_layers_minus1);
    vps->vps_sub_layer_ordering_info_present_flag = READ_FLAG("vps_sub_layer_ordering_info_present_flag");

    vps->vps_max_dec_pic_buffering_minus1.resize(vps->vps_max_sub_layers_minus1 + 1);
    vps->vps_max_num_reorder_pics.resize(vps->vps_max_sub_layers_minus1 + 1);
    vps->vps_max_latency_increase_plus1.resize(vps->vps_max_sub_layers_minus1 + 1);
    int i = vps->vps_sub_layer_ordering_info_present_flag ? 0 : vps->vps_max_sub_layers_minus1;
    for (; i <= vps->vps_max_sub_layers_minus1; i++) {
        vps->vps_max_dec_pic_buffering_minus1[i] = READ_UVLC("vps_max_dec_pic_buffering_minus1[i]");
        vps->vps_max_num_reorder_pics[i] = READ_UVLC("vps_max_num_reorder_pics[i]");
        vps->vps_max_latency_increase_plus1[i] = READ_UVLC("vps_max_latency_increase_plus1[i]");
    }

    vps->vps_max_layer_id = READ_CODE(6, "vps_max_layer_id");
    vps->vps_num_layer_sets_minus1 = READ_UVLC("vps_num_layer_sets_minus1");

    vps->layer_id_included_flag.resize(vps->vps_num_layer_sets_minus1 + 1);
    for (int i = 1; i <= (int)vps->vps_num_layer_sets_minus1; i++) {
        vps->layer_id_included_flag[i].resize(vps->vps_max_layer_id + 1);
        for (int j = 0; j <= vps->vps_max_layer_id; j++) {
            vps->layer_id_included_flag[i][j] = READ_FLAG("layer_id_included_flag[i][j]");
        }
    }

    vps->vps_timing_info_present_flag = READ_FLAG("vps_timing_info_present_flag");
    if (vps->vps_timing_info_present_flag) {
        vps->vps_num_units_in_tick = READ_CODE(32, "vps_num_units_in_tick");
        vps->vps_time_scale = READ_CODE(32, "vps_time_scale");
        vps->vps_poc_proportional_to_timing_flag = READ_FLAG("vps_poc_proportional_to_timing_flag");
        if (vps->vps_poc_proportional_to_timing_flag)
            vps->vps_num_ticks_poc_diff_one_minus1 = READ_UVLC("vps_num_ticks_poc_diff_one_minus1");
        vps->vps_num_hrd_parameters = READ_UVLC("vps_num_hrd_parameters");

        // TODO: hrd data
    }

    vps->vps_extension_flag = READ_FLAG("vps_extension_flag");

    ps->setVps(vps);
}

void CavlcReader::parseSps(shared_ptr<ParameterSet> ps) {
    auto sps = std::make_shared<HevcSps>();

    sps->sps_video_parameter_set_id = READ_CODE(4, "sps_video_parameter_set_id");
    sps->sps_max_sub_layers_minus1 = READ_CODE(3, "sps_max_sub_layers_minus1");
    sps->sps_temporal_id_nesting_flag = READ_FLAG("sps_temporal_id_nesting_flag");
    parsePtl(&sps->profile_tier_level, true, sps->sps_max_sub_layers_minus1);
    sps->sps_seq_parameter_set_id = READ_UVLC("sps_seq_parameter_set_id");
    sps->chroma_format_idc = READ_UVLC("chroma_format_idc");
    if (sps->chroma_format_idc == 3) sps->separate_colour_plane_flag = READ_FLAG("separate_colour_plane_flag");

    sps->pic_width_in_luma_samples = READ_UVLC("pic_width_in_luma_samples");
    sps->pic_height_in_luma_samples = READ_UVLC("pic_height_in_luma_samples");
    sps->conformance_window_flag = READ_FLAG("conformance_window_flag");
    if (sps->conformance_window_flag) {
        sps->conf_win_left_offset = READ_UVLC("conf_win_left_offset");
        sps->conf_win_right_offset = READ_UVLC("conf_win_right_offset");
        sps->conf_win_top_offset = READ_UVLC("conf_win_top_offset");
        sps->conf_win_bottom_offset = READ_UVLC("conf_win_bottom_offset");
    }
    sps->bit_depth_luma_minus8 = READ_UVLC("bit_depth_luma_minus8");
    sps->bit_depth_chroma_minus8 = READ_UVLC("bit_depth_chroma_minus8");
    sps->log2_max_pic_order_cnt_lsb_minus4 = READ_UVLC("log2_max_pic_order_cnt_lsb_minus4");
    sps->sps_sub_layer_ordering_info_present_flag = READ_FLAG("sps_sub_layer_ordering_info_present_flag");

    sps->sps_max_dec_pic_buffering_minus1.resize(sps->sps_max_sub_layers_minus1 + 1);
    sps->sps_max_num_reorder_pics.resize(sps->sps_max_sub_layers_minus1 + 1);
    sps->sps_max_latency_increase_plus1.resize(sps->sps_max_sub_layers_minus1 + 1);
    int i = sps->sps_sub_layer_ordering_info_present_flag ? 0 : sps->sps_max_sub_layers_minus1;
    for (; i <= sps->sps_max_sub_layers_minus1; i++) {
        sps->sps_max_dec_pic_buffering_minus1[i] = READ_UVLC("sps_max_dec_pic_buffering_minus1[i]");
        sps->sps_max_num_reorder_pics[i] = READ_UVLC("sps_max_num_reorder_pics[i]");
        sps->sps_max_latency_increase_plus1[i] = READ_UVLC("sps_max_latency_increase_plus1[i]");
    }

    sps->log2_min_luma_coding_block_size_minus3 = READ_UVLC("log2_min_luma_coding_block_size_minus3");
    sps->log2_diff_max_min_luma_coding_block_size = READ_UVLC("log2_diff_max_min_luma_coding_block_size");
    sps->log2_min_luma_transform_block_size_minus2 = READ_UVLC("log2_min_luma_transform_block_size_minus2");
    sps->log2_diff_max_min_luma_transform_block_size = READ_UVLC("log2_diff_max_min_luma_transform_block_size");
    sps->max_transform_hierarchy_depth_inter = READ_UVLC("max_transform_hierarchy_depth_inter");
    sps->max_transform_hierarchy_depth_intra = READ_UVLC("max_transform_hierarchy_depth_intra");
    sps->scaling_list_enabled_flag = READ_FLAG("scaling_list_enabled_flag");
    if (sps->scaling_list_enabled_flag) {
        sps->sps_scaling_list_data_present_flag = READ_FLAG("sps_scaling_list_data_present_flag");
        if (sps->sps_scaling_list_data_present_flag) {
            // TODO: parse scaling_list_data
        }
    }

    sps->amp_enabled_flag = READ_FLAG("amp_enabled_flag");
    sps->sample_adaptive_offset_enabled_flag = READ_FLAG("sample_adaptive_offset_enabled_flag");
    sps->pcm_enabled_flag = READ_FLAG("pcm_enabled_flag");
    if (sps->pcm_enabled_flag) {
        sps->pcm_sample_bit_depth_luma_minus1 = READ_CODE(4, "pcm_sample_bit_depth_luma_minus1");
        sps->pcm_sample_bit_depth_chroma_minus1 = READ_CODE(4, "pcm_sample_bit_depth_chroma_minus1");
        sps->log2_min_pcm_luma_coding_block_size_minus3 = READ_UVLC("log2_min_pcm_luma_coding_block_size_minus3");
        sps->log2_diff_max_min_pcm_luma_coding_block_size = READ_UVLC("log2_diff_max_min_pcm_luma_coding_block_size");
        sps->pcm_loop_filter_disabled_flag = READ_FLAG("pcm_loop_filter_disabled_flag");
    }

    sps->num_short_term_ref_pic_sets = READ_UVLC("num_short_term_ref_pic_sets");
    sps->st_ref_pic_set.resize(sps->num_short_term_ref_pic_sets);
    for (int i = 0; i < sps->num_short_term_ref_pic_sets; i++)
        parseStRps(&sps->st_ref_pic_set[i], i, sps->num_short_term_ref_pic_sets);

    sps->long_term_ref_pics_present_flag = READ_FLAG("long_term_ref_pics_present_flag");
    if (sps->long_term_ref_pics_present_flag) {
        sps->num_long_term_ref_pics_sps = READ_UVLC("num_long_term_ref_pics_sps");
        sps->lt_ref_pic_poc_lsb_sps.resize(sps->num_long_term_ref_pics_sps);
        sps->used_by_curr_pic_lt_sps_flag.resize(sps->num_long_term_ref_pics_sps);
        for (int i = 0; i < sps->num_long_term_ref_pics_sps; i++) {
            uint32_t bits = sps->log2_max_pic_order_cnt_lsb_minus4 + 4;
            sps->lt_ref_pic_poc_lsb_sps[i] = READ_CODE(bits, "lt_ref_pic_poc_lsb_sps[i]");
            sps->used_by_curr_pic_lt_sps_flag[i] = READ_FLAG("used_by_curr_pic_lt_sps_flag[i]");
        }
    }

    sps->sps_temporal_mvp_enabled_flag = READ_FLAG("sps_temporal_mvp_enabled_flag");
    sps->strong_intra_smoothing_enabled_flag = READ_FLAG("strong_intra_smoothing_enabled_flag");
    sps->vui_parameters_present_flag = READ_FLAG("vui_parameters_present_flag");
    if (sps->vui_parameters_present_flag) {
        // TODO: parse vui parameters
    }

    sps->sps_extension_present_flag = READ_FLAG("sps_extension_present_flag");
    if (sps->sps_extension_present_flag) {
        sps->sps_range_extension_flag = READ_FLAG("sps_range_extension_flag");
        sps->sps_multilayer_extension_flag = READ_FLAG("sps_multilayer_extension_flag");
        sps->sps_3d_extension_flag = READ_FLAG("sps_3d_extension_flag");
        sps->sps_scc_extension_flag = READ_FLAG("sps_scc_extension_flag");
        sps->sps_extension_4bits = READ_FLAG("sps_extension_4bits");
    }

    if (sps->sps_range_extension_flag) {
        // TODO: parse sps_range_extension
    }
    if (sps->sps_multilayer_extension_flag) {
        // TODO: parse sps_multilayer_extension
    }
    if (sps->sps_3d_extension_flag) {
        // TODO: parse sps_3d_extension
    }
    if (sps->sps_scc_extension_flag) {
        // TODO: parse sps_scc_extension
    }

    // calculated value
    sps->bitDepthY = sps->bit_depth_luma_minus8 + 8;
    sps->bitDepthC = sps->bit_depth_chroma_minus8 + 8;
    sps->qpBdOffsetY = sps->bit_depth_luma_minus8 * 6;
    sps->qpBdOffsetC = sps->bit_depth_chroma_minus8 * 6;
    sps->chromaArrayType = sps->chroma_format_idc;
    if (sps->separate_colour_plane_flag) sps->chromaArrayType = 0;
    sps->subWidthC = sps->subHeightC = 1;
    if (sps->chroma_format_idc == 1)
        sps->subWidthC = sps->subHeightC = 2;
    else if (sps->chroma_format_idc == 2)
        sps->subWidthC = 2;
    sps->width = sps->pic_width_in_luma_samples;
    sps->height = sps->pic_height_in_luma_samples;
    sps->log2MinCtbSize = sps->log2_min_luma_coding_block_size_minus3 + 3;
    sps->log2CtbSize = sps->log2MinCtbSize + sps->log2_diff_max_min_luma_coding_block_size;
    sps->ctbWidth = (sps->width + (1 << sps->log2CtbSize) - 1) >> sps->log2CtbSize;
    sps->ctbHeight = (sps->height + (1 << sps->log2CtbSize) - 1) >> sps->log2CtbSize;
    sps->ctbCount = sps->ctbWidth * sps->ctbHeight;
    sps->log2MinTbSize = sps->log2_min_luma_transform_block_size_minus2 + 2;
    sps->log2MaxTbSize = sps->log2MinTbSize + sps->log2_diff_max_min_luma_transform_block_size;

    auto calculateDiagonalScanOrder = [](int blkSize, vector<PointAddr> &scanOrder) {
        int i = 0, x = 0, y = 0;
        while (i < blkSize * blkSize) {
            while (y >= 0) {
                if (x < blkSize && y < blkSize) {
                    scanOrder.push_back({x, y});
                    i++;
                }
                y--;
                x++;
            }
            y = x;
            x = 0;
        }
    };

    auto calculateHorizontalScanOrder = [](int blkSize, vector<PointAddr> &scanOrder) {
        int i = 0;
        for (int y = 0; y < blkSize; y++) {
            for (int x = 0; x < blkSize; x++) {
                scanOrder.push_back({x, y});
                i++;
            }
        }
    };

    auto calculateVerticalScanOrder = [](int blkSize, vector<PointAddr> &scanOrder) {
        int i = 0;
        for (int x = 0; x < blkSize; x++) {
            for (int y = 0; y < blkSize; y++) {
                scanOrder.push_back({x, y});
                i++;
            }
        }
    };

    auto calculateTraverseScanOrder = [](int blkSize, vector<PointAddr> &scanOrder) {
        int i = 0;
        for (int y = 0; y < blkSize; y++) {
            if (y & 1) {
                for (int x = (blkSize - 1); x >= 0; x--) {
                    scanOrder.push_back({x, y});
                    i++;
                }
            } else {
                for (int x = 0; x < blkSize; x++) {
                    scanOrder.push_back({x, y});
                    i++;
                }
            }
        }
    };

    calculateDiagonalScanOrder(1, sps->scanOrderDiagonal1x1);
    calculateDiagonalScanOrder(2, sps->scanOrderDiagonal2x2);
    calculateDiagonalScanOrder(4, sps->scanOrderDiagonal4x4);
    calculateDiagonalScanOrder(8, sps->scanOrderDiagonal8x8);
    calculateHorizontalScanOrder(1, sps->scanOrderHorizontal1x1);
    calculateHorizontalScanOrder(2, sps->scanOrderHorizontal2x2);
    calculateHorizontalScanOrder(4, sps->scanOrderHorizontal4x4);
    calculateHorizontalScanOrder(8, sps->scanOrderHorizontal8x8);
    calculateVerticalScanOrder(1, sps->scanOrderVertical1x1);
    calculateVerticalScanOrder(2, sps->scanOrderVertical2x2);
    calculateVerticalScanOrder(4, sps->scanOrderVertical4x4);
    calculateVerticalScanOrder(8, sps->scanOrderVertical8x8);
    calculateTraverseScanOrder(4, sps->scanOrderTraverse4x4);
    calculateTraverseScanOrder(8, sps->scanOrderTraverse8x8);
    calculateTraverseScanOrder(16, sps->scanOrderTraverse16x16);
    calculateTraverseScanOrder(32, sps->scanOrderTraverse32x32);

    ps->setSps(sps);
}

void CavlcReader::parsePps(shared_ptr<ParameterSet> ps) {
    auto pps = std::make_shared<HevcPps>();

    pps->pps_pic_parameter_set_id = READ_UVLC("pps_pic_parameter_set_id");
    pps->pps_seq_parameter_set_id = READ_UVLC("pps_seq_parameter_set_id");
    pps->dependent_slice_segments_enabled_flag = READ_FLAG("dependent_slice_segments_enabled_flag");
    pps->output_flag_present_flag = READ_FLAG("output_flag_present_flag");
    pps->num_extra_slice_header_bits = READ_CODE(3, "num_extra_slice_header_bits");
    pps->sign_data_hiding_enabled_flag = READ_FLAG("sign_data_hiding_enabled_flag");
    pps->cabac_init_present_flag = READ_FLAG("cabac_init_present_flag");
    pps->num_ref_idx_l0_default_active_minus1 = READ_UVLC("num_ref_idx_l0_default_active_minus1");
    pps->num_ref_idx_l1_default_active_minus1 = READ_UVLC("num_ref_idx_l1_default_active_minus1");
    pps->init_qp_minus26 = READ_SVLC("init_qp_minus26");
    pps->constrained_intra_pred_flag = READ_FLAG("constrained_intra_pred_flag");
    pps->transform_skip_enabled_flag = READ_FLAG("transform_skip_enabled_flag");
    pps->cu_qp_delta_enabled_flag = READ_FLAG("cu_qp_delta_enabled_flag");
    if (pps->cu_qp_delta_enabled_flag) pps->diff_cu_qp_delta_depth = READ_UVLC("diff_cu_qp_delta_depth");
    pps->pps_cb_qp_offset = READ_SVLC("pps_cb_qp_offset");
    pps->pps_cr_qp_offset = READ_SVLC("pps_cr_qp_offset");
    pps->pps_slice_chroma_qp_offsets_present_flag = READ_FLAG("pps_slice_chroma_qp_offsets_present_flag");
    pps->weighted_pred_flag = READ_FLAG("weighted_pred_flag");
    pps->weighted_bipred_flag = READ_FLAG("weighted_bipred_flag");
    pps->transquant_bypass_enabled_flag = READ_FLAG("transquant_bypass_enabled_flag");
    pps->tiles_enabled_flag = READ_FLAG("tiles_enabled_flag");
    pps->entropy_coding_sync_enabled_flag = READ_FLAG("entropy_coding_sync_enabled_flag");
    if (pps->tiles_enabled_flag) {
        pps->num_tile_columns_minus1 = READ_UVLC("num_tile_columns_minus1");
        pps->num_tile_rows_minus1 = READ_UVLC("num_tile_rows_minus1");
        pps->uniform_spacing_flag = READ_FLAG("uniform_spacing_flag");
        if (!pps->uniform_spacing_flag) {
            pps->column_width_minus1.resize(pps->num_tile_columns_minus1);
            for (int i = 0; i < (int)pps->num_tile_columns_minus1; i++)
                pps->column_width_minus1[i] = READ_UVLC("column_width_minus1");
            pps->row_height_minus1.resize(pps->num_tile_rows_minus1);
            for (int i = 0; i < (int)pps->num_tile_rows_minus1; i++)
                pps->row_height_minus1[i] = READ_UVLC("row_height_minus1");
        }
        pps->loop_filter_across_tiles_enabled_flag = READ_FLAG("loop_filter_across_tiles_enabled_flag");
    }

    pps->pps_loop_filter_across_slices_enabled_flag = READ_FLAG("pps_loop_filter_across_slices_enabled_flag");
    pps->deblocking_filter_control_present_flag = READ_FLAG("deblocking_filter_control_present_flag");
    if (pps->deblocking_filter_control_present_flag) {
        pps->deblocking_filter_override_enabled_flag = READ_FLAG("deblocking_filter_override_enabled_flag");
        pps->pps_deblocking_filter_disabled_flag = READ_FLAG("pps_deblocking_filter_disabled_flag");
        if (!pps->pps_deblocking_filter_disabled_flag) {
            pps->pps_beta_offset_div2 = READ_SVLC("pps_beta_offset_div2");
            pps->pps_tc_offset_div2 = READ_SVLC("pps_tc_offset_div2");
        }
    }

    pps->pps_scaling_list_data_present_flag = READ_FLAG("pps_scaling_list_data_present_flag");
    if (pps->pps_scaling_list_data_present_flag) {
        // TODO: parse scaling_list_data
    }

    pps->lists_modification_present_flag = READ_FLAG("lists_modification_present_flag");
    pps->log2_parallel_merge_level_minus2 = READ_UVLC("log2_parallel_merge_level_minus2");
    pps->slice_segment_header_extension_present_flag = READ_FLAG("slice_segment_header_extension_present_flag");
    pps->pps_extension_present_flag = READ_FLAG("pps_extension_present_flag");
    if (pps->pps_extension_present_flag) {
        pps->pps_range_extension_flag = READ_FLAG("pps_range_extension_flag");
        pps->pps_multilayer_extension_flag = READ_FLAG("pps_multilayer_extension_flag");
        pps->pps_3d_extension_flag = READ_FLAG("pps_3d_extension_flag");
        pps->pps_scc_extension_flag = READ_FLAG("pps_scc_extension_flag");
        pps->pps_extension_4bits = READ_CODE(4, "pps_extension_4bits");
    }

    if (pps->pps_range_extension_flag) {
        // TODO: parse pps_range_extension_flag
    }
    if (pps->pps_multilayer_extension_flag) {
        // TODO: parse pps_multilayer_extension_flag
    }
    if (pps->pps_3d_extension_flag) {
        // TODO: parse pps_3d_extension_flag
    }
    if (pps->pps_scc_extension_flag) {
        // TODO: parse pps_scc_extension_flag
    }

    // calculated value
    auto sps = ps->getSps(pps->pps_seq_parameter_set_id);

    pps->log2MinCuQpDeltaSize = sps->log2CtbSize - pps->diff_cu_qp_delta_depth;
    pps->log2ParMrgLevel = pps->log2_parallel_merge_level_minus2 + 2;

    int tileColumnCount = (int)pps->num_tile_columns_minus1 + 1;
    int tileRowCount = (int)pps->num_tile_rows_minus1 + 1;
    if (pps->uniform_spacing_flag) {
        for (int i = 0; i < tileColumnCount; i++) {
            uint32_t width = (i + 1) * sps->ctbWidth / tileColumnCount - i * sps->ctbWidth / tileColumnCount;
            pps->tileColumnWidth.push_back(width);
        }

        for (int i = 0; i < tileRowCount; i++) {
            uint32_t height = (i + 1) * sps->ctbHeight / tileRowCount - i * sps->ctbHeight / tileRowCount;
            pps->tileRowHeight.push_back(height);
        }
    } else {
        uint32_t width = sps->ctbWidth;
        for (int i = 0; i < (tileColumnCount - 1); i++) {
            pps->tileColumnWidth.push_back(pps->column_width_minus1[i] + 1);
            width -= pps->tileColumnWidth[i];
        }
        pps->tileColumnWidth.push_back(width);

        uint32_t height = sps->ctbHeight;
        for (int i = 0; i < (tileRowCount - 1); i++) {
            pps->tileRowHeight.push_back(pps->row_height_minus1[i] + 1);
            height -= pps->tileRowHeight[i];
        }
        pps->tileRowHeight.push_back(height);
    }

    pps->tileColumnBoundary.resize(tileColumnCount + 1);
    pps->tileColumnBoundary[0] = 0;
    for (int i = 0; i < tileColumnCount; i++)
        pps->tileColumnBoundary[i + 1] = pps->tileColumnBoundary[i] + pps->tileColumnWidth[i];

    pps->tileRowBoundary.resize(tileRowCount + 1);
    pps->tileRowBoundary[0] = 0;
    for (int i = 0; i < tileRowCount; i++)
        pps->tileRowBoundary[i + 1] = pps->tileRowBoundary[i] + pps->tileRowHeight[i];

    pps->ctbAddrRsToTs.resize(sps->ctbCount);
    pps->ctbAddrTsToRs.resize(sps->ctbCount);
    for (int rsAddr = 0; rsAddr < (int)sps->ctbCount; rsAddr++) {
        uint32_t tbX = rsAddr % sps->ctbWidth;
        uint32_t tbY = rsAddr / sps->ctbWidth;
        int tileX = tileColumnCount - 1;
        int tileY = tileRowCount - 1;
        while (tileX >= 0 && tbX < pps->tileColumnBoundary[tileX]) tileX--;
        while (tileY >= 0 && tbY < pps->tileRowBoundary[tileY]) tileY--;
        pps->ctbAddrRsToTs[rsAddr] = 0;
        for (int i = 0; i < tileX; i++)
            pps->ctbAddrRsToTs[rsAddr] += pps->tileRowHeight[tileY] * pps->tileColumnWidth[i];
        for (int i = 0; i < tileY; i++) pps->ctbAddrRsToTs[rsAddr] += sps->ctbWidth * pps->tileRowHeight[i];
        pps->ctbAddrRsToTs[rsAddr] +=
            (tbY - pps->tileRowBoundary[tileY]) * pps->tileColumnWidth[tileX] + tbX - pps->tileColumnBoundary[tileX];
        pps->ctbAddrTsToRs[pps->ctbAddrRsToTs[rsAddr]] = rsAddr;
    }

    pps->ctbTileId.resize(sps->ctbCount);
    uint32_t tileIdx = 0;
    for (int j = 0; j < tileRowCount; j++) {
        for (int i = 0; i < tileColumnCount; i++, tileIdx++) {
            for (int y = pps->tileRowBoundary[j]; y < (int)pps->tileRowBoundary[j + 1]; y++) {
                for (int x = pps->tileColumnBoundary[i]; x < (int)pps->tileColumnBoundary[i + 1]; x++) {
                    pps->ctbTileId[pps->ctbAddrRsToTs[y * sps->ctbWidth + x]] = tileIdx;
                }
            }
        }
    }

    pps->minTbWidth = sps->ctbWidth << (sps->log2CtbSize - sps->log2MinTbSize);
    pps->minTbHeight = sps->ctbHeight << (sps->log2CtbSize - sps->log2MinTbSize);
    pps->minTbAddrZs.resize(pps->minTbHeight, vector<uint32_t>(pps->minTbWidth, 0));
    for (int y = 0; y < (int)pps->minTbHeight; y++) {
        for (int x = 0; x < (int)pps->minTbWidth; x++) {
            uint32_t tbX = (x << sps->log2MinTbSize) >> sps->log2CtbSize;
            uint32_t tbY = (y << sps->log2MinTbSize) >> sps->log2CtbSize;
            pps->minTbAddrZs[y][x] = pps->ctbAddrRsToTs[sps->ctbWidth * tbY + tbX]
                                     << (2 * (sps->log2CtbSize - sps->log2MinTbSize));
            for (int i = 0; i < (int)(sps->log2CtbSize - sps->log2MinTbSize); i++) {
                uint32_t m = 1 << i;
                pps->minTbAddrZs[y][x] += m & x ? m * m : 0;
                pps->minTbAddrZs[y][x] += m & y ? 2 * m * m : 0;
            }
        }
    }

    ps->setPps(pps);
}

void CavlcReader::parseSei(shared_ptr<ParameterSet> ps) {}

void CavlcReader::parseSliceHeader(shared_ptr<ParameterSet> ps, int prevTid0Poc, vector<shared_ptr<HevcFrame>> &dpb) {
    auto sliceHeader = std::make_shared<HevcSliceHeader>();
    shared_ptr<HevcSps> sps;
    shared_ptr<HevcPps> pps;

    auto parsePoc = [this, sliceHeader](int prevTid0Poc, int maxPocLsb) {
        if (mNaluType == NAL_UNIT_CODED_SLICE_IDR_W_RADL || mNaluType == NAL_UNIT_CODED_SLICE_IDR_N_LP)
            sliceHeader->poc = 0;
        else {
            int prevPocLsb = prevTid0Poc & (maxPocLsb - 1);
            int prevPocMsb = prevTid0Poc - prevPocLsb;
            int pocLsb = sliceHeader->slice_pic_order_cnt_lsb;
            int pocMsb = 0;

            if ((pocLsb < prevPocLsb) && ((prevPocLsb - pocLsb) >= (maxPocLsb / 2)))
                pocMsb = prevPocMsb + maxPocLsb;
            else if ((pocLsb > prevPocLsb) && ((pocLsb - prevPocLsb) > (maxPocLsb / 2)))
                pocMsb = prevPocMsb - maxPocLsb;
            else
                pocMsb = prevPocMsb;

            if (mNaluType == NAL_UNIT_CODED_SLICE_BLA_W_LP || mNaluType == NAL_UNIT_CODED_SLICE_BLA_W_RADL ||
                mNaluType == NAL_UNIT_CODED_SLICE_BLA_N_LP)
                pocMsb = 0;

            sliceHeader->poc = pocMsb + pocLsb;
        }
    };

    auto parseRpsList = [this, sliceHeader](HevcStRefPicSet *stRps, vector<shared_ptr<HevcFrame>> &dpb) {
        vector<int> pocStCurrBefore, pocStCurrAfter, pocStFoll;
        int numPocStCurrBefore, numPocStCurrAfter, numPocStFoll;
        int numPicTotalCurr = 0;
        int poc = sliceHeader->poc, tid = mTid;

        for (int i = 0; i < (int)stRps->numNegativePics; i++) {
            if (stRps->usedByCurrPicS0[i]) {
                pocStCurrBefore.push_back(poc + stRps->deltaPocS0[i]);
                numPicTotalCurr++;
            } else
                pocStFoll.push_back(poc + stRps->deltaPocS0[i]);
        }
        for (int i = 0; i < (int)stRps->numPositivePics; i++) {
            if (stRps->usedByCurrPicS1[i]) {
                pocStCurrAfter.push_back(poc + stRps->deltaPocS1[i]);
                numPicTotalCurr++;
            } else
                pocStFoll.push_back(poc + stRps->deltaPocS1[i]);
        }

        numPocStCurrBefore = (int)pocStCurrBefore.size();
        numPocStCurrAfter = (int)pocStCurrAfter.size();
        numPocStFoll = (int)pocStFoll.size();

        vector<shared_ptr<HevcFrame>> refPicSetStCurrBefore(numPocStCurrBefore);
        vector<shared_ptr<HevcFrame>> refPicSetStCurrAfter(numPocStCurrAfter);
        vector<shared_ptr<HevcFrame>> refPicSetStFoll(numPocStFoll);

        for (int i = 0; i < dpb.size(); i++) dpb[i]->setUsedForReference(false);

        for (int i = 0; i < numPocStCurrBefore; i++) {
            int refPoc = pocStCurrBefore[i];
            auto pic = std::find_if(dpb.begin(), dpb.end(), [refPoc, tid](shared_ptr<HevcFrame> frame) {
                return (frame->getPoc() == refPoc) && (frame->getTid() == tid);
            });
            if (pic != dpb.end()) {
                refPicSetStCurrBefore[i] = *pic;
                refPicSetStCurrBefore[i]->setUsedForReference(true);
            } else
                refPicSetStCurrBefore[i] = nullptr;
        }
        for (int i = 0; i < numPocStCurrAfter; i++) {
            int refPoc = pocStCurrAfter[i];
            auto pic = std::find_if(dpb.begin(), dpb.end(), [refPoc, tid](shared_ptr<HevcFrame> frame) {
                return (frame->getPoc() == refPoc) && (frame->getTid() == tid);
            });
            if (pic != dpb.end()) {
                refPicSetStCurrAfter[i] = *pic;
                refPicSetStCurrAfter[i]->setUsedForReference(true);
            } else
                refPicSetStCurrAfter[i] = nullptr;
        }
        for (int i = 0; i < numPocStFoll; i++) {
            int refPoc = pocStFoll[i];
            auto pic = std::find_if(dpb.begin(), dpb.end(), [refPoc, tid](shared_ptr<HevcFrame> frame) {
                return (frame->getPoc() == refPoc) && (frame->getTid() == tid);
            });
            if (pic != dpb.end()) {
                refPicSetStFoll[i] = *pic;
                refPicSetStFoll[i]->setUsedForReference(true);
            } else
                refPicSetStFoll[i] = nullptr;
        }

        {
            int numRpsCurrTempList0 = std::max<int>(sliceHeader->num_ref_idx_l0_active_minus1 + 1, numPicTotalCurr);
            int rIdx = 0;
            vector<shared_ptr<HevcFrame>> refPicListTemp0;

            while (rIdx < numRpsCurrTempList0) {
                for (int i = 0; i < numPocStCurrBefore && rIdx < numRpsCurrTempList0; i++) {
                    refPicListTemp0.push_back(refPicSetStCurrBefore[i]);
                    rIdx++;
                }
                for (int i = 0; i < numPocStCurrAfter && rIdx < numRpsCurrTempList0; i++) {
                    refPicListTemp0.push_back(refPicSetStCurrAfter[i]);
                    rIdx++;
                }
            }

            for (int i = 0; i <= (int)sliceHeader->num_ref_idx_l0_active_minus1; i++)
                sliceHeader->refPicList[0].push_back(refPicListTemp0[i]);
        }

        if (sliceHeader->slice_type == B_SLICE) {
            int numRpsCurrTempList1 = std::max<int>(sliceHeader->num_ref_idx_l1_active_minus1 + 1, numPicTotalCurr);
            int rIdx = 0;
            vector<shared_ptr<HevcFrame>> refPicListTemp1;

            while (rIdx < numRpsCurrTempList1) {
                for (int i = 0; i < numPocStCurrAfter && rIdx < numRpsCurrTempList1; i++) {
                    refPicListTemp1.push_back(refPicSetStCurrAfter[i]);
                    rIdx++;
                }
                for (int i = 0; i < numPocStCurrBefore && rIdx < numRpsCurrTempList1; i++) {
                    refPicListTemp1.push_back(refPicSetStCurrBefore[i]);
                    rIdx++;
                }
            }

            for (int i = 0; i <= (int)sliceHeader->num_ref_idx_l1_active_minus1; i++)
                sliceHeader->refPicList[1].push_back(refPicListTemp1[i]);
        }
    };

    sliceHeader->first_slice_segment_in_pic_flag = READ_FLAG("first_slice_segment_in_pic_flag");
    if (mNaluType >= NAL_UNIT_CODED_SLICE_BLA_W_LP && mNaluType <= NAL_UNIT_RESERVED_IRAP_VCL23)
        sliceHeader->no_output_of_prior_pics_flag = READ_FLAG("no_output_of_prior_pics_flag");
    sliceHeader->slice_pic_parameter_set_id = READ_UVLC("slice_pic_parameter_set_id");

    pps = ps->getPps(sliceHeader->slice_pic_parameter_set_id);
    sps = ps->getSps(pps->pps_seq_parameter_set_id);

    if (!sliceHeader->first_slice_segment_in_pic_flag) {
        if (pps->dependent_slice_segments_enabled_flag)
            sliceHeader->dependent_slice_segment_flag = READ_FLAG("dependent_slice_segment_flag");
        sliceHeader->slice_segment_address = READ_CODE(ceilLog2(sps->ctbCount), "slice_segment_address");
    }

    sliceHeader->slice_beta_offset_div2 = pps->pps_beta_offset_div2;
    sliceHeader->slice_tc_offset_div2 = pps->pps_tc_offset_div2;
    sliceHeader->slice_loop_filter_across_slices_enabled_flag = pps->pps_loop_filter_across_slices_enabled_flag;

    if (!sliceHeader->dependent_slice_segment_flag) {
        sliceHeader->slice_reserved_flag.resize(pps->num_extra_slice_header_bits);
        for (int i = 0; i < pps->num_extra_slice_header_bits; i++)
            sliceHeader->slice_reserved_flag[i] = READ_FLAG("slice_reserved_flag[i]");
        sliceHeader->slice_type = READ_UVLC("slice_type");
        if (pps->output_flag_present_flag) sliceHeader->pic_output_flag = READ_FLAG("pic_output_flag");
        if (sps->separate_colour_plane_flag == 1) sliceHeader->colour_plane_id = READ_CODE(2, "colour_plane_id");

        if (mNaluType != NAL_UNIT_CODED_SLICE_IDR_W_RADL && mNaluType != NAL_UNIT_CODED_SLICE_IDR_N_LP) {
            uint8_t bits = sps->log2_max_pic_order_cnt_lsb_minus4 + 4;
            sliceHeader->slice_pic_order_cnt_lsb = READ_CODE(bits, "slice_pic_order_cnt_lsb");
            parsePoc(prevTid0Poc, 1 << (sps->log2_max_pic_order_cnt_lsb_minus4 + 4));

            sliceHeader->short_term_ref_pic_set_sps_flag = READ_FLAG("short_term_ref_pic_set_sps_flag");

            if (!sliceHeader->short_term_ref_pic_set_sps_flag)
                parseStRps(&sliceHeader->st_ref_pic_set, sps->num_short_term_ref_pic_sets,
                           sps->num_long_term_ref_pics_sps);
            else if (sps->num_short_term_ref_pic_sets > 1)
                sliceHeader->short_term_ref_pic_set_idx =
                    READ_CODE(ceilLog2(sps->num_short_term_ref_pic_sets), "short_term_ref_pic_set_idx");

            if (sps->long_term_ref_pics_present_flag) {
                // TODO: parse long term ref
            }

            if (sps->sps_temporal_mvp_enabled_flag)
                sliceHeader->slice_temporal_mvp_enabled_flag = READ_FLAG("slice_temporal_mvp_enabled_flag");
        }

        if (sps->sample_adaptive_offset_enabled_flag) {
            sliceHeader->slice_sao_luma_flag = READ_FLAG("slice_sao_luma_flag");
            if (sps->chromaArrayType) sliceHeader->slice_sao_chroma_flag = READ_FLAG("slice_sao_chroma_flag");
        }

        if (sliceHeader->slice_type == P_SLICE || sliceHeader->slice_type == B_SLICE) {
            sliceHeader->num_ref_idx_active_override_flag = READ_FLAG("num_ref_idx_active_override_flag");
            if (sliceHeader->num_ref_idx_active_override_flag) {
                sliceHeader->num_ref_idx_l0_active_minus1 = READ_UVLC("num_ref_idx_l0_active_minus1");
                if (sliceHeader->slice_type == B_SLICE)
                    sliceHeader->num_ref_idx_l1_active_minus1 = READ_UVLC("num_ref_idx_l1_active_minus1");
            }

            if (!sliceHeader->short_term_ref_pic_set_sps_flag)
                parseRpsList(&sliceHeader->st_ref_pic_set, dpb);
            else if (sps->num_short_term_ref_pic_sets > 1)
                parseRpsList(&sps->st_ref_pic_set[sliceHeader->short_term_ref_pic_set_idx], dpb);

            if (pps->lists_modification_present_flag /* && NumPicTotalCurr > 1*/) {
                // TODO: parse ref_pic_lists_modification
            }

            if (sliceHeader->slice_type == B_SLICE) sliceHeader->mvd_l1_zero_flag = READ_FLAG("mvd_l1_zero_flag");
            if (pps->cabac_init_present_flag) sliceHeader->cabac_init_flag = READ_FLAG("cabac_init_flag");
            if (sliceHeader->slice_temporal_mvp_enabled_flag) {
                if (sliceHeader->slice_type == B_SLICE)
                    sliceHeader->collocated_from_l0_flag = READ_FLAG("collocated_from_l0_flag");
                if ((sliceHeader->collocated_from_l0_flag && sliceHeader->num_ref_idx_l0_active_minus1 > 0) ||
                    (!sliceHeader->collocated_from_l0_flag && sliceHeader->num_ref_idx_l1_active_minus1 > 0))
                    sliceHeader->collocated_ref_idx = READ_UVLC("collocated_ref_idx");
            }

            if ((pps->weighted_pred_flag && sliceHeader->slice_type == P_SLICE) ||
                (pps->weighted_bipred_flag && sliceHeader->slice_type == B_SLICE)) {
                auto parsePwt = [this, sps, sliceHeader](HevcPredWeightTable *table) {
                    table->luma_log2_weight_denom = READ_UVLC("luma_log2_weight_denom");
                    if (sps->chromaArrayType)
                        table->delta_chroma_log2_weight_denom = READ_SVLC("delta_chroma_log2_weight_denom");

                    table->luma_weight_l0_flag.resize(sliceHeader->num_ref_idx_l0_active_minus1 + 1, 0);
                    for (int i = 0; i <= (int)sliceHeader->num_ref_idx_l0_active_minus1; i++) {
                        if (sliceHeader->refPicList[0][i]->getTid() != mTid ||
                            sliceHeader->refPicList[0][i]->getPoc() != sliceHeader->poc)
                            table->luma_weight_l0_flag[i] = READ_FLAG("luma_weight_l0_flag[i]");
                    }

                    if (sps->chromaArrayType) {
                        table->chroma_weight_l0_flag.resize(sliceHeader->num_ref_idx_l0_active_minus1 + 1, 0);
                        for (int i = 0; i <= (int)sliceHeader->num_ref_idx_l0_active_minus1; i++) {
                            if (sliceHeader->refPicList[0][i]->getTid() != mTid ||
                                sliceHeader->refPicList[0][i]->getPoc() != sliceHeader->poc)
                                table->chroma_weight_l0_flag[i] = READ_FLAG("chroma_weight_l0_flag[i]");
                        }
                    }

                    // TODO: skip delta_luma_weight_l0 luma_offset_l0 delta_chroma_weight_l0 delta_chroma_offset_l0

                    if (sliceHeader->slice_type == B_SLICE) {
                        table->luma_weight_l1_flag.resize(sliceHeader->num_ref_idx_l1_active_minus1 + 1, 0);
                        for (int i = 0; i <= (int)sliceHeader->num_ref_idx_l1_active_minus1; i++) {
                            if (sliceHeader->refPicList[0][i]->getTid() != mTid ||
                                sliceHeader->refPicList[1][i]->getPoc() != sliceHeader->poc)
                                table->luma_weight_l1_flag[i] = READ_FLAG("luma_weight_l1_flag[i]");
                        }

                        if (sps->chromaArrayType) {
                            table->chroma_weight_l1_flag.resize(sliceHeader->num_ref_idx_l1_active_minus1 + 1, 0);
                            for (int i = 0; i <= (int)sliceHeader->num_ref_idx_l1_active_minus1; i++) {
                                if (sliceHeader->refPicList[0][i]->getTid() != mTid ||
                                    sliceHeader->refPicList[1][i]->getPoc() != sliceHeader->poc)
                                    table->chroma_weight_l1_flag[i] = READ_FLAG("chroma_weight_l1_flag[i]");
                            }
                        }
                    }

                    // TODO: skip delta_luma_weight_l1 luma_offset_l1 delta_chroma_weight_l1 delta_chroma_offset_l1
                };
                parsePwt(&sliceHeader->pred_weight_table);
            }

            sliceHeader->five_minus_max_num_merge_cand = READ_UVLC("five_minus_max_num_merge_cand");
            // TODO: motion_vector_resolution_control_idc in sps_scc_extension
        }

        sliceHeader->slice_qp_delta = READ_SVLC("slice_qp_delta");
        if (pps->pps_slice_chroma_qp_offsets_present_flag) {
            sliceHeader->slice_cb_qp_offset = READ_SVLC("slice_cb_qp_offset");
            sliceHeader->slice_cr_qp_offset = READ_SVLC("slice_cr_qp_offset");
        }
        // TODO: pps_slice_act_qp_offsets_present_flag in pps_scc_extension
        // TODO: chroma_qp_offset_list_enabled_flag in pps_range_extension

        if (pps->deblocking_filter_override_enabled_flag)
            sliceHeader->deblocking_filter_override_flag = READ_FLAG("deblocking_filter_override_flag");
        if (sliceHeader->deblocking_filter_override_flag) {
            sliceHeader->slice_deblocking_filter_disabled_flag = READ_FLAG("slice_deblocking_filter_disabled_flag");
            if (!sliceHeader->slice_deblocking_filter_disabled_flag) {
                sliceHeader->slice_beta_offset_div2 = READ_SVLC("slice_beta_offset_div2");
                sliceHeader->slice_tc_offset_div2 = READ_SVLC("slice_tc_offset_div2");
            }
        }

        if (pps->pps_loop_filter_across_slices_enabled_flag &&
            (sliceHeader->slice_sao_luma_flag || sliceHeader->slice_sao_chroma_flag ||
             !sliceHeader->slice_deblocking_filter_disabled_flag))
            sliceHeader->slice_loop_filter_across_slices_enabled_flag =
                READ_FLAG("slice_loop_filter_across_slices_enabled_flag");
    }

    if (pps->tiles_enabled_flag || pps->entropy_coding_sync_enabled_flag) {
        sliceHeader->num_entry_point_offsets = READ_UVLC("num_entry_point_offsets");
        if (sliceHeader->num_entry_point_offsets > 0) {
            sliceHeader->offset_len_minus1 = READ_UVLC("offset_len_minus1");
            uint8_t bits = sliceHeader->offset_len_minus1 + 1;
            sliceHeader->entry_point_offset_minus1.resize(sliceHeader->num_entry_point_offsets);
            for (int i = 0; i < (int)sliceHeader->num_entry_point_offsets; i++)
                sliceHeader->entry_point_offset_minus1[i] = READ_CODE(bits, "entry_point_offset_minus1[i]");
        }
    }

    if (pps->slice_segment_header_extension_present_flag) {
        sliceHeader->slice_segment_header_extension_length = READ_UVLC("slice_segment_header_extension_length");
        sliceHeader->slice_segment_header_extension_data_byte.resize(
            sliceHeader->slice_segment_header_extension_length);
        for (int i = 0; i < (int)sliceHeader->slice_segment_header_extension_length; i++)
            sliceHeader->slice_segment_header_extension_data_byte[i] =
                READ_CODE(8, "slice_segment_header_extension_data_byte[i]");
    }

    // byte alignment
    READ_FLAG("alignment_bit_equal_to_one");
    uint8_t bits = numBitsUntilByteAligned();
    skipBits(bits);

    // calculated value
    sliceHeader->sliceQp = pps->init_qp_minus26 + sliceHeader->slice_qp_delta + 26;
    sliceHeader->maxNumMergeCand = 5 - sliceHeader->five_minus_max_num_merge_cand;
    ps->setSliceHeader(sliceHeader);
}

void CavlcReader::parsePtl(HevcProfileTierLevel *ptl, bool profilePresentFlag, uint8_t maxNumSubLayersMinus1) {
#define PTL_TXT(txt) (isSubLayer ? ("sub_layer_" txt) : ("general_" txt))

    auto parseCommonPtl = [this](HevcCommonProfileTierLevel *commonPtl, bool isSubLayer) {
        commonPtl->profile_space = READ_CODE(2, PTL_TEXT("profile_space"));
        commonPtl->tier_flag = READ_FLAG(PTL_TEXT("tier_flag"));
        commonPtl->profile_idc = READ_CODE(5, PTL_TEXT("profile_idc"));
        for (int j = 0; j < 32; j++)
            commonPtl->profile_compatibility_flag[j] = READ_FLAG(PTL_TEXT("profile_compatibility_flag[j]"));
        commonPtl->progressive_source_flag = READ_FLAG(PTL_TEXT("progressive_source_flag"));
        commonPtl->interlaced_source_flag = READ_FLAG(PTL_TEXT("interlaced_source_flag"));
        commonPtl->non_packed_constraint_flag = READ_FLAG(PTL_TEXT("non_packed_constraint_flag"));
        commonPtl->frame_only_constraint_flag = READ_FLAG(PTL_TEXT("frame_only_constraint_flag"));
        if (commonPtl->profile_idc == Profile::MAINREXT || commonPtl->profile_compatibility_flag[Profile::MAINREXT] ||
            commonPtl->profile_idc == Profile::HIGHTHROUGHPUTREXT ||
            commonPtl->profile_compatibility_flag[Profile::HIGHTHROUGHPUTREXT]) {
            commonPtl->max_12bit_constraint_flag = READ_FLAG(PTL_TEXT("max_12bit_constraint_flag"));
            commonPtl->max_10bit_constraint_flag = READ_FLAG(PTL_TEXT("max_10bit_constraint_flag"));
            commonPtl->max_8bit_constraint_flag = READ_FLAG(PTL_TEXT("max_8bit_constraint_flag"));
            commonPtl->max_422chroma_constraint_flag = READ_FLAG(PTL_TEXT("max_422chroma_constraint_flag"));
            commonPtl->max_420chroma_constraint_flag = READ_FLAG(PTL_TEXT("max_420chroma_constraint_flag"));
            commonPtl->max_monochrome_constraint_flag = READ_FLAG(PTL_TEXT("max_monochrome_constraint_flag"));
            commonPtl->intra_constraint_flag = READ_FLAG(PTL_TEXT("intra_constraint_flag"));
            commonPtl->one_picture_only_constraint_flag = READ_FLAG(PTL_TEXT("one_picture_only_constraint_flag"));
            commonPtl->lower_bit_rate_constraint_flag = READ_FLAG(PTL_TEXT("lower_bit_rate_constraint_flag"));
            READ_CODE(16, PTL_TEXT("reserved_zero_34bits[0..15]"));
            READ_CODE(16, PTL_TEXT("reserved_zero_34bits[16..31]"));
            READ_CODE(2, PTL_TEXT("reserved_zero_34bits[32..33]"));
        } else if (commonPtl->profile_idc == Profile::MAIN10 ||
                   commonPtl->profile_compatibility_flag[Profile::MAIN10]) {
            READ_CODE(7, PTL_TEXT("reserved_zero_7bits"));
            commonPtl->one_picture_only_constraint_flag = READ_FLAG(PTL_TEXT("one_picture_only_constraint_flag"));
            READ_CODE(16, PTL_TEXT("reserved_zero_35bits[0..15]"));
            READ_CODE(16, PTL_TEXT("reserved_zero_35bits[16..31]"));
            READ_CODE(3, PTL_TEXT("reserved_zero_35bits[32..34]"));
        } else {
            READ_CODE(16, PTL_TEXT("reserved_zero_43bits[0..15]"));
            READ_CODE(16, PTL_TEXT("reserved_zero_43bits[16..31]"));
            READ_CODE(11, PTL_TEXT("reserved_zero_43bits[32..42]"));
        }

        if (commonPtl->profile_idc == Profile::MAIN || commonPtl->profile_compatibility_flag[Profile::MAIN] ||
            commonPtl->profile_idc == Profile::MAIN10 || commonPtl->profile_compatibility_flag[Profile::MAIN10] ||
            commonPtl->profile_idc == Profile::MAINSTILLPICTURE ||
            commonPtl->profile_compatibility_flag[Profile::MAINSTILLPICTURE] ||
            commonPtl->profile_idc == Profile::MAINREXT || commonPtl->profile_compatibility_flag[Profile::MAINREXT] ||
            commonPtl->profile_idc == Profile::HIGHTHROUGHPUTREXT ||
            commonPtl->profile_compatibility_flag[Profile::HIGHTHROUGHPUTREXT]) {
            commonPtl->inbld_flag = READ_FLAG(PTL_TEXT("inbld_flag"));
        } else
            READ_FLAG(PTL_TEXT("reserved_zero_bit"));
    };

    if (profilePresentFlag) parseCommonPtl(&ptl->general_profile_tier_level, false);
    ptl->general_level_idc = READ_CODE(8, "general_level_idc");

    ptl->sub_layer_profile_present_flag.resize(maxNumSubLayersMinus1);
    ptl->sub_layer_level_present_flag.resize(maxNumSubLayersMinus1);
    for (int i = 0; i < maxNumSubLayersMinus1; i++) {
        ptl->sub_layer_profile_present_flag[i] = READ_FLAG("sub_layer_profile_present_flag[i]");
        ptl->sub_layer_level_present_flag[i] = READ_FLAG("sub_layer_level_present_flag[i]");
    }

    if (maxNumSubLayersMinus1 > 0)
        for (int i = maxNumSubLayersMinus1; i < 8; i++) READ_CODE(2, "reserved_zero_2bits[i]");

    ptl->sub_layer_profile_tier_level.resize(maxNumSubLayersMinus1);
    ptl->sub_layer_level_idc.resize(maxNumSubLayersMinus1);
    for (int i = 0; i < maxNumSubLayersMinus1; i++) {
        if (ptl->sub_layer_profile_present_flag[i]) parseCommonPtl(&ptl->sub_layer_profile_tier_level[i], true);
        if (ptl->sub_layer_level_present_flag[i]) ptl->sub_layer_level_idc[i] = READ_CODE(8, "sub_layer_level_idc[i]");
    }
}

void CavlcReader::parseStRps(HevcStRefPicSet *stRps, int stRpsIdx, int numStRps) {
    if (stRpsIdx) stRps->inter_ref_pic_set_prediction_flag = READ_FLAG("inter_ref_pic_set_prediction_flag");

    if (stRps->inter_ref_pic_set_prediction_flag) {
        if (stRpsIdx == numStRps) stRps->delta_idx_minus1 = READ_UVLC("delta_idx_minus1");
        stRps->delta_rps_sign = READ_FLAG("delta_rps_sign");
        stRps->abs_delta_rps_minus1 = READ_UVLC("abs_delta_rps_minus1");

        // TODO: parse short term rps
    } else {
        stRps->num_negative_pics = READ_UVLC("num_negative_pics");
        stRps->num_positive_pics = READ_UVLC("num_positive_pics");

        stRps->delta_poc_s0_minus1.resize(stRps->num_negative_pics);
        stRps->used_by_curr_pic_s0_flag.resize(stRps->num_negative_pics);
        for (int i = 0; i < (int)stRps->num_negative_pics; i++) {
            stRps->delta_poc_s0_minus1[i] = READ_UVLC("delta_poc_s0_minus1[i]");
            stRps->used_by_curr_pic_s0_flag[i] = READ_FLAG("used_by_curr_pic_s0_flag[i]");
        }

        stRps->delta_poc_s1_minus1.resize(stRps->num_positive_pics);
        stRps->used_by_curr_pic_s1_flag.resize(stRps->num_positive_pics);
        for (int i = 0; i < (int)stRps->num_positive_pics; i++) {
            stRps->delta_poc_s1_minus1[i] = READ_UVLC("delta_poc_s1_minus1[i]");
            stRps->used_by_curr_pic_s1_flag[i] = READ_FLAG("used_by_curr_pic_s1_flag[i]");
        }
    }

    // calculated
    stRps->numNegativePics = stRps->num_negative_pics;
    stRps->numPositivePics = stRps->num_positive_pics;
    stRps->numDeltaPocs = stRps->numNegativePics + stRps->numPositivePics;

    stRps->usedByCurrPicS0.resize(stRps->numNegativePics);
    stRps->deltaPocS0.resize(stRps->numNegativePics);
    for (int i = 0; i < (int)stRps->numNegativePics; i++) {
        stRps->usedByCurrPicS0[i] = stRps->used_by_curr_pic_s0_flag[i];
        if (i == 0)
            stRps->deltaPocS0[i] = -(int)(stRps->delta_poc_s0_minus1[i] + 1);
        else
            stRps->deltaPocS0[i] = stRps->deltaPocS0[i - 1] - (stRps->delta_poc_s0_minus1[i] + 1);
    }

    stRps->usedByCurrPicS1.resize(stRps->numPositivePics);
    stRps->deltaPocS1.resize(stRps->numPositivePics);
    for (int i = 0; i < (int)stRps->numPositivePics; i++) {
        stRps->usedByCurrPicS1[i] = stRps->used_by_curr_pic_s1_flag[i];
        if (i == 0)
            stRps->deltaPocS1[i] = (stRps->delta_poc_s1_minus1[i] + 1);
        else
            stRps->deltaPocS1[i] = stRps->deltaPocS1[i - 1] + (stRps->delta_poc_s1_minus1[i] + 1);
    }
}
