#include "ParameterSet.h"

#include <assert.h>

ParameterSet::ParameterSet() {
    mVpsList.resize(HEVC_MAX_VPS_COUNT);
    mSpsList.resize(HEVC_MAX_SPS_COUNT);
    mPpsList.resize(HEVC_MAX_PPS_COUNT);
}

ParameterSet::~ParameterSet() {}

void ParameterSet::setVps(std::shared_ptr<HevcVps> vps) {
    assert(vps->vps_video_parameter_set_id < HEVC_MAX_VPS_COUNT);
    mVpsList[vps->vps_video_parameter_set_id] = vps;
}

void ParameterSet::setSps(std::shared_ptr<HevcSps> sps) {
    assert(sps->sps_seq_parameter_set_id < HEVC_MAX_SPS_COUNT);
    mSpsList[sps->sps_seq_parameter_set_id] = sps;
}

void ParameterSet::setPps(std::shared_ptr<HevcPps> pps) {
    assert(pps->pps_pic_parameter_set_id < HEVC_MAX_PPS_COUNT);
    mPpsList[pps->pps_pic_parameter_set_id] = pps;
}

void ParameterSet::setSliceHeader(std::shared_ptr<HevcSliceHeader> sliceHeader) {
    mSliceHeader = sliceHeader;
}

std::shared_ptr<HevcVps> ParameterSet::getVps(uint8_t id) {
    assert(id < HEVC_MAX_VPS_COUNT);
    return mVpsList[id];
}

std::shared_ptr<HevcSps> ParameterSet::getSps(uint8_t id) {
    assert(id < HEVC_MAX_SPS_COUNT);
    return mSpsList[id];
}

std::shared_ptr<HevcPps> ParameterSet::getPps(uint8_t id) {
    assert(id < HEVC_MAX_PPS_COUNT);
    return mPpsList[id];
}

std::shared_ptr<HevcSliceHeader> ParameterSet::getSliceHeader() {
    return mSliceHeader;
}
