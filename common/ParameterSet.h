#ifndef PARAMETER_SET_H
#define PARAMETER_SET_H

#include "Typedef.h"

#include <memory>
#include <vector>

class ParameterSet {
private:
    std::vector<std::shared_ptr<HevcVps>> mVpsList;
    std::vector<std::shared_ptr<HevcSps>> mSpsList;
    std::vector<std::shared_ptr<HevcPps>> mPpsList;

    std::shared_ptr<HevcSliceHeader> mSliceHeader;

public:
    ParameterSet();
    ~ParameterSet();

    void setVps(std::shared_ptr<HevcVps> vps);
    void setSps(std::shared_ptr<HevcSps> sps);
    void setPps(std::shared_ptr<HevcPps> pps);
    void setSliceHeader(std::shared_ptr<HevcSliceHeader> sliceHeader);

    std::shared_ptr<HevcVps> getVps(uint8_t id);
    std::shared_ptr<HevcSps> getSps(uint8_t id);
    std::shared_ptr<HevcPps> getPps(uint8_t id);
    std::shared_ptr<HevcSliceHeader> getSliceHeader();
};

#endif
