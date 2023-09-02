#ifndef _GNSSCONFIG_H
#define _GNSSCONFIG_H

#include <yaml-cpp/yaml.h>

#include "rtklib.h"
#include <memory>

#ifdef __cplusplus
extern "C" {
#endif

#undef lock
#undef unlock

#ifdef __cplusplus
}
#endif

class CGNSSConfig 
{
public:
    typedef std::shared_ptr<CGNSSConfig> Ptr;
    prcopt_t prcopt = prcopt_default;
    solopt_t solopt = solopt_default;
    filopt_t filopt = filopt_default;
    sta_t sta = {};

    CGNSSConfig();
    int LoadYAML(const char *file);
};

#endif