#ifndef _IMUCONFIG_H
#define _IMUCONFIG_H

#include <yaml-cpp/yaml.h>
#include "PSINS.h"

using namespace std;
class CIMUConfig 
{
public:
    imuopt_t imuopt;

    CIMUConfig();
    int LoadYAML(const char *file);
};

#endif