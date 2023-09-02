#ifndef _GNSSMANAGE_H
#define _GNSSMANAGE_H

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

struct sol {
    double ti;
    double rr[6];
    double qr[6];
};
struct sols {
    int n, nmax;
    sol *data;
};

class CGNSSManage      
{
public:
    typedef std::shared_ptr<CGNSSManage> Ptr;
    uint8_t *stat;
    gtime_t *ti;
    char   outdir[1024] = "";
    double pos_blh[3], vel_enu[3], rk_blh[3], rk_vel[3];       
               
    CGNSSManage(); 
    ~CGNSSManage();
    int  Init(const prcopt_t &popt, const solopt_t &sopt, const filopt_t &fopt, const sta_t &sta); 
    void Update(obsd_t *obs, int n, int mode = 0, const double *ins_xyz = NULL, const double *ins_Pxyz = NULL);
    void TraceOpen(const char *file);
    void Sol2kml();
    void OutGPSResult();
    void Update(gtime_t ti);

private:
    sols  solss = {};
    rtk_t rtk = {};                  
    nav_t nav = {};
    solopt_t solopt = solopt_default;
    FILE* fp_outs[2] = {};
};

void removeUnusedData(gtime_t &ts, gtime_t &te, nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr);

#endif