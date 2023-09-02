
#ifndef _KFAPP_H
#define _KFAPP_H

#include "rtklib.h"
#include "CGNSSManage.h"
#include "PSINS.h"

class CMyAlignment : public CSINSGNSS       
{
public: 
    enum eMode {
        STATIC = 1,
        STATIC_FINE = 2,
        KINEMA = 3
    };

    BOOL alnOK;
    CQuat  qnb;                             

    CMyAlignment(int mode, double T1 = 100.0, double T2 = 100.0); 
    ~CMyAlignment();
    void SetParam(double vel00, double wz00, double dyaw00);            
    void Init(const CSINS &sins0);  
    void SetMeasGNSS(const CVect3 &pos, const CVect3 &vn);         
    int  Update(const CVect3 *pwm, const CVect3 *pvm, int nn, double ts, int nSteps = 5);  

private:
    BOOL posInit, alnkfInit;
    CAligni0 alni0;
    int alnMode, cntYaw;
    double T1, T2, velPre, yawPre, vel0, wz0, dyaw0;       
};

class CMyNavFilter : public CSINSGNSS        
{
public:              
    CVect3 pos_xyz;
    CMat3  Rk_xyz;

    CMyNavFilter(int nq0, int nr0, double ts, int yawHkRow0 = 6);
    ~CMyNavFilter();
    void Loadopt(const imuopt_t &opt);
    void Init(const CSINS &sins0);
    void ZUPTtest(void);
    void ZIHRtest(void);
    void NHCtest(void);         /* 未启用 */
    void SetMeasGNSS(const CVect3 &pos, const CVect3 &vn, double qfactor = 1.0);   
    int  Update(const CVect3 *pwm, const CVect3 *pvm, int nn, double ts, int nSteps = 5);
	void operator>>(CFileRdWt &f);
    void OutINSResult(double dt);

private:
    CWzhold wzhd;
    CVect3 P0_phi, P0_vel, P0_pos, P0_gyr, P0_acc;         
    CVect3 Qt_phi, Qt_vel, Qt_pos, Qt_gyr, Qt_acc;   
};

class CMyAutoDrive
{  
public:
    enum eState {
        INITIALIZE = 0,
        ALIGNMENT  = 1,
        NAVIGATION = 2
    };
    int bIMU2GPS;   
    double tk, tgps;
    eState mState;

    CMyAutoDrive(); 
    ~CMyAutoDrive();                       
    int  Init(const prcopt_t &popt, const solopt_t &sopt, const filopt_t &fopt, const sta_t &sta, const imuopt_t &iopt);
    void GPSProcess(obsd_t *obs, int nobs); 
    void SetMeasGNSS(const CVect3 &pos, const CVect3 &vn, double dt);      
    void IMUProcess(imud_t *imu); 

    void TraceOpen(const char *filename);
    void Sol2kml();

    friend void SetInterrupt(CMyAutoDrive &app, double T1, double T2);
    friend int  getSamples(CMyAutoDrive &app); 
 
private:
    int week, nn, tdStep;         
    double ts, T1, T2;  

    CFileRdWt    *mpFile;
    CGNSSManage  *mpGPS;            
    CMyAlignment *mpAln;   
    CMyNavFilter *mpFIR;                   
};  

int  getSamples(CMyAutoDrive &app);
void SetInterrupt(CMyAutoDrive &app, double T1, double T2);
void printbrk(gtime_t ti, int gnss_stat, int gins_stat);

#endif