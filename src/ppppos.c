/****************************************
GNSS?????????||???????PPP???????||???
*****************************************/
#include "rtklib.h"
#include "ppp_state.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define MAX_UINT8  1000				 /* max count of uint8 */
#define MAX_ITER     8               /* max number of iterations */
#define MIN_NSAT_SOL 7               // 7

#define VAR_GRA     SQR(0.01)        /* init variance gradient (m^2) */    	 		

typedef struct {
	unsigned char n;
	int flg[256];
} exc_t;
static char id[8];

static int sys2num(int sys) {
	switch (sys) {
		case SYS_GPS: return 0;
		case SYS_GLO: return 1;
		case SYS_GAL: return 2;
		case SYS_CMP: return 3;
		default: return -1;
	}
}
static void readvflg(int vflg, int *type, int *freq, int *sat1, int *sat2) 
{
	*type = vflg & 0xF;	
	if (*type & 3) {
		*sat1 = (vflg >> 8) & 0xFF;		*freq = (vflg >> 4) & 0xF;
	}
	if (*type & 8) {
		*sat1 = (vflg >> 8) & 0xFF;		*sat2 = (vflg >> 16) & 0xFF;
	}	
}
static void setvflg(int *vflg, int type, int freq, int sat1, int sat2) 
{
	*vflg = 0x0000;
	if (type & 3) {
		*vflg = (sat1 << 8) | (freq << 4) | type;
	}
	if (type & 4) {
		*vflg = type;
	}
	if (type & 8) {
		*vflg = (sat1 << 16) | (sat2 << 8) | type;
	}
}
/* number of estimated states ------------------------------------------------*/
extern int pppnx(const prcopt_t *opt)
{
	return NX(opt);
}

/* standard deviation of state -----------------------------------------------*/
static double STD(const rtk_t *rtk, int i)
{
    if (rtk->sol.stat==SOLQ_FIX) return SQRT(rtk->Pa[i+i*rtk->nx]);
    return SQRT(rtk->P[i+i*rtk->nx]);
}

/* write solution status ---------------------------------------------------- */
static void pppoutstat(const rtk_t *rtk, FILE *fp) 
{
	int i, week;
    double tow, pos[3], vel[3], *x, zhd, zazel[] = { 0.0, PI / 2.0 };
    
    tow=time2gpst(rtk->sol.time,&week);
    
    x = rtk->sol.stat==SOLQ_FIX?rtk->xa:rtk->x;

    /* receiver position */
	ecef2pos(rtk->sol.rr,pos);
    fprintf(fp,"$POS,%d,%.0f,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%.4f\n",week,tow,
            rtk->sol.stat,rtk->sol.ns,x[0],x[1],x[2],STD(rtk,0),STD(rtk,1),STD(rtk,2),
			rtk->sol.nsat, rtk->sol.pdop);		/* 2022-08-17 add by zh */
    
    /* receiver velocity and acceleration */
	if (rtk->opt.mode == PMODE_PPP_KINEMA) {
		ecef2enu(pos,rtk->sol.rr+3,vel);
		fprintf(fp,"$VELACC,%d,%.0f,%d,%.4f,%.4f,%.4f\n",week,tow,
			rtk->sol.stat,vel[0],vel[1],vel[2]);
	}
    
    /* receiver clocks */
    i=IC(0,&rtk->opt);
    fprintf(fp,"$CLK,%d,%.0f,%d,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
            week,tow,rtk->sol.stat,x[i]*1E9/CLIGHT,x[i+2]*1E9/CLIGHT,x[i+3]*1E9/CLIGHT,
            STD(rtk,i)*1E9/CLIGHT, STD(rtk,i+2)*1E9/CLIGHT, STD(rtk,i+3)*1E9/CLIGHT);
    
    /* tropospheric parameters */
    if (rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG) {
        i = IT(&rtk->opt);
		zhd = tropmodel(rtk->sol.time, pos, zazel, 1, REL_HUMI, 1);
        fprintf(fp,"$TRP,%d,%.0f,%d,%.4f,%.4f,%.4f,%.4f\n",week,tow,
		        rtk->sol.stat,x[i],STD(rtk,i),zhd,SQR(0.1));
    }
}
static void pppoutiono(const rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, FILE *fp) 
{
	const ssat_t *ssat;
	int i, ll, week;
    double tow, *x;

	tow=time2gpst(rtk->sol.time,&week);

	x = rtk->sol.stat==SOLQ_FIX?rtk->xa:rtk->x;

	if (rtk->opt.ionoopt==IONOOPT_EST) {
		for (i=0; i<n; i++) {
			ssat=rtk->ssat + obs[i].sat - 1;
			ll = II(obs[i].sat, &rtk->opt);
			if (!ssat->vs || rtk->x[ll]==0.0) continue; 
	
			if ((rtk->opt.modear==6 && ssat->fix[0]==2    && ssat->fix[1]==2) || 
				(rtk->opt.modear==0 && ssat->vsat[0]&0x2 && ssat->vsat[1]&0x2)) {
				satno2id(obs[i].sat, id);
				fprintf(fp,"$ION,%d,%.0f,%d,%s,%d,%d,%.1f,%.1f,%.4f,%.4f\n",week,tow,
                    rtk->sol.stat,id,ssat->lock[0],ssat->lock[1],ssat->azel[0]*R2D,
                    ssat->azel[1]*R2D,x[ll],STD(rtk,ll));
			}
		}
    }
}
static void pppoutambi(const rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, FILE *fp) 
{
	const ssat_t *ssat;
	int i, j, ll, mm, nn, week;
	double tow, a, b, lam1, lam2, lam3, EL, WL, NL, IF, Pe, Pw, Pl, Pf, freq[NFREQ], *x, *P;

	tow=time2gpst(rtk->sol.time,&week);
    
    x = rtk->sol.stat==SOLQ_FIX?rtk->xa:rtk->x;
    P = rtk->sol.stat==SOLQ_FIX?rtk->Pa:rtk->P;

	for (i = 0; i < n; i++) {						
		ssat = rtk->ssat + obs[i].sat - 1;
		if (!ssat->vs) continue;

		EL = WL = NL = IF = 0.0;			Pe = Pw = Pl = Pf = 1e6;
		ll = IB(obs[i].sat, 0, &rtk->opt);
		mm = IB(obs[i].sat, 1, &rtk->opt);

		for (j = 0; j < NFREQ; j++) freq[j] = sat2freq(obs[i].sat, obs[i].code[j], nav);
		lam1 = CLIGHT / freq[0];
		lam2 = CLIGHT / freq[1];
		lam3 = CLIGHT / freq[2];

		//el ambiguity
		if ((ssat->vsat[1] & 0x2) && (ssat->vsat[2] & 0x2)) {
			nn = IB(obs[i].sat, 2, &rtk->opt);
			EL = x[mm] / lam2 - x[nn] / lam3;
			Pe = SQR(1 / lam2)*P[mm + mm*rtk->nx] + SQR(1 / lam3)*P[nn + nn*rtk->nx] - 
				 2 / (lam2*lam3)*P[nn + mm*rtk->nx];
		}
		//wl&if ambiguity
		if ((ssat->vsat[0] & 0x2) && (ssat->vsat[1] & 0x2)) {
			WL = x[ll] / lam1 - x[mm] / lam2;
			Pw = SQR(1 / lam1)*P[ll + ll*rtk->nx] + SQR(1 / lam2)*P[mm + mm*rtk->nx] - 
				 2 / (lam1*lam2)*P[mm + ll*rtk->nx];

			a = SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));
			b = SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));
			IF = a*x[ll] - b*x[mm];
			Pf = a*a*P[ll + ll*rtk->nx] + b*b*P[mm + mm*rtk->nx] - 2*a*b*P[mm + ll*rtk->nx];
		}
		satno2id(obs[i].sat, id);
		fprintf(fp,"$AMB,%d,%.0f,%d,%s,%d,%.2f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d,%d,%d\n",week,tow, 
			rtk->sol.stat,id,(ssat->ini>>4)|(ssat->reset>>4),ssat->azel[1]*R2D,EL,sqrt(Pe),WL,sqrt(Pw),
			IF,sqrt(Pf),ssat->lock[0],ssat->lock[1],ssat->lock[2]);
	}
}
static void pppoutresi(const rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, FILE *fp) 		
{
	const ssat_t *ssat;
	int i, week;	
	double tow;		
	
	tow=time2gpst(rtk->sol.time,&week);
	for (i = 0; i < n; i++) {
		ssat = rtk->ssat + obs[i].sat - 1;
		if (!ssat->vs) continue;
		satno2id(obs[i].sat, id);
		if (0) {
			fprintf(fp,"$RES,%d,%.0f,%d,%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",week,tow,rtk->sol.stat,id, 
				ssat->vsat[0]&0x1?ssat->resp[0]:1E+03,ssat->vsat[1]&0x1?ssat->resp[1]:1E+03,ssat->vsat[2]&0x1?ssat->resp[2]:1E+03,
				ssat->vsat[0]&0x2?ssat->resc[0]:1E+03,ssat->vsat[1]&0x2?ssat->resc[1]:1E+03,ssat->vsat[2]&0x2?ssat->resc[2]:1E+03);
		} else {
			fprintf(fp,"$RES,%d,%.0f,%d,%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",week,tow,rtk->sol.stat,id, 
				ssat->resp[0],ssat->resp[1],ssat->resp[2],ssat->resc[0],ssat->resc[1],ssat->resc[2]);
		}
	}
}
static void pppoutqual(const rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, FILE *fp) 
{
	const ssat_t *ssat;
	int i, week;	
	double tow;	

	tow=time2gpst(rtk->sol.time,&week);
	for (i = 0; i < n; i++) {
		ssat = rtk->ssat + obs[i].sat - 1;
		if (!ssat->vs) continue;
		satno2id(obs[i].sat, id);						/* 25 */
		fprintf(fp,"$QUL,%d,%.0f,%d,%s,%.2f,%d,"		/* 周、秒、解状态、卫星、高度角、周跳 6 */
			"%.4f,%.4f,%.4f,"							/* SNR1、SNR2、SNR3 3 */
			"%.3f,%.3f,%.3f,"							/* p1-p2 3 */
			"%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,"			/* gf12、gf23、gf13 6 */
			"%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,"			/* mw12、mw23、mw13 6 */
			"%.3f,%.3f,%.3f,"							/* tdcp resi 3 */
			"%.3f,%.3f,%.3f,"							/* MP1、MP2、MP3 */
			"%.3f,%.3f\n",								/* spp后验、多普勒辅助伪距探测 */						
			week,tow,rtk->sol.stat,id,ssat->azel[1]*R2D,(ssat->ini>>4)|(ssat->reset>>4),
			obs[i].SNR[0]*SNR_UNIT,obs[i].SNR[1]*SNR_UNIT,obs[i].SNR[2]*SNR_UNIT,
			ssat->diff[9],ssat->diff[10],ssat->diff[11],
			ssat->LC[0],ssat->diff[0],ssat->LC[1],ssat->diff[1],ssat->LC[2],ssat->diff[2],
			ssat->LC[3],ssat->diff[3],ssat->LC[4],ssat->diff[4],ssat->LC[5],ssat->diff[5],
			ssat->diff[6],ssat->diff[7],ssat->diff[8],  
			ssat->LC[6],ssat->LC[7],ssat->LC[8],
			ssat->respp, ssat->diff[12]);
	}
}
static void pppoutddif(const rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, FILE *fp) 
{
	static double  DIF0[MAXSAT] = { 0 };
	static uint8_t flag[MAXSAT] = { 0 };

	const ssat_t *ssat;
	int i, week, isat;
	uint8_t vf[MAXSAT] = {0};
	double tow;
	tow=time2gpst(rtk->sol.time,&week);
	
	for (i = 0; i < n; i++) {
		isat = obs[i].sat - 1;
		ssat = rtk->ssat + isat;
		if ((ssat->vs) && !(ssat->slip[0]|ssat->slip[1]|ssat->slip[2]) && !ISZERO(ssat->LC[9])) {
			if (flag[isat]) {
				satno2id(obs[i].sat, id);
				fprintf(fp,"$DIF,%d,%.0f,%d,%s,%.2f,%d,%.6f\n",		/* 周、秒、解状态、卫星、高度角、周跳 DDIF 7 */
						week,tow,rtk->sol.stat,id,ssat->azel[1]*R2D,(ssat->slip[0]|ssat->slip[1]|ssat->slip[2]), ssat->LC[9] - DIF0[isat]);
			}
			DIF0[isat] = ssat->LC[9];
			vf[isat] = 1;
		}
	}
	for (i = 0; i < MAXSAT; i++) flag[i] = vf[i];
}
/* output solution status ----------------------------------------------------*/
extern void outppp(FILE *fp, const rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    if (!fp || !rtk->opt.outopt) return;			
	if (rtk->opt.outopt & 0x01) pppoutstat(rtk, fp); 
	if (rtk->opt.outopt & 0x02) pppoutiono(rtk, obs, n, nav, fp);
	if (rtk->opt.outopt & 0x04) pppoutambi(rtk, obs, n, nav, fp);
	if (rtk->opt.outopt & 0x08) pppoutresi(rtk, obs, n, nav, fp);
	if (rtk->opt.outopt & 0x10) pppoutqual(rtk, obs, n, nav, fp);
	if (rtk->opt.outopt & 0x20) pppoutddif(rtk, obs, n, nav, fp);
}
/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
	int j;
	rtk->x[i] = xi;
	for (j = 0; j<rtk->nx; j++) {
		rtk->P[i + j*rtk->nx] = rtk->P[j + i*rtk->nx] = i == j ? var : 0.0;
	}
}
/* exclude meas of eclipsing satellite (block IIA) ---------------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
	double rsun[3], esun[3], r, ang, erpv[5] = { 0 }, cosa;
	int i, j;
	const char *type;

	trace(3, "testeclipse:\n");

	/* unit vector of sun direction (ecef) */
	sunmoonpos(gpst2utc(obs[0].time), erpv, rsun, NULL, NULL);
	normv3(rsun, esun);

	for (i = 0; i<n; i++) {
		type = nav->pcvs[obs[i].sat - 1].type;

		if ((r = norm(rs + i * 6, 3)) <= 0.0) continue;

		/* only block IIA */
		if (*type&&!strstr(type, "BLOCK IIA")) continue;

		/* sun-earth-satellite angle */
		cosa = dot(rs + i * 6, esun, 3) / r;
		cosa = cosa<-1.0 ? -1.0 : (cosa>1.0 ? 1.0 : cosa);
		ang = acos(cosa);

		/* test eclipse */
		if (ang<PI / 2.0 || r*sin(ang)>RE_WGS84) continue;

		trace(3, "eclipsing sat excluded %s sat=%2d\n", time_str(obs[0].time, 0),
			obs[i].sat);

		for (j = 0; j<3; j++) rs[j + i * 6] = 0.0;
	}
}
/* initial rtk ssat_t ------------------------------------------------------ */
static void initial_ssat(const obsd_t *obs, int n, const prcopt_t *opt, double dt, 
	uint8_t *reset, ssat_t *ssat) 
{																																	
	ssat_t *psat;
    int i, k, nf = opt->nf;

	/* reinit pos only kinematic mode */
	if (fabs(dt) > opt->maxtdiff && opt->dynamics <= 1) (*reset) = 1;		

	for (i = 0; i < MAXSAT; i++) {			
		psat = ssat + i; 	

		for (k = 0; k < nf; k++) {
			psat->outc[k] += (int)floor(fabs(dt/opt->ti)+0.5);
			psat->outc[k]  = (psat->outc[k] > MAX_UINT8) ? MAX_UINT8 : psat->outc[k]; 	
			psat->resp[k] = psat->resc[k] = 1E03;	
		}

		/* ion */
		if (psat->outc[0] >= opt->maxout[1]) {
			psat->reset |= 0x71;
		}
		
		/* amb */
		for (k = 0; k < nf; k++) 
		if (psat->outc[k] >= opt->maxout[3]) {			
			psat->reset |= 0x10 << k;
		}
	}

	for (i = 0; i < n && i < MAXOBS; i++) {
		psat = ssat + obs[i].sat - 1;

		psat->ini = 0;	
		// psat->ini = psat->eobs = 0;	
		for (k = 0; k < nf; k++) { 
			psat->vsat[k] &= 0x00; 			
			psat->slip[k] &= 0xFC;
			psat->fix[k]   = 0;	
		}
		
		/* ion */
		if (psat->rejp[0] >= opt->maxrej[1] && psat->rejp[1] >= opt->maxrej[1]) {    // modified by zh 2022-09-12
			psat->reset |= 0x71;
		}  
		/* wl */
		if (psat->lock[0] >= 10 && psat->rejw >= 3) {		/* 初始化后要坚持一段时间 */
			psat->reset |= 0x71;
		}

		/* ifb */
		if (psat->rejp[2] >= opt->maxrej[2])   psat->reset |= 0x42;	

		/* amb */		
		if (psat->rejc[0] >= opt->maxrej[3])   psat->reset |= 0x10;	
		if (psat->rejc[1] >= opt->maxrej[3])   psat->reset |= 0x20;
		if (psat->rejc[2] >= opt->maxrej[3])   psat->reset |= 0x40;
	}
}
/* pseudo & doppler error detect */
static void pseudo_detect(const obsd_t *obs, int n, const prcopt_t *opt, ssat_t *ssat) 
{ 
    ssat_t *psat;
	double thresp12 = opt->threscheck[3];
	
	for (int i = 0; i < n && i < MAXOBS; i++) {	

		psat = ssat + obs[i].sat - 1;

		psat->eobs = 0;

		psat->diff[9] = psat->diff[10] = psat->diff[11] = psat->diff[12] = 0.0;

		if (fabs(timediff(obs[i].time, psat->prev_T[0])) > opt->acctime) continue;

		if (!ISZERO(obs[i].P[0]) && !ISZERO(obs[i].P[1]) && !ISZERO(psat->prev_P[0][0]) && !ISZERO(psat->prev_P[0][1])) {			
			psat->diff[9] =  (obs[i].P[0] - obs[i].P[1]) - (psat->prev_P[0][0] - psat->prev_P[0][1]);   
		}
		if (!ISZERO(obs[i].P[1]) && !ISZERO(obs[i].P[2]) && !ISZERO(psat->prev_P[0][1]) && !ISZERO(psat->prev_P[0][2])) {			
			psat->diff[10] = (obs[i].P[1] - obs[i].P[2]) - (psat->prev_P[0][1] - psat->prev_P[0][2]);
		}
		if (!ISZERO(obs[i].P[0]) && !ISZERO(obs[i].P[2]) && !ISZERO(psat->prev_P[0][0]) && !ISZERO(psat->prev_P[0][2])) {		
			psat->diff[11] = (obs[i].P[0] - obs[i].P[2]) - (psat->prev_P[0][0] - psat->prev_P[0][2]);
		}
		if (!ISZERO(obs[i].P[0]) && !ISZERO(obs[i].D[0]) && !ISZERO(psat->prev_P[0][0]) && !ISZERO(psat->prev_D[0][0])) {
			double dt = timediff(obs[i].time, psat->prev_T[0]);
			double lambda;
			switch (psat->sys) {
				case SYS_GPS: lambda = CLIGHT/FREQ1; break;
				case SYS_GAL: lambda = CLIGHT/FREQ1; break;
				case SYS_CMP: lambda = CLIGHT/FREQ1_CMP; break;
			}
			if (fabs(obs[i].D[0]-psat->prev_D[0][0]) > 5*lambda || dt > 3) continue;
			psat->diff[12] = (obs[i].P[0] - psat->prev_P[0][0]) + (obs[i].D[0]+psat->prev_D[0][0])*lambda*dt/2;
		}
		if (fabs(psat->diff[9]) > thresp12 || fabs(psat->diff[10]) > thresp12 || fabs(psat->diff[11]) > thresp12) {
			psat->eobs = 0x7;
			satno2id(obs[i].sat, id);
			trace(2, "%s %15.4f %15.4f %15.4f\n", id, psat->diff[9], psat->diff[10], psat->diff[11], psat->eobs);
		}
	}	
}
/* clock detect and repair */
static void clkJump_detect(const obsd_t *obs, int n, const nav_t *nav, const prcopt_t *opt, ssat_t *ssat)
{
	int i, sat;
	double k1, k2, freq, dp, dh, Ns, Na, S, M, clkJump;

	k1 = 1E-6*CLIGHT - 3 * 5;
	k2 = 1E-4;
	
	// clock detect
	Ns = Na = S = M = clkJump = 0;
	for (i = 0; i < n && i < MAXOBS; i++) {			
		sat = obs[i].sat;
		freq = sat2freq(sat, obs[i].code[0], nav);

		if (!ssat[sat - 1].vs) continue;
		if ( ssat[sat - 1].eobs&1) continue;
		if (ISZERO(obs[i].P[0]) || ISZERO(obs[i].L[0])) continue;
		if (fabs(fabs(timediff(obs[i].time, ssat[sat - 1].prev_T[0])) - opt->ti) > DTTOL) continue;	

		Na++;
		dp =  obs[i].P[0] - ssat[sat - 1].prev_P[0][0];
		dh = (obs[i].L[0] - ssat[sat - 1].prev_L[0][0]) * (CLIGHT / freq);

		if (fabs(dp - dh) < k1) continue;		
		S += dp - dh;
		Ns++;
	}

	// clock repair
	if (Ns > Na / 2) {
		M = S / Ns;
		if (fabs(M) > 1E-3*CLIGHT) {
			M = 1E3*M / CLIGHT;
			if (fabs(M - ROUND(M)) < k2) {
				clkJump = ROUND(M)*(1E-3*CLIGHT);
			}
		}
		else if (fabs(M) > 1E-6*CLIGHT && fabs(M) < 1E-3*CLIGHT) {
			clkJump = M;
		}
		trace(1, "detect clock jump: detected time=%s clkJump=%8.6f\n", 
				  time_str(obs[0].time, 0), clkJump);
	}

	for (i = 0; i < n && i < MAXOBS; i++) {			
		sat = obs[i].sat;
		if (ssat[sat - 1].vs) {		
			ssat[sat - 1].php += clkJump;
		}
	}
}
/* LC observation ---------------------------------------------------------------------- */
static void linear_combination(const obsd_t *obs, int n, const nav_t *nav, const prcopt_t *opt, ssat_t *ssat) 
{
	ssat_t *psat;
	int i, k;
	double L[NFREQ], P[NFREQ], freq[NFREQ], alfa, lambda;

	for (i = 0; i < n && i < MAXOBS; i++) {
		psat = ssat + obs[i].sat - 1;

		for (k = 0; k < 9; k++) psat->LC[k] = 0.0;
		for (k = 0; k < NFREQ; k++) P[k] = L[k] = freq[k] = 0;
		for (k = 0; k < opt->nf; k++) {								
			P[k] = obs[i].P[k];
			L[k] = obs[i].L[k];
			freq[k] = sat2freq(obs[i].sat, obs[i].code[k], nav);
			if (psat->eobs&(0x1<<k)) P[k] = 0.0;
		}
		
		// GF:L1-L2
		if (!ISZERO(freq[0]) && !ISZERO(freq[1]) && !ISZERO(L[0]) && !ISZERO(L[1])) {
			psat->LC[0] = CLIGHT * (L[0] / freq[0] - L[1] / freq[1]);
		} 
		// GF:L2-L3
		if (!ISZERO(freq[1]) && !ISZERO(freq[2]) && !ISZERO(L[1]) && !ISZERO(L[2])) {
			psat->LC[1] = CLIGHT * (L[1] / freq[1] - L[2] / freq[2]);
		}
		// GF:L1-L3
		if (!ISZERO(freq[0]) && !ISZERO(freq[2]) && !ISZERO(L[0]) && !ISZERO(L[2])) {
			psat->LC[2] = CLIGHT * (L[0] / freq[0] - L[2] / freq[2]);
		}
		// MW:L1L2-P1P2
		if (!ISZERO(freq[0]) && !ISZERO(freq[1]) && !ISZERO(L[0]) && !ISZERO(L[1]) && !ISZERO(P[0]) && !ISZERO(P[1])) {
			lambda = CLIGHT / (freq[0] - freq[1]);
			psat->LC[3] = (L[0] - L[1]) - (freq[0] * P[0] + freq[1] * P[1]) / (freq[0] + freq[1]) / lambda;			
		}
		// MW:L2L3-P1P2
		if (!ISZERO(freq[1]) && !ISZERO(freq[2]) && !ISZERO(L[1]) && !ISZERO(L[2]) && !ISZERO(P[0]) && !ISZERO(P[1])) {
			alfa   = (SQR(freq[0])*SQR(freq[1]) - SQR(freq[1])*freq[1] * freq[2]) / 
				     (freq[1] * freq[2] * (SQR(freq[0]) - SQR(freq[1])));
			lambda = CLIGHT / (freq[1] - freq[2]);
			psat->LC[4] = (L[1] - L[2]) - ((1 - alfa) * P[0] + alfa* P[1]) / lambda;					
		}
		// MW:L1L3-P1P2 
		if (!ISZERO(freq[0]) && !ISZERO(freq[2]) && !ISZERO(L[0]) && !ISZERO(L[2]) && !ISZERO(P[0]) && !ISZERO(P[1])) {
			alfa   = (SQR(freq[0])*SQR(freq[1]) - SQR(freq[1])*freq[0] * freq[2]) / 
				     (freq[0] * freq[2] * (SQR(freq[0]) - SQR(freq[1])));
			lambda = CLIGHT / (freq[0] - freq[2]);
			psat->LC[5] = (L[0] - L[2]) - ((1 - alfa) * P[0] + alfa* P[1]) / lambda;					
		}

		// MP1 = Pi-(fi^2+fj^2)/(fi^2-fj^2)*Li+(2*fj^2)/(fi^2-fj^2)*Lj    L:cycle
		if (!ISZERO(freq[0]) && !ISZERO(freq[1]) && !ISZERO(L[0]) && !ISZERO(L[1]) && !ISZERO(P[0])) {
			alfa   = (SQR(freq[0]) + SQR(freq[1])) / (SQR(freq[0]) - SQR(freq[1])) * CLIGHT / freq[0];
			lambda = 2 * SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1])) * CLIGHT / freq[1];
			psat->LC[6] = P[0] - alfa * L[0] + lambda * L[1];
		}
		// MP2
		if (!ISZERO(P[1]) && !ISZERO(freq[0]) && !ISZERO(freq[1]) && !ISZERO(L[0]) && !ISZERO(L[1])) {
			alfa   = (SQR(freq[1]) + SQR(freq[0])) / (SQR(freq[1]) - SQR(freq[0])) * CLIGHT / freq[1];
			lambda = 2 * SQR(freq[0]) / (SQR(freq[1]) - SQR(freq[0])) * CLIGHT / freq[0];
			psat->LC[7] = P[1] - alfa * L[1] + lambda * L[0];
		}
		// MP3 (改进)
		if (!ISZERO(P[2]) && !ISZERO(freq[0]) && !ISZERO(freq[2]) && !ISZERO(L[0]) && !ISZERO(L[2])) {	
			alfa   = (SQR(freq[2]) + SQR(freq[0])) / (SQR(freq[2]) - SQR(freq[0])) * CLIGHT / freq[2];
			lambda = 2 * SQR(freq[0]) / (SQR(freq[2]) - SQR(freq[0])) * CLIGHT / freq[0];
			psat->LC[8] = P[2] - alfa * L[2] + lambda * L[0];
		}
		// DIF
		if (!ISZERO(freq[0]) && !ISZERO(freq[1]) && !ISZERO(freq[2]) && !ISZERO(L[0]) && !ISZERO(L[1]) && !ISZERO(L[2])) {
			alfa 	= CLIGHT / (SQR(freq[0]) - SQR(freq[1]));
			lambda  = CLIGHT / (SQR(freq[0]) - SQR(freq[2]));
			psat->LC[9] = (alfa * (freq[0]*L[0] - freq[1]*L[1])) - (lambda * (freq[0]*L[0] - freq[2]*L[2]));
		}
	}
}
/* detect cycle slip by LLI ------------------------------------------------ */
static uint8_t detslp_ll(const obsd_t *obs, const prcopt_t *opt, ssat_t *psat, int rcv)
{
    uint32_t LLI, slip = 0;

    for (int f=0;f<NF(opt);f++) {
        
        if (ISZERO(obs->L[f]) || fabs(timediff(obs->time, psat->prev_T[rcv-1])) < DTTOL) {	
            continue;
        }
        /* restore previous LLI */
        if (rcv==1) LLI=getbitu(&psat->slip[f],0,2); /* rover */
        else        LLI=getbitu(&psat->slip[f],2,2); /* base  */
        
        /* detect slip by cycle slip flag in LLI */
		if(obs->LLI[f] & 1) slip |= 0x1<<f;
        
        /* detect slip by parity unknown flag transition in LLI */
        if (((LLI&2)&&!(obs->LLI[f]&2))||(!(LLI&2)&&(obs->LLI[f]&2))) {		
            slip |= 0x1<<f;
        }
        /* save current LLI */
        if (rcv==1) setbitu(&psat->slip[f],0,2,obs->LLI[f]);		
        else        setbitu(&psat->slip[f],2,2,obs->LLI[f]);
    }
	return slip;
}
/* cycle Jump detect ------------------------------------------------------- */
static void cycJump_detect(const obsd_t *obs, int n, const nav_t *nav, const prcopt_t *opt, ssat_t *ssat)		
{
	int i, k, slip, lli, ltt, lcc, vflg;
	double el, thresgf[2], thresmw[2];
	ssat_t *psat;

	for (i = 0; i < n && i < MAXOBS; i++) {

		psat = ssat + obs[i].sat - 1;

		slip = lli = ltt = lcc = 0;

		// LLI
		lli = detslp_ll(obs + i, opt, psat, 1);	

		ltt = (fabs(timediff(obs[i].time, psat->prev_T[0])>opt->acctime)) ? 7 : 0;

		// TurboEdit	
		for (k = 0, vflg = 0; k < 6; k++) {
			psat->diff[k] = 0.0;
			if (!ISZERO(psat->LC[k]) && !ISZERO(psat->prev_LC[k])) {
				psat->diff[k] = psat->LC[k] - psat->prev_LC[k];
				vflg |= (0x1)<<k;
			}
		}
		el = psat->azel[1]*R2D;
		thresgf[0] = (el >= 15) ? opt->thresslipgf[0] : (7-0.4*el)*opt->thresslipgf[0];
		thresgf[1] = (el >= 15) ? opt->thresslipgf[1] : (7-0.4*el)*opt->thresslipgf[1];
		thresmw[0] = (el >= 15) ? opt->thresslipmw[0] : (4-0.2*el)*opt->thresslipmw[0];
		thresmw[1] = (el >= 15) ? opt->thresslipmw[1] : (4-0.2*el)*opt->thresslipmw[1];
		
		if (vflg&0011 && (fabs(psat->diff[0])>thresgf[0] || fabs(psat->diff[3])>thresmw[0])) lcc |= vflg&0066 ? 03 : 07;	
		if (vflg&0022 && (fabs(psat->diff[1])>thresgf[1] || fabs(psat->diff[4])>thresmw[1])) lcc |= 04;     
		if (vflg&0044 && (fabs(psat->diff[2])>thresgf[0] || fabs(psat->diff[5])>thresmw[0])) lcc |= 04;   

		// if ((vflg&0011)==9  && (fabs(psat->diff[0])<thresgf[0] && fabs(psat->diff[3])<thresmw[0])) lcc &= 04;
		// if ((vflg&0022)==18 && (fabs(psat->diff[1])<thresgf[1] && fabs(psat->diff[4])<thresmw[1])) lcc &= 01;
		// if ((vflg&0044)==36 && (fabs(psat->diff[2])<thresgf[0] && fabs(psat->diff[5])<thresmw[0])) lcc &= 02;
		
		// slip
		slip = lli | ltt | lcc;
		for (k = 0; k < opt->nf; k++) {
			psat->slip[k] |= (slip >> k)&1;
			psat->half[k]  = !(obs[i].LLI[k]&2);
		}
		
		psat->reset |= slip<<4;		/* 2022-08-17 add by zh */
	}
}
/* update satellite information ------------------------------------------------ */
static void update_ssat(const obsd_t *obs, int n, const prcopt_t *opt, double *rs, double *dts, 
	const int *svh, int *stat, exc_t *exc, uint8_t *ns, ssat_t *ssat) 
{
	int i, k, type, freq, sat1, sat2, nsat[4] = {0}, isok = 0, nn = 0;
	double zeros[6] = { 0 };
	ssat_t *psat;

	for (i = 0; i < exc->n; i++) {
		readvflg(exc->flg[i], &type, &freq, &sat1, &sat2);
		switch(type) {
		case 1: ssat[sat1-1].rejc[freq-1] += (ssat[sat1-1].rejc[freq-1]>100) ? 0 : 1; break;
		case 2: ssat[sat1-1].rejp[freq-1] += (ssat[sat1-1].rejp[freq-1]>100) ? 0 : 1; break;
		}
	}

	for (i = (*ns) = 0; i < n && i < MAXOBS; i++) {			
		psat = ssat + obs[i].sat - 1;

		if (psat->ini&0x01) {	// ION
			if (psat->vsat[0]==3 && psat->vsat[1]==3) setbitu(&psat->reset,7,1,0);
			// if (psat->vsat[0]&0x1 || psat->vsat[1]&0x1) setbitu(&psat->reset,7,1,0);
		}
		if (psat->ini&0x02) {	// IFB
			if (psat->vsat[0]==3 && psat->vsat[2]==3) setbitu(&psat->reset,6,1,0);  
			// if (psat->vsat[0]&0x1 && psat->vsat[2]&0x1) setbitu(&psat->reset,6,1,0);		
		}
		for (k = 0; k < opt->nf; k++) if (psat->ini&(0x10<<k)) {	//AMB
			if (psat->vsat[k]==3) setbitu(&psat->reset,3-k,1,0);
			// if (psat->vsat[k]&0x2) setbitu(&psat->reset,3-k,1,0);
		}

		for (k = 0; k < opt->nf; k++) {
			if (psat->vsat[k]&0x1) {
				psat->rejp[k] = 0;
			}
			if (psat->vsat[k]&0x2) {
				psat->rejc[k]  = psat->outc[k] = 0; 
				psat->lock[k] += (psat->lock[k]>MAX_UINT8) ? 0 : 1;
				psat->fix[k]   = 1;
			}
		}
		
		if (ISZERO(obs[i].P[1]) || ISZERO(obs[i].L[1])) continue;	// add by zh 2022-09-10 
		if (!(psat->vsat[0]&1) || (!(psat->vsat[1]&1))) continue;	

		switch(psat->sys) {
			case SYS_GPS: nsat[0]++; break;
			case SYS_GAL: nsat[1]++; break;
			case SYS_CMP: nsat[2]++; break;
		}
		
		psat->prev_T[0] = obs[i].time;
		for (k = 0; k < opt->nf; k++) {
			psat->prev_P[0][k] = obs[i].P[k];		
			psat->prev_L[0][k] = obs[i].L[k];
			psat->prev_D[0][k] = obs[i].D[k];
		}
		for (k = 0; k < 10; k++) psat->prev_LC[k] = psat->LC[k];
		if (svh[i] < 0) {
			matcpy(psat->prev_rs,   zeros, 6, 1);
			matcpy(psat->prev_dts,  zeros, 2, 1);
		}
		else {
			matcpy(psat->prev_rs,   rs + i*6, 6, 1);
			matcpy(psat->prev_dts, dts + i*2, 2, 1);
		}

		if ((psat->vsat[0]&2) && (psat->vsat[1]&2)) (*ns)++;	
	}
	// double tow=time2gpst(obs->time,NULL);
	// trace(1,"%.1f %5d %5d %5d\n", tow, nsat[0], nsat[1], nsat[2]);

	/* update solution status */
	for (i = 0; i < 4; i++) if (nsat[i] > 0) nn += 1;
	*stat = (*ns) < 5 + nn ? SOLQ_NONE : *stat;

	// *stat = (*ns) < (MIN_NSAT_SOL) ? SOLQ_NONE : *stat;

	// for (i = 0; i < 4; i++) if (nsat[i] >= 4) isok = 1;
	// *stat = !isok ? SOLQ_NONE : *stat;
}
/* temporal update of position --------------------------------------------- */
static void udpos_ppp(rtk_t *rtk, const double *xg, const double *Pg)			
{
	double *F, *P, *FP, *x, *xp, pos[3], Q[9] = { 0 }, Qv[9], var = 0.0, std = 1E+010;			
	int i, j, *ix, nx;

	trace(3, "udpos_ppp:\n");

	/* reinitialize */
	if (rtk->reset&1) {
		for (i = 0; i<3; i++) initx(rtk, 0.0, 0.0, i);
		setbitu(&rtk->reset, 7, 1, 0);								
	}

	/* fixed mode */
	if (rtk->opt.mode == PMODE_PPP_FIXED) {
		for (i = 0; i<3; i++) initx(rtk, rtk->opt.ru[i], 1E-8, i);
		return;
	}
	/* initialize position for first epoch */									
	if (norm(rtk->x, 3) <= 0.0) {				
		for (i = 0; i<3; i++) initx(rtk, rtk->sol.rr[i], SQR(rtk->opt.P0[0]), i);	//SPP位置精度														
		for (i = 3; i<6; i++) initx(rtk, rtk->sol.rr[i], SQR(rtk->opt.P0[6]), i);
		for (i = 6; i<9; i++) initx(rtk, 1E-6, SQR(rtk->opt.P0[7]), i);
		return;
	}
	/* static ppp mode */
	if (rtk->opt.mode == PMODE_PPP_STATIC) {
		for (i = 0; i<3; i++) {
			rtk->P[i*(1 + rtk->nx)] += SQR(rtk->opt.Qt[0])*fabs(rtk->tt);				
		}
		return;
	}

	/* kinmatic mode */
	/* generate valid state index */
	ix = imat(rtk->nx, 1);
	for (i = nx = 0; i<rtk->nx; i++) {
		if (rtk->x[i] != 0.0&&rtk->P[i + i*rtk->nx]>0.0) ix[nx++] = i;
	}
	if (nx<9) {
		free(ix);
		return;
	}

	/* state transition of position/velocity/acceleration */
	F = eye(nx); P = mat(nx, nx); FP = mat(nx, nx); x = mat(nx, 1); xp = mat(nx, 1);

	for (i = 0; i<6; i++) {
		F[i + (i + 3)*nx] = rtk->tt;				
	}
	for (i = 0; i<3; i++) {
		F[i + (i + 6)*nx] = SQR(rtk->tt) / 2.0;
	}
	for (i = 0; i<nx; i++) {
		x[i] = rtk->x[ix[i]];								
		for (j = 0; j<nx; j++) {
			P[i + j*nx] = rtk->P[ix[i] + ix[j] * rtk->nx];
		}
	}
	/* x=F*x, P=F*P*F+Q */
	matmul("NN", nx,  1, nx, 1.0, F,  x, 0.0, xp);					
	matmul("NN", nx, nx, nx, 1.0, F,  P, 0.0, FP);
	matmul("NT", nx, nx, nx, 1.0, FP, F, 0.0, P);			

	for (i = 0; i<nx; i++) {
		rtk->x[ix[i]] = xp[i];
		for (j = 0; j<nx; j++) {
			rtk->P[ix[i] + ix[j] * rtk->nx] = P[i + j*nx];
		}
	}
	/* process noise added to only acceleration */
	Q[0] = Q[4] = SQR(rtk->opt.Qt[6])*fabs(rtk->tt);	
	Q[8] = SQR(rtk->opt.Qt[7])*fabs(rtk->tt);			
	ecef2pos(rtk->x, pos);
	covecef(pos, Q, Qv);
	for (i = 0; i<3; i++) for (j = 0; j<3; j++) {
		rtk->P[i + 6 + (j + 6)*rtk->nx] += Qv[i + j * 3];   
	}
	free(ix); free(F); free(P); free(FP); free(x); free(xp);

	/* check variance of estimated position */
    for (i=0,var=0.0;i<3;i++) var+=rtk->P[i+i*rtk->nx]/3.0;

    if (var>SQR(30.0)) {
		/* reset position with large variance */
		if (Pg) for (i=0,std=0.0;i<3;i++) std += sqrt(Pg[i+i*3])/3.0;
		if (xg && std < 10.0) {
			for (i=0;i<3;i++) initx(rtk, xg[i], SQR(rtk->opt.P0[0]), i);
			for (i=3;i<6;i++) initx(rtk, rtk->sol.rr[i], SQR(rtk->opt.P0[6]), i);
        	for (i=6;i<9;i++) initx(rtk, 1E-6, SQR(rtk->opt.P0[7]), i);
		} else {
			for (i=0;i<3;i++) initx(rtk, rtk->sol.rr[i], SQR(rtk->opt.P0[0]), i);
        	for (i=3;i<6;i++) initx(rtk, rtk->sol.rr[i], SQR(rtk->opt.P0[6]), i);
        	for (i=6;i<9;i++) initx(rtk, 1E-6, SQR(rtk->opt.P0[7]), i);
		}
        trace(2,"reset rtk position due to large variance: var=%.3f\n",var);
        return;
    }
}
/* temporal update of clock --------------------------------------------------*/
static void udclk_ppp(rtk_t *rtk)
{
	double dtr; 

	trace(3, "udclk_ppp:\n");

	/* initialize every epoch for clock (white noise) */
	for (int i = 0; i < 4; i++) {
		dtr = rtk->sol.dtr[i];
		initx(rtk, CLIGHT*dtr, SQR(rtk->opt.P0[1]), IC(i, &rtk->opt));
	}
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(rtk_t *rtk)
{
	const double azel[] = { 0.0, PI / 2.0 };
	double pos[3], zwd;
	int i = IT(&rtk->opt), j;

	trace(3, "udtrop_ppp:\n");

	if (rtk->reset&2) {
		initx(rtk, 0.0, 0.0, i);
		setbitu(&rtk->reset, 6, 1, 0); 
	}

	if (!ISZERO(rtk->x[i])) {
		rtk->P[i + i*rtk->nx] += SQR(rtk->opt.Qt[2])*fabs(rtk->tt);			

		if (rtk->opt.tropopt >= TROPOPT_ESTG) {
			for (j = i + 1; j<i + 3; j++) {
				rtk->P[j + j*rtk->nx] += SQR(rtk->opt.Qt[2] * 0.1)*fabs(rtk->tt);
			}
		}	
	}
	else {											
		ecef2pos(rtk->sol.rr, pos);
		zwd = tropmodel(rtk->sol.time, pos, azel, 1, REL_HUMI, 2);				
		initx(rtk, zwd, SQR(rtk->opt.P0[2]), i);

		if (rtk->opt.tropopt >= TROPOPT_ESTG) {
			for (j = i + 1; j<i + 3; j++) initx(rtk, 1E-6, VAR_GRA, j);
		}
	}
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void udiono_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)					
{
	int i, k, ii, sat;
	double freq[NFREQ], ion, sinel, pos[3], dantr[3] = { 0 }, dants[3] = { 0 };
	double L[NFREQ], P[NFREQ], Lc, Pc;
	ssat_t *ssat;

	trace(3, "udiono_ppp:\n");

	for (i = 0; i < MAXSAT; i++) {
		ii = II(i + 1, &rtk->opt);
		if (!ISZERO(rtk->x[ii]) && rtk->ssat[i].reset&0x01) initx(rtk, 0.0, 0.0, ii);	
		if (!ISZERO(rtk->x[ii])) {
			sinel = sin(MAX(rtk->ssat[i].azel[1], 15.0*D2R));
			rtk->P[ii + ii*rtk->nx] += SQR(rtk->opt.Qt[3] / sinel)*fabs(rtk->tt);		
		}	
	}
	ecef2pos(rtk->sol.rr, pos);
	for (i = 0; i<n; i++) {
		sat = obs[i].sat;
		ssat = rtk->ssat + sat - 1;
		if (!ssat->vs) continue;
		ii = II(sat, &rtk->opt);

		for (k = 0; k < NFREQ; k++) freq[k] = sat2freq(sat, obs[i].code[k], nav);	
		
		if (ISZERO(rtk->x[ii])) {

			corr_meas(obs + i, nav, ssat->azel, &rtk->opt, dantr, dants, 
					  ssat->phw, ssat->php, L, P, &Lc, &Pc);	
			ion = 0.0;
			if (!ISZERO(P[0]) && !ISZERO(P[1]) && !(ssat->eobs&3) && !ISZERO(freq[0]) && !ISZERO(freq[1])) {		
				ion = (P[0] - P[1]) / (1.0 - SQR(freq[0] / freq[1]));				
			}

			if (!ISZERO(ion)) {
				rtk->ssat[sat - 1].ini |= 0x01;			
				initx(rtk, ion, SQR(rtk->opt.P0[3]), ii);
				satno2id(sat, id);
				trace(2,"$ION %s %s %10.4f\n", id, time_str(rtk->sol.time, 0), ion);
			}
		}
	}
}
/* temporal update of inter frequency bias ---------------------------------- */
static void udifb_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	int i, k, ii, sat;
	double ion, ifb, gamma, freq[3];
	for (i = 0; i < MAXSAT; i++) {
		ii = IS(i + 1, &rtk->opt);

		if (!ISZERO(rtk->x[ii]) && rtk->ssat[i].reset&0x02) initx(rtk, 0.0, 0.0, ii);
		if (!ISZERO(rtk->x[ii])) {
			rtk->P[ii + ii*rtk->nx] += SQR(rtk->opt.Qt[4])*fabs(rtk->tt);		
		}
	}
	for (i = 0; i < n; i++) {
		sat = obs[i].sat;
		if (!rtk->ssat[sat - 1].vs) continue;

		ii = IS(sat, &rtk->opt);
		for(k = 0; k < NFREQ; k++) freq[k] = sat2freq(sat, obs[i].code[k], nav);
		
		if (ISZERO(rtk->x[ii])) {

			ifb = 0.0;
			ion = rtk->x[II(sat, &rtk->opt)];

			if (rtk->ssat[sat - 1].eobs&5) continue;
			if (!ISZERO(ion) && !ISZERO(freq[0]) && !ISZERO(freq[2]) && !ISZERO(obs[i].P[0]) && !ISZERO(obs[i].P[2])) {
				gamma = SQR(freq[0]/freq[2]) - 1;
				ifb = (obs[i].P[2] - obs[i].P[0] - gamma*ion);  		
			}
			if (!ISZERO(ifb)) {
				rtk->ssat[sat - 1].ini |= 0x02;
				initx(rtk, ifb, SQR(rtk->opt.P0[4]), ii);
				satno2id(sat, id);
				trace(2,"$IFB %s %s %10.4f\n", id, time_str(rtk->sol.time, 0), ifb);
			}
		}	
	}
}
/* temporal update of phase biases ------------------------------------------- */
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)											
{
	int i, ii, f, sat;
	double L[NFREQ], P[NFREQ], Lc, Pc, freq[NFREQ], bias, gamma;
	double ion, dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };

	trace(3, "udbias  : n=%d\n", n);

	for (f = 0; f<NF(&rtk->opt); f++) {	

		/* reset phase-bias if expire obs outage counter */
		for (i = 0; i<MAXSAT; i++) {
			ii = IB(i + 1, f, &rtk->opt);
			if (!ISZERO(rtk->x[ii]) && (rtk->ssat[i].reset&(0x10<<f))) {  
				rtk->ssat[i].lock[f] = 0;
				initx(rtk, 0.0, 0.0, ii);
			}
			if (!ISZERO(rtk->x[ii])) {
				rtk->P[ii + ii*rtk->nx] += SQR(rtk->opt.Qt[5])*fabs(rtk->tt);		
			}
		}
		for (i = 0; i<n&&i<MAXOBS; i++) {
			
			sat = obs[i].sat;
			ii = IB(sat, f, &rtk->opt);
			if (!rtk->ssat[sat-1].vs || !ISZERO(rtk->x[ii])) continue;

			corr_meas(obs + i, nav, rtk->ssat[sat - 1].azel, &rtk->opt, dantr, dants, 
					  rtk->ssat[sat - 1].phw, rtk->ssat[sat - 1].php, L, P, &Lc, &Pc);

			bias = 0.0;
			if (rtk->opt.ionoopt == IONOOPT_IFLC) {
				if (ISZERO(Lc) || ISZERO(Pc)) continue;
				bias = Lc - Pc;
			}
			else if (rtk->opt.ionoopt == IONOOPT_EST) {
				freq[0] = sat2freq(sat, obs[i].code[0], nav);
				freq[f] = sat2freq(sat, obs[i].code[f], nav);
				ion = rtk->x[II(sat, &rtk->opt)];

				if (rtk->ssat[sat - 1].eobs&0x1) continue;
				if (ISZERO(L[f]) || ISZERO(P[0]) || ISZERO(freq[f]) || ISZERO(ion)) continue;

				gamma = SQR(freq[0] / freq[f]) + 1;	  	
				bias = L[f] - P[0] + ion*gamma;
			}
			if (!ISZERO(bias)) {
				rtk->ssat[sat - 1].ini |= 0x10 << f; 				// ini
				initx(rtk, bias, SQR(rtk->opt.P0[5]), ii); 
				satno2id(sat, id);
				trace(2,"$AM%d %10.0f %s %s %10.4f %10.4f %10.4f\n", f+1, time2gpst(rtk->sol.time, NULL), id, time_str(rtk->sol.time, 0), bias, ion, 
					sqrt(rtk->P[II(sat, &rtk->opt)+II(sat, &rtk->opt)*rtk->nx]));
			}
		}
	}//end frequency
}
/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
	const double *xg, const double *Pg)
{
	trace(3, "udstate_ppp: n=%d\n", n);

	/* temporal update of position */
	udpos_ppp(rtk, xg, Pg);

	/* temporal update of clock */
	udclk_ppp(rtk);										

	/* temporal update of tropospheric parameters */
	if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {					
		udtrop_ppp(rtk);
	}
	/* temporal update of ionospheric parameters */
	if (rtk->opt.ionoopt == IONOOPT_EST) {			
		udiono_ppp(rtk, obs, n, nav);
	}
	/* temporal update of inter frequency bias */
	if (rtk->opt.nf > 2) {
		udifb_ppp(rtk, obs, n, nav);
	}
	/* temporal update of phase-bias */
	udbias_ppp(rtk, obs, n, nav);						
}
/* satellite antenna phase center variation ----------------------------------*/
static void satantpcv(const double *rs, const double *rr, const pcv_t *pcv, double *dant)
{
	double ru[3], rz[3], eu[3], ez[3], nadir, cosa;
	int i;

	for (i = 0; i<3; i++) {
		ru[i] = rr[i] - rs[i];
		rz[i] = -rs[i];
	}
	if (!normv3(ru, eu) || !normv3(rz, ez)) return;

	cosa = dot(eu, ez, 3);
	cosa = cosa<-1.0 ? -1.0 : (cosa>1.0 ? 1.0 : cosa);
	nadir = acos(cosa);

	antmodel_s(pcv, nadir, dant);							
}
/* exclude observation -------------------------------------------------------*/
static int exclude(int itype, int ifreq, int isat1, int isat2, const exc_t *exc) 
{
	int i, type, freq, sat1, sat2, stat = 0;
	for (i = 0; i < exc->n; i++) {
		readvflg(exc->flg[i], &type, &freq, &sat1, &sat2);
		if (type == itype) {
			switch (type){
			case 1: stat = (freq == ifreq && sat1 == isat1) ? 1 : 0; break;
			case 2: stat = (freq == ifreq && sat1 == isat1) ? 1 : 0; break;
			case 4: stat = 1; break;
			case 8: stat = (sat1 + sat2 == isat1 + isat2) ? 1 : 0; break;
			}
		}
		if (stat) return stat;
	}
	return stat;
}
/* phase and code residuals --------------------------------------------------*/
static int ppp_res(const obsd_t *obs, int n, const double *rs, const double *dts, 
	const double *var_rs, const int *svh, const nav_t *nav, const double *x, 
	rtk_t *rtk, double *v, double *H, double *var, int *vflg, int *iref, const exc_t *exc, int cc)
{
	int i, j, k, sat, sys, nv = 0, nx = rtk->nx;
	char str[32], id[32];
	prcopt_t *opt = &rtk->opt;
	double freq[NFREQ] = { 0 }, y, r, cdtr, bias, C, rr[3], pos[3], e[3], dtdx[3], L[NFREQ], P[NFREQ], Lc, Pc;
	double dtrp = 0.0, dion = 0.0, vart = 0.0, vari = 0.0, ifb;
	double dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 }, mAzel[NSYS] = { 0 }, *azel, *Hi;

	time2str(obs[0].time, str, 0);

	if (iref) for(i = 0; i<NSYS; i++) iref[i] = -1;
	for (i = 0; i<3; i++) rr[i] = x[i] + rtk->dr[i];
	ecef2pos(rr, pos);
	for (i = 0; i<n&&i<MAXOBS; i++) {
		sat = obs[i].sat;
		satno2id(sat, id);
		sys  = rtk->ssat[sat - 1].sys;
		azel = rtk->ssat[sat - 1].azel;

		for (j = 0; j < NFREQ; j++) freq[j] = sat2freq(sat, obs[i].code[j], nav);			

		if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 || satazel(pos, e, azel)<opt->elmin) {			
			continue;
		}
		
		if (!(sys) || !rtk->ssat[sat - 1].vs || satexclude(obs[i].sat, svh[i], opt)) {							
			continue;
		}
		
		/* tropospheric and ionospheric model */
		if (!tropcorr(obs[i].time, pos, azel, opt, x, dtdx, nav, &dtrp, &vart) ||
			!ionocorr(obs[i].time, pos, azel, opt, sat, x, nav, &dion, &vari)) {
			continue;
		}

		/* satellite and receiver antenna model */
		if (opt->posopt[0]) satantpcv(rs + i * 6, rr, nav->pcvs + sat - 1, dants);
		antmodel(opt->pcvr, opt->antdel[0], azel, opt->posopt[1], dantr);
		
		/* phase windup model */
		if (!model_phw(rtk->sol.time, sat, nav->pcvs[sat - 1].type, opt->posopt[2] ? 2 : 0, 
					   rs + i * 6, rr, &rtk->ssat[sat - 1].phw)) {	
			continue;
		}

		/* corrected phase and code measurements */
		corr_meas(obs + i, nav, azel, &rtk->opt, dantr, dants, rtk->ssat[sat - 1].phw, 
				  rtk->ssat[sat - 1].php, L, P, &Lc, &Pc);

		/* stack phase and code residuals {L1,P1,L2,P2,...} */
		for (j = 0; j < 2 * NF(opt); j++) {										

			bias = ifb = 0.0;	

			if (opt->ionoopt == IONOOPT_IFLC) {						
				if ((y = j % 2 == 0 ? Lc : Pc) == 0.0) continue;			
			}
			else {
				if ((y = j % 2 == 0 ? L[j / 2] : P[j / 2]) == 0.0) continue;
				if (freq[j / 2] == 0.0 || freq[0] == 0.0) continue;
			}
			
			/* exclude obs */
			if (exclude((j % 2) + 1, (j / 2) + 1, sat, NULL, exc)) continue;
			if (j % 2 == 1 && (rtk->ssat[sat - 1].eobs&(0x1<<(j/2)))) continue; 

			/* coordinate */
			if (H) {
				Hi = H + nv*nx;
				for (k = 0; k < nx; k++) Hi[k] = k < 3 ? -e[k] : 0.0;					
			}

			/* receiver clock */
			if ((k = sys2num(sys)) < 0) continue;			
			cdtr = x[IC(k, opt)];
			if (H) {
				Hi[IC(k, opt)] = 1.0;
			}

			/* Trop */
			if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
				if (H) {
					for (k = 0; k<(opt->tropopt >= TROPOPT_ESTG ? 3 : 1); k++) {
						Hi[IT(opt) + k] = dtdx[k];
					}
				}
			}

			/* Iono */
			C = SQR(freq[0] / freq[j / 2])*(j % 2 == 0 ? -1.0 : 1.0);	
			if (opt->ionoopt == IONOOPT_EST) {
				if (rtk->x[II(sat, opt)] == 0.0) continue;
				if (H) {
					Hi[II(sat, opt)] = C;
				}
			}

			/* phase bias */
			if (j % 2 == 0) {
				if ((bias = x[IB(sat, j / 2, opt)]) == 0.0) continue;
				if (H) {
					Hi[IB(sat, j / 2, opt)] = 1.0;
				}
			}

			/* Inter frequency bias */
			if (j % 2 != 0 && j / 2 == 2) {
				if ((ifb = x[IS(sat, opt)]) == 0.0) continue;
				if (H) {
					Hi[IS(sat, opt)] = 1.0;
				}
			}

			/* residual */
			v[nv] = y - (r + cdtr - CLIGHT*dts[i * 2] + dtrp + C*dion + ifb + bias);

			/* variance */
			if (j % 2 == 0) var[nv] = varerr(obs[i].sat, azel[1], j / 2, j % 2, opt) + vart + SQR(C)*vari;
			if (j % 2 == 1)	var[nv] = varerr(obs[i].sat, azel[1], j / 2, j % 2, opt) + vart + SQR(C)*vari;
			if (sys == SYS_GLO && j % 2 == 1) var[nv] += SQR(opt->noise[1]);

			if (cc) {
				satno2id(sat, id);
				trace(1,"%s %s%d %s %10.4f ion=%10.4f bias=%12.4f %10.4f %12.5f\n", time_str(rtk->sol.time,0), (j % 2)==0?"L":"P", 
					j / 2 + 1, id, v[nv], dion, (j % 2)==0?bias:0.0, azel[1]*R2D, sqrt(var[nv]));
			}
			/* refence sate */
			if (iref) {
				k = sys2num(sys);
				if (mAzel[k] < azel[1]) {
					if (pppcorr_stat(obs[i].time, pos, sat, nav)) {		
						mAzel[k] = azel[1];		
						iref[k] = sat;
					}
				}
			}
			setvflg(&vflg[nv++], (j % 2) + 1, (j / 2 + 1), sat, 0);
			trace(2,"%s %s%d %10.6f %10.6f \n", id, (j % 2)==0?"L":"P", j / 2 + 1, sqrt(var[nv - 1]), v[nv - 1]);
		}	//end styple
	}	//end obs
	return nv;
}
/* constraint to local correction --------------------------------------------*/
static int cor_res(const obsd_t *obs, const int n, const int *refs, const nav_t *nav, const rtk_t *rtk, 
	const double *x, double *v, double *H, double *var, int *vflg, const exc_t *exc)
{
	int i, k, ii, jj, sat, sys, nx, nv = 0;
	double trop[3], std_trop[3], iono, std_iono;
	double *azel, rr[3], pos[3], *Hi;

	nx = rtk->nx;			
	for (i = 0; i < 3; i++) { rr[i] = x[i] + rtk->dr[i]; }	ecef2pos(rr, pos);

	/* constraint to external troposphere correction */
	if (rtk->opt.posopt[4]) {
		if (!exclude(4, NULL, NULL, NULL, exc)) {
			if (pppcorr_trop(obs[0].time, pos, trop, std_trop, nav)) {
				ii = IT(&rtk->opt);
				v[nv]   = trop[0] - x[ii];
				var[nv] = SQR(rtk->opt.err[5]);  	// SQR(0.01);		
				// var[nv] = SQR(0.05);
				if (H) {
					Hi = H + nv*nx;
					for (k = 0; k < nx; k++) Hi[k] = k == ii ? 1.0 : 0.0;
				}
#if 0
				trace(1,"%s TROP: %10.5f res=%10.5f\n", time_str(rtk->sol.time,0), trop[0], v[nv]);
#endif
				setvflg(&vflg[nv++], 4, 0, 0, 0);
			}
		}
	}
	/* constraint to external ionosphere correction */
	if (rtk->opt.posopt[5]) {			
		for (i = 0; i < n; i++) {
			sat  = obs[i].sat;
			sys  = rtk->ssat[sat - 1].sys;
			azel = rtk->ssat[sat - 1].azel;

			k = sys2num(sys);
			if (k < 0 || refs[k] < 0 || refs[k] == sat) continue;
			if (!rtk->ssat[sat - 1].vs || azel[1] < 15 * D2R) continue;
			if (exclude(8, NULL, sat, refs[k], exc)) continue;

			if (!pppcorr_stec(obs[i].time, pos, sat, refs[k], &iono, &std_iono, nav, 2)) continue;		
						
			ii = II(sat,     &rtk->opt);
			jj = II(refs[k], &rtk->opt);

			v[nv] = iono - (x[ii] - x[jj]);
#if 0
			char idref[32];
			satno2id(sat,id);
			satno2id(refs[k], idref);
			trace(1,"%s ION: %s %s %10.4f %15.6f %15.6f\n",time_str(obs[i].time,0), id, idref, v[nv], 
				iono, (x[ii] - x[jj]));
			// var[nv] = SQR(std_iono);
#endif			
			// var[nv] = SQR(0.02);  
			var[nv] = SQR(rtk->opt.err[6]);	// SQR(0.015); 
			if (H) {
				Hi = H + nv*nx;
				for (k = 0; k < nx; k++) Hi[k] = (k == ii) ? 1 : (k == jj) ? -1 : 0.0;
			}
			setvflg(&vflg[nv++], 8, 0, sat, refs[k]);
		}
	}
	return nv;
}
/* Add prior cooridition */
static int pri_res(rtk_t *rtk, const double *x, const double *Pp, const double *xg, 
	const double *Pg, double *v, double *H, double *var) 
{
	int nx = rtk->nx, nv = 0;
	double *Hi, var1 = 0.0, var2 = 0.0;

	if (xg && !ISZERO(xg[0]) && Pg) {
		for (int i = 0; i < 3; i++) {
			var1 += sqrt(Pg[i + i * 3]) / 3;
			var2 += sqrt(Pp[i + i * nx]) / 3;
		}
		if (var1 > var2) return nv;
		for (int i = 0; i < 3; i++) {
			v[nv] = xg[i] - x[i];

			var[nv] = sqrt(Pg[i + i * 3]) > 0.05 ? sqrt(Pg[i + i * 3]) : 0.05;

			if (H) {
				Hi = H + nv * nx;
				for (int j = 0; j < nx; j++) Hi[j] = j == i ? 1.0 : 0.0; 
			}
			// char id[32];
			// time2str(rtk->sol.time, id, 2);
			// trace(2, "%s gnss %12.4f %12.8f %12.8f %12.8f\n", id, v[nv], var[nv], 
			// 	sqrt(Pg[i + i*3]), sqrt(Pp[i+i*nx]));
			nv++;
		}
	}
	return nv;
}
/* reject obs by pre-fit residuals */
static int valpre(rtk_t *rtk, double *v, double *H, double *r, int *vflg, int nv, int nx, exc_t *exc) 
{
	int i, k, type, freq, sat1, sat2, nn, sys, num[NSYS] = { 0 }, *ind, *stat;
	double mid, *D;	
	char str[32], id[8];

	if (!rtk->opt.posopt[6] || rtk->reset&1) return nv;

	time2str(rtk->sol.time, str, 0);
	D = mat(nv, NSYS);			ind = imat(nv, NSYS);		stat = imat(nv, 1);

	for (i = 0; i < nv; i++) stat[i] = 1;
	for (i = 0; i < nv; i++) {
		readvflg(vflg[i], &type, &freq, &sat1, &sat2);
		if (type&3) { 
			sys = rtk->ssat[sat1 - 1].sys;
			switch (sys) {
				case SYS_GPS: D[0*nv+num[0]] = v[i]; ind[0*nv+num[0]++] = i; break;
				case SYS_GLO: D[1*nv+num[1]] = v[i]; ind[1*nv+num[1]++] = i; break;
				case SYS_GAL: D[2*nv+num[2]] = v[i]; ind[2*nv+num[2]++] = i; break;
				case SYS_CMP: D[3*nv+num[3]] = v[i]; ind[3*nv+num[3]++] = i; break;
			}
		}
	}
	for (i = 0; i < NSYS; i++) {
		mid = median(D+i*nv, num[i]);
		int cc=0;			
		for (k = 0; k < num[i]; k++) {
			readvflg(vflg[ind[i*nv+k]], &type, &freq, &sat1, &sat2);
			if (fabs(D[i*nv+k] - mid) > rtk->opt.threscheck[0]) {
				exc->flg[exc->n++] = vflg[ind[i*nv+k]];
				stat[ind[i*nv+k]] = 0;

				satno2id(sat1, id);
				trace(1,"pri-eli %s %s %s%d\tres=%10.4f\tlock=%4d\tazel=%10.4f\t %10.4f ion=%10.4f\n",str, id, type==1?"L":"P", 
						freq, v[ind[i*nv+k]], rtk->ssat[sat1-1].lock[freq-1], rtk->ssat[sat1-1].azel[1]*R2D, fabs(D[i*nv+k] - mid),
						rtk->x[II(sat1,&rtk->opt)]);
				cc=1;	
			}
		}
		if (cc) {
			for (k = 0; k < num[i]; k++) {
				readvflg(vflg[ind[i*nv+k]], &type, &freq, &sat1, &sat2);

				satno2id(sat1, id);
				trace(2,"%31s %s%d\tres=%10.4f\tlock=%4d\tazel=%10.4f\t %10.4f\n", id, type==1?"L":"P", freq, 
						v[ind[i*nv+k]], rtk->ssat[sat1-1].lock[freq-1], rtk->ssat[sat1-1].azel[1]*R2D, rtk->ssat[sat1-1].azel[0]*R2D);
			}
		}
	}
	for (i = nn = 0; i < nv; i++) {
		if (stat[i]) {
			v[nn] = v[i];		r[nn] = r[i];		vflg[nn] = vflg[i];	
			matcpy(H + nn*nx, H + i*nx, nx, 1);		nn++;	
		}
	}

	free(D);	 free(ind);	   free(stat);
	return nn;
}
/* reject obs by pos-fit residuals */
static int valpos(rtk_t *rtk, double *v, double *r, int *vflg, int na, double thres, exc_t *exc)	
{
	int i, k, type, sat1, sat2, freq, stat = 0;				
	char str[32], id[8];
	time2str(rtk->sol.time, str, 0);

	if (rtk->opt.posopt[7]) {
		for (i = 0; i < na; i++) {	
			readvflg(vflg[i], &type, &freq, &sat1, &sat2);
			if (type == 4) continue;	
			
			// pos_thres = (type == 1) ? 3.5 : (type == 2) ? 4.0 : 5;
			if (fabs(v[i]) > sqrt(r[i])*thres) {
				exc->flg[exc->n++] = vflg[i];
				stat++;

				satno2id(sat1, id);
				trace(2,"POS-eli %7.0f %s %s %s%d\tres=%10.4f\tlock=%4d\tazel=%10.4f\tstd=%10.4f\n",time2gpst(rtk->sol.time,NULL),str,id,
					type==1?"L":"P",freq,v[i],rtk->ssat[sat1-1].lock[freq-1],rtk->ssat[sat1-1].azel[1]*R2D,sqrt(r[i]));
			}
		}
	}
	if (stat == 0) {
		// for (i = 0; i < MAXSAT; i++) {
		// 	for (k = 0; k < NFREQ; k++) {
		// 		rtk->ssat[i].resp[k] = rtk->ssat[i].resc[k] = 1e3;
		// 	}
		// }
		for (i = 0; i < na; i++) {		
			readvflg(vflg[i], &type, &freq, &sat1, &sat2);
			switch (type) {
			case 1: rtk->ssat[sat1 - 1].resc[freq - 1] = v[i]; 
					rtk->ssat[sat1 - 1].vsat[freq - 1] |= 0x2;    break; 		
			case 2: rtk->ssat[sat1 - 1].resp[freq - 1] = v[i];
					rtk->ssat[sat1 - 1].vsat[freq - 1] |= 0x1;    break; 	
			}
		}
	}
	return stat;
}
/* update solution status ----------------------------------------------------*/
static void update_sol(rtk_t *rtk, int stat)
{
	const prcopt_t *opt = &rtk->opt;

	rtk->sol.stat = stat;

	if (rtk->sol.stat == SOLQ_FIX) {
		for (int i = 0; i<3; i++) {
			rtk->sol.rr[i] = rtk->xa[i];
			rtk->sol.pr[i] = rtk->xa[i];
			rtk->sol.qr[i] = (float)rtk->Pa[i + i*rtk->na];
		}
		rtk->sol.qr[3] = (float)rtk->Pa[1];
		rtk->sol.qr[4] = (float)rtk->Pa[1 + 2 * rtk->na];
		rtk->sol.qr[5] = (float)rtk->Pa[2];

		rtk->sol.dtr[0] = rtk->xa[IC(0, opt)] / CLIGHT;
		rtk->sol.dtr[1] = rtk->xa[IC(1, opt)] / CLIGHT;
	}
	else if(rtk->sol.stat == SOLQ_PPP) {	
		for (int i = 0; i<3; i++) {
			rtk->sol.rr[i] = rtk->x[i];
			rtk->sol.pr[i] = rtk->x[i];
			rtk->sol.qr[i] = (float)rtk->P[i + i*rtk->nx];
		}
		rtk->sol.qr[3] = (float)rtk->P[1];
		rtk->sol.qr[4] = (float)rtk->P[2 + rtk->nx];
		rtk->sol.qr[5] = (float)rtk->P[2];

		rtk->sol.dtr[0] = rtk->x[IC(0, opt)] / CLIGHT;
		rtk->sol.dtr[1] = rtk->x[IC(1, opt)] / CLIGHT;
		rtk->sol.dtr[2] = rtk->x[IC(2, opt)] / CLIGHT;
		rtk->sol.dtr[3] = rtk->x[IC(3, opt)] / CLIGHT;
	}
	else {

	}
}
/* precise point positioning -------------------------------------------------*/
extern void ppppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, 
	const double *xg, const double *Pg)
{
	const prcopt_t *opt = &rtk->opt;
	int i, j, k, nv, na, info, *vflg, refs[NSYS], svh[MAXOBS], stat = SOLQ_SINGLE;
	double *rs, *dts, *var, *v, *H, *r, *R, *xp, *Pp;
	exc_t exc = { 0 };
	char str[32];

	time2str(rtk->sol.time, str, 0);

	rs = mat(6, n); dts = mat(2, n); var = mat(1, n);		

	/* initial satellite status */			
	initial_ssat(obs, n, opt, rtk->tt, &rtk->reset, rtk->ssat);							

	/* pseudorange error detect */
	pseudo_detect(obs, n, opt, rtk->ssat);

	/* clock jump detect and repair */
	clkJump_detect(obs, n, nav, opt, rtk->ssat);	

	/* LC observation */  
	linear_combination(obs, n, nav, opt, rtk->ssat);		

	/* cycle jump detect  */					
	cycJump_detect(obs, n, nav, opt, rtk->ssat);

	/* satellite positions and clocks */
	satposs(rtk->sol.time, obs, n, nav, rtk->opt.sateph, rs, dts, var, svh);

	/* exclude measurements of eclipsing satellite (block IIA) */
	if (rtk->opt.posopt[3]) testeclipse(obs, n, nav, rs);			

	/* tdcp check */
	if (rtk->opt.posopt[8]) tdcpos(obs, n, nav, opt, rs, dts, svh, &rtk->sol, rtk->ssat);		
	
	/* temporal update of ekf states */			
	udstate_ppp(rtk, obs, n, nav, xg, Pg);	 
	
	/* earth tides correction */
	if (opt->tidecorr) {
		tidedisp(gpst2utc(obs[0].time), rtk->x, opt->tidecorr == 1 ? 1 : 7, &nav->erp, opt->odisp[0], rtk->dr);
	}

	/* kalman filter */
	nv = n*rtk->opt.nf * 2 + 1 + MAXSAT;					
	xp = mat(rtk->nx, 1); 	Pp = mat(rtk->nx, rtk->nx);		vflg = imat(nv, 1);
	v  = mat(nv, 1);		 H = mat(rtk->nx, nv);			r = mat(nv, 1);			R = mat(nv, nv);

	for (i = 0; i < MAX_ITER; i++) {
		matcpy(xp, rtk->x, rtk->nx, 1);
		matcpy(Pp, rtk->P, rtk->nx, rtk->nx);

		/* reject obs by pre-fit residuals */
		nv  = ppp_res(obs, n, rs, dts, var, svh, nav, xp, rtk, v, H, r, vflg, refs, &exc, 0);				
		nv += cor_res(obs, n, refs, nav, rtk, xp, v + nv, H + nv*rtk->nx, r + nv, vflg + nv, &exc);		
		if (opt->posopt[6]) nv = valpre(rtk, v, H, r, vflg, nv, rtk->nx, &exc);
		nv += pri_res(rtk, xp, Pp, xg, Pg, v + nv, H + nv*rtk->nx, r + nv);

		/* measurement update of ekf states */
		for (j = 0; j < nv; j++) for (k = 0; k < nv; k++) {
			R[j + k*nv] = j == k ? r[j] : 0.0;
		}
		if ((info = filter(xp, Pp, H, v, R, rtk->nx, nv))) {									
			trace(2, "%s ppp (%d) filter error info=%d\n", str, i + 1, info);
			break;
		}

		/* reject obs by pos-fit residuals */
		na  = ppp_res(obs, n, rs, dts, var, svh, nav, xp, rtk, v, NULL, r, vflg, NULL, &exc, 0); 
		// na += cor_res(obs, n, refs, nav, rtk, xp, v + na, NULL, r + na, vflg + na, &exc);
		if (!valpos(rtk, v, r, vflg, na, opt->threscheck[1], &exc)) {					
			stat = SOLQ_PPP;		
			break;
		}
	}

	/* update satellite information */
	update_ssat(obs, n, opt, rs, dts, svh, &stat, &exc, &rtk->sol.ns, rtk->ssat);		

	if (stat == SOLQ_PPP) {
		matcpy(rtk->x, xp, rtk->nx, 1);
		matcpy(rtk->P, Pp, rtk->nx, rtk->nx);

		/* ambiguity resolution (ppp) */
		if ((opt->modear >= ARMODE_PPPAR) && pppamb(rtk, obs, n, nav)) {		

			/* reject obs by pos-fit residuals */
			na = ppp_res(obs, n, rs, dts, var, svh, nav, rtk->xa, rtk, v, NULL, r, vflg, NULL, &exc, 0);
			if (!valpos(rtk, v, r, vflg, na, opt->threscheck[5], &exc)) {				
				stat = SOLQ_FIX;
			}
		}
	}

	double tow = time2gpst(rtk->sol.time,NULL);
	//DF-INS: 443915-443920    DF: 443525-443533
	if (stat==SOLQ_FIX &&(tow>443525 && tow<443533)){        
		// stat = SOLQ_PPP;
	}

    /* update solution status */
	update_sol(rtk, stat);

	free(rs); free(dts); free(var);	
	free(xp); free(Pp);	 free(vflg);
	free(v);  free(H);	 free(r);	 free(R);
}