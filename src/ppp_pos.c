/****************************************
GNSS双频非差非组合||消电离层PPP（固定）||三频
*****************************************/
#include "rtklib.h"
#include "ppp_state.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define REJION   	 30
#define MAX_ITER     8               /* max number of iterations */
#define MAX_STD_FIX  0.15            /* max std-dev (3d) to fix solution */
#define MIN_NSAT_SOL 4               /* min satellite number for solution */
#define THRES_REJECT 4.0             /* reject threshold of posfit-res (sigma) */

#define VAR_POS     SQR(10.0)         /* init variance receiver position (m^2) */			
#define VAR_VEL     SQR(10.0)        /* init variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0)        /* init variance of receiver acc ((m/ss)^2) */
#define VAR_CLK     SQR(10.0)         /* init variance receiver clock (m^2) */   
#define VAR_ZWD     SQR(0.12)		 /* init variance ZWD */ 
#define VAR_GRA     SQR(0.01)        /* init variance gradient (m^2) */    
#define VAR_IONO    SQR(30.0)        /* init variance iono-delay */
#define VAR_IFB		SQR(10.0)		 /* init variance ifb (m^2) */		  
#define VAR_BIAS    SQR(60.0)        /* init variance phase-bias (m^2) */  

#define ERR_GLO_IFB SQR( 0.6)        /* variance of glonass ifb */		  
#define ERR_CBIAS   SQR( 0.3)        /* code bias error std (m) */		 

#define GAP_RESION  120              /* default gap to reset ionos parameters (ep) && ifb */		

typedef struct {
	unsigned char n;
	int flg[256];
} exc_t;

/* number of estimated states ------------------------------------------------*/
extern int pppnx(const prcopt_t *opt)
{
	return NX(opt);
}
/* initialize rtk control ------------------------------------------------------*/
extern void pppinit(rtk_t *rtk, const prcopt_t *opt)
{
	sol_t sol0 = { { 0 } };
	ambc_t ambc0 = { { { 0 } } };
	ssat_t ssat0 = { 0 };

	trace(3, "rtkinit :\n");

	rtk->sol = sol0;
	for (int i = 0; i < 6; i++) rtk->rb[i] = 0.0;
	rtk->nx = NX(opt);
	rtk->na = NX(opt);
	rtk->tt = 0.0;
	rtk->x = zeros(rtk->nx, 1);
	rtk->P = zeros(rtk->nx, rtk->nx);
	rtk->xa = zeros(rtk->na, 1);
	rtk->Pa = zeros(rtk->na, rtk->na);
	for (int i = 0; i<MAXSAT; i++) {
		rtk->ambc[i] = ambc0;
		rtk->ssat[i] = ssat0;
	}
	rtk->opt = *opt;
}
/* free rtk control ------------------------------------------------------------*/
extern void pppfree(rtk_t *rtk)
{
	trace(3, "rtkfree :\n");

	rtk->nx = rtk->na = 0;
	free(rtk->x); rtk->x = NULL;
	free(rtk->P); rtk->P = NULL;
	free(rtk->xa); rtk->xa = NULL;
	free(rtk->Pa); rtk->Pa = NULL;
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
/* clock detect and repair */
static void clkJump_detect(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)		
{
	int i, sat;
	gtime_t time = rtk->sol.time;
	double k1, k2, freq, dp, dh, Ns, Na, S, M, clkJump;

	k1 = 1E-6*CLIGHT - 3 * 5;
	k2 = 1E-4;
	
	// clock detect
	Ns = Na = S = M = clkJump = 0;
	for (i = 0; i < n && i < MAXOBS; i++) {			
		sat = obs[i].sat;
		freq = sat2freq(sat, obs[i].code[0], nav);

		if (!rtk->ssat[sat - 1].vs) continue;
		if (obs[i].P[0] == 0.0 || obs[i].L[0] == 0.0) continue;
		if (timediff(time, rtk->ssat[sat - 1].pt[0][0]) > 120 || timediff(time, rtk->ssat[sat - 1].pt[1][0]) > 120) continue;

		Na++;
		dp =  obs[i].P[0] - rtk->ssat[sat - 1].ph[0][0];
		dh = (obs[i].L[0] - rtk->ssat[sat - 1].ph[1][0]) * (CLIGHT / freq);

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
		trace(1, "detect clock jump: detected time=%s clkJump=%8.3f\n", time_str(time, 0), clkJump);
	}

	for (i = 0; i < n && i < MAXOBS; i++) {			
		sat = obs[i].sat;
		if (rtk->ssat[sat - 1].vs) {		
			rtk->ssat[sat - 1].php += clkJump;
		}
	}
}
/* LC observation */
static void lc_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav) 
{
	int i, k, sat;
	double dantr[3] = { 0 }, dants[3] = { 0 }, *azel, L[NFREQ], P[NFREQ], Lc, Pc, freq[NFREQ], alfa, lambda;
	for (i = 0; i < n && i < MAXOBS; i++) {
		sat = obs[i].sat;
		azel = rtk->ssat[sat - 1].azel;
		if (!rtk->ssat[sat - 1].vs) continue;

		corr_meas(obs + i, nav, azel, &rtk->opt, dantr, dants, rtk->ssat[sat - 1].phw, rtk->ssat[sat - 1].php, L, P, &Lc, &Pc);		
		for (k = 0; k < NFREQ; k++) freq[k] = sat2freq(sat, obs[i].code[k], nav);

		for (k = 0; k < 7; k++) rtk->ssat[sat - 1].LC[k] = 0.0;

		// GF:L1-L2
		if (freq[0] != 0.0 && freq[1] != 0.0 && L[0] != 0.0 && L[1] != 0.0) {
			rtk->ssat[sat - 1].LC[0] = L[0] - L[1];
		} 

		// GF:L1-L3
		if (freq[0] != 0.0 && freq[2] != 0.0 && L[0] != 0.0 && L[2] != 0.0) {
			rtk->ssat[sat - 1].LC[1] = L[0] - L[2];
		}

		// MW:L1L2-P1P2
		if (freq[0] != 0.0 && freq[1] != 0.0 && L[0] != 0.0 && L[1] != 0.0 && P[0] != 0.0 && P[1] != 0.0) {
			lambda = CLIGHT / (freq[0] - freq[1]);
			rtk->ssat[sat - 1].LC[2] = (freq[0] * L[0] - freq[1] * L[1]) / CLIGHT - (freq[0] * P[0] + freq[1] * P[1]) / (freq[0] + freq[1]) / lambda;
		}

		// MW:L2L3-P2P3 
		if (freq[1] != 0.0 && freq[2] != 0.0 && L[1] != 0.0 && L[2] != 0.0 && P[1] != 0.0 && P[2] != 0.0) {
			lambda = CLIGHT / (freq[1] - freq[2]);
			rtk->ssat[sat - 1].LC[3] = (freq[1] * L[1] - freq[2] * L[2]) / CLIGHT - (freq[1] * P[1] + freq[2] * P[2]) / (freq[1] + freq[2]) / lambda;
		}

		// MW:L2L3-P1P2
		if (freq[1] != 0.0 && freq[2] != 0.0 && L[1] != 0.0 && L[2] != 0.0 && P[0] != 0.0 && P[1] != 0.0) {
			alfa = (SQR(freq[0])*SQR(freq[1]) - SQR(freq[1])*freq[1] * freq[2]) / (freq[1] * freq[2] * (SQR(freq[0]) - SQR(freq[1])));
			lambda = CLIGHT / (freq[1] - freq[2]);
			rtk->ssat[sat - 1].LC[4] = (freq[1] * L[1] - freq[2] * L[2]) / CLIGHT - ((1 - alfa) * P[0] + alfa* P[1]) / lambda;				//一个错误找了多久，兄弟		
		}

		// gfPL
		if (freq[0] != 0.0 && freq[1] != 0.0 && L[0] != 0.0 && L[1] != 0.0 && P[0] != 0.0 && P[1] != 0.0) {
			rtk->ssat[sat - 1].LC[5] = L[0] - L[1] + P[0] - P[1];
		}

		// MP
		if (P[0] != 0.0 && freq[0] != 0.0 && freq[1] != 0.0 && L[0] != 0.0 && L[1] != 0.0) {
			rtk->ssat[sat - 1].LC[6] = P[0] - (SQR(freq[0]) + SQR(freq[1])) / (SQR(freq[0]) - SQR(freq[1]))*L[0] + 2 * SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]))*L[1];
		}
	}
}
/* cycle Jump detect ------------------------------------------------------- */
static void cycJump_detect(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	int i, j, k, sat, slip, lli, mw12, mw23, gf12, gf23;
	double el, thresgf, thresmw, diffamb;

	for (i = 0; i < MAXSAT; i++) for (j = 0; j < rtk->opt.nf; j++) rtk->ssat[i].slip[j] = 0;

	/* LC组合观测值 */
	for (i = 0; i < n && i < MAXOBS; i++) {
		slip = lli = mw12 = mw23 = gf12 = gf23 = 0;
		sat = obs[i].sat;
		el = rtk->ssat[sat - 1].azel[1] / PI * 180;
		if (!rtk->ssat[sat - 1].vs) continue;

		if (el >= 15) {
			thresgf = rtk->opt.thresslipgf;
			thresmw = rtk->opt.thresslipmw;
		}
		else {
			thresgf = (7 - 0.4*el)*rtk->opt.thresslipgf;
			thresmw = (4 - 0.2*el)*rtk->opt.thresslipmw;
		}
		for (k = 0; k < 6; k++) {
			if (fabs(timediff(rtk->sol.time, rtk->ambc[sat - 1].epoch[k])) + DTTOL > 300) {
				rtk->ambc[sat - 1].n[k] = 0;
				rtk->ambc[sat - 1].LC[k] = rtk->ambc[sat - 1].LCv[k] = 0.0;
			}
		}

		if (rtk->ssat[sat - 1].LC[0] != 0.0 && rtk->ssat[sat - 1].LC[1] != 0.0	&& rtk->ssat[sat - 1].LC[3] != 0.0) {	// 三频周跳探测
			// LLI
			lli = (obs[i].LLI[0] & 3 || obs[i].LLI[1] & 3 || obs[i].LLI[2] & 3) ? 1 : 0;
			
			// GF L1L2
			diffamb = rtk->ssat[sat - 1].LC[0] - rtk->ambc[sat - 1].LC[0];
			gf12 = (rtk->ambc[sat - 1].n[0] == 0.0 || (rtk->ambc[sat - 1].n[0] != 0.0 && fabs(diffamb)>thresgf)) ? 1 : 0;

			// GF L2L3
			diffamb = rtk->ssat[sat - 1].LC[1] - rtk->ambc[sat - 1].LC[1];
			gf23 = (rtk->ambc[sat - 1].n[1] == 0.0 || (rtk->ambc[sat - 1].n[1] != 0.0 && fabs(diffamb)>thresgf)) ? 1 : 0;

			// MW L2L3 P2P3
			diffamb = rtk->ssat[sat - 1].LC[3] - rtk->ambc[sat - 1].LC[3];
			mw23 = (rtk->ambc[sat - 1].n[3] == 0.0 || (rtk->ambc[sat - 1].n[3] > 0 && fabs(diffamb) > thresmw)) ? 1 : 0;
		}
		else if (rtk->ssat[sat - 1].LC[0] != 0.0 && rtk->ssat[sat - 1].LC[2] != 0.0) {									// 双频周跳探测
			// LLI
			lli = (obs[i].LLI[0] & 3 || obs[i].LLI[1] & 3) ? 1 : 0;

			// GF L1L2
			diffamb = rtk->ssat[sat - 1].LC[0] - rtk->ambc[sat - 1].LC[0];
			gf12 = (rtk->ambc[sat - 1].n[0] == 0.0 || (rtk->ambc[sat - 1].n[0] != 0.0 && fabs(diffamb)>thresgf)) ? 1 : 0;

			// MW L1L2 P1P2
			diffamb = rtk->ssat[sat - 1].LC[2] - rtk->ambc[sat - 1].LC[2];
			mw12 = (rtk->ambc[sat - 1].n[2] == 0.0 || (rtk->ambc[sat - 1].n[2] > 0 && fabs(diffamb) > thresmw)) ? 1 : 0;
			
		}
		else if (rtk->ssat[sat - 1].LC[0] != 0.0) {																		// 单频周跳探测
			// LLI
			lli = (obs[i].LLI[0] & 3) ? 1 : 0;

			// L(t)-L(t-1) - (P(t)-P(t-1))

			// L(t)-L(t-1)- deltaT*(D(t)+D(t-1))/2			

		}

		if (lli)  trace(2, "LLI		    : slip detected time=%s sat=%2d \n", time_str(obs[i].time, 0), sat);
		if (gf12) trace(2, "detslip_gf12: slip detected time=%s sat=%2d gf=%8.3f->%8.3f\n", time_str(obs[i].time, 0), obs[i].sat, rtk->ambc[sat - 1].LC[0], rtk->ssat[sat - 1].LC[0]);
		if (gf23) trace(2, "detslip_gf23: slip detected time=%s sat=%2d gf=%8.3f->%8.3f\n", time_str(obs[i].time, 0), obs[i].sat, rtk->ambc[sat - 1].LC[1], rtk->ssat[sat - 1].LC[1]);
		if (mw12) trace(2, "detslip_mw12: slip detected time=%s sat=%2d num=%4d mw=%8.3f->%8.3f\n", time_str(obs[i].time, 0), obs[i].sat, rtk->ambc[sat - 1].n[2], 
			rtk->ambc[sat - 1].LC[2], rtk->ssat[sat - 1].LC[2]);
		if (mw23) trace(2, "detslip_mw23: slip detected time=%s sat=%2d num=%4d mw=%8.3f->%8.3f\n", time_str(obs[i].time, 0), obs[i].sat, rtk->ambc[sat - 1].n[3], 
			rtk->ambc[sat - 1].LC[3], rtk->ssat[sat - 1].LC[3]);

		// 周跳标记
		slip = lli || mw12 || mw23 || gf12 || gf23;
		for (k = 0; k < rtk->opt.nf; k++) rtk->ssat[sat - 1].slip[k] = slip;
	}
}
/* temporal update of position --------------------------------------------- */
static void udpos_ppp(rtk_t *rtk)			
{
	double *F, *P, *FP, *x, *xp, pos[3], Q[9] = { 0 }, Qv[9];
	int i, j, *ix, nx;

	trace(3, "udpos_ppp:\n");

	/* fixed mode */
	if (rtk->opt.mode == PMODE_PPP_FIXED) {
		for (i = 0; i<3; i++) initx(rtk, rtk->opt.ru[i], 1E-8, i);
		return;
	}
	/* initialize position for first epoch */									
	if (norm(rtk->x, 3) <= 0.0) {				
		for (i = 0; i<3; i++) initx(rtk, rtk->sol.rr[i], VAR_POS, i);			//initx(rtk, rtk->sol.rr[i], VAR_POS, i);			
		if (rtk->opt.dynamics) {												//动态模型（匀加速模型）
			for (i = 3; i<6; i++) initx(rtk, rtk->sol.rr[i], VAR_VEL, i);
			for (i = 6; i<9; i++) initx(rtk, 1E-6, VAR_ACC, i);
		}
	}
	/* static ppp mode */
	if (rtk->opt.mode == PMODE_PPP_STATIC) {
		for (i = 0; i<3; i++) {
			rtk->P[i*(1 + rtk->nx)] += SQR(rtk->opt.prn[5])*fabs(rtk->tt);
		}
		return;
	}
	/* kinmatic mode without dynamics */
	if (!rtk->opt.dynamics) {
		for (i = 0; i<3; i++) {
			initx(rtk, rtk->sol.rr[i], VAR_POS, i);								//动态ppp
		}
		return;
	}
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
	matmul("NN", nx, 1, nx, 1.0, F, x, 0.0, xp);
	matmul("NN", nx, nx, nx, 1.0, F, P, 0.0, FP);
	matmul("NT", nx, nx, nx, 1.0, FP, F, 0.0, P);

	for (i = 0; i<nx; i++) {
		rtk->x[ix[i]] = xp[i];
		for (j = 0; j<nx; j++) {
			rtk->P[ix[i] + ix[j] * rtk->nx] = P[i + j*nx];
		}
	}
	/* process noise added to only acceleration */
	Q[0] = Q[4] = SQR(rtk->opt.prn[3])*fabs(rtk->tt);
	Q[8] = SQR(rtk->opt.prn[4])*fabs(rtk->tt);
	ecef2pos(rtk->x, pos);
	covecef(pos, Q, Qv);
	for (i = 0; i<3; i++) for (j = 0; j<3; j++) {
		rtk->P[i + 6 + (j + 6)*rtk->nx] += Qv[i + j * 3];
	}
	free(ix); free(F); free(P); free(FP); free(x); free(xp);
}
/* temporal update of clock --------------------------------------------------*/
static void udclk_ppp(rtk_t *rtk)
{
	double dtr, isb; 
	trace(3, "udclk_ppp:\n");
	/* initialize every epoch for clock (white noise) */
	for (int i = 0; i < 4; i++) {
		double dtr = rtk->sol.dtr[i];
		rtk->sol.dsp[i] = CLIGHT*rtk->sol.dtr[i];
		initx(rtk, CLIGHT*dtr, VAR_CLK, IC(i, &rtk->opt));
	}
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(rtk_t *rtk)
{
	const double azel[] = { 0.0, PI / 2.0 };
	double pos[3], zwd, var;
	int i = IT(&rtk->opt), j;

	trace(3, "udtrop_ppp:\n");

	if (rtk->x[i] == 0.0) {
		ecef2pos(rtk->sol.rr, pos);
		zwd = tropmodel(rtk->sol.time, pos, azel, 1, REL_HUMI, 2);					//湿延迟初值,约0.1m
		initx(rtk, zwd, VAR_ZWD, i);

		if (rtk->opt.tropopt >= TROPOPT_ESTG) {
			for (j = i + 1; j<i + 3; j++) initx(rtk, 1E-6, VAR_GRA, j);
		}
	}
	else {											
		rtk->P[i + i*rtk->nx] += SQR(rtk->opt.prn[2])*fabs(rtk->tt);

		if (rtk->opt.tropopt >= TROPOPT_ESTG) {
			for (j = i + 1; j<i + 3; j++) {
				rtk->P[j + j*rtk->nx] += SQR(rtk->opt.prn[2] * 0.1)*fabs(rtk->tt);
			}
		}
	}
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void udiono_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)					
{
	int i, k, ii, sat;
	double freq[NFREQ], ion, sinel, pos[3], dantr[3] = { 0 }, dants[3] = { 0 }, *azel, L[NFREQ], P[NFREQ], Lc, Pc;
	
	trace(3, "udiono_ppp:\n");

	for (i = 0; i<MAXSAT; i++) {
		ii = II(i + 1, &rtk->opt);
		if (rtk->x[ii] != 0.0) {	
			sinel = sin(MAX(rtk->ssat[i].azel[1], 15.0*D2R));
			rtk->P[ii + ii*rtk->nx] += SQR(rtk->opt.prn[1] / sinel)*fabs(rtk->tt);
			if (rtk->ssat[i].outc[0] > GAP_RESION) initx(rtk, 0.0, 0.0, ii);
		}
	}
	ecef2pos(rtk->sol.rr, pos);
	for (i = 0; i<n; i++) {
		sat = obs[i].sat;
		if (!rtk->ssat[sat - 1].vs) continue;
		ii = II(sat, &rtk->opt);
		azel = rtk->ssat[sat - 1].azel;

		for (k = 0; k < NFREQ; k++) freq[k] = sat2freq(sat, obs[i].code[k], nav);	

		if (rtk->x[ii] == 0.0) {
			corr_meas(obs + i, nav, azel, &rtk->opt, dantr, dants, rtk->ssat[sat - 1].phw, rtk->ssat[sat - 1].php, L, P, &Lc, &Pc);				
			if (obs[i].P[0] != 0.0 && obs[i].P[1] != 0.0 && freq[0] != 0.0 && freq[1] != 0.0) {		// 双频IONO初始化
				ion = (obs[i].P[0] - obs[i].P[1]) / (1.0 - SQR(freq[0] / freq[1]));
			}
			else if (obs[i].P[0] != 0.0 && freq[0] != 0.0) {										// 单频IONO初始化
				ion = ionmodel(obs[i].time, nav->ion_gps, pos, azel);			
				ion *= SQR(FREQ1 / freq[0]);
			}
			if (ion != 0.0) initx(rtk, ion, VAR_IONO, ii);
		}
	}
}
/* temporal update of inter frequency bias ---------------------------------- */
static void udifb_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)		
{
	int i, ii, sat;
	for (i = 0; i < MAXSAT; i++) {
		ii = IS(i + 1, &rtk->opt);
		if (rtk->x[ii] != 0.0) {
			if ((int)rtk->ssat[i].outc[2] > GAP_RESION) initx(rtk, 0.0, 0.0, ii);
		}
	}
	for (i = 0; i < n; i++) {
		sat = obs[i].sat;
		if (!rtk->ssat[sat - 1].vs) continue;
		ii = IS(sat, &rtk->opt);
		if (rtk->x[ii] == 0.0) initx(rtk, 1e-6, VAR_IFB, ii);    //1e-6(这里的初始化需要考虑下)
	}
}
/* temporal update of phase biases ------------------------------------------- */
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)											
{
	int i, j, f, sat, reset;
	double L[NFREQ], P[NFREQ], Lc, Pc, freq[NFREQ], bias[MAXOBS], gamma;
	double ion, dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };

	trace(3, "udbias  : n=%d\n", n);

	for (f = 0; f<NF(&rtk->opt); f++) {												
		/* reset phase-bias if expire obs outage counter */
		for (i = 0; i<MAXSAT; i++) {
			j = IB(i + 1, f, &rtk->opt);
			reset = rtk->ssat[i].outc[f] > (unsigned int)rtk->opt.maxout;							//连续失锁5历元重置	
			if (rtk->x[j] != 0.0) {	
				rtk->P[j + j*rtk->nx] += SQR(rtk->opt.prn[0])*fabs(rtk->tt);
				if (reset) initx(rtk, 0.0, 0.0, j);
			}
		}
		for (i = 0; i<n&&i<MAXOBS; i++) {
			sat = obs[i].sat;
			j = IB(sat, f, &rtk->opt);

			corr_meas(obs + i, nav, rtk->ssat[sat - 1].azel, &rtk->opt, dantr, dants, rtk->ssat[sat - 1].phw, rtk->ssat[sat - 1].php, L, P, &Lc, &Pc);

			bias[i] = 0.0;
			if (rtk->x[j] != 0 && !rtk->ssat[sat - 1].slip[f]) continue;

			if (rtk->opt.ionoopt == IONOOPT_IFLC) {
				if (Lc == 0.0 || Pc == 0.0) continue;
				bias[i] = Lc - Pc;
			}
			else if (rtk->opt.ionoopt == IONOOPT_EST) {
				freq[0] = sat2freq(sat, obs[i].code[0], nav);
				freq[f] = sat2freq(sat, obs[i].code[f], nav);
				ion = rtk->x[II(sat, &rtk->opt)];
				if (L[f] == 0.0 || P[f] == 0.0 || freq[f]==0.0 || ion == 0.0) continue;
		
				gamma = SQR(freq[0] / freq[f]);					
				bias[i] = L[f] - P[f] + 2.0*ion*gamma;
			}
			rtk->x[j] = 0.0;
			if (bias[i] != 0.0) initx(rtk, bias[i], VAR_BIAS, j);
		}
	}//end frequency
}
/* temporal update of states --------------------------------------------------*/
static void udstate_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	trace(3, "udstate_ppp: n=%d\n", n);

	/* temporal update of position */
	udpos_ppp(rtk);

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
static int exclude(int itype, int ifreq, int isat1, int isat2, exc_t *exc) 
{
	int i, type, freq, sat1, sat2, stat = 0;
	for (i = 0; i < exc->n; i++) {
		type = exc->flg[i] & 0xF;
		freq = (exc->flg[i] >> 4) & 0xF;
		sat1 = (exc->flg[i] >> 8) & 0xFF;
		sat2 = (exc->flg[i] >> 16) & 0xFF;
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
static int ppp_res(const obsd_t *obs, int n, const double *rs, const double *dts, const double *var_rs, const int *svh, 
	const nav_t *nav, const double *x, rtk_t *rtk, double *v, double *H, double *var, int *vflg, const exc_t *exc)
{
	int i, j, k, sat, sys, flg, nv = 0, nsat = 0, nx = rtk->nx;
	char str[32], id[32];
	prcopt_t *opt = &rtk->opt;
	double freq[NFREQ] = { 0 }, y, r, cdtr, bias, C, rr[3], pos[3], e[3], dtdx[3], L[NFREQ], P[NFREQ], Lc, Pc;
	double dtrp = 0.0, dion = 0.0, vart = 0.0, vari = 0.0, ifb;
	double dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 }, *azel, *hh, vv, qq;

	time2str(obs[0].time, str, 0);

	hh = mat(nx, 1);

	for (i = 0; i<3; i++) rr[i] = x[i] + rtk->dr[i];
	ecef2pos(rr, pos);
	for (i = 0; i<n&&i<MAXOBS; i++) {
		sat = obs[i].sat;
		satno2id(sat, id);
		azel = rtk->ssat[sat - 1].azel;

		for (j = 0; j < NFREQ; j++) freq[j] = sat2freq(sat, obs[i].code[j], nav);			

		if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 || satazel(pos, e, azel)<opt->elmin) {			
			continue;
		}
		if (!(sys = satsys(sat, NULL)) || !rtk->ssat[sat - 1].vs || satexclude(obs[i].sat, svh[i], opt)) {							
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
		if (!model_phw(rtk->sol.time, sat, nav->pcvs[sat - 1].type, opt->posopt[2] ? 2 : 0, rs + i * 6, rr, &rtk->ssat[sat - 1].phw)) {	
			continue;
		}

		/* corrected phase and code measurements */
		corr_meas(obs + i, nav, azel, &rtk->opt, dantr, dants, rtk->ssat[sat - 1].phw, rtk->ssat[sat - 1].php, L, P, &Lc, &Pc);

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

			/* coordinate */
			for (k = 0; k < nx; k++) hh[k] = k < 3 ? -e[k] : 0.0;

			/* receiver clock */
			k = sys == SYS_GLO ? 1 : (sys == SYS_GAL ? 2 : (sys == SYS_CMP ? 3 : 0));				
			cdtr = x[IC(k, opt)];
			hh[IC(k, opt)] = 1.0;

			/* Trop */
			if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
				for (k = 0; k<(opt->tropopt >= TROPOPT_ESTG ? 3 : 1); k++) {
					hh[IT(opt) + k] = dtdx[k];
				}
			}

			/* Iono */
			C = SQR(freq[0] / freq[j / 2])*(j % 2 == 0 ? -1.0 : 1.0);	
			if (opt->ionoopt == IONOOPT_EST) {
				if (rtk->x[II(sat, opt)] == 0.0) continue;
				hh[II(sat, opt)] = C;
			}

			/* phase bias */
			if (j % 2 == 0) {
				if ((bias = x[IB(sat, j / 2, opt)]) == 0.0) continue;
				hh[IB(sat, j / 2, opt)] = 1.0;
			}

			/* Inter frequency bias */
			if (j % 2 != 0 && j / 2 > 1) {
				if ((ifb = x[IS(sat, opt)]) == 0.0) continue;
				hh[IS(sat, opt)] = 1.0;
			}

			/* residual */
			vv = y - (r + cdtr - CLIGHT*dts[i * 2] + dtrp + C*dion + ifb + bias);

			/* variance */
			if (j % 2 == 0) qq = varerr(obs[i].sat, azel[1], j / 2, j % 2, opt) + vart + SQR(C)*vari + var_rs[i];
			if (j % 2 == 1)	qq = varerr(obs[i].sat, azel[1], j / 2, j % 2, opt) + vart + SQR(C)*vari + var_rs[i] + ERR_CBIAS;
			if (sys == SYS_GLO && j % 2 == 1) qq += ERR_GLO_IFB;

			if (v) v[nv] = vv;
			if (var) var[nv] = qq;
			if (H) matcpy(H + nv*nx, hh, nx, 1);
			vflg[nv++] = (sat << 8) | ((j / 2 + 1) << 4) | (j % 2) + 1;
		}	//end styple
	}	//end obs
	free(hh);
	return nv;
}
/* select reference satellite ------------------------------------------------*/
static void selsat(const ssat_t *ssat, const int *vflg, const int nv, int *iref)
{
	int i, k, sys, sat, type;
	double maxAzel[NSYS], *azel;

	for (i = 0; i < NSYS; i++) { maxAzel[i] = 0; iref[i] = -1; }
	for (i = 0; i < nv; i++) {
		type = vflg[i] & 0xF;
		if (type == 1) {
			sat = (vflg[i] >> 8) & 0xFF;
			sys = satsys(sat, NULL);
			azel = ssat[sat - 1].azel;
			k = (sys == SYS_GPS) ? 0 : (sys == SYS_GLO) ? 1 : (sys == SYS_GAL) ? 2 : (sys == SYS_CMP) ? 3 : -1;
			if (k > 0 && maxAzel[k] < azel[1]) {
				maxAzel[k] = azel[1];
				iref[k] = sat;				
			}
		}
	}
}
/* constraint to local correction --------------------------------------------*/
static int cor_res(const obsd_t *obs, const int n, const int *refs, const nav_t *nav, const rtk_t *rtk, 
	const double *x, double *v, double *H, double *var, int *vflg, const exc_t *exc)
{
	int i, k, ii, jj, sat1, sat2, sys, nx, nv = 0;
	double trop[3], std_trop[3], zhd, iono, std_iono, fact, freq;
	double *azel, rr[3], pos[3], vv, qq, *hh;
	const double zazel[] = { 0.0, PI / 2.0 };

	nx = rtk->nx;			hh = mat(nx, 1);
	for (i = 0; i < 3; i++) rr[i] = x[i] + rtk->dr[i];	ecef2pos(rr, pos);
	/* constraint to external troposphere correction */
	if (rtk->opt.posopt[4]) {
		if (!exclude(4, NULL, NULL, NULL, exc)) {
			if (pppcorr_trop(obs[0].time, pos, trop, std_trop)) {
				ii = IT(&rtk->opt);
				zhd = tropmodel(obs[0].time, pos, zazel, 1, 0, 1);

				vv = (trop[0] - zhd) - x[ii];
				qq = SQR(std_trop[0]);
				for (k = 0; k < rtk->nx; k++) hh[k] = k == ii ? 1.0 : 0.0;

				if (v) v[nv] = vv;
				if (var) var[nv] = qq;
				if (H) matcpy(H + nv*nx, hh, nx, 1);
				vflg[nv++] = 4;
			}
		}
	}
	/* constraint to external ionosphere correction */
	if (rtk->opt.posopt[5]) {
		for (i = 0; i < n; i++) {
			sat1 = obs[i].sat;
			sys = satsys(sat1, NULL);
			azel = rtk->ssat[sat1 - 1].azel;
			freq = sat2freq(obs[i].sat, obs[i].code[0], nav);
			fact = 40.30E16 / SQR(freq);

			k = (sys == SYS_GPS) ? 0 : (sys == SYS_GLO) ? 1 : (sys == SYS_GAL) ? 2 : (sys == SYS_CMP) ? 3 : -1;
			if (k < 0 || refs[k] < 0 || !rtk->ssat[sat1 - 1].vs || azel[1] < 15 * D2R) continue;

			sat2 = refs[k]; 

			if (exclude(8, NULL, sat1, sat2, exc)) continue;

			ii = II(sat1, &rtk->opt);
			jj = II(sat2, &rtk->opt);
			if (!pppcorr_stec(obs[i].time, pos, sat1, sat2, &iono, &std_iono)) continue;

			vv = iono - (x[ii] - x[jj]) / fact;
			qq = SQR(std_iono);
			for (k = 0; k < rtk->nx; k++) hh[k] = (k == ii) ? 1 / fact : (k == jj) ? -1 / fact : 0.0;

			if (v) v[nv] = vv;
			if (var) var[nv] = qq;
			if (H) matcpy(H + nv*nx, hh, nx, 1);
			vflg[nv++] = (sat1 << 16) | (sat2 << 8) | 8;
		}
	}
	free(hh);
	return nv;
}
/* reject obs by pre-fit residuals */
static int valpre(double *v, double *H, double *r, int *vflg, int nv, int nx, exc_t *exc) 
{
	int i, type, stat, nn = 0;
	for (i = 0; i < nv; i++) {
		stat = 0;
		type = vflg[i] & 0xF;
		switch (type) {
		case 1: if (fabs(v[i]) > REJION) { exc->flg[exc->n++] = vflg[i]; stat = 1; } break;
		case 2: if (fabs(v[i]) > REJION) { exc->flg[exc->n++] = vflg[i]; stat = 1; } break;
		case 4: if (fabs(v[i]) > REJION) { exc->flg[exc->n++] = vflg[i]; stat = 1; } break;
		case 8: if (fabs(v[i]) > REJION) { exc->flg[exc->n++] = vflg[i]; stat = 1; } break;
		}
		if (!stat) {
			v[nn] = v[i];		matcpy(H + nn*nx, H + i*nx, nx, 1);		
			r[nn] = r[i];		vflg[nn++] = vflg[i];
		}
	}
	return nn;
}
/* reject obs by pos-fit residuals */
static int valpos(rtk_t *rtk, double *v, double *r, int *vflg, int na, double thres, exc_t *exc)
{
	int i, k, type, sat, freq, stat = 0;
	for (i = 0; i < na; i++) {
		type = vflg[i] & 0xF;
		if (fabs(v[i]) > sqrt(r[i])*thres) {
			exc->flg[exc->n++] = vflg[i];
			stat++;
		}
	}
	if (stat == 0) {
		for (i = 0; i < MAXSAT; i++) {
			for (k = 0; k < rtk->opt.nf; k++) {
				rtk->ssat[i].resp[k] = rtk->ssat[i].resc[k] = rtk->ssat[i].vsat[k] = 0;
			}
		}
		for (i = 0; i < na; i++) {
			type = vflg[i] & 0xF;
			if (type & 3) {
				sat  = (vflg[i] >> 8) & 0xFF;
				freq = (vflg[i] >> 4) & 0xF;
				switch (type) {
				case 1: rtk->ssat[sat - 1].resc[freq - 1] = v[i]; rtk->ssat[sat - 1].vsat[freq - 1] = 1;  break; 
				case 2: rtk->ssat[sat - 1].resp[freq - 1] = v[i]; break;
				}
			}
		}
	}
	return stat;
}
/* update information --------------------------------------------------------*/
static void update_satinfo(rtk_t *rtk, const obsd_t *obs, int n) 
{
	int i, j, k, sat, slip, num;
	double sig0, diffamb, LC;

	for (i = 0; i < n && i<MAXOBS; i++) {
		sat = obs[i].sat;
		slip = rtk->ssat[obs[i].sat - 1].slip[0];
		if (!rtk->ssat[sat - 1].vs) continue;

		/* update satellite snr */
		rtk->ssat[obs[i].sat - 1].snr[0] = MIN(obs[i].SNR[0], obs[i].SNR[1]);

		/* update saatellite slip cout */
		if (slip) rtk->ssat[sat - 1].slipc[0]++;

		/* update ambiguity message */
		for (j = 0; j < 6; j++) {
			LC = rtk->ssat[sat - 1].LC[j];
			if (LC) {
				rtk->ambc[sat - 1].epoch[j] = rtk->sol.time;
			}
			if (slip || LC == 0.0 || j < 2) {
				rtk->ambc[sat - 1].n[j] = LC == 0.0 ? 0 : 1;
				rtk->ambc[sat - 1].LC[j] = LC;
				rtk->ambc[sat - 1].LCv[j] = LC == 0.0 ? 0 : SQR(0.5);
			}
			else {
				num = ++rtk->ambc[sat - 1].n[j];
				diffamb = LC - rtk->ambc[sat - 1].LC[j];
				rtk->ambc[sat - 1].LC[j] += (diffamb / num);
				sig0 = rtk->ambc[sat - 1].LCv[j];
				rtk->ambc[sat - 1].LCv[j] = num == 1 ? SQR(0.5) : sig0 + (SQR(diffamb) - sig0) / num;
			}
		}

		/* update previous pseudorange & carrier-phase */
		for (k = 0; k < NFREQ; k++) {
			if (obs[i].P[k] != 0) {
				rtk->ssat[sat - 1].ph[0][k] = obs[i].P[k];
				rtk->ssat[sat - 1].pt[0][k] = obs[i].time;
			}
			if (obs[i].L[k] != 0) {
				rtk->ssat[sat - 1].ph[1][k] = obs[i].L[k];
				rtk->ssat[sat - 1].pt[1][k] = obs[i].time;
			}
		}
	}
}
/* update solution status ----------------------------------------------------*/
static int update_stat(rtk_t *rtk, const obsd_t *obs, int n, int stat)
{
	int i, k, sat;
	const prcopt_t *opt = &rtk->opt;

	/* test # of valid satellites */
	rtk->sol.ns = 0;
	for (i = 0; i < MAXSAT; i++) {
		for (k = 0; k < NF(opt); k++) {
			if (rtk->ssat[i].vsat[k]) {				
				rtk->ssat[i].lock[k]++;		
				rtk->ssat[i].outc[k] = 0;
			}
			else {
				rtk->ssat[i].lock[k] = 0;
				rtk->ssat[i].outc[k]++;				
			}
		}
		if (rtk->ssat[i].vsat[0]) rtk->sol.ns++;
	}

	rtk->sol.stat = rtk->sol.ns<MIN_NSAT_SOL ? SOLQ_NONE : stat;

	if (rtk->sol.stat == SOLQ_FIX) {
		for (i = 0; i<3; i++) {
			rtk->sol.rr[i] = rtk->xa[i];
			rtk->sol.qr[i] = (float)rtk->Pa[i + i*rtk->na];
		}
		rtk->sol.qr[3] = (float)rtk->Pa[1];
		rtk->sol.qr[4] = (float)rtk->Pa[1 + 2 * rtk->na];
		rtk->sol.qr[5] = (float)rtk->Pa[2];

		rtk->sol.dtr[0] = rtk->xa[IC(0, opt)] / CLIGHT;
		rtk->sol.dtr[1] = rtk->xa[IC(1, opt)] / CLIGHT;
	}
	else if(rtk->sol.stat == SOLQ_PPP){
		for (i = 0; i<3; i++) {
			rtk->sol.rr[i] = rtk->x[i];
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
extern void pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	const prcopt_t *opt = &rtk->opt;
	char str[32];
	int i, j, k, f, nv, na, info, *vflg, refs[NSYS], svh[MAXOBS], stat = SOLQ_SINGLE;
	double *rs, *dts, *var, *v, *H, *r, *R, *xp, *Pp;
	exc_t exc = { { 0 } }, exc1;
	
	time2str(obs[0].time, str, 2);
	trace(2, "pppos   : time=%s nx=%d n=%d\n", str, rtk->nx, n);

	rs = mat(6, n); dts = mat(2, n); var = mat(1, n);
	for (i = 0; i < MAXSAT; i++) for (j = 0; j < opt->nf; j++) rtk->ssat[i].fix[j] = 0;			

	/* clock jump detect and repair */
	clkJump_detect(rtk, obs, n, nav);

	/* LC observation */
	lc_ppp(rtk, obs, n, nav);	

	/* cycle jump detect  */					
	cycJump_detect(rtk, obs, n, nav);			

	/* temporal update of ekf states */			
	udstate_ppp(rtk, obs, n, nav);			

	/* update ambiguity information */
	update_satinfo(rtk, obs, n);		

	/* satellite positions and clocks */
	satposs(obs[0].time, obs, n, nav, rtk->opt.sateph, rs, dts, var, svh);	

	/* exclude measurements of eclipsing satellite (block IIA) */
	if (rtk->opt.posopt[3]) testeclipse(obs, n, nav, rs);

	/* earth tides correction */
	if (opt->tidecorr) tidedisp(gpst2utc(obs[0].time), rtk->x, opt->tidecorr == 1 ? 1 : 7, &nav->erp, opt->odisp[0], rtk->dr);

	/* kalman filter */
	nv = n*rtk->opt.nf * 2 + 1 + MAXSAT;					
	xp = mat(rtk->nx, 1); Pp = mat(rtk->nx, rtk->nx);	vflg = imat(nv, 1);
	v = mat(nv, 1);		  H = mat(rtk->nx, nv);			r = mat(nv, 1);			R = mat(nv, nv);

	for (i = 0; i < MAX_ITER; i++) {
		matcpy(xp, rtk->x, rtk->nx, 1);
		matcpy(Pp, rtk->P, rtk->nx, rtk->nx);

		/* reject obs by pre-fit residuals */
		nv = ppp_res(obs, n, rs, dts, var, svh, nav, xp, rtk, v, H, r, vflg, &exc);
		selsat(rtk->ssat, vflg, nv, refs);
		nv += cor_res(obs, n, refs, nav, rtk, xp, v + nv, H + nv*rtk->nx, r + nv, vflg + nv, &exc);
		nv = valpre(v, H, r, vflg, nv, rtk->nx, &exc);

		/* measurement update of ekf states */
		for (j = 0; j < nv; j++) for (k = 0; k < nv; k++) {
			R[j + k*nv] = j == k ? r[j] : 0.0;
		}
		if ((info = filter(xp, Pp, H, v, R, rtk->nx, nv))) {									
			trace(2, "%s ppp (%d) filter error info=%d\n", str, i + 1, info);
			break;
		}

		/* reject obs by pos-fit residuals */
		na  = ppp_res(obs, n, rs, dts, var, svh, nav, xp, rtk, v, NULL, r, vflg, &exc); 
		na += cor_res(obs, n, refs, nav, rtk, xp, v + na, NULL, r + na, vflg + na, &exc);
		if (!valpos(rtk, v, r, vflg, na, 4.0, &exc)) {					
			stat = SOLQ_PPP;		
			break;
		}
	}
	
	if (stat == SOLQ_PPP) {
		matcpy(rtk->x, xp, rtk->nx, 1);
		matcpy(rtk->P, Pp, rtk->nx, rtk->nx);

		/* ambiguity resolution (ppp) */
		if ((opt->modear >= ARMODE_PPPAR) && pppamb(rtk, obs, n, nav, refs)) {

			/* reject obs by pos-fit residuals */
			na = ppp_res(obs, n, rs, dts, var, svh, nav, rtk->xa, rtk, v, NULL, r, vflg, &exc);
			na += cor_res(obs, n, refs, nav, rtk, rtk->xa, v + na, NULL, r + na, vflg + na, &exc);
			if (!valpos(rtk, v, r, vflg, na, 4.0, &exc)) {				
				stat = SOLQ_FIX;
			}
		}
	}

    /* update solution status */
	update_stat(rtk, obs, n, stat);

	free(rs); free(dts); free(var);	
	free(xp); free(Pp);	 free(vflg);
	free(v);  free(H);	 free(r);	 free(R);
}
extern int ppprocess(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)								
{
	prcopt_t *opt = &rtk->opt;
	gtime_t time;
	char msg[128] = "";
	trace(3, "rtkpos  : time=%s n=%d\n", time_str(obs[0].time, 3), n);

	time = rtk->sol.time;
	/* rover position by single point positioning */
	if (!pntpos(obs, n, nav, &rtk->opt, &rtk->sol, NULL, rtk->ssat, msg)) {				
		showmsg(msg);
		return 0;
	}
	if (time.time != 0) rtk->tt = timediff(rtk->sol.time, time);

	/* precise point positioning */
	if (opt->mode >= PMODE_PPP_KINEMA) {	
		pppos(rtk, obs, n, nav);			
		return 1;
	}
}