/****************************************
GPS双频非差非组合||消电离层PPP（固定）
*****************************************/
#include "rtklib.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define MAX_ITER    8               /* max number of iterations */
#define MAX_STD_FIX 0.15            /* max std-dev (3d) to fix solution */
#define MIN_NSAT_SOL 4              /* min satellite number for solution */
#define THRES_REJECT 4.0            /* reject threshold of posfit-res (sigma) */

#define THRES_MW_JUMP 10.0

#define VAR_POS     SQR(10.0)       /* init variance receiver position (m^2) */			//60--10
#define VAR_VEL     SQR(10.0)       /* init variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0)       /* init variance of receiver acc ((m/ss)^2) */
#define VAR_CLK     SQR(60.0)       /* init variance receiver clock (m^2) */   //60---100
#define VAR_ZTD     SQR( 0.6)       /* init variance ztd (m^2) */
#define VAR_ZWD     SQR(0.15)		//0.15
#define VAR_GRA     SQR(0.01)       /* init variance gradient (m^2) */
#define VAR_DCB     SQR(30.0)       /* init variance dcb (m^2) */
#define VAR_IFB		SQR(30.0)			/* init variance ifb (m^2) */
#define VAR_BIAS    SQR(100.0)      /* init variance phase-bias (m^2) */		//60---100
#define VAR_BIAS_W  SQR(2E7)		/* init variance phase-bias (m^2) */
#define VAR_IONO    SQR(30.0)       /* init variance iono-delay */
#define VAR_GLO_IFB SQR( 0.6)       /* variance of glonass ifb */

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.5             /* relative humidity for saastamoinen model */
#define GAP_RESION  120             /* default gap to reset ionos parameters (ep) */

#define EFACT_GPS_L5 10.0           /* error factor of GPS/QZS L5 */

#define MUDOT_GPS   (0.00836*D2R)   /* average angular velocity GPS (rad/s) */
#define MUDOT_GLO   (0.00888*D2R)   /* average angular velocity GLO (rad/s) */
#define EPS0_GPS    (13.5*D2R)      /* max shadow crossing angle GPS (rad) */
#define EPS0_GLO    (14.2*D2R)      /* max shadow crossing angle GLO (rad) */
#define T_POSTSHADOW 1800.0         /* post-shadow recovery time (s) */
#define QZS_EC_BETA 20.0            /* max beta angle for qzss Ec (deg) */

#define MIN_ARC_GAP     300.0       /* min arc gap (s) */
/* number and index of states */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)				//频率
#define NP(opt)     ((opt)->dynamics?9:3)									//位置
#define NC(opt)     (NSYS)													//系统
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))//对流层
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)					//电离层（每颗卫星一量）
#define ND(opt)		(((opt)->ionoopt==IONOOPT_EST&&(opt)->posopt[5]!=0)?NSYS:0)	//接收机DCB
#define NS(opt)		((opt)->nf>2?MAXSAT:0)									//IFB
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt)+NS(opt))		//P+C+trop+Iono+rDCB+IFB
#define NB(opt)     (NF(opt)*MAXSAT)										//模糊度
#define NX(opt)     (NR(opt)+NB(opt))										//总状态
#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)
#define ID(s,opt)	(NP(opt)+NC(opt)+NT(opt)+NI(opt)+s)
#define IS(s,opt)	(NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt)+(s)-1)
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)


/* standard deviation of state -----------------------------------------------*/
static double STD(rtk_t *rtk, int i)
{
	if (rtk->sol.stat == SOLQ_FIX) return SQRT(rtk->Pa[i + i*rtk->nx]);
	return SQRT(rtk->P[i + i*rtk->nx]);
}
/* global variables ----------------------------------------------------------*/
static int statlevel = 0;          /* rtk status output level (0:off) */		
static FILE *fp_stat = NULL;       /* rtk status file pointer */

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
	int i;

	trace(3, "rtkinit :\n");

	rtk->sol = sol0;
	for (i = 0; i<6; i++) rtk->rb[i] = 0.0;
	rtk->nx = pppnx(opt);
	rtk->na = pppnx(opt);
	rtk->tt = 0.0;
	rtk->x = zeros(rtk->nx, 1);
	rtk->P = zeros(rtk->nx, rtk->nx);
	rtk->xa = zeros(rtk->na, 1);
	rtk->Pa = zeros(rtk->na, rtk->na);
	//rtk->nfix = rtk->neb = 0;
	for (i = 0; i<MAXSAT; i++) {
		rtk->ambc[i] = ambc0;
		rtk->ssat[i] = ssat0;
	}
	/*rtk->holdamb = 0;
	rtk->excsat = 0;
	for (i = 0; i<MAXERRMSG; i++) rtk->errbuf[i] = 0;*/
	rtk->opt = *opt;
	/*rtk->initial_mode = rtk->opt.mode;
	rtk->com_bias = 0;
	rtk->sol.thres = (float)opt->thresar[0];*/
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
/* open solution status file--------------------------------------------------*/
extern int pppopenstat(const char *file, int level)
{
	if (!(fp_stat = fopen(file, "w"))) {
		trace(1, "pppopenstat: file open error path=%s\n", file);
		return 0;
	}
	statlevel = level;
	return 1;
}
/* close solution status file --------------------------------------------------*/
extern void pppclosestat(void)
{
	trace(3, "rtkclosestat:\n");

	if (fp_stat) fclose(fp_stat);
	fp_stat = NULL;
	statlevel = 0;
}
/* output solution status ----------------------------------------------------*/
static void outsolstat(rtk_t *rtk, obsd_t *obs, int *nu, const nav_t *nav)
{
	ssat_t *ssat;
	ambc_t *amb;
	double tow, pos[3], freq[NFREQ] = { 0.0 };
	char buff[MAXSOLMSG + 1], id[32];
	int i, j, k, m, n, week, nf = NF(&rtk->opt), sys;
	double *x;

	if (statlevel <= 0 || !fp_stat || !rtk->sol.stat || rtk->sol.stat == SOLQ_SINGLE) return;

	tow = time2gpst(rtk->sol.time, &week);
	x = rtk->sol.stat == SOLQ_FIX ? rtk->xa : rtk->x;

	/* fcb */    
	for (i = 0; i < nu; i++) {
		amb = rtk->ambc + obs[i].sat - 1;
		ssat = rtk->ssat + obs[i].sat - 1;
		n = amb->n[1];
		satno2id(obs[i].sat, id);
		sys = satsys(obs[i].sat, NULL);
		if (!ssat->vs || sys != SYS_GPS) continue;

		/*k = IB(obs[i].sat, 0, &rtk->opt);
		fprintf(fp_stat, "%s %s %4d %10.4f %4d %10.4f %10.4f\n", time_str(rtk->sol.time, 0), id, ssat->slip[0], ssat->mw, ssat->vsat[0], x[k], rtk->P[k + k*rtk->nx]);
		continue;*/

		/* Ione */
		if (!ssat->vsat[0] || !ssat->vsat[1]) continue;
		k = II(obs[i].sat, &rtk->opt);
		ecef2pos(rtk->sol.rr, pos);
		for (j = 0; j < NFREQ; j++) freq[j] = sat2freq(obs[i].sat, obs[i].code[j], nav);
		double gamma = SQR(freq[0]) / SQR(freq[1]);
		double Ion = (CLIGHT / freq[0] * obs[i].L[0] - CLIGHT / freq[1] * obs[i].L[1] - amb->LC[2])/(gamma-1);			
		fprintf(fp_stat, "%s %s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n", time_str(rtk->sol.time, 0), id, pos[0] * R2D, pos[1] * R2D,
			ssat->azel[0] * R2D, ssat->azel[1] * R2D, x[k], Ion);
	}

	return;

	/* receiver position */
	fprintf(fp_stat, "POS\t%.3f\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t", tow,
		rtk->sol.stat, x[0], x[1], x[2], STD(rtk, 0), STD(rtk, 1), STD(rtk, 2));

	/* receiver clocks */
	i = IC(0, &rtk->opt);
	fprintf(fp_stat, "C\t%.3f\t%.3f\t", x[i] * 1E9 / CLIGHT, STD(rtk, i)*1E9 / CLIGHT);

	/* tropospheric parameters */
	if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {
		i = IT(&rtk->opt);
		fprintf(fp_stat, "T\t%.4f\t%.4f\t", x[i], STD(rtk, i));
	}
	if (rtk->opt.tropopt == TROPOPT_ESTG) {
		i = IT(&rtk->opt);
		fprintf(fp_stat, "%.5f\t%.5f\t%.5f\t%.5f", x[i + 1], x[i + 2], STD(rtk, i + 1), STD(rtk, i + 2));
	}
	fprintf(fp_stat, "\n");
	/* ionosphere parameters */
	if (rtk->opt.ionoopt == IONOOPT_EST) {					//伪距计算

		for (i = 0; i<MAXSAT; i++) {
			ssat = rtk->ssat + i;
			if (!ssat->vs) continue;
			j = II(i + 1, &rtk->opt);
			if (rtk->x[j] == 0.0) continue;
			satno2id(i + 1, id);
			fprintf(fp_stat, "ION\t%.3f\t%s\t%.4f\t%.4f\n", tow, id, x[j], STD(rtk, j));
		}
	}
	/* ambiguity parameters */
	for (i = 0; i<MAXSAT; i++) {							//伪距-相位（顾及电离层）
		int stat = 0;
		for (j = 0; j<nf; j++){
			int k = IB(i + 1, j, &rtk->opt);
			if (rtk->x[k] == 0.0) continue;
			else stat = 1;
		}
		if (stat == 0) continue;
		satno2id(i + 1, id);
		fprintf(fp_stat, "AMB\t%.3f\t%s\t", tow, id);
		for (j = 0; j < nf; j++){
			int k = IB(i + 1, j, &rtk->opt);
			fprintf(fp_stat, "%.4f\t%.8f\t", x[k], STD(rtk, k));
		}
		fprintf(fp_stat, "\n");
	}
	/* write residuals and status */
	for (i = 0; i<MAXSAT; i++) {
		ssat = rtk->ssat + i;
		if (!ssat->vs) continue;
		satno2id(i + 1, id);
		fprintf(fp_stat, "SAT\t%.3f\t%s\t%.1f\t%.1f\t", tow, id, rtk->ssat[i].azel[0] * R2D, rtk->ssat[i].azel[1] * R2D);
		for (j = 0; j<nf; j++) {
			k = IB(i + 1, j, &rtk->opt);
			//fprintf(fp_stat, "%.1f\t%.1f\t%d\t%d\t%d\t",ssat->resp[j], ssat->resc[j], ssat->slip[j] & 3, ssat->lock[j], ssat->outc[j]);
			fprintf(fp_stat, "%d\t%d\t%d\t", ssat->slip[j] & 3, ssat->lock[j], ssat->outc[j]);
		}
		fprintf(fp_stat, "\n");
		//ssat->azel[0] * R2D, ssat->azel[1] * R2D, ssat->vsat[j], ssat->snr[j] * 0.25, ssat->slipc[j], ssat->rejc[j]
	}
}
/* exclude meas of eclipsing satellite (block IIA)不包括日食卫星 --------------------*/
static void testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs)
{
	double rsun[3], esun[3], r, ang, erpv[5] = { 0 }, cosa;
	int i, j, sat;
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
		//sat = obs[i].sat;

		//if (cosa<0 && r*sin(ang)<RE_WGS84){
		//	memcpy(&rtk->ssat[sat - 1].lastEclipse, &obs[0].time, sizeof(gtime_t));
		//	satno2id(sat, id);
		//	trace(2, "eclipsing sat excluded %s sat=%s\n", time_str(obs[0].time, 0), id);
		//	for (j = 0; j<3; j++) rs[j + i * 6] = 0.0;
		//}
		//else if (rtk->ssat[sat - 1].lastEclipse.time != 0){
		//	diff = fabs(timediff(obs[0].time, rtk->ssat[sat - 1].lastEclipse));
		//	if (diff <= 1800) { // 30 minutes
		//		satno2id(sat, id);
		//		trace(2, "under eclipse %4.0lf seconds ago %s sat=%s\n", diff, time_str(obs[0].time, 0), id);
		//		for (j = 0; j<3; j++) rs[j + i * 6] = 0.0;
		//	}
		//}

		/* test eclipse */
		if (ang<PI / 2.0 || r*sin(ang)>RE_WGS84) continue;			//判断是否日食（将卫星位置置0）

		trace(3, "eclipsing sat excluded %s sat=%2d\n", time_str(obs[0].time, 0),
			obs[i].sat);

		for (j = 0; j<3; j++) rs[j + i * 6] = 0.0;
	}
}
/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
	if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
	return atan2(-tan(beta), sin(mu)) + PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
extern int yaw_angle(int sat, const char *type, int opt, double beta, double mu,
	double *yaw)
{
	*yaw = yaw_nominal(beta, mu);
	return 1;
}
/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char *type, int opt,
	const double *rs, double *exs, double *eys)
{
	double rsun[3], ri[6], es[3], esun[3], n[3], p[3], en[3], ep[3], ex[3], E, beta, mu;
	double yaw, cosy, siny, erpv[5] = { 0 };
	int i;

	sunmoonpos(gpst2utc(time), erpv, rsun, NULL, NULL);

	/* beta and orbit angle */
	matcpy(ri, rs, 6, 1);
	ri[3] -= OMGE*ri[1];
	ri[4] += OMGE*ri[0];
	cross3(ri, ri + 3, n);
	cross3(rsun, n, p);
	if (!normv3(rs, es) || !normv3(rsun, esun) || !normv3(n, en) ||
		!normv3(p, ep)) return 0;
	beta = PI / 2.0 - acos(dot(esun, en, 3));
	E = acos(dot(es, ep, 3));
	mu = PI / 2.0 + (dot(es, esun, 3) <= 0 ? -E : E);
	if (mu<-PI / 2.0) mu += 2.0*PI;
	else if (mu >= PI / 2.0) mu -= 2.0*PI;

	/* yaw-angle of satellite */
	if (!yaw_angle(sat, type, opt, beta, mu, &yaw)) return 0;

	/* satellite fixed x,y-vector */
	cross3(en, es, ex);
	cosy = cos(yaw);
	siny = sin(yaw);
	for (i = 0; i<3; i++) {
		exs[i] = -siny*en[i] + cosy*ex[i];
		eys[i] = -cosy*en[i] - siny*ex[i];
	}
	return 1;
}
/* phase windup model --------------------------------------------------------*/
static int model_phw(gtime_t time, int sat, const char *type, int opt,
	const double *rs, const double *rr, double *phw)
{
	double exs[3], eys[3], ek[3], exr[3], eyr[3], eks[3], ekr[3], E[9];
	double dr[3], ds[3], drs[3], r[3], pos[3], cosp, ph;
	int i;

	if (opt <= 0) return 1; /* no phase windup */

	/* satellite yaw attitude model */
	if (!sat_yaw(time, sat, type, opt, rs, exs, eys)) return 0;

	/* unit vector satellite to receiver */
	for (i = 0; i<3; i++) r[i] = rr[i] - rs[i];
	if (!normv3(r, ek)) return 0;

	/* unit vectors of receiver antenna */
	ecef2pos(rr, pos);
	xyz2enu(pos, E);
	exr[0] = E[1]; exr[1] = E[4]; exr[2] = E[7]; /* x = north */
	eyr[0] = -E[0]; eyr[1] = -E[3]; eyr[2] = -E[6]; /* y = west  */

	/* phase windup effect */
	cross3(ek, eys, eks);
	cross3(ek, eyr, ekr);
	for (i = 0; i<3; i++) {
		ds[i] = exs[i] - ek[i] * dot(ek, exs, 3) - eks[i];			//卫星天线有效偶极向量
		dr[i] = exr[i] - ek[i] * dot(ek, exr, 3) + ekr[i];			//接收机天线有效偶极向量
	}
	cosp = dot(ds, dr, 3) / norm(ds, 3) / norm(dr, 3);
	if (cosp<-1.0) cosp = -1.0;
	else if (cosp> 1.0) cosp = 1.0;
	ph = acos(cosp) / 2.0 / PI;
	cross3(ds, dr, drs);
	if (dot(ek, drs, 3)<0.0) ph = -ph;

	*phw = ph + floor(*phw - ph + 0.5); /* in cycle */		//floor(*phw - ph + 0.5)
	return 1;
}

/* measurement error variance ------------------------------------------------*/
static double varerr(int sat, double el, int freq, int type, const prcopt_t *opt)
{
	int sys, prn=0;
	double fact = 1.0, sinel = sin(el), varr;

	sys = satsys(sat,&prn);
	varr = SQR(opt->err[1]) + SQR(opt->err[2] / sinel);

	if (type == 1) fact *= opt->eratio[freq == 0 ? 0 : 1];					
	fact *= sys == SYS_GLO ? EFACT_GLO : (sys == SYS_SBS ? EFACT_SBS : EFACT_GPS);				

	if (sys == SYS_GPS || sys == SYS_QZS) {
		if (freq == 2) fact *= EFACT_GPS_L5; /* GPS/QZS L5 error factor */
	}
	if (sys == SYS_CMP && prn < 19) {
		if (type == 1) fact *= 2;
	}
	if (opt->ionoopt == IONOOPT_IFLC) fact *= 3.0;
	return SQR(fact)*varr;
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
/* geometry-free phase measurement -------------------------------------------*/
static double gfmeas(const obsd_t *obs, const nav_t *nav, int type)			
{
	int i, j;
	double freq1, freq2;
	if (!type) { i = 0;	j = 1;}
	else { i = 1; j = 2;}

	freq1 = sat2freq(obs->sat, obs->code[i], nav);
	freq2 = sat2freq(obs->sat, obs->code[j], nav);

	if (freq1 == 0.0 || freq2 == 0.0 || obs->L[i] == 0.0 || obs->L[j] == 0.0) return 0.0;
	return (obs->L[i] / freq1 - obs->L[j] / freq2)*CLIGHT;
}
/* Melbourne-Wubbena linear combination --------------------------------------*/
static double mwmeas(const obsd_t *obs, const nav_t *nav, double *azel, int type)										
{
	int i, j, prn, sys;
	double freq1, freq2;
	double temp = 0, P1, P2, P1_C1, P2_C2, lam1, lam2, lam4, lam5;

	if (!type) { i = 0; j = 1; }
	else { i = 1; j = 2; }
	freq1 = sat2freq(obs->sat, obs->code[i], nav);
	freq2 = sat2freq(obs->sat, obs->code[j], nav);
	sys = satsys(obs->sat, &prn);
	if (freq1 == 0.0 || freq2 == 0.0 || obs->L[i] == 0.0 || obs->L[j] == 0.0 ||
		obs->P[i] == 0.0 || obs->P[j] == 0.0) return 0.0;

	P1 = obs->P[i];
	P2 = obs->P[j];
	P1_C1 = nav->cbias[obs->sat - 1][1];
	P2_C2 = nav->cbias[obs->sat - 1][2];
	/* code with dcb correction */					
	if (sys == SYS_GPS || sys == SYS_GLO) {
		if (obs->code[0] == CODE_L1C) P1 += P1_C1; /* C1->P1 */				//2021-06-09 + -
		if (obs->code[1] == CODE_L2C) P2 += P2_C2; /* C2->P2 */
	}
	else if (sys == SYS_CMP) {
		P1 += bd2smp(obs, prn, azel, i, NULL);
		P2 += bd2smp(obs, prn, azel, j, NULL);
	}

	lam1 = CLIGHT / freq1; lam2 = CLIGHT / freq2;

	lam4 = 1 / (1 / lam1 - 1 / lam2); //宽巷组合观测值波长
	lam5 = 1 / (1 / lam1 + 1 / lam2); //窄巷组合观测值波长
	return (obs->L[i] - obs->L[j]) - (lam5*(P1 / lam1 + P2 / lam2) / lam4);
}
/* geometry-free phase measurement -------------------------------------------*/
static double gfPLmeas(const obsd_t *obs, const nav_t *nav, double *azel)
{
	double freq1, freq2;
	double P1, P2, P1_C1, P2_C2;
	int i = (satsys(obs->sat, NULL)&(SYS_SBS)) ? 2 : 1, prn;
	int sys;

	freq1 = sat2freq(obs->sat, obs->code[0], nav);
	freq2 = sat2freq(obs->sat, obs->code[i], nav);
	sys = satsys(obs->sat, &prn);
	if (freq1 == 0.0 || freq2 == 0.0 || obs->L[0] == 0.0 || obs->L[i] == 0.0
		|| obs->P[0] == 0.0 || obs->P[i] == 0.0) return 0.0;
	P1 = obs->P[0];
	P2 = obs->P[i];
	P1_C1 = nav->cbias[obs->sat - 1][1];
	P2_C2 = nav->cbias[obs->sat - 1][2];
	/* code with dcb correction */					
	if (sys == SYS_GPS || sys == SYS_GLO){									//卫星端多路径需要改正
		if (obs->code[0] == CODE_L1C) P1 += P1_C1; /* C1->P1 */				//2021-06-09 + -
		if (obs->code[1] == CODE_L2C) P2 += P2_C2; /* C2->P2 */
	}
	else if (sys == SYS_CMP) {
		P1 += bd2smp(obs, prn, azel, 0, NULL);
		P2 += bd2smp(obs, prn, azel, i, NULL);
	}
	return (obs->L[0] / freq1 - obs->L[1] / freq2)*CLIGHT + P1 - P2;		//相位平滑伪距计算的电离层吗？
}
/* antenna corrected measurements --------------------------------------------*/
static void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,
	const prcopt_t *opt, const double *dantr, const double *dants, double phw,
	double *L, double *P, double *Lc, double *Pc)
{
	double freq[NFREQ], c1, c2, C1, C2, gamma;
	int i, k, sys, prn, ix;

	for (i = 0; i < NFREQ; i++) freq[i] = sat2freq(obs->sat, obs->code[i], nav);
	sys = satsys(obs->sat, &prn);

	for (i = 0; i<NFREQ; i++) {					
		L[i] = P[i] = 0.0;
		if (freq[i] == 0.0 || obs->L[i] == 0.0 || obs->P[i] == 0.0) continue;
		if (testsnr(0, i, azel[1], obs->SNR[i] * 0.25, &opt->snrmask)) continue;

		/* antenna phase center and phase windup correction */
	
		L[i] = CLIGHT * obs->L[i] / freq[i] - dants[i] - dantr[i] - CLIGHT * phw / freq[i];		
		P[i] = obs->P[i] - dants[i] - dantr[i];

		if (opt->sateph == EPHOPT_SSRAPC || opt->sateph == EPHOPT_SSRCOM) {
			/* use SSR code correction */
			if (sys == SYS_GPS) {
				ix = (i == 0 ? CODE_L1W - 1 : CODE_L2W - 1);
			}
			else if (sys == SYS_GLO) {
				ix = (i == 0 ? CODE_L1P - 1 : CODE_L2P - 1);
			}
			P[i] += (nav->ssr[obs->sat - 1].cbias[obs->code[i] - 1] - nav->ssr[obs->sat - 1].cbias[ix]);
		}
		else {
			/* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
			if (sys == SYS_GPS || sys == SYS_GLO) {				
				if (obs->code[i] == CODE_L1C) P[i] += nav->cbias[obs->sat - 1][1];				
				if (obs->code[i] == CODE_L2C) P[i] += nav->cbias[obs->sat - 1][2];
			}
			else if (sys == SYS_CMP) {
				P[i] += bd2smp(obs, prn, azel, i, NULL);
			}
			else if (sys == SYS_GAL) {
				if (obs->code[i] == CODE_L1X) P[i] += nav->cbias[obs->sat - 1][1];
				if (obs->code[i] == CODE_L5X) P[i] += nav->cbias[obs->sat - 1][2];
			}
		}
		if (opt->ionoopt == IONOOPT_EST && opt->posopt[5]) {
			gamma = SQR(freq[0]) / SQR(freq[1]);			
			c1 = -1.0 / (gamma - 1.0);				c2 = -gamma / (gamma - 1.0);
			if (i == 0) P[i] -= c1*nav->cbias[obs->sat - 1][0];						//2021-06-09 - +
			else if (i == 1) P[i] -= c2*nav->cbias[obs->sat - 1][0];
		}
	}
	/* iono-free LC */
	*Lc = *Pc = 0.0;
	if (freq[0] == 0.0 || freq[1] == 0.0) return;

	C1 =  SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));			//SQR(lam[i]) / (SQR(lam[i]) - SQR(lam[0]));
	C2 = -SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));			//-SQR(lam[0]) / (SQR(lam[i]) - SQR(lam[0]));

#if 0
	/* P1-P2 dcb correction (P1->Pc,P2->Pc) */
	if (sys&(SYS_GPS | SYS_GLO | SYS_QZS)) {
		if (P[0] != 0.0) P[0] -= C2*nav->cbias[obs->sat - 1][0];
		if (P[1] != 0.0) P[1] += C1*nav->cbias[obs->sat - 1][0];
	}
#endif
	if (L[0] != 0.0&&L[1] != 0.0) *Lc = C1*L[0] + C2*L[1];
	if (P[0] != 0.0&&P[1] != 0.0) *Pc = C1*P[0] + C2*P[1];
}
/* LC observation */
static void lc_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav) 
{
	int i, sat,sys;
	double el;
	for (i = 0; i < n && i < MAXOBS; i++) {
		sat = obs[i].sat;
		el = rtk->ssat[sat - 1].azel[1] / PI * 180;
		if (!rtk->ssat[sat - 1].vs) continue;
		if (!(sys = satsys(sat, NULL))) continue;

		rtk->ssat[sat - 1].gf = gfmeas(obs + i, nav, 0);
		rtk->ssat[sat - 1].mw = mwmeas(obs + i, nav, rtk->ssat[sat - 1].azel, 0);
		rtk->ssat[sat - 1].gfPL = gfPLmeas(obs + i, nav, rtk->ssat[sat - 1].azel);
		rtk->ssat[sat - 1].gf2 = gfmeas(obs + i, nav, 1);
	}
}
/* temporal update of position -----------------------------------------------*/
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
	/* initialize position for first epoch */									//首历元初始化
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
	double dtr;
	int i;

	trace(3, "udclk_ppp:\n");

	/* initialize every epoch for clock (white noise) */
	for (i = 0; i<NSYS; i++) {										
		dtr = i == 0 ? rtk->sol.dtr[0] : rtk->sol.dtr[0] + rtk->sol.dtr[i];			//spp钟差 dtr[i]:ISB
		initx(rtk, CLIGHT*dtr, VAR_CLK, IC(i, &rtk->opt));
	}
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop_ppp(rtk_t *rtk)
{
	double pos[3], azel[] = { 0.0, PI / 2.0 }, ztd, var;
	int i = IT(&rtk->opt), j;

	trace(3, "udtrop_ppp:\n");

	if (rtk->x[i] == 0.0) {
		ecef2pos(rtk->sol.rr, pos);
		ztd = tropmodel1(rtk->sol.time, pos, azel, 1, REL_HUMI, 2);					//湿延迟初值,约0.1m
		initx(rtk, ztd, VAR_ZWD, i);

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
	double freq[NFREQ], ion, sinel, pos[3], dantr[3] = { 0 }, dants[3] = { 0 }, *azel, vion, L[NFREQ], P[NFREQ], Lc, Pc;
	char *p;
	int i, j, k, sat, gap_resion = GAP_RESION;

	trace(3, "udiono_ppp:\n");

	if ((p = strstr(rtk->opt.pppopt, "-GAP_RESION="))) {
		sscanf(p, "-GAP_RESION=%d", &gap_resion);
	}
	for (i = 0; i<MAXSAT; i++) {
		j = II(i + 1, &rtk->opt);
		if (rtk->x[j] != 0.0 && (int)rtk->ssat[i].outc[0]>gap_resion) {				
			initx(rtk, 0.0, 0.0, j);
		}
	}
	ecef2pos(rtk->sol.rr, pos);
	for (i = 0; i<n; i++) {
		sat = obs[i].sat;
		if (!rtk->ssat[sat - 1].vs) continue;
		j = II(sat, &rtk->opt);
		azel = rtk->ssat[sat - 1].azel;
		for (k = 0; k < NFREQ; k++) freq[k] = sat2freq(sat, obs[i].code[k], nav);				
		if (rtk->x[j] == 0.0) {
			if (rtk->opt.posopt[5]) {												//外部电离层约束
				if (!iontec(obs[i].time, nav, pos, azel, 3, &ion, &vion)) continue;
				ion *= SQR(FREQ1 / freq[0]);
				//vion *= SQR(FREQ1 / freq[0]) * SQR(FREQ1 / freq[0]);
				trace(2, "UPPP ION:sat%3d %10.4f %10.4f\n", sat, ion, SQR(FREQ1 / freq[0]));
			}
			else {
				k = 1;
				corr_meas(obs + i, nav, azel, &rtk->opt, dantr, dants, 0.0, L, P, &Lc, &Pc);					
				if (obs[i].P[0] != 0.0 && obs[i].P[k] != 0.0 && freq[0] != 0.0 && freq[k] != 0.0) {
					ion = (obs[i].P[0] - obs[i].P[k]) / (1.0 - SQR(freq[0] / freq[k]));
				}
				else {
					ion = ionmodel(obs[i].time, nav->ion_gps, pos, azel);			
					ion *= SQR(FREQ1 / freq[0]);
				}
			}
			if (ion != 0.0) initx(rtk, ion, VAR_IONO, j);
		}
		else {
			sinel = sin(MAX(azel[1], 5.0*D2R));
			rtk->P[j + j*rtk->nx] += SQR(rtk->opt.prn[1] / sinel)*fabs(rtk->tt);	//电离层过程噪声（斜路径）
		}
	}
}
/* temporal update of receiver dcb parameters --------------------------------*/
static void uddcb_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	int i, j;
	double rdcb, C;
	trace(3, "uddcb_ppp:\n");

	C = 1.0 / (1.0 - SQR(lam_carr[1] / lam_carr[0]));
	for (i = 0; i<NSYS; i++) {
		j = ID(i, &rtk->opt);
		if (rtk->opt.posopt[4]) {						
			/*rdcb = -11.461;
			initx(rtk, C*rdcb*1E-9*CLIGHT, 1e-8, j);*/
			continue;
		}
		else {
			if (rtk->x[j] == 0.0) initx(rtk, 1e-6, VAR_DCB, j);    //1e-6
			else rtk->P[j + j*rtk->nx] += SQR(rtk->opt.prn[5])*fabs(rtk->tt)/*/3600*/;
		}
	}
}
static void udifb_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav) 
{
	int i, j, sat;
	for (i = 0; i < n; i++) {
		sat = obs[i].sat;
		if (!rtk->ssat[sat - 1].vs) continue;
		j = IS(sat, &rtk->opt);
		if (rtk->x[j] == 0.0) initx(rtk, 1e-6, VAR_IFB, j);    //1e-6
	}
}
/* temporal update of phase biases ------------------------------------------- */
static void udbias_ppp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)							
{
	gtime_t t0;
	double L[NFREQ], P[NFREQ], Lc, Pc, freq[NFREQ], bias[MAXOBS], offset = 0.0, gamma;
	double ion, dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 }, gfPL, gf0, el, fact, thresgf, thresmw, diffamb, sig0;
	int i, j, k, l, f, sat, sys, clk_jump = 0, flag, slip, num;

	trace(3, "udbias  : n=%d\n", n);

	/* handle day-boundary clock jump */
	if (rtk->opt.posopt[5]) {
		clk_jump = ROUND(time2gpst(obs[0].time, NULL) * 10) % 864000 == 0;
	}
	for (i = 0; i<MAXSAT; i++) for (j = 0; j<rtk->opt.nf; j++) {
		rtk->ssat[i].slip[j] = 0;
	}

	/* LC组合观测值 */
	for (i = 0; i < n && i < MAXOBS; i++) {					
		sat = obs[i].sat;
		el = rtk->ssat[sat - 1].azel[1] / PI * 180;
		if (!rtk->ssat[sat - 1].vs) continue;	
		if (!(sys = satsys(sat, NULL))) continue;
		
		if (rtk->ssat[sat - 1].gf == 0.0 || rtk->ssat[sat - 1].mw == 0.0
			|| (sys&SYS_GAL && rtk->opt.nf>2 && rtk->ssat[sat - 1].gf2 == 0.0)) {
			for (j = 0; j < rtk->opt.nf; j++) rtk->ssat[sat - 1].slip[j] = 1;			//周跳标记
			for (j = 0; j < 4; j++) {
				rtk->ambc[sat - 1].n[j] = 0;
				rtk->ambc[sat - 1].LC[j] = rtk->ambc[sat - 1].LCv[j] = 0.0;
			}	
		}
		else {
			slip = 0;
			for (j = 0; j < 4; j++) {
				t0 = rtk->ambc[sat - 1].epoch[j];
				if (fabs(timediff(rtk->sol.time, t0)) + DTTOL>300) {				
					rtk->ambc[sat - 1].n[j] = 0;
					rtk->ambc[sat - 1].LC[j] = rtk->ambc[sat - 1].LCv[j] = 0.0;
				}
			}
			/*if (rtk->ssat[sat - 1].gf != 0.0) rtk->ambc[sat - 1].epoch[0] = rtk->sol.time;
			if (rtk->ssat[sat - 1].mw != 0.0) rtk->ambc[sat - 1].epoch[1] = rtk->sol.time;
			if (rtk->ssat[sat - 1].gfPL != 0.0) rtk->ambc[sat - 1].epoch[2] = rtk->sol.time;
			if (rtk->ssat[sat - 1].gf2 != 0.0) rtk->ambc[sat - 1].epoch[3] = rtk->sol.time;*/

			// LLI周跳探测
			for (j = 0; j < rtk->opt.nf; j++) {
				flag = (obs[i].L[j] == 0.0 || (obs[i].LLI[j] & 3)) ? 1 : 0;
				if (flag) trace(2, "%s sat=%2d LL cycle slip f=%d\n", time_str(obs[i].time, 0), sat, j + 1);
				if (j < 2) slip |= flag;
			}

			// gf周跳探测
			if (el >= 15) fact = 1.0;
			else		  fact = (7 - 0.4*el);
			thresgf = fact*rtk->opt.thresslipgf;
			diffamb = rtk->ssat[sat - 1].gf - rtk->ambc[sat - 1].LC[0];					
			flag = (rtk->ambc[sat - 1].LC[0] == 0.0 || (rtk->ambc[sat - 1].LC[0] != 0.0 && fabs(diffamb)>thresgf)) ? 1 : 0;
			slip |= flag;
			if (flag) trace(2, "detslip_gf: slip detected time=%s sat=%2d gf=%8.3f->%8.3f\n", time_str(obs[i].time, 0), obs[i].sat, rtk->ambc[sat - 1].LC[0], rtk->ssat[sat - 1].gf);

			// gf2周跳探测
			if (sys&SYS_GAL && rtk->opt.nf > 2) {
				if (el >= 15) fact = 1.0;
				else		  fact = (7 - 0.4*el);
				thresgf = fact*rtk->opt.thresslipgf;
				diffamb = rtk->ssat[sat - 1].gf2 - rtk->ambc[sat - 1].LC[3];
				flag = (rtk->ambc[sat - 1].LC[3] == 0.0 || (rtk->ambc[sat - 1].LC[3] != 0.0 && fabs(diffamb)>thresgf)) ? 1 : 0;
				slip |= flag;
				if (flag) trace(2, "detslip_gf: slip detected time=%s sat=%2d gf=%8.3f->%8.3f\n", time_str(obs[i].time, 0), obs[i].sat, rtk->ambc[sat - 1].LC[3], rtk->ssat[sat - 1].gf2);
			}

			// mw周跳探测
			if (el >= 15) fact = 1;
			else		  fact = (4 - 0.2*el);				
			thresmw = fact*rtk->opt.thresslipmw;
			diffamb = rtk->ssat[sat - 1].mw - rtk->ambc[sat - 1].LC[1];
			flag = (rtk->ambc[sat - 1].n[1] == 0 ||(rtk->ambc[sat - 1].n[1] >0 && fabs(diffamb) > thresmw)) ? 1 : 0;
			slip |= flag;
			if (flag) trace(2, "detslip_mw: slip detected time=%s sat=%2d num=%4d mw=%8.3f->%8.3f\n", time_str(obs[i].time, 0), obs[i].sat, rtk->ambc[sat - 1].n[1], 
				rtk->ambc[sat - 1].LC[1], rtk->ssat[sat - 1].mw);
			
			// 周跳标记
			for (j = 0; j < rtk->opt.nf; j++) rtk->ssat[sat - 1].slip[j] = slip;

			// ambc更新
			/*for (j = 0; j < 4; j++) {
				double LC = j == 0 ? rtk->ssat[sat - 1].gf : j == 1 ? rtk->ssat[sat - 1].mw : j == 2 ? rtk->ssat[sat - 1].gfPL : rtk->ssat[sat - 1].gf2;
				if (LC) {
					rtk->ambc[sat - 1].epoch[j] = rtk->sol.time;
				}
				if (slip || j == 0 || j == 3 || LC == 0.0) {
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
			}*/
		}
	}

	for (f = 0; f<NF(&rtk->opt); f++) {												
		/* reset phase-bias if expire obs outage counter */
		for (i = 0; i<MAXSAT; i++) {
			if (++rtk->ssat[i].outc[f]>(unsigned int)rtk->opt.maxout && rtk->x[IB(i + 1, f, &rtk->opt)] != 0.0) {
				initx(rtk, 0.0, 0.0, IB(i + 1, f, &rtk->opt));
			}
		}
		for (i = 0; i<n&&i<MAXOBS; i++) {
			sat = obs[i].sat;
			j = IB(sat, f, &rtk->opt);

			corr_meas(obs + i, nav, rtk->ssat[sat - 1].azel, &rtk->opt, dantr, dants, 0.0, L, P, &Lc, &Pc);				

			bias[i] = 0.0;
			rtk->P[j + j*rtk->nx] += SQR(rtk->opt.prn[0])*fabs(rtk->tt);
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

				if (rtk->opt.posopt[5]) {
					bias[i] += rtk->x[ID(0, &rtk->opt)] * gamma;
				}
			}
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
	if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {					//湿延迟
		udtrop_ppp(rtk);
	}
	/* temporal update of ionospheric parameters */
	if (rtk->opt.ionoopt == IONOOPT_EST) {
		udiono_ppp(rtk, obs, n, nav);					
		/* temporal update of receiver DCB */
		if (rtk->opt.posopt[5]) {						
			uddcb_ppp(rtk, obs, n, nav);
		}
	}
	/* temporal update of inter frequency bias */
	if (rtk->opt.nf > 2){
		udifb_ppp(rtk, obs, n, nav);
	}
	/* temporal update of phase-bias */
	udbias_ppp(rtk, obs, n, nav);						
}
/* satellite antenna phase center variation ----------------------------------*/
static void satantpcv(const double *rs, const double *rr, const pcv_t *pcv,
	double *dant)
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

	antmodel_s(pcv, nadir, dant);						//rtkcmn
}
/* precise tropospheric model ------------------------------------------------*/
static double trop_model_prec(gtime_t time, const double *pos,
	const double *azel, const double *x, double *dtdx, double *var)
{
	const double zazel[] = { 0.0, PI / 2.0 };
	double zhd, m_h, m_w, cotz, grad_n, grad_e, grade = 0.0;

	/* zenith hydrostatic delay */
	//zhd = tropmodel(time, pos, zazel, 0.0);
	zhd = tropmodel1(time, pos, zazel, 1, .0, 1);				//这里求得应该是干延迟

	/* mapping function */
	m_h = tropmapf(time, pos, azel, 1, &m_w);					//m_h：静力学投影函数		m_w：动力学投影函数

	if (azel[1]>0.0) {

		/* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
		cotz = 1.0 / tan(azel[1]);
		grad_n = m_w*cotz*cos(azel[0]);
		grad_e = m_w*cotz*sin(azel[0]);
		//m_w += grad_n*x[1] + grad_e*x[2];
		//dtdx[1] = grad_n*(x[0] - zhd);						//x:ZTD   ZTD-ZHD=ZWD（估梯度）
		//dtdx[2] = grad_e*(x[0] - zhd);
		grade = grad_n*x[1] + grad_e*x[2];
		dtdx[1] = grad_n*x[0];										//x:ZTD   ZTD-ZHD=ZWD（估梯度）
		dtdx[2] = grad_e*x[0];
	}
	dtdx[0] = m_w;
	*var = SQR(0.01);
	return  m_h*zhd + (m_w + grade)*x[0];						//m_h*zhd + m_w*(x[0]) + grade;		//m_h*zhd + m_w*(x[0] - zhd);
}
/* tropospheric model ---------------------------------------------------------*/
static int model_trop(gtime_t time, const double *pos, const double *azel,
	const prcopt_t *opt, const double *x, double *dtdx,
	const nav_t *nav, double *dtrp, double *var)						//dtrp:延迟项	||   var:方差项
{
	double trp[3] = { 0 }, std[3];

	if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
		matcpy(trp, x + IT(opt), opt->tropopt == TROPOPT_EST ? 1 : 3, 1);
		*dtrp = trop_model_prec(time, pos, azel, trp, dtdx, var);
		*var = 0.0;
		return 1;
	}
	if (opt->tropopt == TROPOPT_ZTD) {									//对流层ZTD改正
		if (pppcorr_trop(&nav->pppcorr, time, pos, trp, std)) {
			*dtrp = trop_model_prec(time, pos, azel, trp, dtdx, var);
			*var = SQR(dtdx[0] * std[0]);
			return 1;
		}
	}
	return 0;
}
/* ionospheric model ---------------------------------------------------------*/
static int model_iono(gtime_t time, const double *pos, const double *azel,
	const prcopt_t *opt, int sat, const double *x,
	const nav_t *nav, double *dion, double *var)
{
	static double iono_p[MAXSAT] = { 0 }, std_p[MAXSAT] = { 0 };
	static gtime_t time_p;

	if (opt->ionoopt == IONOOPT_TEC) {
		return iontec(time, nav, pos, azel, 1, dion, var);				//斜延迟(L1频段)
	}
	if (opt->ionoopt == IONOOPT_BRDC) {
		*dion = ionmodel(time, nav->ion_gps, pos, azel);
		*var = SQR(*dion*ERR_BRDCI);
		return 1;
	}
	if (opt->ionoopt == IONOOPT_EST) {									//待估参数斜延迟	
		*dion = x[II(sat, opt)];
		*var = 0.0;
		return 1;
	}
	if (opt->ionoopt == IONOOPT_IFLC) {
		*dion = *var = 0.0;
		return 1;
	}
	if (opt->ionoopt == IONOOPT_STEC) {									//斜延迟信息改正(只改正，不估计)
		if (timediff(time, time_p) != 0.0&&
			!pppcorr_stec(&nav->pppcorr, time, pos, iono_p, std_p)) return 0;
		if (iono_p[sat - 1] == 0.0 || std_p[sat - 1]>0.1) return 0;
		time_p = time;
		*dion = iono_p[sat - 1];
		*var = SQR(std_p[sat - 1]);
		return 1;
	}
	return 0;
}
/* constraint to local correction --------------------------------------------*/
static int const_corr(const obsd_t *obs, int n, const nav_t *nav,
	const double *dr, const double *x, rtk_t *rtk, double *v, double *H, double *var)
{
	gtime_t time = obs[0].time;
	double trop[3], std_trop[3], iono, std_iono, rr[3], pos[3], tempv, *azel, freq[NFREQ], fact, maxvc = 4.0;
	int i, j, k, sat, nv = 0, sys;
	for (i = 0; i<3; i++) rr[i] = x[i] + dr[i];
	char id[4];
	ecef2pos(rr, pos);
	/* constraint to external troposphere correction */
	if ((rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) &&
		pppcorr_trop(&nav->pppcorr, time, pos, trop, std_trop)) {

		for (i = 0; i<(rtk->opt.tropopt == TROPOPT_EST ? 1 : 3); i++) {
			if (std_trop[i] == 0.0) continue;
			j = IT(&rtk->opt) + i;
			v[nv] = trop[i] - x[j];
			for (k = 0; k<rtk->nx; k++) H[k + nv*rtk->nx] = k == j ? 1.0 : 0.0;
			var[nv++] = SQR(0.13);
			//var[nv++] = SQR(std_trop[i]);
		}
	}
	/* constraint to external ionosphere correction */
	if (rtk->opt.ionoopt == IONOOPT_EST && rtk->opt.posopt[5]) {					
		for (i = 0; i<n; i++) {
			sat = obs[i].sat;
			azel = rtk->ssat[sat - 1].azel;
			for (j = 0; j < NFREQ; j++) freq[j] = sat2freq(sat, obs[i].code[j], nav);
			if (azel[1] < 15 * D2R) continue;
			if (!iontec(obs[i].time, nav, pos, azel, 3, &iono, &std_iono)) continue;
			fact = SQR(FREQ1/freq[0]);
			iono *= fact;
			std_iono = 0.09 + 0.16 / SQR(sin(azel[1]));			//elevation dependent weighting
			//std_iono = azel[i * 2 + 1]<30 * D2R ? SQR(0.4) / (2 * sin(azel[1])) : SQR(0.4);
			j = II(sat, &rtk->opt);
			tempv = iono - x[j];
			v[nv] = tempv;
			for (k = 0; k<rtk->nx; k++) H[k + nv*rtk->nx] = k == j ? 1.0 : 0.0;
			var[nv++] = fabs(tempv)>maxvc?1.0e8:std_iono;

			/*satno2id(sat, id);
			trace(2, "%s %s %6.2f %10.4f\n", time_str(obs[0].time, 0), id, azel[1] / PI * 180, iono);*/
		}
	}
	return nv;
}
/* phase and code residuals --------------------------------------------------*/
static int ppp_res(const obsd_t *obs, int n, const double *rs, const double *dts,
	const double *var_rs, const int *svh, const double *dr, const nav_t *nav,
	const double *x, rtk_t *rtk, double *v, double *H, double *var)
{
	prcopt_t *opt = &rtk->opt;
	double freq[NFREQ], y, r, cdtr, bias, C, rr[3], pos[3], e[3], dtdx[3], L[NFREQ], P[NFREQ], Lc, Pc;				
	double dtrp = 0.0, dion = 0.0, vart = 0.0, vari = 0.0, rDCB = 0.0, ifb,gamma12;
	double dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };
	double *azel;
	char str[32];
	int i, j, k, sat, sys, prn=0, nv = 0, nx = rtk->nx;

	time2str(obs[0].time, str, 0);

	for (i = 0; i<MAXSAT; i++) for (j = 0; j<opt->nf; j++) rtk->ssat[i].vsat[j] = 0;

	for (i = 0; i<3; i++) rr[i] = x[i] + dr[i];
	ecef2pos(rr, pos);
	for (i = 0; i<n&&i<MAXOBS; i++) {
		sat = obs[i].sat;
		for (j = 0; j < NFREQ; j++) freq[j] = sat2freq(sat, obs[i].code[j], nav);			
		azel = rtk->ssat[sat - 1].azel;

		if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 || satazel(pos, e, azel)<opt->elmin) {
			continue;
		}
		if (!(sys = satsys(sat, &prn)) || !rtk->ssat[sat - 1].vs ||
			satexclude(obs[i].sat, svh[i], opt)) {							
			continue;
		}
		
		/* tropospheric and ionospheric model */
		if (!model_trop(obs[i].time, pos, azel, opt, x, dtdx, nav, &dtrp, &vart) ||
			!model_iono(obs[i].time, pos, azel, opt, sat, x, nav, &dion, &vari)) {
			continue;
		}
		/* satellite and receiver antenna model */
		if (opt->posopt[0]) satantpcv(rs + i * 6, rr, nav->pcvs + sat - 1, dants);
		antmodel(opt->pcvr, opt->antdel[0], azel, opt->posopt[1], dantr);

		/* phase windup model */
		/*if (!model_phw(rtk->sol.time, sat, nav->pcvs[sat - 1].type,
			opt->posopt[2] ? 2 : 0, rs + i * 6, rr, &rtk->ssat[sat - 1].phw)) {
			continue;
		}*/

		if (opt->posopt[2]) {
			windupcorr(rtk->sol.time, rs + i * 6, rr, &rtk->ssat[sat - 1].phw);
		}
		/* corrected phase and code measurements */
		corr_meas(obs + i, nav, azel, &rtk->opt, dantr, dants, rtk->ssat[sat - 1].phw, L, P, &Lc, &Pc);			
		if (0) {
			char id[32];
			satno2id(obs[i].sat, id);
			//printf("skip=%s %3d %3d %10.4f %10.4f\n", id, prn, rtk->ssat[sat - 1].vs, r, obs[i].P[0]);
			printf("skip=%s %3d %3d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n", id, prn, rtk->ssat[sat - 1].vs, L[0], L[1], L[2], P[0], P[1], P[2]);
		}
		
		/* stack phase and code residuals {L1,P1,L2,P2,...} */
		for (j = 0; j<2 * NF(opt); j++) {										

			bias = rDCB = ifb = 0.0;

			if (opt->ionoopt == IONOOPT_IFLC) {						
				if ((y = j % 2 == 0 ? Lc : Pc) == 0.0) continue;			
			}
			else {
				if ((y = j % 2 == 0 ? L[j / 2] : P[j / 2]) == 0.0) continue;
				if (freq[j / 2] == 0.0 || freq[0] == 0.0) continue;
			}

			/* coordinate */
			for (k = 0; k<nx; k++) H[k + nx*nv] = k<3 ? -e[k] : 0.0;
			/* receiver clock */
			k = sys == SYS_GLO ? 1 : (sys == SYS_GAL ? 2 : (sys == SYS_CMP ? 3 : 0));				
			
			cdtr = x[IC(k, opt)];
			H[IC(k, opt) + nx*nv] = 1.0;
			/* Trop */
			if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
				for (k = 0; k<(opt->tropopt >= TROPOPT_ESTG ? 3 : 1); k++) {
					H[IT(opt) + k + nx*nv] = dtdx[k];
				}
			}
			/* Iono */
			C = SQR(freq[0] / freq[j / 2])*(j % 2 == 0 ? -1.0 : 1.0);			
			if (opt->ionoopt == IONOOPT_EST) {
				if (rtk->x[II(sat, opt)] == 0.0) continue;
				H[II(sat, opt) + nx*nv] = C;
				/* rDCB */
				if ((j % 2 == 1) && opt->posopt[5]) {
					k = sys == SYS_GLO ? 1 : (sys == SYS_GAL ? 2 : (sys == SYS_CMP ? 3 : 0));
					if ((rDCB = x[ID(k, opt)]) == 0.0) continue;
					H[ID(k, opt) + nx*nv] = C;
				}
			}
			/* phase bias */
			if (j % 2 == 0) {
				if ((bias = x[IB(sat, j / 2, opt)]) == 0.0) continue;
				H[IB(sat, j / 2, opt) + nx*nv] = 1.0;
			}
			/* Inter frequency bias */
			if (j % 2 != 0 && j / 2 == 2) {
				if ((ifb = x[IS(sat, opt)]) == 0.0) continue;
				H[IS(sat, opt) + nx*nv] = 1.0;
			}
			/* residual */
			v[nv] = y - (r + cdtr - CLIGHT*dts[i * 2] + dtrp + C*dion + C*+rDCB + ifb + bias);
			if (j % 2 == 0) rtk->ssat[sat - 1].resc[j / 2] = v[nv];
			else        rtk->ssat[sat - 1].resp[j / 2] = v[nv];
			/* variance */
			if (j % 2 == 0) var[nv] = varerr(obs[i].sat, azel[1], j / 2, j % 2, opt) + vart + SQR(C)*vari + var_rs[i];							
			else		var[nv] = varerr(obs[i].sat, azel[1], j / 2, j % 2, opt) + vart + SQR(C)*vari + var_rs[i] + SQR(ERR_CBIAS);			
			if (sys == SYS_GLO&&j % 2 == 1) var[nv] += VAR_GLO_IFB;	
			char id[32];
			satno2id(sat,id);
			trace(2, "%s %s %s%d res=%9.4f sig=%9.4f el=%4.1f\n", str, id, j % 2 ? "P" : "L", j / 2 + 1, v[nv], sqrt(var[nv]), azel[1] * R2D);
			//trace(1, "%s sat=%2d %s%d res=%9.4f sig=%9.4f el=%4.1f\n", str, sat, j % 2 ? "P" : "L", j / 2 + 1, v[nv], sqrt(var[nv]), azel[1] * R2D);

			/* reject satellite by pre-fit residuals */
			if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) {				
				rtk->ssat[sat - 1].rejc[j % 2]++;
				trace(2, "outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n", 0, str, sat, j % 2 ? "P" : "L", j / 2 + 1, v[nv], azel[1] * R2D);
				continue;
			}
			if (j % 2 == 0) rtk->ssat[sat - 1].vsat[j / 2] = 1;
			nv++;
		}	//end styple
	}	//end obs
	return nv;
}
static int ppp_pos(const obsd_t *obs, int n, const double *rs, const double *dts, 
	const double *var_rs, const int *svh, const double *dr,
	const nav_t *nav, const double *x, rtk_t *rtk)
{
	double freq[NFREQ];
	const double *azel;
	prcopt_t *opt = &rtk->opt;
	double y, r, cdtr, bias, C, rr[3], pos[3], e[3], dtdx[3], L[NFREQ], P[NFREQ], Lc, Pc;
	double dtrp = 0.0, dion = 0.0, vart = 0.0, vari = 0.0;
	double dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };
	char str[32];
	double ve[MAXOBS * 2 * NFREQ] = { 0 }, vmax = 0;
	int ne = 0, obsi[MAXOBS * 2 * NFREQ] = { 0 }, frqi[MAXOBS * 2 * NFREQ], maxobs, maxfrq;
	int i, j, k, sat, sys, prn, nx = rtk->nx, stat = 1;
	double v, var;
	time2str(obs[0].time, str, 2);

	for (i = 0; i<3; i++) rr[i] = x[i] + dr[i];
	ecef2pos(rr, pos);

	for (i = 0; i<n&&i<MAXOBS; i++) {
		sat = obs[i].sat;
		for (j = 0; j < NFREQ; j++) freq[j] = sat2freq(sat, obs[i].code[j], nav);
		azel = rtk->ssat[sat - 1].azel;

		if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 || satazel(pos, e, azel)<opt->elmin) {
			continue;
		}
		if (!(sys = satsys(sat, &prn)) || !rtk->ssat[sat - 1].vs ||
			satexclude(obs[i].sat, svh[i], opt)) {
			continue;
		}

		/* tropospheric and ionospheric model */
		if (!model_trop(obs[i].time, pos, azel, opt, x, dtdx, nav, &dtrp, &vart) ||
			!model_iono(obs[i].time, pos, azel, opt, sat, x, nav, &dion, &vari)) {
			continue;
		}
		/* satellite and receiver antenna model */
		if (opt->posopt[0]) satantpcv(rs + i * 6, rr, nav->pcvs + sat - 1, dants);
		antmodel(opt->pcvr, opt->antdel[0], azel, opt->posopt[1], dantr);

		/* phase windup model */
		/*if (!model_phw(rtk->sol.time, sat, nav->pcvs[sat - 1].type, opt->posopt[2] ? 2 : 0, rs + i * 6, rr, &rtk->ssat[sat - 1].phw)) {
			continue;
		}*/
		if (opt->posopt[2]) {
			windupcorr(rtk->sol.time, rs + i * 6, rr, &rtk->ssat[sat - 1].phw);
		}
		/* corrected phase and code measurements */
		corr_meas(obs + i, nav, azel, &rtk->opt, dantr, dants, rtk->ssat[sat - 1].phw, L, P, &Lc, &Pc);

		/* stack phase and code residuals {L1,P1,L2,P2,...} */
		for (j = 0; j<2 * NF(opt); j++) {

			bias = 0.0;

			if (opt->ionoopt == IONOOPT_IFLC) {
				if ((y = j % 2 == 0 ? Lc : Pc) == 0.0) continue;
			}
			else {
				if ((y = j % 2 == 0 ? L[j / 2] : P[j / 2]) == 0.0) continue;
				if (freq[j / 2] == 0.0 || freq[0] == 0.0) continue;
			}
			C = SQR(freq[0] / freq[j / 2])*(j % 2 == 0 ? -1.0 : 1.0);

			/* receiver clock */
			sys = satsys(sat, NULL);
			k = sys == SYS_GLO ? 1 : (sys == SYS_GAL ? 2 : (sys == SYS_CMP ? 3 : 0));
			cdtr = x[IC(k, opt)];

			/* phase bias */
			if (j % 2 == 0) {
				if ((bias = x[IB(sat, j / 2, opt)]) == 0.0) continue;
			}

			/* residual */
			v = y - (r + cdtr - CLIGHT*dts[i * 2] + dtrp + C*dion + bias);

			if (j % 2 == 0) rtk->ssat[sat - 1].resc[j / 2] = v;
			else        rtk->ssat[sat - 1].resp[j / 2] = v;

			/* variance */
			if (j % 2 == 0) var = varerr(obs[i].sat, azel[1], j / 2, j % 2, opt) + vart + SQR(C)*vari + var_rs[i];
			else		var = varerr(obs[i].sat, azel[1], j / 2, j % 2, opt) + vart + SQR(C)*vari + var_rs[i] + SQR(ERR_CBIAS);
			if (sys == SYS_GLO&&j % 2 == 1) var += VAR_GLO_IFB;

			/* record large post-fit residuals */
			if (fabs(v)>sqrt(var)*THRES_REJECT) {
				obsi[ne] = i; frqi[ne] = j; ve[ne] = v; ne++;
			}
			char id[32];
			satno2id(sat, id);
			trace(1, "%s %s %s%d res=%9.4f sig=%9.4f el=%4.1f\n", str, id, j % 2 ? "P" : "L", j / 2 + 1, v, sqrt(var), azel[1] * R2D);
		}	//end styple
	}	//end obs
	/* reject satellite with large and max post-fit residual */
	if (ne>0) {
		vmax = ve[0]; maxobs = obsi[0]; maxfrq = frqi[0];
		for (j = 1; j<ne; j++) {
			if (fabs(vmax) >= fabs(ve[j])) continue;
			vmax = ve[j]; maxobs = obsi[j]; maxfrq = frqi[j];
		}
		sat = obs[maxobs].sat;
		trace(2, "outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n", 1, str, sat, maxfrq % 2 ? "P" : "L", maxfrq / 2 + 1, vmax, azel[1] * R2D);
		rtk->ssat[sat - 1].rejc[maxfrq % 2]++; stat = 0;
	}
	return stat;
}
/* update solution status ----------------------------------------------------*/
static void update_stat(rtk_t *rtk, const obsd_t *obs, int n, int stat)
{
	const prcopt_t *opt = &rtk->opt;
	int i, j, sat, slip, num;
	double sig0, diffamb;

	/* test # of valid satellites */
	rtk->sol.ns = 0;

	for (i = 0; i<n&&i<MAXOBS; i++) {
		for (j = 0; j < NF(opt); j++) {
			if (!rtk->ssat[obs[i].sat - 1].vsat[j]) continue;
			rtk->ssat[obs[i].sat - 1].lock[j]++;
			rtk->ssat[obs[i].sat - 1].outc[j] = 0;
		}
		if (!rtk->ssat[obs[i].sat - 1].vsat[0]) continue;
		rtk->sol.ns++;
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
	}
	else {
		for (i = 0; i<3; i++) {
			rtk->sol.rr[i] = rtk->x[i];
			rtk->sol.qr[i] = (float)rtk->P[i + i*rtk->nx];
		}
		rtk->sol.qr[3] = (float)rtk->P[1];
		rtk->sol.qr[4] = (float)rtk->P[2 + rtk->nx];
		rtk->sol.qr[5] = (float)rtk->P[2];
	}
	rtk->sol.dtr[0] = rtk->x[IC(0, opt)] / CLIGHT;
	rtk->sol.dtr[1] = (rtk->x[IC(1, opt)] - rtk->x[IC(0, opt)]) / CLIGHT;

	for (i = 0; i<n&&i<MAXOBS; i++) {
		rtk->ssat[obs[i].sat - 1].snr[0] = MIN(obs[i].SNR[0], obs[i].SNR[1]); //obs[i].SNR[j];
	}
	for (i = 0; i<MAXSAT; i++) {
		for (j = 0; j < NF(opt); j++) {
			if (rtk->ssat[i].slip[j] & 3) rtk->ssat[i].slipc[j]++;					//弧段数
		}
	}
	for (i = 0; i < n; i++) {
		sat = obs[i].sat;
		slip = rtk->ssat[obs[i].sat - 1].slip[0];
		for (j = 0; j < 4; j++) {
			double LC = j == 0 ? rtk->ssat[sat - 1].gf : j == 1 ? rtk->ssat[sat - 1].mw : j == 2 ? rtk->ssat[sat - 1].gfPL : rtk->ssat[sat - 1].gf2;
			if (LC) {
				rtk->ambc[sat - 1].epoch[j] = rtk->sol.time;
			}
			if (slip || j == 0 || j == 3 || LC == 0.0) {
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
				rtk->ambc[sat - 1].epoch[j] = rtk->sol.time;
			}
		}
	}
}
/* precise point positioning -------------------------------------------------*/
static void pppos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	const prcopt_t *opt = &rtk->opt;
	double *rs, *dts, *var, *v, *H, *r, *R, *xp, *Pp, dr[3] = { 0 }, std[3];
	char str[32];
	int i, j, k, f, nv, info, svh[MAXOBS], stat = SOLQ_SINGLE;

	time2str(obs[0].time, str, 2);
	trace(2, "pppos   : time=%s nx=%d n=%d\n", str, rtk->nx, n);

	rs = mat(6, n); dts = mat(2, n); var = mat(1, n);

	for (i = 0; i<MAXSAT; i++) for (j = 0; j<opt->nf; j++) rtk->ssat[i].fix[j] = 0;
	/* LC observation */
	lc_ppp(rtk, obs, n, nav);

	/* temporal update of ekf states */
	udstate_ppp(rtk, obs, n, nav);

	/* satellite positions and clocks */
	satposs(obs[0].time, obs, n, nav, rtk->opt.sateph, rs, dts, var, svh);

	/* exclude measurements of eclipsing satellite (block IIA) */
	if (rtk->opt.posopt[3]) {
		testeclipse(obs, n, nav, rs);
	}
	/* earth tides correction */
	if (opt->tidecorr) {
		tidedisp(gpst2utc(obs[0].time), rtk->x, opt->tidecorr == 1 ? 1 : 7, &nav->erp, opt->odisp[0], dr);
	}
	nv = n*rtk->opt.nf * 2 + MAXSAT + 3;
	xp = mat(rtk->nx, 1); Pp = zeros(rtk->nx, rtk->nx);
	v = mat(nv, 1); H = mat(rtk->nx, nv); r = mat(nv, 1); R = mat(nv, nv);

	matcpy(xp, rtk->x, rtk->nx, 1);
	matcpy(Pp, rtk->P, rtk->nx, rtk->nx);
	for (i = 0; i < 1; i++) {
		/* prefit residuals */
		if (!(nv = ppp_res(obs, n, rs, dts, var, svh, dr, nav, xp, rtk, v, H, r))) {		//v H R	exc	
			trace(2, "%s ppp (%d) no valid obs data\n", str, i + 1);
			break;
		}
		if (rtk->opt.ionoopt == IONOOPT_IFLC && nv < 2 * 4) break;
		if (rtk->opt.ionoopt == IONOOPT_EST && nv < 4 * 4) break;
		/* constraint to local correction */
		nv += const_corr(obs, n, nav, dr, xp, rtk, v + nv, H + nv*rtk->nx, r + nv);					
		for (j = 0; j < nv; j++) for (k = 0; k < nv; k++) {
			R[j + k*nv] = j == k ? r[j] : 0.0;
		}
		/* measurement update of ekf states */
		if ((info = filter(xp, Pp, H, v, R, rtk->nx, nv))) {								//xp	Pp(k/k-1）
			trace(2, "%s ppp (%d) filter error info=%d\n", str, i + 1, info);
			break;
		}
		//ppp_pos(obs, n, rs, dts, var, svh, dr, nav, xp, rtk);
		stat = SOLQ_PPP;
	}

	double rr[3];
	for (i = 0; i < 3; i++) rr[i] = rtk->sol.rr[i] - rtk->opt.ru[i];
	if (stat == SOLQ_PPP) {
		matcpy(rtk->x, xp, rtk->nx, 1);
		matcpy(rtk->P, Pp, rtk->nx, rtk->nx);

		/* ambiguity resolution in ppp */
		/*if ((opt->modear >= ARMODE_PPPAR)) {
		if (pppamb(rtk, obs, n, nav)) {
		stat = SOLQ_FIX;
		}
		}*/
		/* update solution status */
		update_stat(rtk, obs, n, stat);
	}
	else {
		rtk->sol.ns = 0;
		for (i = 0; i<n&&i<MAXOBS; i++) {
			if (!rtk->ssat[obs[i].sat - 1].vsat[0]) continue;
			rtk->ssat[obs[i].sat - 1].lock[0]++;
			rtk->ssat[obs[i].sat - 1].outc[0] = 0;
			rtk->sol.ns++;
		}
		rtk->sol.stat = stat;
	}

	free(rs); free(dts); free(var);
	free(xp); free(Pp);	 free(v); free(H); free(R);
}
extern int ppprocess(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)				
{
	prcopt_t *opt = &rtk->opt;
	gtime_t time;
	int nu = n;
	char msg[128] = "";
	trace(3, "rtkpos  : time=%s n=%d\n", time_str(obs[0].time, 3), n);

	time = rtk->sol.time;

	/* rover position by single point positioning */
	if (!pntpos(obs, nu, nav, &rtk->opt, &rtk->sol, NULL, rtk->ssat, msg)) {					
		showmsg(msg);
		return 0;
	}

	if (time.time != 0) rtk->tt = timediff(rtk->sol.time, time);
	/* precise point positioning */
	if (opt->mode >= PMODE_PPP_KINEMA) {
		pppos(rtk, obs, nu, nav);
		outsolstat(rtk, obs, nu, nav);
		return 1;
	}
}
