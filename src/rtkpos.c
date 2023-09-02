/* RTK */
#include <stdarg.h>
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/
#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))

#define MAX_ITER     8               /* max number of iterations */
#define MIN_NSAT_SOL 4               /* min satellite number for solution */
#define REJION		30 

#define VAR_POS     SQR(30.0) /* initial variance of receiver pos (m^2) */
#define VAR_VEL     SQR(10.0) /* initial variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0) /* initial variance of receiver acc ((m/ss)^2) */
#define VAR_ZWD     SQR(0.12)		 /* init variance ZWD */ 
#define VAR_GRA     SQR(0.001) /* initial variance of gradient (m^2) */
#define VAR_BIAS    SQR(60.0)        /* init variance phase-bias (m^2) */

#define INIT_ZWD    0.15     /* initial zwd (m) */
#define GAP_RESION  120      /* gap to reset ionosphere parameters (epochs) */
#define MAXACC      30.0     /* max accel for doppler slip detection (m/s^2) */
#define NONLIEAR	0.1 	 /* threshold for nonliearity (v.2.3.0) */
#define TTOL_MOVEB  (1.0+2*DTTOL)	/* time sync tolerance for moving-baseline (s) */

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:9)
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt))
#define NX(opt)     (NR(opt)+NB(opt))

/* state variable index */
#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1) /* phase bias (s:satno,f:freq) */

typedef struct {
	unsigned char n;
	int flg[256];
} exc_t;

/* number of estimated states ------------------------------------------------*/
extern int relnx(const prcopt_t *opt)
{
	return NX(opt);
}
/* standard deviation of state -----------------------------------------------*/
static double STD(rtk_t *rtk, int i)
{
    if (rtk->sol.stat==SOLQ_FIX) return SQRT(rtk->Pa[i+i*rtk->nx]);
    return SQRT(rtk->P[i+i*rtk->nx]);
}
/* write solution status to buffer -------------------------------------------*/
extern int reloutstat(rtk_t *rtk, char *buff)
{
    ssat_t *ssat;
    double tow,pos[3], *x;
    int i, j, k, week;
    char id[32], *p=buff;
    
    if (!rtk->sol.stat) return 0;

	tow=time2gpst(rtk->sol.time,&week);

	x=rtk->sol.stat==SOLQ_FIX?rtk->xa:rtk->x;
    
    /* receiver position */
   	p+=sprintf(p,"$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",week,tow,
                rtk->sol.stat,x[0],x[1],x[2],STD(rtk,0),STD(rtk,1),STD(rtk,2));

    /* receiver clocks */
    p+=sprintf(p,"$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n",
               week,tow,rtk->sol.stat,1,rtk->sol.dtr[0]*1E9,rtk->sol.dtr[1]*1E9,
               rtk->sol.dtr[2]*1E9,rtk->sol.dtr[3]*1E9);
    
    // /* ionospheric parameters */
    // if (rtk->opt.ionoopt==IONOOPT_EST) {
    //     for (i=0;i<MAXSAT;i++) {
    //         ssat=rtk->ssat+i;
    //         if (!ssat->vs) continue;
    //         satno2id(i+1,id);
    //         j=II(i+1,&rtk->opt);
    //         p+=sprintf(p,"$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f\n",week,tow,rtk->sol.stat,
    //                    id,ssat->azel[0]*R2D,ssat->azel[1]*R2D,x[j]);
    //     }
    // }
    // /* tropospheric parameters */
    // if ((rtk->opt.tropopt==TROPOPT_EST||rtk->opt.tropopt==TROPOPT_ESTG)) {
    //     for (i=0;i<2;i++) {
    //         j=IT(i,&rtk->opt);
    //         p+=sprintf(p,"$TROP,%d,%.3f,%d,%d,%.4f\n",week,tow,
    //                    rtk->sol.stat,i+1,rtk->x[j],x[j]);
    //     }
    // }

	/* ambiguity parameters */
    // for (i = 0; i<MAXSAT; i++) {
	// 	ssat=rtk->ssat+i;
    //     if (!ssat->vs) continue;
	// 	for (j = 0; j<NF(&rtk->opt); j++) {
    //     	k = IB(i+1,j,&rtk->opt);
    //     	if (rtk->x[k]==0.0) continue;
    //     	satno2id(i+1, id);
    //     	p+=sprintf(p,"$AMB,%d,%.3f,%d,%s,%d,%.4f,%.4f,%d,%d,%d\n",week,tow,
    //                	   rtk->sol.stat,id,j+1,x[k],STD(rtk,k), 
	// 				   ssat->lock[j], ssat->slip[j]&1, ssat->outc[j]);
	// 	}
    // }

    return (int)(p-buff);
}
/* output solution status ----------------------------------------------------*/
extern void outrel(FILE *fp, const rtk_t *rtk)
{
	ssat_t *ssat;
	double tow;
    char buff[MAXSOLMSG+1], id[32];
    int i, j, n, week, nfreq, nf=NF(&rtk->opt);
    
    if (!fp||rtk->sol.stat==SOLQ_NONE) return;
    
    /* write solution status */
    n=reloutstat(rtk, buff); buff[n]='\0';
    
    fputs(buff, fp);

	tow=time2gpst(rtk->sol.time,&week);
    nfreq=rtk->opt.mode>=PMODE_DGPS?nf:1;

	/* write residuals and status */
    for (i = 0; i<MAXSAT; i++) {
        ssat=rtk->ssat+i;
        if (!ssat->vs) continue;
        satno2id(i+1, id);
        for (j=0;j<nfreq;j++) {
            fprintf(fp,"$SAT,%d,%.3f,%s,%d,%.1f,%.1f,%.4f,%.4f,%d,%d,%d,%d\n",
                    week,tow,id,j+1,ssat->azel[0]*R2D,ssat->azel[1]*R2D,
                    ssat->resp[j],ssat->resc[j],ssat->vsat[j],
                    ssat->fix[j],ssat->slip[j]&3,ssat->lock[j],ssat->outc[j]);
        }
    }
}
/* select common satellites between rover and reference station --------------*/
static int selsat(const obsd_t *obs, double *azel, int nu, int nr,
	const prcopt_t *opt, int *sat, int *iu, int *ir)
{
	int i, j, k = 0;

	trace(3, "selsat  : nu=%d nr=%d\n", nu, nr);

	for (i = 0, j = nu; i<nu&&j<nu + nr; i++, j++) {
		if (obs[i].sat<obs[j].sat) j--;
		else if (obs[i].sat>obs[j].sat) i--;
		else if (azel[1 + j * 2] >= opt->elmin) { /* elevation at base station */
			sat[k] = obs[i].sat; iu[k] = i; ir[k++] = j;
			trace(4, "(%2d) sat=%3d iu=%2d ir=%2d\n", k - 1, obs[i].sat, i, j);
		}
	}
	return k;
}
/* single-differenced observable ---------------------------------------------*/
static void sdobs(const nav_t *nav, const obsd_t *obs, int iu, int ir, double *L, 
	double *P, double *Lc, double *Pc)
{
	int i, j, k;
	double Pi, Pj, Li, Lj, freq[NFREQ], C1, C2;

	for (k = 0; k < NFREQ; k++) {
		freq[k] = sat2freq(obs[iu].sat, obs[iu].code[k], nav);
		L[k] = P[k] = 0.0;
	}

	for (k = 0; k < NFREQ; k++) {
		Pi = obs[iu].P[k];
		Pj = obs[ir].P[k];
		Li = obs[iu].L[k];
		Lj = obs[ir].L[k];
		if(Pi!=0.0 && Pj!=0.0) P[k] = Pi - Pj;
		if(Li!=0.0 && Lj!=0.0 && freq[k]!=0.0) L[k] = CLIGHT * (Li - Lj) / freq[k];
	}

	/* iono-free LC */				
	*Lc = *Pc = 0.0;
	if (freq[0] == 0.0 || freq[1] == 0.0) return;

	C1 =  SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));
	C2 = -SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));

	if (L[0] != 0.0&&L[1] != 0.0) *Lc = C1*L[0] + C2*L[1];
	if (P[0] != 0.0&&P[1] != 0.0) *Pc = C1*P[0] + C2*P[1];
}
/* LC observation --------------------------------------------------------------- */
static void lc_rel(rtk_t *rtk, const obsd_t *obs, const int *sat,
	const int *iu, const int *ir, int ns, const nav_t *nav) 	
{
	int i, k;
	double freq[NFREQ], L[NFREQ], P[NFREQ], Lc, Pc, alfa, lambda;
	for (i = 0; i<MAXSAT; i++) for(k = 0; k<7; k++) rtk->ssat[i].LC[k] = 0.0;
	for (i = 0; i < ns; i++) {

		sdobs(nav, obs, iu[i], ir[i], L, P, &Lc, &Pc);
		for (k = 0; k < NFREQ; k++) freq[k] = sat2freq(sat[i], obs[iu[i]].code[k], nav);

		// GF:L1-L2
		if (freq[0] * freq[1] * L[0] * L[1] != 0.0) {
			rtk->ssat[sat[i] - 1].LC[0] = L[0] - L[1];
		} 

		// GF:L2-L3
		if (freq[1] * freq[2] * L[1] * L[2] != 0.0) {
			rtk->ssat[sat[i] - 1].LC[1] = L[1] - L[2];
		}

		// MW:L1L2-P1P2
		if (freq[0] * freq[1] * L[0] * L[1] * P[0] * P[1] != 0.0) {
			lambda = CLIGHT / (freq[0] - freq[1]);
			rtk->ssat[sat[i] - 1].LC[2] = (freq[0] * L[0] - freq[1] * L[1]) / CLIGHT - 
				(freq[0] * P[0] + freq[1] * P[1]) / (freq[0] + freq[1]) / lambda;
		}

		// MW:L2L3-P2P3 
		if (freq[1] * freq[2] * L[1] * L[2] * P[1] * P[2] != 0.0) {
			lambda = CLIGHT / (freq[1] - freq[2]);
			rtk->ssat[sat[i] - 1].LC[3] = (freq[1] * L[1] - freq[2] * L[2]) / CLIGHT - 
				(freq[1] * P[1] + freq[2] * P[2]) / (freq[1] + freq[2]) / lambda;
		}

		// MW:L2L3-P1P2
		if (freq[1] * freq[2] * L[1] * L[2] * P[0] * P[1] != 0.0) {
			alfa = (SQR(freq[0])*SQR(freq[1]) - SQR(freq[1])*freq[1] * freq[2]) / 
				(freq[1] * freq[2] * (SQR(freq[0]) - SQR(freq[1])));
			lambda = CLIGHT / (freq[1] - freq[2]);
			rtk->ssat[sat[i] - 1].LC[4] = (freq[1] * L[1] - freq[2] * L[2]) / CLIGHT - 
				((1 - alfa) * P[0] + alfa* P[1]) / lambda;				
		}
	}
}
/* detect cycle slip by LLI ------------------------------------------------ */
static uint8_t detslp_ll(rtk_t *rtk, const obsd_t *obs, int i, int rcv)
{
    uint32_t LLI, byte, slip = 0;
    int f, sat = obs[i].sat;
    
    for (f=0;f<NF(&rtk->opt);f++) {
        
        if (obs[i].L[f]==0.0||
            fabs(timediff(obs[i].time,rtk->ssat[sat-1].prev_T[rcv-1]))<DTTOL) {
            continue;
        }
        /* restore previous LLI */
        if (rcv==1) LLI=getbitu(&rtk->ssat[sat-1].slip[f],0,2); /* rover */
        else        LLI=getbitu(&rtk->ssat[sat-1].slip[f],2,2); /* base  */
        
        /* detect slip by cycle slip flag in LLI */
		byte = (rtk->tt>=0.0) ? obs[i].LLI[f] : LLI; 
		if(byte & 1) slip |= 0x1<<f;
        
        /* detect slip by parity unknown flag transition in LLI */
        if (((LLI&2)&&!(obs[i].LLI[f]&2))||(!(LLI&2)&&(obs[i].LLI[f]&2))) {
            slip |= 0x1<<f;
        }
        /* save current LLI */
        if (rcv==1) setbitu(&rtk->ssat[sat-1].slip[f],0,2,obs[i].LLI[f]);
        else        setbitu(&rtk->ssat[sat-1].slip[f],2,2,obs[i].LLI[f]);
    }
	return slip;
}
/* cycle Jump detect ------------------------------------------------------- */
static void cycJump_detect(rtk_t *rtk, const obsd_t *obs, const int *sat,
	const int *iu, const int *ir, int ns, const nav_t *nav) 	
{
	int i, k, nf, slip, lli, mw12, mw23, gf12, gf13, dopp;
	double thresgf, thresmw, diffamb;

	nf = rtk->opt.nf;

	thresgf = rtk->opt.thresslipgf[0];
	thresmw = rtk->opt.thresslipmw[0];

	for (i = 0; i<ns; i++) {
		slip = lli = gf12 = gf13 = mw12 = mw23 = dopp = 0;

		// LLI
		lli |= detslp_ll(rtk, obs, iu[i], 1);
		lli |= detslp_ll(rtk, obs, ir[i], 2);

		// GF12
		if (rtk->ssat[sat[i]-1].LC[0] != 0.0 && rtk->ssat[sat[i]-1].prev_LC[0] != 0.0) {
			diffamb = rtk->ssat[sat[i] - 1].LC[0] - rtk->ssat[sat[i]-1].prev_LC[0];
			gf12 = (fabs(diffamb) > thresgf) ? 3 : 0;
		}

		// GF23
		if (rtk->ssat[sat[i]-1].LC[1] != 0.0 && rtk->ssat[sat[i]-1].prev_LC[1] != 0.0) {
			diffamb = rtk->ssat[sat[i] - 1].LC[1] - rtk->ssat[sat[i]-1].prev_LC[1];
			gf13 = (fabs(diffamb) > thresgf) ? 5 : 0;
		}

		// MW L1L2 P1P2
		if (rtk->ssat[sat[i]-1].LC[2] != 0.0 && rtk->ssat[sat[i]-1].prev_LC[2] != 0.0) {
			diffamb = rtk->ssat[sat[i] - 1].LC[2] - rtk->ssat[sat[i]-1].prev_LC[2];
			mw12 = (fabs(diffamb) > thresmw) ? 3 : 0;
		}

		// MW L2L3 P2P3
		if (rtk->ssat[sat[i]-1].LC[3] != 0.0 && rtk->ssat[sat[i]-1].prev_LC[3] != 0.0) {
			diffamb = rtk->ssat[sat[i] - 1].LC[3] - rtk->ssat[sat[i]-1].prev_LC[3];
			mw23 = (fabs(diffamb) > thresmw) ? 6 : 0;
		}

		slip = lli | mw12 | mw23 | gf12 | gf13 | dopp;
		for (k = 0; k < nf; k++) {
			rtk->ssat[sat[i] - 1].slip[k] |= (slip>>k) & 0x1;
			rtk->ssat[sat[i] - 1].half[k]  = !((obs[iu[i]].LLI[k]&2) || 
											   (obs[ir[i]].LLI[k]&2)); 
		}
	}
}
/* update information -------------------------------------------------------- */
static void update_sat_amb(rtk_t *rtk, int *stat, const obsd_t *obs, const int *sat,
	const int *iu, const int *ir, int ns) 
{
	int i, k, slip, num, nf = NF(&rtk->opt);
	double sig0, diffamb, LC;

	rtk->sol.ns = 0;
	for (i = 0; i < MAXSAT; i++) {

		/* update satellite status */
		for (k = 0; k < nf; k++) {
			if (rtk->ssat[i].vsat[k]) {
				if(rtk->ssat[i].slip[k]&1) {
					rtk->ssat[i].lock[k] = 1;		
					rtk->ssat[i].outc[k] = 0;
				}
				else {
					rtk->ssat[i].lock[k]++;		
					rtk->ssat[i].outc[k] = 0;
				}
			}
			else {
				rtk->ssat[i].lock[k] = 0;
				rtk->ssat[i].outc[k]++;				
			}
		}
		if (rtk->ssat[i].vsat[0]) rtk->sol.ns++;
	}

	for (i = 0; i<ns; i++) {

		/* update satellite status */
		rtk->ssat[sat[i] - 1].prev_T[0] = obs[iu[0]].time;
		rtk->ssat[sat[i] - 1].prev_T[1] = obs[ir[0]].time;
		for (k = 0; k < nf; k++) {
			slip = rtk->ssat[sat[i] - 1].slip[k];	

			rtk->ssat[sat[i] - 1].prev_P[0][k] = obs[iu[i]].P[k];
			rtk->ssat[sat[i] - 1].prev_L[0][k] = obs[iu[i]].L[k];
			rtk->ssat[sat[i] - 1].prev_D[0][k] = obs[iu[i]].D[k];
			rtk->ssat[sat[i] - 1].prev_P[1][k] = obs[ir[i]].P[k];
			rtk->ssat[sat[i] - 1].prev_L[1][k] = obs[ir[i]].L[k];
			rtk->ssat[sat[i] - 1].prev_D[1][k] = obs[ir[i]].D[k];
		}
		for (k = 0; k < 4; k++) {
			rtk->ssat[sat[i] - 1].prev_LC[k] = rtk->ssat[sat[i] - 1].LC[k];
		}
	}

	/* update solution status */
	*stat = rtk->sol.ns<MIN_NSAT_SOL ? SOLQ_NONE : *stat;
}
/* baseline length -----------------------------------------------------------*/
static double baseline(const double *ru, const double *rb, double *dr)
{
	int i;
	for (i = 0; i<3; i++) dr[i] = ru[i] - rb[i];
	return norm(dr, 3);
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
/* temporal update of position/velocity/acceleration -------------------------*/
static void udpos(rtk_t *rtk, double tt)
{
	double *F, *P, *FP, *x, *xp, pos[3], Q[9] = { 0 }, Qv[9], var = 0.0;
	int i, j, *ix, nx;

	trace(3, "udpos   : tt=%.3f\n", tt);

	/* fixed mode */
	if (rtk->opt.mode == PMODE_FIXED) {
		for (i = 0; i<3; i++) initx(rtk, rtk->opt.ru[i], 1E-8, i);
		return;
	}
	/* initialize position for first epoch */
	if (norm(rtk->x, 3) <= 0.0) {
		for (i = 0; i<3; i++) initx(rtk, rtk->sol.rr[i], SQR(rtk->opt.P0[0]), i);
		if (rtk->opt.dynamics) {
			for (i = 3; i<6; i++) initx(rtk, rtk->sol.rr[i], SQR(rtk->opt.P0[6]), i);
			for (i = 6; i<9; i++) initx(rtk, 1E-6, SQR(rtk->opt.P0[7]), i);
		}
	}
	/* static mode */
	if (rtk->opt.mode == PMODE_STATIC) return;

	/* kinmatic mode without dynamics */
	if (!rtk->opt.dynamics) {
		for (i = 0; i<3; i++) initx(rtk, rtk->sol.rr[i], SQR(rtk->opt.P0[0]), i);  
		return;
	}
	/* check variance of estimated postion */
	for (i = 0; i<3; i++) var += rtk->P[i + i*rtk->nx];
	var /= 3.0;

	if (var>SQR(rtk->opt.P0[0])) {
		/* reset position with large variance */
		for (i = 0; i<3; i++) initx(rtk, rtk->sol.rr[i], SQR(rtk->opt.P0[0]), i);
		for (i = 3; i<6; i++) initx(rtk, rtk->sol.rr[i], SQR(rtk->opt.P0[6]), i);
		for (i = 6; i<9; i++) initx(rtk, 1E-6, SQR(rtk->opt.P0[7]), i);
		trace(2, "reset rtk position due to large variance: var=%.3f\n", var);
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
		F[i + (i + 3)*nx] = tt;
	}
	for (i = 0; i<3; i++) {
		F[i + (i + 6)*nx] = SQR(tt) / 2.0;
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
	Q[0] = Q[4] = SQR(rtk->opt.Qt[6])*fabs(tt);
	Q[8] = SQR(rtk->opt.Qt[7])*fabs(tt);
	ecef2pos(rtk->x, pos);
	covecef(pos, Q, Qv);
	for (i = 0; i<3; i++) for (j = 0; j<3; j++) {
		rtk->P[i + 6 + (j + 6)*rtk->nx] += Qv[i + j * 3];
	}
	free(ix); free(F); free(P); free(FP); free(x); free(xp);
}
/* temporal update of ionospheric parameters ---------------------------------*/
static void udion(rtk_t *rtk, double tt, double bl, const int *sat, int ns)
{
	int i, j;
	double el, fact;

	for (i = 1; i <= MAXSAT; i++) {
		j = II(i, &rtk->opt);
		if (rtk->x[j] != 0.0&&
			rtk->ssat[i - 1].outc[0]>GAP_RESION&&rtk->ssat[i - 1].outc[1]>GAP_RESION)
			rtk->x[j] = 0.0;
	}
	for (i = 0; i<ns; i++) {
		j = II(sat[i], &rtk->opt);

		if (rtk->x[j] != 0.0) {
			/* elevation dependent factor of process noise */
			el = rtk->ssat[sat[i] - 1].azel[1];
			fact = cos(el);
			rtk->P[j + j*rtk->nx] += SQR(rtk->opt.Qt[3] * bl / 1E4*fact)*fabs(tt);	
		}
		else {
			initx(rtk, 1E-6, SQR(rtk->opt.P0[3] * bl / 1E4), j);	
		}
	}
}
/* temporal update of tropospheric parameters --------------------------------*/
static void udtrop(rtk_t *rtk, double tt)
{
	int i, j;

	for (i = 0; i<2; i++) {			
		j = IT(i, &rtk->opt);

		if (rtk->x[j] != 0.0) {
			rtk->P[j + j*rtk->nx] += SQR(rtk->opt.Qt[3])*fabs(tt);

			if (rtk->opt.tropopt >= TROPOPT_ESTG) {
				for (j = i + 1; j<i + 3; j++) {
					rtk->P[j + j*rtk->nx] += SQR(rtk->opt.Qt[2] * 0.3)*fabs(rtk->tt);
				}
			}
		}
		else {
			initx(rtk, INIT_ZWD, SQR(rtk->opt.P0[2]), j); /* initial zwd */

			if (rtk->opt.tropopt >= TROPOPT_ESTG) {
				for (j = i + 1; j<i + 3; j++) initx(rtk, 1E-6, VAR_GRA, j);
			}
		}
	}
}
/* temporal update of phase biases -------------------------------------------*/
static void udbias(rtk_t *rtk, double tt, const obsd_t *obs, const int *sat,
	const int *iu, const int *ir, int ns, const nav_t *nav)
{
	int i, j, f, slip, reset;
	double L[NFREQ], P[NFREQ], Lc, Pc, bias, freqi;

	for (f = 0; f<NF(&rtk->opt); f++) {												
		/* reset phase-bias if expire obs outage counter */
		for (i = 0; i<MAXSAT; i++) {
			j = IB(i + 1, f, &rtk->opt);
			reset = rtk->ssat[i].outc[f] > (unsigned int)rtk->opt.maxout;
			if (rtk->x[j] != 0.0) {	
				rtk->P[j + j*rtk->nx] += SQR(rtk->opt.Qt[5])*fabs(rtk->tt);
				if (reset) initx(rtk, 0.0, 0.0, j);
			}
		}

		for (i = 0; i<ns; i++) {
			bias = 0.0;
			j = IB(sat[i], f, &rtk->opt);
			slip = rtk->ssat[sat[i] - 1].slip[f];
			freqi = sat2freq(sat[i], obs[iu[i]].code[f], nav);

			sdobs(nav, obs, iu[i], ir[i], L, P, &Lc, &Pc);		
			
			if (rtk->x[j] != 0 && !(slip&1)) continue;

			if (rtk->opt.ionoopt == IONOOPT_IFLC) {
				if (Lc == 0.0 || Pc == 0.0) continue;
				bias = Lc - Pc;
			}
			else {
				if (L[f] == 0.0 || P[f] == 0.0) continue;	
				bias = (L[f] - P[f])/CLIGHT*freqi;
			}
			rtk->x[j] = 0.0;
			if (bias != 0.0) initx(rtk, bias, SQR(rtk->opt.P0[5]), j);
		}
	}
}
/* temporal update of states --------------------------------------------------*/
static void udstate(rtk_t *rtk, const obsd_t *obs, const int *sat,
	const int *iu, const int *ir, int ns, const nav_t *nav)
{
	double tt = rtk->tt, bl, dr[3];

	trace(3, "udstate : ns=%d\n", ns);

	/* temporal update of position/velocity/acceleration */
	udpos(rtk, tt);

	/* temporal update of ionospheric parameters */
	if (rtk->opt.ionoopt >= IONOOPT_EST) {
		bl = baseline(rtk->x, rtk->rb, dr);
		udion(rtk, tt, bl, sat, ns);				
	}
	/* temporal update of tropospheric parameters */
	if (rtk->opt.tropopt >= TROPOPT_EST) {
		udtrop(rtk, tt);
	}
	/* temporal update of phase-bias */
	if (rtk->opt.mode>PMODE_DGPS) {
		udbias(rtk, tt, obs, sat, iu, ir, ns, nav);
	}
}
/* UD (undifferenced) phase/code residual for satellite ----------------------*/
static void zdres_sat(int base, double r, const obsd_t *obs, const nav_t *nav,
	const double *azel, const double *dant, const prcopt_t *opt, double *y, double *freq)
{
	double freq1, freq2, C1, C2, dant_if;
	int i, nf = NF(opt);

	if (opt->ionoopt == IONOOPT_IFLC) { /* iono-free linear combination */
		freq1 = sat2freq(obs->sat, obs->code[0], nav);
		freq2 = sat2freq(obs->sat, obs->code[1], nav);
		if (freq1 == 0.0 || freq2 == 0.0) return;

		if (testsnr(base, 0, azel[1], obs->SNR[0] * SNR_UNIT, &opt->snrmask) ||
			testsnr(base, 1, azel[1], obs->SNR[1] * SNR_UNIT, &opt->snrmask)) return;

		C1 =  SQR(freq1) / (SQR(freq1) - SQR(freq2));
		C2 = -SQR(freq2) / (SQR(freq1) - SQR(freq2));
		dant_if = C1*dant[0] + C2*dant[1];

		if (obs->L[0] != 0.0&&obs->L[1] != 0.0) {
			y[0] = C1*obs->L[0] * CLIGHT / freq1 + C2*obs->L[1] * CLIGHT / freq2 - r - dant_if;		
		}
		if (obs->P[0] != 0.0&&obs->P[1] != 0.0) {
			y[1] = C1*obs->P[0] + C2*obs->P[1] - r - dant_if;
		}
		freq[0] = 1.0;
	}
	else {
		for (i = 0; i<nf; i++) {
			if ((freq[i] = sat2freq(obs->sat, obs->code[i], nav)) == 0.0) continue;

			/* check SNR mask */
			if (testsnr(base, i, azel[1], obs->SNR[i] * SNR_UNIT, &opt->snrmask)) {
				continue;
			}
			/* residuals = observable - pseudorange */
			if (obs->L[i] != 0.0) y[i] = obs->L[i] * CLIGHT / freq[i] - r - dant[i];
			if (obs->P[i] != 0.0) y[i + nf] = obs->P[i] - r - dant[i];
		}
	}
}
/* UD (undifferenced) phase/code residuals -----------------------------------*/
static void zdres(int base, const obsd_t *obs, int n, const double *rs,
	const double *dts, const double *var, const int *svh,
	const nav_t *nav, const double *rr, const prcopt_t *opt,
	int index, double *y, double *e, double *azel, double *freq)
{
	double r, rr_[3], pos[3], dant[NFREQ] = { 0 }, disp[3];
	double zhd, zazel[] = { 0.0, 90.0*D2R };
	int i, nf = NF(opt);

	trace(3, "zdres   : n=%d\n", n);

	for (i = 0; i<n*nf * 2; i++) y[i] = 0.0;

	if (norm(rr, 3) <= 0.0) return ; /* no receiver position */

	for (i = 0; i<3; i++) rr_[i] = rr[i];

	/* earth tide correction */
	if (opt->tidecorr) {
		tidedisp(gpst2utc(obs[0].time), rr_, opt->tidecorr, &nav->erp,
			opt->odisp[base], disp);
		for (i = 0; i<3; i++) rr_[i] += disp[i];			
	}
	ecef2pos(rr_, pos);

	for (i = 0; i<n; i++) {
		/* compute geometric-range and azimuth/elevation angle */
		if ((r = geodist(rs + i * 6, rr_, e + i * 3)) <= 0.0) continue;
		if (satazel(pos, e + i * 3, azel + i * 2)<opt->elmin) continue;

		/* excluded satellite? */
		if (satexclude(obs[i].sat, svh[i], opt)) continue;

		/* satellite clock-bias */
		r += -CLIGHT*dts[i * 2];

		/* troposphere delay model (hydrostatic) */
		zhd = tropmodel(obs[0].time, pos, zazel, 1, 0.0, 1);
		r += tropmapf(obs[i].time, pos, azel + i * 2, 1, NULL)*zhd;

		/* receiver antenna phase center correction */
		antmodel(opt->pcvr + index, opt->antdel[index], azel + i * 2, opt->posopt[1], dant);
		
		/* UD phase/code residual for satellite */
		zdres_sat(base, r, obs + i, nav, azel + i * 2, dant, opt, y + i*nf * 2, freq + i*nf);		
	}
	trace(3, "rr_=%.3f %.3f %.3f\n", rr_[0], rr_[1], rr_[2]);
	trace(3, "pos=%.9f %.9f %.3f\n", pos[0] * R2D, pos[1] * R2D, pos[2]);
	for (i = 0; i<n; i++) {
		trace(3, "sat=%2d %13.3f %13.3f %13.3f %13.10f %6.1f %5.1f\n",
			obs[i].sat, rs[i * 6], rs[1 + i * 6], rs[2 + i * 6], dts[i * 2], azel[i * 2] * R2D,
			azel[1 + i * 2] * R2D);
	}
	trace(3, "y=\n"); tracemat(4, y, nf * 2, n, 13, 3);

	return 1;
}
/* time-interpolation of residuals (for post-processing) ---------------------*/
static double intpres(gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
	rtk_t *rtk, double *y)
{
	static obsd_t obsb[MAXOBS];
	static double yb[MAXOBS*NFREQ * 2], rs[MAXOBS * 6], dts[MAXOBS * 2], var[MAXOBS];
	static double e[MAXOBS * 3], azel[MAXOBS * 2], freq[MAXOBS*NFREQ];
	static int nb = 0, svh[MAXOBS * 2];				
	prcopt_t *opt = &rtk->opt;
	double tt = timediff(time, obs[0].time), ttb, *p, *q;
	int i, j, k, nf = NF(opt);

	trace(3, "intpres : n=%d tt=%.1f\n", n, tt);

	if (nb == 0 || fabs(tt)<DTTOL) {
		nb = n; for (i = 0; i<n; i++) obsb[i] = obs[i];
		return tt;
	}
	ttb = timediff(time, obsb[0].time);
	if (fabs(ttb)>opt->maxtdiff*2.0 || ttb == tt) return tt;			

	satposs(time, obsb, nb, nav, opt->sateph, rs, dts, var, svh);

	zdres(1, obsb, nb, rs, dts, var, svh, nav, rtk->rb, opt, 1, yb, e, azel, freq);

	for (i = 0; i<n; i++) {
		for (j = 0; j<nb; j++) if (obsb[j].sat == obs[i].sat) break;
		if (j >= nb) continue;
		for (k = 0, p = y + i*nf * 2, q = yb + j*nf * 2; k<nf * 2; k++, p++, q++) {
			if (*p == 0.0 || *q == 0.0) *p = 0.0; else *p = (ttb*(*p) - tt*(*q)) / (ttb - tt);			
		}
	}
	return fabs(ttb)>fabs(tt) ? ttb : tt;
}
/* test valid observation data -----------------------------------------------*/
static int validobs(int i, int j, int f, int nf, double *y)
{
	/* if no phase observable, psudorange is also unusable */
	return y[f + i*nf * 2] != 0.0&&y[f + j*nf * 2] != 0.0 &&
		(f<nf || (y[f - nf + i*nf * 2] != 0.0&&y[f - nf + j*nf * 2] != 0.0));
}
/* test satellite system (m=0:GPS/SBS,1:GLO,2:GAL,3:BDS,4:QZS,5:IRN) ---------*/
static int test_sys(int sys, int m)
{
	switch (sys) {
	case SYS_GPS: return m == 0;
	case SYS_GLO: return m == 1;
	case SYS_GAL: return m == 2;
	case SYS_CMP: return m == 3;
	}
	return 0;
}
/* exclude observation -------------------------------------------------------*/
static int exclude(int isat1, int isat2, int f, int nf, exc_t *exc) 
{
	int i, itype, ifreq, type, freq, sat1, sat2, stat = 0;

	itype = (f < nf ? 1 : 2);
	ifreq = (f % nf) + 1;

	for (i = 0; i < exc->n; i++) {
		type =  exc->flg[i] & 0xF;
		freq = (exc->flg[i] >> 4) & 0xF;
		sat1 = (exc->flg[i] >> 8) & 0xFF;
		sat2 = (exc->flg[i] >> 16) & 0xFF;
		if (type == itype){
			if (type <= 2 && (freq == ifreq)) {
				if((isat1+isat2) == (sat1+sat2) && (isat1 == sat1 || isat1 == sat2)) {
					stat = 1;     break;
				}; 
			}
		} 
	}
	return stat;
}
/* baseline length constraint ------------------------------------------------*/
static int constbl(rtk_t *rtk, const double *x, const double *P, double *v,
				   double *H, double *Ri, double *Rj, int index)
{
	int i;
	double xb[3], b[3], bb, var = 0.0;

	/* no constraint */
	if (rtk->opt.baseline[0] <= 0.0) return 0;

	/* time-adjusted baseline vector and length */
	for (i = 0; i<3; i++) {
		xb[i] = rtk->rb[i];
		b[i] = x[i] - xb[i];
	}
	bb = norm(b, 3);

	/* approximate variance of solution */
	if (P) {
		for (i = 0; i<3; i++) var += P[i + i*rtk->nx];
		var /= 3.0;
	}
	/* check nonlinearity */
	if (var>SQR(NONLIEAR*bb)) {
		trace(3, "constbl : equation nonlinear (bb=%.3f var=%.3f)\n", bb, var);
		return 0;
	}
	/* constraint to baseline length */
	v[index] = rtk->opt.baseline[0] - bb;
	if (H) {
		for (i = 0; i<3; i++) H[i + index*rtk->nx] = b[i] / bb;			
	}
	Ri[index] = 0.0;
	Rj[index] = SQR(rtk->opt.baseline[1]);

	return 1;
}
/* DD (double-differenced) measurement error covariance ----------------------*/
static void ddcov(const int *nb, int n, const double *Ri, const double *Rj,
				  int nv, double *R)
{
	int i, j, k = 0, b;

	for (i = 0; i<nv*nv; i++) R[i] = 0.0;
	for (b = 0; b<n; k += nb[b++]) {

		for (i = 0; i<nb[b]; i++) for (j = 0; j<nb[b]; j++) {
			R[k + i + (k + j)*nv] = Ri[k + i] + (i == j ? Rj[k + i] : 0.0);				
		}
	}
	trace(5, "R=\n"); tracemat(5, R, nv, nv, 8, 6);
}
/* DD (double-differenced) phase/code residuals ------------------------------*/
static int ddres(rtk_t *rtk, const nav_t *nav, double dt, const double *x, 
				 const int *sat, double *y, double *e, double *azel, double *freq, 
				 const int *iu, const int *ir, int ns, double *v, double *H, 
				 double *R, int *vflg, exc_t *exc)
{
	prcopt_t *opt = &rtk->opt;
	double bl, dr[3], posu[3], posr[3], didxi = 0.0, didxj = 0.0, *im, vartu, vartr;
	double *tropr, *tropu, *dtdxr, *dtdxu, *Ri, *Rj, freqi, freqj, *Hi = NULL;
	int i, j, k, m, f, nv = 0, nb[NFREQ * 4 * 2 + 2] = { 0 }, b = 0, sysi, sysj, nf = NF(opt);

	bl = baseline(x, rtk->rb, dr);		ecef2pos(x, posu); 		ecef2pos(rtk->rb, posr);

	Ri = mat(ns*nf * 2 + 2, 1);   Rj = mat(ns*nf * 2 + 2, 1); 	im = mat(ns, 1);
	tropu = mat(ns, 1);   tropr = mat(ns, 1); 	dtdxu = mat(ns, 3);   dtdxr = mat(ns, 3);

	/* compute factors of ionospheric and tropospheric delay */
	for (i = 0; i<ns; i++) {
		if (opt->ionoopt >= IONOOPT_EST) {
			im[i] = (ionmapf(posu, azel + iu[i] * 2) + ionmapf(posr, azel + ir[i] * 2)) / 2.0;		
		}
		if (opt->tropopt >= TROPOPT_EST) {
			tropcorr(rtk->sol.time, posu, azel + iu[i] * 2, opt, x, dtdxu + i * 3, nav, tropu+i, &vartu);
			tropcorr(rtk->sol.time, posr, azel + ir[i] * 2, opt, x, dtdxr + i * 3, nav, tropr+i, &vartr);
		}
	}
	for (m = 0; m<4; m++) {/* m=0:GPS,1:GLO,2:GAL,3:BDS */

		for (f = 0; f<nf * 2; f++) {				

			/* search reference satellite with highest elevation */
			for (i = -1, j = 0; j<ns; j++) {
				sysi = rtk->ssat[sat[j] - 1].sys;
				if (!test_sys(sysi, m)) continue;
				if (!validobs(iu[j], ir[j], f, nf, y)) continue;
				if (i<0 || azel[1 + iu[j] * 2] >= azel[1 + iu[i] * 2]) i = j;			
			}
			if (i<0) continue;

			/* make DD (double difference) */
			for (j = 0; j<ns; j++) {				
				if (i == j) continue;
				sysi = rtk->ssat[sat[i] - 1].sys;
				sysj = rtk->ssat[sat[j] - 1].sys;
				freqi = freq[f%nf + iu[i] * nf];
				freqj = freq[f%nf + iu[j] * nf];
				if (!test_sys(sysj, m)) continue;
				if (!validobs(iu[j], ir[j], f, nf, y)) continue;
				if (exclude(sat[i], sat[j], f, nf, exc)) continue;

				if (H) {
					Hi = H + nv*rtk->nx;
					for (k = 0; k<rtk->nx; k++) Hi[k] = 0.0;
				}
				/* DD residual */
				v[nv] = (y[f + iu[i] * nf * 2] - y[f + ir[i] * nf * 2]) -
					(y[f + iu[j] * nf * 2] - y[f + ir[j] * nf * 2]);

				/* partial derivatives by rover position */
				if (H) {
					for (k = 0; k<3; k++) {
						Hi[k] = -e[k + iu[i] * 3] + e[k + iu[j] * 3];
					}
				}
				/* DD ionospheric delay term */
				if (opt->ionoopt == IONOOPT_EST) {
					didxi = (f<nf ? -1.0 : 1.0)*im[i] * SQR(FREQ1 / freqi);
					didxj = (f<nf ? -1.0 : 1.0)*im[j] * SQR(FREQ1 / freqj);
					v[nv] -= didxi*x[II(sat[i], opt)] - didxj*x[II(sat[j], opt)];				
					if (H) {
						Hi[II(sat[i], opt)] =  didxi;
						Hi[II(sat[j], opt)] = -didxj;
					}
				}
				/* DD tropospheric delay term */
				if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
					v[nv] -= (tropu[i] - tropu[j]) - (tropr[i] - tropr[j]);
					for (k = 0; k<(opt->tropopt<TROPOPT_ESTG ? 1 : 3); k++) {
						if (!H) continue;
						Hi[IT(0, opt) + k] =  (dtdxu[k + i * 3] - dtdxu[k + j * 3]);
						Hi[IT(1, opt) + k] = -(dtdxr[k + i * 3] - dtdxr[k + j * 3]);
					}
				}
				/* DD phase-bias term */
				if (f<nf) {
					if (opt->ionoopt!=IONOOPT_IFLC) {
						v[nv]-=CLIGHT/freqi*x[IB(sat[i],f,opt)] - CLIGHT/freqj*x[IB(sat[j],f,opt)];
						if (H) {
                        	Hi[IB(sat[i],f,opt)]= CLIGHT/freqi;
                        	Hi[IB(sat[j],f,opt)]=-CLIGHT/freqj;
                    	}
					}
					else {
						v[nv]-=x[IB(sat[i],f,opt)]-x[IB(sat[j],f,opt)];
                    	if (H) {
                        	Hi[IB(sat[i],f,opt)]= 1.0;
                        	Hi[IB(sat[j],f,opt)]=-1.0;
                    	}
					}
				}

				/* SD (single-differenced) measurement error variances */
				Ri[nv] = varerrSD(sat[i], sysi, azel[1 + iu[i] * 2], bl, dt, f, opt);
				Rj[nv] = varerrSD(sat[j], sysj, azel[1 + iu[j] * 2], bl, dt, f, opt);

				trace(2, "sat=%3d-%3d %s%d v=%13.3f R=%8.6f %8.6f\n", sat[i],
					sat[j], f<nf ? "L" : "P", f%nf + 1, v[nv], Ri[nv], Rj[nv]);

				vflg[nv++] = (sat[i] << 16) | (sat[j] << 8) | ((f%nf+1) << 4 ) | (f<nf ? 1 : 2);			
				nb[b]++;
			}
			b++;
		}
	}
	/* baseline length constraint for moving baseline */
    if (opt->mode==PMODE_MOVEB && constbl(rtk, x, NULL, v, H, Ri, Rj, nv)) {
        vflg[nv++] = 3;
        nb[b++]++;
    }
	
	/* DD measurement error covariance */
	ddcov(nb, b, Ri, Rj, nv, R);			

	free(Ri); 		free(Rj); 		free(im);
	free(tropu); 	free(tropr); 	free(dtdxu); 	free(dtdxr);

	return nv;
}
/* index for SD to DD transformation matrix D --------------------------------*/
static int ddidx(rtk_t *rtk, int *ix)
{
	int i, j, k, m, f, nb = 0, na = NR(&rtk->opt), nf = NF(&rtk->opt), nofix;

	for (m = 0; m<4; m++) { /* m=0:GPS,1:GLO,2:GAL,3:BDS */

		nofix = (m == 1);

		for (f = 0, k = na; f<nf; f++, k += MAXSAT) {

			for (i = k; i<k + MAXSAT; i++) {
				if (rtk->x[i] == 0.0 || !test_sys(rtk->ssat[i - k].sys, m) ||
					!rtk->ssat[i - k].vsat[f] || !rtk->ssat[i - k].half[f]) {
					continue;
				}
				if (rtk->ssat[i - k].lock[f] > rtk->opt.minlock && !(rtk->ssat[i - k].slip[f] & 2) &&
					rtk->ssat[i - k].azel[1] >= rtk->opt.elmaskar&&!nofix) {		
					rtk->ssat[i - k].fix[f] = 2; /* fix */
					break;
				}
				else rtk->ssat[i - k].fix[f] = 1;
			}
			for (j = k; j<k + MAXSAT; j++) {
				if (i == j || rtk->x[j] == 0.0 || !test_sys(rtk->ssat[j - k].sys, m) ||
					!rtk->ssat[j - k].vsat[f]) {
					continue;
				}
				if (rtk->ssat[j - k].lock[f] > rtk->opt.minlock && !(rtk->ssat[j - k].slip[f] & 2) &&
					rtk->ssat[i - k].vsat[f] &&
					rtk->ssat[j - k].azel[1] >= rtk->opt.elmaskar&&!nofix) {
					ix[nb * 2] = i;		/* state index of ref bias */
					ix[nb * 2 + 1] = j; /* state index of target bias */
					nb++;
					rtk->ssat[j - k].fix[f] = 2; /* fix */
				}
				else rtk->ssat[j - k].fix[f] = 1;
			}
		}
	}
	return nb;
}
/* resolve integer ambiguity by LAMBDA ---------------------------------------*/
static int resamb_LAMBDA(rtk_t *rtk)
{
	int *iD, *ix, *ii, i, j, m, n, nx, na, nb, fix;
	double *a, *b, *Qaa, *Qbb, *Qba;

	nx = rtk->nx;
	a = mat(nx, 1); 		b = mat(nx, 1);			
	iD = imat(nx, 2);		ix = imat(nx, 1);		ii = imat(nx, 1);

	/* index of SD to DD transformation matrix D */
	na = ddidx(rtk, iD);
	for (i = nb = 0; i < rtk->nx; i++) {
		ii[i] = -1;
		if (rtk->x[i] != 0.0 && rtk->P[i + i*rtk->nx] > 0.0) {
			ii[i] = nb;
			b[nb] = rtk->x[i];
			ix[nb++] = i;	
		}
	}
	
	Qbb = mat(nb, nb);		Qba = mat(nb, na); 		Qaa = mat(na, na); 

	for (i = 0; i<na; i++) {
		a[i] = rtk->x[iD[i * 2]] - rtk->x[iD[i * 2 + 1]];
	}
	for (j = 0; j<nb; j++) for (i = 0; i<na; i++) {
		Qba[j + i*nb] = rtk->P[ix[j] + iD[i*2] * nx] - rtk->P[ix[j] + iD[i*2 + 1] * nx];
	}
	for (j = 0; j<na; j++) {
		m = ii[iD[j*2]]; 
		n = ii[iD[j*2 + 1]];
		for (i = 0; i<na; i++) {
			Qaa[j + i*na] = Qba[m + i*nb] - Qba[n + i*nb];
		}
	}
	for (j = 0; j<nb; j++) for (i = 0; i<nb; i++) {
		Qbb[i + j*nb] = rtk->P[ix[i] + ix[j] * nx];
	}
	/* LAMBDA/MLAMBDA ILS (integer least-square) estimation */
	if ((fix = lambda_PAR(rtk, na, nb, 3, a, b, Qaa, Qbb, Qba))) {
		for (i = 0; i < nb; i++) {
			rtk->xa[ix[i]] = b[i];
			for (j = 0; j < nb; j++) {
				rtk->Pa[ix[j] + ix[i] * nx] = Qbb[j + i*nb];
			}
		}
	}

	free(iD); free(ix); free(ii); free(a); free(b); free(Qaa); free(Qbb); free(Qba); 

	return fix; /* number of ambiguities */
}
/* reject obs by pre-fit residuals */
static int valpre(double *v, double *H, double *R, int *vflg, int nv, int nx, exc_t *exc) 
{
	int i, j, type, stat, nn = 0;
	int *ix = imat(nv, 1);

	for (i = 0; i < nv; i++) ix[i] = i;
	for (i = 0; i < nv; i++) {
		stat = 0;
		type = vflg[i] & 0xF;
		switch (type) {
		case 1: if (fabs(v[i]) > REJION) { exc->flg[exc->n++] = vflg[i]; stat = 1; } break;
		case 2: if (fabs(v[i]) > REJION) { exc->flg[exc->n++] = vflg[i]; stat = 1; } break;
		}
		if (!stat) {
			v[nn]  = v[i];		matcpy(H + nn*nx, H + i*nx, nx, 1);		
			ix[nn] = ix[i];		vflg[nn++] = vflg[i];
		}
	}
	if(nn != nv) {
		for (i = 0; i < nn; i++) {
			for (j = 0; j < nn; j++) { 
				R[j + i*nn] = R[ix[j] + ix[i]*nv];
			}
		}
	}
	
	free(ix);
	return nn;
}
/* validation of solution ----------------------------------------------------*/
static int valpos(rtk_t *rtk, double *v, double *R, int *vflg, int na, double thres, exc_t *exc)
{
	int i, k, type, sat1, sat2, freq, stat = 0;
		
	for (i = 0; i < na; i++) {
		type = vflg[i] & 0xF;
		if (fabs(v[i]) > sqrt(R[i + i*na])*thres) {
			char str[215];
			time2str(rtk->sol.time,str,0);
			sat1  = (vflg[i] >> 16) & 0xFF;
			sat2  = (vflg[i] >> 8) & 0xFF;
			freq = (vflg[i] >> 4) & 0xF;
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
				sat1  = (vflg[i] >> 16) & 0xFF;
				sat2  = (vflg[i] >> 8) & 0xFF;
				freq = (vflg[i] >> 4) & 0xF;
				switch (type) {
				case 1: rtk->ssat[sat2 - 1].resc[freq - 1] = v[i]; 
						rtk->ssat[sat1 - 1].vsat[freq - 1] = 1;
						rtk->ssat[sat2 - 1].vsat[freq - 1] = 1;  break; 
				case 2: rtk->ssat[sat2 - 1].resp[freq - 1] = v[i]; break;
				}
			}
		}
	}
	return stat;
}
/* update solution status ----------------------------------------------------*/
static int update_sol(rtk_t *rtk, int stat)
{
	rtk->sol.stat = stat;

	/* save solution status */
	if (rtk->sol.stat == SOLQ_FIX) {
		for (int i = 0; i<3; i++) {
			rtk->sol.rr[i] = rtk->xa[i];
			rtk->sol.qr[i] = (float)rtk->Pa[i + i*rtk->na];
		}
		rtk->sol.qr[3] = (float)rtk->Pa[1];
		rtk->sol.qr[4] = (float)rtk->Pa[1 + 2 * rtk->na];
		rtk->sol.qr[5] = (float)rtk->Pa[2];
	}
	else if (rtk->sol.stat == SOLQ_FLOAT){
		for (int i = 0; i<3; i++) {
			rtk->sol.rr[i] = rtk->x[i];
			rtk->sol.qr[i] = (float)rtk->P[i + i*rtk->nx];
		}
		rtk->sol.qr[3] = (float)rtk->P[1];
		rtk->sol.qr[4] = (float)rtk->P[1 + 2 * rtk->nx];
		rtk->sol.qr[5] = (float)rtk->P[2];
		rtk->nfix = 0;
	}
}
/* relative positioning ------------------------------------------------------*/
extern void relpos(rtk_t *rtk, const obsd_t *obs, int nu, int nr, const nav_t *nav)
{
	prcopt_t *opt = &rtk->opt;
	gtime_t time = obs[0].time;
	exc_t exc = { 0 };
	double *rs, *dts, *var, *y, *e, *azel, *freq, *v, *H, *R, *xp, *Pp, dt;
	int n, nf, ns, ny, nv, sat[MAXSAT], iu[MAXSAT], ir[MAXSAT];
	int i, j, info, *vflg, *svh, stat = SOLQ_NONE;
	
	trace(2, "relpos  : nx=%d nu=%d nr=%d\n", rtk->nx, nu, nr);

	dt = timediff(time, obs[nu].time);

	n = nu + nr; 		nf = NF(opt);
	rs = mat(6, n); 	dts = mat(2, n); 	var  = mat(1, n);	svh  = imat(1, n); 	
	y = mat(nf * 2, n); e = mat(3, n);		azel = zeros(2, n); freq = zeros(nf, n);

	for (i = 0; i<MAXSAT; i++) {
		rtk->ssat[i].sys = satsys(i + 1, NULL);
		for (j = 0; j<nf; j++) { rtk->ssat[i].vsat[j] = rtk->ssat[i].fix[j] = 0;
								 rtk->ssat[i].slip[j] &= 0xFC; }
	}

	/* satellite positions/clocks */
	satposs(time, obs, n, nav, opt->sateph, rs, dts, var, svh);

	/* UD (undifferenced) residuals for base station */
	zdres(1, obs + nu, nr, rs + nu * 6, dts + nu * 2, var + nu, svh + nu, nav, rtk->rb, opt, 1,
		y + nu*nf * 2, e + nu * 3, azel + nu * 2, freq + nu*nf);
	
	/* time-interpolation of residuals (for post-processing) */
	if (opt->intpref) {
		dt = intpres(time, obs + nu, nr, nav, rtk, y + nu*nf * 2);		
	}

	/* select common satellites between rover and base-station */
	if ((ns = selsat(obs, azel, nu, nr, opt, sat, iu, ir)) <= 0) {				
		free(rs); free(dts); free(var); free(y); free(e); free(azel); free(freq);
		return ;
	}

	/* LC observation */
	lc_rel(rtk, obs, sat, iu, ir, ns, nav);			

	/* cycle jump detect  */
	cycJump_detect(rtk, obs, sat, iu, ir, ns, nav);						

	/* temporal update of states */
	udstate(rtk, obs, sat, iu, ir, ns, nav);
	
	ny = ns*nf * 2 + 2;
	xp = mat(rtk->nx, 1); 	Pp = zeros(rtk->nx, rtk->nx); 	vflg = imat(ny, 1);
	v = mat(ny, 1); 		H = zeros(rtk->nx, ny); 		R = mat(ny, ny); 

	for (i = 0; i < MAX_ITER; i++) {	
		matcpy(xp, rtk->x, rtk->nx, 1);
		matcpy(Pp, rtk->P, rtk->nx, rtk->nx);

		/* reject obs by pre-fit residuals */
		zdres(0, obs, nu, rs, dts, var, svh, nav, xp, opt, 0, y, e, azel, freq);					
		nv = ddres(rtk, nav, dt, xp, sat, y, e, azel, freq, iu, ir, ns, v, H, R, vflg, &exc);	
		nv = valpre(v, H, R, vflg, nv, rtk->nx, &exc);

		/* Kalman filter measurement update */
		if ((info = filter(xp, Pp, H, v, R, rtk->nx, nv))) {						
			trace(2, "filter error (info=%d)\n", info);
			break;
		}

		/* reject obs by pos-fit residuals */	
		zdres(0, obs, nu, rs, dts, var, svh, nav, xp, opt, 0, y, e, azel, freq);				
		nv = ddres(rtk, nav, dt, xp, sat, y, e, azel, freq, iu, ir, ns, v, NULL, R, vflg, &exc);		
		if (!valpos(rtk, v, R, vflg, nv, 4.0, &exc)) {				
			stat = SOLQ_FLOAT;
			break;
		}
	}

	/* update satellite && ambiguity information */
	update_sat_amb(rtk, &stat, obs, sat, iu, ir, ns);

	/* resolve integer ambiguity by LAMBDA */
	if (stat == SOLQ_FLOAT) {
		matcpy(rtk->x, xp, rtk->nx, 1);
		matcpy(rtk->P, Pp, rtk->nx, rtk->nx);

		if(resamb_LAMBDA(rtk) > 4) {
			/* validation of fixed solution */
			zdres(0, obs, nu, rs, dts, var, svh, nav, rtk->xa, opt, 0, y, e, azel, freq);
			nv = ddres(rtk, nav, dt, rtk->xa, sat, y, e, azel, freq, iu, ir, ns, v, NULL, R, vflg, &exc);
			if (!valpos(rtk, v, R, vflg, nv, 4.0, &exc)) {				
				stat = SOLQ_FIX;
			}
		}
	}

	 /* update solution status */
	update_sol(rtk, stat);
	
	free(rs); free(dts); free(var); free(svh);	free(y); free(e); free(azel); free(freq);
	free(xp); free(Pp);  free(v);   free(H); 	free(R); free(vflg);

}
