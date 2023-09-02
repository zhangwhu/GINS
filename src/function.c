
#include "rtklib.h"
#include "ppp_state.h"

#define SQR(x)      ((x)*(x))
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define ERR_ION     5.0         /* ionospheric delay std (m) */
#define ERR_TROP    3.0         /* tropspheric delay std (m) */
#define ERR_BRDCI   0.5         /* broadcast iono model error factor */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */

#define SWAP_I(x,y)     do {int _z=x; x=y; y=_z;} while (0)
#define SWAP_D(x,y)     do {double _z=x; x=y; y=_z;} while (0)
/* select median */
extern double median(const double *data, int n) 
{
	if (n < 1) return 0;

	int i, k, j, *ind = imat(n, 1);
	double max, mid, *D = mat(n, 1);

	for (i = 0; i<n; i++) { ind[i] = i; D[i] = data[i]; }
	for (i = 0; i<n-1; i++) {
		max = D[i];
		for(k = i+1; k<n; k++) {
			if(D[k] > max) {
				j = k; 		max = D[k];	
			}
		}
		if (max > D[i]) {
			SWAP_D(D[i], D[j]);		SWAP_I(ind[i], ind[j]);
		}
	}

	mid = n % 2 ? D[n/2] : (D[n/2] + D[n/2-1])/2;
	free(ind);		free(D);
	return mid;
}
/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
	if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
	return atan2(-tan(beta), sin(mu)) + PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
static int yaw_angle(int sat, const char *type, int opt, double beta, double mu,
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
/* initialize RTK control ------------------------------------------------------
* initialize RTK control struct
* args   : rtk_t    *rtk    IO  TKk control/result struct
*          prcopt_t *opt    I   positioning options (see rtklib.h)
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkinit(rtk_t *rtk, const prcopt_t *opt)
{
	sol_t sol0 = { { 0 } };
	ssat_t ssat0 = { 0 };
	int i;

	trace(3, "rtkinit :\n");

	rtk->sol = sol0;
	for (i = 0; i<6; i++) rtk->rb[i] = 0.0;
	rtk->nx = opt->mode <= PMODE_FIXED ? relnx(opt) : pppnx(opt);
	rtk->na = opt->mode <= PMODE_FIXED ? relnx(opt) : pppnx(opt);
	rtk->tt = opt->ti;
	rtk->x = zeros(rtk->nx, 1);
	rtk->P = zeros(rtk->nx, rtk->nx);
	rtk->xa = zeros(rtk->na, 1);
	rtk->Pa = zeros(rtk->na, rtk->na);
	rtk->nfix = rtk->neb = 0;
	for (i = 0; i<MAXSAT; i++) {
		rtk->ssat[i] = ssat0;
		rtk->ssat[i].sys = satsys(i + 1, NULL);
	}
	for (i = 0; i<MAXERRMSG; i++) rtk->errbuf[i] = 0;
	rtk->reset = 0;
	rtk->opt = *opt;
}

/* free rtk control ------------------------------------------------------------
* free memory for rtk control struct
* args   : rtk_t    *rtk    IO  rtk control/result struct
* return : none
*-----------------------------------------------------------------------------*/
extern void rtkfree(rtk_t *rtk)
{
	trace(3, "rtkfree :\n");

	rtk->nx = rtk->na = 0;
	free(rtk->x); rtk->x = NULL;
	free(rtk->P); rtk->P = NULL;
	free(rtk->xa); rtk->xa = NULL;
	free(rtk->Pa); rtk->Pa = NULL;
}

/* bd-2 satellite code multipath--------------------------------------------*/
extern double bd2smp(int orb, double *azel, int nq) {
	int i, el;
	double v[3][2], mp = 0.0;

	el = (int)((azel[1] * R2D) / 10) * 10;
	for (i = 0; i < 20; i++) {
		if (bdsmptable[i][0] == orb && bdsmptable[i][1] == el) {
			v[nq][0] = bdsmptable[i][nq + 2];
			v[nq][1] = bdsmptable[i + 1][nq + 2];
			mp = (azel[1] * R2D - el) / 10 * v[nq][1] + ((el + 10) - azel[1] * R2D) / 10 * v[nq][0];
			break;
		}
	}
	//trace(1, "prn=%4d freq=%4d azel=%10.4f el=%4d %6.3f %6.3f %6.3f\n", prn, nq, azel[1] * R2D, el, mp, v[nq][0], v[nq][1]);
	return mp;
}
/* set antenna parameters ----------------------------------------------------*/
extern void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs, 
	const pcvs_t *pcvr, const sta_t *sta)
{
	int i, j, mode=PMODE_DGPS<=popt->mode&&popt->mode<=PMODE_FIXED;
	char id[64];
	double pos[3], del[3];
	pcv_t *pcv, pcv0={ 0 };

	/* set satellite antenna parameters */
	for (i = 0; i<MAXSAT; i++) {
		if (!(satsys(i + 1, NULL)&popt->navsys)) continue;
		if (!(pcv = searchpcv(i + 1, "", time, pcvs))) {
			satno2id(i + 1, id);
			trace(3, "no satellite antenna pcv: %s\n", id);
			continue;
		}
		nav->pcvs[i] = *pcv;
	}

	for (i = 0; i<(mode?2:1); i++) {
		// PCO
		if (sta[i].deltype == 1) { /* xyz */
			if (norm(sta[i].pos,3) > 0.0) {
				ecef2pos(sta[i].pos, pos);
				ecef2enu(pos,sta[i].del, del);
				for (j=0;j<3;j++) popt->antdel[i][j] = del[j];
			}
		}
		else { /* enu */
			for (j = 0; j<3; j++) popt->antdel[i][j] = sta[i].del[j];
		}

		// PCV
		popt->pcvr[i] = pcv0;
		strcpy(popt->anttype[i], sta[i].antdes);
		if (!(pcv = searchpcv(0, popt->anttype[i], time, pcvr))) {
			trace(1,"no receiver antenna pcv: %s\n",popt->anttype[i]);
			*popt->anttype[i]='\0';
			continue;
		}
		trace(1, "%s >> %s\n", popt->anttype[i], pcv->type);
		strcpy(popt->anttype[i], pcv->type);
		popt->pcvr[i]=*pcv;
	}
}
/* phase windup model -------------------------------------------------------- */
extern int model_phw(gtime_t time, int sat, const char *type, int opt,
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
		ds[i] = exs[i] - ek[i] * dot(ek, exs, 3) - eks[i];			
		dr[i] = exr[i] - ek[i] * dot(ek, exr, 3) + ekr[i];			
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
/* troposphere correction ----------------------------------------------------- */
extern int tropcorr(gtime_t time, const double *pos, const double *azel, const prcopt_t *opt, 
	const double *x, double *dtdx, const nav_t *nav, double *dtrp, double *var)
{
	const double zazel[] = { 0.0, PI / 2.0 };
	double trp[3] = { 0 };
	double zhd, m_h, m_w, cotz, grad_n, grad_e;

	if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
		matcpy(trp, x + IT(opt), opt->tropopt == TROPOPT_EST ? 1 : 3, 1);

		/* zenith hydrostatic delay */
		zhd = tropmodel(time, pos, zazel, 1, REL_HUMI, 1);

		/* mapping function */
		m_h = tropmapf(time, pos, azel, 1, &m_w);					

		/* gradient */
		if (opt->tropopt == TROPOPT_ESTG && azel[1]>0.0) {
			cotz = 1.0 / tan(azel[1]);
			grad_n = m_w*cotz*cos(azel[0]);
			grad_e = m_w*cotz*sin(azel[0]);
			m_w += grad_n*trp[1] + grad_e*trp[2];
			dtdx[1] = grad_n*trp[0];								
			dtdx[2] = grad_e*trp[0];
		}
		dtdx[0] = m_w;
		*dtrp = m_h*zhd + m_w*trp[0];
		*var = 0.0;
		return 1;
	}
	/* troposphere ZTD correction*/
	if (opt->tropopt == TROPOPT_ZTD) {									

	}
	/* saastamoinen model */
	if (opt->tropopt == TROPOPT_SAAS) {
		*dtrp = tropmodel(time, pos, azel, 1, REL_HUMI, 0);
		*var = SQR(ERR_SAAS / (sin(azel[1]) + 0.1));
		return 1;
	}
	/* no correction */
	if (opt->tropopt == TROPOPT_OFF) {				
		*dtrp = 0.0;
		*var = SQR(ERR_TROP);
		return 1;
	}
	
	return 0;
}
/* ionospheric correction ----------------------------------------------------- */
extern int ionocorr(gtime_t time, const double *pos, const double *azel, const prcopt_t *opt, int sat, const double *x,
	const nav_t *nav, double *dion, double *var)
{
	if (opt->ionoopt == IONOOPT_TEC) {
		return iontec(time, nav, pos, azel, 1, dion, var);				
	}
	if (opt->ionoopt == IONOOPT_BRDC) {
		*dion = ionmodel(time, nav->ion_gps, pos, azel);
		*var = SQR(*dion*ERR_BRDCI);
		return 1;
	}
	if (opt->ionoopt == IONOOPT_EST) {									
		*dion = x[II(sat, opt)];
		*var = 0.0;
		return 1;
	}
	if (opt->ionoopt == IONOOPT_IFLC) {
		*dion = *var = 0.0;
		return 1;
	}
	if (opt->ionoopt == IONOOPT_STEC) {

	}
	return 0;
}
/* observation correction ----------------------------------------------------- */
extern void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel, const prcopt_t *opt, const double *dantr,
	const double *dants, double phw, const double php, double *L, double *P, double *Lc, double *Pc)
{
	int i, k, sys, prn, orb, ix, nf = opt->nf;
	double freq[NFREQ], c1, c2, C1, C2, gamma;
	ifcb_t *ifcb = NULL;

	for (i = 0; i < NFREQ; i++) {
		freq[i] = sat2freq(obs->sat, obs->code[i], nav);
		L[i] = P[i] = 0.0;
	}
	sys = satsys(obs->sat, &prn);

	for (i = 0; i < nav->ni - 1; i++) {
		if (timediff(obs[0].time, nav->ifcb[i].ti) >= 0 && timediff(obs[0].time, nav->ifcb[i+1].ti) < 0) {
			ifcb = &nav->ifcb[i];
			continue;
		}
	}

	for (i = 0; i< nf; i++) {
		if (testsnr(0, i, azel[1], obs->SNR[i] * SNR_UNIT, &opt->snrmask)) continue;

		/* antenna phase center and phase windup correction */
		if (freq[i] != 0.0 && obs->L[i] != 0.0) {
			L[i] = CLIGHT * obs->L[i] / freq[i] - dants[i] - dantr[i] - CLIGHT * phw / freq[i] + php;
		}
		if (freq[i] != 0.0 && obs->L[i] != 0.0 && freq[0] != 0.0 && i == 2 && ifcb) {
			if (fabs(ifcb->bias[obs->sat - 1]) < 10) {
				L[i] -= (SQR(freq[0]/freq[2]) - 1) * ifcb->bias[obs->sat - 1];
			}
			// printf("IFCB %.6f\n", (SQR(freq[0]/freq[2]) - 1) * ifcb->bias[obs->sat - 1]);
		}
		
		if (obs->P[i] != 0.0) {
			P[i] = obs->P[i] - dants[i] - dantr[i];

			/* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
			if (sys == SYS_GPS || sys == SYS_GLO) {
				if (obs->code[i] == CODE_L1C) P[i] += nav->cbias[obs->sat - 1][1];
				if (obs->code[i] == CODE_L2C) P[i] += nav->cbias[obs->sat - 1][2];
			}
			else if (sys == SYS_CMP) {
				if (prn < 19) {
					satorb(obs->sat, &orb);
					P[i] += bd2smp(orb, azel, i);
				}
			}
			else if (sys == SYS_GAL) {
				if (obs->code[i] == CODE_L1X) P[i] += nav->cbias[obs->sat - 1][1];
				if (obs->code[i] == CODE_L5X) P[i] += nav->cbias[obs->sat - 1][2];
			}
		}
	}

	/* iono-free LC */				
	*Lc = *Pc = 0.0;
	if (freq[0] == 0.0 || freq[1] == 0.0) return;

	C1 =  SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));
	C2 = -SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));

	if (L[0] != 0.0&&L[1] != 0.0) *Lc = C1*L[0] + C2*L[1];
	if (P[0] != 0.0&&P[1] != 0.0) *Pc = C1*P[0] + C2*P[1];
}

/* measurement error variance ------------------------------------------------*/
extern double varerr(unsigned char sat, double el, int freq, int type, const prcopt_t *opt)
{
	int sys, prn, orb;
	double fact = 1.0, sinel, a, b;				

	sys = satsys(sat, &prn);	satorb(sat, &orb);
	el = (el < MIN_EL) ? MIN_EL : el; 	sinel = sin(el);

	fact = (type == 1) ? opt->eratio[freq] : 1;

	switch (sys) {
		case SYS_GPS: fact *= (freq == 2 /* && type == 1 */) ? EFACT_GPS_L5 : 1; break;		// GPS的第三频点设置不太多
		case SYS_GLO: fact *= EFACT_GPS_L5; break;
		case SYS_GAL: fact *= 1; break;
		case SYS_CMP: fact *= (prn < 19) ? 2 : 1;  break;
		default: fact *= 100;
	}
	a = fact*opt->err[1];
	b = fact*opt->err[2];

	return (opt->ionoopt == IONOOPT_IFLC ? 9.0 : 1.0) * (a*a + b*b / sinel / sinel);
}

/* single-differenced measurement error variance -----------------------------*/
extern double varerrSD(int sat, int sys, double el, double bl, double dt, int f,
	const prcopt_t *opt)
{
	double a, b, c = opt->err[3] * bl / 1E4, d = CLIGHT*opt->sclkstab*dt, fact = 1.0;
	double sinel = sin(el);
	int nf = NF(opt);

	if (f >= nf) fact = opt->eratio[f - nf];
	if (fact <= 0.0) fact = opt->eratio[0];
	fact *= sys == SYS_GLO ? EFACT_GLO : (sys == SYS_SBS ? EFACT_SBS : EFACT_GPS);
	a = fact*opt->err[1];
	b = fact*opt->err[2];

	return 2.0*(opt->ionoopt == IONOOPT_IFLC ? 3.0 : 1.0)*(a*a + b*b / sinel / sinel + c*c) + d*d;
}
