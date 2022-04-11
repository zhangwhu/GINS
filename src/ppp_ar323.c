/*------------------------------------------------------------------------------
* ppp_ar.c : ppp ambiguity resolution
*
* options : -DREV_WL_FCB reversed polarity of WL FCB
*
* reference :
*    [1] H.Okumura, C-gengo niyoru saishin algorithm jiten (in Japanese),
*        Software Technology, 1991
*
*          Copyright (C) 2012-2013 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2014-08-20 12:56:44 $
* history : 2013/03/11 1.0  new
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#include "ppp_state.h"

/* constants/macros ----------------------------------------------------------*/
#define SQR(x)          ((x)*(x))
#define ROUND(x)        (int)floor((x)+0.5)
#define FROUND(x)		((x)-(int)floor((x)+0.5))
#define NROUND(x)		(x-(int)floor(x))

#define SWAP_I(x,y)     do {int _z=x; x=y; y=_z;} while (0)
#define SWAP_D(x,y)     do {double _z=x; x=y; y=_z;} while (0)
#define MIN(x,y)        ((x)<(y)?(x):(y))

#define CONST_AMB       0.0001       /* constraint to fixed ambiguity */
#define CONST_FCB		0.0008		/* constraint to FCB (cycle) */
#define CONST_AFIF      0.3         /* constraint to AFIF (m) */
#define LOG_PI          1.14472988584940017 /* log(pi) */
#define SQRT2           1.41421356237309510 /* sqrt(2) */
#define MAXCAMB			10000
#define MAXVAR			1.2


static FILE *fp_stat = NULL;;
/* open solution status file--------------------------------------------------*/
extern int ppparopenstat(const char *file)
{
	if (!(fp_stat = fopen(file, "w"))) {
		trace(1, "pppopenstat: file open error path=%s\n", file);
		return 0;
	}
	return 1;
}
/* close solution status file --------------------------------------------------*/
extern void ppparclosestat(void)
{
	trace(3, "rtkclosestat:\n");

	if (fp_stat) fclose(fp_stat);
	fp_stat = NULL;
}
/* linear dependency check ---------------------------------------------------*/
static int is_depend(int sat1, int sat2, int *flgs, int *max_flg)
{
	int i;

	if (flgs[sat1 - 1] == 0 && flgs[sat2 - 1] == 0) {
		flgs[sat1 - 1] = flgs[sat2 - 1] = ++(*max_flg);
	}
	else if (flgs[sat1 - 1] == 0 && flgs[sat2 - 1] != 0) {
		flgs[sat1 - 1] = flgs[sat2 - 1];
	}
	else if (flgs[sat1 - 1] != 0 && flgs[sat2 - 1] == 0) {
		flgs[sat2 - 1] = flgs[sat1 - 1];
	}
	else if (flgs[sat1 - 1]>flgs[sat2 - 1]) {
		for (i = 0; i<MAXSAT; i++) if (flgs[i] == flgs[sat2 - 1]) flgs[i] = flgs[sat1 - 1];
	}
	else if (flgs[sat1 - 1]<flgs[sat2 - 1]) {
		for (i = 0; i<MAXSAT; i++) if (flgs[i] == flgs[sat1 - 1]) flgs[i] = flgs[sat2 - 1];
	}
	else return 0; /* linear depenent */
	return 1;
}
static int ranking(double *Q, int *sat1, int *sat2, double *a, int *NW, int na, int *ia) {
	int i, j, m, ind, flgs[MAXSAT] = { 0 }, max_flg = 0;
	double *q, *T, *F, min;

	q = mat(na, 1);		T = mat(na, na);		F = mat(na, na);
	for (i = 0; i < na; i++) { q[i] = Q[i + i*na];	 ia[i] = i; }
	for (i = 0; i < na - 1; i++) {
		ind = i;		min = q[i];
		for (j = i + 1; j < na; j++) {
			if (q[j] < min) {
				ind = j;   min = q[j];
			}
		}
		if (ind != i) {
			SWAP_I(ia[i], ia[ind]);
			SWAP_I(NW[i], NW[ind]);
			SWAP_I(sat1[i], sat1[ind]);
			SWAP_I(sat2[i], sat2[ind]);
			SWAP_D(a[i], a[ind]);
			SWAP_D(q[i], q[ind]);
		}
	}
	for (i = 0, m = 0; i < na; i++) {		
		if (!is_depend(sat1[i], sat2[i], flgs, &max_flg)) continue;
		a[m]  = a[i];
		NW[m] = NW[i];
		ia[m++] = ia[i];	
	}
	for (i = 0; i < m; i++) {
		for (j = 0; j < na; j++) {
			T[j + i * na] = (j == ia[i]) ? 1 : 0;
		}
	}
	matmul("TN", m, na, na, 1.0, T, Q, 0.0, F);
	matmul("NN", m, m, na, 1.0, F, T, 0.0, Q);

	free(q);	free(T);	free(F);
	return m;
}
/* select fcb data struct   --------------------------------------------------*/
static fcbd_t* selfcb(const nav_t *nav, gtime_t tt)
{
	for (int i = 0; i < nav->nf; i++) {
		if ((timediff(nav->fcb[i].ts, tt) <= 0) && (timediff(nav->fcb[i].te, tt) > 0)) {
			return &nav->fcb[i];
		}
	}
	return NULL;
}
static int sel_sat_SD(const rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, 
	const int *refs,int *isat1, int *isat2)			
{
	int i, j, k, m = 0, sys, sat, stat, isat[4][20] = { 0 }, sn[4] = { 0 }, iref[4];

	for (k = 0; k < 4; k++) {
		stat = -1;
		for (i = 0; i < n; i++) {
			if (refs[k] == obs[i].sat) {
				stat = i;  break;
			}
		}
		iref[k] = stat > 0 ? stat : -1;
	}
	for (i = 0; i < n; i++) {
		sat = obs[i].sat;
		if (!rtk->ssat[sat - 1].vsat[0] || !rtk->ssat[sat - 1].vsat[1] || rtk->ssat[sat - 1].azel[1] < rtk->opt.elmin) continue;
		if (!(sys = satsys(sat, NULL))) continue;	

		switch (sys){
		case SYS_GPS: isat[0][sn[0]++] = i; break;
		case SYS_GLO: isat[1][sn[1]++] = i; break;
		case SYS_GAL: isat[2][sn[2]++] = i; break;
		case SYS_CMP: isat[3][sn[3]++] = i; break;
		}
	}
	/* method 1 */
	for (k = 0; k < 4; k++) {
		if (sn[k] < 2 || iref[k] < 0) continue;
		for (i = 0; i < sn[k]; i++) {
			if (iref[k] != isat[k][i]) {
				isat1[m] = iref[k];
				isat2[m++] = isat[k][i];
			}
		}
	}

	/* method 2 */
	/*for (k = 0; k < 4; k++) {
		if (sn[k] < 2) continue;
		for (i = 0; i < sn[k] - 1; i++) {
			for (j = i + 1; j < sn[k]; j++) {
				sat1[m] = ssat[k][i];		sat2[m++] = ssat[k][j];
			}
		}
	}*/

#if 0
	for (i = 0; i < m; i++) {
		char id1[32], id2[32];
		
		satno2id(obs[isat1[i]].sat, id1);
		satno2id(obs[isat2[i]].sat, id2);
		trace(1, "%s--%s\n", id1, id2);
	}
#endif
	return m;
}
static void fix_EWL(const rtk_t *rtk, const obsd_t *obs, const nav_t *nav, const int *isat1, const int *isat2, int n, int *NE)				
{
	int i, j, k, f, sat1, sat2, info, nv = 0, stat = 1;
	double LC, BE, vE, freq[NFREQ], lam2, lam3, el, *xp, *Pp, *v, *H, *R, *var;
	const prcopt_t *opt = &rtk->opt;

	xp = mat(rtk->nx, 1); Pp = mat(rtk->nx, rtk->nx);
	v = mat(n, 1);		  H = zeros(rtk->nx, n);		R = mat(n, n);		var = mat(n, 1);

	matcpy(xp, rtk->xa, rtk->nx, 1);
	matcpy(Pp, rtk->Pa, rtk->nx, rtk->nx);
	for (i = 0; i < n; i++) {
		sat1 = obs[isat1[i]].sat;
		sat2 = obs[isat2[i]].sat;

		for (f = 0; f < NFREQ; f++) freq[f] = sat2freq(sat1, obs[isat1[i]].code[f], nav);
		lam2 = CLIGHT / freq[1];		
		lam3 = CLIGHT / freq[2];

		if (nav->elbias[sat1 - 1] > 999 || nav->elbias[sat2 - 1] > 999) continue;
		if (rtk->ssat[sat1 - 1].azel[1] < 15 * D2R || !rtk->ssat[sat1 - 1].vsat[1] || !rtk->ssat[sat1 - 1].vsat[2]) continue;
		if (rtk->ssat[sat2 - 1].azel[1] < 15 * D2R || !rtk->ssat[sat2 - 1].vsat[1] || !rtk->ssat[sat2 - 1].vsat[2]) continue;

		/* wide-lane ambiguity */
		j = IB(sat1, 1, opt);
		k = IB(sat2, 1, opt);
		LC = (xp[j] - xp[k]) / lam2 - (xp[j + MAXSAT] - xp[k + MAXSAT]) / lam3;
		el = (nav->elbias[sat1 - 1] - nav->elbias[sat2 - 1]);	//ÐÞ¸Ä

		/* validation of integer wide-lane ambigyity */
		BE = LC - el;
		if (fabs(FROUND(BE)) <= opt->thresar[2]) {			
			NE[i] = ROUND(BE);	

			v[nv] = (NE[i] + el) - LC;
			H[j + nv*rtk->nx] = 1.0 / lam2;
			H[k + nv*rtk->nx] = -1.0 / lam2;
			H[j + MAXSAT + nv*rtk->nx] = -1.0 / lam3;
			H[k + MAXSAT + nv*rtk->nx] = 1.0 / lam3;
			var[nv++] = SQR(0.01);
#if 0
			char id1[32], id2[32];
			satno2id(obs[isat1[i]].sat, id1);
			satno2id(obs[isat2[i]].sat, id2);
			trace(1, "EL-KF:%s--%s %10.4f %10.4f\n", id1, id2, LC, NE[i] + el);
#endif
		}
	}

	for (j = 0; j < nv; j++) for (k = 0; k < nv; k++) {
		R[j + k*nv] = j == k ? var[j] : 0.0;
	}

	/* update states with constraints */
	/*if (nv > 0 && (info = filter(xp, Pp, H, v, R, rtk->nx, nv))) {
		trace(1, "filter error (info=%d)\n", info);
		stat = 0;
	}*/

	if (stat) {
		matcpy(rtk->xa, xp, rtk->nx, 1);
		matcpy(rtk->Pa, Pp, rtk->nx, rtk->nx);
	}
	free(xp); free(Pp);  free(v);  free(H); free(R); free(var);
}
static void fix_WL(const rtk_t *rtk, const obsd_t *obs, const nav_t *nav, const int *isat1, const int *isat2, int n, int *NW, int *NE)							
{
	int i, j, k, f, sat1, sat2, info, nv = 0, stat = 1;
	double LC, BW, vW, freq[NFREQ], lam1, lam2, wl, p1, p2, *xp, *Pp, *v, *H, *R, *var;
	const prcopt_t *opt = &rtk->opt;

	xp = mat(rtk->nx, 1); Pp = mat(rtk->nx, rtk->nx);
	v = mat(n, 1);		  H = zeros(rtk->nx, n);		R = mat(n, n);		var = mat(n, 1);

	matcpy(xp, rtk->xa, rtk->nx, 1);
	matcpy(Pp, rtk->Pa, rtk->nx, rtk->nx);
	for (i = 0; i < n; i++) {
		sat1 = obs[isat1[i]].sat;
		sat2 = obs[isat2[i]].sat;

		for (f = 0; f < NFREQ; f++) freq[f] = sat2freq(sat1, obs[isat1[i]].code[f], nav);
		lam1 = CLIGHT / freq[0];
		lam2 = CLIGHT / freq[1];

		if (nav->wlbias[sat1 - 1] > 999 || nav->wlbias[sat2 - 1] > 999) continue;
		if (rtk->ssat[sat1 - 1].azel[1] < 15 * D2R || !rtk->ssat[sat1 - 1].vsat[0] || !rtk->ssat[sat1 - 1].vsat[1]) continue;
		if (rtk->ssat[sat2 - 1].azel[1] < 15 * D2R || !rtk->ssat[sat2 - 1].vsat[0] || !rtk->ssat[sat2 - 1].vsat[1]) continue;
		
		/* wide-lane ambiguity */
		j = IB(sat1, 0, opt);
		k = IB(sat2, 0, opt);
		LC = (xp[j] - xp[k]) / lam1 - (xp[j + MAXSAT] - xp[k + MAXSAT]) / lam2;
		wl = (nav->wlbias[sat1 - 1] - nav->wlbias[sat2 - 1]);		

		/* validation of integer wide-lane ambigyity */
		if (opt->modear == ARMODE_PPPAR) BW = LC - wl;
		else if (opt->modear == ARMODE_PPPAR_ILS) BW = LC + wl;			//×¢Òâ£¡£¡£¡

		if (fabs(FROUND(BW)) <= opt->thresar[2]) {			
			NW[i] = ROUND(BW);

			v[nv] = (NW[i] + wl) - LC;
			H[j + nv*rtk->nx] = 1.0 / lam1;
			H[k + nv*rtk->nx] = -1.0 / lam1;
			H[j + MAXSAT + nv*rtk->nx] = -1.0 / lam2;
			H[k + MAXSAT + nv*rtk->nx] = 1.0 / lam2;
			var[nv++] = SQR(0.02);
#if 0
			char id1[32], id2[32];
			satno2id(obs[isat1[i]].sat, id1);
			satno2id(obs[isat2[i]].sat, id2);
			trace(1, "WL-KF:%s--%s %10.4f %10.4f %10.4f\n", id1, id2, LC, NW[i] + wl, v[nv-1]);
#endif

		} 
	}
	for (j = 0; j < nv; j++) for (k = 0; k < nv; k++) {
		R[j + k*nv] = j == k ? var[j] : 0.0;
	}

	/* update states with constraints */
	/*if (nv>0 && (info = filter(xp, Pp, H, v, R, rtk->nx, nv))) {
		trace(1, "filter error (info=%d)\n", info);
		stat = 0;
	}*/

	if (stat) {
		matcpy(rtk->xa, xp, rtk->nx, 1);
		matcpy(rtk->Pa, Pp, rtk->nx, rtk->nx);
	}
	free(xp); free(Pp); free(v);  free(H); free(R); free(var);
}
static int fix_NL(rtk_t *rtk, const obsd_t *obs, const nav_t *nav, int *isat1, int *isat2, int n, int *NW)
{
	int i, j, k, f, sat1, sat2, na, naa, nb, *ix, *ia, fix = 0;
	double *a, *b, *D, *E, *Qaa, *Qba, *Qbb, *Z;
	double LC, nl, lam_NL, lam1, lam2, C1, C2, freq[NFREQ] = { 0.0 };
	fcbd_t *fcb;

	if (!(fcb = selfcb(nav, rtk->sol.time))) return 0;

	a = zeros(n, 1);		b = mat(rtk->nx, 1);	D = zeros(rtk->nx, n);		ix = imat(rtk->nx, 1);	

	for (i = na = 0; i<n; i++) {
		sat1 = obs[isat1[i]].sat;
		sat2 = obs[isat2[i]].sat;

		if (NW[i] > MAXCAMB - 1 || fcb->bias[sat1 - 1][0] > 100 || fcb->bias[sat2 - 1][0] > 100) continue;

		/* float narrow-lane ambiguity (cycle) */
		j = IB(sat1, 0, &rtk->opt);
		k = IB(sat2, 0, &rtk->opt);
		for (f = 0; f < NFREQ; f++) freq[f] = sat2freq(sat1, obs[isat1[i]].code[f], nav);
		lam1 = CLIGHT / freq[0];
		lam2 = CLIGHT / freq[1];
		C1 = freq[0] / (freq[0] - freq[1]);
		C2 = freq[1] / (freq[0] - freq[1]);
		LC = C1*(rtk->xa[j] - rtk->xa[k]) / lam1 - C2*(rtk->xa[j + MAXSAT] - rtk->xa[k + MAXSAT]) / lam2 - C2*NW[i];	
		nl = fcb->bias[sat1 - 1][0] - fcb->bias[sat2 - 1][0];
		a[na] = LC - nl;

		/* validation of narrow-lane ambiguity */
		if (fabs(FROUND(a[na])) < 0.2/*rtk->opt.thresar[2]*/) {
			isat1[na] = isat1[i];
			isat2[na] = isat2[i];
			D[j + na*rtk->nx] =  C1*1.0 / lam1;
			D[k + na*rtk->nx] = -C1*1.0 / lam1;
			D[j + MAXSAT + na*rtk->nx] = -C2*1.0 / lam2;
			D[k + MAXSAT + na*rtk->nx] =  C2*1.0 / lam2;
			NW[na++] = NW[i];
#if 0
			char id1[32], id2[32];
			satno2id(obs[isat1[i]].sat, id1);
			satno2id(obs[isat2[i]].sat, id2);
			trace(1, "NL-MM:%s--%s %6.3f\n", id1, id2, FROUND(B1[m-1]));
#endif
		}
	}

	if (na >= 4) {
		for (i = nb = 0; i < rtk->nx; i++) {
			if (rtk->xa[i] != 0.0&&rtk->Pa[i + i*rtk->nx] > 0.0) {
				b[nb] = rtk->xa[i];
				ix[nb++] = i;
			}
		}
		E = mat(rtk->nx, na);	    ia = imat(na, 1);	
		Qaa = mat(na, na);			Qba = mat(nb, na);		Qbb = mat(nb, nb);

		/* covariance of narrow-lane ambiguities */
		matmul("NN", rtk->nx, na, rtk->nx, 1.0, rtk->Pa, D, 0.0, E);
		matmul("TN", na, na, rtk->nx, 1.0, D, E, 0.0, Qaa);

		naa = ranking(Qaa, isat1, isat2, a, NW, na, ia);

		for (i = 0; i < nb; i++) {
			for (k = 0; k < nb; k++) {
				Qbb[k + i*nb] = rtk->Pa[ix[k] + ix[i] * rtk->nx];
			}
			for (k = 0; k < naa; k++) {
				Qba[i + k*nb] = E[ix[i] + ia[k]*rtk->nx];
			}
		}

		if ((fix = lambda_PAR(rtk, naa, nb, 3, a, b, Qaa, Qbb, Qba))) {
			for (i = 0; i < nb; i++) {
				rtk->xa[ix[i]] = b[i];
				for (j = 0; j < nb; j++) {
					rtk->Pa[ix[j] + ix[i] * rtk->nx] = Qbb[j + i*nb];
				}
			}
		}

		free(E);	free(ia);	 free(Qaa);		free(Qba);	 free(Qbb);
	}
	free(a);	free(b);	free(D);	free(ix);
	return 	fix;
}
/* resolve integer ambiguity for DF/TF-PPP -----------------------------------------*/
extern int pppamb(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, const int *refs)				
{
	int m, stat=0, *isat1, *isat2, *NE, *NW;

	NE = imat(n*n, 1);			NW = imat(n*n, 1);				
	isat1 = imat(n*n, 1);		isat2 = imat(n*n, 1);

	trace(2, "%s\n", time_str(rtk->sol.time, 0));
	matcpy(rtk->xa, rtk->x, rtk->nx, 1);
	matcpy(rtk->Pa, rtk->P, rtk->nx, rtk->nx);
	for (int i = 0; i < n*n; i++) NE[i] = NW[i] = MAXCAMB;

	/* single-diff sat */
	m = sel_sat_SD(rtk, obs, n, nav, refs, isat1, isat2);					
	
	/* fix extra wide-lane ambiguity */
	fix_EWL(rtk, obs, nav, isat1, isat2, m, NE);								

	/* fix wide-lane ambiguity */
	fix_WL(rtk, obs, nav, isat1, isat2, m, NW, NE);										

	/* fix narrow-lane ambiguity */
	stat = fix_NL(rtk, obs, nav, isat1, isat2, m, NW);				

	free(isat1);		free(isat2);		 free(NE);		free(NW);	
	return stat;
}