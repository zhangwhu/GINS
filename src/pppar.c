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

#define MAXCAMB			1E+08

static FILE *fp_stat = NULL;
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
	double dt;
	for (int i = 0; i < nav->nf; i++) {
		dt = fabs(timediff(nav->fcb[i].ti, tt));
		if (dt < nav->fcb[i].iod) {
			return &nav->fcb[i];
		}
	}
	return NULL;
}
static int sel_sat_SD(const rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav, int *isat1, int *isat2)			
{
	int i, j, k, m = 0, sys, sat, isat[4][20] = { 0 }, sn[4] = { 0 }, ind, tempi;
	double azel[4][20], tempd, max, azel0, azel1;
	const ssat_t *ssat;

	for (i = 0; i < n; i++) {
		sat = obs[i].sat;
		ssat = rtk->ssat + sat - 1;
		
		if (ssat->azel[1] < rtk->opt.elmaskar) continue;
		if (!(ssat->vsat[0] & 0x2) || !(ssat->vsat[1] & 0x2)) continue;						
		if (!(sys = satsys(sat, NULL))) continue;
		if (nav->wlbias[sat - 1] > 999) continue;	

		switch (sys){
		case SYS_GPS: isat[0][sn[0]] = i; k = 0; azel[0][sn[0]++] = ssat->azel[1]; break;
		case SYS_GLO: isat[1][sn[1]] = i; k = 1; azel[1][sn[1]++] = ssat->azel[1]; break;
		case SYS_GAL: isat[2][sn[2]] = i; k = 2; azel[2][sn[2]++] = ssat->azel[1]; break;
		case SYS_CMP: isat[3][sn[3]] = i; k = 3; azel[3][sn[3]++] = ssat->azel[1]; break;
		default : k = 0;
		}
	}
	for (i = 0; i<NSYS; i++) {
		for(j = 0; j < sn[i]-1; j++) {
			max = azel[i][j];   ind = j;
			for (k = j+1; k < sn[i]; k++) {
				if (max < azel[i][k]) {
					max = azel[i][k]; 	ind = k; 
				}
			}
			if (ind != j) {
				tempd = azel[i][j]; azel[i][j] = azel[i][ind]; azel[i][ind] = tempd;
				tempi = isat[i][j]; isat[i][j] = isat[i][ind]; isat[i][ind] = tempi;
			}
		}
	}

	/* method 2 */
	for (k = 0; k < 4; k++) {
		if (sn[k] < 2) continue;
		azel0 = rtk->ssat[obs[isat[k][0]].sat-1].azel[1];
		for (i = 0; i < sn[k]-1; i++) {
			// find ref sat
			if (rtk->ssat[obs[isat[k][i]].sat-1].lock[0] < 5 && 
				rtk->ssat[obs[isat[k][i+1]].sat-1].lock[0] > rtk->ssat[obs[isat[k][i]].sat-1].lock[0]) continue;

			// SD 
			azel1 = rtk->ssat[obs[isat[k][i]].sat-1].azel[1];
			ind = (fabs(azel0 - azel1) > 10*D2R) ? i - 1 : i;
			for (j = 0; j < sn[k]; j++) {
				if (ind == j) continue;
				isat1[m] = isat[k][ind];		isat2[m++] = isat[k][j];
			}
			break;
		}
	}

#if 0
	for (i = 0; i < m; i++) {
		char id1[32], id2[32];
		
		satno2id(obs[isat1[i]].sat, id1);
		satno2id(obs[isat2[i]].sat, id2);
		trace(1, "%s--%s %10.4f %10.4f\n", id1, id2, rtk->ssat[obs[isat1[i]].sat-1].azel[1]*R2D, 
			  rtk->ssat[obs[isat2[i]].sat-1].azel[1]*R2D);
	}
#endif
	return m;
}
static int fix_EWL(const rtk_t *rtk, const obsd_t *obs, const nav_t *nav, const int *isat1, const int *isat2, int n, int *NE)				
{
	int i, j, k, f, sat1, sat2, info, nv = 0;
	double LC, BE, freq[NFREQ], lam2, lam3, el, cov[2], temp, *xp, *Pp, *v, *H, *R, *var;
	const prcopt_t *opt = &rtk->opt;

	xp = rtk->xa; 		  Pp = rtk->Pa;
	v = mat(n, 1);		  H = zeros(rtk->nx, n);		R = mat(n, n);		var = mat(n, 1);

	for (i = 0; i < n; i++) {
		sat1 = obs[isat1[i]].sat;
		sat2 = obs[isat2[i]].sat;

		if (rtk->ssat[sat1 - 1].sys == SYS_GPS) continue;		// GPS不参与处理
		if (nav->elbias[sat1 - 1] > 999 || nav->elbias[sat2 - 1] > 999) continue;
		if (!(rtk->ssat[sat1-1].vsat[2] & 0x2) || !(rtk->ssat[sat2-1].vsat[2] & 0x2)) continue;
		
		/* wide-lane ambiguity */
		for (f = 0; f < NFREQ; f++) freq[f] = sat2freq(sat1, obs[isat1[i]].code[f], nav);

		j = IB(sat1, 1, opt);
		k = IB(sat2, 1, opt);
		lam2 = CLIGHT / freq[1];		
		lam3 = CLIGHT / freq[2];
		LC = (xp[j] - xp[k]) / lam2 - (xp[j + MAXSAT] - xp[k + MAXSAT]) / lam3;
		el = (nav->elbias[sat1 - 1] - nav->elbias[sat2 - 1]);

		/* validation of integer wide-lane ambigyity */
		BE = LC - el;
		if (fabs(FROUND(BE)) <= opt->thresar[0]) {			
			NE[i] = ROUND(BE);	

			v[nv] = (NE[i] + el) - LC;
			H[j + nv*rtk->nx] =  1.0 / lam2;
			H[k + nv*rtk->nx] = -1.0 / lam2;
			H[j + MAXSAT + nv*rtk->nx] = -1.0 / lam3;
			H[k + MAXSAT + nv*rtk->nx] =  1.0 / lam3;

			cov[0] = SQR(1 / lam2)*Pp[j + j*rtk->nx] + SQR(1 / lam3)*Pp[(j+MAXSAT) + (j+MAXSAT)*rtk->nx] -
					 2 / (lam2*lam3)*Pp[(j+MAXSAT) + j*rtk->nx];
			cov[1] = SQR(1 / lam2)*Pp[k + k*rtk->nx] + SQR(1 / lam3)*Pp[(k+MAXSAT) + (k+MAXSAT)*rtk->nx] -
					 2 / (lam2*lam3)*Pp[(k+MAXSAT) + k*rtk->nx];
			temp = (cov[1] > cov[0]) ?  cov[0] : cov[1];
			
			var[nv++] = (temp / 100 < SQR(1e-2)) ? temp / 100 : SQR(1e-2);							
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
	if (nv > 0 && (info = filter(xp, Pp, H, v, R, rtk->nx, nv))) {
		trace(1, "EL-filter error (info=%d)\n", info);
		nv = -1;
	}

	free(v);  free(H); free(R); free(var);
	return nv;
}
static int fix_WL(rtk_t *rtk, const obsd_t *obs, const nav_t *nav, const int *isat1, const int *isat2, int n, int *NW)	
{	
	int i, j, k, f, sat1, sat2, info, nv = 0;
	double LC, BW, freq[NFREQ], lam1, lam2, wl, cov[2], temp, *xp, *Pp, *v, *H, *R, *var;
	const prcopt_t *opt = &rtk->opt;

	xp = rtk->xa; 		  Pp = rtk->Pa;
	v = mat(n, 1);		  H = zeros(rtk->nx, n);		R = mat(n, n);		var = mat(n, 1);

	for (i = 0; i < n; i++) {
		sat1 = obs[isat1[i]].sat;
		sat2 = obs[isat2[i]].sat;
		
		/* wide-lane ambiguity */
		for (f = 0; f < NFREQ; f++) freq[f] = sat2freq(sat1, obs[isat1[i]].code[f], nav);
		j = IB(sat1, 0, opt);
		k = IB(sat2, 0, opt);
		lam1 = CLIGHT / freq[0];
		lam2 = CLIGHT / freq[1];
		LC = (xp[j] - xp[k]) / lam1 - (xp[j + MAXSAT] - xp[k + MAXSAT]) / lam2;
		wl = (nav->wlbias[sat1 - 1] - nav->wlbias[sat2 - 1]);		

		/* validation of integer wide-lane ambigyity */
		if (opt->modear == ARMODE_PPPAR) BW = LC - wl;
		else if (opt->modear == ARMODE_PPPAR_ILS) BW = LC + wl;	

		if (fabs(FROUND(BW)) <= opt->thresar[1]) {		
			NW[i] = ROUND(BW);
			rtk->ssat[sat2 - 1].rejw = 0;

			v[nv] = (NW[i] + wl) - LC;
			H[j + nv*rtk->nx] =  1.0 / lam1;
			H[k + nv*rtk->nx] = -1.0 / lam1;
			H[j + MAXSAT + nv*rtk->nx] = -1.0 / lam2;
			H[k + MAXSAT + nv*rtk->nx] =  1.0 / lam2;

			cov[0] = SQR(1 / lam1)*Pp[j + j*rtk->nx] + SQR(1 / lam2)*Pp[(j+MAXSAT) + (j+MAXSAT)*rtk->nx] -
					 2 / (lam1*lam2)*Pp[(j+MAXSAT) + j*rtk->nx];
			cov[1] = SQR(1 / lam1)*Pp[k + k*rtk->nx] + SQR(1 / lam2)*Pp[(k+MAXSAT) + (k+MAXSAT)*rtk->nx] -
					 2 / (lam1*lam2)*Pp[(k+MAXSAT) + k*rtk->nx];
			temp = (cov[1] > cov[0]) ?  cov[0] : cov[1];

			var[nv++] = (temp / 100 < SQR(1e-2)) ? temp / 100 : SQR(1e-2);

			
#if 0
			char id1[32], id2[32];
			satno2id(obs[isat1[i]].sat, id1);
			satno2id(obs[isat2[i]].sat, id2);
			trace(1, "WL-KF %s--%s %10.4f %15.4f %15.4f\n", id1, id2, v[nv-1], LC, NW[i] + wl);
#endif

		} 
		else {
			rtk->ssat[sat2 - 1].rejw++;
#if 0		
			char id1[32], id2[32];
			satno2id(sat1, id1);
			satno2id(sat2, id2);	
			trace(1, "%s %8.0f WL-KF %s--%s %10.4f %15.4f %15.4f %4d %4d %4d\n", time_str(rtk->sol.time,0), 
				time2gpst(rtk->sol.time,NULL),
				id1,id2,(ROUND(BW) + wl)-LC,rtk->ssat[sat1-1].azel[1]*R2D,rtk->ssat[sat2-1].azel[1]*R2D,
				rtk->ssat[sat1-1].lock[0], rtk->ssat[sat2-1].lock[0], rtk->ssat[sat2 - 1].rejw);
#endif
		}
	}
	for (j = 0; j < nv; j++) for (k = 0; k < nv; k++) {
		R[j + k*nv] = j == k ? var[j] : 0.0;
	}

	/* update states with constraints */
	if (nv > 0 && (info = filter(xp, Pp, H, v, R, rtk->nx, nv))) {
		trace(1, "WL-filter error (info=%d)\n", info);
		nv = -1;
	}

	free(v);  free(H); free(R); free(var);
	return nv;
}
static int fix_sol(rtk_t *rtk, const obsd_t *obs, const nav_t *nav, const int *isat1, const int *isat2, 
		   const double *F, int n);
static int fix_NL(rtk_t *rtk, const obsd_t *obs, const nav_t *nav, int *isat1, int *isat2, int n, int *NW)
{
	int i, j, k, f, sat1, sat2, na, naa, *ia, info, fix, minfixsats, sys, stat = 0;
	double *xp, *Pp, *a, *D, *E, *F, *Qaa;
	double LC, nl, lam1, lam2, C1, C2, freq[NFREQ] = { 0.0 }, s[2];
	fcbd_t *fcb;

	if (!(fcb = selfcb(nav, rtk->sol.time))) return 0;

	xp = rtk->xa; 		Pp = rtk->Pa;		
	a = zeros(n, 1);	D = zeros(rtk->nx, n);	 Qaa = E = F = ia = NULL;

	for (i = na = naa = sys = 0; i<n; i++) {
		sat1 = obs[isat1[i]].sat;
		sat2 = obs[isat2[i]].sat;

		if (NW[i] > MAXCAMB - 1) continue;
		if (fcb->bias[sat1 - 1] > 100 || fcb->bias[sat2 - 1] > 100) continue;
		if (rtk->ssat[sat1 - 1].lock[0] < rtk->opt.minlock) continue;		
		if (rtk->ssat[sat1 - 1].lock[1] < rtk->opt.minlock) continue;

		/* float narrow-lane ambiguity (cycle) */
		for (f = 0; f < NFREQ; f++) freq[f] = sat2freq(sat1, obs[isat1[i]].code[f], nav);
		j = IB(sat1, 0, &rtk->opt);
		k = IB(sat2, 0, &rtk->opt);
		lam1 = CLIGHT / freq[0];
		lam2 = CLIGHT / freq[1];
		C1 = freq[0] / (freq[0] - freq[1]);
		C2 = freq[1] / (freq[0] - freq[1]);
		LC = C1*(xp[j] - xp[k]) / lam1 - C2*(xp[j + MAXSAT] - xp[k + MAXSAT]) / lam2 - C2*NW[i];	
		nl = fcb->bias[sat1 - 1] - fcb->bias[sat2 - 1];
		a[na] = LC - nl;

		/* validation of narrow-lane ambiguity */
		if (fabs(FROUND(a[na])) < rtk->opt.thresar[2]) {
			sys |= satsys(sat1, NULL);
			isat1[na] = isat1[i];
			isat2[na] = isat2[i];
			D[j + na*rtk->nx] =  C1 / lam1;
			D[k + na*rtk->nx] = -C1 / lam1;
			D[j + MAXSAT + na*rtk->nx] = -C2 / lam2;
			D[k + MAXSAT + na*rtk->nx] =  C2 / lam2;	
			NW[na++] = NW[i];
#if 0
			char id1[32], id2[32];
			satno2id(obs[isat1[i]].sat, id1);
			satno2id(obs[isat2[i]].sat, id2);
			trace(1, "NL-MM:%s--%s %6.3f IF=%10.4f NW=%6d\n", id1, id2, FROUND(a[na-1]), LC, NW[na-1]);
#endif
		}
		else{
#if 0
			char id1[32], id2[32];
			satno2id(obs[isat1[i]].sat, id1);
			satno2id(obs[isat2[i]].sat, id2);
			trace(1, "%s NL-MM %s--%s %6.3f\n", time_str(rtk->sol.time,0), id1, id2, FROUND(a[na]));		
#endif
		}
	}

	minfixsats = rtk->opt.minfixsats - 1;
	for (i = 0; i<8; i++) minfixsats += (sys&(0x1<<i)) ? 1 : 0;

	if (na >= minfixsats) {				
		/* allocate memory */
		Qaa = mat(na, na);		E = mat(rtk->nx, na);		ia = imat(na, 1);	
		
		/* covariance of narrow-lane ambiguities */
		matmul("NN", rtk->nx, na, rtk->nx, 1.0, Pp, D, 0.0, E);			
		matmul("TN", na, na, rtk->nx, 1.0, D, E, 0.0, Qaa);

		/* decorrelation and sorting  */
		naa = ranking(Qaa, isat1, isat2, a, NW, na, ia);
	}

	if (naa >= minfixsats) {
		F = zeros(naa, 2);	
		for (fix = naa; fix >= minfixsats; fix--) {			
			for (i = 0; i < fix; i++) for (j = 0; j < fix; j++) 
				E[j + i*fix] = Qaa[j + i*naa];
			
			if ((info = lambda(fix, 2, a, E, F, s))) {
				trace(1, "lambda error: info=%d\n", info);  break;
			}

			/* varidation by ratio-test */
			if (s[0] <= 0.0) break;
			rtk->sol.ratio = MIN(s[1] / s[0], 999.9);
			if (rtk->sol.ratio > rtk->opt.thresar[3]) {
				stat = 1;	 break;
			}
		}		
	}

	if (stat) {
		for(i = 0; i < fix; i++) {
			sat1 = obs[isat1[i]].sat;
			sat2 = obs[isat2[i]].sat;
			nl = fcb->bias[sat1 - 1] - fcb->bias[sat2 - 1];
			for (f = 0; f < NFREQ; f++) freq[f] = sat2freq(sat1, obs[isat1[i]].code[f], nav);
			F[i] = F[i] + nl + (freq[1] / (freq[0] - freq[1]))*NW[i];
		}
		stat = fix_sol(rtk, obs, nav, isat1, isat2, F, fix);
	}

	free(a); free(D); free(Qaa); free(E); free(F); free(ia);
	return 	stat ? fix : 0;
}
static int fix_sol(rtk_t *rtk, const obsd_t *obs, const nav_t *nav, const int *isat1, const int *isat2, 
		   const double *F, int n) 
{
	int i, j, k, f, sat1, sat2, info, stat = 1;
	double lam1, lam2, C1, C2, freq[NFREQ], cov[2], temp;
	double *xp, *Pp, *v, *H, *R, *var;

	xp = rtk->xa;	   Pp = rtk->Pa;
	v = zeros(n, 1);   H = zeros(rtk->nx, n);    R = zeros(n, n);     var = zeros(n, 1);

	for (i = 0; i < n; i++) {
		sat1 = obs[isat1[i]].sat;
		sat2 = obs[isat2[i]].sat;
		j = IB(sat1, 0, &rtk->opt);
		k = IB(sat2, 0, &rtk->opt);

		for (f = 0; f < NFREQ; f++) freq[f] = sat2freq(sat1, obs[isat1[i]].code[f], nav);
		lam1 = CLIGHT / freq[0];
		lam2 = CLIGHT / freq[1];
		C1 = freq[0] / (freq[0] - freq[1]);
		C2 = freq[1] / (freq[0] - freq[1]);

		v[i] = F[i] - (C1*(xp[j] - xp[k]) / lam1 - C2*(xp[j + MAXSAT] - xp[k + MAXSAT]) / lam2);
		H[j + i*rtk->nx] =  C1 / lam1;
		H[k + i*rtk->nx] = -C1 / lam1;
		H[j + MAXSAT + i*rtk->nx] = -C2 / lam2;
		H[k + MAXSAT + i*rtk->nx] =  C2 / lam2;
		
		cov[0] = SQR(C1 / lam1)*Pp[j + j*rtk->nx] + SQR(C2 / lam2)*Pp[(j+MAXSAT) + (j+MAXSAT)*rtk->nx] -
					 2*(C1*C2)/(lam1*lam2) * Pp[(j+MAXSAT) + j*rtk->nx];

		cov[1] = SQR(C1 / lam1)*Pp[k + k*rtk->nx] + SQR(C2 / lam2)*Pp[(k+MAXSAT) + (k+MAXSAT)*rtk->nx] -
					 2*(C1*C2)/(lam1*lam2) * Pp[(k+MAXSAT) + k*rtk->nx];

		temp = (cov[1] > cov[0]) ?  cov[0] : cov[1];
		var[i] = (temp / 100 < SQR(1e-2)) ? temp / 100 : SQR(1e-2);

#if 0
		char id1[8], id2[8];
		satno2id(sat1, id1);
		satno2id(sat2, id2);
		trace(1,"SOL: %s %s %10.4f %15.6f %15.6f\n", id1, id2, v[i], F[i], 
			(C1*(xp[j] - xp[k]) / lam1 - C2*(xp[j + MAXSAT] - xp[k + MAXSAT]) / lam2));
#endif
	}

	for (j = 0; j < n; j++) for (k = 0; k < n; k++) {
		R[j + k*n] = j == k ? var[j] : 0.0;
	}

	/* update states with constraints */
	if (n > 0 && (info = filter(xp, Pp, H, v, R, rtk->nx, n))) {
		trace(1, "NL-filter error (info=%d)\n", info);
		stat = 0;
	}

	/* set flags */
	if (stat) {
		for (i = 0; i<n; i++) {
			sat1 = obs[isat1[i]].sat;
			sat2 = obs[isat2[i]].sat;
			rtk->ssat[sat1 - 1].fix[0] = rtk->ssat[sat1 - 1].fix[1] = 2;
			rtk->ssat[sat2 - 1].fix[0] = rtk->ssat[sat2 - 1].fix[1] = 2;
			// rtk->ssat[sat1 - 1].vsat[0] |= 0x10;
			// rtk->ssat[sat2 - 1].vsat[0] |= 0x10;
		}
	}

	free(v); free(H); free(R); free(var);
	return stat;
}
/* resolve integer ambiguity for DF/TF-PPP -----------------------------------------*/
extern int pppamb(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)			
{
	int *isat1, *isat2, *NE, *NW, m, stat = 0;

	NE = imat(n*n, 1);			NW = imat(n*n, 1);				
	isat1 = imat(n*n, 1);		isat2 = imat(n*n, 1);

	trace(2, "%s %10.0f\n", time_str(rtk->sol.time, 0), time2gpst(rtk->sol.time, NULL));
	matcpy(rtk->xa, rtk->x, rtk->nx, 1);			
	matcpy(rtk->Pa, rtk->P, rtk->nx, rtk->nx);
	for (int i = 0; i < n*n; i++) NE[i] = NW[i] = MAXCAMB;

	/* single-diff sat */
	m = sel_sat_SD(rtk, obs, n, nav, isat1, isat2);	

	/* fix extra wide-lane ambiguity */
	if (fix_EWL(rtk, obs, nav, isat1, isat2, m, NE) >= 0) {

		/* fix wide-lane ambiguity */
		if (fix_WL(rtk, obs, nav, isat1, isat2, m, NW) >= 0) {	
			
			/* fix narrow-lane ambiguity */
			stat = fix_NL(rtk, obs, nav, isat1, isat2, m, NW);			
		}	
	}
	trace(2, "ratio=%5.2f\n", rtk->sol.ratio);
	free(isat1);		free(isat2);		 free(NE);		free(NW);	
	return stat;
}