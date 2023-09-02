
#include "rtklib.h"

#define MAXITR      10          /* max number of iteration for tdcp pos */

static int  restdc(const obsd_t *obs, int n, const nav_t *nav, const double *rs, const double *dts, 
	const int *svh, const ssat_t *ssat, const double *prr, int nf, const prcopt_t *opt, double *x, double *v, 
	double *H, int *vsat, double *resp, int *exc) 
{
    int i, j, k, nv = 0;
    double freq, sclk, deltaD, deltaG, xcorr, dsagnac, var;
    double r1, r2, e1[3], e2[3], rr[3], vv;
	const ssat_t *psat;

    for (i = 0; i<n; i++) {
		for (k = 0; k < nf; k++) vsat[k + nf*i] = 0;

		psat = ssat + obs[i].sat - 1;

        if (!psat->vs || svh[i] < 0) continue;
		
        if (fabs(fabs(timediff(obs[i].time, psat->prev_T[0])) - opt->ti) > DTTOL) continue;

        if ((r1 = geodist(psat->prev_rs, prr, e1)) <= 0) continue;		
		for (k = 0; k < 3; k++) rr[k] = prr[k] + x[k];
		if ((r2 = geodist(rs + 6*i, rr, e2)) <= 0) continue;				

		dsagnac = OMGE*((rs[6*i+0]*rr[1]-rs[6*i+1]*rr[0]) - (psat->prev_rs[0]*prr[1]-psat->prev_rs[1]*prr[0]))/CLIGHT;

		// char id[32], time[32];
		// satno2id(obs[i].sat, id);
		// time2str(obs[i].time, time, 0);
		// if (psat->slip[0]&01) trace(1, "TDCP %s %s \n", id, time);

		xcorr = -(e2[0]*x[0] + e2[1]*x[1] + e2[2]*x[2]) + x[3];
        sclk   = CLIGHT*(dts[i * 2] - psat->prev_dts[0]);
        deltaD = dot(e2, rs + 6*i, 3) - dot(e1, psat->prev_rs, 3);
        deltaG = dot(e2, prr, 3) - dot(e1, prr, 3);

        for (k = 0; k < nf; k++) {
			
			freq = sat2freq(obs[i].sat, obs[i].code[k], nav);

            if (ISZERO(obs[i].L[k]) || ISZERO(freq) || ISZERO(psat->prev_L[0][k]) || psat->slip[k]&1) continue;
			
			vv = CLIGHT*(obs[i].L[k] - psat->prev_L[0][k])/freq + sclk - deltaD + deltaG - xcorr - dsagnac;

			resp[k + nf*i] = vv / (CLIGHT / freq);

			if (exc[k+i*nf]) continue;
			
			var = varerr(obs[i].sat, psat->azel[1], k, 0, opt);

            v[nv]  = vv / sqrt(var);

            /* design matrix */
            for (j = 0; j<4; j++) {
			    H[j + 4*nv] = ((j<3) ? -e2[j] : 1.0) / sqrt(var);
		    }

			vsat[k + nf*i] = 1;
			nv++;
        }
    }
	return nv;
}
extern int tdcpos(const obsd_t *obs, int n, const nav_t *nav, const prcopt_t *opt, const double *rs, 
	const double *dts, const int *svh, sol_t *sol, ssat_t *ssat) 
{
    int i, j, k, f, nv, sat, info, *vsat, *exc, maxi, nf = opt->nf, stat = 0;
    double x[4] = { 0 }, dx[4] = { 0 }, Q[16], *v, *H, *resp, vmax, *pr;

	vsat = imat(nf*n, 1);	exc = imat(nf*n, 1); 
    v = mat(nf*n, 1);    	H = mat(4, nf*n);		resp = mat(nf*n, 1); 		
	
    /* tdcp residuals */
	pr = ISZERO(sol->pr[0]) ? sol->rr : sol->pr;
	for (i = 0; i < nf*n; i++) resp[i] = exc[i] = 0;
	for (i = 0; i < MAXITR; i++) {	
		nv = restdc(obs, n, nav, rs, dts, svh, ssat, pr, nf, opt, x, v, H, vsat, resp, exc);	
		if (nv<8) {
			trace(2, "lack of valid sats ns=%d", nv);
			break;
		}

		/* least square estimation */
		if ((info = lsq(H, v, 4, nv, dx, Q))) {
			trace(1, "lsq error info=%d", info);
			break;
		}
		for (j = 0; j<4; j++) x[j] += dx[j];
		if (norm(dx,4)>1E-4) continue;	

		vmax = maxi = 0;
		restdc(obs, n, nav, rs, dts, svh, ssat, pr, nf, opt, x, v, H, vsat, resp, exc);
		for (k = 0; k < n; k++) {
			for (f = 0; f < nf; f++) {
				if (vsat[f+k*nf] && fabs(resp[f+k*nf]) > vmax) {
					vmax = fabs(resp[f+k*nf]);
					maxi = f + k*nf;
				}
			}
		}
		if (vmax < opt->threscheck[2]) {	stat = 1; break; }
		exc[maxi] = 1;
	}

    if (stat) {
		matcpy(sol->pr+3, x, 3, 1);			//修改
		for (k = 0; k < n; k++) {
			sat = obs[k].sat;
			ssat[sat-1].diff[6] = ssat[sat-1].diff[7] = ssat[sat-1].diff[8] = 0.0;

			for (f = 0; f < nf; f++) {
				if (exc[f + nf*k]) {
					ssat[sat-1].slip[f] |= 1;			//需要修改2022-09-15
				}
				ssat[sat-1].diff[6+f] = resp[f + nf*k];

				// char id[32], time[32];
				// satno2id(sat, id);
				// time2str(sol->time, time,0);
				// trace(1,"%s %s %2d %10.3f\n",time, id, f, resp[f + nf*k]);
			}
		}
    }

    free(v);    free(H);	free(resp);	  free(vsat);	free(exc);
    return stat;
}