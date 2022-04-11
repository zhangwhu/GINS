#include "rtklib.h"
#include "ppp_state.h"

#define SQR(x)      ((x)*(x))
#define ionoThread 0.08
#define tropThread 0.1

static outtrop_header(FILE *fp)
{

}
static outiono_header(FILE *fp)
{
	
}
static outfcbs_header(FILE *fp)
{
	
}
extern void outppp(const char *outfile, const prcopt_t *popt, const char *staname, const double *pos, char type)
{
	FILE *fp = NULL;

	trace(3, "outppp: outfile=%s\n", outfile);

	if (*outfile) {
		createdir(outfile);
		if (!(fp = fopen(outfile, "w"))) {
			showmsg("error : open output file %s", outfile);
			return ;
		}
	}
	if (fp == NULL) return;

	/* output header */
	fprintf(fp, "%s staname   : %s\n", COMMENTH, staname);
	fprintf(fp, "%s position  :%14.5f %14.5f %14.5f\n", COMMENTH, pos[0], pos[1], pos[2]);
	outprcopt(fp, popt);
	fprintf(fp, "%s\n", COMMENTH);
	switch (type) {
		case 'T': outtrop_header(fp); break;
		case 'I': outiono_header(fp); break; 
		case 'F': outfcbs_header(fp); break; 
		default: break;
	}
	if (*outfile) fclose(fp);
}

/* write outbody file ---------------------------------------------------- */
extern void outtrop(const FILE *fp, const rtk_t *rtk, const obsd_t *obs, const int *n, const nav_t *nav)
{
	if (!fp || rtk->sol.stat == SOLQ_NONE || rtk->sol.stat == SOLQ_SINGLE || rtk->opt.tropopt != TROPOPT_EST) return;

	int ii;
	char id[32];
	double *x, *P, pos[3], zazel[] = { 0.0, PI / 2.0 }, zhd;

	x = rtk->sol.stat == SOLQ_FIX ? rtk->xa : rtk->x;
	P = rtk->sol.stat == SOLQ_FIX ? rtk->Pa : rtk->P;

	ii = IT(&rtk->opt);
	if (sqrt(P[ii + ii*rtk->nx]) > tropThread) return;

	ecef2pos(rtk->sol.rr, pos);
	zhd = tropmodel(rtk->sol.time, pos, zazel, 1, 0, 1);
	fprintf(fp, "%s %10.4f %10.4f\n", time_str(rtk->sol.time, 0), x[ii] + zhd, sqrt(P[ii + ii*rtk->nx]));
}
extern void outiono(const FILE *fp, const rtk_t *rtk, const obsd_t *obs, const int *n, const nav_t *nav)
{
	if (!fp || rtk->sol.stat == SOLQ_NONE || rtk->sol.stat == SOLQ_SINGLE || rtk->opt.ionoopt != IONOOPT_EST) return;

	ambc_t *amb;
	ssat_t *ssat;

	char id[32];
	int i, k, ii;
	double *x, *P, freq1, fact, stec;

	x = rtk->sol.stat == SOLQ_FIX ? rtk->xa : rtk->x;
	P = rtk->sol.stat == SOLQ_FIX ? rtk->Pa : rtk->P;

	for (i = 0; i < n; i++) {
		amb = rtk->ambc + obs[i].sat - 1;
		ssat = rtk->ssat + obs[i].sat - 1;
		satno2id(obs[i].sat, id);

		ii = II(obs[i].sat, &rtk->opt);
	
		if (!ssat->vs || !ssat->vsat[0] || !ssat->vsat[1] ||  ssat->azel[1] < 15 * D2R /*|| sqrt(P[ii + ii*rtk->nx]) > ionoThread*/) continue;

		freq1 = sat2freq(obs[i].sat, obs[i].code[0], nav);
		fact = 40.30E16 / SQR(freq1);
		stec = x[ii] / fact;
		fprintf(fp, "%s %s %10.3f %10.3f %10.4f %10.4f\n", time_str(rtk->sol.time, 0), id,
			ssat->azel[0] * R2D, ssat->azel[1] * R2D, stec, sqrt(P[ii + ii*rtk->nx] / SQR(fact)));
	}
}
extern void outfcbs(const FILE *fp, const rtk_t *rtk, const obsd_t *obs, const int *n, const nav_t *nav)
{
	if (!fp || rtk->sol.stat == SOLQ_NONE || rtk->sol.stat == SOLQ_SINGLE) return;

	ssat_t *ssat;
	double freq[NFREQ], pos[3];
	char id[32];
	int i, j, ll, mm, nn, vsat;
	double *x, *P, a, b, lam1, lam2, lam3, EL, WL, IF, Pe, Pw, Pf;

	x = rtk->sol.stat == SOLQ_FIX ? rtk->xa : rtk->x;
	P = rtk->sol.stat == SOLQ_FIX ? rtk->Pa : rtk->P;

	/* fcb */
	for (i = 0; i < n; i++) {
		ssat = rtk->ssat + obs[i].sat - 1;
		if (!ssat->vs) continue;

		EL = WL = IF = 0.0;			Pe = Pw = Pf = 1e6;
		ll = IB(obs[i].sat, 0, &rtk->opt);
		mm = IB(obs[i].sat, 1, &rtk->opt);

		for (j = 0; j < NFREQ; j++) freq[j] = sat2freq(obs[i].sat, obs[i].code[j], nav);
		lam1 = CLIGHT / freq[0];
		lam2 = CLIGHT / freq[1];
		lam3 = CLIGHT / freq[2];

		//el ambiguity
		if (ssat->vsat[1] && ssat->vsat[2]) {
			nn = IB(obs[i].sat, 2, &rtk->opt);
			EL = x[mm] / lam2 - x[nn] / lam3;
			Pe = SQR(1 / lam2)*rtk->P[mm + mm*rtk->nx] + SQR(1 / lam3)*rtk->P[nn + nn*rtk->nx] - 2 / (lam2*lam3)*rtk->P[nn + mm*rtk->nx];
		}
		//wl&if ambiguity
		if (ssat->vsat[0] && ssat->vsat[1]) {
			WL = x[ll] / lam1 - x[mm] / lam2;
			Pw = SQR(1 / lam1)*rtk->P[ll + ll*rtk->nx] + SQR(1 / lam2) *rtk->P[mm + mm*rtk->nx] - 2 / (lam1*lam2)*rtk->P[mm + ll*rtk->nx];

			a = SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));
			b = SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));
			IF = a*x[ll] - b*x[mm];
			Pf = a*a*rtk->P[ll + ll*rtk->nx] + b*b*rtk->P[mm + mm*rtk->nx] - 2 * a*b*rtk->P[mm + ll*rtk->nx];
		}

		satno2id(obs[i].sat, id);
		fprintf(fp, "%s %s %4d %6.3f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n", time_str(rtk->sol.time, 0), id, 
			ssat->slip[0], ssat->azel[1]*R2D, EL, sqrt(Pe), WL, sqrt(Pw), IF, sqrt(Pf));
	}
	return;
}
extern void outelse(const FILE *fp, const rtk_t *rtk, const obsd_t *obs, const int *n, const nav_t *nav) 
{
	int ii, jj, kk;
	double *x, *P;
	if (!fp || rtk->sol.stat == SOLQ_NONE || rtk->sol.stat == SOLQ_SINGLE) return;

	x = rtk->sol.stat == SOLQ_FIX ? rtk->xa : rtk->x;
	P = rtk->sol.stat == SOLQ_FIX ? rtk->Pa : rtk->P;

	ii = IC(0, &rtk->opt);
	jj = IC(2, &rtk->opt);
	kk = IC(3, &rtk->opt);

	fprintf(fp, "%s %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n", time_str(rtk->sol.time, 0), x[ii]-x[jj],
		x[jj] - x[kk], x[ii] - x[kk], x[ii] * 1E6 / CLIGHT, x[jj] * 1E6 / CLIGHT, x[kk] * 1E6 / CLIGHT, rtk->sol.dsp[0], rtk->sol.dsp[2], rtk->sol.dsp[3]);
	/*int sat;
	char id[32];
	double da;
	for (int i = 0; i < n; i++) {
		sat = obs[i].sat;
		if (rtk->ssat[sat - 1].vs) {
			satno2id(sat,id);
			fprintf(fp, "%s %s %10.4f ", time_str(rtk->sol.time, 0), id, rtk->ssat[sat-1].azel[1]*R2D);
			for (int k = 0; k < 3; k++) {
				da = rtk->ssat[sat - 1].resp[k];
				fprintf(fp, "%12.4f", da);
			}
			for (int k = 0; k < 3; k++) {
				da = 0.0;
				if (rtk->ssat[sat - 1].vsat[k]) da = rtk->ssat[sat - 1].resc[k];
				fprintf(fp, "%12.4f", da);
			}
			fprintf(fp, "\n");
		}
	}*/
}