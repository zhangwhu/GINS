#include "rtklib.h"


/* add satellite fcb ---------------------------------------------------------*/
static int addfcb(nav_t *nav, gtime_t ts, gtime_t te, int sat,
	const double *bias, const double *std)
{
	fcbd_t *nav_fcb;
	int i, j;

	if (nav->nf>0 && fabs(timediff(ts, nav->fcb[nav->nf - 1].ts)) <= 1e-3) {
		for (i = 0; i<3; i++) {
			nav->fcb[nav->nf - 1].bias[sat - 1][i] = bias[i];
			nav->fcb[nav->nf - 1].std[sat - 1][i] = std[i];
		}
		return 1;
	}
	if (nav->nf >= nav->nfmax) {
		nav->nfmax = nav->nfmax <= 0 ? 2048 : nav->nfmax * 2;
		if (!(nav_fcb = (fcbd_t *)realloc(nav->fcb, sizeof(fcbd_t)*nav->nfmax))) {
			free(nav->fcb); nav->nf = nav->nfmax = 0;
			return 0;
		}
		nav->fcb = nav_fcb;
	}
	for (i = 0; i<MAXSAT; i++) for (j = 0; j<3; j++) {
		nav->fcb[nav->nf].bias[i][j] = nav->fcb[nav->nf].std[i][j] = 0.0;
	}
	for (i = 0; i<3; i++) {
		nav->fcb[nav->nf].bias[sat - 1][i] = bias[i];
		nav->fcb[nav->nf].std[sat - 1][i] = std[i];
	}
	nav->fcb[nav->nf].ts = ts;
	nav->fcb[nav->nf++].te = te;
	return 1;
}
/* read satellite fcb file ---------------------------------------------------*/
static int readfcbf(const char *file, nav_t *nav)
{
	FILE *fp;
	gtime_t ts, te;
	double ep1[6], ep2[6], bias[3] = { 0 }, std[3] = { 0 };
	char buff[1024], str[32], *p;
	int sat;

	trace(3, "readfcbf: file=%s\n", file);

	if (!(fp = fopen(file, "r"))) {
		trace(2, "fcb parameters file open error: %s\n", file);
		return 0;
	}
	while (fgets(buff, sizeof(buff), fp)) {
		if ((p = strchr(buff, '#'))) *p = '\0';
		if (sscanf(buff, "%lf/%lf/%lf %lf:%lf:%lf %lf/%lf/%lf %lf:%lf:%lf %s"
			"%lf %lf %lf %lf %lf %lf", ep1, ep1 + 1, ep1 + 2, ep1 + 3, ep1 + 4, ep1 + 5,
			ep2, ep2 + 1, ep2 + 2, ep2 + 3, ep2 + 4, ep2 + 5, str, bias, std, bias + 1, std + 1,
			bias + 2, std + 2)<17) continue;
		if (!(sat = satid2no(str))) continue;
		ts = epoch2time(ep1);
		te = epoch2time(ep2);
		if (!addfcb(nav, ts, te, sat, bias, std)) return 0;
	}
	fclose(fp);
	return 1;
}
/* compare satellite fcb -----------------------------------------------------*/
static int cmpfcb(const void *p1, const void *p2)
{
	fcbd_t *q1 = (fcbd_t *)p1, *q2 = (fcbd_t *)p2;
	double tt = timediff(q1->ts, q2->ts);
	return tt<-1E-3 ? -1 : (tt>1E-3 ? 1 : 0);
}
/* read satellite fcb data -----------------------------------------------------
* read satellite fractional cycle bias (dcb) parameters
* args   : char   *file       I   fcb parameters file (wild-card * expanded)
*          nav_t  *nav        IO  navigation data
* return : status (1:ok,0:error)
* notes  : fcb data appended to navigation data
*-----------------------------------------------------------------------------*/
extern int readfcb(const char *file, nav_t *nav)
{
	char *efiles[MAXEXFILE] = { 0 };
	int i, n;

	trace(3, "readfcb : file=%s\n", file);

	for (i = 0; i<MAXEXFILE; i++) {
		if (!(efiles[i] = (char *)malloc(1024))) {
			for (i--; i >= 0; i--) free(efiles[i]);
			return 0;
		}
	}
	n = expath(file, efiles, MAXEXFILE);

	for (i = 0; i<n; i++) {
		readfcbf(efiles[i], nav);
	}
	for (i = 0; i<MAXEXFILE; i++) free(efiles[i]);

	if (nav->nf>1) {
		qsort(nav->fcb, nav->nf, sizeof(fcbd_t), cmpfcb);
	}
	return 1;
}

extern int file2nav(const prcopt_t *prcopt, const filopt_t *fopt, nav_t *nav)
{
	int i, j;
	char* ext;
	gtime_t ts = { 0 }, te = { 0 };

	/* read brdc data */
	if (*fopt->brdc) {
		free(nav->eph);	 nav->eph = NULL;  nav->n = nav->nmax = 0;
		free(nav->geph); nav->geph = NULL; nav->ng = nav->ngmax = 0;
		free(nav->seph); nav->seph = NULL; nav->ns = nav->nsmax = 0;
		if (readrnxt(fopt->brdc, 0, ts, te, 0, prcopt->rnxopt[0], NULL, nav, NULL) < 0) {
			checkbrk("error : insufficient memory");		return 0;
		};
		/* delete duplicated ephemeris */
		uniqnav(nav);
	}

	/* read procise clk */
	if (*fopt->clk && (ext = strrchr(fopt->clk, '.')) && (!strcmp(ext, ".clk") || !strcmp(ext, ".CLK"))) {
		free(nav->pclk); nav->pclk = NULL; nav->nc = nav->ncmax = 0;
		readrnxc(fopt->clk, nav);
	}

	/* read procise orb */
	if (*fopt->sp3 && (ext = strrchr(fopt->sp3, '.')) && (!strcmp(ext, ".sp3") || !strcmp(ext, ".SP3"))) {
		free(nav->peph); nav->peph = NULL; nav->ne = nav->nemax = 0;
		readsp3(fopt->sp3, nav, 0);
	}

	/* read erp data */
	if (*fopt->eop) {
		free(nav->erp.data); nav->erp.data = NULL; nav->erp.n = nav->erp.nmax = 0;
		if (!readerp(fopt->eop, &nav->erp)) {
			checkbrk("no erp data\n");
		}
	}

	/* read dcb parameters */
	if (*fopt->dcb && (ext = strrchr(fopt->dcb, '.')) && ((!strcmp(ext, ".BIA")) || !strcmp(ext, ".DCB"))) {
		readdcb(fopt->dcb, nav);
	}

	/* read GIM inon */
	if (*fopt->iono && (ext = strrchr(fopt->iono, '.')) && (ext[3] == 'i' || ext[3] == 'I')) {
		for (i = 0; i<nav->nt; i++) {
			free(nav->tec[i].data);
			free(nav->tec[i].rms);
		}
		free(nav->tec); nav->tec = NULL; nav->nt = nav->ntmax = 0;
		readtec(fopt->iono, nav, 1);
	}

	/* read fcb file */
	if (*fopt->fcb && (ext = strrchr(fopt->fcb, '.')) /*&& (!strcmp(ext, ".fcb"))*/) {
		free(nav->fcb);	    nav->nf = nav->nfmax = 0;
		readfcb_sgg(fopt->fcb, nav);
	}

	return 1;
}