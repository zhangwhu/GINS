
#include "rtklib.h"

#define MIN(x,y)    ((x)<(y)?(x):(y))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))

#define MAXPRCDAYS  100				/* max days of continuous processing */
#define MAXINFILE   1000			/* max number of input files */ 

/* constants/global variables ------------------------------------------------*/
static pcvs_t pcvss = { 0 };        /* receiver antenna parameters */
static pcvs_t pcvsr = { 0 };        /* satellite antenna parameters */
static obs_t obss = { 0 };          /* observation data */
static nav_t navs = { 0 };          /* navigation data */			
static stas_t stas = { 0 };			/* station list */
static int nepoch = 0;				/* number of observation epochs */
static int iobsu = 0;				/* current rover observation data index */
static int iobsr = 0;				/* current reference observation data index */
static int revs = 0;				/* analysis direction (0:forward,1:backward) */
static sol_t *solf;					/* forward solutions */
static sol_t *solb;					/* backward solutions */
static double *rbf;					/* forward base positions */
static double *rbb;					/* backward base positions */
static int isolf = 0;				/* current forward solutions index */
static int isolb = 0;				/* current backward solutions index */

/* show message --------------------------------------------------------------*/
extern int showmsg(char *format, ...)
{
	va_list arg;
	va_start(arg, format); vfprintf(stderr, format, arg); va_end(arg);
	fprintf(stderr, "\r");				/* \r���������� */
	return 0;
}
extern void settspan(gtime_t ts, gtime_t te) {}
extern void settime(gtime_t time) {}

/* show message and check break ----------------------------------------------*/
extern int checkbrk(const char *format, ...)
{
	va_list arg;
	char buff[1024], *p = buff;
	if (!*format) return showmsg("");
	va_start(arg, format);
	p += vsprintf(p, format, arg);
	va_end(arg);
	return showmsg(buff);
}
/* open output file for append -----------------------------------------------*/
static FILE *openfile(const char *outfile)
{
	trace(3, "openfile: outfile=%s\n", outfile);

	return !*outfile ? stdout : fopen(outfile, "a");
}
/* free prec ephemeris and sbas data -----------------------------------------*/
static void freedata(nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr, stas_t *stas)
{
	free(nav->eph);		nav->eph = NULL;	nav->n = nav->nmax = 0;
	free(nav->geph);	nav->geph = NULL;	nav->ng = nav->ngmax = 0;
	free(nav->seph);	nav->seph = NULL;	nav->ns = nav->nsmax = 0;
	free(nav->peph);	nav->peph = NULL;	nav->ne = nav->nemax = 0;
	free(nav->pclk);	nav->pclk = NULL;	nav->nc = nav->ncmax = 0;
	free(nav->fcb);		nav->fcb = NULL;	nav->nf = nav->nfmax = 0;

	for (int i = 0; i<nav->nt; i++) {
		free(nav->tec[i].data);
		free(nav->tec[i].rms);
	}
	free(nav->tec);		nav->tec = NULL;	nav->nt = nav->ntmax = 0;

	/* free erp data */
	free(nav->erp.data); nav->erp.data = NULL; nav->erp.n = nav->erp.nmax = 0;

	/* free antenna parameters */
	free(pcvs->pcv);	pcvs->pcv = NULL;	pcvs->n = pcvs->nmax = 0;
	free(pcvr->pcv);	pcvr->pcv = NULL;	pcvr->n = pcvr->nmax = 0;

	/* free station list */
	free(stas->data);   stas->data = NULL;  stas->n = stas->nmax = 0;
}
/* output header -------------------------------------------------------------*/
static void outheader(FILE *fp, const prcopt_t *popt, const solopt_t *sopt)
{
	const char *s1[] = { "GPST", "UTC", "JST" };
	gtime_t ts, te;
	double t1, t2;
	int i, j, w1, w2;
	char s2[32], s3[32];
	const char *sep = sopt->sep;

	trace(3, "outheader: n=%d\n");

	if (sopt->posf == SOLF_NMEA || sopt->posf == SOLF_STAT)  return;
	if (sopt->outhead) {
		fprintf(fp, "%s program   : RTKLIB ver.%s %s\n", COMMENTH, VER_RTKLIB, PATCH_LEVEL);

		for (i = 0; i<obss.n; i++)    if (obss.data[i].rcv == 1) break;
		for (j = obss.n - 1; j >= 0; j--) if (obss.data[j].rcv == 1) break;
		if (j<i) { fprintf(fp, "\n%s no rover obs data\n", COMMENTH); return; }
		ts = obss.data[i].time;
		te = obss.data[j].time;
		t1 = time2gpst(ts, &w1);
		t2 = time2gpst(te, &w2);
		if (sopt->times >= 1) ts = gpst2utc(ts);
		if (sopt->times >= 1) te = gpst2utc(te);
		if (sopt->times == 2) ts = timeadd(ts, 9 * 3600.0);
		if (sopt->times == 2) te = timeadd(te, 9 * 3600.0);
		time2str(ts, s2, 1);
		time2str(te, s3, 1);
		fprintf(fp, "%s obs start : %s %s (week%04d %8.1fs)\n", COMMENTH, s2, s1[sopt->times], w1, t1);
		fprintf(fp, "%s obs end   : %s %s (week%04d %8.1fs)\n", COMMENTH, s3, s1[sopt->times], w2, t2);
	}
	if (sopt->outopt) {
		outprcopt(fp, popt);				
	}

	if (sopt->posf == SOLF_ENU) {
		fprintf(fp, "%s ref pos   :%14.4f%s%14.4f%s%14.4f\n", COMMENTH, popt->ru[0], sep, popt->ru[1], sep, popt->ru[2]);
	}
	if (sopt->outhead || sopt->outopt) fprintf(fp, "%s\n", COMMENTH);
	outsolhead(fp, sopt);
}
/* write header to output file -----------------------------------------------*/
static int outhead(const char *outfile, const prcopt_t *popt, const solopt_t *sopt)		
{
	FILE *fp = stdout;

	trace(3, "outhead: outfile=%s\n", outfile);

	if (*outfile) {
		createdir(outfile);						
		if (!(fp = fopen(outfile, "w"))) {
			showmsg("error : open output file %s", outfile);
			return 0;
		}
	}
	/* output header */
	outheader(fp, popt, sopt);

	if (*outfile) fclose(fp);

	return 1;
}
/* read obs data ----------------------------------------------------- */
static int readobs(const char *infile, int rcv, const prcopt_t *prcopt, obs_t *obs, sta_t *sta)
{
	int i, j, ind = 0, nobs = 0;
	gtime_t ts = { 0 }, te = { 0 };
	double ti = 0;
	trace(3, "readobs:  n=%d\n");

	obs->data = NULL; obs->n = obs->nmax = 0;
	nepoch = 0;

	if (checkbrk("")) return 0;

	/* read rinex obs and nav file */
	if (readrnxt(infile, rcv, ts, te, ti, prcopt->rnxopt[0], obs, NULL, sta)<0) {
		checkbrk("error : insufficient memory");
		trace(1, "insufficient memory\n");
		return 0;
	}
	if (obs->n <= 0) {
		checkbrk("error : no obs data");
		trace(1, "\n");
		return 0;
	}
	/* sort observation data */
	nepoch = sortobs(obs);

	return 1;
}
/* read nav data ----------------------------------------------------- */
static int readfile(const prcopt_t *popt, const filopt_t *fopt)
{
	/* generate stas */		
	if (!readstas(fopt->staname, &stas)) return 0;

	/* generate nav */
	if (!file2nav(popt, fopt, &navs)) return 0;

	/* generate ppp corr */
	if (!file2corr(fopt)) return 0;

	/* read satellite antenna parameters */
	if (*fopt->satantp) {
		free(pcvss.pcv);	pcvss.n = pcvss.nmax = 0;
		if (!(readpcv(fopt->satantp, &pcvss))) {
			showmsg("error : no sat ant pcv in %s", fopt->satantp);
			return 0;
		}
		/* use satellite L2 offset if L5 offset does not exists */
		for (int i = 0; i<pcvss.n; i++) {
			if (norm(pcvss.pcv[i].off[2], 3)>0.0) continue;
			matcpy(pcvss.pcv[i].off[2], pcvss.pcv[i].off[1], 3, 1);
			matcpy(pcvss.pcv[i].var[2], pcvss.pcv[i].var[1], 19, 1);
		}
	}

	/* read receiver antenna parameters */
	if (*fopt->rcvantp) {
		free(pcvsr.pcv);	pcvsr.n = pcvsr.nmax = 0;
		if (!(readpcv(fopt->rcvantp, &pcvsr))) {
			showmsg("error : no rec ant pcv in %s", fopt->rcvantp);
			return 0;
		}
		/* use receive L2 offset if L5 offset does not exists */
		for (int i = 0; i<pcvsr.n; i++) {
			if (norm(pcvsr.pcv[i].off[2], 3)>0.0) continue;
			matcpy(pcvsr.pcv[i].off[2], pcvsr.pcv[i].off[1], 3, 1);
			matcpy(pcvsr.pcv[i].var[2], pcvsr.pcv[i].var[1], 19, 1);
		}
	}
	return 1;
}
/* search next observation data index ----------------------------------------*/
static int nextobsf(const obs_t *obs, int *i, int rcv)			
{
	int n;
	double tt;

	for (; *i<obs->n; (*i)++) if (obs->data[*i].rcv == rcv) break;
	for (n = 0; *i + n<obs->n; n++) {
		tt = timediff(obs->data[*i + n].time, obs->data[*i].time);
		if (obs->data[*i + n].rcv != rcv || tt>DTTOL) break;
	}
	return n;   
}
static int nextobsb(const obs_t *obs, int *i, int rcv)			
{
	double tt;
	int n;

	for (; *i >= 0; (*i)--) if (obs->data[*i].rcv == rcv) break;
	for (n = 0; *i - n >= 0; n++) {
		tt = timediff(obs->data[*i - n].time, obs->data[*i].time);
		if (obs->data[*i - n].rcv != rcv || tt<-DTTOL) break;				
	}
	return n;
}
/* input obs data, navigation messages and sbas correction -------------------*/
static int inputobs(obsd_t *obs, int solq, const prcopt_t *popt)
{
	gtime_t time = { 0 };
	int i, nu, n = 0;

	trace(3, "infunc  : revs=%d iobsu=%d\n", revs, iobsu);					

	if (0 <= iobsu&&iobsu<obss.n) {
		settime((time = obss.data[iobsu].time));					
		if (checkbrk("processing : %s Q=%d", time_str(time, 0), solq)) {
			showmsg("aborted"); return -1;
		}
	}
	if (!revs) { /* input forward data */
		if ((nu = nextobsf(&obss, &iobsu, 1)) <= 0) return -1;
		for (i = 0; i<nu&&n<MAXOBS * 2; i++) obs[n++] = obss.data[iobsu + i];
		iobsu += nu;
	}
	else { /* input backward data */
		if ((nu = nextobsb(&obss, &iobsu, 1)) <= 0) return -1;
		for (i = 0; i<nu&&n<MAXOBS * 2; i++) obs[n++] = obss.data[iobsu - nu + 1 + i];
		iobsu -= nu;
	}
	return n;
}
/* process positioning ----------------------------------------------------*/
static void procpos(const prcopt_t *popt, const solopt_t *sopt, const filopt_t *fopt, rtk_t *rtk, int mode)				
{
	obsd_t obs[MAXOBS];
	FILE *fp_poss, *fp_iono, *fp_trop, *fp_fcbs, *fp_else;
	int i, nobs, n;
	gtime_t ts;
	char stime[32] = "2021 08 11 04 00 00.00";
	str2time(stime, 0, 32, &ts);

	if (mode == 0) iobsu = 0; else iobsu = nepoch - 1; 
	fp_poss = fp_iono = fp_trop = fp_fcbs = fp_else = NULL;

	if (*fopt->outposs) fp_poss = openfile(fopt->outposs);
	if (*fopt->outtrop) fp_trop = openfile(fopt->outtrop);
	if (*fopt->outiono) fp_iono = openfile(fopt->outiono);
	if (*fopt->outfcbs) fp_fcbs = openfile(fopt->outfcbs);
	if (*fopt->outelse) fp_else = openfile(fopt->outelse);

	pppinit(rtk, popt);

	while ((nobs = inputobs(obs, rtk->sol.stat, popt)) >= 0) {			
		/* exclude satellites */
		for (i = n = 0; i<nobs; i++) {
			if ((satsys(obs[i].sat, NULL)&popt->navsys) &&
				popt->exsats[obs[i].sat - 1] != 1) obs[n++] = obs[i];
		}
		if (n <= 0) continue;
		if (timediff(obs[0].time, ts) < 1E-2) continue;

		if (!ppprocess(rtk, obs, n, &navs)) continue;															
		if (mode == 0) { 
			if (rtk->sol.stat == SOLQ_PPP || rtk->sol.stat == SOLQ_FIX || rtk->sol.stat == SOLQ_SINGLE)	{
				if (fp_poss) outsol (fp_poss, &rtk->sol, rtk->opt.ru, sopt);
				if (fp_trop) outtrop(fp_trop, rtk, obs, n, &navs);
				if (fp_iono) outiono(fp_iono, rtk, obs, n, &navs);
				if (fp_fcbs) outfcbs(fp_fcbs, rtk, obs, n, &navs);
				if (fp_else) outelse(fp_else, rtk, obs, n, &navs);
			}
		}
		else {
		
		}
	}

	if (fp_poss) fclose(fp_poss);	if (fp_iono) fclose(fp_iono);
	if (fp_trop) fclose(fp_trop);	if (fp_fcbs) fclose(fp_fcbs);	if (fp_else) fclose(fp_else);

	checkbrk("                                        ");
}
/* post-processing positioning -------------------------------------------------*/
extern int postpos(const prcopt_t *popt, const solopt_t *sopt, filopt_t *fopt)		
{
	rtk_t rtk;
	prcopt_t popt_ = *popt;
	char obsfile[1024], posfile[512], staname[32];

	/* read nav & atx & corr data */					
	if (!readfile(popt, fopt)) return 0;			
	
	for (int i = 0; i < stas.n; i++) {	
		unsigned int tick = tickget();
		strcpy(staname, stas.data[i].name);
		trace(1, "(%d)station: %s\n", i + 1, staname);

		fopt->outposs[0] = fopt->outtrop[0] = fopt->outiono[0] = fopt->outfcbs[0] = fopt->outelse[0] = '\0';

		sprintf(obsfile, "%s%s%s", fopt->indir,  staname, "*.*o");	
		sprintf(posfile, "%s%s%s", fopt->outdir, staname, "-COD-par.pos");		strcpy(fopt->outposs, posfile);
		//sprintf(posfile, "%s%s%s", fopt->outdir, staname, ".trop");		strcpy(fopt->outtrop, posfile);
		//sprintf(posfile, "%s%s%s", fopt->outdir, staname, ".stec");		strcpy(fopt->outiono, posfile);
		//sprintf(posfile, "%s%s%s", fopt->outdir, staname, "-COD.fcbs");		strcpy(fopt->outfcbs, posfile);
		//sprintf(posfile, "%s%s%s", fopt->outdir, staname, "-CLOCK.else");		strcpy(fopt->outelse, posfile);

		/* read obs &&  set antenna paramters &&  set ocean tide loading parameters */
		if (!readobs(obsfile, 1, &popt_, &obss, &stas.data[i])) {			
			continue;
		}
		strcpy(stas.data[i].name, staname);

		setpcv(obss.n > 0 ? obss.data[0].time : timeget(), &popt_, &navs, &pcvss, &pcvsr, &stas.data[i]);
		readblq(fopt->blq, stas.data[i].name, popt->odisp[0]);

		matcpy(popt_.ru, stas.data[i].pos, 3, 1);
		pppcorr_stat(stas.data[i].name, popt_.ru);
		
		/* write outfile header */	
		if (!outhead(fopt->outposs, &popt_, sopt)) {		
			freeobs(&obss);
			continue;
		}
		outppp(fopt->outtrop, &popt_, &stas.data[i].name, &stas.data[i].pos, 'T');
		outppp(fopt->outiono, &popt_, &stas.data[i].name, &stas.data[i].pos, 'I');
		outppp(fopt->outfcbs, &popt_, &stas.data[i].name, &stas.data[i].pos, 'F');
		outppp(fopt->outelse, &popt_, &stas.data[i].name, &stas.data[i].pos, 'E');

		/* process */
		if (popt_.soltype == 0) {			
			procpos(&popt_, sopt, fopt, &rtk, 0);				
			pppfree(&rtk);
		}

		freeobs(&obss);
		checkbrk("Time=%.1f s\n", (tickget() - tick)*0.001);
	}	
	pppcorr_free();
	freedata(&navs, &pcvss, &pcvsr, &stas);	

	return 1;
}

