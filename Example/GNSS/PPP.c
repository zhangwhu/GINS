
#include "rtklib.h"

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define SQRT(x) ((x) <= 0.0 || (x) != (x) ? 0.0 : sqrt(x))

/* constants/global variables ------------------------------------------------*/
static pcvs_t pcvss = {0};	/* receiver antenna parameters */
static pcvs_t pcvsr = {0};	/* satellite antenna parameters */
static obs_t obss = {0};	/* observation data */
static nav_t navs = {0};	/* navigation data */
static stas_t stas = {{0}}; /* station list */
static int nepoch = 0;		/* number of observation epochs */
static int iobsu  = 0;		/* current rover observation data index */
static int iobsr  = 0;		/* current reference observation data index */
static int revs   = 0;		/* analysis direction (0:forward,1:backward) */

static FILE *openfile(const char *file, gtime_t time, char *name)
{
	FILE *fp;
	char path[MAXEXFILE];

	reppath(file, path, time, name, "");
	createdir(path);
	if (!(fp = fopen(path, "w")))
	{
		trace(1, "openstat: file open error path=%s\n", path);
		return NULL;
	}
	return fp;
}
static void closefile(FILE *fp[], int n)
{
	for (int i = 0; i < n; i++)
	{
		if (fp[i])
			fclose(fp[i]);
		fp[i] = NULL;
	}
}

/* search next observation data index ----------------------------------------*/
static int nextobsf(const obs_t *obs, int *i, int rcv)
{
	int n;
	double tt;

	for (; *i < obs->n; (*i)++)
		if (obs->data[*i].rcv == rcv)
			break;
	for (n = 0; *i + n < obs->n; n++)
	{
		tt = timediff(obs->data[*i + n].time, obs->data[*i].time);
		if (obs->data[*i + n].rcv != rcv || tt > DTTOL)
			break;
	}
	return n;
}
static int nextobsb(const obs_t *obs, int *i, int rcv)
{
	double tt;
	int n;

	for (; *i >= 0; (*i)--)
		if (obs->data[*i].rcv == rcv)
			break;
	for (n = 0; *i - n >= 0; n++)
	{
		tt = timediff(obs->data[*i - n].time, obs->data[*i].time);
		if (obs->data[*i - n].rcv != rcv || tt < -DTTOL)
			break;
	}
	return n;
}

/* input obs data, navigation messages and sbas correction -------------------*/
static int inputobs(obsd_t *obs, int solq, const prcopt_t *popt)
{
	gtime_t time = {0};
	int i, nu, nr, n = 0;

	trace(3, "infunc  : revs=%d iobsu=%d\n", revs, iobsu);

	if (0 <= iobsu && iobsu < obss.n)
	{
		settime((time = obss.data[iobsu].time));
		if (checkbrk("processing : %s Q=%d", time_str(time, 0), solq))
		{
			showmsg("aborted");
			return -1;
		}
	}
	if (!revs)
	{ /* input forward data */
		if ((nu = nextobsf(&obss, &iobsu, 1)) <= 0)
			return -1;
		if (popt->intpref)
		{
			for (; (nr = nextobsf(&obss, &iobsr, 2)) > 0; iobsr += nr)
				if (timediff(obss.data[iobsr].time, obss.data[iobsu].time) > -DTTOL)
					break;
		}
		else
		{
			for (i = iobsr; (nr = nextobsf(&obss, &i, 2)) > 0; iobsr = i, i += nr)
				if (timediff(obss.data[i].time, obss.data[iobsu].time) > DTTOL)
					break;
		}
		nr = nextobsf(&obss, &iobsr, 2);
		if (nr <= 0)
		{
			nr = nextobsf(&obss, &iobsr, 2);
		}
		for (i = 0; i < nu && n < MAXOBS * 2; i++)
			obs[n++] = obss.data[iobsu + i];
		for (i = 0; i < nr && n < MAXOBS * 2; i++)
			obs[n++] = obss.data[iobsr + i];
		iobsu += nu;
	}
	else
	{ /* input backward data */
		if ((nu = nextobsb(&obss, &iobsu, 1)) <= 0)
			return -1;
		if (popt->intpref)
		{
			for (; (nr = nextobsb(&obss, &iobsr, 2)) > 0; iobsr -= nr)
				if (timediff(obss.data[iobsr].time, obss.data[iobsu].time) < DTTOL)
					break;
		}
		else
		{
			for (i = iobsr; (nr = nextobsb(&obss, &i, 2)) > 0; iobsr = i, i -= nr)
				if (timediff(obss.data[i].time, obss.data[iobsu].time) < -DTTOL)
					break;
		}
		nr = nextobsb(&obss, &iobsr, 2);
		for (i = 0; i < nu && n < MAXOBS * 2; i++)
			obs[n++] = obss.data[iobsu - nu + 1 + i];
		for (i = 0; i < nr && n < MAXOBS * 2; i++)
			obs[n++] = obss.data[iobsr - nr + 1 + i];
		iobsu -= nu;
	}
	return n;
}

static int process(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	prcopt_t *opt = &rtk->opt;
	gtime_t time;
	char msg[128] = "";
	trace(3, "rtkpos  : time=%s n=%d\n", time_str(obs[0].time, 3), n);

	time = rtk->sol.time;			

	/* rover position by single point positioning */
	if (!pntpos(obs, n, nav, &rtk->opt, &rtk->sol, NULL, rtk->ssat, msg))			
	{
		rtk->sol.time = time;			
		showmsg(msg);
		return 0;
	}
	if (time.time != 0) rtk->tt = fabs(timediff(rtk->sol.time, time));		

	/* precise point positioning */
	if (opt->mode >= PMODE_PPP_KINEMA)
	{
		ppppos(rtk, obs, n, nav, NULL, NULL);
		return 1;
	}
}

/* process positioning ----------------------------------------------------*/
static void procpos(const prcopt_t *popt, const solopt_t *sopt, rtk_t *rtk,
					int mode, FILE *fp_outs[])
{
	int i, nobs, n;
	obsd_t obs[MAXOBS];

	iobsu = (mode == 0) ? 0 : obss.n - 1;
	rtkinit(rtk, popt);

	while ((nobs = inputobs(obs, rtk->sol.stat, popt)) >= 0)	
	{
		/* exclude satellites */
		for (i = n = 0; i < nobs; i++)									
		{
			if (!(rtk->ssat[obs[i].sat - 1].sys&popt->navsys)) continue;
			if (popt->exsats[obs[i].sat - 1]) continue;
			if (ISZERO(obs[i].P[0]) || ISZERO(obs[i].L[0])) continue;
			obs[n++] = obs[i];
		}
		if (n <= 0) continue;
		
		if (!process(rtk, obs, n, &navs)) continue;				
		if (mode == 0) {	
			outppp(fp_outs[1], rtk, obs, n, &navs);
			outsol(fp_outs[0], &rtk->sol, rtk->opt.ru, sopt);
		}
	}
	checkbrk("                                        ");
}

int main(int argc, char *argv[])
{
	sta_t sta = {0};
	rtk_t rtk = {0};
	prcopt_t prcopt = prcopt_default;
	solopt_t solopt = solopt_default;
	filopt_t filopt = {""};
	FILE *fp_outs[2] = {NULL}, *fp_stas = NULL;
	char path[MAXEXFILE], outfile1[MAXEXFILE], outfile2[MAXEXFILE], outfile3[MAXEXFILE];

	if(!(fp_stas = fopen("fcb-enu-0522.txt","w"))) return 0;

	/* load conf */
	resetsysopts();
	if (!loadopts(argv[1], sysopts)) return -1;
	getsysopts(&prcopt, &solopt, &filopt);
	/* read product */
	if (!readproduct(&prcopt, &filopt, &navs, &pcvss, &pcvsr, &stas)) return 0;

	for (int i = 0; i < stas.n; i++)	
	{
		unsigned int tick = tickget();
		trace(1, "%03d:%s\n", i + 1, stas.data[i].name);

		/* read obs */
		reppath(filopt.rovobs, path, prcopt.ts, stas.data[i].name, "");
		if (!readobs(path, 1, &prcopt, &obss, &sta, &nepoch)) continue;

		/* read ppp corrections */
		reppath(filopt.corr, path, prcopt.ts, "", "");
		pppcorr_read(path, &navs);

		/* set antenna &&  ocean tide && ref position */
		setpcv(obss.n > 0 ? obss.data[0].time : timeget(), &prcopt, &navs, &pcvss, &pcvsr, &sta);		
		readblq(filopt.blq, sta.name, prcopt.odisp[0]);
		matcpy(prcopt.ru, (norm(stas.data[i].pos, 3) > 0.0) ? stas.data[i].pos : sta.pos, 3, 1);

		/* open outfile &&  write file header*/
		sprintf(outfile1, "%s%s", filopt.outdir, filopt.outfile1);
		sprintf(outfile2, "%s%s", filopt.outdir, filopt.outfile2);
		sprintf(outfile3, "%s%s%s", filopt.outdir, stas.data[i].name,filopt.trace);			
		if (!(fp_outs[0] = openfile(outfile1, prcopt.ts, stas.data[i].name)) ||
			!(fp_outs[1] = openfile(outfile2, prcopt.ts, stas.data[i].name)))
		{
			freeobs(&obss);		closefile(fp_outs, 2);		continue;
		}
		traceopen(outfile3);	tracelevel(1);
		outheader(fp_outs[0], &prcopt, &solopt);
		fprintf(fp_outs[1], "$NAME,%s\n", stas.data[i].name);

		/* process */
		if (prcopt.soltype == 0)
		{
			procpos(&prcopt, &solopt, &rtk, 0, fp_outs);
#if 1
			char time[32];
			double rr[3], enu[3], pos[3];
			ecef2pos(prcopt.ru, pos);
			for (int ii=0; ii<3; ii++) rr[ii] = rtk.sol.rr[ii] - prcopt.ru[ii];
			ecef2enu(pos, rr, enu);
			time2str(rtk.sol.time, time, 0);
			fprintf(fp_stas, "%s %25s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %8s%s\n", stas.data[i].name, sta.antdes, 
				sta.del[0], sta.del[1],	sta.del[2], enu[0], enu[1], enu[2], " ", time);
			fflush(fp_stas);
#else
			fprintf(fp_stas, "%s %15.4f %15.4f %15.4f\n", stas.data[i].name, rtk.sol.rr[0], rtk.sol.rr[1], rtk.sol.rr[2]);
			fflush(fp_stas);
#endif
			rtkfree(&rtk);
		}

		freeobs(&obss);
		closefile(fp_outs, 2);	traceclose();
		checkbrk("Time=%.1f s\n", (tickget() - tick) * 0.001);
	}
	freeproduct(&navs, &pcvss, &pcvsr, &stas);
	// sol2kml(filopt.outdir);
	getchar();
	return 1;
}
