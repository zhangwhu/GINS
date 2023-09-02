
#include "rtklib.h"

#include "rtklib.h"

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define SQRT(x) ((x) <= 0.0 || (x) != (x) ? 0.0 : sqrt(x))

#define MAXINFILE 1000				 /* max number of input files */
#define TTOL_MOVEB (1.0 + 2 * DTTOL) /* time sync tolerance for moving-baseline (s) */

/* constants/global variables ------------------------------------------------*/
static pcvs_t pcvss = {0}; /* receiver antenna parameters */
static pcvs_t pcvsr = {0}; /* satellite antenna parameters */
static obs_t obss = {0};   /* observation data */
static nav_t navs = {0};   /* navigation data */
static stas_t stas = {0};  /* station list */
static int nepoch = 0;	   /* number of observation epochs */
static int iobsu = 0;	   /* current rover observation data index */
static int iobsr = 0;	   /* current reference observation data index */
static int revs = 0;	   /* analysis direction (0:forward,1:backward) */
static FILE *fp_sol = NULL;
static FILE *fp_rel = NULL;

static FILE *openfile(const char *file)
{
	gtime_t time = utc2gpst(timeget());
	FILE *fp;
	char path[1024];

	trace(3, "openstat: file=%s\n", file);

	reppath(file, path, time, "", "");

	createdir(path);
	if (!(fp = fopen(path, "w")))
	{
		trace(1, "openstat: file open error path=%s\n", path);
		return NULL;
	}
	return fp;
}
static void closefile(FILE *fp)
{
	if (fp)
		fclose(fp);
	fp = NULL;
}

static int setsta(const stas_t *stas, sta_t *sta)
{
	int i, ret = 0;
	for (i = 0; i < stas->n; i++)
	{
		if (!strnicmp(stas->data[i].name, sta->name, 4))
		{
			matcpy(sta->pos, stas->data[i].pos, 3, 1);
			ret = 1;
			break;
		}
	}
	return ret;
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

static int process(rtk_t *rel, const obsd_t *obs, int n, const nav_t *nav)
{
	prcopt_t *opt = &rel->opt;
	sol_t solb = {{0}};
	gtime_t time = rel->sol.time;
	int i, nu, nr;
	char msg[128] = "";

	trace(3, "process  : time=%s n=%d\n", time_str(obs[0].time, 3), n);

	/* count rover/base station observations */
	for (nu = 0; nu < n && obs[nu].rcv == 1; nu++)
		;
	for (nr = 0; nu + nr < n && obs[nu + nr].rcv == 2; nr++)
		;

	/* set base staion position */
	if (opt->mode != PMODE_MOVEB)
	{
		for (i = 0; i < 6; i++)
			rel->rb[i] = i < 3 ? opt->rb[i] : 0.0;
	}

	/* rover position by single point positioning */
	if (!pntpos(obs, nu, nav, &rel->opt, &rel->sol, NULL, rel->ssat, msg))
	{
		showmsg(msg);
		return 0;
	}
	if (time.time != 0)
		rel->tt = timediff(rel->sol.time, time);

	/* check number of data of base station and age of differential */
	if (nr == 0)
	{
		trace(1, "no base station observation data for rtk\n");
		return 1;
	}

	if (opt->mode == PMODE_MOVEB)
	{ /*  moving baseline */

		/* estimate position/velocity of base station */
		if (!pntpos(obs + nu, nr, nav, &rel->opt, &solb, NULL, NULL, msg))
			return 0;
		rel->sol.age = (float)timediff(rel->sol.time, solb.time);

		if (fabs(rel->sol.age) > TTOL_MOVEB)
			return 0;
		for (i = 0; i < 6; i++)
			rel->rb[i] = solb.rr[i];
	}
	else
	{
		rel->sol.age = (float)timediff(obs[0].time, obs[nu].time);

		if (fabs(rel->sol.age) > opt->maxtdiff)
			return 1;
	}

	/* relative potitioning */
	relpos(rel, obs, nu, nr, nav);
	return 1;
}

/* process positioning ----------------------------------------------------*/
static void procpos(const prcopt_t *popt, const solopt_t *sopt, rtk_t *rel, int mode)
{
	int i, nobs, n;
	obsd_t obs[MAXOBS];

	if (mode == 0)
		iobsu = 0;
	else
		iobsu = obss.n - 1;
	rtkinit(rel, popt);

	while ((nobs = inputobs(obs, rel->sol.stat, popt)) >= 0)
	{
		/* exclude satellites */
		for (i = n = 0; i < nobs; i++)
		{
			if ((satsys(obs[i].sat, NULL) & popt->navsys) &&
				popt->exsats[obs[i].sat - 1] != 1)
				obs[n++] = obs[i];
		}
		if (n <= 0)
			continue;

		if (!process(rel, obs, n, &navs))
			continue;
		if (mode == 0)
		{
			if (rel->sol.stat == SOLQ_FLOAT || rel->sol.stat == SOLQ_FIX || rel->sol.stat == SOLQ_SINGLE)
			{
				outsol(fp_sol, &rel->sol, rel->rb, sopt);
				outrel(fp_rel, rel);
			}
		}
	}
	checkbrk("                                        ");
}

int main(int argc, char *argv[])
{
	sta_t stas[2] = {0};
	rtk_t rel = {0};
	prcopt_t prcopt = prcopt_default;
	solopt_t solopt = solopt_default;
	filopt_t filopt = {""};
	unsigned int tick = tickget();
	char solfile[MAXINFILE], relfile[MAXINFILE];

	/* load conf */
	if (!loadopts(argv[1], sysopts))
		return -1;
	getsysopts(&prcopt, &solopt, &filopt);

	/* read product */
	if (!readproduct(&prcopt, &filopt, &navs, &pcvss, &pcvsr, NULL))
		return 0;

	/* read obs */
	if (!readobs(filopt.rovobs, 1, &prcopt, &obss, stas, &nepoch))
		return 0;
	if (!readobs(filopt.refobs, 2, &prcopt, &obss, stas + 1, &nepoch))
		return 0;

	/* set antenna &&  ocean tide && ref position */
	setpcv(obss.n > 0 ? obss.data[0].time : timeget(), &prcopt, &navs, &pcvss, &pcvsr, stas);
	readblq(filopt.blq, stas[0].name, prcopt.odisp[0]);
	readblq(filopt.blq, stas[1].name, prcopt.odisp[1]);

	/* open outfile */
	// sprintf(solfile, "%s%s%s", filopt.outdir, stas[0].name, "-%y-%m-%d.sol");
	// sprintf(relfile, "%s%s%s", filopt.outdir, stas[0].name, "-%y-%m-%d.rtk");
	if (!(fp_sol = openfile(solfile)))
		return 0;
	if (!(fp_rel = openfile(relfile)))
		return 0;

	/* write outfile header */
	outheader(fp_sol, &prcopt, &solopt);

	/* process */
	if (prcopt.soltype == 0)
	{
		procpos(&prcopt, &solopt, &rel, 0);
		rtkfree(&rel);
	}
	closefile(fp_sol);
	closefile(fp_rel);
	freeobs(&obss);

	freeproduct(&navs, &pcvss, &pcvsr, NULL);

	checkbrk("Time=%.1f s\n", (tickget() - tick) * 0.001);
	getchar();
	return 1;
}