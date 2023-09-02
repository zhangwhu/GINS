#include "rtklib.h"

#define MAXRCV 300
#define MAXARC 15
#define MAXITR 10
#define MAXIND 192
#define MAXOBS MAXRCV * 40

#define THRESWL 0.05 // 0.1>0.05
#define THRESEL 0.02
#define THRESCN 50
#define THRESIN 0.25

#define NS 5
#define NR 4

#define SQR(x) ((x) * (x))
#define SQRT(x) ((x) <= 0.0 || (x) != (x) ? 0.0 : sqrt(x))
#define ROUND(x) (int)floor((x) + 0.5)
#define FROUND(x) ((x) - (int)floor((x) + 0.5))
#define IFROUND(x) ((x) < 0 ? (x + 1) : (x))

#define IS(x) (x)
#define IR(sys, x) (x + MAXSAT + sys * MAXRCV)
#define NX (MAXSAT + 4 * MAXRCV)

typedef struct
{				/* ambiguity control type */
	gtime_t ts; /* start time */
	int n;		/* number of epochs */
	double LC;	/* linear combination average(m) */
	double LCv; /* linear combination variance */
} ambd_t;

typedef struct
{
	gtime_t ts, te;
	int rcv, sat, lenArc, NN, flag;
	double val, std;
} arcInf;
typedef struct
{
	int n, nmax, nsta;
	unsigned int narc[MAXRCV][MAXSAT], index[MAXRCV][MAXSAT][MAXARC];
	char name[MAXRCV][MAXANT];
	uint8_t flag[MAXRCV];
	arcInf *data;
} arc_t;
typedef struct
{
	int ind, rcv, sat, len;
	double val, std;
} ambInf;
typedef struct
{
	int n, nmax;
	ambInf *data;
} amb_t;
typedef struct
{
	unsigned char flag[NX];
	double ini[NX], sol[NX];
} fcb_t;
static FILE *fp;
static ambd_t ambd[MAXSAT] = {0};
static arc_t arc_el = {0}, arc_wl = {0};
static amb_t amb_el = {0}, amb_wl = {0}, amb_nl = {0};
static fcb_t fcb_el = {0}, fcb_wl = {0}, fcb_nl = {0};

static int sat2sys(int sat)
{
	static int listsys[4] = {SYS_GPS, SYS_GLO, SYS_GAL, SYS_CMP};
	int sys, id = -1;
	sys = satsys(sat, NULL);
	for (int k = 0; k < 4; k++)
	{
		if (sys & listsys[k])
			id = k;
	}
	return id;
}
static void outarc(arc_t *arcs, char *file)
{
	FILE *fp;
	int i, j, k, rcv, sat, sys;
	double BN, NN;
	char id[32];
	if (!(fp = fopen(file, "w")))
		return;
	for (i = 0; i < arcs->n; i++)
	{
		rcv = arcs->data[i].rcv;
		sat = arcs->data[i].sat;
		satno2id(sat + 1, id);
		if (!(sys = satsys(sat + 1, NULL)))
			continue;
		fprintf(fp, "%s\t", time_str(arcs->data[i].ts, 0));
		fprintf(fp, "%s %s %s %10.4f %10.4f %4d\n", time_str(arcs->data[i].te, 0), arcs->name[rcv], id, arcs->data[i].val,
				arcs->data[i].std, arcs->data[i].lenArc);
		if (arcs->data[i].std > THRESWL)
		{
			trace(1, "%s\t", time_str(arcs->data[i].ts, 0));
			trace(1, "%s %s %s %10.4f %10.4f %4d\n", time_str(arcs->data[i].te, 0), arcs->name[rcv], id, arcs->data[i].val,
				  arcs->data[i].std, arcs->data[i].lenArc);
		}
	}
	fclose(fp);
}
static void outamb(amb_t *ambs, arc_t *arcs, char *file)
{
	FILE *fp;
	int i, j, k, rcv, sat, sys, n;
	double BN, NN;
	char id[32];
	if (!(fp = fopen(file, "w")))
		return;
	for (i = 0; i < ambs->n; i++)
	{
		rcv = ambs->data[i].rcv;
		sat = ambs->data[i].sat;
		satno2id(sat + 1, id);
		if (!(sys = satsys(sat + 1, NULL)))
			continue;
		fprintf(fp, "%s %s %4d %10.4f %10.4f\n", arcs->name[rcv], id, ambs->data[i].ind, ambs->data[i].val, ambs->data[i].std);
	}
	fclose(fp);
}
static void outarcpos(arc_t *arcs, fcb_t *fcb, char *file)
{
	FILE *fp;
	int i, j, k, rcv, sat;
	double BN, NN;
	char id[32];
	if (!(fp = fopen(file, "w")))
		return;
	for (i = 0; i < arcs->n; i++)
	{
		rcv = arcs->data[i].rcv;
		sat = arcs->data[i].sat;
		if ((k = sat2sys(sat + 1)) < 0)
			continue;
		if (fcb->flag[sat] && fcb->flag[IR(k, rcv)])
		{
			BN = arcs->data[i].val - (fcb->ini[sat] - fcb->ini[IR(k, rcv)]);
			NN = arcs->data[i].val - (fcb->sol[sat] - fcb->sol[IR(k, rcv)]);
			satno2id(sat + 1, id);
			fprintf(fp, "%s\t", time_str(arcs->data[i].ts, 0));
			fprintf(fp, "%s %s %d %s %10.4f %10.4f\n", time_str(arcs->data[i].te, 0), arcs->name[rcv], k, id, FROUND(BN), FROUND(NN));
		}
	}
	fclose(fp);
}
static void outambpos(ambInf *ambs, int n, fcb_t *fcb, arc_t *arcs, char *file)
{
	FILE *fp;
	int i, j, k, rcv, sat;
	double BN, NN;
	char id[32];
	if (!(fp = fopen(file, "a")))
		return;
	for (i = 0; i < n; i++)
	{
		rcv = ambs[i].rcv;
		sat = ambs[i].sat;
		if ((k = sat2sys(sat + 1)) < 0)
			continue;
		// if (k != 3)
		// 	continue;
		if (fcb->flag[sat] && fcb->flag[IR(k, rcv)])
		{
			BN = ambs[i].val - (fcb->ini[sat] - fcb->ini[IR(k, rcv)]);
			NN = ambs[i].val - (fcb->sol[sat] - fcb->sol[IR(k, rcv)]);
			satno2id(sat + 1, id);
			fprintf(fp, "%4d %s %s %d %10.4f %10.4f\n", ambs[i].ind, arcs->name[rcv], id, k, FROUND(BN), FROUND(NN));
		}
	}
	fclose(fp);
}
static int time2ind(const gtime_t ts, gtime_t t, double nhour)
{
	double n;
	if (timediff(t, ts) >= 0.0f)
	{
		n = (int)(timediff(t, ts) / (nhour * 3600));
		return n;
	}
	return -1;
}
static gtime_t ind2time(const gtime_t ts, int ind, double nhour)
{
	return timeadd(ts, ind * nhour * 3600);
}
static int rejectMaxRes(double *a, int *flag, double b, int n, double thread)
{
	int pos = 0;
	double temp = 0.0;
	for (int i = 0; i < n; i++)
	{
		if (a[i] == 0)
			continue;
		if (temp <= fabs(a[i] - b))
		{
			temp = fabs(a[i] - b);
			pos = i;
		}
	}
	if (temp <= thread)
		return 0.0;
	else
	{
		a[pos] = 0.0;
		flag[pos] = 1;
		return temp;
	}
}
static double getAverage(double *a, int n, double *std, int ext)
{
	if (n == 0)
		return 0.0;
	int i, j, k = n, *flag;
	double *aTmp, *aTmp1, *aTmp2, aveTmp[2] = {0}, stdTmp[2] = {0}, ave, sum, cov = 0;

	aTmp1 = mat(n, 1);
	aTmp2 = mat(n, 1);
	flag = imat(n, 1);
	for (i = 0; i < n; i++)
	{
		flag[i] = 0;
		aTmp1[i] = FROUND(a[i]);
		aTmp2[i] = IFROUND(FROUND(a[i]));
		aveTmp[0] += (aTmp1[i] - aveTmp[0]) / (i + 1);
		aveTmp[1] += (aTmp2[i] - aveTmp[1]) / (i + 1);
	}

	for (i = 0; i < n; i++)
	{
		stdTmp[0] += pow(aTmp1[i] - aveTmp[0], 2);
		stdTmp[1] += pow(aTmp2[i] - aveTmp[1], 2);
	}
	if (stdTmp[0] > stdTmp[1])
	{
		ave = aveTmp[1];
		aTmp = aTmp2;
	}
	else
	{
		ave = aveTmp[0];
		aTmp = aTmp1;
	}

	if (ext)
	{
		while (1)
		{
			if (k <= 15 || rejectMaxRes(aTmp, flag, ave, n, 0.15) == 0)
			{
				break;
			}
			else
			{
				k--;
				sum = 0.0;
				for (i = 0; i < n; i++)
					sum += aTmp[i];
				ave = sum / k;
			}
		}
	}
	if (std)
	{
		for (i = 0; i < n; i++)
		{
			if (!flag[i])
			{
				cov += pow(aTmp[i] - ave, 2);
			}
		}
		*std = sqrt(cov / (k - 1));
	}
	free(aTmp1);
	free(aTmp2);
	free(flag);
	return FROUND(ave);
}
static void initarc(arcInf *arc, int rcv, int sat, gtime_t ts, gtime_t te, int lenArc, double val, double std)
{
	arc->rcv = rcv;
	arc->sat = sat;
	arc->ts = ts;
	arc->te = te;
	arc->val = val;
	arc->std = std;
	arc->lenArc = lenArc;
}
static int addarc(arc_t *arcs, arcInf *data)
{
	int rcv, sat;
	arcInf *arc_data;
	if (arcs->nmax <= arcs->n)
	{
		if (arcs->nmax <= 0)
			arcs->nmax = 9000;
		else
			arcs->nmax *= 2;
		if (!(arc_data = (arcInf *)realloc(arcs->data, sizeof(arcInf) * arcs->nmax)))
		{
			trace(1, "addarcdata: memalloc error n=%dx%d\n", sizeof(arcInf), arcs->nmax);
			free(arcs->data);
			arcs->data = NULL;
			arcs->n = arcs->nmax = 0;
			return -1;
		}
		arcs->data = arc_data;
	}
	rcv = data->rcv;
	sat = data->sat;
	arcs->index[rcv][sat][arcs->narc[rcv][sat]++] = arcs->n;
	arcs->data[arcs->n++] = *data;
	return 1;
}
static void freearc(arc_t *arcs)
{
	free(arcs->data);
	arcs->data = NULL;
	arcs->n = arcs->nmax = 0;
}
static void initamb(ambInf *amb, int ind, int rcv, int sat, double val, double std, int len)
{
	amb->ind = ind;
	amb->rcv = rcv;
	amb->sat = sat;
	amb->val = val;
	amb->std = std;
	amb->len = len;
}
static int addamb(amb_t *amb, ambInf *data)
{
	ambInf *amb_data;
	if (amb->nmax <= amb->n)
	{
		if (amb->nmax <= 0)
			amb->nmax = 5000;
		else
			amb->nmax *= 2;
		if (!(amb_data = (ambInf *)realloc(amb->data, sizeof(ambInf) * amb->nmax)))
		{
			trace(1, "addambdata: memalloc error n=%dx%d\n", sizeof(ambInf), amb->nmax);
			free(amb->data);
			amb->data = NULL;
			amb->n = amb->nmax = 0;
			return -1;
		}
		amb->data = amb_data;
	}
	amb->data[amb->n++] = *data;
	return 1;
}
static void freeamb(amb_t *ambs)
{
	free(ambs->data);
	ambs->data = NULL;
	ambs->n = ambs->nmax = 0;
}
static int isRefSat(int sat, int *refSat, int nRefSat)
{
	for (int i = 0; i < nRefSat; i++)
	{
		if (sat == refSat[i])
			return 1;
	}
	return 0;
}
static void arc2amb(arc_t *arcs, amb_t *ambs)
{
	ambInf amb0 = {0};
	int i, j, k, narc, len;
	double std, val[200];
	for (i = 0; i < arcs->nsta; i++)
	{
		for (j = 0; j < MAXSAT; j++)
		{
			narc = arcs->narc[i][j];
			if (narc > 0)
			{
				for (k = len = 0, std = 1E4; k < narc && k < 200; k++)
				{
					if (std > arcs->data[arcs->index[i][j][k]].std)
					{
						std = arcs->data[arcs->index[i][j][k]].std;
					}
					if (len < arcs->data[arcs->index[i][j][k]].lenArc)
					{
						len = arcs->data[arcs->index[i][j][k]].lenArc;
					}
					val[k] = arcs->data[arcs->index[i][j][k]].val;
				}
				amb0.rcv = i;
				amb0.sat = j;
				amb0.len = len;
				amb0.val = getAverage(val, k, NULL, 0);
				amb0.std = std;
				if (addamb(ambs, &amb0) < 0)
					return;
			}
		}
	}
}
static int sNum[MAXSAT], rNum[MAXRCV], len[MAXRCV][MAXSAT];
static double iniSatFcbSum[MAXSAT][MAXRCV], iniRcvFcbSum[MAXRCV][MAXSAT];
static double val[MAXRCV][MAXSAT], std[MAXRCV][MAXSAT];
static int inifcb(const ambInf *amb, int num, int *flag, double *ini)
{
	int i, j, k, exc, rcv, sat, ref, temp, nfix = 0, maxSat[MAXSAT], ind[MAXSAT], is[MAXSAT];
	double tempd, std0;

	for (i = 0; i < MAXRCV; i++)
		rNum[i] = 0;
	for (i = 0; i < MAXSAT; i++)
	{
		maxSat[i] = 0;
		ind[i] = 1E4;
	}
	for (i = 0; i < MAXSAT + MAXRCV; i++)
		flag[i] = 0;
	for (i = 0; i < MAXRCV; i++)
		for (j = 0; j < MAXSAT; j++)
			len[i][j] = 0;
	for (i = 0, temp = 0; i < num; i++)
	{
		rcv = amb[i].rcv;
		sat = amb[i].sat;
		val[rcv][sat] = amb[i].val;
		len[rcv][sat] = amb[i].len;
		std[rcv][sat] = amb[i].std;
		maxSat[sat]++;
		if (temp < maxSat[sat])
		{
			temp = maxSat[sat];
			ref = sat;
		}
	}
	flag[ref] = 1;
	ini[ref] = 0;
	ind[ref] = nfix;

	while (1)
	{
		for (sat = 0; sat < MAXSAT; sat++)
		{
			sNum[sat] = 0;
			is[sat] = sat;
		}
		for (rcv = 0; rcv < MAXRCV; rcv++)
		{
			for (sat = 0, temp = 1E2, ref = -1; sat < MAXSAT; sat++)
			{
				if (flag[sat] && len[rcv][sat] > 0 && ind[sat] < temp)
				{
					temp = ind[sat];
					ref = sat;
				}
			}
			if (ref < 0)
				continue;
			for (sat = 0; sat < MAXSAT; sat++)
			{
				if (!flag[sat] && len[rcv][sat] > 0)
				{
					iniSatFcbSum[sat][sNum[sat]++] = FROUND(val[rcv][sat] - val[rcv][ref] + ini[ref]);
				}
			}
		}

		for (i = 0; i < MAXSAT - 1; i++)
			for (j = i + 1; j < MAXSAT; j++)
			{
				if (sNum[is[j]] > sNum[is[i]])
				{
					temp = is[i];
					is[i] = is[j];
					is[j] = temp;
				}
			}
		for (i = exc = 0; i < MAXSAT; i++)
		{
			sat = is[i];
			if (sNum[sat] >= NS)
			{
				tempd = getAverage(iniSatFcbSum[sat], sNum[sat], &std0, 0); // sFCB
				if (std0 < THRESIN)
				{
					exc = 1;
					flag[sat] = 1;
					ini[sat] = tempd;
					ind[sat] = ++nfix;
					trace(2, "sat%4d=%10.4f %4d\n", sat + 1, tempd, nfix);
					break;
				}
			}
		}
		if (!exc)
			break;
	}

	for (rcv = 0; rcv < MAXRCV; rcv++)
	{
		for (sat = 0; sat < MAXSAT; sat++)
		{
			if (flag[sat] && len[rcv][sat] > 0)
			{
				iniRcvFcbSum[rcv][rNum[rcv]++] = FROUND(ini[sat] - val[rcv][sat]); // NN+sFCB-rFCB
			}
		}
	}
	for (rcv = 0; rcv < MAXRCV; rcv++)
	{
		if (rNum[rcv] >= NR)
		{
			tempd = getAverage(iniRcvFcbSum[rcv], rNum[rcv], &std0, 0); // rFCB
			if (1 /*std0 < THRESIN*/)
			{
				flag[MAXSAT + rcv] = 1;
				ini[MAXSAT + rcv] = tempd;
				trace(2, "rcv%4d=%10.4f\n", rcv, tempd);
			}
		}
	}
	return 1;
}
static int lsqfcb(const ambInf *amb, int num, int *flag, double *ini, double *sol)
{
	int i, j, k, nrcv, nsat, nx, nv, iter, info, stat, nfix, ircv[MAXRCV], isat[MAXSAT];
	double BN, temp, *x0, *x, *dx, *v, *H, *Q;

	for (i = 0; i < MAXRCV; i++)
		for (j = 0; j < MAXSAT; j++)
			len[i][j] = 0;
	for (i = 0, nsat = 0; i < MAXSAT; i++)
		if (flag[i])
			isat[nsat++] = i;
	for (i = 0, nrcv = 0; i < MAXRCV; i++)
		if (flag[MAXSAT + i])
			ircv[nrcv++] = i;

	for (i = 0; i < num; i++)
	{
		val[amb[i].rcv][amb[i].sat] = amb[i].val;
		len[amb[i].rcv][amb[i].sat] = amb[i].len;
		std[amb[i].rcv][amb[i].sat] = amb[i].std;
	}
	nx = nsat + nrcv;
	x0 = mat(nx, 1);
	x = zeros(nx, 1);
	dx = zeros(nx, 1);
	v = zeros(num + 1, 1);
	H = zeros(num + 1, nx);
	Q = zeros(nx, nx);

	for (i = 0; i < nsat; i++)
		x0[i] = ini[isat[i]];
	for (i = 0; i < nrcv; i++)
		x0[nsat + i] = ini[MAXSAT + ircv[i]];

	for (iter = nfix = stat = 0; iter < MAXITR; iter++)
	{
		nv = 0;
		for (i = 0; i < nrcv; i++)
		{
			for (j = 0; j < nsat; j++)
			{
				if (len[ircv[i]][isat[j]] > 0)
				{
					BN = val[ircv[i]][isat[j]] - (x0[j] - x0[nsat + i]);
					if (fabs(FROUND(BN) < 0.25))
					{
						temp = (val[ircv[i]][isat[j]] - ROUND(BN) - (x[j] - x[nsat + i]));
						v[nv] = temp;
						H[j + nx * nv] = 1;
						H[nsat + i + nx * nv] = -1;
						nv++;
					}
				}
			}
		}
		for (i = 0, temp = 0; i < nsat; i++)
		{
			H[i + nx * nv] = 1 / 1e-4;
			temp += x[i];
		}
		v[nv++] = (0.0 - temp) / 1e-4;

		if (nfix >= nv && nv > 1)
		{
			stat = 1;
			break;
		}
		nfix = nv;
		/* least square estimation */
		if ((info = lsq(H, v, nx, nv, dx, Q)))
		{
			trace(1, "lsq error info=%d  nv=%d  nsat=%d nrcv=%d\n", info, nv, nsat, nrcv);
			break;
		}
		for (i = 0; i < nx; i++)
		{
			x[i] = x[i] + dx[i];
			x0[i] = x[i];
		}
		trace(2, "(%d) %10.4f obsN=%4d nrcv=%4d nsat=%4d\n", iter, norm(dx, nx), nv, nrcv, nsat);
	}
	for (i = 0; i < nsat; i++)
		sol[isat[i]] = x[i];
	for (i = 0; i < nrcv; i++)
		sol[MAXSAT + ircv[i]] = x[nsat + i];

	free(x0);
	free(x);
	free(dx);
	free(v);
	free(H);
	free(Q);
	return stat;
}
static void estfcb(const ambInf *amb, int *num, fcb_t *fcb)
{
	int i, k, n, sys, flag[MAXSAT + MAXRCV] = {0};
	double ifcb[MAXSAT + MAXRCV], sfcb[MAXSAT + MAXRCV];
	ambInf amb0[MAXOBS];

	for (i = 0; i < NX; i++)
		fcb->flag[i] = 0;
	for (k = 0; k < 4; k++)
	{
		for (i = 0, n = 0; i < num; i++)
		{
			if (sat2sys(amb[i].sat + 1) == k)
			{
				amb0[n++] = amb[i];
			}
		}
		if (n > 0)
		{
			if (inifcb(&amb0, n, flag, ifcb) && lsqfcb(&amb0, n, flag, ifcb, sfcb))
			{
				for (i = 0; i < MAXSAT; i++)
					if (flag[i])
					{
						fcb->flag[IS(i)] = flag[i];
						fcb->ini[IS(i)] = ifcb[i];
						fcb->sol[IS(i)] = sfcb[i];
					}
				for (i = 0; i < MAXRCV; i++)
					if (flag[MAXSAT + i])
					{
						fcb->flag[IR(k, i)] = flag[MAXSAT + i];
						fcb->ini[IR(k, i)] = ifcb[MAXSAT + i];
						fcb->sol[IR(k, i)] = sfcb[MAXSAT + i];
					}
			}
		}
	}
}
static void corr_fcb(arc_t *arc, const fcb_t *fcb)
{
	arcInf *arcd;
	int i, k, rcv, sat;
	double BN;

	for (i = 0; i < arc->n; i++)
	{
		rcv = arc->data[i].rcv;
		sat = arc->data[i].sat;
		arcd = &arc->data[i];
		if ((k = sat2sys(sat + 1)) < 0)
			continue;
		if (fcb->flag[IS(sat)] && fcb->flag[IR(k, rcv)])
		{
			BN = arcd->val - (fcb->sol[IS(sat)] - fcb->sol[IR(k, rcv)]);
			if (fabs(FROUND(BN)) < 0.25)
			{
				arcd->NN = ROUND(BN);
				arcd->flag = 1;
			}
		}
	}
}
/* compare observation data -------------------------------------------------*/
static int cmpamb(const void *p1, const void *p2)
{
	ambInf *q1 = (ambInf *)p1, *q2 = (ambInf *)p2;
	if (q1->ind != q2->ind)
		return (int)q1->ind - (int)q2->ind;
	if (q1->rcv != q2->rcv)
		return (int)q1->rcv - (int)q2->rcv;
	if (q1->sat != q2->sat)
		return (int)q1->sat - (int)q2->sat;
}
static int sortamb(amb_t *ambs)
{
	if (ambs->n <= 0)
		return 0;
	qsort(ambs->data, ambs->n, sizeof(ambInf), cmpamb);
	return 1;
}
static int inputamb(const amb_t *ambs, ambInf *data)
{
	static int iamb = 0;
	int i, n, ind;

	ind = ambs->data[iamb].ind;
	for (i = iamb, n = 0; i < ambs->n; i++)
	{
		if (ambs->data[i].ind != ind)
			break;
		data[n++] = ambs->data[i];
	}
	iamb += n;
	return n;
}
/* open output file for append -----------------------------------------------*/
static FILE *openfile(const char *outfile)
{
	trace(3, "openfile: outfile=%s\n", outfile);

	return !*outfile ? stdout : fopen(outfile, "a");
}
static int outhead(char *outfile, fcb_t *el, fcb_t *wl)
{
	FILE *fp;
	char id[32];
	if (outfile == NULL)
		return 0;

	createdir(outfile);
	if (!(fp = fopen(outfile, "w")))
	{
		showmsg("error : open output file %s", outfile);
		return 0;
	}
	/* output header */
	fprintf(fp, "Widelane Satellite Fractional Cycle Biases                  COMMENT\n");
	if (el)
	{
		for (int i = 0; i < MAXSAT; i++)
		{
			if (!el->flag[i])
				continue;
			satno2id(i + 1, id);
			fprintf(fp, "EL  %s%13.3f                                        COMMENT\n", id, el->sol[i]);
		}
	}
	if (wl)
	{
		for (int i = 0; i < MAXSAT; i++)
		{
			if (!wl->flag[i])
				continue;
			satno2id(i + 1, id);
			fprintf(fp, "WL  %s%13.3f                                        COMMENT\n", id, wl->sol[i]);
		}
	}
	fprintf(fp, "                                                            END OF HEADER\n");
	fclose(fp);
	return 1;
}
static void outbody(FILE *fp, const fcb_t *fcb, gtime_t ts, int ind)
{
	char id[32];
	gtime_t ti;
	ti = ind2time(ts, ind, 0.25);
	fprintf(fp, "* %s\n", time_str(ti, 5));
	for (int i = 0; i < MAXSAT; i++)
	{
		if (!fcb->flag[i])
			continue;
		satno2id(i + 1, id);
		fprintf(fp, "P%s%26.3f\n", id, FROUND(fcb->sol[i]));
	}
}
static void read_arc(gtime_t ts, gtime_t te, const char *file, arc_t *arcs, int type)
{
	gtime_t ti, t_1 = {0}, t_2;
	ambd_t ambd0 = {{0}, 0, 0, 1e4};
	arcInf arc0 = {0};

	int i, week, solq, slip, flg, sat, rcv, ns, stat = 0;
	int sflag[MAXSAT] = {0}, sflag_1[MAXSAT] = {0}, lock[3], num;
	double tow, azel, EL, Pe, WL, Pw, IF, Pf, NN, std, thread;
	char *p, buff[256] = "", satid[8] = "";

	rcv = arcs->nsta;
	ns = arcs->n;
	for (i = 0; i < MAXSAT; i++)
	{
		ambd[i] = ambd0;
		arcs->narc[rcv][i] = 0;
	}

	if (!(fp = fopen(file, "r")) || !(p = strrchr(file, FILEPATHSEP)))
	{
		trace(2, "input file open error: file=%s\n", file);
		return;
	}
	strncpy(arcs->name[rcv], p + 1, 4);
	arcs->name[rcv][4] = '\0';

	while (fgets(buff, sizeof(buff), fp))
	{
		if (buff[0] != '$')
			continue;
		if (!strncmp(buff + 1, "AMB", 3))
		{
			for (p = buff; *p; p++)
				if (*p == ',')
					*p = ' ';
			if (sscanf(buff, "$AMB%d%lf%d%s%d%lf%lf%lf%lf%lf%lf%lf%d%d%d\n", &week,
					   &tow, &solq, satid, &slip, &azel, &EL, &Pe, &WL, &Pw, &IF, &Pf,
					   &lock[0], &lock[1], &lock[2]) < 15)
			{
				continue;
			}
			ti = gpst2time(week, tow);
			sat = (unsigned char)satid2no(satid) - 1;

			if (timediff(ts, ti) > 0 || timediff(ti, te) > 0)
				continue;
			switch (type)
			{
			case 0:
				NN = WL;
				std = Pw;
				thread = THRESWL;
				flg = 0x3;
				num = (lock[0] > lock[1]) ? lock[1] : lock[0];
				break;
			case 1:
				NN = EL;
				std = Pe;
				thread = THRESEL;
				flg = 0x5;
				num = (lock[1] > lock[2]) ? lock[2] : lock[1];
				break;
			default:
				NN = 0;
			}

			if (timediff(ti, t_1) != 0.0)
			{
				for (i = 0; i < MAXSAT; i++)
				{
					if (sflag_1[i] && !sflag[i])
					{
						if (ambd[i].LCv < thread && ambd[i].n > THRESCN)
						{
							initarc(&arc0, rcv, i, ambd[i].ts, t_2, ambd[i].n, ambd[i].LC, ambd[i].LCv);
							if ((stat = addarc(arcs, &arc0)) < 0)
								break;
						}
						ambd[i].n = ambd[i].LC = 0;
						ambd[i].LCv = 1e4;
					}
					sflag_1[i] = sflag[i];
					sflag[i] = 0;
				}
				t_2 = t_1;
				t_1 = ti;
			}
			sflag[sat] = 1;
			if (slip & flg)
			{
				if (ambd[sat].LCv < thread && ambd[sat].n > THRESCN)
				{
					initarc(&arc0, rcv, sat, ambd[sat].ts, t_2, ambd[sat].n, ambd[sat].LC, ambd[sat].LCv);
					if ((stat = addarc(arcs, &arc0)) < 0)
						break;
				}
				ambd[sat].n = ambd[sat].LC = 0;
				ambd[sat].LCv = 1e4;
			}
			if (NN)
			{
				if (ambd[sat].n == 0)
				{
					ambd[sat].ts = ti;
					ambd[sat].LC = NN;
					ambd[sat].LCv = std;
				}
				if (ambd[sat].LCv > std)
				{
					ambd[sat].LC = NN;
					ambd[sat].LCv = std;
				}
				ambd[sat].n = num;
			}
		}
	}
	fclose(fp);
	for (i = 0; i < MAXSAT; i++)
	{
		if (sflag[i])
		{
			if (ambd[i].LCv < thread && ambd[i].n > THRESCN)
			{
				initarc(&arc0, rcv, i, ambd[i].ts, t_1, ambd[i].n, ambd[i].LC, ambd[i].LCv);
				if ((stat = addarc(arcs, &arc0)) < 0)
					break;
			}
		}
	}
	arcs->nsta++;
	if (arcs->n - ns < 10)
	{
		arcs->n = ns;
		arcs->nsta--;
	}
}
static void read_amb(gtime_t ts, gtime_t te, const char *file, arc_t *arc_wl, amb_t *amb_nl)
{
	gtime_t ti;
	ambInf amb0;
	char *p, buff[256], ymd[12] = "", dms[10] = "", satid[8] = "", time[32] = "", staname[32];
	int i, j, k, rcv, week, solq, slip, ind, sat, sys, id, stat, num[MAXIND][MAXSAT], exc[MAXIND][MAXSAT], lock[3], qq;
	double tow, azel, EL, Pe, WL, Pw, IF, Pf, NL, temp, NN, val[MAXIND][MAXSAT], std[MAXIND][MAXSAT], minstd[MAXIND][MAXSAT];

	for (i = 0; i < MAXIND; i++)
		for (j = 0; j < MAXSAT; j++)
		{
			minstd[i][j] = 1E4;
			exc[i][j] = 0;
		}

	if (!(fp = fopen(file, "r")) || !(p = strrchr(file, FILEPATHSEP)))
		return;
	for (i = 0, rcv = -1; i < arc_wl->nsta; i++)
	{
		if (!strncmp(arc_wl->name[i], p + 1, 4))
		{
			rcv = i;
			break;
		}
	}
	if (rcv < 0)
		return;

	while (fgets(buff, sizeof(buff), fp))
	{
		if (buff[0] != '$')
			continue;
		if (!strncmp(buff + 1, "AMB", 3))
		{
			for (p = buff; *p; p++)
				if (*p == ',')
					*p = ' ';
			if (sscanf(buff, "$AMB%d%lf%d%s%d%lf%lf%lf%lf%lf%lf%lf%d%d%d\n", &week,
					   &tow, &solq, satid, &slip, &azel, &EL, &Pe, &WL, &Pw, &IF, &Pf,
					   &lock[0], &lock[1], &lock[2]) < 15)
			{
				continue;
			}
			ti = gpst2time(week, tow);
			sat = (unsigned char)satid2no(satid) - 1;
			qq = (lock[0] > lock[1]) ? lock[1] : lock[0];
			if (timediff(ts, ti) > 0 || timediff(ti, te) > 0)
				continue;

			if (!(sys = satsys(sat + 1, NULL)) || qq < THRESCN || azel < 20)
				continue;
			if ((ind = time2ind(ts, ti, 0.25)) == -1)
				continue;
			for (k = 0; k < arc_wl->narc[rcv][sat]; k++)
			{
				id = arc_wl->index[rcv][sat][k];
				if (arc_wl->data[id].flag && (timediff(arc_wl->data[id].ts, ti) < 0) && (timediff(arc_wl->data[id].te, ti) >= 0))
				{
					switch (sys)
					{
					case SYS_GPS:
						NL = IF * (FREQ1 + FREQ2) / CLIGHT - arc_wl->data[id].NN * FREQ2 / (FREQ1 - FREQ2);
						break;
					case SYS_GAL:
						NL = IF * (FREQ1 + FREQ5) / CLIGHT - arc_wl->data[id].NN * FREQ5 / (FREQ1 - FREQ5);
						break;
					case SYS_CMP:
						NL = IF * (FREQ1_CMP + FREQ3_CMP) / CLIGHT - arc_wl->data[id].NN * FREQ3_CMP / (FREQ1_CMP - FREQ3_CMP);
						break;
					default:
						continue;
					}
					if (minstd[ind][sat] > Pf)
					{
						val[ind][sat] = NL;
						std[ind][sat] = minstd[ind][sat] = Pf;
						num[ind][sat] = qq;
						exc[ind][sat] = 1;
					}
				}
			}
		}
	}
	fclose(fp);
	for (i = 0; i < MAXIND; i++)
	{
		for (j = 0; j < MAXSAT; j++)
		{
			if (exc[i][j])
			{
				initamb(&amb0, i, rcv, j, val[i][j], std[i][j], num[i][j]);
				if ((stat = addamb(amb_nl, &amb0)) < 0)
					return;
			}
		}
	}
}

int main(int argc, char *argv[])
{
	int i, n;
	gtime_t ts, te;
	ambInf ambs[4000];
	prcopt_t prcopt = prcopt_default;
	filopt_t filopt = {""};
	char infile[1024], fcbfile[1024], nlfile[1024], *files[MAXEXFILE] = {0};
	unsigned int tick = tickget();

	/* load conf */
	if (!loadopts(argv[1], sysopts))
		return -1;
	getsysopts(&prcopt, NULL, &filopt);

	/* pro processing */
	ts = timeadd(prcopt.ts, 7200);
	te = prcopt.te;
	sprintf(infile, "%s%s", filopt.outdir, "*-23-03-23-WUM.ppp");
	sprintf(fcbfile, "%s%s", filopt.outdir, "FCB_230323-utl.txt");
	sprintf(nlfile, "%s%s", filopt.outdir, "FCB_DRAW_230323-utl.txt");

	for (i = 0; i < MAXEXFILE; i++)
	{
		if (!(files[i] = (char *)malloc(1024)))
		{
			for (i--; i >= 0; i--)
				free(files[i]);
			return 0;
		}
	}
	n = expath(infile, files, MAXEXFILE);

	/* generate el wl amb */
	printf("(1)read ambiguity\n");
	for (i = 0; i < n; i++)
	{
		read_arc(ts, te, files[i], &arc_el, 1);
		read_arc(ts, te, files[i], &arc_wl, 0);
		checkbrk("\t%d/%d", i, n);
	}
	// outarc(&arc_el,"arc_el.txt");
	// outarc(&arc_wl,"arc_wl.txt");
	arc2amb(&arc_el, &amb_el);
	arc2amb(&arc_wl, &amb_wl);
	estfcb(amb_el.data, amb_el.n, &fcb_el);
	estfcb(amb_wl.data, amb_wl.n, &fcb_wl);

	outarcpos(&arc_el, &fcb_el, "arc_el-pos.txt");
	outarcpos(&arc_wl, &fcb_wl, "arc_wl-pos.txt");

	if (!outhead(fcbfile, &fcb_el, &fcb_wl))
		return -1;

	/* generate nl amb */
	printf("(2)generate nl ambiguity\n");
	corr_fcb(&arc_wl, &fcb_wl);

	for (i = 0; i < n; i++)
	{
		read_amb(ts, te, files[i], &arc_wl, &amb_nl);
		checkbrk("\t%d/%d", i, n);
	}

	/* process nl fcb */
	printf("(3)process nl fcb\n");
	if (!sortamb(&amb_nl))
		return 0;
	// outamb(&amb_nl, &arc_wl, nlfile);
	if ((fp = openfile(fcbfile)))
	{
		while ((n = inputamb(&amb_nl, ambs)) > 0)
		{
			estfcb(ambs, n, &fcb_nl);
			outambpos(ambs, n, &fcb_nl, &arc_wl, "arc_nl-pos.txt");
			outbody(fp, &fcb_nl, ts, ambs[0].ind);
		}
		fclose(fp);
	}
	drawnlfcb(fcbfile, nlfile);

	for (i = 0; i < MAXEXFILE; i++)
		free(files[i]);
	freearc(&arc_el);
	freearc(&arc_wl);
	freeamb(&amb_el);
	freeamb(&amb_wl);
	freeamb(&amb_nl);

	checkbrk("Time=%.1f s\n", (tickget() - tick) * 0.001);

	getchar();
}