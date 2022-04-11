/*------------------------------------------------------------------------------
* ppp_corr.c : ppp corrections functions
*
*          Copyright (C) 2015 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2015/05/20 1.0 new
*           2016/05/10 1.1 delete codes
*-----------------------------------------------------------------------------*/
#include "rtklib.h"
#define MAXSTA 50
#define MAXTIM 2880
#define DOTTLT 10

typedef struct {
	char name[8];		/* marker name */
	double pos[3];      /* station position (ecef) (m) */
} stat_t;
typedef struct {
	int n, nmax;        /* number of data/allocated */
	stat_t *data;
} corrstas_t;

typedef struct {        /* trop data type */
	gtime_t time;       /* time (GPST) */
	double trp[3];      /* zenith tropos delay/gradient (m) */
	double std[3];       /* std-dev (m) */
} trop_t;
typedef struct {        /* ppp corrections type */
	char stas[8];		/* station names */
	double rr[3];		/* station ecef positions (m) */
	int nt, ntmax;		/* number of trop data */
	trop_t *data;		/* trop data */
} corrtrop_t;

typedef struct {        /* stec data type */
	uint8_t sat;		/* satellite number */
	double	ion;        /* slant ionos delay (m) */
	float	std;        /* std-dev (m) */
} stec_t;
typedef struct {        /* ppp corrections type */
	char stas[8];		/* station names */
	double rr[3];		/* station ecef positions (m) */
	gtime_t time[MAXTIM];			/* time list (GPST) */
	int ns, nsmax, nn, ii[MAXTIM];	/* number of stec data */
	stec_t *data;		/* stec data */		
} corrstec_t;

static int ntrop = 0; /* number of stations */
static int nstec = 0; /* number of stations */
static corrstas_t stas = { 0 };
static corrtrop_t corrtrop[MAXSTA] = { 0 };
static corrstec_t corrstec[MAXSTA] = { 0 };

/* free ppp corrections --------------------------------------------------------*/
extern void pppcorr_free()
{
	for (int i = 0; i < MAXSTA; i++) {
		free(corrtrop[i].data);		corrtrop[i].data = NULL;	corrtrop[i].nt = corrtrop[i].ntmax = 0;
		free(corrstec[i].data);		corrstec[i].data = NULL;	corrstec[i].ns = corrstec[i].nsmax = corrstec[i].nn =0;
	}
	free(stas.data);	stas.data = NULL;	stas.n = stas.nmax = 0;
}
/* decode stas corrections -------------------------------------------------------- */
static int decode_stas(char *file, corrstas_t *corr)
{
	FILE *fp;
	char buff[256];
	stat_t *stas_data;

	if (!(fp = fopen(file, "r"))) return 0;
	while (fgets(buff, sizeof(buff), fp)) {
		if (corr->nmax <= corr->n) {
			corr->nmax += 100;
			if (!(stas_data = (stat_t *)realloc(corr->data, sizeof(stat_t)*corr->nmax))) {
				free(corr->data);	corr->data = NULL;		corr->n = corr->nmax = 0;
				return -1;
			}
			corr->data = stas_data;
		}
		if (sscanf(buff, "%s %lf %lf %lf\n", corr->data[corr->n].name, &corr->data[corr->n].pos[0], &corr->data[corr->n].pos[1], &corr->data[corr->n].pos[2]) < 4) {
			continue;
		}
		corr->data[corr->n++].name[4] = '\0';
	}
	return 1;
}
/* decode trop corrections -------------------------------------------------------- */
static int decode_trop(char *file, corrtrop_t *corr)			
{
	FILE *fp;
	gtime_t time;
	trop_t *ptrop;
	int i, n, ii, stat = 0;
	char buff[256];

	ii = corr->nt = 0;
	/* read trop files */
	if (!(fp = fopen(file, "r"))) {
		trace(2, "trop file open error %s\n", file);
		return stat;
	}
	/*while (fgets(buff, sizeof(buff), fp)) {
		if (!strncmp(buff, "+TROP/STA_COORDINATES", 21)) {
			while (fgets(buff, sizeof(buff), fp)) {
				if (buff[0] == '*') continue;
				if (!strncmp(buff, "-TROP/STA_COORDINATES", 21)) break;
				strcpy(corr->stas, buff, 4);
				corr->rr[0] = str2num(buff, 16, 13);
				corr->rr[1] = str2num(buff, 29, 13);
				corr->rr[2] = str2num(buff, 42, 13);
			}
		}
		if (!strncmp(buff, "+TROP/SOLUTION", 14)) {
			while (fgets(buff, sizeof(buff), fp)) {
				if (corr->nt >= corr->ntmax) {
					corr->ntmax += 300;
					if (!(ptrop = (trop_t *)realloc(corr->data, sizeof(trop_t)*(corr->ntmax)))) {
						trace(1, "readpppcorr trop malloc error: nmax=%d\n", corr->ntmax);
						free(corr->data); corr->data = NULL; corr->nt = corr->ntmax = 0;
						return -1;
					}
					corr->data = ptrop;
				}
				if (buff[0] == '*') continue;
				if (!strncmp(buff, "-TROP/SOLUTION", 14)) { stat = 1;  break; }
				str2time(buff, 6, 18, &time);

				corr->data[ii].time = time;
				corr->data[ii].trp[0] = str2num(buff, 19, 6)*0.001;
				corr->data[ii].std[0] = str2num(buff, 26, 6)*0.001;
				corr->data[ii].trp[1] = str2num(buff, 34, 6)*0.001;
				corr->data[ii].std[1] = str2num(buff, 41, 6)*0.001;
				corr->data[ii].trp[2] = str2num(buff, 49, 6)*0.001;
				corr->data[ii].std[2] = str2num(buff, 56, 6)*0.001;
				ii = corr->nt++;
			}
		}
	}*/

	char *p, tt[32], ymd[32], dms[32];
	double trop, std;
	while (fgets(buff, sizeof(buff), fp)) {
		if ((p = strstr(buff, "staname"))) {
			strncpy(corr->stas, buff + 14, 4);
			corr->stas[4] = '\0';
		}
		else if ((p = strstr(buff, "position"))) {
			corr->rr[0] = str2num(buff, 13, 14);
			corr->rr[1] = str2num(buff, 27, 14);
			corr->rr[2] = str2num(buff, 42, 14);
		}
		else if (buff[0] != '%') {
			if (sscanf(buff, "%s %s %lf %lf", ymd, dms, &trop, &std) < 4) {
				continue;
			}
			if (corr->nt >= corr->ntmax) {
				corr->ntmax += 300;
				if (!(ptrop = (trop_t *)realloc(corr->data, sizeof(trop_t)*(corr->ntmax)))) {
					trace(1, "readpppcorr trop malloc error: nmax=%d\n", corr->ntmax);
					free(corr->data); corr->data = NULL; corr->nt = corr->ntmax = 0;
					return -1;
				}
				corr->data = ptrop;
			}

			sprintf(tt, "%s %s", ymd, dms);
			str2time(tt, 0, 32, &time);
			corr->data[ii].time = time;
			corr->data[ii].trp[0] = trop;
			corr->data[ii].std[0] = std;
			ii = corr->nt++;
		}
	}
	if (corr->nt > 10) stat = 1;
	return stat;
}
/* decode stec corrections -------------------------------------------------------- */
static int decode_stec(char *file, corrstec_t *corr)
{
	FILE *fp;
	gtime_t ti, t0 = { 0 };
	stec_t *stec_data;
	int i, n, ii, sat, stat = 0;
	double azel[2], ion, std;
	char *p, buff[256], ymd[12] = "", dms[10] = "", satid[8] = "", time[32] = "";

	/* read stec files */
	if (!(fp = fopen(file, "r"))) {
		trace(2, "stec file open error %s\n", file);
		return -1;
	}

	ii = corr->ns = corr->nsmax = corr->nn = 0;
	while (fgets(buff, sizeof(buff), fp)) {
		if ((p = strstr(buff, "staname"))) {
			strncpy(corr->stas, buff + 14, 4);
			corr->stas[4] = '\0';
		}
		else if ((p = strstr(buff, "position"))) {
			corr->rr[0] = str2num(buff, 13, 14);
			corr->rr[1] = str2num(buff, 27, 14);
			corr->rr[2] = str2num(buff, 42, 14);
		}
		else if (buff[0] != '%') {
			if (sscanf(buff, "%s %s %s %lf %lf %lf %lf", &ymd, &dms, satid, &azel[0], &azel[1], &ion, &std) < 7) {			
				continue;
			}
			if (corr->ns >= corr->nsmax) {
				corr->nsmax = corr->nsmax <= 0 ? 70000 : corr->nsmax * 2;
				stec_data = (stec_t *)realloc(corr->data, sizeof(stec_t)*corr->nsmax);
				if (!stec_data) {
					free(corr->data);
					fclose(fp);
					return stat;
				}
				corr->data = stec_data;
			}
			sprintf(time, "%s %s", ymd, dms);
			str2time(time, 0, 32, &ti);
			sat = (unsigned char)satid2no(satid);

			if (timediff(ti, t0) > DOTTLT) {
				corr->time[corr->nn] = ti;
				corr->ii[corr->nn] = ii;
				if (++corr->nn >= MAXTIM) break;
				t0 = ti;
			}
			corr->data[ii].sat = sat;
			corr->data[ii].ion = ion;
			corr->data[ii].std = std;
			ii = corr->ns++;
		}
	}
	if (corr->ns > 10) stat = 1;
	return stat;
}
/* read ppp corrections -------------------------------------------------------- */
extern int file2corr(const filopt_t *fopt) 
{
	int i, n = 0;
	char *ext, *files[MAXEXFILE] = { 0 };

	for (i = 0; i < MAXEXFILE; i++) {
		if (!(files[i] = (char *)malloc(1024))) {
			for (i--; i >= 0; i--) free(files[i]);
			return 0;
		}
	}
	/* free ppp corrections */
	pppcorr_free();

	/* read snx position */
	readsnx(fopt->snx, fopt->staname, fopt->prosta);	

	/* read local products */
	if (fopt->prosta) n += expath(fopt->prosta, files + n, MAXEXFILE);
	if (fopt->trop)   n += expath(fopt->trop,   files + n, MAXEXFILE);
	if (fopt->stec)   n += expath(fopt->stec,   files + n, MAXEXFILE);

	for (i = 0; i < n; i++) {
		if ((ext = strrchr(files[i], '.')) && (!strcmp(ext, ".stas"))) {			/* read station position */
			if (decode_stas(files[i], &stas) > 0);
		}
		else if ((ext = strrchr(files[i], '.')) && (!strcmp(ext, ".trop"))) {		/* read trop file */
			if (decode_trop(files[i], corrtrop + ntrop) > 0) ntrop++;
		}
		else if ((ext = strrchr(files[i], '.')) && (!strcmp(ext, ".stec"))) {		/* read stec file */
			if (decode_stec(files[i], corrstec + nstec) > 0) nstec++;
		}
	}

	for (i = 0; i < MAXEXFILE; i++) free(files[i]);
	return 1;
}
/* get sta-position correction ------------------------------------------------- */
extern int pppcorr_stat(const char *staname, double *pos) 
{
	int stat = 0;
	for (int i = 0; i < stas.n; i++) {
		if (!strncmp(stas.data[i].name, staname, 4)) {
			pos[0] = stas.data[i].pos[0];
			pos[1] = stas.data[i].pos[1];
			pos[2] = stas.data[i].pos[2];
			stat = 1;   break;
		}
	}
	return stat;
}
/* get tropospheric correction ------------------------------------------------- */
extern int pppcorr_trop(gtime_t time, const double *pos, double *trp, double *std)			
{
	int i, stat;
	corrtrop_t corr = corrtrop[0];

	stat = 0;
	for (i = 0; i < corr.nt; i++) {
		if (fabs(timediff(time, corr.data[i].time)) < 5) {
			*trp = corr.data[i].trp[0];
			*std = corr.data[i].std[0];
			stat = 1;
			break;
		}
	}
    return stat;
}
/* get ionospherec correction -------------------------------------------------- */
extern int pppcorr_stec(gtime_t time, const double *pos, const int sat1, const int sat2, double *dstec, double *std)				
{
	int i, is, ie, stat;
	double ion[2], std_ion[2];
	corrstec_t corr = corrstec[0];

	is = ie = stat = 0;
	
	for (i = 0; i < corr.nn; i++) {			
		if (fabs(timediff(time, corr.time[i])) < 5) {
			is = corr.ii[i];
			ie = i == corr.nn - 1 ? corr.ns : corr.ii[i + 1];
			break;
		}
	}
	for (i = is; i < ie; i++) {
		if (corr.data[i].sat == sat1) { ion[0] = corr.data[i].ion; std_ion[0] = corr.data[i].std; stat++; }
		if (corr.data[i].sat == sat2) { ion[1] = corr.data[i].ion; std_ion[1] = corr.data[i].std; stat++; }
	}
	if (stat == 2) {
		*dstec = ion[0] - ion[1];
		*std = sqrt(std_ion[0] * std_ion[0] + std_ion[1] * std_ion[1]);
	}
	return stat > 1;
}
