#include "rtklib.h"

#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max rinex record length */
#define MAXPOSHEAD	1024

/* read fcb header & body -----------------------------------------------------------*/
static int readfcbh(FILE *fp, nav_t *nav)
{
	double bias;
	char buff[MAXRNXLEN], *label = buff + 60;
	int i=0, k, sat;
	for (k = 0; k < MAXSAT; k++) {
		nav->wlbias[k] = 10000;
		nav->elbias[k] = 10000;
	}
	while (fgets(buff, MAXRNXLEN, fp)) {
		if (strlen(buff) <= 60) continue;
		else if (strstr(label, "COMMENT")) {
			if (!strncmp(buff, "WL", 2) && (sat = satid2no(buff + 4)) &&
				sscanf(buff + 14, "%lf", &bias) == 1) {
				nav->wlbias[sat - 1] = bias;
			}
			else if (!strncmp(buff, "EL", 2) && (sat = satid2no(buff + 4)) &&
				sscanf(buff + 14, "%lf", &bias) == 1) {
				nav->elbias[sat - 1] = bias;
			}
		}
		if (strstr(label, "END OF HEADER")) return 1;
		if (++i >= MAXPOSHEAD) break; /* no fcb file */
	}
	return 0;
}
static int readfcbb(FILE *fp, nav_t *nav) 
{
	gtime_t time;
	fcbd_t * nav_fcb;
	char buff[1024], satid[8] = "";;
	int num = 0, sat,i;
	while (fgets(buff, sizeof(buff), fp)) {
		if (!strncmp(buff, "*", 1)) {
			if (str2time(buff, 2, 21, &time)) {
				trace(2, "fcb invalid epoch: %34.34s\n", buff);
				continue;
			}
			nav->nf++;	
			if (nav->nf >= nav->nfmax) {
				nav->nfmax += 1024;
				if (!(nav_fcb = (fcbd_t *)realloc(nav->fcb, sizeof(fcbd_t)*(nav->nfmax)))) {
					trace(1, "readrnxclk malloc error: nmax=%d\n", nav->nfmax);
					free(nav->fcb); nav->fcb = NULL; nav->nf = nav->nfmax = 0;
					return -1;
				}
				nav->fcb = nav_fcb;
			}
			
			nav->fcb[nav->nf - 1].ti = timeadd(time, 450);
			nav->fcb[nav->nf - 1].iod = 450;
			for (i = 0; i<MAXSAT; i++) {
				nav->fcb[nav->nf - 1].bias[i] = 999.9;
				nav->fcb[nav->nf - 1].std[i] = 0.0f;
			}
		}
		
		if (!strncmp(buff, "P", 1)) {
			strncpy(satid, buff + 1, 4);
			sat = (unsigned char)satid2no(satid) - 1;
			nav->fcb[nav->nf - 1].bias[sat] = str2num(buff, 24, 6);
			nav->fcb[nav->nf - 1].std[sat] = (float)str2num(buff, 55, 6);		
		}
	}
	return 1;
}

extern int readfcb_sgg(const char *file, nav_t *nav) {
	FILE *fp;

	if (!(fp = fopen(file, "r"))) {
		trace(1, "rinex file open error: %s\n", file);
		return 0;
	}
	/* read rinex header */
	if (!readfcbh(fp, nav)) return 0;
	if (!readfcbb(fp, nav)) return 0;
	
	fclose(fp);
	return 1;
}

static int readifcbb(FILE *fp, nav_t *nav) 
{
	gtime_t time;
	ifcb_t * nav_ifcb;
	char buff[1024], satid[8] = "";;
	int num = 0, sat, i;
	while (fgets(buff, sizeof(buff), fp)) {
		if (!strncmp(buff, "*", 1)) {
			if (str2time(buff, 2, 21, &time)) {
				trace(2, "ifcb invalid epoch: %34.34s\n", buff);
				continue;
			}
			nav->ni++;	
			if (nav->ni >= nav->nimax) {
				nav->nimax += 1024;
				if (!(nav_ifcb = (ifcb_t *)realloc(nav->ifcb, sizeof(ifcb_t)*(nav->nimax)))) {
					trace(1, "readifcb malloc error: nmax=%d\n", nav->nimax);
					free(nav->ifcb); nav->ifcb = NULL; nav->ni = nav->nimax = 0;
					return -1;
				}
				nav->ifcb = nav_ifcb;
			}
			
			nav->ifcb[nav->ni - 1].ti = time;
			for (i = 0; i<MAXSAT; i++) {
				nav->ifcb[nav->ni - 1].bias[i] = 999.9;
			}
		}
		
		if (!strncmp(buff, "P", 1)) {
			strncpy(satid, buff + 1, 4);
			sat = (unsigned char)satid2no(satid) - 1;
			nav->ifcb[nav->ni - 1].bias[sat] = str2num(buff, 7, 19);
		}
	}
	return 1;
}

extern int readifcb_sgg(const char *file, nav_t *nav) {
	FILE *fp;

	if (!(fp = fopen(file, "r"))) {
		trace(1, "rinex file open error: %s\n", file);
		return 0;
	}
	/* read rinex header */
	if (!readifcbb(fp, nav)) return 0;
	
	fclose(fp);
	return 1;
}

static int uniqsnx(const char *file, stas_t *stas) {
	FILE *fp;
	sta_t sta;
	char *p, *q, name[1024], l_name[1024], buff[1024];
	int n = 0;

	if (!(fp = fopen(file, "r"))) {
		fprintf(stderr, "station list file read error %s\n", file);
		return 0;
	}
	while (fgets(buff, sizeof(buff), fp)) {
		if ((p = strchr(buff, '#'))) *p = '\0';
		for (p = strtok(buff, " \r\n"); p; p = strtok(NULL, " \r\n")) {
			strcpy(name, p);
			for (int i = n; i < stas->n; i++) {
				for (p = l_name, q = (char *)name; (*p = (char)toupper(*q)); p++, q++);		
				if (!strncmp(stas->data[i].name, l_name, 4)) {
					sta = stas->data[i];
					stas->data[i] = stas->data[n];
					stas->data[n] = sta;
					n++;	 break;
				}
			}
		}
	}
	stas->n = n;
	fclose(fp);
	return 1;
}

/* free stas data ------------------------------------------------------------*/
extern void freestas(stas_t *stas) 
{
	/* free station list */
	free(stas->data);   stas->data = NULL;  stas->n = stas->nmax = 0;
}

/* generate stationlist ------------------------------------------------------*/
extern int readstas(const char *file, stas_t *stas)
{
	FILE *fp;
	char buff[256];
	sta_t *stas_data;

	if (!(fp = fopen(file, "r"))) return 0;
	while (fgets(buff, sizeof(buff), fp)) {
		if (buff[0] == '%' || buff[0] == '#') continue;			
		if (stas->nmax <= stas->n) {
			stas->nmax += 100;
			if (!(stas_data = (sta_t *)realloc(stas->data, sizeof(sta_t)*stas->nmax))) {
				free(stas->data);	stas->data = NULL;		stas->n = stas->nmax = 0;
				return  0;
			}
			stas->data = stas_data;
		}
		if (sscanf(buff, "%s %lf %lf %lf\n", stas->data[stas->n].name, &stas->data[stas->n].pos[0], 
			&stas->data[stas->n].pos[1], &stas->data[stas->n].pos[2]) < 4) {
			if (sscanf(buff, "%s\n", stas->data[stas->n].name) < 1) {
				continue;
			}
		}
		stas->data[stas->n++].name[4] = '\0';
	}
	return stas->n;
}

static nav_t nav;
extern int drawnlfcb(const char *infile, const char *outfile) {
	FILE *fp;
	int week;
	double tow;

	if (!readfcb_sgg(infile, &nav)) return 0;
	if (!(fp = fopen(outfile, "w"))) {
		trace(1, "rinex file open error: %s\n", outfile);
		return 0;
	}
	for (int i = 0; i < nav.nf; i++) {
		tow=time2gpst(nav.fcb[i].ti,&week);
		fprintf(fp, "\n%d %.0f\t", week, tow);
		for (int j = 0; j < MAXSAT; j++) {
			fprintf(fp,"%10.4f", nav.fcb[i].bias[j]);
		}
	}
	fclose(fp);
	return 1;
}

/* read snxfile ---------------------------------------------------------------- */
extern void readsnx(const char *snxfile, const char *infile, const char *outfile)
{
	FILE *fp, *fp_out;
	char buff[256];
	int index;
	stas_t stas = { 0 };

	if (!(fp = fopen(snxfile, "r"))) return;
	while (fgets(buff, sizeof(buff), fp)) {
		if (!stas.nmax && !strncmp(buff, "%=SNX", 5)) {
			stas.nmax = str2num(buff, 60, 65);
			if (!(stas.data = (sta_t*)malloc(sizeof(sta_t)*stas.nmax))) {
				stas.nmax = stas.n = 0;
				return 0;
			}
		}
		else if (!strncmp(buff, "+SOLUTION/ESTIMATE", 18)) {
			while (fgets(buff, sizeof(buff), fp)) {
				if (!strncmp(buff, "-SOLUTION/ESTIMATE", 18)) break;
				index = (int)str2num(buff, 0, 6);
				if (!index) continue;
				if (index % 3 == 1) {
					strncpy(stas.data[stas.n].name, buff + 14, 4);
					stas.data[stas.n].name[4] = '\0';
					stas.data[stas.n].pos[0] = str2num(buff, 47, 68);
				}
				if (index % 3 == 2) {
					stas.data[stas.n].pos[1] = str2num(buff, 47, 68);
				}
				if (index % 3 == 0) {
					stas.data[stas.n].pos[2] = str2num(buff, 47, 68);
					stas.n++;
				}
			}
		}
	}
	if (!(fp_out = fopen(outfile, "w"))) {
		fclose(fp);		freestas(&stas);	  return;
	}
	if (infile) uniqsnx(infile, &stas);

	for (int i = 0; i < stas.n; i++) {
		if (stas.data[i].name[0] != '\0') {
			// double blh[3];
			// ecef2pos(stas.data[i].pos, blh);
			// fprintf(fp_out, "%s %15.6f %15.6f %15.6f\n", stas.data[i].name, 
			// blh[0] * 180 / PI, blh[1] * 180 / PI, blh[2]);
			fprintf(fp_out, "%s %15.6f %15.6f %15.6f\n", stas.data[i].name, stas.data[i].pos[0], 
					stas.data[i].pos[1], stas.data[i].pos[2]);
		}
	}

	fclose(fp);			fclose(fp_out);			freestas(&stas);
}

/* free ephemeris and pcv data -----------------------------------------*/
extern void freeproduct(nav_t *nav, pcvs_t *pcvs, pcvs_t *pcvr, stas_t *stas)
{
	/* free nav data */
	if (nav) freenav(nav,0x288);

	/* free antenna parameters */
	if (pcvs) { free(pcvs->pcv);	pcvs->pcv = NULL;	pcvs->n = pcvs->nmax = 0;}
	if (pcvr) { free(pcvr->pcv);	pcvr->pcv = NULL;	pcvr->n = pcvr->nmax = 0;}

	/* free stas data */
	if (stas) { freestas(stas);}
}

/* file to nav&pcv data ------------------------------------------------ */
extern int readproduct(const prcopt_t *prcopt, const filopt_t *fopt, nav_t *nav, 
	pcvs_t *pcvss, pcvs_t *pcvsr, stas_t *stas)
{
	int i, j;
	char *ext, path[1024];
	gtime_t ts = prcopt->ts, te = prcopt->te;

	/* read brdc data */
	if (*fopt->brdc) {
		free(nav->eph);	 nav->eph = NULL;  nav->n = nav->nmax = 0;
		free(nav->geph); nav->geph = NULL; nav->ng = nav->ngmax = 0;
		free(nav->seph); nav->seph = NULL; nav->ns = nav->nsmax = 0;
		reppath(fopt->brdc, path, ts, "", "");
		if (readrnxt(path, 0, ts, te, 0, prcopt->rnxopt[0], NULL, nav, NULL) <= 0) {
			checkbrk("error : insufficient memory");		return 0;
		};
		/* delete duplicated ephemeris */
		uniqnav(nav);
	}

	/* read procise clk */
	if (*fopt->clk && (ext = strrchr(fopt->clk, '.')) && (!strcmp(ext, ".clk") || !strcmp(ext, ".CLK"))) {
		free(nav->pclk); nav->pclk = NULL; nav->nc = nav->ncmax = 0;
		reppath(fopt->clk, path, ts, "", "");
		readrnxc(path, nav);
	}

	/* read procise orb */
	if (*fopt->sp3 && (ext = strrchr(fopt->sp3, '.')) && (!strcmp(ext, ".sp3") || !strcmp(ext, ".SP3"))) {
		free(nav->peph); nav->peph = NULL; nav->ne = nav->nemax = 0;
		reppath(fopt->sp3, path, ts, "", "");
		readsp3(path, nav, 0);
	}
	
	/* read erp data */
	if (*fopt->eop) {
		free(nav->erp.data); nav->erp.data = NULL; nav->erp.n = nav->erp.nmax = 0;
		reppath(fopt->eop, path, ts, "", "");
		if (!readerp(path, &nav->erp)) {
			checkbrk("no erp data\n");
		}
	}

	/* read dcb parameters */
	if (*fopt->dcb && (ext = strrchr(fopt->dcb, '.')) && ((!strcmp(ext, ".BIA")) || !strcmp(ext, ".DCB"))) {
		reppath(fopt->dcb, path, ts, "", "");
		readdcb(path, nav);
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
		reppath(fopt->fcb, path, ts, "", "");
		readfcb_sgg(path, nav);
	}
	/* read ifcb file */
	if (*fopt->ifcb && (ext = strrchr(fopt->ifcb, '.')) && (!strcmp(ext, ".ifcb"))) {
		free(nav->ifcb);	    nav->ni = nav->nimax = 0;
		reppath(fopt->ifcb, path, ts, "", "");
		readifcb_sgg(path, nav);
	}

	/* read satellite antenna parameters */
	if (*fopt->satantp && pcvss) {
		free(pcvss->pcv);	pcvss->n = pcvss->nmax = 0;			
		if (!(readpcv(fopt->satantp, pcvss))) {
			showmsg("error : no sat ant pcv in %s", fopt->satantp);
			return 0;
		}		
	}

	/* read receiver antenna parameters */
	if (*fopt->rcvantp && pcvsr) {
		free(pcvsr->pcv);	pcvsr->n = pcvsr->nmax = 0;
		if (!(readpcv(fopt->rcvantp, pcvsr))) {
			showmsg("error : no rec ant pcv in %s", fopt->rcvantp);
			return 0;
		}
	}

	/* read station position */
	if (*fopt->stapos && stas) {
		free(stas->data);stas->data = NULL; stas->n = stas->nmax = 0;
		readstas(fopt->stapos, stas);
	}

	return 1;
}

/* read obs data ----------------------------------------------------- */
extern int readobs(const char *infile, int rcv, const prcopt_t *popt, obs_t *obs, 
	sta_t *sta, int *nepoch)
{
	trace(3, "readobs:  n=%d\n");

	if(nepoch) *nepoch = 0;

	if (checkbrk("")) return 0;

	/* read rinex obs and nav file */
	if (readrnxt(infile, rcv, popt->ts, popt->te, popt->ti, popt->rnxopt[0], obs, NULL, sta)<0) {
		checkbrk("error : insufficient memory");
		trace(1, "insufficient memory\n");
		return 0;
	}
	if (obs->n <= 0) {
		checkbrk("error : no obs data");
		return 0;
	}
	/* sort observation data */
	if(nepoch) *nepoch = sortobs(obs);

	return 1;
}

extern void decode_corr(const char *file, nav_t *nav, int dt, int opt) 
{
	FILE *fp;
	gtime_t ttrop = { 0 }, tstec = { 0 },ti = { 0 };
	int init = 0, flag = 0, week, solq, nsat, lock1, lock2;
	double tow, xyz[3], azel[2], blh[3], zwd, zhd, stec, std;
	char *p, buff[256], name[8] = "", satid[8] = "", time[32] = "";
	corrtrop_t *ptrop = nav->corrtrop + nav->ntrop;
	corrstec_t *pstec = nav->corrstec + nav->nstec;

	ptrop->nt = ptrop->ntmax = pstec->ns = pstec->nsmax = pstec->nn = 0;
	if (!(fp = fopen(file, "r"))) {
		trace(2, "stec file open error %s\n", file);
		return;
	}
	while (fgets(buff, sizeof(buff), fp)) {
		if (buff[0]!='$') continue;
		for (p=buff;*p;p++) if (*p==',') *p=' ';
		if (!init) {
			if (!strncmp(buff+1,"NAME"   ,4)) {
				if (sscanf(buff,"$NAME%s\n",name)<1) continue;
				strcpy(ptrop->stas, name);
				strcpy(pstec->stas, name);
				flag += 1;
				trace(2,"%4d %4d %s\n", nav->ntrop, nav->nstec, pstec->stas);
			}
			if (!strncmp(buff+1,"POS"   ,3)) {
				if (sscanf(buff,"$POS%d%lf%d%d%lf%lf%lf\n", &week, &tow, &solq, &nsat,   
					&xyz[0], &xyz[1], &xyz[2]) < 7) continue;
				ecef2pos(xyz, blh);
				matcpy(pstec->rr, blh, 3, 1);	
				matcpy(ptrop->rr, blh, 3, 1);
				flag += 2;
				trace(2,"%4d %4d %.4f %.4f %.4f\n", nav->ntrop, nav->nstec, blh[0], blh[1], blh[2]);
			}
		}
		if (init) {
			if (!strncmp(buff+1,"TRP"   ,3)) {
				if (sscanf(buff,"$TRP%d%lf%d%lf%lf%lf\n", &week, &tow, &solq, &zwd, &std, &zhd)<6) continue;
				
				if (ptrop->nt >= ptrop->ntmax) {
					trop_t *pd;
					ptrop->ntmax += ptrop->ntmax <= 0 ? 300 : 50;
					if (!(pd = (trop_t *)realloc(ptrop->data, sizeof(trop_t)*(ptrop->ntmax)))) {
						trace(1, "readpppcorr trop malloc error: nmax=%d\n", ptrop->ntmax);
						free(ptrop->data);	ptrop->data = NULL; ptrop->nt = ptrop->ntmax = 0;
						break;
				  	}
					ptrop->data = pd;
				}
				ti = gpst2time(week, tow);
				if (timediff(ti, ttrop) >= dt) ttrop = ti;
				if (timediff(ti, ttrop) < DTTOL) {
					ptrop->data[ptrop->nt].time =  ti;
					ptrop->data[ptrop->nt].trp[0] = zwd;
					ptrop->data[ptrop->nt++].std[0] = std;
					time2str(ti,time,0);
					trace(2,"trop:%s %4d %.4f %.4f\n", time, ptrop->nt-1, zwd, std);
				}				
			}
			if (!strncmp(buff+1,"ION"   ,3)) {
				if (sscanf(buff,"$ION%d%lf%d%s%d%d%lf%lf%lf%lf\n", &week, &tow, &solq, satid, &lock1, &lock2, &azel[0], 
					&azel[1], &stec, &std)<10) continue;

				if (opt == 1 && solq != 1) continue;
				if (std > 0.07) continue;		// 基准站不行时，需要控制
				// if (std > 0.1) continue;		// by zh 23-05-24
				
				if (pstec->ns >= pstec->nsmax) {
					pstec->nsmax += pstec->nsmax <= 0 ? 90000 : 4000;
					stec_t *pd;
					if (!(pd = (stec_t *)realloc(pstec->data, sizeof(stec_t)*(pstec->nsmax)))) {
						trace(1, "readpppcorr stec malloc error: nmax=%d\n", pstec->nsmax);
						free(pstec->data);	pstec->data = NULL; pstec->ns = pstec->nsmax = 0;
						break;
				  	}
					pstec->data = pd;
				}
				ti = gpst2time(week, tow);

				if (timediff(ti, tstec) >= dt) {
					tstec = ti;
					pstec->time[pstec->nn] = ti;
					pstec->ii[pstec->nn++] = pstec->ns;
					if (pstec->nn >= MAXTIM) break;
				}
				if (timediff(ti, tstec) < DTTOL) {
					pstec->data[pstec->ns].sat = (unsigned char)satid2no(satid);
					pstec->data[pstec->ns].ion = stec;
					pstec->data[pstec->ns].el  = azel[1];
					pstec->data[pstec->ns].lock1  = lock1;
					pstec->data[pstec->ns++].std = std;	
					time2str(ti,time,0);
					trace(2,"stec:%s %4d %.4f %.4f\n", time, pstec->nn-1, stec, std);				
				}			
			}
		}	
		init = (flag == 3) ? 1 : 0;
	}
	fclose(fp);
	if (ptrop->nt >= 10) nav->ntrop++;
	if (pstec->ns >= 10) nav->nstec++;
}

extern int pppcorr_read(const char *infile, nav_t *nav) 
{
	int i, n = 0;
	char *ext, *files[MAXEXFILE] = { 0 };

	for (i = 0; i < MAXEXFILE; i++) {
		if (!(files[i] = (char *)malloc(1024))) {
			for (i--; i >= 0; i--) free(files[i]);
			return 0;
		}
	}
	freenav(nav, 0x200);

	/* read local products */
	if (infile) n = expath(infile, files, MAXEXFILE);
	for (i = 0; i < n; i++) {
		if ((ext = strrchr(files[i], '.')) && !strcmp(ext, ".ppp")) {	
			decode_corr(files[i], nav, 30, 1);			
		}
	}
	for (i = 0; i < MAXEXFILE; i++) free(files[i]);

	return 1;
}