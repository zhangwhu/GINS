#include "rtklib.h"

#define MAXRNXLEN   (16*MAXOBSTYPE+4)   /* max rinex record length */
#define MAXPOSHEAD	1024

/* compare satellite fcb -----------------------------------------------------*/
static int cmpfcb(const void *p1, const void *p2)
{
	fcbd_t *q1 = (fcbd_t *)p1, *q2 = (fcbd_t *)p2;
	double tt = timediff(q1->ts, q2->ts);
	return tt<-1E-3 ? -1 : (tt>1E-3 ? 1 : 0);
}

/* read rinex file -----------------------------------------------------------*/
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
			
			nav->fcb[nav->nf - 1].ts = time;
			nav->fcb[nav->nf - 1].te = timeadd(time, 900);
			for (i = 0; i<MAXSAT; i++) {
				nav->fcb[nav->nf - 1].bias[i][0] = 999.9;
				nav->fcb[nav->nf - 1].std[i][0] = 0.0f;
			}
		}
		
		if (!strncmp(buff, "P", 1)) {
			strncpy(satid, buff + 1, 4);
			sat = (unsigned char)satid2no(satid) - 1;
			nav->fcb[nav->nf - 1].bias[sat][0] = str2num(buff, 24, 6);
			nav->fcb[nav->nf - 1].std[sat][0] = (float)str2num(buff, 55, 6);		
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
static nav_t nav;
extern int drawnlfcb(const char *infile, const char *outfile) {
	FILE *fp;
	if (!readfcb_sgg(infile, &nav)) return 0;
	if (!(fp = fopen(outfile, "w"))) {
		trace(1, "rinex file open error: %s\n", outfile);
		return 0;
	}
	for (int i = 0; i < nav.nf; i++) {
		fprintf(fp, "\n%s\t", time_str(nav.fcb[i].ts, 0));
		for (int j = 0; j < MAXSAT; j++) {
			fprintf(fp,"%10.4f", nav.fcb[i].bias[j][0]);
		}
	}
	fclose(fp);
	return 1;
}