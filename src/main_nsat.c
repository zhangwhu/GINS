#include <stdarg.h>
#include "rtklib.h"

#define NT 24
static nav_t nav = { 0 };

static double outazel(double *rs, double *blh) 
{
	double rr[3], e[3], r;
	pos2ecef(blh, rr);
	if (norm(rs, 3)<RE_WGS84) return -1.0;
	for (int i = 0; i < 3; i++) e[i] = rs[i] - rr[i];
	r = norm(e, 3);             
	for (int i = 0; i<3; i++) e[i] /= r;
	return satazel(blh, e, NULL);
}
static int readbrdc(char *file, nav_t *nav) 
{
	gtime_t ts = { 0 }, te = { 0 };
	char opt[2] = "";

	free(nav->eph);	 nav->eph = NULL;  nav->n = nav->nmax = 0;
	free(nav->geph); nav->geph = NULL; nav->ng = nav->ngmax = 0;
	free(nav->seph); nav->seph = NULL; nav->ns = nav->nsmax = 0;
	if (readrnxt(file, 0, ts, te, 0, opt, NULL, nav, NULL) < 0) {
		checkbrk("error : insufficient memory");		return 0;
	};
	/* delete duplicated ephemeris */
	uniqnav(nav);
	return 1;
}
static int process(gtime_t ti, nav_t *nav, double *lat, double *lon, int *nsat)			
{
	FILE *fp;
	int flag[MAXSAT] = { 0 }, mm, nn, sys, prn;
	double blh[3];
	double *rs, *dts, *var, *svh;

	mm = (lat[1] - lat[0]) / lat[2] + 1;
	nn = (lon[1] - lon[0]) / lon[2] + 1;
	rs = mat(MAXSAT, 6);		dts = mat(MAXSAT, 2);		var = mat(MAXSAT, 1);		svh = mat(MAXSAT, 1);

	for (int i = 0; i < MAXSAT; i++) {
		if (satpos(ti, ti, i + 1, EPHOPT_BRDC, nav, rs + i * 6, dts + i * 2, var + i, svh + i)) {
			flag[i] = 1;
		}
	}

	blh[2] = 1000.0;
	for (int m = 0; m < mm; m++) {	
		blh[0] = (lat[0] + lat[2] * m)*D2R;
		for (int n = 0; n < nn; n++) {	
			blh[1] = (lon[0] + lon[2] * n)*D2R;
			for (int i = 0; i < MAXSAT; i++) {
				sys = satsys(i + 1, &prn);
				if (!flag[i]) continue;
				//if ((sys == SYS_CMP && prn > 18) || sys == SYS_GAL) {
				if (sys == SYS_CMP && prn > 18 /*&& prn<59*/) {
				//if (sys == SYS_CMP && prn < 19) {
				//if (sys == SYS_GAL) {
				//if (sys == SYS_GLO) {
					if (outazel(rs + i * 6, blh) > 10 * D2R) {
						nsat[n + m*nn] += 1;
					}
				}
			}
		}
	}

	free(rs);	free(dts);	 free(var);	  free(svh);
	return 1;
}
int main2() 
{
	FILE *fp;
	gtime_t ti;
	int *nsat, mm, nn;
	char time[32] = "2021 08 10 00 00 0.00";
	char brdc[1024] = "E:\\GNSS_DATA\\products\\2021-08-10\\brdc.21p";
	/*char time[32] = "2017 07 22 00 00 0.00";
	char brdc[1024] = "E:\\GNSS_DATA\\brdc203.17p";*/
	char outfile[1024] = "E:\\GNSS_DATA\\result\\2021-08-10\\nsat.txt";
	double lat[3] = { 89, -89, -1.0 }, lon[3] = { -180.0, 180.0, 2.0 };

	str2time(time, 0, 32, &ti);

	if (!readbrdc(brdc, &nav)) return 0;

	mm = (lat[1] - lat[0]) / lat[2] + 1;
	nn = (lon[1] - lon[0]) / lon[2] + 1;
	nsat = imat(mm, nn);
	for (int m = 0; m < mm; m++) for (int n = 0; n < nn; n++) nsat[n + m*nn] = 0;
	for (int i = 0; i < NT; i++) {		
		ti = timeadd(ti, 3600);
		if (!process(ti, &nav, lat, lon, nsat)) continue;
	}

	if (!(fp = fopen(outfile, "w"))) {
		trace(1, "error open file %s\n", outfile);
		return 0;
	}

	for (int m = 0; m < mm; m++) {
		for (int n = 0; n < nn; n++) {
			fprintf(fp, "%10f", (double)nsat[n + m*nn] / NT);
		}
		fprintf(fp,"\n");
	}

	free(nsat);
	fclose(fp);
}