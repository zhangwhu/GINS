#include "rtklib.h"

#define MAXIND       3000
#define MAXRCV       300

typedef struct {
    uint8_t num[MAXIND][MAXSAT];
	float   azl[MAXIND][MAXSAT];
	double  dif[MAXIND][MAXSAT];
} dif_t;

static dif_t difs[MAXRCV] = { 0 };

static int time2ind(const gtime_t ts, gtime_t t, double nsec)
{
	double n;
	if (timediff(t, ts) >= 0.0f)
	{
		n = (int)(timediff(t, ts) / (nsec));
		return n;
	}
	return -1;
}

static gtime_t ind2time(const gtime_t ts, int ind, double nsec)
{
	return timeadd(ts, ind * nsec);
}

void read_dif(gtime_t ts, gtime_t te, const char *file, dif_t *difs) 
{
    FILE *fp;
    gtime_t ti;
    char *p, buff[256] = "", satid[8] = "";
    int week, solq, slip, sat, ind;
    double ddif, tow, azel;

    if (!(fp = fopen(file, "r")) || !(p = strrchr(file, FILEPATHSEP)))
	{
		trace(2, "input file open error: file=%s\n", file);
		return;
	}

	while (fgets(buff, sizeof(buff), fp)) {
        if (buff[0] != '$')
			continue;
		if (!strncmp(buff + 1, "DIF", 3))
		{
			for (p = buff; *p; p++)
				if (*p == ',')
					*p = ' ';
			if (sscanf(buff, "$DIF%d%lf%d%s%lf%d%lf\n", &week,
					   &tow, &solq, satid, &azel, &slip, &ddif) < 7)
			{
				continue;
			}
			ti = gpst2time(week, tow);
			sat = (unsigned char)satid2no(satid) - 1;

			if (timediff(ts, ti) > 0 || timediff(ti, te) > 0)
				continue;
            ind = time2ind(ts, ti, 30);
			if (difs->num[ind][sat] > 1) {
				double mm = difs->dif[ind][sat] / difs->num[ind][sat];
				// if (fabs(ddif - mm) > 0.02) {
				// 	printf("Skip: %s %4d %.4f\n", satid, ind, fabs(ddif - mm));
				// 	continue;
				// }
			}
			// printf("%4d %4d\n", ind, sat);
            difs->dif[ind][sat] = ddif;
			difs->azl[ind][sat] = azel;
            difs->num[ind][sat]++;
        }
    }
}

static void outdifpos(gtime_t ts, dif_t *dif, char *file)
{
	FILE *fp;
    gtime_t ti; 
	int ii, is,flag[MAXSAT] = {0};
    double IFCB, temp[MAXSAT], tow;
	char id[32], str[32];
	if (!(fp = fopen(file, "w")))
		return;
	for (ii = 0; ii < MAXIND; ii++)
	{   
        ti  = ind2time(ts, ii, 30);
		tow = time2gpst(ti, NULL);
        time2str(ti, str, 0);
		fprintf(fp,"* %s\n", str);
        for (is = 0; is < MAXSAT; is++) {
            // if ((k = sat2sys(is + 1)) < 0)
			//     continue;
		    // if (k != 3)
			//     continue;
            if (dif->num[ii][is] < 1) continue;
            satno2id(is + 1, id);
            IFCB = dif->dif[ii][is] / dif->num[ii][is];
			if (!flag[is]) {
				temp[is] = -IFCB;
				flag[is] = 1;
				continue;
			} else {
				if (fabs(IFCB) > 0.01) continue;
				temp[is] += IFCB;
			}
			fprintf(fp, "P%s \t \t %10.6f\n", id, temp[is]);
        }
	}
	fclose(fp);
}

int main(int argc, char *argv[]) {
    gtime_t ts, te;
    prcopt_t prcopt = prcopt_default;		  
	filopt_t filopt = { "" };
    dif_t data = { 0 };
    int i, n, is, ii, ir;
    char infile[1024], outfile[1024], *files[MAXEXFILE] = { 0 };

    /* load conf */
	if (!loadopts(argv[1], sysopts)) return -1;			
	getsysopts(&prcopt, NULL, &filopt);

    ts = prcopt.ts;
	te = prcopt.te;
	
    sprintf( infile, "%s%s", filopt.outdir,"*-23-03-23-WUM.ppp");
    sprintf(outfile, "%s%s", filopt.outdir, "IFCB_230323-utl.ifcb");

    for (i = 0; i < MAXEXFILE; i++) {
		if (!(files[i] = (char *)malloc(1024))) {
			for (i--; i >= 0; i--) free(files[i]);
			return 0;
		}
	}
	n = expath(infile, files, MAXEXFILE);

    /* generate el wl amb */
	printf("(1)read dif\n");
	for (i = 0; i < n; i++) {
		read_dif(ts, te, files[i], difs + i);
		checkbrk("\t%d/%d", i, n);
	}

	for (ii = 0; ii < MAXIND; ii++) {
		for (is = 0; is < MAXSAT; is++) {
			int num = 0, dd = 0;
			double temp[MAXRCV] = { 0 };
			double result = 0;
			for (ir = 0; ir < n; ir++) {
				if (difs[ir].num[ii][is] < 1) continue; 
				temp[num++] = difs[ir].dif[ii][is];
			}
			/* 中位数 */
			double mid = median(temp, num);
			/* 踢掉数据 */
			for (int c = 0; c < num; c++) {
				if (fabs(temp[c] - mid) < 0.01) {
					result += temp[c];
					dd++;
				}
			}
			/* 求平均值 */
			if (dd > 0.1) {
				data.dif[ii][is] = result;
				data.num[ii][is] = dd;
			}
		}
	}

    outdifpos(ts, &data, outfile);

}