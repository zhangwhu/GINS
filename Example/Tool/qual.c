
#include "rtklib.h"

#define MAXARC	2000


typedef struct {
	gtime_t ts, te;
	int sat, len;
	double mp[3];
} arcInf;

typedef struct {
	int n, nmax;
	arcInf *data;
} arc_t;

static int addarc(arc_t *arcs, arcInf *data)			
{
	arcInf *arc_data;

	if (arcs->nmax <= arcs->n) {
		if (arcs->nmax <= 0) arcs->nmax = 9000; else arcs->nmax *= 2;
		if (!(arc_data = (arcInf*)realloc(arcs->data,sizeof(arcInf)*arcs->nmax))) {
			trace(1, "addarcdata: memalloc error n=%dx%d\n", sizeof(arcInf), arcs->nmax);
			free(arcs->data); arcs->data = NULL; arcs->n = arcs->nmax = 0;
			return -1;
		}
		arcs->data = arc_data;
	}
	
	arcs->data[arcs->n++] = *data;

	return 1;
}

void decodeQual(char *infile) 
{	
	gtime_t ti, t_1 = { 0 }, t_2;
	FILE *fp_in, *fp_out;
	int i, n, k, ii, week, sat, solq, ini;
	int sflag[MAXSAT] = { 0 }, sflag_1[MAXSAT] = { 0 }, stat = 0;
	double tow, temp, azel, snr[3], mp[3], dGF[3], dMW[3], p12[3], tdcp[3], dmp[3], spppos, dpdopper;
	char outfile[1024], buff[256]="", *p, satid[8]="";
	arcInf arcd[MAXSAT] = { 0 };
	arc_t  arcs = {0};

	if ((p=strrchr(infile,'.'))) {
        strncpy(outfile, infile, p-infile);
        strcpy(outfile+(p-infile),".mp");
    }
	else return 0;

	if (!(fp_in = fopen(infile, "r")) || !(fp_out = fopen(outfile, "w"))) {
		trace(2, "file open error: file=%s\n", infile);
		return ;
	}

	// 提取弧段信息
	for (i = 0; i < MAXSAT; i++) arcd[i].len = 0; 
	while (fgets(buff, sizeof(buff), fp_in)) {
		if (buff[0]!='$') continue;
		if (!strncmp(buff+1,"QUL"   ,3)) {
			for (p=buff;*p;p++) if (*p==',') *p=' ';
			if (sscanf(buff,"$QUL%d%lf%d%s%lf%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf"
				"%lf%lf%lf%lf%lf%lf%lf%lf\n", 
				&week,   &tow,    &solq,    satid, &azel, &ini, 
				&snr[0], &snr[1], &snr[2], 
				&temp,   &temp,   &temp,	
				&temp,   &temp,   &temp,    &temp, &temp, &temp, 
				&temp,   &temp,   &temp,    &temp, &temp, &temp, 
				&temp,   &temp,   &temp,    &mp[0],&mp[1],&mp[2]) < 30) {
        		continue;
    		}
			ti = gpst2time(week,tow);
			sat = (unsigned char)satid2no(satid) - 1;

			if (timediff(ti, t_1) != 0.0) {
				for (i = 0; i < MAXSAT; i++) {
					if (sflag_1[i] && !sflag[i]) {						
						if (arcd[i].len > 1) {
							arcd[i].te = t_2;
							if ((stat = addarc(&arcs, &arcd[i])) < 0) return 0;
						}
						arcd[i].len = 0;		
					}
					sflag_1[i] = sflag[i];
					sflag[i] = 0;
				}
				t_2 = t_1;
				t_1 = ti;
			}
			sflag[sat] = 1;
			if (ini&0x3) {												
				if (arcd[sat].len > 1) {
					arcd[sat].te = t_2;
					if ((stat = addarc(&arcs, &arcd[sat])) < 0) return 0;
				}
				arcd[sat].len = 0;	
			}
			if (mp[0] && mp[1]) {
				if (arcd[sat].len == 0) {
					arcd[sat].sat = sat;
					arcd[sat].ts = ti;
					arcd[sat].mp[0] = mp[0];
					arcd[sat].mp[1] = mp[1];
					arcd[sat].mp[2] = mp[2];
					arcd[sat].len = 1;
				}
				else {
					n = ++arcd[sat].len;
					arcd[sat].mp[0] = arcd[sat].mp[0] * (n-1)/n + mp[0]/n;
					arcd[sat].mp[1] = arcd[sat].mp[1] * (n-1)/n + mp[1]/n;
					arcd[sat].mp[2] = arcd[sat].mp[2] * (n-1)/n + mp[2]/n;
				}
			}
		}
	}
	for (i = 0; i < MAXSAT; i++) {									
		if (sflag[i]) {
			if (arcd[i].len > 1) {
				arcd[i].te = t_1;
				if ((stat = addarc(&arcs, &arcd[i])) < 0) break;
			}
		}
	}

	// 获取MP信息并输出
	rewind(fp_in);
	while (fgets(buff, sizeof(buff), fp_in)) {
		if (buff[0]!='$') continue;
		if (!strncmp(buff+1,"QUL"   ,3)) {
			for (p=buff;*p;p++) if (*p==',') *p=' ';
			if (sscanf(buff,"$QUL%d%lf%d%s%lf%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf"
				"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n", 
				&week,   &tow,    &solq,   satid,   &azel,   &ini, 
				&snr[0], &snr[1], &snr[2], 
				&p12[0], &p12[1], &p12[2], 
				&temp,   &dGF[0], &temp,   &dGF[1], &temp,   &dGF[2],
				&temp,   &dMW[0], &temp,   &dMW[1], &temp,   &dMW[2],
				&tdcp[0],&tdcp[1],&tdcp[2],&mp[0],  &mp[1],  &mp[2], &spppos, &dpdopper) < 32) {
        		continue;
    		}
			ti = gpst2time(week, tow);
			sat = (unsigned char)satid2no(satid) - 1;

			dmp[0] = dmp[1] = dmp[2] = 0.0;
			for (k = 0; k < arcs.n; k++) {
				if (arcs.data[k].sat != sat) continue;
				if ((timediff(ti, arcs.data[k].ts) >= 0) && (timediff(ti, arcs.data[k].te) <= 0)) {
					dmp[0] = mp[0] ? mp[0] - arcs.data[k].mp[0] : 0.0;
					dmp[1] = mp[1] ? mp[1] - arcs.data[k].mp[1] : 0.0;
					dmp[2] = mp[2] ? mp[2] - arcs.data[k].mp[2] : 0.0;
					
					if (fabs(dmp[0]) > 7)  {
						trace(1,"%s %s %.f %10.4f\t", time_str(ti,0),satid,time2gpst(ti, NULL),mp[0]);
						trace(1,"%.f %.f %10.4f\n", time2gpst(arcs.data[k].ts,NULL),time2gpst(arcs.data[k].te, NULL),arcs.data[k].mp[0]);
					}
					break;
				}
			}
			// if (fabs(dMW[0])<1E-07)   dMW[0]   = 1E+03;
			// if (fabs(dGF[0])<1E-07)   dGF[0]   = 1E+03;
			// if (fabs(tdcp[0])<1E-07)  tdcp[0]  = 1E+03;
			// if (fabs(p12[0])<1E-07)   p12[0]   = 1E+03;
			// if (fabs(spppos)<1E-07)   spppos   = 1E+03;
			// if (fabs(dpdopper)<1E-07) dpdopper = 1E+03;

			fprintf(fp_out, "%s %s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f"
				" %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n", 
				time_str(ti,0),satid,snr[0],snr[1],snr[2],p12[0],p12[1],p12[2],dGF[0],dGF[1],dGF[2],
				dMW[0],dMW[1],dMW[2],tdcp[0],tdcp[1],tdcp[2],dmp[0],dmp[1],dmp[2], spppos, dpdopper);
		}
	}
	fclose(fp_in);	fclose(fp_out);
}
void decodePos(char *infile) 
{	
	gtime_t ti;
	FILE *fp_in, *fp_out;
	int i, week, solq, nsat;
	double tow, pdop, temp, tow_1 = 0.0, dt; 
	char outfile[1024], buff[256]="", *p;

	if ((p=strrchr(infile,'.'))) {
        strncpy(outfile, infile, p-infile);
        strcpy(outfile+(p-infile),".dop");
    }
	else return ;

	if (!(fp_in = fopen(infile, "r")) || !(fp_out = fopen(outfile, "w"))) {
		trace(2, "file open error: file=%s\n", infile);
		return ;
	}

	while (fgets(buff, sizeof(buff), fp_in)) {
		if (buff[0]!='$') continue;
		if (!strncmp(buff+1,"POS"   ,3)) { 
			for (p=buff;*p;p++) if (*p==',') *p=' ';
			if (sscanf(buff,"$POS%d%lf%d%d%lf%lf%lf%lf%lf%lf%d%lf\n",&week,&tow, 
				&solq,&nsat,&temp,&temp,&temp,&temp,&temp,&temp,&nsat,&pdop) < 12) {
        		continue;
    		}

			dt = tow - tow_1;		
		
			if (tow_1==0.0 || dt<1.5) {
				ti = gpst2time(week, tow);
				fprintf(fp_out, "%s %6d %10.4f\n",time_str(ti,0),nsat,pdop);
			}
			else {
				for (i = 0; i<round(dt); i++) {
					ti = gpst2time(week, tow_1+(i+1));
					fprintf(fp_out, "%s %6d %10.4f\n",time_str(ti,0),0,15.0);
				}
			}
			tow_1 = tow;
		}
	}
	fclose(fp_in);	fclose(fp_out);
}

void decodeRes(char *infile) 
{	
	gtime_t ti;
	FILE *fp_in, *fp_out;
	int week, solq;
	double tow, temp, resp[3], resc[3]; 
	char outfile[1024], buff[256]="", *p, satid[8]="";

	if ((p=strrchr(infile,'.'))) {
        strncpy(outfile, infile, p-infile);
        strcpy(outfile+(p-infile),".res");
    }
	else return 0;

	if (!(fp_in = fopen(infile, "r")) || !(fp_out = fopen(outfile, "w"))) {
		trace(2, "file open error: file=%s\n", infile);
		return ;
	}

	while (fgets(buff, sizeof(buff), fp_in)) {
		if (buff[0]!='$') continue;
		if (!strncmp(buff+1,"RES"   ,3)) { 
			for (p=buff;*p;p++) if (*p==',') *p=' ';
			if (sscanf(buff,"$RES%d%lf%d%s%lf%lf%lf%lf%lf%lf\n",&week,&tow,&solq, 
				satid,&resp[0],&resp[1],&resp[2],&resc[0],&resc[1],&resc[2]) < 10) {
        		continue;
    		}
			ti = gpst2time(week, tow);
			fprintf(fp_out, "%s %s %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",time_str(ti,0), 
				satid,resp[0],resp[1],resp[2],resc[0],resc[1],resc[2]);
		}
	}
	fclose(fp_in);	fclose(fp_out);
}
	

int main(int argc, char *argv[]) 
{
    int i, n;
    char infile[MAXEXFILE], *files[MAXEXFILE] = { 0 };
    prcopt_t prcopt = prcopt_default;
	solopt_t solopt = solopt_default;
	filopt_t filopt = {""};

    /* load conf */
	// if (!loadopts(argv[1], sysopts)) return -1;
	// getsysopts(&prcopt, &solopt, &filopt);
    // sprintf(infile, "%s*.ppp", filopt.outdir);

	sprintf(infile, "/media/zhuhang/D/Data/PPPRTK_INS-Wuh-220519/Result/Sept/*f2-urban-pp.ppp");
	// sprintf(infile, "/media/zhuhang/D/Data/PPPRTK_INS-Wuh-220519/Result/Sept/*05-24.ppp");
	// sprintf(infile, "/media/zhuhang/D/Data/PPPRTK_INS-Wuh-220519/Result/*-19.ppp");
	// sprintf(infile, "/media/zhuhang/D/Data/PPPRTK_INS_VISUAL_Wuh-230416/Result/PP7/*.ppp");
    /* expand wild-card */
    for (i = 0; i<MAXEXFILE; i++) {
		if (!(files[i] = (char *)malloc(1024))) {
			for (i--; i >= 0; i--) free(files[i]);
			return -1;
		}
	}
	if ((n = expath(infile, files, MAXEXFILE)) <= 0) {
		for (i = 0; i<MAXEXFILE; i++) free(files[i]);
		return 0;
	}
	
    for (i = 0; i<n; i++) {
        decodeQual(files[i]);
		decodePos(files[i]);
		decodeRes(files[i]);
    }

    for (i = 0; i<MAXEXFILE; i++) free(files[i]);
    return 1;
}