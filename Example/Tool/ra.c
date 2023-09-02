
#include "rtklib.h"

#define SQR(x)      ((x)*(x))

static nav_t nav = { 0 };
static nav_t nav0 = { 0 };

int main(int argc, char *argv[]) 
{
    prcopt_t prcopt = prcopt_default;
	solopt_t solopt = solopt_default;
	filopt_t filopt = {""};
	FILE *fp_stec, *fp_trop;
	char id1[8], id2[8], time[32], path[1024];
	int i, j, k, m, n, sat1, sat2, iref[4], isat[4][20], sn[4], sys, week;
	double tow, trp, ion, var, value, maxel[4], el, dd, qq;
	
    /* load conf */
	if (!loadopts(argv[1], sysopts)) return -1;
	getsysopts(&prcopt, &solopt, &filopt);

	reppath(filopt.corr, path, prcopt.ts, "", "");

	pppcorr_read(path, &nav);
	decode_corr("/media/zhuhang/E/additional/Results/RA/WUH2-23-03-23-GEC-new.ppp", &nav0, 30, 1);		
	
	char out1[1024], out2[1024], out3[1024];
	sprintf(out1, "%s%s",filopt.outdir, "WUH2-ion-0323-1-new.txt");
	fp_stec = fopen(out1,"w");
	sprintf(out2, "%s%s",filopt.outdir, "WUH2-trp-0323-1-new.txt");
	fp_trop = fopen(out2,"w");
	sprintf(out3, "%s%s",filopt.outdir, "WUH2-trace-exc-new.txt");
	traceopen(out3);	tracelevel(1);

	corrstec_t *stec = nav0.corrstec;
	corrtrop_t *trop = nav0.corrtrop;
	gtime_t tss, tee;
	char time1[32] = "2023 03 23 08 00 00.00";		// "2023 04 17 03 00 00.00"
	char time2[32] = "2023 03 24 00 00 00.00";		// "2023 04 17 05 00 00.00"
	str2time(time1, 0, 32, &tss);
	str2time(time2, 0, 32, &tee);
	
	for (i = 0; i<stec->nn-1; i++) {
		for (k = 0; k<4; k++) { maxel[k] = sn[k] = 0; }
		for (k = stec->ii[i]; k<stec->ii[i+1]; k++) {	
			if (!(sys = satsys(stec->data[k].sat, NULL))) continue; 
			if (!(pppcorr_stat(stec->time[i], stec->rr, stec->data[k].sat, &nav))) continue; 
			el = stec->data[k].el;
			switch (sys) {
			case SYS_GPS: isat[0][sn[0]++] = k; if (el > maxel[0]) { iref[0] = k; maxel[0] = el; } break;
			case SYS_GAL: isat[1][sn[1]++] = k; if (el > maxel[1]) { iref[1] = k; maxel[1] = el; } break;
			case SYS_CMP: isat[2][sn[2]++] = k; if (el > maxel[2]) { iref[2] = k; maxel[2] = el; } break;
			}
		} 
		
		for (k = 0; k < 3; k++) {
			n = iref[k];
			for (j = 0; j < sn[k]; j++) {
				m = isat[k][j]; 
				if (m == n) continue;
				sat1 = stec->data[m].sat;
				sat2 = stec->data[n].sat;
				if (sat1 == sat2) continue;

				dd = (stec->data[m].ion - stec->data[n].ion);
				qq = (SQR(stec->data[m].std) + SQR(stec->data[n].std));
				if (timediff(stec->time[i], tss) < 0 || timediff(stec->time[i], tee) > 0) continue;

				if (pppcorr_stec(stec->time[i], stec->rr, sat1, sat2, &ion, &var, &nav, 2)) {	

					satno2id(sat1, id1);
					satno2id(sat2, id2);
					time2str(stec->time[i], time, 0);

					value =  dd - ion;

					if (fabs(value)<=0.2) {
					// if (1/* fabs(value)>0.1 */) {
					tow = time2gpst(stec->time[i], &week);
					fprintf(fp_stec, "%d %.1f %d %8s %10.3f %8s %10.3f %10.4f %10.4f %10.4f %10.4f %10.4f\n", 
						week, tow, k, id1, stec->data[m].el, id2, stec->data[n].el, dd, sqrt(qq), ion, sqrt(var), value);
					pppcorr_stec(stec->time[i], stec->rr, sat1, sat2, &ion, &var, &nav, 1);
					}
				}	
			}
		}
	}

	for (i = 0; i<trop->nt-1; i++) {
		if (timediff(trop->data[i].time, tss) < 0 || timediff(trop->data[i].time, tee) > 0) continue;
		if (pppcorr_trop(trop->data[i].time, trop->rr, &trp, &var, &nav)) {		

			value = trop->data[i].trp[0] - trp;
			tow = time2gpst(trop->data[i].time, &week);
			fprintf(fp_trop, "%d %.1f %7.4f %7.4f %7.4f %7.4f %7.4f\n", week, tow, trop->data[i].trp[0], 
					trop->data[i].std[0], trp, sqrt(var), value);
		}
	}
	/* fclose(fp_stec); */		/* fclose(fp_trop); */ 		/* traceclose(); */
}
