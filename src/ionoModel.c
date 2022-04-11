#include "rtklib.h"

#define SQR(x)      ((x)*(x))

/* tec data record */
typedef struct {
	gtime_t time;
	unsigned char sat, rcv, sys;
	double ipp[3];						// 穿刺点经纬度、天顶距
	double tec, std;
} itecd_t;

typedef struct {			
	int n, nmax;
	itecd_t *data;
} itec_t;

static itec_t tec = { 0 };	
static stas_t stas = { 0 };

static const int    optsys = 1;						// "1:gps+2:sbas+4:glo+8:gal+16:qzs+32:comp"
static const int    order  = 2;
static const int    nhour  = 2;
static const double Re     = 6371.0;	//6371.0  6371.137
static const double hion   = 450;	    //450 506.7

static const double BM = 80.33*D2R;				//  79.91				//2008年地磁参数（刘长建）
static const double LM = -72.67*D2R;			    // -72.05

static double BL0[2] = { 0.0 };

static char start_time[32] = "2021 08 11 02 00 0.0";			//4h
static char   end_time[32] = "2021 08 12 00 00 0.0";

static int listsys[4] = { SYS_GPS, SYS_GLO, SYS_GAL, SYS_CMP };

static double *x = NULL;
static int ns = 0, nr = 0;
static int sat2ind[MAXSAT] = { 0 }, ind2sat[MAXSAT] = { 0 };
static int rcv2sta[MAXRCV] = { 0 }, rcv2ind[MAXRCV][4] = { 0 }, 
		   ind2rcv[MAXRCV*4] = { 0 }, ind2sys[MAXRCV*4] = { 0 };
static int nexp = -1;

/* free amb data ------------------------------------------------------------*/
static void freetec(itec_t *stec) 
{
	free(stec->data);	stec->data = NULL;	stec->n = stec->nmax = 0;
}
static void freestas(stas_t *stas) 
{
	free(stas->data);   stas->data = NULL;  stas->n = stas->nmax = 0;
}
static void geo2mag(const double *geo, double *mag)				
{
	double sinL, cosL;
	mag[0] = asin(sin(BM)*sin(geo[0])+cos(BM)*cos(geo[0])*cos(geo[1]-LM));
	sinL =  cos(geo[0])*sin(geo[1] - LM) / cos(mag[0]);
	cosL = -(sin(geo[0]) - sin(BM)*sin(mag[0])) / (cos(BM)*cos(mag[0]));
	mag[1] = atan2(sinL,cosL);
}
static void geo2solarmag(const double *ipp_geo, double tt, double *solar_mag) 
{
	double ipp_mag[2], sun_geo[2], sun_mag[2];
	geo2mag(ipp_geo,ipp_mag);
	sun_geo[0] = 0;				sun_geo[1] = PI - tt*PI / 43200;
	geo2mag(sun_geo, sun_mag);
	solar_mag[0] = ipp_mag[0];	solar_mag[1] = ipp_mag[1] - sun_mag[1];

	if (solar_mag[1] < -PI) solar_mag[1] += 2 * PI;
	else if (solar_mag[1] > PI) solar_mag[1] -= 2 * PI;
}
static int time2ind(const gtime_t ts, gtime_t t, int nhour, double *dt)		
{
	int i;
	double tt;
	if (timediff(t, ts) >= 0.0f) {
		tt = timediff(t, ts);
		i = (int)(tt / (nhour * 3600));
		if (dt) *dt = (i + 1) - tt / (nhour * 3600);
		return i;
	}
	return -1;
}
static void P_nm(double x,double *P, int N) 
{
	int n, m;
	double a = 0, b = 0, c = 0;
	P[0 + 0] = 1;		P[0 + 1 * N] = sqrt(3)*cos(x);			P[1 + 1 * N] = -sqrt(3)*sin(x);				//P[1 + 1 * N] = sqrt(3)*sin(x);
	for (n = 2; n < N; n++) {
		c = sqrt((2.0*n + 1) / (2.0*n));
		P[n + n*N] = -c*sin(x)*P[n - 1 + (n - 1)*N];			//P[n + n*N] = c*sin(x)*P[n - 1 + (n - 1)*N]							
		for (m = 0; m < n; m++) {
			a = sqrt(((2.0 * n - 1)*(2.0 * n + 1)) / ((n - m)*(n + m)));
			b = sqrt(((2.0 * n + 1)*(n + m - 1)*(n - m - 1)) / ((2.0 * n - 3)*(n + m)*(n - m)));	
			P[m + n*N] = a*cos(x)*P[m + (n - 1)*N] - b*P[m + (n - 2)*N];
		}
	}
}
/* 分段线性 */
static void Harmony(int order, int ndcb, const double *solar_mag, double *b, double dt, int iperiod) 
{
	int i = 0, n, m, ind, N;
	double *P, s, p_sin = 0, p_cos = 0;
	
	N = order + 1;
	P = zeros(N, N);
	P_nm(PI / 2 - solar_mag[0], P, N);		
	s = solar_mag[1];						
	ind = ndcb + N*N*iperiod;
	for (n = 0; n <= order; n++) {
		for (m = 0; m <= n; m++) {
			p_sin = P[m + n*N] * sin(m*s);
			p_cos = P[m + n*N] * cos(m*s);
			if (m == 0) {
				b[ind + i] = P[m + n*N] * dt;
				b[ind + i + N*N] = P[m + n*N] * (1 - dt);
			}
			else {
				b[ind + i] = p_cos*dt;
				b[ind + i + N*N] = p_cos*(1 - dt);
				i++;
				b[ind + i] = p_sin*dt;
				b[ind + i + N*N] = p_sin*(1 - dt);
			}
			i++;
		}
	}
	free(P);
}
static double result_Harmony(int iperiod, int order, int ndcb, double *solar_mag, double dt, double *x)			
{
	int m, n, ind, N, kk = 0;
	double *P, s, a, b,vtec = 0;

	N = order + 1;
	P = zeros(N, N);
	P_nm(PI / 2 - solar_mag[0], P, N);
	s = solar_mag[1];
	ind = ndcb + N*N*iperiod;
	for (n = 0; n <= order; n++) {
		for (m = 0; m <= n; m++) {
			if (m == 0) {
				a = x[ind + kk] * dt + x[ind + N*N + kk] * (1 - dt);
				vtec += P[m + n*N] * a;
			}
			else {
				a = x[ind + kk] * dt + x[ind + N*N + kk] * (1 - dt);
				kk++;
				b = x[ind + kk] * dt + x[ind + N*N + kk] * (1 - dt);
				vtec += P[m + n*N] * (a*cos(m*s) + b*sin(m*s));
			}
			kk++;
		}
	}
	free(P);
	return vtec;
}
/* data index (i:lat,j:lon,k:hgt) --------------------------------------------*/
static int dataindex(int i, int j, int k, const int *ndata)
{
	if (i<0 || ndata[0] <= i || j<0 || ndata[1] <= j || k<0 || ndata[2] <= k) return -1;
	return i + ndata[0] * (j + ndata[1] * k);
}
static nav_t nav = { 0 };

static int file2tec(const char *infile, const sta_t *sta, int ircv, unsigned char *sflag, unsigned char *syslist)
{
	FILE *fp;
	gtime_t ts, te, ti;
	itecd_t *tec_data;
	char buff[256], ymd[12] = "", dms[10] = "", satid[8] = "", time[32] = "";
	double pos[3], azel[2], ipp[3], stec1, stec2, temp, std;
	unsigned char sys, sat, flag;

	str2time(start_time, 0, 32, &ts);
	str2time(  end_time, 0, 32, &te);
	*syslist = 0;

	ecef2pos(sta->pos, pos);
	if (!(fp = fopen(infile, "r"))) return 0;
	while (fgets(buff, sizeof(buff), fp)) {
		if (sscanf(buff, "%s %s %s %lf %lf %lf %lf %lf", &ymd, &dms, satid, &azel[0], &azel[1], &stec1, &stec2, &std) < 8) {
			continue;
		}
		if (tec.n >= tec.nmax) {
			tec.nmax = tec.nmax <= 0 ? 131072 : tec.nmax * 2;
			tec_data = (itecd_t *)realloc(tec.data, sizeof(itecd_t)*tec.nmax);
			if (!tec_data) {
				freetec(&tec);
				fclose(fp);
				return 0;
			}
			tec.data = tec_data;
		}
		sprintf(time, "%s %s", ymd, dms);
		str2time(time, 0, 32, &ti);
		sat = (unsigned char)satid2no(satid);
		sys = (unsigned char)satsys(sat, NULL);
		azel[0] *= D2R;		 azel[1] *= D2R;

		if (timediff(ts, ti) > 0.0f || timediff(ti, te) > 0.0f || azel[1] < 15 * D2R || !(sys & optsys)) continue;
		
		sflag[sat - 1] = 1;
		*syslist |= sys;

		tec.data[tec.n].time = ti;
		tec.data[tec.n].rcv = ircv;
		tec.data[tec.n].sat = sat;
		tec.data[tec.n].sys = sys;
		
		ipp[2] = ionppp(pos, azel, Re, hion, ipp);										//ipp[2]:投影函数   pos:纬度、经度

		for (int i = 0; i < 3; i++) tec.data[tec.n].ipp[i] = ipp[i];
		tec.data[tec.n].std = std;
		tec.data[tec.n++].tec = stec1;	
	}
	fclose(fp);

	BL0[0] += pos[0] / (ircv + 1);
	BL0[1] += pos[1] / (ircv + 1);
	return 1;
}
static int generateData(const filopt_t *fopt)
{
	char stafile[256];
	unsigned char sflag[MAXSAT] = { 0 }, syslist;			
	if (!readstas(fopt->prosta, &stas)) return 0;

	trace(1,"Generate data....\n");
	for (int i = 0, ircv=0; i < stas.n; i++) {
		checkbrk("(%4d) %s", i, stas.data[i].name);
		sprintf(stafile, "%s%s%s", fopt->outdir, stas.data[i].name, "-gec-3.ppp");

		if (!file2tec(stafile, &stas.data[i], ircv, sflag, &syslist)) continue;			
	
		for (int k = 0; k < 4; k++) {
			if (syslist & listsys[k]) {
				rcv2ind[ircv][k] = nr;
				ind2rcv[nr] = ircv;
				ind2sys[nr++] = listsys[k];
			}
			else rcv2ind[ircv][k] = -1;
		}	
		rcv2sta[ircv++] = i;
	}

	/* init tec struct */
	for (int i = 0; i < MAXSAT; i++) {				
		if (sflag[i] != 0) {
			sat2ind[i] = ns;
			ind2sat[ns++] = i;
		}
		else sat2ind[i] = -1;
	}
	return 1;
}
static int estiono()		//这是建模的核心地带			
{
	gtime_t ts, te;
	double *N, *W, *b, *ipp;
	int sat, sys, flg, jj, kk, need, is, ir, iperiod, nperiod, info;
	double l_tec, tt, dt, fact, gamma, solar_mag[2], std;
	unsigned int tick = tickget();

	trace(1,"Estimate iono....\n");
	str2time(start_time, 0, 32, &ts);
	str2time(  end_time, 0, 32, &te);
	
	nperiod = time2ind(ts, te, nhour, NULL) + 1;
	flg = (order + 1)*(order + 1);
	need = ns + nr + flg * nperiod;
	N = zeros(need, need);		W = zeros(need, 1);		
	x = zeros(need,1);			b = zeros(need, 1);

	for (int i = 0; i < tec.n; i++) {				
		sat = tec.data[i].sat;
		is = sat2ind[sat - 1];
		for (int k = 0; k < 4; k++) {
			if (tec.data[i].sys & listsys[k]) {
				ir = rcv2ind[tec.data[i].rcv][k];
			}
		}
		ipp = tec.data[i].ipp;
		l_tec = tec.data[i].tec;
		std = tec.data[i].std;
		if (std > 0.1) continue;

		iperiod = time2ind(ts, tec.data[i].time, nhour, &dt);
		tt = (int)(tec.data[i].time.time % 86400) + tec.data[i].time.sec;

		/* 生成系数阵 */
		sys = satsys(sat, NULL);
		if (sys&SYS_GPS) {
			fact = 40.30E16 / SQR(FREQ1);
			gamma = SQR(FREQ2) / (SQR(FREQ2) - SQR(FREQ1));
		}
		else if (sys&SYS_GLO) {
			fact = 40.30E16 / SQR(FREQ1_GLO);
			gamma = SQR(FREQ2_GLO) / (SQR(FREQ2_GLO) - SQR(FREQ1_GLO));
		}
		else if (sys&SYS_GAL) {
			fact = 40.30E16 / SQR(FREQ1);
			gamma = SQR(FREQ5) / (SQR(FREQ5) - SQR(FREQ1));
		}
		else if (sys&SYS_CMP) {
			fact = 40.30E16 / SQR(FREQ1_CMP);
			gamma = SQR(FREQ3_CMP) / (SQR(FREQ3_CMP) - SQR(FREQ1_CMP));
		}
		else {
			continue;
		}

		b[is] = -1.0 / fact / ipp[2] * (-gamma);					
		b[ns + ir] = -1.0 / fact / ipp[2] * (-gamma);
		double p = 1 / (ipp[2] * ipp[2]);
		l_tec = l_tec / fact / ipp[2];

		geo2solarmag(ipp, tt, solar_mag);	
		Harmony(order, ns + nr, solar_mag, b, dt, iperiod);

		/* 法方程叠加（上三角阵） */
		N[is + is * need] += b[is] * b[is] * p;
		N[(ns + ir) + is*need] += b[is] * b[ns + ir] * p;
		N[(ns + ir) + (ns + ir)*need] += b[ns + ir] * b[ns + ir] * p;
		W[is] += b[is] * l_tec * p;
		W[ns + ir] += b[ns + ir] * l_tec * p;
		for (int j = 0; j < 2 * flg; j++) {
			jj = ns + nr + j + flg*iperiod;
			N[jj + is*need] += b[is] * b[jj]*p;
			N[jj + (ns + ir)*need] += b[ns + ir] * b[jj]*p;
			for (int k = 2 * flg - 1; k >= j; k--) {
				kk = ns + nr + k + flg*iperiod;
				N[kk + jj*need] += b[jj] * b[kk] * p;
			}
			W[jj] += b[jj] * l_tec*p;
		}
		checkbrk("(%3d)%s %s %4d %6.3f", ind2rcv[ir], stas.data[rcv2sta[ind2rcv[ir]]].name, time_str(tec.data[i].time, 0), iperiod, dt);
	}
	checkbrk("LS...\n");

	/* sDCB constraint */
	for (int i = 0; i < ns; i++) {								
		for (int j = ns - 1; j >= i; j--) {
			if (satsys(ind2sat[i], NULL) & satsys(ind2sat[j], NULL)) {
				N[j + i*need] += 1 / 1E-8;				
			}
		}
		W[i] += 0;
	}

	/* 下三角阵赋值 */
	for (int i = 1; i < need; i++) {
		for (int j = 0; j < i; j++) {
			N[j + i*need] = N[i + j*need];
		}
	}

	/* 法方程解算 */
	if ((info = matinv(N, need))) {
		trace(1, "lsq error info=%d", info);
		free(N); free(W); free(b); free(x);
		return 0;
	}
	matmul("NN", need, 1, need, 1.0, N, W, 0.0, x);							/* x=N^-1*W */
	trace(1, "estimate iono time=%.1f s\n", (tickget() - tick)*0.001);

	free(N); free(W); free(b);
	return 1;
}
static int outionogrid(const filopt_t *fopt, double *lat, double *lon, int nh)
{
	FILE *fp;
	gtime_t ts, te, ti;
	char outfile[256], id[32];
	double ep[6], BL[2], solar_mag[2], tt, dt, vtec_cal;
	int mm, nn, iperiod, nperiod;

	str2time(start_time, 0, 32, &ts);
	str2time(  end_time, 0, 32, &te);
	nperiod = time2ind(ts, te, nh, NULL) - 1;

	sprintf(outfile, "%s%s", fopt->outdir, "iono-g.txt");
	if (!(fp = fopen(outfile, "w"))) return 0;

	fprintf(fp, "     1.0            IONOSPHERE MAPS     GNSS                IONEX VERSION / TYPE\n");
	fprintf(fp, "%s\n", "Global ionospheric map produced by SGG                      DESCRIPTION         ");
	fprintf(fp, "%s\n", "Web page:       http://ybyao.sgg.whu.edu.cn/index.aspx      DESCRIPTION         ");
	fprintf(fp, "%6d%74.20s\n", nperiod, "# OF MAPS IN FILE   ");
	fprintf(fp, "%8.1f%63s\n", 6371.0, "BASE RADIUS");
	fprintf(fp, "%8.1f%6.1f%6.1f%60s\n", hion, hion, 0.0, "HGT1 / HGT2 / DHGT  ");
	fprintf(fp, "%8.1f%6.1f%6.1f%60s\n", lat[0], lat[1], lat[2], "LAT1 / LAT2 / DLAT  ");
	fprintf(fp, "%8.1f%6.1f%6.1f%60s\n", lon[0], lon[1], lon[2], "LON1 / LON2 / DLON  ");
	fprintf(fp, "%6d%74s\n", nexp, "EXPONENT            ");
	fprintf(fp, "%s\n", "DIFFERENTIAL CODE BIASES                                    START OF AUX DATA   ");

	/****************** sDCB ******************/
	for (int i = 0; i < ns; i++) {
		satno2id(ind2sat[i] + 1, id);
		fprintf(fp, "%6s%10.3f%10.3f%54s\n", id, x[i]/CLIGHT*1E9, 0.0, "PRN / BIAS / RMS    ");			
	}
	fprintf(fp, "%s\n", "DCB values in ns; zero-mean condition wrt satellite values  COMMENT             ");
	fprintf(fp, "%s\n", "DIFFERENTIAL CODE BIASES                                    END OF AUX DATA");
	fprintf(fp, "%73s\n", "END OF HEADER");

	/****************** vtec ******************/
	char startof[] = "START OF TEC MAP";
	char epochof[] = "EPOCH OF CURRENT MAP";
	char endof[] = "END OF TEC MAP";
	char title[] = "LAT/LON1/LON2/DLON/H";

	for (int i = 0; i < nperiod; i++) {
		ti = timeadd(ts, 3600 * nh * (i + 1));
		time2epoch(ti, ep);
		iperiod = time2ind(ts, ti, nhour, &dt);
		tt = (int)(ti.time % 86400) + ti.sec;
		fprintf(fp, "%6d%70s\n", i + 1, startof);
		fprintf(fp, "%6d%6d%6d%6d%6d%6d%44s\n", (int)ep[0], (int)ep[1], (int)ep[2], (int)ep[3], (int)ep[4], (int)ep[5], epochof);
		mm = (lat[1] - lat[0]) / lat[2] + 1;
		nn = (lon[1] - lon[0]) / lon[2] + 1;
		for (int m = 0; m < mm; m++) {
			BL[0] = (lat[1] - lat[2] * m)*D2R;
			fprintf(fp, "%8.1f%6.1f%6.1f%6.1f%6.1f%48s\n", BL[0]*R2D, lon[0], lon[1], lon[2], hion, title);
			for (int n = 0; n < nn; n++) {
				BL[1] = (lon[0] + lon[2] * n)*D2R;
				if (n % 16 == 0 && n != 0) {
					fprintf(fp, "\n");
				}
				geo2solarmag(BL, tt, solar_mag);
				vtec_cal = result_Harmony(iperiod, order, ns + nr, solar_mag, dt, x)*pow(10.0, -nexp);
				if (vtec_cal < 0) vtec_cal = 0;
				fprintf(fp, "%5d", (int)vtec_cal);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "%6d%68s\n", i + 1, endof);
	}
	fprintf(fp, "%80s\n", "END OF FILE         ");
	fclose(fp);

	return 1;
}
static int postiono(const filopt_t *fopt)   
{
	FILE *fp = NULL;
	gtime_t ts;
	char outfile[256], id[32];
	int sat, is, ir, iperiod, sys, ircv_1 = -1;
	double *ipp, l_tec, tt, dt, fact, gamma, solar_mag[2], vtec, vtec_cal, stec;

	unsigned int tick = tickget();

	str2time(start_time, 0, 32, &ts);

	for (int i = 0; i < tec.n; i++) {
		sat = tec.data[i].sat;
		is = sat2ind[sat - 1];
		for (int k = 0; k < 4; k++) {
			if (tec.data[i].sys & listsys[k]){
				ir = rcv2ind[tec.data[i].rcv][k];
			}
		}
		ipp = tec.data[i].ipp;
		l_tec = tec.data[i].tec;
		iperiod = time2ind(ts, tec.data[i].time, nhour, &dt);
		tt = (int)(tec.data[i].time.time % 86400) + tec.data[i].time.sec;

		/* 生成系数阵 */
		sys = satsys(sat, NULL);
		if (sys&SYS_GPS) {
			fact = 40.30E16 / SQR(FREQ1);
			gamma = SQR(FREQ2) / (SQR(FREQ2) - SQR(FREQ1));
		}
		else if (sys&SYS_GLO) {
			fact = 40.30E16 / SQR(FREQ1_GLO);
			gamma = SQR(FREQ2_GLO) / (SQR(FREQ2_GLO) - SQR(FREQ1_GLO));
		}
		else if (sys&SYS_GAL) {
			fact = 40.30E16 / SQR(FREQ1);
			gamma = SQR(FREQ5) / (SQR(FREQ5) - SQR(FREQ1));
		}
		else if (sys&SYS_CMP) {
			fact = 40.30E16 / SQR(FREQ1_CMP);
			gamma = SQR(FREQ3_CMP) / (SQR(FREQ3_CMP) - SQR(FREQ1_CMP));
		}
		else {
			continue;
		}
		vtec = (l_tec - gamma*(x[is] + x[ns + ir])) / fact / ipp[2];	/* 单位tecu */
		geo2solarmag(ipp, tt, solar_mag);
		vtec_cal = result_Harmony(iperiod, order, ns + nr, solar_mag, dt, x);
		stec = fact*ipp[2] * vtec_cal + gamma*(x[is] + x[ns + ir]);
		if (ircv_1 != ind2rcv[ir]) {
			sprintf(outfile, "%s%s-g1.iono", fopt->outdir, stas.data[rcv2sta[ind2rcv[ir]]].name);
			if (fp == NULL) {
				if (!(fp = fopen(outfile, "w"))) return 0;
			}
			else {
				fclose(fp);
				if (!(fp = fopen(outfile, "w"))) return 0;
			}
			ircv_1 = ind2rcv[ir];
		}
		satno2id(ind2sat[is] + 1, id);
		fprintf(fp, "%s %6s %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n", time_str(tec.data[i].time, 0), id, vtec, vtec_cal, fabs(vtec - vtec_cal), l_tec, stec, fabs(l_tec - stec));
	}
	trace(1, "post iono time=%.1f s\n", (tickget() - tick)*0.001);
	fclose(fp);
	return 1;
}
static int drawipp(const filopt_t *fopt)
{
	FILE *fp = NULL;
	gtime_t ts;
	char outfile[256];
	int iperiod, ir;
	double *ipp, tt, solar_mag[2], dt;
	sprintf(outfile, "%sipp.txt", fopt->outdir);
	if (!(fp = fopen(outfile, "w"))) return 0;

	str2time(start_time, 0, 32, &ts);
	for (int i = 0; i < tec.n; i++) {
		ipp = tec.data[i].ipp;
		ir = tec.data[i].rcv;
		tt = (int)(tec.data[i].time.time % 86400) + tec.data[i].time.sec;
		iperiod = time2ind(ts, tec.data[i].time, nhour, &dt);
		geo2solarmag(ipp, tt, solar_mag);
		if (ir == 20) 
			fprintf(fp, "%10.4f %10.4f\n", ipp[0] * R2D, ipp[1] * R2D);
	}
	fclose(fp);
	return 1;
}
extern int ionoModel(const filopt_t *fopt)	
{			
	double lat[3] = { 37, 57, 0.5 }, lon[3] = { -6, 30, 1.0 };
	//double lat[3] = { 87.5, -87.5, -2.5 }, lon[3] = { -180.0, 180.0, 5.0 };

	/* read station list */
	if (!generateData(fopt)) return 0;				//	

	/* process */
	if (!estiono()) return 0;				

	/* output */
	if (!outionogrid(fopt, lat, lon, 1)) return 0;			

	/* post residuals */
	if (!postiono(fopt)) return 0;

	/* draw ipp */
	if (!drawipp(fopt)) return 0;

	free(x);	freetec(&tec);	  freestas(&stas);
	return 1;
}

//FILE *fp;
//char infile[256] = "E:\\GNSS_DATA\\Result\\Iono-11\\iono-g.txt";
////char infile[256] = "E:\\GNSS_DATA\\products\\2021-07-10\\codg1910.21i";
//readtec(infile, &nav, 0);
//char outfile[256];
//double BL[2] = { 0.0 }, vtec_cal;
//int index;
//for (int i = 0; i < nav.nt; i++) {
//	sprintf(outfile, "%s%d%s", fopt->outdir, 3 + i, "vtec.txt");
//	fp = fopen(outfile, "w");
//	for (int m = 0; m < nav.tec[i].ndata[0]; m++) {
//		BL[0] = (nav.tec[i].lats[0] + nav.tec[i].lats[2] * m)*D2R;
//		for (int n = 0; n < nav.tec->ndata[1]; n++) {
//			BL[1] = (nav.tec[i].lons[0] + nav.tec[i].lons[2] * n)*D2R;
//			if ((index = dataindex(m, n, 0, nav.tec[i].ndata)) < 0) continue;
//			vtec_cal = nav.tec[i].data[index];
//			fprintf(fp, "%6.3f ", vtec_cal);
//		}
//		fprintf(fp, "\n");
//	}
//	fclose(fp);
//}
//return 0;
