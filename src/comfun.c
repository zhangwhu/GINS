
#include "rtklib.h"
#include "ppp_state.h"

#define SQR(x)      ((x)*(x))
#define ERR_SAAS    0.3         /* saastamoinen model error std (m) */
#define ERR_ION     5.0         /* ionospheric delay std (m) */
#define ERR_TROP    3.0         /* tropspheric delay std (m) */
#define ERR_BRDCI   0.5         /* broadcast iono model error factor */
#define MIN_EL      (5.0*D2R)   /* min elevation for measurement error (rad) */

/* nominal yaw-angle ---------------------------------------------------------*/
static double yaw_nominal(double beta, double mu)
{
	if (fabs(beta)<1E-12&&fabs(mu)<1E-12) return PI;
	return atan2(-tan(beta), sin(mu)) + PI;
}
/* yaw-angle of satellite ----------------------------------------------------*/
static int yaw_angle(int sat, const char *type, int opt, double beta, double mu,
	double *yaw)
{
	*yaw = yaw_nominal(beta, mu);
	return 1;
}
/* satellite attitude model --------------------------------------------------*/
static int sat_yaw(gtime_t time, int sat, const char *type, int opt,
	const double *rs, double *exs, double *eys)
{
	double rsun[3], ri[6], es[3], esun[3], n[3], p[3], en[3], ep[3], ex[3], E, beta, mu;
	double yaw, cosy, siny, erpv[5] = { 0 };
	int i;

	sunmoonpos(gpst2utc(time), erpv, rsun, NULL, NULL);

	/* beta and orbit angle */
	matcpy(ri, rs, 6, 1);
	ri[3] -= OMGE*ri[1];
	ri[4] += OMGE*ri[0];
	cross3(ri, ri + 3, n);
	cross3(rsun, n, p);
	if (!normv3(rs, es) || !normv3(rsun, esun) || !normv3(n, en) ||
		!normv3(p, ep)) return 0;
	beta = PI / 2.0 - acos(dot(esun, en, 3));
	E = acos(dot(es, ep, 3));
	mu = PI / 2.0 + (dot(es, esun, 3) <= 0 ? -E : E);
	if (mu<-PI / 2.0) mu += 2.0*PI;
	else if (mu >= PI / 2.0) mu -= 2.0*PI;

	/* yaw-angle of satellite */
	if (!yaw_angle(sat, type, opt, beta, mu, &yaw)) return 0;

	/* satellite fixed x,y-vector */
	cross3(en, es, ex);
	cosy = cos(yaw);
	siny = sin(yaw);
	for (i = 0; i<3; i++) {
		exs[i] = -siny*en[i] + cosy*ex[i];
		eys[i] = -cosy*en[i] - siny*ex[i];
	}
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
				for (p = l_name, q = (char *)name; (*p = (char)toupper(*q)); p++, q++);					//字符串->大写字母
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
		fclose(fp);
		free(stas.data);	stas.data = NULL;	stas.n = stas.nmax = 0;
		return;
	}
	if (infile) uniqsnx(infile, &stas);

	for (int i = 0; i < stas.n; i++) {
		if (stas.data[i].name[0] != '\0') {
			/*double blh[3];
			ecef2pos(stas.data[i].pos, blh);
			fprintf(fp_out, "%s %15.6f %15.6f\n", stas.data[i].name, blh[0] * 180 / PI, blh[1] * 180 / PI);*/
			fprintf(fp_out, "%s %15.6f %15.6f %15.6f\n", stas.data[i].name, stas.data[i].pos[0], stas.data[i].pos[1], stas.data[i].pos[2]);
		}
	}

	fclose(fp);			fclose(fp_out);
	free(stas.data);	stas.data = NULL;	stas.n = stas.nmax = 0;
}
/* generate stationlist */
extern int readstas(const char *file, stas_t *stas)
{
	FILE *fp;
	char buff[256];
	stas_t *stas_data;

	if (!(fp = fopen(file, "r"))) return 0;
	while (fgets(buff, sizeof(buff), fp)) {
		if (stas->nmax <= stas->n) {
			stas->nmax += 100;
			if (!(stas_data = (sta_t *)realloc(stas->data, sizeof(sta_t)*stas->nmax))) {
				free(stas->data);	stas->data = NULL;		stas->n = stas->nmax = 0;
				return  0;
			}
			stas->data = stas_data;
		}
		if (sscanf(buff, "%s %lf %lf %lf\n", stas->data[stas->n].name, &stas->data[stas->n].pos[0], &stas->data[stas->n].pos[1], &stas->data[stas->n].pos[2]) < 4) {
			if (sscanf(buff, "%s\n", stas->data[stas->n].name) < 1) {
				continue;
			}
		}
		stas->data[stas->n++].name[4] = '\0';
	}
	return stas->n;
}
extern double bd2smp(int orb, double *azel, int nq) {
	int i, el;
	double v[3][2], mp = 0.0;

	el = (int)((azel[1] * R2D) / 10) * 10;
	for (i = 0; i < 20; i++) {
		if (bdsmptable[i][0] == orb && bdsmptable[i][1] == el) {
			v[nq][0] = bdsmptable[i][nq + 2];
			v[nq][1] = bdsmptable[i + 1][nq + 2];
			mp = (azel[1] * R2D - el) / 10 * v[nq][1] + ((el + 10) - azel[1] * R2D) / 10 * v[nq][0];
			break;
		}
	}
	//trace(1, "prn=%4d freq=%4d azel=%10.4f el=%4d %6.3f %6.3f %6.3f\n", prn, nq, azel[1] * R2D, el, mp, v[nq][0], v[nq][1]);
	return mp;
}
/* set antenna parameters ----------------------------------------------------*/
extern void setpcv(gtime_t time, prcopt_t *popt, nav_t *nav, const pcvs_t *pcvs, const pcvs_t *pcvr, const sta_t *sta)
{
	int i, j;
	char id[64];
	double pos[3], del[3];
	pcv_t *pcv;

	/* set satellite antenna parameters */
	for (i = 0; i<MAXSAT; i++) {
		if (!(satsys(i + 1, NULL)&popt->navsys)) continue;
		if (!(pcv = searchpcv(i + 1, "", time, pcvs))) {
			satno2id(i + 1, id);
			trace(3, "no satellite antenna pcv: %s\n", id);
			continue;
		}
		nav->pcvs[i] = *pcv;
	}

	*popt->anttype[0] = '\0';
	for (j = 0; j<3; j++) popt->antdel[0][j] = 0;

	strcpy(popt->anttype[0], sta->antdes);
	if (sta->deltype == 1) {							/* xyz */
		if (norm(sta->pos, 3)>0.0) {
			ecef2pos(sta->pos, pos);
			ecef2enu(pos, sta->del, del);
			for (j = 0; j<3; j++) popt->antdel[0][j] = del[j];
		}
	}
	else {												/* enu */
		for (j = 0; j<3; j++) popt->antdel[0][j] = sta->del[j];
	}
	if (!(pcv = searchpcv(0, popt->anttype[0], time, pcvr))) {
		trace(2, "no receiver antenna pcv: %s\n", popt->anttype[i]);
		*popt->anttype[i] = '\0';
		return;
	}
	strcpy(popt->anttype[0], pcv->type);
	popt->pcvr[0] = *pcv;
}
/* phase windup model -------------------------------------------------------- */
extern int model_phw(gtime_t time, int sat, const char *type, int opt,
	const double *rs, const double *rr, double *phw)
{
	double exs[3], eys[3], ek[3], exr[3], eyr[3], eks[3], ekr[3], E[9];
	double dr[3], ds[3], drs[3], r[3], pos[3], cosp, ph;
	int i;

	if (opt <= 0) return 1; /* no phase windup */

	/* satellite yaw attitude model */
	if (!sat_yaw(time, sat, type, opt, rs, exs, eys)) return 0;		//卫星姿态引起的相位缠绕影响

	/* unit vector satellite to receiver */
	for (i = 0; i<3; i++) r[i] = rr[i] - rs[i];
	if (!normv3(r, ek)) return 0;

	/* unit vectors of receiver antenna */
	ecef2pos(rr, pos);
	xyz2enu(pos, E);
	exr[0] = E[1]; exr[1] = E[4]; exr[2] = E[7]; /* x = north */
	eyr[0] = -E[0]; eyr[1] = -E[3]; eyr[2] = -E[6]; /* y = west  */

	/* phase windup effect */
	cross3(ek, eys, eks);
	cross3(ek, eyr, ekr);
	for (i = 0; i<3; i++) {
		ds[i] = exs[i] - ek[i] * dot(ek, exs, 3) - eks[i];			//卫星天线有效偶极向量
		dr[i] = exr[i] - ek[i] * dot(ek, exr, 3) + ekr[i];			//接收机天线有效偶极向量
	}
	cosp = dot(ds, dr, 3) / norm(ds, 3) / norm(dr, 3);
	if (cosp<-1.0) cosp = -1.0;
	else if (cosp> 1.0) cosp = 1.0;
	ph = acos(cosp) / 2.0 / PI;
	cross3(ds, dr, drs);
	if (dot(ek, drs, 3)<0.0) ph = -ph;

	*phw = ph + floor(*phw - ph + 0.5); /* in cycle */		//floor(*phw - ph + 0.5)
	return 1;
}
/* troposphere correction ----------------------------------------------------- */
extern int tropcorr(gtime_t time, const double *pos, const double *azel, const prcopt_t *opt, 
	const double *x, double *dtdx, const nav_t *nav, double *dtrp, double *var)
{
	const double zazel[] = { 0.0, PI / 2.0 };
	double trp[3] = { 0 };
	double zhd, m_h, m_w, cotz, grad_n, grad_e;

	if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG) {
		matcpy(trp, x + IT(opt), opt->tropopt == TROPOPT_EST ? 1 : 3, 1);

		/* zenith hydrostatic delay */
		zhd = tropmodel(time, pos, zazel, 1, REL_HUMI, 1);

		/* mapping function */
		m_h = tropmapf(time, pos, azel, 1, &m_w);					

		/* gradient */
		if (opt->tropopt == TROPOPT_ESTG && azel[1]>0.0) {
			cotz = 1.0 / tan(azel[1]);
			grad_n = m_w*cotz*cos(azel[0]);
			grad_e = m_w*cotz*sin(azel[0]);
			m_w += grad_n*trp[1] + grad_e*trp[2];
			dtdx[1] = grad_n*trp[0];								
			dtdx[2] = grad_e*trp[0];
		}
		dtdx[0] = m_w;
		*dtrp = m_h*zhd + m_w*trp[0];
		*var = 0.0;
		return 1;
	}
	/* troposphere ZTD correction*/
	if (opt->tropopt == TROPOPT_ZTD) {									

	}
	/* saastamoinen model */
	if (opt->tropopt == TROPOPT_SAAS) {
		*dtrp = tropmodel(time, pos, azel, 1, REL_HUMI, 0);
		*var = SQR(ERR_SAAS / (sin(azel[1]) + 0.1));
		return 1;
	}
	/* no correction */
	if (opt->tropopt == TROPOPT_OFF) {				
		*dtrp = 0.0;
		*var = SQR(ERR_TROP);
		return 1;
	}
	
	return 0;
}
/* ionospheric correction ----------------------------------------------------- */
extern int ionocorr(gtime_t time, const double *pos, const double *azel, const prcopt_t *opt, int sat, const double *x,
	const nav_t *nav, double *dion, double *var)
{
	if (opt->ionoopt == IONOOPT_TEC) {
		return iontec(time, nav, pos, azel, 1, dion, var);				//斜延迟(L1频段)
	}
	if (opt->ionoopt == IONOOPT_BRDC) {
		*dion = ionmodel(time, nav->ion_gps, pos, azel);
		*var = SQR(*dion*ERR_BRDCI);
		return 1;
	}
	if (opt->ionoopt == IONOOPT_EST) {									//待估参数斜延迟	
		*dion = x[II(sat, opt)];
		*var = 0.0;
		return 1;
	}
	if (opt->ionoopt == IONOOPT_IFLC) {
		*dion = *var = 0.0;
		return 1;
	}
	if (opt->ionoopt == IONOOPT_STEC) {

	}
	return 0;
}
/* observation correction ----------------------------------------------------- */
extern void corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel, const prcopt_t *opt, const double *dantr,
	const double *dants, double phw, const double php, double *L, double *P, double *Lc, double *Pc)
{
	int i, k, sys, prn, orb, ix, nf;
	double freq[NFREQ], c1, c2, C1, C2, gamma;
	nf = opt->ionoopt == IONOOPT_IFLC ? 1 : opt->nf;
	for (i = 0; i < NFREQ; i++) {
		freq[i] = sat2freq(obs->sat, obs->code[i], nav);
		L[i] = P[i] = 0.0;
	}
	sys = satsys(obs->sat, &prn);

	for (i = 0; i< nf; i++) {
		if (freq[i] == 0.0 || obs->L[i] == 0.0 || obs->P[i] == 0.0) continue;
		if (testsnr(0, i, azel[1], obs->SNR[i] * SNR_UNIT, &opt->snrmask)) continue;

		/* antenna phase center and phase windup correction */
		L[i] = CLIGHT * obs->L[i] / freq[i] - dants[i] - dantr[i] - CLIGHT * phw / freq[i] + php;
		P[i] = obs->P[i] - dants[i] - dantr[i];

		/* P1-C1,P2-C2 dcb correction (C1->P1,C2->P2) */
		if (sys == SYS_GPS || sys == SYS_GLO) {
			if (obs->code[i] == CODE_L1C) P[i] += nav->cbias[obs->sat - 1][1];
			if (obs->code[i] == CODE_L2C) P[i] += nav->cbias[obs->sat - 1][2];
		}
		else if (sys == SYS_CMP) {
			if (prn < 19) {
				satorb(obs->sat, &orb);
				P[i] += bd2smp(orb, azel, i);
			}
		}
		else if (sys == SYS_GAL) {
			if (obs->code[i] == CODE_L1X) P[i] += nav->cbias[obs->sat - 1][1];
			if (obs->code[i] == CODE_L5X) P[i] += nav->cbias[obs->sat - 1][2];
		}
	}

	/* iono-free LC */				
	*Lc = *Pc = 0.0;
	if (freq[0] == 0.0 || freq[1] == 0.0) return;

	C1 =  SQR(freq[0]) / (SQR(freq[0]) - SQR(freq[1]));
	C2 = -SQR(freq[1]) / (SQR(freq[0]) - SQR(freq[1]));

	if (L[0] != 0.0&&L[1] != 0.0) *Lc = C1*L[0] + C2*L[1];
	if (P[0] != 0.0&&P[1] != 0.0) *Pc = C1*P[0] + C2*P[1];
}

/* measurement error variance ------------------------------------------------*/
extern double varerr(unsigned char sat, double el, int freq, int type, const prcopt_t *opt)
{
	int sys, prn, orb;
	double fact = 1.0, sinel, varr, a;

	if (el<MIN_EL) el = MIN_EL;

	sys = satsys(sat, &prn);
	satorb(sat, &orb);
	sinel = sin(el);

	if (sys == SYS_GPS) {
		a = (type == 1) ? 0.3 : 0.003;
		if (freq == 2) fact *= EFACT_GPS_L5;		/* GPS/QZS L5 error factor */ 
	}
	if (sys == SYS_GLO) {
		a = (type == 1) ? 0.3 : 0.003;
		fact *= EFACT_GLO;
	}
	if (sys == SYS_GAL) {
		a = (type == 1) ? 0.3 : 0.003;
	}
	if (sys == SYS_CMP) {
		if (prn < 19) {
			if (orb == SAT_GEO) a = (type == 1) ? 0.6 : 0.009;
			if (orb == SAT_MEO || orb == SAT_IGSO) a = (type == 1) ? 0.6 : 0.003;
		}
		if (prn > 18) {
			a = (type == 1) ? 0.3 : 0.003;
		}
	}

	varr = SQR(a);
	if (el < 30 * D2R) varr /= SQR(2 * sinel);
	if (opt->ionoopt == IONOOPT_IFLC) fact *= 3.0;

	return SQR(fact)*varr;
}
