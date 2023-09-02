
#include "CGNSSManage.h"

/***************************  class CGNSSManage  *********************************/
CGNSSManage::CGNSSManage() 
{	
	stat = &rtk.sol.stat;
	ti   = &rtk.sol.time;

	FILE *fp;
	char buff[256];
	int q;
	double ti, pos[3], vel[3], att[3], eb[3], db[3];
	if (!(fp = fopen("/media/zhuhang/D/Data/PPPRTK_INS_VISUAL_Wuh-230416/Result/TC/gvins-adis-1.txt", "r"))) {
		trace(2, "gvins file open error\n");
		return;
	}
	while (fgets(buff, sizeof(buff), fp)) {
		if (sscanf(buff,"%lf%d%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf\n", 
			&ti, &q, &pos[0], &pos[1], &pos[2], &vel[0], &vel[1], &vel[2],
			&att[0], &att[1], &att[2], &eb[0], &eb[1], &eb[2], &db[0], &db[1], &db[2]) < 17) continue;
		if (solss.n >= solss.nmax) {
			sol *pd;
			solss.nmax += solss.nmax <= 0 ? 300 : 50;
			if (!(pd = (sol *)realloc(solss.data, sizeof(sol)*(solss.nmax)))) {
				trace(1, "readsol sol malloc error: nmax=%d\n", solss.nmax);
				free(solss.data);	solss.data = NULL; solss.n = solss.nmax = 0;
				break;
			}
			solss.data = pd;
		}
		sol data;
		data.ti = ti;
		data.rr[0] = pos[0];  data.rr[1] = pos[1];  data.rr[2] = pos[2]; 
		data.qr[0] = 0.05/RE_WGS84; 	data.qr[1] = 0.05/RE_WGS84;   data.qr[2] = 0.1; 
		solss.data[solss.n++] = data;
	}
	fclose(fp);
}

CGNSSManage::~CGNSSManage() 
{	
	rtkfree(&rtk);
	freeproduct(&nav, NULL, NULL, NULL);
	if (fp_outs[0])  { fflush(fp_outs[0]); fclose(fp_outs[0]); }  else { fp_outs[0] = NULL; }
	if (fp_outs[1])  { fflush(fp_outs[1]); fclose(fp_outs[1]); }  else { fp_outs[1] = NULL; }
	traceclose();
}

int CGNSSManage::Init(const prcopt_t &popt, const solopt_t &sopt, const filopt_t &fopt, 
	const sta_t &sta) 
{	
	prcopt_t popt_ = popt;
	pcvs_t pcvss = { 0 }, pcvsr = { 0 };
	char outfile1[1024], outfile2[1024], outfile[1024], path[1024];
	
	/* read product */
	if (!readproduct(&popt_, &fopt, &nav, &pcvss, &pcvsr, NULL)) return 0;
	
	/* read ppp corrections */
	reppath(fopt.corr, path, popt_.ts, "", "");
	pppcorr_read(path, &nav);

	/* set antenna &&  ocean tide && ref position */
	setpcv(popt_.ts, &popt_, &nav, &pcvss, &pcvsr, &sta);		// sta正常情况应该是数组

	/* remove unused product data */
	removeUnusedData(popt_.ts, popt_.te, &nav, &pcvss, &pcvsr);

	/* rtk init */
	rtkinit(&rtk, &popt_);
	
	/* output file */
	solopt = sopt;
	strcpy(outdir, fopt.outdir);
	sprintf(outfile1, "%s%s", outdir, fopt.outfile1);
	sprintf(outfile2, "%s%s", outdir, fopt.outfile2);
	reppath(outfile1, outfile, popt_.ts, sta.name, "");
	if (fp_outs[0] = fopen(outfile, "w")) {
		outheader(fp_outs[0], &popt, &sopt);	
	}
	reppath(outfile2, outfile, popt_.ts, sta.name, "");
	if (fp_outs[1] = fopen(outfile, "w")) {
		fprintf(fp_outs[1], "$NAME,%s\n", sta.name);
	}
	if (!fp_outs[0] || !fp_outs[1]) return 0;

	return 1;
}

void CGNSSManage::Update(obsd_t *obs, int nobs, int bimu2gps, const double *ins_xyz, const double *ins_Pxyz) 
{
	const prcopt_t *prcopt_ = &rtk.opt;
	static gtime_t time0 = { 0 };
	int i, n, sat;
	double blh[3], venu[3], xg[3] = { 0.0 }, Pg[3*3] = { 0.0 }, P[3*3] = { 0.0 };
	char msg[128] = "";
	
	/* exclude satellites */
	for (i = n = 0; i < nobs; i++)		
	{
		sat = obs[i].sat;
		if (!(rtk.ssat[sat - 1].sys&prcopt_->navsys))continue;
		if (prcopt_->exsats[sat - 1]) continue;
		if (ISZERO(obs[i].P[0]) || ISZERO(obs[i].L[0])) continue;
		// if (ISZERO(obs[i].P[1]) && !ISZERO(obs[i].P[2])) {
		// 	char id[8];
		// 	satno2id(obs[i].sat, id);
		// 	trace(1,"P %s %s\n",time_str(obs[i].time,0), id);
		// }
		// if (ISZERO(obs[i].L[1]) && !ISZERO(obs[i].L[2])) {
		// 	char id[8];
		// 	satno2id(obs[i].sat, id);
		// 	trace(1,"L %s %s\n",time_str(obs[i].time,0), id);
		// }

		if (rtk.ssat[sat - 1].sys&SYS_GPS) {
			obs[i].P[2] = obs[i].L[2] = 0.0;
		}
		obs[n++] = obs[i];
	}
	if (n <= 0) return;
	
	/* rover position by single point positioning */
	if (!pntpos(obs, n, &nav, &rtk.opt, &rtk.sol, NULL, rtk.ssat, msg))
	{	
		rtk.sol.time = obs[0].time;			
		showmsg(msg);
		return;
	}
	
	// if (time.time != 0) rtk.tt = fabs(timediff(rtk.sol.time, time));
	if (time0.time != 0) rtk.tt = fabs(timediff(rtk.sol.time, time0));
	time0 = rtk.sol.time;				// 上解决方案的时刻

	/* precise point positioning */
	if (prcopt_->mode >= PMODE_PPP_KINEMA)
	{			
		if (bimu2gps && (ins_xyz) && (ins_Pxyz)) 
		{	
			// pos2ecef(blh, xg);	 covecef(blh, P, Pg);  
            matcpy(xg,  ins_xyz, 1, 3); 
            matcpy(Pg, ins_Pxyz, 3, 3);
		}
		ppppos(&rtk, obs, n, &nav, xg, Pg);			
	}

	/* output result */
	outsol(fp_outs[0], &rtk.sol, NULL, &solopt);	
	outppp(fp_outs[1], &rtk, obs, n, &nav);			
	ecef2pos(rtk.sol.rr, blh);			matcpy(pos_blh,  blh, 1, 3);
	ecef2enu(blh,rtk.sol.rr+3, venu);	matcpy(vel_enu, venu, 1, 3);
	// if (*stat == SOLQ_FIX || *stat == SOLQ_PPP) {
	// 	double dr[3], enu[3];
	// 	double ti = time2gpst(rtk.sol.time, NULL); 
	// 	for (int i = 0; i < 3; i++) dr[i] = xg[i] - rtk.sol.rr[i];
	// 	ecef2enu(blh, dr, enu);
	// 	trace(1,"%10.2f %4d %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n", 
	// 		ti, *stat, enu[0], enu[1], enu[2], 
	// 		xg[0], xg[1], xg[2], rtk.sol.rr[0], rtk.sol.rr[1], rtk.sol.rr[2]);
	// }
}

void CGNSSManage::TraceOpen(const char *filename) 
{	
	char file[1024];

	sprintf(file, "%s%s", outdir, filename);
	traceopen(file);	
	tracelevel(1);
}

void CGNSSManage::Sol2kml()
{	
	sol2kml(outdir);
}

void CGNSSManage::OutGPSResult() 			
{	
	if (*stat == SOLQ_FIX) {
        float std[] = { 0.05, 0.1};						
		rk_blh[0] = rk_blh[1] = std[0]/RE_WGS84;    rk_blh[2] = std[1];
	}
	else if(*stat == SOLQ_PPP) {
        float std[] = { 0.3,  0.6};     // 0.3 0.4
		rk_blh[0] = rk_blh[1] = std[0]/RE_WGS84;    rk_blh[2] = std[1];   
	}
	else if (*stat == SOLQ_SINGLE) {
        float std[] = { 1, 10};
		rk_blh[0] = rk_blh[1] = std[0]/RE_WGS84;    rk_blh[2] = std[1];
	}
	else {
        float std[] = { 1000, 2000};
		rk_blh[0] = rk_blh[1] = std[0]/RE_WGS84;    rk_blh[2] = std[1];
	}
}

void CGNSSManage::Update(gtime_t time) 
{
	double ti = time2gpst(time, NULL);
	for (int k = 0; k < solss.n; k++) {
		if (fabs(solss.data[k].ti - ti) < 0.1) {
			pos_blh[0] = solss.data[k].rr[0];
			pos_blh[1] = solss.data[k].rr[1];
			pos_blh[2] = solss.data[k].rr[2];
			rk_blh[0]  = solss.data[k].qr[0]; 
			rk_blh[1]  = solss.data[k].qr[1]; 
			rk_blh[2]  = solss.data[k].qr[2]; 
			return ;
		}
	}
}

static int screent(gtime_t time, gtime_t ts, gtime_t te)
{
	return ((ts.time == 0 || timediff(time, ts) >= -DTTOL) &&
		    (te.time == 0 || timediff(time, te) <  DTTOL));
}

void removeUnusedData(gtime_t &ts, gtime_t &te, nav_t *nav, pcvs_t *pcvss, pcvs_t *pcvsr) 
{
	int i, k, ind, ii;
	peph_t *nav_peph;
	pclk_t *nav_pclk;
	stec_t *nav_stec;
	trop_t *nav_trop;
	gtime_t tss;
	gtime_t tee;

	/* 删除冗余的卫星和接收机天线数据 */
	freeproduct(NULL, pcvss, pcvsr, NULL);
	
	/* 删除冗余的精密星历数据 */
	tss = timeadd(ts, -1500);
	tee = timeadd(te,  1500);
	for (i = 0, k = 0; k < nav->ne; k++) {
		if (screent(nav->peph[k].time, tss, tee)) {			
			nav->peph[i++] = nav->peph[k];
		}
	}
	nav->ne = i;
	if (!(nav_peph = (peph_t *)realloc(nav->peph, sizeof(peph_t)*nav->ne))) {
		free(nav->peph); nav->peph = NULL; nav->ne = nav->nemax = 0;
		trace(1, "combpeph malloc error ne=%d\n", nav->ne);
	} else {
		nav->peph = nav_peph;		nav->nemax = nav->ne;
	}
	
	/* 删除冗余的精密钟差数据 */
	tss = timeadd(ts, -900);
	tee = timeadd(te,  900);
	for (i = 0, k = 0; k < nav->nc; k++) {
		if (screent(nav->pclk[k].time, tss, tee)) {			
			nav->pclk[i++] = nav->pclk[k];
		}
	}
	nav->nc = i;
	if (!(nav_pclk = (pclk_t *)realloc(nav->pclk, sizeof(pclk_t)*nav->nc))) {
		free(nav->pclk); nav->pclk = NULL; nav->nc = nav->ncmax = 0;
		trace(1, "combpclk malloc error nc=%d\n", nav->nc);
	} else {
		nav->pclk = nav_pclk;		nav->ncmax = nav->nc;
	}

	/* 删除冗余的区域增强-对流层数据 */
	tss = timeadd(ts, -60);
	tee = timeadd(te,  60);
	for (int n = 0; n < nav->ntrop; n++) {
		corrtrop_t *ptrop = nav->corrtrop + n;
		for (i = 0, k = 0; k < ptrop->nt; k++) {
			if (screent(ptrop->data[k].time, tss, tee)) {			
				ptrop->data[i++] = ptrop->data[k];
			}
		}
		ptrop->nt = i;
		if (!(nav_trop = (trop_t *)realloc(ptrop->data, sizeof(trop_t)*ptrop->nt))) {
			free(ptrop->data); 		ptrop->data = NULL; 	ptrop->nt = ptrop->ntmax = 0;
			trace(1, "combtrop malloc error nt=%d\n", ptrop->nt);
		} else {
			ptrop->data = nav_trop;			ptrop->ntmax = ptrop->nt;
		}
	}
	
	/* 删除冗余的区域增强-电离层数据 */
	for (int n = 0; n < nav->nstec; n++) {
		corrstec_t *pstec = nav->corrstec + n;
		for (i = ind = 0, k = 0; k < pstec->nn; k++) {
			if (screent(pstec->time[k], tss, tee)) {		  	
				pstec->time[i] = pstec->time[k];
				pstec->ii[i++] = ind;
				int endind = (k + 1 < pstec->nn) ? pstec->ii[k + 1] : pstec->ns;
				for (ii = pstec->ii[k]; ii < endind; ii++) {
					pstec->data[ind++] = pstec->data[ii];
				}
			}
		}
		pstec->nn = i;				pstec->ns = ind;
		if (!(nav_stec = (stec_t *)realloc(pstec->data, sizeof(stec_t)*pstec->ns))) {
			free(pstec->data); 		pstec->data = NULL; 	pstec->ns = pstec->nsmax = 0;
			trace(1, "combstec malloc error ns=%d\n", pstec->ns);
		} else {
			pstec->data  = nav_stec;		pstec->nsmax = pstec->ns;
		}
	}
}