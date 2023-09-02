
#include "CKFApp.h"

/***************************  class CMyAutoDrive  *********************************/
CMyAlignment::CMyAlignment(int mode, double T1, double T2):CSINSGNSS(15, 4, 0.02)
{	
	this->alnMode = mode;		this->T1 = T1;		this->T2 = T2;

	alnOK = 0;		qnb = qI;	

	posInit = alnkfInit = 0;

	vel0 = 1.0;	  wz0 = 5.0 * DPS; 	 dyaw0 = 5.0 * DEG;

	Hk(0,3) = Hk(1,4) = Hk(2,5) = 1.0;					// vn 
	Hk(3,6) = Hk(4,7) = Hk(5,8) = 0.0;	Hk(3,2) = 1.0;	// dyaw
	SetMeasMask(1, 017);   					
}

CMyAlignment::~CMyAlignment() 
{

}

void CMyAlignment::SetParam(double vel00, double wz00, double dyaw00) 
{
	vel0 = vel00;		wz0 = wz00; 	dyaw0 = dyaw00;
}

void CMyAlignment::Init(const CSINS &sins0) 		
{	
	alni0.Init(sins0.pos);	 CSINSTDKF::Init(sins0); 	

	posInit = 1;
	velPre = yawPre = cntYaw = alnkfInit = 0; 		sins.mvnT = 0.1;

	/* Pk */				
	Pmax.Set2(fDEG3(30.0), fXXX(50.0), fdPOS(1.0e4), fDPH3(10.0), fMG3(10.0));
	Pmin.Set2(fXXZU(1.0,10.0, SEC), fXXX(0.001), fXXX(0.01), fDPH3(0.001), fUG3(1.0));
	Pk.SetDiag2(fXXZU(1.10,10.0, DEG), fXXX(1.0), fXXX(100.0), fDPH3(0.05), fUG3(100.0));

	/* Qt */
	Qt.Set2(fDPSH3(0.001), fUGPSHZ3(1.0), fXXX(0.0), fXXX(0.0), fXXX(0.0));

	/* Xk */
	Xmax.Set(fINF9, fDPH3(0.1), fMG3(1.0));

	/* Rt */
	Rt.Set2(fXXX(0.1), fdPOS(100.0));
	Rmax = Rt*100;  Rmin = Rt*0.01;  Rb = 0.5;

	/* FB */
	FBTau.Set(fXXX(0.1), fXXX(0.1), fINF3, fXXX(1.1), fXXX(1.1));
}

void CMyAlignment::SetMeasGNSS(const CVect3 &pos, const CVect3 &vnr) 
{	
	if (!posInit) Init(pos);

	if (alnMode == STATIC_FINE && alnkfInit) {
		if(norm(vnr)>0)	{
			*(CVect3*)&Zk.dd[0] = sins.vn-vnr;			
			if(norm(sins.wnb)<1.1*glv.dps) {
				SetMeasFlag(0007);
			}
		}
	}
	if (alnMode == KINEMA) {		
		double vel=normXY(vnr), wz, yaw, dyaw;
		if(vel>1.0e-3) {
			wz = fabs(sins.wnb.k);
			if(wz<wz0) {  // vel meas
				*(CVect3*)&Zk.dd[0] = sins.vn-vnr;
				SetMeasFlag(0007);
			}
			yaw = atan2(-vnr.i,vnr.j);	dyaw = fabs(diffYaw(yaw,yawPre));
			if(velPre>vel0 && vel>vel0 && wz<wz0 && dyaw<dyaw0) {  // yaw meas
				Zk.dd[3] = -diffYaw(sins.att.k, yaw);			
				if(fabs(Zk.dd[3])>20*DEG) {  // yaw init
					cntYaw = 0;
					SetYaw(yaw); sins.vn = vnr;
					Pset.dd[0] = Pset.dd[1] = pow2(10*DEG); Pset.dd[2] = pow2(20*DEG);
					Pset.dd[3] = Pset.dd[4] = Pset.dd[5] = pow2(10.0);
				}
				else {
					cntYaw++;
					SetMeasFlag(0010);
				}
			}
			velPre = vel;
			yawPre = yaw;
		}
	}
}

int  CMyAlignment::Update(const CVect3 *pwm, const CVect3 *pvm, int nn, double ts, int nSteps) 		
{	
	if (!posInit) return 0;

	if (alnMode == STATIC) {
		alni0.Update(pwm, pvm, nn, ts);
		if (alni0.tk > T1) {
			qnb = alni0.qnb;	alnOK = 1;
		}
	}
	if (alnMode == STATIC_FINE) {
		if (alni0.tk < T1) {
			alni0.Update(pwm, pvm, nn, ts);
		}
		else {
			if (!alnkfInit) {			
				sins = CSINS(alni0.qnb, O31, alni0.pos0);
				alnkfInit = 1;
			}
			TDUpdate(pwm, pvm, nn, ts, nSteps);
			if(sins.mvnk==0 && norm(sins.wnb)<0.1*glv.dps) {
				*(CVect3*)&Zk.dd[0] = sins.mvn;
				SetMeasFlag(0007);				
			}
			if (sins.tk > T2) {
				qnb = sins.qnb;		alnOK = 1;
			}
		}	
	}
	if (alnMode == KINEMA) {
		TDUpdate(pwm, pvm, nn, ts, nSteps);
		if (cntYaw > 5) {
			qnb = sins.qnb;		alnOK = 1;
		}
	}

	return 1;
}

/***************************  class CMySINSGNSS  **********************************/
CMyNavFilter::CMyNavFilter(int nq0, int nr0, double ts, int yawHkRow0):CSINSGNSS(nq0, nr0, ts, yawHkRow0)
{	
	Hk.SetMat3(6, 3, I33);  // 6-8: ZUPT
	Hk(9, 11) = 1.0;		// 9: WzHold (ZIHR)
	Hk(10, 3) = 1.0; 		// 10-11: NHC
	Hk(11, 4) = 1.0; 
}

CMyNavFilter::~CMyNavFilter() 
{

}

void CMyNavFilter::Loadopt(const imuopt_t &opt)                                                  
{
	/* lvGNSS */
	lvGNSS = opt.antlever;

	/* P0 */
	P0_phi = opt.P0_phi;	P0_vel = opt.P0_vel;	P0_pos = opt.P0_pos;
	P0_gyr = opt.P0_gyr;	P0_acc = opt.P0_acc;

	/* Qt */
	Qt_phi = opt.Qt_phi;	Qt_vel = opt.Qt_vel;	Qt_pos = opt.Qt_pos;
	Qt_gyr = opt.Qt_gyr;	Qt_acc = opt.Qt_acc;

	/* Out */
	pos_xyz = O31;			Rk_xyz = I33;
}	

void CMyNavFilter::Init(const CSINS &sins0) 
{
	// 0-14 phi,dvn,dpos,eb,db
	CSINSTDKF::Init(sins0);
	sins.lever(-lvGNSS);  sins.pos = sins.posL; 	
	sins.SetTauGA(CVect3(1800),CVect3(1800)); 				// by zh in 23-08-07
	sins.mmwb = sins.mmfb = CMaxMin(50);					// by zh in 23-08-07
	avpi.Init(sins, kfts, 1);	
	wzhd.Init(200.0*DPH, kfts, 10.0);		

	/* P0 */									
	Pk.SetDiagVect3( 0,  0, pow(P0_phi, 2));	
	Pk.SetDiagVect3( 3,  3, pow(P0_vel, 2));	
	Pk.SetDiagVect3( 6,  6, pow(P0_pos, 2));		
	Pk.SetDiagVect3( 9,  9, pow(P0_gyr, 2));	
	Pk.SetDiagVect3(12, 12, pow(P0_acc, 2));
	/* Qt */
	Qt.Set2Vect3( 0, Qt_phi);
	Qt.Set2Vect3( 3, Qt_vel);
	Qt.Set2Vect3( 6, Qt_pos);
	Qt.Set2Vect3( 9, Qt_gyr);
	Qt.Set2Vect3(12, Qt_acc);	

	// // Pmax.Set2(fDEG3(50.0),  fXXX(100.0),  fdPOS(1.0e6), fDPH3(5000.0),  fMG3(30.0),
	// // 	fXXX(100.0), 0.5, fdKG9(1000.0,900.0), fdKA6(100.0,100.0));
	// Pmin.Set2(fPHI(1,1),  fXXZ(0.01,0.01),  fdLLH(.05,0.1),  fDPH3(0.01),  fXYZU(20,20,500,UG));
		
	// /* Rt */
	// Rt.Set2(fXXZ(0.2,0.2), fdLLH(0.1,0.3), fXXZ(0.05,0.05), 10.0*DPH);
	// // Rmax = Rt*100;  Rmin = Rt;  Rb = 0.6;		
	// // Rt.Set2(fXXZ(0.1,0.2), fdLLH(0.5,1.0), fXXZ(0.05,0.05), 20.0*DPH);
	// FBTau.Set(fXX9(1.0),  fXX6(INF),  fINF3, INF,  fINF9,  fINF6);

	// Pmax.Set2(fDEG3(1.0),  fXXX(100.0),  fdPOS(1.0e6), fDPH3(0.5),  fMG3(1.0),
	// 	fXXX(100.0), 0.5, fdKG9(1000.0,900.0), fdKA6(100.0,100.0));
	// Pmin.Set2(fPHI(0.1,1.0),  fXXX(0.001),  fdLLH(0.005,0.01),  fDPH3(0.001),  fUG3(10.0),
	// 	fXXX(0.001), 0.0001, fdKG9(1.0,1.0), fdKA6(1.0,1.0));
	// Pk.SetDiag2(fPHI(500,600),  fXXX(1.0),  fdPOS(5.0),  fDPH3(0.1),  fMG3(1.0),
	// 	fXXX(1.0),  0.01,  fdKG9(100.0,90.0), fdKA6(10.0,10.0));
	// Qt.Set2(fDPSH3(0.001),  fUGPSHZ3(1.0),  fOOO,  fOO6,
	// 	fOOO, 0.0,  fOO9,  fOO6);
	// Rt.Set2(fXXZ(0.5,1.0),   fdLLH(0.05,0.1), fXXZ(0.005,0.005), 2.0*DPH);
	// FBTau.Set(fXX9(0.5),  fXX6(INF),  fINF3, INF,  fINF9,  fINF6);
	// SetMeasMask(1, 0000070);

	/* INSD */                                                                
	// Pmin.Set2(fPHI(1,1), fXXX(0.001), fdPOS(.01), fDPH3(1.0), fXYZU(250,250,150,UG));
	// Pk.SetDiag2(fPHI(600,600), fXXX(1.0), fdPOS(5.0), fDPH3(1000.0), fMG3(10.0));             
	// Qt.Set2(fDPSH3(1.0), fUGPSHZ3(10.0), fOOO, fOO6, fOOO, 0.0, fOO9, fOO6);
	// Qt.Set2(fDPSH3(1.0), fUGPSHZ3(10.0), fOOO, fOOO, fOOO, fOOO, 0.0, fOO9, fOO6);
	Rt.Set2(fXXZ(0.5,1.0),   fdLLH(0.2,0.5), fXXZ(0.005,0.005), 10.0*DPH); 
	FBTau.Set(fXX9(0.5),  fXX6(1),  fINF3, INF,  fINF9,  fINF6);               
	SetMeasMask(1, 0001770);
	
}

void CMyNavFilter::ZUPTtest(void)										
{
	if(sins.mmwb.flag==1) {
		double maxwb = fabs(sins.mmwb.maxRes), 	minwb = fabs(sins.mmwb.minRes);
		double maxfb = fabs(sins.mmfb.maxRes),  minfb = fabs(sins.mmfb.minRes);
		double dwb  = fabs(maxwb - minwb),		dfb = fabs(maxfb - minfb);
		double mmwb = maxwb > minwb ? maxwb : minwb;
		double mmfb = maxfb > minfb ? maxfb : minfb;
		if(dwb < 0.2*DPS && mmwb < 0.2*DPS && dfb < 20*MG && mmfb < 25*MG && norm(sins.vn) < 0.03) {
			*(CVect3*)&Zk.dd[6] = sins.vn;			
			SetMeasFlag(000700);
			trace(1,"ZUPTtest:%10.4f %10.6f %10.6f %10.6f\n", sins.tk, sins.vn.i, sins.vn.j, sins.vn.k);
		}		
		// trace(1, "%10.2f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n", sins.tk, 
		// 	dwb/glv.dps, dfb/glv.mg, mmwb/glv.dps, mmfb/glv.mg,
		// 	maxwb/glv.dps, minwb/glv.mg, maxfb/glv.dps, minfb/glv.mg);
	}
}

void CMyNavFilter::ZIHRtest(void)								// 参考：https://zhuanlan.zhihu.com/p/115529319
{
	wzhd.Update(sins.imu.phim.k/sins.nts);
	if(wzhd.retres==3 && sins.fn.k>9.5) {						
		Zk.dd[9] = wzhd.meanwpre-sins.eth.wnie.k;			
		SetMeasFlag(001000);
		Xk.dd[2] -= Zk.dd[9]*wzhd.T;									
		trace(2,"ZIHRtest:%10.4f\n", sins.tk);					
	}
}

void CMyNavFilter::NHCtest(void)
{

}

void CMyNavFilter::SetMeasGNSS(const CVect3 &posgnss, const CVect3 &vngnss, double qfactor) 
{
	if(!IsZero(vngnss) && avpi.Interp(vnGNSSdelay+dtGNSSdelay,0x2)) {
		*(CVect3*)&Zk.dd[0] = avpi.vn - vngnss;
		SetMeasFlag(00007);
	}
	if(!IsZero(posgnss) && avpi.Interp(posGNSSdelay+dtGNSSdelay,0x4)) {
		*(CVect3*)&Zk.dd[3] = avpi.pos - posgnss;	
		trace(2,"%10.2f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n", 
			sins.tk, Zk.dd[3]*glv.Re, Zk.dd[4]*glv.Re, Zk.dd[5], 
			sqrt(Rt.dd[3])*glv.Re, sqrt(Rt.dd[4])*glv.Re, sqrt(Rt.dd[5]), 
			sqrt(Pk(6,6))*glv.Re,  sqrt(Pk(7,7))*glv.Re, sqrt(Pk(8,8)),
			sqrt(Pk(3,3)),  sqrt(Pk(4,4)), sqrt(Pk(5,5)));					
		SetMeasFlag(00070);
	}
}

int  CMyNavFilter::Update(const CVect3 *pwm, const CVect3 *pvm, int nn, double ts, int nSteps) 
{
	int res=TDUpdate(pwm, pvm, nn, ts, nSteps);					
	sins.lever(lvGNSS);
	avpi.Push(sins, 1);	
	// /* if(*gnsslost>5) */ ZUPTtest();
	// ZIHRtest();	
	return res;
}

void CMyNavFilter::operator>>(CFileRdWt &f)
{	
	static double dt_1 = 3.0;

	double dt = kftk-(int)(kftk);	

	if (dt < dt_1 && dt_1 < 1.0) {
		avpi.Interp(-dt);
		f << avpi.att << avpi.vn << avpi.pos;					//By zh in 23-08-08
		for (int i=0; i<9; i++) {
			f << Xk(i);
		}
		f << sins.eb.i << sins.eb.j << sins.eb.k;
		f << sins.db.i << sins.db.j << sins.db.k;	
		for (int i=0; i<15; i++) {
			f << sqrt(Pk.dd[i+i*nq]);
		}
		f << (int)(kftk);
	}
	dt_1 = dt;
}

void CMyNavFilter::OutINSResult(double dt) 
{	
	CMat3 T(glv.Re, glv.Re, 1.0), Cen, PP(0.0);

	/* xyz */
	avpi.Interp(dt, 0x4);	
	pos_xyz = blh2xyz(avpi.pos);
	
	/* Rk_neu2xyz */
	Cen = pos2Cen(avpi.pos);
	PP  = T*Pk.GetMat3(6, 6)*T;
	Rk_xyz = Cen*PP*~Cen;	
}

/***************************  class CMyAutoDrive  *********************************/
CMyAutoDrive::CMyAutoDrive()		
{	
	time2gpst(timeget(),&week);	
	ts = 0.01;	 		tk = 0.0;			nn = 1;				tdStep = 10; 
	T1 = 0.0; 			T2 = 0.0;			bIMU2GPS = 0;	

	mState = INITIALIZE;
	mpFile = NULL;		mpGPS  = NULL;		mpAln  = NULL;		mpFIR  = NULL;	
}

CMyAutoDrive::~CMyAutoDrive() 
{	
	delete(mpFile);		delete(mpGPS);		delete(mpAln);		delete(mpFIR);
}

int CMyAutoDrive::Init(const prcopt_t &popt, const solopt_t &sopt, const filopt_t &fopt,
	const sta_t &sta, const imuopt_t &iopt)	
{	
	/* Basic Initialize */	
	time2gpst(popt.ts, &week);	
	ts = 1.0/iopt.freq;		nn = iopt.nn;		tdStep = iopt.tdStep;	

	CFileRdWt::DirO(iopt.outpath);	
	mpFile = new CFileRdWt("kfdebug-urban-f2.bin");

	/* GPS Initialize */
	mpGPS = new CGNSSManage();
	if (!(mpGPS->Init(popt, sopt, fopt, sta))) {
		showmsg("Failed Initialize GNSS!\n");
		return 0;
	}

	/* Aln Initialize */
	mpAln = new CMyAlignment(iopt.aln_mode, iopt.aln_t1, iopt.aln_t2);

	/* FIR Initialize */
	mpFIR = new CMyNavFilter(15, 10, nn*ts);
	mpFIR->Loadopt(iopt);

	/* Else */
	bIMU2GPS = iopt.imu2gps;

	return 1;
}

void CMyAutoDrive::SetMeasGNSS(const CVect3 &pos, const CVect3 &vn, double dt) 		
{	
	if (mState == INITIALIZE) {
		mpAln->Init(CSINS(O31, vn, pos, this->tk));			
		mState = ALIGNMENT;	
	}
	else if (mState == ALIGNMENT) {
		mpAln->SetMeasGNSS(pos, vn);
		if (mpAln->alnOK) {
			mpFIR->Init(CSINS(mpAln->qnb, vn, pos, this->tk));	
			mState = NAVIGATION;
		}
	}
	else if (mState == NAVIGATION) {
		mpFIR->posGNSSdelay = mpFIR->vnGNSSdelay = dt;	
		mpFIR->SetMeasGNSS(pos, vn);	
	}
}

void CMyAutoDrive::GPSProcess(obsd_t *obs, int n) 
{	
	double X[3], P[9];
	/* GNSS Process */
	if (mState == NAVIGATION) {
		mpFIR->OutINSResult(time2gpst(obs[0].time, NULL) - tk);
		X[0] = mpFIR->pos_xyz.i;	X[1] = mpFIR->pos_xyz.j;	X[2] = mpFIR->pos_xyz.k;
		P[0] = mpFIR->Rk_xyz.e00;	P[1] = mpFIR->Rk_xyz.e01;	P[2] = mpFIR->Rk_xyz.e02;
		P[3] = mpFIR->Rk_xyz.e10;	P[4] = mpFIR->Rk_xyz.e11;	P[5] = mpFIR->Rk_xyz.e12;
		P[6] = mpFIR->Rk_xyz.e20;	P[7] = mpFIR->Rk_xyz.e21;	P[8] = mpFIR->Rk_xyz.e22;	
		mpGPS->Update(obs, n, bIMU2GPS, X, P);

		// mpGPS->Update(obs[0].time);			// By zh in 23-08-08			
	}
	else {
		mpGPS->Update(obs, n);
	}
	tgps = time2gpst(*mpGPS->ti, NULL);

	if ((T1 < tgps) && (tgps < T2)) return;	
	if (*mpGPS->stat == SOLQ_NONE) return;

	mpGPS->OutGPSResult();
	mpFIR->Rt.Set2Vect3(3, mpGPS->rk_blh);
	SetMeasGNSS(mpGPS->pos_blh, mpGPS->vel_enu, tgps - tk);	
}

void CMyAutoDrive::IMUProcess(imud_t *pimu)		 	
{	
	CVect3 wm[10], vm[10];
	for (int i = 0; i < nn; i++) {
		wm[i] = pimu[i].wm; 
		vm[i] = pimu[i].vm;
		if (tk != 0.0 && (pimu[i].t-tk) > 0.011) {
			trace(1,"%10.4f----%10.4f\n",pimu[i].t, tk);
		}
		tk = pimu[i].t;
	}

	if (mState == ALIGNMENT) {
		mpAln->Update(wm, vm, nn, ts, tdStep);	
	}
	if (mState == NAVIGATION) {
		mpFIR->Update(wm, vm, nn, ts, tdStep);
		*mpFIR >> *mpFile;
	}

	printbrk(gpst2time(week,tk), *mpGPS->stat, mState);	
}

void CMyAutoDrive::TraceOpen(const char *filename) 
{
	mpGPS->TraceOpen(filename);
}

void CMyAutoDrive::Sol2kml() 
{
	mpGPS->Sol2kml();
}

/***************************  Extern Function  *********************************/
int getSamples(CMyAutoDrive &app) 
{
	return app.nn;
}
void SetInterrupt(CMyAutoDrive &app, double T1, double T2)
{
	app.T1 = T1;		app.T2 = T2;
}
void printbrk(gtime_t ti, int gnss_stat, int gins_stat) 
{	
	char str[32] = "processing : %s GINS=%d GNSS=%d";
	checkbrk(str, time_str(ti,2), gins_stat, gnss_stat);	
}