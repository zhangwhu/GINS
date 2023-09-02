
#include "CGNSSConfig.h"
#include "CIMUConfig.h"
#include "CKFApp.h"

/* constants/global variables ------------------------------------------------*/
static obs_t  obss = { 0 };	/* observation data */
static imu_t  imus = { 0 };   /* imu data */
static sta_t   sta = { 0 };
static int    iobs = 0;		/* current rover observation data index */
static int    iimu = 0;		/* current imu data index */
static int  nepoch;

typedef struct {
	double t;
	CVect3 vngps, posgps;
} gps_t;

static int data_align(CRMemory &memimu, CRMemory &memobs, int &iimu, int &iobsu) 
{
	double tow, refs;
	imud_t *pimu = (imud_t *)memimu.get(0);
	obsd_t *pobs = (obsd_t *)memobs.get(0);

	refs = (pimu->t > time2gpst(pobs->time, NULL)) ? pimu->t : time2gpst(pobs->time, NULL);

	for (iobsu = 0; iobsu<memobs.recordNum; iobsu++) {
		tow = time2gpst(pobs[iobsu].time, NULL);
		if (tow >= refs) break;
	}
	for ( iimu = 0;  iimu<memimu.recordNum; iimu++) {
		tow = pimu[iimu].t;
		if (tow >= refs) break;
	}

	return (iimu<memimu.recordNum && iobsu<memobs.recordNum) ? 1 : 0;
}

static int inputobs(CRMemory &memobs, int &iobsu, obsd_t *obs) 
{
	int n;
	obsd_t *pobs = (obsd_t *)memobs.get(iobsu);

	for (n = 0; n<memobs.recordNum-iobsu; n++) {
		if (timediff(pobs[n].time, pobs[0].time)>DTTOL) {
			break;
		}
		obs[n] = pobs[n];
	}
	iobsu += n;

	return n;
}

static void Forward(CRMemory &memimu, CRMemory &memobs, CMyAutoDrive &app) 
{
	imud_t *pimu;
	obsd_t obs[MAXOBS];
	int    nobs = inputobs(memobs, iobs, obs), nn = getSamples(app);			
	double tobs = time2gpst(obs[0].time, NULL);
	for (; iimu + nn < imus.n; iimu += nn) 
	{
		pimu = (imud_t *)memimu.get(iimu);	
		app.IMUProcess(pimu);
		if (app.tk >= tobs) {		
			app.GPSProcess(obs, nobs);											
			nobs = inputobs(memobs, iobs, obs);			
			tobs = time2gpst(obs[0].time, NULL);
		}
	}
}	

int main(int argc, char *argv[])	
{
	unsigned int tick = tickget();

	/* Read YAML */
	CGNSSConfig Gyaml;				
	CIMUConfig  Iyaml;

	if (!Gyaml.LoadYAML(argv[1]) || !Iyaml.LoadYAML(argv[1])) return -1;
	
	/* Read Data */
	if (!readobs(Gyaml.filopt.rovobs, 1, &Gyaml.prcopt, &obss, &sta, &nepoch)) return 0;
	if (!readimu(&Iyaml.imuopt, &imus)) return 0;		

	CRMemory memimu((BYTE*)imus.data, imus.n*sizeof(imud_t), sizeof(imud_t));							
	CRMemory memobs((BYTE*)obss.data, obss.n*sizeof(obsd_t), sizeof(obsd_t));
	if (!data_align(memimu, memobs, iimu, iobs)) return 0;
	
	/* Processing instantiation */		
	CMyAutoDrive app;		
	if (!(app.Init(Gyaml.prcopt, Gyaml.solopt, Gyaml.filopt, sta, Iyaml.imuopt))) return 0;		

	app.TraceOpen(Gyaml.filopt.trace);
	
	// SetInterrupt(app, 371287, 371327);
	/* Forward Processing */
	Forward(memimu, memobs, app);	  	

	app.Sol2kml();
	
	checkbrk("%40s","");
	printf("Time=%.1f\n", (tickget() - tick) * 0.001);   
	return 1;
}