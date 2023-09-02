/*------------------------------------------------------------------------------
* ionex.c : ionex functions������������Ʒ��
*
*          Copyright (C) 2011-2013 by T.TAKASU, All rights reserved.
*
* references:
*     [1] S.Schear, W.Gurtner and J.Feltens, IONEX: The IONosphere Map EXchange
*         Format Version 1, February 25, 1998
*     [2] S.Schaer, R.Markus, B.Gerhard and A.S.Timon, Daily Global Ionosphere
*         Maps based on GPS Carrier Phase Data Routinely producted by CODE
*         Analysis Center, Proceeding of the IGS Analysis Center Workshop, 1996
*
* version : $Revision:$ $Date:$
* history : 2011/03/29 1.0 new
*           2013/03/05 1.1 change api readtec()
*                          fix problem in case of lat>85deg or lat<-85deg
*           2014/02/22 1.2 fix problem on compiled as C++
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define SQR(x)      ((x)*(x))
#define VAR_NOTEC   SQR(30.0)   /* variance of no tec */
#define MIN_EL      0.0         /* min elevation angle (rad) */
#define MIN_HGT     -1000.0     /* min user height (m) */

/* get index -----------------------------------------------------------------*/
static int getindex(double value, const double *range)
{
    if (range[2]==0.0) return 0;
    if (range[1]>0.0&&(value<range[0]||range[1]<value)) return -1;
    if (range[1]<0.0&&(value<range[1]||range[0]<value)) return -1;
    return (int)floor((value-range[0])/range[2]+0.5);
}
/* get number of items -------------------------------------------------------*/
static int nitem(const double *range)
{
    return getindex(range[1],range)+1;
}
/* data index (i:lat,j:lon,k:hgt) --------------------------------------------*/
static int dataindex(int i, int j, int k, const int *ndata)
{
    if (i<0||ndata[0]<=i||j<0||ndata[1]<=j||k<0||ndata[2]<=k) return -1;
    return i+ndata[0]*(j+ndata[1]*k);
}
/* add tec data to navigation data -------------------------------------------*/
static tec_t *addtec(const double *lats, const double *lons, const double *hgts,
                     double rb, nav_t *nav)
{
    tec_t *p,*nav_tec;
    gtime_t time0={0};
    int i,n,ndata[3];
    
    trace(3,"addtec  :\n");
    
    ndata[0]=nitem(lats);       
    ndata[1]=nitem(lons);
    ndata[2]=nitem(hgts);
    if (ndata[0]<=1||ndata[1]<=1||ndata[2]<=0) return NULL;
    
    if (nav->nt>=nav->ntmax) {
        nav->ntmax+=256;
        if (!(nav_tec=(tec_t *)realloc(nav->tec,sizeof(tec_t)*nav->ntmax))) {
            trace(1,"readionex malloc error ntmax=%d\n",nav->ntmax);
            free(nav->tec); nav->tec=NULL; nav->nt=nav->ntmax=0;
            return NULL;
        }
        nav->tec=nav_tec;
    }
    p=nav->tec+nav->nt;
    p->time=time0;
    p->rb=rb;
    for (i=0;i<3;i++) {
        p->ndata[i]=ndata[i];
        p->lats[i]=lats[i];
        p->lons[i]=lons[i];
        p->hgts[i]=hgts[i];
    }
    n=ndata[0]*ndata[1]*ndata[2];
    
    if (!(p->data=(double *)malloc(sizeof(double)*n))||
        !(p->rms =(float  *)malloc(sizeof(float )*n))) {
        return NULL;
    }
    for (i=0;i<n;i++) {
        p->data[i]=0.0;
        p->rms [i]=0.0f;
    }
    nav->nt++;
    return p;
}
/* read ionex dcb aux data ----------------------------------------------------*/
static void readionexdcb(FILE *fp, double *dcb, double *rms)
{
    int i,sat;
    char buff[1024],id[32],*label;
    
    trace(3,"readionexdcb:\n");
    
    for (i=0;i<MAXSAT;i++) dcb[i]=rms[i]=0.0;
    
    while (fgets(buff,sizeof(buff),fp)) {
        if (strlen(buff)<60) continue;
        label=buff+60;
        
        if (strstr(label,"PRN / BIAS / RMS")==label) {			
            
            strncpy(id,buff+3,3); id[3]='\0';
            
            if (!(sat=satid2no(id))) {
                trace(2,"ionex invalid satellite: %s\n",id);
                continue;
            }
            dcb[sat-1]=str2num(buff, 6,10);
            rms[sat-1]=str2num(buff,16,10);
        }
        else if (strstr(label,"END OF AUX DATA")==label) break;
    }
}
/* read ionex header ---------------------------------------------------------*/
static double readionexh(FILE *fp, double *lats, double *lons, double *hgts,
                         double *rb, double *nexp, double *dcb, double *rms)
{
    double ver=0.0;
    char buff[1024],*label;
    
    trace(3,"readionexh:\n");
    
    while (fgets(buff,sizeof(buff),fp)) {
        
        if (strlen(buff)<60) continue;
        label=buff+60;
        
        if (strstr(label,"IONEX VERSION / TYPE")==label) {
            if (buff[20]=='I') ver=str2num(buff,0,8);
        }
        else if (strstr(label,"BASE RADIUS")==label) {
            *rb=str2num(buff,0,8);
        }
        else if (strstr(label,"HGT1 / HGT2 / DHGT")==label) {
            hgts[0]=str2num(buff, 2,6);
            hgts[1]=str2num(buff, 8,6);
            hgts[2]=str2num(buff,14,6);
        }
        else if (strstr(label,"LAT1 / LAT2 / DLAT")==label) {
            lats[0]=str2num(buff, 2,6);
            lats[1]=str2num(buff, 8,6);
            lats[2]=str2num(buff,14,6);
        }
        else if (strstr(label,"LON1 / LON2 / DLON")==label) {
            lons[0]=str2num(buff, 2,6);
            lons[1]=str2num(buff, 8,6);
            lons[2]=str2num(buff,14,6);
        }
        else if (strstr(label,"EXPONENT")==label) {          //ָ��
            *nexp=str2num(buff,0,6);
        }
        else if (strstr(label,"START OF AUX DATA")==label&&
                 strstr(buff,"DIFFERENTIAL CODE BIASES")) {
            readionexdcb(fp,dcb,rms);
        }
        else if (strstr(label,"END OF HEADER")==label) {
            return ver;
        }
    }
    return 0.0;
}
/* read ionex body -----------------------------------------------------------*/
static int readionexb(FILE *fp, const double *lats, const double *lons,
                      const double *hgts, double rb, double nexp, nav_t *nav)
{
    tec_t *p=NULL;
    gtime_t time={0};
    double lat,lon[3],hgt,x;
    int i,j,k,n,m,index,type=0;
    char buff[1024],*label=buff+60;
    
    trace(3,"readionexb:\n");
    
    while (fgets(buff,sizeof(buff),fp)) {
        
        if (strlen(buff)<60) continue;
        
        if (strstr(label,"START OF TEC MAP")==label) {
            if ((p=addtec(lats,lons,hgts,rb,nav))) type=1;
        }
        else if (strstr(label,"END OF TEC MAP")==label) {
            type=0;
            p=NULL;
        }
        else if (strstr(label,"START OF RMS MAP")==label) {
            type=2;
            p=NULL;
        }
        else if (strstr(label,"END OF RMS MAP")==label) {
            type=0;
            p=NULL;
        }
        else if (strstr(label,"EPOCH OF CURRENT MAP")==label) {
            if (str2time(buff,0,36,&time)) {
                trace(2,"ionex epoch invalid: %-36.36s\n",buff);
                continue;
            }
            if (type==2) {
                for (i=nav->nt-1;i>=0;i--) {
                    if (fabs(timediff(time,nav->tec[i].time))>=1.0) continue;
                    p=nav->tec+i;
                    break;
                }
            }
            else if (p) p->time=time;
        }
        else if (strstr(label,"LAT/LON1/LON2/DLON/H")==label&&p) {
            lat   =str2num(buff, 2,6);
            lon[0]=str2num(buff, 8,6);
            lon[1]=str2num(buff,14,6);
            lon[2]=str2num(buff,20,6);
            hgt   =str2num(buff,26,6);
            
            i=getindex(lat,p->lats);
            k=getindex(hgt,p->hgts);
            n=nitem(lon);
            
            for (m=0;m<n;m++) {
                if (m%16==0&&!fgets(buff,sizeof(buff),fp)) break;
                
                j=getindex(lon[0]+lon[2]*m,p->lons);
                if ((index=dataindex(i,j,k,p->ndata))<0) continue;
                
                if ((x=str2num(buff,m%16*5,5))==9999.0) continue;
                
                if (type==1) p->data[index]=x*pow(10.0,nexp);				//pow(10.nexp):����λת��ΪTECU
                else p->rms[index]=(float)(x*pow(10.0,nexp));
            }
        }
    }
    return 1;
}
/* combine tec grid data -----------------------------------------------------*/
static void combtec(nav_t *nav)
{
    tec_t tmp;
    int i,j,n=0;
    
    trace(3,"combtec : nav->nt=%d\n",nav->nt);
    
    for (i=0;i<nav->nt-1;i++) {
        for (j=i+1;j<nav->nt;j++) {
            if (timediff(nav->tec[j].time,nav->tec[i].time)<0.0) {      //����
                tmp=nav->tec[i];
                nav->tec[i]=nav->tec[j];
                nav->tec[j]=tmp;
            }
        }
    }
    for (i=0;i<nav->nt;i++) {
        if (i>0&&timediff(nav->tec[i].time,nav->tec[n-1].time)==0.0) {   //�ų�����
            free(nav->tec[n-1].data);
            free(nav->tec[n-1].rms );
            nav->tec[n-1]=nav->tec[i];
            continue;
        }
        nav->tec[n++]=nav->tec[i];
    }
    nav->nt=n;
    
    trace(4,"combtec : nav->nt=%d\n",nav->nt);
}
/* read ionex tec grid file ----------------------------------------------------
* read ionex ionospheric tec grid file
* args   : char   *file       I   ionex tec grid file
*                                 (wind-card * is expanded)
*          nav_t  *nav        IO  navigation data
*                                 nav->nt, nav->ntmax and nav->tec are modified
*          int    opt         I   read option (1: no clear of tec data,0:clear)
* return : none
* notes  : see ref [1]
*-----------------------------------------------------------------------------*/
extern void readtec(const char *file, nav_t *nav, int opt)      
{
    FILE *fp;
    double lats[3]={0},lons[3]={0},hgts[3]={0},rb=0.0,nexp=-1.0;
    double dcb[MAXSAT]={0},rms[MAXSAT]={0};
    int i,n;
    char *efiles[MAXEXFILE];
    
    trace(3,"readtec : file=%s\n",file);
    
    /* clear of tec grid data option */
    if (!opt) {
        free(nav->tec); nav->tec=NULL; nav->nt=nav->ntmax=0;
    }
    for (i=0;i<MAXEXFILE;i++) {
        if (!(efiles[i]=(char *)malloc(1024))) {
            for (i--;i>=0;i--) free(efiles[i]);
            return;
        }
    }
    /* expand wild card in file path */
    n=expath(file,efiles,MAXEXFILE);
    
    for (i=0;i<n;i++) {
        if (!(fp=fopen(efiles[i],"r"))) {
            trace(2,"ionex file open error %s\n",efiles[i]);
            continue;
        }
        /* read ionex header */
        if (readionexh(fp,lats,lons,hgts,&rb,&nexp,dcb,rms)<=0.0) {
            trace(2,"ionex file format error %s\n",efiles[i]);
            continue;
        }
        /* read ionex body */
		readionexb(fp, lats, lons, hgts, rb, nexp, nav);
        
        fclose(fp);
    }
    for (i=0;i<MAXEXFILE;i++) free(efiles[i]);
    
    /* combine tec grid data */
    if (nav->nt>0) combtec(nav);			
    
    /* P1-P2 dcb */
    for (i=0;i<MAXSAT;i++) {
        nav->cbias[i][0]=CLIGHT*dcb[i]*1E-9; /* ns->m */   
    }
}
/* interpolate tec grid data -------------------------------------------------*/   //��ֵ��������
static int interptec(const tec_t *tec, int k, const double *posp, double *value,
                     double *rms)
{
    double dlat,dlon,a,b,d[4]={0},r[4]={0};
    int i,j,n,index;
    
    trace(3,"interptec: k=%d posp=%.2f %.2f\n",k,posp[0]*R2D,posp[1]*R2D);
    *value=*rms=0.0;
    
    if (tec->lats[2]==0.0||tec->lons[2]==0.0) return 0;
    
    dlat=posp[0]*R2D-tec->lats[0];
    dlon=posp[1]*R2D-tec->lons[0];
    if (tec->lons[2]>0.0) dlon-=floor( dlon/360)*360.0; /*  0<=dlon<360 */
    else                  dlon+=floor(-dlon/360)*360.0; /* -360<dlon<=0 */
    
    a=dlat/tec->lats[2];
    b=dlon/tec->lons[2];
    i=(int)floor(a); a-=i;
    j=(int)floor(b); b-=j;
    
    /* get gridded tec data */
    for (n=0;n<4;n++) {
        if ((index=dataindex(i+(n%2),j+(n<2?0:1),k,tec->ndata))<0) continue;				
        d[n]=tec->data[index];
        r[n]=tec->rms [index];
    }
    if (d[0]>0.0&&d[1]>0.0&&d[2]>0.0&&d[3]>0.0) {
        
        /* bilinear interpolation (inside of grid) */
        *value=(1.0-a)*(1.0-b)*d[0]+a*(1.0-b)*d[1]+(1.0-a)*b*d[2]+a*b*d[3];
        *rms  =(1.0-a)*(1.0-b)*r[0]+a*(1.0-b)*r[1]+(1.0-a)*b*r[2]+a*b*r[3];
    }
    /* nearest-neighbour extrapolation (outside of grid) */
    else if (a<=0.5&&b<=0.5&&d[0]>0.0) {*value=d[0]; *rms=r[0];}
    else if (a> 0.5&&b<=0.5&&d[1]>0.0) {*value=d[1]; *rms=r[1];}
    else if (a<=0.5&&b> 0.5&&d[2]>0.0) {*value=d[2]; *rms=r[2];}
    else if (a> 0.5&&b> 0.5&&d[3]>0.0) {*value=d[3]; *rms=r[3];}
    else {
        i=0;
        for (n=0;n<4;n++) if (d[n]>0.0) {i++; *value+=d[n]; *rms+=r[n];}
        if(i==0) return 0;
        *value/=i; *rms/=i;
    }
    return 1;
}
/* ionosphere delay by tec grid data -----------------------------------------*/
static int iondelay(gtime_t time, const tec_t *tec, const double *pos,
                    const double *azel, int opt, double *delay, double *var)
{
    const double fact=40.30E16/FREQ1/FREQ1; /* tecu->L1 iono (m) */
    double fs,posp[3]={0},vtec,rms,hion,rp;
    int i;
    
    trace(3,"iondelay: time=%s pos=%.1f %.1f azel=%.1f %.1f\n",time_str(time,0),
          pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);
    
    *delay=*var=0.0;
    
    for (i=0;i<tec->ndata[2];i++) { /* for a layer */					//һ���ǵ���ģ��
        
        hion=tec->hgts[0]+tec->hgts[2]*i;
        
        /* ionospheric pierce point position */
        fs=ionppp(pos,azel,tec->rb,hion,posp);
        
        if (opt&2) {
            /* modified single layer mapping function (M-SLM) ref [2] */
			//hion = 506.7;												//2021-06-5����
            rp=tec->rb/(tec->rb+hion)*sin(0.9782*(PI/2.0-azel[1]));
            fs=1.0/sqrt(1.0-rp*rp);
        }
        if (opt&1) {
            /* earth rotation correction (sun-fixed coordinate) */
            posp[1]+=2.0*PI*timediff(time,tec->time)/86400.0;
        }
        /* interpolate tec grid data */
        if (!interptec(tec,i,posp,&vtec,&rms)) return 0;				//���ݴ��̵����vtec
        
        *delay+=fact*fs*vtec;											//vtec->tecu->m
        *var+=fact*fact*fs*fs*rms*rms;
    }
	//trace(1, "%s %6.3f %10.4f %10.4f\n", time_str(time, 0), azel[1] / PI * 180, fs, *delay);
    trace(4,"iondelay: delay=%7.2f std=%6.2f\n",*delay,sqrt(*var));
    
    return 1;
}
/* ionosphere model by tec grid data -------------------------------------------
* compute ionospheric delay by tec grid data
* args   : gtime_t time     I   time (gpst)
*          nav_t  *nav      I   navigation data
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          int    opt       I   model option
*                                bit0: 0:earth-fixed,1:sun-fixed
*                                bit1: 0:single-layer,1:modified single-layer
*          double *delay    O   ionospheric delay (L1) (m)
*          double *var      O   ionospheric dealy (L1) variance (m^2)
* return : status (1:ok,0:error)
* notes  : before calling the function, read tec grid data by calling readtec()
*          return ok with delay=0 and var=VAR_NOTEC if el<MIN_EL or h<MIN_HGT
*-----------------------------------------------------------------------------*/
extern int iontec(gtime_t time, const nav_t *nav, const double *pos,
                  const double *azel, int opt, double *delay, double *var)
{
    double dels[2],vars[2],a,tt;
    int i,stat[2];
    
    trace(3,"iontec  : time=%s pos=%.1f %.1f azel=%.1f %.1f\n",time_str(time,0),
          pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);
    
    if (azel[1]<MIN_EL||pos[2]<MIN_HGT) {
        *delay=0.0;
        *var=VAR_NOTEC;
        return 1;
    }
    for (i=0;i<nav->nt;i++) {
        if (timediff(nav->tec[i].time,time)>0.0) break;
    }
    if (i==0||i>=nav->nt) {
        trace(2,"%s: tec grid out of period\n",time_str(time,0));
        return 0;
    }
    if ((tt=timediff(nav->tec[i].time,nav->tec[i-1].time))==0.0) {				
        trace(2,"tec grid time interval error\n");
        return 0;
    }
    /* ionospheric delay by tec grid data */
    stat[0]=iondelay(time,nav->tec+i-1,pos,azel,opt,dels  ,vars  );       
    stat[1]=iondelay(time,nav->tec+i  ,pos,azel,opt,dels+1,vars+1);
    
    if (!stat[0]&&!stat[1]) {
        trace(2,"%s: tec grid out of area pos=%6.2f %7.2f azel=%6.1f %5.1f\n",
              time_str(time,0),pos[0]*R2D,pos[1]*R2D,azel[0]*R2D,azel[1]*R2D);
        return 0;
    }
    if (stat[0]&&stat[1]) { /* linear interpolation by time */
        a=timediff(time,nav->tec[i-1].time)/tt;
        *delay=dels[0]*(1.0-a)+dels[1]*a;
        *var  =vars[0]*(1.0-a)+vars[1]*a;
    }
    else if (stat[0]) { /* nearest-neighbour extrapolation by time */
        *delay=dels[0];
        *var  =vars[0];
    }
    else {
        *delay=dels[1];
        *var  =vars[1];
    }
    trace(3,"iontec  : delay=%5.2f std=%5.2f\n",*delay,sqrt(*var));
    return 1;
}

/* get ionospherec correction -------------------------------------------------- */
extern int pppcorr_stat(gtime_t time, const double *blh, const int sat, const nav_t *nav)			
{
    uint8_t flag;
	int i, k, is, im, ie, n;
	double dt1, dt2, c1, c2; 
    const corrstec_t *corr;    

	for (k = n = 0; k < nav->nstec; k++) {                     
        corr = nav->corrstec + k;
        for (i = is = ie = 0, flag = 0; i < corr->nn - 1; i++) {
            dt1 = timediff(time, corr->time[i]);	
            dt2 = timediff(time, corr->time[i+1]);
		    if (dt1 > DTTOL && dt2 < DTTOL) {
			    is = corr->ii[i];       im = corr->ii[i+1];
			    ie = (i+2 < corr->nn) ? corr->ii[i+2] : corr->ns; 
                c1 = fabs(dt2) / (fabs(dt1)+fabs(dt2));
                c2 = fabs(dt1) / (fabs(dt1)+fabs(dt2));	
                flag = 1;    
                break;
		    }
	    }
        if (!flag) continue;
        
	    for (i = is, flag = 0; i < ie; i++) {
		    if (corr->data[i].sat == sat) { 
                flag += (i < im) ? 1 : 4; 
            }
	    }
        if ((flag == 5) || (flag == 1 && c1 > 0.5) || (flag == 4 && c2 > 0.5)) {
            if(++n >= 3) return 1;
        } 
    } 
	
	return 0;
}

extern int pppcorr_stec(gtime_t time, const double *blh, const int sat1, const int sat2, 
	double *stec, double *var, const nav_t *nav, int opt)			
{
    uint8_t flag;
	int i, k, is, im, ie, n, naa, info;
	double dt1, dt2, c1, c2, ion[2*10], std[2*10], mid, D[10]; 
    double A[3*10], L[10], Naa[3*3], W[3] = { -1.0, 0.0, 0.0 }, K[10];
    const corrstec_t *corr;    

	for (k = 0, n = 0; k < nav->nstec; k++) {
        corr = nav->corrstec + k;
        ion[0+2*n] = ion[1+2*n] = std[0+2*n] = std[1+2*n] = 0.0;
                
        for (i = is = ie = 0, flag = 0; i < corr->nn - 1; i++) {
            dt1 = timediff(time, corr->time[i]);	
            dt2 = timediff(time, corr->time[i+1]);
            
		    if (dt1 > DTTOL && dt2 < DTTOL) {
			    is = corr->ii[i];       im = corr->ii[i+1];
			    ie = (i+2 < corr->nn) ? corr->ii[i+2] : corr->ns;    
                c1 = fabs(dt2) / (fabs(dt1)+fabs(dt2));
                c2 = fabs(dt1) / (fabs(dt1)+fabs(dt2));	 
                flag = 1;
                break;
		    }
	    }
        if (!flag) continue;
        
	    for (i = is, flag = 0; i < ie; i++) {
		    if (corr->data[i].sat == sat1) { 
                ion[0+2*n] += (i < im) ? c1*corr->data[i].ion : c2*corr->data[i].ion;
                std[0+2*n] += (i < im) ? c1*corr->data[i].std : c2*corr->data[i].std; 
                flag += (i < im) ? 1 : 4; 
            }
		    if (corr->data[i].sat == sat2) { 
                ion[1+2*n] += (i < im) ? c1*corr->data[i].ion : c2*corr->data[i].ion;  
                std[1+2*n] += (i < im) ? c1*corr->data[i].std : c2*corr->data[i].std; 
                flag += (i < im) ? 2 : 8; 
            }
	    }
        if ((flag == 15) || (flag == 3 && c1 > 0.5) || (flag == 12 && c2 > 0.5)) {
            ion[0+2*n] /= (flag == 3) ? c1 : (flag == 12) ? c2 : 1;
            ion[1+2*n] /= (flag == 3) ? c1 : (flag == 12) ? c2 : 1;
            std[0+2*n] /= (flag == 3) ? c1 : (flag == 12) ? c2 : 1;
            std[1+2*n] /= (flag == 3) ? c1 : (flag == 12) ? c2 : 1;
            A[0+3*n] = 1.0;       
            A[1+3*n] = (corr->rr[0] - blh[0])*RE_WGS84;     
            A[2+3*n] = (corr->rr[1] - blh[1])*RE_WGS84;
            D[n++]   = (ion[0+2*n] - ion[1+2*n]);

            char T[32], id1[8], id2[8];
            time2str(time, T, 0);
            satno2id(sat1, id1);
			satno2id(sat2, id2);
            trace(opt,"%s %.0f %s %s %s %10.4f %10.4f %10.4f\n", T, time2gpst(time,NULL),corr->stas, id1, id2, 
                ion[0+2*(n-1)], ion[1+2*(n-1)], ion[0+2*(n-1)]-ion[1+2*(n-1)]);
        } 
    }

    mid = median(D, n);
    for (i = 0, naa = 0; i < n; i++) {
        // if (fabs(D[i]-mid)<=0.6) {
        if (fabs(D[i]-mid)<=0.28) {          //小网原始值
        // if (fabs(D[i]-mid)<=1.5) {          //大网时需要改
            matcpy(A+3*naa, A+3*i, 1, 3);
            matcpy(ion+2*naa, ion+2*i, 1, 3);
            matcpy(std+2*naa, std+2*i, 1, 3);
            naa++;
        }
    }  

    if (naa >= 3) {
        matmul("NT", 3, 3, naa, 1.0, A, A, 0.0, Naa);
        if (!(info = matinv(Naa, 3))) {
            matmul("NN", 3, 1, 3, -1.0, Naa, W, 0.0, K);
            matmul("TN", naa, 1, 3,  1.0,   A, K, 0.0, L);
            for (i = *stec = *var = 0; i < naa; i++) {   
                if (fabs(L[i])>1) return 0;        
                *stec += L[i]*(ion[0+2*i] - ion[1+2*i]);
                *var  += SQR(L[i]*std[0+2*i]) + SQR(L[i]*std[1+2*i]);
                trace(opt, "%4d %.4f %.4f\n", i + 1, L[i], ion[0+2*i] - ion[1+2*i]); 
            }
            return 1;
        }
    }
	
	return 0;
}

