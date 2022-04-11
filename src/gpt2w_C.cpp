//Converted to C by Haroldo A. Marques (haroldoh2o@gmail.com)
//Date: 28/06/2017

#include "gpt2w_C.h"

GPT2_Data Data[64800];

  const double pi = 3.1415926535;
//---------------------------------------------------------------------------
int main11()    // program sample_mainprogram
{
/*
   This programme serves to better understand GPT2w and the modeling of total delays
   Johannes Boehm, 30 August 2013
   Johannes Boehm, 24 December 2014, converted to Fortran
*/

      double cdlat[10]={0.0},cdlon[10]={0.0},chell[10]={0.0},
             cp[10]={0.0},cT[10]={0.0},cdT[10]={0.0},
             cTm[10]={0.0},ce[10]={0.0},cah[10]={0.0},caw[10]={0.0},
             clambda[10]={0.0},cundu[10]={0.0};

      double zhd=0.0, zwd=0.0, dmfh=0.0, dmfw=0.0;

      //Example input data
      double cdmjd = 56141.0;
      cdlat[0] = 48.20*pi/180.0;
      cdlon[0] = 16.37*pi/180.0;
      chell[0] = 156.0;

      //elevation angle shall be 10 degree
      double zd = (90.0 - 10.0)*pi/180.0; //zenith distance in radians

      //the following two parameters are GPT2-specific; just fix them to:
      int nstat = 1;   //we do the calculations for nstat stations
      int it = 0;       //0 means GPT2 with time variation, 1 means static

      //We call GPT2
	  gpt2_1w(cdmjd, cdlat, cdlon, chell, nstat, it, cp, cT, cdT, cTm, ce, cah, caw, clambda, cundu);

      printf("%lf %lf %lf %lf %lf %lf %lf\n",cp[0],cT[0],ce[0],cTm[0],clambda[0],cah[0],caw[0]);

      //The zenith hydrostatic delay is always calculated with the equation by
      //Saastamoinen (1972) as refined by Davis et al. (1985)
	  saasthyd(cp[0], cdlat[0], chell[0], zhd);

      //The zenith wet delay can be calculated with (2) the equation 18 by Askne and Nordius (1987)
	  asknewet(ce[0], cTm[0], clambda[0], zwd);

      printf("%lf %lf\n",zhd,zwd);

      //Calculate the mapping functions
	  vmf1_ht(cah[0], caw[0], cdmjd, cdlat[0], chell[0], zd, dmfh, dmfw);

	  printf("%lf %lf\n", dmfh, dmfw);

      //total delay = 13.698 m in this example
      double delay = zhd*dmfh + zwd*dmfw;

      printf("%lf", delay);

      getchar();

    return 0;
}
//---------------------------------------------------------------------------
double sign(double x, double y)
{  //similar function 'sign' in fortran

    if (y >= 0)
        return fabs(x);

    if (y<0) return (-1.0*fabs(x));
}
//---------------------------------------------------------------------------
void gpt2_1w (double dmjd,double dlat[],double dlon[],double hell[],int nstat,int it,
                   double p[],double T[],double dT[],double Tm[],double e[],
                   double ah[],double aw[],double la[],double undu[])
{
/*--------------------------------------------------------------------
 (c) Department of Geodesy and Geoinformation, Vienna University of
 Technology, 2013

 The copyright in this document is vested in the Department of Geodesy and
 Geoinformation (GEO), Vienna University of Technology, Austria. This document
 may only be reproduced in whole or in part, or stored in a retrieval
 system, or transmitted in any form, or by any means electronic,
 mechanical, photocopying or otherwise, either with the prior permission
 of GEO or in accordance with the terms of ESTEC Contract No.
 4000107329/12/NL/LvH.
 ---

 This subroutine determines pressure, temperature, temperature lapse rate,
 mean temperature of the water vapor, water vapor pressure, hydrostatic
 and wet mapping function coefficients ah and aw, water vapour decrease
 factor and geoid undulation for specific sites near the Earth surface.
 It is based on a 1 x 1 degree external grid file ('gpt2_1wA.grd') with mean
 values as well as sine and cosine amplitudes for the annual and
 semiannual variation of the coefficients.

 c Reference:
 J. Böhm, G. Möller, M. Schindelegger, G. Pain, R. Weber, Development of an
 improved blind model for slant delays in the troposphere (GPT2w),
 GPS Solutions, 2015, doi:10.1007/s10291-014-0403-7

 input parameters:

 dmjd:  modified Julian date (scalar, only one epoch per call is possible)
 dlat:  ellipsoidal latitude in radians [-pi/2:+pi/2] (vector)
 dlon:  longitude in radians [-pi:pi] or [0:2pi] (vector)
 hell:  ellipsoidal height in m (vector)
 nstat: number of stations in dlat, dlon, and hell
        maximum possible: not relevant for Matlab version
 it:    case 1: no time variation but static quantities
        case 0: with time variation (annual and semiannual terms)

 output parameters:

 p:    pressure in hPa (vector of length nstat)
 T:    temperature in degrees Celsius (vector of length nstat)
 dT:   temperature lapse rate in degrees per km (vector of length nstat)
 Tm:   mean temperature of the water vapor in degrees Kelvin (vector of length nstat)
 e:    water vapor pressure in hPa (vector of length nstat)
 ah:   hydrostatic mapping function coefficient at zero height (VMF1)
       (vector of length nstat)
 aw:   wet mapping function coefficient (VMF1) (vector of length nstat)
 la:   water vapor decrease factor (vector of length nstat)
 undu: geoid undulation in m (vector of length nstat)

 The hydrostatic mapping function coefficients have to be used with the
 height dependent Vienna Mapping Function 1 (vmf_ht.f) because the
 coefficients refer to zero height.

 Example 1 (Vienna, 2 August 2012, with time variation):

 dmjd = 56141.d0
 dlat(1) = 48.20d0*pi/180.d0
 dlon(1) = 16.37d0*pi/180.d0
 hell(1) = 156.d0
 nstat = 1
 it = 0

 output:
 p = 1002.79 hPa
 T = 22.06 deg Celsius
 dT = -6.23 deg / km
 Tm = 281.30 K
 e = 16.65 hPa
 ah = 0.0012646
 aw = 0.0005752
 la = 2.6530
 undu = 45.76 m

 Example 2 (Vienna, 2 August 2012, without time variation, i.e. constant values):

 dmjd = 56141.d0
 dlat(1) = 48.20d0*pi/180.d0
 dlon(1) = 16.37d0*pi/180.d0
 hell(1) = 156.d0
 nstat = 1
 it = 1

 output:
 p = 1003.71 hPa
 T = 11.79 deg Celsius
 dT = -5.49 deg / km
 Tm = 273.22 K
 e = 10.22 hPa
 ah = 0.0012396
 aw = 0.0005753
 la = 2.6358
 undu = 45.76 m

 Klemens Lagler, 2 August 2012
 Johannes Boehm, 6 August 2012, revision
 Klemens Lagler, 21 August 2012, epoch change to January 1 2000
 Johannes Boehm, 23 August 2012, adding possibility to determine constant field
 Johannes Boehm, 27 December 2012, reference added
 Johannes Boehm, 10 January 2013, correction for dlat = -90 degrees
                                  (problem found by Changyong He)
 Johannes Boehm, 21 May 2013, bug with dmjd removed (input parameter dmjd was replaced
                 unintentionally; problem found by Dennis Ferguson)
 Gregory Pain,   17 June 2013, adding water vapor decrease factor la
 Gregory Pain,   21 June 2013, using the 1 degree grid : better for calculating zenith wet delays (la)
 Gregory Pain,   01 July 2013, adding mean temperature of the water vapor Tm
 Gregory Pain,   30 July 2013, changing the method to calculate the water vapor partial pressure (e)
 Gregory Pain,   31 July 2013, correction for (dlat = -90 degrees, dlon = 360 degrees)
 Johannes Boehm, 27 December 2013, copyright notice added
 Johannes Boehm, 25 August 2014, default input file changed to
                 gpt2_1wA.grd (slightly different humidity values)
 Johannes Boehm, 25 August 2014, reference changed to Boehm et al. in GPS
                 Solutions
 Johannes Boehm, 24 December 2014, converted to Fortran
 --------------------------------------------------------------------*/

    double lal[4]={0.0};

    double vec[44]={0.0}; //dimension vec(44)

    double cosfy=0.0, coshy=0.0, sinfy=0.0, sinhy=0.0;

    double undul[4]={0.0},Ql[4]={0.0},dTl[4]={0.0},Tl[4]={0.0},pl[4]={0.0},
           ahl[4]={0.0},awl[4]={0.0},Tml[4]={0.0},el[4]={0.0};

    double plon=0.0, ppod=0.0, diffpod=0.0, difflon=0.0, hgt=0.0, T0=0.0, p0=0.0,
           Q=0.0, redh=0.0, Tv=0.0, c=0.0, e0=0.0, Hs1=0.0, dnpod1=0.0,dnpod2=0.0,
           dnlon1=0.0, dnlon2=0.0, R1=0.0,R2=0.0;

    int ipod=0, ilon=0, ibilinear=0, ix=0, ipod1=0, ilon1=0;

    int indx[4]={0};

    char line[350]; //  character line*80

    FILE *myFile;

    //FILE * new_file; //new_file can be used to check the reading of the grid

    //mean gravity in m/s**2
    double gm = 9.80665;

    //molar mass of dry air in kg/mol
    double dMtr = 28.965e-3;

    //universal gas constant in J/K/mol
    double Rg = 8.3143;

    //change the reference epoch to January 1 2000
    double dmjd1 = dmjd-51544.5;

    //factors for amplitudes
    if (it==1) { //then ! constant parameters
        cosfy = 0.0;
        coshy = 0.0;
        sinfy = 0.0;
        sinhy = 0.0;
      }else {
        cosfy = cos(dmjd1/365.25*2.0*pi);
        coshy = cos(dmjd1/365.25*4.0*pi);
        sinfy = sin(dmjd1/365.25*2.0*pi);
        sinhy = sin(dmjd1/365.25*4.0*pi);
      } //end if

   //File to read grid
    myFile = fopen("gpt2_1wA.grd", "r");

    //Checks file
    if( myFile == NULL ) {
        printf("Cannot open file gpt2_1wA.grd\n");
        return;
    }

  /* new_file can be used to check the reading of the grid
    new_file = fopen("gpt2_out.grd", "w+");
    if( new_file == NULL ) {  //Checks file
        printf("Cannot open file gpt2_out.grd\n");
        return;
    }
*/

    rewind(myFile);  //force to pointer to begin of the file
    memset(line,'\0',sizeof(line));
    fgets(line,350,myFile); // read first comment line

    int n = 0;

    memset(Data,0,sizeof(Data)); //start Data

    do{
         memset(vec,0.0,sizeof(vec)); //start vec
         for (int i = 0; i < 44; i++) fscanf(myFile, "%lf", &vec[i]);

         memset(line,'\0',sizeof(line));
         fgets(line,350,myFile);       //read the remaining characters until find '\n'

         //if you want to check the read of the grid, please use the new_file
         //for (int i = 0; i < 44; i++) fprintf(new_file,"%.2lf ",vec[i]);
         //fprintf(new_file,"\n");

         for (int i = 0; i < 5; i++) Data[n].pgrid[i] = vec[i+2];            //pgrid(n,1:5)  = vec(3:7) -  pressure in Pascal
         for (int i = 0; i < 5; i++) Data[n].Tgrid[i] = vec[i+7];            //Tgrid(n,1:5)  = vec (8:12) - temperature in Kelvin
         for (int i = 0; i < 5; i++) Data[n].Qgrid[i] = vec[i+12]/1000.0;    //Qgrid(n,1:5)  = vec(13:17)/1000.d0 // specific humidity in kg/kg
         for (int i = 0; i < 5; i++) Data[n].dTgrid[i] = vec[i+17]/1000.0;   //dTgrid(n,1:5) = vec(18:22)/1000.d0 // temperature lapse rate in Kelvin/m
         u[n] = vec[22];                                                     // u(n) = vec(23)            // geoid undulation in m
         Hs[n] = vec[23];                                                    //Hs(n) = vec(24)            // orthometric grid height in m
         for (int i = 0; i < 5; i++) Data[n].ahgrid[i] = vec[i+24]/1000.0;   //ahgrid(n,1:5) = vec(25:29)/1000.d0 // hydrostatic mapping function coefficient, dimensionless
         for (int i = 0; i < 5; i++) Data[n].awgrid[i] = vec[i+29]/1000.0;   //awgrid(n,1:5) = vec(30:34)/1000.d0 // wet mapping function coefficient, dimensionless
         for (int i = 0; i < 5; i++) Data[n].lagrid[i] = vec[i+34];          //lagrid(n,1:5) = vec(35:39)         // water vapour decrease factor, dimensionless
         for (int i = 0; i < 5; i++) Data[n].Tmgrid[i] = vec[i+39];          //Tmgrid(n,1:5) = vec(40:44)         // weighted mean temperature, Kelvin

         n++;
   }while(n < 64800);

    fclose(myFile);

    //loop over stations

    int k=0;

    do{  //do  k = 1,nstat

        //only positive longitude in degrees
        if (dlon[k] < 0.0)
          plon = (dlon[k] + 2.0*pi)*180.0/pi;
        else
          plon = dlon[k]*180.0/pi;


        //transform to polar distance in degrees
        ppod = (-dlat[k] + pi/2.0)*180.0/pi;

        //find the index (line in the grid file) of the nearest point
        ipod = floor(ppod+1.0);
        ilon = floor(plon+1.0);

        // normalized (to one) differences, can be positive or negative
        diffpod = ppod - (ipod - 0.5);
        difflon = plon - (ilon - 0.5);

        //added by HCY
        if (ipod == 181) ipod = 180;

	    //added by GP
        if (ilon == 361) ilon = 1;

        if (ilon == 0) ilon = 360;

        // get the number of the corresponding line
        indx[0] = (ipod - 1)*360 + ilon;
        indx[0] -= 1;  //-1 to use in c

        // near the poles: nearest neighbour interpolation, otherwise: bilinear
	    // with the 1 degree grid the limits are lower and upper (GP)
        ibilinear = 0;

        if( (ppod > 0.5) && (ppod < 179.5) ) ibilinear = 1;

        // case of nearest neighborhood
        if (ibilinear == 0) { // then

          ix = indx[0] - 1;  //-1 to use in c

          // transforming ellipsoidal height to orthometric height
          undu[k] = u[ix];
          hgt = hell[k]-undu[k];

          // pressure, temperature at the height of the grid
          T0 = Data[ix].Tgrid[0] +
               Data[ix].Tgrid[1]*cosfy +Data[ix]. Tgrid[2]*sinfy +
               Data[ix].Tgrid[3]*coshy + Data[ix].Tgrid[4]*sinhy;

          p0 = Data[ix].pgrid[0] +
               Data[ix].pgrid[1]*cosfy + Data[ix].pgrid[2]*sinfy+
               Data[ix].pgrid[3]*coshy + Data[ix].pgrid[4]*sinhy;

          // specific humidity
          Q = Data[ix].Qgrid[0] +
              Data[ix].Qgrid[1]*cosfy + Data[ix].Qgrid[2]*sinfy+
              Data[ix].Qgrid[3]*coshy + Data[ix].Qgrid[4]*sinhy;

          // lapse rate of the temperature
          dT[k] = Data[ix].dTgrid[0] +
                  Data[ix].dTgrid[1]*cosfy + Data[ix].dTgrid[2]*sinfy+
                  Data[ix].dTgrid[3]*coshy + Data[ix].dTgrid[4]*sinhy;

          // station height - grid height
          redh = hgt - Hs[ix];

          // temperature at station height in Celsius
          T[k] = T0 + dT[k]*redh - 273.150;

          // temperature lapse rate in degrees / km
          dT[k] = dT[k]*1000.0;

          // virtual temperature in Kelvin
          Tv = T0*(1+0.6077*Q);

          c = gm*dMtr/(Rg*Tv);

          // pressure in hPa
          p[k] = (p0*exp(-c*redh))/100.0;

          // hydrostatic coefficient ah
          ah[k] = Data[ix].ahgrid[0] +
                  Data[ix].ahgrid[1]*cosfy + Data[ix].ahgrid[2]*sinfy+
                  Data[ix].ahgrid[3]*coshy + Data[ix].ahgrid[4]*sinhy;

          // wet coefficient aw
          aw[k] = Data[ix].awgrid[0] +
                  Data[ix].awgrid[1]*cosfy + Data[ix].awgrid[2]*sinfy +
                  Data[ix].awgrid[3]*coshy + Data[ix].awgrid[4]*sinhy;

		  // water vapour decrease factor la - added by GP
          la[k] = Data[ix].lagrid[0] +
                  Data[ix].lagrid[1]*cosfy + Data[ix].lagrid[2]*sinfy +
                  Data[ix].lagrid[3]*coshy + Data[ix].lagrid[4]*sinhy;

		  // mean temperature of the water vapor Tm - added by GP
          Tm[k] = Data[ix].Tmgrid[0] +
                  Data[ix].Tmgrid[1]*cosfy + Data[ix].Tmgrid[2]*sinfy +
                  Data[ix].Tmgrid[3]*coshy + Data[ix].Tmgrid[4]*sinhy;

		  // water vapor pressure in hPa - changed by GP
		  e0 = Q*p0/(0.622 + 0.378*Q)/100.0;  // on the grid

		  e[k] = e0* pow( (100.0*p[k]/p0),(la[k]+1) );  //on the station height - (14) Askne and Nordius, 1987
        }
        else { // bilinear interpolation


          ipod1 = ipod + int(sign(1.0,diffpod));
          ilon1 = ilon + int(sign(1.0,difflon));

		  // changed for the 1 degree grid (GP)
          if (ilon1 == 361) ilon1 = 1;

          if (ilon1 == 0) ilon1 = 360;

          // get the number of the line
          indx[1] = (ipod1 - 1)*360 + ilon;   // along same longitude
          indx[2] = (ipod  - 1)*360 + ilon1;  // along same polar distance
          indx[3] = (ipod1 - 1)*360 + ilon1;  // diagonal

          indx[1] -= 1;  //-1 to use in c
          indx[2] -= 1;  //-1 to use in c
          indx[3] -= 1;  //-1 to use in c

         int l = 0;

         do{ // do l = 1,4

            //transforming ellipsoidal height to orthometric height:
            //Hortho = -N + Hell
            int ix = indx[l];

            undul[l] = u[ix];
            hgt = hell[k]-undul[l];



            //pressure, temperature at the height of the grid
            T0 = Data[ix].Tgrid[0] +
                 Data[ix].Tgrid[1]*cosfy + Data[ix].Tgrid[2]*sinfy +
                 Data[ix].Tgrid[3]*coshy + Data[ix].Tgrid[4]*sinhy;

            p0 = Data[ix].pgrid[0] +
                 Data[ix].pgrid[1]*cosfy + Data[ix].pgrid[2]*sinfy +
                 Data[ix].pgrid[3]*coshy + Data[ix].pgrid[4]*sinhy;

            //humidity
            Ql[l] = Data[ix].Qgrid[0] +
                    Data[ix].Qgrid[1]*cosfy + Data[ix].Qgrid[2]*sinfy +
                    Data[ix].Qgrid[3]*coshy + Data[ix].Qgrid[4]*sinhy;

            //reduction = stationheight - gridheight
            Hs1 = Hs[ix];
            redh = hgt - Hs1;

            //lapse rate of the temperature in degree / m
            dTl[l] = Data[ix].dTgrid[0] +
                     Data[ix].dTgrid[1]*cosfy + Data[ix].dTgrid[2]*sinfy +
                     Data[ix].dTgrid[3]*coshy + Data[ix].dTgrid[4]*sinhy;

            //temperature reduction to station height
            Tl[l] = T0 + dTl[l]*redh - 273.150;

            //virtual temperature
            Tv = T0*(1+0.6077*Ql[l]);
            c = gm*dMtr/(Rg*Tv);

            //pressure in hPa
            pl[l] = (p0*exp(-c*redh))/100.0;

            //hydrostatic coefficient ah
            ahl[l] = Data[ix].ahgrid[0] +
                     Data[ix].ahgrid[1]*cosfy + Data[ix].ahgrid[2]*sinfy +
                     Data[ix].ahgrid[3]*coshy + Data[ix].ahgrid[4]*sinhy;

            //wet coefficient aw
            awl[l] = Data[ix].awgrid[0] +
                     Data[ix].awgrid[1]*cosfy + Data[ix].awgrid[2]*sinfy +
                     Data[ix].awgrid[3]*coshy + Data[ix].awgrid[4]*sinhy;

            //water vapor decrease factor la - added by GP
            lal[l] = Data[ix].lagrid[0] +
                     Data[ix].lagrid[1]*cosfy + Data[ix].lagrid[2]*sinfy +
                     Data[ix].lagrid[3]*coshy + Data[ix].lagrid[4]*sinhy;

            //mean temperature of the water vapor Tm - added by GP
            Tml[l] = Data[ix].Tmgrid[0] +
                     Data[ix].Tmgrid[1]*cosfy + Data[ix].Tmgrid[2]*sinfy +
                     Data[ix].Tmgrid[3]*coshy + Data[ix].Tmgrid[4]*sinhy;

            //water vapor pressure in hPa - changed by GP
            e0 = Ql[l]*p0/(0.622 + 0.378*Ql[l])/100.0; // on the grid
            el[l] = e0 * pow( (100.0*pl[l]/p0), (lal[l]+1.0) );  // on the station height  (14) Askne and Nordius, 1987

            l++;

          }while (l< 4);//end do


          dnpod1 = fabs(diffpod); // distance nearer point
          dnpod2 = 1.0 - dnpod1;   // distance to distant point
          dnlon1 = fabs(difflon);
          dnlon2 = 1.0 - dnlon1;

          //pressure
          R1 = dnpod2*pl[0]+dnpod1*pl[1];
          R2 = dnpod2*pl[2]+dnpod1*pl[3];
          p[k] = dnlon2*R1+dnlon1*R2;

          //temperature
          R1 = dnpod2*Tl[0]+dnpod1*Tl[1];
          R2 = dnpod2*Tl[2]+dnpod1*Tl[3];
          T[k] = dnlon2*R1+dnlon1*R2;

          //temperature in degree per km
          R1 = dnpod2*dTl[0]+dnpod1*dTl[1];
          R2 = dnpod2*dTl[2]+dnpod1*dTl[3];
          dT[k] = (dnlon2*R1+dnlon1*R2)*1000.0;

          //water vapor pressure in hPa - changed by GP
		  R1 = dnpod2*el[0]+dnpod1*el[1];
          R2 = dnpod2*el[2]+dnpod1*el[3];
          e[k] = dnlon2*R1+dnlon1*R2;

          //hydrostatic
          R1 = dnpod2*ahl[0]+dnpod1*ahl[1];
          R2 = dnpod2*ahl[2]+dnpod1*ahl[3];
          ah[k] = dnlon2*R1+dnlon1*R2;

          //wet
          R1 = dnpod2*awl[0]+dnpod1*awl[1];
          R2 = dnpod2*awl[2]+dnpod1*awl[3];
          aw[k] = dnlon2*R1+dnlon1*R2;

          //undulation
          R1 = dnpod2*undul[0]+dnpod1*undul[1];
          R2 = dnpod2*undul[2]+dnpod1*undul[3];
          undu[k] = dnlon2*R1+dnlon1*R2;

          //water vapor decrease factor la - added by GP
          R1 = dnpod2*lal[0]+dnpod1*lal[1];
          R2 = dnpod2*lal[2]+dnpod1*lal[3];
          la[k] = dnlon2*R1+dnlon1*R2;

          //mean temperature of the water vapor Tm - added by GP
          R1 = dnpod2*Tml[0]+dnpod1*Tml[1];
          R2 = dnpod2*Tml[2]+dnpod1*Tml[3];
          Tm[k] = dnlon2*R1+dnlon1*R2;

        } //else

        k++;

      }while(k < nstat); // end do


}  //end subroutine
//---------------------------------------------------------------------------
void asknewet (double e,double Tm,double lambda,double &zwd)
{
/* --------------------------------------------------------------------
 This subroutine determines the zenith wet delay based on the
 equation 22 by Aske and Nordius (1987)

 c Reference:
 Askne and Nordius, Estimation of tropospheric delay for microwaves from
 surface weather data, Radio Science, Vol 22(3): 379-386, 1987.

 input parameters:

 e:      water vapor pressure in hPa
 Tm:     mean temperature in Kelvin
 lambda: water vapor lapse rate (see definition in Askne and Nordius 1987)

 output parameters:

 zwd:  zenith wet delay in meter

 Example 1 :

 e =  10.9621 hPa
 Tm = 273.8720
 lambda = 2.8071

 output:
 zwd = 0.1176 m

 Johannes Boehm, 3 August 2013
 Johannes Boehm, 24 December 2014, converted to Fortran
 --------------------------------------------------------------------*/

  // e =  10.9621; // hPa
  // Tm = 273.8720;
  // lambda = 2.8071;

      double k1=0.0,k2=0.0,k2p=0.0,k3=0.0;

      //coefficients
      k1  = 77.604;                     //K/hPa
      k2 = 64.79;                       //K/hPa
      k2p = k2 - k1*18.0152/28.9644;    //K/hPa
      k3  = 377600.0;                   //KK/hPa

      //mean gravity in m/s**2
      double gm = 9.80665;

      //molar mass of dry air in kg/mol
      double dMtr = 28.965e-3;

      //universal gas constant in J/K/mol
      double R = 8.3143;

      //specific gas constant for dry consituents
      double Rd = R/dMtr;

      zwd = 1.0e-6*(k2p + k3/Tm)*Rd/(lambda + 1.0)/gm*e;

}
//---------------------------------------------------------------------------
void saasthyd (double p,double dlat,double hell,double &zhd)
{
/*
 This subroutine determines the zenith hydrostatic delay based on the
 equation by Saastamoinen (1972) as refined by Davis et al. (1985)

 c Reference:
 Saastamoinen, J., Atmospheric correction for the troposphere and
 stratosphere in radio ranging of satellites. The use of artificial
 satellites for geodesy, Geophys. Monogr. Ser. 15, Amer. Geophys. Union,
 pp. 274-251, 1972.
 Davis, J.L, T.A. Herring, I.I. Shapiro, A.E.E. Rogers, and G. Elgered,
 Geodesy by Radio Interferometry: Effects of Atmospheric Modeling Errors
 on Estimates of Baseline Length, Radio Science, Vol. 20, No. 6,
 pp. 1593-1607, 1985.

 input parameters:

 p:     pressure in hPa
 dlat:  ellipsoidal latitude in radians
 dlon:  longitude in radians
 hell:  ellipsoidal height in m

 output parameters:

 zhd:  zenith hydrostatic delay in meter

 Example 1 :

 p = 1000;
 dlat = 48d0*pi/180.d0
 hell = 200.d0

 output:
 zhd = 2.2695 m

 Johannes Boehm, 8 May 2013
 Johannes Boehm, 24 December 2014, converted to Fortran
 ---
*/

      //calculate denominator f
      double f = 1.0-0.00266*cos(2.0*dlat) - 0.00000028*hell;

      //calculate the zenith hydrostatic delay
      zhd = 0.0022768*p/f;

} // end subroutine
//---------------------------------------------------------------------------
 void vmf1_ht(double ah, double aw, double dmjd, double dlat,
	         double ht, double zd, double &vmf1h, double &vmf1w)
{
/*---------------------------------------------------------------------------
	      !!! This is the version with height correction !!!
	      !!! It has to be used with the grid !!!

	      This subroutine determines the VMF1 (Vienna Mapping Functions 1)
	      Reference: Boehm, J., B. Werl, H. Schuh (2006),
	      Troposphere mapping functions for GPS and very long baseline interferometry
	      from European Centre for Medium-Range Weather Forecasts operational analysis data,
	      J. Geoph. Res., Vol. 111, B02406, doi:10.1029/2005JB003629.

	      input data
	      ----------
	      ah:   hydrostatic coefficient a (www.hg.tuwien.ac.at/~ecmwf1)
	      aw:   wet coefficient a         (www.hg.tuwien.ac.at/~ecmwf1)
	      dmjd: modified julian date
	      dlat: latitude in radians
	      ht:   ellipsoidal height in meter
	      zd:   zenith distance in radians

	      output data
	      -----------
	      vmf1h: hydrostatic mapping function
	      vmf1w: wet mapping function

	      Johannes Boehm, 2005 October 2
	 */

	double a_ht, b_ht, beta, bh, bw, c0h, c10h, c11h, c_ht, ch, cw,
	      doy, gamma, hs_km, ht_corr, ht_corr_coef, phh, sine, topcon;

	//     reference day is 28 January
	//     this is taken from Niell (1996) to be consistent

	doy = dmjd - 44239.e0 + 1 - 28;

	bh = 0.0029;
	c0h = 0.062;
	if( dlat < 0.e0 ){ /* southern hemisphere*/
		phh = pi;
		c11h = 0.007;
		c10h = 0.002;
		}
	else{ /* northern hemisphere*/
		phh = 0.e0;
		c11h = 0.005;
		c10h = 0.001;
		}
	ch = c0h + ((cos(doy/365.25e0*2.e0*pi+phh) + 1.0)*c11h/2.e0 + c10h)*(1.0 - cos(dlat));

	sine = sin(pi/2.e0-zd);
	beta = bh/(sine + ch);
	gamma = ah/(sine + beta);
	topcon = 1.0 + ah/(1.0 + bh/(1.0 + ch));
	vmf1h = topcon/(sine + gamma);

	//  height correction [Niell, 1996]
	a_ht = 2.53e-5;
	b_ht = 5.49e-3;
	c_ht = 1.14e-3;
	hs_km = ht/1000.e0;
	beta = b_ht/(sine + c_ht);
	gamma = a_ht/(sine + beta);
	topcon = 1.0 + a_ht/(1.0 + b_ht/(1.0 + c_ht));
	ht_corr_coef = 1.0/sine - topcon/(sine + gamma);
	ht_corr = ht_corr_coef*hs_km;
	vmf1h += ht_corr;

	bw = 0.00146;
	cw = 0.04391;
	beta = bw/(sine + cw);
	gamma = aw/(sine + beta);
	topcon = 1.0 + aw/(1.0 + bw/(1.0 + cw));
	vmf1w = topcon/(sine + gamma);

	return;
}


