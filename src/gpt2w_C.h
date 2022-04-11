#ifndef gpt2h
#define gpt2h
//---------------------------------------------------------------------------
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
//---------------------------------------------------------------------------
using namespace std;

//------------- Header of Functions ----------------------------------------
void vmf1_ht(double ah, double aw, double dmjd, double dlat,
	         double ht, double zd, double &vmf1h, double &vmf1w);

void saasthyd (double p,double dlat,double hell,double &zhd);

void asknewet (double e,double Tm,double lambda,double &zwd);

void gpt2_1w (double dmjd,double dlat[],double dlon[],double hell[],int nstat,int it,
                   double p[],double T[],double dT[],double Tm[],double e[],
                   double ah[],double aw[],double la[],double undu[]);
//---------------------------------------------------------------------------
//global variables
int const Max_Dimension =  64800;
double u[Max_Dimension]={0.0},Hs[Max_Dimension]={0.0};

/*double pgrid[Max_Dimension][5],Tgrid[Max_Dimension][5],Qgrid[Max_Dimension][5],
       dTgrid[Max_Dimension][5],u[Max_Dimension],Hs[Max_Dimension],ahgrid[Max_Dimension][5],
       awgrid[Max_Dimension][5],Tmgrid[Max_Dimension][5],lagrid[Max_Dimension][5]={0.0};
 */

typedef struct data{

  double pgrid[5],
         Tgrid[5],
         Qgrid[5],
         dTgrid[5],
         ahgrid[5],
         awgrid[5],
         Tmgrid[5],
         lagrid[5];

}GPT2_Data;

#endif
