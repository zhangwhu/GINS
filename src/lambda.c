/*------------------------------------------------------------------------------
* lambda.c : integer ambiguity resolution
*
*          Copyright (C) 2007-2008 by T.TAKASU, All rights reserved.
*
* reference :
*     [1] P.J.G.Teunissen, The least-square ambiguity decorrelation adjustment:
*         a method for fast GPS ambiguity estimation, J.Geodesy, Vol.70, 65-82,
*         1995
*     [2] X.-W.Chang, X.Yang, T.Zhou, MLAMBDA: A modified LAMBDA method for
*         integer least-squares estimation, J.Geodesy, Vol.79, 552-565, 2005
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
* history : 2007/01/13 1.0 new
*           2015/05/31 1.1 add api lambda_reduction(), lambda_search()
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/

#define LOOPMAX     10000           /* maximum count of search loop */

#define SGN(x)      ((x)<=0.0?-1.0:1.0)
#define ROUND(x)    (floor((x)+0.5))
#define SWAP(x,y)   do {double tmp_; tmp_=x; x=y; y=tmp_;} while (0)
#define MIN(x,y)        ((x)<(y)?(x):(y))

/* LD factorization (Q=L'*diag(D)*L) -----------------------------------------*/
extern int LD(int n, const double *Q, double *L, double *D)
{
    int i,j,k,info=0;
    double a,*A=mat(n,n);
    
    memcpy(A,Q,sizeof(double)*n*n);
    for (i=n-1;i>=0;i--) {
		A[i + i*n] += 0.000001;
        if ((D[i]=A[i+i*n])<=0.0) {
			info=-1; break;
		}
        a=sqrt(D[i]);
        for (j=0;j<=i;j++) L[i+j*n]=A[i+j*n]/a;
        for (j=0;j<=i-1;j++) for (k=0;k<=j;k++) A[j+k*n]-=L[i+k*n]*L[i+j*n];				
        for (j=0;j<=i;j++) L[i+j*n]/=L[i+i*n];
    }
    free(A);
    if (info) fprintf(stderr,"%s : LD factorization error\n",__FILE__);
    return info;
}
/* integer gauss transformation ----------------------------------------------*/
static void gauss(int n, double *L, double *Z, int i, int j)
{
    int k,mu;
    
    if ((mu=(int)ROUND(L[i+j*n]))!=0) {
        for (k=i;k<n;k++) L[k+n*j]-=(double)mu*L[k+i*n];
        for (k=0;k<n;k++) Z[k+n*j]-=(double)mu*Z[k+i*n];
    }
}
/* permutations --------------------------------------------------------------*/
static void perm(int n, double *L, double *D, int j, double del, double *Z)
{
    int k;
    double eta,lam,a0,a1;
    
    eta=D[j]/del;
    lam=D[j+1]*L[j+1+j*n]/del;
    D[j]=eta*D[j+1]; D[j+1]=del;
    for (k=0;k<=j-1;k++) {
        a0=L[j+k*n]; a1=L[j+1+k*n];
        L[j+k*n]=-L[j+1+j*n]*a0+a1;
        L[j+1+k*n]=eta*a0+lam*a1;
    }
    L[j+1+j*n]=lam;
    for (k=j+2;k<n;k++) SWAP(L[k+j*n],L[k+(j+1)*n]);
    for (k=0;k<n;k++) SWAP(Z[k+j*n],Z[k+(j+1)*n]);
}
/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) (ref.[1]) ---------------*/
static void reduction(int n, double *L, double *D, double *Z)
{
    int i,j,k;
    double del;
    
    j=n-2; k=n-2;
    while (j>=0) {
        if (j<=k) for (i=j+1;i<n;i++) gauss(n,L,Z,i,j);
        del=D[j]+L[j+1+j*n]*L[j+1+j*n]*D[j+1];
        if (del+1E-6<D[j+1]) { /* compared considering numerical error */
            perm(n,L,D,j,del,Z);
            k=j; j=n-2;
        }
        else j--;
    }
}
/* modified lambda (mlambda) search (ref. [2]) -------------------------------
* args   : n      I  number of float parameters
*          m      I  number of fixed solution
           L,D    I  transformed covariance matrix
           zs     I  transformed double-diff phase biases
           zn     O  fixed solutions
           s      O  sum of residuals for fixed solutions                    */
static int search(int n, int m, const double *L, const double *D,
                  const double *zs, double *zn, double *s)
{
    int i,j,k,c,nn=0,imax=0;
    double newdist,maxdist=1E99,y;
    double *S=zeros(n,n),*dist=mat(n,1),*zb=mat(n,1),*z=mat(n,1),*step=mat(n,1);
    
    k=n-1; dist[k]=0.0;						
	zb[k] = zs[k];			//�����
	z[k] = ROUND(zb[k]);	//�ͽ��̶�
	y = zb[k] - z[k];
    step[k]=SGN(y);  /* step towards closest integer */
    for (c=0;c<LOOPMAX;c++) {			
        newdist=dist[k]+y*y/D[k];  /* newdist=sum(((z(j)-zb(j))^2/d(j))) */			
        if (newdist<maxdist) {
            /* Case 1: move down */
            if (k!=0) {
                dist[--k]=newdist;				
                for (i=0;i<=k;i++)
                    S[k+i*n]=S[k+1+i*n]+(z[k+1]-zb[k+1])*L[k+1+i*n];
                zb[k]=zs[k]+S[k+k*n];
                z[k]=ROUND(zb[k]); /* next valid integer */
                y=zb[k]-z[k];
                step[k]=SGN(y);
            }
            /* Case 2: store the found candidate and try next valid integer */
            else {     
                if (nn<m) {  /* store the first m initial points */
                    if (nn==0||newdist>s[imax]) imax=nn;					
                    for (i=0;i<n;i++) zn[i+nn*n]=z[i];			
                    s[nn++]=newdist;		
                }
                else {
                    if (newdist<s[imax]) {
                        for (i=0;i<n;i++) zn[i+imax*n]=z[i];
                        s[imax]=newdist;
                        for (i=imax=0;i<m;i++) if (s[imax]<s[i]) imax=i;
                    }
                    maxdist=s[imax];				
                }
                z[0]+=step[0]; /* next valid integer */
                y=zb[0]-z[0];
                step[0]=-step[0]-SGN(step[0]);		
            }
        }
        /* Case 3: exit or move up */
        else {					
            if (k==n-1) break;
            else {
                k++;  /* move up */
                z[k]+=step[k];  /* next valid integer */
                y=zb[k]-z[k];
                step[k]=-step[k]-SGN(step[k]);
            }
        }
    }
    for (i=0;i<m-1;i++) { /* sort by s */
        for (j=i+1;j<m;j++) {								
            if (s[i]<s[j]) continue;				//����
            SWAP(s[i],s[j]);    
            for (k=0;k<n;k++) SWAP(zn[k+i*n],zn[k+j*n]);   
        }
    }
    free(S); free(dist); free(zb); free(z); free(step);
    
    if (c>=LOOPMAX) {
        fprintf(stderr,"%s : search loop count overflow\n",__FILE__);
        return -1;
    }
    return 0;
}
/*Cumulative density function of normal distribution (use Taylor series)��̬�ֲ��ۼ��ܶȺ���*/
extern double normcdf(double x)
{
	long double y = 0.0, z = 1.0, sum, b = 1.0;
	int n = 1, a = 1, m = 1;

	if (x >= 4)
		return 1.0;

	while (fabs(z)>1e-15) {
		z = a*pow(x, n) / (b*n);
		y += z;
		a *= -1;
		b *= 2 * m;
		m++;
		n += 2;
	}

	sum = 0.5 + y / sqrt(2 * PI);
	/*double sum1 = 0.5*erfc(-x*sqrt(0.5));
	trace(1, "%f %10.4f %10.4f\n", x, sum, sum1);*/
	return sum;
}
/* lambda/mlambda integer least-square estimation ------------------------------
* integer least-square estimation. reduction is performed by lambda (ref.[1]),
* and search by mlambda (ref.[2]).		����ģ���ȹ̶�
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1) (double-diff phase biases)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
* notes  : matrix stored by column-major order (fortran convension)
*-----------------------------------------------------------------------------*/
extern int lambda(int n, int m, const double *a, const double *Q, double *F, double *s)          
{
    int i,info;
    double *L,*D,*Z,*z,*E, P=1.0, pP=0.0;
    
    if (n<=0||m<=0) return -1;
    L=zeros(n,n); D=mat(n,1); Z=eye(n); z=mat(n,1); E=mat(n,m);
    
    /* LD (lower diaganol) factorization (Q=L'*diag(D)*L) */
    if (!(info=LD(n,Q,L,D))) {      
		
        /* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) */
		reduction(n, L, D, Z);			
        matmul("TN",n,1,n,1.0,Z,a,0.0,z); /* z=Z'*a */
		for (i = 0; i < n; i++)
			P = P*(2 * (normcdf(1 / (2 * sqrt(D[i])))) - 1);
        /* mlambda search 
            z = transformed double-diff phase biases
            L,D = transformed covariance matrix */
        if (!(info=search(n,m,L,D,z,E,s))) {  /* returns 0 if no error */   
            
            info=solve("T",Z,E,n,m,F); /* F=Z'\E */					
        }
		//pP = exp(-0.5*s[0])/exp(-0.5*(s[0]+s[1]+s[2]));
    }
    free(L); free(D); free(Z); free(z); free(E);
	//if (pP<0.99) s[0] = s[1] = 1;
	if (P < 0.3 /*|| pP<0.99*/) s[0] = s[1] = 1;
    return info;
}
/* lambda reduction ------------------------------------------------------------
* reduction by lambda (ref [1]) for integer least square
* args   : int    n      I  number of float parameters
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *Z     O  lambda reduction matrix (n x n)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_reduction(int n, const double *Q, double *Z)
{
    double *L,*D;
    int i,j,info;
    
    if (n<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    for (i=0;i<n;i++) for (j=0;j<n;j++) {
        Z[i+j*n]=i==j?1.0:0.0;
    }
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {
        free(L); free(D);
        return info;
    }
    /* lambda reduction */
    reduction(n,L,D,Z);				
     
    free(L); free(D);
    return 0;
}
/* mlambda search --------------------------------------------------------------
* search by  mlambda (ref [2]) for integer least square
* args   : int    n      I  number of float parameters
*          int    m      I  number of fixed solutions
*          double *a     I  float parameters (n x 1)
*          double *Q     I  covariance matrix of float parameters (n x n)
*          double *F     O  fixed solutions (n x m)
*          double *s     O  sum of squared residulas of fixed solutions (1 x m)
* return : status (0:ok,other:error)
*-----------------------------------------------------------------------------*/
extern int lambda_search(int n, int m, const double *a, const double *Q,
                         double *F, double *s)			
{
    double *L,*D;
    int info;
    
    if (n<=0||m<=0) return -1;
    
    L=zeros(n,n); D=mat(n,1);
    
    /* LD factorization */
    if ((info=LD(n,Q,L,D))) {		
        free(L); free(D);
        return info;
    }
    /* mlambda search */
    info=search(n,m,L,D,a,F,s);
    
    free(L); free(D);
    return info;
}
/* lambda par */
extern int lambda_PAR(rtk_t *rtk, int na, int nb, int m, const double *a, double *b, 
    const double *Qaa, double *Qbb, const double *Qba)			
{
	int i, j, k, info, ufix, fixN = 0;
	double *F, *E, *Z, *Z2, *Qz2z2, *Qbz2, *Qbb_, *s, P, trj1, trj2;
	double *L, *D, *z, *Lp, *Dp, *zp, *dz2;

	if (na <= 0 || m <= 0) return 0;

	F = mat(na, m);			E = mat(nb, nb);		
	Z = eye(na);			Z2 = mat(na, na);	
	s = mat(m, 1);			dz2 = mat(na, 1);
	Qz2z2 = mat(na, na);	Qbz2 = mat(nb, na);		Qbb_ = mat(nb, nb);		
	L  = zeros(na, na);		D  = mat(na, 1);	    z = mat(na, 1);			
	Lp = mat(na, na);		Dp = mat(na, 1);		zp = mat(na, 1);		

	/* LD (lower diaganol) factorization (Q=L'*diag(D)*L) */
	if (!(info = LD(na, Qaa, L, D))) {

		/* lambda reduction (z=Z'*a, Qz=Z'*Q*Z=L'*diag(D)*L) */
		reduction(na, L, D, Z);

		matmul("TN", na, 1, na, 1.0, Z, a, 0.0, z);				/* z=Z'*a */
		
		for (k = na; k >= 4; k--) {			
			ufix = na - k;
			for (i = ufix, P = 1.0; i < na; i++) {
				P = P*(2 * (normcdf(0.5 / sqrt(D[i]))) - 1);
				for (j = 0; j < na; j++) {
					Z2[j + (i - ufix)*na] = Z[j + i*na];
				}
			}
			if (P < 0.99) continue;

			matmul("NN", nb, k, na, 1.0, Qba, Z2, 0.0, Qbz2);
			matmul("TN", k, na, na, 1.0, Z2, Qaa, 0.0, E);
			matmul("NN", k, k, na, 1.0, E, Z2, 0.0, Qz2z2);

			if (info = matinv(Qz2z2, k)) continue;				/*Qz2z2-1 */
			
			matmul("NN", nb, k, k, 1.0, Qbz2, Qz2z2, 0.0, E);
			matmul("NT", nb, nb, k, 1.0, E, Qbz2, 0.0, Qbb_);

			for (i = trj1 = trj2 = 0; i < 3; i++) {		
				trj1 += Qbb[i + i*nb];
				trj2 += Qbb[i + i*nb] - Qbb_[i + i*nb];
			}
			//if (sqrt(trj1 / trj2) < 50) continue;			

			for (i = 0; i < k; i++) {
				for (j = 0; j < k; j++) {
					Lp[j + i*k] = L[(ufix + j) + (ufix + i)*na];
				}
				Dp[i] = D[ufix + i];
				zp[i] = z[ufix + i];
			}
			/* mlambda search */
			if (!(info = search(k, m, Lp, Dp, zp, F, s))) {			  /* returns 0 if no error */					
				rtk->sol.ratio = MIN(s[1] / s[0], 999.9);
				if (rtk->opt.thresar[0] > 0.0 && rtk->sol.ratio >= rtk->opt.thresar[0]) {
					for (i = 0; i < k; i++) {
						dz2[i] = zp[i] - F[i];
					}
					matmul("NN", nb, 1, k, -1.0, E, dz2, 1.0, b);
					for (i = 0; i < nb; i++) for (j = 0; j < nb; j++) {
						Qbb[j + i*nb] -= Qbb_[j + i*nb];
					}
					fixN = k;	break;
				}
			}
		}
	}

	free(F);		free(E);	
	free(Z);		free(Z2);
	free(s);		free(dz2);
	free(Qz2z2);	free(Qbz2);		free(Qbb_);  
	free(L);		free(D);		free(z);		
	free(Lp);		free(Dp);		free(zp);
	return fixN;
}