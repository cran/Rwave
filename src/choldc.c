/***************************************************************
*    $Log: cholde.c,v $                                        *
****************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
****************************************************************/
/* Numeric Recipe of C */

#include "Swave.h"
#include "dyadic.h"

#ifndef _NR_UTILS_H_
#define _NR_UTILS_H_

static float sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static float maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static float minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b),(minarg1) < (minarg2) ?\
        (minarg1) : (minarg2))

static long lmaxarg1,lmaxarg2;
#define LMAX(a,b) (lmaxarg1=(a),lmaxarg2=(b),(lmaxarg1) > (lmaxarg2) ?\
        (lmaxarg1) : (lmaxarg2))

static long lminarg1,lminarg2;
#define LMIN(a,b) (lminarg1=(a),lminarg2=(b),(lminarg1) < (lminarg2) ?\
        (lminarg1) : (lminarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b),(imaxarg1) > (imaxarg2) ?\
        (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void nrerror(char error_text[]);
float *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

#else /* ANSI */
/* traditional - K&R */

void nrerror();
float *vector();
float **matrix();
float **submatrix();
float **convert_matrix();
float ***f3tensor();
double *dvector();
double **dmatrix();
int *ivector();
int **imatrix();
unsigned char *cvector();
unsigned long *lvector();
void free_vector();
void free_dvector();
void free_ivector();
void free_cvector();
void free_lvector();
void free_matrix();
void free_submatrix();
void free_convert_matrix();
void free_dmatrix();
void free_imatrix();
void free_f3tensor();

#endif /* ANSI */

#endif /* _NR_UTILS_H_ */

/*****************************************************************************/
/* choldc, cholsl */
/* ludcmp, lubksb */
/*****************************************************************************/

void float_choldc(float **a, int n, float p[])
{
  int i,j,k;
  float sum;
  
  for (i=1;i<=n;i++) {
    for (j=i;j<=n;j++) {
      for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
      if (i == j) {
	if (sum <= 0.0)
	  printf("choldc failed");
	p[i]=sqrt(sum);
      } else a[j][i]=sum/p[i];
    }
  }
}

/*****************************************************************************/
/* cholde : double precision with C language array format                    */
/*****************************************************************************/

void choldc(float **a, int n, float p[])
{
  float** A;
  float *P;
  int i,j,k;

/*  if(!(A = (float **)(malloc(sizeof(float *) * (n+1)))))
    error("Memory allocation failed for A in choldc.c \n");
*/
  if(!(P = (float *)(malloc(sizeof(float) * (n+1)))))
    error("Memory allocation failed for P in choldc.c \n");
/*		  
  for(i = 0; i <= n; i++) {
    if(!(A[i] = (float *)(malloc(sizeof(float) * (n+1)))))
      error("Memory allocation failed for A in choldc.c \n");
  }
*/  
  for(i = 0; i < n; i++) {
    P[i+1] = (float)(p[i]);
/*    for(j = 0; j <n ; j++)
      A[i+1][j+1] = (float)(a[i][j]); */
  }
    
/*  float_choldc(A,n,P); */
  float_choldc(a,n,P);


  for(i = 0; i < n; i++) {
    p[i] = (float)(P[i+1]);
/*    for(j = 0; j <n ; j++)
      a[i][j] = (float)(A[i+1][j+1]); */
  }
  
/*  for(i = 0; i <= n; i++) {
    free(A[i]);
  }
*/
  free(P);
/*  free(A); */
}



/***************************************/
/* cholsl with C language array format */
/***************************************/
void float_cholsl(float **a, int n, float p[], float b[], float x[])
{
  int i,k;
  float sum;
  
  for (i=1;i<=n;i++) {
    for (sum=b[i],k=i-1;k>=1;k--) sum -= a[i][k]*x[k];
    x[i]=sum/p[i];
  }
  for (i=n;i>=1;i--) {
    for (sum=x[i],k=i+1;k<=n;k++) sum -= a[k][i]*x[k];
    x[i]=sum/p[i];
  }
}

void cholsl(float **a, int n, float p[], float b[], float x[])
{
  float** A;
  float *P;
  int i,j,k;
  float *B, *X;
  
/*  if(!(A = (float **)(malloc(sizeof(float *) * (n+1)))))
    error("Memory allocation failed for A in choldc.c \n");
*/
  if(!(P = (float *)(malloc(sizeof(float) * (n+1)))))
    error("Memory allocation failed for P in choldc.c \n");
  if(!(B = (float *)(malloc(sizeof(float) * (n+1)))))
    error("Memory allocation failed for B in choldc.c \n");
  if(!(X = (float *)(malloc(sizeof(float) * (n+1)))))
    error("Memory allocation failed for X in choldc.c \n");
		  
/*  for(i = 0; i <= n; i++) {
    if(!(A[i] = (float *)(malloc(sizeof(float) * (n+1)))))
      error("Memory allocation failed for A in choldc.c \n");
  }
*/
  for(i = 0; i < n; i++) {
    P[i+1] = (float)(p[i]);
    X[i+1] = (float)(x[i]);
    B[i+1] = (float)(b[i]);
/*    for(j = 0; j <n ; j++)
      A[i+1][j+1] = (float)(a[i][j]);
*/
  }  

/*  float_cholsl(A, n, P, B, X); */
  float_cholsl(a, n, P, B, X);

  for(i = 0; i < n; i++) {
    p[i] = (float)(P[i+1]);
    b[i] = (float)(B[i+1]);
    x[i] = (float)(X[i+1]);
/*    for(j = 0; j <n ; j++)
      a[i][j] = (float)(A[i+1][j+1]); 
*/
  }  


/*  for(i = 0; i <= n; i++) {
    free(A[i]);
  }
*/
  free(P);
  free(B);
  free(X);
/*  free(A); */
}



/*****************************************************************************/

#define NRANSI
#define TINY 1.0e-20;

void ludcmp(a, n, indx, d)
float **a;
int n;
int *indx;
float *d;
{
  int i,imax,j,k;
  float big,dum,sum,temp;
  float *vv;
  
  if(!(vv=(float *) malloc( (n+1) * sizeof(float) )))
    error("Memory allocation failed for vv in choldc.c \n");
  
  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) printf("Singular matrix in routine ludcmp");
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
  free( vv );
}
#undef TINY
#undef NRANSI

/*****************************************************************************/

void lubksb(float **a, int n, int *indx, float b[])
{
  int i,ii=0,ip,j;
  float sum;
  
  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}
