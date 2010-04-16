#include <stdlib.h>

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

static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
        (dmaxarg1) : (dmaxarg2))

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static double minarg1,minarg2;
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
double *vector(long nl, long nh);
int *ivector(long nl, long nh);
unsigned char *cvector(long nl, long nh);
unsigned long *lvector(long nl, long nh);
double *dvector(long nl, long nh);
double **matrix(long nrl, long nrh, long ncl, long nch);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
int **imatrix(long nrl, long nrh, long ncl, long nch);
double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl);
double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch);
double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_vector(double *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_cvector(unsigned char *v, long nl, long nh);
void free_lvector(unsigned long *v, long nl, long nh);
void free_dvector(double *v, long nl, long nh);
void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch);
void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh);

#else /* ANSI */
/* traditional - K&R */

void nrerror();
double *vector();
double **matrix();
double **submatrix();
double **convert_matrix();
double ***f3tensor();
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

void double_choldc(double **a, int n, double p[])
{
  int i,j,k;
  double sum;
  
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

void choldc(double **a, int n, double p[])
{
  double** A;
  double *P;
  int i,j,k;

/*  if(!(A = (double **)(malloc(sizeof(double *) * (n+1)))))
    error("Memory allocation failed for A in choldc.c \n");
*/
  if(!(P = (double *)(R_alloc((n+1), sizeof(double) ))))
    error("Memory allocation failed for P in choldc.c \n");
/*		  
  for(i = 0; i <= n; i++) {
    if(!(A[i] = (double *)(malloc(sizeof(double) * (n+1)))))
      error("Memory allocation failed for A in choldc.c \n");
  }
*/  
  for(i = 0; i < n; i++) {
    P[i+1] = (double)(p[i]);
/*    for(j = 0; j <n ; j++)
      A[i+1][j+1] = (double)(a[i][j]); */
  }
    
/*  double_choldc(A,n,P); */
  double_choldc(a,n,P);


  for(i = 0; i < n; i++) {
    p[i] = (double)(P[i+1]);
/*    for(j = 0; j <n ; j++)
      a[i][j] = (double)(A[i+1][j+1]); */
  }
  
}



/***************************************/
/* cholsl with C language array format */
/***************************************/
void double_cholsl(double **a, int n, double p[], double b[], double x[])
{
  int i,k;
  double sum;
  
  for (i=1;i<=n;i++) {
    for (sum=b[i],k=i-1;k>=1;k--) sum -= a[i][k]*x[k];
    x[i]=sum/p[i];
  }
  for (i=n;i>=1;i--) {
    for (sum=x[i],k=i+1;k<=n;k++) sum -= a[k][i]*x[k];
    x[i]=sum/p[i];
  }
}

void cholsl(double **a, int n, double p[], double b[], double x[])
{
  double** A;
  double *P;
  int i,j,k;
  double *B, *X;
  
/*  if(!(A = (double **)(malloc(sizeof(double *) * (n+1)))))
    error("Memory allocation failed for A in choldc.c \n");
*/
  if(!(P = (double *)(R_alloc((n+1), sizeof(double) ))))
    error("Memory allocation failed for P in choldc.c \n");
  if(!(B = (double *)(R_alloc((n+1), sizeof(double) ))))
    error("Memory allocation failed for B in choldc.c \n");
  if(!(X = (double *)(R_alloc((n+1), sizeof(double) ))))
    error("Memory allocation failed for X in choldc.c \n");
		  
/*  for(i = 0; i <= n; i++) {
    if(!(A[i] = (double *)(malloc(sizeof(double) * (n+1)))))
      error("Memory allocation failed for A in choldc.c \n");
  }
*/
  for(i = 0; i < n; i++) {
    P[i+1] = (double)(p[i]);
    X[i+1] = (double)(x[i]);
    B[i+1] = (double)(b[i]);
/*    for(j = 0; j <n ; j++)
      A[i+1][j+1] = (double)(a[i][j]);
*/
  }  

/*  double_cholsl(A, n, P, B, X); */
  double_cholsl(a, n, P, B, X);

  for(i = 0; i < n; i++) {
    p[i] = (double)(P[i+1]);
    b[i] = (double)(B[i+1]);
    x[i] = (double)(X[i+1]);
/*    for(j = 0; j <n ; j++)
      a[i][j] = (double)(A[i+1][j+1]); 
*/
  }  


}



/*****************************************************************************/

#define NRANSI
#define TINY 1.0e-20;

void ludcmp(a, n, indx, d)
double **a;
int n;
int *indx;
double *d;
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;
  
  if(!(vv=(double *) R_alloc( (n+1) ,  sizeof(double))))
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
}
#undef TINY
#undef NRANSI

/*****************************************************************************/

void lubksb(double **a, int n, int *indx, double b[])
{
  int i,ii=0,ip,j;
  double sum;
  
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
