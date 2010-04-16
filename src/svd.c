#include <stdlib.h>


/* NUMERICAL RECIPE */

/****************************************************************/
/*              (c) Copyright  1997
/*                         by                                   */
/*     Author: Rene Carmona, Andrea Wang, Wen-Liang Hwang       */
/*                 Princeton University                         */
/*                 All right reserved                           */
/****************************************************************/

/* This routine computes the singular value decomposition of a matrix    */
/* a= u w v^t. The maxtix U is replaces a on output. The diagonal matrix */
/* of singular values w is output as a vector w[1..n]. The matrix v(not  */
/* the transpost v^t) is output as v[1..n][1..n]. m  must be greater or  */
/* equal to n; if it is smaller, than a should be filled up to square    */
/* with zero rows.                                                       */

#include "Swave.h"
#include "denoise.h"


#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
static double sqrarg;
#define SQR(a) ((sqrarg = (a)) == 0.0 ? 0.0 : sqrarg * sqrarg)
static double maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a), maxarg2=(b), (maxarg1) > (maxarg2) ? \
(maxarg1) : (maxarg2))
static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a), iminarg2=(b),(iminarg1) < (iminarg2) ? \
(iminarg1) : (iminarg2))

double pythag(double a, double b)
{
  double absa, absb;

  absa = fabs(a);
  absb = fabs(b);
  if(absa > absb) 
    return (absa * sqrt(1.0 + SQR(absb/absa)));
  else 
    return (absb == 0.0 ? 0.0 : absb * sqrt(1.0+SQR(absa/absb)));
}
  
void svdcmp(double **a, int m, int n, double *w, double **v)
{
  int flag, i, its, j, jj, k, l, nm;
  double anorm, c, f, g, h, s, scale,  x, y, z, *rv1;
  
  if(!(rv1 = (double *)R_alloc((n + 1), sizeof(double) )))
    error("Memory allocation failed for rv1 in svd.c \n");
  g = scale = anorm = 0.0;
  
  /* Householder reduction to bidiagonal form */

  for (i = 1; i <= n; i++) {
    l = i + 1;
    rv1[i] = scale * g;
    g = s = scale = 0.0;
    if(i <= m) {
      for(k=i;k <= m; k++) scale += fabs(a[k][i]);
      if(scale) {
	for(k=i; k <=m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i] * a[k][i];
	}
	f = a[i][i];
	g = -SIGN(sqrt(s),f);
	h = f * g - s;
	a[i][i] = f-g;
	for ( j= l; j <= n; j++) {
	  for(s= 0.0, k= i; k<=m; k++) s+= a[k][i] * a[k][j];
	  f = s/h;
	  for(k=i; k<= m;k++) a[k][j] += f * a[k][i];
	}
	for(k=i;k<=m;k++) a[k][i] *= scale;
      }
    }
    w[i] = scale * g;
    g = s = scale = 0.0;
    if(i<= m && i != n) {
      for (k = l; k <= n; k++) scale += fabs(a[i][k]);
      if(scale) {
	for(k= l; k<=n; k++) {
	  a[i][k] /= scale;
	  s += a[i][k] * a[i][k];
	}
	f = a[i][l];
	g = -SIGN(sqrt(s),f);
	h= f*g-s;
	a[i][l] = f-g;
	for(k = l; k <= n; k++) rv1[k] = a[i][k]/h;
	for ( j = l; j <= m; j++) {
	  for (s = 0.0, k= l; k <= n; k++) s += a[j][k] * a[i][k];
	  for (k = l; k <= n; k++) a[j][k] += s * rv1[k];
	}
	for(k = l; k <= n; k++) a[i][k] *= scale;
      }
    }
    anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
  }

/* Accumulation of right-hand transformations */

  for (i = n; i >= 1; i--) {
    if (i < n) {
      if(g) {
	for (j = l; j <= n; j++)
	  v[j][i] = (a[i][j]/a[i][l])/g;
	for(j = l; j <= n; j++) {
	  for (s = 0.0, k = l; k <= n; k++) s += a[i][k] * v[k][j];
	  for(k = l; k <= n; k++)  v[k][j] += s * v[k][i];
	}
      }
      for(j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g= rv1[i];
    l = i;
  }

/* Accumulation of left-hand transformations */

  for(i = IMIN(m,n); i >= 1; i--) {
    l = i + 1;
    g = w[i];
    for(j = l; j <= n; j++) a[i][j] = 0.0;
    if(g) {
      g = 1.0/g;
      for ( j = l; j <= n; j++) {
	for ( s= 0.0, k = l; k <= m; k++) s+= a[k][i] * a[k][j];
	f = (s/a[i][i]) * g;
	for( k = i; k <= m; k++) a[k][j] += f * a[k][i];
      }
      for(j = i; j <= m; j++) a[j][i] *= g;
    }
    else for (j = i; j <= m; j++) a[j][i] = 0.0;
    ++a[i][i];
  }

/* Diagonalization of bidiagonal form */

  for(k = n; k >= 1; k--)   {
    for(its= 1; its <= 30; its++) {
      flag = 1;
      for(l=k; l >= 1; l--) {
	nm = l-1;
	if((double)(fabs(rv1[l]) + anorm) == (double)anorm) {
	  flag = 0;
	  break;
	}
	if((double)(fabs(w[nm]) + anorm) == (double)anorm) break; 
      }
      if(flag) {
	c = 0.0;
	s = 1.0;
	for(i = l; i <= k; i++) {
	  f = s * rv1[i];
	  rv1[i] = c * rv1[i];
	  if ((double)(fabs(f) + anorm) == (double)anorm) break;
	  g = w[i];
	  h = pythag(f,g);
	  w[i] = h;
	  h = 1.0/h;
	  c = g * h;
	  s = -f * h;
	  for( j = 1; j <= m; j++) {
	    y = a[j][nm];
	    z = a[j][i];
	    a[j][nm] = y * c + z * s;
	    a[j][i] = z * c - y * s;
	  }
	}
      }
      z = w[k];
      if( l == k) {
	if ( z < 0.0) {
	  w[k] = -z;
	  for ( j = 1; j <= n; j++) v[j][k] = (-v[j][k]);
	}
	break;
      }
      if(its == 30) error("No convergence in 30 SVDCMP iterations");
      x = w[l];
      nm = k - 1;
      y = w[nm];
      g = rv1[nm];
      h = rv1[k];
      f = ((y-z) * (y + z) + (g - h) * (g + h))/(2.0*h*y);
      g = pythag(f,1.0);
      f= ((x-z) * (x+z) + h * ((y/(f + SIGN(g,f))) - h))/x;
      c=s = 1.0;
      for(j = l; j <= nm; j++) {
	i = j + 1;
	g = rv1[i];
	y = w[i];
	h = s * g;
	g = c * g;
	z = pythag(f,h);
	rv1[j] = z;
	c = f/z;
	s = h/z;
	f = x * c + g * s;
	g = g * c - x * s;
	h = y * s;
	y *= c;
	for ( jj= 1; jj <= n; jj++) {
	  x = v[jj][j];
	  z = v[jj][i];
	  v[jj][j] = x * c + z * s;
	  v[jj][i] = z * c - x * s;
	}
	z = pythag(f,h);
	w[j] = z;
	if(z) {
	  z = 1.0/z;
	  c = f * z;
	  s = h * z;
	}
	f = (c * g) + (s * y);
	x = (c * y) - (s * g);
	for(jj = 1; jj <= m; jj++) {
	  y = a[jj][j];
	  z = a[jj][i];
	  a[jj][j] = y * c + z * s;
	  a[jj][i] = z * c - y * s;
	}
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }

}


/**********************************************************************/
/* Solve AX = B for a vector X, where A is specified by the array     */
/* U, w and v as returned by svdcmp. b is the input right-hand side   */
/* x is the output solution vector. No input quantities are destroyed */
/* so the routin may be called sequentially with different b's.       */
/**********************************************************************/
    

void svbksb(double **U, double *W, double **V, int m, int n,
  double *B, double *X)
{
  int jj, j, i;
  double s, *tmp;
  
  if(!(tmp = (double *)R_alloc((n+1), sizeof(double) )))
    error("Memory allocation failed for tmp in svd.c \n");
  for (j = 1; j <= n; j++) {
    s = 0.0;
    if (W[j]) {
      for(i = 1; i <= m; i++) s += U[i][j] * B[i];
	s = s/W[j];
    }
    tmp[j] = s;
  }
  for(j=1; j <= n; j++) {
    s = 0.0;
    for(jj = 1; jj <= n; jj++) s += V[j][jj] * tmp[jj];
    X[j] = s;
  }
}
  
/* solve linear Ax = B by single value decomposition */
	   
void svdecomp_solve(double **a, double *b, double *x, int m,
  int n, double **w, double ***v)
{
  int i, j;
  double **A, *W, **V, *B, *X;
  void double_residue();

  if(!((*w) = (double *) R_alloc(n, sizeof(double) )))
    error("Memory allocation failed for (*w) in svd.c \n");
  if(!((*v) = (double **)R_alloc(  n, sizeof(double *))))
    error("Memory allocation failed for (*v) in svd.c \n");
  for(i = 0; i < n; i++) 
    if(!((*v)[i] = (double *) R_alloc(n, sizeof(double) )))
      error("Memory allocation failed for (*v)[] in svd.c \n");
  
  if(!(W = (double *)R_alloc((n+1), sizeof(double) )))
    error("Memory allocation failed for W in svd.c \n");
  if(!(V = (double **)R_alloc( (n+1), sizeof(double *))))
    error("Memory allocation failed for V in svd.c \n");
  if(!(A = (double **)R_alloc( (m+1), sizeof(double *))))
    error("Memory allocation failed for A in svd.c \n");
  if(!(B = (double *)R_alloc((m+1), sizeof(double) )))
    error("Memory allocation failed for B in svd.c \n");
  if(!(X = (double *)R_alloc((n+1), sizeof(double) )))
    error("Memory allocation failed for X in svd.c \n");
  for(i = 0; i <= n; i++) 
    if(!(V[i] = (double *)R_alloc((n+1), sizeof(double) )))
      error("Memory allocation failed for V[] in svd.c \n");
  for(i = 0; i <= m; i++) 
    if(!(A[i] = (double *)R_alloc((n+1), sizeof(double) )))
      error("Memory allocation failed for A[] in svd.c \n");
  
  for( i = 0; i < m; i++) {
    B[i+1] = (double)(b[i]);
    for(j = 0; j < n; j++) 
      A[i + 1][j + 1] = (double)(a[i][j]);
  }

  svdcmp(A,m,n,W,V); 
  svbksb(A,W,V,m,n,B,X);
  double_residue(A,W,V,m,n,B,X); 

  for( i = 0; i < m; i++) 
    for(j = 0; j < n; j++) 
      a[i][j] = (double)(A[i+1][j+1]);
  
  for( i = 0; i < n; i++)
    for(j = 0; j < n; j++) 
      (*v)[i][j] = (double)(V[i+1][j+1]);
  
  for(i = 0; i < n; i++) {
    (*w)[i] = (double)(W[i + 1]);
    x[i] = (double)(X[i+1]);
  }
  
/*  residue(a,*w,*v,m,n,b,x);  */
/*
  output_array(a,m,n,"U");
  output_array((*v),n,n,"V");
  output_signal(b,m,"B");
  output_signal(x,n,"X");
  output_signal((*w),n,"W");
*/
}
    
/* compute the L2 Norm */

void residue(double **a, double *w, double **v, int m, int n,
  double *b, double *x)
{
  double **tmp, *tmp1;
  double sum;
  int i, j, k;

  if(!(tmp = (double **)R_alloc( m, sizeof(double *))))
    error("Memory allocation failed for tmp in svd.c \n");
  if(!(tmp1 = (double *) R_alloc(m, sizeof(double) )))
    error("Memory allocation failed for tmp1 in svd.c \n");
  for(i = 0 ; i < m; i++)
    if(!(tmp[i] = (double *) R_alloc(n, sizeof(double) )))
      error("Memory allocation failed for tmp[] in svd.c \n");
  
  for(i = 0; i < m; i++)
    for(j = 0; j < n; j++) {
      tmp[i][j] = 0.0;
      for(k = 0; k < n; k++)
	  tmp[i][j] = tmp[i][j] + w[k] * a[i][k] * v[j][k];
    }


  for(i = 0; i < m; i++) {
    tmp1[i] = 0.0;
    for(k = 0; k < n; k++)
      tmp1[i] = tmp1[i] + tmp[i][k] * x[k];
  }
  
  for(i = 0; i < m; i++)
    tmp1[i] = tmp1[i] - b[i];
  
  sum = 0.0;

  for(i = 0; i < m; i++)
    sum = sum + tmp1[i] * tmp1[i];
  
/*  printf("Residule is %g \n",sum); */
}


void double_residue(double **a, double *w, double **v, int m,
  int n, double *b, double *x)
{
  double **tmp, *tmp1;
  double sum;
  int i, j, k;

  if(!(tmp = (double **)R_alloc( (m+1), sizeof(double *))))
    error("Memory allocation failed for tmp in svd.c \n");
  if(!(tmp1 = (double *)R_alloc((m+1), sizeof(double) )))
    error("Memory allocation failed for tmp1 in svd.c \n");
  for(i = 1 ; i <= m; i++)
    if(!(tmp[i] = (double *)R_alloc((n+1), sizeof(double) )))
      error("Memory allocation failed for tmp[] in svd.c \n");
  
  for(i = 1; i <= m; i++)
    for(j = 1; j <= n; j++) {
      tmp[i][j] = 0.0;
      for(k = 1; k <= n; k++)
	tmp[i][j] = tmp[i][j] + w[k] * a[i][k] * v[j][k];
    }


  for(i = 1; i <= m; i++) {
    tmp1[i] = 0.0;
    for(k = 1; k <= n; k++)
      tmp1[i] = tmp1[i] + tmp[i][k] * x[k];
  }
  
  for(i = 1; i <= m; i++)
    tmp1[i] = tmp1[i] - b[i];
  
  sum = 0.0;

  for(i = 1; i <= m; i++)
    sum = sum + tmp1[i] * tmp1[i];
  
/*   printf("Residule is %g \n",sum); */
}
    
/* compute the L2 Norm */
/*
void Sresidue(double *a, double *w, double *v, int m, int n,
  double *b, double *x)
{
  double *tmp, *tmp1;
  double sum;
  int i, j, k;
  int t;

  if(!(tmp = (double *)malloc(sizeof(double) * m * n)))
    error("Memory allocation failed for tmp in svd.c \n");
  if(!(tmp1 = (double *)malloc(sizeof(double) * m)))
    error("Memory allocation failed for tmp1 in svd.c \n");
  
  for(i = 0; i < m; i++)
    for(j = 0; j < n; j++) {
      t = i * n + j;
      tmp[t] = 0.0;
      for(k = 0; k < n; k++)
	tmp[t] = tmp[t] + w[k] * a[i*n+k] * v[j*n+k];
    }


  for(i = 0; i < m; i++) {
    tmp1[i] = 0.0;
    for(k = 0; k < n; k++)
      tmp1[i] = tmp1[i] + tmp[i*n+k] * x[k];
  }
  
  for(i = 0; i < m; i++)
    tmp1[i] = tmp1[i] - b[i];
  
  sum = 0.0;

  for(i = 0; i < m; i++)
    sum = sum + tmp1[i] * tmp1[i];
  
  printf("Residule is %f \n",sum); 
  free(tmp);
  free(tmp1);
}
*/
/* Called by Splus */
	   
void Ssvdecomp(double *a, int *pm, int *pn, double *w, double *v)
{
  int m, n;
  int i, j, k;
  double **A;
  double **V;
  double *W;

  m = *pm;
  n = *pn;

  if(!(A = (double **)R_alloc( (m+1), sizeof(double *))))
    error("Memory allocation failed for A in svd.c \n");
  if(!(V = (double **)R_alloc( (n+1), sizeof(double *))))
    error("Memory allocation failed for V in svd.c \n");
  if(!(W = (double *)R_alloc((n+1), sizeof(double) )))
    error("Memory allocation failed for W in svd.c \n");

  for(i = 0; i <= m; i++) 
    if(!(A[i] = (double *)R_alloc((n+1), sizeof(double) )))
      error("Memory allocation failed for A[] in svd.c \n");

  for(i = 0; i <= n; i++) 
    if(!(V[i] = (double *)R_alloc((n+1), sizeof(double) )))
      error("Memory allocation failed for V[] in svd.c \n");

  for(j = 0; j < n; j++)   
    for( i = 0; i < m; i++)
      A[i + 1][j + 1] = (double)(a[j * m + i]);

  svdcmp(A,m,n,W,V);  

  for(i = 0,k=0; i < n; i++) 
    for( j = 0; j < m; j++,k++)
      a[k] = A[j+1][i+1];

  for( i = 0; i < n; i++)
    w[i] = W[i+1];

  for( i = 0, k = 0; i < n; i++)
    for(j = 0; j < n; j++, k++) 
      v[k] = V[j+1][i+1];

}
    

