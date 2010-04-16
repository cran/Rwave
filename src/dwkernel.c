#include <stdlib.h>


/*******************************************************************/
/*              (c) Copyright  1997                                */
/*                         by                                      */
/*  Author: Rene Carmona, Bruno Torresani, Wen L. Hwang, A. Wang   */
/*                 Princeton University                            */
/*                 All right reserved                              */
/*******************************************************************/

#include "Swave.h"
#include "dyadic.h"

/****************************************************************
*  Function: fexp2:
*  ----------------
*  returning a double number of power of 2^j
*
*    j: exponent
*
****************************************************************/

/* 2^j  (j can be >= 0 or < 0 ) */
double fexp2(j)
  int j;
{
  int k;
  double s;
  s = 1.0;
  if (j >= 0) {
    return( (double)(1 << j));
  }
  else {
    for (k = j; k < 0 ; k++)
      s /= 2.0;
    return(s);
  }
}

/****************************************************************
*  Function: wavelet_transform_gradient:
*  -------------------------------------
*  Computation of the derivative of the wavelet transform along spatial 
*
*    grad: derivative of wavelet transform
*    s: wavelet transform
*    max_resoln: number of decomposition
*    np: signal size
*
****************************************************************/

void wavelet_transform_gradient( grad, s, max_resoln, np )
     double **grad;
     double **s;
     int max_resoln;
     int np;
{
  int j,t;
  int np_minus1 = np - 1;

  for ( j = 1; j <= max_resoln; j++ )  {
    for ( t = 0 ; t < np_minus1; t++)
      grad[j][t] = s[j][t+1] - s[j][t];
    grad[j][t] = 0.0;
  }
}

/****************************************************************
*  Function: signal_K_compute:
*  ---------------------------
*  Computation of kernel
*
*    K: kernel
*    W: wavelet transform
*    max_resoln: number of decomposition
*    np: signal size
*
****************************************************************/


void signal_K_compute(K,W,max_resoln,np )
     double ***K;
     double **W;
     int max_resoln;
     int np;
{
  int j, z, y, x, t, i, offset;
  double sum;
  double **grad_W, *k_tilda;
  double fexp2();

  if(!(grad_W = (double **) R_alloc( (max_resoln+1) , sizeof(double *) )))
    error("Memory allocation failed for grad_pis in K_compute.c \n");
  if(!(k_tilda = (double *) R_alloc( np , sizeof(double) )))
    error("Memory allocation failed for k_tilda in K_compute.c \n");
  
  for(i = 1; i <= max_resoln; i++)
    if(!(grad_W[i] = (double *)R_alloc(np , sizeof(double))))
      error("Memory allocation failed for grad_W[] in K_compute.c \n");


  wavelet_transform_gradient( grad_W, W, max_resoln, np );

  for ( z = 0; z < np; z++ )  {
    for ( j = 1, sum = 0.0; j <= max_resoln; j++ )    {
      for ( x = 0; x < np; x++ )      {
        y = (x+z)%np;
	sum += W[j][x] * W[j][y] +
	  fexp2(j) * grad_W[j][x] * grad_W[j][y];
      }
    }
    k_tilda[z] = sum;
  }

/*  output_signal( k_tilda, np, "signal_k_tilda");  */

  /**************************************************************/
  /* we compute K in the following ... k is two dimensional and */
  /* k_tilda is one dimensional ...                             */
  /**************************************************************/

  if(!((*K)=(double **) R_alloc( (np+1) , sizeof(double *))))
    error("Memory allocation failed for *k in K_compute.c \n");
  for ( t = 0; t <= np; t++ )
    if(!((*K)[t] = (double *) R_alloc( (np+1) , sizeof(double) )))
      error("Memory allocation failed for (*k)[] in K_compute.c \n");
  
  for ( t = 0; t < np; t++ )    {
    for ( i = t, z = 0; i < np; i++, z++ )
      (*K)[t+1][i+1] = (*K)[i+1][t+1] = k_tilda[z];
  }
  
/*  
  output_array( *K, np,np,"signal_K_matrix" );
*/
  //for(i = 0; i <= max_resoln; i++)
//    free( grad_W[i]);
  
}

/****************************************************************
*  Function: signal_tilda_adjust:
*  ------------------------------
*  Adjustment of the size of (w)k_tilda for a signal (w)k_tilda read
*  from disk ...
*
*    tilda: (w)k_tilda
*    ksize: size of (w)k_tilda required
*    fname: name of (w)k_tilda is stored
*    fsize: size of (w)k_tila in file
*
****************************************************************/

/* please don't write to disk
void signal_tilda_adjust(tilda,ksize,fname,fsize)
     double **tilda;
     char *fname;
     int fsize, ksize;
{
  double *tmp;
  int i, j, k, l, m, n;
  int rsize, middle, rest;
  
  if(!(*tilda = (double *)malloc(sizeof(double) * ksize)))
    error("Memory allocation failed for *tilda in K_op.c \n");
  signal_zero(*tilda, ksize);

  input_signal(fname, &tmp, fsize);
  rsize = min(ksize, fsize)/2;

  for(i= 0; i < rsize; i++) 
    (*tilda)[i] = tmp[i];

  for(i= 1; i <= rsize; i++) {
    k = fsize -  i;
    m = ksize -  i;
    (*tilda)[m] = tmp[k];
  }
  free(tmp);
}

*/ 






