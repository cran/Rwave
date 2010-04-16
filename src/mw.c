#include <stdlib.h>


/*******************************************************************
*              (c) Copyright  1997                                 *
*                         by                                       *
*  Author: Rene Carmona, Bruno Torresani, Wen L. Hwang, A. Wang    *
*                 Princeton University                             *
*                 All right reserved                               *
*******************************************************************/



#include "Swave.h"
#include "dyadic.h"




/****************************************************************
*  Function: Sf_compute:
*  ---------------------
*  Computation of signals from resolution 1 up to resolution
*    2^max_resoln
*
*    Sf: resultant signals
*    f: the original signal
*    max_resoln: number of decomposition
*    np: the signal size
*    filtername: decomposition filter
*
****************************************************************/

void Sf_compute(double *Sf, double *f, int *max_resoln_ptr,
  int *np_ptr, char *filtername)
{
  int max_resoln = *max_resoln_ptr;
  int np = *np_ptr;
  int j, n, k, offset;
  bound *H_bound, *G_bound;
  double **H, sum;
  //char *filtername;


  //filtername = *pfiltername;

  HGfilter_bound(filtername,&H_bound,&G_bound,max_resoln);
  Hfilter_compute(filtername,&H,H_bound,max_resoln);

  for ( j = 0; j <= max_resoln; j++ )  {
    if ( j == 0 )     { /* Sf[0] = original signal f  */
      for ( n = 0; n < np; n++ )
        Sf[n] = f[n];
    }
    else    {
      offset = (j-1)*np;
      for ( n = 0; n < np; n++ )      {
        for ( k = H_bound[j-1].lb, sum = 0.0; k <= H_bound[j-1].ub; k++ ) 
	  sum += H[j-1][k-H_bound[j-1].lb] * Sf[offset+(n-k+np)%np];
        Sf[j*np+n] = sum;
      }
    }
  }
  
}




/****************************************************************
*  Function: Wf_compute:
*  ---------------------
*  Computation of the wavelet transform of a signal
*
*    Wf: the wavelet transform of a signal
*    Sf: the multiresolution representation of the original signal
*    max_resoln: number of decomposition
*    np: the signal size
*    filtername: decomposition filter
*
****************************************************************/

void Wf_compute(double *Wf, double *Sf, int *max_resoln_ptr,
  int *np_ptr, char *filtername)
{
  int max_resoln = *max_resoln_ptr;
  int np = *np_ptr;
  int j, n, k, offset;
  //char *filtername;
  bound *G_bound, *H_bound;
  double **G, sum;

  //filtername = *pfiltername;
  HGfilter_bound(filtername,&H_bound, &G_bound, max_resoln );
  Gfilter_compute(filtername,&G, G_bound, max_resoln );
  	
  for ( j = 1; j <= max_resoln; j++ )  {
    offset = (j-1)*np;
    for ( n = 0; n < np; n++ )    {
      for ( k = G_bound[j-1].lb, sum = 0.0; k <= G_bound[j-1].ub; k++ )
	/* Compute Wf[m] from Sf[m-1] */
        sum += G[j-1][k-G_bound[j-1].lb] * Sf[offset+(n-k+np)%np]; 
      Wf[offset+n] = sum;							
    }					       
  }
}








