
/*******************************************************************/
/*              (c) Copyright  1997                                */
/*                         by                                      */
/*  Author: Rene Carmona, Bruno Torresani, Wen L. Hwang, A. Wang   */
/*                 Princeton University            e               */
/*                 All right reserved                              */
/*******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dau_wave.h"
#include "Swave.h"


/****************************************************************
*  Function: inverse_wavelet_transform:
*  ------------------------------------
*  Computation of the inversion of dyadic wavelet transform using spline
*  wavelet. This commands corresponds to the command dw, which computes
*  dyadic wavelet transform without subsampling at each resolution.
*
*    f_back: reconstructed signal
*    Sf: multiresolution signals
*    Wf: wavelet tranform of signal
*    max_resoln: number of decomposition
*    np: signal size
*    filtername: reconstruction filter
*
****************************************************************/


void inverse_wavelet_transform(f_back,Sf,Wf,max_resoln,np,filtername)
     float *f_back, *Sf, *Wf;
     int max_resoln, np;
     char *filtername;
{
  int i, j, n, k;
  float sum;
  bound *K_bound, *S_bound;
  float **S, **K;
  int offset;
  float *tmp;

  if(!(tmp = (float *) malloc(np * sizeof(float))))
    error("Memory allocation failed for tmp in signal_back.c \n");

  KSfilter_bound(filtername,&K_bound,&S_bound,max_resoln);
  Sfilter_compute(filtername,&S,S_bound,max_resoln);
  Kfilter_compute(filtername,&K,K_bound,max_resoln);

  for( i = 0; i < np; i++)
    f_back[i] = Sf[i];


  for ( j = max_resoln; j >=  1 ; j-- ){

    for ( n = 0; n < np; n++ )    {
      sum = 0.0;
      for ( k = S_bound[j-1].lb; k <= S_bound[j-1].ub; k++ ){
	sum += S[j-1][k-S_bound[j-1].lb] * f_back[(n-k+np)%np];	
      }
      tmp[n] = sum;
    }
    offset = (j - 1) * np;
    for ( n = 0; n < np; n++ )    {
      sum = 0.0;
      for ( k = K_bound[j-1].lb; k <= K_bound[j-1].ub; k++ ) {
	sum += K[j-1][k-K_bound[j-1].lb] * Wf[offset+(n-k+np)%np];
      }
      tmp[n] = tmp[n] + sum;
    }
    signal_copy(tmp,f_back,0,np);
  }

  for ( j = 0; j <= max_resoln; j++ )  {
    free( S[j] );
    free( K[j] );
  }
  free( S );
  free( S_bound );
  free( K );
  free( K_bound );
  free(tmp);

}

/****************************************************************
*  Function: Sinverse_wavelet_transform:
*  -------------------------------------
*  Called by S plus. 
*  Computation of the inversion of dyadic wavelet transform using spline
*  wavelet. This command is corresponding to the command dw, which computes
*  dyadic wavelet transform without subsampling at each resolution.
*
*    f_back: reconstructed signal
*    Sf: multiresolution signals
*    Wf: wavelet tranform of signal
*    max_resoln: number of decomposition
*    np: signal size
*    filtername: reconstruction filter
*
****************************************************************/

Sinverse_wavelet_transform(f_back,Sf,Wf,pmax_resoln,pnp,pfiltername)
     float *f_back, *Sf, *Wf;
     int *pmax_resoln, *pnp;
     char **pfiltername;
{
  int np, max_resoln;
  char *filtername;

  filtername = *pfiltername;
  np = *pnp;
  max_resoln = *pmax_resoln;

  inverse_wavelet_transform(f_back,Sf,Wf,max_resoln,np,filtername); 
}

