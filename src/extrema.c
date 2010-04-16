#include <stdlib.h>


/******************************************************************
*              (c) Copyright  1997                                *
*                         by                                      *
*  Author: Rene Carmona, Bruno Torresani, Wen L. Hwang, A. Wang   *
*                 Princeton University                            *
*                 All right reserved                              *
*******************************************************************/


#include "Swave.h"
#include "dyadic.h"

/****************************************************************
*  Function: modulus_maxima:
*  -------------------------
*  Computation of modulus local maxima of wavelet transform
*
*    extrema: modulus local maxima of wavelet transform
*    wt: wavelet transform
*    resoln_ptr: number of decomposition
*    np_ptr: the signal size
*
****************************************************************/


void modulus_maxima(double *extrema, double *wt, int *resoln_ptr,
  int *np_ptr )
{
  int resoln = *resoln_ptr;
  int np = *np_ptr;
  double *abs;

  int x, j;
  int offset;


  if(!(abs  = (double *) R_alloc( np , sizeof(double) )))
    error("Memory allocation failed for abs in extrema.c");

  for (j = 0; j < resoln; j++)  {
    offset = j * np;
    for (x = 0; x < np; x++)
      abs[x] = fabs( (double) wt[offset+x] );

    extrema[offset] = 0.0;
    extrema[offset+np-1] = 0.0;

    for ( x = 1; x < (np-1); x++ )    {
      if (((abs[x] > abs[x-1]) && (abs[x] >= abs[x+1])) ||
	  ((abs[x] > abs[x+1]) && (abs[x] >= abs[x-1])))
	extrema[offset+x] = wt[offset+x];
      else
	extrema[offset+x] = 0.0;
    }
  }
}

/****************************************************************
*  Function: extrema_input:
*  ------------------------
*  Converting extrema representation from array image_ext structure
*
*    extrema: modulus local maxima of wavelet transform
*    max_resoln: number of decomposition
*    np: signal size
*    ext: structure of image_ext
*    num_of_extrema: number of extrema
*
****************************************************************/

void extrema_input(double *extrema, int max_resoln, int np,
  image_ext **ext, int *num_of_extrema)
{
  int j, k, t, length, offset;
  
  length = max_resoln * np;

  *num_of_extrema = 0;
  for ( t = 0; t < length; t++ )
    if ( extrema[t] != 0.0 )
      (*num_of_extrema)++;

  if(!(*ext = (image_ext *) R_alloc( (*num_of_extrema) , sizeof(image_ext) )))
    error("Memory allocation failed for *ext in point_input.c \n");

  k = 0; 
  for ( j = 1; j <= max_resoln; j++ )  {
    offset = (j-1) * np;
    for ( t = 0; t < np; t++ )
      if ( extrema[offset+t] != 0.0 ) { /* detect one extremum */
        (*ext)[k].resoln= j;
        (*ext)[k].x = t;
        (*ext)[k].W1f = extrema[offset+t];
        k++;
      }
  }   
}



