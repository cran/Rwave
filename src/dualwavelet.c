
/******************************************************************
*              (c) Copyright  1997                                *
*                         by                                      *
*  Author: Rene Carmona, Bruno Torresani, Wen L. Hwang, A. Wang   *
*                 Princeton University                            *
*                 All right reserved                              *
******************************************************************/

#include "Swave.h"
#include "dyadic.h"

/****************************************************************
*  Function: signal_W_tilda:
*  -------------------------
*  Computation of W tilda
*
*    W_tilda: the dual wavelet 
*    W: wavelet 
*    K: kernel
*    max_resoln: number of decomposition
*    np: signal size
*
****************************************************************/


void signal_W_tilda(float ***W_tilda, float **W, float **K,
		    int max_resoln, int np)
{
  float *p, *b;
  int t, i, j;
  char filename[STRING_SIZE];
  
  if(!(p = (float *) malloc( np * sizeof(float) )))
    error("Memory allocation failed for p in image_W_tilda \n");
  if(!(b = (float *) malloc( np * sizeof(float) )))
    error("Memory allocation failed for b in image_W_tilda \n");
  if(!(*W_tilda = (float **) malloc( (max_resoln+1) * sizeof(float *) )))
    error("Memory allocation failed for *W_tilda in image_W_tilda \n");

  for(j = 1; j <= max_resoln; j++) {
    if(!((*W_tilda)[j] = (float *) malloc( np * sizeof(float) )))
      error("Memory allocation failed for (*W_tilda)[] in image_W_tilda \n");
  }

  for ( j = 1; j <= max_resoln; j++ )    {
    /*  printf("computing W_tilda[%d]\n", j ); */
    for ( t = 0; t < np; t++)
      b[t] = W[j][t];
    choldc(K, np, p );
    cholsl(K, np, p, b, (*W_tilda)[j]);
    
    filename_given(filename,"sig_W_tilda");
    filename_inc(filename,j);
    output_signal((*W_tilda)[j], np, filename);

  }
  free( p );
  free( b );
  for(j = 0; j <= np; j++)
    free(K[j]);
  free(K);

  return;
}


/****************************************************************
*  Function: signal_W_tilda_input:
*  ------------------------------
*  Read W tilda from disk
*
*    W_tilda: dual wavelet
*    max_resoln: number of decomposition
*    np: signal size
*
****************************************************************/

void signal_W_tilda_input(float ***W_tilda, int max_resoln, int np)
{
  int j;
  char filename[STRING_SIZE];

  if(!(*W_tilda = (float **) malloc( (max_resoln+1) * sizeof(float *) )))
    error("Memory allocation failed for *W_tilda in signal_W_tilda \n");


  for(j = 1; j <= max_resoln; j++) {
    filename_given(filename, "signal_W_tilda");
    filename_inc(filename, j);
    signal_tilda_adjust(&((*W_tilda)[j]), np, filename, 4096);

    filename_given(filename, "W_tilda");
    filename_inc(filename, j);
    output_signal((*W_tilda)[j], np, filename);
  }
}





