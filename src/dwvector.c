
/*******************************************************************/
/*              (c) Copyright  1997                                */
/*                         by                                      */
/*  Author: Rene Carmona, Bruno Torresani, Wen L. Hwang, A. Wang   */
/*                 Princeton University                            */
/*                 All right reserved                              */
/*******************************************************************/

#include "Swave.h"

/****************************************************************
*  Function: compute_convolution:
*  ------------------------------
*  Computation of convolution product
*
*    s: resultant of convolution 
*    f: signal 1  
*    g: signal 2
*    np: signal size
*
****************************************************************/

/* convolution product s[m] = (f * g)[m] */


void compute_convolution( s, f, g, np )
     float *s, *f, *g;
     int np;
{
  int m, n;
  float sum;

  for ( m = 0; m < np; m++ )  {
    for ( n = 0, sum = 0.0; n < np; n++ )
      sum += f[(m-n+np)%np] * g[n];
    s[m] = sum;
  }
}

/****************************************************************
*  Function: product:
*  ------------------
*  Pointwise product of two real arrays
*
*    image: resultant array
*    f: array 1  
*    g: array 2
*    np: array size
*
****************************************************************/

void product( image, f, g, np )
float ***image, *f, *g;
int np;
{
  int x, y;

  if(!(*image = (float **) malloc( np * sizeof(float *) )))
    error("Memory allocation failed for *image in vector_op.c \n");

  for ( x = 0; x < np; x++ )
  {
    if(!((*image)[x] = (float *) malloc( np * sizeof(float) )))
      error("Memory allocation failed for *image in vector_op.c \n");
    for ( y = 0; y < np; y++ )
      (*image)[x][y] = f[x] * g[y];
  }
}

/****************************************************************
*  Function: complex_product:
*  --------------------------
*  Pointwise product of two complex vectors
*
*    product: resultant vector
*    s1: vector 1 (length of 2*np)  
*    s2: vector 2
*    np: vector size
*
****************************************************************/


void complex_product( product, s1, s2, np )
     float *product;
     float *s1;  /* length of 2*np */
     float *s2;  /* length of 2*np */
     int np;
{
  int m, x, y;
  float a, b, c, d;

  for ( m = 0; m < np; m++ )  {
    x = 2*m;
    y = 2*m+1;

    a = s1[x];     /* (a + bi) * (c + di) */
    b = s1[y];
    c = s2[x];
    d = s2[y];

    product[x] = a*c - b*d;
    product[y] = b*c + a*d;
  }
}


/****************************************************************
*  Function: signal_copy:
*  ----------------------
*  Copy signal
*
*    input: signal to be copied
*    output: resultant signal
*    offset: copy starts at ...
*    np: number of elements to be copied
*
****************************************************************/


signal_copy(input,output,offset,size)
     float *input;
     float *output;
     int size, offset;
{
  int i;

  for(i = 0; i < size; i++)
    output[i] = input[offset+i];
}

/****************************************************************
*  Function: signal_zero:
*  ----------------------
*  zero a signal
*
*    input: signal
*    size: number of signals to be set to zero
*
****************************************************************/

signal_zero(input,size)
     float *input;
     int size;
{
  int i;

  for(i = 0; i < size; i++)
    input[i] = 0.0;
}



