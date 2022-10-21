#include <stdlib.h>

/***************************************************************
*    $Log: cwt_maxima.c,v $                                    *
****************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
***************************************************************/


#define MAX(x,y) ((x)>(y) ? (x) : (y))
#include "Swave.h"
#include "denoise.h"



/****************************************************************
*  function Scwt_gmax
*    compute the global maximum of cwt for fixed position
*
*  input: cwt
*  output: cwt global maxima at fixed b
*  nrow,ncol: parameters of the cwt 
*****************************************************************/

void Scwt_gmax(double *input, double *output,  int *pnrow,  int *pncol,  int *posvector)
{
  int nrow, ncol, i, j;
  int pos;
  double tmp;

  nrow = *pnrow;
  ncol = *pncol;

  for(i = 0; i < nrow; i++) {
    tmp = -99999999.0;
    pos = -1;
    for(j = 0; j < ncol; j++) {
      tmp = MAX(tmp, input[j * nrow + i]);
      if(tmp == input[j * nrow + i]) pos = j;
    }
    posvector[i] = pos;
    output[pos * nrow + i] = tmp;
  }
}


/****************************************************************
*  function Scwt_mridge
*    compute the local maxima of cwt for fixed position
*
*  input: cwt
*  output: cwt global maxima at fixed b
*  nrow,ncol: parameters of the cwt 
*****************************************************************/

void Scwt_mridge(double *input, double *output, int *pnrow, int *pncol)
{
  int nrow, ncol, i, j;
  

  nrow = *pnrow;
  ncol = *pncol;

  for(i = 0; i < nrow; i++) {
    if(input[i] > input[nrow + i]) 
      output[i] = input[i];
    if(input[(ncol-1) * nrow + i] > input[(ncol-2) * nrow + 1])
      output[(ncol-1) * nrow + i] = input[(ncol-1) * nrow + i];

    for(j = 1; j < ncol-1; j++) {
      if(((input[j * nrow + i] > input[(j+1) * nrow + i]) &&
	 (input[j * nrow + i] >= input[(j-1) * nrow + i])) ||
	 ((input[j * nrow + i] > input[(j-1) * nrow + i]) &&
	  (input[j * nrow + i] >= input[(j+1) * nrow + i])))
	output[j * nrow + i] = input[j * nrow + i];
    }
  }
}

	

