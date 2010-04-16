#include <stdlib.h>


/***************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
***************************************************************/

#include "Swave.h"
#include "denoise.h"

#define PRECISION 1.e-16


/***************************************************************
*  Function: Lpnorm
*  ---------
*     L^p norm of a matrix
*       Computes the L^p norm of a (double) complex valued matrix.
*
*   norm: L^p norm
*   Rmat, Imat: real and imag. parts of the matrix
*   p: exponent for the L^p norm
*   length: number of rows
*   width: number of columns
***************************************************************/
void Lpnorm(double *norm, double *p, double *Rmat, double *Imat,
  int *length, int *width)
{
  int i,j;
  double tmp, rtmp, itmp;
  double ntmp = 0.0;

  for(i=0;i<(*length);i++){
    for(j=0;j<(*width);j++){
      rtmp = fabs(*Rmat);
      itmp = fabs(*Imat);
      if ((rtmp >= PRECISION)&&(itmp >= PRECISION)){
	tmp = pow(rtmp,*p) + pow(itmp,*p);
	ntmp += tmp;
      }
      Rmat++; Imat++;
    }
  }
  *norm = pow(ntmp, 1/(*p));
  return;
}


/***************************************************************
*  Function: entropy
*  ---------
*     Entropy of a matrix
*       Computes the entropy of a (double) complex valued matrix.
*
*   entr: entropy
*   Rmat, Imat: real and imag. parts of the matrix
*   length: number of rows
*   width: number of columns
***************************************************************/
void entropy(double *entr, double *Rmat, double *Imat,
  int *length, int *width)
{
  int i,j;
  double tmp, ntmp=0.0;

  for(i=0;i<(*length);i++){
    for(j=0;j<(*width);j++){
      tmp = (*Rmat)*(*Rmat) + (*Imat)*(*Imat);
      if((tmp >= PRECISION))
	ntmp -= tmp * log(tmp);
      Rmat++; Imat++;
    }
  }
  *entr = ntmp;
  return;
}
  
