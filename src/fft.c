#include <stdlib.h>

/****************************************************************
*               (c) Copyright  1997                             *
*                          by                                   *
*      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                  Princeton University                         *
*                  All right reserved                           *
****************************************************************/

#include "Swave.h"
#include "denoise.h"

/* Fast Fourier transform (from numerical recipes routine)
   -------------------------------------------------------
*/

void double_fft(double *Or,double *Oi,double *Ir,double *Ii,
  int isize,int isign)
{
  double *tmp;
  int nt, find2power(int  n), newsize, i;

  nt = find2power(isize);
  newsize = 1 << nt;

  if(!(tmp = (double *)R_alloc(newsize * 2, sizeof(double))))
     Rf_error("Memory allocation failed for tmp in cwt_morlet.c \n");

  for(i = 0; i < isize; i++) {
    tmp[2 * i] = Ir[i];
    tmp[2 * i + 1] = Ii[i];
  }
  four1(tmp-1,newsize,isign);  

  
  for(i = 0; i < isize; i++) {
    if(isign == -1) {
      Or[i] = tmp[2 * i]/newsize;
      Oi[i] = tmp[2 * i + 1]/newsize;
    } 
    else {
      Or[i] = tmp[2 * i];
      Oi[i] = tmp[2 * i + 1];
    }
  }
}
