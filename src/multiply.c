#include <stdlib.h>

/****************************************************************
*               (c) Copyright  1997                             *
*                          by                                   *
*      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                  Princeton University                         *
*                  All right reserved                           *
****************************************************************/

#include "Swave.h"





/* Multiplication of complex-valued signals:
   -----------------------------------------
   Ri1,Ii1: Real and imaginary parts of signal 1
   Ri2,Ii2: Real and imaginary parts of signal 2
   Or,Oi: Real and imaginary parts of output signal 
*/

void multiply(Ri1,Ii1,Ri2,Ii2,Or,Oi,isize)
     double *Ri1, *Ri2, *Ii1, *Ii2, *Or, *Oi;
     int isize;
{
  int i;

  for(i = 0; i < isize; i++) {
    *Or = (*Ri1)*(*Ri2) - (*Ii1)*(*Ii2);
    *Oi = (*Ii1)*(*Ri2) + (*Ri1)*(*Ii2);
    Or++; Oi++;
    Ri1++; Ii1++;
    Ri2++; Ii2++;
  }
  //return;
}
