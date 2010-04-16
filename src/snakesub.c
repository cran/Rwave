#include <stdlib.h>


/***************************************************************
*    $Log: snakesub2.c,v $                                     *
****************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
****************************************************************
*	 Reduction of the size of the snake                    *
***************************************************************/


#include "Swave.h"

void snakesub(rho,rate,snakesize)
     int rate, snakesize;
     double *rho;
{
  int i;

  for(i=0; i< snakesize;i++){
    *rho /= (double)rate;
    /* *rho = floorf(*rho); */
    *rho = (double)floor((double)(*rho));
/*    printf("rho=%f\n",*rho);*/
    rho++;
  }
return;
}

void snakexpand(rho,rate,snakesize)
     int rate, snakesize;
     double *rho;
{
  int i;

  for(i=0; i< snakesize;i++){
    *rho *= (double)rate;
    /* *rho = floorf(*rho); */
    *rho = (double)floor((double)(*rho));
/*    printf("rho=%f\n",*rho);*/
    rho++;
  }
return;
}
      
