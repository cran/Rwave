
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
     float *rho;
{
  int i;

  for(i=0; i< snakesize;i++){
    *rho /= (float)rate;
    /* *rho = floorf(*rho); */
    *rho = (float)floor((double)(*rho));
/*    printf("rho=%f\n",*rho);*/
    rho++;
  }
return;
}

void snakexpand(rho,rate,snakesize)
     int rate, snakesize;
     float *rho;
{
  int i;

  for(i=0; i< snakesize;i++){
    *rho *= (float)rate;
    /* *rho = floorf(*rho); */
    *rho = (float)floor((double)(*rho));
/*    printf("rho=%f\n",*rho);*/
    rho++;
  }
return;
}
      
