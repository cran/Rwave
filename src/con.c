


/***************************************************************
 *    $Log: con.c,v $                                           *
 ****************************************************************
 *              (c) Copyright  1997                             *
 *                         by                                   *
 *     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
 *                 Princeton University                         *
 *                 All right reserved                           *
 ***************************************************************/




#include <stdio.h>
#include <stdlib.h>

void main()
{
  FILE *fp,*fp2;
  float x;
  int i;


  fp = fopen("dumpdata","r");
  fp2 = fopen("../NEW/HI.EXT","w");
  for(i=1;i<=8192;i++){
    fscanf(fp,"%f",&x);
    fprintf(fp2,"%f\n",x);
  }
  fclose(fp);
  fclose(fp2);
}
