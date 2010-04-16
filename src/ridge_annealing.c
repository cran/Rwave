#include <stdlib.h>


/***************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
****************************************************************/

#include "Swave.h"
#include "denoise.h"






/****************************************************************
*  Function: Sridge_annealing:
*  --------------------------
*  Ridge caracterisation from annealing
*
*    smodulus: smoothed modulus of the wavelet transform
*    cost:   cost function
*    phi: ridge
*    lambda: coefficient in front of phi' in the cost function
*    mu: coefficient in front of phi'' in the cost function
*    c: constant on the temperature schedule
*    sigsize: signal size
*    nscale: total number of scales for CWT
*    iteration: maximal number of iterations for the annealing
*    stagnant: allowed number of consecutive steps without
*              move (stopping criterion)
*    seed: seed for random number generator
*    count: number of iterations
*    sub: subsampling rate for ridge extraction
*    blocksize: subsampling of the cost function in cost
*               blocksize=1 means that the whole cost function
*               is returned
*    smodsize: the size of subsampled signal
****************************************************************/

void Sridge_annealing(double *cost, double *smodulus,
  double *phi, double *plambda, double *pmu, double *pc, int *psigsize,
  int *pnscale, int *piteration, int *pstagnant, int *pseed,
  int *pcount, int *psub, int *pblocksize, int *psmodsize)
{
  int i,sigsize,ncount,iteration,up,pos,num,a,count, costcount,sub;
  long idum=-9;
  int again, tbox, blocksize,smodsize;
  int nscale, stagnant, recal;
  double lambda, c, mu;
  double *bcost, *phi2;
  double ran, gibbs;
  double cost1;
  double temperature, tmp=0.0;
  FILE *fp;


  /* Generalities; initializations
     -----------------------------*/
  mu = *pmu;
  stagnant = *pstagnant;
  nscale = *pnscale;
  iteration = *piteration;
  lambda = *plambda;
  c = *pc;
  sigsize = *psigsize;
  idum = *pseed;
  sub = *psub;
  blocksize = *pblocksize;
  smodsize = *psmodsize;

  recal = 1000000; /* recompute cost function every 'recal' iterations */

  if(!(bcost = (double *) R_alloc(blocksize, sizeof(double) )))
    error("Memory allocation failed for bcost at ridge_annealing.c \n");

  if(!(phi2 = (double *)S_alloc((smodsize+1)*sub,sizeof(double))))
    error("Memory allocation failed for phi2 at ridge_annealing.c \n");

/*
  if(blocksize != 1) {
    if((fp = fopen("annealing.cost","w")) == NULL)
      error("can't open file at ridge_annealing.c \n");
  }
*/
  tbox = 0;
  ncount = 0; /* count for cost */
  count = 0; /* total count */
  temperature = c/log(2. + (double)count); /* Initial temperature */
  cost1 = 0;


  /* Smooth and subsample the wavelet transform modulus
     --------------------------------------------------*/
/*  smoothwt(modulus,smodulus,sigsize,nscale,sub); */
/*  smoothwt2(modulus,smodulus,sigsize,nscale,sub, &smodsize); */
/*   printf("smodsize=%d\n",smodsize); */
/*  for(i=0;i<smodsize;i++){
    phi[i] = phi[sub*i];
  } */

  for(i=0;i<smodsize;i++){
    phi[i] = phi[(int)((sigsize-1)/(smodsize-1)*i)];
  }


  /* Iterations:
     -----------*/

  while(1) {
    for(costcount = 0; costcount < blocksize; costcount++) {


      /* Initialize the cost function
	 ----------------------------*/
      if(count == 0) {
	for(i = 1; i < smodsize-1; i++) {
	  tmp = (double)((phi[i-1]+ phi[i+1]-2 * phi[i]));
	  cost1 += (double)((lambda * tmp * tmp));

	  tmp = (double)((phi[i] - phi[i+1]));
	  cost1 += (double)((mu * tmp * tmp));

	  a = (int)phi[i];
          tmp = smodulus[smodsize * a + i];
/*	  cost1 -= (tmp * tmp - noise[a]); */
	  cost1 -= tmp;
	}

        tmp = (double)((phi[0] - phi[1]));
	cost1 += (double) ((mu * tmp * tmp));
	a = (int)phi[0];
        tmp = smodulus[smodsize * a];
/*	cost1 -= (tmp * tmp - noise[a]); */
	cost1 -= tmp;

	a = (int)phi[smodsize-1];
        tmp = smodulus[smodsize * a + smodsize-1];
/*	cost1 -= (tmp * tmp - noise[a]); */
	cost1 -= tmp;

	cost[ncount++] = (double)cost1;
	bcost[0] = (double)cost1;
	count ++;
	costcount = 1;
	if(costcount == blocksize) break;
      }

      /* Generate potential random move
	 ------------------------------*/
      again = YES;
      while(again) {
	randomwalker2(smodsize,&num,&idum); 
	/* returns between 0 and 2 * smodsize - 1*/
	pos = num/2;
	up = -1;
	if(num%2 == 0) up = 1;
	again = NO;
	if((((int)phi[pos] == 0) && (up == -1)) ||
	   (((int)phi[pos] == (nscale-1) && (up == 1)))) again = YES; 
	   /* boundary effects */
      }

      /* Compute corresponding update of the cost function
	 -------------------------------------------------*/
      if(inrange(2,pos,smodsize-3)) {
	tmp = (double)(lambda*up);
	tmp *=(double)((6*up+(12*phi[pos]-8*(phi[pos-1]+phi[pos+1])
		+2*(phi[pos-2]+phi[pos+2]))));

	tmp += (double)(mu*up*(4.0*phi[pos]
		-2.0*(phi[pos-1]+phi[pos+1])+2.0*up));

	a = (int)phi[pos];
/*	tmp += ((smodulus[smodsize*a+pos]*smodulus[smodsize*a+pos])
		-noise[a]);  */
	tmp += smodulus[smodsize*a+pos];
	a = (int)phi[pos] +  up;
/*	tmp -= ((smodulus[smodsize*a+pos]*smodulus[smodsize * a + pos])
		-noise[a]); */
	tmp -= smodulus[smodsize*a+pos];
      }

      if(inrange(2,pos,smodsize-3) == NO) {
	tmp = (double)(lambda*up);
	if(pos == 0) {
	  tmp *= (double)((up+2.0*(phi[0]-2*phi[1]+phi[2])));
	  tmp += (double)(mu*up*((2.0*phi[pos]-2.0*phi[pos+1]) + up));
	}
	else if(pos == 1) {
	  tmp *= (double)((5*up+2.0*(-2*phi[0]+5*phi[1]-4*phi[2]+phi[3])));
	  tmp += (double)(mu*up*(4.0*phi[pos]-2.0*(phi[pos-1]+phi[pos+1])
				 +2.0*up));
	}
	else if(pos == (smodsize-2)) {
	  tmp *= (double)((5*up+2.0*(phi[pos-2]-4*phi[pos-1]+5*phi[pos]
				     -2*phi[pos+1])));
	  tmp += (double)(mu*up*(4.0*phi[pos]-2.0*(phi[pos-1]+phi[pos+1])
				 +2.0*up));
	}
	else if(pos == (smodsize-1)) {
	  tmp *= (double)((up+2.0*(phi[pos-2]-2*phi[pos-1]+phi[pos])));
	  tmp += (double)(mu*up*((2.0*phi[pos]-2.0*phi[pos-1]) + up));
	}
	a = (int)phi[pos];
/*	tmp +=((smodulus[smodsize*a+pos]*smodulus[smodsize*a+pos])-noise[a]); */
	tmp += smodulus[smodsize*a+pos];
	a = (int)phi[pos] +  up;
/*	tmp -=((smodulus[smodsize*a+pos]*smodulus[smodsize*a+pos])-noise[a]); */
	tmp -= smodulus[smodsize*a+pos];
      }

      /* To move or not to move: that's the question
	 -------------------------------------------*/
      if(tmp < (double)0.0) {
	phi[pos] = phi[pos] + up; /* good move */
	if(phi[pos] < 0) printf("Error \n");
	cost1 += tmp;
	tbox = 0;
      }
      else {
	gibbs = exp(-tmp/temperature);
	ran = ran1(&idum); 
	if(ran < gibbs) {      
	  phi[pos] = phi[pos] + up; /* adverse move */
	  cost1 += tmp;
	  tbox = 0;
	}
	tbox ++;
	if(tbox >= stagnant)  {
	  cost[ncount++] = (double)cost1;
	  *pcount = ncount;
/*	  if((blocksize != 1)){
	    for(i = 0; i < costcount+1; i++)
	      fprintf(fp, "%f ", bcost[i]);
	    fclose(fp);
	  }
*/
	  /* Interpolate from subsampled ridge
	     --------------------------------*/
	  splridge(sub, phi, smodsize, phi2);
	  for(i=0;i<sigsize;i++) phi[i]=phi2[i];
	  /* splridge(1, phi, smodsize, phi2);
	  for(i=0;i<sigsize;i++) phi[i]=phi2[i];*/
	  return;
	}
      }
      bcost[costcount] = (double)cost1;

      count ++;
      if(count >=  iteration) 	{
	cost[ncount++] = (double)cost1;
	*pcount = ncount;

	/* Write cost function to a file
	   -----------------------------*/
/*	if((blocksize != 1)){
	  for(i = 0; i < costcount+1; i++)
	    fprintf(fp, "%f ", bcost[i]);
	  fclose(fp);
	}
*/
	/* Interpolate from subsampled ridge
	   --------------------------------*/
	splridge(sub, phi, smodsize, phi2);
	for(i=0;i<sigsize;i++) phi[i]=phi2[i];
	/* splridge(1, phi, smodsize, phi2);
	for(i=0;i<sigsize;i++) phi[i]=phi2[i];*/
	//printf("Done !\n");
	return;
      }
      temperature = c/log(1. + (double)count);
    }

    bcost[blocksize-1] = (double)cost1;
    if((blocksize != 1)){
/*      for(i = 0; i < blocksize; i++)
	fprintf(fp, "%f ", bcost[i]); */
      for(i = 0; i < blocksize; i++)
	bcost[i] = 0.0;
    }

    /* recalculate cost to prevent error propagation 
       ---------------------------------------------*/
    if(count % recal == 0) {
      cost1 = 0.0;
      for(i = 1; i < smodsize-1; i++) {
	tmp = (double)((phi[i-1]+ phi[i+1]-2 * phi[i]));
	cost1 += (double)((lambda * tmp * tmp));
	tmp = (double)((phi[i]-phi[i+1]));
	cost1 += (double)((mu * tmp * tmp));

	a = (int)phi[i];
	cost1 -= smodulus[smodsize * a + i];
/*	cost1 -= (smodulus[smodsize * a + i] * smodulus[smodsize * a + i]
		  -noise[a]); */
      }
      a = (int)phi[0];
      tmp = (double)((phi[0]-phi[1]));
      cost1 += (double)((mu * tmp * tmp));
      cost1 -= smodulus[smodsize * a];
/*      cost1 -= (smodulus[smodsize * a] * smodulus[smodsize * a]
		-noise[a]); */
      a = (int)phi[smodsize-1];
/*      cost1 -= (smodulus[smodsize * a + smodsize-1] *
		smodulus[smodsize * a + smodsize-1] -noise[a]); */
      cost1 -= smodulus[smodsize * a + smodsize-1];
    }
    cost[ncount++] = (double)cost1;
  }
  /* return; */
}

       




    
    
  
  

  

