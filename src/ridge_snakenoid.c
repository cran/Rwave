#include <stdlib.h>

/***************************************************************
*	$Log: ridge_snakoid.c,v	$			       *
****************************************************************
*							       *
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
****************************************************************/

#include "Swave.h"
#include "denoise.h"




/****************************************************************
*  Function: Ssnakoid:					        *
*  -------------------					        *
*  Ridge extraction from annealing, using the snake method	*
*								*
*    smodulus: modulus smoothed by a window 			*
*    cost:   cost function					*
*    phi: ridge (scale coordinate)				*
*    rho: ridge(position coordinate)				*
*    lambda: coefficient in front of phi'' in the cost function	*
*    mu: coefficient in front of phi' in the cost function	*
*    lambda2: coefficient in front of rho'' in the cost function*
*    mu2: coefficient in front of rho' in the cost function	*
*    c: constant on the temperature schedule			*
*    sigsize: signal size					*
*    snakesize: size of the snake (number of points)		*
*    nscale: total number of scales for CWT			*
*    iteration: maximal number of iterations for the annealing	*
*    stagnant: allowed number of consecutive steps without	*
*              move (stopping criterion)			*
*    seed: seed for random number generator			*
*    count: number of iterations				*
*    sub: subsampling rate for ridge extraction			*
*    blocksize: subsampling of the cost function in cost	*
*               blocksize=1 means that the whole cost function	*
*               is returned					*
*								*
****************************************************************/

void Ssnakenoid_annealing(double *cost, double *smodulus,
  double *phi, double *rho, double *plambda, double *pmu,
  double *plambda2, double *pmu2, double *pc, int *psigsize,
  int *psnakesize, int *pnscale, int *piteration,
  int *pstagnant, int *pseed, int *pcount, int *psub,
  int *pblocksize, int *psmodsize)
{
  int sigsize,snakesize,ncount,iteration;
  int i,k,up,right,pos,num,a,b,count,costcount,sub;
  int count1=0;
  long idum=-9;
  int again, tbox, blocksize,smodsize;
  int nscale, stagnant, recal;
  double c, lambda, mu, lambda2, mu2;
  double *bcost, *phi2;
  double ran, gibbs;
  double cost1;
  double temperature, tmp=0, tmp2;
  FILE *fp;
  int *posmap;
  double tmp_cost;
  double der_plus,der_minus,der_sec,der_sec_plus,der_sec_minus;
  double der_plusB,der_minusB,der_secB,der_sec_plusB,der_sec_minusB;

  /* Generalities; initializations
     -----------------------------*/
  mu = *pmu;
  mu2 = *pmu2;
  lambda = *plambda;
  lambda2 = *plambda2;
  stagnant = *pstagnant;
  nscale = *pnscale;
  iteration = *piteration;
  c = *pc;
  sigsize = *psigsize;
  snakesize = *psnakesize;
  idum = (long)(*pseed);
  sub = *psub;
  blocksize = *pblocksize;
  smodsize = *psmodsize;

  recal = 100000; /* recompute cost function every 'recal' iterations */

  if(!(bcost = (double *)S_alloc(blocksize,sizeof(double))))
    error("Memory allocation failed for bcost at snake_annealing.c \n");

  if(!(phi2 = (double *)S_alloc(sigsize,sizeof(double))))
    error("Memory allocation failed for phi2 at snake_annealing.c \n");

  if(!(posmap = (int *)S_alloc(smodsize * nscale,sizeof(int))))
    error("Memory allocation failed for posmap at snake_annealing.c \n");

  tbox = 0;
  ncount = 0; /* count for cost */
  count = 0; /* total count */
  temperature = c/log(2. + (double)count); /* Initial temperature */
  cost1 = 0;
  tmp_cost = 0;


  /* mark the initial positions of snakes */
  for(i = 0; i < snakesize; i++) {
    k = (int)(rho[i]) + smodsize * (int)(phi[i]);
    posmap[k] = 1;
  }

  /* Smooth and subsample the wavelet transform modulus
     --------------------------------------------------*/
    snakesub(rho,sub,snakesize);

  /* Iterations:
     -----------*/
  while(1) {
    for(costcount = 0; costcount < blocksize; costcount++) {


      /* Initialize the cost function
	 ----------------------------*/
      if(count == 0) {
	for(i = 1; i < snakesize-1; i++) {
	  tmp = (double)((phi[i-1]+ phi[i+1] - 2 * phi[i]));
	  tmp_cost = (double)((lambda * tmp * tmp));

	  tmp = (double)((phi[i] - phi[i+1]));
	  tmp_cost += (double)((mu * tmp * tmp));

	  tmp = (double)((rho[i-1]+ rho[i+1] - 2 * rho[i]));
	  tmp_cost += (double)((lambda2 * tmp * tmp));

	  tmp = (double)((rho[i] - rho[i+1]));
	  tmp_cost += (double)((mu2 * tmp * tmp));

	  a = (int)phi[i]; b= (int)rho[i];
          tmp = smodulus[smodsize * a + b];
	  cost1 -= (tmp * (1.-tmp_cost));
	}

        tmp = (double)((phi[0] - phi[1]));
	tmp_cost = (double) ((mu * tmp * tmp));
        tmp = (double)((rho[0] - rho[1]));
	tmp_cost += (double) ((mu2 * tmp * tmp));

	a = (int)phi[0]; b= (int)rho[0];
        tmp = smodulus[smodsize * a + b];
	cost1 -= (tmp * (1.-tmp_cost));

	a = (int)phi[snakesize-1]; b= (int)rho[snakesize-1];
        tmp = smodulus[smodsize * a + b];
	cost1 -= tmp;

	cost[ncount++] = (double)cost1;
	bcost[0] = (double)cost1;
	count ++;
	costcount = 1;

	printf("Initialisation of cost function done\n");
	if(costcount == blocksize) break;
      }


      /* Generate potential random move
	 ------------------------------*/
      again = YES;
      while(again) {
	randomsnaker(snakesize,&num); /* returns between 0 and 4 * snakesize - 1*/
	pos = num/4;
	up = 0; right=0;
	if(num%4 == 0) up = 1;	  
	if(num%4 == 1) up = -1;
	if(num%4 == 2) right = 1;
	if(num%4 == 3) right = -1;	
	again = NO;
	if((((int)phi[pos] == 0) && (up == -1)) ||
	   (((int)phi[pos] == (nscale-1)) && (up == 1)) ||
	   (((int)rho[pos] == (smodsize-1)) && (right == 1)) ||
	   (((int)rho[pos] == 0) && (right == -1)) ||
	   (posmap[(int)(rho[pos]+right + (phi[pos]+up) * smodsize)] == 1))
	  again = YES; /* boundary effects */

      }

      /* Compute corresponding update of the cost function
	 -------------------------------------------------*/
      if(inrange(2,pos,snakesize-3)) {
	der_plus = (double)(phi[pos +1]-phi[pos]);
	der_minus = (double)(phi[pos] -phi[pos-1]);
	der_sec = der_plus - der_minus;
	der_sec_plus = (double)(phi[pos+2] - 2.*phi[pos+1] + phi[pos]);
	der_sec_minus = (double)(phi[pos-2] - 2.*phi[pos-1] + phi[pos]);

	der_plusB = (double)(rho[pos +1]-rho[pos]);
	der_minusB = (double)(rho[pos] -rho[pos-1]);
	der_secB = der_plusB - der_minusB;
	der_sec_plusB = (double)(rho[pos+2] - 2.*rho[pos+1] + rho[pos]);
	der_sec_minusB = (double)(rho[pos-2] - 2.*rho[pos-1] + rho[pos]);

  	a = (int)phi[pos]; b = (int)rho[pos];
	tmp_cost = (smodulus[smodsize*a+b]);

	a = (int)phi[pos] +  up;  b = (int)rho[pos] + right;
	tmp_cost -= (smodulus[smodsize*a+b]);
  
	tmp = tmp_cost*(-1.+mu*der_plus*der_plus+lambda*der_sec*der_sec
			   +mu2*der_plusB*der_plusB+lambda2*der_secB*der_secB);

	tmp_cost = mu*(1. - 2.*up*der_plus);
	tmp_cost +=  4.*lambda*(1.-up*der_sec);
	tmp_cost += mu2*(1. - 2.*right*der_plusB);
	tmp_cost +=  4.*lambda2*(1.-right*der_secB);
	a = (int)phi[pos] +  up;  b = (int)rho[pos] + right;
	tmp += (tmp_cost*smodulus[smodsize*a+b]);

	tmp_cost = mu*(1. + 2.*up*der_minus);
	tmp_cost += 2.*lambda*(up*der_sec_minus +1.);
	tmp_cost += mu2*(1. + 2.*right*der_minusB);
	tmp_cost += 2.*lambda2*(right*der_sec_minusB +1.);
	a = (int)phi[pos-1];  b = (int)rho[pos-1];
	tmp += (tmp_cost*smodulus[smodsize*a+b]);

	tmp_cost = 2.*lambda*(up*der_sec_plus +1.);
	tmp_cost += 2.*lambda2*(right*der_sec_plusB +1.);
	a = (int)phi[pos+1];  b = (int)rho[pos+1];
	tmp += (tmp_cost*smodulus[smodsize*a+b]);
      }

      if(inrange(2,pos,snakesize-3) == NO) {

	if(pos == 0) {

	  der_plus = (double)(phi[pos +1]-phi[pos]);
	  der_sec_plus = (double)(phi[pos+2] - 2.*phi[pos+1] + phi[pos]);
	  der_plusB = (double)(rho[pos +1]-rho[pos]);
	  der_sec_plusB = (double)(rho[pos+2] - 2.*rho[pos+1] + rho[pos]);

	  a = (int)phi[pos] +  up;  b = (int)rho[pos] + right;
	  tmp_cost = smodulus[smodsize*a+b];

	  a = (int)phi[pos];  b = (int)rho[pos];
	  tmp_cost -= smodulus[smodsize*a+b];

	  tmp = tmp_cost*(-1.+mu*der_plus*der_plus
			  + mu2*der_plusB*der_plusB);

	  a = (int)phi[pos] + up;  b = (int)rho[pos] + right;
	  tmp_cost = mu*(1. - 2.*up*der_plus);
	  tmp_cost += mu2*(1. - 2.*right*der_plusB);
	  tmp += (tmp_cost*smodulus[smodsize*a+b]);

	  a = (int)(phi[pos+1]);  b = (int)(rho[pos+1]);
	  tmp_cost = lambda*(1. + 2.*up*der_sec_plus);
	  tmp_cost += lambda2*(1. + 2.*right*der_sec_plusB);
	  tmp += (tmp_cost*smodulus[smodsize*a+b]);
	}
	else if(pos == (snakesize-1)) {
	  der_minus = (double)(phi[pos] -phi[pos-1]);
	  der_sec_minus = (double)(phi[pos-2] - 2.*phi[pos-1] + phi[pos]);
	  der_minusB = (double)(rho[pos] -rho[pos-1]);
	  der_sec_minusB = (double)(rho[pos-2] - 2.*rho[pos-1] + rho[pos]);

	  a = (int)phi[pos] +  up;  b = (int)rho[pos] + right;
	  tmp_cost = smodulus[smodsize*a+b];

	  a = (int)phi[pos];  b = (int)rho[pos];
	  tmp_cost -= smodulus[smodsize*a+b];

	  tmp = -tmp_cost;

	  a = (int)(phi[pos-1]); b=(int)(rho[pos-1]);
	  tmp_cost = mu*(1. +2.*up*der_minus);
	  tmp_cost += mu2*(1. +2.*right*der_minusB);
	  tmp_cost += lambda*(1.+2*up*der_sec_minus);
	  tmp_cost += lambda2*(1.+2*right*der_sec_minusB);
	  tmp += (tmp_cost*smodulus[smodsize*a+b]);

	}
      }
      

      /* To move or not to move: that's the question
	 -------------------------------------------*/
      if(tmp < (double)0.0) {
	posmap[(int)(rho[pos] + smodsize * phi[pos])] = 0;
	phi[pos] += up; 
	rho[pos] += right; /* good move */
	posmap[(int)(rho[pos] + smodsize * phi[pos])] = 1;
	cost1 += tmp;
	tbox = 0;
      }
      else {
	gibbs = exp(-tmp/temperature);
	ran = ran1(&idum);  

	if(ran < gibbs) {      
	  posmap[(int)(rho[pos] + smodsize * phi[pos])] = 0;
	  phi[pos] += up;
	  rho[pos] += right; /* adverse move */
	  posmap[(int)(rho[pos] + smodsize * phi[pos])] = 1;
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
/*	  if(sub != 1){
	    splsnake(sub, rho-1, phi-1, snakesize, phi2-1);
	    printf("interpolation done\n");
	    for(i=0;i<sigsize;i++) phi[i]=phi2[i];
	  }*/
	  snakexpand(rho, sub, snakesize);
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
	   ---------------------------------*/
/*	if(sub !=1){
	  splsnake(sub, rho-1, phi-1, snakesize, phi2-1);
	  printf("interpolation done\n");
	  for(i=0;i<sigsize;i++) phi[i]=phi2[i];
	}*/
	snakexpand(rho, sub, snakesize);
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
      for(i = 1; i < snakesize-1; i++) {

	der_sec = (double)((phi[i-1]+ phi[i+1]-2 * phi[i]));
	tmp_cost = (double)((lambda * der_sec * der_sec));
	der_plus = (double)((phi[i] - phi[i+1]));
	tmp_cost += (double)((mu * der_plus * der_plus));
	der_secB = (double)((rho[i-1]+ rho[i+1]-2 * rho[i]));
	tmp_cost += (double)((lambda2 * der_secB * der_secB));
	der_plusB = (double)((rho[i] - rho[i+1]));
	tmp_cost += (double)((mu2 * der_plusB * der_plusB));

	a = (int)(phi[i]); b = (int)(rho[i]);
	cost1 -= (smodulus[smodsize*a+b]*(1. - tmp_cost));
      }

      tmp = (double)((phi[0] - phi[1]));
      tmp_cost = (double) ((mu * tmp * tmp));
      tmp = (double)((rho[0] - rho[1]));
      tmp_cost += (double) ((mu2 * tmp * tmp));
      a = (int)(phi[0]); b = (int)(rho[0]);
      cost1 -= (smodulus[smodsize*a+b]*(1.- tmp_cost));

      a = (int)phi[snakesize-1]; b= (int)rho[snakesize-1];
      cost1 -= smodulus[smodsize * a + b];
    }
    cost[ncount++] = (double)cost1;
  }
  /* return; */
}



