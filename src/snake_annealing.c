/***************************************************************
*	$Log: s_annealing.c,v	$			       *
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
*  Function: Ssnake_annealing:					*
*  --------------------------					*
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

void Ssnake_annealing(float *cost, double *smodulus,
  float *phi, float *rho, float *plambda, float *pmu,
  float *plambda2, float *pmu2, float *pc, int *psigsize,
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
  float c, lambda, mu, lambda2, mu2;
  float *bcost, *phi2;
  double ran, gibbs;
  double cost1;
  double temperature, tmp=0, tmp2;
  FILE *fp;
  int *posmap;

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

  if(!(bcost = (float *)calloc(blocksize,sizeof(float))))
    error("Memory allocation failed for bcost at snake_annealing.c \n");

  if(!(phi2 = (float *)calloc(sigsize,sizeof(float))))
    error("Memory allocation failed for phi2 at snake_annealing.c \n");

  if(!(posmap = (int *)calloc(smodsize * nscale,sizeof(int))))
    error("Memory allocation failed for posmap at snake_annealing.c \n");

/*  if(blocksize != 1) {
    if((fp = fopen("annealing.cost","w")) == NULL)
      error("can't open file at snake_annealing.c \n");
  }
*/  
  tbox = 0;
  ncount = 0; /* count for cost */
  count = 0; /* total count */
     temperature = c/log(2. + (double)count); /* Initial temperature */
  cost1 = 0;


  /* mark the initial positions of snakes */
  for(i = 0; i < snakesize; i++) {
    /* k = (int)(rho[i] + smodsize * phi[i]); this is probably a bug */
    k = (int)(rho[i]) + smodsize * (int)(phi[i]);
    posmap[k] = 1;
  }

  /* Smooth and subsample the wavelet transform modulus
     --------------------------------------------------*/
/*  if(sub !=1){*/
/*    smoothwt(modulus,smodulus,sigsize,nscale,sub); */
    snakesub(rho,sub,snakesize);
/*  }
  if(sub == 1)
    for(i=0; i<sigsize*nscale;i++)
      smodulus[i]=modulus[i];
*/
  /* Iterations:
     -----------*/
  while(1) {
    for(costcount = 0; costcount < blocksize; costcount++) {


      /* Initialize the cost function
	 ----------------------------*/
      if(count == 0) {
	for(i = 1; i < snakesize-1; i++) {
	  tmp = (double)((phi[i-1]+ phi[i+1] - 2 * phi[i]));
	  cost1 += (double)((lambda * tmp * tmp));

	  tmp = (double)((phi[i] - phi[i+1]));
	  cost1 += (double)((mu * tmp * tmp));

	  tmp = (double)((rho[i-1]+ rho[i+1] - 2 * rho[i]));
	  cost1 += (double)((lambda2 * tmp * tmp));

	  tmp = (double)((rho[i] - rho[i+1]));
	  cost1 += (double)((mu2 * tmp * tmp));

	  a = (int)phi[i]; b= (int)rho[i];
          tmp = smodulus[smodsize * a + b];
/*	  cost1 -= (tmp * tmp - noise[a]); */
	  cost1 -= tmp;
	}

        tmp = (double)((phi[0] - phi[1]));
	cost1 += (double) ((mu * tmp * tmp));
        tmp = (double)((rho[0] - rho[1]));
	cost1 += (double) ((mu2 * tmp * tmp));

	a = (int)phi[0]; b= (int)rho[0];
        tmp = smodulus[smodsize * a + b];
/*	cost1 -= (tmp * tmp - noise[a]); */
	cost1 -= tmp;

	a = (int)phi[snakesize-1]; b= (int)rho[snakesize-1];
        tmp = smodulus[smodsize * a + b];
/*	cost1 -= (tmp * tmp - noise[a]); */
	cost1 -= tmp;

	cost[ncount++] = (float)cost1;
	bcost[0] = (float)cost1;
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
	tmp = (double)(lambda*up);
	tmp *=(double)((6*up+(12*phi[pos]-8*(phi[pos-1]+phi[pos+1])+2*(phi[pos-2]+phi[pos+2]))));

	tmp += (double)(mu*up*(4.0*phi[pos]-2.0*(phi[pos-1]+phi[pos+1])+2.0*up));

	tmp2 = (double)(lambda2*right);
	tmp2 *=(double)((6*right+(12*rho[pos]-8*(rho[pos-1]+rho[pos+1])+2*(rho[pos-2]+rho[pos+2]))));

	tmp2 += (double)(mu2*right*(4.0*rho[pos]-2.0*(rho[pos-1]+rho[pos+1])+2.0*right));
	tmp += tmp2;

	a = (int)phi[pos]; b = (int)rho[pos];
	tmp += (smodulus[smodsize*a+b]);
/*	tmp += ((smodulus[smodsize*a+b]*smodulus[smodsize*a+b])-noise[a]); */
	a = (int)phi[pos] +  up;  b = (int)rho[pos] + right;
	tmp -= (smodulus[smodsize*a+b]);
/*	tmp -= ((smodulus[smodsize*a+b]*smodulus[smodsize * a + b])-noise[a]); */
      }

      if(inrange(2,pos,snakesize-3) == NO) {
	tmp = (double)(lambda*up);
	tmp2 = (double)(lambda2*right);
	if(pos == 0) {
	  tmp *= (double)((up+2.0*(phi[0]-2*phi[1]+phi[2])));
	  tmp += (double)(mu*up*((2.0*phi[pos]-2.0*phi[pos+1]) + up));
	  tmp2 *= (double)((right+2.0*(rho[0]-2*rho[1]+rho[2])));
	  tmp2 += (double)(mu2*right*((2.0*rho[pos]-2.0*rho[pos+1]) + right));
	}
	else if(pos == 1) {
	  tmp *= (double)((5*up+2.0*(-2*phi[0]+5*phi[1]-4*phi[2]+phi[3])));
	  tmp += (double)(mu*up*(4.0*phi[pos]-2.0*(phi[pos-1]+phi[pos+1])+2.0*up));
	  tmp2 *= (double)((5*right+2.0*(-2*rho[0]+5*rho[1]-4*rho[2]+rho[3])));
	  tmp2 += (double)(mu2*right*(4.0*rho[pos]-2.0*(rho[pos-1]+rho[pos+1])+2.0*right));
	}
	else if(pos == (snakesize-2)) {
	  tmp *= (double)((5*up+2.0*(phi[pos-2]-4*phi[pos-1]+5*phi[pos]-2*phi[pos+1])));
	  tmp += (double)(mu*up*(4.0*phi[pos]-2.0*(phi[pos-1]+phi[pos+1])+2.0*up));
	  tmp2 *= (double)((5*right+2.0*(rho[pos-2]-4*rho[pos-1]+5*rho[pos]-2*rho[pos+1])));
	  tmp2 += (double)(mu2*right*(4.0*rho[pos]-2.0*(rho[pos-1]+rho[pos+1])+2.0*right));
	}
	else if(pos == (snakesize-1)) {
	  tmp *= (double)((up+2.0*(phi[pos-2]-2*phi[pos-1]+phi[pos])));
	  tmp += (double)(mu*up*((2.0*phi[pos]-2.0*phi[pos-1]) + up));
	  tmp2 *= (double)((right+2.0*(rho[pos-2]-2*rho[pos-1]+rho[pos])));
	  tmp2 += (double)(mu2*right*((2.0*rho[pos]-2.0*rho[pos-1]) + right));
	}
	tmp += tmp2;

	a = (int)phi[pos]; b = (int)rho[pos];
	tmp +=(smodulus[smodsize*a + b]);
/*	tmp +=((smodulus[smodsize*a + b]*smodulus[smodsize*a+b])-noise[a]); */
	a = (int)phi[pos] +  up; b = (int)rho[pos] + right;
	tmp -=(smodulus[smodsize*a + b]);
/*	tmp -=((smodulus[smodsize*a + b]*smodulus[smodsize*a+b])-noise[a]); */
      
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
	  cost[ncount++] = (float)cost1;
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
	  free((char *)bcost);
	  free((char *)posmap);
	  return;
	}
      }
      bcost[costcount] = (float)cost1;

      count ++;
      if(count >=  iteration) 	{
	cost[ncount++] = (float)cost1;
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
	free((char *)bcost);
	free((char *)posmap);
	return;
      }
      temperature = c/log(1. + (double)count);
    }

    bcost[blocksize-1] = (float)cost1;
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

	tmp = (double)((phi[i-1]+ phi[i+1] - 2 * phi[i]));
	cost1 += (double)((lambda * tmp * tmp));
	
	tmp = (double)((phi[i] - phi[i+1]));
	cost1 += (double)((mu * tmp * tmp));
	
	tmp = (double)((rho[i-1]+ rho[i+1] - 2 * rho[i]));
	cost1 += (double)((lambda2 * tmp * tmp));
	
	tmp = (double)((rho[i] - rho[i+1]));
	cost1 += (double)((mu2 * tmp * tmp));

	a = (int)phi[i]; b= (int)rho[i];
	tmp = smodulus[smodsize * a + b];
	cost1 -= tmp;
/*	cost1 -= (tmp * tmp - noise[a]); */


      }

      tmp = (double)((phi[0] - phi[1]));
      cost1 += (double) ((mu * tmp * tmp));
      tmp = (double)((rho[0] - rho[1]));
      cost1 += (double) ((mu2 * tmp * tmp));
      
      a = (int)phi[0]; b= (int)rho[0];
      tmp = smodulus[smodsize * a + b];
      cost1 -= tmp;
/*      cost1 -= (tmp * tmp - noise[a]); */
      
      a = (int)phi[snakesize-1]; b= (int)rho[snakesize-1];
      tmp = smodulus[smodsize * a + b];
      cost1 -= tmp;
/*      cost1 -= (tmp * tmp - noise[a]); */
    }
    cost[ncount++] = (float)cost1;
  }
  /* return; */
}



