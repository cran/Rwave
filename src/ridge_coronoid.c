
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
*  Function: Sridge_coronoid:
*  --------------------------
*  Ridge caracterisation from annealing;
*  Modification of Sridge_annealing, the cost function is
*  replaced with another one, in which the smoothness penalty
*  is proportional to the transform modulus.
*
*    smodulus: smoothed modulus of the wavelet transform
*    cost:   cost function
*    noise: additional potential (coming from learning the noise)
*    phi: ridge
*    lambda: coefficient in front of phi'' in the cost function
*    mu: coefficient in front of phi' in the cost function
*    c: constant on the temperature schedule
*    sigsize: signal size
*    nscale: total number of scales for CWT
*    iteration: maximal number of iterations for the annealing
*    stagnant: allowed number of consecutive steps without
*              move (stopping criterion)
*    seed: seed for random number generator
*    count: number of iterations
*    sub: subsampling rate for ridge extraction
*    costsub: subsampling of the cost function in cost
*               costsub=1 means that the whole cost function
*               is returned
*    smodsize: the subsampling size of signal size
****************************************************************/

void Sridge_coronoid(float *cost, double *smodulus,
  float *phi, float *plambda, float *pmu, float *pc, int *psigsize,
  int *pnscale, int *piteration, int *pstagnant, int *pseed,
  int *pcount, int *psub, int *pblocksize, int *psmodsize)
{
  int i,sigsize,ncount,iteration,up,pos,num,a,count,costcount,sub;
  long idum=-9;
  int again, tbox, blocksize,smodsize;
  int nscale, stagnant, recal;
  double lambda, c, mu;
  float *bcost, *phi2;
  double ran, gibbs;
  double cost1;
  double temperature, tmp=0.0, tmp_cost =0.0, tmp_mod;
  double der_plus,der_minus,der_sec,der_sec_plus,der_sec_minus;
  FILE *fp;


  /* Generalities; initializations
     -----------------------------*/
  mu = (double)(*pmu);
  stagnant = *pstagnant;
  nscale = *pnscale;
  iteration = *piteration;
  lambda = (double)(*plambda);
  c = (double)(*pc);
  sigsize = *psigsize;
  idum = (long) (*pseed);
  sub = *psub;
  blocksize = *pblocksize;
  smodsize = *psmodsize;

  recal = 1000; /* recompute cost function every 'recal' iterations */

  if(!(bcost = (float *)malloc(sizeof(float) * blocksize)))
    error("Memory allocation failed for bcost at ridge_annealing.c \n");

  if(!(phi2 = (float *)calloc((smodsize+1)*sub,sizeof(float))))
    error("Memory allocation failed for phi2 at ridge_annealing.c \n");

/*  if(blocksize != 1) {
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
/*  smoothwt2(modulus,smodulus,sigsize,nscale,sub,&smodsize); */

  /* Subsample the initial guess for the ridge
     -----------------------------------------*/
  for(i=0;i<smodsize;i++){
    phi[i] = phi[sub*i];
  }

  /* Take into account discretization
     --------------------------------*/
  tmp = (double)(sub*sub);
  mu /= tmp;
  lambda /= (tmp*tmp);


  /* Iterations:
     -----------*/
  while(1) {
    for(costcount = 0; costcount < blocksize; costcount++) {


      /* Initialize the cost function
	 ----------------------------*/
      if(count == 0) {
	for(i = 1; i < smodsize-1; i++) {
	  der_sec = (double)(phi[i-1]+ phi[i+1]-2 * phi[i]);
	  tmp_cost = (double)(lambda * der_sec * der_sec);

	  der_plus = (double)(phi[i] - phi[i+1]);
	  tmp_cost += (double)(mu * der_plus * der_plus);

	  a = (int)(phi[i]);
          tmp_mod = smodulus[smodsize * a + i];
/*	  cost1 -= (tmp_mod * tmp_mod - noise[a])*(1. - tmp_cost); */
	  cost1 -= (tmp_mod*(1. - tmp_cost)); 

	}

        tmp = (double)(phi[0] - phi[1]);
	tmp_cost = mu * tmp * tmp;
	a = (int)(phi[0]);
        tmp_mod = smodulus[smodsize * a];
	cost1 -= (tmp_mod*(1.- tmp_cost)); 
/*	cost1 -= (tmp_mod * tmp_mod - noise[a])*(1.- tmp_cost); */

	a = (int)(phi[smodsize-1]);
        tmp_mod = smodulus[smodsize * a + smodsize-1];
/*	cost1 -= (tmp_mod * tmp_mod - noise[a]); */
	cost1 -= tmp_mod;

	cost[ncount++] = (float)cost1;
	bcost[0] = (float)cost1;
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
	if((((int)(phi[pos]) == 0) && (up == -1)) ||
	   (((int)(phi[pos]) == (nscale-1) && (up == 1)))) again = YES; 
	   /* boundary effects */
      }

      /* Compute corresponding update of the cost function
	 -------------------------------------------------*/
      /* tmp: update of cost function */

      if(inrange(2,pos,smodsize-3)) {
	der_plus = (double)(phi[pos +1]-phi[pos]);
	der_minus = (double)(phi[pos] -phi[pos-1]);
	der_sec = der_plus - der_minus;
	der_sec_plus = (double)(phi[pos+2] - 2.*phi[pos+1] + phi[pos]);
	der_sec_minus = (double)(phi[pos-2] - 2.*phi[pos-1] + phi[pos]);

	a = (int)(phi[pos]) +  up;
	tmp_mod = smodulus[smodsize*a+pos];
/*	tmp_cost = (tmp_mod * tmp_mod - noise[a]); */
	tmp_cost = tmp_mod;

	a = (int)(phi[pos]);
	tmp_mod = smodulus[smodsize*a+pos];
/*	tmp_cost -= (tmp_mod * tmp_mod - noise[a]); */
	tmp_cost -= tmp_mod;

	tmp = tmp_cost*(-1.+mu*der_plus*der_plus+lambda*der_sec*der_sec);


	tmp_cost = mu*(1. - 2.*up*der_plus);
	tmp_cost +=  4.*lambda*(1.-up*der_sec);
	a = (int)(phi[pos]) +up;
	tmp_mod = smodulus[smodsize*a+pos];
/*	tmp += tmp_cost*(tmp_mod*tmp_mod - noise[a]); */
	tmp += (tmp_cost*tmp_mod);


	tmp_cost = mu*(1. + 2.*up*der_minus);
	tmp_cost += 2.*lambda*(up*der_sec_minus +1.);
	a = (int)(phi[pos-1]);
	tmp_mod = smodulus[smodsize*a+pos-1];
/*	tmp += tmp_cost*(tmp_mod*tmp_mod - noise[a]); */
	tmp += (tmp_cost*tmp_mod);


	tmp_cost = 2.*lambda*(up*der_sec_plus +1.);
	a = (int)(phi[pos+1]);
	tmp_mod = smodulus[smodsize*a+pos+1];
/*	tmp += tmp_cost*(tmp_mod*tmp_mod - noise[a]); */
	tmp += (tmp_cost * tmp_mod);
      }

      if(inrange(2,pos,smodsize-3) == NO) {

	if(pos == 0) {
	  der_plus = (double)(phi[pos +1]-phi[pos]);
	  der_sec_plus = (double)(phi[pos+2] - 2.*phi[pos+1] + phi[pos]);

	  a = (int)(phi[pos]) +  up;
	  tmp_mod = smodulus[smodsize*a+pos];
/*	  tmp_cost = (tmp_mod * tmp_mod - noise[a]); */
	  tmp_cost = tmp_mod;

	  a = (int)(phi[pos]);
	  tmp_mod = smodulus[smodsize*a+pos];
/*	  tmp_cost -= (tmp_mod * tmp_mod - noise[a]); */
	  tmp_cost -= tmp_mod;

	  tmp = tmp_cost*(-1. + mu*der_plus*der_plus);

	  a = (int)(phi[pos]) + up;
	  tmp_mod = smodulus[smodsize*a+pos];
	  tmp_cost = mu*(1. - 2.*up*der_plus);
/*	  tmp += tmp_cost*(tmp_mod*tmp_mod -noise[a]); */
	  tmp += (tmp_cost*tmp_mod);

	  a = (int)(phi[pos+1]);
	  tmp_mod = smodulus[smodsize*a+pos+1];
	  tmp_cost = lambda*(1. + 2.*up*der_sec_plus);
/*	  tmp += tmp_cost*(tmp_mod*tmp_mod -noise[a]); */
	  tmp += (tmp_cost*tmp_mod);
	}
	else if(pos == (smodsize-1)) {
	  der_minus = (double)(phi[pos] -phi[pos-1]);
	  der_sec_minus = (double)(phi[pos-2] - 2.*phi[pos-1] + phi[pos]);

	  a = (int)(phi[pos]) +  up;
	  tmp_mod = smodulus[smodsize*a+pos];
/*	  tmp_cost = (tmp_mod * tmp_mod - noise[a]); */
	  tmp_cost = tmp_mod;

	  a = (int)(phi[pos]);
	  tmp_mod = smodulus[smodsize*a+pos];
/*	  tmp_cost -= (tmp_mod * tmp_mod - noise[a]); */
	  tmp_cost -= tmp_mod;

	  tmp = -tmp_cost;

	  a = (int)(phi[pos-1]);
	  tmp_mod = smodulus[smodsize*a+pos-1];
	  tmp_cost = mu*(1. +2.*up*der_minus);
	  tmp_cost += lambda*(1.+2*up*der_sec_minus);
/*	  tmp += tmp_cost*(tmp_mod*tmp_mod -noise[a]); */
	  tmp += (tmp_cost*tmp_mod);
	}

      }

      /* To move or not to move: that's the question
	 -------------------------------------------*/
      if(tmp < (double)0.0) {
	phi[pos] = phi[pos] + up; /* good move */
	cost1 += tmp;
	tbox = 0;
      }
      else {
	gibbs = exp(-tmp/temperature);
	ran = ran1(&idum); 
	if(ran < gibbs) {      
	  phi[pos] = phi[pos] + (float)up; /* adverse move */
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
	  } */

	  /* Interpolate from subsampled ridge
	     --------------------------------*/
	  splridge(sub, phi, smodsize, phi2);
	  for(i=0;i<sigsize;i++) phi[i]=phi2[i];
	  /* splridge(1, phi, smodsize, phi2);
	  for(i=0;i<sigsize;i++) phi[i]=phi2[i];*/
	  free((char *)bcost);
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
	   --------------------------------*/
	splridge(sub, phi, smodsize, phi2);
	for(i=0;i<sigsize;i++) phi[i]=phi2[i];
	/* splridge(1, phi, smodsize, phi2);
	for(i=0;i<sigsize;i++) phi[i]=phi2[i];*/
	free((char *)bcost);
	printf("Done !\n");
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
	for(i = 1; i < smodsize-1; i++) {
	  der_sec = (double)((phi[i-1]+ phi[i+1]-2 * phi[i]));
	  tmp_cost = (double)((lambda * der_sec * der_sec));

	  der_plus = (double)((phi[i] - phi[i+1]));
	  tmp_cost += (double)((mu * der_plus * der_plus));

	  a = (int)(phi[i]);
          tmp_mod = smodulus[smodsize * a + i];
	  cost1 -= (tmp_mod*(1. - tmp_cost));
/*	  cost1 -= (tmp_mod * tmp_mod - noise[a])*(1. - tmp_cost); */
	}
      
      tmp = (double)((phi[0] - phi[1]));
      tmp_cost = (double) ((mu * tmp * tmp));
      a = (int)(phi[0]);
      tmp_mod = smodulus[smodsize * a];
/*      cost1 -= (tmp_mod * tmp_mod - noise[a])*(1.- tmp_cost); */
      cost1 -= (tmp_mod*(1.- tmp_cost));

      a = (int)(phi[smodsize-1]);
      tmp_mod = smodulus[smodsize * a + smodsize-1];
/*      cost1 -= (tmp_mod * tmp_mod - noise[a]); */
      cost1 -= tmp_mod;
      
    }
    cost[ncount++] = (float)cost1;
  }
  /* return; */
}

       




    
    
  
  

  

