
/***************************************************************
*    $Log: bee_annealing.c,v $                                 *
****************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
****************************************************************/



/****************************************************************
*  Function: Sbee_annealing:
*  --------------------------
*  Ridges caracterisation with crazy_climber algorithm
*
*    modulus: modulus of the wavelet transform
*     (coming from learning the noise)
*    beemap: ridges
*    c: constant on the temperature schedule
*    sigsize: signal size
*    nscale: total number of scales for CWT
*     (or frequencies for Gabor)
*    nbbee: number of crazy climbers
*    iteration: number of moves for each climber
*    seed: seed for random number generator
*    bstep: stepsize for moves in b direction
*    integral: flag; if set to FALSE, the function computes
*      the "characteristic function" of the ridge; if set to
*      TRUE, computes the characteristic function, weighted by
*      the value of the modulus.
*    chain: flag; if set to TRUE, returns "chained ridges"
*
****************************************************************/


#include "Swave.h"


void Sbee_annealing(double *smodulus, double *beemap,
	       float *pc,
	       int *psigsize, int *pnscale, int *piteration,
               int *pseed, int *pbstep, int *pnbbee,
               int *pintegral, int *pchain, int *flag)
{
  double r, dd, ee;
  int i, bstep, k, k1, k2, bee;
  int *a, *b, integral, chain, s, count;
  int seed, nscale, iteration, sigsize, nbbee, tstep;
  long idum;
  float c;
  double ran1();

  chain = *pchain;
  integral = *pintegral;
  nbbee = *pnbbee;
  bstep = *pbstep;
  nscale = *pnscale;
  iteration = *piteration;
  c = *pc;
  sigsize = *psigsize;
  seed = *pseed;
  idum = (long)seed;

  if(!(a = (int *)malloc(sizeof(int) * iteration)))
    error("Memory allocation failed for a in bee_annealing.c \n");
  if(!(b = (int *)malloc(sizeof(int) * iteration)))
    error("Memory allocation failed for b in bee_annealing.c \n");
  

  /* generation of the random moves of the climbers
     ----------------------------------------------*/
  for(bee = 0; bee < nbbee; bee++) {

    /* Initialisation */
    a[0] = (int)((nscale-1) * ran1(&idum));
    b[0] = (int)((sigsize-1) * ran1(&idum));
    a[0] = min(nscale-1, a[0]);
    b[0] = min(sigsize-1, b[0]);
    a[0] = max(0, a[0]);
    b[0] = max(0, b[0]);
    k = b[0] + a[0] * sigsize;
    /* By Wen 6/5/97 */
    if(integral) 
      beemap[k] = beemap[k] + smodulus[k];
    else 
      beemap[k] = beemap[k] +1;
    /*    beemap[k] = beemap[k] +1;*/

    /* Iterations */
    for(i =1; i < iteration; i++) {
      r = (double)ran1(&idum);
      if (r >= 0.5) 
	b[i] = min(sigsize-1,b[i-1]+bstep);
      else 
	b[i] = max(0,b[i-1]-bstep);

      r = (double)ran1(&idum);

      if (r >= 0.5) 
	a[i] = min(nscale-1,a[i-1]+1);
      else 
	a[i] = max(0,a[i-1]-1);

      k1 = b[i] + sigsize * a[i];
      k2 = b[i] + sigsize * a[i-1];
      dd = smodulus[k1] - smodulus[k2];
/*      dd = smodulus[k1] - smodulus[k2] - noise[a[i]] + noise[a[i-1]]; */
      if (dd<0) {
	r = (double)ran1(&idum);
	ee = exp(dd * log((double)(3.0 + i))/c);
	if((*flag)==1) ee=exp(dd * log((double)(3.0))/c); /* constant temperature */
	if (r>ee) a[i] = a[i-1];
      }

      /* Chaining of the ridges */
      if(chain) {
	count = 1;

        /* The real move By Wen 6/7/97 */
        tstep = b[i] - b[i-1];
        if(tstep  < 0) tstep = -tstep;
	while(count < tstep) {
	  /*	while(count < bstep) { By Wen 6/7/97 */
	  if(b[i] - b[i-1] > 0) { 
	    k1 = b[i-1] + count + sigsize * a[i];
	    k2 = b[i-1] + count + sigsize * a[i-1];
/*	    if(smodulus[k1] - noise[a[i]] > smodulus[k2] - noise[a[i-1]]){ */
	    if(smodulus[k1] > smodulus[k2]){
	      k = k1;
	      s = a[i];
	    }
	    else {
	      k = k2;
	      s = a[i-1];
	    }
	  }
	  if(b[i] - b[i-1] < 0) {
	    k1 = b[i-1] - count + sigsize * a[i];
	    k2 = b[i-1] - count + sigsize * a[i-1];
/*	    if(smodulus[k1] - noise[a[i]] > smodulus[k2] - noise[a[i-1]])  { */
	    if(smodulus[k1] > smodulus[k2])  {
	      k = k1;
	      s = a[i];
	    }
	    else {
	      k = k2;
	      s = a[i-1];
	    }
	  }
	  /* update of the ridges */
	  if(integral) 
	    beemap[k] = beemap[k] + smodulus[k];
	  else 
	    beemap[k] = beemap[k] +1;
	  count ++;
	}
      }

      k = b[i] + a[i] * sigsize;
      if(integral) 
	beemap[k] = beemap[k] + smodulus[k];
      else 
	beemap[k] = beemap[k] +1;
    }
  }
}


       



    

  


    
    
  
  

  

