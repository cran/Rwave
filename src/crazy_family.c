#include <stdlib.h>


/***************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
****************************************************************/


/***************************************************************
*                                                              *
*                                                              *
*         Bug fix by Tom Price 2007, t0mpr1c3@gmail.com        *
*                                                              *
*                                                              *
****************************************************************/

#include "Swave.h"
#include "denoise.h"



/****************************************************************
*   Function: Scrazy_family
*   ---------
*      Organize the output of the crazy climber algorithm into
*      a series of connected ridges.
*
*   ridgemap: output of the crazy climber algorithm
*   sigsize: input signal size
*   nscale: number of scales (or frequencies) of the transform
*   nbchain: number of ridges
*   chain: array containing the ridges; the first element
*          contains the starting point of the ridge, the
*          second contains the ridge length, and the other
*          are the values of the ridge
*   id: pointer containing the indes of the considered ridge
*   threshold: minimal value of transform modulus considered
*              on a ridge
*****************************************************************/

void Scrazy_family(double *ridgemap,double *orderedmap,int *chain,
		   int *pnbchain,int *psigsize,int *pnscale,
		   int *pbstep,double *pthreshold)
{
  int bstep, nscale, sigsize, nbchain;
  int i, j, k, id, count, a, b, found, k1;
  double *mridge;
  double threshold;

  threshold = *pthreshold;
  bstep = *pbstep;
  nscale = *pnscale;
  sigsize = *psigsize;
  nbchain = *pnbchain;

  if(!(mridge = (double *)S_alloc(sigsize * nscale, sizeof(double))))
    error("Memory allocation failed for mridge in crazy_family.c \n");


  /* Compute local maxima of modulus (for fixed b)
     -------------------------------------------- */
  Scwt_mridge(ridgemap,mridge,psigsize,pnscale);


  id = 0;

  /* Start looking for ridges as connected curves
     -------------------------------------------- */
  for(i = 0; i < sigsize; i+= bstep) {
    for(j = 0; j < nscale; j++) {
      b = i;
      a = j;
      k = b + sigsize * a;
      if((mridge[k] > 0.000001) && (orderedmap[k] == 0.0)) {
	found = YES;

	/* backwarding: looking for previous points of the chain
	   ----------------------------------------------------- */
	while(found && (b > 0)) {
	  found = NO;
	  b = b-1;
	  k1 = b + sigsize * max(a-1, 0);
	  if((mridge[k1] > 0.000001) && (orderedmap[k1] == 0.0)) {
	    found = YES;
	    a = max(a-1,0);
	  }
	  else {
	    k1 = b + sigsize * max(a, 0);
	    if((mridge[k1] > 0.000001) && (orderedmap[k1] == 0.0)) {
	      found = YES;
	    }
	    else {
	      k1 = b + sigsize * min(a+1, nscale-1);
	      if((mridge[k1] > 0.000001) && (orderedmap[k1] == 0.0)) {
		found = YES;
		a = min(a+1,nscale-1);
	      }
	    }
	  }
	}

	/* forwarding
	   ---------- */
	id ++;
	if(id > nbchain) {
	  printf("Nb of chains > reserved number. Increase the nbchain. \n");
	  return;
	}
	b = b+1;
	k = b + sigsize * a;
	chain[(id-1)] = b;
	chain[(id-1) + nbchain]= a;
	count = 2;
	found = YES;

	while(found) {
	  orderedmap[k] =(double)(id);
	  found = NO;

	  b = min(b+1, sigsize-1);
	  k1 = b + sigsize * max(a-1, 0);
	  if((mridge[k1] > 0.000001) && (orderedmap[k1] == 0.0)) {
	    found = YES;
	    a = max(a-1,0);
	  }
	  else {
	    k1 = b + sigsize * max(a, 0);
	    if((mridge[k1] > 0.000001) && (orderedmap[k1] == 0.0)) {
	      found = YES;
	    }
	    else {
	      k1 = b + sigsize * min(a+1, nscale-1);
	      if((mridge[k1] > 0.000001) && (orderedmap[k1] == 0.0)) {
		found = YES;
		a = min(a+1,nscale-1);
	      }
	    }
	  }
	  if(found) {
	    k = b + sigsize * a;
	    chain[(id -1) + nbchain * count] = a;
	    count ++;
	  }
	}
	/* Threshold and chain the ridges
	   ------------------------------ */
	chain_thresholded(mridge,sigsize,chain,&id,nbchain,threshold,bstep);
      }
    }
  }

  /* Generate the image of the chains
     -------------------------------- */
  orderedmap_thresholded(orderedmap,sigsize,nscale,chain,nbchain);


  /* Order the output according to the chosen data structure
     ------------------------------------------------------- */
  reordering(chain, sigsize, nbchain);


  printf("There are %d chains. \n", id);
  *pnbchain = id;

  return;
}


/****************************************************************
*   Function: orderedmap_thresholded
*   ---------
*      Organizes the thresholded chains into an image whose
*      size is the same as that of the time-frequency
*      transform
*
*   orderedmap: image containing the ridges
*   sigsize: input signal size
*   nscale: number of scales (or frequencies) of the transform
*   nbchain: number of ridges
*   chain: array containing the ridges
*****************************************************************/

void orderedmap_thresholded(double *orderedmap,int sigsize,int nscale,
			    int *chain,int nbchain)
{
  int id, i, j, k;
  int a, b;

  for(i = 0; i < sigsize; i++)
    for(j = 0; j < nscale; j++) {
      k = i + j * sigsize;
      orderedmap[k] = 0.0;
    }

  for(id = 0; id < nbchain; id++) {
    k = id;
    b = chain[k];
    k += nbchain;
    a = chain[k];
    while((a != -1)) {
      orderedmap[b + a * sigsize] = (double)(id + 1);
      b++;
      k += nbchain;
      a = chain[k];
    }
  }
  return;
}


/****************************************************************
*   Function: chain_thresholded
*   ---------
*      Check the length of chained ridges, and threshold them
*
*   mridge: fixed time local maxima of transform (output of
*           Scwt_mridge
*   sigsize: input signal size
*   nscale: number of scales (or frequencies) of the transform
*   nbchain: number of ridges
*   chain: array containing the ridges; the first element
*          contains the starting point of the ridge, the
*          second contains the ridge length, and the other
*          are the values of the ridge
*   id: pointer containing the indes of the considered ridge
*   threshold: minimal value of transform modulus considered
*              on a ridge
*****************************************************************/

void chain_thresholded(double *mridge,int sigsize,int *chain,int *id,
		       int nbchain,double threshold, int bstep)
{
  int i, j, k, k1, a, b;
  int count, found, lng;
  int bstart, astart, bend, aend, oldbstart;
  int chnlng;

  /* look for actual beginning of the chain */

  k = (*id)-1;
  b = chain[k];
  k +=  nbchain;
  a = chain[k];
  k1 = b + sigsize * a;
  while((chain[k] != -1) && (mridge[k1] < (double)threshold)) {
    k += nbchain;
    a = chain[k];
    b ++;
    k1 = b + sigsize * a;
  }
  /* Invalid chain (not enough energy) */
  if(chain[k] == -1) {
    chnlng = sigsize + 2;
    k = (*id)-1;
    for(i = 0; i < chnlng; i++) {
      chain[k] = -1;
      k += nbchain;
    }
    (*id)--;
    return;
  }

  astart = a;
  bstart = b;

  /* move along the non-thresholded chain */
  while((b<sigsize)&&(chain[k] != -1)){
    k+= nbchain;
    b++;
  }
  /* This code removed for bug fix */
  /*
  b--;
  k-= nbchain;
  */
  /* end of code removed for bug fix */
  /* This code inserted for bug fix */
  if ( b > bstart ) {
    b--;
    k-= nbchain;
  }
  /* end of code inserted for bug fix */

  a=chain[k];
  k1 = b + sigsize * a;

  /* look for actual end of the chain */
  while(mridge[k1] < threshold) {
    k -= nbchain;
    a = chain[k];
    b--;
    k1 = b + sigsize * a;
  }
  aend=a;
  bend=b;

  /* shift */
  oldbstart = chain[(*id)-1];
  chain[(*id)-1]=bstart;
  lng = bend-bstart+1;
  if(lng <= bstep) {
    (*id) --;
    return;
  }
  b = (bstart-oldbstart);
  for(count = 1; count < lng; count++){
    b ++;
    k = (*id)-1 + b * nbchain;
    chain[(*id)-1 + count *  nbchain] = chain[k];
  }
  b++;
  k = (*id) - 1 + count * nbchain;
  while((b < (sigsize+bstart-oldbstart)) && ((int)(chain[k]) != -1)) {
    chain[k] = -1;
    b++;
    k += nbchain;
  }

  return;
}

/****************************************************************
*   Function: reordering
*   ---------
*      Organizes chain[] according to the chosen data structure
*****************************************************************/

void reordering(int *chain,int sigsize,int nbchain)
{
  int i,j,cnt;
  for(i=0; i<nbchain-1;i++){
    j = sigsize; cnt=0;
    while((j>0)&&(chain[i+j*nbchain]==-1)) j--;
    while((j>0)&&(chain[i+j*nbchain]!=-1)){
      chain[i+(j+1)*nbchain]=chain[i+j *nbchain];
      j--;
      cnt++;
    }
    chain[i+nbchain] = cnt;
  }
  return;
}

