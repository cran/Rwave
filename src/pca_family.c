
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
*  function Stf_pcaridge
*    compute the local maxima along the minor direction of pca
*
*  input: tf
*  output: tf local maxima
*  orientmap : orientation of pca 
*  nrow,ncol: parameters of the cwt 
*****************************************************************/

Stf_pcaridge(input, output, pnrow, pncol, orientmap)
     double *input, *output;
     int *pnrow, *pncol;
     int *orientmap;
{
  int nrow, ncol, i, j,k;
  int pos, dir;
  float tmp;

  nrow = *pnrow;
  ncol = *pncol;
  for(i = 1; i < nrow-1; i++) {
    for(j = 1; j < ncol-1; j++) {
      k = j * nrow + i;
      dir = orientmap[k];
      if(dir == 1) {
	/* minor direction at dir == 3 */
	if(((input[k] > input[(j+1) * nrow + i]) &&
	    (input[k] >= input[(j-1) * nrow + i])) ||
	   ((input[k] > input[(j-1) * nrow + i]) &&
	    (input[k] >= input[(j+1) * nrow + i])))
	  output[k] = input[k];
      }
      if(dir == 3) {
	/* minor direction at dir == 1 */
	if(((input[k] > input[j * nrow + (i+1)]) &&
	    (input[k] >= input[j * nrow + (i-1)])) ||
	   ((input[k] > input[j * nrow + (i-1)]) &&
	    (input[k] >= input[j * nrow + (i+1)])))
	  output[k] = input[k];
      }
      if(dir == 2) {
	/* minor direction at dir == 4 */
	if(((input[k] > input[(j+1) * nrow + (i-1)]) &&
	    (input[k] >= input[(j-1) * nrow + (i+1)])) ||
	   ((input[k] > input[(j-1) * nrow + (i+1)]) &&
	    (input[k] >= input[(j+1) * nrow + (i-1)])))
	  output[k] = input[k];
      }
      if(dir == 4) {
	/* minor direction at dir == 2 */
	if(((input[k] > input[(j+1) * nrow + (i+1)]) &&
	    (input[k] >= input[(j-1) * nrow + (i-1)])) ||
	   ((input[k] > input[(j-1) * nrow + (i-1)]) &&
	    (input[k] >= input[(j+1) * nrow + (i+1)])))
	  output[k] = input[k];
      }
    }
  }
}



/*******************************************************************
 * dir == 1 along x, dir == 2, along x=y, dir == 3 along y, dir ==4 
 * along x=-y 
 *******************************************************************/
/* Three majors direction of dir going downwards and left */

previous_a_b(int a,int b,int dir,
	     int *a0,int *b0,int *a1,int *b1,int *a2,int *b2) 
{
  if(dir == 1) {
    *a0 = a;
    *b0 = b-1;
    *a1 = a-1;
    *b1 = b-1;
    *a2 = a+1;
    *b2 = b-1;
  }
  if(dir == 3) {
    *a0 = a-1;
    *b0 = b;
    *a1 = a-1;
    *b1 = b-1;
    *a2 = a-1;
    *b2 = b+1;
  }
  if(dir == 2) {
    *a0 = a-1;
    *b0 = b-1;
    *a1 = a-1;
    *b1 = b;
    *a2 = a;
    *b2 = b-1;
  }
  if(dir == 4) {
    *a0 = a-1;
    *b0 = b+1;
    *a1 = a-1;
    *b1 = b;
    *a2 = a;
    *b2 = b+1;
  }
}

/* Three major directions of dir going upwards or right */
next_a_b(int a,int b,int dir,
	     int *a0,int *b0,int *a1,int *b1,int *a2,int *b2) 
{
  if(dir == 1) {
    *a0 = a;
    *b0 = b+1;
    *a1 = a-1;
    *b1 = b+1;
    *a2 = a+1;
    *b2 = b+1;
  }
  if(dir == 3) {
    *a0 = a+1;
    *b0 = b;
    *a1 = a+1;
    *b1 = b-1;
    *a2 = a+1;
    *b2 = b+1;
  }
  if(dir == 2) {
    *a0 = a+1;
    *b0 = b+1;
    *a1 = a+1;
    *b1 = b;
    *a2 = a;
    *b2 = b+1;
  }
  if(dir == 4) {
    *a0 = a+1;
    *b0 = b-1;
    *a1 = a+1;
    *b1 = b;
    *a2 = a;
    *b2 = b-1;
  }
}












/****************************************************************
*   Function: pca_orderedmap_thresholded
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

void pca_orderedmap_thresholded(float *orderedmap,int sigsize,int nscale,
			    int *chain,int nbchain)
{
  int id, i, j, k;
  int a, b, count, chnlng;
  
  for(i = 0; i < sigsize; i++)
    for(j = 0; j < nscale; j++) {
      k = i + j * sigsize;
      orderedmap[k] = 0.0;
    }

  for(id = 0; id < nbchain; id++) {
    k = id;
    chnlng = chain[k];
    k += nbchain;
    a = chain[k];
    k += nbchain;
    b = chain[k];
    count = 1;
    while((count <= chnlng)) {
      orderedmap[b + a * sigsize] = (float)(id + 1);
      k += nbchain;
      a = chain[k];
      k += nbchain;
      b = chain[k];
      count++;
    }
  }
  return;
}
    

/****************************************************************
*   Function: pca_chain_thresholded
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

void pca_chain_thresholded(double *mridge,int sigsize,int *chain,int *id,
		       int nbchain,float threshold, int bstep)
{
  int i, j, k, k1, a, b;
  int count, found, lng;
  int kstart, kend;
  int chnlng;



  /* look for actual beginning of the chain */

  k = (*id)-1;
  chnlng = chain[k];
  k += nbchain;
  kstart = k; /* beginning of the chain */
  a = chain[k];
  k +=  nbchain;
  b = chain[k];
  k1 = b + sigsize * a;
  count = 1;
  while((count <= chnlng)&&(mridge[k1]<(double)threshold)) {
    k += nbchain;
    kstart = k;
    a = chain[k];
    k += nbchain;
    b = chain[k];
    k1 = b + sigsize * a;
    count++;
  }
  /* Invalid chain (not enough energy) */
  if(count > chnlng) {
    k = (*id)-1;
    chain[k] = -1;
    count = 0;
    while(count <= chnlng) {
      k += nbchain;
      chain[k] = -1;
      k += nbchain;
      chain[k] = -1;
      count++;
    }
    (*id)--;
    return;
  }


  /* move to the end of the chain */
  while(count < chnlng){
    k+= nbchain;
    k+= nbchain;
    count ++;
  }

  kend = k;
  b = chain[k];
  k-= nbchain;
  a=chain[k];
  k1 = b + sigsize * a;
  /*   look for actual end of the chain  */
  while(mridge[k1]<threshold) {
    k -= nbchain;
    kend = k;
    b = chain[k];
    k -= nbchain;
    a = chain[k];
    k1 = b + sigsize * a;
  }
  /* shift the chain */
  count = 1;
  k = kstart;
  chain[(*id)-1+count*nbchain]=chain[k];
  while(k != kend) {
    count++;
    k += nbchain;
    chain[(*id)-1+count*nbchain]=chain[k];
  }
  count++;
  k += nbchain;
  chain[(*id)-1+count*nbchain]=chain[k];
  chain[(*id)-1] = count/2;/* the size of the chain */

  /* remove the chain that are shorter than bstep */
  if(chain[(*id)-1] < bstep) {
    chnlng = chain[(*id)-1];
    k = (*id)-1;
    chain[k] = -1;
    count = 0;
    while(count <= chnlng) {
      k += nbchain;
      chain[k] = -1;
      k += nbchain;
      chain[k] = -1;
      count++;
    }
    (*id)--;

    return;
  }

  return;
}


    


/****************************************************************
*   Function: Spca_family
*   ---------
*      Organize the output of the pca climber algorithm into
*      a series of connected ridges.
*
*   ridgemap: output of the crazy climber algorithm
*   orientmap: direction of 1st eigen vector 
*   sigsize: input signal size
*   nscale: number of scales (or frequencies) of the transform
*   nbchain: number of ridges
*   chain: array containing the ridges; the first element
*          contains the length of the ridge, the following in
*          order of (a,b). The maximal length of chain is given.
*   id: pointer containing the indes of the considered ridge
*   threshold: minimal value of transform modulus considered
*              on a ridge
*   maxchnlng: maximal chain length
*****************************************************************/

void Spca_family(double *ridgemap,int *orientmap, float *orderedmap,int *chain,
		   int *pnbchain, int *psigsize,int *pnscale,
		   int *pbstep,float *pthreshold, int* pmaxchnlng)
{
  int bstep, nscale, sigsize, nbchain;
  int i, j, k, id, count, a, b, found, k1, maxchnlng;
  int a0,a1,a2,b0,b1,b2;
  double *mridge;
  float threshold;
  int dir;

  threshold = *pthreshold;
  bstep = *pbstep;
  nscale = *pnscale;
  sigsize = *psigsize;
  nbchain = *pnbchain;
  maxchnlng = *pmaxchnlng;

  if(!(mridge = (double *)calloc(sigsize * nscale, sizeof(double))))
    error("Memory allocation failed for mridge in crazy_family.c \n");


  /* Compute local maxima of modulus (along minor axis of pca
     ------------------------------------------------------- */
  Stf_pcaridge(ridgemap,mridge,psigsize,pnscale,orientmap);  


  id = 0;

  /* Start looking for ridges as connected curves
     -------------------------------------------- */
  for(i = 0; i < sigsize; i+= bstep) {
    for(j = 0; j < nscale; j++) {
      b = i;
      a = j;
      k = b + sigsize * a;
      dir = orientmap[k];

      if((mridge[k] > 0.000001) && (orderedmap[k] == 0.0)) {
	found = YES;

	/* backwarding: looking for previous points of the chain 
	   ----------------------------------------------------- */

        while(found){
	  found = NO;
          previous_a_b(a,b,dir,&a0,&b0,&a1,&b1,&a2,&b2);
          if(inrange(0,a0,nscale-1) &&inrange(0,b0,sigsize-1)) {
            k1 = b0 + sigsize * a0;
            dir = orientmap[k1];
	    if((mridge[k1]>0.000001)&&(orderedmap[k1]==0.0)) {
	      found = YES;	    
              b = b0;
              a = a0;
	    }
	  }
	  /*          if (found == NO) {
	    if(inrange(0,a1,nscale-1) &&inrange(0,b1,sigsize-1)) {
	      k1 = b1 + sigsize * a1;
	      dir = orientmap[k1];
	      if((mridge[k1]>0.000001)&&(orderedmap[k1]==0.0)) {
		found = YES;	    
                b = b1;
                a = a1;
	      }
	    }
	  }
          if(found == NO) {
	    if(inrange(0,a2,nscale-1) &&inrange(0,b2,sigsize-1)) {
	      k1 = b2 + sigsize * a2;
	      dir = orientmap[k1];
	      if((mridge[k1]>0.000001)&&(orderedmap[k1]==0.0)) {
		found = YES;	    
		b= b2;
		a= a2;
	      }
	    }
	  } */
	}

	/* forwarding
	   ---------- */
	id ++;
	if(id >= nbchain) {
	  printf("Nb of chains > reserved number %d. Returned. \n",nbchain); 
	  return;
	}
	count = 1;
	found = YES;
	while(found) {
	  chain[(id-1)+count*nbchain] = a;
          count++;
          if(count/2 > maxchnlng) 
	    error("Longer than max chain length. Returned. \n");
	  chain[(id-1)+count*nbchain]= b;
          count++;

	  k= b + sigsize * a;
	  dir = orientmap[k];
	  next_a_b(a,b,dir,&a0,&b0,&a1,&b1,&a2,&b2);
	  orderedmap[k] =(float)(id);
	  found = NO;
          if(inrange(0,a0,nscale-1) && inrange(0,b0,sigsize-1)) {
	    k1 = b0 + sigsize * a0;
	    if((mridge[k1]>0.000001)&&(orderedmap[k1]==0.0)) {
	      found = YES;
	      orderedmap[k1] =(float)(id);
	      a = a0;
	      b = b0;
	    }
	  }
          if(found == NO) {
	    if(inrange(0,a1,nscale-1) && inrange(0,b1,sigsize-1)) {
	      k1 = b1 + sigsize * a1;
	      if((mridge[k1]>0.000001)&&(orderedmap[k1]==0.0)) {
		found = YES;
		orderedmap[k1] =(float)(id);
		a = a1;
		b = b1;
	      }
	    }
	  }
	  if(found == NO) {
	      if(inrange(0,a2,nscale-1) && inrange(0,b2,sigsize-1)) {
		k1 = b2 + sigsize * a2;
		if((mridge[k1]>0.000001)&&(orderedmap[k1]==0.0)) {
		  found = YES;
		  orderedmap[k1] =(float)(id);
		  a = a2;
		  b = b2;
		}
	      }
	  }
	}
	chain[(id-1)]= (count-1)/2;
	/*	printf("number of chain %d with lng %d \n",id,chain[(id-1)]); */
	/* Threshold and chain the ridges
	   ------------------------------ */
	pca_chain_thresholded(mridge,sigsize,chain,&id,
		  nbchain,threshold,bstep);
      }
    }
  }

  /* Generate the image of the chains
     -------------------------------- */
  pca_orderedmap_thresholded(orderedmap,sigsize,nscale,chain,nbchain); 


  free((char *)mridge);
  printf("There are %d chains. \n", id);
  *pnbchain = id;

  return;
}




