
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
*  Function: Spointmap:
*  --------------------
*  Road map designed using principle component analysis
*
*    smodulus: squared-modulus of the wavelet transform
*     (coming from learning the noise)
*    sigsize: signal size
*    nscale: total number of scales for CWT
*    gridx : the window size for pca along x
*    gridy : the window size for pac along y
*    nbblock: number of window used by pca
*    nbpoint: number of sampler each window by pca
*    pointmap: map of sampling points
*    tst: the first 4 locations containing the coordinates of left
*         and right corner, followed with a sequence of the coor-
*         dinates of sampled points in each block.
*    tstsize: equal to 4 + 2 * number of sampled point each pca block
*    count: the maximal number of repetition if location is chosen in 
*           a block
*    seed: seed for random number generator
****************************************************************/

void Spointmap(double *sqmodulus, int *psigsize, int *pnscale,
	       int *pgridx, int *pgridy, int *pnbblock, int *pnbpoint, 
	       int *pointmap, float *tst, int *ptstsize, int *pcount, int *pseed)
{
  long l;
  int sigsize, nscale, i, j, k, tstsize;
  int gridx, gridy, nbpoint;
  long u1, u2, u3;
  double um, pm, largest;
  long seed;
  float p1, p2, p3;
  int a, b, t, up, down, left, right;
  int count1, count2, nbblock, bnumber, lnb;
  int lefto, righto, upo, downo;

  /* seed for random number
     ---------------------- */
  seed = (*pseed);

  sigsize = *psigsize;
  nscale = *pnscale;


  /* size of grid along b and a
     -------------------- */
  gridx = *pgridx;
  gridy = *pgridy;
  tstsize = *ptstsize;
/*
  for(l = 0; l < sigsize * nscale; l++)
    pointmap[l] = 0;
*/
  /* number of points to be in a block  and the total number of block
     --------------------------------- */
  nbblock = *pnbblock;
  nbpoint = *pnbpoint;

  bnumber = 0;

  for(a = 0; a < nscale; a+= gridy) {
    down = a;
    up = min(nscale-1, a + gridy);

    /* overlapping with neighborhood block 
       -----------------------------------*/
    downo = max(0,a - gridy/2);
    upo = min(nscale-1, a + gridy + gridy/2);

    for(b = 0; b  < sigsize; b +=  gridx) {

      left = b;
      right = min(sigsize-1,b + gridx);

      /* overlapping with neighborhood block 
	 -----------------------------------*/
      lefto = max(0,b - gridx/2);
      righto = min(sigsize-1, b + gridx + gridx/2);

      /* largest sqared modulus in a block 
         --------------------------------- */
      largest = 0.0;

/*      for(i = left; i < right; i++) {
	for(j = down; j < up; j++) {
	  k = j * sigsize + i;
	  largest = max(sqmodulus[k], largest);
	}
      }
*/

      for(i = lefto; i < righto; i++) {
	for(j = downo; j < upo; j++) {
	  k = j * sigsize + i;
	  pointmap[k] = 0;
	  largest = max(sqmodulus[k], largest);
	}
      }

      /* locations of the block; Add 1 to compensate difference between S and c
         ---------------------------------------------------------------------- */
      tst[tstsize*bnumber] = (float)(left + 1.0);
      tst[tstsize*bnumber+1] = (float)(down + 1.0);
      tst[tstsize*bnumber+2] = (float)(right + 1.0);
      tst[tstsize*bnumber+3] = (float)(up + 1.0);

      for(t = 1; t <= nbpoint; t++) {

	count2 = 0;
	while(1) {
/*
	  p1 = (float)(ran1(&seed));
	  u1 = (long)min((right-left-1) * p1 + left, sigsize-1);
	  p2 = (float)(ran1(&seed));
	  u2 = (long)min((up-down-1) * p2 + down, nscale-1);
*/

	  p1 = (float)(ran1(&seed));
	  u1 = (long)min((righto-lefto-1) * p1 + lefto, sigsize-1);
	  p2 = (float)(ran1(&seed));
	  u2 = (long)min((upo-downo-1) * p2 + downo, nscale-1);

	  count1 = 0;
	  while(1) {

	    k = u1 + u2 * sigsize;
            /* does the position chosen ? 
               -------------------------- */
	    if((pointmap[k] == 0) || (count1 > (*pcount))) break;
/*	    
	    p1 = (float)(ran1(&seed));
	    u1 = (long)((right-left-1) * p1 + left);
	    u1 = min(u1, sigsize-1);
	    p2 = (float)(ran1(&seed));
	    u2 = (long)((up-down-1) * p2 + down);
	    u2 = min(u2, nscale-1);
*/

	    p1 = (float)(ran1(&seed));
	    u1 = (long)((righto-lefto-1) * p1 + lefto);
	    u1 = min(u1, sigsize-1);
	    p2 = (float)(ran1(&seed));
	    u2 = (long)((upo-downo-1) * p2 + downo);
	    u2 = min(u2, nscale-1);

	    count1 ++; 
	  }

	  k = u2 * sigsize + u1;
	  pm = ran1(&seed);
	  um = pm * largest;

	  if((sqmodulus[k] >= um) || (count2 > (*pcount))) {
	    
	    /* choose only when modulus at (u1, u2) is at least *
	     * as large as a random reference.                  *
               ------------------------------------------------ */
	    pointmap[k] = 1;

            /* the locations of selected points, normalized as prob. 
               --------------------------------------------------- */
            tst[tstsize*bnumber+ 2*(t+1)] =  p1;
            tst[tstsize*bnumber+ 2*(t+1)+1] = p2;

 	    break;
	  }
	  count2 ++;
	}
      }
      /* advance to next block 
	 --------------------- */
      bnumber ++;
    }
  }
}


/****************************************************************
*  Function: Spca_annealing:
*  --------------------------
*  Ridges characterisation with pca climber algorithm
*
*    smodulus: squared-modulus of the wavelet transform
*     (coming from learning the noise)
*    beemap: ridges
*    crazymap: roadmap for pca
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

void Spca_annealing(double *smodulus, double *beemap,
               int *crazymap,
	       float *pc,
	       int *psigsize, int *pnscale, int *piteration,
               int *pseed, int *pbstep, int *pnbbee,
               int *pintegral, int *pchain, int *flag)
{
  double r, dd, ee;
  int i, bstep, k, k1, k2, bee;
  int *a, *b, integral, chain, s, count, dir;
  int arest, brest, afree, bfree;
  int seed, nscale, iteration, sigsize, nbbee;
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
    if(integral) 
      beemap[k] = beemap[k] + smodulus[k];
    else 
      beemap[k] = beemap[k] +1;

    /* Iterations */
    for(i =1; i < iteration; i++) {
      
      dir = crazymap[k];

      /* free along a or along b
	 -----------------------*/

      if((dir == 1) || (dir == 3)) {  

	if(dir == 1) { /* more freely along b */
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
	}

	if(dir == 3) { /* mover freely along b */
	  r = (double)ran1(&idum);
	  if (r >= 0.5) 
	    b[i] = min(sigsize-1,b[i-1]+1);
	  else 
	    b[i] = max(0,b[i-1]-1);

	  r = (double)ran1(&idum);

	  if (r >= 0.5) 
	    a[i] = min(nscale-1,a[i-1]+bstep);
	  else 
	    a[i] = max(0,a[i-1]-bstep);
	  
	  k1 = b[i] + sigsize * a[i];
	  k2 = b[i-1] + sigsize * a[i];
	}

	dd = smodulus[k1] - smodulus[k2];

	if (dd<0) {
	  r = (double)ran1(&idum);
	  ee = exp(dd * log((double)(3.0 + i))/c);
	  if((*flag)==1) ee=exp(dd * log((double)(3.0))/c); /* constant temperature */

	  /* freeze 
	     ------ */
	  if (r>ee) {
	    if(dir == 1) a[i] = a[i-1]; 
	    else b[i] = b[i-1];
	  }
	}
      }
 
      /* free along a=b or a=-b 
	 ----------------------*/

      if((dir == 2) || (dir == 4)) { 

	if(dir == 2) { /* move freely along a = -b */
	  r = (double)ran1(&idum);

	  /* choose free move 
	     ----------------*/
	  if (r >= 0.5) {
	    bfree = min(sigsize-1,b[i-1]+bstep);
	    afree = max(0,a[i-1]-bstep);
	  }
	  else {
	    bfree = max(0,b[i-1]-bstep);
	    afree = min(nscale-1,a[i-1]+bstep);
	  }
	  /* choose restricted move 
	     ----------------------*/
	  r = (double)ran1(&idum);
	  if (r >= 0.5) {
	    arest = max(0,afree - 1);
	    brest = max(0,bfree - 1);
	  }
	  else {
	    arest = min(nscale-1,afree + 1);
	    brest = min(sigsize-1,bfree + 1);
	  }
	}

	if(dir == 4) { /* more freely along a = b */
	  r = (double)ran1(&idum);

	  /* choose free move 
	     ----------------*/
	  if (r >= 0.5) {
	    bfree = min(sigsize-1,b[i-1]+bstep);
	    afree = min(nscale-1,a[i-1]+bstep);
	  }
	  else {
	    bfree = max(0,b[i-1]-bstep);
	    afree = max(0,a[i-1]-bstep);
	  }

	  /* choose restricted move 
	     ----------------------*/
	  r = (double)ran1(&idum);
	  if (r >= 0.5) {
	    arest = min(nscale-1,afree+1);
	    brest = max(0,bfree-1);
	  }
	  else {
	    arest = max(0,afree-1);
	    brest = min(sigsize-1,bfree+1);
	  }
	}

	k1 = brest + arest * sigsize;
	k2 = bfree + afree * sigsize;

	dd = smodulus[k1] - smodulus[k2];
	b[i] = brest;
	a[i] = arest;
	if (dd<0) {
	  r = (double)ran1(&idum);
	  
	  if((*flag) == 1)
	    ee=exp(dd * log((double)(3.0))/c); 
	  else 
	    ee = exp(dd * log((double)(3.0 + i))/c);

	  /* freeze 
	     ------ */
	  if (r>ee) {
	    a[i] = afree;
	    b[i] = bfree;
	  }
	}
      }

      /* Chaining a curve between (a[i-1],b[i-1]) and (a[i],b[i]) with larger modulus 
	 ---------------------------------------------------------------------------- */

      if(chain) {

	if(dir == 1) { 

	  count = 1;
	  while(count < bstep) {
	    if(b[i] - b[i-1] > 0) {
	      if(b[i-1] + count >= sigsize) break;

	      k1 = b[i-1] + count + sigsize * a[i];
	      k2 = b[i-1] + count + sigsize * a[i-1];
	      if(smodulus[k1] > smodulus[k2])
		k = k1;
	      else 
		k = k2;
	    }
	    else {
	      if(b[i-1] - count < 0) break;

	      k1 = b[i-1] - count + sigsize * a[i];
	      k2 = b[i-1] - count + sigsize * a[i-1];
	      if(smodulus[k1] > smodulus[k2]) 
		k = k1;
	      else 
		k = k2;
	    }

	    /* update of the ridges */
	    if(integral) 
	      beemap[k] = beemap[k] + smodulus[k];
	    else 
	      beemap[k] = beemap[k] +1;
	    
	    count ++;
	  }
	}


	if(dir == 3) { 

	  count = 1;
	  while(count < bstep) {
	    if(a[i] - a[i-1] > 0) {

	      if(a[i-1] + count >= nscale) break;

	      k1 = b[i] + sigsize * (a[i-1] + count);
	      k2 = b[i-1] + sigsize * (a[i-1] + count);
	      if(smodulus[k1] > smodulus[k2])
		k = k1;
	      else 
		k = k2;
	    }
	    else {
	      if((a[i-1] - count) < 0 ) break;

	      k1 = b[i] + sigsize * (a[i-1] - count);
	      k2 = b[i-1] + sigsize * (a[i-1] - count);
	      if(smodulus[k1] > smodulus[k2]) 
		k = k1;
	      else 
		k = k2;
	    }

	    /* update of the ridges */
	    if(integral) 
	      beemap[k] = beemap[k] + smodulus[k];
	    else 
	      beemap[k] = beemap[k] +1;
	    
	    count ++;
	  }
	}


	if(dir == 2) { 
	  if((a[i] < a[i-1]) || (b[i] > b[i-1])) {

	    count = 1;
	    while(count < bstep) {

	      if( ((b[i-1] + count + 1) >= sigsize) ||
		  ((a[i-1] - count - 1) < 0))
		break;

	      if((bfree == brest) && (afree == arest)) {
		k = b[i-1] + count + sigsize * (a[i-1] - count);
	      }
	      if(arest > afree) {

		k1 = (b[i-1] + count) + sigsize * (a[i-1] - count);
		k2 = (b[i-1] + count + 1) + sigsize * (a[i-1] - count + 1);
		if(smodulus[k1] > smodulus[k2]) 
		  k = k1;
		else 
		  k = k2;
	      }
	      if(arest < afree) {
		k1 = (b[i-1] + count) + sigsize * (a[i-1] - count);
		k2 = (b[i-1] + count - 1) + sigsize * (a[i-1] - count - 1);
		if(smodulus[k1] > smodulus[k2]) 
		  k = k1;
		else 
		  k = k2;
	      }
	    
	      /* update of the ridges */
	      if(integral) 
		beemap[k] = beemap[k] + smodulus[k];
	      else 
		beemap[k] = beemap[k] +1;
	      
	      count ++;
	    }
	  }
	  if((a[i] > a[i-1]) || (b[i] < b[i-1])) {
	    
	    count = 1;
	    while(count < bstep) {

	      if( ((b[i-1] - count - 1) < 0) ||
		  ((a[i-1] + count + 1) >= nscale) ||
		  ((a[i-1] - count - 1) < 0)) 
		break;

	      if((bfree == brest) && (afree == arest)) {
		k = b[i-1] - count + sigsize * (a[i-1] + count);
	      }
	      if(arest > afree) {
		k1 = (b[i-1] - count) + sigsize * (a[i-1] + count);
		k2 = (b[i-1] - count + 1) + sigsize * (a[i-1] + count + 1);
		if(smodulus[k1] > smodulus[k2]) 
		  k = k1;
		else 
		  k = k2;
	      }
	      if(arest < afree) {
		k1 = (b[i-1] - count) + sigsize * (a[i-1] + count);
		k2 = (b[i-1] - count - 1) + sigsize * (a[i-1] - count - 1);
		if(smodulus[k1] > smodulus[k2]) 
		  k = k1;
		else 
		  k = k2;
	      }
	    
	      /* update of the ridges */
	      if(integral) 
		beemap[k] = beemap[k] + smodulus[k];
	      else 
		beemap[k] = beemap[k] +1;
	      
	      count ++;
	    }
	  }
	}


	if(dir == 4) {  /* free along a = b */
	  if((a[i] > a[i-1]) || (b[i] > b[i-1])) {

	    count = 1;
	    while(count < bstep) {

	      if( ((b[i-1] + count + 1) >= sigsize) ||
		  (a[i-1] + count + 1>= nscale))
		break;


	      if((bfree == brest) && (afree == arest)) {
		k = b[i-1] + count + sigsize * (a[i-1] + count);
	      }
	      if(arest > afree) {
		k1 = (b[i-1] + count) + sigsize * (a[i-1] + count);
		k2 = (b[i-1] + count - 1) + sigsize * (a[i-1] + count + 1);
		if(smodulus[k1] > smodulus[k2]) 
		  k = k1;
		else 
		  k = k2;
	      }
	      if(arest < afree) {
		k1 = (b[i-1] + count) + sigsize * (a[i-1] + count);
		k2 = (b[i-1] + count + 1) + sigsize * (a[i-1] + count - 1);
		if(smodulus[k1] > smodulus[k2]) 
		  k = k1;
		else 
		  k = k2;
	      }
	    
	      /* update of the ridges */
	      if(integral) 
		beemap[k] = beemap[k] + smodulus[k];
	      else 
		beemap[k] = beemap[k] +1;
	      
	      count ++;
	    }
	  }
	  if((a[i] < a[i-1]) || (b[i] < b[i-1])) {
	    
	    count = 1;
	    while(count < bstep) {

	      if( ((b[i-1] - count - 1) < 0) ||
		  (a[i-1] - count - 1 < 0))
		break;

	      if((bfree == brest) && (afree == arest)) {
		k = b[i-1] - count + sigsize * (a[i-1] - count);
	      }
	      if(arest > afree) {
		k1 = (b[i-1] - count) + sigsize * (a[i-1] - count);
		k2 = (b[i-1] - count - 1) + sigsize * (a[i-1] - count + 1);
		if(smodulus[k1] > smodulus[k2]) 
		  k = k1;
		else 
		  k = k2;
	      }
	      if(arest < afree) {
		k1 = (b[i-1] - count) + sigsize * (a[i-1] - count);
		k2 = (b[i-1] - count + 1) + sigsize * (a[i-1] - count - 1);
		if(smodulus[k1] > smodulus[k2]) 
		  k = k1;
		else 
		  k = k2;
	      }
	    
	      /* update of the ridges */
	      if(integral) 
		beemap[k] = beemap[k] + smodulus[k];
	      else 
		beemap[k] = beemap[k] +1;
	      
	      count ++;
	    }
	  }
	}
      }

      /* advance to next postion 
	 -----------------------*/

      k = b[i] + a[i] * sigsize;
      if(integral) 
	beemap[k] = beemap[k] + smodulus[k];
      else 
	beemap[k] = beemap[k] +1;
    }
  }
}


       



    

  


    
    
  
  

  


       



    

  


    
    
  
  

  




    

  


    
    
  
  

  

