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
*  Function: Sridge_newicm:
*  --------------------------
*  Ridge characterization with Besag's ICM algorithm
*
*    smodulus: smoothed modulus of the wavelet transform
*    cost:   cost function
*    phi: ridge
*    lambda: coefficient in front of phi' in the cost function
*    mu: coefficient in front of phi'' in the cost function
*    sigsize: signal size
*    nscale: total number of scales for CWT
*    iteration: maximal number of iterations for the annealing
*    stagnant: allowed number of consecutive steps without
*              move (stopping criterion)
*    count: number of iterations
*    sub: subsampling rate for ridge extraction
*    smodsize: the size of sub-sampled signal
*
****************************************************************/

void Sridge_icm(double *cost, double *smodulus, double *phi,
  double *plambda, double *pmu, int *psigsize, int *pnscale,
  int *piteration,int *pcount, int *psub, int *psmodsize)
{
  int i,sigsize,iteration,up, best_up,pos,a,count,sub;
  int smodsize, tbox=0, ttbox=1000;
  int nscale;
  double lambda, mu;
  double *phi2;
   double cost1;
  double tmp=0.0, best_tmp;


  /* Generalities; initializations
     -----------------------------*/
  mu = *pmu;
  nscale = *pnscale;
  iteration = *piteration;
  lambda = *plambda;
  sigsize = *psigsize;
  sub = *psub;
  smodsize = *psmodsize;

  if(!(phi2 = (double *)S_alloc((smodsize+1)*sub,sizeof(double))))
    error("Memory allocation failed for phi2 at icm.c \n");

  count = 0; /* total count */
  cost1 = 0;

  for(i=0;i<smodsize;i++){
    phi[i] = phi[sub*i];
  }
/*  printf("smodsize=%d\n",smodsize);  */



  /* Iterations:
     -----------*/
  while((ttbox > 1)&&(count < iteration)) {

    /* Initialize the cost function
     ----------------------------*/
    tbox = 0;

    if(count == 0) {
      for(i = 1; i < smodsize-1; i++) {
        tmp = (double)((phi[i-1]+ phi[i+1]-2 * phi[i]));
        cost1 += (double)((lambda * tmp * tmp));
	
        tmp = (double)((phi[i] - phi[i+1]));
        cost1 += (double)((mu * tmp * tmp));
	
        a = (int)phi[i];
        tmp = smodulus[smodsize * a + i];
/*      cost1 -= (tmp * tmp - noise[a]); */
        cost1 -= tmp;
      }
      
      tmp = (double)((phi[0] - phi[1]));
      cost1 += (double) ((mu * tmp * tmp));
      a = (int)phi[0];
      tmp = smodulus[smodsize * a];
/*      cost1 -= (tmp * tmp - noise[a]); */
      cost1 -= tmp;
      
      a = (int)phi[smodsize-1];
      tmp = smodulus[smodsize * a + smodsize-1];
/*      cost1 -= (tmp * tmp - noise[a]); */
      cost1 -= tmp;
      
    }				   


    /* Generate moves
       --------------*/
    for(pos=0; pos < smodsize; pos++){

      best_tmp = (double)0.0; 
      best_up = 0;

      for(up = -(int)phi[pos]; up < nscale- (int)phi[pos]; up++){

        /* Compute corresponding update of the cost function
           -------------------------------------------------*/
        if(inrange(2,pos,smodsize-3)) {
          tmp = (double)(lambda*up);
          tmp *=(double)((6*up+(12*phi[pos]
            -8*(phi[pos-1]+phi[pos+1])
            +2*(phi[pos-2]+phi[pos+2]))));

          tmp += (double)(mu*up*(4.0*phi[pos]
            -2.0*(phi[pos-1]+phi[pos+1])+2.0*up));

          a = (int)phi[pos];
/*          tmp += (smodulus[smodsize*a+pos]
            *smodulus[smodsize*a+pos]);
          tmp -= noise[a]; */
          tmp += smodulus[smodsize*a+pos];

          a = (int)phi[pos] +  up;
/*          tmp -= (smodulus[smodsize*a+pos]
            *smodulus[smodsize*a+pos]);
          tmp += noise[a]; */
	  tmp -= smodulus[smodsize*a+pos];
        }

        if(inrange(2,pos,smodsize-3) == NO) {
          tmp = (double)(lambda*up);
          if(pos == 0) {
            tmp *= (double)((up+2.0*(phi[0]
              -2*phi[1]+phi[2])));
            tmp += (double)(mu*up*((2.0*phi[pos]
              -2.0*phi[pos+1])+up));
          }
          else if(pos == 1) {
            tmp *=(double)((5*up+2.0*(-2*phi[0]+5*phi[1]
              -4*phi[2]+phi[3])));
            tmp += (double)(mu*up*(4.0*phi[pos]
              -2.0*(phi[pos-1]+phi[pos+1]-up)));
          }
	  else if(pos == (smodsize-2)) {
            tmp *= (double)((5*up+2.0*(phi[pos-2]-4*phi[pos-1]
              +5*phi[pos]-2*phi[pos+1])));
            tmp += (double)(mu*up*(4.0*phi[pos]-2.0*(phi[pos-1]
              +phi[pos+1])+2.0*up));
	  }
	  else if(pos == (smodsize-1)) {
            tmp *= (double)((up+2.0*(phi[pos-2]
              -2*phi[pos-1]+phi[pos])));
            tmp += (double)(mu*up*((2.0*phi[pos]
              -2.0*phi[pos-1])+up));
          }

          a = (int)phi[pos];
/*          tmp +=(smodulus[smodsize*a+pos]*smodulus[smodsize*a+pos]);
          tmp -= noise[a]; */
          tmp +=smodulus[smodsize*a+pos];

          a = (int)phi[pos] +  up;
/*          tmp -=(smodulus[smodsize*a+pos]*smodulus[smodsize*a+pos]);
          tmp += noise[a]; */
          tmp -=smodulus[smodsize*a+pos];
	}


        /* Compare with other moves
	   ------------------------*/
        if(tmp < best_tmp) {
	  best_tmp = tmp;
          best_up = up;
        }
      }

      /* Best move
         ---------*/
      if (best_up != 0) {
        cost1 += best_tmp;
        phi[pos] += (double)best_up;
        tbox++;
      }
    }
    ttbox = tbox;
    cost[count++] = cost1;
  }


  /* Interpolate from subsampled ridge
     --------------------------------*/
  if (sub != 1){
  splridge(sub, phi, smodsize, phi2);
  for(i=0;i<sigsize;i++)
    phi[i]=phi2[i];
}
  *pcount = count;

  return;
}

