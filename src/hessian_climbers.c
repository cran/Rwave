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
*  Function: Shessianmap:
*  ---------------------
*  Road map designed using hessian matrix
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

void Shessianmap(double *sqmodulus, int *psigsize, int *pnscale, int *pnbblock,
		 int *pgridx, int *pgridy, double *tst)
{
  int a, b, sigsize, nscale, gridx, gridy, nbblock;
  int left, right, down, up;
  double mxx, mxy, myx, myy;
  int bnumber, k;

  sigsize = *psigsize;
  nscale = *pnscale;
  gridx = *pgridx;
  gridy = *pgridy;

  bnumber = 0;

  for(a = 2; a < nscale-2; a += gridy) {
    for(b = 2; b  < sigsize-2; b += gridx) {

      down = a;
      up = min(a + gridy, nscale-1);
      left = b;
      right = min(b + gridx, sigsize-1);

      k = b + a * sigsize;
      mxx = 0.25*(sqmodulus[k+2]+sqmodulus[k-2]-2*sqmodulus[k]);
      myy = 0.25*(sqmodulus[k+2*sigsize]+sqmodulus[k-2*sigsize]-2 * sqmodulus[k]);
      mxy = 0.25*(sqmodulus[k+1+sigsize]+sqmodulus[k-1-sigsize]
		  -sqmodulus[k+1-sigsize]-sqmodulus[k-1+sigsize]);
      myx = mxy;
      
      /* first four containing the coordinate of left low and right up corner
         --------------------------------------------------------------------*/

      tst[8*bnumber] = (double)(left+1);
      tst[8*bnumber+1] = (double)(down+1);
      tst[8*bnumber+2] = (double)(right+1);
      tst[8*bnumber+3] = (double)(up+1);
      
      /* negative Hessian as Gaussian 
         ----------------------------*/

      tst[8*bnumber+4] = -mxx;
      tst[8*bnumber+5] = -mxy;
      tst[8*bnumber+6] = -myx;
      tst[8*bnumber+7] = -myy;

      bnumber++;
    }
  }

  *pnbblock = bnumber;
}



       



    

  


    
    
  
  

  


       



    

  


    
    
  
  

  




    

  


    
    
  
  

  

