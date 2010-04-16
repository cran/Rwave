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



/***************************************************************
*  Function: smoothwt
*  ---------
*     Compute a smoothed and subsampled (in the b direction)
*      version of the wavelet transform modulus.
*      Here: rectangular window of size windowsize
*
*   wt: original wavelet transform modulus.
*   swt: smoothed wavelet transform modulus.
***************************************************************/

void smoothwt(double *wt, double *swt, int sigsize, int nbscale, int windowlength)
{
  int a,b,k,adr, cnt=0;
  double normal;

  normal = (double)(2*windowlength -1);
  for(a = 0; a < nbscale; a++){
    for(b = 0; b < sigsize; b += windowlength){
      for(k = 1-windowlength; k< windowlength; k++){
	adr = (b-k+sigsize)%sigsize;
	adr += a*sigsize;
	*swt += wt[adr];
      }
      *swt /= normal;
      swt++; cnt++;
    }
  }
  printf("smoothing done\n");
  /* printf("%d coefficients computed\n",cnt);*/
  return;
}



/***************************************************************
*  Function: smoothwt1
*  ---------
*     Compute a smoothed (but not subsampled)
*       version of the wavelet transform modulus.
*       Here: rectangular window of size windowsize
*
*  wt: original wavelet transform modulus.
*  swt: smoothed wavelet transform modulus.
***************************************************************/

void smoothwt1(double *wt, double *swt, int sigsize, int nbscale, int windowlength)
{
  int a,b,k,adr,cnt=0;
  double normal;

  normal = (double)(2*windowlength -1);
  for(a = 0; a < nbscale; a++){
    for(b = 0; b < sigsize; b ++){
      for(k = 1-windowlength; k< windowlength; k++){
	adr = (b-k+sigsize)%sigsize;
	adr += a*sigsize;
	*swt += wt[adr];
      }
      *swt /= normal;
      swt++;cnt++;
    }
  }
  printf("smoothing done\n");
  printf("%d coefficients computed\n",cnt);
  return;
}



/***************************************************************
*  Function: smoothwt2
*  ---------
*     Compute a smoothed and subsampled (in the b direction)
*      version of the wavelet transform modulus.
*      Here: rectangular window of size windowsize
*
*   wt: original wavelet transform modulus.
*   swt: smoothed wavelet transform modulus.
*   smodsize: actual length of the smoothed transform.
***************************************************************/

void smoothwt2(double *wt, double *swt, int sigsize, int nbscale,
  int windowlength, int *smodsize)
{
  int a,b,k,adr,kmin,kmax;
  double normal;
  int cnt = 0;
  printf("smodsize %d \n",*smodsize);
  printf("number of scales %d \n",nbscale);
  printf("windowlength %d \n",windowlength);
  for(a = 0; a < nbscale; a++){
    for(b = 0; b < sigsize; b += windowlength){
      kmin = max(0,1-windowlength+b);
      kmax = min(sigsize-1,windowlength+b);
      normal = (double)(kmax-kmin+1);
      for(k = kmin; k<= kmax; k++){
/*	adr = (b-k+sigsize)%sigsize;
	adr += a*sigsize; */
	adr = k + a*sigsize; 
	*swt += wt[adr];
      }
      *swt /= normal;
      swt++; cnt++;
    }
  }
  if((cnt%nbscale) != 0){
    printf("Error in smoothwt2\n");
    exit(0);
  }
  *smodsize = cnt/nbscale;
  printf("smoothing done\n");
  printf("%d coefficients computed\n",cnt);
  return;
}


/****************************************************************
*  Function: Smodulus_smoothing:
*  ----------------------------
*  Smoothing and subsampling the modulus of transform
*
*    modulus: modulus of the wavelet transform
*    smodulus: modulus smoothed by a box function
*    sigsize: signal size
*    modsize: signal size after subsampling
*    nscale: total number of scales for CWT
*    subrate: subsampling rate for ridge extraction
*
****************************************************************/


void Smodulus_smoothing(double *modulus, double *smodulus, 
  int *psigsize, int *psmodsize, int *pnscale, int *psubrate)
{
  int sigsize,sub;
  int smodsize;
  int nscale;

  /* Generalities; initializations
     -----------------------------*/
  nscale = *pnscale;
  sigsize = *psigsize;
  smodsize = *psmodsize;
  sub = *psubrate;  

  /* Smooth and subsample the wavelet transform modulus
     --------------------------------------------------*/
  smoothwt2(modulus,smodulus,sigsize,nscale,sub, &smodsize);
  *psmodsize = smodsize;

  return;
}


/***************************************************************
*  Function: Ssmoothwt
*  ---------
*     Interface between the 2 previous functions and S code.
***************************************************************/

void Ssmoothwt(double *smodulus,double * modulus, int *psigsize, int *pnscale, int *psubrate, int *pflag)
{
  int sigsize, nscale, subrate;
  int flag;

  flag = *pflag;

  sigsize = *psigsize;
  nscale = *pnscale;
  subrate = *psubrate;
  if(flag)
    smoothwt1(modulus,smodulus,sigsize,nscale,subrate);
  else
    smoothwt(modulus,smodulus,sigsize,nscale,subrate);
  return;
}


