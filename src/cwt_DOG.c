
/***************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
***************************************************************/

#include "Swave.h"
#include "denoise.h"



/***************************************************************
*  Function: DOG_frequency:
*  ---------
*     Generates a DOG wavelet in the frequency domain.
*     The wavelet is centered at the origin, and normalized so
*     that psi(0) =1/sqrt(2 pi)
*
*   w: wavelet
*   scale: scale of the wavelet
*   isize: window size
*   M: number of vanishing moments
***************************************************************/

void DOG_frequency(int M,float scale,double *w,int isize)
{
  double tmp, cst;
  int i;

  cst = exp(-(double)M*(1.-log((double)M)))/2.;
  for(i = 0; i < isize; i++) {
    tmp = (double)(scale * (double)i * sqrt((double)M)/(double)isize);
    *w = exp(-tmp*tmp/2.)*pow(tmp,(double)M)/cst;
    w++;
  }
  return;
}


/*****************************************************************
*  function:  Scwt_DOG_r
*    Continuous wavelet transform :
*
*   input: (real-valued) input signal
*   Ri1, Ii1: Fourier transform of input signal (real and
*      imaginary parts).
*   Ri2: Real part of Fourier transform of Morlet wavelet
*   Oreal,Oimage: real and imaginary parts of CWT
*   pinputsize: signal size
*   pnboctave: number of scales (powers of 2)
*   pnvoice: number of scales between 2 consecutive powers of 2
*   pcenterfrequency: centralfrequency of Morlet wavelet
******************************************************************/

void Scwt_DOG_r(float *input, double *Oreal, double *Oimage,
   int *pnboctave, int *pnbvoice, int *pinputsize, int *pM)
{	
  int nboctave, nbvoice, i, j, inputsize, M;
  float a;
  double *Ri2, *Ri1, *Ii1, *Ii, *Ri;


  M= *pM;
  nboctave = *pnboctave;
  nbvoice = *pnbvoice;
  inputsize = *pinputsize;
  if(!(Ri2 = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ri2 in cwt_DOG.c \n");
  if(!(Ri1 = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ri1 in cwt_DOG.c \n");
  if(!(Ii1 = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ii1 in cwt_DOG.c \n");
  if(!(Ri = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ri in cwt_DOG.c \n");
  if(!(Ii = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ii in cwt_DOG.c \n");

  for(i = 0; i < inputsize; i++) {
    Ri[i] = (double)input[i]; 
    Ii[i] = 0.0;
  }

  double_fft(Ri1,Ii1,Ri,Ii,inputsize,-1);   
  
  for(i = 1; i <= nboctave; i++) {
    for(j=0; j < nbvoice; j++) {
      a = (float)(pow((double)2,(double)(i+j/((double)nbvoice))));
      DOG_frequency(M,a,Ri2,inputsize); 
      multi(Ri1,Ii1,Ri2,Oreal,Oimage,inputsize);
      double_fft(Oreal,Oimage,Oreal,Oimage,inputsize,1); 
      Oreal = Oreal + inputsize;
      Oimage = Oimage + inputsize;  
    }
  }

  free((char *)Ri2);
  free((char *)Ri1);
  free((char *)Ii1);
  free((char *)Ri);
  free((char *)Ii);
  return;
}




/*****************************************************************
*  function:  Scwt_DOG
*    Continuous wavelet transform :
*
*   input: (a priori complex-valued) input signal
*   Ri1, Ii1: Fourier transform of input signal (real and
*      imaginary parts).
*   Ri2: Real part of Fourier transform of Morlet wavelet
*   Oreal,Oimage: real and imaginary parts of CWT
*   pinputsize: signal size
*   pnboctave: number of scales (powers of 2)
*   pnvoice: number of scales between 2 consecutive powers of 2
*   pM: number of vanishing moments
******************************************************************/

void Scwt_DOG(float *Rinput,float *Iinput,double *Oreal,
   double *Oimage,int *pnboctave,int *pnbvoice,
   int *pinputsize,int *pM)
{	
  int nboctave, nbvoice, i, j, inputsize, M;
  float a;
  double *Ri2, *Ri1, *Ii1, *Ii, *Ri;


  M = *pM;
  nboctave = *pnboctave;
  nbvoice = *pnbvoice;
  inputsize = *pinputsize;
  if(!(Ri2 = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ri2 in cwt_DOG.c \n");
  if(!(Ri1 = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ri1 in cwt_DOG.c \n");
  if(!(Ii1 = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ii1 in cwt_DOG.c \n");
  if(!(Ri = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ri in cwt_DOG.c \n");
  if(!(Ii = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ii in cwt_DOG.c \n");

  for(i = 0; i < inputsize; i++) {
    Ri[i] = (double)Rinput[i]; 
    Ii[i] = (double)Iinput[i];
  }

  double_fft(Ri1,Ii1,Ri,Ii,inputsize,-1);   
  
  for(i = 1; i <= nboctave; i++) {
    for(j=0; j < nbvoice; j++) {
      a = (float)(pow((double)2,(double)(i+j/((double)nbvoice))));
      DOG_frequency(M,a,Ri2,inputsize); 
      multi(Ri1,Ii1,Ri2,Oreal,Oimage,inputsize);
      double_fft(Oreal,Oimage,Oreal,Oimage,inputsize,1); 
      Oreal = Oreal + inputsize;
      Oimage = Oimage + inputsize;  
    }
  }

  free((char *)Ri2);
  free((char *)Ri1);
  free((char *)Ii1);
  free((char *)Ri);
  free((char *)Ii);
}



/***************************************************************
*  Function: Svwt_DOG:
*  ---------
*     Computes the continuous wavelet transform of input
*     signal with Morlet wavelet, a fixed scale a
*
*     Rinput, Iinput: real and imaginary parts of the signal
*     Oreal, Oimage: real and imaginary parts of the cwt.
***************************************************************/

void Svwt_DOG(float *Rinput,float *Iinput,double *Oreal,
   double *Oimage,float *pa,int *pinputsize,
   int *pM)
{	
  int i, inputsize, M;
  float a;
  double *Ri2, *Ri1, *Ii1, *Ii, *Ri;


  M = *pM;
  a = *pa;
  inputsize = *pinputsize;
  if(!(Ri2 = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ri2 in cwt_DOG.c \n");
  if(!(Ri1 = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ri1 in cwt_DOG.c \n");
  if(!(Ii1 = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ii1 in cwt_DOG.c \n");
  if(!(Ri = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ri in cwt_DOG.c \n");
  if(!(Ii = (double *)malloc(sizeof(double) * inputsize)))
    error("Memory allocation failed for Ii in cwt_DOG.c \n");

  for(i = 0; i < inputsize; i++) {
    Ri[i] = (double)Rinput[i]; 
    Ii[i] = (double)Iinput[i];
  }

  double_fft(Ri1,Ii1,Ri,Ii,inputsize,-1);   
  
  DOG_frequency(M,a,Ri2,inputsize); 
  multi(Ri1,Ii1,Ri2,Oreal,Oimage,inputsize);
  double_fft(Oreal,Oimage,Oreal,Oimage,inputsize,1); 
  
  free((char *)Ri2);
  free((char *)Ri1);
  free((char *)Ii1);
  free((char *)Ri);
  free((char *)Ii);
}












