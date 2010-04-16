#include <stdlib.h>


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
*  Function: multi:
*  ---------
*     Multiplication of 2 vectors in the Fourier domain; the
*     first one is complex-valued, as well as the output; the
*     second one is real-valued.
*
*   Ri1, Ii1: first input vector.
*   Ri2: second input vector (real).
*   Or, Oi: output vector.
*   isize: length of the vectors.
***************************************************************/

void multi(double *Ri1, double *Ii1, double *Ri2, double *Or,
  double *Oi, int isize)
{
  int i;

  for(i = 0; i < isize; i++) {
    Or[i] = Ri1[i] * Ri2[i];
    Oi[i] = Ii1[i] * Ri2[i];
  }
  return;
}


/***************************************************************
*  Function: morlet_frequency:
*  ---------
*     Generates a Morlet wavelet in the frequency domain.
*     The wavelet is centered at the origin, and normalized so
*     that psi(0) =1/sqrt(2 pi)
*
*   w: wavelet
*   scale: scale of the wavelet
*   isize: window size
*   cf: central frequency of the wavelet
***************************************************************/

void morlet_frequency(double cf,double scale,double *w,int isize)
{
  double tmp, tmp1=0;
  int i;
  double twopi;

  twopi = 6.28318530717959;
  
/*  tmp1 = exp(-(cf  * cf)/2); */
  for(i = 0; i < isize; i++) {
    tmp = (double)(scale * i * twopi/isize - cf);
    tmp = -(tmp * tmp)/2;
/*    w[i] = exp(tmp) - tmp1; */
    w[i] = exp(tmp);
  }
  return;
}

/***************************************************************
*  Function: morlet_time:
*  ---------
*     Generates a Morlet wavelet in the time domain.
*     The wavelet is centered at the origin, and normalized so
*     that psi(0) = 1
*
*   w: wavelet
*   scale: scale of the wavelet
*   isize: window size
*   cf: central frequency of the wavelet
*
* remark: unlike the other functions, this one generates an 
*         array starting at 1, for compatibility with S. 
***************************************************************/

void morlet_time(double *pcf,double *pscale, int *pb, 
		 double *w_r, double *w_i,int *pisize)
{
  double tmp, tmp2;
  double cf = *pcf, scale = *pscale;
  int b = *pb, isize = *pisize;
  int i;

  for(i = 1; i <= isize; i++) {
    tmp = (double)((double)(i-b)/scale); 
    tmp2 = exp(-(tmp * tmp)/2.);
    w_r[i-1] = tmp2*cos(tmp*cf)/scale;
    w_i[i-1] = tmp2*sin(tmp*cf)/scale;
  }
  return;
}



/***************************************************************
*  Function: vmorlet_time:
*  ---------
*     Generates Morlet wavelets located on nodes of the ridge,
*     in the time domain.
*     The mother wavelet is centered at the origin, and normalized so
*     that psi(0) = 1/sqrt(2 pi)
*
*   w_r, w_i: wavelet (real and imaginary parts)
*   scale: scale of the wavelets on the ridge samples
*   b: position of the wavelets on the ridge samples
*   isize: window size
*   cf: central frequency of the wavelet
*   nbnodes: number of ridge samples
*
* remark: unlike the other functions, this one generates an 
*         array starting at 1, for compatibility with S. 
***************************************************************/
void vmorlet_time(double *pcf,double *pscale, int *b, 
		 double *w_r, double *w_i,int *pisize, int *pnbnode)
{
  double tmp, tmp2;
  double cf = (double)(*pcf);
  int isize = *pisize;
  int i, j, nbnode = *pnbnode;
  int position;
  double sqtwopi;

  sqtwopi = sqrt(6.28318530717959);

  for(j = 0; j < nbnode; j++) {
         position = b[j];
     for(i = 1; i <= isize; i++) {
         tmp = (double)((double)(i-position)/pscale[j]); 
         tmp2 = exp(-(tmp * tmp)/2.)/pscale[j]/sqtwopi;
         w_r[j * isize + i-1] = tmp2*cos(tmp*cf);
         w_i[j * isize + i-1] = tmp2*sin(tmp*cf);
      }
   }
  return;
}


/*****************************************************************
*  function:  Scwt_morlet_r
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

void Scwt_morlet_r(double *input, double *Oreal, double *Oimage,
   int *pnboctave, int *pnbvoice, int *pinputsize, double *pcenterfrequency)
{	
  int nboctave, nbvoice, i, j, inputsize;
  double centerfrequency, a;
  double *Ri2, *Ri1, *Ii1, *Ii, *Ri;


  centerfrequency = *pcenterfrequency;
  nboctave = *pnboctave;
  nbvoice = *pnbvoice;
  inputsize = *pinputsize;
  if(!(Ri2 = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ri2 in cwt_morlet.c \n");
  if(!(Ri1 = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ri1 in cwt_morlet.c \n");
  if(!(Ii1 = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ii1 in cwt_morlet.c \n");
  if(!(Ri = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ri in cwt_morlet.c \n");
  if(!(Ii = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ii in cwt_morlet.c \n");

  for(i = 0; i < inputsize; i++) {
    Ri[i] = (double)input[i]; 
    Ii[i] = 0.0;
  }

  double_fft(Ri1,Ii1,Ri,Ii,inputsize,-1);   
  
  for(i = 1; i <= nboctave; i++) {
    for(j=0; j < nbvoice; j++) {
      a = (double)(pow((double)2,(double)(i+j/((double)nbvoice))));
      morlet_frequency(centerfrequency,a,Ri2,inputsize); 
      multi(Ri1,Ii1,Ri2,Oreal,Oimage,inputsize);
      double_fft(Oreal,Oimage,Oreal,Oimage,inputsize,1); 
      Oreal = Oreal + inputsize;
      Oimage = Oimage + inputsize;  
    }
  }

}




/*****************************************************************
*  function:  Scwt_morlet
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
*   pcenterfrequency: centralfrequency of Morlet wavelet
******************************************************************/

void Scwt_morlet(double *Rinput,double *Iinput,double *Oreal,
   double *Oimage,int *pnboctave,int *pnbvoice,
   int *pinputsize,double *pcenterfrequency)
{	
  int nboctave, nbvoice, i, j, k, inputsize;
  double centerfrequency, a;
  double *Ri2, *Ri1, *Ii1, *Ii, *Ri;


  centerfrequency = *pcenterfrequency;
  nboctave = *pnboctave;
  nbvoice = *pnbvoice;
  inputsize = *pinputsize;
  if(!(Ri2 = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ri2 in cwt_morlet.c \n");
  if(!(Ri1 = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ri1 in cwt_morlet.c \n");
  if(!(Ii1 = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ii1 in cwt_morlet.c \n");
  if(!(Ri = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ri in cwt_morlet.c \n");
  if(!(Ii = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ii in cwt_morlet.c \n");

  for(i = 0; i < inputsize; i++) {
    Ri[i] = (double)Rinput[i]; 
    Ii[i] = (double)Iinput[i];
  }

  double_fft(Ri1,Ii1,Ri,Ii,inputsize,-1);   
  
  for(i = 1; i <= nboctave; i++) {
    for(j=0; j < nbvoice; j++) {
      a = (double)(pow((double)2,(double)(i+j/((double)nbvoice))));
      morlet_frequency(centerfrequency,a,Ri2,inputsize); 
      multi(Ri1,Ii1,Ri2,Oreal,Oimage,inputsize);
      double_fft(Oreal,Oimage,Oreal,Oimage,inputsize,1); 
      Oreal = Oreal + inputsize;
      Oimage = Oimage + inputsize;  
    }
  }

}



/***************************************************************
*  Function: Svwt_morlet:
*  ---------
*     Computes the continuous wavelet transform of input
*     signal with Morlet wavelet, a fixed scale a
*
*     Rinput, Iinput: real and imaginary parts of the signal
*     Oreal, Oimage: real and imaginary parts of the cwt.
***************************************************************/

void Svwt_morlet(double *Rinput,double *Iinput,double *Oreal,
   double *Oimage,double *pa,int *pinputsize,
   double *pcenterfrequency)
{	
  int octave, voice, nbvoice, i, j, k, inputsize;
  double centerfrequency, a;
  double *Ri2, *Ri1, *Ii1, *Ii, *Ri;


  centerfrequency = *pcenterfrequency;
  a = *pa;
  inputsize = *pinputsize;
  if(!(Ri2 = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ri2 in cwt_morlet.c \n");
  if(!(Ri1 = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ri1 in cwt_morlet.c \n");
  if(!(Ii1 = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ii1 in cwt_morlet.c \n");
  if(!(Ri = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ri in cwt_morlet.c \n");
  if(!(Ii = (double *) R_alloc(inputsize, sizeof(double) )))
    error("Memory allocation failed for Ii in cwt_morlet.c \n");

  for(i = 0; i < inputsize; i++) {
    Ri[i] = (double)Rinput[i]; 
    Ii[i] = (double)Iinput[i];
  }

  double_fft(Ri1,Ii1,Ri,Ii,inputsize,-1);   
  
  morlet_frequency(centerfrequency,a,Ri2,inputsize); 
  multi(Ri1,Ii1,Ri2,Oreal,Oimage,inputsize);
  double_fft(Oreal,Oimage,Oreal,Oimage,inputsize,1); 
  
}


