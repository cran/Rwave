/****************************************************************
*               (c) Copyright  1997                             *
*                          by                                   *
*      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                  Princeton University                         *
*                  All right reserved                           *
****************************************************************/

#include "Swave.h"
#include "denoise.h"



/*****************************************************************
*  Function: morlet_frequencyph:
*    Generation of Morlet wavelet and derivative (frequency domain)
*
*   cf: central frequency
*   w: scaled wavelet (real part)
*   wd: scaled derivative of wavelet (imaginary part,up to a
*       - sign and scale factor)
*   isize: signal size
******************************************************************/

void morlet_frequencyph(float cf,float scale,double *w,
  double *wd,int isize)
{
  double tmp, tmp1, tmp0;
  int i;
  double twopi;

  twopi = 6.28318530717959;
  
  tmp1 = exp(-(cf  * cf)/2);
  for(i = 0; i < isize; i++) {
    tmp0 = (double)(scale * i * twopi/isize);
    tmp = (double)(tmp0 - cf);
    tmp = -(tmp * tmp)/2;
    w[i] = exp(tmp) - tmp1;
    wd[i]= w[i]*tmp0/(double)scale;
  }
  return;
}

/*****************************************************************
*  function normalization
*   Normalize the derivative of CWT by the square-modulus of CWT
*
*  Oreal, Oimage: real and imaginary parts of wavelet transform.
*  Odreal, Odimage: real and imaginary parts of wavelet transform
*      derivative.
******************************************************************/

void normalization(double *Oreal, double *Oimage, double *Odreal,
  double *Odimage, int cwtlength)
{
  double tmp;
  int i;
  
  for(i=0;i<cwtlength;i++){
    tmp = (*Oreal)*(*Oreal) + (*Oimage)*(*Oimage);
    *Odreal /=tmp; *Odimage /= tmp;
    Oreal++; Oimage++;
    Odreal ++; Odimage++;
  }
  return;
}


/*****************************************************************
*  function: f_function
*   compute the b-derivative of the phase of the CWT at 
*   scale a, and substract the "natural frequency" of the
*   wavelet at the same scale.
*
*  Oreal, Oimage: real and imaginary parts of wavelet transform.
*  Odreal, Odimage: real and imaginary parts of wavelet transform
*      derivative.
*  f:  f function
*  cf: central frequency of thye wavelet
*  inputsize,nbvoice,nboctave: parameters of wavelet transform.
******************************************************************/

void f_function(double *Oreal, double *Oimage, double *Odreal,
  double *Odimage, double *f, float cf,int inputsize,int nbvoice,
  int nboctave)
{
  int i, j, k;
  float scale;
  
  for(i = 1; i <= nboctave; i++) {
    for(j=0; j < nbvoice; j++) {
      scale =(float)(pow((double)2,(double)(i+j/((double)nbvoice))));
      for(k=0;k<inputsize;k++){
	*f = -(*Odreal)*(*Oimage)+ (*Odimage)*(*Oreal);
	*f -= cf/scale;
	Oreal++; Oimage++;
	Odreal ++; Odimage++;
	f++;
      }
    }
  }
  return;
}


/*****************************************************************
*  function: w_reassign
*   compute the b-derivative of the phase of the CWT at 
*   scale a, and use it to reassign vertically the energy in the
*   time-frequency plane.
*
*  Oreal, Oimage: real and imaginary parts of wavelet transform.
*  Odreal, Odimage: real and imaginary parts of wavelet transform
*      derivative.
*  squeezed:  squeezed (reassigned) time-frequency representation
*  cf: central frequency of the wavelet
*  inputsize,nbvoice,nboctave: parameters of wavelet transform.
******************************************************************/

void w_reassign(double *Oreal, double *Oimage, double *Odreal,
  double *Odimage, double *squeezed_r, double *squeezed_i, float cf,
  int inputsize,int nbvoice,int nboctave)
{
  int i, j, k, scale2;
  float scale, tmp;
  
  for(i = 1; i <= nboctave; i++) {
    for(j=0; j < nbvoice; j++) {
      scale =(float)(pow((double)2,(double)(i+j/((double)nbvoice))));
      for(k=0;k<inputsize;k++){
	tmp = (float)(-(*Odreal)*(*Oimage)+ (*Odimage)*(*Oreal));
	scale2 =  (int)(nbvoice*(log(cf/tmp/2.)/log(2.))+.5);
	if (inrange(0,scale2,nbvoice*nboctave-1)){
	  squeezed_r[scale2*inputsize + k] += (*Oreal);
	  squeezed_i[scale2*inputsize + k] += (*Oimage);
	}
	Oreal++; Oimage++;
	Odreal ++; Odimage++;
      }
    }
  }
  return;
}


/*****************************************************************
*  function:  Scwt_phase
*    Continuous wavelet transform (and derivative):
*
*   input: input signal
*   Ri1, Ii1: Fourier transform of input signal (real and
*      imaginary parts).
*   Ri2: Real part of Fourier transform of Morlet wavelet
*   Idi2: Imaginary part of Fourier transform of Derivative
*      of Morlet wavelet
*   Oreal,Oimage: real and imaginary parts of CWT
*   Odreal,Odimage: real and imaginary parts of CWT derivative
*   pinputsize: signal size
*   pnboctave: number of scales (powers of 2)
*   pnvoice: number of scales between 2 consecutive powers of 2
*   pcenterfrequency: centralfrequency of Morlet wavelet
******************************************************************/

void Scwt_phase(float *input, double *Oreal, double *Oimage,
  double *f, int *pnboctave, int *pnbvoice, int *pinputsize,
  float *pcenterfrequency)
{
  int nboctave, nbvoice, i, j, inputsize;
  float centerfrequency, a;
  double *Ri2, *Ri1, *Ii1, *Ii2, *Rdi2, *Idi2, *Ii, *Ri;
  double *Odreal, *Odimage;


  centerfrequency = *pcenterfrequency;
  nboctave = *pnboctave;
  nbvoice = *pnbvoice;
  inputsize = *pinputsize;

  /* Memory allocations
     ------------------*/
  if(!(Odreal = (double *)calloc(inputsize*nbvoice*nboctave, sizeof(double))))
    error("Memory allocation failed for Ri1 in cwt_phase.c \n");
  if(!(Odimage = (double *)calloc(inputsize*nbvoice*nboctave, sizeof(double))))
    error("Memory allocation failed for Ii1 in cwt_phase.c \n");

  if(!(Ri1 = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ri1 in cwt_phase.c \n");
  if(!(Ii1 = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ii1 in cwt_phase.c \n");

  if(!(Ii2 = (double *)calloc(inputsize,sizeof(double))))
    error("Memory allocation failed for Ri2 in cwt_phase.c \n");
  if(!(Ri2 = (double *)calloc(inputsize,sizeof(double))))
    error("Memory allocation failed for Ri2 in cwt_phase.c \n");

  if(!(Idi2 = (double *)calloc(inputsize,sizeof(double))))
    error("Memory allocation failed for Ri2 in cwt_phase.c \n");
  if(!(Rdi2 = (double *)calloc(inputsize,sizeof(double))))
    error("Memory allocation failed for Ri2 in cwt_phase.c \n");

  if(!(Ri = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ri in cwt_phase.c \n");
  if(!(Ii = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ii in cwt_phase.c \n");

  for(i = 0; i < inputsize; i++){
    *Ri = (double)(*input);
    Ri++; input++;
  }
  Ri -= inputsize;
  input -= inputsize;
  
  /* Compute fft of the signal
     -------------------------*/
  double_fft(Ri1,Ii1,Ri,Ii,inputsize,-1);   
  
  /* Multiply signal and wavelets in the Fourier space
     -------------------------------------------------*/
  for(i = 1; i <= nboctave; i++) {
    for(j=0; j < nbvoice; j++) {
      a = (float)(pow((double)2,(double)(i+j/((double)nbvoice))));
      morlet_frequencyph(centerfrequency,a,Ri2,Idi2,inputsize); 
      multiply(Ri1,Ii1,Ri2,Ii2,Oreal,Oimage,inputsize);
      multiply(Ri1,Ii1,Rdi2,Idi2,Odreal,Odimage,inputsize);
      double_fft(Oreal,Oimage,Oreal,Oimage,inputsize,1); 
      double_fft(Odreal,Odimage,Odreal,Odimage,inputsize,1); 
      Oreal += inputsize;
      Oimage += inputsize;  
      Odreal += inputsize;
      Odimage += inputsize; 
    }
  }

  Oreal -= inputsize*nbvoice*nboctave;
  Odreal -= inputsize*nbvoice*nboctave;
  Oimage -= inputsize*nbvoice*nboctave;
  Odimage -= inputsize*nbvoice*nboctave;

  free((char *)Ri2);
  free((char *)Ri1);
  free((char *)Ii1);
  free((char *)Ii2);
  free((char *)Ri);
  free((char *)Ii);

  /* Normalize the cwt and compute the f function
     --------------------------------------------*/
  normalization(Oreal, Oimage, Odreal, Odimage,
    inputsize*nbvoice*nboctave);

  f_function(Oreal, Oimage, Odreal, Odimage, f,
    centerfrequency,inputsize,nbvoice,nboctave);

  return;
}



/*****************************************************************
*  function:  Scwt_squeezed
*    Squeezed continuous wavelet transform:
*
*   input: input signal
*   Ri1, Ii1: Fourier transform of input signal (real and
*      imaginary parts).
*   Ri2: Real part of Fourier transform of Morlet wavelet
*   Idi2: Imaginary part of Fourier transform of Derivative
*      of Morlet wavelet
*   Oreal,Oimage: real and imaginary parts of CWT
*   Odreal,Odimage: real and imaginary parts of CWT derivative
*   pinputsize: signal size
*   pnboctave: number of scales (powers of 2)
*   pnvoice: number of scales between 2 consecutive powers of 2
*   pcenterfrequency: centralfrequency of Morlet wavelet
******************************************************************/

void Scwt_squeezed(float *input, double *squeezed_r,
  double *squeezed_i, int *pnboctave, int *pnbvoice,
  int *pinputsize, float *pcenterfrequency)
{
  int nboctave, nbvoice, i, j, inputsize, bigsize;
  float centerfrequency, a;
  double *Ri2, *Ri1, *Ii1, *Ii2, *Rdi2, *Idi2, *Ii, *Ri;
  double *Oreal, *Oimage, *Odreal, *Odimage;


  centerfrequency = *pcenterfrequency;
  nboctave = *pnboctave;
  nbvoice = *pnbvoice;
  inputsize = *pinputsize;
  bigsize = inputsize*nbvoice*nboctave;

  /* Memory allocations
     ------------------*/
  if(!(Oreal = (double *)calloc(bigsize, sizeof(double))))
    error("Memory allocation failed for Ri1 in cwt_phase.c \n");
  if(!(Oimage = (double *)calloc(bigsize, sizeof(double))))
    error("Memory allocation failed for Ii1 in cwt_phase.c \n");

  if(!(Odreal = (double *)calloc(bigsize, sizeof(double))))
    error("Memory allocation failed for Ri1 in cwt_phase.c \n");
  if(!(Odimage = (double *)calloc(bigsize, sizeof(double))))
    error("Memory allocation failed for Ii1 in cwt_phase.c \n");

  if(!(Ri1 = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ri1 in cwt_phase.c \n");
  if(!(Ii1 = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ii1 in cwt_phase.c \n");

  if(!(Ii2 = (double *)calloc(inputsize,sizeof(double))))
    error("Memory allocation failed for Ri2 in cwt_phase.c \n");
  if(!(Ri2 = (double *)calloc(inputsize,sizeof(double))))
    error("Memory allocation failed for Ri2 in cwt_phase.c \n");

  if(!(Idi2 = (double *)calloc(inputsize,sizeof(double))))
    error("Memory allocation failed for Ri2 in cwt_phase.c \n");
  if(!(Rdi2 = (double *)calloc(inputsize,sizeof(double))))
    error("Memory allocation failed for Ri2 in cwt_phase.c \n");

  if(!(Ri = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ri in cwt_phase.c \n");
  if(!(Ii = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ii in cwt_phase.c \n");

  for(i = 0; i < inputsize; i++){
    *Ri = (double)(*input);
    Ri++; input++;
  }
  Ri -= inputsize;
  input -= inputsize;
  
  /* Compute fft of the signal
     -------------------------*/
  double_fft(Ri1,Ii1,Ri,Ii,inputsize,-1);   
  
  /* Multiply signal and wavelets in the Fourier space
     -------------------------------------------------*/
  for(i = 1; i <= nboctave; i++) {
    for(j=0; j < nbvoice; j++) {
      a = (float)(pow((double)2,(double)(i+j/((double)nbvoice))));
      morlet_frequencyph(centerfrequency,a,Ri2,Idi2,inputsize); 
      multiply(Ri1,Ii1,Ri2,Ii2,Oreal,Oimage,inputsize);
      multiply(Ri1,Ii1,Rdi2,Idi2,Odreal,Odimage,inputsize);
      double_fft(Oreal,Oimage,Oreal,Oimage,inputsize,1); 
      double_fft(Odreal,Odimage,Odreal,Odimage,inputsize,1); 
      Oreal += inputsize;
      Oimage += inputsize;  
      Odreal += inputsize;
      Odimage += inputsize; 
    }
  }

  Oreal -= bigsize;
  Odreal -= bigsize;
  Oimage -= bigsize;
  Odimage -= bigsize;

  free((char *)Ri2);
  free((char *)Ri1);
  free((char *)Ii1);
  free((char *)Ii2);
  free((char *)Ri);
  free((char *)Ii);

  /* Normalize the cwt and compute the squeezed transform
     ----------------------------------------------------*/
  normalization(Oreal, Oimage, Odreal, Odimage, bigsize);

  w_reassign(Oreal, Oimage, Odreal, Odimage, squeezed_r,
    squeezed_i, centerfrequency,inputsize, nbvoice, nboctave);

  return;

}



