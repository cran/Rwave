/****************************************************************
*               (c) Copyright  1997                             *
*                          by                                   *
*      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                  Princeton University                         *
*                  All right reserved                           *
****************************************************************/

#include "Swave.h"
#include "denoise.h"



/***********************************************************
*  Function: gabor_frequency:
*  ---------
*          Generate Gabor function and derivative in 
*            frequency domain:
*
*   sigma: scale of the Gabor function
*   w: modulated Gabor function (real part)
*   isize: signal size
*
************************************************************/

void gabor_frequency(float sigma,float frequency,double *w,int isize)
{
  double tmp, tmp1;
  int i;
  double twopi;

  twopi = 6.28318530717959;
  
  for(i = 0; i < isize; i++) {
    tmp = (double)(sigma * ((i-frequency*(double)isize/2.) * twopi/isize));
    tmp = -(tmp * tmp)/2;
    w[i] = exp(tmp);
  }
  return;
}



/***********************************************************
*  Function: Sgabor:
*  ---------
*      Continuous Gabor transform.
*
*   input: input signal
*   Ri1, Ii1: Fourier transform of input signal (re. and imag. parts)
*   Ri2: Real part of Fourier transform of Gabor function
*   Oreal: real part of WFT
*   Oimage: Imaginary part of WFT
*   pinputsize: signal size
*   pnbfreq: Number of values for the frequency
*   pfreqstep: frequency step
*   pscale: scale of Gabor function
*
***********************************************************/

void Sgabor(float *input,double *Oreal,double *Oimage,int *pnbfreq,
       float *pfreqstep,int *pinputsize,float *pscale)
{
/*  void multiply();
  void FFT();*/
  int nbfreq, i, j, k, inputsize;
  float scale, freqstep, frequency;
  double *Ri2, *Ri1, *Ii1, *Ii2, *Ii, *Ri;


  scale = *pscale;
  nbfreq = *pnbfreq;
  freqstep = *pfreqstep;
  inputsize = *pinputsize;



  if(!(Ri1 = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ri1 in gabor.c \n");
  if(!(Ii1 = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ii1 in gabor.c \n");

  if(!(Ii2 = (double *)calloc(inputsize,sizeof(double))))
    error("Memory allocation failed for Ri2 in gabor.c \n");
  if(!(Ri2 = (double *)calloc(inputsize,sizeof(double))))
    error("Memory allocation failed for Ri2 in gabor.c \n");


  if(!(Ri = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ri in gabor.c \n");
  if(!(Ii = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ii in gabor.c \n");

  for(i = 0; i < inputsize; i++){
    *Ri = (double)(*input);
    Ri++; input++;
  }
  Ri -= inputsize;
  input -= inputsize;
  
  /* Compute fft of the signal */
  double_fft(Ri1,Ii1,Ri,Ii,inputsize,-1);   
  
  /* Multiply signal and wavelets in the Fourier space */
  frequency = 0;
  for(i = 1; i <= nbfreq; i++) {
    frequency += freqstep;
    gabor_frequency(scale,frequency,Ri2,inputsize); 
    multiply(Ri1,Ii1,Ri2,Ii2,Oreal,Oimage,inputsize);
    double_fft(Oreal,Oimage,Oreal,Oimage,inputsize,1); 
    Oreal += inputsize;
    Oimage += inputsize;  
  }

  Oreal -= inputsize*nbfreq;
  Oimage -= inputsize*nbfreq;

  free((char *)Ri2);
  free((char *)Ri1);
  free((char *)Ii1);
  free((char *)Ii2);
  free((char *)Ri);
  free((char *)Ii);


}




/***********************************************************
*  Function: Svgabor:
*  ---------
*      Continuous Gabor transform for one frequency
*
*   input: input signal
*   Ri1, Ii1: Fourier transform of input signal (re. and imag. parts)
*   Ri2: Real part of Fourier transform of Gabor function
*   Oreal: real part of WFT
*   Oimage: Imaginary part of WFT
*   pinputsize: signal size
*   pfreq: Value of the frequency
*   pscale: scale of Gabor function
*
************************************************************/

void Svgabor(float *input,double *Oreal,double *Oimage,float *pfreq,
	int *pinputsize,float *pscale)
{
/*  void multiply();
  void FFT();*/
  int i, j, k, inputsize;
  float scale, freqstep, frequency;
  double *Ri2, *Ri1, *Ii1, *Ii2, *Ii, *Ri;


  scale = *pscale;
  frequency = *pfreq;
  inputsize = *pinputsize;



  if(!(Ri1 = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ri1 in gabor.c \n");
  if(!(Ii1 = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ii1 in gabor.c \n");

  if(!(Ii2 = (double *)calloc(inputsize,sizeof(double))))
    error("Memory allocation failed for Ri2 in gabor.c \n");
  if(!(Ri2 = (double *)calloc(inputsize,sizeof(double))))
    error("Memory allocation failed for Ri2 in gabor.c \n");


  if(!(Ri = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ri in gabor.c \n");
  if(!(Ii = (double *)calloc(inputsize, sizeof(double))))
    error("Memory allocation failed for Ii in gabor.c \n");

  for(i = 0; i < inputsize; i++){
    *Ri = (double)(*input);
    Ri++; input++;
  }
  Ri -= inputsize;
  input -= inputsize;
  
  /* Compute fft of the signal */
  double_fft(Ri1,Ii1,Ri,Ii,inputsize,-1);   
  
  /* Multiply signal and wavelets in the Fourier space */

  gabor_frequency(scale,frequency,Ri2,inputsize); 
  multiply(Ri1,Ii1,Ri2,Ii2,Oreal,Oimage,inputsize);
  double_fft(Oreal,Oimage,Oreal,Oimage,inputsize,1); 


  free((char *)Ri2);
  free((char *)Ri1);
  free((char *)Ii1);
  free((char *)Ii2);
  free((char *)Ri);
  free((char *)Ii);


}


/***************************************************************
*  Function: gabor_time:
*  ---------
*     Generates a Gabor function in the time domain.
*     The Gabor is centered at the b, and normalized so
*     that psi(0) = 1
*
*   g: Gabor
*   scale: scale of the Gabor
*   isize: window size
*   frequency: frequency of the Gabor
*
* remark: unlike the other functions, this one generates an 
*         array starting at 1, for compatibility with S. 
***************************************************************/

void gabor_time(float *pfrequency,float *pscale, int *pb, 
		 double *g_r, double *g_i,int *pisize)
{
  double tmp, tmp2;
  float frequency = *pfrequency, scale = *pscale;
  int b = *pb, isize = *pisize;
  int i;
  double pi;

  pi = 3.141593; 
  for(i = 1; i <= isize; i++) {
    tmp = (double)((double)(i-b)/scale); 
    tmp2 = exp(-(tmp * tmp)/2.)/scale/sqrt(2.0*pi);
    g_r[i-1] = tmp2*cos(((double)(i-b)) * pi * frequency);
    g_i[i-1] = tmp2*sin(((double)(i-b)) * pi * frequency);
  }
  return;
}


/***************************************************************
*  Function: vgabor_time:
*  ---------
*     Generates many Gabor functions in the time domain.
*     The Gabors are centered on the node of the ridge, and normalized so
*     that psi(0) = 1
*
*   frequency: frequencies along the ridge
*   scale: scale of the Gabor
*   b: locations along the ridge
*   g: Gabor
*   isize: window size
*   nbnode: number of samples at the ridge
*
* remark: unlike the other functions, this one generates an 
*         array starting at 1, for compatibility with S. 
***************************************************************/

void vgabor_time(float *frequency,float *pscale, int *b, 
		 double *g_r, double *g_i,int *pisize, int *pnbnode)
{
  double tmp, tmp2;
  float scale = *pscale;
  int isize = *pisize;
  int i, j, nbnode = *pnbnode;
  double pi;

  int position;
  pi = 3.141593; 

  for(j = 0; j < nbnode; j++) {
         position = b[j];
     for(i = 1; i <= isize; i++) {
         tmp = (double)((double)(i-position)/scale); 
         tmp2 = exp(-(tmp * tmp)/2.)/scale/sqrt(2.0*pi);
         g_r[j * isize + i-1] = tmp2*cos(((double)(i-position)) * pi * frequency[j]);
         g_i[j * isize + i-1] = tmp2*sin(((double)(i-position)) * pi * frequency[j]);
      }
   }
  return;
}




