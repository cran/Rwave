/****************************************************************
*	$ Log: denoise.h,v	$				*
*               (c) Copyright  1994                             *
*                          by                                   *
*      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                  University of California, Irvine             *
*                  All right reserved                           *
****************************************************************/






/**************************************************************************

Type definitions:
	wave
	WT
**************************************************************************/

#include "complex.h"



/**************************************************************************
*  Function declarations:
**************************************************************************/


/* In FFT.c
   --------*/
void double_fft(double *Or,double *Oi,double *Ir,double *Ii,
  int isize,int isign);



/* In four1.c
   ----------*/
void four1(double data[],int nn,int isign);



/* In cwt_maxima.c
   ---------------*/
void Scwt_gmax(double *input, double *output, int *pnrow, int *pncol);

void Scwt_mridge(double *input, double *output, int *pnrow, int *pncol);



/* In cwt_phase.c
   --------------*/
void morlet_frequencyph(float cf,float scale,double *w,
  double *wd,int isize);

void normalization(double *Oreal, double *Oimage, double *Odreal,
  double *Odimage, int cwtlength);

void f_function(double *Oreal, double *Oimage, double *Odreal,
  double *Odimage, double *f, float cf,int inputsize,int nbvoice,
  int nboctave);

void Scwt_phase(float *input, double *Oreal, double *Oimage,
  double *f, int *pnboctave, int *pnbvoice, int *pinputsize,
  float *pcenterfrequency);

void w_reassign(double *Oreal, double *Oimage, double *Odreal,
  double *Odimage, double *squeezed_r, double *squeezed_i, float cf,
  int inputsize,int nbvoice,int nboctave);

void Scwt_squeezed(float *input, double *squeezed_r,
  double *squeezed_i, int *pnboctave, int *pnbvoice,
  int *pinputsize, float *pcenterfrequency);



/* In cwt_morlet.c
   ---------------*/
void multi(double *Ri1, double *Ii1, double *Ri2, double *Or,
   double *Oi, int isize);

void morlet_frequency(float cf,float scale,double *w,int isize);

void morlet_time(float *pcf,float *pscale, int *pb, 
		 fcomplex *w,int *pisize);

void vmorlet_time(float *pcf,float *pscale, int *b, 
		 double *w_r, double *w_i,int *pisize, int *pnbnode);

void Scwt_morlet_r(float *input, double *Oreal, double *Oimage,
   int *pnboctave, int *pnbvoice, int *pinputsize, float *pcenterfrequency);

void Scwt_morlet(float *Rinput,float *Iinput,double *Oreal,
   double *Oimage,int *pnboctave,int *pnbvoice,
   int *pinputsize,float *pcenterfrequency);

void Svwt_morlet(float *Rinput,float *Iinput,double *Oreal,
   double *Oimage,float *pa,int *pinputsize,
   float *pcenterfrequency);



/* In cwt_thierry.c
   ----------------*/
void thierry_frequency(int M,float scale,double *w,int isize);

void Scwt_thierry_r(float *input, double *Oreal, double *Oimage,
   int *pnboctave, int *pnbvoice, int *pinputsize, int *pM);

void Scwt_thierry(float *Rinput,float *Iinput,double *Oreal,
   double *Oimage,int *pnboctave,int *pnbvoice,
   int *pinputsize,int *pM);

void Svwt_thierry(float *Rinput,float *Iinput,double *Oreal,
   double *Oimage,float *pa,int *pinputsize,
   int *pM);



/* In cwt_DOG.c
   ------------*/
void DOG_frequency(int M,float scale,double *w,int isize);

void Scwt_DOG_r(float *input, double *Oreal, double *Oimage,
   int *pnboctave, int *pnbvoice, int *pinputsize, int *pM);

void Scwt_DOG(float *Rinput,float *Iinput,double *Oreal,
   double *Oimage,int *pnboctave,int *pnbvoice,
   int *pinputsize,int *pM);

void Svwt_DOG(float *Rinput,float *Iinput,double *Oreal,
   double *Oimage,float *pa,int *pinputsize,
   int *pM);



/* In icm.c
   --------*/
void Sridge_icm(float *cost, double *smodulus, float *phi,
  float *plambda, float *pmu, int *psigsize, int *pnscale,
  int *piteration,int *pcount, int *psub, int *psmodsize);




/* In gabor.c
   ----------*/
void gabor_frequency(float sigma,float frequency,double *w,int isize);

void Sgabor(float *input, double *Oreal, double *Oimage, int *pnbfreq,
   float *pfreqstep, int *pinputsize, float *pscale);

void Svgabor(float *input,double *Oreal,double *Oimage,float *pfreq,
	int *pinputsize,float *pscale);

void vgabor_time(float *frequency,float *pscale, int *b, 
		 double *g_r, double *g_i,int *pisize, int *pnbnode);



/* In randomwalker.c
   -----------------*/
double ran1(long *idum);

void randomwalker(int sigsize,int *num);

void randomsnaker(int sigsize,int *num);

void randomwalker2(int sigsize,int *num, long *seed);



/* In randomwalker.c
   -----------------*/
double oldran1(long *idum);



/* In smoothwt.c
   -------------*/
void smoothwt(double *wt, double *swt, int sigsize, int nbscale,
   int windowlength);

void smoothwt1(double *wt, double *swt, int sigsize, int nbscale,
   int windowlength);

void smoothwt2(double *wt, double *swt, int sigsize, int nbscale,
   int windowlength, int *smodsize);

void Smodulus_smoothing(double *modulus, double *smodulus, 
  int *psigsize, int *psmodsize, int *pnscale, int *psubrate);

void Ssmoothwt(double *smodulus,double * modulus, int *psigsize,
   int *pnscale, int *psubrate, int *pflag);



/* In splridge.c
   -------------*/
void splridge(int rate, float *y, int n, float *yy);



/* In splsnake.c
   -------------*/
void splsnake(int rate, float *x, float *y, int n, float *yy);    



/* In snakesub.c
   -------------*/
void snakesub(float *rho,int rate,int snakesize);

void snakexpand(float *rho,int rate,int snakesize);



/* In ridge_annealing.c
   --------------------*/
void Sridge_annealing(float *cost, double *smodulus,
  float *phi, float *plambda, float *pmu, float *pc, int *psigsize,
  int *pnscale, int *piteration, int *pstagnant, int *pseed,
  int *pcount, int *psub, int *pblocksize, int *psmodsize);



/* In ridge_coronoid.c
   -------------------*/
void Sridge_coronoid(float *cost, double *smodulus,
  float *phi, float *plambda, float *pmu, float *pc, int *psigsize,
  int *pnscale, int *piteration, int *pstagnant, int *pseed,
  int *pcount, int *psub, int *pblocksize, int *psmodsize);


/* In snake_annealing.c
   --------------------*/
void Ssnake_annealing(float *cost, double *smodulus,
  float *phi, float *rho, float *plambda, float *pmu,
  float *plambda2, float *pmu2, float *pc, int *psigsize,
  int *psnakesize, int *pnscale, int *piteration,
  int *pstagnant, int *pseed, int *pcount, int *psub,
  int *pblocksize, int *psmodsize);


/* In multiply.c
   -------------*/
void multiply(double *Ri1, double *Ii1, double *Ri2, double *Ii2,
   double *Or,double *Oi, int isize);



/* In spline.c
   -----------*/
void spline(double x[], double y[], int n, double yp1, double ypn,
   double y2[]);



/* In splint.c
   -----------*/
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);



/* In splint2.c
   ------------*/
void splint2(double xa[], double ya[], double y2a[], int n, double x,
   double *y, double *yp);



/* In Util_denoising.c
   -------------------*/
int find2power(int n);

void error(char *s);



/* In bee_annealing.c
   ------------------*/
void Sbee_annealing(double *smodulus, double *beemap,
	       float *pc,
	       int *psigsize, int *pnscale, int *piteration,
               int *pseed, int *pbstep, int *pnbbee,
               int *pintegral, int *pchain, int *flag);

/* In ridrep.c
   -----------*/
void ridrep(float *signal,double *transform,float *phi,int length,
   int width, int bmin, int bmax, int amin, int amax,char *filename);

void snakerep(float *signal,double *transform,float *phi,float *rho,
   int nb_nodes,int length,int width, int bmin, int bmax, int amin, 
   int amax,char *filename);

void marsrep(float *signal,double *transform,float *phi, int length,
   int width, int b_start, int b_end, int bmin, int bmax, int amin,
   int amax,char *filename);



/* In delog.c
   ----------*/
void delog(float *phi, float *phi2, int A, int nvoice, int B);

void delog_inv(float *phi, float *phi2, int A, int nvoice, int B);



/* In initialization.c
   -------------------*/
void initialization(int *b_start , int *b_end,int *nitermax,
   float *a_0, int *rate);


/* In ridrecon.c
   -------------
void ridrecon(float *ridge, float *skel,float *signal, wave *W,
   int start, int end, float omega);
   */


/* In variance.c
   -------------*/
float variance(float *signal, int length);


/* In normalize.c
   --------------*/
void normalize(float *signal, float norm, int length);

void dnormalize(double *signal, double norm, int length);


/* In V_pot.c
   ----------*/
void V_pot(double *modulus, float *V, int B, int A);


/* In clean.c
   ----------*/
void fclean(float *array,  int length);

void dclean(double *array,  int length);

void iclean(int *array,  int length);


/* In power_law.c
   --------------*/
void power_law(float *V, float alpha, float cst, int nscale, int nvoice);


/* In ridge_rec.c
   --------------*/
void Sridge_rec(float *chain, float *recsig, double *gabor, int *pnbchain,
   int *psigsize, int *pnfreq);


/* In crazy_family.c
   -----------------*/
void Scrazy_family(double *ridgemap,float *orderedmap,int *chain,
   int *pnbchain,int *psigsize,int *pnscale,int *pbstep,float *pthreshold);

void orderedmap_thresholded(float *orderedmap,int sigsize,int nscale,
   int *chain,int nbchain);

void chain_thresholded(double *mridge,int sigsize,int *chain,int *id,
   int nbchain,float threshold, int bstep);

void reordering(int *chain, int sigsize, int nbchain);


/* In transpose.c
   --------------*/
void transpose(float *inmat,  float *outmat, int length1, int length2);

void itranspose(int *inmat,  int *outmat, int length1, int length2);

void dtranspose(double *inmat,  double *outmat, int length1, int length2);


/* in simul.c
   --------- */
void local_mean(float *mean, float *s, int np );

float gasdev(long *idum);

float variance(float *s, int np );

float denominator(float *Wf, int np );

float numerator(float *Wf, int resoln, int np );

void normal_histo( float ***histo, int max_resoln, int sample_size );

void bootstrap_histo(float ***histo, float *s, int max_resoln,
  int sample_size );

void normal_pval_compute(float *pval, float *s, int *max_resoln_ptr,
  int *np_ptr, int *num_of_windows_ptr, int *window_size_ptr );

void bootstrap_pval_compute(float *pval, float *s, int *max_resoln_ptr,
  int *np_ptr, int *num_of_windows_ptr, int *window_size_ptr );

void nthresh_compute(float *nthresh, float *s, int *maxresoln_ptr,
  int *sample_size_ptr, float prct );

void bthresh_compute(float *bthresh, float *s, int *maxresoln_ptr,
  int *sample_size_ptr, float prct );

float p_value(float T, float **histo, int resoln, int histo_size );

void compute_pval_average(float *pval, float **p, int max_resoln, int np,
  int num_of_windows, int window_size );



/* in qcksrt.c
   ----------- */
void qcksrt(int n,float arr[]);


/* in WV.c
   ------- */
void freq_parity(float frequency,double *win,double *wout,
   int isize,int sign);

void WV_freq_mult(float frequency,double *Ri,double *Ii,
  double *Ro, double *Io,int isize);

void WV(float *input, double *Oreal,double *Oimage,int *pnbfreq,
   float *pfreqstep,int *pinputsize);



/* in optimize.c
   ------------- */
void Lpnorm(double *norm, double *p, double *Rmat, double *Imat,
   int *length, int *width);

void entropy(double *entr, double *Rmat, double *Imat,
   int *length, int *width);



