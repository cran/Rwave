

/***************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
****************************************************************/




/***************************************************************************/
/*                                                                         */
/* Functions for computing p-values for mallat wavelets                    */
/*                                                                         */
/***************************************************************************/


#include "denoise.h"
#include "Swave.h"

#define HISTO_SIZE    500

extern long idum;




/**************************************************************************/
/* local_mean */
/**************************************************************************/
#define LOCAL_LENGTH    17

void local_mean(float *mean, float *s, int np )
{
  float *s_sym;
  int i, j, t;
  float sum;

  if(!(s_sym = (float *) malloc( 2*np*sizeof(float))))
    error("Memory allocation failed in s_sym at simul.c \n");

  for ( t = 0; t < np; t++ )
    s_sym[t] = s_sym[2*np-t-1] = s[t];

  for ( t = 0; t < np; t++ )
  {
    sum = 0.0;
    for ( j = t-8, i = 0; i < LOCAL_LENGTH; i++, j++ )
      sum += (j < 0) ? s_sym[-j-1] : s_sym[j];
    mean[t] = sum / (float) LOCAL_LENGTH;
  }

  free( s_sym );
}
#undef LOCAL_LENGTH



/**************************************************************************/
/* variance */
/**************************************************************************/

float variance(float *s, int np )
{
  float mean, sum;
  float *temp;
  int i;

  if(!(temp = (float *)malloc(np * sizeof(float))))
    error("Memory allocation failed for temp at simul.c ");

  for ( i = 0, mean = 0.0; i < np; i++ )
    mean += s[i];
  mean /= (float) np;

  for ( i = 0; i < np; i++ )
    temp[i] = s[i] - mean;

  for ( i = 0, sum = 0.0; i < np; i++ )
    sum += temp[i] * temp[i];

  free( temp );
  return sum / (float) np;
}

/****************************************************************************/
/* Compute p-value                                                          */
/* Def:  p_value = P{ T > T_obs }                                           */
/****************************************************************************/

float p_value(float T, float **histo, int resoln, int histo_size )
{
  int count = 0;
  int i;

  for ( i = 0; i < histo_size; i++ )
    if ( histo[resoln][i] > T )
    {
      count = histo_size - i;  /* num of values that > T */
      break;
    }
  return (float) count / histo_size;
}

/****************************************************************************/
/* Compute p-value average                                                  */
/****************************************************************************/

void compute_pval_average(float *pval, float **p, int max_resoln, int np,
  int num_of_windows, int window_size )

{
  int interval_length = window_size / 4;
  int num_of_intervals = np / interval_length;
  float *temp;
  int m, i, j, k;

  if(!(temp = (float *)malloc(num_of_intervals * sizeof(float))))
    error("Memory allocation failed for temp at simul.c \n");

  for ( j = 1; j <= max_resoln; j++ )
  {
    /* Compute the average for each interval */
    
    temp[0] = p[j][0];
    temp[1] = (p[j][0] + p[j][1]) / 2;
    temp[2] = (p[j][0] + p[j][1] + p[j][2]) / 3;

    for ( i = 3, k = i; i < num_of_intervals-3; i++, k++ )
      temp[i] = (p[j][k-3] + p[j][k-2] + p[j][k-1] + p[j][k])/ 4;

    temp[num_of_intervals-1] = p[j][num_of_windows-1];
    temp[num_of_intervals-2] = (p[j][num_of_windows-1] +
				p[j][num_of_windows-2]) / 2;
    temp[num_of_intervals-3] = (p[j][num_of_windows-1] +
				p[j][num_of_windows-2] +
				p[j][num_of_windows-3]) / 3;

    for ( i = 0 ; i < num_of_intervals; i++ )
    {
      for ( m = (j-1)*np+i*interval_length, k = 0; k < interval_length; k++, m++ )
	pval[m] = temp[i];
    }
  }
  free( temp );
}



/****************************************************************************/
/* Gaussian random variables generator                                      */
/****************************************************************************/
float gasdev(long *idum)
{
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;

  if  (iset == 0) {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
}




/****************************************************************************/
/* denominator of the statistics = (Wf[1]^2) + (Wf[2]^2)                    */
/****************************************************************************/

float denominator(float *Wf, int np )
{
  float den = 0.0;
  int t;

  for ( t = 0; t < 2*np; t++ )  /* resoln 1 and 2 */
    den += Wf[t] * Wf[t];
  return den;
}

/****************************************************************************/
/* numerator of the statistics = sqrt of (Wf[j]^4)                          */
/****************************************************************************/

float numerator(float *Wf, int resoln, int np )
{
  float sum = 0.0;
  int t, i;

  for ( t = (resoln-1)*np, i = 0; i < np; i++, t++ )
    sum += Wf[t] * Wf[t] * Wf[t] * Wf[t];
  return (float) sqrt( (double) sum );
}

/***************************************************************************/
/* Compute mallat normal histogram for computing p-values */
/***************************************************************************/

void normal_histo( float ***histo, int max_resoln, int sample_size )
{
  float *Sf = (float *) malloc( (max_resoln+1) * sample_size * sizeof(float) );
  float *Wf = (float *) malloc( max_resoln * sample_size * sizeof(float) );
  float *sample = (float *) malloc( sample_size * sizeof(float) );

  float den;
  int b, i, j;


  *histo = (float **) malloc( (max_resoln+1) * sizeof(float *) );
  for ( j = 1; j <= max_resoln; j++ )
    (*histo)[j] = (float *) malloc( HISTO_SIZE * sizeof(float) );

  for ( b = 0; b < HISTO_SIZE; b++ )
  {
    for ( i = 0; i < sample_size; i++ )
      sample[i] = gasdev( &idum );

    Sf_compute( Sf, sample, &max_resoln, &sample_size,"Gaussian1" );
    Wf_compute( Wf, Sf, &max_resoln, &sample_size,"Gaussian1" );
    den = denominator( Wf, sample_size );

    for ( j = 1; j <= max_resoln; j++ )
      (*histo)[j][b] = numerator( Wf, j, sample_size ) / den;
  }

  for ( j = 1; j <= max_resoln; j++ )
    qcksrt( HISTO_SIZE, (*histo)[j]-1 );

  free( Sf );
  free( Wf );
  free( sample );
}

/***************************************************************************/
/* Compute mallat bootstrap histogram for computing p-values */
/***************************************************************************/

void bootstrap_histo(float ***histo, float *s, int max_resoln,
  int sample_size )
{
  float *Sf = (float *) malloc( (max_resoln+1) * sample_size * sizeof(float) );
  float *Wf = (float *) malloc( max_resoln * sample_size * sizeof(float) );

  float *sample = (float *) malloc( sample_size * sizeof(float) );
  float *bsample = (float *) malloc( sample_size * sizeof(float) );
  float *mean = (float *) malloc( sample_size * sizeof(float) );

  float den;
  int b, i, j;
  int k = sample_size - 16;  

  for ( i = 0; i < sample_size; i++ )
    bsample[i] = s[i];
  local_mean( mean, bsample, sample_size );
  for ( i = 0; i < sample_size; i++ )
    bsample[i] -= mean[i];

  *histo = (float **) malloc( (max_resoln+1) * sizeof(float *) );
  for ( j = 1; j <= max_resoln; j++ )
    (*histo)[j] = (float *) malloc( HISTO_SIZE * sizeof(float) );

  for ( b = 0; b < HISTO_SIZE; b++ )
  {
    for ( i = 0; i < sample_size; i++ )
      sample[i] = bsample[8+ (int) (k * ran1(&idum))];

    Sf_compute( Sf, sample, &max_resoln, &sample_size, "Gaussian1" );
    Wf_compute( Wf, Sf, &max_resoln, &sample_size, "Gaussian1" );
    den = denominator( Wf, sample_size );

    for ( j = 1; j <= max_resoln; j++ )
      (*histo)[j][b] = numerator( Wf, j, sample_size ) / den;
  }

  for ( j = 1; j <= max_resoln; j++ )
    qcksrt( HISTO_SIZE, (*histo)[j]-1 );

  free( Sf );
  free( Wf );
  free( bsample );
  free( mean );
  free( sample );
}

/***************************************************************************/
/* Compute mallat pvalue */
/***************************************************************************/

void normal_pval_compute(float *pval, float *s, int *max_resoln_ptr,
  int *np_ptr, int *num_of_windows_ptr, int *window_size_ptr )
{
  int max_resoln = *max_resoln_ptr;
  int np = *np_ptr;
  int num_of_windows = *num_of_windows_ptr;
  int window_size = *window_size_ptr;
  char **pfiltername = (char **) malloc(sizeof(char *));

  float *window_data;
  float *Sf;
  float *Wf;
  float **p;

  int step = window_size / 4;

  float T, den, **histo;
  int w, i, j, t;

  if(!(window_data = (float *) malloc( window_size * sizeof(float) )))
    error("Memory allocation failed for window_data in simul.c \n");
  if(!(histo = (float **)malloc((max_resoln + 1) * sizeof(float *))))
    error("Memory allocation failed for histo in simul.c \n");
  if(!(pfiltername = (char **) malloc(sizeof(char *))))
    error("Memory allocation failed for pfiltername in simul.c \n");
  if(!(Sf = (float *) malloc((max_resoln+1)* window_size * sizeof(float))))
    error("Memory allocation failed for *Sf in simul.c \n");
  if(!(Wf = (float *) malloc(max_resoln* window_size * sizeof(float))))
    error("Memory allocation failed for *Wf in simul.c \n");
  if(!(p = (float **) malloc( (max_resoln+1) * sizeof(float *) )))
    error("Memory allocation failed for p in simul.c \n");



  normal_histo( &histo, max_resoln, window_size );
  filename_given(*pfiltername,"Gaussian1");
  for ( j = 1; j <= max_resoln; j++ )
    if(!(p[j] = (float *) malloc( num_of_windows * sizeof(float) )))
      error("Memory failed for p[j] in simul.c ");

  for ( w = 0; w < num_of_windows; w++ )
  {
    for ( i = 0, t = step*w; i < window_size; i++, t++ )
      window_data[i] = s[t];

    Sf_compute(Sf, window_data, &max_resoln, &window_size,pfiltername);
    Wf_compute(Wf, Sf, &max_resoln, &window_size,pfiltername);
    den = denominator( Wf, window_size );

    for ( j = 1; j <= max_resoln; j++ )
    {
      T = numerator( Wf, j, window_size ) / den;
      p[j][w] = p_value( T, histo, j, HISTO_SIZE );
    }
  }

  compute_pval_average( pval, p, max_resoln, np, num_of_windows, window_size );

  free( window_data );
  free( Sf );
  free( Wf );
  for ( j = 1; j <= max_resoln; j++ )
  {
    free( histo[j] );
    free( p[j] );
  }
  free( histo );
  free( p );
}

/***************************************************************************/
/* Compute mallat bootstrap pvalue */
/***************************************************************************/

void bootstrap_pval_compute(float *pval, float *s, int *max_resoln_ptr,
  int *np_ptr, int *num_of_windows_ptr, int *window_size_ptr )
{
  int max_resoln = *max_resoln_ptr;
  int np = *np_ptr;
  int num_of_windows = *num_of_windows_ptr;
  int window_size = *window_size_ptr;

  float *window_data;
  float *Sf;
  float *Wf;
  float **p;
  char **pfiltername;

  int step = window_size / 4;

  float T, den, **histo;
  int w, i, j, offset;

  if(!(window_data = (float *) malloc( window_size * sizeof(float) )))
    error("Memory allocation failed for window_data in simul.c \n");
  if(!(histo = (float **)malloc((max_resoln + 1) * sizeof(float *))))
    error("Memory allocation failed for histo in simul.c \n");
  if(!(pfiltername = (char **) malloc(sizeof(char *))))
    error("Memory allocation failed for pfiltername in simul.c \n");
  if(!(Sf = (float *) malloc((max_resoln+1)* window_size * sizeof(float))))
    error("Memory allocation failed for *Sf in simul.c \n");
  if(!(Wf = (float *) malloc(max_resoln* window_size * sizeof(float))))
    error("Memory allocation failed for *Wf in simul.c \n");
  if(!(p = (float **) malloc( (max_resoln+1) * sizeof(float *) )))
    error("Memory allocation failed for p in simul.c \n");

  bootstrap_histo( &histo, s, max_resoln, window_size );

  for ( j = 1; j <= max_resoln; j++ )
    if(!(p[j] = (float *) malloc( num_of_windows * sizeof(float) )))
      error("Memory allocation failed for p[j] in simul.c \n ");

  filename_given(*pfiltername,"Gaussian1");

  for ( w = 0; w < num_of_windows; w++ )
  {
    for ( i = 0, offset = step*w; i < window_size; i++ )
      window_data[i] = s[offset+i];

    Sf_compute(Sf, window_data, &max_resoln, &window_size,pfiltername);
    Wf_compute(Wf, Sf, &max_resoln, &window_size,pfiltername);

    den = denominator( Wf, window_size );

    for ( j = 1; j <= max_resoln; j++ )
    {
      T = numerator( Wf, j, window_size ) / den;
      p[j][w] = p_value( T, histo, j, HISTO_SIZE );
    }
  }

  compute_pval_average( pval, p, max_resoln, np, num_of_windows, window_size );

  free( window_data );
  free( Sf );
  free( Wf );
  for ( j = 1; j <= max_resoln; j++ )
  {
    free( histo[j] );
    free( p[j] );
  }
  free( histo );
  free( p );
}

/***************************************************************************/
/* Compute mallat normal threshold for trimming */
/***************************************************************************/

void nthresh_compute(float *nthresh, float *s, int *maxresoln_ptr,
  int *sample_size_ptr, float prct )
{
  int max_resoln = *maxresoln_ptr;
  int sample_size = *sample_size_ptr;

  float *mean;
  float *sample;
  float **histo;
  float *Sf;
  float *Wf;
  char **pfiltername;
  
  float var, std;
  int j, b, i, t;

  if(!(histo = (float **)malloc((max_resoln + 1) * sizeof(float *))))
    error("Memory allocation failed for histo in simul.c \n");
  if(!(pfiltername = (char **) malloc(sizeof(char *))))
    error("Memory allocation failed for pfiltername in simul.c \n");
  if(!(mean = (float *) malloc( sample_size * sizeof(float))))
    error("Memory allocation failed for *mean in simul.c \n");
  if(!(sample = (float *) malloc( sample_size * sizeof(float))))
    error("Memory allocation failed for *sample in simul.c \n");
  if(!(Sf = (float *) malloc((max_resoln+1)* sample_size * sizeof(float))))
    error("Memory allocation failed for *Sf in simul.c \n");
  if(!(Wf = (float *) malloc(max_resoln* sample_size * sizeof(float))))
    error("Memory allocation failed for *Wf in simul.c \n");

  /*  printf("Idum = %d\n",idum);  */

  for ( i = 0; i < sample_size; i++ )
    sample[i] = s[i];

  local_mean( mean, sample, sample_size );
  for ( i = 0; i < sample_size; i++ )
    sample[i] -= mean[i];

  var = variance( sample, sample_size );
  std = (float) sqrt( var );

  for ( j = 1; j <= max_resoln; j++ )
    if(!(histo[j] = (float *) malloc( HISTO_SIZE * sizeof(float) )))
      error("Memory allocation failed for histo[i] in simul.c \n");

  if(!(*pfiltername = (char *)malloc(STRING_SIZE * sizeof(char))))
    error("Memory allocation failed for *pfilename in simul.c \n");


  filename_given(*pfiltername,"Gaussian1");

  for ( b = 0; b < HISTO_SIZE; b++ )  {
    for ( i = 0; i < sample_size; i++ )
      sample[i] = std * gasdev( &idum ); 

    Sf_compute(Sf, sample, &max_resoln, &sample_size,pfiltername);
    Wf_compute(Wf, Sf, &max_resoln, &sample_size,pfiltername);

     for ( j = 1; j <= max_resoln; j++ )
    {
      for ( i = 0, t = (j-1)*sample_size; i < sample_size; i++, t++ )
	sample[i] = Wf[t];

      qcksrt( sample_size, sample-1 );
     /* Take the max of the absolute value of psi  */

      histo[j][b] = max( fabs(sample[0]), fabs(sample[sample_size-1]) );
     } 
  }
 
  for ( j = 1; j <= max_resoln; j++ )
  {
    qcksrt( HISTO_SIZE, histo[j]-1 );
    /*    nthresh[j-1] = histo[j][(int)(HISTO_SIZE * 0.95) - 1];   */
    nthresh[j-1] = histo[j][(int)(HISTO_SIZE * prct) - 1];
    free( histo[j] );
  }

  free( histo );
  free( mean );
  free( sample );
  free( Sf );
  free( Wf );
}

/***************************************************************************/
/* Compute mallat bootstrap threshold for trimming */
/***************************************************************************/

void bthresh_compute(float *bthresh, float *s, int *maxresoln_ptr,
  int *sample_size_ptr, float prct )
{
  int max_resoln = *maxresoln_ptr;
  int sample_size = *sample_size_ptr;

  float *mean;
  float *sample;
  float *bsample;
  float **histo;
  float *Sf;
  float *Wf;
  char **pfiltername;
  int k = sample_size - 16;  /* k depends on LOCAL_LENGTH in local_mean */
  int j, b, i, t;

  if(!(histo = (float **)malloc((max_resoln + 1) * sizeof(float *))))
    error("Memory allocation failed for histo in simul.c \n");
  if(!(pfiltername = (char **) malloc(sizeof(char *))))
    error("Memory allocation failed for pfiltername in simul.c \n");
  if(!(mean = (float *) malloc( sample_size * sizeof(float))))
    error("Memory allocation failed for *mean in simul.c \n");
  if(!(sample = (float *) malloc( sample_size * sizeof(float))))
    error("Memory allocation failed for *sample in simul.c \n");
  if(!(bsample = (float *) malloc( sample_size * sizeof(float))))
    error("Memory allocation failed for *bample in simul.c \n");
  if(!(Sf = (float *) malloc((max_resoln+1)* sample_size * sizeof(float))))
    error("Memory allocation failed for *Sf in simul.c \n");
  if(!(Wf = (float *) malloc(max_resoln* sample_size * sizeof(float))))
    error("Memory allocation failed for *Wf in simul.c \n");

  for ( i = 0; i < sample_size; i++ )
    bsample[i] = s[i];
  local_mean( mean, bsample, sample_size );
  for ( i = 0; i < sample_size; i++ )
    bsample[i] -= mean[i];

  for ( j = 1; j <= max_resoln; j++ )
    if(!(histo[j] = (float *) malloc( HISTO_SIZE * sizeof(float) )))
      error("Memory allocation failed for histo[i] in simul.c \n");
  *pfiltername = (char *)malloc(STRING_SIZE * sizeof(char));

  filename_given(*pfiltername,"Gaussian1");

  for ( b = 0; b < HISTO_SIZE; b++ )
  {
    for ( i = 0; i < sample_size; i++ )
      sample[i] = bsample[8+ (int) (k * ran1(&idum))];
    
    Sf_compute(Sf, sample, &max_resoln, &sample_size,pfiltername);
    Wf_compute(Wf, Sf, &max_resoln, &sample_size,pfiltername);

    for ( j = 1; j <= max_resoln; j++ )
    {
      for ( i = 0, t = (j-1)*sample_size; i < sample_size; i++, t++ )
	sample[i] = Wf[t];

      qcksrt( sample_size, sample-1 );
      /* Take the max of the absolute value of psi */
      histo[j][b] = max( fabs(sample[0]), fabs(sample[sample_size-1]) );
    }
  }  

  for ( j = 1; j <= max_resoln; j++ )
  {
    qcksrt( HISTO_SIZE, histo[j]-1 );
    /*    bthresh[j-1] = histo[j][(int)(HISTO_SIZE * 0.95) - 1];   */
    bthresh[j-1] = histo[j][(int)(HISTO_SIZE * prct) - 1];
    free( histo[j] );
  }
  
  free( histo );
  free( mean );
  free( sample );
  free( bsample );
  free( Sf );
  free( Wf );
}
