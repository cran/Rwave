#include <stdlib.h>



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
#include "dyadic.h"
#include "Swave.h"

#define HISTO_SIZE    500

extern long idum;




/**************************************************************************/
/* local_mean */
/**************************************************************************/
#define LOCAL_LENGTH    17

void local_mean(double *mean, double *s, int np )
{
  double *s_sym;
  int i, j, t;
  double sum;

  if(!(s_sym = (double *) R_alloc( 2*np,sizeof(double))))
    error("Memory allocation failed in s_sym at simul.c \n");

  for ( t = 0; t < np; t++ )
    s_sym[t] = s_sym[2*np-t-1] = s[t];

  for ( t = 0; t < np; t++ )
  {
    sum = 0.0;
    for ( j = t-8, i = 0; i < LOCAL_LENGTH; i++, j++ )
      sum += (j < 0) ? s_sym[-j-1] : s_sym[j];
    mean[t] = sum / (double) LOCAL_LENGTH;
  }

  //free( s_sym );
}
#undef LOCAL_LENGTH



/**************************************************************************/
/* variance */
/**************************************************************************/

double variance(double *s, int np )
{
  double mean, sum;
  double *temp;
  int i;

  if(!(temp = (double *)R_alloc(np , sizeof(double))))
    error("Memory allocation failed for temp at simul.c ");

  for ( i = 0, mean = 0.0; i < np; i++ )
    mean += s[i];
  mean /= (double) np;

  for ( i = 0; i < np; i++ )
    temp[i] = s[i] - mean;

  for ( i = 0, sum = 0.0; i < np; i++ )
    sum += temp[i] * temp[i];

  //free( temp );
  return sum / (double) np;
}

/****************************************************************************/
/* Compute p-value                                                          */
/* Def:  p_value = P{ T > T_obs }                                           */
/****************************************************************************/

double p_value(double T, double **histo, int resoln, int histo_size )
{
  int count = 0;
  int i;

  for ( i = 0; i < histo_size; i++ )
    if ( histo[resoln][i] > T )
    {
      count = histo_size - i;  /* num of values that > T */
      break;
    }
  return (double) count / histo_size;
}

/****************************************************************************/
/* Compute p-value average                                                  */
/****************************************************************************/

void compute_pval_average(double *pval, double **p, int max_resoln, int np,
  int num_of_windows, int window_size )

{
  int interval_length = window_size / 4;
  int num_of_intervals = np / interval_length;
  double *temp;
  int m, i, j, k;

  if(!(temp = (double *)R_alloc(num_of_intervals , sizeof(double))))
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
  //free( temp );
}



/****************************************************************************/
/* Gaussian random variables generator                                      */
/****************************************************************************/
double gasdev(long *idum)
{
  static int iset=0;
  static double gset;
  double fac,rsq,v1,v2;

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

double denominator(double *Wf, int np )
{
  double den = 0.0;
  int t;

  for ( t = 0; t < 2*np; t++ )  /* resoln 1 and 2 */
    den += Wf[t] * Wf[t];
  return den;
}

/****************************************************************************/
/* numerator of the statistics = sqrt of (Wf[j]^4)                          */
/****************************************************************************/

double numerator(double *Wf, int resoln, int np )
{
  double sum = 0.0;
  int t, i;

  for ( t = (resoln-1)*np, i = 0; i < np; i++, t++ )
    sum += Wf[t] * Wf[t] * Wf[t] * Wf[t];
  return (double) sqrt( (double) sum );
}

/***************************************************************************/
/* Compute mallat normal histogram for computing p-values */
/***************************************************************************/

void normal_histo( double ***histo, int max_resoln, int sample_size )
{
  double *Sf = (double *) R_alloc( (max_resoln+1) * sample_size , sizeof(double) );
  double *Wf = (double *) R_alloc( max_resoln * sample_size , sizeof(double) );
  double *sample = (double *) R_alloc( sample_size , sizeof(double) );

  double den;
  int b, i, j;


  *histo = (double **) R_alloc( (max_resoln+1) , sizeof(double *) );
  for ( j = 1; j <= max_resoln; j++ )
    (*histo)[j] = (double *) R_alloc( HISTO_SIZE , sizeof(double) );

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

/*
  free( Sf );
  free( Wf );
  free( sample );
*/
}

/***************************************************************************/
/* Compute mallat bootstrap histogram for computing p-values */
/***************************************************************************/

void bootstrap_histo(double ***histo, double *s, int max_resoln,
  int sample_size )
{
  double *Sf = (double *) R_alloc( (max_resoln+1) * sample_size , sizeof(double) );
  double *Wf = (double *) R_alloc( max_resoln * sample_size , sizeof(double) );

  double *sample = (double *) R_alloc( sample_size , sizeof(double) );
  double *bsample = (double *) R_alloc( sample_size , sizeof(double) );
  double *mean = (double *) R_alloc( sample_size , sizeof(double) );

  double den;
  int b, i, j;
  int k = sample_size - 16;  

  for ( i = 0; i < sample_size; i++ )
    bsample[i] = s[i];
  local_mean( mean, bsample, sample_size );
  for ( i = 0; i < sample_size; i++ )
    bsample[i] -= mean[i];

  *histo = (double **) R_alloc( (max_resoln+1) , sizeof(double *) );
  for ( j = 1; j <= max_resoln; j++ )
    (*histo)[j] = (double *) R_alloc( HISTO_SIZE , sizeof(double) );

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

/*
  free( Sf );
  free( Wf );
  free( bsample );
  free( mean );
  free( sample );
*/
}

/***************************************************************************/
/* Compute mallat pvalue */
/***************************************************************************/

void normal_pval_compute(double *pval, double *s, int *max_resoln_ptr,
  int *np_ptr, int *num_of_windows_ptr, int *window_size_ptr )
{
  int max_resoln = *max_resoln_ptr;
  int np = *np_ptr;
  int num_of_windows = *num_of_windows_ptr;
  int window_size = *window_size_ptr;
  char *pfiltername = "Gaussian1";
  //char *pfiltername = (char *) R_alloc(1, sizeof(char *));

  double *window_data;
  double *Sf;
  double *Wf;
  double **p;

  int step = window_size / 4;

  double T, den;
  double **histo;
  int w, i, j, t;

  if(!(window_data = (double *) R_alloc( window_size , sizeof(double) )))
    error("Memory allocation failed for window_data in simul.c \n");
  if(!(histo = (double **)R_alloc((max_resoln + 1) , sizeof(double *))))
    error("Memory allocation failed for histo in simul.c \n");
  //if(!(pfiltername = (char **) R_alloc(1, sizeof(char *))))
   // error("Memory allocation failed for pfiltername in simul.c \n");
  if(!(Sf = (double *) R_alloc((max_resoln+1)* window_size , sizeof(double))))
    error("Memory allocation failed for *Sf in simul.c \n");
  if(!(Wf = (double *) R_alloc(max_resoln* window_size , sizeof(double))))
    error("Memory allocation failed for *Wf in simul.c \n");
  if(!(p = (double **) R_alloc( (max_resoln+1) , sizeof(double *) )))
    error("Memory allocation failed for p in simul.c \n");



  normal_histo( &histo, max_resoln, window_size );
  //filename_given(*pfiltername,"Gaussian1");
  for ( j = 1; j <= max_resoln; j++ )
    if(!(p[j] = (double *) R_alloc( num_of_windows , sizeof(double) )))
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

/*
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
*/
}

/***************************************************************************/
/* Compute mallat bootstrap pvalue */
/***************************************************************************/

void bootstrap_pval_compute(double *pval, double *s, int *max_resoln_ptr,
  int *np_ptr, int *num_of_windows_ptr, int *window_size_ptr )
{
  int max_resoln = *max_resoln_ptr;
  int np = *np_ptr;
  int num_of_windows = *num_of_windows_ptr;
  int window_size = *window_size_ptr;

  double *window_data;
  double *Sf;
  double *Wf;
  double **p;
  char *pfiltername = "Gaussian1";

  int step = window_size / 4;

  double T, den;
  double **histo;
  int w, i, j, offset;

  if(!(window_data = (double *) R_alloc( window_size , sizeof(double) )))
    error("Memory allocation failed for window_data in simul.c \n");
  if(!(histo = (double **)R_alloc((max_resoln + 1) , sizeof(double *))))
    error("Memory allocation failed for histo in simul.c \n");
  //if(!(pfiltername = (char **) R_alloc(1, sizeof(char *))))
    error("Memory allocation failed for pfiltername in simul.c \n");
  if(!(Sf = (double *) R_alloc((max_resoln+1)* window_size , sizeof(double))))
    error("Memory allocation failed for *Sf in simul.c \n");
  if(!(Wf = (double *) R_alloc(max_resoln* window_size , sizeof(double))))
    error("Memory allocation failed for *Wf in simul.c \n");
  if(!(p = (double **) R_alloc( (max_resoln+1) , sizeof(double *) )))
    error("Memory allocation failed for p in simul.c \n");

  bootstrap_histo( &histo, s, max_resoln, window_size );

  for ( j = 1; j <= max_resoln; j++ )
    if(!(p[j] = (double *) R_alloc( num_of_windows , sizeof(double) )))
      error("Memory allocation failed for p[j] in simul.c \n ");

  //filename_given(*pfiltername,"Gaussian1");

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

/*
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
*/
}

/***************************************************************************/
/* Compute mallat normal threshold for trimming */
/***************************************************************************/

void nthresh_compute(double *nthresh, double *s, int *maxresoln_ptr,
  int *sample_size_ptr, double prct )
{
  int max_resoln = *maxresoln_ptr;
  int sample_size = *sample_size_ptr;

  double *mean;
  double *sample;
  double **histo;
  double *Sf;
  double *Wf;
  char *pfiltername = "Gaussian1";
  
  double var, std;
  int j, b, i, t;

  if(!(histo = (double **)R_alloc((max_resoln + 1) , sizeof(double *))))
    error("Memory allocation failed for histo in simul.c \n");
  //if(!(pfiltername = (char **) R_alloc(1, sizeof(char *))))
    //error("Memory allocation failed for pfiltername in simul.c \n");
  if(!(mean = (double *) R_alloc( sample_size , sizeof(double))))
    error("Memory allocation failed for *mean in simul.c \n");
  if(!(sample = (double *) R_alloc( sample_size , sizeof(double))))
    error("Memory allocation failed for *sample in simul.c \n");
  if(!(Sf = (double *) R_alloc((max_resoln+1)* sample_size , sizeof(double))))
    error("Memory allocation failed for *Sf in simul.c \n");
  if(!(Wf = (double *) R_alloc(max_resoln* sample_size , sizeof(double))))
    error("Memory allocation failed for *Wf in simul.c \n");

  /*  printf("Idum = %d\n",idum);  */

  for ( i = 0; i < sample_size; i++ )
    sample[i] = s[i];

  local_mean( mean, sample, sample_size );
  for ( i = 0; i < sample_size; i++ )
    sample[i] -= mean[i];

  var = variance( sample, sample_size );
  std = (double) sqrt( var );

  for ( j = 1; j <= max_resoln; j++ )
    if(!(histo[j] = (double *) R_alloc( HISTO_SIZE , sizeof(double) )))
      error("Memory allocation failed for histo[i] in simul.c \n");

  //if(!(*pfiltername = (char *)R_alloc(STRING_SIZE , sizeof(char))))
   // error("Memory allocation failed for *pfilename in simul.c \n");


  //filename_given(*pfiltername,"Gaussian1");

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
    //free( histo[j] );
  }

/*
  free( histo );
  free( mean );
  free( sample );
  free( Sf );
  free( Wf );
*/
}

/***************************************************************************/
/* Compute mallat bootstrap threshold for trimming */
/***************************************************************************/

void bthresh_compute(double *bthresh, double *s, int *maxresoln_ptr,
  int *sample_size_ptr, double prct )
{
  int max_resoln = *maxresoln_ptr;
  int sample_size = *sample_size_ptr;

  double *mean;
  double *sample;
  double *bsample;
  double **histo;
  double *Sf;
  double *Wf;
  char *pfiltername = "Gaussian1";
  int k = sample_size - 16;  /* k depends on LOCAL_LENGTH in local_mean */
  int j, b, i, t;

  histo = (double **) R_alloc((max_resoln + 1), sizeof(double *));
  //if(!(pfiltername = (char **) R_alloc(1, sizeof(char *))))
    //error("Memory allocation failed for pfiltername in simul.c \n");
  if(!(mean = (double *) R_alloc( sample_size, sizeof(double))))
    error("Memory allocation failed for *mean in simul.c \n");
  if(!(sample = (double *) R_alloc( sample_size , sizeof(double))))
    error("Memory allocation failed for *sample in simul.c \n");
  if(!(bsample = (double *) R_alloc( sample_size , sizeof(double))))
    error("Memory allocation failed for *bample in simul.c \n");
  if(!(Sf = (double *) R_alloc((max_resoln+1)* sample_size , sizeof(double))))
    error("Memory allocation failed for *Sf in simul.c \n");
  if(!(Wf = (double *) R_alloc(max_resoln* sample_size , sizeof(double))))
    error("Memory allocation failed for *Wf in simul.c \n");

  for ( i = 0; i < sample_size; i++ )
    bsample[i] = s[i];
  local_mean( mean, bsample, sample_size );
  for ( i = 0; i < sample_size; i++ )
    bsample[i] -= mean[i];

  for ( j = 1; j <= max_resoln; j++ )
    if(!(histo[j] = (double *) R_alloc( HISTO_SIZE , sizeof(double) )))
      error("Memory allocation failed for histo[i] in simul.c \n");
  //*pfiltername = (char *)R_alloc(STRING_SIZE , sizeof(char));

  //filename_given(*pfiltername,"Gaussian1");

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
    //free( histo[j] );
  }
  
  /*free( histo );
  free( mean );
  free( sample );
  free( bsample );
  free( Sf );
  free( Wf );
  */
}
