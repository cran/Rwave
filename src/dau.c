
/* #include "wavelet.h" */
#include "dau_wave.h"
#include "pvalue.h"


/******************************************************************
*              (c) Copyright  1997                                *
*                         by                                      *
*  Author: Rene Carmona, Bruno Torresani, Wen L. Hwang, A. Wang   *
*                 Princeton University                            *
*                 All right reserved                              *
******************************************************************/




/****************************************************************************/

#define HISTO_SIZE    500

/****************************************************************************/



double *a, **c;
int NW, *twoto;

/****************************************************************************/
/*		COMPUTE WAVELET COEFFICIENT RANGE FOR ALL RESOLUTIONS	    */
/****************************************************************************/

void compute_d_phi_range_for_all_resoln( d_phi_range, max_resoln, np )
bound **d_phi_range;
int max_resoln, np;
{
  int j;

  *d_phi_range = (bound *) malloc( (max_resoln+1) * sizeof(bound) );
  for ( j = 0; j <= max_resoln; j++ )
  {
    (*d_phi_range)[j].lb = (int) ceil((1 - 1.0 / twoto[j]) * (1 - 2*NW));
    (*d_phi_range)[j].ub = (int) ((np -1) / twoto[j]);
    (*d_phi_range)[j].size = (*d_phi_range)[j].ub - (*d_phi_range)[j].lb +1;
  }
}

/****************************************************************************/
/*         COMPUTE WAVELET COEFFICIENT RANGE FOR ALL RESOLUTIONS	    */
/****************************************************************************/

void compute_d_psi_range_for_all_resoln( d_psi_range, d_phi_range, max_resoln, np )
bound **d_psi_range;
bound *d_phi_range;
int max_resoln;
int np;
{
  int j;

  *d_psi_range = (bound *) malloc( (max_resoln+1) * sizeof(bound) );
  for ( j = 1; j <= max_resoln; j++ )
  {
    (*d_psi_range)[j].lb = (int) ceil( (double) ((d_phi_range[j-1].lb -1) / 2) );
    (*d_psi_range)[j].ub = (int) ( d_phi_range[j-1].ub / 2 + NW -1 );
    (*d_psi_range)[j].size = (*d_psi_range)[j].ub - (*d_psi_range)[j].lb +1;
  }
}

/****************************************************************************/
/*		COMPUTE PHI COEFFICIENTS FOR ALL RESOLUTIONS                */
/****************************************************************************/

void compute_d_phi_for_all_resoln( d_phi, d_phi_range, s, max_resoln )
float **d_phi;
bound *d_phi_range;
float *s;
int max_resoln;
{
  int j, k, n_min, n_max, n;
  double sum;

  for ( j = 0; j <= max_resoln; j++ )
  {
    d_phi[j] = (float *) malloc( d_phi_range[j].size * sizeof(float) );
		
    if ( j == 0 )
      for ( k = d_phi_range[j].lb; k <= d_phi_range[j].ub; k++ )
	d_phi[j][k] = s[k];
    else
      for ( k = d_phi_range[j].lb; k <= d_phi_range[j].ub; k++ )
      {
	n_min = max( d_phi_range[j-1].lb, 2*k );
	n_max = min( d_phi_range[j-1].ub, 2*k+2*NW-1 );
	
	sum = 0.0;
	for ( n = n_min; n <= n_max; n++ )
	  sum += c[NW][n - 2*k] * d_phi[j-1][n - d_phi_range[j-1].lb];
	d_phi[j][k - d_phi_range[j].lb] = (float) sum;
      }
  }   
}

/****************************************************************************/
/*		COMPUTE PSI COEFFICIENTS FOR ALL RESOLUTIONS       	    */
/****************************************************************************/

void compute_d_psi_for_all_resoln( d_psi, d_psi_range, d_phi, d_phi_range,
				  max_resoln )
float **d_psi;
bound *d_psi_range;
float **d_phi;
bound *d_phi_range;
int max_resoln;
{
  int j, k, n_min, n_max, n;
  double sum;

  for ( j = 1; j <= max_resoln; j++ )
  {
    d_psi[j] = (float *) malloc( d_psi_range[j].size * sizeof( float ) );
		
    for ( k = d_psi_range[j].lb; k <= d_psi_range[j].ub; k++ )
    {
      n_min = max( d_phi_range[j-1].lb, 2*(k-NW+1) );
      n_max = min( d_phi_range[j-1].ub, 2*k+1 );

      sum = 0.0;
      for ( n = n_min; n <= n_max; n++ )
	sum += minus1to(n) * c[NW][2*k+1-n] * d_phi[j-1][n - d_phi_range[j-1].lb];
      d_psi[j][k - d_psi_range[j].lb] = (float) sum;
    }
  }
}

/****************************************************************************/
/*		PHI RECONSTRUCTION					    */
/****************************************************************************/

void phi_reconstruction( phi, d_phi, phi_array, d_phi_range, max_resoln, np )
float *phi, **d_phi, *phi_array;
bound *d_phi_range;
int max_resoln, np;
{
  int j, t, k_min, k_max, k;
  float sum, two_to_j, two_to_j_half, two_to_j_times_t;
	
  for ( j = 0; j <= max_resoln; j++ )
  {
    two_to_j = 1.0 / (float)pow( 2.0, (double) j );	  /* j = -m */
    two_to_j_half = 1.0 / (float)pow( 2.0, ((double)j)/2 );
		
    for( t = 0; t < np; t++ )
    {
      two_to_j_times_t = two_to_j * t;
      k_min = max((int) ceil((double)(two_to_j_times_t - 2*NW +1)), d_phi_range[j].lb);
      k_max = (int) floor( (double) two_to_j_times_t );
			
      /* NOTE: k_max <= d_phi_range[j].ub */
      sum = 0.0;
      for ( k = k_min; k <= k_max; k++ )
	sum += d_phi[j][k - d_phi_range[j].lb] * 
	       phi_array[(int) ((two_to_j_times_t - k) * twoto[max_resoln])];

      phi[j*np+t] = two_to_j_half * sum;
    }
  }
}

/****************************************************************************/
/*		PSI RECONSTRUCTION					    */
/****************************************************************************/

void psi_reconstruction( psi, d_psi, psi_array, d_psi_range, max_resoln, np )
float *psi, **d_psi, *psi_array;
bound *d_psi_range;
int max_resoln, np;
{
  int j, t, k_min, k_max, k;
  float sum, two_to_j, two_to_j_half, two_to_j_times_t;
	
  for ( j = 1; j <= max_resoln; j++ )
  {
    two_to_j = 1.0 / (float) pow( 2.0, (float) j );     /* j = -m */
    two_to_j_half = 1.0 / (float) pow( 2.0, ((float) j)/2 );
		
    for( t = 0; t < np; t++ )
    {
      two_to_j_times_t = two_to_j * t;
      k_min = max( (int) ceil( (double)(two_to_j_times_t - NW +1) ), d_psi_range[j].lb );
      k_max = min( (int) floor((double)(two_to_j_times_t + NW) ), d_psi_range[j].ub );

      sum = 0.0;
      for ( k = k_min; k <= k_max; k++ )
	sum += d_psi[j][k - d_psi_range[j].lb] * 
               psi_array[(int) ((two_to_j_times_t - k + NW) * twoto[max_resoln])];

      psi[(j-1)*np+t] = two_to_j_half * sum;
    }
  }
}

/****************************************************************************/
/*  daubechies_reconst, called by Splus                                     */
/****************************************************************************/

void daubechies_wt( phi, psi, s, NW_ptr, maxresoln_ptr, np_ptr )
float *phi;  /* (maxresoln+1) by np, where np is a power of 2 */
float *psi;  /* maxresoln by np, where np is a power of 2 */
float *s;
int *NW_ptr;
int *maxresoln_ptr;
int *np_ptr;
{
  int max_resoln = *maxresoln_ptr;
  int np = *np_ptr;
  int num_of_resoln = max_resoln + 1;

  bound *d_phi_range, *d_psi_range;
  float **d_phi, **d_psi, *phi_array, *psi_array;
  int j;
  
  NW = *NW_ptr;
  open_read();
  compute_a();
  init_twoto( max_resoln );

  d_psi_range = (bound *) malloc( num_of_resoln * sizeof(bound) );
  d_phi = (float **) malloc( num_of_resoln * sizeof(float *) );
  d_psi = (float **) malloc( num_of_resoln * sizeof(float *) );

  init_phi_array( &phi_array, max_resoln );
  init_psi_array( &psi_array, max_resoln );

  compute_d_phi_range_for_all_resoln( &d_phi_range, max_resoln, np );
  compute_d_psi_range_for_all_resoln( &d_psi_range, d_phi_range, max_resoln, np );

  compute_d_phi_for_all_resoln( d_phi, d_phi_range, s, max_resoln );
  compute_d_psi_for_all_resoln( d_psi, d_psi_range, d_phi, d_phi_range, max_resoln );

  phi_reconstruction( phi, d_phi, phi_array, d_phi_range, max_resoln, np );
  psi_reconstruction( psi, d_psi, psi_array, d_psi_range, max_resoln, np );

  free( twoto );
  free( phi_array );
  free( psi_array );
  free( d_phi_range );
  free( d_psi_range );
  for ( j = 0; j <= max_resoln; j++ )
    free( d_phi[j] );
  free( d_phi );
  for ( j = 1; j <= max_resoln; j++ )
    free( d_psi[j] );
  free( d_psi );

  free( a );
  for ( j = NMIN; j < NMAX + 1; j++ )
    free( c[j] );
  free( c );
}

/****************************************************************************/
/*                                                                          */
/*          Discrete Daubechies Wavelet Transform                           */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/* Compute dH filter */
/****************************************************************************/

void compute_dH_bound( dH_bound, max_resoln )
bound **dH_bound;
int max_resoln;
{
  int j;
  int temp = 2*NW-1;

  *dH_bound = (bound *) malloc( max_resoln * sizeof(bound) );

  for ( j = 0; j < max_resoln; j++ )
  {
    (*dH_bound)[j].lb = 0; 
    (*dH_bound)[j].ub = twoto[j]*temp;
    (*dH_bound)[j].size = (*dH_bound)[j].ub - (*dH_bound)[j].lb + 1;
  }
/*
  for ( j = 0; j < max_resoln; j++ )
    printf("dH_bound[%d] = [%d, %d]\n", j, (*dH_bound)[j].lb, (*dH_bound)[j].ub );
  printf("\n");
*/
}

/****************************************************************************/
/* Compute dG filter */
/****************************************************************************/

void compute_dG_bound( dG_bound, max_resoln )
bound **dG_bound;
int max_resoln;
{
  int j;
  int temp = 2 - 2*NW;

  *dG_bound = (bound *) malloc( max_resoln * sizeof(bound) );

  for ( j = 0; j < max_resoln; j++ )
  {
    (*dG_bound)[j].lb = twoto[j] * temp;
    (*dG_bound)[j].ub = twoto[j];
    (*dG_bound)[j].size = (*dG_bound)[j].ub - (*dG_bound)[j].lb + 1;
  }
/*
  for ( j = 0; j < max_resoln; j++ )
    printf("dG_bound[%d] = [%d, %d]\n", j, (*dG_bound)[j].lb, (*dG_bound)[j].ub );
  printf("\n");
*/
}

/****************************************************************************/
/* Compute dH filter for each resolution */
/****************************************************************************/

void compute_dH( dH, dH_bound, max_resoln )
float ***dH;
bound *dH_bound;
int max_resoln;
{
  int j, i;

  *dH = (float **) malloc( max_resoln * sizeof(float *) );
  for ( j = 0; j < max_resoln; j++ )
  {
    (*dH)[j] = (float *) malloc( dH_bound[j].size * sizeof(float) );
    if ( j == 0 )
    {
      for ( i = 0; i < dH_bound[j].size; i++ )
	(*dH)[j][i] = c[NW][i];
    }  
    else
    {
      for ( i = 0; i < dH_bound[j].size; i++ )  /* insert zeros */
	(*dH)[j][i] = (i % 2) == 0 ? (*dH)[j-1][i/2] : 0.0;
    }

  }
/*
  for ( j = 0; j < max_resoln; j++ )
  {
    for ( i = 0; i < dH_bound[j].size; i++ )
      printf("%f\n", (*dH)[j][i] );
    printf("\n");
  }
*/
}

/****************************************************************************/
/* Compute dG filter for each resolution */
/****************************************************************************/

void compute_dG( dG, dG_bound, max_resoln )
float ***dG;
bound *dG_bound;
int max_resoln;
{
  int j, i, n;

  *dG = (float **) malloc( max_resoln * sizeof(float *) );
  for ( j = 0; j < max_resoln; j++ )
  {
    (*dG)[j] = (float *) malloc( dG_bound[j].size * sizeof(float) );
    if ( j == 0 )
    { 
      for ( i = 0, n = 2-2*NW; i < dG_bound[j].size; i++, n++ )
	(*dG)[j][i] = minus1to(n) * c[NW][-n+1];
    }  
    else
    {
      for ( i = 0; i < dG_bound[j].size; i++ )  /* insert zeros */
	(*dG)[j][i] = (i % 2) == 0 ? (*dG)[j-1][i/2] : 0.0;
    }
  }
/*
  for ( j = 0; j < max_resoln; j++ )
  {
    for ( i = 0; i < dG_bound[j].size; i++ )
      printf("%f\n", (*dG)[j][i] );
    printf("\n");
  }
*/
}

/****************************************************************************/
/* COMPUTE DDWAVE (Discrete Daubechies wavelet) */
/****************************************************************************/

void compute_ddwave( phi, psi, s, max_resoln_ptr, np_ptr, NW_ptr )
float *phi;
float *psi;
float *s;
int *max_resoln_ptr;
int *np_ptr;
int *NW_ptr;
{
  int max_resoln = *max_resoln_ptr;
  int np = *np_ptr;

  bound *dH_bound, *dG_bound;
  float **dH, **dG;
  float *sym = (float *) malloc( 2*np * sizeof(float) );
  int j, n, k, t;
  float sum;

  NW = *NW_ptr;
  open_read();
  init_twoto( max_resoln );

  compute_dH_bound( &dH_bound, max_resoln );
  compute_dG_bound( &dG_bound, max_resoln );

  compute_dH( &dH, dH_bound, max_resoln );
  compute_dG( &dG, dG_bound, max_resoln );

  for ( j = 0; j <= max_resoln; j++ )
  {
    if ( j == 0 )
      for ( n = 0; n < np; n++ )
	phi[n] = s[n];
    else
    {
      t = (j-1)*np;
      for ( n = 0; n < np; n++ )
      {
	for ( k = dH_bound[j-1].lb, sum = 0.0; k <= dH_bound[j-1].ub; k++ )
	  sum += dH[j-1][k] * phi[t+(n-k+np)%np];
	phi[j*np+n] = sum;
      }
    }
  }

  for ( j = 1; j <= max_resoln; j++ )
  {
    t = (j-1)*np;
    for ( n = 0; n < np; n++ )
    {
      for ( k = dG_bound[j-1].lb, sum = 0.0; k <= dG_bound[j-1].ub; k++ )
	sum += dG[j-1][k-dG_bound[j-1].lb] * phi[t+(n-k+np)%np];
      psi[t+n] = sum;
    }
  }

  free( twoto );
  free( sym );
  free( dH_bound );
  free( dG_bound );
  for ( j = 0; j < max_resoln; j++ )
  {
    free( dH[j] );
    free( dG[j] );
  }
  free( dH );
  free( dG );
  for ( j = NMIN; j < NMAX + 1; j++ )
    free( c[j] );
  free( c );
}

