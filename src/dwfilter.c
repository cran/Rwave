
/* *****************************************************************
*              (c) Copyright  1997                                 *
*                         by                                       *
* Author: Rene Carmona, Bruno Torresani, Wen L. Hwang and A. Wang  *
*                 Princeton University                             *
*                 All right reserved                               *
*******************************************************************/


#include "dau_wave.h"

/****************************************************************
*  Function: iexp2:
*  ----------------
*  returning interger which is the jth power of 2
*
*    j: interger
*
****************************************************************/

/* 2^j  (j > 0) */
iexp2(j)
  int j;
{
  if(j == 0) return(1);

  return( 1 << j);
}

/****************************************************************
*  Function: Hfilter_compute:
*  --------------------------
*  Computation of H filter
*
*    filename: filter name
*    H: H filters
*    H_bound: size of H filters
*    max_resoln: number of decomposition
*
****************************************************************/

void Hfilter_compute(filtername,H,H_bound,max_resoln )
     char *filtername;
     float ***H;
     bound *H_bound;
     int max_resoln;
{
  int j, i;

  if(!(*H = (float **) malloc( (max_resoln+1) * sizeof(float *))))
    error("Memory allocation failed for *H in filter.c \n");

  for ( j = 0; j <= max_resoln; j++ )  {
    if(!((*H)[j] = (float *) malloc( H_bound[j].size * sizeof(float))))
      error("Memory allocation failed for H[] in filter.c \n");
    signal_zero((*H)[j],H_bound[j].size);

    if ( j == 0 )    {
      if(strcmp(filtername,"Haar") == 0) {
	(*H)[0][0] = 0.5; (*H)[0][1] = 0.5; 
      }
      else {
	(*H)[0][0] = 0.125; (*H)[0][1] = 0.375; 
	(*H)[0][2] = 0.375; (*H)[0][3] = 0.125;
      }
    }
    else    {
      for ( i = 0; i < H_bound[j-1].size; i++ )
	(*H)[j][2 * i] = (*H)[j-1][i];
    }
  }
}

/****************************************************************
*  Function: Gfilter_compute:
*  --------------------------
*  Computation of G filter
*
*    filename: filter name
*    G: G filters
*    G_bound: size of G filters
*    max_resoln: number of decomposition
*
****************************************************************/

void Gfilter_compute(filtername,G,G_bound,max_resoln)
     char *filtername;
     float ***G;
     bound *G_bound;
     int max_resoln;
{
  int j, i;

  if(!(*G = (float **) malloc( (max_resoln+1) * sizeof(float *))))
    error("Memory allocation failed for G in filter.c \n");    

  for ( j = 0; j <= max_resoln; j++ )  {

    if(!((*G)[j] = (float *) malloc( G_bound[j].size * sizeof(float))))
      error("Memory allocation failed for G[] in filter.c \n");    
    signal_zero((*G)[j],G_bound[j].size);

    if ( j == 0 )    {
      if(strcmp(filtername,"Haar") == 0) {
	(*G)[0][0] = 0.5; 
	(*G)[0][1] = -0.5;  
      }
      else {
	(*G)[0][0] = 0.5; 
	(*G)[0][1] = -0.5;  
      }
    }
    else      {
      for ( i = 0; i < G_bound[j-1].size; i++ )
        (*G)[j][2 * i] = (*G)[j-1][i];
    }
  }
}


/****************************************************************
*  Function: HGfilter_bound:
*  -------------------------
*  Computation of filter bound for each resolution
*
*    filterename: filter name
*    H_bound: size of H filters
*    G_bound: size of G filters
*    max_resoln: number of decomposition
*
****************************************************************/


void HGfilter_bound(filtername,H_bound,G_bound,max_resoln )
     char *filtername;
     bound **H_bound, **G_bound;
     int max_resoln;
{
  int j;
  int iexp2();

  if(!(*H_bound = (bound *) malloc( (max_resoln+1) * sizeof(bound) )))
    error("Memory allocation failed for *H_bound in filter.c \n");    
  if(!(*G_bound = (bound *) malloc( (max_resoln+1) * sizeof(bound) )))
    error("Memory allocation failed for *G_bound in filter.c \n");    


  for ( j = 0; j <= max_resoln; j++ )  {
    if(strcmp(filtername,"Haar") == 0) {
      if(j == 0)  {
	(*H_bound)[j].lb = 0;
	(*H_bound)[j].ub = 1; 
	(*H_bound)[j].size = (*H_bound)[j].ub - (*H_bound)[j].lb + 1;

	(*G_bound)[j].lb =  0;
	(*G_bound)[j].ub =  1; 
	(*G_bound)[j].size = (*G_bound)[j].ub - (*G_bound)[j].lb + 1;
      }
      else {
	(*H_bound)[j].lb = -iexp2(j-1);
	(*H_bound)[j].ub = iexp2(j-1);
	(*H_bound)[j].size = (*H_bound)[j].ub - (*H_bound)[j].lb + 1;

	(*G_bound)[j].lb = -iexp2(j-1);
	(*G_bound)[j].ub = iexp2(j-1);
	(*G_bound)[j].size = (*G_bound)[j].ub - (*G_bound)[j].lb + 1;
      }
    }
    else { 
      if(j == 0)  {
	(*H_bound)[j].lb = -1;
	(*H_bound)[j].ub = 2; 

	(*H_bound)[j].size = (*H_bound)[j].ub - (*H_bound)[j].lb + 1;

	(*G_bound)[j].lb =  0;
	(*G_bound)[j].ub =  1; 
	(*G_bound)[j].size = (*G_bound)[j].ub - (*G_bound)[j].lb + 1;
      }
      else {
	(*H_bound)[j].lb = -iexp2(j-1) * 3;
	(*H_bound)[j].ub = 3 * iexp2(j-1); /*  changed 2/1/94 j to j-1  */
	(*H_bound)[j].size = (*H_bound)[j].ub - (*H_bound)[j].lb + 1;

	(*G_bound)[j].lb = -iexp2(j-1);
	(*G_bound)[j].ub = iexp2(j-1); /*  changed 2/1/94 j to j-1  */
	(*G_bound)[j].size = (*G_bound)[j].ub - (*G_bound)[j].lb + 1;
      }
    }
  }

}



/****************************************************************
*  Function: HG_hat_compute:
*  -------------------------
*  Computation of the Fourier transform of H and G
*  H_hat and G_hat are obtained from analytical functions
*
*    filterename: filter name
*    H_hat: the Fourier transform of H filters
*    G_hat: the Fourier transform of G filters
*    max_resoln: number of decomposition
*    np: signal size
*
****************************************************************/

#define pi 3.141592653589793

void HG_hat_compute(filtername,H_hat,G_hat,max_resoln,np)
     float ***H_hat;
     float ***G_hat;
     int max_resoln;
     int np;
     char *filtername;
{
  double temp;
  double arg;
  int m, j, t;
  int iexp2();

  if(strcmp(filtername,"Gaussian1") != 0) {
    printf("Need Gaussian1 filter \n");
    return;
  }

/*  printf("computing H_hat & G_hat with Gaussian1 filter\n");  */
  if(!(*H_hat = (float **) malloc( (max_resoln+1) * sizeof(float) )))
    error("Memory allocation failed for *H_hat in filter.c \n");
  if(!(*G_hat = (float **) malloc( (max_resoln+1) * sizeof(float) )))
    error("Memory allocation failed for *G_hat in filter.c \n");
  
  for ( j = 0; j <= max_resoln; j++ )  {
    if(!((*H_hat)[j] = (float *) malloc( 2*(np+1) * sizeof(float) )))
      error("Memory allocation failed for *H_hat[] in filter.c \n");
    if(!((*G_hat)[j] = (float *) malloc( 2*(np+1) * sizeof(float) )))
      error("Memory allocation failed for *G_hat[] in filter.c \n");
    
    if ( j == 0 )    {
      temp = pi / (double) np;
      
      for ( m = 0; m < np; m++ )      {
	arg = (double) m * temp;
	(*H_hat)[j][2*m]    = (float)(cos(arg) * cos(arg) * cos(arg) * cos(arg));
	(*H_hat)[j][2*m+1] = (float)(cos(arg) * cos(arg) * cos(arg) * sin(arg));
	
	/* remove 4 in the following due to normalize G in L1 */
	(*G_hat)[j][2*m]   = (float)(sin(arg) * sin(arg));
	(*G_hat)[j][2*m+1] = -(float)(sin(arg) * cos(arg));
      }
    }
    else {
      temp = iexp2(j) * pi / (double) np;
      
      for ( m = 0; m < np; m++ )      {
	arg = (double) m * temp;
	(*H_hat)[j][2*m]    = (float)(cos(arg) * cos(arg) * cos(arg));
	(*H_hat)[j][2*m+1] = 0.0;
	  
	/* remove 4 in the following due to normalize G in L1 */
	(*G_hat)[j][2*m]   = 0.0;
	(*G_hat)[j][2*m+1] = (float)(-sin(arg));
      }
    }
  }
}

#undef pi

/****************************************************************
*  Function: Sfilter_compute:
*  --------------------------
*  Computation of the S filters
*
*    filterename: filter name
*    S: S filters
*    S_bound: the size of S filters
*    max_resoln: number of decomposition
*
****************************************************************/

void Sfilter_compute(filtername,S,S_bound,max_resoln)
     char *filtername;
     float ***S;
     bound *S_bound;
     int max_resoln;
{
  int j, i;

  if(!(*S = (float **) malloc( (max_resoln+1) * sizeof(float *))))
    error("Memory allocation failed for *S in filter.c \n");

  for ( j = 0; j <= max_resoln; j++ )  {
    if(!((*S)[j] = (float *) malloc( S_bound[j].size * sizeof(float))))
      error("Memory allocation failed for S[] in filter.c \n");
    signal_zero((*S)[j], S_bound[j].size);

    if ( j == 0 )    {
      if(strcmp(filtername,"Haar") == 0) {
	(*S)[0][0] = 0.5; 
	(*S)[0][1] = 0.5; 
      }
      else {
	(*S)[0][0] = 0.125; 
	(*S)[0][1] = 0.375; 
	(*S)[0][2] = 0.375; 
	(*S)[0][3] = 0.125; 
      }
    }
    else    {
      for ( i = 0; i < S_bound[j-1].size; i++ )
        (*S)[j][2 * i] = (*S)[j-1][i];
    }
  }

  /* begugging ..... */
/*
  {
    FILE *fp = fopen( "Sfilter", "w" );
    int j, i;
    
    for ( j = 0; j <= max_resoln; j++ )    {
      fprintf( fp, "j = %d\n", j );
      for ( i = 0; i < S_bound[j].size; i++ )
	fprintf( fp, "%.3f\n", (*S)[j][i] );
    }
    fclose( fp );
  }
*/  
}

/****************************************************************
*  Function: Kfilter_compute:
*  --------------------------
*  Computation of K filters
*
*    filterename: filter name
*    K: K filters
*    K_bound: the size of K filters
*    max_resoln: number of decomposition
*
****************************************************************/


void Kfilter_compute(filtername,K,K_bound,max_resoln)
  char *filtername;
  float ***K;
  bound *K_bound;
  int max_resoln;
{
  int j, i;

  if(!(*K = (float **) malloc( (max_resoln+1) * sizeof(float *))))
    error("Memory allocation failed for K in filter.c \n");    

  for ( j = 0; j <= max_resoln; j++ )  {

    if(!((*K)[j] = (float *) malloc( K_bound[j].size * sizeof(float))))
      error("Memory allocation failed for K[] in filter.c \n");    
    signal_zero((*K)[j], K_bound[j].size);

    if ( j == 0 )    {
      if(strcmp(filtername,"Haar") == 0) {
	(*K)[0][0] = -0.5; 
	(*K)[0][1] = 0.5; 
      }
      else {
	(*K)[0][0] = -0.03125; 
	(*K)[0][1] = -0.21875; 
	(*K)[0][2] = -0.6875;
	(*K)[0][3] = 0.6875;  
	(*K)[0][4] = 0.21875; 
	(*K)[0][5] = 0.03125; 
      }
    }
    else      {
      for ( i = 0; i < K_bound[j-1].size; i++ )
        (*K)[j][2 * i] = (*K)[j-1][i];
    }
  }	

  /* degugging .... */
/*
  {
    FILE *fp = fopen( "Kfilter", "w" );
    int j, i;
    
    for ( j = 0; j <= max_resoln; j++ )    {
      fprintf( fp, "j = %d\n", j );
      for ( i = 0; i < K_bound[j].size; i++ )
	fprintf( fp, "%f\n", (*K)[j][i] );
    }
    fclose( fp );
  }
*/
}

/****************************************************************
*  Function: Lfilter_compute:
*  --------------------------
*  Computation of L filters
*
*    filterename: filter name
*    L: L filters
*    L_bound: the size of L filters
*    max_resoln: number of decomposition
*
****************************************************************/

void Lfilter_compute(filtername,L,L_bound,max_resoln)
     char *filtername;
     float ***L;
     bound *L_bound;
     int max_resoln;
{
  int j, i;

  if(!(*L = (float **) malloc( (max_resoln+1) * sizeof(float *))))
    error("Memory allocation failed for L in filter.c \n");    

  for ( j = 0; j <= max_resoln; j++ )  {
    if(!((*L)[j] = (float *) malloc( L_bound[j].size * sizeof(float))))
      error("Memory allocation failed for L[] in filter.c \n");    
    signal_zero((*L)[j], L_bound[j].size);

    if ( j == 0 )    {
      if(strcmp(filtername,"Haar") == 0) {
	(*L)[0][0] = 0.125; 
	(*L)[0][1] = 0.75;
	(*L)[0][2] = 0.125;
      }
      else {
	(*L)[0][0] = 0.0078125; 
	(*L)[0][1] = 0.046875; 
	(*L)[0][2] = 0.1171875;
	(*L)[0][3] = 0.65625;  
	(*L)[0][4] = 0.1171875; 
	(*L)[0][5] = 0.046875; 
	(*L)[0][6] = 0.0078125;
      }
    }
    else      {
      for ( i = 0; i < L_bound[j-1].size; i++ )
        (*L)[j][2 * i] = (*L)[j-1][i];
    }
  }

  /* degugging .... */
/*
  {
    FILE *fp = fopen( "Lfilter", "w" );
    int j, i;
    
    for ( j = 0; j <= max_resoln; j++ )    {
      fprintf( fp, "j = %d\n", j );
      for ( i = 0; i < L_bound[j].size; i++ )
	fprintf( fp, "%f\n", (*L)[j][i] );
    }
    fclose( fp );
  }
*/
}


/****************************************************************
*  Function: KSfilter_bound:
*  -------------------------
*  Computation of the size of K and S filters
*
*    filterename: filter name
*    K_bound: the size of K filters
*    S_bound: the size of L filters
*    max_resoln: number of decomposition
*
****************************************************************/

void KSfilter_bound(filtername,K_bound,S_bound,max_resoln)
     char *filtername;
     bound **K_bound, **S_bound;
     int max_resoln;
{
  int j;
  int iexp2();

  if(!(*K_bound = (bound *) malloc( (max_resoln+1) * sizeof(bound) )))
    error("Memory allocation failed for *K_bound in signal_back.c \n");
  if(!(*S_bound = (bound *) malloc( (max_resoln+1) * sizeof(bound) )))
    error("Memory allocation failed for *S_bound in filter.c \n");
  
  for ( j = 0; j <= max_resoln; j++ )    {
    if(strcmp(filtername,"Haar") == 0) {
      if(j == 0) {
	(*S_bound)[0].lb = -1;
	(*S_bound)[0].ub = 0;  
	(*S_bound)[0].size = (*S_bound)[0].ub - (*S_bound)[0].lb + 1;
	
	(*K_bound)[0].lb = -1;
	(*K_bound)[0].ub = 0; 
	(*K_bound)[0].size = (*K_bound)[0].ub - (*K_bound)[0].lb + 1;
      }
      else {
	(*S_bound)[j].lb = -iexp2(j-1);
	(*S_bound)[j].ub = iexp2(j-1); 
	(*S_bound)[j].size = (*S_bound)[j].ub - (*S_bound)[j].lb + 1;

	(*K_bound)[j].lb = -iexp2(j-1);
	(*K_bound)[j].ub = iexp2(j-1);
	(*K_bound)[j].size = (*K_bound)[j].ub - (*K_bound)[j].lb + 1;	
      }
    }
    else {
      if(j == 0) {
	(*S_bound)[0].lb = -2;
	(*S_bound)[0].ub = 1;  
	(*S_bound)[0].size = (*S_bound)[0].ub - (*S_bound)[0].lb + 1;
      
	(*K_bound)[0].lb = -3;
	(*K_bound)[0].ub = 2; 
	(*K_bound)[0].size = (*K_bound)[0].ub - (*K_bound)[0].lb + 1;
      }
      else {
	(*S_bound)[j].lb = -3 * iexp2(j-1);
	(*S_bound)[j].ub = 3 * iexp2(j-1); /* changed 2/1/94 j to j-1 */
	(*S_bound)[j].size = (*S_bound)[j].ub - (*S_bound)[j].lb + 1;

	(*K_bound)[j].lb = -5 * iexp2(j-1);
	(*K_bound)[j].ub = 5 * iexp2(j-1);/* changed 2/1/94 j to j-1 */
	(*K_bound)[j].size = (*K_bound)[j].ub - (*K_bound)[j].lb + 1;	
      }
    }
  }

/*
  {
    int j;
    for ( j = 0; j <= max_resoln; j++ )
      fprint("S_bound[%d] = [%d, %d], size = %d\n", j,
	      (*S_bound)[j].lb, (*S_bound)[j].ub, (*S_bound)[j].size );
    fprint( fp, "\n" );
    for ( j = 0; j <= max_resoln; j++ )
      fprintf( fp, "K_bound[%d] = [%d, %d], size = %d\n", j,
	      (*K_bound)[j].lb, (*K_bound)[j].ub, (*K_bound)[j].size );
    fprint( fp, "\n" );
  }
*/
}

/****************************************************************
*  Function: Lfilter_bound:
*  -------------------------
*  Computation of the size of L filters
*
*    filterename: filter name
*    L_bound: the size of L filters
*    max_resoln: number of decomposition
*
****************************************************************/

void Lfilter_bound(filtername,L_bound,max_resoln)
     char *filtername;
     bound **L_bound;
     int max_resoln;
{
  int j;
  int iexp2();

  if(!(*L_bound = (bound *) malloc( (max_resoln+1) * sizeof(bound) )))
    error("Memory allocation failed for *L_bound in filter.c \n");
  
  for ( j = 0; j <= max_resoln; j++ )    {
    if(strcmp(filtername,"Haar") == 0) {
      if(j == 0) {
	(*L_bound)[0].lb = -1;
	(*L_bound)[0].ub = 1;  
	(*L_bound)[0].size = (*L_bound)[0].ub-(*L_bound)[0].lb + 1;
      }
      else {
	(*L_bound)[j].lb = -iexp2(j);
	(*L_bound)[j].ub = iexp2(j); 
	(*L_bound)[j].size = (*L_bound)[j].ub-(*L_bound)[j].lb + 1;
      }
    }
    else {
      if(j == 0) {
	(*L_bound)[0].lb = -3;
	(*L_bound)[0].ub = 3;  
	(*L_bound)[0].size = (*L_bound)[0].ub-(*L_bound)[0].lb + 1;
      }
      else {
	(*L_bound)[j].lb = -3 * iexp2(j);
	(*L_bound)[j].ub = 3 * iexp2(j); 
	(*L_bound)[j].size = (*L_bound)[j].ub-(*L_bound)[j].lb + 1;
      }
    }
  }
/*
  {
    int j;
    for ( j = 0; j <= max_resoln; j++ )
      fprint("L_bound[%d] = [%d, %d], size = %d\n", j,
	      (*L_bound)[j].lb,(*L_bound)[j].ub,(*L_bound)[j].size);
    fprint( fp, "\n" );
  }
*/
}

/****************************************************************
*  Function: PsiPhifilter_bound:
*  -----------------------------
*  Computation of the size of Psi and Phi filters
*
*    psi: the size of psi filters
*    phi: the size of phi filters
*    H_bound: the size of H filters
*    G_bound: the size of G filters   
*    max_resoln: number of decomposition
*
****************************************************************/

void PsiPhifilter_bound(psi,phi,H_bound,G_bound,max_resoln)
     bound **psi, **phi;
     bound *G_bound;
     bound *H_bound;
     int max_resoln;
{
  int j;

  if(!(*psi = (bound *) malloc( (max_resoln+1) * sizeof(bound) )))
    error("Memory allocation failed for *psi in K_compute.c \n");

  if(!(*phi = (bound *) malloc( (max_resoln+1) * sizeof(bound) )))
    error("Memory allocation failed for *phi in K_compute.c \n");
  
  (*phi)[0].lb = (*phi)[0].ub = 0;
  (*phi)[0].size = 1;

  
  for ( j = 1; j <= max_resoln; j++ )  {
    if( j == 1 ) {
      (*psi)[j].lb =  G_bound[j].lb;
      (*psi)[j].ub =  G_bound[j].ub;
      (*phi)[j].lb =  H_bound[j].lb;
      (*phi)[j].ub =  H_bound[j].ub;
    }
    else {
      (*psi)[j].lb = (*psi)[j-1].lb +  G_bound[j].lb;
      (*psi)[j].ub = (*psi)[j-1].ub +  G_bound[j].ub;
      (*phi)[j].lb = (*phi)[j-1].lb +  H_bound[j].lb;
      (*phi)[j].ub = (*phi)[j-1].ub +  H_bound[j].ub;
    }
    (*psi)[j].size = (*psi)[j].ub - (*psi)[j].lb + 1;
    (*phi)[j].size = (*phi)[j].ub - (*phi)[j].lb + 1;

/*    printf("<%d> [%d,%d] %d  ====  [%d,%d] %d\n", j,
	(*psi)[j].lb,(*psi)[j].ub,(*psi)[j].size,
	(*phi)[j].lb,(*phi)[j].ub,(*phi)[j].size ); */
  
  }
}

/****************************************************/
/*                                                  */
/*           H0    H1    H2   .....     HJ-1        */
/*      o-----o-----o-----o---------o-----o S[J]    */
/*       \     \     \     \ .....   \              */
/*     G0 \  G1 \  G2 \  G3 \    GJ-1 \             */
/*	     o     o     o     o       o            */
/*        W[1]  W[2]  W[3]  W[4]      W[J]          */  
/*                                                  */
/*                                                  */
/****************************************************/


/****************************************************************
*  Function: signal_W_S:
*  ---------------------
*  Computation of the W and S filters
*
*    W: W filters
*    S: S filters
*    max_resoln: number of decomposition
*    np: signal size
*
****************************************************************/

void signal_W_S(W,S,max_resoln,np)
     float ***W, ***S;
     int max_resoln, np;
{
  int j, m, n, t;
  char filename1[STRING_SIZE],filename2[STRING_SIZE],filtername[STRING_SIZE];
  bound *H_bound,*G_bound;
  float **H_filter,**G_filter;
  float **H;
  float **G;
  float *prev,*curr,*temp,*normalize_factor;
  
  if(!(H = (float **) malloc( max_resoln * sizeof(float *) )))
    error("Memory allocation failed for H in oneD_filter.c \n");
  if(!(G = (float **) malloc( max_resoln * sizeof(float *) )))
    error("Memory allocation failed for G in oneD_filter.c \n");
  if(!(prev = (float *) malloc( np * sizeof(float) )))
    error("Memory allocation failed for prev in oneD_filter.c \n");
  if(!(curr = (float *) malloc( np * sizeof(float) )))
    error("Memory allocation failed for curr in oneD_filter.c \n");
  if(!(temp = (float *) malloc( np * sizeof(float) )))
    error("Memory allocation failed for temp in oneD_filter.c \n");

  filename_given(filtername,"Gaussian1");
  HGfilter_bound(filtername,&H_bound,&G_bound,max_resoln );  
  Hfilter_compute(filtername,&H_filter, H_bound, max_resoln );
  Gfilter_compute(filtername,&G_filter, G_bound, max_resoln);
/*  printf("Using Gaussian1 filter \n"); */

  for ( j = 0; j < max_resoln; j++ )   {
    if(!(H[j] = (float *) malloc( np * sizeof(float) )))
      error("Memory allocation failed for H[] in oneD_filter.c \n");
    if(!(G[j] = (float *) malloc( np * sizeof(float) )))
      error("Memory allocation failed for G[] in oneD_filter.c \n");
    
    for ( m = 0; m < np; m++ ) 
      H[j][m] = G[j][m] = 0.0;


    for ( m = H_bound[j].lb , t = 0; t < H_bound[j].size; t++, m++ )
      H[j][(m+np)%np] = H_filter[j][t];
    for ( m = G_bound[j].lb, t = 0; t < G_bound[j].size; t++, m++ )
      G[j][(m+np)%np] = G_filter[j][t]; 
/*    
    filename_given(filename1,"G");
    filename_given(filename2,"H");
    filename_inc(filename1,j);
    output_signal(G[j],np,filename1);
    filename_inc(filename2,j);
    output_signal(H[j],np,filename2);   
*/
  }


  if(!(*W = (float **) malloc( (max_resoln+1) * sizeof(float *))))
    error("Memory allocation failed for *W in oneD_filter.c \n");
  if(!(*S = (float **) malloc( (max_resoln+1) * sizeof(float *) )))
    error("Memory allocation failed for *S in oneD_filter.c \n");

  for ( j = 1; j <= max_resoln; j++ )  {
    if(!((*W)[j] = (float *) malloc( np * sizeof(float))))
      error("Memory allocation failed for (*W)[] in oneD_filter.c \n");
    if(!((*S)[j] = (float *) malloc( np * sizeof(float) )))
      error("Memory allocation failed for (*S)[] in oneD_filter.c \n");

    if ( j == 1 )    {
      for ( m = 0; m < np; m++ )      {
	(*W)[j][m] = G[0][m];
	(*S)[j][m] = H[0][m];
      }
    }
    else if ( j == 2 )    {
      compute_convolution( (*W)[j], G[j-1], H[j-2], np );
      compute_convolution( (*S)[j], H[j-1], H[j-2], np );
      for ( m = 0; m < np; m++ )
	prev[m] = H[0][m];
    }
    else    {
      compute_convolution( curr, H[j-2], prev, np );
      compute_convolution( (*W)[j], G[j-1], curr, np );
      compute_convolution( (*S)[j], H[j-1], curr, np );


      if ( j < max_resoln )      {
	for ( m = 0; m < np; m++ )
	  prev[m] = curr[m];
      }
    }
/*
    filename_given(filename1,"W");
    filename_given(filename2,"S");
    filename_inc(filename1,j);
    output_signal((*W)[j],np,filename1);
    filename_inc(filename2,j);
    output_signal((*S)[j],np,filename2);
*/
  }


  free( H_bound );
  free( G_bound );
  free( prev );
  free( curr );
  for ( j = 0; j < max_resoln; j++ )
  {
    free( H_filter[j] );
    free( G_filter[j] );
    free( H[j] );
    free( G[j] );
  }
  free( H_filter );
  free( G_filter );
  free( H );
  free( G );

}

/****************************************************************
*  Function: signal_W_hat_S_hat:
*  -----------------------------
*  Computation of the Fourier transform of W and S filters
*
*    W_hat: Fourier transform of W filters
*    S_hat: Fourier transform of S filters
*    max_resoln: number of decomposition
*    np: signal size
*
****************************************************************/


/* 
   Note:
   W_hat[1] = G_hat[0]                     S_hat[1] = H_hat[0]          
   W_hat[2] = G_hat[1]*H_hat[0]            S_hat[2] = H_hat[1]*H_hat[0]
   W_hat[3] = G_hat[2]*H_hat[1]*H_hat[0]   S_hat[3] = H_hat[2]*H_hat[1]*H_hat[0]
    :
    :
   W_hat[J] = G_hat[J-1]*H_hat[J-2]*......*H[0]
*/

void signal_W_hat_S_hat(W_hat,S_hat,max_resoln,np)
     float ***W_hat, ***S_hat;
     int max_resoln; 
     int np;
{
  char filename1[STRING_SIZE],filename2[STRING_SIZE],filtername[STRING_SIZE];
  int two_np;
  float *prev, *curr, **H_hat, **G_hat;
  int j, m;

  two_np = 2 * np;   /* real and imaginary */
  if(!(prev = (float *) malloc( two_np * sizeof(float))))
     error("Memory allocation failed for prev in oneD_filter.c \n");
  if(!(curr = (float *) malloc( two_np * sizeof(float) )))
     error("Memory allocation failed for curr in oneD_filter.c \n");

  filename_given(filtername,"Gaussian1");
  HG_hat_compute(filtername,&H_hat,&G_hat,max_resoln,np);
/*  printf("computing W_hat & S_hat with Gaussian1 filter\n");  */

  if(!(*W_hat = (float **) malloc( (max_resoln+1) * sizeof(float) )))
     error("Memory allocation failed for *W_hat in oneD_filter.c \n");
  if(!(*S_hat = (float **) malloc( (max_resoln+1) * sizeof(float) )))
    error("Memory allocation failed for *S_hat in oneD_filter.c \n");

  if(!((*S_hat)[0] = (float *) malloc( two_np * sizeof(float) )))
     error("Memory allocation failed for *S_hat in oneD_filter.c \n");

  for ( m = 0; m < np; m++ )  {
    (*S_hat)[0][2*m] = 1.0;
    (*S_hat)[0][2*m+1] = 0.0;
  }
  for ( j = 1; j <= max_resoln; j++ )  {
    if(!((*W_hat)[j] = (float *) malloc( two_np * sizeof(float) )))
      error("Memory allocation failed for (*W_hat)[] in oneD_filter.c \n");
    if(!((*S_hat)[j] = (float *) malloc( two_np * sizeof(float) )))
      error("Memory allocation failed for (*S_hat)[] in oneD_filter.c \n");

    if ( j == 1 ) {    /* H_hat & G_hat:  j = 0 to (max_reaoln-1) */
      for ( m = 0; m < two_np; m++ )  {
	(*W_hat)[j][m] = G_hat[j-1][m];
	(*S_hat)[j][m] = H_hat[j-1][m];
      }
    }
    else if ( j == 2 )    {
      complex_product( (*W_hat)[j], G_hat[j-1], H_hat[j-2], np );
      complex_product( (*S_hat)[j], H_hat[j-1], H_hat[j-2], np );
      for ( m = 0; m < two_np; m++ )
	prev[m] = H_hat[0][m];
    }
    else    {
      complex_product( curr, H_hat[j-2], prev, np );      
      complex_product( (*W_hat)[j], G_hat[j-1], curr, np );
      complex_product( (*S_hat)[j], H_hat[j-1], curr, np );
      for ( m = 0; m < two_np; m++ )
	prev[m] = curr[m];      
    }

/*
    filename_given(filename1,"WhatR");
    filename_given(filename2,"WhatI");
    filename_inc(filename1,j);
    filename_inc(filename2,j);
    output_complex((*W_hat)[j],np,filename1,filename2);

    filename_given(filename1,"ShatR");
    filename_given(filename2,"ShatI");

    filename_inc(filename1,j);
    filename_inc(filename2,j);
    output_complex((*S_hat)[j],np,filename1,filename2);
*/
  }

  free( prev );
  free( curr );
  for ( j = 0; j < max_resoln; j++ )  {
    free( H_hat[j] );
    free( G_hat[j] );
  }
  free( H_hat );
  free( G_hat );
}





