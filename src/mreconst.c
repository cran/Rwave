#include <stdlib.h>


/******************************************************************
*              (c) Copyright  1997                                *
*                         by                                      *
*  Author: Rene Carmona, Bruno Torresani, Wen L. Hwang, A. Wang   *
*                 Princeton University                            *
*                 All right reserved                              *
******************************************************************/


#include "Swave.h"
#include "dyadic.h"

/****************************************************************
*  Function: signal_penalty_function:
*  ----------------------------------
*  Reconstruction of signal which preserves the locations and 
*  values of extrema
*
*    f: reconstructed signal
*    lambda: Lagragian multiplier
*    W_tilda: dual wavelet 
*    ext: extrema representation
*    num_of_extrema: number of extrema
*    np: signal size
*
****************************************************************/

void signal_penalty_function(double *f, double *lambda,
  double **W_tilda, image_ext *ext, int num_of_extrema, int np)
{
  int s, t;
  
  for ( t = 0; t < np; t++ )
  {
    for ( s = 0, f[t] = 0.0; s < num_of_extrema; s++ )
      f[t] += lambda[s] * W_tilda[ext[s].resoln][(ext[s].x-t+np)%np];
  }
  return;
}


/****************************************************************
*  Function: signal_position:
*  --------------------------
*  Position matrix of the extrema of a signal
*
*    filtername: filter name
*    lambda: Lagragian multiplier
*    ext: extrema representation
*    Wtilda:
*    W:
*    num_of_extrema: number of extrema
*    max_resoln: number of decomposition
*    np: signal size
*
****************************************************************/

void signal_position(char *filtername, double **lambda,
  image_ext *ext, double **Wtilda, double **W, int num_of_extrema,
  int max_resoln, int np)
{
  int s, t, r, i, j, diff;
  bound *psi, *phi;
  bound *H_bound, *G_bound;
  char filename[STRING_SIZE];
  double **position_matrix, *w, **v;
  double sum;
  int p, x;
  double *b;

  double d;
  int *indx;

  if(!(indx = (int *) R_alloc(num_of_extrema , sizeof(int))))
    error("Memory allocation failed for indx in signal_position.c \n");

  HGfilter_bound(filtername,&H_bound,&G_bound,max_resoln); 
  PsiPhifilter_bound( &psi, &phi, H_bound, G_bound, max_resoln );

  if(!(position_matrix = (double **) R_alloc(num_of_extrema , sizeof(double *) )))
    error("Memory allocation failed for position matrix in image_lambda \n");
  for ( r = 0; r < num_of_extrema; r++ )
    if(!(position_matrix[r] = (double *) R_alloc((num_of_extrema) , sizeof(double) )))
      error("Memory allocation failed for position_matrix[] in image_lambda \n");

  for ( r = 0; r < num_of_extrema; r++ ) {
    i = ext[r].resoln;
    for ( s = 0; s < num_of_extrema; s++ )    {
      j = ext[s].resoln;
      diff = ext[s].x - ext[r].x;
      sum = 0.0;

      for (x = psi[i].lb; x <= psi[i].ub; x++ ) {
	p = (diff+x+2 * np)%np;
	sum +=  W[i][(x+np)%np] * Wtilda[j][p];
      }
      position_matrix[r][s] = sum;
    }
  }

/*   init_filename(filename); */
/*  filename_given(filename,"signal_POS"); */
/*  printf("num_of_extrema = %d \n",num_of_extrema); */
/*  output_array(position_matrix,num_of_extrema,num_of_extrema,filename); */


  /**************************************************/
  /* solve lambda from position_matrix (lambda) = b */
  /**************************************************/

  if(!(*lambda = (double *) R_alloc(num_of_extrema , sizeof(double) )))
    error("Memory allocation failed for lambda in image_position.c \n");

  if(!(b = (double *) R_alloc( num_of_extrema , sizeof(double))))
    error("Memory allocation failed for b in image_position.c \n");  


  for ( i = 0; i < num_of_extrema; i++ )  {
    b[i] = ext[i].W1f;
  }
 

/*
  for ( i = 0; i < num_of_extrema; i++ )  {
    (*lambda)[i] = ext[i].W1f;
  }
  ludcmp(position_matrix-1,num_of_extrema,indx-1,&d);
  printf("here");
  lubksb(position_matrix-1,num_of_extrema,indx-1,(*lambda)-1);
  printf("There");
*/
  svdecomp_solve(position_matrix,b,*lambda,num_of_extrema,
		 num_of_extrema,&w,&v);
}


/****************************************************************
*  Function: extrema reconstruction:
*  ---------------------------------
*  Called by Splus. 
*  Reconstruction of a signal which preserves the location and value
*  at the extrema representation
*
*    filtername: filter name
*    f: reconstructed signal
*    extrema: extrema representation
*    max_resoln_ptr: pointer to number of decomposition
*    np_pr: pointer to signal size
*    preadfile: if set to TRUE, the kernel is read from files that 
*    are precomputed and stored in the local directory. OW, Compute
*    the kernel, this usually takes quite some time if the size of 
*    signal is large.
*    
****************************************************************/

void extrema_reconst(char *filtername, double *f, double *extrema,
  int *max_resoln_ptr, int *np_ptr, int *preadfile)
{
  double **W, **S;
  double **W_tilda;
  double **K;
  int max_resoln = *max_resoln_ptr;
  int np = *np_ptr;
  int readfileflag = *preadfile;
  image_ext *ext;                  
  double *lambda;
  int num_of_extrema, j;
  
  signal_W_S(&W, &S, max_resoln, np); 
  
  if(!readfileflag) {
    signal_K_compute(&K, W, max_resoln, np); 
   // printf("K "); 
    signal_W_tilda(&W_tilda, W, K, max_resoln,np);    
    //printf("W ");
  }
  else {
    signal_W_tilda_input(&W_tilda, max_resoln, np);   
  }
  
  extrema_input(extrema,max_resoln,np,&ext,&num_of_extrema); 
  signal_position(filtername,&lambda,ext,W_tilda,W,num_of_extrema,
		  max_resoln,np);   
  
  signal_penalty_function(f,lambda,W_tilda,ext,num_of_extrema,np);    

  //free( lambda );
  //free( ext );
  /*
    for ( j = 0; j <= max_resoln; j++ )  {
    free( W[j] );
    free S[j] );
    free( W_tilda[j] );
    }
    free( W );
    free( S );
    free( W_tilda );
    return;
  */
}



















