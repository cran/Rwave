#include <stdlib.h>
#include "Swave.h"


/***************************************************************
*    $log: gkernel.c,v $                                       *
****************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
***************************************************************/





#include "kernel.h"



/*********************************************************
*  Function:  ghermite_sym
*  ---------
*      Complete a matrix by symmetry
*
*   ker: matrix to be filled in.
*   lng: number of rows (and columns) of the matrix.
**********************************************************/

void ghermite_sym(double *ker,int lng)
{
  int i,j;

  for (i=0;i<lng;i++){
    for (j=lng-1;j>i;j--){
      ker[i*lng +j] = ker[j*lng +i];
      }
  }
  return;
}



/******************************************************************
*  Function: gintegrand
*  ---------
*     Evaluates the integrand for Romberg integration.
*      (case of the Gabor transform)
*
*    b: center position of wavelets and derivatives.
*    x,y: Integration variables.
*    nodes, phi_nodes: position and scale of the nodes of the ridge.
*    nb_nodes: number of nodes of the ridge.
*******************************************************************/

double gintegrand(double b,int x,int y,double *p2,double *nodes,
		   double *phi_nodes,int nb_nodes,double scale)

{
  double xpos,ypos,tmp,tmp2,phi,phi_b,phip_b, u;
  double g_x,g_y,gprime_x,gprime_y;
  double integ=0;


  /* Spline interpolation of the ridge; computation of its derivative
     ---------------------------------------------------------------*/
  splint2(nodes,phi_nodes,p2,nb_nodes,b,&phi_b,&phip_b);


  /* Evaluation of the integrand
     --------------------------*/
  xpos = (x-b);
  g_x = gfunc(xpos,scale);
  gprime_x = gprime(xpos,scale);
  ypos = (y-b);
  g_y = gfunc(ypos,scale);
  gprime_y = gprime(ypos,scale);

  u = phi_b*(x-y);

  /* First part */

  tmp= (phip_b*phip_b*xpos*ypos - phi_b*phip_b*(xpos+ypos));
  tmp *= g_x*g_y;
  integ = tmp;

  tmp=gprime_x*gprime_y;
  integ += tmp;

  integ *= cos(u);

  /* Second part */

  tmp2 = phip_b*xpos -phi_b;
  tmp2 *= (gprime_y*g_x);

  tmp = phip_b*ypos -phi_b;
  tmp *= (gprime_x*g_y);
  tmp2 -= tmp;

  tmp2 *= sin(u);

  integ += tmp2;

  return integ;
}





/***********************************************************
*  Function: gkernel
*  ---------
*     Computation of the kernel
*
*   ker_r, ker_i: real and imaginary parts of the kernel.
*   x_min, x_max: limiting values for the integration.
*   x_inc: distance between 2 consecutive x,y values.
*   nodes, phi_nodes: position and scale of the samples
*       of the ridge.
*   nb_nodes: number of nodes of the sampled ridge.
*   b_start, b_end: integration bounds.
************************************************************/

void gkernel(double *ker, int *px_min,int *px_max,
	    int *px_inc, int *plng, double *nodes,double *phi_nodes,
	    int *pnb_nodes,double *pscale,double *pb_start,double *pb_end)
{
  double *p2, b_start=*pb_start, b_end=*pb_end, scale=*pscale;
  double phimax, b_lo,b_hi;
  int x,y,yy;
  int x_min=*px_min,x_max=*px_max,x_inc=*px_inc,lng=*plng,nb_nodes=*pnb_nodes;
  int i=0,up_bound,gamma_min,lng2;
  double *p_tmp;
  double tmp;


  /* printf("xmin=%d, xmax=%d\n",x_min,x_max); */
  p2 = (double *)S_alloc(nb_nodes,sizeof(double));
  p_tmp=ker; /* mark the first element of ker */

  phimax = scale;
  up_bound = (int)(phimax * sqrt(-2.0 * log(EPS))+1);
  lng2=lng*lng;

  /* Compute second derivative of the ridge for spline interpolation
     ---------------------------------------------------------------*/
  spline(nodes-1,phi_nodes-1,nb_nodes,(double)0,(double)0,p2-1);
  /* printf("spline done\n"); */

  /* Integrate
     --------*/
  for(x=x_min;x<=x_max;x+=x_inc){
    /* fprintf(stderr,"x = %d;  ", x);
    fflush(stderr); */

    /* Evaluate the range of computation of the kernel */
    gamma_min = MAX(x_min,(x-2*up_bound) -(x-x_min-2*up_bound)%x_inc);
    ker += (gamma_min-x_min)/x_inc; i = (gamma_min-x_min)/x_inc;

    for(y = gamma_min; y <= x; y+=x_inc){
      
      /* Estimation of integration bounds */
      b_lo = MAX(MAX(x-2*up_bound,y-2*up_bound),b_start); 
      b_hi = MIN(MIN(x+2*up_bound,y+2*up_bound),b_end);
      
     /* Evaluation of the kernel */
      *ker = gqrombmod(x,y,p2-1,nodes,phi_nodes,nb_nodes,scale,b_lo,b_hi); 
      ker++; i++;
    }
    ker -= (i - lng) ;
  }
  ker = p_tmp;


  /* Finish to fill in the kernel by Hermite symmetry
     ------------------------------------------------ */
  ghermite_sym(ker,lng);

}





/***********************************************************
*  Function: fastgkernel
*  ---------
*     Computation of the kernel (Gabor case); 
*      the integral is approximated by a Riemann sum
*
*   ker_r, ker_i: real and imaginary parts of the kernel.
*   x_min, x_max: limiting values for the integration.
*   x_inc: distance between 2 consecutive x,y values.
*   nodes, phi_nodes: position and scale of the samples
*       of the ridge.
*   nb_nodes: number of nodes of the sampled ridge.
*   b_start, b_end: integration bounds.
************************************************************/


void fastgkernel(double *ker, int *px_min,int *px_max,
	    int *px_inc, int *plng, double *nodes,double *phi_nodes,
	    int *pnb_nodes,double *pscale,double *pb_start,double *pb_end)
{
  double *p2, b_start=*pb_start, b_end=*pb_end, scale=*pscale;
  double phimax, b_lo,b_hi;
  int x,y,yy,b;
  int x_min=*px_min,x_max=*px_max,x_inc=*px_inc,lng=*plng,nb_nodes=*pnb_nodes;
  int i=0,up_bound,gamma_min,lng2;
  double *p_tmp;


  p2 = (double *)S_alloc(nb_nodes,sizeof(double));
  p_tmp=ker; /* mark the first element of ker */

  phimax = scale;
  up_bound = (int)(phimax * sqrt(-2.0 * log(EPS))+1);
  lng2=lng*lng;

  /* Compute second derivative of the ridge for spline interpolation
     ---------------------------------------------------------------*/
  spline(nodes-1,phi_nodes-1,nb_nodes,(double)0,(double)0,p2-1);


  /* Integrate
     --------*/
  for(x=x_min;x<=x_max;x+=x_inc){
/*    fprintf(stderr,"x = %d;  ", x);
    fflush(stderr); */

    /* Evaluate the range of computation of the kernel */
    gamma_min = MAX(x_min,(x-2*up_bound) -(x-x_min-2*up_bound)%x_inc);
/*    fprintf(stderr,"gamma_min = %d; \n ", gamma_min);
    fflush(stderr); */
    ker += (gamma_min-x_min)/x_inc; i = (gamma_min-x_min)/x_inc;

    for(y = gamma_min; y <= x; y+=x_inc){
      
      /* Estimation of integration bounds */
      b_lo = MAX(MAX(x-2*up_bound,y-2*up_bound),b_start); 
      b_hi = MIN(MIN(x+2*up_bound,y+2*up_bound),b_end);

      /* Evaluation of the kernel */
      for(b = (int)b_lo; b<= (int)b_hi;b++){
	*ker += gintegrand(b,x,y,p2-1,nodes,phi_nodes,nb_nodes,scale);
	/* to be checked */
      }
      ker++; i++;
    }
    ker -= (i - lng) ;
  }
  ker = p_tmp;


  /* Finish to fill in the kernel by Hermite symmetry
     -----------------------------------------------*/
  ghermite_sym(ker,lng);

}




/*********************************************************
*  Computation of the gabor window
**********************************************************/

double gfunc(double x, double scale)
{
  double u,v;

  v = x/scale;
  u=exp(-v*v/(double)2)/scale;
  return u; 

}

/*********************************************************
*  Computation of the window derivative
**********************************************************/

double  gprime(double x,double scale)
{
  double u, v;
  
  v = x/scale;
  u= - v*exp(-v*v/(double)2.)/scale/scale;
  return u;
}


