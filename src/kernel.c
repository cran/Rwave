/***************************************************************
*    $Log: kernel.c,v $                                        *
****************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
***************************************************************/





#include "kernel.h"





/******************************************************************
*  Function: integrand
*  ---------
*     Evaluates the integrand for Romberg integration.
*
*    b: center position of wavelets and derivatives.
*    x,y: Integration variables.
*    p2: second derivative of ridge function.
*    nodes, phi_nodes: position and scale of the nodes of the ridge.
*              (warning: this is the real scale, and not its log).
*    nb_nodes: number of nodes of the ridge.
*    w0: central frequency of Morlet's wavelet.
*******************************************************************/

fcomplex integrand(double b,int x,int y,double *p2,double *nodes,
		   double *phi_nodes,int nb_nodes,double w0)
{
  fcomplex ctmp,psi_x,psi_y,psiprime_x,psiprime_y;
  double xpos,ypos,xpos0,ypos0,phi_b2,tmp,phi,phi_b,phip_b,w02;
  fcomplex integ;

  integ.r =0.0; integ.i = 0.0;
  w02 = w0 * w0;


  /* Spline interpolation of the ridge; computation of its derivative
     ---------------------------------------------------------------*/
  splint2(nodes,phi_nodes,p2,nb_nodes,b,&phi_b,&phip_b);


  /* Evaluation of the integrand
     --------------------------*/
  xpos0 = x-b;
  xpos = xpos0/phi_b;
  ypos0 = y-b;
  ypos = ypos0/phi_b;

  phi_b2 = phi_b * phi_b;

  psi_x = psi(xpos,w0);
  psiprime_x = psiprime(xpos,w0);
  psi_y = psi(ypos,w0);
  psiprime_y = psiprime(ypos,w0);

  /* psi_psi term */
  integ = Cmul(psi_y,Conjg(psi_x)) ;
  tmp= (phip_b*phip_b - w02);
  ctmp=Complex(tmp,0);
  integ = Cmul(integ,ctmp);

  /* psiprime_psiprime term */
  ctmp=Complex(1.0 + xpos + ypos + xpos*ypos,0);
  ctmp = Cmul(ctmp,Conjg(psiprime_x));
  ctmp = Cmul(ctmp,psiprime_y);
  integ = Cadd(integ,ctmp);

  /* psi_psiprime terms */
  tmp = 1. + ypos;
  ctmp = Complex(tmp*phip_b,0.);
  ctmp = Cmul(ctmp,Conjg(psi_x));
  ctmp = Cmul(ctmp,psiprime_y);
  integ = Cadd(integ,ctmp);

  tmp = 1. + xpos;
  ctmp = Complex(tmp*phip_b,0.);
  ctmp = Cmul(ctmp,Conjg(psiprime_x));
  ctmp = Cmul(ctmp,psi_y);
  integ = Cadd(integ,ctmp);

  /* normalisation */
  tmp =phi_b2*phi_b2;
  ctmp = Complex(1/tmp,0);
  integ = Cmul(integ,ctmp);

  return (integ);
}




/******************************************************************
*  Function: rintegrand
*  ---------
*     Evaluates the integrand for Romberg integration, in the
*     case of a real valued wavelet kernel.
*
*    b: center position of wavelets and derivatives.
*    x,y: Integration variables.
*    p2: second derivative of ridge function.
*    nodes, phi_nodes: position and scale of the nodes of the ridge.
*              (warning: this is the real scale, and not its log).
*    nb_nodes: number of nodes of the ridge.
*    w0: central frequency of Morlet's wavelet.
*******************************************************************/

double rintegrand(double b,int x,int y,double *p2,double *nodes,
		   double *phi_nodes,int nb_nodes,double w0)
{
  fcomplex psi_x,psi_y,psiprime_x,psiprime_y;
  double xpos,ypos,xpos0,ypos0,phi_b2,tmp,phi,phi_b,phip_b,w02;
  double integ;

  integ =0.0;
  w02 = w0 * w0;


  /* Spline interpolation of the ridge; computation of its derivative
     ---------------------------------------------------------------*/
  splint2(nodes,phi_nodes,p2,nb_nodes,b,&phi_b,&phip_b);


  /* Evaluation of the integrand
     --------------------------- */
  xpos0 = x-b;
  xpos = xpos0/phi_b;
  ypos0 = y-b;
  ypos = ypos0/phi_b;

  phi_b2 = phi_b * phi_b;

  psi_x = psi(xpos,w0);
  psiprime_x = psiprime(xpos,w0);
  psi_y = psi(ypos,w0);
  psiprime_y = psiprime(ypos,w0);

  /* psi_psi term */
  integ = (psi_y).r * (psi_x).r + (psi_y).i * (psi_x).i ;
  tmp= (phip_b*phip_b - w02);
  integ *= tmp;

  /* psiprime_psiprime term */
  tmp=1.0 + xpos + ypos + xpos*ypos;
  tmp *= ((psiprime_y).r*(psiprime_x).r + (psiprime_y).i*(psiprime_x).i);
  integ += tmp;

  /* psi_psiprime terms */
  tmp = 1. + ypos;
  tmp *= phip_b;
  tmp *= ((psiprime_y).r*(psi_x).r + (psiprime_y).i*(psi_x).i);
  integ += tmp;

  tmp = 1. + xpos;
  tmp *= phip_b;
  tmp *= ((psi_y).r*(psiprime_x).r + (psi_y).i*(psiprime_x).i);
  integ += tmp;

  /* normalisation */
  tmp =phi_b2*phi_b2;
  integ /= tmp;

  return (integ);
}





/******************************************************************
*  Function:  maxvalue
*  ---------
*      Maximal element of an array
*
*     vector: array whose maximal element is seeked.
*     length: length of the array.
*******************************************************************/

double maxvalue(double *vector, int length)
{
  double maxi;
  int i;
  
  maxi=0;
  for(i=0; i< length; i++){
    maxi = MAX(maxi,*vector);
    vector++;
  }
  return maxi;
}




/*********************************************************
*  Function:  hermite_sym
*  ---------
*      Complete a matrix by Hermite symmetry
*
*   ker: matrix to be filled in.
*   lng: number of rows (and columns) of the matrix.
**********************************************************/

void hermite_sym(fcomplex *ker,int lng)
{
  int i,j;

  for (i=0;i<lng;i++){
    for (j=lng-1;j>i;j--){
/*	ker[i*lng +j] = Conjg(ker[j*lng +i]);*/
      (ker[i*lng +j]).r = (ker[j*lng +i]).r;
      (ker[i*lng +j]).i = -(ker[j*lng +i]).i;
      }
  }
  return;
}



/***********************************************************
*  Function: kernel
*  ---------
*     Computation of the kernel
*
*   ker_r, ker_i: real and imaginary parts of the kernel.
*   px_min, px_max: limiting values for the integration.
*   px_inc: distance between 2 consecutive x,y values.
*   plng:length
*   nodes, phi_nodes: position and (true) scale of the 
*       samples of the ridge.
*   pnb_nodes: number of nodes of the sampled ridge.
*   pw0: central frequency of Morlet'x wavelet.
*   pb_start, pb_end: integration bounds.
************************************************************/


void kernel(double *ker_r, double *ker_i,int *px_min,int *px_max,
	    int *px_inc, int *plng, double *nodes,double *phi_nodes,
	    int *pnb_nodes,double *pw0,double *pb_start,double *pb_end)
{
  double *p2, b_start=*pb_start, b_end=*pb_end, w0=*pw0;
  double phimax, b_lo,b_hi;
  int x,y,yy;
  int x_min=*px_min,x_max=*px_max,x_inc=*px_inc;
  int lng=*plng,nb_nodes=*pnb_nodes;
  int i=0,up_bound,gamma_min,lng2;
  fcomplex *ker,*p_tmp;
  fcomplex tmp;


  p2 = (double *)calloc(nb_nodes,sizeof(double));
  ker = (fcomplex *)calloc(lng*lng,sizeof(fcomplex));
  p_tmp=ker; /* mark the first element of ker */

  phimax = maxvalue(phi_nodes,nb_nodes);
  up_bound = (int)(phimax * sqrt(-2.0 * log(EPS))+1);
  lng2=lng*lng;

  /* Compute second derivative of the ridge for spline interpolation
     ---------------------------------------------------------------*/
  spline(nodes-1,phi_nodes-1,nb_nodes,(double)0,(double)0,p2-1);


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
      *ker = qrombmod(x,y,p2-1,nodes,phi_nodes,nb_nodes,w0,b_lo,b_hi); 
      ker++; i++;
    }
    ker -= (i - lng) ;
  }
  ker = p_tmp;


  /* Finish to fill in the kernel by Hermite symmetry
     -----------------------------------------------*/
/*  printf("Now going to hermite_sym\n");*/
  hermite_sym(ker,lng);

/*
  i=0;
  for(x=x_min;x<=x_max;x+=x_inc){
    for(y=x_min;y<=x_max;y+=x_inc){
      *ker_r = ker->r; ker_r++;
      *ker_i = ker->i; ker_i++;
      ker++; i++;
    }
  }

  ker -= i;
  ker_r -= i; ker_i -= i;
*/

/*  printf("Now loading in ker_r and ker_i\n");*/
  for(i=0;i<lng2;i++){
    *ker_r = ker->r; ker_r++;
    *ker_i = ker->i; ker_i++;
    ker++;
  }

  ker -= lng2;
  ker_r -= lng2; ker_i -= lng2;

  free((char *)ker);
  free(p2);
}



/***********************************************************
*  Function: rkernel
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


void rkernel(double *ker,int *px_min,int *px_max,int *px_inc,
	    int *plng, double *nodes,double *phi_nodes,int *pnb_nodes,
	    double *pw0,double *pb_start,double *pb_end)
{
  double *p2, b_start=*pb_start, b_end=*pb_end, w0=*pw0;
  double phimax, b_lo,b_hi;
  int x,y,yy;
  int x_min=*px_min,x_max=*px_max,x_inc=*px_inc,lng=*plng,nb_nodes=*pnb_nodes;
  int i=0,j,up_bound,gamma_min,lng2;
  double *p_tmp;
  fcomplex tmp;


  p2 = (double *)calloc(nb_nodes,sizeof(double));
  p_tmp=ker; /* mark the first element of ker */

  phimax = maxvalue(phi_nodes,nb_nodes);
  up_bound = (int)(phimax * sqrt(-2.0 * log(EPS))+1);
  lng2=lng*lng;

  /* Compute second derivative of the ridge for spline interpolation
     ---------------------------------------------------------------*/
  spline(nodes-1,phi_nodes-1,nb_nodes,(double)0,(double)0,p2-1);


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
      *ker = rqrombmod(x,y,p2-1,nodes,phi_nodes,nb_nodes,w0,b_lo,b_hi); 
      ker++; i++;
    }
    ker -= (i - lng) ;
  }
  ker = p_tmp;


  /* Finish to fill in the kernel by Hermite symmetry
     -----------------------------------------------*/
/*  printf("Now going to hermite_sym\n");*/
  ghermite_sym(ker,lng);


  /*  for(i=0;i<lng;i++)
    for(j=0;j<lng;j++)
      printf("%f; ",ker[i+lng*j]); */


  free(p2);
}



/***********************************************************
*  Function: fastkernel
*  ---------
*     Computation of the kernel; the integral is
*      approximated by a Riemann sum
*
*   ker_r, ker_i: real and imaginary parts of the kernel.
*   x_min, x_max: limiting values for the integration.
*   x_inc: distance between 2 consecutive x,y values.
*   nodes, phi_nodes: position and scale of the samples
*       of the ridge.
*   nb_nodes: number of nodes of the sampled ridge.
*   b_start, b_end: integration bounds.
************************************************************/


void fastkernel(double *ker_r, double *ker_i,int *px_min,int *px_max,
	    int *px_inc, int *plng, double *nodes,double *phi_nodes,
	    int *pnb_nodes,double *pw0,double *pb_start,double *pb_end)
{
  double *p2, b_start=*pb_start, b_end=*pb_end, w0=*pw0;
  double phimax, b_lo,b_hi;
  int x,y,yy,b;
  int x_min=*px_min,x_max=*px_max,x_inc=*px_inc,lng=*plng,nb_nodes=*pnb_nodes;
  int i=0,up_bound,gamma_min,lng2;
  fcomplex *ker,*p_tmp;
  fcomplex tmp;


  p2 = (double *)calloc(nb_nodes,sizeof(double));
  ker = (fcomplex *)calloc(lng*lng,sizeof(fcomplex));
  p_tmp=ker; /* mark the first element of ker */

  phimax = maxvalue(phi_nodes,nb_nodes);
  up_bound = (int)(phimax * sqrt(-2.0 * log(EPS))+1);
  lng2=lng*lng;

  /* Compute second derivative of the ridge for spline interpolation
     ---------------------------------------------------------------*/
  spline(nodes-1,phi_nodes-1,nb_nodes,(double)0,(double)0,p2-1);


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
      for(b = (int)b_lo; b<= (int)b_hi;b++){
	*ker = Cadd(*ker,integrand(b,x,y,p2-1,nodes,phi_nodes,nb_nodes,w0));
      }
      ker++; i++;
    }
    ker -= (i - lng) ;
  }
  ker = p_tmp;


  /* Finish to fill in the kernel by Hermite symmetry
     -----------------------------------------------*/
/*  printf("Now going to hermite_sym\n");*/
  hermite_sym(ker,lng);

/*
  i=0;
  for(x=x_min;x<=x_max;x+=x_inc){
    for(y=x_min;y<=x_max;y+=x_inc){
      *ker_r = ker->r; ker_r++;
      *ker_i = ker->i; ker_i++;
      ker++; i++;
    }
  }

  ker -= i;
  ker_r -= i; ker_i -= i;
*/

/*  printf("Now loading in ker_r and ker_i\n");*/
  for(i=0;i<lng2;i++){
    *ker_r = ker->r; ker_r++;
    *ker_i = ker->i; ker_i++;
    ker++;
  }

  ker -= lng2;
  ker_r -= lng2; ker_i -= lng2;

  free((char *)ker);
  free(p2);
}


/*********************************************************
*  Function: psi
*  ---------
*    Computation of the wavelet
*
*   x: location of the wavelet
*   w0: center frequency of Morlet's wavelet
**********************************************************/

fcomplex psi(double x,double w0)
{
  double u;
  
  u=exp(-x*x/2.);
  return(Cmul(Complex(u,0.0),Complex(cos(w0*x),sin(w0*x)))); 

}

/*********************************************************
*  Function: psiprime
*  ---------
*  Computation of the wavelet derivative
*
*   x: location of the wavelet
*   w0: center frequency of Morlet's wavelet
**********************************************************/

fcomplex psiprime(double x,double w0)
{
  double u;
  fcomplex v;
  
  u=exp(-x*x/2.);
  v = Complex(-x,w0);
  v = Cmul(v,Complex(u,0.0));
  return(Cmul(v,Complex(cos(w0*x),sin(w0*x))));
}




























