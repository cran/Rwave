


/****************************************************************/
/*              (c) Copyright  1997
/*                         by                                   */
/*     Author: Rene Carmona, Andrea Wang, Wen-Liang Hwang       */
/*                 Princeton University                         */
/*                 All right reserved                           */
/****************************************************************/

/****************************************************************************
*	$Log: splsnake.c,v $
*****************************************************************************
*                                                                           *
*        This file  contains proprietary information                        *          
*                                                                           *
*****************************************************************************
*	 Cubic spline interpolation of the ridge of wavelet transform	    *
*	 of amplitude and frequency modulated signals			    *
*	 Modification of the routines spline.c and splint.c		    *
*	 (numerical Recipes)				                    *
*                                                                           *
*	y: input vector                                         	    *
*       yy: output (interpolated) vector                                    *
*       rate: subsampling rate                                              *
****************************************************************************/




#include "Swave.h"

/***************************************************************************
*   n: number of elements of the snake
*   rate: subsampling rate for the wavelet transform (b direction)
*   
***************************************************************************/

void splsnake(rate, x, y, n, yy)
     int rate,n;
     float *x, *y, *yy;
     
{
  int i,k, khi, klo, ilo, ihi;
  float p,qn,sig,un,*u,yp1,ypn,a,b,h;
  float *y2;
  
  u=(float *)calloc(n,sizeof(float));
  y2=(float *)calloc(n+1,sizeof(float));
  yp1 = ypn =0;
  
  if (yp1 > 0.99e30)
    y2[1]=u[1]=0.0;
  else {
    y2[1] = -0.5;
    u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
  }
  
  for (i=2;i<=n-1;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
  }
  y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
  for (k=n-1;k>=1;k--){
    y2[k]=y2[k]*y2[k+1]+u[k];
  }


/* Interpolation */
  
  ilo = (int)(x[1])*rate;
  ihi = (int)(x[n])*rate;

  for(i=ilo;i<ihi;i++){
    klo=1;
    khi=n;
    
    while (khi-klo > 1) {
      k=(khi+klo) >> 1;
      if (x[k]*rate > (float)i) khi=k;
      else klo=k;
    }
    h=(x[khi]-x[klo])*rate;
    if (h == 0.0) error("Impossible interpolation");
    a=(rate*x[khi]-i)/h;
    b=(i-x[klo]*rate)/h;
    yy[i]=a*y[klo]+b*y[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0;
  }
  free((float *)u);
  free((float *)y2);
  return;
}



