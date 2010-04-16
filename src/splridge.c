#include <stdlib.h>



/****************************************************************/
/*              (c) Copyright  1997
/*                         by                                   */
/*     Author: Rene Carmona, Andrea Wang, Wen-Liang Hwang       */
/*                 Princeton University                         */
/*                 All right reserved                           */
/****************************************************************/



/****************************************************************************
*	$Log: splridge.c,v $
 * Revision 1.1  1994/05/22  01:44:51  bruno
 * Initial revision
 *
 * Revision 1.1  1994/05/22  00:59:26  bruno
 * Initial revision
 *
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


/* #include <stdlib.h>
#include <stdio.h>
#include <math.h>
*/

#include "Swave.h"

void splridge(int rate, double *y, int n, double *yy)
     
{
  int i,k, x, khi, klo;
  double p,qn,sig,un,*u,yp1,ypn,a,b,h;
  double *y2;
  
  u=(double *)S_alloc(n-1,sizeof(double));
  y2=(double *)S_alloc(n,sizeof(double));
  yp1 = ypn =0;
  
  if (yp1 > 0.99e30)
    y2[0]=u[0]=0.0;
  else {
    y2[0] = -0.5;
    u[0]=(3.0/rate)*((y[1]-y[0])/rate-yp1);
  }
  
  for (i=1;i<=n-2;i++) {
    sig=2.0;
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/rate - (y[i]-y[i-1])/rate;
    u[i]=(6.0*u[i]/rate/2.0-sig*u[i-1])/p;
  }
  if (ypn > 0.99e30)
    qn=un=0.0;
  else {
    qn=0.5;
    un=(3.0/rate)*(ypn-(y[n-1]-y[n-2])/rate);
  }
  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
  for (k=n-2;k>=0;k--)
    y2[k]=y2[k]*y2[k+1]+u[k];
  
  
  for(x=0;x<n*rate;x++){
    klo=1;
    khi=n;
    
    while (khi-klo > 1) {
      k=(khi+klo) >> 1;
      if (k*rate > x) khi=k;
      else klo=k;
    }
    h=(khi-klo)*rate;
    if (h == 0.0) error("Impossible interpolation");
    a=(rate*khi-x)/h;
    b=(x-klo*rate)/h;
    *yy=a*y[klo]+b*y[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0;
    yy++;
  }
  
  
}
