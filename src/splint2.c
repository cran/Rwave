/***************************************************************
*    $Log: splint2.c,v $                                       *
****************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
***************************************************************/


#include "kernel.h"


/***************************************************************
*  Function: splint2
*  ---------
*      Cubic spline interpolation (modified from Numerical
*       Recipes, in order to incorporate the computation
*       of the first derivative).
*
*    xa,xb: arrays containing the x and y values at the nodes 
*    ya2: second derivative.
*    x: value at which the function is to be estimated.
*    y: value of the function at point x.
*    yp: derivative of the function at point x.
***************************************************************/

void splint2(double xa[], double ya[], double y2a[], int n, double x, double *y, double *yp)
{
	int klo,khi,k;
	double h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) {
	  printf("Bad xa input to routine splint2 \n");
	  exit(1);
	}
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
	*yp = (ya[khi]-ya[klo] - ((3*a*a -1)*y2a[klo] - (3*b*b -1)*y2a[khi])*h*h/6)/h;

}








