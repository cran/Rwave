/****************************************************************
*	$ Log: kernel.h,v	$				*
*               (c) Copyright  1995                             *
*                          by                                   *
*      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                  University of California, Irvine             *
*                  All right reserved                           *
****************************************************************/

#include "complex.h"
#define EPS 1.0e-3

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#ifndef Macintosh
#include <sys/file.h>
#include <sys/types.h>
#endif

#include <time.h>

#ifndef Macintosh
#include <sys/time.h>
#endif

#define YES 1
#define NO 0

#define MAX( a, b ) 	( (a) > (b) ? (a) : (b) )
#define MIN( a, b ) 	( (a) < (b) ? (a) : (b) )

#define STRING_SIZE 256

/**************************************************************************
Function declarations:
**************************************************************************/

/* In kernel.c
   ----------*/
fcomplex integrand(double b,int x,int y,double *p2,double *nodes,
		   double *phi_nodes,int nb_nodes,double w0);

void kernel(double *ker_r, double *ker_i,int *px_min,int *px_max,
	    int *px_inc, int *plng,double *nodes,double *phi_nodes, 
	    int *pnb_nodes,double *pw0,double *pb_start,double *pb_end);

void fastkernel(double *ker_r, double *ker_i,int *px_min,int *px_max,
	    int *px_inc, int *plng, double *nodes,double *phi_nodes,
	    int *pnb_nodes,double *pw0,double *pb_start,double *pb_end);

double rintegrand(double b,int x,int y,double *p2,double *nodes,
		   double *phi_nodes,int nb_nodes,double w0);

void rkernel(double *ker,int *px_min,int *px_max,int *px_inc,
	    int *plng, double *nodes,double *phi_nodes,int *pnb_nodes,
	    double *pw0,double *pb_start,double *pb_end);

double maxvalue(double *vector, int length);

void hermite_sym(fcomplex *ker,int lng);

fcomplex psi(double x,double w0);

fcomplex psiprime(double x,double w0);

/* In gkernel.c
   -----------*/
double gintegrand(double b,int x,int y,double *p2,double *nodes,
		   double *phi_nodes,int nb_nodes,double w0);

void gkernel(double *ker, int *px_min,int *px_max,
	    int *px_inc, int *plng, double *nodes,double *phi_nodes,
	    int *pnb_nodes,double *pw0,double *pb_start,double *pb_end);

void fastgkernel(double *ker,int *px_min,int *px_max,
	    int *px_inc, int *plng, double *nodes,double *phi_nodes,
	    int *pnb_nodes,double *pscale,double *pb_start,double *pb_end);

double gfunc(double x, double scale);

void ghermite_sym(double *ker,int lng);

double  gprime(double x,double scale);


/* In splint2.c
   ------------*/
void splint2(double xa[], double ya[], double y2a[], int n, 
	     double x, double *y, double *yp);

/* In spline.c
   -----------*/
void spline(double x[], double y[], int n, double yp1, double ypn,
	    double y2[]);

/* In compinteg.c
   --------------*/
fcomplex qrombmod(int x, int y, double *p2, double *nodes, double *phi_nodes, 
		  int nb_nodes,double cent_freq,double b_start, double b_end);

fcomplex trapzdmod(int x, int y, double *p2, double *nodes, double *phi_nodes,
		   int nb_nodes,double cent_freq,double b_start, double b_end,
		   int n);

double rqrombmod(int x, int y, double *p2, double *nodes, double *phi_nodes, 
		  int nb_nodes,double cent_freq,double b_start, double b_end);

double rtrapzdmod(int x, int y, double *p2, double *nodes, double *phi_nodes, 
		   int nb_nodes,double cent_freq,double b_start, double b_end,
		   int n);

double gqrombmod(int x, int y, double *p2, double *nodes, double *phi_nodes, 
		  int nb_nodes,double scale,double b_start, double b_end);

double gtrapzdmod(int x, int y, double *p2, double *nodes, double *phi_nodes, 
		   int nb_nodes,double scale,double b_start, double b_end, 
		  int n);

void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
