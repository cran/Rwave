/***************************************************************
*    $Log: compinteg.c,v $                                     *
****************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
***************************************************************/


#include "kernel.h"
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5


/***************************************************************
*  Function: qrombmod
*  ---------
*       Romberg quadrature for integration; 
*     complex-valued functions
*
*   x,y: position of the wavelets in the integrand.
*   b_start, b_end: integration domain.
*   nodes, phi_nodes: position and scale of the nodes of the ridge.
*   nb_nodes: number of nodes of the ridge.
****************************************************************/

fcomplex qrombmod(int x, int y, double *p2, double *nodes, double *phi_nodes, 
		  int nb_nodes,double cent_freq,double b_start, double b_end)
{
  double h[JMAXP+1];
  fcomplex ss,dss, tmp3;
  fcomplex *s;
  int j;
  double tmpr[JMAX+1],tmpi[JMAX+1], tmp1, tmp2;


  s = (fcomplex *)calloc(JMAXP+1,sizeof(fcomplex));

  /* Initialisation of tmp arrays */
  for(j = 0; j < JMAX+1; j++)
    tmpi[j] = tmpr[j] = 0.0;

  h[1]=1.0;

  /* Loop and test accuracy
     ----------------------*/
  for (j=1;j<=JMAX;j++) {

    /* Trapezoidal rule */
    s[j] = trapzdmod(x,y,p2,nodes,phi_nodes,nb_nodes,cent_freq,b_start,b_end,j);
    tmpr[j] = (s[j]).r;    tmpi[j] = (s[j]).i;
/*    printf("step=%d ,integral.r=%g12 ,integral.i=%g\n",j,(s[j]).r,(s[j]).i); */

    if (j >= K) {

      /* Polynomial (Lagrange) extrapolation to estimate the limit of
	 the integral as step goes to 0 */
      polint(&h[j-K],&tmpr[j-K],(int)K,0.0,&(ss.r),&(dss.r));
      polint(&h[j-K],&tmpi[j-K],(int)K,0.0,&(ss.i),&(dss.i));

      /* Test accuracy */
      if (((fabs(dss.r) < EPS*fabs(ss.r))&&(EPS * fabs(ss.r) > fabs(ss.i))) ||
           ((fabs(dss.i) < EPS*fabs(ss.i))&&(EPS * fabs(ss.i) > fabs(ss.r))) ||
           ((fabs(dss.r) < EPS*fabs(ss.r))&&(fabs(dss.i) < EPS*fabs(ss.i)))) {
	free((char *)s);
	return ss;
      }
    }
    (s[j+1]).r=(s[j]).r;
    (s[j+1]).i=(s[j]).i;
    h[j+1]=0.25*h[j];
  }
  printf("Too many steps in routine qrombmod (x=%d, y=%d) \n",x,y);
  free((char *)s);
  return(ss);
}




/***************************************************************
*  Function: rqrombmod
*  ---------
*       Romberg quadrature for integration of WT kernel; 
*     complex-valued functions
*
*   x,y: position of the wavelets in the integrand.
*   b_start, b_end: integration domain.
*   nodes, phi_nodes: position and scale of the nodes of the ridge.
*   nb_nodes: number of nodes of the ridge.
****************************************************************/

double rqrombmod(int x, int y, double *p2, double *nodes, double *phi_nodes, 
		  int nb_nodes,double cent_freq,double b_start, double b_end)
{
  double h[JMAXP+1];
  double ss,dss, tmp3;
  double *s;
  int j;
  double tmpr[JMAX+1], tmp1, tmp2;


  s = (double *)calloc(JMAXP+1,sizeof(double));

  /* Initialisation of tmp arrays */
  for(j = 0; j < JMAX+1; j++)
    tmpr[j] = 0.0;

  h[1]=1.0;

  /* Loop and test accuracy
     ----------------------*/
  for (j=1;j<=JMAX;j++) {

    /* Trapezoidal rule */
    s[j] = rtrapzdmod(x,y,p2,nodes,phi_nodes,nb_nodes,cent_freq,b_start,b_end,j);
    tmpr[j] = s[j];
/*    printf("step=%d ,integral.r=%g12 ,integral.i=%g\n",j,(s[j]).r,(s[j]).i); */

    if (j >= K) {

      /* Polynomial (Lagrange) extrapolation to estimate the limit of
	 the integral as step goes to 0 */
      polint(&h[j-K],&tmpr[j-K],(int)K,0.0,&ss,&dss);

      /* Test accuracy */
      if (fabs(dss) < EPS*fabs(ss)) {
	free((char *)s);
	return ss;
      }
    }
    s[j+1]=s[j];
    h[j+1]=0.25*h[j];
  }
  printf("Too many steps in routine rqrombmod (x=%d, y=%d) \n",x,y);
  free((char *)s);
  return(ss);
}



/***************************************************************
*  Function: gqrombmod
*  ---------
*       Romberg quadrature for integration; 
*     real-valued functions
*
*   x,y: position of the wavelets in the integrand.
*   b_start, b_end: integration domain.
*   nodes, phi_nodes: position and scale of the nodes of the ridge.
*   nb_nodes: number of nodes of the ridge.
****************************************************************/

double gqrombmod(int x, int y, double *p2, double *nodes, double *phi_nodes, 
		  int nb_nodes,double scale,double b_start, double b_end)
{
  double h[JMAXP+1];
  double ss,dss, tmp3;
  double *s;
  int j;
  double tmpr[JMAX+1], tmp1, tmp2;


  s = (double *)calloc(JMAXP+1,sizeof(double));

  /* Initialisation of tmp arrays */
  for(j = 0; j < JMAX+1; j++)
    tmpr[j] = 0.0;

  h[1]=1.0;

  /* Loop and test accuracy
     ----------------------*/
  for (j=1;j<=JMAX;j++) {

    /* Trapezoidal rule */
    s[j] = gtrapzdmod(x,y,p2,nodes,phi_nodes,nb_nodes,scale,b_start,b_end,j);
    tmpr[j] = s[j];

    if (j >= K) {

      /* Polynomial (Lagrange) extrapolation to estimate the limit of
	 the integral as step goes to 0 */
      polint(&h[j-K],&tmpr[j-K],(int)K,0.0,&ss,&dss);

      /* Test accuracy */
      if (fabs(dss) < EPS*fabs(ss)) {
	free((char *)s);
	return ss;
      }
    }
    (s[j+1])=(s[j]);
    h[j+1]=0.25*h[j];
  }
  printf("Too many steps in routine gqrombmod (x=%d, y=%d) \n",x,y);
  free((char *)s);
  return(ss);
}



#undef EPS
#undef JMAX
#undef JMAXP
#undef K




/***************************************************************
*  Function: trapzdmod
*  ---------
*       Integration with trapezoidal rule; 
*     complex-valued functions
*
*   x,y: position of the wavelets in the integrand.
*   b_start, b_end: integration domain.
*   nodes, phi_nodes: position and scale of the nodes of the ridge.
*   nb_nodes: number of nodes of the ridge.
*   n: order of the integration.
****************************************************************/


fcomplex trapzdmod(int x, int y, double *p2, double *nodes, double *phi_nodes, 
		   int nb_nodes,double cent_freq,double b_start, double b_end,
		   int n)
{
  double xx,tnm,del;
  static fcomplex s;
  fcomplex ctmp,sum;
  int it,j;
  
  if (n == 1) {
    ctmp = integrand(b_start,x,y,p2,nodes,phi_nodes,nb_nodes,cent_freq);
    ctmp = Cadd(ctmp,integrand(b_end,x,y,p2,nodes,phi_nodes,nb_nodes,cent_freq));
    return (s=Cmul(Complex(0.5*(b_end-b_start),0),ctmp));
  } 
  else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b_end-b_start)/tnm;
    xx=b_start+0.5*del;
    for (sum=Complex(0.0,0.0),j=1;j<=it;j++,xx+=del) 
      sum = Cadd(sum,integrand(xx,x,y,p2,nodes,phi_nodes,nb_nodes,cent_freq));
    ctmp = Cmul(Complex((b_end-b_start)/tnm,0),sum);
    ctmp = Cadd(ctmp,s);
    s = Cmul(Complex(0.5,0),ctmp);
    return s; 
  } 
}



/***************************************************************
*  Function: rtrapzdmod
*  ---------
*       Integration with trapezoidal rule; 
*     real-valued functions for wavelet transform
*
*   x,y: position of the wavelets in the integrand.
*   b_start, b_end: integration domain.
*   nodes, phi_nodes: position and scale of the nodes of the ridge.
*   nb_nodes: number of nodes of the ridge.
*   n: order of the integration.
****************************************************************/


double rtrapzdmod(int x, int y, double *p2, double *nodes, double *phi_nodes, 
		   int nb_nodes,double cent_freq,double b_start, double b_end,
		   int n)
{
  double xx,tnm,del;
  static double s;
  double ctmp,sum;
  int it,j;
  
  if (n == 1) {
    ctmp = rintegrand(b_start,x,y,p2,nodes,phi_nodes,nb_nodes,cent_freq);
    ctmp += rintegrand(b_end,x,y,p2,nodes,phi_nodes,nb_nodes,cent_freq);
    return (s=0.5*(b_end-b_start)*ctmp);
  } 
  else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b_end-b_start)/tnm;
    xx=b_start+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,xx+=del) 
      sum += rintegrand(xx,x,y,p2,nodes,phi_nodes,nb_nodes,cent_freq);
    ctmp = sum *(b_end-b_start)/tnm;
    ctmp += s;
    s = 0.5*ctmp;
    return s; 
  } 
}



/***************************************************************
*  Function: gtrapzdmod
*  ---------
*       Integration with trapezoidal rule (case of 
*       Gabor transform); 
*     real-valued functions
*
*   x,y: position of the wavelets in the integrand.
*   b_start, b_end: integration domain.
*   nodes, phi_nodes: position and scale of the nodes of the ridge.
*   nb_nodes: number of nodes of the ridge.
*   n: order of the integration.
****************************************************************/


double gtrapzdmod(int x, int y, double *p2, double *nodes, double *phi_nodes, 
		  int nb_nodes,double scale,double b_start, double b_end, 
		  int n)
{
  double xx,tnm,del;
  static double s;
  double ctmp,sum;
  int it,j;
  
  if (n == 1) {
    ctmp = gintegrand(b_start,x,y,p2,nodes,phi_nodes,nb_nodes,scale);
    ctmp += gintegrand(b_end,x,y,p2,nodes,phi_nodes,nb_nodes,scale);
    return (s=0.5*(b_end-b_start)*ctmp);
  } 
  else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b_end-b_start)/tnm;
    xx=b_start+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,xx+=del) 
      sum += gintegrand(xx,x,y,p2,nodes,phi_nodes,nb_nodes,scale);
    ctmp = sum *(b_end-b_start)/tnm;
    ctmp += s;
    s = 0.5*ctmp;
    return s; 
  } 
}






