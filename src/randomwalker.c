#include <stdlib.h>

/***************************************************************
*              (c) Copyright  1997                             *
*                         by                                   *
*     Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                 Princeton University                         *
*                 All right reserved                           *
***************************************************************/

#include "Swave.h"
#include "denoise.h"


/* int idum = -7; */

extern long idum;



/******************************************************************
*		UNIFORM RANDOM NUMBER GENERATOR	
*******************************************************************/

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

/* ran1 (from Numerical Recipes): uniform random numbers between 0 and 1 */
double oldran1(long *idum)
{
    static long ix1,ix2,ix3;
    static double r[98];
    double temp;
    static int iff=0;
    int j;

    if (*idum < 0 || iff == 0) 
    {
	iff=1;
	ix1=(IC1-(*idum)) % M1;
	ix1=(IA1*ix1+IC1) % M1;
	ix2=ix1 % M2;
	ix1=(IA1*ix1+IC1) % M1;
	ix3=ix1 % M3;
	for (j=1;j<=97;j++) 
	{
            ix1=(IA1*ix1+IC1) % M1;
            ix2=(IA2*ix2+IC2) % M2;
            r[j]=(ix1+ix2*RM2)*RM1;
	}
	*idum=1;
    }
    ix1=(IA1*ix1+IC1) % M1;
    ix2=(IA2*ix2+IC2) % M2;
    ix3=(IA3*ix3+IC3) % M3;
    j=1 + ((97*ix3)/M3);
    if (j > 97 || j < 1)
    {
	fprintf(stderr,"RAN1: This cannot happen.\n");
	fprintf(stderr,"Exiting now.\n");
	exit(1);
    }	
    temp=r[j];
    r[j]=(ix1+ix2*RM2)*RM1;
    return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3



/***************************************************************
* Function: randomwalker
* ---------
*    random number between 0 and 2* sigsize 
***************************************************************/
void randomwalker(int sigsize,int *num)
{
  int conf;
  double tmp;
  
  conf = 2 * sigsize;
  tmp = (double)ran1(&idum);
  *num = (int)((double)conf * ran1(&idum));
  
  if(*num >= 2 * sigsize) *num = 2 * sigsize -1;
  return;
}  


/***************************************************************
* Function: randomwalker2
* ---------
*    random number between 0 and 2* sigsize 
***************************************************************/
void randomwalker2(int sigsize,int *num, long *seed)
{
  int conf;
  
  conf = 2 * sigsize;
  *num = floor((double)conf * ran1(seed));
  
  if(*num >= 2 * sigsize) *num = 2 * sigsize -1;
  return;
}  



/***************************************************************
* Function: randomsnaker
* ---------
*    random number between 0 and 4* sigsize 
***************************************************************/
void randomsnaker(int sigsize,int *num)
{
  int conf;
  double tmp;
  
  conf = 4 * sigsize;
  tmp = (double)ran1(&idum);
  *num = (int)((double)conf * ran1(&idum));
  
  if(*num >= 4 * sigsize) *num = 4 * sigsize -1;
  return;
}  


