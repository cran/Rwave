/***********************************/
/* Swave.h Some Basic include file */
/***********************************/

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "R.h"

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
#define ZERO 0.000001

#define STRING_SIZE 256

#define max( a, b ) 	( (a) > (b) ? (a) : (b) )
#define min( a, b ) 	( (a) < (b) ? (a) : (b) )
#define inrange(inf,x,sup) ((inf) <= (x) && (x) <= (sup))

/*****************************/
/* Type Definition           */
/*****************************/
/*
typedef struct
{
  int lb;			 lower_bound 
  int ub;			 upper_bound 
  int size;
} bound;

*/

/* structure of an image extrema used in point reconst */

/* typedef struct
{
  int resoln;
  int x; 
  int y; 
  double W1f; 
  double W2f; 
} image_ext;
*/

extern double my_exp2();
extern double my_log2();
extern int nint();
extern double log2();
extern double exp2();
extern int iexp2();
extern double fexp2();
extern int find2power();

/* uniform random number generator */
extern double urand_();
