
#define NMIN       2
#define NMAX       10

/****************************************************************************/

#define max( a, b ) 	( (a) > (b) ? (a) : (b) )
#define min( a, b ) 	( (a) < (b) ? (a) : (b) )
#define minus1to( n ) 	( ( (n) % 2 == 0 ) ? (1) : -1 )
#define STRING_SIZE 256




/****************************
* Global Variables          *
****************************/

extern double *a, **c;
extern int NW;
extern int *twoto;

/****************************
* Structure Definition      *
****************************/
typedef struct
{
  int lb;
  int ub;
  int size;
} bound;

/****************************************************************************/

int open_read( void );
int compute_a( void );

double phi( double x );
double psi( double x );

void init_twoto( int max_resoln );
void init_phi_array( float **phi_array, int max_resoln );
void init_psi_array( float **psi_array, int max_resoln );

