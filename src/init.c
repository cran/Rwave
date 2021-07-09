
#include <R_ext/Rdynload.h>
#include <R.h>
#include <Rinternals.h>
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL

#include "denoise.h"
#include "dyadic.h"

void Shessianmap(double *sqmodulus, int *psigsize, int *pnscale, int *pnbblock,
                 int *pgridx, int *pgridy, double *tst);


void Spca_annealing(double *smodulus, double *beemap,
               int *crazymap,
               double *pc,
               int *psigsize, int *pnscale, int *piteration,
               int *pseed, int *pbstep, int *pnbbee,
		    int *pintegral, int *pchain, int *flag);

void Spca_family(double *ridgemap,int *orientmap, double *orderedmap,int *chain,
		   int *pnbchain, int *psigsize,int *pnscale,
		 int *pbstep,double *pthreshold, int* pmaxchnlng) ;

void Ssnakenoid_annealing(double *cost, double *smodulus,
  double *phi, double *rho, double *plambda, double *pmu,
  double *plambda2, double *pmu2, double *pc, int *psigsize,
  int *psnakesize, int *pnscale, int *piteration,
  int *pstagnant, int *pseed, int *pcount, int *psub,
			  int *pblocksize, int *psmodsize);
  
void daubechies_wt(void *,void *,void *,void *,void *,void *);

		

void compute_ddwave( double *phi, double *psi, double *s, int *max_resoln_ptr, int *np_ptr, int *NW_ptr );

void Sinverse_wavelet_transform(double *f_back,double *Sf,double *Wf,int *pmax_resoln,int *pnp,char **pfiltername);

    
void Spointmap(double *sqmodulus, int *psigsize, int *pnscale,
	       int *pgridx, int *pgridy, int *pnbblock, int *pnbpoint, 
	       int *pointmap, double *tst, int *ptstsize, int *pcount, int *pseed);

void Stf_pcaridge(double *input, double *output, int *pnrow, int *pncol, int *orientmap);
    
void gabor_time(double *pfrequency,double *pscale, int *pb, 
		double *g_r, double *g_i,int *pisize);

void modulus_maxima(double *extrema, double *wt, int *resoln_ptr,
		    int *np_ptr );


void nthresh_compute(double *nthresh, double *s, int *maxresoln_ptr,
		     int *sample_size_ptr, double prct );

 void Svwt_thierry(double *Rinput,double *Iinput,double *Oreal,
   double *Oimage,double *pa,int *pinputsize,
		   int *pM);

   void normal_pval_compute(double *pval, double *s, int *max_resoln_ptr,
			    int *np_ptr, int *num_of_windows_ptr, int *window_size_ptr );

void bootstrap_pval_compute(double *pval, double *s, int *max_resoln_ptr,
			    int *np_ptr, int *num_of_windows_ptr, int *window_size_ptr );
void bthresh_compute(double *bthresh, double *s, int *maxresoln_ptr,
		     int *sample_size_ptr, double prct );

void rkernel(double *ker,int *px_min,int *px_max,int *px_inc,
	    int *plng, double *nodes,double *phi_nodes,int *pnb_nodes,
	     double *pw0,double *pb_start,double *pb_end);

void fastkernel(double *ker_r, double *ker_i,int *px_min,int *px_max,
	    int *px_inc, int *plng, double *nodes,double *phi_nodes,
		int *pnb_nodes,double *pw0,double *pb_start,double *pb_end);

  void gkernel(double *ker, int *px_min,int *px_max,
	    int *px_inc, int *plng, double *nodes,double *phi_nodes,
	       int *pnb_nodes,double *pscale,double *pb_start,double *pb_end);

void fastgkernel(double *ker, int *px_min,int *px_max,
	    int *px_inc, int *plng, double *nodes,double *phi_nodes,
		 int *pnb_nodes,double *pscale,double *pb_start,double *pb_end);
 



#define CDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const  R_CMethodDef CEntries[] = {
  CDEF(Svwt_thierry , 7),
  CDEF(normal_pval_compute , 6),
  CDEF(bootstrap_pval_compute ,6 ),
  CDEF(bthresh_compute ,5 ),
  CDEF( rkernel,11 ),
  CDEF(fastkernel ,12 ),
  CDEF( gkernel , 11),
  CDEF(fastgkernel ,11 ),
  CDEF(nthresh_compute, 5),
  CDEF(modulus_maxima, 4),
  CDEF(Ssmoothwt,  6),
CDEF(Smodulus_smoothing,  6),
CDEF(Sbee_annealing,  12),
CDEF(Scrazy_family, 8 ),
CDEF(Scwt_morlet,  8),
 CDEF(Svwt_morlet, 7 ),
CDEF(morlet_time, 6 ),
CDEF(vmorlet_time,  7),
CDEF(Scwt_squeezed, 7 ),
CDEF(Scwt_thierry, 8 ),
CDEF(Scwt_phase, 8 ),
CDEF(Shessianmap,  7),
CDEF(Spca_annealing,  13),
CDEF(Spointmap,  12),
CDEF(Spca_family, 10 ),
CDEF(Stf_pcaridge,  5),
CDEF(Sridge_annealing, 15 ),
CDEF(Sridge_coronoid, 15 ),
CDEF(Sridge_icm,  11),
CDEF(Ssnake_annealing,  19),
CDEF(Ssnakenoid_annealing,  19),
CDEF(Scwt_gmax,  5),
CDEF(Scwt_mridge,  4),
CDEF(Scwt_DOG, 8 ),
CDEF(Svwt_DOG,  7),
CDEF(Svgabor,  6),
CDEF(gabor_time, 6),
CDEF(vgabor_time,  7),
CDEF(daubechies_wt,  6),
CDEF(compute_ddwave, 6 ),
CDEF(Sgabor,  7),
CDEF(entropy, 5 ),
CDEF(Lpnorm,  6),
CDEF(extrema_reconst, 6 ),
CDEF(Sinverse_wavelet_transform,  6),
CDEF(Sf_compute,  5),
CDEF(Wf_compute,  5),
CDEF(Ssvdecomp, 5 ),
  CDEF(WV,6),
    {NULL, NULL, 0}
     };


 void R_init_Rwave(DllInfo *dll)
     {
        R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
        R_useDynamicSymbols(dll, FALSE);  
     }


