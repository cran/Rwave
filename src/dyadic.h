/****************************************************************
*	$ Log: dyadic.h,v	$				*
*               (c) Copyright  1995                             *
*                          by                                   *
*      Author: Rene Carmona, Bruno Torresani, Wen-Liang Hwang   *
*                  University of California, Irvine             *
*                  All right reserved                           *
****************************************************************/


/****************************************************************
* Type Definition           
****************************************************************/

typedef struct
{
  int lb;                       /* lower_bound */
  int ub;                       /* upper_bound */
  int size;
} bound;


typedef struct
{
  int resoln;
  int x; /* coordinate col from 0 ... */
  int y; /* coordinate row from 0 ... */
  double W1f; /* W1f */
  double W2f; /* W2f */
} image_ext;



/* in choldc.c
   ------- */
void double_choldc(double **a, int n, double p[]);
void cholsl(double **a, int n, double p[], double b[], double x[]);
void choldc(double **a, int n, double p[]);

/* in mw.c
   ------- */
void Sf_compute(double *Sf,double *f,int *max_resoln_ptr,
  int *np_ptr,char *filtername);

void Wf_compute(double *Wf, double *Sf, int *max_resoln_ptr,
  int *np_ptr, char *filtername);



/* in dwfilter.c
   ------------- */
int iexp2(int j);

void Hfilter_compute(char *filtername, double ***H,
  bound *H_bound, int max_resoln );

void Gfilter_compute(char *filtername, double ***G,
  bound *G_bound, int max_resoln);

void HGfilter_bound(char *filtername, bound **H_bound,
  bound **G_bound, int max_resoln );

void HG_hat_compute(char *filtername, double ***H_hat,
  double ***G_hat, int max_resoln, int np);

void Sfilter_compute(char *filtername, double ***S, bound *S_bound,
  int max_resoln);

void Kfilter_compute(char *filtername, double ***K, bound *K_bound,
  int max_resoln);

void Lfilter_compute(char *filtername, double ***L, bound *L_bound,
  int max_resoln);

void KSfilter_bound(char *filtername, bound **K_bound,
  bound **S_bound, int max_resoln);

void Lfilter_bound(char *filtername, bound **L_bound,
  int max_resoln);

void PsiPhifilter_bound(bound **psi, bound **phi, bound *H_bound,
  bound *G_bound, int max_resoln);

void signal_W_S(double ***W, double ***S, int max_resoln, int np);

void signal_W_hat_S_hat(double ***W_hat, double ***S_hat,
  int max_resoln, int np);



/* in dwvector.c
   ------------- */
void compute_convolution(double *s, double *f, double *g, int np );
void signal_zero( double *input, int size);
void signal_copy( double *input, double *output, int size, int offset);
void complex_product( double *product, double *s1, double *s2, int np);




/* in dwfileio.c
please don't write to disk
   ------------- 
void input_signal(char *fname, double **Pic, int size);

void init_filename(char filename[]);

void filename_given(char filename[], char *name);

void strconcate(char s[] , char t[], char buff[]);

void filename_inc(char filename[], int inc);

void output_signal(double *s, int np, char *fname);

void output_array(double **array, int nrow, int ncol, char *file_name );
*/


/* in extrema.c
   ------------ */
void modulus_maxima(double *extrema, double *wt, int *resoln_ptr,
  int *np_ptr );

void extrema_input(double *extrema, int max_resoln, int np,
  image_ext **ext, int *num_of_extrema);



/* in m_reconst.c
   -------------- */
void signal_penalty_function(double *f, double *lambda,
  double **W_tilda, image_ext *ext, int num_of_extrema, int np);

void signal_position(char *filtername, double **lambda,
  image_ext *ext, double **Wtilda, double **W, int num_of_extrema,
  int max_resoln, int np);

void extrema_reconst(char *filtername, double *f, double *extrema,
  int *max_resoln_ptr, int *np_ptr, int *preadfile);



/* in dualwavelet.c
   ---------------- */
void signal_W_tilda(double ***W_tilda, double **W, double **K,
		    int max_resoln, int np);

void signal_W_tilda_input(double ***W_tilda, int max_resoln, int np);



/* in svd.c
   -------- */
double pythag(double a, double b);

void svdcmp(double **a, int m, int n, double *w, double **v);

void svbksb(double **U, double *W, double **V, int m, int n,
  double *B, double *X);

void svdecomp_solve(double **a, double *b, double *x, int m,
  int n, double **w, double ***v);

void residue(double **a, double *w, double **v, int m, int n,
  double *b, double *x);

void double_residue(double **a, double *w, double **v, int m,
  int n, double *b, double *x);

void Sresidue(double *a, double *w, double *v, int m, int n,
  double *b, double *x);

void Ssvdecomp(double *a, int *pm, int *pn, double *w, double *v,
  double *b, double *x);



/* in dwkernel.c
   ------------- */
double fexp2(int j);

void wavelet_transform_gradient(double **grad, double **s,
  int max_resoln, int np );

void signal_K_compute(double ***K, double **W, int max_resoln,
  int np );

/* please don't write to disk
void signal_tilda_adjust(double **tilda, int ksize, char *fname,
  int fsize);
*/



