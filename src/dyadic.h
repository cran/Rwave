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
  float W1f; /* W1f */
  float W2f; /* W2f */
} image_ext;


/* in mw.c
   ------- */
void Sf_compute(float *Sf,float *f,int *max_resoln_ptr,
  int *np_ptr,char **pfiltername);

void Wf_compute(float *Wf, float *Sf, int *max_resoln_ptr,
  int *np_ptr, char **pfiltername);



/* in dwfilter.c
   ------------- */
int iexp2(int j);

void Hfilter_compute(char *filtername, float ***H,
  bound *H_bound, int max_resoln );

void Gfilter_compute(char *filtername, float ***G,
  bound *G_bound, int max_resoln);

void HGfilter_bound(char *filtername, bound **H_bound,
  bound **G_bound, int max_resoln );

void HG_hat_compute(char *filtername, float ***H_hat,
  float ***G_hat, int max_resoln, int np);

void Sfilter_compute(char *filtername, float ***S, bound *S_bound,
  int max_resoln);

void Kfilter_compute(char *filtername, float ***K, bound *K_bound,
  int max_resoln);

void Lfilter_compute(char *filtername, float ***L, bound *L_bound,
  int max_resoln);

void KSfilter_bound(char *filtername, bound **K_bound,
  bound **S_bound, int max_resoln);

void Lfilter_bound(char *filtername, bound **L_bound,
  int max_resoln);

void PsiPhifilter_bound(bound **psi, bound **phi, bound *H_bound,
  bound *G_bound, int max_resoln);

void signal_W_S(float ***W, float ***S, int max_resoln, int np);

void signal_W_hat_S_hat(float ***W_hat, float ***S_hat,
  int max_resoln, int np);



/* in dwvector.c
   ------------- */
void compute_convolution(float *s, float *f, float *g, int np );


/* in dwfileio.c
   ------------- */
void input_signal(char *fname, float **Pic, int size);

void init_filename(char filename[]);

void filename_given(char filename[], char *name);

void strconcate(char s[] , char t[], char buff[]);

void filename_inc(char filename[], int inc);

void output_signal(float *s, int np, char *fname);

void output_array(float **array, int nrow, int ncol, char *file_name );



/* in extrema.c
   ------------ */
void modulus_maxima(float *extrema, float *wt, int *resoln_ptr,
  int *np_ptr );

void extrema_input(float *extrema, int max_resoln, int np,
  image_ext **ext, int *num_of_extrema);



/* in m_reconst.c
   -------------- */
void signal_penalty_function(float *f, float *lambda,
  float **W_tilda, image_ext *ext, int num_of_extrema, int np);

void signal_position(char *filtername, float **lambda,
  image_ext *ext, float **Wtilda, float **W, int num_of_extrema,
  int max_resoln, int np);

void extrema_reconst(char *filtername, float *f, float *extrema,
  int *max_resoln_ptr, int *np_ptr, int *preadfile);



/* in dualwavelet.c
   ---------------- */
void signal_W_tilda(float ***W_tilda, float **W, float **K,
		    int max_resoln, int np);

void signal_W_tilda_input(float ***W_tilda, int max_resoln, int np);



/* in svd.c
   -------- */
double pythag(double a, double b);

void svdcmp(double **a, int m, int n, double *w, double **v);

void svbksb(double **U, double *W, double **V, int m, int n,
  double *B, double *X);

void svdecomp_solve(float **a, float *b, float *x, int m,
  int n, float **w, float ***v);

void residue(float **a, float *w, float **v, int m, int n,
  float *b, float *x);

void double_residue(double **a, double *w, double **v, int m,
  int n, double *b, double *x);

void Sresidue(float *a, float *w, float *v, int m, int n,
  float *b, float *x);

void Ssvdecomp(float *a, int *pm, int *pn, float *w, float *v,
  float *b, float *x);



/* in dwkernel.c
   ------------- */
float fexp2(int j);

void wavelet_transform_gradient(float **grad, float **s,
  int max_resoln, int np );

void signal_K_compute(float ***K, float **W, int max_resoln,
  int np );

void signal_tilda_adjust(float **tilda, int ksize, char *fname,
  int fsize);


