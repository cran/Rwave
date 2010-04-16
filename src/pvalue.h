
double p_value( double T, double **histo, int resoln, int histo_size );

void compute_pval_average( double *pval, double **p, int max_resoln, int np, 
			   int num_of_windows, int window_size );

void local_mean( double *mean, double *s, int np );

double variance( double *s, int np );

void local_var( double *var, double *s, int *resoln_ptr, int *np_ptr );
