
float p_value( float T, float **histo, int resoln, int histo_size );

void compute_pval_average( float *pval, float **p, int max_resoln, int np, 
			   int num_of_windows, int window_size );

void local_mean( float *mean, float *s, int np );

float variance( float *s, int np );

void local_var( float *var, float *s, int *resoln_ptr, int *np_ptr );
