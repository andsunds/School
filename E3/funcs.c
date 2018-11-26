#include "funcs.h"

double get_f1(double x){
	return x*(1-x);
}

double get_f2(double x){
  double f2;
  double epsilon = 1e-10;
  if ( (x > epsilon) && ( (1-x) > epsilon) ) { 
    f2 =  2 * x*(1-x)/( PI*sin(PI*x) );
	}else{ // the limits where x<<1 and 1-x << 1 are tricky to resolve numerically
	  printf("x = %f\n",x); 
	  f2 = 2/(PI*PI);
	}
	return f2;
}

void set_uniform_random(int N, double vec[N])
{
  const  gsl_rng_type *T;    // static  info  about  rngs 
  gsl_rng *q;                // rng  instance 
  gsl_rng_env_setup ();      // setup  the  rngs 
  T = gsl_rng_default;       // specify  default  rng 
  q = gsl_rng_alloc(T);      // allocate  default  rng 
  gsl_rng_set(q,time(NULL)); // Initialize  rng 
  
  // Loops over all the elemtens in the matrix, to which we want to add noise
  for (int i=0; i<N; i++){
    vec[i] = gsl_rng_uniform(q); 
  }
  gsl_rng_free(q); // deallocate  rng 
} 

void set_zero (int M, int N, double mat[M][N]){
  /* Sets all the elements of matrix `mat` to zero */
  //loops over all indices
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      mat[i][j] = 0;
    } 
  }
} 

