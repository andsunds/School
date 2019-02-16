#include "funcs.h"

gsl_rng* initialize_rng()
{
  gsl_rng *q;
  const  gsl_rng_type *T;    // static  info  about  rngs  
  gsl_rng_env_setup ();      // setup  the  rngs 
  T = gsl_rng_default;       // specify  default  rng 
  q = gsl_rng_alloc(T);      // allocate  default  rng 
  gsl_rng_set(q,time(NULL)); // Initialize  rng 
  return q;
} 

void take_BD_step(double *a, double *v, double *x, 
    double c0, gsl_rng *q, double v_th, double dt, int N_part, double omega0){
  double sqrt_c0 = sqrt(c0);
  double sqrt_1c0 = sqrt(1-c0);
  double omega0_sq = omega0 * omega0;
  
  for (int i = 0; i < N_part; i++) {
    double G1 = gsl_ran_gaussian(q,1);
    double G2 = gsl_ran_gaussian(q,1);
    v[i] = dt*0.5*a[i] + sqrt_c0 * v[i] + v_th*sqrt_1c0*G1;
    x[i] += v[i]* dt;
    a[i] = -x[i] * omega0_sq;
    v[i] = dt*0.5*sqrt_c0*a[i] + sqrt_c0*v[i] + v_th*sqrt_1c0*G2;
  }
}

