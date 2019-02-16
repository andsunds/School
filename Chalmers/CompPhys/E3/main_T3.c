#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "funcs.h"

#define PI 3.141592653589
#define PI_To_minus_ThreeHalfs 0.179587122125167

/* Main program */
int main()
{
  gsl_rng *q;  // rng  instance 
  const  gsl_rng_type *T;    // static  info  about  rngs  
  gsl_rng_env_setup ();      // setup  the  rngs 
  T = gsl_rng_default;       // specify  default  rng 
  q = gsl_rng_alloc(T);      // allocate  default  rng 
  gsl_rng_set(q,time(NULL)); // Initialize  rng 
  
  int N = 1e6;
  int N_equilibrate = 1e3;
  int N_accept = 0;
 
  double mean_value = 0; 
  //double sigma_s = 0;
  double exact_value = 7/((double) 8);
  clock_t startTime;
  double spentTime;
  
  printf("true value: \t I = %f \n", exact_value); 
  
  double x1 = 0; double y1 = 0; double z1 = 0;
  double x2 = 0; double y2 = 0; double z2 = 0;
  double xi; // new random
  
  double delta = 2; 
  double weight_function_ratio;
  double f;
  
  startTime = clock();
  
  //equilibrate
  for (int i = 0; i < N_equilibrate; i++){
    x2 = x1 + delta * ( gsl_rng_uniform(q) - 0.5 );
    y2 = y1 + delta * ( gsl_rng_uniform(q) - 0.5 );
    z2 = z1 + delta * ( gsl_rng_uniform(q) - 0.5 );
    xi = gsl_rng_uniform(q);
    
    weight_function_ratio = get_weightfunction(x2, y2, z2)/
      get_weightfunction(x1, y1, z1);
    if (weight_function_ratio > xi){
      x1 = x2;
      y1 = y2;
      z1 = z2;
    }  
  }
  
  // real run
  for (int i = 0; i < N; i++){
    x2 = x1 + delta * ( gsl_rng_uniform(q) - 0.5 );
    y2 = y1 + delta * ( gsl_rng_uniform(q) - 0.5 );
    z2 = z1 + delta * ( gsl_rng_uniform(q) - 0.5 );
    xi = gsl_rng_uniform(q);
    
    weight_function_ratio = get_weightfunction(x2, y2, z2)/
      get_weightfunction(x1, y1, z1);
    if (weight_function_ratio > xi){
      x1 = x2;
      y1 = y2;
      z1 = z2;
      N_accept++;
      //printf("accepted!! x = %.e, y = %.e, z = %.e\n", x1,y1,z1);
    }
    f =  get_normalized_integrand(x1, y1, z1);
  
    mean_value += f;
    //sigma_s += f * f;   
  }

  mean_value *= 1/(double)N;
  //sigma_s = sigma_s/(double)N - mean_value * mean_value;
  //sigma_s = sqrt(sigma_s/(N-1));

  spentTime = (double)(clock() - startTime)/CLOCKS_PER_SEC*1e3;
    
  printf("I_N = %.4f, \t I-I_N = %.2e, \t t = %.1f ms\n", 
    mean_value, exact_value-mean_value, spentTime);
  
  printf("N_accept = %d (%.1f %%) \n", N_accept, (double)N_accept/N*100); 
  gsl_rng_free(q); // deallocate  rng 
  
  
  return 0;    
}
