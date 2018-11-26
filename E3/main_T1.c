#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "funcs.h"

#define PI 3.141592653589

/* Main program */
int main()
{
  double Ns[] = {1e1,1e2,1e3,1e4};
  int N_Ns = 4;
  double (*function_pointer)(double);
   
  double mean_value; 
  double sigma_s;
  double exact_value = 1/((double) 6);
  double y;
  clock_t startTime;
  double spentTime;
  
  double *x = malloc(Ns[0]*sizeof(double));
  printf("true value: \t I = %f \n", exact_value); 
  
  
  printf("Monte Carlo:  \n"); 
  function_pointer = &get_f1;
  for (int i = 0; i < N_Ns; i++){
    mean_value = 0;
    sigma_s = 0;
    startTime = clock();
    
    if (i>0){ // need to change size of x and f
      x = realloc(x, Ns[i]*sizeof(double));
      
    }
    set_uniform_random(Ns[i], x);
  
    for( int j = 0; j<Ns[i]; j++){
      mean_value += (*function_pointer)(x[j]);
      sigma_s += (*function_pointer)(x[j]) * (*function_pointer)(x[j]);
    }
    mean_value *= 1/Ns[i];
    sigma_s = sigma_s*1/Ns[i] - mean_value * mean_value;
    sigma_s = sqrt(sigma_s/(Ns[i]-1));
    
    spentTime = (double)(clock() - startTime)/CLOCKS_PER_SEC*1e3;
    
    printf("N= %.0e, \t I_N = %.4f +- %.4f, \t I-I_N = %.1e, \t t = %.4f ms\n", (double)Ns[i], 
      mean_value, sigma_s, exact_value-mean_value, spentTime);
  }
  
  // with importance sampling
  printf("With importance sampling:  \n"); 
  function_pointer = &get_f2; // f(x)/p(x)
  for (int i = 0; i < N_Ns; i++){
    mean_value = 0;
    sigma_s = 0;
    startTime = clock();
    
    x = realloc(x, Ns[i]*sizeof(double));
    set_uniform_random(Ns[i], x);
    
    for( int j = 0; j<Ns[i]; j++){
      y = 1/PI * acos( 1 - 2*x[j]); // F^-1 (x)
      mean_value += (*function_pointer)(y);
      sigma_s += (*function_pointer)(y) * (*function_pointer)(y);
    }
    mean_value *= 1/Ns[i];
    sigma_s = sigma_s*1/Ns[i] - mean_value * mean_value;
    if (sigma_s<0){
        printf("NOOOO! %f \n", sigma_s );
    }
    sigma_s = sqrt(sigma_s/Ns[i]);
    
    spentTime = (double)(clock() - startTime)/CLOCKS_PER_SEC*1e3;
    printf("N= %.0e, \t I_N = %.4f +- %.4f, \t I-I_N = %.4f, \t t = %.1e ms\n", (double)Ns[i], 
      mean_value, sigma_s, exact_value-mean_value, spentTime);
  }
  
  
  free(x);
  return 0;    
}
