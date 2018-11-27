#ifndef  _FUNCS_H
#define  _FUNCS_H

#define PI 3.141592653589
#define PI_To_minus_ThreeHalfs 0.179587122125167
//  includes
#include  <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>  // also needed for random stuff
#include <stdio.h>
#include <stdlib.h>

//  prototypes

double get_f1(double x);

double get_f2(double x);

double get_weightfunction(double x, double y, double z);

double get_normalized_integrand(double x, double y, double z);

void initialize_rng(gsl_rng *q);

void set_uniform_random(int N, double vec[N]);



//void set_zero (int M, int N, double mat[M][N]);


#endif
