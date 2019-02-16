#ifndef  _FUNCS_H
#define  _FUNCS_H

#define PI 3.141592653589
#define PI_To_minus_ThreeHalfs 0.179587122125167
//  includes
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>  // also needed for random stuff
#include <stdio.h>
#include <stdlib.h>

//  prototypes

gsl_rng* initialize_rng();

void take_BD_step(double *a, double *v, double *x, 
    double c0, gsl_rng *q, double v_th, double dt, int N_particles, double omega0);


#endif
