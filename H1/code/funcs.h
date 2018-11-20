#ifndef  _FUNCS_H
#define  _FUNCS_H

//  includes
#include  <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>  // also needed for random stuff
#include "alpotential.h"

//  prototypes

extern void add_noise(int N, int M, double mat[M][N], double noise_amplitude );

extern void timestep_Verlet ( int N_atoms, double (*pos)[3],
			      double (*momentum)[3],
			      double (*forces)[3], double m, double dt,
			      double cell_length);

extern double get_kin_energy ( int N_atoms,  double (*momentum)[3], double m );

extern void get_displacements ( int N_atoms,  double (*positions)[3],
				double (*initial_positions)[3], double disp[]);

extern void set_zero (int M, int N, double mat[M][N]);

extern void scale_mat (int M, int N, double mat[M][N], double alpha);
#endif
