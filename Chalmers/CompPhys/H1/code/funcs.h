#ifndef  _FUNCS_H
#define  _FUNCS_H

//  includes
#include  <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>  // also needed for random stuff
#include "alpotential.h"
#include "fft_func.h"

//  prototypes

void add_noise(int N, int M, double mat[M][N], double noise_amplitude );

void timestep_Verlet ( int N_atoms, double (*pos)[3],
			      double (*momentum)[3],
			      double (*forces)[3], double m, double dt,
			      double cell_length);

double get_kin_energy ( int N_atoms,  double (*momentum)[3], double m );

void get_displacements ( int N_atoms,  double (*positions)[3],
				double (*initial_positions)[3], double disp[]);
				
void get_MSD ( int N_atoms,  int N_times, double all_pos[N_times][N_atoms][3], 
                 double MSD[N_times]);				

void get_vel_corr ( int N_atoms,  int N_times, double all_mom[N_times][N_atoms][3], 
                 double vel_corr[N_times]);
                 
void get_powerspectrum ( int N_atoms,  int N_times, double all_vel[N_times][N_atoms][3], 
                 double pow_spec[N_times]);                 

void copy_mat (int M, int N, double mat_from[M][N], double mat_to[M][N]);

void set_zero (int M, int N, double mat[M][N]);

void scale_mat (int M, int N, double mat[M][N], double alpha);
#endif
