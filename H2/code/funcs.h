#ifndef  _FUNCS_H
#define  _FUNCS_H

//  includes
#include  <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>  // also needed for random stuff

#define N_neigh 8

//  prototypes

/************************ get functions **************************************/
double get_bond_E(int site_1, int site_2);

double get_order_parameter(int *lattice, int N_Cu);

double get_short_range_order_parameter(int *lattice, int(*nearest)[N_neigh],
				       int N_Cu);

double get_Etot(int *lattice, int N_atoms, int (*nearest)[N_neigh]);

void get_phi (double *phi, int N_times, double f_mean, double f_var,
	      double *data,int N_k, int N_skip);

void get_varF_block_average(double *var_F, int N_times, double f_mean,
			    double f_var, double *data, int N_k, int N_skip);


/************** Monte Carlo step functions ************************************/
void MC_step( double *Etot, double *r, double *P, gsl_rng *q,
              int *lattice, int (*nearest)[N_neigh], double beta, int N_Cu);

void update_E_P_r(int iT, double E_dev, double *E_mean, double *E_sq_mean,
		  double P, double *P_mean, double *P_sq_mean,
		  double r, double *r_mean,double *r_sq_mean,
		  int *lattice, int (*nearest)[N_neigh], int N_Cu);


/************************* initializing functions******************************/
void* init_temps(int *nT, double dT_small, double dT_large,
		 double T_start, double T_end, double T_start_fine,
		 double T_end_fine);

void init_nearestneighbor(int Nc, int (*nearest)[N_neigh]);

void init_ordered_lattice(int N_atoms, int N_Cu, int *lattice);

void init_random_lattice(int N_atoms, int N_Cu, int *lattice, gsl_rng *q);

void* init_random();


/************************* file I/O functions *********************************/
void write_equil_to_file(double T_degC, double *E_equilibration, int N_bonds,
			 double *P, int N_eq);

void write_production(double *T_degC, int nT, double *E_mean_approx,
		      double *E_mean, double *E_sq_mean,
		      double *P_mean, double *P_sq_mean,
		      double *r_mean, double *r_sq_mean);

void write_stat_inefficiency_to_file(double T_degC, double *phi, double *var_F,
				     int N_k,int N_skip);

#endif
