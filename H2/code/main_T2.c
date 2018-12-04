/*
  H2a, Task 2
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "funcs.h"

#define Nc 10 //number of cells
#define N_neigh 8
#define degC_to_K 273.15
#define kB 8.61733e-5

/* Main program */
int main()
{   
  int N_Cu = Nc*Nc*Nc;
  int N_atoms = 2*N_Cu;
  int N_bonds = 8*N_Cu;
  double Etot, E_Var, r, P; 
  gsl_rng *q = init_random();

  // done for all saved steps
  int N_timesteps = 1e7; 
  int N_eq = 1e6; 
  int N_eq_short = 5e5;
  double *E_equilibration = malloc(sizeof(double[N_eq]));
  double *P_equilibration = malloc(sizeof(double[N_eq]));
  double *E_production = malloc(sizeof(double[N_timesteps]));
  
  // statistical inefficiency
  int N_k = 500;
  int N_skip = 1000; // k_Max = N_k * N_skip;
  double *phi = malloc(sizeof(double[N_k]));
  double *var_F = malloc(sizeof(double[N_k]));
  
  /* set Temperature steps */
  double dT_small = 2;
  double dT_large = 10;
  double T_start = -200;
  double T_end = 600;
  double T_start_fine = 410;
  double T_end_fine = 460;
  int nT;
  double *T_degC = init_temps(&nT, dT_small, dT_large, T_start, T_end, 
    T_start_fine, T_end_fine);   
  double beta; 
  /* save equilibration data and stat inefficiency at T%20 =0*/
  int T_save_step = 20;  
   // done for all temps
  double *E_mean = malloc(sizeof(double[nT]));
  double *E_mean_approx = malloc(sizeof(double[nT]));
  double *E_sq_mean = malloc(sizeof(double[nT]));
  double *P_mean = malloc(sizeof(double[nT]));
  double *P_sq_mean = malloc(sizeof(double[nT]));
  double *r_mean = malloc(sizeof(double[nT]));
  double *r_sq_mean = malloc(sizeof(double[nT]));
  
  // initialize lattice 
  int (*nearest)[N_neigh] = malloc(sizeof(int[N_atoms][N_neigh])); // nearest neighbors
  int *lattice = malloc(sizeof(int[N_atoms]));
  init_nearestneighbor(Nc, nearest);
  init_ordered_lattice(N_atoms, N_Cu, lattice);
  Etot = get_Etot(lattice, N_atoms, nearest);
  P = get_order_parameter(lattice, N_Cu);
  r = get_short_range_order_parameter(lattice, nearest, N_Cu);
  
  // start simulation
  for (int iT=0; iT<nT; iT++){ // loop over all temps
    printf("Now running T = %.0f degC\n",T_degC[iT]);    
    beta = 1/(kB*(T_degC[iT] + degC_to_K));
    
    // equilibration run
    if (iT!=0){// First run needs longer equlibration 
      N_eq=N_eq_short; 
    }
    for( int i=0; i<N_eq; i++){
      //take Monte Carlo step.
      MC_step( &Etot, &r, &P, q, lattice, nearest, beta, N_Cu);     
      E_equilibration[i] = Etot;
      P_equilibration[i]= P;
    }
    //Print to file
    if ( ((int)T_degC[iT]) % T_save_step==0){
      write_equil_to_file(T_degC[iT], E_equilibration,N_bonds,P_equilibration,N_eq);
    }
    
    // initialize at temperature[iT]
    E_mean_approx[iT] = Etot; // shift to get higher accuracy in variance
    E_mean[iT]        = 0;
    E_sq_mean[iT]     = 0;
    P_mean[iT]        = 0;
    P_sq_mean[iT]     = 0;
    r_mean[iT]        = 0;
    r_sq_mean[iT]     = 0;
    // production run
    for( int i=0; i<N_timesteps; i++){
      MC_step( &Etot, &r, &P, q, lattice, nearest, beta, N_Cu);
      E_production[i] = Etot- E_mean_approx[iT];
      update_E_P_r(iT, Etot-E_mean_approx[iT], E_mean, E_sq_mean, P, P_mean, 
        P_sq_mean, r, r_mean,r_sq_mean, lattice, nearest, N_Cu);
    }
    E_mean[iT]    *= 1/(double)N_timesteps;
    E_sq_mean[iT] *= 1/(double)N_timesteps; 
    P_mean[iT]    *= 1/(double)N_timesteps;
    P_sq_mean[iT] *= 1/(double)N_timesteps;
    r_mean[iT]    *= 1/(double)N_timesteps;
    r_sq_mean[iT] *= 1/(double)N_timesteps;

    if ( ((int)T_degC[iT]) %T_save_step==0){ // calculate stat inefficiency
      E_Var = E_sq_mean[iT] - E_mean[iT]*E_mean[iT];
      printf("Calculating statistical inefficiencies \n");
      get_phi (phi, N_timesteps, E_mean[iT], E_Var, E_production,N_k,N_skip);
      get_varF_block_average(var_F, N_timesteps, E_mean[iT], E_Var, 
        E_production, N_k, N_skip);
      write_stat_inefficiency_to_file(T_degC[iT], phi, var_F, N_k, N_skip);
    }
  }//END temp for
  
  
  //PRINT TO FILE
  write_production(T_degC, nT, E_mean_approx, E_mean, E_sq_mean, 
    P_mean, P_sq_mean, r_mean, r_sq_mean);
  
  // DON'T FORGET TO FREE ALL malloc's. 
  free(nearest);  nearest = NULL;
  free(lattice);  lattice = NULL;
  free(E_equilibration); E_equilibration = NULL;
  free(P_equilibration); P_equilibration = NULL;
  free(E_mean); E_mean = NULL;
  free(E_mean_approx); E_mean_approx = NULL;
  free(E_sq_mean);  E_sq_mean = NULL;
  free(P_mean);  P_mean = NULL;
  free(P_sq_mean);  P_sq_mean = NULL;
  free(r_mean); r_mean = NULL;
  free(r_sq_mean); r_sq_mean = NULL;
  free(E_production); E_production = NULL;
  free(phi); phi = NULL;
  free(var_F); var_F = NULL;
  free(T_degC); T_degC = NULL;
  
  gsl_rng_free(q); // deallocate  rng
  return 0;
}
