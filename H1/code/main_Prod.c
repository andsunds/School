/*
 MD_main.c
 
 Created by Anders Lindman on 2013-10-31.
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "initfcc.h"
#include "alpotential.h"
#include "funcs.h"

#define N_cells 4
/* define constants in atomic units: eV, Å, ps, K */
#define AMU 1.0364e-4
#define degC_to_K 273.15
#define bar 6.2415e-07
#define kB 8.61733e-5

/* Main program */
int main()
{
  char file_name[100];
  
  int N_atoms = 4*N_cells*N_cells*N_cells;
  double m_Al = 27*AMU;
  /*
    Values of Young's and shear modulus, Y and G resp., taken from 
    Physics Handbook, table T 1.1. Bulk mudulus then calculated as
    B = Y*G / (9*G - 3*Y)   [F 1.15, Physics Handbook]
    kappa = 1/B
  */
//  double kappa_Al = 100/(6.6444e+05 * bar); // STRANGE FACTOR 100 OFF !!!
  double cell_length = 0;
  double inv_volume;
  

  double T_eq_C   = 500;
  double P_eq_bar = 1;
//  double T_eq     = T_eq_C + degC_to_K;
//  double P_eq     = P_eq_bar*bar;
  double dt       = 5e-3;
  double t_end    = 10;
//  double tau_T = 100*dt;
//  double tau_P = 100*dt; 
  
  int N_timesteps = t_end/dt;
  int N_save_timesteps = N_timesteps / 100; //for the displacements
  
//  double alpha_T, alpha_P,alpha_P_cube_root;  
  double t, E_kin, virial;
    
  double (*pos)[3]      = malloc(sizeof(double[N_atoms][3]));
  double (*pos_0)[3]    = malloc(sizeof(double[N_atoms][3]));
  double (*momentum)[3] = malloc(sizeof(double[N_atoms][3]));
  double (*mom_0)[3]    = malloc(sizeof(double[N_atoms][3]));
  double (*forces)[3]   = malloc(sizeof(double[N_atoms][3]));
  double (*displacements)[N_atoms] = malloc(sizeof(double[N_save_timesteps][N_atoms]));
  double *temperature   = malloc(sizeof(double[N_timesteps]));
  double *pressure      = malloc(sizeof(double[N_timesteps]));
  //double *msd           = malloc(sizeof(double[N_timesteps]));

  
    
  FILE *file_pointer;
    
  /* ----------------------------- TASK 3 ----------------------------------*/
  
  // Read in the positions
  sprintf(file_name,"../data/pos_temp-%d_pres-%d.bin",
	  (int) T_eq_C, (int) P_eq_bar);
  file_pointer = fopen(file_name, "rb");
  for (int i=0; i<N_atoms; i++){
    fread(pos[i], sizeof(double), 3, file_pointer);
    for (int j=0; j<3; j++){
      pos_0[i][j]=pos[i][j];
    }
  }
  fclose(file_pointer);
  // Read in the momentums
  sprintf(file_name,"../data/mom_temp-%d_pres-%d.bin",
	  (int) T_eq_C, (int) P_eq_bar);
  file_pointer = fopen(file_name, "rb");
  for (int i=0; i<N_atoms; i++){
    fread(momentum[i], sizeof(double), 3, file_pointer);
    for (int j=0; j<3; j++){
      mom_0[i][j]=momentum[i][j];
    }
  }
  fclose(file_pointer);
  
  // read cell_length
  sprintf(file_name,"../data/cell-length_temp-%d_pres-%d.bin",
	  (int) T_eq_C, (int) P_eq_bar);
  file_pointer = fopen(file_name, "rb");
  fread(&cell_length, sizeof(double), 1, file_pointer);
  fclose(file_pointer);
  inv_volume = pow(N_cells*cell_length, -3);
  
  
  get_forces_AL( forces, pos, cell_length, N_atoms); //initial cond forces
  
  for (int i=0; i<N_timesteps; i++){
    /* 
       The loop over the timesteps first takes a timestep according to the 
       Verlet algorithm, then calculates the energies and temeperature.
    */
    timestep_Verlet(N_atoms, pos,  momentum, forces, m_Al, dt, cell_length);
    
    E_kin  = get_kin_energy(N_atoms, momentum, m_Al );
    virial = get_virial_AL(pos, cell_length, N_atoms);
    
    /* PV = NkT + virial */
    pressure[i] = inv_volume * (1.5*E_kin + virial);
    /* 3N*kB*T/2 = 1/(2m) * \sum_{i=1}^{N} p_i^2  = p_sq/(2m) */
    temperature[i] =  E_kin * 1/(1.5*N_atoms*kB);
    

    /*
    // Re-equlibrate temperature 
    alpha_T = 1 + 2*dt*(T_eq - temperature[i]) / (tau_T * temperature[i]);
    scale_mat(N_atoms, 3, momentum, sqrt(alpha_T));
    
    // Equlibrate pressure by scaling the posistions by a factor of alpha_P^(1/3)
    alpha_P = 1 - kappa_Al* dt*(P_eq - pressure[i])/tau_P;
    alpha_P_cube_root = pow(alpha_P, 1.0/3.0);
    scale_mat(N_atoms, 3, pos, alpha_P_cube_root);

    cell_length*=alpha_P_cube_root;
    inv_volume*=1/alpha_P;
    pressure[i] *= alpha_P;
    temperature[i] *= alpha_T;
    
    //*/
  }
 

  /* Write tempertaure to file */
  //*
  sprintf(file_name,"../data/temp-%d_pres-%d_Prod-test.tsv",
	  (int) T_eq_C, (int) P_eq_bar);
  file_pointer = fopen(file_name, "w");
  for (int i=0; i<N_timesteps; i++){
    t = i*dt; // time at step i
    fprintf(file_pointer, "%.4f \t %.8f \t %.8f \n",
	    t, temperature[i],pressure[i]);
  }
  fclose(file_pointer);
  //*/
       
  free(pos);           pos = NULL;
  free(momentum);      momentum = NULL;
  free(forces);        forces = NULL;
  free(temperature);   temperature = NULL;
  free(pressure);      pressure = NULL;
  free(displacements); displacements = NULL;
  return 0;
}
