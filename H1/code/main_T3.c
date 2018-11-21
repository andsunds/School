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
  double kappa_Al = 100/(6.6444e+05 * bar); // STRANGE FACTOR 100 OFF !!!
  double a_eq = 4.03;
  double cell_length = a_eq*N_cells;
  double inv_volume = pow(N_cells*cell_length, -3);
  double noise_amplitude = 6.5e-2 * a_eq;

  double T_eq_C= 500;
  double P_eq_bar= 1;
  double T_eq  = T_eq_C + degC_to_K;
  double P_eq  = P_eq_bar*bar;
  double dt    = 5e-3;
  double tau_T = 100*dt;
  double tau_P = 100*dt;
  double t_T_eq= 10*tau_T; //equlibration times
  double t_P_eq= 15*tau_P; //equlibration times
  int N_timesteps_T_eq = t_T_eq/dt;
  int N_timesteps_P_eq = t_P_eq/dt;
  int N_timesteps = N_timesteps_T_eq + N_timesteps_P_eq;
  //int N_save_timesteps = N_timesteps / 100; //for the displacements
  
  double alpha_T, alpha_P,alpha_P_cube_root;
  double t, E_kin, virial;

    
  double (*pos)[3] = malloc(sizeof(double[N_atoms][3]));
  double (*momentum)[3] = malloc(sizeof(double[N_atoms][3]));
  double (*forces)[3] = malloc(sizeof(double[N_atoms][3]));
  double *temperature = malloc(sizeof(double[N_timesteps]));
  double *pressure = malloc(sizeof(double[N_timesteps]));
  //double *displacements = malloc(sizeof(double[N_atoms][N_save_timesteps]));
  
    
  FILE *file_pointer;
    
  /* ----------------------------- TASK 3 ----------------------------------*/
  
  init_fcc(pos, N_cells, a_eq); // initialize fcc lattice
  add_noise( N_atoms, 3, pos, noise_amplitude ); // adds random noise to pos
  set_zero( N_atoms, 3, momentum); // set momentum to 0
  get_forces_AL( forces, pos, cell_length, N_atoms); //initial cond forces

  for (int i=0; i<N_timesteps_T_eq; i++){
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

    /* Equlibrate temperature by scaling momentum by a factor sqrt(alpha_T).
       N.B. It is equally valid to scale the momentum instead of the velocity,
       since they only differ by a constant factor m.
    */
    alpha_T = 1 + 2*dt*(T_eq - temperature[i]) / (tau_T * temperature[i]);
    scale_mat(N_atoms, 3, momentum, sqrt(alpha_T));
    temperature[i]*=alpha_T;
  }

  for (int i=N_timesteps_T_eq; i<N_timesteps; i++){
    /* 
       The loop over the timesteps first takes a timestep according to the 
       Verlet algorithm, then calculates the energies and temeperature.
    */
    timestep_Verlet(N_atoms, pos,  momentum, forces, m_Al, dt, cell_length);
    
    //volume[i] = cell_length;//1/inv_volume;
    
    E_kin  = get_kin_energy(N_atoms, momentum, m_Al );
    virial = get_virial_AL(pos, cell_length, N_atoms);

    /* 3N*kB*T/2 = 1/(2m) * \sum_{i=1}^{N} p_i^2  = p_sq/(2m) */
    temperature[i] =  E_kin * 1/(1.5*N_atoms*kB);
    /* PV = NkT + virial */
    pressure[i] = inv_volume * (1.5*E_kin + virial);

    // Re-equlibrate temperature 
    alpha_T = 1 + 2*dt*(T_eq - temperature[i]) / (tau_T * temperature[i]);
    scale_mat(N_atoms, 3, momentum, sqrt(alpha_T));
    
    // Equlibrate pressure by scaling the posistions by a factor of alpha_P^(1/3)
    
    alpha_P = 1 - kappa_Al* dt*(P_eq - pressure[i])/tau_P;
    alpha_P_cube_root = pow(alpha_P, 1.0/3.0);
    scale_mat(N_atoms, 3, pos, alpha_P_cube_root);

    cell_length*=alpha_P_cube_root;
    inv_volume*=1/alpha_P;

    temperature[i]*=alpha_T;
    pressure[i]*=alpha_P;
  }


  /* Write tempertaure to file */
  sprintf(file_name,"../data/temp-%d_pres-%d_Task3.tsv",
	  (int) T_eq_C, (int) P_eq_bar);
  file_pointer = fopen(file_name, "w");
  for (int i=0; i<N_timesteps; i++){
    t = i*dt; // time at step i
    fprintf(file_pointer, "%.4f \t %.8f \t %.8f \n",
	    t, temperature[i],pressure[i]);
  }
  fclose(file_pointer);

  /* Write phase space coordinates to file */
  sprintf(file_name,"../data/phase-space_temp-%d_pres-%d.tsv",
	  (int) T_eq_C, (int) P_eq_bar);
  file_pointer = fopen(file_name, "w");
  for (int i=0; i<N_atoms; i++){
    for (int j=0;j<3;j++){
      fprintf(file_pointer, " %.16e \t", pos[i][j]);
    }
    for (int j=0;j<3;j++){
      fprintf(file_pointer, " %.16e \t", momentum[i][j]);
    }
    fprintf(file_pointer,"\n");
  }
  fclose(file_pointer);

  /* save equlibrated position and momentum as a binary file */  
  sprintf(file_name,"../data/INIDATA_temp-%d_pres-%d.bin",
	  (int) T_eq_C, (int) P_eq_bar);
  file_pointer = fopen(file_name, "wb");
  fwrite(pos, sizeof(double), 3*N_atoms, file_pointer);
  fwrite(momentum, sizeof(double), 3*N_atoms, file_pointer);
  fwrite(&cell_length, sizeof(double), 1, file_pointer);
  fclose(file_pointer);


  /*  
  printf("T=%0.2f\tP=%0.2e\n",
	 temperature[N_timesteps-1],pressure[N_timesteps-1]);
  */
  
  free(pos); pos = NULL;
  free(momentum); momentum = NULL;
  free(forces); forces = NULL;
  free(temperature); temperature = NULL;
  free(pressure); pressure = NULL;
  //free(volume); volume = NULL;
  return 0;
}
