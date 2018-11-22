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
  

  double T_eq_C   = 700;
  double P_eq_bar = 1;
//  double T_eq     = T_eq_C + degC_to_K;
//  double P_eq     = P_eq_bar*bar;
  double dt       = 5e-4; // higher res for spectral function
  double t_end    = 5;
//  double tau_T = 100*dt;
//  double tau_P = 100*dt; 
  
  int N_timesteps = t_end/dt;
  
  int N_between_steps = 1;
  int N_save_timesteps = N_timesteps / N_between_steps; //for the displacements
  int N_save_atoms = 5;
  
//  double alpha_T, alpha_P,alpha_P_cube_root;  
  double t, E_kin, virial;
    
  double (*pos)[3]      = malloc(sizeof(double[N_atoms][3]));
  double (*pos_0)[3]    = malloc(sizeof(double[N_atoms][3]));
  double (*momentum)[3] = malloc(sizeof(double[N_atoms][3]));
  double (*forces)[3]   = malloc(sizeof(double[N_atoms][3]));
  double (*displacements)[N_save_atoms] = 
  			malloc(sizeof(double[N_save_timesteps][N_save_atoms]));
  double (*pos_all)[N_atoms][3] = 
  			malloc(sizeof(double[N_save_timesteps][N_atoms][3]));
  double (*vel_all)[N_atoms][3] = 
  			malloc(sizeof(double[N_save_timesteps][N_atoms][3]));
  double *temperature   = malloc(sizeof(double[N_timesteps]));
  double *pressure      = malloc(sizeof(double[N_timesteps]));
  double *msd           = malloc(sizeof(double[N_save_timesteps])); 
  double *vel_corr      = malloc(sizeof(double[N_save_timesteps])); 
  
  for (int i = 0; i<N_save_timesteps; i++){
    msd[i] = 0;
  }
  FILE *file_pointer;
    
  /* ----------------------------- TASK 3 ----------------------------------*/
  
  // read positions, momenta and cell_length
  sprintf(file_name,"../data/INIDATA_temp-%d_pres-%d.bin",
	  (int) T_eq_C, (int) P_eq_bar);
  file_pointer = fopen(file_name, "rb");
  fread(pos, sizeof(double), 3*N_atoms, file_pointer);
  fread(momentum, sizeof(double), 3*N_atoms, file_pointer);
  fread(&cell_length, sizeof(double), 1, file_pointer);
  fclose(file_pointer);
  
  for (int i=0; i<N_atoms; i++){
    for (int j=0; j<3; j++){
      pos_0[i][j]=pos[i][j];
    }
  }
  inv_volume = pow(N_cells*cell_length, -3);
  get_forces_AL( forces, pos, cell_length, N_atoms); //initial cond forces
  
  printf("Initialized. Starting with Verlet timestepping.\n");
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
    
    if (i % N_between_steps == 0){
	     int k = i/N_between_steps; // number of saved timesteps so far
		  get_displacements (N_save_atoms,  pos, pos_0, displacements[k]);
	     copy_mat(N_atoms, 3, pos, pos_all[k]);
	     
	     copy_mat(N_atoms, 3, momentum, vel_all[k]);
	     scale_mat(N_atoms, 3, vel_all[k], 1/m_Al);
	 }
    if ((i*10) % N_timesteps == 0){
       printf("done %d0 %% of Verlet timestepping\n", (i*10)/N_timesteps);
    }
  }
  printf("calculating MSD\n");
  get_MSD(N_atoms, N_save_timesteps, pos_all, msd);
  printf("calculating velocity correlation\n");
  get_vel_corr(N_atoms, N_save_timesteps, vel_all, vel_corr);

  printf("writing to file\n");
  /* Write tempertaure to file */
  
  sprintf(file_name,"../data/temp-%d_pres-%d_Prod-test.tsv",
	  (int) T_eq_C, (int) P_eq_bar);
  file_pointer = fopen(file_name, "w");
  for (int i=0; i<N_timesteps; i++){
    t = i*dt; // time at step i
    fprintf(file_pointer, "%.4f \t %.8f \t %.8f \n",
	    t, temperature[i],pressure[i]);
  }
  fclose(file_pointer);
  
  /* Write displacements to file */
  sprintf(file_name,"../data/temp-%d_pres-%d_displacements.tsv",
	  (int) T_eq_C, (int) P_eq_bar);
  file_pointer = fopen(file_name, "w");
  for (int i=0; i<N_save_timesteps; i++){
    t = i*dt*N_between_steps; // time at step i
    fprintf(file_pointer, "%.4f", t);
    for (int j=0; j<N_save_atoms; j++){
    	fprintf(file_pointer, "\t %.8f", displacements[i][j]);
	 }
	 fprintf(file_pointer, "\n");
  }
  fclose(file_pointer);
  
   /* Write MSD to file */
  sprintf(file_name,"../data/temp-%d_pres-%d_dynamicProperties.tsv",
	  (int) T_eq_C, (int) P_eq_bar);
  file_pointer = fopen(file_name, "w");
  // write header 
  fprintf(file_pointer, "%% t[ps] \t MSD[A^2] \t vel_corr [A/ps]^2 \n"); 
  for (int i=0; i<N_save_timesteps; i++){
     t = i*dt*N_between_steps; // time at step i
     fprintf(file_pointer, "%.4f \t %.8f \t %.8f \n", t, msd[i], vel_corr[i]);
  }
  fclose(file_pointer);
       
  free(pos);           pos = NULL;
  free(pos_0);         pos_0 = NULL;
  free(momentum);      momentum = NULL;
  free(forces);        forces = NULL;
  free(temperature);   temperature = NULL;
  free(pressure);      pressure = NULL;
  free(displacements); displacements = NULL;
  free(pos_all); pos_all = NULL;
  free(vel_all); vel_all = NULL;
  return 0;
}
