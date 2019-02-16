/*
  main_Prod.c, Production runs, H1b
  In this program, we use the equlibrated micro-states from Tasks 3-4 to study 
  dynamical properties, such as mean squared displacement (MSD), velocity 
  auto-correlation function, and the power spectral density of the atom 
  movements.
 
  System of units:
  Energy   - eV
  Time     - ps
  Length   - Angstrom
  Temp     - K
  Mass     - eV (ps)^2 A^(-2)
  Pressure - eV A^(-3) 
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
  double cell_length;
  double inv_volume;
  
  double T_eq_C   = 500;
  double P_eq_bar = 1;

  double dt       = 5e-4; // higher res for spectral function
  double t_end    = 5;
  int N_timesteps = t_end/dt;  
  int N_between_steps = 1; // save all steps for max res in spectral function
  int N_save_timesteps = N_timesteps / N_between_steps; 
  int N_save_atoms = 5;
  
  double t, E_kin, virial;
    
  double (*pos)[3]      = malloc(sizeof(double[N_atoms][3]));
  double (*pos_0)[3]    = malloc(sizeof(double[N_atoms][3]));//for displacements
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
  double *pow_spec      = malloc(sizeof(double[N_save_timesteps])); 
  double *freq		= malloc(sizeof(double[N_save_timesteps])); 

  // Initialize to 0
  for (int i = 0; i<N_save_timesteps; i++){
    msd[i] = 0;
    pow_spec[i] = 0;
    vel_corr[i] = 0;
  }
  FILE *file_pointer;
  
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
      // Saves the displacements of some atoms into `displacements`
      get_displacements (N_save_atoms,  pos, pos_0, displacements[k]);

      // Saves all the positions
      copy_mat(N_atoms, 3, pos, pos_all[k]);
      
      // Saves all the velocities
      copy_mat(N_atoms, 3, momentum, vel_all[k]);
      //But we need to scale the momenta to get the velocities
      scale_mat(N_atoms, 3, vel_all[k], 1/m_Al);
    }
    
    if ((i*10) % N_timesteps == 0){ //Print out progress at every 10%
      printf("done %d0%% of Verlet timestepping\n", (i*10)/N_timesteps);
    }
  }
  printf("done 100%% of Verlet timestepping\n");
  
  //Calculating MSD
  printf("calculating MSD\n");
  get_MSD(N_atoms, N_save_timesteps, pos_all, msd);
  
  //Calculating the velocity correlation function
  printf("calculating velocity correlation\n");
  get_vel_corr(N_atoms, N_save_timesteps, vel_all, vel_corr);
  
  //Calculating the velocity power spectrum
  printf("calculating power spectrum\n");
  get_powerspectrum(N_atoms, N_save_timesteps, vel_all, pow_spec);
  fft_freq(freq, dt, N_save_timesteps);

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
  
  /* Write power spectrum to file */
  sprintf(file_name,"../data/temp-%d_pres-%d_power-spectrum.tsv",
	  (int) T_eq_C, (int) P_eq_bar);
  file_pointer = fopen(file_name, "w");
  // write header 
  fprintf(file_pointer, "%% f[1/ps] \t P[A/ps]^2 \n"); 
  for (int i=0; i<N_save_timesteps/2; i++){ // only save from f=0 to f_crit
    fprintf(file_pointer, "%.4f \t %.8f \n", freq[i], pow_spec[i]);
  }
  fclose(file_pointer);

  // Freeing all the memory
  free(pos);           pos = NULL;
  free(pos_0);         pos_0 = NULL;
  free(momentum);      momentum = NULL;
  free(forces);        forces = NULL;
  free(temperature);   temperature = NULL;
  free(pressure);      pressure = NULL;
  free(displacements); displacements = NULL;
  free(pos_all);       pos_all = NULL;
  free(vel_all);       vel_all = NULL;
  free(msd);           msd = NULL;
  free(vel_corr);      vel_corr = NULL;
  free(pow_spec);      pow_spec = NULL;
  free(freq);          freq = NULL;
  
  return 0;
}
