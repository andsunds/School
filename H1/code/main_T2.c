/*
  main_T2.c, Task 2, H1b
  In this task, we add random noise to the particle positions and see how the 
  system evolves in time. Using the kinetic energy of the particles, we can 
  derive an instantaneous temperature of the system.

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
#define AMU 1.0364e-4
#define kB 8.6173303e-5

/* Main program */
int main()
{
  int N_atoms = 4*N_cells*N_cells*N_cells;
  double m_Al = 27*AMU;
    
  double a_eq = 4.03; // Min potential energy lattice constant
    
  double noise_amplitude = 6.5e-2 * a_eq;
  double t_max=10; //
  double dt = 1e-3;
  int N_timesteps = t_max/dt;
  double t, E_kin;
    
  double (*pos)[3] = malloc(sizeof(double[N_atoms][3]));
  double (*momentum)[3] = malloc(sizeof(double[N_atoms][3]));
  double (*forces)[3] = malloc(sizeof(double[N_atoms][3]));
  double *temperature = malloc(sizeof(double[N_timesteps]));
  double *E_tot = malloc(sizeof(double[N_timesteps]));
    
  FILE *file_pointer;

    
  /* ----------------------------- TASK 2 ----------------------------------*/
  
  init_fcc(pos, N_cells, a_eq); // initialize fcc lattice
  add_noise( N_atoms, 3, pos, noise_amplitude ); // adds random noise to pos
  set_zero( N_atoms, 3, momentum); // set momentum to 0
  get_forces_AL( forces, pos, a_eq*N_cells, N_atoms); //initial cond forces
  
  for (int i=0; i<N_timesteps; i++){
    /* 
       The loop over the timesteps first takes a timestep according to the 
       Verlet algorithm, then calculates the energies and temeperature.
    */
    timestep_Verlet (N_atoms, pos,  momentum, forces, m_Al, dt, a_eq*N_cells);

    E_kin    = get_kin_energy(N_atoms, momentum, m_Al );
    E_tot[i] = (E_kin + get_energy_AL(pos, a_eq*N_cells, N_atoms))*4/N_atoms;
    
    /* 3N*kB*T/2 = 1/(2m) * \sum_{i=1}^{N} p_i^2  = p_sq/(2m) */
    temperature[i] =  E_kin * 2/(3*N_atoms*kB);
  }

  /* Write tempertaure to file */
  char file_name[100];
  sprintf(file_name,"../data/temperature_dt-%0.0e_Task2.tsv", dt);
  file_pointer = fopen(file_name, "w");
  for (int i=0; i<N_timesteps; i++){
    t = i*dt; // time at step i
    fprintf(file_pointer, "%.4f \t %.8f \n", t, temperature[i]);
  }
  fclose(file_pointer);

  /* Write total energy to file */
  sprintf(file_name,"../data/total_energy_dt-%0.0e_Task2.tsv", dt);
  file_pointer = fopen(file_name, "w");
  for (int i=0; i<N_timesteps; i++){
    t = i*dt; // time at step i
    fprintf(file_pointer, "%.4f \t %.8f \n", t, E_tot[i]);
  }
  fclose(file_pointer);
     
  free(pos);         pos = NULL;
  free(momentum);    momentum = NULL;
  free(forces);      forces = NULL;
  free(temperature); temperature = NULL;
  free(E_tot);       E_tot = NULL;
  return 0;
}
