/*
  main_T1.c Task 1 H1b
  In this task, we scan over a range of lattice parameters, a0, to determine 
  which results in the lowest potetntial energy stored in the lattice.
  
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

#include "initfcc.h"
#include "alpotential.h"

#define N_cells 4
#define N_lattice_params 25

/* Main program */
int main()
{   
  int N_atoms = 4*N_cells*N_cells*N_cells;
  double a0;
  double a0_min = 4.0;
  double a0_max = 4.2;
  double da0 = (a0_max - a0_min)/N_lattice_params;

  double (*pos)[3] = malloc(sizeof(double[N_atoms][3]));
  double *energy = malloc(sizeof(double[N_lattice_params]));
  
  FILE *file_pointer;    
    
  for (int i=0; i<N_lattice_params; i++){
    a0 = a0_min + i*da0;        // The lattice constant of this iteration	
    init_fcc(pos, N_cells, a0); // Init, FCC cells with lattice constant `a0`
    // energy per unit cell
    energy[i] = get_energy_AL(pos, N_cells*a0, N_atoms )*4/N_atoms;    
  }

  // Write to files
  file_pointer = fopen("../data/lattice_energies.tsv", "w");
  for (int i=0; i<N_lattice_params; i++){
    a0 = a0_min + i*da0; 
    fprintf(file_pointer, "%.8f \t %.8f \n", a0, energy[i]);
  }
  fclose(file_pointer);    
  
  free(pos);    pos = NULL;
  free(energy); energy = NULL;
  return 0;
}
