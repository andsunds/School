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
#define N_lattice_params 20
#define N_timesteps 10000
#define AMU 1.0364e-4

/* Main program */
int main()
{
    
    int N_atoms = 4*N_cells*N_cells*N_cells;
    double m_Al = 27*AMU;
    
    double a0;
    double a0_min = 4.0;
    double a0_max = 4.2;
    double da0 = (a0_max - a0_min)/N_lattice_params;
    double a_eq = 4.03; 
    
    double noise_amplitude = 6.5e-2;
    double dt = 2e-3; 
    double t;
    
    double (*pos)[3] = malloc(sizeof(double[N_atoms][3]));
    double (*momentum)[3] = malloc(sizeof(double[N_atoms][3]));
    double (*forces)[3] = malloc(sizeof(double[N_atoms][3]));
    double *energy = malloc(sizeof(double[N_lattice_params]));
    
    
	FILE *file_pointer;
	 
    
    
    /*
     Code for generating a uniform random number between 0 and 1. srand should only
     be called once.
    */
    /*
     srand(time(NULL));
     double random_value;
     random_value = (double) rand() / (double) RAND_MAX;
    */
    
    /*
     Descriptions of the different functions in the files initfcc.c and
     alpotential.c are listed below.
    */
    
    /* 
     Function that generates a fcc lattice in units of [Å]. Nc is the number of 
     primitive cells in each direction and a0 is the lattice parameter. The
     positions of all the atoms are stored in pos which should be a matrix of the
     size N x 3, where N is the number of atoms. The first, second and third column
     correspond to the x,y and z coordinate respectively.
    */
    
    
    
    /* 
     Function that calculates the potential energy in units of [eV]. pos should be
     a matrix containing the positions of all the atoms, L is the length of the 
     supercell and N is the number of atoms.
    */
    
    /* ----------------------------- TASK 1 ----------------------------------*/    
    
	for (int i=0; i<N_lattice_params; i++){
		a0 = a0_min + i*da0; 	
 		init_fcc(pos, N_cells, a0);
 		
 		// energy per unit cell 	
 		energy[i] = get_energy_AL(pos, N_cells*a0, N_atoms )*4/N_atoms;    
 	}

    file_pointer = fopen("../data/lattice_energies.tsv", "w");
	for (int i=0; i<N_lattice_params; i++){
		a0 = a0_min + i*da0; 
 		fprintf(file_pointer, "%.8f \t %.8f \n", a0, energy[i]);
 	}
    fclose(file_pointer);    
    
    /* ----------------------------- TASK 2 ----------------------------------*/
    
    energy = realloc(energy, sizeof(double[N_timesteps]));
    init_fcc(pos, N_cells, a_eq);
    add_noise( N_atoms, 3, pos, noise_amplitude );
    set_zero( N_atoms, 3, momentum);
    get_forces_AL( forces, pos, a_eq*N_cells, N_atoms);
    
    for (int i=0; i<N_timesteps; i++){
        timestep_Verlet (N_atoms, pos,  momentum, forces, m_Al, dt, a_eq*N_cells);
        energy[i] = get_energy_AL(pos, N_cells*a_eq, N_atoms )/N_atoms;
    }
    
    file_pointer = fopen("../data/atom_energies_tmp.tsv", "w");
	for (int i=0; i<N_lattice_params; i++){
		t = i*dt; 
 		fprintf(file_pointer, "%.4f \t %.8f \n", t, energy[i]);
 	}
    fclose(file_pointer);
    
    /* 
     Function that calculates the virial in units of [eV]. pos should be a matrix
     containing the positions of all the atoms, L is the length of the supercell 
     and N is the number of atoms.
    */
    /*
     double virial;
     virial = get_virial_AL(pos, L, N);
    */
    
    /*
     Function that calculates the forces on all atoms in units of [eV/Å]. the 
     forces are stored in f which should be a matrix of size N x 3, where N is the
     number of atoms and column 1,2 and 3 correspond to the x,y and z component of
     the force resepctively . pos should be a matrix containing the positions of 
     all the atoms, L is the length of the supercell and N is the number of atoms.
    */
    /*
     get_forces_AL(f,pos, L, N);
    */
    
    free(pos); pos = NULL;
    free(momentum); momentum = NULL;
    free(forces); forces = NULL;
    free(energy); energy = NULL;
    return 0;
}
