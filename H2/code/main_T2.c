/*
  main_T1.c Task 1 H1b
  In this task, we scan over a range of lattice parameters, a0, to determine 
  which one results in the lowest potential energy stored in the lattice.
  
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

#include "funcs.h"

#define Nc 10 //number of cells
#define N_neigh 8
#define degC_to_K 273.15
#define kB 8.61733e-5

/* Main program */
int main()
{   
  int N_Cu = Nc*Nc*Nc;
  int N_CuCu, N_ZnZn, N_CuZn;
  double E_CuCu = -0.436;
  double E_ZnZn = -0.113;
  double E_CuZn = -0.294;
  
  double Etot, dE;
  double T_min_degC = 400;
  double dT = 10; 
  int nT = 20;
  double beta,T_degC; 
  
  int N_atoms = 2*N_Cu;
  int N_bonds = 8*N_Cu;
  int N_timesteps = 1e5; 
  int N_eq = 1e5; 
  
  int test_1, test_2, i1,i2;
  int (*nearest)[N_neigh] = malloc(sizeof(double[N_atoms][N_neigh])); // nearest neighbors
  int *lattice = malloc(sizeof(double[N_atoms]));
  double *E_equilibration = malloc(sizeof(double[N_eq]));
    
  
  double *E_mean = malloc(sizeof(double[nT]));
  double *E_sq_mean = malloc(sizeof(double[nT]));
  
  
  const  gsl_rng_type *rng_T;    // static  info  about  rngs 
  gsl_rng *q;                // rng  instance 
  gsl_rng_env_setup ();      // setup  the  rngs 
  rng_T = gsl_rng_default;       // specify  default  rng 
  q = gsl_rng_alloc(rng_T);      // allocate  default  rng 
  gsl_rng_set(q,time(NULL)); // Initialize  rng 
  
  FILE *file_pointer;
  char file_name[100];
  
  
  init_nearestneighbor(Nc, nearest);

  // initialize lattice with Cu atoms (1) in Cu lattice and Zn (0) in Zn lattice
  for( int i=0; i<N_Cu; i++){
    lattice[i] = 1; 
  }
  for( int i=N_Cu; i<N_atoms; i++){
    lattice[i] = 0; 
  }
  
  Etot = N_bonds * E_CuZn;
  
  for (int iT=0; iT<nT; iT++){
  T_degC = T_min_degC + dT*iT;
  beta = 1/(kB*(T_degC+ degC_to_K));
  E_mean[iT] = 0;
  E_sq_mean[iT] = 0;
  // equilibration run
  for( int i=0; i<N_eq; i++){
    i1 = (int)(N_atoms*gsl_rng_uniform(q)); 
    i2 = (int)(N_atoms*gsl_rng_uniform(q)); 
    test_1 = lattice[i2];
    test_2 = lattice[i1];
     
    dE = 0;
    if (test_1 != test_2){
      for( int j=0; j<N_neigh; j++){
        dE+= get_bond_E(test_1, lattice[nearest[i1][j]], E_ZnZn, E_CuZn, E_CuCu) 
            +get_bond_E(test_2, lattice[nearest[i2][j]], E_ZnZn, E_CuZn, E_CuCu)
            -get_bond_E(test_1, lattice[nearest[i2][j]], E_ZnZn, E_CuZn, E_CuCu)
           -get_bond_E(test_2, lattice[nearest[i1][j]], E_ZnZn, E_CuZn, E_CuCu);    
      }
      
      if ( (dE<=0)|| ( exp(-beta * dE) >  gsl_rng_uniform(q) ) ){
        lattice[i1] = test_1;
        lattice[i2] = test_2;
      }else{
        dE = 0;
      }
    }
    Etot += dE;
    E_equilibration[i] = Etot;
  }

  // production run
  for( int i=0; i<N_timesteps; i++){
    i1 = (int)(N_atoms*gsl_rng_uniform(q)); 
    i2 = (int)(N_atoms*gsl_rng_uniform(q)); 
    
    test_1 = lattice[i2];
    test_2 = lattice[i1];
    dE = 0;
    
    if (test_1 != test_2){
      for( int j=0; j<N_neigh; j++){
        dE+= get_bond_E(test_1, lattice[nearest[i1][j]], E_ZnZn, E_CuZn, E_CuCu) 
            +get_bond_E(test_2, lattice[nearest[i2][j]], E_ZnZn, E_CuZn, E_CuCu)
            -get_bond_E(test_1, lattice[nearest[i2][j]], E_ZnZn, E_CuZn, E_CuCu)
           -get_bond_E(test_2, lattice[nearest[i1][j]], E_ZnZn, E_CuZn, E_CuCu);
      }
      
      if (dE<0 || ( exp(-beta * dE) >  gsl_rng_uniform(q)) ){
        lattice[i1] = test_1;
        lattice[i2] = test_2;
      }else{
        dE = 0;
      }
    }
    Etot += dE;
    E_mean[iT] += Etot;
    E_sq_mean[iT] += Etot * Etot;
  }
  E_mean[iT] *= 1/(double)N_timesteps;
  E_sq_mean[iT] *= 1/(double)N_timesteps;
  
  
  
  //PRINT TO FILE
  
  sprintf(file_name,"../data/E_equilibration-T%d.tsv",
	  (int) T_degC);    
  file_pointer = fopen(file_name, "w");
  for (int i=0; i<N_eq; i++){
    fprintf(file_pointer, "%.8f \n", E_equilibration[i]/N_bonds);
  }
  fclose(file_pointer);    
  
  }
  
  //PRINT TO FILE
  sprintf(file_name,"../data/E_mean.tsv");    
  file_pointer = fopen(file_name, "w");
  for (int iT=0; iT<nT; iT++){
    T_degC = T_min_degC + dT*iT;
    fprintf(file_pointer, "%.2f \t %.8f \t %.8f \n", T_degC, E_mean[iT], E_sq_mean[iT]);
  }
  fclose(file_pointer); 
  
  
  free(nearest);  nearest = NULL;
  free(lattice); lattice = NULL;
  free(E_equilibration); E_equilibration = NULL;
  gsl_rng_free(q); // deallocate  rng
  return 0;
}
