#ifndef  _FUNCS_H
#define  _FUNCS_H

//  includes
#include  <math.h>
#include <gsl/gsl_rng.h>
#include <time.h>  // also needed for random stuff

//  prototypes

double get_bond_E(int site_1, int site_2, double E_ZnZn, double E_CuZn, 
  double E_CuCu);
  
void init_nearestneighbor(int Nc, int (*nearest)[8]);

double get_P(int *lattice, int N_Cu);

double get_Etot(int *lattice, int N_atoms, int (*nearest)[8],
		double E_ZnZn, double E_CuZn, double E_CuCu);
#endif
