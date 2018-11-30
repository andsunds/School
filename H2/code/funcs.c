#include "funcs.h"

double get_bond_E(int site_1, int site_2,
		  double E_ZnZn, double E_CuZn, double E_CuCu){
  double tmp=0;
  switch(site_1 + site_2 ) {
    case 0 :
      //return E_ZnZn;
      tmp=E_ZnZn;
      break;
    case 1 :
      //return E_CuZn;
      tmp= E_CuZn;
      break;
    case 2 :
      //return E_CuCu;
      tmp=E_CuCu;
      break;
  }
  //printf("%0.3f ", tmp);
  return tmp;
}

double get_P(int *lattice, int N_Cu){
  int N_Cu_in_Cu_lattice=0;
  for(int i=0;i<N_Cu;i++){
    N_Cu_in_Cu_lattice+=lattice[i];
  }
  return (double)N_Cu_in_Cu_lattice/N_Cu *2 -1;
}

double get_Etot(int *lattice, int N_atoms, int (*nearest)[8],
		double E_ZnZn, double E_CuZn, double E_CuCu){
  double Etot=0;
  for(int i=0; i<N_atoms; i++){
    for( int j=0; j<8; j++){
      Etot+= get_bond_E(lattice[i], lattice[nearest[i][j]], E_ZnZn, E_CuZn, E_CuCu);
    }
  }
  return Etot/2;
}

void init_nearestneighbor(int Nc, int (*nearest)[8]){
    // create nearest neighbor matrix 
  int i_atom;
  int N_Cu = Nc*Nc*Nc;
  for( int i=0; i<Nc; i++){
    for( int j=0; j<Nc; j++){
      for( int k=0; k<Nc; k++){
        i_atom = k + Nc*j + Nc*Nc*i;
        // k i j in one lattice <=> "k-0.5" "i-0.5" "j-0.5" in the other lattice
        // use mod to handle periodic boundary conditions 
        nearest[i_atom][0] = k        + Nc*j          + Nc*Nc*i           +N_Cu;
        nearest[i_atom][1] = k        + Nc*j          + Nc*Nc*((i+1)%Nc)  +N_Cu;
        nearest[i_atom][2] = k        + Nc*((j+1)%Nc) + Nc*Nc*i           +N_Cu;
        nearest[i_atom][3] = k        + Nc*((j+1)%Nc) + Nc*Nc*((i+1)%Nc)  +N_Cu;
        nearest[i_atom][4] = (k+1)%Nc + Nc*j          + Nc*Nc*i           +N_Cu;
        nearest[i_atom][5] = (k+1)%Nc + Nc*j          + Nc*Nc*((i+1)%Nc)  +N_Cu;
        nearest[i_atom][6] = (k+1)%Nc + Nc*((j+1)%Nc) + Nc*Nc*i           +N_Cu;
        nearest[i_atom][7] = (k+1)%Nc + Nc*((j+1)%Nc) + Nc*Nc*((i+1)%Nc)  +N_Cu;
        
        // k i j in one lattice <=> "k+0.5" "i+0.5" "j+0.5" in the other lattice
        // use mod to handle periodic boundary conditions 
        // note that mod([negative])<0 :/
        i_atom += N_Cu;
        nearest[i_atom][0] = k        + Nc*j          + Nc*Nc*i;
        nearest[i_atom][1] = k        + Nc*j          + Nc*Nc*((i-1+Nc)%Nc);
        nearest[i_atom][2] = k        + Nc*((j-1+Nc)%Nc) + Nc*Nc*i;
        nearest[i_atom][3] = k        + Nc*((j-1+Nc)%Nc) + Nc*Nc*((i-1+Nc)%Nc);
        nearest[i_atom][4] = (k-1+Nc)%Nc + Nc*j          + Nc*Nc*i;
        nearest[i_atom][5] = (k-1+Nc)%Nc + Nc*j          + Nc*Nc*((i-1+Nc)%Nc);
        nearest[i_atom][6] = (k-1+Nc)%Nc + Nc*((j-1+Nc)%Nc) + Nc*Nc*i;
        nearest[i_atom][7] = (k-1+Nc)%Nc + Nc*((j-1+Nc)%Nc) + Nc*Nc*((i-1+Nc)%Nc);
      }
    }
  }
}
