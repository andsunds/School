#include "funcs.h"

/************************ get functions **************************************/
double get_bond_E(int site_1, int site_2){
  double tmp=0;
  switch(site_1 + site_2 ) {
  case 0 :
    //return E_ZnZn;
    tmp=-0.113;
    break;
  case 1 :
    //return E_CuZn;
    tmp= -0.294;
    break;
  case 2 :
    //return E_CuCu;
    tmp=-0.436;
    break;
  }
  return tmp;
}

double get_order_parameter(int *lattice, int N_Cu){
  int N_Cu_in_Cu_lattice=0;
  for(int i=0;i<N_Cu;i++){
    N_Cu_in_Cu_lattice+=lattice[i];
  }
  return (double)N_Cu_in_Cu_lattice/N_Cu *2 -1;
}

double get_short_range_order_parameter(int *lattice, int(*nearest)[N_neigh],
				       int N_Cu){
  int N_CuZnBonds=0;
  for(int i=0;i<N_Cu;i++){
    for( int j=0; j<N_neigh; j++){
      N_CuZnBonds+= (lattice[i] + lattice[nearest[i][j]]) == 1 ;
    }
  }
  return (double) N_CuZnBonds/(4*N_Cu)-1;
}

double get_Etot(int *lattice, int N_atoms, int (*nearest)[N_neigh]){
  double Etot=0;
  for(int i=0; i<N_atoms; i++){
    for( int j=0; j<N_neigh; j++){
      Etot+= get_bond_E(lattice[i], lattice[nearest[i][j]]);
    }
  }
  return Etot/2;
}

void get_phi (double *phi, int N_times, double f_mean,
	      double f_var, double *data, int N_k, int N_skip){
  for (int k=0; k<N_k; k++) {
    phi[k] = 0;
    for (int i=0; (i+k)*N_skip<N_times; i++) {
      phi[k] += data[i*N_skip]*data[(i+k)*N_skip];
    }
    phi[k] = (phi[k]/(N_times/N_skip - k) - f_mean*f_mean)/f_var;
  }
}

void get_varF_block_average(double *var_F, int N_times, double f_mean,
			    double f_var, double *data, int N_k, int N_skip){
  // block average
  int block_size;
  double Fj;
  int number_of_blocks;
  for (int k=0; k<N_k; k++) { // block size loop
    block_size = N_skip * (k+1);
    number_of_blocks = N_times/block_size;
    var_F[k] = 0;
    for (int j=0; j<number_of_blocks; j++) {// loop over all blocks
      Fj = 0;
      for (int i=0; i<block_size; i++) {// internal block loop
        Fj += data[j*block_size + i];
      }
      Fj *= 1/(double)block_size; // these are the values we need the variance of F
      var_F[k] += Fj*Fj; // will become the variance soon
    }
    var_F[k] = var_F[k]/number_of_blocks - f_mean*f_mean;
    var_F[k] *= block_size/f_var;
  }
}

/************** Monte Carlo step functions ************************************/
void MC_step( double *Etot, double *r, double *P, gsl_rng *q,
              int *lattice, int (*nearest)[N_neigh], double beta, int N_Cu){
  /* 
     takes a Monte Carlo step. Updates the lattice and `Etot`
  */
  // Picks two random sites in the whole lattice.
  int i1 = (int)(2*N_Cu*gsl_rng_uniform(q));
  int i2 = (int)(2*N_Cu*gsl_rng_uniform(q));
  // saves the original values
  int old_1 = lattice[i1];
  int old_2 = lattice[i2];
  // Used to clacluate the change in `Etot` and `r`
  double dr = 0;
  double dE = 0;
  // We only need to do something if the two atoms aer different
  if (old_1 != old_2){
    for( int j=0; j<N_neigh; j++){
      /* 
	 The change in `Etot` and `r` are first _minus_ the old energies and `r` 
	 contributtions.
      */
      dE-= get_bond_E(lattice[i1], lattice[nearest[i1][j]])
	  +get_bond_E(lattice[i2], lattice[nearest[i2][j]]);

      dr -= ((lattice[i1] + lattice[nearest[i1][j]]) == 1)
	   +((lattice[i2] + lattice[nearest[i2][j]]) == 1);
    }
    /* Then we do the change of the two atoms */
    lattice[i1] = old_2;
    lattice[i2] = old_1;
    for( int j=0; j<N_neigh; j++){
      /* 
	 And _add_ the contribtions to `Etot` and `r` from the updated lattice.
      */
      dE+= +get_bond_E(lattice[i1], lattice[nearest[i1][j]])
	+get_bond_E(lattice[i2], lattice[nearest[i2][j]]);

      dr += ((lattice[i1] + lattice[nearest[i1][j]]) == 1)
	   +((lattice[i2] + lattice[nearest[i2][j]]) == 1);
    }

    if ( (dE<=0)|| (exp(-beta * dE) >  gsl_rng_uniform(q)) ){
      /* 
	 The test is accepted if dE < 0 (accept immediately), OR
	 otherwise it's accepted with a probability of `exp(-beta * dE)`
      */
      // Updates P
      if (i1 < N_Cu) *P += (double)(lattice[i1] - old_1 )/N_Cu *2;
      if (i2 < N_Cu) *P += (double)(lattice[i2] - old_2 )/N_Cu *2;
    }else{
      /* 
	 If the test failed, we change back to the old lattice configuration 
	 and no change happes to `Etot` or `r`
      */
      lattice[i1] = old_1;
      lattice[i2] = old_2;
      dE = 0;
      dr = 0;
    }// end if step is accepted
    *Etot += dE;
    *r += dr/(4*N_Cu);
  }// end if atoms are different
}

void update_E_P_r(int iT, double E_dev, double *E_mean, double *E_sq_mean,
		  double P, double *P_mean, double *P_sq_mean,
		  double r, double *r_mean, double *r_sq_mean,
		  int *lattice, int (*nearest)[N_neigh], int N_Cu){
  /*
    Updates the macro parameters `E`, `P`, and `r`, as well as their squares.
    Runs in every Monte Carlo step during the producction run.
  */
  E_mean[iT] += E_dev;
  E_sq_mean[iT] += E_dev * E_dev;

  P_mean[iT] += P;
  P_sq_mean[iT] += P*P;

  r_mean[iT] += r;
  r_sq_mean[iT] += r*r;
}

/************************* initializing functions******************************/
void * init_temps( int *nT, double dT_small, double dT_large,
		   double T_start, double T_end, double T_start_fine,
		   double T_end_fine){
  /*
    Creates an array `T_degC` with the temperatures to loop over in the main 
    function, given the fine temperature step range and the sizes of the 
    temperature steps.
  */
  *nT = (int) ((T_end_fine - T_start_fine)/dT_small
	       +(T_start_fine-T_start + T_end-T_end_fine)/dT_large +1);
  double *T_degC = malloc(sizeof(double[*nT]));
  T_degC[0] = T_start;
  for (int iT=1; iT<*nT; iT++){ // loop over all temps
    if (T_degC[iT-1]>=T_start_fine && T_degC[iT-1]<T_end_fine){
      T_degC[iT] = T_degC[iT-1] + dT_small;
    }else{
      T_degC[iT] = T_degC[iT-1] + dT_large;
    }
  }
  return T_degC;
}


void init_ordered_lattice(int N_atoms, int N_Cu, int *lattice){
  /* 
     Initialize lattice with Cu atoms (1) in Cu lattice (i=0:N_Cu-1)
     and Zn (0) in Zn lattice (i=N_cu:N_atoms-1):
  */
  for( int i=0; i<N_Cu; i++){
    lattice[i] = 1;
  }
  for( int i=N_Cu; i<N_atoms; i++){
    lattice[i] = 0;
  }
}

void init_random_lattice(int N_atoms, int N_Cu, int *lattice, gsl_rng *q){
  /*
    Initialize lattice with Cu and Zn atoms randomly distributed:
  */
  for( int i=0; i<N_Cu; i++){
    lattice[i] = (int)(gsl_rng_uniform(q)+0.5);
    lattice[i+N_Cu] = 1-lattice[i];
  }
}


void init_nearestneighbor(int Nc, int (*nearest)[N_neigh]){
  /*
    Create a matrix `nearest[i][j]` with the index of the `j`th neares 
    neighbors to site `i`.
    N.B. Each site has `N_neigh` (8) nearest neighbors.
  */
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
        // note that mod([negative])<0 :
        i_atom += N_Cu;
        nearest[i_atom][0] =k           + Nc*j             + Nc*Nc*i;
        nearest[i_atom][1] =k           + Nc*j             + Nc*Nc*((i-1+Nc)%Nc);
        nearest[i_atom][2] =k           + Nc*((j-1+Nc)%Nc) + Nc*Nc*i;
        nearest[i_atom][3] =k           + Nc*((j-1+Nc)%Nc) + Nc*Nc*((i-1+Nc)%Nc);
        nearest[i_atom][4] =(k-1+Nc)%Nc + Nc*j             + Nc*Nc*i;
        nearest[i_atom][5] =(k-1+Nc)%Nc + Nc*j             + Nc*Nc*((i-1+Nc)%Nc);
        nearest[i_atom][6] =(k-1+Nc)%Nc + Nc*((j-1+Nc)%Nc) + Nc*Nc*i;
        nearest[i_atom][7] =(k-1+Nc)%Nc + Nc*((j-1+Nc)%Nc) + Nc*Nc*((i-1+Nc)%Nc);
      }
    }
  }
}

void* init_random(){
  /*
    Initializes a GSL random nuber generator, and returns the pointer.
  */
  gsl_rng *q;
  const  gsl_rng_type *rng_T;    // static  info  about  rngs
  gsl_rng_env_setup ();          // setup  the  rngs
  rng_T = gsl_rng_default;       // specify  default  rng
  q = gsl_rng_alloc(rng_T);      // allocate  default  rng
  gsl_rng_set(q,time(NULL));     // Initialize  rng
  return q;
}


/************************* file I/O functions *********************************/
void write_equil_to_file(double T_degC, double *E_equilibration, int N_bonds,
			 double *P, int N_eq){
  /*
    Writes the energy per bond `E_equilibration`/`N_bonds` and order 
    parameter `P`, at each Monte Carlo step during the equlibration runs.
  */
  FILE *file_pointer;
  char file_name[256];
  sprintf(file_name,"../data/E_equilibration-T%d.tsv", (int) T_degC);
  file_pointer = fopen(file_name, "w");
  for (int i=0; i<N_eq; i++){
    fprintf(file_pointer, "%.8f\t%.8f \n", E_equilibration[i]/N_bonds,P[i]);
  }
  fclose(file_pointer);
}

void write_production(double *T_degC, int nT, double *E_mean_approx,
		      double *E_mean, double *E_sq_mean,
		      double *P_mean, double *P_sq_mean,
		      double *r_mean, double *r_sq_mean){
  /*
    Writes the macro parameters `E_mean_approx`, `E_mean`, `E_sq_mean`, 
    `P_mean`, `P_sq_mean`, `r_mean`, and `r_sq_mean` for each temperature 
    to file.
  */
  void* init_random();
  FILE *file_pointer;
  char file_name[256];
  sprintf(file_name,"../data/E_production.tsv");
  file_pointer = fopen(file_name, "w");
  fprintf(file_pointer, "%% T[degC]\t E_approx\t<E-E_approx>\t<(E-E_approx)^2>\tP\tr\n");
  for (int iT=0; iT<nT; iT++){
    fprintf(file_pointer, "%.2f\t%.8e\t%.8e\t%.8e\t%.8f\t%.8f\t %.8f\t%.8f \n",
	    T_degC[iT], E_mean_approx[iT], E_mean[iT], E_sq_mean[iT], P_mean[iT],
	    P_sq_mean[iT], r_mean[iT], r_sq_mean[iT]);
  }
  fclose(file_pointer);
}

void write_stat_inefficiency_to_file(double T_degC, double *phi, double *var_F,
				     int N_k, int N_skip){
  /*
    Writes the auto-correlation `phi` and block varaiances `var_F` for each
    tested temperature to file.
  */
  FILE *file_pointer;
  char file_name[256];
  sprintf(file_name,"../data/stat_inefficiency-T%d.tsv", (int) T_degC);
  file_pointer = fopen(file_name, "w");
  for (int i=0; i<N_k; i++){
    fprintf(file_pointer, "%d\t%.8f\t%.8f \n", i*N_skip, phi[i],var_F[i]);
  }
  fclose(file_pointer);
}
