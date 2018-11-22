#include "funcs.h"

void add_noise(int M, int N, double mat[M][N], double noise_amplitude )
{
  const  gsl_rng_type *T; /*  static  info  about  rngs */
  gsl_rng *q; /* rng  instance  */
  gsl_rng_env_setup (); /*  setup  the  rngs */
  T = gsl_rng_default; /*  specify  default  rng */
  q = gsl_rng_alloc(T); /*  allocate  default  rng */
  gsl_rng_set(q,time(NULL)); /*  Initialize  rng */
	
  for (int i=0; i<N; i++){
    for (int j=0; j<M; j++){
      // adds uniformly distributed random noise in range +-`noise_amplitude`
      mat[i][j] += noise_amplitude * (2*gsl_rng_uniform(q)-1); 
    }
  }
  gsl_rng_free(q); /*  deallocate  rng */
} 

void timestep_Verlet ( int N_atoms, double (*pos)[3],  double (*momentum)[3],
		       double (*forces)[3], double m, double dt, 
		       double cell_length){
  for (int i = 0; i < N_atoms; i++) {
    for (int j = 0; j < 3; j++) {
      /* p(t+dt/2) */
      momentum[i][j] += dt * 0.5 * forces[i][j];
      /* q(t+dt) */
      pos[i][j] += dt * momentum[i][j] / m;
    } 
  }
  /* F(t+dt) */
  get_forces_AL( forces, pos, cell_length, N_atoms);
  for (int i = 0; i < N_atoms; i++) {
    for (int j = 0; j < 3; j++) {
      /* p(t+dt/2) */
      momentum[i][j] += dt * 0.5 * forces[i][j];
    } 
  }
}

double get_kin_energy ( int N_atoms,  double (*momentum)[3], double m ) {
  double p_sq=0; // momentum squared
  for (int i = 0; i < N_atoms; i++) {
    for (int j = 0; j < 3; j++) {
      p_sq += momentum[i][j] * momentum[i][j];
    } 
  }
  return p_sq / (2*m); 
}

void get_displacements ( int N_atoms,  double (*positions)[3],
			 double (*initial_positions)[3], double disp[]) {
  for (int i = 0; i < N_atoms; i++) {
    for (int j = 0; j < 3; j++) {
      disp[i] += (positions[i][j] - initial_positions[i][j])
	        *(positions[i][j] - initial_positions[i][j]);
    }
    disp[i] = sqrt(disp[i]);
  }
}
								   

void get_MSD ( int N_atoms,  int N_times, double all_pos[N_times][N_atoms][3], 
                 double MSD[N_times]) {
   /* all_pos = positions of all particles at all (saved) times */
   /* outer time index it starts at outer it = 1, since MSD[0] = 0*/        
   for (int it = 1; it < N_times; it++) { // 	 
      for (int jt = 0; jt < N_times-it; jt++) { // summed time index
         for (int kn = 0; kn < N_atoms; kn++) { // particle index
            for (int kd = 0; kd < 3; kd++) { // three dimensions
               MSD[it] += (all_pos[it+jt][kn][kd] - all_pos[jt][kn][kd])
	                      *(all_pos[it+jt][kn][kd] - all_pos[jt][kn][kd]);
	         }
	      }
      }
      MSD[it] *= 1/( (double)N_atoms * (N_times-it));
   } 
}

void get_vel_corr ( int N_atoms,  int N_times, double all_vel[N_times][N_atoms][3], 
                 double vel_corr[N_times]) {
   /* all_vel = velocity of all particles at all (saved) times */        
   for (int it = 0; it < N_times; it++) { // 	 
      for (int jt = 0; jt < N_times-it; jt++) { // summed time index
         for (int kn = 0; kn < N_atoms; kn++) { // particle index
            for (int kd = 0; kd < 3; kd++) { // three dimensions
               vel_corr[it] += (all_vel[it+jt][kn][kd] * all_vel[jt][kn][kd]);
	         }
	      }
      }
      vel_corr[it] *= 1/( (double)N_atoms * (N_times-it));
   } 
}


void copy_mat (int M, int N, double mat_from[M][N], double mat_to[M][N]){
   /* Copies matrix `mat_from` to `mat_to` */
   for (int i = 0; i < M; i++) {
      for (int j = 0; j < N; j++) {
      mat_to[i][j] = mat_from[i][j];
      } 
   }
}

void set_zero (int M, int N, double mat[M][N]){
  /* Sets the matrix `mat` to zero */
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      mat[i][j] = 0;
    } 
  }
} 

void scale_mat (int M, int N, double mat[M][N], double alpha){
  /* Scales the matrix `mat` by factor `alpha` */
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      mat[i][j] *= alpha;
    } 
  }
}
