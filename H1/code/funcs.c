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
		    mat[i][j] += noise_amplitude * (2*gsl_rng_uniform(q)-1); /*  generate  random  number 
		(repeatable) */
		}
	}
	gsl_rng_free(q); /*  deallocate  rng */
} 

void timestep_Verlet ( int N_atoms, double (*pos)[3],  double (*momentum)[3],
    double (*forces)[3], double m, double dt, double cell_length){
    
    
	for (int i = 0; i < N_atoms; i++) {
	    for (int j = 0; j < 3; j++) {
	        /* p(t+dt/2) */
            momentum[i][j] += dt * 0.5 * forces[i][j];
             /* q(t+dt) */
            pos[i][j] += dt * momentum[i][j] / m;
        } 
    }

	/* a(t+dt) */
	get_forces_AL( forces, pos, cell_length, N_atoms);

	/* v(t+dt) */
	for (int i = 0; i < N_atoms; i++) {
	    for (int j = 0; j < 3; j++) {
	        /* p(t+dt/2) */
            momentum[i][j] += dt * 0.5 * forces[i][j];
        } 
    }
}

void set_zero (int M, int N, double mat[M][N]){
    for (int i = 0; i < M; i++) {
	    for (int j = 0; j < N; j++) {
	        mat[i][j] = 0;
        } 
    }
} 
