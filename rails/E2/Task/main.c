/*
main.c

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "func.h"
#define PI 3.141592653589
#define N_time 1e7
#define N_between_steps 1000
#define N_part 32 /* The number of particles */


/* Main program */
int main()
{
	/* Declartion of variables */
  	double timestep = 0.1;
	int i,j,k;
	int N_saved = N_time/N_between_steps;
	double timestep_sq,current_time;//timestep_sq = timestep * timestep;
	double m  = 1;
	double kappa = 1;
	double omega0 = sqrt(kappa/m);
	double alpha = 0.1;//0.01;
	double E0 = N_part;  
	
	/* declare file variable */
	FILE *file;

	/* displacement, velocity and acceleration */
	double q[N_part] = {0};
	double p[N_part] = {0};
	double a[N_part] = {0}; 
   double P[N_part] = {0};
   double Q[N_part] = {0};
   double Esum[N_part] = {0};
   
   double trans_matrix[N_part][N_part];
   construct_trans_matrix( N_part, trans_matrix); 

	/* Allocating memory for large vectors */
	/* displacements for writing to file */
	double (*Emat)[N_part] = malloc((N_saved*N_part) * sizeof (double));
	double (*Eavg)[N_part] = malloc((N_saved*N_part) * sizeof (double));
	//double (*E)[2] = malloc( (N_time) * 2 * sizeof (double));
	//double *p_spectrum = malloc((N_time) * sizeof (double)); deallocate!!
	//double *freq = malloc((N_time) * sizeof (double));

	/* Initial conditions */
	/* Set initial displacements and velocites */
	P[0] = sqrt(2*E0);
	get_normal_modes( p, P, N_part, trans_matrix);
	calc_mode_energy(Emat[0], Q, P, N_part, omega0);
	// Qinit = 0 => q[t=0] = 0

	/* Calculate initial accelerations based on initial displacements */
	calc_acc(a, q, m, kappa, alpha, N_part);

	/* timesteps according to velocity Verlet algorithm */
	for (i = 1; i < N_time; i++) {
		/* v(t+dt/2) */
		for (j = 0; j < N_part; j++) {
		    p[j] += timestep * 0.5 * a[j]*m;
		} 

		/* q(t+dt) */
		for (j = 0; j < N_part; j++) {
		    q[j] += timestep * p[j]/m;
		}

		/* a(t+dt) */
		calc_acc(a, q, m, kappa, alpha, N_part);

		/* v(t+dt) */
		for (j = 0; j < N_part; j++) {
		    p[j] += timestep * 0.5 * a[j]*m;
		} 
		
		if (i % N_between_steps == 0){
			/* calculate normal modes */
			k = i/N_between_steps; // number of saved timesteps so far
			get_normal_modes( Q, q, N_part, trans_matrix);
			get_normal_modes( P, p, N_part, trans_matrix);
			calc_mode_energy(Emat[k], Q, P, N_part, omega0);
			
			for (j = 0; j < N_part; j++) {
		   	Esum[j] += Emat[k][j];
		   	Eavg[k][j] = Esum[j]/k;
		   }
		} 
			
		
		/* Save the displacement of the three atoms */
				
		//E[i][0] = calc_pe(q, kappa, alpha, N_part);
		//E[i][1] = calc_ke(p, m, N_part);
		
	}

	/* Print displacement data to output file */
	
	/*
	file = fopen("energy.dat","w");
	for (i = 0; i < N_time; i++) {
		current_time = i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e", current_time, 
			E[i][0] + E[i][1], E[i][0], E[i][1]);	//E[i][0] + E[i][1] is not initiated?
		fprintf(file, "\n");
	}
	fclose(file);	
	*/
	
	file = fopen("Emat.dat","w");
	for (i = 0; i < N_saved; i++) {
		current_time = i * timestep * N_between_steps;
		fprintf(file, "%.4f \t", current_time);
		for (j = 0; j < N_part; j++) { 
			fprintf(file, "%e \t", Emat[i][j]);
		}
		fprintf(file, "\n");
	}
	fclose(file);	
	
	
	file = fopen("Eavg.dat","w");
	for (i = 0; i < N_saved; i++) {
		current_time = i * timestep * N_between_steps;
		fprintf(file, "%.4f \t", current_time);
		for (j = 0; j < N_part; j++) { 
			fprintf(file, "%e \t", Eavg[i][j]);
		}
		fprintf(file, "\n");
	}
	fclose(file);	
		
	
	/* make FFT (powerspectrum) */
	/* Frequency given in inverse time units:
	   [ps^(-1)]=[THz]
	   
	   
	 */
	/*
	double *data_to_plot = q_3;
	powerspectrum(data_to_plot, p_spectrum, N_time);
	powerspectrum_shift(p_spectrum, N_time);
	fft_freq_shift(freq, timestep, N_time);
*/
	/*Save powerspectrum data in file */
 /* 	file = fopen("powerspectrum.dat","w");
	for (i = 0; i < N_time; i++)	{
		fprintf (file,"%e \t %e\n", freq[i], p_spectrum[i]);
	}
	fclose(file);
*/

	/* Free allocated memory */ 
	free(Emat); Emat = NULL;
	free(Eavg); Eavg = NULL;
	//free(E); E = NULL; 
	return 0;    
}
