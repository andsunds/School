/*
E1_main.c

NIST:
* σg+	1	Sym str	1333	 C	 ia		1388.15	gas	FR(2ν2)
* σg+	1	Sym str	1333	 C	 ia		1285.40	gas	FR(2ν2)
  πu	2	Bend	667	 A	667.38 S	gas	 ia		
* σu+	3	Anti str	2349	 A	2349.16 VS	gas	 ia		
*=our modes
This corresponds to:
   39.962 (sym)    70.421 (anti-sym) THz



MATHEMATICA for calcualting eigen values:

M2 = {{1, -1,  0}, {-1 mRatio , 2 mRatio, -1 mRatio}, {0, -1,  1}};
Simplify[Eigenvectors[M2] /. mRatio -> (16/12)] // TraditionalForm
Simplify[Eigenvalues[M2] /. mRatio -> (16/12)] // TraditionalForm

{{-1, 0, 1}, {1, 1, 1}, {1, -(8/3), 1}}

{1,0,11/3}


Theoretical base frequencies: 39.058747 (sym)   74.791807 (anti-sym) THz
 
Created by AL on 2013-10-24.
Further developed by Martin Gren on 2014-10-20.
*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "func.h"
#include "fft_func.h"
#define PI 3.141592653589
#define N_time 32768 /* N_time+1 = power of 2, for best speed */
#define nbr_of_particles 3 /* The number of particles is 3 */
#define AMU 1.0364e-4
#define kappa_unit 6.242e-2


/* Main program */
int main()
{
	/* Declartion of variables */
  	double timestep;
	int i,j;
	double timestep_sq,current_time;
	double m[nbr_of_particles]={16*AMU,12*AMU,16*AMU};
	double kappa;

	/* declare file variable */
	FILE *file;

	/* displacement, velocity and acceleration */
	double q[nbr_of_particles];
	double v[nbr_of_particles];
	double a[nbr_of_particles]; 


	/* Allocating memory for large vectors */
	/* displacements for writing to file */
	double *q_1 = malloc((N_time) * sizeof (double));
	double *q_2 = malloc((N_time) * sizeof (double));
	double *q_3 = malloc((N_time) * sizeof (double));
	double (*E)[2] = malloc( (N_time) * 2 * sizeof (double));
	double *p_spectrum = malloc((N_time) * sizeof (double));
	double *freq = malloc((N_time) * sizeof (double));
	/* Set variables */
	kappa = 1.6e3*kappa_unit;
	timestep = 0.01*sqrt(m[0]/kappa);
	timestep_sq = timestep * timestep;

	printf("%f\n",timestep);
	printf("Base frequencies: %f\t%f\n",sqrt(kappa/m[0])/2/PI,sqrt(kappa/m[0]*11/3)/2/PI);
	
	/* Initial conditions */
	/* Set initial displacements and velocites */
	q[0] = 0.20;
	q[1] = -0.15;
	q[2] = -0.05;
	v[0] = 0;

	for (i = 1; i < nbr_of_particles; i++) {
	  //q[i] = 0;
		v[i] = 0;
	}
	q_1[0] = q[0];
	q_2[0] = q[1];
	q_3[0] = q[2];

	/* Calculate initial accelerations based on initial displacements */
	calc_acc(a, q, m, kappa, nbr_of_particles);

	/* timesteps according to velocity Verlet algorithm */
	for (i = 1; i < N_time; i++) {
		/* v(t+dt/2) */
		for (j = 0; j < nbr_of_particles; j++) {
		    v[j] += timestep * 0.5 * a[j];
		} 

		/* q(t+dt) */
		for (j = 0; j < nbr_of_particles; j++) {
		    q[j] += timestep * v[j];
		}

		/* a(t+dt) */
		calc_acc(a, q, m, kappa, nbr_of_particles);

		/* v(t+dt) */
		for (j = 0; j < nbr_of_particles; j++) {
		    v[j] += timestep * 0.5 * a[j];
		} 

		/* Save the displacement of the three atoms */
		q_1[i] = q[0];
		q_2[i] = q[1];
		q_3[i] = q[2];
		
		
		E[i][0] = calc_pe(q, kappa, nbr_of_particles);
		E[i][1] = calc_ke(v, m, nbr_of_particles);
		
	}

	/* Print displacement data to output file */
	file = fopen("disp.dat","w");

	for (i = 0; i < N_time; i++) {
		current_time = i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e", current_time, q_1[i], q_2[i], q_3[i] );	
		fprintf(file, "\n");
	}
	fclose(file);
	
	file = fopen("energy.dat","w");
	for (i = 0; i < N_time; i++) {
		current_time = i * timestep;
		fprintf(file, "%.4f \t %e \t %e \t %e", current_time, 
			E[i][0] + E[i][1], E[i][0], E[i][1]);	
		fprintf(file, "\n");
	}
	fclose(file);	
	
	
	/* make FFT (powerspectrum) */
	/* Frequency given in inverse time units:
	   [ps^(-1)]=[THz]
	   
	   
	 */
	double *data_to_plot = q_3;
	powerspectrum(data_to_plot, p_spectrum, N_time);
	powerspectrum_shift(p_spectrum, N_time);
	fft_freq_shift(freq, timestep, N_time);

	/*Save powerspectrum data in file */
  	file = fopen("powerspectrum.dat","w");
	for (i = 0; i < N_time; i++)	{
		fprintf (file,"%e \t %e\n", freq[i], p_spectrum[i]);
	}
	fclose(file);

	/* Free allocated memory */ 
	free(q_1); q_1 = NULL;
	free(q_2); q_2 = NULL;
	free(q_3); q_3 = NULL;
	free(E); E = NULL; 
	return 0;    
}
