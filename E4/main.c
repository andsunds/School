#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "funcs.h"

#define PI 3.141592653589
#define N_part 30000
#define N_timesteps 1000
#define N_dist_times 5
#define N_dist_points 100

#define kB 1.38064852e-23


/* Main program */
int main()
{

  /* UNITS: kg, µs, µm  */

  double dt = 3; 
  
  double xlim = 0.15;
  double vlim = 3e-3;
  double dx = 2*xlim / N_dist_points;
  double dv = 2*vlim / N_dist_points;
  

  double omega0 = 2 * PI * 3e-3; 
  double T = 297;
  double m = 3.0134e-14;
  double v_th = sqrt(kB*T / m);  
  

  double tau =  48.5; //147.3/  48.5; // Time unit: µs, length unit µm 
  double eta = 1/tau;
  
  double c0 = exp( - eta * dt);
  double time;
  
  gsl_rng *q = initialize_rng();
  
  double *a = malloc(N_part*sizeof(double));
  double *v = malloc(N_part*sizeof(double));
  double *x = malloc(N_part*sizeof(double));
  
  double (*x_5)[5] = malloc(5*N_timesteps*sizeof(double));
  double (*v_5)[5] = malloc(5*N_timesteps*sizeof(double));
  
  double *mu_x = malloc(N_timesteps*sizeof(double));
  double *sigma_x = malloc(N_timesteps*sizeof(double));
  double *mu_v = malloc(N_timesteps*sizeof(double));
  double *sigma_v = malloc(N_timesteps*sizeof(double));
  
  double (*rho_x)[N_dist_times]= malloc(N_dist_times*N_dist_points*sizeof(double));
  double (*rho_v)[N_dist_times] = malloc(N_dist_times*N_dist_points*sizeof(double));
  int save_ind[N_dist_times] = {0, 33, 100, 300, N_timesteps-1};
  
  for (int i = 0; i<N_part; i++){
    x[i] = 0.1;
    v[i] = 2e-3;
    a[i] = -omega0*omega0 * x[i];
  }
  
  // step in time
  int k = 0;
  for (int i = 0; i<N_timesteps; i++){
    if ( i==save_ind[0] || i==save_ind[1] || i==save_ind[2] || i==save_ind[3] 
      || i==save_ind[4] ){
      for (int j = 0; j<N_part; j++){
        rho_x[(int)( (x[j] + xlim)/dx )][k] +=1;
        rho_v[(int)( (v[j] + vlim)/dv )][k] +=1;
      }  
      k++;
    }
  
  
    take_BD_step(a, v,x, c0, q, v_th, dt, N_part, omega0);
    for (int j = 0; j<N_part; j++){
      mu_x[i] += x[j];
      sigma_x[i] += x[j] * x[j];
      mu_v[i] += v[j];
      sigma_v[i] += v[j] * v[j];
    }
    

    for (int j = 0; j<5; j++){
      x_5[i][j] = x[j];
      v_5[i][j] = v[j];
    }
    mu_x[i] *= 1/(double)N_part;
    sigma_x[i] = sqrt( sigma_x[i]/N_part - mu_x[i]*mu_x[i] );
    mu_v[i] *= 1/(double)N_part;
    sigma_v[i] = sqrt( sigma_v[i]/N_part - mu_v[i]*mu_v[i] );
  }
  
  char file_name[256];
  
  /* Write phase space coordinates to file */
  
  sprintf(file_name,"xv_tau-%d.tsv", (int) tau);
	  
  FILE *file_pointer = fopen(file_name, "w");
  for (int i=0; i<N_timesteps; i++){
    time = i*dt;
    fprintf(file_pointer, "%.2f\t", time);
    for (int j=0;j<5;j++){
      fprintf(file_pointer, "%.8e \t%.8e\t ", x_5[i][j],v_5[i][j]);
    }
    fprintf(file_pointer,"\n");
  }
  fclose(file_pointer);
  
   /* write mu and sigma */
  sprintf(file_name,"mu_sigma_tau-%d.tsv", (int) tau);
  file_pointer = fopen(file_name, "w");
  
  for (int i=0; i<N_timesteps; i++){
    time = i*dt;
    fprintf(file_pointer, "%.2f\t %.8e\t%.8e\t%.8e\t%.8e\n ",
      time, mu_x[i], sigma_x[i], mu_v[i], sigma_v[i]);
  }
  fclose(file_pointer);
  
  
  /* write hist */
  sprintf(file_name,"dist_tau-%d.tsv", (int) tau);
  file_pointer = fopen(file_name, "w");
  for (int j=0; j<N_dist_times; j++){
    time = dt * save_ind[j];
    fprintf(file_pointer, "%.2f\t%.2f\t", time,time);
  }
  fprintf(file_pointer, "\n");
  for (int i=0; i<N_dist_points; i++){
    for (int j=0; j<N_dist_times; j++){
      fprintf(file_pointer, "%.8e\t%.8e\t ",rho_x[i][j]/N_part/dx,rho_v[i][j]/N_part/dv);
      }
    fprintf(file_pointer, "\n");
  }
  fclose(file_pointer);
  
  
  gsl_rng_free(q);
  free(a);
  free(v);
  free(x);
  return 0;    
}
