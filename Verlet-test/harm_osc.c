/***
1D Harmonic oscillator
***/
#include <stdio.h>

double acc(double x){
  return -x;
}

/*
  This is not quite the explicit Euler forward method, sice the velocity,
  XV_new[1], is updated using the acceleration at the new position,
  acc(k, m, XV_new[0]). This makes this algorithm semi-implicit and also more 
  stable.
 */
void euler2F_step(double XV_old[], double XV_new[], double dt){
  /* XV arrays are structured as XV[0]=x, XV[1]=v. */
  XV_new[0]=XV_old[0]+XV_old[1]*dt;
  XV_new[1]=XV_old[1]+acc(XV_new[0],alpha)*dt;
}
void stepper_euler2F(double XV0[2], double *XV[][2], double dt, int N_steps){
  // /*
  *XV[0][0]=XV0[0];
  *XV[0][1]=XV0[1];
  // */
  //XV[0][0]=&XV0;
  /* Starting at the second time step */
  for (int i=1; i<N_steps; i++)
    //                v-- old  v-- new phase-space vector.
    euler2F_step(*XV[i-1], *XV[i], dt);
}

/*
  This Verlet velocity algorithm is not explicit for the same reson as the above
  Euler algorithm. The only difference is that the velocity is updated in two half
  steps, and the position is updated inbetween. This results in much smaller energy 
  fluctuations.
 */
void verlet_step(double XV_old[], double XV_new[], double dt){
  /* XV arrays are structured as XV[0]=x, XV[1]=v. */
  double dt_half=0.5*dt; //pre-calculate size of the half step.
  XV_new[1]=XV_old[1]+acc(XV_old[0], alpha)*dt_half;
  XV_new[0]=XV_old[0]+XV_new[1]*dt;
  XV_new[1]=XV_new[1]+acc(XV_new[0], alpha)*dt_half;
}
void stepper_verlet(double XV0[2], double XV[][2], double dt, int N_steps){
  XV[0][0]=XV0[0];
  XV[0][1]=XV0[1];
  /* starting at the second time step */
  for (int i=1; i<N_steps; i++)
    //                v-- old  v-- new phase-space vector.
    verlet_step(XV[i-1], XV[i], dt);
}

void save_to_file(char *filename, double XV[][2], double dt, int N_steps){
  double t=0;
  FILE *savefile=fopen(filename, "w"); // open the file
  for (int i=0; i<N_steps; i++){
    fprintf(savefile, "%0.8f\t",t); // prints the time
    for (int j=0; j<2; j++) // prints the position and velocity
      fprintf(savefile, "%0.8f\t", XV[i][j]);
    fprintf(savefile, "\n");// newline
    t+=dt; // updates the time
  } 
  fclose(savefile);
}
