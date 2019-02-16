/***
1D Coupled Harmonic oscillator
***/
#include <stdio.h>
#include "harm_osc.h"


void stepper_verlet(double XV0[2], double XV[][2], double dt, int N_steps, int){
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
