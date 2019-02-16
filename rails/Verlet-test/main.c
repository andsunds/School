#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "harm_osc.h"

int main() {
  double alpha=0.9;
  double dt=1e-2;
  int N_steps=100000;
  double XV0[2]={1,0};
  double (*XV)[]=malloc(2*N_steps*sizeof(double));
  //double (*XV)[N_steps]=malloc(sizeof(double[N_steps][2]));
  //double XV[N_steps][2];

  char filename_verlet[]="verlet.tsv";

  stepper_verlet(XV0, XV, dt, N_steps);
  printf("Saving to %s\n",filename_verlet);
  save_to_file(filename_verlet, XV, dt, N_steps);

  /*
  char filename_euler[]="euler2F.tsv";
  stepper_euler2F(XV0, XV, dt, N_steps);
  printf("Saving to %s\n",filename_euler);
  save_to_file(filename_euler, XV, dt, N_steps);
  */

  free(XV); XV=NULL;
  return(0);
}

