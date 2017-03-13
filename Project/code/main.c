/* N.B.
   In this implementation all matrices will be
   represented as a 1D array of length rows*cols.
*/
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"

void play_w_print_matrix(int L,int N, int index){
  int mtrx[N];
  for (int a=0; a<N; a++ )
    mtrx[a]=a; // initializing
  
  print_matrix(mtrx, L, L);
  int NN[5]; // VERY important that it's _5_ here!
  get_NN(NN, index, L, L);
  for (int b=0; b<NN[0]; b++)
    printf("NN%d: %2d\n", b+1, NN[b+1]);
}

void play_w_ising_init(int L, double J){
  int *pmtrx=ising_init(L, L);
  /* //creates non-random marix
  int pmtrx[N];
  for (int c=0; c<N; c++ ){
    //pmtrx[c] = 1; // const, min energy
    pmtrx[c] = 1 - 2*(c%2); // alternating, max E
  } */
  print_matrix_sign(pmtrx, L, L);
  double E=totE(J, pmtrx, L,L);
  printf("Total energy: %3.1f\n",E);
  free(pmtrx); pmtrx=NULL;
}


static int simulate_ising_write_all_to_bin(){
  /* As the name suggests, simulates a bunch
     of differnt times and writes all data to
     a binary file.
  */
  srand(time(NULL)); //seeds rand, should only be called once
  int L=16;

  double J=1; // energy factor in tha Hamlitonian
  double beta = 0.0002; // invers temperature
  double factor = 1.096; // this ^ 100 = 1e4
  int Nsims = 1; // # simulations, set to 0 if not active
  int Nsteps = (int)1e4; // # Monte Carlo steps
  int write_chunk = 256; // # doubles to be written at once
  
  clock_t begin2, end2; //init for time tracking
  float time_spent;
  char save_directory[] = "../data/bin/";

  /* Pre-setting a beta value.
    for (int p=0; p<26; ++p )
    beta = beta/factor;
  */
  printf("\n   Starting %d simulation(s) of %2g steps each.\n\n",
	 Nsims, (float)Nsteps);
  for ( int i=0; i<Nsims; ++i ){
    //    beta = beta/factor; // the new/next value 
    
    /* A timer for how long the execution took */
    begin2 = clock();

    /* Here is the actual calculation */
    int isOK = montecarlo_ising_full(L,L, J, beta, Nsteps,
				     write_chunk, save_directory);

    end2 = clock();
    time_spent = (float)(end2 - begin2) / CLOCKS_PER_SEC;
    printf("Execution time (beta=%0.4f): %0.3f s.\n", 
	   beta, time_spent);

    if ( isOK != 0 ){
      /* Just a small check to see that nothing went wrong
	 in the last simulation.
      */
      printf(" ERROR! \n\n");
      return 1;
    }
  } // end for
  return 0;
}








int main(){
  /* Timings */
  clock_t begin1, end1; //init for time tracking
  float time_spent;
  begin1 = clock();
  /* Computation goes in here: */

  int isOK = simulate_ising_write_all_to_bin();
  //int N=L*L;int index=4; play_w_print_matrix(L, N, index);
  //play_w_ising_init(L, J);

  end1 = clock();
  time_spent = (float)(end1 - begin1) / CLOCKS_PER_SEC / 60; //min
  printf("Total execution time: %0.2f min.\n", time_spent);
  
  return isOK;
}
