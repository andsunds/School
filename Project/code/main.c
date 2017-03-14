/* N.B.
   In this implementation all matrices will be
   represented as a 1D array of length rows*cols.
*/
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
  double J=1;            // energy factor in tha Hamlitonian
  double beta = 0.0002;  // invers temperature
  double factor = 1.096; // this ^ 100 = 1e4
  int Nsims = 1;         // # simulations, set to 0 if not active
  int Nsteps = (int)1e7; // # Monte Carlo steps
  int write_chunk = 256; // # doubles to be written at once
  
  clock_t begin2, end2; //init, for time tracking
  float time_spent;
  char save_directory[] = "../data/bin/";

  /* Pre-setting a beta value.
    for (int p=0; p<26; ++p )
    beta = beta/factor;
  */
    printf("\n   Starting %d simulation(s) of %2g steps each (write all).\n\n", Nsims, (float)Nsteps);
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


static int simulate_ising_write_avg_to_tsv(){
  /* As the name suggests, simulates a bunch of
     differnt times and then oly writes the averages
     and std's of E and M to a tsv-file.
  */
  srand(time(NULL));// seeds rand, should only be called once
  
  int no_of_values = 5;        // do NOT change this! 
  double TEMarr[no_of_values]; // array with mean and std of E & M.
  int    L=16;                 // side length og the grid
  double J=1;                  // energy factor in the Hamlitonian
  double beta0 = 0.99;         // invers temperature
  double d_beta = 1e-3;        // step size in beta
  double beta;                 // init
  int Nsims      = 256;        // # sims, set to 0 if not active
  int Nsteps     = (int)2e7;   // # Monte Carlo steps
  int disc_first = (int)5e4;   // discard values from warm-up period

  int isOK;                    //init, for error checking
  clock_t begin2,end2, begin3; //init, for time tracking
  float runtime2, runtime3;    //init, for time tracking

  /*    I/O    */
  char filename[] = "../data/TEMstdEstdM_beta_.99-1.246_256.tsv"; 
  FILE *filePTR;
  filePTR=fopen(filename,"w");
  if ( !filePTR ){ //check if the file opened.
    printf(" ERROR in simulate_ising_write_avg_to_tsv(): unable to open file: %s\n",filename);
    return 1; // 1 means that something went wrong.
  }

  /* Starting the computations */
  printf("\n   Starting %d simulation(s) of %2g steps each (write avg).\n \n", Nsims, (float)Nsteps);
  begin3 = clock();

  for ( int i=0; i<Nsims; ++i ){
    /* A timer for how long the execution took */
    begin2 = clock();
    /* Here is the actual calculation */
    beta = beta0 + i*d_beta;
    isOK = montecarlo_ising_average
            (L,L, J, beta, TEMarr, Nsteps, disc_first);
    if ( isOK != 0 ){
      /* Just a small check to see that nothing went wrong
	 in the last simulation.
      */
      printf(" ERROR! \n\n");
      return 1;
    }
    /* Write to file */
    for ( int j=0; j<no_of_values; ++j ){
      fprintf(filePTR, "%16e",TEMarr[j]);
    }
    fprintf(filePTR, "\n"); //newline after the simulation

    end2 = clock();
    runtime2= (float)(end2 - begin2) / CLOCKS_PER_SEC;
    runtime3= (float)(end2 - begin3) / CLOCKS_PER_SEC;
    printf("Execution time (beta=%0.4f): %0.3f s. ( %0.1f s )\n", 
	   beta, runtime2, runtime3);
  } // end for

  fclose(filePTR);
  return 0;
}





int main(){
  /* Timings */
  clock_t begin1, end1; //init for time tracking
  float time_spent;
  begin1 = clock();
  /* Computation goes in here: */

  //int isOK = simulate_ising_write_all_to_bin();
  int isOK = simulate_ising_write_avg_to_tsv();
  //int N=L*L;int index=4; play_w_print_matrix(L, N, index);
  //play_w_ising_init(L, J);

  end1 = clock();
  time_spent = (float)(end1 - begin1) / CLOCKS_PER_SEC; 
  int h   = floor(time_spent/3600);
  int m   = floor(time_spent/60) - 60*h; 
  float s = time_spent-3600*h-60*m;

  printf("Total execution time: %dh %02dm %0.3fs. \n",
	 h, m, s);
  
  return isOK;
}
