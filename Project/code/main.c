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










int montecarlo(int rows, int cols, 
		double J, double beta,
		int Nsteps, int write_chunk){
  /* This function Monte Carlo simulates a 2D ising model
     and writes the energy, E, and order paramter, M, to a
     binary file. 

     The Ising model is that of a grid of size [rows*cols]
     with Hamiltonian:
               H = -J * \sum_{<i,j> NN} s_i*s_j
     at temperature 1/[beta].
     [Nsteps] of Monte Carlo simulations are performed.


     - The output file is named ("beta_%1.2f.bin",[beta]) 
       and formated:
             {E(1),M(1),E(2),M(2), ... ,E(end),M(end)}.
     - It contains doubles (8 bytes).
     - The size of the file is:
        (1 + [Nsteps]/[write_chunk])*[write_chunk]*8 bytes,
       i.e. up to
             8*([Nsteps] + [write_chunk]) bytes.

     Output is written to the file in chunks of size
                 [write_chunk] doubles
     at a time. This will supposedly make the file I/O
     faster. At least this method is save memory compared
     to writing everything to the file at the end. 
   */



  /* Initializations */
  int N=rows*cols; //total # sites
  double E, dE; //energy and its change
  //  double threshold;
  int random_index; // random index that decides which site to flip
 
  if ( write_chunk % 2 == 1 )
    write_chunk++; //makes sure write_chunk is even
  double arr_EM[write_chunk]; //values to be written to file
  int *ising=ising_init(rows, cols); // init of random ising grid

  int loop1 = Nsteps*2/write_chunk + 1; //length of outer loop
  int loop2 = write_chunk/2; //length of inner loop

  double dbl_rand; //a random dbl, to be used in a check
  double dbl_RAND_MAX = (double)RAND_MAX; // a dbl const of RAND_MAX


  /* I/O */
  char filename[64]; //longer then the filename
  sprintf(filename,"../data/beta_%0.2f.bin",beta);
  FILE *filePTR;
  filePTR=fopen(filename,"wb");
  if ( !filePTR ){ //check if the file opened.
    printf(" ERROR: unable to open file: %s\n",filename);
    return 1; // 1 means that something went wrong.
  }


  /* It's cheeper to calculate deltaE each iteration
     and just add that to E, to get the next energy
     value.

     This has a minor flaw in that, we don't get the
     very first value of E, but that should not be any
     problems.
  */
  E = totE(J, ising, rows, cols);
  for (int a=0; a<loop1; ++a ){ // for #1
    /* In each iteration of this loop we write data
       to file in a chunk of size 2*write_chunk.
    */
    for (int b=0; b<loop2; ++b ){ // for #2
      /* Flip a random site and check if you want to
	 keep it. 
      */
      random_index = rand() % N; // random index.
      ising[random_index] = -ising[random_index]; //flip
      dE = deltaE(J, ising, random_index, rows, cols); // deltaE
      /* If deltaE is positive, we have to do an extra
	 check to see whether or not to keep the change.

	 If deltaE is negative, the change is keept and
	 we just continue on.
      */
      if ( dE>0 ){
	/* Generate a random number in [0,1] */
	dbl_rand = (double)rand()/dbl_RAND_MAX;
	if ( dbl_rand > exp(-beta*dE) ){
	  /* If r is too big, then flip back and the
	     change in energy is 0.
	  */
	  ising[random_index]=-ising[random_index];
	  dE=0;
	}
      }
      E = E + dE;
      arr_EM[2*b]   = E; //totE(J, ising, rows, cols);
      arr_EM[2*b+1] = order_param(ising, N);
    } // end for #2
    // printf("%2.2f\n",arr_EM[write_chunk-1]); //DEBUG
    /* Writes the contents of this chuck to file. */
    fwrite(&arr_EM, sizeof(double), write_chunk, filePTR);
  }//end for #1
  /* The loop structure is as it is, so that we can
     write data to the file in smaller chunks. This
     way it's supposed to be faster, and this will
     definetly save memory usage.
  */

  fclose(filePTR);
  free(ising); ising=NULL;
  return 0;
}







int main(){
  // seeds rand, should only be called once
  srand(time(NULL)); 

  int L=16;

  double J=1.0; // energy factor in tha Hamlitonian
  double beta = 1.0; // temp
  int Nsteps = 10000; // # Monte Carlo steps
  int write_chunk = 256; // # doubles to be written at once


int isOK = montecarlo(L,L, J, beta, Nsteps, write_chunk);
  if ( isOK != 0 )
    printf(" ERROR in montecarlo(). \n\n");


  //int N=L*L;int index=4; play_w_print_matrix(L, N, index);
  //play_w_ising_init(L, J);

  return 0;
}
