#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"

int *ising_init(int rows, int cols){
  /* Initializes a matrix, in the form of an
     array, filled with +-1 randomly.

     The array is allocated as a pointer inside
     this function, so it must be freed after use
     in the implementation. 
  */
  // Don't forget the sizeof!
  int N=rows*cols;
  int *ising_mtrx= malloc(N*sizeof(int)); 
  for( int i=0; i<N; i++ ){
    if( rand()<RAND_MAX/2 ) //random # to decide +-1.
      ising_mtrx[i]=-1;
    else
      ising_mtrx[i]= 1;
  }
  /* Do not forget to free ising_mtrx after use
     in the implementeation! */
  return ising_mtrx;
}


double totE(double J, int *mtrx_as_arr, int rows, int cols){
  /* Returns the total energy of the system according to
     the Hamiltonian:
             H = -J \sum_{<i,j> is NN} s_i * s_j,
     where s_i and s_j are described in mtrx_as_arr.

     No need (in problem 1) to use double for J and the
     return value, but it might come in handy later. 
   */
  int sum = 0; // sum of spin products
  int NN[5];   // init for use be get_NN
  int N=rows*cols;
  for (int a=0; a<N; a++ ){
    /* Loop/sum over all sites */
    get_NN(NN, a, rows, cols);
    for ( int b=0; b<NN[0]; b++){
      /* For each site, loop/sum over all NN's. */
      sum += mtrx_as_arr[a]*mtrx_as_arr[NN[b+1]];
    }
  }
  return -J*sum/2;
  /* We have to divide by 2 due to double counting
     of NN's. We sum for EVERY site's NN's, which
     means that if A is NN to B, then B is NN to A. 
 */
}

double deltaE
(double J, int *mtrx_as_arr, int index, int rows, int cols){
  /* Calculates the change in energy after _one_ spin
     (at index) has been flipped. The calculation is
     based on the Hamiltonian:
             H = -J \sum_{<i,j> is NN} s_i * s_j,
     where s_i and s_j are described in mtrx_as_arr.
  */
  int sum = 0; // sum of spin products
  int NN[5];   // init for use be get_NN
  get_NN(NN, index, rows, cols);
  for ( int b=0; b<NN[0]; b++){ // loop over all NN's
    sum += mtrx_as_arr[index]*mtrx_as_arr[NN[b+1]];
  }

  return -J*sum;
  /* The return should realy be:
           -J*( sum - (-sum) )/2  =  -J*(2*sum)/2.
     This is because dE= E_f-E_i, where the only differnce
     between E_f and E_i occurs due to the flipped spin at
     site index. I.e. summation would look like
       sum_f += +mtrx_as_arr[a]*mtrx_as_arr[NN[b+1]];
       sum_i += -mtrx_as_arr[a]*mtrx_as_arr[NN[b+1]];
                ^-- Here lies the only differnce,
		    i.e. sum_i = -sum_f.     
  */
}


double order_param(int *mtrx_as_arr, int N){
  /* Returns the order parameter of the sysyem:
        M = (1/N) * \sum_{all i} s_i,
     where the spin values s_i are described in the
     array mtrx_as_arr, and N is the total number of
     spins. 
   */
  int sum = 0; // sum of spins
  for (int a=0; a<N; a++ ){
    /* Loop/sum over all sites */
    sum += mtrx_as_arr[a];
  }

  return ( (double)sum )/( (double)N );
  /* Convert the sum and total # spins to double,
     then calculate the average. 
 */
}


















static void motecarlo_ising_step
(int *ising, int rows, int cols, int N,
 double *arr_EM, double J, double beta, double *E_tot,
 int current_iteration){
  /* Flips a random site and checks whether or not to
     keep it. 
  */
  double dbl_rand; //a random dbl, to be used in a check
  double dbl_RAND_MAX = (double)RAND_MAX; // a dbl const of RAND_MAX

  int random_index = rand() % N; // random index.
  ising[random_index] = -ising[random_index]; //flip
  double dE = deltaE(J, ising, random_index, rows, cols); // deltaE
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
  *E_tot  = (*E_tot) + dE;
  /* Writes the new values to the array passed to this
     function. 

     MAKE SURE TO KEEP TRACK OF <current_iteration> WHEN
     PASSING IT TO THIS FUNCTION!!!
  */
  arr_EM[2*current_iteration]   = *E_tot; 
  arr_EM[2*current_iteration+1] = order_param(ising, N);
}



int montecarlo_ising_full
(int rows, int cols, double J, double beta, int Nsteps,
 int chunk, char *save_directory, FILE *logPTR){
  /* This function Monte Carlo simulates a 2D ising model
     and writes the energy, E, and order paramter, M, to a
     binary file. 

     The Ising model is that of a grid of size <rows>*<cols>
     with Hamiltonian:
               H = -J * \sum_{<i,j> NN} s_i*s_j
     at temperature 1/<beta>.
     <Nsteps> of Monte Carlo simulations are performed.


     - The output file is named ("beta_%1.2f.bin",<beta>) 
       in the directory specified in <save_directory>
       and formated:
             {E(1),M(1),E(2),M(2), ... ,E(end),M(end)}.
     - It contains doubles (8 bytes).
     - The size of the file is:
        ( 1 + <Nsteps>*2/<chunk> )*( <chunk>/2 )*8 bytes,
       i.e. up to
             8*(<Nsteps> + <chunk>) bytes.
       The facor 8 is a consequence of the data type
       beeing double, which is 8 bytes. 

     Output is written to the file in chunks of size
                    <chunk> doubles
     at a time. This will supposedly make the file I/O
     faster. At least this method is save memory compared
     to writing everything to the file at the end.
     
     <logPTR> is a pointer to a log file, to which the
     printpouts will be written.
   */


  /* Initializations */
  int N=rows*cols; //total # sites
  double E; //energy and its change
   
  if ( chunk % 2 == 1 )
    chunk++; //makes sure chunk is even
  double arr_EM[chunk]; //values to be written to file
  int *ising=ising_init(rows, cols); // init of random ising grid

  // The length of the loops. These values are such
  // that loop1*loop2 should be just over Nsteps.
  int loop1 = Nsteps*2/chunk + 1; //length of outer loop
  int loop2 = chunk/2; //length of inner loop


  /* Lots of differents messages */
  char file_errorMSG[] = " ERROR in montecarlo_ising_full(): unable to open file: ";
  char endMSG[256]; // Tells us which file we've written to


  /*   I/O   */
  char filename[64]; //longer then the filename
  sprintf(filename,"%sbeta_%0.5f.bin",save_directory,beta);
  FILE *filePTR;
  filePTR=fopen(filename,"wb");
  if ( !filePTR ){ //check if the file opened.
    printf(        "%s%s\n",file_errorMSG,filename);
    fprintf(logPTR,"%s%s\n",file_errorMSG,filename);
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
       to file in a chunk of size 2*chunk.
    */
    for (int b=0; b<loop2; ++b ){ // for #2
      /* Performs an Ising Monte Carlo step. 

	 This step requires the total energy, E, to be
	 passed to it, because the total energy is only
	 updated via the calculation of deltaE.
 */
      motecarlo_ising_step
	(ising, rows, cols, N, arr_EM, J, beta, &E, b);
    } // end for #2
    // printf("%2.2f\n",arr_EM[chunk-1]); //DEBUG
    /* Writes the contents of this chuck to file. */
    fwrite(&arr_EM, sizeof(double), chunk, filePTR);
  }//end for #1
  /* The loop structure is as it is, so that we can
     write data to the file in smaller chunks. This
     way it's supposed to be faster, and this will
     definetly save memory usage.
  */

  sprintf(endMSG,"\nSimulation done, data written to:\n   %s\n",
	  filename);
  printf(        "%s",endMSG);
  fprintf(logPTR,"%s",endMSG);

  fclose(filePTR);
  free(ising); ising=NULL;
  return 0;
}







int montecarlo_ising_average
(int rows, int cols, 
 double J, double beta, double *return_values,
 int Nsteps, int discard_first){
  /* This function Monte Carlo simulates a 2D ising model
     and retruns the average and std of the energy, E, and
     order paramter, M. The values are retuned in the array
     <return_values> in the fashion:
                   {T, E, sdtE, M, stdM},
     where T=1/<beta>.

     The Ising model is that of a grid of size <rows>*<cols>
     with Hamiltonian:
               H = -J * \sum_{<i,j> NN} s_i*s_j
     at temperature 1/<beta>.
     <Nsteps> of Monte Carlo simulations are performed.

     MAKE SURE THAT <return_values> IS AN ARRAY WITH 
                             5
     ELEMENTS.
  */


  /* Initializations */
  int N=rows*cols; //total # sites
  double E, M; 

  // Giving the pointers (human) understandable names.
  return_values[0] = 1/beta;
  double *meanE = return_values +1; *meanE = 0;
  double *stdE  = return_values +2; *stdE  = 0;
  double *meanM = return_values +3; *meanM = 0;
  double *stdM  = return_values +4; *stdM  = 0;
  
  double arr_EM[2]; //array of current E and M values
  int *ising=ising_init(rows, cols); // init of random ising grid


  /* It's cheeper to calculate deltaE each iteration
     and just add that to E, to get the next energy
     value.

     This has a minor flaw in that, we don't get the
     very first value of E, but that should not be any
     problems.
  */
  E = totE(J, ising, rows, cols);

  for (int a=0; a<discard_first; ++a ){ //for #1
    /* In this loop we just perform the Monte Carlo
       step <discard_first> number of times, to get
       the initial warm-up oscillations out.
    */
    motecarlo_ising_step
      (ising, rows, cols, N, arr_EM, J, beta, &E, 0);
     //Passing 0 as the <current_index> because
     //arr_EM only has 2 elements.
  }//end for #1

  for (int b=0; b<Nsteps; ++b ){ //for #2
    /* This is one ising step */
    motecarlo_ising_step
      (ising, rows, cols, N, arr_EM, J, beta, &E, 0);
    //Passing 0 as the <current_index> because
    //arr_EM only has 2 elements.
    M = order_param(ising, N);

    /* These expressions are not fully the mean and sdt,
       but they will be modified after the loop.
    */
    *meanE += E; 
    *stdE  += E*E;
    *meanM += M; 
    *stdM  += M*M;
  }//end for #2

  /* Calculating the means */
  *meanE /= Nsteps;
  *meanM /= Nsteps;

  /*        Calculating the std's.
     In this stage the *stdX's are just sums of X*X (X^2).

     var = <x^2> - <x>^2 = sum x^2/N - sum (x/N)^2
     std = sqrt( var )
     In the estimates for var one should use N-1, however
     in this case N>>1, so it makes no difference.
  */
  *stdE = sqrt( *stdE/Nsteps - ( *meanE )*( *meanE ) );
  *stdM = sqrt( *stdM/Nsteps - ( *meanM )*( *meanM ) );

  free(ising); ising=NULL;
  return 0;
  /* Remember that the important values returned by this
     fuction are in <return_values>.
  */
}
