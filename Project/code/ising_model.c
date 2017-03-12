#include <stdlib.h>
#include "main.h"

int *ising_init(int rows, int cols){
  /* Initializes a matrix, in the form of an
     array, filled with +-1 randomly.

     The array is allocated as a pointer inside
     this function, so it must be freed after use
     in the implementation. 
  */
  // Don't forget the sizeof!
  int *ising_mtrx= malloc(rows*cols*sizeof(int)); 
  for( int i=0; i<rows*cols; i++ ){
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
  for (int a=0; a<rows*cols; a++ ){
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
