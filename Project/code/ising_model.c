#include <stdlib.h>
#include "main.h"

int *ising_init(int rows, int cols){
  /* Initializes an Ising matrix */
  // Don't forget the sizeof!
  int *ising= malloc(rows*cols*sizeof(int)); 
  for( int i=0; i<rows*cols; i++ ){
    if( rand()<RAND_MAX/2 ) //random # to decide +-1.
      ising[i]=-1;
    else
      ising[i]= 1;
  }
  /* Do not forget to free ising after use
     in implementeation! */
  return ising;
}

double hamiltonian(double J, int *mtrx_as_arr, int rows, int cols){
  
  int sum = 0; 
  int NN[5];
  for (int a=0; a<rows*cols; a++ ){
    get_NN(NN, a, rows, cols);
    for ( int b=0; b<NN[0]; b++){
      sum += mtrx_as_arr[a]*mtrx_as_arr[NN[b+1]];
    }
  }

  return -J*sum/2;
}
