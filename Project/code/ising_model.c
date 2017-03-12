#include <stdlib.h>
//#include "grid2D.h"
//#include "nearest_neighbour.h"
#include "main.h"

int *ising_init(int rows, int cols){
  /* Initializes an Ising matrix */
  // Don't forget the sizeof!
  int *ising= malloc(rows*cols*sizeof(int)); 
  for( int i=0; i<rows*cols; i++ ){
    if( rand()<RAND_MAX/2 )
      ising[i]=-1;
    else
      ising[i]= 1;
  }
  return ising;
}
