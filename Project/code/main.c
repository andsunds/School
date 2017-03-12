/* N.B.
   In this implementation all matrices will be
   represented as an 1D array of length rows*cols.
 */
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


int main(){
  int L=9;
  int N=L*L;
  // seeds rand, should only be called once
  srand(time(NULL)); 

  int mtrx[N];
  for (int a=0; a<N; a++ )
    mtrx[a]=a;
  
  print_matrix(mtrx, L, L);
  int index = 8;
  int *NN=get_NN(index, L, L);
  for (int b=0; b<*(NN); b++)
    printf("NN%d: %2d\n", b+1, *(NN+b+1));
  /* int *pmtrx=ising_init(L, L); */
  /* print_matrix_sign(pmtrx, L, L); */
  /* free(pmtrx); */
  return 0;
}
