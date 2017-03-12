/* N.B.
   In this implementation all matrices will be
   represented as an 1D array of length rows*cols.
 */
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include "main.h"


int main(){
  int L=3;
  int N=L*L;
  double J=1.0;
  // seeds rand, should only be called once
  srand(time(NULL)); 

  /*
  int mtrx[N];
  for (int a=0; a<N; a++ )
    mtrx[a]=a;
  
  print_matrix(mtrx, L, L);
  int index = 49;
  int NN[5]; // VERY important that it's _5_ here!
  get_NN(NN, index, L, L);
  for (int b=0; b<NN[0]; b++)
    printf("NN%d: %2d\n", b+1, NN[b+1]);
  */

  int *pmtrx=ising_init(L, L);
  /* int pmtrx[N];
  for (int c=0; c<N; c++ )
  pmtrx[c] = 1 - 2*(c%2); */
  print_matrix_sign(pmtrx, L, L);
  double E=hamiltonian(J, pmtrx, L,L);
  printf("Total energy: %3.2f\n",E);
  free(pmtrx);pmtrx=NULL;
  return 0;
}
