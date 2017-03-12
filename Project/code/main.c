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

int main(){
  // seeds rand, should only be called once
  srand(time(NULL)); 

  int L=3;
  int N=L*L;
  double J=1.0;
  int index = 4;

  play_w_print_matrix(L, N, index);

  play_w_ising_init(L, J);

  return 0;
}
