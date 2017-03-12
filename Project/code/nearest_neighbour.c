#include <stdlib.h>
#include <stdio.h>
#include "main.h"


static void horizontal_check
(int index, int col, int cols, int *arr_len,int *NNa, int *NNb){
  /* The second level of the if statements in get_NN.
     Here, the horizontal check is perforemed.*/
  if ( col==0 ){ //............2
    //printf("DEBUG 2a\n"); //DEBUG
    *NNa=index+1;
    *arr_len=*arr_len-1;
    //printf("#NN = %d\n", *arr_len); //DEBUG
  }else if ( col==cols-1 ){ //...2
    *NNa=index-1;
    *arr_len=*arr_len-1;
  }else{ //....................2
    *NNa=index-1;
    *NNb=index+1;
  }
}



int *get_NN(int index, int rows, int cols){
  /* Returns a pointer to an array of the form
     { #NN, NN1, NN2, ?NN3, ?NN4 }

            NN = Nearest Neighbour 
  */
  int row=get_row(index, cols);
  int col=get_col(index, cols);

  int NN1, NN2, NN3, NN4; // inti
  int nbr_of_NN=4; /* there should be 4 NN
		      if there isn't, decreas by 1 (--) */
  
  /* There are 2 levels of if statements:
     - Level 1 checks vertically, if we're at the top or bottom.
     - Level 2 checks horizontally, implemented in horizontal_check.
     If we're at the end in any of the directions, nbr_of_NN
     is decreased by 1, and the last (not already killed NN)
     is killed off by setting it to NULL.
  */
  if ( row==0 ){ //__________1
    NN1=index+cols;
    nbr_of_NN--;
    horizontal_check(index, col, cols, &nbr_of_NN, &NN2, &NN3);
  }else if ( row==rows-1 ){ //_1
    NN1=index-cols;
    nbr_of_NN--;
    horizontal_check(index, col, cols, &nbr_of_NN, &NN2, &NN3);
  }else{//___________________1
    NN1=index-cols;
    NN2=index+cols;
    horizontal_check(index, col, cols, &nbr_of_NN, &NN3, &NN4);
  }
  /*   DEBUG
    printf("%2d, %2d, %2d, %2d \n",NN1, NN2, NN3, NN4); */

  /* Creating an array, of NN's, and a pointer to return */
  int NN_arr[5]={nbr_of_NN, NN1, NN2, NN3, NN4};
  int *ret = malloc(5*sizeof(int));
  ret=NN_arr;
  return ret;
}
