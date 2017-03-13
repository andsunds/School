#include <stdlib.h>
#include <stdio.h>
#include "main.h"


static void horizontal_check
(int index, int col, int cols, int *arr_len,int *NNa, int *NNb){
  /* The second level of the if statements in get_NN.
     Here, the horizontal check is perforemed.*/
  if ( col==0 ){
    /* If first column, then no NN to the left */
    *NNa=index+1;
    *arr_len=*arr_len-1;
  }else if ( col==cols-1 ){
    /* If last column, then no NN to the right */
    *NNa=index-1;
    *arr_len=*arr_len-1;
  }else{
    /* If inner column, then NN both to the left AND right */
    *NNa=index-1;
    *NNb=index+1;
  }
}


void get_NN(int *NN_arr, int index, int rows, int cols){
  /* Finds the nearest neighbours (NN) to the site 
     described by index. This method will take boundaries
     into account when calculating the NN's. 
       I.e. regular boundary sites will have 3 NN's,
       while corner sites only have 2 NN's, and inner
       sites have 4 NN's.

     The argument NN_arr has to be a pointer to an ARRAY
     of size _5_, otherwise expect big fuckups!
  */
  int row=get_row(index, cols);
  int col=get_col(index, cols);

  if ( index>=rows*cols || index<0 ){
    for (int a=0; a<5; a++ )
      NN_arr[a]=-1; 
    printf("\n ERROR in get_NN. Index out of bounds. \n   Got index=%d, should be 0<=index<%d.\n\n", index, rows*cols);
  }else{
    NN_arr[0]=4; /* there should be 4 NN
		    if there isn't, decreas by 1 (--) for each
		    case where a NN won't exist.
		 */
  }
  
  /* There are 2 levels of if statements:
     - Level 1 checks vertically, if we're at the top or bottom.
     - Level 2 checks horizontally, implemented in horizontal_check.
     If we're at the end in any of the directions, nbr_of_NN
     is decreased by 1, and the last (not already killed NN)
     is killed off by setting it to NULL.
  */
  if ( row==0 ){
    /* If first row, then no NN above */
    NN_arr[0]=NN_arr[0]-1; // decreases #NN
    NN_arr[1]=index+cols; //sets 1st NN
    horizontal_check
      (index, col, cols, &NN_arr[0], &NN_arr[2], &NN_arr[3]);
  }else if ( row==rows-1 ){ 
    /* If last row, then no NN below */
    NN_arr[0]=NN_arr[0]-1; // decreases #NN
    NN_arr[1]=index-cols; //sets 1st NN
    horizontal_check
      (index, col, cols, &NN_arr[0], &NN_arr[2], &NN_arr[3]);
  }else{
    /* If inner row, then NN both above AND below */
    NN_arr[1]=index-cols; //sets 1st NN
    NN_arr[2]=index+cols; //sets 2nd NN
    horizontal_check
      (index, col, cols, &NN_arr[0], &NN_arr[3], &NN_arr[4]);
  }  
}


