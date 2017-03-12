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
  /* The argument NN_arr has to be a pointer to an ARRAY
     of size _5_, otherwise expect big fuckups!

            NN = Nearest Neighbour 
  */
  int row=get_row(index, cols);
  int col=get_col(index, cols);

  
  NN_arr[0]=4; /* there should be 4 NN
		  if there isn't, decreas by 1 (--) for each
		  case where a NN won't exist.
	       */
  
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






/*    Still based on passing an array and modifying it,
      however, this is probably somewhat less efficient
      than the above version.
*/
/* void get_NN(int *NN_arr, int index, int rows, int cols){ */
/*   /\* The argument NN_arr has to be a pointer to an ARRAY */
/*      of size _5_, otherwise expect big fuckups! */

/*             NN = Nearest Neighbour  */
/*   *\/ */
/*   int row=get_row(index, cols); */
/*   int col=get_col(index, cols); */

/*   int NN1, NN2, NN3, NN4; // inti */
/*   int nbr_of_NN=4; /\* there should be 4 NN */
/* 		      if there isn't, decreas by 1 (--) *\/ */
  
/*   /\* There are 2 levels of if statements: */
/*      - Level 1 checks vertically, if we're at the top or bottom. */
/*      - Level 2 checks horizontally, implemented in horizontal_check. */
/*      If we're at the end in any of the directions, nbr_of_NN */
/*      is decreased by 1, and the last (not already killed NN) */
/*      is killed off by setting it to NULL. */
/*   *\/ */
/*   if ( row==0 ){ //__________1 */
/*     NN1=index+cols; */
/*     nbr_of_NN--; */
/*     horizontal_check(index, col, cols, &nbr_of_NN, &NN2, &NN3); */
/*   }else if ( row==rows-1 ){ //_1 */
/*     NN1=index-cols; */
/*     nbr_of_NN--; */
/*     horizontal_check(index, col, cols, &nbr_of_NN, &NN2, &NN3); */
/*   }else{//___________________1 */
/*     NN1=index-cols; */
/*     NN2=index+cols; */
/*     horizontal_check(index, col, cols, &nbr_of_NN, &NN3, &NN4); */
/*   } */
/*   /\*   DEBUG */
/*     printf("%2d, %2d, %2d, %2d \n",NN1, NN2, NN3, NN4); *\/ */

/*   /\* Creating an array, of NN's, and a pointer to return *\/ */
/*   //int NN_arr[5]={nbr_of_NN, NN1, NN2, NN3, NN4}; */
/*   //int *ret = malloc(5*sizeof(int)); */
/*   //ret=NN_arr; */
/*   //return ret; */
/*   int arr[5]={nbr_of_NN, NN1, NN2, NN3, NN4}; */
/*   for (int a=0; a<5; a++ ) */
/*     NN_arr[a] = arr[a]; */
/* } */
