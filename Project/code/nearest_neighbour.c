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


void get_NN(int *arr_NN_index, int index, int rows, int cols){
  /* Finds the nearest neighbours (NN) to the site 
     described by <index>. This method will take boundaries
     into account when calculating the NN's. 
       I.e. regular boundary sites will have 3 NN's,
       while corner sites only have 2 NN's, and inner
       sites have 4 NN's.

     <arr_NN_index> is set up as:
        {#NN's, NN1, NN2, (NN3), (NN4)}.
     

     The argument <arr_NN_index> has to be a pointer to an
     array of size _5_, otherwise expect big fuckups!
  */
  int row=get_row(index, cols);
  int col=get_col(index, cols);

  if ( index>=rows*cols || index<0 ){// Index out of bounds?
    for (int a=0; a<5; a++ )
      arr_NN_index[a]=-1; 
    printf("\n ERROR in get_NN. Index out of bounds. \n   Got index=%d, should be 0<=index<%d.\n\n", index, rows*cols);
  }else{
    /* Initially there are 4 NN's, but in the edge cases,
       we decrese this (--) by the appropriate amount.
    */
    arr_NN_index[0]=4;
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
    arr_NN_index[0]=arr_NN_index[0]-1; // decreases #NN
    arr_NN_index[1]=index+cols; //sets 1st NN
    horizontal_check
      (index, col, cols,
       &arr_NN_index[0], &arr_NN_index[2], &arr_NN_index[3]);
  }else if ( row==rows-1 ){ 
    /* If last row, then no NN below */
    arr_NN_index[0]=arr_NN_index[0]-1; // decreases #NN
    arr_NN_index[1]=index-cols; //sets 1st NN
    horizontal_check
      (index, col, cols,
       &arr_NN_index[0], &arr_NN_index[2], &arr_NN_index[3]);
  }else{
    /* If inner row, then NN both above AND below */
    arr_NN_index[1]=index-cols; //sets 1st NN
    arr_NN_index[2]=index+cols; //sets 2nd NN
    horizontal_check
      (index, col, cols,
       &arr_NN_index[0], &arr_NN_index[3], &arr_NN_index[4]);
  }  
}



void get_all_NN(int *arr_all_NN, int rows, int cols){
  /* This funtion returns an array <arr_all_NN> with
     information about all NN's.

     It is important that <arr_all_NN> is
                5*<rows>*<cols> = 5*N
     ints long. The factor 5 comes in because there
     are 4 NN's, and for each index there is also a
     value telling us 

     <arr_all_NN> is structured acording to:
      {#NN(i=0), NN1(i=0), NN2(i=0), (NN3(i=0)), (NN4(i=0)),
       #NN(i=1), NN1(i=1), NN2(i=1), (NN3(i=1)), (NN4(i=1)),
       ...,
       #NN(i=N), NN1(i=N), NN2(i=N), (NN3(i=N)), (NN4(i=N))}
     So to acces the relevant infrmation for index i, just
     use:
             for ( int j=0; j<5; ++j )
               ... = arr_all_NN[i*5 + j];
  */

  int N = rows * cols;
  //int arr_NN_index[5];

  for (int index = 0; index<N; ++index ){
    /* get_NN() takes a pointer and writes the information
       to the 5 following addresses. What we want is to
       get those pieces of information into the right spots
       in <arr_all_NN>.
    */
    get_NN( (arr_all_NN+5*index), index, rows, cols);
  }
}



/* /////////////// Periodic boundary conditions /////////////// */

static void horizontal_check_periodic
(int index, int col, int cols, int *NNa, int *NNb){
  /* The second level of the if statements in get_NN.
     Here, the horizontal check is perforemed.*/
  if ( col==0 ){
    /* If first column, then the NN left is on last column */
    *NNa = index + 1;        //sets 3rd NN
    *NNb = index + (cols-1); //sets 4th NN
  }else if ( col==cols-1 ){
    /* If first column, then the NN left is on last column */
    *NNa = index - 1;        //sets 3rd NN
    *NNb = index - (cols-1); //sets 4th NN
  }else{
    /* If inner column, then NN both to the left AND right as ususal */
    *NNa = index - 1;
    *NNb = index + 1;
  }
}



void get_NN_periodic(int *arr_NN_index, int index, int rows, int cols){
  /* Finds the nearest neighbours (NN) to the site 
     described by <index>. This method applies periodic boundary
     conditions to the NN's.

     <arr_NN_index> is set up as:
        {#NN's, NN1, NN2, (NN3), (NN4)}.
     (This is to keep backwards compatibility with the non-periodic
     version of this functino.)

     The argument <arr_NN_index> has to be a pointer to an
     array of size _5_, otherwise expect big fuckups!
  */
  int row=get_row(index, cols);
  int col=get_col(index, cols);

  if ( index>=rows*cols || index<0 ){// Index out of bounds?
    for (int a=0; a<5; a++ )
      arr_NN_index[a]=-1; 
    printf("\n ERROR in get_NN. Index out of bounds. \n   Got index=%d, should be 0<=index<%d.\n\n", index, rows*cols);
  }else{
    /* With periodic BC's there are always 4 NN's */
    arr_NN_index[0]=4; 
  }
  
  /* There are 2 levels of if statements:
     - Level 1 checks vertically, if we're at the top or bottom.
     - Level 2 checks horizontally, implemented in horizontal_check.
     If we're at the end in any of the directions, nbr_of_NN
     is decreased by 1, and the last (not already killed NN)
     is killed off by setting it to NULL.
  */
  if ( row==0 ){
    /* If first row, then the NN above is on last row */
    arr_NN_index[1] = index + cols;          //sets 1st NN
    arr_NN_index[2] = index + (rows-1)*cols; //sets 2nd NN
    horizontal_check_periodic
      (index, col, cols, &arr_NN_index[3], &arr_NN_index[4]);
  }else if ( row==rows-1 ){ 
    /* If last row, then the NN below is on 0th row */
    arr_NN_index[1] = index - cols;          //sets 1st NN
    arr_NN_index[2] = index - (rows-1)*cols; //sets 2nd NN
    horizontal_check_periodic
      (index, col, cols, &arr_NN_index[3], &arr_NN_index[4]);
  }else{
    /* If inner row, then NN both above AND below as usual */
    arr_NN_index[1] = index - cols; //sets 1st NN
    arr_NN_index[2] = index + cols; //sets 2nd NN
    horizontal_check_periodic
      (index, col, cols, &arr_NN_index[3], &arr_NN_index[4]);
  }  
}


void get_all_NN_periodic(int *arr_all_NN, int rows, int cols){
  /* This funtion returns an array <arr_all_NN> with
     information about all NN's, with _periodic_ boundary
     conditinos.

     It is important that <arr_all_NN> is
                5*<rows>*<cols> = 5*N
     ints long. The factor 5 comes in because there
     are 4 NN's, and for each index there is also a
     value telling us 

     <arr_all_NN> is structured acording to:
      {#NN(i=0), NN1(i=0), NN2(i=0), (NN3(i=0)), (NN4(i=0)),
       #NN(i=1), NN1(i=1), NN2(i=1), (NN3(i=1)), (NN4(i=1)),
       ...,
       #NN(i=N), NN1(i=N), NN2(i=N), (NN3(i=N)), (NN4(i=N))}
     So to acces the relevant infrmation for index i, just
     use:
             for ( int j=0; j<5; ++j )
               ... = arr_all_NN[i*5 + j];
  */

  int N = rows * cols;
  //int arr_NN_index[5];

  for (int index = 0; index<N; ++index ){
    /* get_NN() takes a pointer and writes the information
       to the 5 following addresses. What we want is to
       get those pieces of information into the right spots
       in <arr_all_NN>.
    */
    get_NN_periodic( (arr_all_NN+5*index), index, rows, cols);
  }
}
