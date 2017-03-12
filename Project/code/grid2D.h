int index2d(int row, int col, int cols){
  /* Returns the 1D index corresponding to 
     the 2D index of (row,col) in a matrix 
     with cols # columns. 

     This is for when storing a 2D matrix a 
     a 1D array.
  */
  return row*cols+col;
}

int get_row(int index, int cols){
  /* Returns the row corresponding to the
     index provided.*/
  return index / cols;
}

int get_col(int index, int cols){
  /* Returns the row corresponding to the
     index provided.*/
  return index % cols;
}

void print_matrix(int *matrix_as_arr, int rows, int cols){
  /* Prints out the matrix to stdout. */
  for( int i=0; i<rows; i++ ){
    for( int j=0; j<cols; j++ ){
      printf(" %2d ", *(matrix_as_arr+index2d(i,j,cols)));
    }
    printf("\n");
  }
}

void print_matrix_sign(int *matrix_as_arr, int rows, int cols){
  /* Prints out the sign of each element to stdout. */
  for( int i=0; i<rows; i++ ){
    for( int j=0; j<cols; j++ ){
      if(*(matrix_as_arr+index2d(i,j,cols))<0)
	printf(" - ");
      else
	printf(" + ");
    }
    printf("\n");
  }
}

/*
int * zeros(int rows, int cols){
  // Initializes a zero matrix 
  static int z[rows*cols];
  for( int i=0; i<rows*cols; i++ ){
    z[i]=0;
  }
  return z;
}
*/
