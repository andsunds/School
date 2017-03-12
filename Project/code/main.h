#ifndef NEAREST_NEIGHBOUR_H_
#define NEAREST_NEIGHBOUR_H_

static void horizontal_check
(int index, int col, int cols, int *arr_len,int *NNa, int *NNb);

int *get_NN(int index, int rows, int cols);


#endif

#ifndef ISING_MODEL_H_
#define ISING_MODEL_H_

int *ising_init(int rows, int cols);

#endif


#ifndef GRID2d_H_
#define GRID2d_H_

int index2d(int row, int col, int cols);

int get_row(int index, int cols);

int get_col(int index, int cols);

void print_matrix(int *matrix_as_arr, int rows, int cols);

void print_matrix_sign(int *matrix_as_arr, int rows, int cols);

#endif
