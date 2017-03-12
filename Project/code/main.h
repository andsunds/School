#ifndef MAIN_H_
#define MAIN_H_

/* nearest_neighbour.c */

void get_NN(int *NN_arr, int index, int rows, int cols);


/* ising_model.c */
int *ising_init(int rows, int cols);

double hamiltonian(double J, int *mtrx_as_arr, int rows, int cols);


/* grid2D.c */
int index2d(int row, int col, int cols);

int get_row(int index, int cols);

int get_col(int index, int cols);

void print_matrix(int *matrix_as_arr, int rows, int cols);

void print_matrix_sign(int *matrix_as_arr, int rows, int cols);

#endif
