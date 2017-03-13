/* N.B.
   In this implementation all matrices will be
   represented as a 1D array of length rows*cols.
*/
#ifndef MAIN_H_ //A check so to not rerun all the definitions below.
#define MAIN_H_


/******************** grid2D.c ********************/
int index2d(int row, int col, int cols);
  /* Returns the 1D index corresponding to the 2D
     index of (row,col) in a matrix with cols beeing
      the number of columns. 

     This is for when storing a 2D matrix a 1D array.
  */

int get_row(int index, int cols);
  /* Returns the row corresponding to the
     index provided.
  */

int get_col(int index, int cols);
  /* Returns the column corresponding to the
     index provided.
  */

void print_matrix(int *matrix_as_arr, int rows, int cols);
  /* Prints out the matrix element's value to stdout. 
  */

void print_matrix_sign(int *matrix_as_arr, int rows, int cols);
  /* Prints out the sign of each element to stdout.
  */



/************** nearest_neighbour.c **************/
void get_NN(int *NN_arr, int index, int rows, int cols);
  /* Finds the nearest neighbours (NN) to the site 
     described by index. This method will take boundaries
     into account when calculating the NN's. 
       I.e. regular boundary sites will have 3 NN's,
       while corner sites only have 2 NN's, and inner
       sites have 4 NN's.

     The argument NN_arr has to be a pointer to an ARRAY
     of size _5_, otherwise expect big fuckups!
  */



/***************** ising_model.c *****************/
int *ising_init(int rows, int cols);
  /* Initializes a matrix, in the form of an
     array, filled with +-1 randomly.

     The array is allocated as a pointer inside
     this function, so it must be freed after use
     in the implementation. 
  */

double totE
   (double J, int *mtrx_as_arr, int rows, int cols);
  /* Returns the total energy of the system according to
     the Hamiltonian:
             H = -J \sum_{<i,j> is NN} s_i * s_j,
     where s_i and s_j are described in mtrx_as_arr.

     No need (in problem 1) to use double for J and the
     return value, but it might come in handy later. 
   */

double deltaE
(double J, int *mtrx_as_arr, int index, int rows, int cols);
  /* Calculates the change in energy after _one_ spin
     (at index) has been flipped. The calculation is
     based on the Hamiltonian:
             H = -J \sum_{<i,j> is NN} s_i * s_j,
     where s_i and s_j are described in mtrx_as_arr.
  */

double order_param(int *mtrx_as_arr, int N);
  /* Returns the order parameter of the sysyem:
        M = (1/N) * \sum_{all i} s_i,
     where the spin values s_i are described in the
     array mtrx_as_arr, and N is the total number of
     spins. 
   */

int montecarlo_ising(int rows, int cols, 
		     double J, double beta,
		     int Nsteps, int write_chunk);
  /* This function Monte Carlo simulates a 2D ising model
     and writes the energy, E, and order paramter, M, to a
     binary file. 

     The Ising model is that of a grid of size [rows*cols]
     with Hamiltonian:
               H = -J * \sum_{<i,j> NN} s_i*s_j
     at temperature 1/[beta].
     [Nsteps] of Monte Carlo simulations are performed.


     - The output file is named ("beta_%1.2f.bin",[beta]) 
       and formated:
             {E(1),M(1),E(2),M(2), ... ,E(end),M(end)}.
     - It contains doubles (8 bytes).
     - The size of the file is:
 ( 1 + [Nsteps]*2/[write_chunk] )*( [write_chunk]/2 )*8 bytes,
       i.e. up to
             8*([Nsteps] + [write_chunk]) bytes.
       The facor 8 is a consequence of the data type
       beeing double, which is 8 bytes. 

     Output is written to the file in chunks of size
                 [write_chunk] doubles
     at a time. This will supposedly make the file I/O
     faster. At least this method is save memory compared
     to writing everything to the file at the end. 
   */



#endif
