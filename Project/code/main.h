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

void get_all_NN(int *arr_all_NN, int rows, int cols);
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


void get_NN_periodic(int *arr_NN_index, int index, int rows, int cols);
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


void get_all_NN_periodic(int *arr_all_NN, int rows, int cols);
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

/***************** ising_model.c *****************/


int montecarlo_ising_full
(int rows, int cols, double J, double beta, int Nsteps, int isP,
 int chunk, char *save_directory, FILE *logPTR);
  /* This function Monte Carlo simulates a 2D ising model
     and writes the energy, E, and order paramter, M, to a
     binary file. 

     The Ising model is that of a grid of size <rows>*<cols>
     with Hamiltonian:
               H = -J * \sum_{<i,j> NN} s_i*s_j
     at temperature 1/<beta>.
     <Nsteps> of Monte Carlo simulations are performed.


     - The output file is named ("beta_%1.2f.bin",<beta>) 
       and formated:
             {E(1),M(1),E(2),M(2), ... ,E(end),M(end)}.
     - It contains doubles (8 bytes).
     - The size of the file is:
        ( 1 + <Nsteps>*2/<chunk> )*( <chunk>/2 )*8 bytes,
       i.e. up to
             8*(<Nsteps> + <chunk>) bytes.
       The facor 8 is a consequence of the data type
       beeing double, which is 8 bytes. 

     Output is written to the file in chunks of size
                    <chunk> doubles
     at a time. This will supposedly make the file I/O
     faster. At least this method is save memory compared
     to writing everything to the file at the end. 

     <logPTR> is a pointer to a log file, to which the
     printpouts will be written.
   */

int montecarlo_ising_average
(int rows, int cols, 
 double J, double beta, double *return_values,
 int Nsteps, int discard_first, int isP);
  /* This function Monte Carlo simulates a 2D ising model
     and retruns the average and std of the energy, E, and
     order paramter, M. The values are retuned in the array
     <return_values> in the fashion:
                   {T, E, sdtE, M, stdM},
     where T=1/<beta>.

     The Ising model is that of a grid of size <rows>*<cols>
     with Hamiltonian:
               H = -J * \sum_{<i,j> NN} s_i*s_j
     at temperature 1/<beta>.
     <Nsteps> of Monte Carlo simulations are performed.

     MAKE SURE THAT <return_values> IS AN ARRAY WITH 
                             5
     ELEMENTS.
  */



/***************** XY_model.c *****************/
int montecarlo_XY_average
(int rows, int cols, 
 double J, double beta, double *return_values,
 int Nsteps, int discard_first, int isP);
  /* This function Monte Carlo simulates a 2D ising model
     and retruns the average and std of the energy, E, and
     order paramter, M. The values are retuned in the array
     <return_values> in the fashion:
                   {T, E, sdtE, rho_sX, rho_sY},
     where T=1/<beta>.

     The Ising model is that of a grid of size <rows>*<cols>
     with Hamiltonian:
         H = -J * \sum_{<i,j> NN} cos(theta_i - theta_j),
     at temperature 1/<beta>.
     <Nsteps> of Monte Carlo simulations are performed.

     MAKE SURE THAT <return_values> IS AN ARRAY WITH 
                             5
     ELEMENTS.
  */

#endif
