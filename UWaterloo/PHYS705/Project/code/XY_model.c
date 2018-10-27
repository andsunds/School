#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "main.h"


static void theta_init(double *theta_mtrx, int rows, int cols){
  /* Initializes a matrix, passed to here, in the form of an
     array, filled with random angles between 0 and 2*pi.

     The array is allocated as a pointer inside
     this function, so it must be freed after use
     in the implementation. 
  */
  
  int N=rows*cols;
  double TWO_PI = 6.283185307;
  for( int i=0; i<N; i++ )
    theta_mtrx[i]= (rand()*TWO_PI)/RAND_MAX;
}

static double totE(double J, double *theta_mtrx, int rows, int cols,
		   int *arr_all_NN){
  /* Returns the total energy of the system according to
     the Hamiltonian:
      H = -J \sum_{<i,j> is NN} cos(theta_i * theta_j),
     where s_i and s_j are described in mtrx_as_arr.

     No need (in problem 1) to use double for J and the
     return value, but it might come in handy later. 

     NOT to be used inside simulation loops, since this
     function is not optimized. INSTEAD use this _once_
     before the loop, and then use deltaE() inside the
     loop to then get the new value. 
   */
  double sum = 0; // sum of cosines
  int N=rows*cols;
  int start;
  for (int a=0; a<N; a++ ){
    /* Loop/sum over all sites */
    start = a*5;
    for ( int b=1; b<=arr_all_NN[start]; b++){
      /* Loop over all NN's, and sum the products.
	 To access the NN's for <index> start at <index>*5,
	 and then check as many slots as are indicated in
	 that location.
      */
      sum += cos( theta_mtrx[a] - theta_mtrx[arr_all_NN[start+b]] );
    }
  }

  return -J*sum/2;
  /* We have to divide by 2 due to double counting
     of NN's. We sum for EVERY site's NN's, which
     means that if A is NN to B, then B is NN to A. 
  */
}

static double deltaE
(double J, double *theta_mtrx, int index, double old_theta,
 int rows, int cols, int *arr_all_NN){
  /* Calculates the change in energy after _one_ spin
     (at index) has been flipped. The calculation is
     based on the Hamiltonian:
             H = -J \sum_{<i,j> is NN} s_i * s_j,
     where s_i and s_j are described in mtrx_as_arr.
  */
  double new_sum = 0; // sum of new spin products
  double old_sum = 0; // sum of old spin products
  int start = index*5;
  for ( int b=1; b<=arr_all_NN[start]; b++){
    /* Loop over all NN's, and sum the cosines.
       To access the NN's for <index> start at <index>*5,
       and then check as many slots as are indicated in
       that location.
     */
    new_sum += cos( theta_mtrx[index] - theta_mtrx[ arr_all_NN[start+b] ] );
    old_sum += cos( old_theta         - theta_mtrx[ arr_all_NN[start+b] ] );
  }

  return J*(old_sum - new_sum); //should be -J*(new_sum - old_sum)
}



static double rho_s
(double J, double beta, double *theta_mtrx, int rows, int cols,
 int *arr_all_NN, int XorY){
  /* Calculates the "instantaneous spin stiffness" in the X (XorY==2)
     or Y (XorY==1) direction. 
 */
  int N = rows*cols;

  double cos_sum = 0;
  double sin_sum = 0;
  double theta1, theta2;

  int start_index;
  switch( XorY ){
  case 1: //The Y direction is represented by 1
    start_index = 1;  break;
  case 2: //The X direction is represented by 2
    start_index = 3;  break;
  default:
    printf("Error in choosing direction. Must be 1 or 2.\n");
    return -1;
  }
  
  int odd, site, site5;  //inits, for the "pick every other" algorithm
  int cols_by2 = cols/2; //don't calculate this each iteration.
  for (int a=0; a<rows; ++a){ 
    /* Since sin(.) is an odd funtion, we can not just sum over
       ALL sites and add up  sin(theta_i - theta_j)  for all sites.
       Then we would end up with a sum that is 0, since we would
       end up getting
          sin(theta_i - theta_j) + sin(theta_j - theta_i) = 0.
       We therefore have to only pick every other site in such a way
       that none of the picked sites are mutual NN. 

       This is done by picking every other column-index, but shifting
       it by one if the row-index is odd. 
    */
    odd=a%2;
    for (int b=0; b<cols_by2; ++b){
      site  = index2d(a, odd+2*b, cols); //col and row to 1D index
      site5 = site*5; //To access the NN's of <site> we need to go to 5*<site>
      /* These are the two NN's in the chosen direction. */
      theta1 = theta_mtrx[site] - theta_mtrx[arr_all_NN[site5+start_index   ]];
      theta2 = theta_mtrx[site] - theta_mtrx[arr_all_NN[site5+start_index +1]];
    
      cos_sum += cos(theta1) + cos(theta2);
      sin_sum += sin(theta1) - sin(theta2);
      /* Note the minus sig in the sin sum. This is becuase it is
	 not mearely enough to keep track of not double suming the
	 sines, but we also have to keep track in which direction
	 the NN is.

	 The NN function is set up in such a way so that 
	                arr_all_NN[site5+start_index] 
	 will always be the site in the positive direction, while 
                       arr_all_NN[site5+start_index +1]
	 will be the site in the negative direction relative <site>.
	 This is why there is a minus sign on one of the two sines. 
      */
    }
  }

  /* <rho_s> is the spin stiffness in the chosen direction. It is defined
     as:
       rho_s = (1/N)*sum_{<i,j>_a}( cos(theta_i-theta_j) )
               -(beta*J/N)*( sum_{<i,j>_a}(sin(theta_i-theta_j)) )^2,
     where a is the driection.
  */

  return ( cos_sum - beta*J*sin_sum*sin_sum )/N;
  /* No divisions by 2 here, since we're actually only summing over
     each 'bond', and not just over all sites.
  */
}

/*
static double rho_s_OLD //NOT WORKING
(double J, double beta, double *theta_mtrx, int rows, int cols,
 int *arr_all_NN, int XorY){
  /* Calculates the "instantaneous spin stiffness" in the X (XorY==2)
     or Y (XorY==1) direction. 
  * /
  int N = rows*cols;

  double cos_sum = 0;
  double sin_sum = 0;
  double theta1, theta2;

  int start_index;
  switch( XorY ){
  case 1: //The Y direction is represented by 1
    start_index = 1;  break;
  case 2: //The X direction is represented by 2
    start_index = 3;  break;
  default:
    printf("Error in choosing direction. Must be 1 or 2.\n");
    return -1;
  }
  
  /* Think of some way to not have to all this in each step! * /
  for (int a=0; a<N; ++a){ 
    /* These are the two NN's in the chosen direction. * /
    theta1 = theta_mtrx[a] - theta_mtrx[arr_all_NN[5*a+start_index   ]];
    theta2 = theta_mtrx[a] - theta_mtrx[arr_all_NN[5*a+start_index +1]];
    
    cos_sum += cos(theta1) + cos(theta2);
    sin_sum += sin(theta1) + sin(theta2);
  }

  /* <rho_s> is the spin stiffness in the chosen direction. It is defined
     as:
       rho_s = (1/N)*sum_{<i,j>_a}( cos(theta_i-theta_j) )
               -(beta*J/N)*( sum_{<i,j>_a}(sin(theta_i-theta_j)) )^2,
     where a is the driection.
  * /

  return ( cos_sum/2 - beta*J*sin_sum*sin_sum/4 )/N;
}
  */




/////////////////////////////////////////////////////////////////
static void motecarlo_XY_step
(double *theta_mtrx, int rows, int cols, int N,
 double J, double beta, double max_flip,
 double *E_tot, double *rhoX_sum, double *rhoY_sum,
 int *arr_all_NN 
 ){
  /* Makes a Monte Carlo step in the XY model. 
  */
  int random_index = rand() % N; // random index.
  double old_theta = theta_mtrx[random_index];
  /* Produces a random angle between -<max_flip> to +<max_flip. */
  double d_theta   = max_flip*2*((double)(rand()-RAND_MAX/2))/RAND_MAX;
  theta_mtrx[random_index] += d_theta; //updating theta_mtrx 

  double dE = deltaE(J, theta_mtrx, random_index,
		     old_theta, rows, cols, arr_all_NN);

  /* If dE is positive, we make anonther check to see wheter or not
     to keep the change. This extra check is based on the temperature
     (beta). 
  */
  if ( dE>0 ){ 
    if ( ((double)rand())/RAND_MAX > exp( -beta*dE ) ){
      /* If r is too big, then flip back
	 and the change in energy is 0.
      */
      theta_mtrx[random_index] = old_theta;
      dE=0;
    }
  }

  *E_tot    += dE;
  *rhoX_sum += rho_s(J, beta, theta_mtrx, rows, cols, arr_all_NN, 2);
  *rhoY_sum += rho_s(J, beta, theta_mtrx, rows, cols, arr_all_NN, 1);
}


static void motecarlo_XY_step_no_data
(double *theta_mtrx, int rows, int cols, int N,
 double J, double beta, double max_flip,
 double *E_tot, 
 int *arr_all_NN 
 ){
  /* Makes a Monte Carlo step in the XY model with out recording
     any data on the spin stiffness.
  */
  int random_index = rand() % N; // random index.
  double old_theta = theta_mtrx[random_index];
  /* Produces a random angle between -<max_flip> to +<max_flip. */
  double d_theta   = max_flip*2*((double)(rand()-RAND_MAX/2))/RAND_MAX;
  theta_mtrx[random_index] += d_theta; //updating theta_mtrx 

  double dE = deltaE(J, theta_mtrx, random_index,
		     old_theta, rows, cols, arr_all_NN);

  /* If dE is positive, we make anonther check to see wheter or not
     to keep the change. This extra check is based on the temperature
     (beta). 
  */
  if ( dE>0 ){ 
    if ( ((double)rand())/RAND_MAX > exp( -beta*dE ) ){
      /* If r is too big, then flip back
	 and the change in energy is 0.
      */
      theta_mtrx[random_index] = old_theta;
      dE=0;
    }
  }

  *E_tot    += dE;
}


/////////////////////////////////////////////////////////////////


int montecarlo_XY_average
(int rows, int cols, 
 double J, double beta, double *return_values,
 int Nsteps, int discard_first, int isP){
  /* This function Monte Carlo simulates a 2D XY model
     and retruns the average and std of the energy, E. 
     The values are retuned in the array <return_values> 
     according to:
                   {T, E, sdtE, rho_sX, rho_sY},
     where T=1/<beta>.

     The XY model is that of a grid of size <rows>*<cols>
     with Hamiltonian:
         H = -J * \sum_{<i,j> NN} cos(theta_i - theta_j),
     at temperature 1/<beta>.
     <Nsteps> of Monte Carlo simulations are performed.

     MAKE SURE THAT <return_values> IS AN ARRAY WITH 
                             5
     ELEMENTS.
  */


  /* Initializations */
  int N           = rows*cols;              //total # sites
  double *Ept     = malloc(sizeof(double));
  double max_flip = 3.14159265358979;       // pi 
  max_flip *=0.1+.5/(beta*J); //Adjusting the maximum spin flip by temperature


  int  arr_all_NN[5*N]; // An array with info on all the NN's.
  // Only periodic BC's will do here.
  get_all_NN_periodic(arr_all_NN, rows, cols);
  


  // Giving the pointers (human) understandable names.
  return_values[0] = 1/beta;
  double *meanE;    meanE    = return_values +1; *meanE    = 0;
  double *stdE;     stdE     = return_values +2; *stdE     = 0;
  double *rhoX_sum; rhoX_sum = return_values +3; *rhoX_sum = 0;
  double *rhoY_sum; rhoY_sum = return_values +4; *rhoY_sum = 0;
 
  
  
  double theta_mtrx[N];                  // init, of random ising grid
  theta_init(theta_mtrx, rows, cols);    // init with random angles


  *Ept = totE(J, theta_mtrx, rows, cols, arr_all_NN);

  for (int a=0; a<discard_first; ++a ){ //for #1
    /* In this loop we just perform the Monte Carlo
       step <discard_first> number of times, to get
       the initial warm-up oscillations out.
    */
    motecarlo_XY_step_no_data
      (theta_mtrx, rows, cols, N,
       J, beta, max_flip,
       Ept, 
       arr_all_NN );

     //Passing 0 as the <current_index> because
     //arr_EM only has 2 elements.
  }//end for #1

  *rhoX_sum = 0;
  *rhoY_sum = 0;
  int Ndata=0;
  for (int b=0; b<Nsteps; ++b ){ //for #2
    /* This is one MC step */

    if (b%N == 0){ //only collect data each Nth step.
      Ndata++;
      motecarlo_XY_step
	(theta_mtrx, rows, cols, N,
	 J, beta, max_flip,
	 Ept, rhoX_sum, rhoY_sum,
	 arr_all_NN );
      
      /* These expressions are not fully the mean and sdt,
	 but they will be modified after the loop.
      */
      *meanE += *Ept; //E;
      *stdE  += (*Ept) * (*Ept);  //E*E;
    }else{//The rest of the turn, we just trod along. 
      motecarlo_XY_step_no_data
	(theta_mtrx, rows, cols, N,
	 J, beta, max_flip,
	 Ept, // E still needs to be updated
	 arr_all_NN );
    }
    
  }//end for #2

  /* Calculating the means */
  *meanE    /= Ndata;
  *rhoX_sum /= Ndata;
  *rhoY_sum /= Ndata;
  /*        Calculating the std's.
     In this stage the *stdX's are just sums of X*X (X^2).

     var = <x^2> - <x>^2 = sum x^2/N - (sum x/N)^2
     std = sqrt( var )
     In the estimates of variance, one should use N-1, however
     in this case N>>1, so it makes no difference.
  */
  *stdE = sqrt( *stdE/Ndata - ( *meanE )*( *meanE ) );

  free(Ept); Ept = NULL;
  return 0;
}