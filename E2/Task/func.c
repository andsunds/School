/*
func.c
 
Created by AL on 2013-10-24.
Further developed by Martin Gren on 2015-10-23.
*/

#define PI 3.141592653589
#include <math.h>
/*
Function that calculates the acceleration based on the Hamiltonian.
The acceleration is calculated based on the displacements u and then stored in a.
u and a should be vectors of the same size, N_u
*/
void calc_acc(double *a, double *u, double m, double kappa, double alpha, int N_u)
{
    /* Declaration of variables */
    int i;
    
    /* Calculating the acceleration on the boundaries */
    /* EOM : m * a[i] = kappa * (u[i+1] - 2u[i] + u[i-1] ) + ...
    		alpha * [ (u[i] - u[i+1])^2 - (u[i] - u[i-1])^2 ] */
    
    a[0] = kappa*(- 2*u[0] + u[1])/m 
    	+ alpha/m * ( - u[0]*u[0] + (u[1]-u[0])*(u[1]-u[0]) );
    	
    a[N_u - 1] = kappa*(u[N_u - 2] - 2*u[N_u - 1])/m +
    	+ alpha/m * ( -(u[N_u-1] - u[N_u-2]) * (u[N_u-1] - u[N_u-2]) 
    	+ u[N_u-1]*u[N_u-1] );
    
    /* Calculating the acceleration of the inner points */
    for (i = 1; i < N_u - 1; i++){
        a[i] = kappa*(u[i - 1] - 2*u[i] + u[i + 1])/m 
        + alpha/m * ( - (u[i] - u[i-1]) * (u[i] - u[i-1]) 
    	   + (u[i+1] - u[i]) * (u[i+1] - u[i]) );
    }
}

/* Function that calculates the potential energy based on the displacements */
double calc_pe(double *u, double kappa, double alpha, int N_u)
{
    /* Declaration of variables */
    int i;
    double e = 0;
    double tmp_du;
    
    /* Calculating the energy on the boundaries */
    
    tmp_du = u[0];
    e += kappa*tmp_du*tmp_du/2 + alpha/3 * tmp_du * tmp_du * tmp_du;
    
    tmp_du = - u[N_u -1]; //u_N = 0;
    e += kappa*tmp_du*tmp_du/2 + alpha/3 * tmp_du * tmp_du * tmp_du;
    
    /* Calculating the energy of the inner points */
    
    for (i = 0; i < N_u - 1; i++){
    	tmp_du = u[i+1] - u[i];
        e += kappa*tmp_du*tmp_du/2 + alpha/3 * tmp_du * tmp_du * tmp_du;
    }
    
    return e;	
}


/* Function that calculates and returns the kinetic energy based on the velocities and masses */
double calc_ke(double *p, double m, int size_of_p)
{
    /* Declaration of variables */
    int i;
    double e = 0; 
    /* Calculating the energy of the inner points */
    for (i = 0; i < size_of_p; i++){
        e += (p[i])*(p[i])/(2*m);
    }
    return e;	
}

// note that if m=/=1, we need some normalization factor!
void construct_trans_matrix(int N_part, double trans_matrix[N_part][N_part])
{
    double factor = 1 / ((double) N_part + 1);

    for (int i=0; i < N_part; i++) {
        for (int j=0; j < N_part; j++) {
            trans_matrix[i][j] = sqrt(2 * factor) * 
                sin((j + 1) * (i + 1) * PI * factor);
        }
    }
}    


    /* Transformation to normal modes Q from displacements q.  */
void get_normal_modes(double *Q, double *q, int N_part, 
    double trans_matrix[N_part][N_part])
{
    double sum;
    for (int i = 0; i < N_part; i++){
        sum = 0;
        for (int j = 0; j < N_part; j++){
            sum += q[j] * trans_matrix[i][j];
        }
        Q[i] = sum;
    }
}

void calc_mode_energy(double *E, double  *Q, double *P, int N_part, 
	double omega0)
{
	double omega_k;// = 0;
	double factor = PI * 0.5 / ((double) N_part + 1);
	
	for (int k=1; k<N_part+1; k++){
		omega_k = 2 * omega0 * sin ( k * factor);
		E[k-1] = 0.5 * ( P[k-1]*P[k-1] + omega_k*omega_k * Q[k-1]*Q[k-1]);
	}
}	
	
