/*
E1_func.c
 
Created by AL on 2013-10-24.
Further developed by Martin Gren on 2015-10-23.
*/


/*
Function that calculates the acceleration based on the Hamiltonian.
The acceleration is calculated based on the displacements u and then stored in a.
u and a should be vectors of the same size, size_u
*/
void calc_acc(double *a, double *u, double *m, double kappa, int size_u)
{
    /* Declaration of variables */
    int i;
    
    /* Calculating the acceleration on the boundaries */
    /* No factors of 2 here since now, the ends are free. */
    a[0] = kappa*(-u[0] + u[1])/m[0]; 
    a[size_u - 1] = kappa*(u[size_u-2] - u[size_u-1])/m[size_u-1];
    
    /* Calculating the acceleration of the inner points */
    for (i = 1; i<size_u-1; i++){
        a[i] = kappa*(u[i-1] - 2*u[i] + u[i+1])/m[i];
    }
}

/* Function that calculates the potential energy based on the displacements */
double calc_pe(double *u, double kappa, int size_u)
{
    /* Declaration of variables */
    int i;
    double e = 0;
    /* No need for boundary cases. */
    
    /* Calculating the energy of the inner springs. */
    for (i = 0; i < size_u - 1; i++){
        e += kappa*(u[i] - u[i + 1])*(u[i] - u[i + 1])/2;
    }
    return e;	
}


/* Function that calculates and returns the kinetic energy based on the velocities and masses */
double calc_ke(double *v, double *m, int size_v)
{
    /* Declaration of variables */
    int i;
    double e = 0; 
    /* Calculating the energy of the inner points */
    for (i = 0; i < size_v; i++){
        e += m[i]*(v[i])*(v[i])/2;
    }
    return e;	
}
