/*
func.h
 
Created by AL on 2013-10-24.
Further developed by Martin Gren on 2015-10-23.
*/

#ifndef _func_h
#define _func_h

extern void calc_acc(double *, double *, double, double, double alpha, int);

extern double calc_pe(double *, double, double alpha, int);

extern double calc_ke(double *, double, int);

extern void get_normal_modes(double *Q, double *q, int N_part, 
	double trans_matrix[N_part][N_part]);

extern void construct_trans_matrix( int N_part, double trans_matrix[N_part][N_part]);

extern void calc_mode_energy(double *E, double  *Q, double *P, int N_part, 
	double omega0);
#endif
