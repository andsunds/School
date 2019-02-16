
#include <stdio.h>
#include <stdlib.h>

int main()
{
  int i,k,j, nbr_of_lines;
  FILE *in_file;
  
  nbr_of_lines = 1e6; /* The number of lines in MC.txt. */
  double *data = malloc((nbr_of_lines) * sizeof (double));
  
  double f_mean = 0;
  double f_squared_mean = 0;
  double f_var;
  int N_k = 1000;
  double *phi = malloc((N_k) * sizeof (double));
  
  /* Read data from file. */
  in_file = fopen("MC.txt","r");
  for (i=0; i<nbr_of_lines; i++) {
    fscanf(in_file,"%lf",&data[i]);
    f_mean += data[i];
    
  }
  fclose(in_file);
  
  f_mean *= (double)1/nbr_of_lines;
  
  //printf("f_mean = %.2e, Var[f], %.2e \n", f_mean, f_var);
  
  // shift the data with mean around zero.
  for (i=0; i<nbr_of_lines; i++) {
    data[i] -= f_mean;
    f_squared_mean += data[i]*data[i];
  }
  f_mean = 0;
  f_squared_mean *= (double)1/nbr_of_lines;
  f_var = f_squared_mean - f_mean*f_mean;
  
  for (k=0; k<N_k; k++) {
    phi[k] = 0;
    for (i=0; i<nbr_of_lines-k; i++) {
      phi[k] += data[i]*data[i+k];
    }
    phi[k] = (phi[k]/(nbr_of_lines-k) - f_mean*f_mean)/
      (f_squared_mean- f_mean*f_mean);
    //printf("phi = %.2f, k = %d \n", phi[k], k);
  }
  
  // block average
  int block_size;
  double *var_F = malloc((N_k) * sizeof (double));
  double Fj;
  int number_of_blocks;
  for (k=0; k<N_k; k++) { // block size look
    block_size = k+1;
    number_of_blocks = nbr_of_lines/block_size;
    for (j=0; j<number_of_blocks; j++) {// loop over all blocks  
      Fj = 0;
      //Fj_squared = 0;
      for (i=0; i<block_size; i++) {// internal block loop
        Fj += data[j*block_size + i];
        //Fj_squared += data[j*block_size + i]*data[j*block_size + i];
      }
      Fj *= 1/(double)block_size; // these are the values we need the variance of
      var_F[k] += Fj*Fj; // will become the variance soon
    }
    var_F[k] = var_F[k]/number_of_blocks - f_mean*f_mean;
    var_F[k] *= block_size/f_var;
  }
  
  
  FILE *out_file;
  out_file = fopen("MC_analysis.tsv","w");
  for (k=0; k<N_k; k++){
    fprintf(out_file, "%d\t%.8f\t%.8f\n",k,phi[k], var_F[k]);
  }
  fclose(out_file);
  
  
  
  free(data); data = NULL;
  free(phi); phi = NULL;
}



