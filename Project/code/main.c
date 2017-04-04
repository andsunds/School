/* N.B.
   In this implementation all matrices will be
   represented as a 1D array of length rows*cols.
*/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "main.h"


static int simulate_ising_write_all_to_bin
( int L, double J, double beta0, double d_beta, 
  int Nsims, int Nsteps, 
  int isP, FILE *logPTR, char *save_directory ){
  /* As the name suggests, simulates a bunch
     of differnt times and writes all data to
     a binary file.
  */
  //srand((unsigned) time(NULL)); //seeds rand, should only be called once

  //int    L        = 16;         // side length of the grid
  //double J        = 1.0;        // energy factor in the Hamlitonian
  //double beta0    = 0.25;       // invers temperature
  //double d_beta   = 16e-3;      // step size in beta
  double beta;                  // init
  //int Nsims       = 128;        // # sims, set to 0 if not active
  //int Nsteps      = (int)2e7;   // # Monte Carlo steps
  int write_chunk = 256;        // # doubles to be written at once
  
  clock_t begin2, begin3, end2; //init, for time tracking
  float runtime2, runtime3;  


  /* Lots of differents messages */
  char startMSG[128];
  sprintf(startMSG,"\n   Starting %d simulation(s) of %2g steps each (write all).\n \n", Nsims, (float)Nsteps);
  
  char sim_errorMSG[] = " ERROR in simulate_ising_write_all_to_tbin(): error while simulating in sim number ";

  char timeMSG[128]; //rewritten after each iteration



  /* Pre-setting a beta value.
    for (int p=0; p<26; ++p )
    beta = beta/factor;
  */

  /* A message before we start */
  printf(        "%s",startMSG);
  fprintf(logPTR,"%s",startMSG);


  /* Starting the computations */
  begin3 = clock();
  for ( int i=0; i<Nsims; ++i ){
    beta = beta0 + i*d_beta;
    
    /* A timer for how long the execution took */
    begin2 = clock();

    /* Here is the actual calculation */
    int isOK = montecarlo_ising_full
      (L,L, J, beta, Nsteps, isP,
       write_chunk, save_directory, logPTR);
    if ( isOK != 0 ){
      /* Just a small check to see that nothing went wrong
	 in the last simulation.
      */
      printf(         "%s%d", sim_errorMSG,i);
      fprintf(logPTR, "%s%d", sim_errorMSG,i);
      return 1;
    }

    /* Time tracking */
    end2 = clock();
    runtime2= (float)(end2 - begin2) / CLOCKS_PER_SEC;
    runtime3= (float)(end2 - begin3) / CLOCKS_PER_SEC / 60;//min

    if ( runtime3 < 60 ){ //diplay accumulated time in min's
      sprintf(timeMSG,"Execution time (beta=%0.4f, %d of %d): %0.3f s. ( %0.2f min )\n", beta, i+1, Nsims, runtime2, runtime3);
    }else{ //diplay accumulated time in hours
      sprintf(timeMSG,"Execution time (beta=%0.4f, %d of %d): %0.3f s. ( %0.2f hours )\n", beta, i+1, Nsims, runtime2, runtime3/60);
    }
      printf(        "%s",timeMSG);
      fprintf(logPTR,"%s",timeMSG);
  } // end for
  return 0;
}








static int simulate_ising_write_avg_to_tsv
( int L, double J, double T0, double dT, 
  int Nsims, int Nsteps, int disc_first,
  int isP, FILE *logPTR, char *save_directory ){
  /* As the name suggests, simulates a bunch of
     differnt times and then oly writes the averages
     and std's of E and M to a tsv-file.

     <logPTR> is a pointer to a log file.
  */
  //srand((unsigned) time(NULL)); //seeds rand, should only be called once

  /* Inits */
  int no_of_values = 5;         // do NOT change this! 
  double TEMarr[no_of_values];  // array with mean and std of E & M.
  //int    L       = 16;          // side length og the grid
  //double J       = 1.0;         // energy factor in the Hamlitonian
  //double beta0   = 0.40;        // invers temperature
  //double d_beta  = 2.5e-5;      // step size in beta
  double beta;                  // init
  //int Nsims      = 2048;        // # sims, set to 0 if not active
  //int Nsteps     = (int)2e7;    // # Monte Carlo steps
  //int disc_first = (int)5e4;    // discard values from warm-up period

  int isOK;                     //init, for error checking
  clock_t begin2,end2, begin3;  //init, for time tracking
  float runtime2, runtime3;     //init, for time tracking


  /* Lots of differents messages */
  char startMSG[128];
  sprintf(startMSG,"\n   Starting %d simulation(s) of %2g steps each (write avg).\n \n", Nsims, (float)Nsteps);

  char file_errorMSG[] = " ERROR in simulate_ising_write_avg_to_tsv(): unable to open file: ";

  char sim_errorMSG[] = " ERROR in simulate_ising_write_avg_to_tsv(): error while simulating in sim number ";

  char timeMSG[128]; //rewritten after each iteration

  char endMSG[256]; //written later


  /*    file I/O    */
  char filename[128];
  //sprintf(filename,"%sDEBUG_.tsv",save_directory); //DEBUG
  if ( isP == +1 )
    sprintf(filename, "%sEstdEMstdM_T_%0.3f-%0.3f_%d_L-%d_PERIODIC.tsv",
	    save_directory, T0, T0+(Nsims-1)*dT, Nsims,L);
  else 
    sprintf(filename, "%sEstdEMstdM_T_%0.3f-%0.3f_%d_L-%d.tsv",
	    save_directory, T0, T0+(Nsims-1)*dT, Nsims,L);


  FILE *filePTR;
  filePTR = fopen(filename,"w");
    /* Some error checking if the file opened correctly */
  if ( !filePTR ){
    printf(         "%s%s\n", file_errorMSG,filename);
    fprintf(logPTR, "%s%s\n", file_errorMSG,filename);
    return 1; 
  }


  /* A message before we start */
  printf(        "%s",startMSG);
  fprintf(logPTR,"%s",startMSG);
  printf("Writing data to: %s\n",filename);

  /* Starting the computations */
  begin3 = clock();
  for ( int i=0; i<Nsims; ++i ){
    /* A timer for how long the execution took */
    begin2 = clock();
    /* Here is the actual calculation */
    beta = 1/(T0 + i*dT);
    isOK = montecarlo_ising_average
      (L,L, J, beta, TEMarr, Nsteps, disc_first, isP);
    if ( isOK != 0 ){
      /* Just a small check to see that nothing went wrong
	 in the last simulation.
      */
      printf(         "%s%d", sim_errorMSG,i);
      fprintf(logPTR, "%s%d", sim_errorMSG,i);
      return 1;
    }
    /* Write to file */
    for ( int j=0; j<no_of_values; ++j ){
      fprintf(filePTR, "%16e",TEMarr[j]);
    }
    fprintf(filePTR, "\n"); //newline after the simulation

    /* Time tracking */
    end2 = clock();
    runtime2= (float)(end2 - begin2) / CLOCKS_PER_SEC;
    runtime3= (float)(end2 - begin3) / CLOCKS_PER_SEC / 60;//min

    if ( runtime3 < 60 ){ //diplay accumulated time in min's
      sprintf(timeMSG,"Execution time (beta=%0.4f, %d of %d): %0.3f s. ( %0.2f min )\n", beta, i+1, Nsims, runtime2, runtime3);
    }else{ //diplay accumulated time in hours
      sprintf(timeMSG,"Execution time (beta=%0.4f, %d of %d): %0.3f s. ( %0.2f hours )\n", beta, i+1, Nsims, runtime2, runtime3/60);
    }
    printf(        "%s",timeMSG);
    fprintf(logPTR,"%s",timeMSG);
  } // end for

  sprintf(endMSG,"\nSimulations done, data written to:\n   %s\n",
	  filename);
  printf(        "%s",endMSG);
  fprintf(logPTR,"%s",endMSG);

  fclose(filePTR); 
  return 0;
}






static int simulate_XY_write_avg_to_tsv
( int L, double J, double T0, double dT, 
  int Nsims, int Nsteps, int disc_first,
  int isP, FILE *logPTR, char *save_directory ){
  /* As the name suggests, simulates a bunch of
     differnt times and then oly writes the averages
     and std's of E and M to a tsv-file.

     <logPTR> is a pointer to a log file.
  */
  //srand((unsigned) time(NULL)); //seeds rand, should only be called once

  /* Inits */
  int no_of_values = 5;           // do NOT change this! 
  double TErho_arr[no_of_values]; // array with mean and std of E & M.
  //int    L       = 16;            // side length og the grid
  //double J       = 1.0;           // energy factor in the Hamlitonian
  //double beta0   = 0.72;          // invers temperature
  //double d_beta  = 2e-3;          // step size in beta
  double beta;                    // init
  //int Nsims      = 256;           // # sims, set to 0 if not active
  //int Nsteps     = (int)1e6;      // # Monte Carlo steps
  //int disc_first = (int)5e4;      // discard values from warm-up period

  int isOK;                       //init, for error checking
  clock_t begin2,end2, begin3;    //init, for time tracking
  float runtime2, runtime3;       //init, for time tracking


  /* Lots of differents messages */
  char startMSG[128];
  sprintf(startMSG,"\n   Starting %d simulation(s) of %2g steps each (write avg, XY).\n \n", Nsims, (float)Nsteps);

  char file_errorMSG[] = " ERROR in simulate_XY_write_avg_to_tsv(): unable to open file: ";

  char sim_errorMSG[] = " ERROR in simulate_XY_write_avg_to_tsv(): error while simulating in sim number ";

  char BC_errorMSG[] = "ERROR: Boundary conditions not periodic. For XY model, need periodic BC.";

  char timeMSG[128]; //rewritten after each iteration

  char endMSG[256]; //written later


  /*    file I/O    */
  char filename[128];
  //sprintf(filename,"%sDEBUG_.tsv",save_directory); //DEBUG
  if ( isP == +1 )
    sprintf(filename, "%sTEstdErhoXY_T_%0.3f-%0.3f_%d_L-%d_PERIODIC.tsv",
	    save_directory, T0, T0+(Nsims-1)*dT, Nsims,L);
  else {
    printf(         "%s\n", BC_errorMSG);
    fprintf(logPTR, "%s\n", BC_errorMSG);
    return 1;
  }

  FILE *filePTR;
  filePTR = fopen(filename,"w");
    /* Some error checking if the file opened correctly */
  if ( !filePTR ){
    printf(         "%s%s\n", file_errorMSG,filename);
    fprintf(logPTR, "%s%s\n", file_errorMSG,filename);
    return 1; 
  }


  /* A message before we start */
  printf(        "%s",startMSG);
  fprintf(logPTR,"%s",startMSG);
  printf("Writing data to: %s\n",filename);

  /* Starting the computations */
  begin3 = clock();
  for ( int i=0; i<Nsims; ++i ){
    /* A timer for how long the execution took */
    begin2 = clock();
    /* Here is the actual calculation */
    beta = 1/(T0 + i*dT);
    isOK = montecarlo_XY_average
      (L,L, J, beta, TErho_arr, Nsteps, disc_first, isP);
    if ( isOK != 0 ){
      /* Just a small check to see that nothing went wrong
	 in the last simulation.
      */
      printf(         "%s%d", sim_errorMSG,i);
      fprintf(logPTR, "%s%d", sim_errorMSG,i);
      return 1;
    }
    /* Write to file */
    for ( int j=0; j<no_of_values; ++j ){
      fprintf(filePTR, "%16e",TErho_arr[j]);
    }
    fprintf(filePTR, "\n"); //newline after the simulation

    /* Time tracking */
    end2 = clock();
    runtime2= (float)(end2 - begin2) / CLOCKS_PER_SEC;
    runtime3= (float)(end2 - begin3) / CLOCKS_PER_SEC / 60;//min

    if ( runtime3 < 60 ){ //diplay accumulated time in min's
      sprintf(timeMSG,"Execution time (beta=%0.4f, %d of %d): %0.3f s. ( %0.2f min )\n", beta, i+1, Nsims, runtime2, runtime3);
    }else{ //diplay accumulated time in hours
      sprintf(timeMSG,"Execution time (beta=%0.4f, %d of %d): %0.3f s. ( %0.2f hours )\n", beta, i+1, Nsims, runtime2, runtime3/60);
    }
    printf(        "%s",timeMSG);
    fprintf(logPTR,"%s",timeMSG);
  } // end for

  sprintf(endMSG,"\nSimulations done, data written to:\n   %s\n",
	  filename);
  printf(        "%s",endMSG);
  fprintf(logPTR,"%s",endMSG);

  fclose(filePTR); 
  return 0;
}




int main(){
  int isOK = 0; // init (0 = ok, 1 = not ok)

  srand((unsigned) time(NULL)); //seeds rand, should only be called once

  
  int L, Nsims, Nsteps, disc_first;
  double J, T0, dT, beta0, d_beta;
 
  J          = 1.0;         // energy factor in the Hamlitonian
  disc_first = (int)5e4;    // discard values from warm-up period

  



  //char save_directory[] = "../data/Ising/";
  char save_directory[] = "../data/XY/"; 
  //char save_directory[] = "../data/"; 
  //char save_directory[] = "../data/DEBUG/"; //DEBUG

  /* Variable to indicate wheter or not to use periodic BC's.
     isP == +1, means to use periodic BC's, otherwise not.
  */
  int isP = +1; 


  /* Opens a log file with runtimes and possible errors */
  FILE *logPTR;
  logPTR = fopen("sim.log","w");
  if ( !logPTR )
    return 1; // Sumting went wong...


  /* Timings */
  clock_t begin1, end1; //init for time tracking
  float time_spent, sec;
  int hour, min;
  char timeMSG[128];
  
  begin1 = clock();


  /* Values for the Ising simulations */
  L          = 16;          // side length of the grid
  Nsteps     = (int)1e7;    // # Monte Carlo steps
  
  /*
  // Write all 
  Nsims  = 2;  // # sims, set to 0 if not active
  beta0  = 1;
  d_beta = 2;
  isOK = simulate_ising_write_all_to_bin
    ( L, J, beta0,d_beta, Nsims, Nsteps,
      isP, logPTR, save_directory );
  */
  
  
  /*
  // Write avg 
  T0         = 1.40;        // invers temperature
  dT         = 1e-3;        // step size in beta
  Nsims      = 2048;        // # sims, set to 0 if not active
  for (int i=0; i<4; ++i ){
    // Do the simulations for L=12, 16, 20, 24 
    L = 12 +4*i;
    isOK = simulate_ising_write_avg_to_tsv 
      ( L, J, T0, dT, Nsims, Nsteps, disc_first,
	isP, logPTR, save_directory );
  }
  
  fprintf(logPTR,"%s", "Starting XY.\n");
  printf(        "%s", "Starting XY.\n");
  */
  // Values for the XY simulations 
  T0         = .01;         // temperature
  dT         = 3e-2;        // step size in beta
  Nsims      = 64;          // # sims, set to 0 if not active
  Nsteps     = (int)5e4;    // # Monte Carlo steps
  disc_first = (int)2e3;
  L = 32;
  isOK = simulate_XY_write_avg_to_tsv
    ( L, J, T0, dT, Nsims, Nsteps*L*L+1, disc_first*L*L,
      isP, logPTR, save_directory );
  ///*
  T0         = .9;      // temperature
  dT         = 1e-2;    // step size in beta
  Nsims      = 6;       // # sims, set to 0 if not active
  L=8;
  for (int i=0; i<3;++i){
    L*=2;
    isOK = simulate_XY_write_avg_to_tsv
      ( L, J, T0, dT, Nsims, Nsteps*L*L+1, disc_first*L*L,
  	isP, logPTR, save_directory );
  }
  //*/



  /* Time tracking and logging */
  end1 = clock();
  time_spent = (float)(end1 - begin1) / CLOCKS_PER_SEC; 
  hour = floor(time_spent/3600);
  min  = floor(time_spent/60) - 60*hour; 
  sec  = time_spent-3600*hour-60*min;

  sprintf(timeMSG,"\nTotal execution time: %dh %02dm %0.3fs. \n\n",  hour, min, sec);
  fprintf(logPTR,"%s", timeMSG);
  printf(        "%s", timeMSG);

  fclose(logPTR);
  return isOK;
}
