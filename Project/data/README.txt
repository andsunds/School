
The data files contains data from simulations of the Ising model:
             H = -J * sum_{<i,j> are NN} s_i * s_j.   

This is unless their file names start with XY, then the model used 
were:
       H = -J * sum_{<i,j> are NN} cos( theta_i - theta_j).   

Unless otherwise stated the simulations were done on a 16x16 grid.


======================================================================

The .bin file names are formated 
                  EM_beta_%0.5f[_PERIODIC].bin,
where the value is the beta value; if periodic BC's were used, then 
the file name also has "_PERIODIC" in the end.   

The .bin files contains values of E and M as doubles. The value of 
J used was 1, therefore the units of T and E are in terms of J.


======================================================================
        
The .tsv file names are formated   
EstdEMstdM_beta_%0.3f-%0.3f_%d[_PERIODIC].tsv,   
where the first value is the lowest beta value and the second is the 
higest, followed by the number of steps; if periodic BC's were used, 
then the file name also has "_PERIODIC" in the end.   

The files are formated as columns with T, E, stdE, M, and stdM. The
value of J used was 1, therefore the units of T and E are in terms of 
J.


======================================================================

Up to 2017-03-26, the files with periodic BC's are a litte of in their
energies, E. This is due to a miss in the calculation of the initial 
energy, where hard BC's were accidentally used.


