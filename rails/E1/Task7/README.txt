################################################################################
# FKA121/FIM540 Computational Physics 2015
# C programs for exercise E1
#
# Martin Gren, gmartin@chalmers.se 
################################################################################
Once gcc, make, and gsl are installed you can compile the code. 
Type 'make' to compile the code. 

To run the program type './run' in the terminal. 

The displcement data is dumped into 'displacements.dat'. The programs 
plot_displacements.py and plot_displacements.m plots the displacements in Python 
and MATLAB respectively.   

There are corresponding plotting programs for the energies and powerspectrum 
as well. In Task 6 you need to dump the energies and power spectrum into the 
.data files given in the plotting programs to make use of them. 

Remember to change the prototypes in the header files if you make changes to the
functions.
 
(There is a bug in matplotlib (python) on the StuDAT computers. 
You get the following warning:
	/usr/lib64/python2.6/site-packages/matplotlib/backends/backend_gtk.py:621: 
	DeprecationWarning: Use the new widget gtk.Tooltip
	 self.tooltips = gtk.Tooltips()
This warning can be ignored.
)
