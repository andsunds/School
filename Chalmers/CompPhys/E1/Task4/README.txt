################################################################################
# FKA121/FIM540 Computational Physics 2015
# C programs for exercise E1
#
# Martin Gren, gmartin@chalmers.se 
################################################################################
The following commands will build the executable powerspectrum:

	gcc -c -o fft_func.o fft_func.c -O3 -W
	gcc -c -o main.o main.c -O3 -W
	gcc -o powerspectrum fft_func.o main.o -O3 -W -lm -lgsl -lgslcblas

Already with two .c files the compilation starts to become a bit messy.
One therefore uses makefiles that automatically executes the commands 
needed to compile the program. If the makefile is called Makefile or 
makefile one can type the command 'make' to compile the code. 

To run the program type './powerspectrum' in the terminal. 

The data is dumped into files with suffix '.dat'. The programs plot_*.py 
and plot_*.m plots the powerspectrum in Python and MATLAB respectively.   
 
(There is a bug in matplotlib (python) on the StuDAT computers. 
You get the following warning:
	/usr/lib64/python2.6/site-packages/matplotlib/backends/backend_gtk.py:621: 
	DeprecationWarning: Use the new widget gtk.Tooltip
	 self.tooltips = gtk.Tooltips()
This warning can be ignored.
)
