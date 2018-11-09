################################################################################
# FKA121/FIM540 Computational Physics 2015
# C programs for exercise E1
#
# Martin Gren, gmartin@chalmers.se 
################################################################################
If you have installed the gcc compiler the program can be compiled with the 
following command:

			gcc -o run main.c -O3 -W -lm

which will create an executable file called "run".

To run the program type './run' in the terminal. 

The data is dumped into files with suffix '.dat'. The programs plot_*.py and 
plot_*.m plots the function in Python and MATLAB respectively.   

(There is a bug in matplotlib (python) on the StuDAT computers. 
You get the following warning:
	/usr/lib64/python2.6/site-packages/matplotlib/backends/backend_gtk.py:621: 
	DeprecationWarning: Use the new widget gtk.Tooltip
	 self.tooltips = gtk.Tooltips()
This warning can be ignored.
)
