# -*- coding: utf-8 -*-

# plot the displacements
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# load input file
dataE = np.loadtxt('Emat.dat')

# initial size of plot window

plt.figure(figsize=(8,6))

# plot
for i in range(1,6):
	plt.plot(dataE[:,0], dataE[:,i],'-',label=('mode ' + str(i) ))



# labels
plt.xlabel('Time / [dimensionless]', fontsize=20)
plt.ylabel('E / [dimensionless]', fontsize=20)

# legend
plt.legend()
leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12) 

# axis limits
# plt.ylim([0,100])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# display the plot

plt.savefig('Emat.pdf')

#######################################

# load input file
dataE = np.loadtxt('Eavg.dat')

# initial size of plot window

plt.figure(figsize=(8,6))

# plot

plt.plot(dataE[:,0], dataE[:,1:],'-')



# labels
plt.xlabel('Time / [dimensionless]', fontsize=20)
plt.ylabel('E / [dimensionless]', fontsize=20)

# legend
#plt.legend()
#leg = plt.gca().get_legend()
ltext  = leg.get_texts()
plt.setp(ltext, fontsize=12) 
plt.gca().set_yscale('log') 
plt.gca().set_xscale('log')

# axis limits
plt.ylim([1e-12,40])
plt.xlim([1e2,1e6])

# tick fontsize
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# display the plot

plt.savefig('Eavg.pdf')

plt.show()

