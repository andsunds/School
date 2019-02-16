# plot the powerspectrum
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pylab as plt
import numpy as np

# input file
filename = 'powerspectrum.dat'

# import data
data = np.loadtxt(filename)

# initial size of plot window
plt.figure(figsize=(8,6))

# plot
plt.plot(data[:,0], data[:,1],'-')

# labels
plt.xlabel('Frequency / [arb. unit]', fontsize=20)
plt.ylabel('Power spectrum / [arb. unit]', fontsize=20)

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)
plt.xlim({-0.5, 0.5})

# display the plot

plt.savefig('powerspectrum.pdf')
plt.show()
