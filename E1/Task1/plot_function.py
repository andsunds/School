# plot the h(t)
# Created by Martin Gren 2014-10-25.

# imports
import matplotlib.pyplot as plt
import numpy as np

# input file
filename = 'function.dat'

# import data
data = np.loadtxt(filename)

# initial size of plot window
plt.figure(figsize=(8,6))

# plot
plt.plot(data[:,0], data[:,1],'-')

# labels
plt.xlabel('t / [arb. unit]', fontsize=20)
plt.ylabel('h(t) / [arb. unit]', fontsize=20)

# set tick fontsize
plt.yticks(fontsize=12)
plt.xticks(fontsize=12)

# display the plot

plt.savefig('h(t).pdf')
plt.show()
