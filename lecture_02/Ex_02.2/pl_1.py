import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x1, y1, error1 = np.loadtxt("Plot_lattice.txt", usecols=(0,1,2), delimiter='\t', unpack='true')
x2, y2, error2 = np.loadtxt("Plot_uniform.txt", usecols=(0,1,2), delimiter='\t', unpack='true')

f, axarr = plt.subplots(1,2)

axarr[0].errorbar(x1,y1,yerr=error1)
axarr[0].set_title('Lattice Random Walk')
axarr[1].errorbar(x2,y2,yerr=error2)
axarr[1].set_title('Uniform Random Walk')

#plt.xscale('log')
#plt.yscale('log')

plt.show()
