import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x, f, error = np.loadtxt("Plot_1.txt", usecols=(0,1,2), delimiter='\t', unpack='true')
plt.errorbar(x,f,yerr=error)
plt.xlabel('# blocks')
plt.ylabel(r'$\pi$')
plt.grid(True)

plt.plot([1, 50, 100], [3.1415, 3.1415, 3.1415], 'r')

plt.show()
