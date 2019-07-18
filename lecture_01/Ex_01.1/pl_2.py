import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x, f, error = np.loadtxt("Plot_2.txt", usecols=(0,1,2), delimiter='\t', unpack='true')
plt.errorbar(x,f,yerr=error)
plt.xlabel('# blocks')
plt.ylabel('<(r-1/2)^2> - 1/12')
plt.grid(True)

plt.show()
