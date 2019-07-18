import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x, f, error = np.loadtxt("Plot_1.txt", usecols=(0,1,2), delimiter='\t', unpack='true')
plt.errorbar(x,f,yerr=error)
plt.xlabel('Number of blocks')
plt.ylabel('<r> - 0.5')
plt.grid(True)

plt.show()
