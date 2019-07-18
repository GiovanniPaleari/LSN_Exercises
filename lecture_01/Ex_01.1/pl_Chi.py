import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x, f = np.loadtxt("Plot_Chi.txt", usecols=(0,1), delimiter='\t', unpack=True)
plt.scatter(x,f)
#plt.xlabel('# blocks')
plt.ylabel(r'$\chi^2$')
plt.grid(True)
plt.ylim((0,200))

plt.plot([1, 50, 100], [100, 100, 100], 'r')

plt.show()
