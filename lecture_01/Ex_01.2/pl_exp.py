import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x1, y1 = np.loadtxt("histo_N=1.txt", usecols=(3,4), delimiter='\t', unpack='true')
x2, y2 = np.loadtxt("histo_N=2.txt", usecols=(3,4), delimiter='\t', unpack='true')
x3, y3 = np.loadtxt("histo_N=10.txt", usecols=(3,4), delimiter='\t', unpack='true')
x4, y4 = np.loadtxt("histo_N=100.txt", usecols=(3,4), delimiter='\t', unpack='true')

f, axarr = plt.subplots(2,2)

axarr[0, 0].bar(x1, y1, width=0.03)
axarr[0, 0].set_title(r'$N=1$')
axarr[0, 1].bar(x2, y2, width=0.03)
axarr[0, 1].set_title(r'$N=2$')
axarr[1, 0].bar(x3, y3, width=0.03)
axarr[1, 0].set_title(r'$N=10$')
axarr[1, 1].bar(x4, y4, width=0.03)
axarr[1, 1].set_title(r'$N=100$')



f.suptitle('Exponential dice')

plt.show()
