import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x1, y1, error1 = np.loadtxt("Put_1step.txt", usecols=(0,1,2), delimiter='\t', unpack='true')
x2, y2, error2 = np.loadtxt("Put_100step.txt", usecols=(0,1,2), delimiter='\t', unpack='true')


f, axarr = plt.subplots(1,2,sharey=True)

axarr[0].errorbar(x1,y1,yerr=error1)
axarr[0].set_title('1 step')
axarr[0].plot([1,50,100], [5.4595325819072364, 5.4595325819072364, 5.4595325819072364])
axarr[0].set(xlabel='# blocks')
axarr[0].set(ylabel=r'$P[S(0),0]$')

axarr[1].errorbar(x2,y2,yerr=error2)
axarr[1].plot([1,50,100], [5.4595325819072364, 5.4595325819072364, 5.4595325819072364])
axarr[1].set_title('100 step')
axarr[1].set(xlabel='# blocks')

plt.show()
