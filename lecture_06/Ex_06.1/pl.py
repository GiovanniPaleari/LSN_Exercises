import matplotlib
import matplotlib.pyplot as plt
import numpy as np


temp, u, err1 = np.loadtxt("Gibbs_MofT.txt",usecols=(0,1,2), delimiter='\t', unpack='true')


plt.errorbar(temp, u, yerr=err1)
#plt.plot(mcstep, inst_m, label='Temperature')


#plt.legend(loc='upper right')


#f, axarr = plt.subplots(1,2,sharey=True)

#axarr[0].errorbar(x1,y1,yerr=error1)
#axarr[0].set_title('1 step')
#axarr[0].plot([1,50,100], [14.975790778311286, 14.975790778311286, 14.975790778311286])
#axarr[0].set(xlabel='# blocks')
#axarr[0].set(ylabel=r'$C[S(0),0]$')

#axarr[1].errorbar(x2,y2,yerr=error2)
#axarr[1].plot([1,50,100], [14.975790778311286, 14.975790778311286, 14.975790778311286])
#axarr[1].set_title('100 step')
#axarr[1].set(xlabel='# blocks')



#plt.xscale('log')
#plt.yscale('log')

plt.show()
