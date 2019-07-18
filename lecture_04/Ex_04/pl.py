import matplotlib
import matplotlib.pyplot as plt
import numpy as np


temp = np.loadtxt("output_temp.dat", delimiter='\t', unpack='true')
ekin = np.loadtxt("output_ekin.dat", delimiter='\t', unpack='true')
epot = np.loadtxt("output_epot.dat", delimiter='\t', unpack='true')
etot = np.loadtxt("output_etot.dat", delimiter='\t', unpack='true')
pr = np.loadtxt("output_pr.dat", delimiter='\t', unpack='true')

plt.plot(temp, label='Temperature')
plt.plot(ekin, label='Kinetic Energy')
plt.plot(epot, label='Potential Energy')
plt.plot(etot, label='Total Energy')
plt.plot([1.,1300],[0.8,0.8])
plt.plot(pr, label='Pressure')

plt.legend(loc='upper right')


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
