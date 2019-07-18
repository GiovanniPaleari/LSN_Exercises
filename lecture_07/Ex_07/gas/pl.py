import matplotlib
import matplotlib.pyplot as plt
import numpy as np


step, ene = np.loadtxt("inst_energy_gas.txt",usecols=(0,1), delimiter='\t', unpack='true')
press = np.loadtxt("inst_pressure_gas.txt",usecols=(1), delimiter='\t', unpack='true')


plt.plot(step, ene, label='Instant Energy')
plt.plot(step, press, label='Instant Pressure')
#plt.plot(mcstep, inst_m, label='Temperature')


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
