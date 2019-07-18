import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def f2(mu, sigma, c, x):
    return c*(np.exp(-0.5*((x-mu)/sigma)**2)+np.exp(-0.5*((x+mu)/sigma)**2))**2

#iblk, ene, err = np.loadtxt("Energy_sigma1.0_mu1.0.txt",usecols=(0,1,2), unpack='true')

#plt.errorbar(iblk, ene, yerr=err)

bin, val = np.loadtxt("histo.txt",usecols=(0,1), unpack='true')
plt.figure()
plt.plot(bin, val)

g = f2(0.8,0.6,1.,bin)
plt.figure()
plt.plot(bin,g)



'''
step, ene = np.loadtxt("inst_energy_gas.txt",usecols=(0,1), delimiter='\t', unpack='true')
press = np.loadtxt("inst_pressure_gas.txt",usecols=(1), delimiter='\t', unpack='true')


plt.plot(step, ene, label='Instant Energy')
plt.plot(step, press, label='Instant Pressure')
#plt.plot(mcstep, inst_m, label='Temperature')


plt.legend(loc='upper right')
'''

#x, gofr = np.loadtxt("output.epot.0",usecols=(0,1), delimiter='\t', unpack='true')



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
